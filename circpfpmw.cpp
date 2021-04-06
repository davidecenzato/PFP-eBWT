#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctime>
#include <map>
#include <set>
#include <assert.h>
#include <errno.h>
#include <zlib.h>
#ifdef GZSTREAM
#include <gzstream.h>
#endif
extern "C" {
#include "utils.h"
#include "xerrors.h"
}
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace __gnu_cxx;

// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;
typedef pair <uint32_t,uint32_t> p;

vector<uint64_t> primes{1999999913, 1999999927, 1999999943, 1999999973, 2000000011, 2000000033, 2000000063, 2000000087, 2000000089, 3000000109};
vector<uint64_t> primelw{999999929, 999999937, 1000000021, 1000000033, 1000000087, 1000000093, 1000000097, 1000000103, 1000000207, 1000000297};

// struct containing command line parameters and other globals
struct Args {
   string inputFileName = "";
   size_t w = 10;            // sliding window size and its default
   size_t p = 100;           // modulus for establishing stopping w-tuples
   uint16_t n = 1;
   int th=0;              // number of helper threads
   int verbose=0;         // verbosity level
};

struct word_stats {
  string str;  // parse phrase
  occ_int_t occ;  // no. of phrases
  word_int_t rank=0; // rank of the phrase
};

void print_help(char** argv, Args &args) {
  cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
        << "\t-w W\tsliding window size, def. " << args.w << endl
        << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
        #ifndef NOTHREADS
        << "\t-t M\tnumber of helper threads, def. none " << endl
        #endif
        << "\t-h  \tshow help and exit" << endl;
  exit(1);
}

void parseArgs( int argc, char** argv, Args& arg ) {
    int c;
    extern char *optarg;
    extern int optind;

    puts("==== Command line:");
    for(int i=0;i<argc;i++)
      printf(" %s",argv[i]);
    puts("");
  
    string sarg;
    while ((c = getopt( argc, argv, "p:w:ht:vn:") ) != -1) {
       switch(c) {
         case 'w':
         sarg.assign( optarg );
         arg.w = stoi( sarg ); break;
         case 'n':
         sarg.assign( optarg );
         arg.n = stoi( sarg ); break;
         case 'p':
         sarg.assign( optarg );
         arg.p = stoi( sarg ); break;
         case 't':
         sarg.assign( optarg );
         arg.th = stoi( sarg ); break;
         case 'v':
            arg.verbose++; break;
         case 'h':
            print_help(argv, arg); exit(1);
         case '?':
         cout << "Unknown option. Use -h for help." << endl;
         exit(1);
       }
    }
    // the only input parameter is the file name
    if (argc == optind+1) {
      arg.inputFileName.assign( argv[optind] );
    }
    else {
       cout << "Invalid number of arguments" << endl;
       print_help(argv,arg);
    }
    // check algorithm parameters
    if(arg.w < 2) {
      cout << "Windows size must be at least 2\n";
      exit(1);
    }
    if(arg.p<10) {
      cout << "Modulus must be at leas 10\n";
      exit(1);
    }
    #ifdef NOTHREADS
    if(arg.th!=0) {
      cout << "The NT version cannot use threads\n";
      exit(1);
    }
    #else
    if(arg.th<0) {
      cout << "Number of threads cannot be negative\n";
      exit(1);
    }
    #endif
}

struct KR_window {
  int wsize;
  int current;
  int *window;
  int asize;
  const uint64_t prime;
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot;   // asize^(wsize-1) mod prime

  KR_window(int w, uint64_t prime_): wsize(w), prime(prime_) {
    asize = 256;
    asize_pot = 1;
    //  asize_pot = (asize_pot*asize)% prime; // ugly linear-time power algorithm
    for(int i=1;i<wsize;i++){
        asize_pot = (asize_pot*asize)% prime;}
    // alloc and clear window
    window = new int[wsize];
    reset();
  }

  // init window, hash, and tot_char
  void reset() {
    for(int i=0;i<wsize;i++){ window[i]=0; }
    // init hash value and related values
    hash=tot_char=current=0;
  }

  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
    ++current;
    current = min(wsize,current);
    // complex expression to avoid negative numbers
    hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution
    hash = (asize*hash + c) % prime;      //  add char i
    window[k]=c;
    return hash;
  }
  // debug only
  string get_window() {
    string w = "";
    int k = (tot_char-1) % wsize;
    for(int i=k+1;i<k+1+wsize;i++)
      w.append(1,window[i%wsize]);
    return w;
  }

  void delete_window() {
     delete[] window;
  }
  
};

static void save_update_word(string& w, unsigned int minsize, map<uint64_t,word_stats>& freq, FILE *tmp_parse_file, bool last_word);

#ifndef NOTHREADS
#include "circpfpmw.hpp"
#endif

// compute 64-bit KR hash of a string
// to avoid overflows in 64 bit arithmetic the prime is taken < 2**55
uint64_t kr_hash(string s) {
    uint64_t hash = 0;
    //const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
    const uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for(size_t k=0;k<s.size();k++) {
      int c = (unsigned char) s[k];
      assert(c>=0 && c< 256);
      hash = (256*hash + c) % prime;    //  add char k
    }
    return hash;
}

// save current word in the freq map and update it leaving only the
// last minsize chars which is the overlap with next word
static void save_update_word(string& w, unsigned int minsize,map<uint64_t,word_stats>&  freq, FILE *tmp_parse_file, bool last_word)
{
  assert(w.size() >= minsize);
  if(w.size() <= minsize) return;
  // get the hash value and write it to the temporary parse file
  uint64_t hash = kr_hash(w);
  if(fwrite(&hash,sizeof(hash),1,tmp_parse_file)!=1) die("parse write error");
  if(last_word){
      string lw(minsize,Dollar);
      uint64_t hash = kr_hash(lw);
      if(fwrite(&hash,sizeof(hash),1,tmp_parse_file)!=1) die("parse write error");
  } 

#ifndef NOTHREADS
  xpthread_mutex_lock(&map_mutex,__LINE__,__FILE__);
#endif
  // update frequency table for current hash
  if(freq.find(hash)==freq.end()) {
      freq[hash].occ = 1; // new hash
      freq[hash].str = w;
  }
  else {
      freq[hash].occ += 1; // known hash
      if(freq[hash].occ <=0) {
        cerr << "Emergency exit! Maximum # of occurence of dictionary word (";
        cerr<< MAX_WORD_OCC << ") exceeded\n";
        exit(1);
      }
      if(freq[hash].str != w) {
        cerr << "Emergency exit! Hash collision for strings:\n";
        cerr << freq[hash].str << "\n  vs\n" <<  w << endl;
        exit(1);
      }
  }
#ifndef NOTHREADS
  xpthread_mutex_unlock(&map_mutex,__LINE__,__FILE__);
#endif
  
  // keep only the overlapping part of the window
  w.erase(0,w.size() - minsize);
}


uint16_t compute_windows_n(Args& arg){
 
    string fnam = arg.inputFileName;   
    // main loop on the chars of the input file
    uint8_t c; uint16_t nw = 1; size_t i=0;
    vector<KR_window> windows;
    for(i=0;i<nw;i++){ windows.push_back(KR_window(arg.w,primes[i]));}
    
    gzFile fp;
    kseq_t *seq;
    size_t nseq = 1;
    //p st;
    fp = gzopen(fnam.c_str(), "r");
    seq = kseq_init(fp);
    long int l = kseq_read(seq);
    while (l >= 0) {
        string word(""); bool trg_f = 0; 
        for (i = 0; i < seq->seq.l; i++) {
            c = std::toupper(seq->seq.s[i]);
            word.append(1, c);
            for(uint16_t j = 0; j<nw; ++j){
                uint64_t hash = windows[j].addchar(c);
                if (hash%arg.p==0 && windows[j].current == arg.w) {trg_f=1;break;}
            }
        }
        if(trg_f==0){
            for(i=0;i<arg.w-1;i++){
                c = word[i];
                for(uint16_t j = 0; j<nw; ++j){
                    uint64_t hash = windows[j].addchar(c);
                    if (hash%arg.p==0 && windows[j].current == arg.w) {trg_f=1;break;}
                }
            }
        }
        if(trg_f){ l =  kseq_read(seq); nseq++; }
        else{ windows.push_back(KR_window(arg.w,primes[nw])); ++nw; }
        for(uint16_t j=0;j<nw-1;j++){ windows[j].reset(); }
    }
    for(uint16_t j=0;j<nw;j++){ windows[j].delete_window(); }
    return nw;
}

// circular prefix free parse of fname, w is the window size, p is the modulus
// use a KR-hash as the word ID that is immediately written to the parse file
uint64_t parse_fasta_reads(Args& arg, map<uint64_t,word_stats>& wordFreq)
{
    //open a, possibly compressed, input file
    string fnam = arg.inputFileName;
  
    // open the 1st pass parsing file
    FILE *parse_file = open_aux_file(arg.inputFileName.c_str(),EXTPARS0,"wb");
    // open the words offset file
    FILE *offset_file = open_aux_file(arg.inputFileName.c_str(),EXTOFF0,"wb");
  
    // main loop on the chars of the input file
    uint8_t c; int nw = arg.n;
    vector<KR_window> windows;
    for(uint16_t i=0;i<nw;i++){ windows.push_back(KR_window(arg.w,primes[i]));}
    uint64_t total_char = 0;
    
    gzFile fp;
    kseq_t *seq;
    long int l;
    fp = gzopen(fnam.c_str(), "r");
    seq = kseq_init(fp);
    while ((l =  kseq_read(seq)) >= 0) {
        bool f_trg = 0, win_f = 0;
        uint64_t start_char=0; size_t i=0;
        string first_word(""); string next_word(""); 
        for (i = 0; i < seq->seq.l; i++) {
            c = std::toupper(seq->seq.s[i]);
            if (c <= Dollar) {cerr << "Invalid char found in input file: no additional chars will be read\n"; break;}
            next_word.append(1, c);
            win_f = 0;
            for(uint16_t j = 0; j<nw; ++j){
                uint64_t hash = windows[j].addchar(c);
                if(win_f == 0){
                    if (hash%arg.p==0 && windows[j].current == arg.w) {
                        start_char = i; f_trg = 1;
                        if(fwrite(&start_char,sizeof(start_char),1,offset_file)!=1) die("offset write error");
                        first_word = string(next_word);
                        next_word.erase(0,next_word.size() - arg.w); win_f = 1;
                    }
                }
            }
            if(f_trg==1) break;
        }
        
        for (i=i+1; i< seq->seq.l; i++){
            c = std::toupper(seq->seq.s[i]);
            next_word.append(1, c);
            win_f = 0;
            for(uint16_t j=0;j<nw;j++){
                uint64_t hash = windows[j].addchar(c);
                if(win_f == 0){
                    if (hash%arg.p==0) {
                        save_update_word(next_word,arg.w,wordFreq,parse_file,0); win_f = 1;
                    }
                }
            }
        }
           
        total_char += windows[0].tot_char;
        if(f_trg) { assert(first_word.size() >= arg.w); }
        // check if exist a trigger string in final word
        if(!f_trg) first_word = next_word.substr(0,arg.w - 1);
        for (i = 0; i < arg.w - 1; i++) {
            c = first_word[i];
            next_word.append(1, c);
            win_f = 0;
            for(int j=0;j<nw;j++){
                uint64_t hash = windows[j].addchar(c);
                if(win_f == 0){
                    if(hash%arg.p==0){
                        if(!f_trg) { start_char = windows[0].tot_char; f_trg = 1; 
                                    if(fwrite(&start_char,sizeof(start_char),1,offset_file)!=1) die("offset write error"); 
                                    first_word = string(next_word);
                                    next_word.erase(0,next_word.size() - arg.w); }
                        else{
                            save_update_word(next_word,arg.w,wordFreq,parse_file,0); 
                            }
                        win_f = 1;
                    }
                }
            }
        }
        if(!f_trg) { cerr << "No trigger strings found, use more windows." << endl; exit(1);}
        // join first and last word
        string final_word = next_word + first_word.erase(0,arg.w-1);
        save_update_word(final_word,arg.w,wordFreq,parse_file,1);
        for(int j=0;j<nw;j++){ windows[j].reset(); }
        if (c <= Dollar) break;
    }
    for(uint16_t j=0;j<nw;j++){ windows[j].delete_window(); }
    kseq_destroy(seq);
    gzclose(fp);

    // close input and output files
    if(fclose(parse_file)!=0) die("Error closing parse file");
    if(fclose(offset_file)!=0) die("Error closing offset file");

    return total_char;
}

// given the sorted dictionary and the frequency map write the dictionary and occ files
// also compute the 1-based rank for each hash
void writeDictOcc(Args &arg, map<uint64_t,word_stats> &wfreq, vector<const string *> &sortedDict)
{
  assert(sortedDict.size() == wfreq.size());
  FILE *fdict;
  fdict = open_aux_file(arg.inputFileName.c_str(),EXTDICT,"wb");
  FILE *focc = open_aux_file(arg.inputFileName.c_str(),EXTOCC,"wb");

  word_int_t wrank = 1; // current word rank (1 based)
  for(auto x: sortedDict) {
    const char *word = (*x).data();       // current dictionary word
    int offset=0; size_t len = (*x).size();  // offset and length of word
    assert(len>(size_t)arg.w);
    size_t s = fwrite(word,1,len, fdict);
    if(s!=len) die("Error writing to DICT file");
    if(fputc(EndOfWord,fdict)==EOF) die("Error writing EndOfWord to DICT file");
    uint64_t hash = kr_hash(*x);
    auto& wf = wfreq.at(hash);
    assert(wf.occ>0);
    s = fwrite(&wf.occ,sizeof(wf.occ),1, focc);
    if(s!=1) die("Error writing to OCC file");
    assert(wf.rank==0);
    wf.rank = wrank++;
  }
  if(fputc(EndOfDict,fdict)==EOF) die("Error writing EndOfDict to DICT file");
  if(fclose(focc)!=0) die("Error closing OCC file");
  if(fclose(fdict)!=0) die("Error closing DICT file");
}

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b)
{
  return *a <= *b;
}

void remapParse(Args &arg, map<uint64_t,word_stats> &wfreq, int th)
{
  // open parse files. the old parse can be stored in a single file or in multiple files
  mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), EXTPARS0, th);
  mFile *moff = mopen_aux_file(arg.inputFileName.c_str(), EXTOFF0, th);
  FILE *newp   = open_aux_file(arg.inputFileName.c_str(), EXTPARSE, "wb");
  FILE *newoff  = open_aux_file(arg.inputFileName.c_str(), EXTOFF, "wb");
  FILE *strt   = open_aux_file(arg.inputFileName.c_str(), EXTSTART, "wb");
  FILE *fchar  = open_aux_file(arg.inputFileName.c_str(), EXTFCHAR, "wb");
  
  // recompute occ as an extra check
  vector<occ_int_t> occ(wfreq.size()+1,0); // ranks are zero based
  uint64_t hash, phash, fc;
  uint64_t start = 0, len = 0;
  string separator(arg.w,Dollar);
  uint64_t hash_sep = kr_hash(separator);
  //map <p,uint32_t> startFreq;
  set<p> startChr;

  while(true) {
    size_t s = mfread(&hash,sizeof(hash),1,moldp);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    if(hash != hash_sep){
        len++;
        word_int_t rank = wfreq.at(hash).rank;
        occ[rank]++;
        phash = hash;
        s = fwrite(&rank,sizeof(rank),1,newp);
        if(s!=1) die("Error writing to new parse file");}
    else{
        s = fwrite(&start,sizeof(start),1,strt);
        if(s!=1) die("Error writing to start file");
        start += len;
        len=0;
        s = mfread(&fc,sizeof(fc),1,moff);
        if(s!=1) die("Unexpected offset EOF");
        word_int_t rank = wfreq.at(phash).rank;
        uint64_t len = wfreq.at(phash).str.length();
        uint32_t off = uint32_t(len-fc-1);
        s = fwrite(&off,sizeof(off),1,newoff);
        if(s!=1) die("Error writing to new offset file");
        p st = p(rank,off);
        if(startChr.find(st)==startChr.end()){
          startChr.insert(st);
        }
        //if(startFreq.find(st)==startFreq.end()){
        //    startFreq[st] = 1;
        //    }else{startFreq[st]+=1;}   

        }
    }

  //for (auto& x: startFreq) {
  //    if(fwrite(&x.first,sizeof(x.first),1,fchar)!=1) die("error writing to first char file");
  //}
  for (auto& x: startChr) {
      if(fwrite(&x,sizeof(x),1,fchar)!=1) die("error writing to first char file");
  }
    
  if(fclose(newp)!=0) die("Error closing new parse file");
  if(fclose(fchar)!=0) die("Error closing first char positions file");
  if(mfclose(moldp)!=0) die("Error closing old parse segment");
  if(mfclose(moff)!=0) die("Error closing offset file");
  if(fclose(strt)!=0) die("Error closing starting positions file");
  if(fclose(newoff)!=0) die("Error closing new offsets file");
  // check old and recomputed occ coincide
  for(auto& x : wfreq)
    assert(x.second.occ == occ[x.second.rank]);
}

int main(int argc, char** argv) {
    
    // translate command line parameters
    Args arg;
    parseArgs(argc, argv, arg);
    cout << "File name: " << arg.inputFileName << endl;
    cout << "Windows size: " << arg.w << endl;
    cout << "Stop word modulus: " << arg.p << endl;
    cout << "Windows number: " << arg.n << endl;
    
    if(arg.w<5){ primes = primelw; }
    // measure elapsed wall clock time
    time_t start_main = time(NULL);
    time_t start_wc = start_main;
    // init sorted map counting the number of occurrences of parse phrases
    map <uint64_t,word_stats> wordFreq;
    uint64_t totChar; int nt = 0; // tot characters seen
    
    // ------------ parse input fasta file
    //exit(1);
    try{
        if(arg.th<=1){
            if(arg.n==0){
                uint16_t nw = compute_windows_n(arg);
                cout << "nw " << nw << endl;
                arg.n = nw; }
            totChar = parse_fasta_reads(arg,wordFreq);         
        }
        else
        {
            #ifdef NOTHREADS
            cerr << "Sorry, this is the no-threads executable and you requested " << arg.th << " threads\n";
            exit(1);
            #else
            if(arg.n==0){
                Res res = parallel_parse_fasta(arg, wordFreq, 0);
                arg.n = res.nw;
            }
            Res res = parallel_parse_fasta(arg, wordFreq, 1);
            totChar = res.tot_char; nt = res.us_th;
            #endif
        }
    }
    catch(const std::bad_alloc&) {
    cout << "Out of memory (parsing phase)... emergency exit\n";
    die("bad alloc exception");
    }
    
    uint64_t totDWord = wordFreq.size();
    cout << "Total input symbols: " << totChar << endl;
    cout << "Found " << totDWord << " distinct words" <<endl;
    cout << "Parsing took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    // check # distinct words
    if(totDWord>MAX_DISTINCT_WORDS) {
      cerr << "Emergency exit! The number of distinc words (" << totDWord << ")\n";
      cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS << ")\n";
      exit(1);
    }
    
    // -------------- second pass
    start_wc = time(NULL);
    // create array of dictionary words
    vector<const string *> dictArray;
    dictArray.reserve(totDWord);
    // fill array
    uint64_t sumLen = 0;
    uint64_t totWord = 0;
    for (auto& x: wordFreq) {
      sumLen += x.second.str.size();
      totWord += x.second.occ;
      dictArray.push_back(&x.second.str);
    }
    assert(dictArray.size()==totDWord);
    cout << "Sum of lenghts of dictionary words: " << sumLen << endl;
    cout << "Total number of words: " << totWord << endl;
    // sort dictionary
    sort(dictArray.begin(), dictArray.end(),pstringCompare);
    // write plain dictionary and occ file, also compute rank for each hash
    cout << "Writing plain dictionary and occ file\n";
    writeDictOcc(arg, wordFreq, dictArray);
    dictArray.clear(); // reclaim memory
    cout << "Dictionary construction took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    
    // remap parse file
    start_wc = time(NULL);
    cout << "Generating remapped parse file\n";
    remapParse(arg, wordFreq, nt);
    cout << "Remapping parse file took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    cout << "==== Elapsed time: " << difftime(time(NULL),start_main) << " wall clock seconds\n";
    
    return 0;
}