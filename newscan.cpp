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
//#include "parse.hpp"

using namespace std;
using namespace __gnu_cxx;

// =============== algorithm limits ===================
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;

// struct containing command line parameters and other globals
struct Args {
   string inputFileName = "";
   size_t w = 10;            // sliding window size and its default
   size_t p = 100;           // modulus for establishing stopping w-tuples
   bool SAinfo = false;   // compute SA information
   bool is_fasta = false;  // read a fasta file
   bool compress = false; // parsing called in compress mode
   bool p_val = false; // validate the parse
   int th=0;              // number of helper threads
   int verbose=0;         // verbosity level
};

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
  string str;
  occ_int_t occ;
  word_int_t rank=0;
};

bool is_gzipped(std::string fname) {
    FILE* fp = fopen(fname.c_str(), "rb");
    int byte1 = 0, byte2 = 0;
    fread(&byte1, sizeof(char), 1, fp);
    fread(&byte2, sizeof(char), 1, fp);
    fclose(fp);
    return (byte1 == 0x1f && byte2 == 0x8b);
}

void print_help(char** argv, Args &args) {
  cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
        << "\t-w W\tsliding window size, def. " << args.w << endl
        << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
        #ifndef NOTHREADS
        << "\t-t M\tnumber of helper threads, def. none " << endl
        #endif
        << "\t-h  \tshow help and exit" << endl
        << "\t-s  \tcompute suffix array info" << endl;
  #ifdef GZSTREAM
  cout << "If the input file is gzipped it is automatically extracted\n";
  #endif
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
    while ((c = getopt( argc, argv, "p:w:fasht:v") ) != -1) {
       switch(c) {
         case 's':
         arg.SAinfo = true; break;
         case 'c':
         arg.compress = true; break;
         case 'w':
         sarg.assign( optarg );
         arg.w = stoi( sarg ); break;
         case 'p':
         sarg.assign( optarg );
         arg.p = stoi( sarg ); break;
         case 'f':
         arg.is_fasta = true; break;
         case 'a':
         arg.p_val = true; break;
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
    if(arg.w <4) {
      cout << "Windows size must be at least 4\n";
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
  const uint64_t prime = 1999999973;
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot;   // asize^(wsize-1) mod prime

  KR_window(int w): wsize(w) {
    asize = 256;
    asize_pot = 1;
    for(int i=1;i<wsize;i++)
      asize_pot = (asize_pot*asize)% prime; // ugly linear-time power algorithm
    // alloc and clear window
    window = new int[wsize];
    reset();
  }

  // init window, hash, and tot_char
  void reset() {
    for(int i=0;i<wsize;i++) window[i]=0;
    // init hash value and related values
    hash=tot_char=current=0;
  }

  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
    current++;
    current = min(wsize,current);
    // complex expression to avoid negative numbers
    hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution
    hash = (asize*hash + c) % prime;      //  add char i
    window[k]=c;
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
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

  ~KR_window() {
    delete[] window;
  }

};

static void save_update_word(string& w, unsigned int minsize, map<uint64_t,word_stats>& freq, FILE *tmp_parse_file, FILE *last, FILE *sa, FILE *o, uint64_t &pos, uint64_t &offset, bool last_word);

#ifndef NOTHREADS
#include "newscan.hpp"
#endif

// compute 64-bit KR hash of a string
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
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
static void save_update_word(string& w, unsigned int minsize,map<uint64_t,word_stats>&  freq, FILE *tmp_parse_file, FILE *last, FILE *sa, FILE *off, uint64_t &pos, uint64_t &offset,  bool last_word)
{
  assert(pos==0 || w.size() > minsize);
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
  
  // output char w+1 from the end
  if(fputc(w[w.size()- minsize-1],last)==EOF) die("Error writing to .last file");
  // compute ending position +1 of current word and write it to sa file
  // pos is the ending position+1 of the previous word and is updated here
  //if(pos==0) pos = w.size(); // -1 is for the initial $ of the first word
  if(!last_word){pos += w.size() - minsize;}
  if(sa) if(fwrite(&pos,IBYTES,1,sa)!=1) die("Error writing to sa info file");
  
  if(fwrite(&offset,sizeof(offset),1,off)!=1) die("Error writing to offset file");
  offset += w.size() - minsize;
  if(last_word){
      offset = 0;
      if(fwrite(&offset,sizeof(offset),1,off)!=1) die("Error writing last element to offset file");
  }
  
  // keep only the overlapping part of the window
  w.erase(0,w.size() - minsize);
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is immediately written to the parse file
uint64_t process_file(Args& arg, map<uint64_t,word_stats>& wordFreq)
{
  //open a, possibly compressed, input file
  string fnam = arg.inputFileName;
  
  // open the 1st pass parsing file
  FILE *g = open_aux_file(arg.inputFileName.c_str(),EXTPARS0,"wb");
  // open output file containing the char at position -(w+1) of each word
  FILE *last_file = open_aux_file(arg.inputFileName.c_str(),EXTLST,"wb");
  // if requested open file containing the ending position+1 of each word
  // open the words offset file
  FILE *offset_file = open_aux_file(arg.inputFileName.c_str(),EXTOFF,"wb");
  FILE *sa_file = NULL;
  if(arg.SAinfo)
    sa_file = open_aux_file(arg.inputFileName.c_str(),EXTSAI,"wb");
  
  // main loop on the chars of the input file
  int c;
  int no_seq = 0;
  KR_window krw(arg.w);
  uint64_t total_char = 0;
  std::string line;
  
  if (arg.is_fasta) {
    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(fnam.c_str(), "r");
    seq = kseq_init(fp);
    while ((l =  kseq_read(seq)) >= 0) {
        no_seq++;
        uint64_t pos = 0, offset = 1, first_offset = 1;
        uint64_t first_pos = 0;
        assert(IBYTES<=sizeof(pos));
        string first_word("");
        string next_word("");
        bool first_trigger = 0;
        for (size_t i = 0; i < seq->seq.l; i++) {
            c = std::toupper(seq->seq.s[i]);
            if (c <= Dollar) {cerr << "Invalid char found in input file: no additional chars will be read\n"; break;}
            next_word.append(1, c);
            if(first_trigger == 0){first_word.append(1, c);}
            uint64_t hash = krw.addchar(c);
            if (hash%arg.p==0 && krw.current == arg.w) {
                if(first_trigger==0){
                    first_trigger = 1;
                    first_pos = pos = first_word.size() - arg.w;
                    first_offset = offset = first_pos + 1;
                    next_word.erase(0,next_word.size() - arg.w);
                }
                else{
                    save_update_word(next_word,arg.w,wordFreq,g,last_file,sa_file,offset_file,pos,offset,0);
                }
            }
        } 
        total_char += krw.tot_char;
        assert(first_word.size() >= arg.w);
        pos = krw.tot_char - arg.w;
        // check if exist a trigger string in final word
        for (size_t i = 0; i < arg.w - 1; i++) {
            c = first_word[i];
            next_word.append(1, c);
            pos++;
            uint64_t hash = krw.addchar(c);
            if(hash%arg.p==0){
                save_update_word(next_word,arg.w,wordFreq,g,last_file,sa_file,offset_file,pos,offset,0);
            }
        }
        // join first and last word
        string final_word = next_word + first_word.erase(0,arg.w-1);
        pos = first_pos;
        offset = first_offset;
        save_update_word(final_word,arg.w,wordFreq,g,last_file,sa_file,offset_file,pos,offset,1);
        krw.reset();
        if (c <= Dollar) break;
    }
    kseq_destroy(seq);
    gzclose(fp);
  }
  else {
      cerr << "ERROR: Only fasta files are supported!" << endl; 
      exit(1);
  }

  // close input and output files
  if(sa_file) if(fclose(sa_file)!=0) die("Error closing SA file");
  if(fclose(last_file)!=0) die("Error closing last file");
  if(fclose(g)!=0) die("Error closing parse file");
  if(fclose(offset_file)!=0) die("Error closing offset file");

  return total_char;
}

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b)
{
  return *a <= *b;
}

// given the sorted dictionary and the frequency map write the dictionary and occ files
// also compute the 1-based rank for each hash
void writeDictOcc(Args &arg, map<uint64_t,word_stats> &wfreq, vector<const string *> &sortedDict)
{
  assert(sortedDict.size() == wfreq.size());
  FILE *fdict;
  // open dictionary and occ files
  if(arg.compress)
    fdict = open_aux_file(arg.inputFileName.c_str(),EXTDICZ,"wb");
  else
    fdict = open_aux_file(arg.inputFileName.c_str(),EXTDICT,"wb");
  FILE *focc = open_aux_file(arg.inputFileName.c_str(),EXTOCC,"wb");

  word_int_t wrank = 1; // current word rank (1 based)
  for(auto x: sortedDict) {
    const char *word = (*x).data();       // current dictionary word
    int offset=0; size_t len = (*x).size();  // offset and length of word
    assert(len>(size_t)arg.w);
    if(arg.compress) {  // if we are compressing remove overlapping and extraneous chars
      len -= arg.w;     // remove the last w chars
    }
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

void remapParse(Args &arg, map<uint64_t,word_stats> &wfreq)
{
  // open parse files. the old parse can be stored in a single file or in multiple files
  mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), EXTPARS0, arg.th);
  FILE *newp = open_aux_file(arg.inputFileName.c_str(), EXTPARSE, "wb");
  FILE *lenf = open_aux_file(arg.inputFileName.c_str(), EXTLEN, "wb");
  FILE *strt = open_aux_file(arg.inputFileName.c_str(), EXTSTART, "wb");
  
  // recompute occ as an extra check
  vector<occ_int_t> occ(wfreq.size()+1,0); // ranks are zero based
  uint64_t hash;
  uint32_t len = 0;
  uint64_t start = 0;
  string separator(arg.w,Dollar);
  uint64_t hash_sep = kr_hash(separator);
  //size_t s = fwrite(&start,sizeof(start),1,strt);
  //if(s!=1) die("Error writing to new parse file");
  while(true) {
    size_t s = mfread(&hash,sizeof(hash),1,moldp);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    if(hash != hash_sep){
        len++;
        word_int_t rank = wfreq.at(hash).rank;
        occ[rank]++;
        s = fwrite(&rank,sizeof(rank),1,newp);
        if(s!=1) die("Error writing to new parse file");}
    else{
        s = fwrite(&start,sizeof(start),1,strt);
        if(s!=1) die("Error writing to start file");
        start += len;
        //word_int_t rank = 0;
        //s = fwrite(&rank,sizeof(rank),1,newp);
        //if(s!=1) die("Error writing to new parse file");
        s = fwrite(&len,sizeof(len),1,lenf);
        if(s!=1) die("Error writing to lengths file");
        len=0;}
    }
  if(fclose(newp)!=0) die("Error closing new parse file");
  if(mfclose(moldp)!=0) die("Error closing old parse segment");
  if(fclose(lenf)!=0) die("Error closing lengths file");
  if(fclose(strt)!=0) die("Error closing starting positions file");
  // check old and recomputed occ coincide
  for(auto& x : wfreq)
    assert(x.second.occ == occ[x.second.rank]);
}

void remapOffset(Args &arg)
{
    FILE *newoff = open_aux_file(arg.inputFileName.c_str(), EXTOFF, "wb");
    uint64_t offset = 0;
    for(int i=0; i<arg.th; i++){
        string filename = arg.inputFileName + string(".") + to_string(i) + string(".offset");
        //cout << "filename: " << filename << endl;
        FILE *off = fopen(filename.c_str(), "r");
        //int fn = fileno(off);
        //cout << fn << endl;
        while(true){
            size_t s = fread(&offset,sizeof(offset),1,off);
            //cout << "off: " << offset << endl;
            if(s==0) break;
            if(s!=1) die("Unexpected parse EOF");
            s = fwrite(&offset,sizeof(offset),1,newoff);
            if(s!=1) die("Error writing to offset file");
        }
        if(fclose(off)!=0) die("Error closing old offset segment");
    }
    if(fclose(newoff)!=0) die("Error closing new offset file");
}


int main(int argc, char** argv) {
    
    // translate command line parameters
    Args arg;
    parseArgs(argc, argv, arg);
    cout << "File name: " << arg.inputFileName << endl;
    cout << "Windows size: " << arg.w << endl;
    cout << "Stop word modulus: " << arg.p << endl;
    
    // measure elapsed wall clock time
    time_t start_main = time(NULL);
    time_t start_wc = start_main;
    // init sorted map counting the number of occurrences of each word
    map <uint64_t,word_stats> wordFreq;
    uint64_t totChar;
  
    // ------------ parsing input file
    // force threads=0 if gzipped
    if (is_gzipped(arg.inputFileName)) {
        cerr << "input is gzipped, forcing single thread!" << endl;
        arg.th = 0;
    }
    try {
      if(arg.th==0){
          totChar = process_file(arg,wordFreq);
        }
      else {
        #ifdef NOTHREADS
        cerr << "Sorry, this is the no-threads executable and you requested " << arg.th << " threads\n";
        exit(EXIT_FAILURE);
        #else
        if (arg.is_fasta){ 
            totChar = cyclic_mt_process_file_fasta(arg, wordFreq);
        }
        else{ cerr << "Error: Only fasta files are supported." << endl; //totChar = mt_process_file(arg,wordFreq);
              exit(1);} 
        #endif
        }
      }
      catch(const std::bad_alloc&) {
      cout << "Out of memory (parsing phase)... emergency exit\n";
      die("bad alloc exception");
      }
    
    // first report

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
    remapParse(arg, wordFreq);
    cout << "Generating remapped offset file\n";
    remapOffset(arg);
    cout << "Remapping parse file took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
    cout << "==== Elapsed time: " << difftime(time(NULL),start_main) << " wall clock seconds\n";
    
    return 0;
}