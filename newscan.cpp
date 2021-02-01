/* ******************************************************************************
 * newscan.cpp
 *
 * parsing algorithm for bwt construction of repetitive sequences based
 * on prefix free parsing. See:
 *   Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini
 *   Prefix-Free Parsing for Building Big BWTs
 *   [Proc. WABI '18](http://drops.dagstuhl.de/opus/volltexte/2018/9304/)
 *
 * Usage:
 *   newscan.x wsize modulus file
 *
 * Accepts any kind of file that does not contain the chars 0x0, 0x1, 0x2
 * which are used internally. If input file is gzipped use cnewscan.x which
 * automatically extracts the content
 *
 * The parameters wsize and modulus are used to define the prefix free parsing
 * using KR-fingerprints (see paper)
 *
 * The algorithm computes the prefix free parsing of
 *     T = (0x2)file_content(0x2)^wsize
 * in a dictionary of words D and a parsing P of T in terms of the
 * dictionary words. Note that the words in the parsing overlap by wsize.
 * Let d denote the number of words in D and p the number of phrases in
 * the parsing P
 *
 * newscan outputs the following files:
 *
 * file.dict
 * containing the dictionary words in lexicographic order with a 0x1 at the end of
 * each word and a 0x0 at the end of the file. Size: |D| + d + 1 where
 * |D| is the sum of the word lengths
 *
 * file.occ
 * the number of occurrences of each word in lexicographic order.
 * We assume the number of occurrences of each word is at most 2^32-1
 * so the size is 4d bytes
 *
 * file.parse
 * containing the parse P with each word identified with its 1-based lexicographic
 * rank (ie its position in D). We assume the number of distinct words
 * is at most 2^32-1, so the size is 4p bytes
 *
 * file.last
 * contaning the charater in positon w+1 from the end for each dictionary word
 * Size: d
 *
 * file.sai (if option -s is given on the command line)
 * containing the ending position +1 of each dictionary word in the original
 * text written using IBYTES bytes for each entry
 * Size: d*IBYTES
 *
 * The output of newscan must be processed by bwtparse, which invoked as
 *
 *    bwtparse file
 *
 * computes the BWT of file.parse and produces file.ilist of size 4p+4 bytes
 * contaning, for each dictionary word in lexicographic order, the list
 * of BWT positions where that word appears (ie i\in ilist(w) <=> BWT[i]=w).
 * There is also an entry for the EOF word which is not in the dictionary
 * but is assumed to be the smallest word.
 *
 * In addition, bwtparse permutes file.last according to
 * the BWT permutation and generates file.bwlast such that file.bwlast[i]
 * is the char from P[SA[i]-2] (if SA[i]==0 , BWT[i]=0 and file.bwlast[i]=0,
 * if SA[i]==1, BWT[i]=P[0] and file.bwlast[i] is taken from P[n-1], the last
 * word in the parsing).
 *
 * If the option -s is given to bwtparse, it permutes file.sai according
 * to the BWT permutation and generate file.bwsai using again IBYTES
 * per entry.  file.bwsai[i] is the ending position+1 ofBWT[i] in the
 * original text
 *
 * The output of bwtparse (the files .ilist .bwlast) together with the
 * dictionary itself (file .dict) and the number of occurrences
 * of each word (file .occ) are used to compute the final BWT by the
 * pfbwt algorithm.
 *
 */
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
#ifdef GZSTREAM
#include <gzstream.h>
#endif
extern "C" {
#include "utils.h"
}
#include <zlib.h>
#include <stdio.h>
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

uint8_t asc2dnacat[] = {
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0,
           /*                                        - */
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
           /*    A  B  C  D        G  H        K     M  N */
    /*  80 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
           /*       R  S  T     V  W  X  Y */
    /*  96 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
           /*    a  b  c  d        g  h        k     m  n */
    /* 112 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
           /*       r  s  t     v  w  x  y */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
  string str;
  occ_int_t occ;
  word_int_t rank=0;
};

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
   string inputFileName = "";
   size_t w = 10;            // sliding window size and its default
   size_t p = 100;           // modulus for establishing stopping w-tuples
   bool SAinfo = false;   // compute SA information
   bool is_fasta = false;  // read a fasta file
   bool compress = false; // parsing called in compress mode
   int th=0;              // number of helper threads
   int verbose=0;         // verbosity level
};


// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
struct KR_window {
  int wsize;
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
    hash=tot_char=0;
  }

  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
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
// -----------------------------------------------------------

static void save_update_word(string& w, unsigned int minsize,map<uint64_t,word_stats>&  freq, FILE *tmp_parse_file, FILE *last, FILE *sa, uint64_t &pos);

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
static void save_update_word(string& w, unsigned int minsize,map<uint64_t,word_stats>&  freq, FILE *tmp_parse_file, FILE *last, FILE *sa, uint64_t &pos)
{
  assert(pos==0 || w.size() > minsize);
  if(w.size() <= minsize) return;
  // get the hash value and write it to the temporary parse file
  uint64_t hash = kr_hash(w);
  if(fwrite(&hash,sizeof(hash),1,tmp_parse_file)!=1) die("parse write error");

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
  if(pos==0) pos = w.size()-1; // -1 is for the initial $ of the first word
  else pos += w.size() -minsize;
  if(sa) if(fwrite(&pos,IBYTES,1,sa)!=1) die("Error writing to sa info file");
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
  FILE *sa_file = NULL;
  if(arg.SAinfo)
    sa_file = open_aux_file(arg.inputFileName.c_str(),EXTSAI,"wb");

  // main loop on the chars of the input file
  int c;
  uint64_t pos = 0; // ending position +1 of previous word in the original text, used for computing sa_info
  assert(IBYTES<=sizeof(pos)); // IBYTES bytes of pos are written to the sa info file
  // init first word in the parsing with a Dollar char
  string word("");
  word.append(1,Dollar);
  KR_window krw(arg.w);
  std::string line;
  if (arg.is_fasta) {
      gzFile fp;
      kseq_t *seq;
      int l;
      fp = gzopen(fnam.c_str(), "r");
      seq = kseq_init(fp);
      while ((l =  kseq_read(seq)) >= 0) {
          for (size_t i = 0; i < seq->seq.l; i++) {
              c = std::toupper(seq->seq.s[i]);
              if (c <= Dollar) {cerr << "Invalid char found in input file: no additional chars will be read\n"; break;}
              word.append(1, c);
              uint64_t hash = krw.addchar(c);
              if (hash%arg.p==0) {
                  save_update_word(word,arg.w,wordFreq,g,last_file,sa_file,pos);
              }
          }
          if (c <= Dollar) break;
      }
      kseq_destroy(seq);
      gzclose(fp);
  }
  else {
      #ifdef GZSTREAM
      igzstream f(fnam.c_str());
      #else
      ifstream f(fnam);
      #endif
      if(!f.rdbuf()->is_open()) {// is_open does not work on igzstreams
        perror(__func__);
        throw new std::runtime_error("Cannot open input file " + fnam);
      }
      while( (c = f.get()) != EOF ) {
        if(c<=Dollar) {cerr << "Invalid char found in input file: no additional chars will be read\n"; break;}
        word.append(1,c);
        uint64_t hash = krw.addchar(c);
        if(hash%arg.p==0) {
          // end of word, save it and write its full hash to the output file
          // cerr << "~"<< c << "~ " << hash << " ~~ <" << word << "> ~~ <" << krw.get_window() << ">" <<  endl;
          save_update_word(word,arg.w,wordFreq,g,last_file,sa_file,pos);
        }
      }
      f.close();
  }
  // virtually add w null chars at the end of the file and add the last word in the dict
  word.append(arg.w,Dollar);
  save_update_word(word,arg.w,wordFreq,g,last_file,sa_file,pos);
  // close input and output files
  if(sa_file) if(fclose(sa_file)!=0) die("Error closing SA file");
  if(fclose(last_file)!=0) die("Error closing last file");
  if(fclose(g)!=0) die("Error closing parse file");
  if(pos!=krw.tot_char+arg.w) cerr << "Pos: " << pos << " tot " << krw.tot_char << endl;
  return krw.tot_char;
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
      if(word[0]==Dollar) {offset=1; len -= 1;} // remove the very first Dollar
    }
    size_t s = fwrite(word+offset,1,len, fdict);
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

  // recompute occ as an extra check
  vector<occ_int_t> occ(wfreq.size()+1,0); // ranks are zero based
  uint64_t hash;
  while(true) {
    size_t s = mfread(&hash,sizeof(hash),1,moldp);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    word_int_t rank = wfreq.at(hash).rank;
    occ[rank]++;
    s = fwrite(&rank,sizeof(rank),1,newp);
    if(s!=1) die("Error writing to new parse file");
  }
  if(fclose(newp)!=0) die("Error closing new parse file");
  if(mfclose(moldp)!=0) die("Error closing old parse segment");
  // check old and recomputed occ coincide
  for(auto& x : wfreq)
    assert(x.second.occ == occ[x.second.rank]);
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
   while ((c = getopt( argc, argv, "p:w:fsht:v") ) != -1) {
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


bool is_gzipped(std::string fname) {
    FILE* fp = fopen(fname.c_str(), "rb");
    int byte1 = 0, byte2 = 0;
    fread(&byte1, sizeof(char), 1, fp);
    fread(&byte2, sizeof(char), 1, fp);
    fclose(fp);
    return (byte1 == 0x1f && byte2 == 0x8b);
}


int main(int argc, char** argv)
{
  // translate command line parameters
  Args arg;
  parseArgs(argc, argv, arg);
  cout << "Windows size: " << arg.w << endl;
  cout << "Stop word modulus: " << arg.p << endl;

  // measure elapsed wall clock time
  time_t start_main = time(NULL);
  time_t start_wc = start_main;
  // init sorted map counting the number of occurrences of each word
  map<uint64_t,word_stats> wordFreq;
  uint64_t totChar;

  // ------------ parsing input file
  // force threads=0 if gzipped
  if (is_gzipped(arg.inputFileName)) {
      cerr << "input is gzipped, forcing single thread!" << endl;
      arg.th = 0;
  }
  try {
      if(arg.th==0)
        totChar = process_file(arg,wordFreq);
      else {
        #ifdef NOTHREADS
        cerr << "Sorry, this is the no-threads executable and you requested " << arg.th << " threads\n";
        exit(EXIT_FAILURE);
        #else
        if (arg.is_fasta) totChar = mt_process_file_fasta(arg, wordFreq);
        else totChar = mt_process_file(arg,wordFreq);
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
  cout << "Remapping parse file took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";
  cout << "==== Elapsed time: " << difftime(time(NULL),start_main) << " wall clock seconds\n";
  return 0;
}

