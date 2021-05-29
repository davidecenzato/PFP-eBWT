/*
 * Multithread Prefix-free parse implementation to compute the circular PFP of sequence collections.
 * 
 * This code is adapted from https://github.com/alshai/Big-BWT/blob/master/newscan.hpp
 *
 */

extern "C" {
#include "xerrors.h"
}
#include <vector>
#include <istream>
pthread_mutex_t map_mutex = PTHREAD_MUTEX_INITIALIZER;

struct Res {
    size_t tot_char = 0;
    int us_th = 0;
};

// struct shared via mt_parse
typedef struct {
  map<uint64_t,word_stats> *wordFreq; // shared dictionary
  Args *arg;       // command line input
  size_t true_start, true_end; // input
  size_t parsed, words;  // output
  FILE *parse, *o;
} mt_data;

// modified from mt_parse to skip newlines and fasta header lines (ie. lines starting with ">")
void *cyclic_mt_parse_fasta(void *dx)
{
  // extract input data
  mt_data *d = (mt_data *) dx;
  Args *arg = d->arg;
  map<uint64_t,word_stats> *wordFreq = d-> wordFreq;

  if(arg->verbose>1)
    printf("Scanning from %ld, size %ld as a FASTA record\n",d->true_start,d->true_end-d->true_start);
  if(d->true_end-d->true_start == 0) return NULL;

  // open input file
  ifstream f(arg->inputFileName);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open file " + arg->inputFileName);
  }
  
  // prepare for parsing
  f.seekg(d->true_start); // move to the beginning of assigned region
  KR_window krw(arg->w);
  uint8_t c, pc = '\n'; string word = ""; string fword = ""; string final_word = "";
  uint64_t start_char = 0;
  bool first_trigger = 0;
  
  // skip the header
  uint64_t current_pos = d->true_start; uint64_t i=0;
  // step when we find the newline
  while((c != '\n')){
      c = f.get();
      current_pos++;
  }
  pc = c;
  //assert(is_valid_base(std::toupper(f.peek())));
  // parse the sequence
  while( (pc != EOF) && current_pos <= d->true_end) {
      c = f.get();
      current_pos++;
      c = std::toupper(c);
      //if(is_valid_base(c)){
      if(c > 64){ // A is 65 in ascii table. 
          if(c<= Dollar || c> 90) die("Invalid char found in input file. Exiting...");
          word.append(1, c);
          if(first_trigger == 0){fword.append(1, c);}
          uint64_t hash = krw.addchar(c);
            if (hash%arg->p==0 && krw.current == arg->w){
                if(first_trigger==0){
                    first_trigger = 1, start_char = i;
                    if(fwrite(&start_char,sizeof(start_char),1,d->o)!=1) die("offset write error");
                    word.erase(0,word.size() - arg->w);
                }
                else{
                    save_update_word(word,arg->w,*wordFreq,d->parse,0);
                    d->words++;
                }
            }
          pc = c; i++;
      }
      else{ 
        if(c == '>' || c == EOF || current_pos >= d->true_end){
            d->parsed += krw.tot_char;
            for (size_t i = 0; i < arg->w - 1; i++) {
                c = fword[i];
                word.append(1, c);
                if(first_trigger == 0){fword.append(1, c);}
                uint64_t hash = krw.addchar(c);
                if(hash%arg->p==0){
                    if(first_trigger==0){
                        first_trigger = 1, start_char = krw.tot_char;
                        if(fwrite(&start_char,sizeof(start_char),1,d->o)!=1) die("offset write error");
                        word.erase(0,word.size() - arg->w);
                    }else{
                        save_update_word(word,arg->w,*wordFreq,d->parse,0);
                        d->words++;}
                }
            }
            if(first_trigger==0) { cerr << "No trigger strings found. Please use '--reads' flag. Exiting..." << endl; exit(1); }
            final_word = word + fword.erase(0,arg->w - 1);
            save_update_word(final_word,arg->w,*wordFreq,d->parse,1);
            d->words++; 
            krw.reset();
            first_trigger = 0; start_char=0; i=0;
            word = fword = final_word = "";
            while(  ((c = f.get()) != EOF) && current_pos <= d->true_end ){
                current_pos++;
                if(c=='\n'){break;}
             }
            pc = c;
        }
      }
    }
  
  f.close(); 
  return NULL;
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
Res parallel_parse_fasta(Args& arg, map<uint64_t,word_stats>& wf)
{
    assert(arg.th>0);
    pthread_t t[arg.th];
    mt_data td[arg.th];
    // scan file for start positions and execute threads
    FILE* fp = fopen(arg.inputFileName.c_str(), "r");
    if (fp == NULL) {
      throw new std::runtime_error("Cannot open input file " +arg.inputFileName);
    }
    fseek(fp, 0L, SEEK_END);
    size_t size = ftell(fp);
    rewind(fp);
    std::vector<size_t> th_sts(arg.th+1);
    th_sts[0] = 1; th_sts[arg.th] = size;
    for (int i = 1; i < arg.th; ++i) {
      th_sts[i] = (size_t) (size / arg.th) * i;
    }
  
    if(arg.verbose) {
    cout << "Thread: " << arg.th << endl;
    cout << "Total size: " << size << endl;
    cout << "------------------------" << endl; }
  
    uint8_t c = fgetc(fp); size_t hp = 1; bool sf = 0;
    size_t tstart = 1, tend = 1;
    size_t nt = 1;
    // this loop scans the Fasta file in order to properly divide it up
    // for the threads, so that they don't accidently start in a ">" header.
    // As soon as a proper start and end position has been found, execute the thread 
    for(int i=1;i<th_sts.size()-1;i++){
        
        size_t start = th_sts[i];
        size_t end = th_sts[i+1];
        fseek(fp, start, SEEK_SET);
        //cout << "scanning " << start << " - " << end << endl;
        for(size_t j=start; j<end; ++j){
            c = fgetc(fp);
            if(c == '>'){hp=j;sf=1;break;}
        }
        
        if(sf==1){
            sf=0;
            //tend = hp-1;
            tend = hp;
            // prepare and execute thread j-1
            td[nt-1].wordFreq = &wf;
            td[nt-1].arg = &arg;
            td[nt-1].true_start = tstart;
            td[nt-1].true_end = tend;
            td[nt-1].words = 0;
            td[nt-1].parsed = 0;
            assert(td[nt-1].true_end <=size);
            td[nt-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,nt-1,"wb");
            td[nt-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,nt-1,"wb");
            xpthread_create(&t[nt-1],NULL,&cyclic_mt_parse_fasta,&td[nt-1],__LINE__,__FILE__);
            nt++;
            tstart = hp+1;
        }
    }
    
    tend = size;
    td[nt-1].wordFreq = &wf;
    td[nt-1].arg = &arg;
    td[nt-1].true_start = tstart;
    td[nt-1].true_end = tend;
    td[nt-1].words = 0;
    td[nt-1].parsed = 0;
    assert(td[nt-1].true_end <=size);
    td[nt-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,nt-1,"wb");
    td[nt-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,nt-1,"wb");
    xpthread_create(&t[nt-1],NULL,&cyclic_mt_parse_fasta,&td[nt-1],__LINE__,__FILE__);
    
    // wait for the threads to finish (in order) and close output files
    size_t tot_char=0;
    for(int i=0;i<nt;i++) {
      xpthread_join(t[i],NULL,__LINE__,__FILE__); 
      if(arg.verbose) {
      cout << "s:" << td[i].true_start << "  e:" << td[i].true_end << "  pa:";
      }
    // close thread-specific output files
    fclose(td[i].parse);
    fclose(td[i].o);
    if(td[i].words>0) {
      // extra check
      assert(td[i].parsed>arg.w);
      tot_char += td[i].parsed;
    }
    else assert(i>0); // the first thread must produce some words 
    }
    fclose(fp);
    Res res; res.tot_char = tot_char; res.us_th = nt;
    return res;
}
