/*
 * Multithread PFP parse implementation to compute the circular PFP of a collection of sequences.
 * This version can use multiple primes to find trigger strings in the sequences.
 * 
 * This code is adapted from https://github.com/alshai/Big-BWT/blob/master/newscan.cpp
 *
 */

#ifndef CIRCPFPMW_H
#define CIRCPFPMW_H

extern "C" {
#include "xerrors.h"
}
#include <vector>
#include <istream>
pthread_mutex_t map_mutex = PTHREAD_MUTEX_INITIALIZER;

struct Res {
    size_t tot_char = 0;
    int us_th = 0;
    uint16_t nw = 0;
};

// struct shared via mt_parse
typedef struct {
  map<uint64_t,word_stats> *wordFreq; // shared dictionary
  Args *arg;       // command line input
  size_t true_start, true_end; // input
  size_t parsed, words;  // output
  FILE *parse, *o;
  uint16_t nw = 0;
  vector<bool> used;
  vector<uint64_t> to_rmv;
  string filtered;
  int th_id;
} mt_data;

// skip filtered sequences
void skip_seq(uint64_t &nseq, uint64_t &current_pos, size_t &next_tr, size_t &tr_i,
              string &filt_seq, uint8_t &c, ifstream &f, mt_data &d){
    
    while(nseq==next_tr){
        while(  ((c = f.get()) != '>') && current_pos < d.true_end){
            ++current_pos; filt_seq.append(1,std::toupper(c));
        }
        ++nseq; ++current_pos;
        if(tr_i+1 < to_remove_mt[d.th_id].size()){ ++tr_i; next_tr = to_remove_mt[d.th_id][tr_i]; }
        else{next_tr = 0;}
    }
    while( ((c = f.get()) != EOF) && current_pos < d.true_end ){
        current_pos++;
        if(c=='\n'){break;}
    }
}

// modified from mt_parse to skip newlines and fasta header lines (ie. lines starting with ">")
void *cyclic_mt_parse_fasta(void *dx)
{
  // extract input data
  mt_data *d = (mt_data *) dx;
  Args *arg = d->arg;
  map<uint64_t,word_stats> *wordFreq = d-> wordFreq;
  string filt_seq = "";

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
  vector<KR_window> windows;
  for(uint16_t i=0;i<arg->n;i++){ windows.push_back(KR_window(arg->w,primes[i]));}
  uint8_t c; string word = ""; string fword = ""; string final_word = "";
  uint64_t start_char = 0, nseq = 1;
  size_t next_tr = 0, tr_i = 0; if(to_remove_mt[d->th_id].size()>0){ next_tr = to_remove_mt[d->th_id][tr_i]; }
  bool first_trigger = 0;
  
  // skip the header
  uint64_t current_pos = d->true_start; uint64_t sp=0; size_t i,j;
  // step when we find the correct starting point
  skip_seq(nseq, current_pos, next_tr, tr_i, filt_seq, c, f, *d);
  // parse the sequence
  while( current_pos < d->true_end ) {
      c = f.get();
      ++current_pos;
      c = std::toupper(c);
      if(c > 64){ // A is 65 in ascii table. 
          if(c<= 64 || c> 90) { die("Invalid char found in input file. Exiting..."); }
          word.append(1, c);
          if(first_trigger == 0){ fword.append(1, c); }
          bool wind_f = 0;
          for(j=0; j<arg->n; ++j){
              uint64_t hash = windows[j].addchar(c);
              if (hash%arg->p==0 && !wind_f && windows[j].current == arg->w){
                  if(first_trigger==0){
                      first_trigger = 1, start_char = sp;
                      if(fwrite(&start_char,sizeof(start_char),1,d->o)!=1) die("offset write error");
                      word.erase(0,word.size() - arg->w); wind_f=1;
                  }
                  else{
                      save_update_word(word,arg->w,*wordFreq,d->parse,0);
                      d->words++; wind_f=1;
                    }
                }
            }
            ++sp;
      }
      else{ 
          if(c == '>' || current_pos >= d->true_end-1){
            d->parsed += windows[0].tot_char;
            for (i = 0; i < arg->w - 1; ++i) {
                c = fword[i];
                word.append(1, c);
                if(first_trigger == 0){fword.append(1, c);}
                bool wind_f = 0;
                for(j=0; j<arg->n; ++j){
                    uint64_t hash = windows[j].addchar(c);
                    if(hash%arg->p==0 && !wind_f){
                        if(first_trigger==0){
                            first_trigger = 1, start_char = windows[0].tot_char;
                            if(fwrite(&start_char,sizeof(start_char),1,d->o)!=1) die("offset write error");
                            word.erase(0,word.size() - arg->w); wind_f=1;
                        }else{
                            save_update_word(word,arg->w,*wordFreq,d->parse,0);
                            d->words++; wind_f=1; }
                    }
                }
            }
            if(first_trigger==0) { cerr << "No trigger strings found, use more windows." << endl; exit(1); }
            final_word = word + fword.erase(0,arg->w - 1);
            save_update_word(final_word,arg->w,*wordFreq,d->parse,1);
            d->words++; 
            for(j=0; j<arg->n; ++j){ windows[j].reset(); }
            first_trigger = 0; start_char=0; sp=0;
            word = fword = final_word = ""; ++nseq;
            
            skip_seq(nseq, current_pos, next_tr, tr_i, filt_seq, c, f, *d);
        }
      }
    }
  
  for(j=0; j<arg->n; ++j){ windows[j].delete_window(); }
  d->filtered = filt_seq;
  f.close(); 
  return NULL;
}

// modified from mt_parse to skip newlines and fasta header lines (ie. lines starting with ">")
void *parallel_windows_n(void *dx)
{
  // extract input data
  mt_data *d = (mt_data *) dx;
  Args *arg = d->arg;
  
  // open input file
  ifstream f(arg->inputFileName);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open file " + arg->inputFileName);
  }
  
  // prepare for parsing
  uint16_t nw = 1; size_t beginning = d->true_start; 
  f.seekg(d->true_start); // move to the beginning of assigned region
  vector<KR_window> windows; vector<bool> used(primes.size(),0); used[0]=1;
  for(uint16_t i=0;i<nw;i++){ windows.push_back(KR_window(arg->w,primes[i]));}
  uint8_t c; string word = ""; bool trg_f = 0, new_w = 0; size_t nseq = 1;
  vector<uint64_t> to_remove; uint16_t cw = 0;
  
  // skip the header
  uint64_t current_pos = d->true_start; uint64_t i=0,j=0;
  // step when we find the newline
  while((c != '\n')){
      c = f.get();
      ++current_pos;
  }
  beginning = current_pos;
  // parse the sequence
  while( current_pos < d->true_end ) {
      c = f.get(); c = std::toupper(c);
      ++current_pos;
      if(c > 64){ // A is 65 in ascii table. 
          if(c<= 64 || c> 90) die("Invalid char found in input file. Exiting...");
          word.append(1, c);
          if(!trg_f){
            for(j=0; j<nw; ++j){
                uint64_t hash = windows[j].addchar(c);
                if (hash%arg->p==0 && windows[j].current == arg->w){ trg_f=1; break; }
              }
          }
      }
      else{ 
        if(c == '>' || current_pos >= d->true_end-1){
            if(!trg_f){
                for (i = 0; i < arg->w - 1; ++i) {
                    c = word[i];
                    bool wind_f = 0;
                    for(j=0; j<nw; ++j){
                        uint64_t hash = windows[j].addchar(c);
                        if(hash%arg->p==0 && !wind_f){ trg_f=1; break; }
                    }
                    if(trg_f){ break; }
                }
            }
            
            bool rem = 0;
            if(!trg_f){
                // if we reached the maximum window number and we are not testing a new window
                if(nw == arg->n && !new_w){
                    trg_f = 1;  to_remove.push_back(nseq); rem = 1;
                }else{
                    if(cw+1 == arg->n){
                        trg_f = 1;  to_remove.push_back(nseq); rem = 1;
                        windows[windows.size()-1].delete_window(); windows.pop_back(); used[cw] = 0; --nw; 
                    }
                }
            }

            if(arg->c && !rem && !new_w){
                vector<uint16_t> lps(word.size(),0); LPS(word,word.size(),lps);
                if(word.size() > 0 && word.size()%(word.size()-lps[lps.size()-1]) == 0 && lps[lps.size()-1] > 0){ trg_f = 1; to_remove.push_back(nseq); }
            }
            
            word.erase(0,word.size());
            
            if(trg_f){
                new_w = 0;  cw = 0; ++nseq;
                while(  ((c = f.get()) != EOF) && current_pos < d->true_end ){
                    ++current_pos;
                    if(c=='\n'){break;}
                }
                beginning = current_pos; trg_f = 0;
            }
            else{ 
                current_pos = beginning; f.seekg(beginning);
                if(new_w){ windows[windows.size()-1].delete_window(); windows.pop_back(); used[cw] = 0; ++cw;}
                else { new_w = 1; ++nw; }
                for(uint16_t j=cw; j<arg->n; j++){ if(!used[j]){windows.push_back(KR_window(arg->w,primes[j])); cw=j; used[j]=1; break;}}
                if(j == arg->n){ --nw; }
            }
            
            for(j=0; j<nw; ++j){ windows[j].reset(); }
        }
     }
  }
  
  for(j=0; j<nw; ++j){ windows[j].delete_window(); }
  d->nw = nw;
  d->used = used;
  d->to_rmv = to_remove;
  f.close(); 
  return NULL;
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
Res parallel_parse_fasta(Args& arg, map<uint64_t,word_stats>& wf, bool mode)
{
    assert(arg.th>0);
    pthread_t t[arg.th];
    mt_data td[arg.th];
    // open filtered reads file
    FILE *filtered = open_aux_file(arg.inputFileName.c_str(),EXTFILT,"wb"); 
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
    size_t tstart = 1, tend = 1, nseq = 1, seq_b = 0;
    size_t nt = 1;
    // this loop scans the Fasta file in order to properly divide it up
    // for the threads, so that they don't accidently start in a ">" header.
    // As soon as a proper start and end position has been found, execute the thread 
    for(int i=1;i<th_sts.size()-1;i++){
        
        size_t start = th_sts[i];
        size_t end = th_sts[i+1];
        fseek(fp, start, SEEK_SET);
        for(size_t j=start; j<end; ++j){
            c = fgetc(fp);
            if(c == '>'){hp=j;sf=1;break;}
        }
        
        if(sf==1){
            sf=0;
            tend = hp;
            // prepare and execute thread j-1
            td[nt-1].wordFreq = &wf;
            td[nt-1].arg = &arg;
            td[nt-1].true_start = tstart;
            td[nt-1].true_end = tend;
            td[nt-1].words = 0;
            td[nt-1].parsed = 0;
            td[nt-1].th_id = nt-1;
            assert(td[nt-1].true_end <=size);
            if(mode){
            td[nt-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,nt-1,"wb");
            td[nt-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,nt-1,"wb");
            }
            if(mode==1){xpthread_create(&t[nt-1],NULL,&cyclic_mt_parse_fasta,&td[nt-1],__LINE__,__FILE__);}
            else{xpthread_create(&t[nt-1],NULL,&parallel_windows_n,&td[nt-1],__LINE__,__FILE__);}
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
    td[nt-1].th_id = nt-1;
    assert(td[nt-1].true_end <=size);
    if(mode){
    td[nt-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,nt-1,"wb");
    td[nt-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,nt-1,"wb");
    }
    if(mode==1){xpthread_create(&t[nt-1],NULL,&cyclic_mt_parse_fasta,&td[nt-1],__LINE__,__FILE__);}
    else{xpthread_create(&t[nt-1],NULL,&parallel_windows_n,&td[nt-1],__LINE__,__FILE__);}
    
    size_t tot_char=0; vector<bool> used(primes.size(),0); size_t tot_rem = 0; uint16_t wind_n;
    if(!mode){
        for(int i=0;i<nt;i++) {
            xpthread_join(t[i],NULL,__LINE__,__FILE__); 
            for(uint16_t j = 0; j < td[i].used.size(); j++){if(td[i].used[j]){used[j]=1;} } 
            tot_rem += td[i].to_rmv.size(); 
            to_remove_mt.push_back(td[i].to_rmv);
        }
        uint16_t cnt = 0; wind_n = count(used.begin(), used.end(), 1);
        for(uint16_t j = 0; j < used.size(); j++){if(used[j]){primes[cnt]=primes[j]; ++cnt;}}
        cout << "no. windows needed: " << wind_n << endl;
        cout << tot_rem << " sequences filtered" << endl;
    }
    else{
        // wait for the threads to finish (in order) and close output files
        for(int i=0;i<nt;i++) {
            xpthread_join(t[i],NULL,__LINE__,__FILE__); 
            if(td[i].filtered.size() > 0){
                size_t s = fwrite(&td[i].filtered[0], sizeof(char), td[i].filtered.size(), filtered);
                if(s!=td[i].filtered.size()) die("Error writing to filtered file"); 
            }
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
        //to_remove_mt.clear();
    }
    fclose(fp);
    if(fclose(filtered)!=0) die("Error closing filtered file");
    Res res; res.tot_char = tot_char; res.us_th = nt; res.nw = wind_n;
    return res;
}

#endif /* CIRCPFPMW_H */
