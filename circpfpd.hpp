/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   circpfpdmt.hpp
 * Author: hejans
 *
 * Created on April 14, 2021, 2:51 PM
 */

#ifndef CIRCPFPDMT_HPP
#define CIRCPFPDMT_HPP

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
  size_t parsed, words; //seqn;  // output
  FILE *parse, *o;
  uint64_t nd = 1;
  vector<bool> used;
  vector<uint64_t> to_rmv, pl;
  string filtered;
  int th_id;
} mt_data;

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
    while(  ((c = f.get()) != EOF) && current_pos < d.true_end ){
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
  KR_window krw(arg->w);
  uint8_t c; string word = ""; string fword = ""; string final_word = "";
  uint64_t start_char = 0, nseq = 1;
  size_t next_tr = 0, tr_i = 0; if(to_remove_mt[d->th_id].size()>0){ next_tr = to_remove_mt[d->th_id][tr_i]; }
  bool first_trigger = 0; int nd = d->pl.size();
  

  // skip the header
  uint64_t current_pos = d->true_start; uint64_t sp=0; size_t i, j;
  // step when we find the correct starting point
  skip_seq(nseq, current_pos, next_tr, tr_i, filt_seq, c, f, *d);
  // parse the sequence
  while( current_pos < d->true_end ) {
      c = f.get();
      ++current_pos;
      c = std::toupper(c);
      if(c > 64){ // A is 65 in ascii table. 
          if(c<= Dollar || c> 90) { die("Invalid char found in input file. Exiting..."); }
          word.append(1, c);
          if(first_trigger == 0){ fword.append(1, c); }
          uint64_t hash = krw.addchar(c);
          for(j=0; j<nd; ++j){
              if (hash%arg->p == d->pl[j] && krw.current == arg->w){
                  if(first_trigger==0){
                      first_trigger = 1, start_char = sp;
                      if(fwrite(&start_char,sizeof(start_char),1,d->o)!=1) die("offset write error");
                      word.erase(0,word.size() - arg->w);
                  }
                  else{
                      save_update_word(word,arg->w,*wordFreq,d->parse,0);
                      d->words++; 
                    }
                  break;
                }
            }
        ++sp;
      }
      else{ 
          if(c == '>' || current_pos >= d->true_end-1){
            d->parsed += krw.tot_char;
            for (i = 0; i < arg->w - 1; ++i) {
                c = fword[i];
                word.append(1, c);
                if(first_trigger == 0){fword.append(1, c);}
                uint64_t hash = krw.addchar(c);
                for(j=0; j<nd; ++j){   
                    if(hash%arg->p == d->pl[j] && krw.current == arg->w){
                        if(first_trigger==0){
                            first_trigger = 1, start_char = krw.tot_char;
                            if(fwrite(&start_char,sizeof(start_char),1,d->o)!=1) die("offset write error");
                            word.erase(0,word.size() - arg->w); 
                        }else{
                            save_update_word(word,arg->w,*wordFreq,d->parse,0);
                            d->words++;
                          }
                        break;
                    }
                }
            }
            if(first_trigger==0) { cerr << "No trigger strings found, use more windows." << endl; /*cout << "th_i "<< d->th_id << " " << nseq << endl;*/ exit(1); }
            final_word = word + fword.erase(0,arg->w - 1);
            save_update_word(final_word,arg->w,*wordFreq,d->parse,1);
            d->words++; 
            krw.reset();
            first_trigger = 0; start_char=0; sp=0;
            word = fword = final_word = ""; ++nseq;
            
            skip_seq(nseq, current_pos, next_tr, tr_i, filt_seq, c, f, *d);
          
        }
      }
    }
  
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
  uint64_t nd = 1; size_t beginning = d->true_start; 
  f.seekg(d->true_start); // move to the beginning of assigned region

  KR_window krw(arg->w);
  uint8_t c; /*pc = '\n';*/ string word(""); bool trg_f = 0, new_w = 0; uint64_t nseq = 1;
  vector<uint64_t> to_remove; uint16_t cw = 0;
  vector<uint64_t> pl; pl.push_back(0);
  vector<bool> used(arg->p,0); used[0]=1;
  vector<uint32_t> pf(arg->p,0);
  
  // skip the header
  uint64_t current_pos = d->true_start; uint64_t i, j;
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
          if(c<= 64 || c> 90){ die("Invalid char found in input file. Exiting...");}
          word.append(1, c);
          uint64_t hash = krw.addchar(c);
          if(krw.current == arg->w) ++pf[hash%arg->p];
          if(!trg_f){
            for(j=0; j<nd; ++j){
                if (hash%arg->p==pl[j] && krw.current == arg->w){ trg_f=1; break; }
              }
          }
      }
      else{ 
        if(c == '>' || current_pos >= d->true_end-1){
            if(!trg_f){
                for (i = 0; i < arg->w - 1; ++i) {
                    c = word[i];
                    uint64_t hash = krw.addchar(c);
                    if(krw.current == arg->w) ++pf[hash%arg->p];
                    for(j=0; j<nd; ++j){
                        if(hash%arg->p==pl[j] && krw.current == arg->w){ trg_f=1; break; }
                    }
                    if(trg_f){ break; }
                }
            }
            
            bool rem = 0;
            if(!trg_f && nd == arg->n){
                // if we reached the maximum window number and we are not testing a new window  
                trg_f = 1;  to_remove.push_back(nseq); rem = 1;
            }

            if(arg->c && !rem){
                vector<uint16_t> lps(word.size(),0); LPS(word,word.size(),lps);
                if(word.size() > 0 && word.size()%(word.size()-lps[lps.size()-1]) == 0 && lps[lps.size()-1] > 0){ trg_f = 1; to_remove.push_back(nseq); }
            }
            
            word.erase(0,word.size());
            
            if(!trg_f){ 
                vector<uint32_t> I(pf.size());
                std::iota(I.begin(),I.end(),0);
                sort(I.begin(),I.end(),[&](uint32_t i, uint32_t j){return pf[i]<pf[j];});
                
                for(i = 0;i<pf.size();++i){
                    size_t ind = I[i];
                    if(pf[ind]>0 && ind < arg->n){ if(!used[ind]){ used[ind]=1; ++nd; pl.push_back(ind); break; } }
                }
                
                if(i==pf.size()){
                    to_remove.push_back(nseq);
                }
            }
            
            ++nseq;
            while(  ((c = f.get()) != EOF) && current_pos < d->true_end ){
                ++current_pos;
                if(c=='\n'){break;}
            }
            beginning = current_pos; trg_f = 0;
            
            krw.reset();
            memset(&pf[0], 0, pf.size() * sizeof(pf[0]));
        }
     }
  }
  
  d->nd = nd;
  d->pl = pl;
  d->used = used;
  d->to_rmv = to_remove;
  f.close(); 
  return NULL;
}

// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
Res parallel_parse_fasta(Args& arg, map<uint64_t,word_stats>& wf, vector<uint64_t>& pl, bool mode)
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
            td[nt-1].pl = pl;
            td[nt-1].th_id = nt-1;
            assert(td[nt-1].true_end <=size);
            td[nt-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,nt-1,"wb");
            td[nt-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,nt-1,"wb");
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
    td[nt-1].pl = pl;
    td[nt-1].th_id = nt-1;
    assert(td[nt-1].true_end <=size);
    td[nt-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,nt-1,"wb");
    td[nt-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,nt-1,"wb");
    if(mode==1){xpthread_create(&t[nt-1],NULL,&cyclic_mt_parse_fasta,&td[nt-1],__LINE__,__FILE__);}
    else{xpthread_create(&t[nt-1],NULL,&parallel_windows_n,&td[nt-1],__LINE__,__FILE__);}
    
    size_t tot_char=0; uint64_t max_nd = 1; vector<bool> used(arg.p,0); size_t tot_rem = 0;
    if(!mode){
        for(int i=0;i<nt;++i) {
            xpthread_join(t[i],NULL,__LINE__,__FILE__); 
            max_nd = max(max_nd,td[i].nd); 
            for(uint16_t j = 0; j < td[i].used.size(); j++){if(td[i].used[j]){used[j]=1;} } 
            tot_rem += td[i].to_rmv.size(); 
            to_remove_mt.push_back(td[i].to_rmv);  
        }
        for(size_t j=0; j<used.size(); ++j){ if(used[j]){ pl.push_back(j); } }
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
        to_remove_mt.clear();
    }
    fclose(fp); 
    if(fclose(filtered)!=0) die("Error closing filtered file");
    Res res; res.tot_char = tot_char; res.us_th = nt; 
    return res;
}


#endif /* CIRCPFPDMT_HPP */

