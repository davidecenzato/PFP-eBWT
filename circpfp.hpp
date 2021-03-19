extern "C" {
#include "xerrors.h"
}
#include <vector>
pthread_mutex_t map_mutex = PTHREAD_MUTEX_INITIALIZER;

// struct shared via mt_parse
typedef struct {
  map<uint64_t,word_stats> *wordFreq; // shared dictionary
  Args *arg;       // command line input
  size_t true_start, true_end; // input
  size_t parsed, words;  // output
  FILE *parse, *o;
} mt_data;

bool is_valid_base(char base) {
    switch (std::toupper(base)) {
        case 'A': case 'C': case 'G': case 'T': return true;
        default: return false;
    }
}

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
  int c, pc = '\n'; string word = ""; string fword=""; string final_word="";
  uint64_t pos = 0, first_pos = 0, offset = 1, first_offset = 1, start_char = 0;
  bool first_trigger = 0;
  
  uint64_t current_pos = d->true_start; uint64_t i=0;
  while( (pc != EOF) && current_pos <= d->true_end + 1) {
      c = f.get();
      current_pos++;
      c = std::toupper(c);
      if(is_valid_base(c)){
          if(c<= Dollar) die("Invalid char found in input file. Exiting...");
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
        if(c == EOF || c == '>'){
            d->parsed += krw.tot_char;
            for (size_t i = 0; i < arg->w - 1; i++) {
                c = fword[i];
                word.append(1, c);
                pos++;
                uint64_t hash = krw.addchar(c);
                if(hash%arg->p==0){
                    save_update_word(word,arg->w,*wordFreq,d->parse,0);
                    d->words++;
                }
            }
            final_word = word + fword.erase(0,arg->w - 1);
            save_update_word(final_word,arg->w,*wordFreq,d->parse,1);
            d->words++; 
            krw.reset();
            first_trigger = 0; start_char=0; i=0;
            word = fword = final_word = "";
            while(  ((c = f.get()) != EOF) ){
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
uint64_t parallel_parse_fasta(Args& arg, map<uint64_t,word_stats>& wf)
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
  std::vector<size_t> th_sts(arg.th);
  th_sts[0] = 0;
  for (int i = 1; i < arg.th; ++i) {
    th_sts[i] = (size_t) (size / arg.th) * i;
  }
  
  cout << "Thread: " << arg.th << endl;
  cout << "Total size: " << size << endl;
  cout << "------------------------" << endl;
  int j = 1, c = 0;
  int header_len = 1;
  
  // this loop scans till he find the starting point of the first sequence
  int current_pos = 0;
  int true_start = 0, true_end = 0; 
  while ( ((c = fgetc(fp)) != EOF) ) {
      current_pos ++;
      if(is_valid_base(c)){break;}
      true_start ++;
  }
  // this loop scans the Fasta file in order to properly divide it up
  // for the threads, so that they don't accidently start in a ">" header.
  // As soon as a proper start and end position has been found, execute the thread
  while ( ((c = fgetc(fp)) != EOF) ) {
    current_pos++;
    if (j == arg.th + 1) break;
    if (c == '>' ){
        true_end = current_pos - 2; 
        header_len = 1;
        while ( ((c = fgetc(fp)) != EOF) ) {
            current_pos ++;
            if(is_valid_base(c)){break;}
            header_len ++;
        }
    }        
    
    if(current_pos == th_sts[j] && j < arg.th + 1 && true_end > true_start){
        
        // prepare and execute thread j-1
        td[j-1].wordFreq = &wf;
        td[j-1].arg = &arg;
        td[j-1].true_start = true_start;
        td[j-1].true_end = true_end;
        td[j-1].words = td[j-1].parsed = 0;
        assert(td[j-1].true_end <=size);
        td[j-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,j-1,"wb");
        td[j-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,j-1,"wb");
        xpthread_create(&t[j-1],NULL,&cyclic_mt_parse_fasta,&td[j-1],__LINE__,__FILE__);
        true_start = true_end + header_len + 1;
        j++;
    }
  }
  
    true_end = current_pos-1;
    
    td[j-1].wordFreq = &wf;
    td[j-1].arg = &arg;
    td[j-1].true_start = true_start;
    td[j-1].true_end = true_end;
    td[j-1].words = td[j-1].parsed = 0;
    assert(td[j-1].true_end<=size);
    td[j-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,j-1,"wb");
    td[j-1].o = open_aux_file_num(arg.inputFileName.c_str(),EXTOFF0,j-1,"wb");
    xpthread_create(&t[j-1],NULL,&cyclic_mt_parse_fasta,&td[j-1],__LINE__,__FILE__);
    fclose(fp);
    
    // wait for the threads to finish (in order) and close output files
    size_t tot_char=0;
    if(j != arg.th){throw new std::runtime_error("Cannot divide input file for threads, check input file format.");}
    for(int i=0;i<arg.th;i++) {
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
    return tot_char;
}

