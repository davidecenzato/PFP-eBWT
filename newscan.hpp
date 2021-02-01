extern "C" {
#include "xerrors.h"
}
#include <vector>
pthread_mutex_t map_mutex = PTHREAD_MUTEX_INITIALIZER;

// struct shared via mt_parse
typedef struct {
  map<uint64_t,word_stats> *wordFreq; // shared dictionary
  Args *arg;       // command line input
  size_t true_start, true_end, start, end; // input
  size_t skipped, parsed, words;  // output
  FILE *parse, *last, *sa;
} mt_data;


void *mt_parse(void *dx)
{
  // extract input data
  mt_data *d = (mt_data *) dx;
  Args *arg = d->arg;
  map<uint64_t,word_stats> *wordFreq = d->wordFreq;

  if(arg->verbose>1)
    printf("Scanning from %ld, size %ld\n",d->start,d->end-d->start);

  // open input file
  ifstream f(arg->inputFileName);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open file " + arg->inputFileName);
  }

  // prepare for parsing
  f.seekg(d->start); // move to the begining of assigned region
  KR_window krw(arg->w);
  int c; string word = "";
  d->skipped = d->parsed = d->words = 0;
  if(d->start==0) {
    word.append(1,Dollar);// no need to reach the next kr-window
  }
  else {   // reach the next breaking window
    while( (c = f.get()) != EOF ) {
      if(c<=Dollar) die("Invalid char found in input file. Exiting...");
      d->skipped++;
      if(d->true_start + d->skipped == d->true_end + arg->w) {f.close(); return NULL;}
      word.append(1,c);
      uint64_t hash = krw.addchar(c);
      if(hash%arg->p==0 && d->skipped >= arg->w) break;
    }
    if(c==EOF) {f.close(); return NULL;} // reached EOF without finding a breaking point nothing to do
    d->parsed = arg->w;   // the kr-window is part of the next word
    d->skipped -= arg->w; // ... so w less chars have been skipped
    word.erase(0,word.size() - arg->w);// keep only the last w chars
  }

  // there is some parsing to do
  uint64_t pos = d->start;             // ending position+1 in text of previous word
  if(pos>0) pos+= d->skipped+ arg->w;  // or 0 for the first word
  assert(IBYTES<=sizeof(pos)); // IBYTES bytes of pos are written to the sa info file
  while( (c = f.get()) != EOF ) {
    if(c<=Dollar) die("Invalid char found in input file. Exiting...");
    word.append(1,c);
    uint64_t hash = krw.addchar(c);
    d->parsed++;
    if(hash%arg->p==0 && d->parsed>arg->w) {
      // end of word, save it and write its full hash to the output file
      // pos is the ending position+1 of previous word and is updated in the next call
      save_update_word(word,arg->w,*wordFreq,d->parse,d->last,d->sa,pos);
      d->words++;
      if(d->true_start+d->skipped+d->parsed>=d->true_end+arg->w) {f.close(); return NULL;}
    }
  }
  // end of file reached
  // virtually add w null chars at the end of the file and add the last word in the dict
  word.append(arg->w,Dollar);
  save_update_word(word,arg->w,*wordFreq,d->parse,d->last,d->sa,pos);
  // close input file and return
  f.close();
  return NULL;
}


// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
uint64_t mt_process_file(Args& arg, map<uint64_t,word_stats>& wf)
{
  // get input file size
  ifstream f(arg.inputFileName, std::ifstream::ate);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " +arg.inputFileName);
  }
  size_t size = f.tellg();
  f.close();

  // prepare and execute threads
  assert(arg.th>0);
  pthread_t t[arg.th];
  mt_data td[arg.th];
  for(int i=0;i<arg.th;i++) {
    td[i].wordFreq = &wf;
    td[i].arg = &arg;
    td[i].start = i*(size/arg.th); // range start
    td[i].end = (i+1==arg.th) ? size : (i+1)*(size/arg.th); // range end
    assert(td[i].end<=size);
    // open the 1st pass parsing file
    td[i].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,i,"wb");
    // open output file containing the char at position -(w+1) of each word
    td[i].last = open_aux_file_num(arg.inputFileName.c_str(),EXTLST,i,"wb");
    // if requested open file containing the ending position+1 of each word
    td[i].sa = arg.SAinfo ?open_aux_file_num(arg.inputFileName.c_str(),EXTSAI,i,"wb") : NULL;
    xpthread_create(&t[i],NULL,&mt_parse,&td[i],__LINE__,__FILE__);
  }

  // wait for the threads to finish (in order) and close output files
  size_t tot_char=0;
  for(int i=0;i<arg.th;i++) {
    xpthread_join(t[i],NULL,__LINE__,__FILE__);
    if(arg.verbose) {
      cout << "s:" << td[i].start << "  e:" << td[i].end << "  pa:";
      cout << td[i].parsed << "  sk:" << td[i].skipped << "  wo:" << td[i].words << endl;
    }
    // close thread-specific output files
    fclose(td[i].parse);
    fclose(td[i].last);
    if(td[i].sa) fclose(td[i].sa);
    if(td[i].words>0) {
      // extra check
      assert(td[i].parsed>arg.w);
      tot_char += td[i].parsed - (i!=0? arg.w: 0); //parsed - overlapping
    }
    else assert(i>0); // the first thread must produce some words
  }
  assert(tot_char==size);
  return size;
}


// modified from mt_parse to skip newlines and fasta header lines (ie. lines starting with ">")
void *mt_parse_fasta(void *dx)
{
  // extract input data
  mt_data *d = (mt_data *) dx;
  Args *arg = d->arg;
  map<uint64_t,word_stats> *wordFreq = d->wordFreq;

  if(arg->verbose>1)
    printf("Scanning from %ld, size %ld as a FASTA record\n",d->start,d->end-d->start);

  // open input file
  ifstream f(arg->inputFileName);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open file " + arg->inputFileName);
  }

  // prepare for parsing
  f.seekg(d->start); // move to the begining of assigned region
  KR_window krw(arg->w);
  int c, pc = '\n'; string word = "";
  d->skipped = d->parsed = d->words = 0;
  int IN_HEADER = 1;
  // for fasta files, figure out if start in a different fashion
  if (d->true_start == 0) {
    word.append(1, Dollar);
  }
  else {   // reach the next breaking window
    while(  ((c = f.get()) != EOF) ) {
      if (pc == '\n') IN_HEADER = (c == '>');
      if (c != '\n' && !IN_HEADER)  {
        c = std::toupper(c);
        if(c<=Dollar) die("Invalid char found in input file. Exiting...");
        d->skipped++;
        if(d->true_start + d->skipped == d->true_end + arg->w) {
          f.close();
          return NULL;
        }
        word.append(1,c);
        uint64_t hash = krw.addchar(c);
        if(hash%arg->p==0 && d->skipped >= arg->w) break;
      }
      pc = c;
    }
    if(c==EOF) {f.close(); return NULL;} // reached EOF without finding a breaking point nothing to do
    d->parsed = arg->w;   // the kr-window is part of the next word
    d->skipped -= arg->w; // ... so w less chars have been skipped
    word.erase(0,word.size() - arg->w);// keep only the last w chars
  }

  // there is some parsing to do
  // uint64_t pos = d->start;             // ending position+1 in text of previous word
  // if(pos>0) pos+= d->skipped+ arg->w;  // or 0 for the first word
  uint64_t pos = d->true_start;
  if (pos>0) pos += d->skipped + arg->w;
  assert(IBYTES<=sizeof(pos)); // IBYTES bytes of pos are written to the sa info file
  // note: IN_HEADER state carries over
  while(  (c = f.get()) != EOF ) {
    if (pc == '\n') IN_HEADER = (c == '>');
    if (c != '\n' && !IN_HEADER)  {
      c = std::toupper(c);
      if(c<=Dollar) die("Invalid char found in input file. Exiting...");
      word.append(1,c);
      uint64_t hash = krw.addchar(c);
      d->parsed++;
      if(hash%arg->p==0 && d->parsed>arg->w) {
        // end of word, save it and write its full hash to the output file
        // pos is the ending position+1 of previous word and is updated in the next call
        save_update_word(word,arg->w,*wordFreq,d->parse,d->last,d->sa,pos);
        d->words++;
        if(d->true_start+d->skipped+d->parsed>=d->true_end+arg->w) {
          f.close(); return NULL;
        }
      }
    }
    pc = c;
  }
  // end of file reached
  // virtually add w null chars at the end of the file and add the last word in the dict
  word.append(arg->w,Dollar);
  save_update_word(word,arg->w,*wordFreq,d->parse,d->last,d->sa,pos);
  // close input file and return
  f.close();
  return NULL;
}


// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
uint64_t mt_process_file_fasta(Args& arg, map<uint64_t,word_stats>& wf)
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
  std::vector<size_t> true_starts(arg.th);
  std::vector<size_t> true_ends(arg.th);
  th_sts[0] = 0;
  for (int i = 1; i < arg.th; ++i) {
    th_sts[i] = (size_t) (size / arg.th) * i;
  }
  int IN_HEADER = 1;
  size_t true_pos = 0, file_pos = 0;
  int j = 0, pc = 0, c = 0;
  // this loop scans the fasta file in order to properly divvy it up
  // for the threads, so they they don't accidently start in a ">" header.
  // As soon as a proper start and end position has been found, execute the thread
  while ( ((c = fgetc(fp)) != EOF) ) {
    if (j == arg.th) break;
    if (pc == '\n') IN_HEADER = (c == '>');
    if (file_pos == th_sts[j]) {
      if (IN_HEADER) th_sts[j]++;
      else {
        true_starts[j] = true_pos;
        if (j) {
          true_ends[j-1] = true_pos;
          // prepare and execute thread j-1
          td[j-1].wordFreq = &wf;
          td[j-1].arg = &arg;
          td[j-1].true_start = true_starts[j-1];
          td[j-1].true_end = true_ends[j-1];
          td[j-1].start = th_sts[j-1]; // range start
          td[j-1].end = (j==arg.th) ? size : th_sts[j];
          assert(td[j-1].end<=size);
          td[j-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,j-1,"wb");
          td[j-1].last = open_aux_file_num(arg.inputFileName.c_str(),EXTLST,j-1,"wb");
          td[j-1].sa = arg.SAinfo ?open_aux_file_num(arg.inputFileName.c_str(),EXTSAI,j-1,"wb") : NULL;
          xpthread_create(&t[j-1],NULL,&mt_parse_fasta,&td[j-1],__LINE__,__FILE__);
        }
        ++j;
        // check if previous thread spilled over computed start position of
        // next thread only possible in rare situations.
        if (j && j < arg.th && th_sts[j-1] >= th_sts[j]) { 
          th_sts[j] = file_pos + 1;
        }
      }
    }
    if (!IN_HEADER && c != '\n') ++true_pos;
    pc = c;
    ++file_pos;
  }
  assert(j == arg.th);
  // execute the last thread
  true_ends[j-1] = size;
  td[j-1].wordFreq = &wf;
  td[j-1].arg = &arg;
  td[j-1].true_start = true_starts[j-1];
  td[j-1].true_end = true_ends[j-1];
  td[j-1].start = th_sts[j-1]; // range start
  td[j-1].end = size;
  td[j-1].parse = open_aux_file_num(arg.inputFileName.c_str(),EXTPARS0,j-1,"wb");
  td[j-1].last = open_aux_file_num(arg.inputFileName.c_str(),EXTLST,j-1,"wb");
  td[j-1].sa = arg.SAinfo ?open_aux_file_num(arg.inputFileName.c_str(),EXTSAI,j-1,"wb") : NULL;
  xpthread_create(&t[j-1],NULL,&mt_parse_fasta,&td[j-1],__LINE__,__FILE__);
  fclose(fp);
  // TODO: we might have a bunch of threads at the end that start at the same
  // position! this can happen in the following situation with >2 threads:
  // >HEADER1
  // AATACBTAC
  // >HEADER2
  // >HEADER3
  // >HEADER4
  // I guess this is fine, since those threads aren't likely to  return
  // anything useful anyway

  // wait for the threads to finish (in order) and close output files
  size_t tot_char=0;
  for(int i=0;i<arg.th;i++) {
    xpthread_join(t[i],NULL,__LINE__,__FILE__);
    if(arg.verbose) {
      cout << "s:" << td[i].start << "  e:" << td[i].end << "  pa:";
      cout << td[i].parsed << "  sk:" << td[i].skipped << "  wo:" << td[i].words << endl;
    }
    // close thread-specific output files
    fclose(td[i].parse);
    fclose(td[i].last);
    if(td[i].sa) fclose(td[i].sa);
    if(td[i].words>0) {
      // extra check
      assert(td[i].parsed>arg.w);
      tot_char += td[i].parsed - (i!=0? arg.w: 0); //parsed - overlapping
    }
    else assert(i>0); // the first thread must produce some words
  }
  cout << "total parsed characters: " << tot_char <<  ", file size: " << size << endl;
  // assert(tot_char==size); // this assert statement doesn't apply to fasta files
  return size;
}
