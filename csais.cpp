/*
 * Circular SAIS and SACA implementations to compute the circular SA of an integer vector.
 * 
 * This code is adapted from https://github.com/kurpicz/saca-bench/blob/master/sa-is/sais.cpp
 * which is the original code of the SA-IS algorithm listed below, and 
 * from https://github.com/felipelouza/gsa-is/blob/master/gsacak.c which is an implementation
 * of the GSACA-K algorithm.
 *
 */
// This is the sample code for the SA-IS algorithm presented in
// our article "Two Efficient Algorithms for Linear Suffix Array Construction"
// (co-authored with Sen Zhang and Wai Hong Chan),
// which can be retrieved at: http://code.google.com/p/ge-nong/


#include "csais.h"

using namespace std;
using namespace sdsl;

#if P64
  const uint_s EMPTY=0xffffffffffffffff; 
#else
  const uint_s EMPTY=0xffffffff; 
#endif

#define chr(i) (cs==sizeof(uint_s)?((uint_s *)s)[i]:((uint_p *)s)[i])
// functions only for SAIS algorithm
unsigned char mask[]={0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};
#define tget(i) ( (t[(i)/8]&mask[(i)%8]) ? 1 : 0 )
#define tset(i, b) t[(i)/8]=(b) ? (mask[(i)%8]|t[(i)/8]) : ((~mask[(i)%8])&t[(i)/8])
#define isLMS(i,j) ((tget(i) && !tget(j)))

// compute the head or end of each bucket
void getBuckets(uint_s *s, uint_s *bkt, size_t n, size_t K, size_t cs, bool end) { 
  size_t i, sum=0;
  for(i=0; i<K; i++) bkt[i]=0; // clear all buckets
  for(i=0; i<n; i++) bkt[chr(i)]++; // compute the size of each bucket
  for(i=0; i<K; i++) { sum+=bkt[i]; bkt[i]= end ? sum-1 : sum-bkt[i]; }
}

// induce L suffixes for SAIS
void induceSAl(unsigned char *t, uint_s *SA, uint_s *s, uint_s *bkt, uint_s *l_bkt, sd_vector<>::rank_1_type& r_s,
               sd_vector<>::select_1_type& s_s, sd_vector<>& b_s, size_t n, size_t K, size_t cs, int level, vector<uint_s>& star) { 
    
    size_t i, j, rank;
    getBuckets(s, bkt, n, K, cs, false); // find heads of buckets
    
    for(i=0;i<star.size();i++){
      j=star[i]; SA[(l_bkt[chr(j)]++)+bkt[chr(j)]]=j;
    }

    for(i=0; i<n; i++){
        if(SA[i]!=EMPTY) {
            j=SA[i];
            if(b_s[j]==1){
                rank=r_s(j+1);j=s_s(rank+1)-1;
                if(!tget(j) && j>=0) SA[bkt[chr(j)]++]=j;    
            }else{
                if(!tget(j-1) && j>=0) SA[bkt[chr(j-1)]++]=j-1; 
            }    
        }
    }
}

// induce L suffixes for SACA
void induceL(uint_s *SA, uint_s *s, uint_s *bkt, sd_vector<>::rank_1_type& r_s, 
             sd_vector<>::select_1_type& s_s, sd_vector<>& b_s, size_t n, size_t K, 
             size_t cs, bool phase, int level) { 
    
    size_t i, j, m, rank;
    getBuckets(s, bkt, n, K, cs, false); // find heads of buckets

    for(i=0; i<n; ++i){
        if(SA[i]!=EMPTY) {
            j=SA[i];
            if(b_s[j]==1){
                rank=r_s(j+1);m=s_s(rank+1)-1;
                if(chr(m) >= chr(j)){ SA[bkt[chr(m)]++]=m; if(phase){ SA[i] = EMPTY; } }
            }else{
                if(chr(j-1) >= chr(j)){ SA[bkt[chr(j-1)]++]=j-1; if(phase){ SA[i] = EMPTY; } } 
            }
        }
    }
}

// induce S suffixes for SAIS
void induceSAs(unsigned char *t, uint_s *SA, uint_s *s, uint_s *bkt, sd_vector<>::rank_1_type& r_s,
               sd_vector<>::select_1_type& s_s, sd_vector<>& b_s, size_t n, size_t K, size_t cs, int level) { 
    
  size_t i, j, rank;
  getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  for(i=0;i<n;i++){
    size_t ni = n-i-1;
    if(SA[ni]!=EMPTY) {
        j=SA[ni];
        if(b_s[j]==1){
            rank=r_s(j+1);j=s_s(rank+1)-1;
            if(tget(j) && j>=0) SA[bkt[chr(j)]--]=j;    
        }else{
            if(tget(j-1) && j>=0) SA[bkt[chr(j-1)]--]=j-1; 
            }    
        }   
    }
}

// induce S suffixes for SACA
void induceS(uint_s *SA, uint_s *s, uint_s *bkt, vector<uint_s>& singletons, sd_vector<>::rank_1_type& r_s,
             sd_vector<>::select_1_type& s_s, sd_vector<>& b_s, size_t n, size_t K, size_t cs, bool phase, int level) { 
    
    size_t i, j, m, rank;
    getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
    for(i=0; i<n; ++i){
        size_t ni = n-i-1;
        if(SA[ni]!=EMPTY) {
            j=SA[ni];
            if(b_s[j]==1){
                rank=r_s(j+1);m=s_s(rank+1)-1;
                if(chr(m) <= chr(j) && bkt[chr(m)]<ni){ SA[bkt[chr(m)]--]=m; if(phase){ SA[ni] = EMPTY; } }    
             }else{
                if(chr(j-1) <= chr(j) && bkt[chr(j-1)]<ni){ SA[bkt[chr(j-1)]--]=j-1; if(phase){ SA[ni] = EMPTY; } } 
            }  
        }   
    }
  // insert singletons
  if(!phase){ for(i = 0;i < singletons.size(); ++i){ SA[bkt[chr(singletons[i])]--]=singletons[i]; } }
}

// compute length of a LMS substring for SACA
uint_s LMSlength(uint_s *s, uint_s sb, uint_s eb, int level, uint_s x, size_t cs) {
  // return the length of the LMS substring
  uint_s len=1, i=1;  
  uint_s prev = x, pos = 0;
  if(x == eb){ pos = sb; }else{ pos = x+1; }
  // S suffixes
  while(true) {
    if(chr(pos)<chr(prev)){ break; }
    ++i;
    if(pos == eb) { pos = sb; } else{ ++pos;  }
    if(prev == eb){ prev = sb; }else{ ++prev; }
  }  
  // L suffixes
  while(true) {
    if(chr(pos)>chr(prev)){ break; }
    if(chr(pos)<chr(prev)){ len=i; }
    ++i;
    if(pos == eb) { pos = sb; } else{ ++pos;  }
    if(prev == eb){ prev = sb; }else{ ++prev; }
  }
  // return the length
  return len+1;
}

// find the circular suffix array SA of s[0..n-1] in {0..K-1}^n
// level starts from 0
/** @brief computes the circular suffix array of string s[0..n-1] using cSAIS algo
 *
 *  @param s	   input string 
 *  @param SA    suffix array 
 *  @param n	   string length
 *  @param K	   alphabet size
 *  @param cs	   integer size
 *  @param level recursion level, debug only
 *  @param b_s   starting positions bit vector
 *  @return      None
 */
void cSAIS(uint_s *s, uint_s *SA, size_t n, size_t K, size_t cs, int level, sd_vector<> &b_s) {
  size_t i, j, nseq, rank;
  size_t sb, eb, fm, len; // maybe uint_s
  sd_vector<>::rank_1_type r_s = sd_vector<>::rank_1_type(&b_s);
  sd_vector<>::select_1_type s_s = sd_vector<>::select_1_type(&b_s);

  nseq = r_s(n);
  unsigned char *t=(unsigned char *)malloc(n/8+1); // LS-type array in bits
  
  // stage 1: reduce the problem by at least 1/2
  
  uint_s *bkt = (uint_s *)malloc(sizeof(uint_s)*K); // bucket counters
  
  uint_s *l_bkt = (uint_s *)malloc(sizeof(uint_s)*K); // l types bucket counters
  for(i=0; i<K; i++) l_bkt[i]=0; 
  
  for(i=0; i<n; i++) SA[i]=EMPTY;
  
  vector<uint_s> sgt; // singletons vector
  for(i=0;i<nseq;i++){
      sb = s_s(i+1), eb = s_s(i+2)-1, fm = eb+1, len=eb-sb+1;
      if(len==1){tset(sb,1);sgt.push_back(sb);}
      else{
        // find first type
        for(j=eb;j>sb;j--){
            if(chr(j)!=chr(j-1)){
              if(chr(j)>chr(j-1)){fm=j;tset(j,0);break;}
              else{fm=j;tset(j,1);break;}
            }else{tset(j,0);} 
        }
        // mismatch not found!
        if(fm==eb+1){cerr << "Error! Sequence without a mismatch detected." << endl; exit(1);}
        else{
          // Classify the type of each character
          // Count the number of L types
          uint_s pos=fm-1, prev=fm;
          bool type = (chr(pos)<chr(prev));
          tset(pos, type); if(!type) l_bkt[chr(pos)]++;
          while(pos != fm){
              if(pos-sb > 0){ pos--; }
              else {pos = eb;}
              if(prev-sb > 0){ prev--; }
              else {prev = eb;}
              
              type=(chr(pos)<chr(prev) || (chr(pos)==chr(prev) && tget(prev)==1));
              if(!type){l_bkt[chr(pos)]++;}
              tset(pos,type);
          }
        }
      }
    }
  
  // sort all the S-substrings
  getBuckets(s, bkt, n, K, cs, true); // find ends of buckets

  // put S star suffixes in their correct positions
  for(i=0; i<n; ++i)
  {
      size_t ni = n-i-1;
      if(b_s[ni]==0){
          if(isLMS(ni,ni-1)) SA[bkt[chr(ni)]--]=ni;
      }else{
          rank=r_s(ni+1), eb=s_s(rank+1)-1;
          if(isLMS(ni,eb)) SA[bkt[chr(ni)]--]=ni;
      }
  }
  
  induceSAl(t, SA, s, bkt, l_bkt, r_s, s_s, b_s, n, K, cs, level, sgt);
  induceSAs(t, SA, s, bkt, r_s, s_s, b_s, n, K, cs, level);  
  
  free(bkt);
  free(l_bkt);
  sgt.clear();
  
  // compact all the sorted substrings into the first n1 items of s
  size_t n1=0;
  uint_s *sts = (uint_s *)malloc(sizeof(uint_s)*(nseq+1)); for(i=0; i<(nseq+1); i++) sts[i]=0;
  for(i=0;i<n;i++){
      uint_s pos=SA[i];
      rank=r_s(pos+1),sb=s_s(rank),eb=s_s(rank+1);
      if(pos==sb){
          if(isLMS(sb,eb-1)){SA[n1++]=pos;sts[rank]++;}
      }else{
          if(isLMS(pos,pos-1)){SA[n1++]=pos;sts[rank]++;}
      }
  }
  
  //create the bit vector of the new string
  sd_vector<> nb_s; vector<size_t> onset; size_t sum=0;
  onset.push_back(0);
  for(size_t i=1;i<(nseq+1);++i){sum+=sts[i]; if(sts[i]>0){onset.push_back(sum);}}
  sd_vector_builder builder(n1+1,onset.size());
  for(auto idx: onset){builder.set(idx);}
  nb_s = sd_vector<>(builder);
  free(sts);
  
  // Init the name array buffer
  for(i=n1; i<n; i++) SA[i]=EMPTY;
  
  // find the lexicographic names of all substrings
  // insert the first substring
  uint_s name=1;
  if(n1 > 0){
    uint_s *names = (uint_s *)malloc(sizeof(uint_s)*n); for(i=0; i<n; i++) names[i]=EMPTY;
    uint_s prev = SA[0]; uint_s pos = SA[0]; 
    names[pos]=name-1;
    size_t rank_prev = r_s(prev+1); size_t sb_prev = s_s(rank_prev);
    for(i=1; i<n1; ++i) {
          size_t pp=0, pv=0;
          pos=SA[i]; bool diff=false;
          rank = r_s(pos+1); sb = s_s(rank);
          
          size_t ind = pos, ind_prev = prev;
          //check starting characters
          if(chr(ind)!=chr(ind_prev) || tget(ind)!=tget(ind_prev)){ diff=true; }
          pp=ind, pv=ind_prev;

          while(!diff){
              //check the other characters in cyclic way
              if(b_s[ind_prev+1]==0){ind_prev++;}
              else{ ind_prev = sb_prev; }
              if(b_s[ind+1]==0){ind++;}
              else{ ind = sb; }
              if(chr(ind)!=chr(ind_prev) || tget(ind)!=tget(ind_prev)) 
              { 
                  diff=true; break;
              }
              else{ // we stop the comparison when we find another s*
                if(isLMS(ind,pp) || isLMS(ind_prev,pv))
                  break;
              }
              pp=ind,pv=ind_prev;
          }
        
        if(diff){ name++; prev=pos; sb_prev = sb; }
        names[pos]=name-1;
    }
  
    j = n-1;
    for(i=0;i<n;i++){
        if(names[n-i-1]!=EMPTY){SA[j--]=names[n-i-1];}
    }
      
    free(names);
  }
  
  // s1 is done now
  uint_s *SA1=SA, *s1=SA+n-n1;  
  
  // stage 2: solve the reduced problem
  // recurse if names are not yet unique
  if(name<n1) {
      cSAIS((uint_s*)s1, SA1, n1, name, sizeof(uint_s), level+1, nb_s);
  } else { // stop the recursion, generate the suffix array of s1 directly
      for(i=0; i<n1; i++){ SA1[s1[i]] = i; }
  }
  
  
  // stage 3: induce the result for the original problem
  
  bkt = (uint_s *)malloc(sizeof(uint_s)*K); // bucket counters
  l_bkt = (uint_s *)malloc(sizeof(uint_s)*K); for(i=0; i<K; i++) l_bkt[i]=0;

  // put all left-most S characters into their buckets
  j=0;
  for(i=0;i<nseq;i++){
      sb=s_s(i+1), eb=s_s(i+2), len= eb-sb;
      if(len==1){sgt.push_back(sb);}
      else
      {
          if(isLMS(sb,eb-1)){ s1[j++]=sb; }
          else{ if(tget(sb)==0){ l_bkt[chr(sb)]++; } }
          for(size_t k=sb+1;k<eb;k++){
              if(isLMS(k,k-1)){ s1[j++]=k; }
              else{ if(tget(k)==0) { l_bkt[chr(k)]++; } } 
          }
      }
  }
  
  assert(j==n1);
    
  for(i=0; i<n1; i++) {SA1[i]=s1[SA1[i]];} // get index in s1
   
  if(n1>0){ // if at least one s* is not a singleton
      for(i=n1; i<n; i++) SA[i]=EMPTY; // init SA[n1..n-1]
      getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  }
  
  for(i=0; i < n1; i++){
      j=SA[n1-i-1]; SA[n1-i-1]=EMPTY;
      SA[bkt[chr(j)]--]=j;
  }

  induceSAl(t, SA, s, bkt, l_bkt, r_s, s_s, b_s, n, K, cs, level, sgt); 
  induceSAs(t, SA, s, bkt, r_s, s_s, b_s, n, K, cs, level);
  
  onset.clear();
  free(bkt); 
  free(l_bkt);
  free(t);
  
}

// find the circular suffix array SA of s[0..n-1] in {0..K-1}^n
// level starts from 0
/** @brief computes the circular suffix array of string s[0..n-1] using cSACA algo
 *
 *  @param s     input string 
 *  @param SA    suffix array 
 *  @param n     string length
 *  @param K     alphabet size
 *  @param cs    integer size
 *  @param level recursion level, debug only
 *  @param b_s   starting positions bit vector
 *  @return      None
 */
void cSACA(uint_s *s, uint_s *SA, size_t n, size_t K, size_t cs, int level, sd_vector<> &b_s) {
  size_t i, j, m, nseq, rank;
  size_t sb, eb, fl, len; 
  // initialize support for rank and select 
  sd_vector<>::rank_1_type r_s = sd_vector<>::rank_1_type(&b_s);
  sd_vector<>::select_1_type s_s = sd_vector<>::select_1_type(&b_s);
  // singletons vector
  vector<uint_s> sn;
  // no. sequences
  nseq = r_s(n);
  // stage 1: reduce the problem by at least 1/2
  uint_s *bkt = (uint_s *)malloc(sizeof(uint_s)*K); // bucket counters
  
  for(i=0; i<n; i++) SA[i]=EMPTY; // initialize SA values to -1
  getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  
  //initialize onset vector
  vector<uint_s> onset; size_t st = 0;
  onset.reserve(nseq+1); onset.push_back(st);
  // place S* suffixes in their buckets
  for(i=0; i<nseq; ++i){
      sb = s_s(i+1), eb = s_s(i+2)-1, fl = eb+1, len=eb-sb+1;
      assert(len > 0);
      bool type;
      if(len > 1){
        // find first type
        for(j=eb; j>sb; --j){
            if(chr(j)!=chr(j-1)){
              if(chr(j)>chr(j-1)){fl=j; type=0; break;}
              else{fl=j-1; type=0; break;}
            } 
        }
        // mismatch not found!
        if(fl==eb+1){cerr << "Error! Sequence without a mismatch detected." << endl; exit(1);}
        else{
          // Classify the S* suffixes
          uint_s pos, prev; prev=fl; 
          if(prev==sb){ pos = eb; } else{ pos = fl-1; }
          while(pos != fl){
              if(chr(pos) > chr(prev)){ if(type){ SA[bkt[chr(prev)]--] = prev; ++st; type = 0; } }
              else{ 
                  if(chr(pos)<chr(prev)){ type = 1; }
              }
              // Skip to next chars
              if(prev==sb){ prev = eb; } else{ --prev; }
              if(pos==sb) { pos = eb;  } else{ --pos;  }
           }
          // check the last suffix
          if(chr(pos) > chr(prev)){ if(type){ SA[bkt[chr(prev)]--] = prev; ++st; type = 0; } }
        }
        onset.push_back(st); 
     }
  }
  
  // Induce L ans S suffixes to sort the LMS substrings
  induceL(SA, s, bkt, r_s, s_s, b_s, n, K, cs, true, level);
  induceS(SA, s, bkt, sn, r_s, s_s, b_s, n, K, cs, true, level); 
  
  free(bkt); // free bucket vector
  
  // compact all the sorted substrings into the first n1 items of s
  size_t n1=0;
  for(i=0; i<n; ++i){
      uint_s pos=SA[i];
      if(pos != EMPTY){ SA[n1++]=pos; }
  }
  
  // build new bit vector
  sd_vector_builder builder(n1+1,onset.size());
  for(auto idx: onset){ builder.set(idx); }
  sd_vector<> nb_s = sd_vector<>(builder);
  
  // Init the name array buffer
  for(i=n1; i<n; ++i) SA[i]=EMPTY;
  
  // find the lexicographic names of all LMS substrings
  uint_s name=1;
  if(n1 > 0){
    // initialize names array
    uint_s *names = (uint_s *)malloc(sizeof(uint_s)*n); for(i=0; i<n; i++) names[i]=EMPTY;
    // insert the first LMS substring
    uint_s pos = SA[0];
    names[pos]=name-1;
    size_t rank = r_s(pos+1), sb = s_s(rank), eb = s_s(rank+1)-1;
    uint_s pre_len = LMSlength(s, sb, eb, level, pos, cs); 
    uint_s prev = pos;  size_t pre_sb = sb, pre_eb = eb;
    // for all S* suffixes
    for(i=1; i<n1; ++i) {
          pos=SA[i]; bool diff=false;
          rank = r_s(pos+1); sb = s_s(rank); eb = s_s(rank+1)-1;
          uint_s len = LMSlength(s, sb, eb, level, pos, cs); 
          // if the LMS length are different skip and increase name counter
          if(len != pre_len){ diff = true; }
          else{
              // if same length compare all the characters till you find a mismatch
              uint_s tp = pos, pre_tp = prev;
              for(j=0;j<len;++j){ 
                  if(chr(pre_tp)!=chr(tp)){ diff=true; break; }
                  if(tp==eb){ tp=sb; }else{ ++tp;}
                  if(pre_tp==pre_eb){ pre_tp=pre_sb;}else{ ++pre_tp;}
              }
          }
        // name the LMS substring
        if(diff){ ++name; prev=pos; pre_len = len; pre_sb = sb; pre_eb = eb; }
        names[pos]=name-1;
    }
    // compact names array
    j = n-1;
    for(i=0; i<n; ++i){
        if(names[n-i-1]!=EMPTY){SA[j--]=names[n-i-1];}
    }
    // free names array
    free(names);
  }
  
  // s1 is done now
  uint_s *SA1=SA, *s1=SA+n-n1;  

  // stage 2: solve the reduced problem
  // recurse if names are not unique yet
  if(name<n1) {
      cSACA((uint_s*)s1, SA1, n1, name, sizeof(uint_s), level+1, nb_s);
  } else { // stop the recursion, generate the suffix array of s1 directly
      for(i=0; i<n1; i++){ SA1[s1[i]] = i; }
  }
  
  // stage 3: induce the result for the original problem
  bkt = (uint_s *)malloc(sizeof(uint_s)*K); // bucket counters
  
  // put all left-most S characters into their buckets
  j=n1-1;
  for(i=0; i<nseq; ++i){
      size_t ni = nseq-i; bool type;
      sb=s_s(ni); eb=s_s(ni+1)-1; len=eb-sb+1;
      if(len==1){ sn.push_back(sb); } // fill singletons vector 
      else{
          // if it is not a singleton find type of the first int
          type = 0;
          if(chr(sb)!=chr(sb+1)){ if(chr(sb)<chr(sb+1)){ type = 1; } }
          else{
              m = 1;
              while(m < len-1){
                  if(chr(sb+m)!=chr(sb+1+m)){ type=(chr(sb+m)>chr(sb+1+m))?0:1; break; }
                  ++m;
              }
          }
          // type of last int
          if(chr(eb)!=chr(sb)){ type = (chr(eb)>chr(sb))?0:1; }
          // find all S* suffixes
          for(m=0; m<len-1; ++m){
              size_t nx = eb-m-1;
              if(chr(nx)>chr(nx+1)){ if(type){ s1[j--]=nx+1; type=0; } }
              else{ 
                  if(chr(nx)<chr(nx+1)){ type=1; }
              }
          }
          // check if first suffix is S*
          if(chr(sb)<chr(eb)){ if(type){ s1[j--]=sb; } }
      }
  }
  
  // we must have the same number of S* of stage 1
  assert(j+1==0);
    
  for(i=0; i<n1; ++i) {SA1[i]=s1[SA1[i]];} // get index in s1
   
  if(n1>0){ // if at least one s* is not a singleton
      for(i=n1; i<n; ++i) SA[i]=EMPTY; // init SA[n1..n-1]
      getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  }
  
  // insert S* suffixes at the end of the buckets
  for(i=0; i < n1; ++i){
      j=SA[n1-i-1]; SA[n1-i-1]=EMPTY;
      SA[bkt[chr(j)]--]=j;
  }
  
  // induce the L and S suffixes to compute the final SA of each level
  induceL(SA, s, bkt, r_s, s_s, b_s, n, K, cs, false, level);
  induceS(SA, s, bkt, sn, r_s, s_s, b_s, n, K, cs, false, level); 
  
  // free onset and bucket vectors
  onset.clear();
  free(bkt); 
}

/**
 * @brief Compute the suffix array of a vector of integers using cSAIS algo
 * 
 * @param s input string
 * @param SA suffix array
 * @param n length of the input string
 * @param K alphabet size
 * @param b_s bitvector of the starting phrases of the parse
 */
void csais_int(uint_p *s, uint_s *SA, size_t n, size_t K, sd_vector<> &b_s){
    if((s == nullptr) || (SA == nullptr) || (n < 0)) {cerr << "Empty input given." << endl; exit(1);}
    cSAIS((uint_s *)s, (uint_s *)SA, n, K, sizeof(uint_p), 0, b_s);
}

/**
 * @brief Compute the suffix array of a vector of integers using cSACA algo
 * 
 * @param s input string
 * @param SA suffix array
 * @param n length of the input string
 * @param K alphabet size
 * @param b_s bitvector of the starting phrases of the parse
 */
void csaca_int(uint_p *s, uint_s *SA, size_t n, size_t K, sd_vector<> &b_s){
    if((s == nullptr) || (SA == nullptr) || (n < 0)) {cerr << "Empty input given." << endl; exit(1);}
    cSACA((uint_s *)s, (uint_s *)SA, n, K, sizeof(uint_p), 0, b_s);
}