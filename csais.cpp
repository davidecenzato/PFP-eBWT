/*
 * Circular SAIS implementation to compute the circular SA of an integer vector.
 * 
 * This code is adapted from https://github.com/kurpicz/saca-bench/blob/master/sa-is/sais.cpp
 * which is the original code of the SA-IS algorithm listed below
 */
// This is the sample code for the SA-IS algorithm presented in
// our article "Two Efficient Algorithms for Linear Suffix Array Construction"
// (co-authored with Sen Zhang and Wai Hong Chan),
// which can be retrieved at: http://code.google.com/p/ge-nong/


#include "csais.h"

using namespace std;
using namespace sdsl;

const int EMPTY=0xffffffff;
unsigned char mask[]={0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};
#define tget(i) ( (t[(i)/8]&mask[(i)%8]) ? 1 : 0 )
#define tset(i, b) t[(i)/8]=(b) ? (mask[(i)%8]|t[(i)/8]) : ((~mask[(i)%8])&t[(i)/8])

#define chr(i) (cs==sizeof(uint_s)?((uint_s *)s)[i]:((uint_p *)s)[i])
#define isLMS(i,j) ((tget(i) && !tget(j)))

// compute the head or end of each bucket
void getBuckets(uint_s *s, uint_s *bkt, size_t n, size_t K, size_t cs, bool end) { 
  size_t i, sum=0;
  for(i=0; i<K; i++) bkt[i]=0; // clear all buckets
  for(i=0; i<n; i++) bkt[chr(i)]++; // compute the size of each bucket
  for(i=0; i<K; i++) { sum+=bkt[i]; bkt[i]= end ? sum-1 : sum-bkt[i]; }
}

// induce L suffixes
void induceSAl(unsigned char *t, uint_s *SA, uint_s *s, uint_s *bkt, uint_s *l_bkt, bit_vector::rank_1_type& r_s,
               bit_vector::select_1_type& s_s, bit_vector& b_s, size_t n, size_t K, size_t cs, int level, vector<uint_s>& star) { 
    
    size_t i, j, rank;
    getBuckets(s, bkt, n, K, cs, false); // find heads of buckets
    
    for(i=0;i<star.size();i++){
      j=star[i]; SA[(l_bkt[chr(j)])+bkt[chr(j)]++]=j;
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

// induce S suffixes
void induceSAs(unsigned char *t, uint_s *SA, uint_s *s, uint_s *bkt, bit_vector::rank_1_type& r_s,
               bit_vector::select_1_type& s_s, bit_vector& b_s, size_t n, size_t K, size_t cs, int level) { 
    
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

// find the circular suffix array SA of s[0..n-1] in {0..K-1}^n
// level starts from 0
/** @brief computes the circular suffix array of string s[0..n-1] 
 *
 *  @param s	 input string 
 *  @param SA    suffix array 
 *  @param n	 string length
 *  @param K	 alphabet size
 *  @param cs	 integer size
 *  @param level recursion level, debug only
 *  @param b_s   starting positions bit vector
 *  @return      None
 */
void cSAIS(uint_s *s, uint_s *SA, size_t n, size_t K, size_t cs, int level, bit_vector &b_s) {
  size_t i, j, nseq, rank;
  size_t sb, eb, fm, len; // maybe uint_s
  bit_vector::rank_1_type r_s = bit_vector::rank_1_type(&b_s);
  bit_vector::select_1_type s_s = bit_vector::select_1_type(&b_s);
  
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
          bool type = (chr(pos)<chr(prev))?1:0;
          tset(pos, type); if(type==0) l_bkt[chr(pos)]++;
          while(pos != fm){
              if(pos-sb > 0){ pos--; }
              else {pos = eb;}
              if(prev-sb > 0){ prev--; }
              else {prev = eb;}
              
              type=(chr(pos)<chr(prev) || (chr(pos)==chr(prev) && tget(prev)==1));
              if(type==0){l_bkt[chr(pos)]++;}
              tset(pos,type?1:0);
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
  bit_vector nb_s(n1+1,0); size_t sum=0;
  for(i=0;i<(nseq+1);i++){ sum+=sts[i]; nb_s[sum]=1;}
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
  
  free(bkt); 
  free(l_bkt);
  free(t);
  
}

/**
 * @brief Compute the suffix array 
 * 
 * @param s input string
 * @param SA suffix array
 * @param n length of the input string
 * @param K alphabet size
 * @param b_s bitvector of the starting phrases of the parse
 */
void csais_int(uint_p *s, uint_s *SA, size_t n, size_t K, bit_vector &b_s){
    if((s == nullptr) || (SA == nullptr) || (n < 0)) {cerr << "Empty input given." << endl; exit(1);}
    cSAIS((uint_s *)s, (uint_s *)SA, n, K, sizeof(uint_p), 0, b_s);
}