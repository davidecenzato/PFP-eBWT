/*
 * Circular SAIS implementation to compute the circular SA of a integer vector.
 */

#include "csais.h"

using namespace std;
using namespace sdsl;

const int EMPTY=0xffffffff;
unsigned char mask[]={0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};
#define tget(i) ( (t[(i)/8]&mask[(i)%8]) ? 1 : 0 )
#define tset(i, b) t[(i)/8]=(b) ? (mask[(i)%8]|t[(i)/8]) : ((~mask[(i)%8])&t[(i)/8])

#define chr(i) (cs==sizeof(char)?((uint8_t *)s)[i]:((unsigned int *)s)[i])
#define isLMSc(i,j) ((tget(i) && !tget(j)))

// compute the head or end of each bucket
void getBuckets(unsigned int *s, int *bkt, int n, int K, int cs, bool end) { 
  int i, sum=0;
  for(i=0; i<K; i++) bkt[i]=0; // clear all buckets
  for(i=0; i<n; i++) bkt[chr(i)]++; // compute the size of each bucket
  for(i=0; i<K; i++) { sum+=bkt[i]; bkt[i]= end ? sum-1 : sum-bkt[i]; }
}

// induce L suffixes
void induceSAl(unsigned char *t, uint_t *SA, unsigned int *s, int *bkt, int *l_bkt, bit_vector::rank_1_type& r_s,
               bit_vector::select_1_type& s_s, bit_vector& b_s, int n, int K, int cs, int level, vector<uint32_t>& star) { 
    
    int i, j, rank;
    getBuckets(s, bkt, n, K, cs, false); // find heads of buckets
    
    for(int i=0;i<star.size();i++){
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
void induceSAs(unsigned char *t, uint_t *SA, unsigned int *s, int *bkt, bit_vector::rank_1_type& r_s,
               bit_vector::select_1_type& s_s, bit_vector& b_s, int n, int K, int cs, int level) { 
    
  int i, j, rank;
  getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  for(i=n-1; i>=0; i--){
    if(SA[i]!=EMPTY) {
        j=SA[i];
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
void cSAIS(unsigned int *s, uint_t *SA, int n, int K, int cs, int level, bit_vector &b_s) {
  int i, j, nseq, rank, sb, eb, fm, len;
  bit_vector::rank_1_type r_s = bit_vector::rank_1_type(&b_s);
  bit_vector::select_1_type s_s = bit_vector::select_1_type(&b_s);
  nseq = r_s(n);
  unsigned char *t=(unsigned char *)malloc(n/8+1); // LS-type array in bits
  
  // stage 1: reduce the problem by at least 1/2
  
  int *bkt = (int *)malloc(sizeof(int)*K); // bucket counters
  
  int *l_bkt = (int *)malloc(sizeof(int)*K); for(i=0; i<K; i++) l_bkt[i]=0; // l types bucket counters
  
  for(i=0; i<n; i++) SA[i]=EMPTY;
  
  vector<uint32_t> sgt; // singletons vector
  for(int i=0;i<nseq;i++){
      sb = s_s(i+1), eb = s_s(i+2)-1, fm = eb+1, len=eb-sb+1;
      if(len==1){tset(sb,1);sgt.push_back(sb);}
      else{
        // find first type
        for(int j=eb;j>sb;j--){
            if(chr(j)!=chr(j-1)){
              if(chr(j)>chr(j-1)){fm=j;tset(j,0);break;}
              else{fm=j;tset(j,1);break;}
            }else{tset(j,0);} 
        }
        // mismatch not found!
        if(fm==eb+1){cerr << "Error! Sequence without a mismatch detected." << endl; exit(1);}
        else{
          // Classify the type of each character
          int pos=fm-1, prev=fm;
          bool type = (chr(pos)<chr(prev))?1:0;
          tset(pos, type); if(type==0) l_bkt[chr(pos)]++;
          while(pos != fm){
              pos=sb+((len+((pos-1-sb)%len))%len); 
              prev=sb+((len+((prev-1-sb)%len))%len); 
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
  for(i=n-1; i>=0; i--)
  {
      if(b_s[i]==0){
          if(isLMSc(i,i-1)) SA[bkt[chr(i)]--]=i;
      }else{
          rank=r_s(i+1), eb=s_s(rank+1)-1;
          if(isLMSc(i,eb)) SA[bkt[chr(i)]--]=i;
      }
  }
   
  induceSAl(t, SA, s, bkt, l_bkt, r_s, s_s, b_s, n, K, cs, level, sgt);  
  induceSAs(t, SA, s, bkt, r_s, s_s, b_s, n, K, cs, level); 

  free(bkt);
  free(l_bkt);
  sgt.clear();
  
  // compact all the sorted substrings into the first n1 items of s
  int n1=0;
  vector<uint_t> sts(nseq+1,0);
  for(int i=0;i<n;i++){
      int pos=SA[i];
      rank=r_s(pos+1),sb=s_s(rank),eb=s_s(rank+1);
      if(pos==sb){
          if(isLMSc(sb,eb-1)){SA[n1++]=pos;sts[rank]++;}
      }else{
          if(isLMSc(pos,pos-1)){SA[n1++]=pos;sts[rank]++;}
      }
  }
  
  //create the bit vector of the new string
  bit_vector nb_s(n1+1,0); int sum=0;
  for(int i=0;i<sts.size();i++){ sum+=sts[i]; nb_s[sum]=1;}
  sts.clear();
  
  // Init the name array buffer
  for(i=n1; i<n; i++) SA[i]=EMPTY;
  
  // find the lexicographic names of all substrings
  // insert the first substring
  int name=1, prev=0; int pos = prev = SA[0]; vector<uint_t> names(n,EMPTY);
  names[pos]=name-1;
  for(i=1; i<n1; i++) {
        int pp=0, pv=0;
        pos=SA[i]; bool diff=false;
        int rank=r_s(pos+1), rank_p=r_s(prev+1);
        int curr=s_s(rank), len=s_s(rank+1)-curr; 
        int curr_p=s_s(rank_p), len_p=s_s(rank_p+1)-curr_p;

        for(int d=0; d<n; d++){
            //compute cyclic index
            int j=curr+((pos-curr+d)%(len));  
            int j_p=curr_p+((prev-curr_p+d)%(len_p)); 
            if(prev==-1 || chr(j)!=chr(j_p) || tget(j)!=tget(j_p))
            { 
                diff=true; break;
            }
            else{
              if(d>0 && (isLMSc(j,pp) || isLMSc(j_p,pv)))
                break;
            }
            pp=j,pv=j_p;
        }

        if(diff){ name++; prev=pos; }
        names[pos]=name-1;
    }
  
  j = n-1;
  for(int i=names.size()-1;i>-1;i--){
      if(names[i]!=EMPTY){SA[j--]=names[i];}
  }
  names.clear();
    
  // s1 is done now
  uint_t *SA1=SA, *s1=SA+n-n1;  
  
  // stage 2: solve the reduced problem
  // recurse if names are not yet unique
  if(name<n1) {
      cSAIS((unsigned int*)s1, SA1, n1, name, sizeof(int), level+1, nb_s);
  } else { // stop the recursion, generate the suffix array of s1 directly
      for(i=0; i<n1; i++){ SA1[s1[i]] = i; }
  }
    
  // stage 3: induce the result for the original problem

  bkt = (int *)malloc(sizeof(int)*K); // bucket counters
  l_bkt = (int *)malloc(sizeof(int)*K); for(i=0; i<K; i++) l_bkt[i]=0;

  // put all left-most S characters into their buckets
  j=0;
  for(i=0;i<nseq;i++){
      sb=s_s(i+1), eb=s_s(i+2), len= eb-sb;
      if(len==1){sgt.push_back(sb);}
      else
      {
          if(isLMSc(sb,eb-1)){ s1[j++]=sb; }
          else{ if(tget(sb)==0){ l_bkt[chr(sb)]++; } }
          for(int x=sb+1;x<eb;x++){
              if(isLMSc(x,x-1)){ s1[j++]=x; }
              else{ if(tget(x)==0) {l_bkt[chr(x)]++; } } 
          }
      }
  }
    
  for(int i=0; i<n1; i++) {SA1[i]=s1[SA1[i]];} // get index in s1
   
  if(n1>0) for(i=n1; i<n; i++) SA[i]=EMPTY; // init SA[n1..n-1]
    
  getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
  
  for(i=n1-1; i>=0; i--) {
      j=SA[i]; SA[i]=EMPTY;
      SA[bkt[chr(j)]--]=j;}

  induceSAl(t, SA, s, bkt, l_bkt, r_s, s_s, b_s, n, K, cs, level, sgt); 
  induceSAs(t, SA, s, bkt, r_s, s_s, b_s, n, K, cs, level); 

  free(bkt); 
  free(l_bkt);
  free(t);
  sgt.clear();
  
}