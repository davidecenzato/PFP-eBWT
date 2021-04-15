/* ******************************************************************
 * invertebwt
 * Simple function that given the eBWT of text and its starting positions S, invert the eBWT to the original text.
 * ****************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <assert.h>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <sdsl/wavelet_trees.hpp>
#include <algorithm>

#include "common.hpp"

using namespace std;

void inverteBWT(std::vector<uint8_t> EBWT, std::vector<uint64_t> ST_P, string I, int alph_size){
    // Function that inverse the eBWT of a int string using sdsl wavelet trees
    string wt_filename = I + std::string(".ebwt");
    cout << "Building the Wavelet tree of the eBWT." << endl;
    sdsl::wt_blcd<> wt; sdsl::construct(wt,wt_filename,1);
    vector<uint32_t> C(alph_size+1,0); vector<uint32_t> bkt(alph_size,0);
    cout << "Building the C vector of the eBWT." << endl;
    for(int i=0;i<EBWT.size();i++){ bkt[EBWT[i]]++; }
    for(int i=1;i<C.size();i++){ C[i]=C[i-1]+bkt[i-1]; }
    
    string tmp_filename = I + std::string(".rfasta");
    FILE* fp = fopen(tmp_filename.c_str(), "w+");

    cout << "Inverting " << ST_P.size() << " sequences." << endl;
    for(size_t i=0;i<ST_P.size();++i){
        std::vector<uint8_t> RP;  
        int index = ST_P[i]; 
        uint8_t p = EBWT[index]; 
        RP.push_back(p);
        int starting = index;
        index = C[p]+wt.rank(index,p);
        while(index != starting){
            p = EBWT[index];
            RP.push_back(p);
            index = C[p]+wt.rank(index,p);
        }
        // write the sequence in the correct order
        reverse(RP.begin(),RP.end()); RP.push_back('\n');
        if((fwrite(&RP[0], sizeof(uint8_t), RP.size(), fp))!=RP.size()) {cerr << "fwrite failed" << endl;}
    }
    fclose(fp);
}

int main(int argc, char *argv[])
{
  vector<uint8_t> eBWT;
  vector<uint64_t> I;
  
  time_t start_wc = time(NULL);
  
    // check input data
  if(argc<2){
    printf("\nUsage:\n\t %s name.fasta\n\n", argv[0]);
    puts("Invert the eBWT of text name and output it to name.rfasta");
    exit(1);
  }
  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");
  
  string tmp_filename = argv[1] + string(".ebwt");
  read_file(tmp_filename.c_str(), eBWT);
  
  tmp_filename = argv[1] + std::string(".I");
  read_file(tmp_filename.c_str(), I);
  
  inverteBWT(eBWT,I,argv[1],256);
  
  return 0;
}