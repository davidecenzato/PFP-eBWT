/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dictionary.hpp
 * Author: hejans
 *
 * Created on February 20, 2021, 1:32 AM
 */

#ifndef DICTIONARY_HPP
#define DICTIONARY_HPP

extern "C" {
    #include "gsa/gsacak.h"
}
#include <vector>
#include <sdsl/bit_vectors.hpp>

// TODO: Extend it to integer alphabets
class dictionary{
private:
  std::vector<uint64_t> fchar;
  
public:
  std::vector<uint8_t> d;
  std::vector<uint32_t> occ;
  std::vector<uint_t> saD;
  std::vector<int_t> lcpD;
  sdsl::bit_vector b_d; // Starting position of each phrase in D  
  sdsl::bit_vector b_s;
  sdsl::bit_vector::rank_1_type rank_b_d;
  sdsl::bit_vector::select_1_type select_b_d;
  sdsl::bit_vector::rank_1_type rank_b_s;

  // default constructor for load.
  dictionary() {}

  dictionary(std::string filename,
             size_t w)
  {
    // Building dictionary from file
    std::string tmp_filename = filename + std::string(".dict");
    read_file(tmp_filename.c_str(), d);
    // Creating bit vector for concatenated dictionary
    b_d = sdsl::bit_vector(d.size(),0);
    b_d[0] = 1;
    for(int i=1;i<d.size()-1;i++)
    {
        if(d[i]==1){b_d[i+1]=1;}
    }

    rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
    select_b_d = sdsl::bit_vector::select_1_type(&b_d);
    // reading occurrences file

    tmp_filename = filename + std::string(".occ");
    read_file(tmp_filename.c_str(), occ);
    
    tmp_filename = filename + std::string(".fchar");
    read_file(tmp_filename.c_str(), fchar);
    
    b_s = sdsl::bit_vector(d.size(),0);
    for(int i=0;i<fchar.size();i+=2){
        b_s[select_b_d(fchar[i])+fchar[i+1]]=1;
    }
    fchar.clear();
    rank_b_s = sdsl::bit_vector::rank_1_type(&b_s);
    
    // build data structures for the dictionary.
    build();
    
    }
  
    void build(){
        
        // initialize data structures to support rank and select on the bitvector
        // of the concatenated dictionary
        
        // resize SA and LCP arrays of the dictionary
        saD.resize(d.size());
        lcpD.resize(d.size());
        verbose("Computing SA and LCP of the dictionary");
        _elapsed_time(
            gsacak(&d[0], &saD[0], &lcpD[0], nullptr, d.size());
        );
    }
};
  

#endif /* DICTIONARY_HPP */

