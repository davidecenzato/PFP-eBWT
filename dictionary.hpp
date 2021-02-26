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

// TODO: Extend it to integer alphabets
class dictionary{
public:
  std::vector<uint8_t> d;
  std::vector<uint32_t> occ;
  std::vector<uint_t> saD;
  std::vector<int_t> lcpD;
  sdsl::bit_vector b_d; // Starting position of each phrase in D  
  sdsl::bit_vector::rank_1_type rank_b_d;
  sdsl::bit_vector::select_1_type select_b_d;

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
    //b_d.resize(d.size(),0);
    b_d[0] = 1;
    size_t numdol=0;
    for(int i=1;i<d.size()-1;i++)
    {
        if(d[i]==1){b_d[i+1]=1;numdol++;}
    }

    rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
    select_b_d = sdsl::bit_vector::select_1_type(&b_d);
    // reading occurrences file

    tmp_filename = filename + std::string(".occ");
    read_file(tmp_filename.c_str(), occ);
    
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
            //for(int i=0;i<saD.size();i++){std::cout << "saD " << int(saD[i]) << std::endl;}
            //for(int i=0;i<lcpD.size();i++){std::cout << "lcpD " << int(lcpD[i]) << std::endl;}
        );
    }
};
  

#endif /* DICTIONARY_HPP */
