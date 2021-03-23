/* 
 * File:   parse.hpp
 * Author: hejans
 *
 * Created on February 12, 2021, 11:11 AM
 */

#ifndef PARSE_HPP
#define PARSE_HPP


#include <sys/stat.h>
#include <algorithm>

#include "common.hpp"
#include "csais.h"

class parse{
private:
    std::vector<uint32_t> ebwtP;
    std::vector<uint32_t> p;
    std::vector<uint64_t> sts;
    std::vector<uint_t> saP;
    sdsl::bit_vector::rank_1_type rank_b_d;
    sdsl::bit_vector::select_1_type select_b_d;
    sdsl::bit_vector b_il;
    size_t size;
    size_t alphabet_size;
    
public:  
    std::vector<uint_t> ilP;
    std::vector<uint64_t> offset;
    sdsl::bit_vector::select_1_type select_ilist;
    sdsl::bit_vector::rank_1_type rank_st;
    sdsl::bit_vector b_d;
    sdsl::bit_vector b_st;
    
    bool saP_flag = false;
    bool ilP_flag = false;

  // Default constructor for load
    parse() {}
    
    parse(std::string filename,
          bool saP_flag_ = true,
          bool ilP_flag_ = true)//:
          //alphabet_size(alphabet_size_)
  {
    // read file
    std::string tmp_filename = filename + std::string(".parse");
    read_file(tmp_filename.c_str(), p);
    alphabet_size = *std::max_element(p.begin(),p.end());
    size = p.size();
  
    checkParseSize();
 
    // read starting positions
    tmp_filename = filename + std::string(".start");
    read_file(tmp_filename.c_str(), sts);
    // create bit vector for starting positions
    sdsl::bit_vector tmp(size+1,0); b_d = tmp; b_d[size]=1;
    for(int i=0;i<sts.size();i++){b_d[sts[i]]=1;}
    rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
    select_b_d = sdsl::bit_vector::select_1_type(&b_d);
    
    std::vector<uint64_t> temp(offset.size(),0);
    tmp_filename = filename + std::string(".offset");
    read_file(tmp_filename.c_str(), temp);
    
    build(saP_flag_, ilP_flag_);
    
    buildBitIl();
    
    // build vector of starting character offsets
    offset.resize(temp.size()); int j=0;
    for(int i=0;i<ilP.size();i++){
        if(b_st[i]==1){
            offset[j] = temp[rank_b_d(saP[ilP[i]]+1)-1]; j++;
        }
    }
    temp.clear();
    
    clearVectors();
    
   }
    
    void build(bool saP_flag_, bool ilP_flag_){
        
        if(saP_flag_){
            saP.resize(p.size());
            // suffix array of the parse.
            verbose("Computing cSA of the parse");
            _elapsed_time(
                // build SA using circular SA-IS algorithm
                cSAIS(&p[0],&saP[0], size, alphabet_size+1, sizeof(int), 0, b_d);
            );
        }
        
        if(ilP_flag_){
            ebwtP.resize(p.size());
            ilP.resize(p.size());
            verbose("Computing Inverted List of the parse");
            _elapsed_time(
                makeEBWT();
                countingSort(ebwtP,ilP);
                );
        }
    }
    
    void buildBitIl(){
        
        // initialize bit vector of the inverted list
        // initialize bit vector of the words at the end of parse phrases
        b_il.resize(ilP.size()+1);
        b_st.resize(ilP.size());
        b_il[0]=1; b_il[ilP.size()]=1;
        size_t prev = ebwtP[ilP[0]];
        if(b_d[saP[ilP[0]]]==1){b_st[0]=1;}else{b_st[0]=0;}
        for(int i=1;i<ilP.size();i++)
        {   
            if(b_d[saP[ilP[i]]]==1){b_st[i]=1;}else{b_st[i]=0;}
            size_t next = ebwtP[ilP[i]];
            if(next!=prev){b_il[i]=1;}
            else{b_il[i]=0;}
            prev = next;
        }
        rank_st = sdsl::bit_vector::rank_1_type(&b_st);
        select_ilist = sdsl::bit_vector::select_1_type(&b_il);
    }
    
    void checkParseSize(){
        if(p.size() > 0x7FFFFFFE) {
            die("Input containing more than 2^31-2 phrases!\n");
            die("Please use 64 bit version\n");
            exit(1); }
    }
    
    void clearVectors(){
        //clear vectors not useful for the next step
        p.clear();
        sts.clear();
        saP.clear();
        ebwtP.clear();
    }
    
    void makeEBWT(){
     // build the eBWT of the parsing using circular SA
     uint32_t pc = 0;
     for (int i=0; i<saP.size();i++){
         if(b_d[saP[i]]==1){ 
             pc = p[select_b_d(rank_b_d(saP[i]+1)+1)-1]-1;
             ebwtP[i] = pc;
         }else{
             pc = p[saP[i]-1]-1;
             ebwtP[i] = pc;
             }
         }
    }
         
    
    void countingSort(std::vector<uint32_t> vec, std::vector<uint_t> &out)
    {
        std::vector<uint32_t> count(alphabet_size,0);
        for(size_t i=0; i<vec.size(); ++i)
        {   
            uint32_t cs = vec[i];
            count[cs]++;
        }  
        std::vector<uint32_t> psum(count.size(),0);
        for(size_t i=1; i<count.size(); ++i)
        {
            psum[i] = psum[i-1] + count[i-1];
        }  
        count.clear();
        for (size_t i = 0; i < vec.size(); ++i)
        {
            uint32_t cs = vec[i];
            size_t index = psum[cs]++;
            out[index] = i;
        }
    }
    
};



#endif /* PARSE_HPP */

