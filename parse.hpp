/*
 * Code to build the inverted list of the eBWT of a circular Prefix Free Parse.
 * 
 * This code is adapted from https://github.com/maxrossi91/pfp-thresholds/blob/master/include/pfp/parse.hpp
 */

#ifndef PARSE_HPP
#define PARSE_HPP


#include <sys/stat.h>
#include <algorithm>

#include "common.hpp"
#include "csais.h"

class parse{
private:
    std::vector<uint_p> ebwtP;
    std::vector<uint_p> p;
    std::vector<uint64_t> sts;
    std::vector<uint_s> saP;
    sdsl::sd_vector<> b_d;
    sdsl::sd_vector<>::rank_1_type rank_b_d;
    sdsl::sd_vector<>::select_1_type select_b_d;
    size_t size;
    size_t alphabet_size;
    
public:  
    std::vector<uint_s> ilP;
    std::vector<uint32_t> offset;
    sdsl::sd_vector<> b_il;
    sdsl::sd_vector<> b_st;
    sdsl::sd_vector<>::select_1_type select_ilist;
    sdsl::sd_vector<>::rank_1_type rank_st;
    
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
    std::string tmp_filename = filename + std::string(".eparse");
    read_file(tmp_filename.c_str(), p);
    alphabet_size = *std::max_element(p.begin(),p.end());
    size = p.size();
  
    #if P64 == 0
        // if we are in 32 bit mode, check that parse has less than 2^32-2 words
        checkParseSize();
    #endif 
 
    // read starting positions
    tmp_filename = filename + std::string(".start");
    read_file(tmp_filename.c_str(), sts);
    // create bit vector for starting positions
    ////std::vector<size_t> onset;
    ////for(size_t i=0;i<sts.size();i++){onset.push_back(sts[i]);}
    ////onset.push_back(size);
    sts.push_back(size);
    sdsl::sd_vector_builder builder(size+1,sts.size());
    for(auto idx: sts){builder.set(idx);}
    b_d = sdsl::sd_vector<>(builder);
    //sts.clear();
    
    build(saP_flag_, ilP_flag_);
    
    buildBitIl();

    std::vector<uint32_t> temp(offset.size(),0);
    tmp_filename = filename + std::string(".offset");
    read_file(tmp_filename.c_str(), temp);
    
    // build vector of starting character offsets
    offset.resize(temp.size()); size_t j=0;
    for(size_t i=0;i<ilP.size();++i){
        if(b_st[i]==1){
            offset[j] = temp[rank_b_d(saP[ilP[i]]+1)-1]; ++j;
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
                std::cout << "Starting computing cSA" << std::endl;
                csais_int(&p[0],&saP[0], size, alphabet_size+1, b_d);
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
    
    /*
    void buildBitIl(){
        
        // initialize bit vector of the inverted list
        // initialize bit vector of the words at the end of parse phrases
        b_il.resize(ilP.size()+1);
        b_st.resize(ilP.size());
        b_il[0]=1; b_il[ilP.size()]=1;
        uint_p prev = ebwtP[ilP[0]];
        if(b_d[saP[ilP[0]]]==1){b_st[0]=1;}else{b_st[0]=0;}
        for(size_t i=1;i<ilP.size();i++)
        {   
            if(b_d[saP[ilP[i]]]==1){b_st[i]=1;}else{b_st[i]=0;}
            uint_p next = ebwtP[ilP[i]];
            if(next!=prev){b_il[i]=1;}
            else{b_il[i]=0;}
            prev = next;
        }
        rank_st = sdsl::bit_vector::rank_1_type(&b_st);
        select_ilist = sdsl::bit_vector::select_1_type(&b_il);
    }*/

    void buildBitIl(){
        
        // initialize bit vector of the inverted list
        // initialize bit vector of the words at the end of parse phrases
        std::vector<size_t> onset_il; onset_il.push_back(0); 
        std::vector<size_t> onset_st;
        uint_p prev = ebwtP[ilP[0]];
        if(b_d[saP[ilP[0]]]==1){ onset_st.push_back(0); }
        for(size_t i=1;i<ilP.size();i++)
        {   
            if(b_d[saP[ilP[i]]]==1){ onset_st.push_back(i); }
            uint_p next = ebwtP[ilP[i]];
            if(next!=prev){ onset_il.push_back(i); }
            prev = next;
        }
        onset_il.push_back(ilP.size());
        // initalize bit vectors
        sdsl::sd_vector_builder builder_il(ilP.size()+1,onset_il.size());
        for(auto idx: onset_il){builder_il.set(idx);}
        b_il = sdsl::sd_vector<>(builder_il);
        sdsl::sd_vector_builder builder_st(ilP.size(),onset_st.size());
        for(auto idx: onset_st){builder_st.set(idx);}
        b_st = sdsl::sd_vector<>(builder_st);
        // initialize support for rank and select
        rank_st = sdsl::sd_vector<>::rank_1_type(&b_st);
        select_ilist = sdsl::sd_vector<>::select_1_type(&b_il);

    }
    
    void checkParseSize(){
        std::cout << "Checking parse size." << std::endl;
        if(p.size() > pow(2,32) - 1) {
            std::cerr << "Input containing more than 2^32 - 1 words" << std::endl;
            std::cout << pow(2,32) - 1 << std::endl;
            std::cout << p.size() << std::endl;
            die("Input containing more than 2^32-1 phrases!\n");
            die("Please use 64 bit version\n");
            exit(1); }
    }
    
    void clearVectors(){
        //clear vectors not useful for the next step
        sts.clear();
        saP.clear();
        ebwtP.clear();
    }
    
    void makeEBWT(){
        // build support for rank and select
        rank_b_d = sdsl::sd_vector<>::rank_1_type(&b_d);
        select_b_d = sdsl::sd_vector<>::select_1_type(&b_d);
        // build the eBWT of the parsing using circular SA
        uint_p pc = 0;
        for (size_t i=0; i<saP.size();i++){
            if(b_d[saP[i]]==1){ 
                pc = p[select_b_d(rank_b_d(saP[i]+1)+1)-1]-1;
                ebwtP[i] = pc;
            }else{
                pc = p[saP[i]-1]-1;
                ebwtP[i] = pc;
            }  
        }
        //rank_b_d.~rank_support();
        //select_b_d.~select_support();
        p.clear();
    }
         
    
    void countingSort(std::vector<uint_p> vec, std::vector<uint_s> &out)
    {
        std::vector<uint_s> count(alphabet_size,0);
        for(size_t i=0; i<vec.size(); ++i)
        {   
            uint_p cs = vec[i];
            count[cs]++;
        }  
        std::vector<uint_s> psum(count.size(),0);
        for(size_t i=1; i<count.size(); ++i)
        {
            psum[i] = psum[i-1] + count[i-1];
        }  
        count.clear();
        for (size_t i = 0; i < vec.size(); ++i)
        {
            uint_p cs = vec[i];
            uint_s index = psum[cs]++;
            out[index] = i;
        }
    }
    
};



#endif /* PARSE_HPP */

