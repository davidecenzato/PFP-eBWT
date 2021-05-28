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
// #include "malloc_count.h"

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
  
    //   std::cout << "Memory peak-eparse: " << malloc_count_peak() << std::endl;

    #if P64 == 0
        // if we are in 32 bit mode, check that parse has less than 2^32-2 words
        checkParseSize();
    #endif 
 
    // read starting positions
    tmp_filename = filename + std::string(".start");
    read_file(tmp_filename.c_str(), sts);
    // create bit vector for starting positions
    sts.push_back(size);
    sdsl::sd_vector_builder builder(size+1,sts.size());
    for(auto idx: sts){builder.set(idx);}
    b_d = sdsl::sd_vector<>(builder);
    build(saP_flag_, ilP_flag_);
    
    buildBitIl(); // build Inverted List

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

    // serialize parse data structures
    std::string output = filename + std::string(".sdsl");
    std::ofstream out(output);
    b_il.serialize(out);
    //select_ilist.serialize(out);
    b_st.serialize(out);
    //rank_st.serialize(out);
    my_serialize(ilP,out);
    my_serialize(offset,out);
    out.close();
    
   }
    
    void build(bool saP_flag_, bool ilP_flag_){
        size_t p_size = p.size();
        if(saP_flag_){
            saP.resize(p_size);
            // std::cout << "Current memory before csaca: " << malloc_count_current() << std::endl;
            // suffix array of the parse.
            verbose("Computing cSA of the parse");
            _elapsed_time(
                std::cout << "Starting computing cSA" << std::endl;
                // build SA using circular SAIS algorithm
                csais_int(&p[0],&saP[0], size, alphabet_size+1, b_d);
            );
        }
        if(ilP_flag_){
            ebwtP.resize(p_size);
            verbose("Computing Inverted List of the parse");
            _elapsed_time(
                makeEBWT();
                p.clear(); p.shrink_to_fit();
                ilP.resize(p_size);
                countingSort(ebwtP,ilP);
                );
        }else{
            p.clear(); p.shrink_to_fit();
        }
    }

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
    }
    
    void checkParseSize(){
        std::cout << "Checking parse size." << std::endl;
        if(p.size() > pow(2,32) - 1) {
            std::cerr << "Input containing more than 2^32 - 1 words" << std::endl;
            std::cerr << pow(2,32) - 1 << std::endl;
            std::cerr << p.size() << std::endl;
            std::cerr << "Please use 64 bit version." << std::endl;
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
        p.clear();
    }
         
    
    void countingSort(std::vector<uint_p> &vec, std::vector<uint_s> &out)
    {
        std::vector<uint_s> count(alphabet_size,0);
        for(size_t i=0; i<vec.size(); ++i)
        {   
            uint_p cs = vec[i];
            count[cs]++;
        }  
        uint_s prev_c = count[0];
        count[0] = 0;
        for(size_t i=1; i<count.size(); ++i)
        {
            uint_s tmp = count[i];
            count[i] = count[i-1] + prev_c;
            prev_c = tmp;
        }  
        for (size_t i = 0; i < vec.size(); ++i)
        {
            uint_p cs = vec[i];
            uint_s index = count[cs]++;
            out[index] = i;
        }
    }
    
};

#endif /* PARSE_HPP */

