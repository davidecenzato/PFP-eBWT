/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   parse.hpp
 * Author: hejans
 *
 * Created on February 12, 2021, 11:11 AM
 */

#ifndef PARSE_HPP
#define PARSE_HPP

#include "common.hpp"

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sys/stat.h>
#include <algorithm>
#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>

typedef std::tuple<int,int,int> triplet;

class parse{
private:
    std::vector<uint32_t> p;
    std::vector<uint64_t> sts;
    std::vector<uint32_t> slen;
    std::vector<uint32_t> saP;
    sdsl::bit_vector::rank_1_type rank_b_d;
    sdsl::bit_vector b_il;
    size_t size;
    size_t alphabet_size;
    
public:  
    std::vector<uint32_t> ilP;
    std::vector<uint32_t> ebwtP;
    sdsl::bit_vector::select_1_type select_ilist;
    sdsl::bit_vector b_d;
    
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
    sdsl::bit_vector b_d(size,0);
    for(uint64_t i=0;i<sts.size();i++){b_d[sts[i]]=1;}
    rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
    
    // read lengths
    tmp_filename = filename + std::string(".len");
    read_file(tmp_filename.c_str(), slen);
    
    build(saP_flag_, ilP_flag_);
    
    buildBitIl();
    
    clearVectors();
    
   }
    
    void build(bool saP_flag_, bool ilP_flag_){
        
        if(saP_flag_){
            saP.resize(p.size());
            // suffix array of the parse.
            verbose("Computing cSA of the parse");
            _elapsed_time(
                buildSA();
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
        b_il.resize(ilP.size()+1);
        b_il[0]=1; b_il[ilP.size()]=1;
        size_t prev = ebwtP[ilP[0]];
        for(int i=1;i<ilP.size();i++)
        {   
            size_t next = ebwtP[ilP[i]];
            if(next!=prev){b_il[i]=1;}
            else{b_il[i]=0;}
            prev = next;
        }
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
        slen.clear();
        saP.clear();
        ebwtP.clear();
    }
   
    void buildSA(){

        std::queue<triplet> buckets;
        //uint32_t size = std::accumulate(slen.begin(),slen.end(),0);
        buckets.push({0,size,0});
        std::iota(saP.begin(),saP.end(),0);
        size_t max_len = *std::max_element(slen.begin(),slen.end());

        while(!buckets.empty()){
           auto bucket = buckets.front(); buckets.pop();
           int start = std::get<0>(bucket);
           int end   = std::get<1>(bucket);
           int depth = std::get<2>(bucket);

           if((start < end)){
                std::vector<uint32_t> count;
                std::map <uint32_t,uint32_t> count_map;
                for(size_t i=start; i<end; ++i)
                {   
                    uint32_t cs = rank_b_d(saP[i]+1)-1;
                    size_t ind = sts[cs]+((saP[i]-sts[cs]+depth)%(slen[cs]));
                    uint32_t symb = p[ind]-1;
                    if(count_map.find(symb)==count_map.end()){
                        count_map.insert(std::pair<uint32_t,uint32_t>(symb,1));
                    }else{
                        count_map.at(symb)++;
                    }
                }
                std::map <uint32_t,uint32_t> psum_map;
                size_t pind = 0;
                for(auto& x: count_map){
                    count.push_back(x.second);
                    psum_map.insert(std::pair<uint32_t,uint32_t>(x.first,pind));
                    pind++;
                }
                count_map.clear();
                std::vector<uint32_t> psum(count.size(),0);
                for(size_t i=1; i<count.size(); ++i)
                {
                    psum[i] = psum[i-1] + count[i-1];
                }   
                std::vector<uint32_t> tmp(end - start, 0);
                for (size_t i = start; i < end; ++i)
                {
                    uint32_t cs = rank_b_d(saP[i]+1)-1;
                    size_t ind = sts[cs] + ((saP[i]-sts[cs]+depth)%(slen[cs]));
                    size_t index = psum[psum_map[p[ind]-1]]++;
                    tmp[index] = saP[i];
                }
                for (uint32_t i=start; i<tmp.size()+start;i++){
                    saP[i] = tmp[i-start];
                }
                tmp.clear();
                for(size_t i=0; i<count.size(); i++)
                {
                    end = start + count[i];
                    if(end - start > 1 && depth < (max_len*2-2)){buckets.push({start,end,depth+1});}
                    start = end;
                }     
            }
        }  
    }
   
    void makeEBWT(){
     // build the eBWT of the parsing using pSA
     size_t bchar = 0;
     for (int i=0; i<saP.size();i++){
         if(saP[i]==0 || rank_b_d(saP[i]+1) != rank_b_d(saP[i])){ 
             bchar = p[(saP[i]-1)+slen[rank_b_d(saP[i]+1)-1]]-1;
             ebwtP[i] = bchar;
         }else{
             bchar = p[saP[i]-1]-1;
             ebwtP[i] = bchar;
             }
         }
    }
    
    void countingSort(std::vector<uint32_t> vec, std::vector<uint32_t> &out)
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
    
  // Serialize to a stream.
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += my_serialize(p, out, child, "parse");
    written_bytes += my_serialize(saP, out, child, "saP");
    written_bytes += my_serialize(ilP, out, child, "ilP");
    //written_bytes += my_serialize(lcpP, out, child, "lcpP");
    //written_bytes += rmq_lcp_P.serialize(out, child, "rmq_lcp_P");
    // written_bytes += b_p.serialize(out, child, "b_p");
    // written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
    // written_bytes += select_b_p.serialize(out, child, "select_b_p");
    written_bytes += sdsl::write_member(alphabet_size, out, child, "alphabet_size");
    // written_bytes += sdsl::serialize(p, out, child, "parse");
    // written_bytes += sdsl::serialize(saP, out, child, "saP");
    // written_bytes += sdsl::serialize(isaP, out, child, "isaP");
    // written_bytes += sdsl::serialize(lcpP, out, child, "lcpP");
    // written_bytes += rmq_lcp_P.serialize(out, child, "rmq_lcp_P");
    // // written_bytes += b_p.serialize(out, child, "b_p");
    // // written_bytes += rank_b_p.serialize(out, child, "rank_b_p");
    // // written_bytes += select_b_p.serialize(out, child, "select_b_p");
    // written_bytes += sdsl::write_member(alphabet_size, out, child, "alphabet_size");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream &in)
  {
    my_load(p, in);
    my_load(saP, in);
    my_load(ilP, in);
    //my_load(lcpP, in);
    //rmq_lcp_P.load(in);
    // b_p.load(in);
    // rank_b_p.load(in);
    // select_b_p.load(in);
    sdsl::read_member(alphabet_size, in);
    // sdsl::load(p, in);
    // sdsl::load(saP, in);
    // sdsl::load(isaP, in);
    // sdsl::load(lcpP, in);
    // rmq_lcp_P.load(in);
    // // b_p.load(in);
    // // rank_b_p.load(in);
    // // select_b_p.load(in);
    // sdsl::read_member(alphabet_size, in);
  }

};



#endif /* PARSE_HPP */

