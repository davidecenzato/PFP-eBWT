/* pfp-parse - prefix free parsing parse
    Copyright (C) 2020 Massimiliano Rossi
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file parse.hpp
   \brief parse.hpp define and build the prefix-free parse data structure.
   \author Massimiliano Rossi
   \date 03/04/2020
*/

#ifndef PARSE_HPP
#define PARSE_HPP

#include "common.hpp"

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sys/stat.h>
#include <algorithm>

class parse{
public:
    std::vector<uint32_t> p;
    std::vector<uint32_t> stp;
    std::vector<uint32_t> len;
    std::vector<uint32_t> saP;
    std::vector<uint32_t> ilP;
    std::vector<uint32_t> ebwtP;
    bool saP_flag = false;
    bool ilP_flag = false;
    size_t alphabet_size;
    typedef size_t size_type;

  // Default constructor for load
    parse() {}
  
   parse(std::vector<uint32_t>& p_,
         std::vector<uint32_t>& sts_,
         std::vector<uint32_t>& len_,
         size_t alphabet_size_,
         bool saP_flag_ = true,
         bool ilP_flag_ = true):
         p(p_),
         alphabet_size(alphabet_size_)
    {
        build(saP_flag_, ilP_flag_);
    } 
    
    parse(std::string filename,
          size_t alphabet_size_,
          bool saP_flag_ = true,
          bool ilP_flag_ = true):
          alphabet_size(alphabet_size_)
  {
    // read file
    std::string tmp_filename = filename + std::string(".parse");
    read_file(tmp_filename.c_str(), p);
    // remove 0s from parse, to change, non efficient
    p.erase(remove(p.begin(),p.end(), 0),p.end());
    // read starting positions
    tmp_filename = filename + std::string(".start");
    read_file(tmp_filename.c_str(), stp);
    // read lengths
    tmp_filename = filename + std::string(".len");
    read_file(tmp_filename.c_str(), len);

    build(saP_flag_, ilP_flag_);
   }
    
    void build(bool saP_flag_, bool ilP_flag_){
        
        if(saP_flag_){
            saP.resize(p.size());
            ilP.resize(p.size());
            ebwtP.resize(p.size());
            // suffix array of the parse.
            verbose("Computing cSA of the parse");
            _elapsed_time(
                sortConjugatesMap(p,saP,ilP,ebwtP,stp,len,alphabet_size);
            );
        }
    }
    
      // Serialize to a stream.
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const
  {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_type written_bytes = 0;

    written_bytes += my_serialize(p, out, child, "parse");
    written_bytes += my_serialize(saP, out, child, "saP");
    written_bytes += my_serialize(ilP, out, child, "isaP");
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