/*
 * Code to build the SA and LCP arrays of the dictionary of a prefix-free parse.
 * 
 * This code is adapted from https://github.com/maxrossi91/pfp-thresholds/blob/master/include/pfp/dictionary.hpp
 */

#ifndef PFP_PARSE_HPP
#define PFP_PARSE_HPP

#ifndef P64
  typedef uint32_t uint_s;
#else
  typedef uint64_t uint_s;
#endif

// TODO: Extend it to integer alphabets
class pfp_parse{
public:  
    std::vector<uint_s> ilP;
    std::vector<uint32_t> offset;
    sdsl::sd_vector<> b_il;
    sdsl::sd_vector<> b_st;
    sdsl::sd_vector<>::select_1_type select_ilist;
    sdsl::sd_vector<>::rank_1_type rank_st;

  pfp_parse() {}

  pfp_parse(std::string filename)
  {

    std::string input = filename + std::string(".sdsl");
    // load serialized parse's data structures.
    std::ifstream in(input);
    b_il.load(in);
    //select_ilist.load(in);
    b_st.load(in);
    //rank_st.load(in);
    my_load(ilP,in);
    my_load(offset,in);
    in.close();

    rank_st = sdsl::sd_vector<>::rank_1_type(&b_st);
    select_ilist = sdsl::sd_vector<>::select_1_type(&b_il);
    
  }
    
};
  

#endif /* PFP_PARSE_HPP */

