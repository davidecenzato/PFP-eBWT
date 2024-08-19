/*
 * Code to build the SA and LCP arrays of the dictionary of a prefix-free parse.
 * 
 * This code is adapted from https://github.com/maxrossi91/pfp-thresholds/blob/master/include/pfp/dictionary.hpp
 */

#ifndef PFP_PARSE_GCA_HPP
#define PFP_PARSE_GCA_HPP

#ifndef P64
  typedef uint32_t uint_s;
#else
  typedef uint64_t uint_s;
#endif

class pfp_parse_gca{
public:  
    std::vector<uint_s> ilP;
    std::vector<uint32_t> offset;
    std::vector<uint8_t> last;
    sdsl::sd_vector<> b_il;
    sdsl::sd_vector<> b_st;
    sdsl::sd_vector<>::select_1_type select_ilist;
    sdsl::sd_vector<>::rank_1_type rank_st;
    std::vector<uint8_t> first;
    sdsl::sd_vector<> b_ns;
    sdsl::sd_vector<>::rank_1_type rank_ns;
    sdsl::sd_vector<>::select_1_type select_ns;
    size_t firstLen, lastLen, bwtLen;

  pfp_parse_gca() {}

  pfp_parse_gca(std::string filename)
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
    //for(int i=0;i<b_st.size();++i){ std::cout << b_st[i] << " "; }
    //std::cout << "\n";
    // compute rank and select data structures for the parse
    rank_st = sdsl::sd_vector<>::rank_1_type(&b_st);
    select_ilist = sdsl::sd_vector<>::select_1_type(&b_il);
    // read last positions file
    input = filename + std::string(".slast"); 
    read_file(input.c_str(), last);
    lastLen = last.size()/IBYTES;

    // read first positions file
    input = filename + std::string(".spos");
    read_file(input.c_str(), first);
    firstLen = first.size()/IBYTES;
    bwtLen = get_myint(&first[0],firstLen,firstLen-1);

    // compute the rank data structure for strings first positions
    sdsl::sd_vector_builder builder(bwtLen+1,firstLen);
    for(size_t i=0;i<firstLen;++i){ builder.set(get_myint(&first[0],firstLen,i)); }
    first.clear();
    b_ns = sdsl::sd_vector<>(builder);
    rank_ns = sdsl::sd_vector<>::rank_1_type(&b_ns);
    select_ns = sdsl::sd_vector<>::select_1_type(&b_ns);
  }

  size_t get_last(size_t i){ return get_myint(&last[0],lastLen,i); }
  size_t get_ebwt_size(){ return bwtLen; }
  size_t get_no_seq(){ return firstLen-1; }
    
};
  

#endif /* PFP_PARSE_GCA_HPP */