/*
 * Code to build the the runlength eBWT and the sampled suffix array of a string collections.
 * 
 * This code is adapted from https://github.com/maxrossi91/pfp-thresholds/blob/master/include/pfp/pfp_thresholds.hpp
 */

#ifndef PFP_EBWT_GCA_HPP
#define PFP_EBWT_GCA_HPP

#include <algorithm>
#include <tuple>
#include <queue>

typedef std::tuple<bool,size_t,size_t> i_tuple;
typedef std::pair<size_t,size_t> il_interval;

class pfp_ebwt_gca{
private:
    typedef struct
    {
        size_t i = 0; 
        size_t phrase = 0;
        size_t suffix_length = 0;
        size_t sn = 0;
        uint8_t bwt_char = 0;
        size_t  st_pos = 0;
        bool is_starting = 0;
    } phrase_suffix_t;

    size_t text_len = 0; // total length of text
    size_t ins_sofar = 0; // Number of characters inserted in the eBWT
    size_t runs = 0; // Number of runs of the eBWT
    int bytes_da, bytes_gca; // Number of bytes to store the GCA entries

    // run length encoded eBWT
    bool rle;
    bool sample_first;
    
    // output files
    FILE *ebwt_file;
    FILE *ebwt_file_len;
    FILE *ebwt_file_heads;
    // for strings tarting positions
    FILE *I_file;
    // gca samples files
    FILE *ebwt_file_ssa;
    FILE *ebwt_file_esa;
    FILE *ebwt_file_dssa;
    FILE *ebwt_file_desa;

public:

    pfp_parse_gca& pars;
    dictionary& dict;
    std::string filename;
    size_t w;
    size_t min_s; 
    size_t pos_s; 

    uint8_t head; // Head of the current run of BWT_T
    size_t start;
    size_t  length = 0; // Length of the current run of BWT_T

    // First occurrence of same suffix phrases in BWT_P
    size_t first_occ = 0;
    // Last occurrence of same suffix phrases in BWT_P 
    size_t last_occ = 0, p_last_occ = 0;   
    // starting sample initialization
    size_t ssa = 0, dssa = 0;
    // ending sample initialization
    size_t esa = 0, desa = 0;

    // initialize struct for current suffix
    phrase_suffix_t curr, prev;

    // function processing next suffix
    // it returns false if the suffix is not valid
    inline bool inc(phrase_suffix_t& s)
    {
        s.i++;
        if (s.i >= dict.saD.size())
            return false;
        s.sn = dict.saD[s.i];
        s.is_starting = dict.b_s[s.sn];
        s.phrase = dict.rank_b_d(s.sn+1);                             //
        if(s.is_starting){ s.st_pos=s.sn-dict.select_b_d(s.phrase); } //10000001000001
        else{ s.st_pos=0; }
        s.suffix_length = dict.select_b_d(dict.rank_b_d(s.sn + 1) + 1) - s.sn -1;
        assert(s.phrase > 0 && s.phrase <= pars.ilP.size());
        
        if(is_valid(s)){
            s.bwt_char = dict.d[s.sn - 1];}
        else {s.bwt_char = 0;}
        return true;
    }

    // function processing the next suffix
    // it returns false if the suffix is not valid
    inline void update_ebwt_sa(uint8_t next_char, size_t length_, phrase_suffix_t &curr, il_interval occ, i_tuple st_pos)
    {
        if (head != next_char)
        {
            // update last gCA sample
            update_sa_sample(prev, p_last_occ, 1);
            p_last_occ = std::get<1>(occ);
            prev = curr;
            // udate first gCA sample
            update_sa_sample(curr, std::get<0>(occ), 0);
            // print a run of the ebwt
            print_ebwt();
            // print SA samples
            print_sa();
            // update BWT head and run length
            head = next_char;
            length = 0;
        }
        //else
        //{
        //    p_last_occ = std::get<1>(occ);
        //    prev = curr;
        //}
        
        ins_sofar += length_;
        
        // print indexes for I file, positions of the strings in gCA
        if(std::get<0>(st_pos) && std::get<2>(st_pos) > 0){
            // check if we are on a starting position
            if(pars.b_st[std::get<1>(st_pos)]==1){
                if(pars.offset[(pars.rank_st(std::get<1>(st_pos)+1)-1)] == std::get<2>(st_pos)){
                    // print the current position in the eBWT
                    start = ins_sofar-1;
                    if(fwrite(&start,sizeof(start),1,I_file)!=1) error("I file write error");
                    // check if this position is already sampled
                    if(sample_first && length > 0)
                    {
                        // update last gCA sample
                        update_sa_sample(prev, p_last_occ, 1);
                        p_last_occ = std::get<1>(occ);
                        prev = curr;
                        // udate first gCA sample
                        update_sa_sample(curr, std::get<0>(occ), 0);
                        // print a run of the ebwt
                        print_ebwt();
                        // print SA samples
                        print_sa();
                        // update BWT head and run lengt
                        length = 0;
                    }
                }
            }
        }

        if ( length > 0 )
        {
            p_last_occ = std::get<1>(occ);
            prev = curr;
        }
        
        length += length_;
        
    }

    inline void print_sa(){

        // skip ending sample of empty run
        if(ins_sofar > 0){
            if (fwrite(&esa, bytes_gca, 1, ebwt_file_esa) != 1)
                error("GCA write error 1");
            if (fwrite(&desa, bytes_da, 1, ebwt_file_desa) != 1)
                error("DA write error 1");
        }
        if (fwrite(&ssa, bytes_gca, 1, ebwt_file_ssa) != 1)
            error("GCA write error 2");
        if (fwrite(&dssa, bytes_da, 1, ebwt_file_dssa) != 1)
            error("DA write error 2");
    }

    // function processing the next suffix
    // it returns false if the suffix is not valid
    inline void print_ebwt()
    {
        if(length > 0)
        {
            if(rle)
            {
                // Write the head
                if (fputc(head, ebwt_file_heads) == EOF)
                    error("BWT write error 1");
                // Write the length
                if (fwrite(&length, BWTBYTES, 1, ebwt_file_len) != 1)
                    error("BWT write error 2");
            }else{
                // write plain ebwt
                for(size_t i = 0; i < length; ++i)
                {
                    if (fputc(head, ebwt_file) == EOF)
                        error("BWT write error 1");
                }
            }
            // one run added
            ++runs;
        }
    }

    void print_last_sa(){
        // compute last sample
        update_sa_sample(curr, p_last_occ, 1);
        // print last sample of the last eBWT run
        if (fwrite(&esa, bytes_gca, 1, ebwt_file_esa) != 1)
            error("GCA write error 1");
        if (fwrite(&desa, bytes_da, 1, ebwt_file_desa) != 1)
            error("DA write error 1");
    }

    // function that returns true if the suffix is valid false otherwise
    inline bool is_valid(phrase_suffix_t& s)
    {
        // Check if the suffix has length at least w and is not the complete phrase.
        if (s.suffix_length < w || dict.b_d[s.sn]==1)
            return false;
        
        return true;
    }
    // function that updates the GCA-samples
    inline void update_sa_sample(phrase_suffix_t &curr, size_t pos, bool lsam)
    {   // get offset of the phrase
        int64_t last_val = pars.get_last(pos);
        size_t rank = pars.rank_ns(last_val+1);
        int64_t first_string_pos = pars.select_ns(rank);
        // compute temp GCA-sample
        int64_t sa_tmp = last_val - curr.suffix_length + 1;
        // circular GCA-sample computation
        if( first_string_pos >  sa_tmp )
            sa_tmp = pars.select_ns(rank+1) - ( curr.suffix_length - (last_val - first_string_pos) ) + 1;
        // update either esa or ssa
        if( lsam ){ desa = rank-1; esa = sa_tmp - first_string_pos; }
        else { dssa = rank-1; ssa = sa_tmp - first_string_pos; }
    }

    // initialize first and last occurence variables
    inline void init_first_last_occ()
    {
        // First occurrence of same suffix phrases in BWT_P
        first_occ = pars.bwtLen; 
        // Last occurrence of same suffix phrases in BWT_P
        last_occ = 0;     
    }

    // compute the first and last position in the inverted list for the sa samples
    inline void update_first_last_occ(phrase_suffix_t &curr) 
    {
        size_t begin = pars.select_ilist(curr.phrase); 
        size_t end = pars.select_ilist(curr.phrase + 1) - 1; 

        size_t first = pars.ilP[begin];
        size_t last  = pars.ilP[end];

        if(first_occ > first){ first_occ = first; }

        if(last_occ < last){ last_occ = last; }
    }


    pfp_ebwt_gca(pfp_parse_gca &p_, dictionary &d_, std::string filename_, size_t w_, bool rle_, bool sample_first_): 
        pars(p_),
        dict(d_),
        filename(filename_),
        w(w_),
        rle(rle_),
        sample_first(sample_first_),
        pos_s(0),
        head(0)
    {
    	// initialize output files
    	// HEADS file
        std::string outfile;
        if( rle ){
    	outfile = filename + std::string(".head"); 
    	if((ebwt_file_heads = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");
        // HEADS file
    	outfile = filename + std::string(".len"); 
    	if((ebwt_file_len = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");
        }
        else{
            outfile = filename + std::string(".ebwt"); 
            if((ebwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");
        }
        // I file
        outfile = filename + std::string(".I");
        if((I_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");
        // starting samples file
        outfile = filename + std::string(".ssam"); 
        if((ebwt_file_ssa = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");
        outfile = filename + std::string(".sda"); 
        if((ebwt_file_dssa = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");
        // ending samples file
        outfile = filename + std::string(".esam"); 
        if((ebwt_file_esa = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");
        outfile = filename + std::string(".eda"); 
        if((ebwt_file_desa = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        bytes_gca = std::ceil(std::log2(pars.get_ebwt_size())/8);
        bytes_da = std::ceil(std::log2(pars.get_no_seq())/8);

        // increment the current suffix
        inc(curr); prev = curr;
        // iterate over SA of dict
        while (curr.i < dict.saD.size())
        {
            // if the current suffix is valid process it
            if(is_valid(curr))
            {
                // initialize vector containing identical suffixes
                std::vector<phrase_suffix_t> same_suffix(1, curr);
                // same_chars == true if all suffixes are preceded by the
                // same BWT characters. st_chars == true if some of the 
                // suffixes are preceded by a character that starts a string
                bool same_chars = true;
                bool st_chars = curr.is_starting;
                // initialize first and last occurrence in the parse inverted list
                init_first_last_occ();
                update_first_last_occ(curr);
                // initialize the variable for the next suffix to process
                phrase_suffix_t next = curr;
                // process all identical suffixes in one block
                while (inc(next) && (dict.lcpD[next.i] >= curr.suffix_length))
                {
                    // if the two suffixes have the same length
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_chars = (same_chars && same_suffix.back().bwt_char == next.bwt_char);
                        st_chars = (st_chars || next.is_starting);
                        same_suffix.push_back(next);
                        update_first_last_occ(next);
                    }                  
                }
                // simple case: all suffixes preceded by the same BWT char
                //              and not starting characters
                if (same_chars && !st_chars)
                {
                    size_t block_length = 0;
                    for (auto curr : same_suffix)
                    {
                        block_length += dict.occ[curr.phrase-1];
                    }
                    update_ebwt_sa(same_suffix[0].bwt_char,block_length,same_suffix[0],
                                   std::make_pair(first_occ,last_occ),std::make_tuple(0, 0, 0));
                }   
                // hard case
                else
                {       
                    //suffix not starting with a character occurring at the beginning of 
                    // a input sequence
                    if(!st_chars){
                        typedef std::pair<uint_s *, std::pair<uint_s *, uint8_t>> pq_t;
                        // using lambda to compare elements.
                        auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
                            return *lhs.first > *rhs.first;
                        };
                        std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
                        
                        for (auto s: same_suffix){
                            size_t begin = pars.select_ilist(s.phrase);
                            size_t end = pars.select_ilist(s.phrase+1);
                            pq.push({&pars.ilP[begin], {&pars.ilP[end], s.bwt_char}});
                        }
                        while(!pq.empty()){
                            auto curr_occ = pq.top();
                            pq.pop();

                            //if(head != curr_occ.second.second) update_ssa(curr, *curr_occ.first);
                            update_ebwt_sa(curr_occ.second.second,1,curr,
                                           std::make_pair(*curr_occ.first,*curr_occ.first),std::make_tuple(0, 0, 0));
                            //update_esa(curr, *curr_occ.first);
                            
                            // Update pq
                            curr_occ.first++;
                            if (curr_occ.first != curr_occ.second.first)
                                pq.push(curr_occ);
                        
                        }                      
                    }else{
                        //check and store the positions at which the start of string characters
                        //are inserted in the ebwt (we need them to invert the ebwt)
                        typedef std::pair<uint_s *, std::tuple<uint_s *, uint8_t, size_t, size_t>> tq_t;
                        
                        auto cmp2 = [](const tq_t &lhs, const tq_t &rhs) {
                            return *lhs.first > *rhs.first;
                        };
                        std::priority_queue<tq_t, std::vector<tq_t>, decltype(cmp2)> tq(cmp2);
                        for (auto s: same_suffix){
                            size_t begin = pars.select_ilist(s.phrase);
                            size_t end = pars.select_ilist(s.phrase+1);
                            tq.push({&pars.ilP[begin], std::make_tuple(&pars.ilP[end], s.bwt_char, begin, s.st_pos)});
                        }
                        while(!tq.empty()){
                            auto curr_occ = tq.top();
                            tq.pop();
                            update_ebwt_sa(std::get<1>(curr_occ.second),1,curr,std::make_pair(*curr_occ.first,*curr_occ.first),
                                           std::make_tuple(1,std::get<2>(curr_occ.second),std::get<3>(curr_occ.second)));
                            
                            // Update pq
                            curr_occ.first++;
                            std::get<2>(curr_occ.second)++;
                            if (curr_occ.first != std::get<0>(curr_occ.second))
                                tq.push(curr_occ);
                        
                        }
                    }
                }
                // the next suffix is now the current suffix
                curr = next;
            }
            // increment the suffix if the previous one was not valid
            else { inc(curr); }
        }
        // print last eBWT run and gCA sample
        print_ebwt();
        print_last_sa();
        // close output files
        fclose(I_file);
        // close rle lengths file if rle was used
        if(rle){ fclose(ebwt_file_len); fclose(ebwt_file_heads);}
        else { fclose(ebwt_file); }   
        // close gCA sample files
        fclose(ebwt_file_ssa);
        fclose(ebwt_file_esa); 
        fclose(ebwt_file_dssa);
        fclose(ebwt_file_desa); 

        std::string info = "Number of runs of the eBWT: " + std::to_string(runs) + "\n"; 
        info += "Length of the eBWT: " + std::to_string(pars.get_ebwt_size()) + "\n"; 
        if(rle)
        {
            info += "Number of bytes per eBWT run character (.heads file): 1\n";
            info += "Number of bytes per eBWT run length (.len file): " + std::to_string(BWTBYTES) + "\n";
        } 
        else
            info += "Number of bytes per eBWT character (.ebwt): 1\n";
        info += "Number of bytes per GCA samples conjugate offset (.ssam and .esam files): " + std::to_string(bytes_gca) + "\n";
        info += "Number of bytes per GCA sequence offset (.sda and .eda files): " + std::to_string(bytes_da) + "\n";

        std::ofstream out(filename + std::string(".info"));
        out << info;
        out.close();
    }
};

#endif /* PFP_EBWT_GCA_HPP */