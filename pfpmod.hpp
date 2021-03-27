/* 
 * File:   pfp.hpp
 * Author: hejans
 *
 * Created on February 20, 2021, 11:12 PM
 */

#ifndef PFP_HPP
#define PFP_HPP

#include "parse.hpp"
#include "dictionary.hpp"

#include <algorithm>
#include <tuple>
#include <queue>

#define BWTBYTES 5

class pfp{
private:
    typedef struct
    {
        size_t i = 0; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        size_t phrase = 0;
        size_t suffix_length = 0;
        size_t sn = 0;
        uint8_t bwt_char = 0;
        size_t  st_pos = 0;
        bool is_starting = 0;
    } phrase_suffix_t;
    
    size_t first_occ = 0; // First occurrence of same suffix phrases in BWT_P
    size_t last_occ = 0;     // Last occurrence of same suffix phrases in BWT_P
    size_t text_len = 0;
    size_t ins_sofar = 0;
    
    bool rle;
    
    // for rle eBWT
    FILE *bwt_file_len;
    FILE *bwt_file_heads;
    // for plain eBWT
    FILE *bwt_file;
    FILE *I_file;
    
public:

    parse& pars;
    dictionary& dict;
    std::string filename;
    size_t w;
    size_t min_s; // Value of the minimum lcp_T in the current run of BWT_T
    size_t pos_s; // Position of the minimum lcp_T in the current run of BWT_T

    uint8_t head; // Head of the current run of BWT_T
    size_t start;
    size_t  length = 0; // Length of the current run of BWT_T
    
    inline bool is_valid(phrase_suffix_t& s)
    {
        // Check if the suffix has length at least w and is not the complete phrase.
        if (s.suffix_length < w || dict.b_d[s.sn]==1)
            return false;
        
        return true;
    }
    
    //static bool sort_pairs (pq_t i, pq_t j) {return(i.first<j.first);}
    
    inline bool inc(phrase_suffix_t& s)
    {
        s.i++;
        if (s.i >= dict.saD.size())
            return false;
        s.sn = dict.saD[s.i];
        s.is_starting = dict.b_s[s.sn];
        s.phrase = dict.rank_b_d(s.sn+1);
        if(s.is_starting){s.st_pos=s.sn-dict.select_b_d(s.phrase);}
        s.suffix_length = dict.select_b_d(dict.rank_b_d(s.sn + 1) + 1) - s.sn -1;
        assert(!is_valid(s) || (s.phrase > 0 && s.phrase < pars.ilP.size()));
        
        if(is_valid(s)){
            s.bwt_char = dict.d[s.sn - 1];}
        else {s.bwt_char = 0;}
        return true;
    }
    
    inline void print_ebwt()
    {
        if(length > 0)
        {
            if(rle)
            {
                // Write the head
                if (fputc(head, bwt_file) == EOF)
                    error("BWT write error 1");
                
                // Write the length
                if (fwrite(&length, BWTBYTES, 1, bwt_file_len) != 1)
                    error("BWT write error 2");
            }else{

                for(size_t i = 0; i < length; ++i)
                {
                    if (fputc(head, bwt_file) == EOF)
                        error("BWT write error 1");
                }
                // uint8_t* run = new uint8_t[length];
                // memset(run,head,length);

                // if (fwrite(run, 1, length, bwt_file) != length)
                //     error("BWT write error 1");

                // delete run;
            }
        }
    }
    
    inline void update_ebwt(uint8_t next_char, size_t length_)
    {
        if (head != next_char)
        {
            // print a run of the ebwt
            print_ebwt();
            ins_sofar += length;

            head = next_char;
            // lengths.push_back(0); // Debug only
            length = 0; // Debug only
        }

        length += length_;
        // lengths.back() += length; // Debug only
    }
    
    inline void update_ebwt_sp(uint8_t next_char, size_t ind, size_t st_pos)
    {
        if (head != next_char)
        {
            // insert a run in the eBWT
            print_ebwt();

            head = next_char;
            // lengths.push_back(0); // Debug only
            length = 0; // Debug only
        }
        
        ++ins_sofar;
        ++length;
        
        if(pars.b_st[ind]==1){
            if(pars.offset[(pars.rank_st(ind+1)-1)] == st_pos){
                start = ins_sofar-1;
                if(fwrite(&start,sizeof(start),1,I_file)!=1) error("I file write error");
            }
        }
        // lengths.back() += length; // Debug only
    }

    pfp(parse &p_, dictionary &d_, std::string filename_, size_t w_, bool rle_): 
            pars(p_),
            dict(d_),
            filename(filename_),
            w(w_),
            rle(rle_),
            pos_s(0),
            head(0)
    {
        assert(dict.d[dict.saD[0]] == EndOfDict);
        
        std::string outfile = "";
        if(rle)
        {
            
            outfile = filename + std::string(".ebwt.len");
            if ((bwt_file_len = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

        }else{

            outfile = filename + std::string(".ebwt");
            if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
                error("open() file " + outfile + " failed");

        }
        
        outfile = filename + std::string(".I");
        if((I_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        phrase_suffix_t curr;
        
        inc(curr);
        while (curr.i < dict.saD.size())
        {
            if(is_valid(curr))
            {
                // Compute the next character of the BWT of T
                std::vector<phrase_suffix_t> same_suffix(1, curr);  // Store the list of all phrase ids with the same suffix.

                bool same_chars = true;
                bool st_chars = curr.is_starting;

                phrase_suffix_t next = curr;

                while (inc(next) && (dict.lcpD[next.i] >= curr.suffix_length))
                {
                    assert(next.suffix_length >= curr.suffix_length);
                    assert((dict.b_d[next.sn] == 0 && next.suffix_length >= w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_chars = (same_chars && same_suffix.back().bwt_char == next.bwt_char);
                        st_chars = (st_chars || next.is_starting);
                        same_suffix.push_back(next);
                    }                  
                }
                // Simple case
                if (same_chars && !st_chars){

                    for (auto curr : same_suffix)
                    {
                        //head = curr.bwt_char;
                        //for(int i=0;i<dict.occ[curr.phrase-1];i++)
                        //{                          
                        //    if(fputc(head, bwt_file) == EOF) error("BWT write error 1"); ins_sofar++;
                        //}
                        update_ebwt(curr.bwt_char, dict.occ[curr.phrase-1]);
                    }
                    
                }   
                // Hard case
                else
                {       
                    //suffix not starting with a character occurring at the beginning of 
                    // a input sequence
                    if(!st_chars){
                        
                        typedef std::pair<int_t *, std::pair<int_t *, uint8_t>> pq_t;
                        // using lambda to compare elements.
                        auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
                            return *lhs.first > *rhs.first;
                        };
                        
                        std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
             
                        //std::vector<pq_t> cand;
                        
                        /*for (auto s: same_suffix){
                            size_t begin = pars.select_ilist(s.phrase);
                            size_t end = pars.select_ilist(s.phrase+1);
                            for(int i=begin;i<end;i++)
                            {
                                cand.push_back({pars.ilP[i],s.bwt_char});
                            }
                        }*/
                        
                        for (auto s: same_suffix){
                            size_t begin = pars.select_ilist(s.phrase);
                            size_t end = pars.select_ilist(s.phrase+1);
                            pq.push({&pars.ilP[begin], {&pars.ilP[end], s.bwt_char}});
                        }
                        exit(1);
                        //size_t prev_occ;
                        while(!pq.empty()){
                            auto curr_occ = pq.top();
                            pq.pop();
                            
                            update_ebwt(curr_occ.second.second, 1);
                            
                            // Update pq
                            curr_occ.first++;
                            if (curr_occ.first != curr_occ.second.first)
                                pq.push(curr_occ);
                        
                        }
                        //std::sort(cand.begin(),cand.end(),sort_pairs);
                        
                        //for(int i=0;i<cand.size();i++){
                        //    head = cand[i].second;
                        //    if(fputc(head, bwt_file) == EOF) error("BWT write error 2"); ins_sofar++; }
                        
                    }else{
                        //check and store the positions at which the start of string characters
                        //are inserted in the ebwt (we need them to invert the ebwt)
                        //std::vector<tq_t> cand2;
                //        typedef std::pair<int_t *, std::tuple<int_t *, uint8_t, size_t, size_t>> tq_t;
                        //typedef std::tuple<size_t, uint8_t, size_t, size_t> tq_t;
               //         auto cmp2 = [](const tq_t &lhs, const tq_t &rhs) {
               //             return *lhs.first > *rhs.first;
               //         };
                        
               //         std::priority_queue<tq_t, std::vector<tq_t>, decltype(cmp2)> tq(cmp2);
                        
               //         for (auto s: same_suffix){
               //             size_t begin = pars.select_ilist(s.phrase);
               //             size_t end = pars.select_ilist(s.phrase+1);
               //             tq.push({&pars.ilP[begin], std::make_tuple(&pars.ilP[end], s.bwt_char, begin, s.st_pos)});
               //         }
               //         
               //         while(!tq.empty()){
               //             auto curr_occ = tq.top();
               //             tq.pop();
               //             
               //             update_ebwt_sp(std::get<1>(curr_occ.second),std::get<2>(curr_occ.second),std::get<3>(curr_occ.second));
               //             
                            // Update pq
               //             curr_occ.first++;
               //             std::get<2>(curr_occ.second)++;
               //             if (curr_occ.first != std::get<0>(curr_occ.second))
               //                 tq.push(curr_occ);
                        
               //         }
                        
                        //for (auto s: same_suffix)
                        //{
                        //    size_t begin = pars.select_ilist(s.phrase);
                        //    size_t end = pars.select_ilist(s.phrase+1);
                        //    for(int i=begin;i<end;i++)
                        //    {
                        //        cand2.push_back(std::make_tuple(pars.ilP[i],s.bwt_char,i,s.st_pos));
                        //    }
                        //}
                        
                        //std::sort(cand2.begin(),cand2.end());
                        
                        
                        /*for(int i=0;i<cand2.size();i++){
                            head = std::get<1>(cand2[i]);
                            if(fputc(head, bwt_file) == EOF) error("BWT write error 3"); ins_sofar++;
                            size_t index = std::get<2>(cand2[i]);
                            //if character in the eBWT of the parse is marked as containing a start of string 
                            // character, store its position in the ebwt.
                            if(pars.b_st[index]==1){
                                if(pars.offset[(pars.rank_st(index+1)-1)] == std::get<3>(cand2[i])){
                                    start = ins_sofar-1;
                                    if(fwrite(&start,sizeof(start),1,I_file)!=1) error("I file write error");
                                }
                            }
                        }*/
                    }
                }
                curr = next;
            }
            else
            {
                inc(curr);
            }      
        }
        
        // print last characters
        print_ebwt();
        // close output files
        fclose(bwt_file);
        fclose(I_file);
        // close rle lengths file if rle was used
        if(rle){ fclose(bwt_file_len); }   
    }        
};


#endif /* PFP_HPP */