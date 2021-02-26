#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include <tuple>
#include <numeric>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <map>

using namespace std;

typedef size_t size_type;
typedef tuple<int,int,int> triplet;

bool sort_int (uint32_t i, uint32_t j) {return(i<j);}

void readParse(vector<uint32_t> &parse, vector<uint32_t> &starting, vector<uint32_t> &lengths, uint32_t &size, uint32_t &max_rank, FILE *parse_file){

    uint32_t rank = 0, w_len=0;
    starting.push_back(0);
    while(true){
        size_t s = fread(&rank,sizeof(rank),1,parse_file);
        if(s!=1){break;}
        if(rank != 0){
            parse.push_back(rank);
            max_rank = max(max_rank,rank);
            size++;
            w_len ++;
        }else{
            lengths.push_back(w_len);
            starting.push_back(w_len+starting.back());
            w_len = 0;
        }
    }
    starting.pop_back();
}

void sortConjugates(vector<uint32_t> &parse, vector<uint32_t> &sts, vector<uint32_t> &slen, uint32_t size, uint32_t max_rank, FILE *bwt){
    
    queue<triplet> buckets;
    buckets.push({0,size,0});
    vector<uint32_t> sa(size,0);
    vector<uint32_t> si(size,0);
    for(uint32_t i=0;i<sa.size();i++){
        sa[i] = i;
    }
    int nseq = 0;
    for(uint32_t i=0;i<slen.size();i++){
        for(uint32_t j=sts[i];j<sts[i]+slen[i];j++){
            si[j] = nseq;
        }
        nseq++;
    }
    // Sort the conjugates
    while(!buckets.empty()){
        auto bucket = buckets.front(); buckets.pop();
        int start = get<0>(bucket);
        int end   = get<1>(bucket);
        int depth = get<2>(bucket);
       
        if((start < end)){
            vector<uint32_t> count(max_rank,0);
            for(size_t i=start; i<end; ++i)
            {   
                uint32_t cs = si[sa[i]];
                size_t ind = sts[cs]+((sa[i]-sts[cs]+depth)%(slen[cs]));
                count[parse[ind]-1]++;
            }              
            vector<uint32_t> psum(count.size(),0);
            for(size_t i=1; i<count.size(); ++i)
            {
                psum[i] = psum[i-1] + count[i-1];
            }           
            vector<uint32_t> tmp(end - start, 0);
            for (size_t i = start; i < end; ++i)
            {
                uint32_t cs = si[sa[i]];
                size_t ind = sts[cs] + ((sa[i]-sts[cs]+depth)%(slen[cs]));
                size_t index = psum[parse[ind]-1]++;
                tmp[index] = sa[i];
            }          
            for (uint32_t i=start; i<tmp.size()+start;i++){
                sa[i] = tmp[i-start];
            }
            tmp.clear();
            for(size_t i=0; i<count.size(); i++)
            {
                
                end = start + count[i];
                if(end - start > 1){buckets.push({start,end,depth+1});}
                start = end;
            }
        }
    }
    // Store eBWT to file.
    uint32_t bchar = 0;
    for (int i=0; i<sa.size();i++){
        if(sa[i]==0 || si[sa[i]] != si[sa[i]-1]){
            bchar = parse[(sa[i]-1)+slen[si[sa[i]]]]-1;
            fwrite(&bchar,sizeof(bchar),1,bwt);
        }else{
            bchar = parse[sa[i]-1]-1;
            fwrite(&bchar,sizeof(bchar),1,bwt);
        }
    }
}

void sortConjugatesMap(vector<uint32_t> &parse, vector<uint32_t> &sts, vector<uint32_t> &slen, uint32_t size, uint32_t max_rank, FILE *bwt){
    
    queue<triplet> buckets;
    buckets.push({0,size,0});
    vector<uint32_t> sa(size,0);
    vector<uint32_t> si(size,0);
    for(uint32_t i=0;i<sa.size();i++){
        sa[i] = i;
    }
    int nseq = 0;
    for(uint32_t i=0;i<slen.size();i++){
        for(uint32_t j=sts[i];j<sts[i]+slen[i];j++){
            si[j] = nseq;
        }
        nseq++;
    }
    // Sort the conjugates
    cout << "Sorting the conjugates..." << endl;
    while(!buckets.empty()){
        auto bucket = buckets.front(); buckets.pop();
        int start = get<0>(bucket);
        int end   = get<1>(bucket);
        int depth = get<2>(bucket);
        
        if((start < end)){
            vector<uint32_t> count;
            map <uint32_t,uint32_t> count_map;
            for(size_t i=start; i<end; ++i)
            {   
                uint32_t cs = si[sa[i]];
                size_t ind = sts[cs]+((sa[i]-sts[cs]+depth)%(slen[cs]));
                uint32_t symb = parse[ind]-1;
                if(count_map.find(symb)==count_map.end()){
                    count_map.insert(pair<uint32_t,uint32_t>(symb,1));
                }else{
                    count_map.at(symb)++;
                }
            }
            map <uint32_t,uint32_t> psum_map;
            uint32_t pind = 0;
            for(auto& x: count_map){
                count.push_back(x.second);
                psum_map.insert(pair<uint32_t,uint32_t>(x.first,pind));
                pind++;
            }
            count_map.clear();
            vector<uint32_t> psum(count.size(),0);
            for(size_t i=1; i<count.size(); ++i)
            {
                psum[i] = psum[i-1] + count[i-1];
            }   
            vector<uint32_t> tmp(end - start, 0);
            for (size_t i = start; i < end; ++i)
            {
                uint32_t cs = si[sa[i]];
                size_t ind = sts[cs] + ((sa[i]-sts[cs]+depth)%(slen[cs]));
                size_t index = psum[psum_map[parse[ind]-1]]++;
                tmp[index] = sa[i];
            }
            for (uint32_t i=start; i<tmp.size()+start;i++){
                sa[i] = tmp[i-start];
            }
            tmp.clear();
            for(size_t i=0; i<count.size(); i++)
            {
                
                end = start + count[i];
                if(end - start > 1){buckets.push({start,end,depth+1});}
                start = end;
            }
        }
    }
    // Store eBWT to file
    uint32_t bchar = 0;
    for (int i=0; i<sa.size();i++){
        if(sa[i]==0 || si[sa[i]] != si[sa[i]-1]){
            bchar = parse[(sa[i]-1)+slen[si[sa[i]]]]-1;
            fwrite(&bchar,sizeof(bchar),1,bwt);
        }else{
            bchar = parse[sa[i]-1]-1;
            fwrite(&bchar,sizeof(bchar),1,bwt);
        }
    }
}

int main(int argc, char** argv) {
    
    cout << "Reading parse file..." << endl;
    string inputFileName = argv[1];
    string parse_fname = inputFileName + ".parse";
    FILE* p = fopen(parse_fname.c_str(), "r");
    string bwt_fname = inputFileName + ".bwt";
    FILE* fb = fopen(bwt_fname.c_str(), "wb");

    //import parse sequences
    vector<uint32_t> parse, sts;
    vector<uint32_t> slen;
    uint32_t size = 0, max_rank = 0;

    readParse(parse,sts,slen,size,max_rank,p);
    
    sortConjugatesMap(parse,sts,slen,size,max_rank,fb);

    return 0;
    
}