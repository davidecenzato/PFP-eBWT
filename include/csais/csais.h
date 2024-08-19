/*
 * Circular SAIS implementations to compute the circular SA of an integer vector.
 * 
 * This code is adapted from https://github.com/kurpicz/saca-bench/blob/master/sa-is/sais.cpp
 * which is the original code of the SA-IS algorithm listed below, and 
 * from https://github.com/felipelouza/gsa-is/blob/master/gsacak.c which is an implementation
 * of the GSACA-K algorithm.
 *
 * csais_int(s, SA, n, K, b_s) // computes circular SA of an integer vector using cSSAIS. 
 *
 */
// This is the sample code for the SA-IS algorithm presented in
// our article "Two Efficient Algorithms for Linear Suffix Array Construction"
// (co-authored with Sen Zhang and Wai Hong Chan),
// which can be retrieved at: http://code.google.com/p/ge-nong/

#ifndef CSAIS_H
#define CSAIS_H

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <inttypes.h>
#include <string.h>
#include <memory.h>

#include <iostream>
#include <ctime>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>


#ifndef P64
	#define P64 0
#else
	#define P64 1
#endif

#if P64 == 1
    typedef uint64_t uint_s;
#else
    typedef uint32_t uint_s;
#endif

typedef uint32_t uint_p;

/**
 * @brief Compute the suffix array of a string of integers using cSACA algorithm 
 * 
 * @param s input string
 * @param SA suffix array
 * @param n length of the input string
 * @param K alphabet size
 * @param b_s bitvector of the starting phrases of the parse
 */
void csais_int(uint_p *s, uint_s *SA, size_t n, size_t K, sdsl::sd_vector<> &b_s);

#endif /* CSAIS_H */

