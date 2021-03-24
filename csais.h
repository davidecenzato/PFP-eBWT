/* 
 * File:   csais.h
 * Author: hejans
 *
 * Created on March 12, 2021, 11:50 AM
 */

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

void csais_int(uint_p *s, uint_s *SA, int n, int K, sdsl::bit_vector &b_s);

/** @brief computes the circular suffix array of string s[0..n-1] 
 *
 *  @param s	 input string 
 *  @param SA    suffix array 
 *  @param n	 string length
 *  @param K	 alphabet size
 *  @param cs	 integer size
 *  @param level recursion level, debug only
 *  @param b_s   starting positions bit vector
 *  @return      None
 */



#endif /* CSAIS_H */

