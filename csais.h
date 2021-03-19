/* 
 * File:   csais.h
 * Author: hejans
 *
 * Created on March 12, 2021, 11:50 AM
 */

#ifndef CSAIS_H
#define CSAIS_H

#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <sdsl/bit_vectors.hpp>

#ifndef M64
	#define M64 0
#endif

#if M64
	typedef uint64_t uint_t;
#else
	typedef uint32_t uint_t;
#endif


//void cSAIS(unsigned int *s, unsigned int *SA, int n, int K, int cs, int level, sdsl::bit_vector &b_s);
void cSAIS(unsigned int *s, uint_t *SA, int n, int K, int cs, int level, sdsl::bit_vector &b_s);

/** @brief computes the suffix array of string s[0..n-1] 
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

