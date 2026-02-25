#ifndef _RANDOM_HPP_
#define _RANDOM_HPP_

#include <random>

typedef std::mt19937 rng_type;
static std::uniform_int_distribution<rng_type::result_type>
    rand_int(0, std::numeric_limits<rng_type::result_type>::max());
static std::uniform_real_distribution<double> rand_double(0.0, 1.0);
// seed rng:
static rng_type rng(0); // Seed 0: deterministic

#endif //_RANDOM_HPP_