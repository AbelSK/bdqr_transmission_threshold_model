#ifndef COMMON_H
#define COMMON_H
#include <cmath>

double get_rand_uniform();
std::int64_t get_rand_poisson(double mean);
std::int64_t get_rand_binom(std::int64_t n, double p);

#endif