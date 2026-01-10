#include "common.h"
#include <cmath>
#include <random>
#include "xoshiro.h" 
#include <cstdint>

namespace {
    thread_local xso::rng gen;

    thread_local std::poisson_distribution<std::int64_t> pois;
    thread_local std::binomial_distribution<std::int64_t> binom;

    inline double get_rand_normal(double mean, double stddev) {
        static thread_local std::normal_distribution<double> norm;
        norm.param(std::normal_distribution<double>::param_type(mean, stddev));
        return norm(gen);
    }

    constexpr std::int64_t BINOM_EXACT_THRESHOLD = 10'000;
    constexpr double POISSON_EXACT_THRESHOLD = 10'000.0;
}

double get_rand_uniform() {
    std::uint64_t raw = gen();
    return (raw >> 11) * (1.0 / (1ULL << 53));
}

std::int64_t get_rand_poisson(double mean) {
    if (mean <= POISSON_EXACT_THRESHOLD) {
        pois.param(std::poisson_distribution<std::int64_t>::param_type(mean));
        return pois(gen);
    } else {
        double sample = get_rand_normal(mean, std::sqrt(mean));
        return std::max<std::int64_t>(0, std::llround(sample));
    }
}

std::int64_t get_rand_binom(std::int64_t n, double p) {
    if (n <= BINOM_EXACT_THRESHOLD) {
        binom.param(std::binomial_distribution<std::int64_t>::param_type(n, p));
        return binom(gen);
    } else {
        double mean = n * p;
        double stddev = std::sqrt(n * p * (1.0 - p));
        double sample = get_rand_normal(mean, stddev);
        return std::clamp<std::int64_t>(std::llround(sample), 0, n);
    }
}


// std::int64_t get_rand_poisson(double mean) {
//     std::poisson_distribution<std::int64_t> dist(mean);
//     return dist(gen);
// }

// std::int64_t get_rand_binom(std::int64_t n, double p) {
//     std::binomial_distribution<std::int64_t> dist(n, p);
//     return dist(gen);
// }