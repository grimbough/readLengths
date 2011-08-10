#include "readLengths.h"
 

int SampleNormal (double mean, double sigma, int seed)
{
    // Create a Mersenne twister random number generator
    // that is seeded once with #seconds since 1970
    static mt19937 rng(static_cast<unsigned> (seed));
 
    // select Gaussian probability distribution
    normal_distribution<double> norm_dist(mean, sigma);
 
    // bind random number generator to distribution, forming a function
    variate_generator<mt19937&, normal_distribution<double> >  normal_sampler(rng, norm_dist);
 
    // sample from the distribution
    return (int) normal_sampler();
}

int SampleUniform (int min, int max, int seed)
{
    // Create a Mersenne twister random number generator
    // that is seeded once with #seconds since 1970
    static mt19937 rng(static_cast<unsigned> (seed));
 
    // select Gaussian probability distribution
    uniform_int<> uni_dist(min, max);
 
    // bind random number generator to distribution, forming a function
    variate_generator<mt19937&, uniform_int<> >  uniform_sampler(rng, uni_dist);
 
    // sample from the distribution
    return uniform_sampler();
}

double SampleBernoulli (double prob, int seed)
{
    // Create a Mersenne twister random number generator
    // that is seeded once with #seconds since 1970
    static mt19937 rng(static_cast<unsigned> (seed));
 
    // select Gaussian probability distribution
    bernoulli_distribution<> b_dist(prob);
 
    // bind random number generator to distribution, forming a function
    variate_generator<mt19937&, bernoulli_distribution<double> >  b_sampler(rng, b_dist);
 
    // sample from the distribution
    return b_sampler();
}
