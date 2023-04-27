// Wrapper for random routines from :
//	1. ROOT
//	2. MT10037
//	3. GSL
//
// Programmer: Xihe Han
//
// Changelog:
//	- 4/26/2023 initial commit
//////////////////////////////////////

// define namespace
#ifndef RANDOM_ROUTINES
#define RANDOM_ROUTINES

// cpp built in random lib
#include <random>
// ROOT MT19937 implementation
#include <TRandom3.h>
// ROOT lorentz vector
#include <TLorentzVector.h>
// GSL random lib
#include <gsl/gsl_rng.h>
// GSL random distribution
#include <gsl/gsl_randist.h>

// function prototypes
namespace Random_routines {
// the following routines requires a seed to be passed in for random number generation 
// if you don't care about seeding, please use 0 for the seed. If so, a random seed will 
// be generated 

// return random integer in the range [lower_bound, upper_bound]
int get_random_int_std(int lower_bound, int upper_bound, unsigned long int seed);
int get_random_int_root(int lower_bound, int upper_bound, unsigned long int seed);
int get_random_int_gsl(int lower_bound, int upper_bound, unsigned long int seed);

// return random double in the range [lower_bound, upper_bound]
double get_random_double_std(double lower_bound, double upper_bound, unsigned long int seed);
double get_random_double_root(double lower_bound, double upper_bound, unsigned long int seed);
double get_random_double_gsl(double lower_bound, double upper_bound, unsigned long int seed);

// return random double based on Gaussian distribution
// with mu = mean, sigma = stdev
double get_random_gaussian_std(double mean, double stdev, unsigned long int seed);
double get_random_gaussian_root(double mean, double stdev, unsigned long int seed);
double get_random_gaussian_gsl(double mean, double stdev, unsigned long int seed);

// return random double based on flat distribution
// in the range [min, max]
double get_random_flat_std(double min, double max, unsigned long int seed);
double get_random_flat_root(double min, double max, unsigned long int seed);
double get_random_flat_gsl(double min, double max, unsigned long int seed);

// return a Lorentz vector ptr with a fixed mass, random transverse momentum
// Pt with a range, random Peta with a range, random phi with a range
TLorentzVector *get_uniform_randomized_lorentz_ptr(double mass, double Pt_min,
                                           double Pt_max, double Peta_min,
                                           double Peta_max, double phi_min,
                                           double phi_max, unsigned long int seed);
}; // namespace Random_routines
#endif