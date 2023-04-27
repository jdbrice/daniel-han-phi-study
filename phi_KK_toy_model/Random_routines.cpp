#include "Random_routines.h"
#include <TLorentzVector.h>
#include <iostream>
#include <random>

// random number generator seed
// This is used when the user enter 0 as seed
using std::endl;

std::random_device rd;
std::mt19937_64 global_rng(rd());

// function definitions see header file "Random_routines.h"
int Random_routines::get_random_int_std(int lower_bound, int upper_bound,
                                        unsigned long int seed) {
  std::uniform_int_distribution<int> dist(lower_bound, upper_bound + 1);
  if (seed == 0) {
    return dist(global_rng);
  }
  std::mt19937_64 local_rng(seed);
  return dist(local_rng);
}

double Random_routines::get_random_double_std(double lower_bound,
                                              double upper_bound,
                                              unsigned long int seed) {
  std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
  if (seed == 0) {
    return dist(global_rng);
  }
  std::mt19937_64 local_rng(seed);
  return dist(local_rng);
}

double Random_routines::get_random_gaussian_std(double mean, double stdev,
                                                unsigned long seed) {
  std::normal_distribution<double> dist(mean, stdev);
  if (seed == 0) {
    return dist(global_rng);
  }
  std::mt19937_64 local_rng(seed);
  return dist(local_rng);
}

double Random_routines::get_random_flat_std(double lower_bound,
                                            double upper_bound,
                                            unsigned long seed) {
  std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
  if (seed == 0) {
    return dist(global_rng);
  }
  std::mt19937_64 local_rng(seed);
  return dist(local_rng);
}

TLorentzVector *Random_routines::get_uniform_randomized_lorentz_ptr(
    double mass, double Pt_min, double Pt_max, double Peta_min, double Peta_max,
    double phi_min, double phi_max, unsigned long seed) {

  if (seed == 0) {
    double pt = get_random_flat_std(Pt_min, Pt_max, 0);
    double peta = get_random_flat_std(Peta_min, Peta_max, 0);
    double phi = get_random_flat_std(phi_min, phi_max, 0);
    TLorentzVector *lorentz_vector_ptr = new TLorentzVector();
    lorentz_vector_ptr->SetPtEtaPhiM(pt, peta, phi, mass);
    return lorentz_vector_ptr;
  }
  double pt = get_random_flat_std(Pt_min, Pt_max, seed);
  double peta = get_random_flat_std(Peta_min, Peta_max, seed);
  double phi = get_random_flat_std(phi_min, phi_max, seed);
  TLorentzVector *lorentz_vector_ptr = new TLorentzVector();
  lorentz_vector_ptr->SetPtEtaPhiM(pt, peta, phi, mass);
  return lorentz_vector_ptr;
}
