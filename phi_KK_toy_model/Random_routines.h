#include "TLorentzVector.h"
#include "TRandom3.h"
#include <vector>

#ifndef RANDOM_ROUTINES
// define namespace
#define RANDOM_ROUTINES
namespace Random_routines {

// return a Lorentz vector with randomized four-vectors with
// bounded Pt, eta, phi and a known mass in GeV
TLorentzVector *get_random_lorentz_vector(double Pt_min, double Pt_max,
                                          double eta_min, double eta_max,
                                          double phi_min, double phi_max,
                                          double mass);
// simulation random decay process from a parent particle with known mass to two
// identical daughte particles with known mass.
//
// Funtion return a list of lorentz vectors with the pointer of the two dauthers
std::vector<TLorentzVector *>
symmetrical_two_body_decay(TLorentzVector *parent_particle,
                           double daughter_mass);

// Add uniform error from zero to max_error 
void add_uniform_pt_error(TLorentzVector* target_vector, double max_error);

// Add gaussian error with mean 0, stdev = stdev 
// to the transverse momentum of 
// a given TLorentzVector. 
void add_gaussian_pt_error(TLorentzVector* target_vector, double stdev);

// Add gaussian error with mean 0, stdev = stdev 
// to the pseudorapidity of 
// a given TLorentzVector. 
void add_gaussian_eta_error(TLorentzVector* target_vector, double stdev);

// Add gaussian error with mean 0, stdev = stdev 
// to the transverse plane angle of 
// a given TLorentzVector. 
void add_gaussian_phi_error(TLorentzVector* target_vector, double stdev);
// set a given target 4-vector's transverse momentum uniformly distributed 
// between pt_min to pt_max
void set_uniform_pt(TLorentzVector* target_vector, double pt_min, double pt_max);

// set a given target 4-vector's transverse momentum with gaussian distribution 
// with mean = pt_mu, stdev = pt_sigma
void set_gauss_pt(TLorentzVector* target_vector, double pt_mu, double pt_sigma);
};

#endif // !RANDOM_ROUTINES
