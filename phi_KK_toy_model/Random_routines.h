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

// Add uniform error from zero to max_percent 
// to the transverse momentum of 
// a given TLorentzVector. 
void add_uniform_pt_error(TLorentzVector* target_vector, double percent_error);

// Add gaussian error with mean 0, stdev = percent_error  
// to the transverse momentum of 
// a given TLorentzVector. 
void add_gaussian_pt_error(TLorentzVector* target_vector);
}; // namespace Random_routines


#endif // !RANDOM_ROUTINES
