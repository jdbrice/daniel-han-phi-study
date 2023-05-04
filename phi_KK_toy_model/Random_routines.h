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
std::vector<TLorentzVector *>
symmetrical_two_body_decay(TLorentzVector *parent_particle,
                           double daughter_mass);
}; // namespace Random_routines

// simulation random decay process from a parent particle with known mass to two
// identical daughte particles with known mass.
//
// Funtion return a list of lorentz vectors with the pointer of the two dauthers

#endif // !RANDOM_ROUTINES
