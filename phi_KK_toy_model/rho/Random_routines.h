#include "TLorentzVector.h"
#include "TRandom3.h"
#include <TUUID.h>
#include <string>
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

// return a Lorentz vectors with four-momentum that follows the distribution
// result of a starlight simulation histogram loaded as TFile Note that
// currently only Kaon as "kaon" in particle_type and Pion as "pion" in
// particle_type are supported. This is to be changed. Note that the root file
// must contain histogram with hardcoded four-momentum histogram names. Namely,
// they are: AzimuthX InvMassX PtX RapX Where X = El for electrons, Mu for
// Muons, Pi for pions, and Ka for kaons
TLorentzVector *get_slight_lorentz_vector(TFile *root_file,
                                          std::string particle_type);

// Funtion return a list of lorentz vectors with the pointer of the two dauthers
// when the two daughters have different masses
std::vector<TLorentzVector *> two_body_decay(TLorentzVector *parent_particle,
                                             double daughter_1_mass,
                                             double daughter_2_mass);

// Add uniform error from zero to max_error
void add_uniform_pt_error(TLorentzVector *target_vector, double max_error);

// Add gaussian error with mean 0, stdev = stdev
// to the transverse momentum of
// a given TLorentzVector.
void add_gaussian_pt_error(TLorentzVector *target_vector, double stdev);

// Add gaussian error with mean 0, stdev = stdev
// to the pseudorapidity of
// a given TLorentzVector.
void add_gaussian_eta_error(TLorentzVector *target_vector, double stdev);

// Add gaussian error with mean 0, stdev = stdev
// to the transverse plane angle of
// a given TLorentzVector.
void add_gaussian_phi_error(TLorentzVector *target_vector, double stdev);
// set a given target 4-vector's transverse momentum uniformly distributed
// between pt_min to pt_max
void set_uniform_pt(TLorentzVector *target_vector, double pt_min,
                    double pt_max);

// set a given target 4-vector's transverse momentum with gaussian distribution
// with mean = pt_mu, stdev = pt_sigma
void set_gauss_pt(TLorentzVector *target_vector, double pt_mu, double pt_sigma);

// add a certain amount of transverse momentum loss
void add_pt_percent_loss(TLorentzVector *target_vector,
                         double energy_lost_percent);

// perform one simulation of two body decay with detector effects
// type includes: kaon, pion, electron, kaon_inc, pion_inc. The mc TLorentzVector list and rc list needs to be provided.
void trigger_two_body_decay(std::string type, int event_num,
                            std::vector<TLorentzVector *> &mc_parent,
                            std::vector<TLorentzVector *> &mc_d1,
                            std::vector<TLorentzVector *> &mc_d2,
                            std::vector<TLorentzVector *> &rc_parent,
                            std::vector<TLorentzVector *> &rc_d1,
                            std::vector<TLorentzVector *> &rc_d2);
}; // namespace Random_routines

#endif // !RANDOM_ROUTINES
