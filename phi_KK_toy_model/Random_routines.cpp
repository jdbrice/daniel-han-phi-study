#include "Random_routines.h"
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <cmath>
#include <iostream>
#include <string>

double pi = M_PI;
thread_local TRandom3 rng(0);

// return a randomly parepared particle lorentz vector ptr with known mass
// and bounded pt-> transverse momentum, eta-> pseudorapidity, phi-> transverse
// azimuthal angle
TLorentzVector *Random_routines::get_random_lorentz_vector(
    double pt_min, double pt_max, double eta_min, double eta_max,
    double phi_min, double phi_max, double mass) {
  TLorentzVector *random_vector_ptr = new TLorentzVector;
  // get random pt from given bounds
  double pt = rng.Uniform(pt_min, pt_max);
  // get random eta from given bounds
  double eta = rng.Uniform(eta_min, eta_max);
  // get random phi from given bounds
  double phi = rng.Uniform(phi_min, phi_max);

  random_vector_ptr->SetPtEtaPhiM(pt, eta, phi, mass);

  return random_vector_ptr;
}

// return a Lorentz vectors with four-momentum that follows the distribution
// result of a starlight simulation histogram loaded as TFile Note that
// currently only Kaon as "kaon" in particle_type and Pion as "pion" in
// particle_type are supported. This is to be changed. Note that the root file
// must contain histogram with hardcoded four-momentum histogram names. Namely,
// they are: AzimuthX InvMassX PtX RapX Where X = El for electrons, Mu for
// Muons, Pi for pions, and Ka for kaons
TLorentzVector *Random_routines::get_slight_lorentz_vector(TFile *root_file,
                                          std::string particle_type) {
  TLorentzVector *rc_vecor = new TLorentzVector();
  if (particle_type == "kaon") {
    TH1F *phi_distr = (TH1F *)root_file->Get("AzimuthKa");
    TH1F *inv_mass_distr = (TH1F *)root_file->Get("InvMassKa");
    TH1F *pt_distr = (TH1F *)root_file->Get("PtKa");
    TH1F *eta_distr = (TH1F *)root_file->Get("RapKa");
    rc_vecor->SetPtEtaPhiM(pt_distr->GetRandom(), eta_distr->GetRandom(),
                           phi_distr->GetRandom(), inv_mass_distr->GetRandom());

  } else if (particle_type == "pion") {
    TH1F *phi_distr = (TH1F *)root_file->Get("AzimuthPi");
    TH1F *inv_mass_distr = (TH1F *)root_file->Get("InvMassPi");
    TH1F *pt_distr = (TH1F *)root_file->Get("PtPi");
    TH1F *eta_distr = (TH1F *)root_file->Get("RapPi");
    rc_vecor->SetPtEtaPhiM(pt_distr->GetRandom(), eta_distr->GetRandom(),
                           phi_distr->GetRandom(), inv_mass_distr->GetRandom());
  } else {
    std::cout << "this particle type is not supported!" << std::endl;
  }
  return rc_vecor;
}

// return a vector that contain two-particle with unequal masses decay process
// from a given parent particle lorentz vector. The address of both dauthers are
// contained in the vector.
std::vector<TLorentzVector *>
Random_routines::two_body_decay(TLorentzVector *parent_particle, double m1,
                                double m2) {
  TLorentzVector *daughter_one = new TLorentzVector;
  TLorentzVector *daughter_two = new TLorentzVector;
  std::vector<TLorentzVector *> daughter_list;
  daughter_list.push_back(daughter_one);
  daughter_list.push_back(daughter_two);

  double M = parent_particle->M();
  // randomly generate polar angle for particle oen
  double daughter1_cos_polar = rng.Uniform(-1., 1.);
  // randomly generate azimuthal angle for particle one
  double daughter1_azimuth = rng.Uniform(0., 2. * pi);

  // calculate the rest frame momentum for one particle using
  // special relativity. Note that this calculation is true due to the
  // conservation of momentum in the rest frame of the parent particle plus
  // energy conservation.
  double daughter1_momentum =
      sqrt((M * M - (m1 + m2) * (m1 + m2)) * (M * M - (m1 - m2) * (m1 - m2))) /
      (2. * M);
  // calculate the four vector of the first daughter particle
  double daughter1_px = daughter1_momentum *
                        sqrt(1. - daughter1_cos_polar * daughter1_cos_polar) *
                        std::cos(daughter1_azimuth);

  double daughter1_py = daughter1_momentum *
                        sqrt(1. - daughter1_cos_polar * daughter1_cos_polar) *
                        std::sin(daughter1_azimuth);

  double daughter1_pz = daughter1_momentum * daughter1_cos_polar;

  double daughter1_E =
      std::sqrt(m1 * m1 + daughter1_momentum * daughter1_momentum);

  double daughter2_E =
      std::sqrt(m2 * m2 + daughter1_momentum * daughter1_momentum);

  daughter_one->SetPxPyPzE(daughter1_px, daughter1_py, daughter1_pz,
                           daughter1_E);
  // due to momentum conservation, the second daughter particle will have the
  // exact opposite signed momentum compared with the first daughter particle.
  daughter_two->SetPxPyPzE(-daughter1_px, -daughter1_py, -daughter1_pz,
                           daughter2_E);

  // we boost the daughter particles to the lab frame, where the parent particle
  // is moving
  TVector3 boost_vector = parent_particle->BoostVector();
  daughter_one->Boost(boost_vector);
  daughter_two->Boost(boost_vector);
  return daughter_list;
}

void Random_routines::add_uniform_pt_error(TLorentzVector *target_vector,
                                           double max_error) {
  double random_error = rng.Uniform(-max_error, max_error);

  target_vector->SetPtEtaPhiM(target_vector->Pt() + random_error,
                              target_vector->Eta(), target_vector->Phi(),
                              target_vector->M());
}
void Random_routines::set_uniform_pt(TLorentzVector *target_vector,
                                     double pt_min, double pt_max) {
  target_vector->SetPtEtaPhiM(rng.Uniform(pt_min, pt_max), target_vector->Eta(),
                              target_vector->Phi(), target_vector->M());
}

void Random_routines::set_gauss_pt(TLorentzVector *target_vector, double pt_mu,
                                   double pt_sigma) {

  target_vector->SetPtEtaPhiM(rng.Gaus(pt_mu, pt_sigma), target_vector->Eta(),
                              target_vector->Phi(), target_vector->M());
}

void Random_routines::add_gaussian_pt_error(TLorentzVector *target_vector,
                                            double stdev) {
  double random_error = rng.Gaus(0., stdev);
  target_vector->SetPtEtaPhiM(target_vector->Pt() + random_error,
                              target_vector->Eta(), target_vector->Phi(),
                              target_vector->M());
}
void Random_routines::add_gaussian_eta_error(TLorentzVector *target_vector,
                                             double stdev) {
  target_vector->SetPtEtaPhiM(target_vector->Pt(),
                              target_vector->Eta() + rng.Gaus(0., stdev),
                              target_vector->Phi(), target_vector->M());
}

void Random_routines::add_gaussian_phi_error(TLorentzVector *target_vector,
                                             double stdev) {
  target_vector->SetPtEtaPhiM(
      target_vector->Pt(), target_vector->Eta(),
      std::fmod((target_vector->Phi() + rng.Gaus(0., stdev)), (2. * pi)),
      target_vector->M());
}
