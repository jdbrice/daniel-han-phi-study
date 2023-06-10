#include "Random_routines.h"
#include <TLorentzVector.h>
#include <TROOT.h>
#include <cmath>
#include <iostream>

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

// return a vector that contain the equal mass two-particle decay process from a
// given parent particle lorentz vector. The address of both dauthers are
// contained in the vector.

std::vector<TLorentzVector *>
Random_routines::symmetrical_two_body_decay(TLorentzVector *parent_particle,
                                            double daughter_mass) {
  TLorentzVector *daughter_one = new TLorentzVector;
  TLorentzVector *daughter_two = new TLorentzVector;
  std::vector<TLorentzVector *> daughter_list;
  daughter_list.push_back(daughter_one);
  daughter_list.push_back(daughter_two);

  // randomly generate polar angle for particle oen
  double kaon1_polar = rng.Uniform(0., pi);
  // randomly generate azimuthal angle for particle one
  double kaon1_azimuth = rng.Uniform(0., 2. * pi);

  // calculate the rest frame momentum for one particle using
  // special relativity. Note that this calculation is true due to the
  // conservation of momentum in the rest frame of the parent particle.
  double kaon1_momentum =
      std::sqrt((parent_particle->M() * parent_particle->M() -
                 4 * daughter_mass * daughter_mass) /
                4.);

  // calculate the four vector of the first daughter particle
  double kaon1_px =
      kaon1_momentum * std::sin(kaon1_polar) * std::cos(kaon1_azimuth);

  double kaon1_py =
      kaon1_momentum * std::sin(kaon1_polar) * std::sin(kaon1_azimuth);

  double kaon1_pz = kaon1_momentum * std::cos(kaon1_polar);

  double kaon1_E =
      std::sqrt(daughter_mass * daughter_mass + kaon1_px * kaon1_px +
                kaon1_py * kaon1_py + kaon1_pz * kaon1_pz);

  daughter_one->SetPxPyPzE(kaon1_px, kaon1_py, kaon1_pz, kaon1_E);
  // due to momentum conservation, the second daughter particle will have the
  // exact opposite signed momentum compared with the first daughter particle
  // with the same amount of energy.
  daughter_two->SetPxPyPzE(-kaon1_px, -kaon1_py, -kaon1_pz, kaon1_E);

  // we boost the daughter particles to the lab frame, where the parent particle
  // is moving
  TVector3 boost_vector = parent_particle->BoostVector();
  daughter_one->Boost(boost_vector);
  daughter_two->Boost(boost_vector);
  return daughter_list;
}

// return a vector that contain two-particle with unequal masses decay process from a
// given parent particle lorentz vector. The address of both dauthers are
// contained in the vector.
std::vector<TLorentzVector *>
Random_routines::two_body_decay(TLorentzVector *parent_particle,
                                            double m1, double m2) {
  TLorentzVector *daughter_one = new TLorentzVector;
  TLorentzVector *daughter_two = new TLorentzVector;
  std::vector<TLorentzVector *> daughter_list;
  daughter_list.push_back(daughter_one);
  daughter_list.push_back(daughter_two);

  double M = parent_particle->M();
  // randomly generate polar angle for particle oen
  double kaon1_polar = rng.Uniform(0., pi);
  // randomly generate azimuthal angle for particle one
  double kaon1_azimuth = rng.Uniform(0., 2. * pi);

  // calculate the rest frame momentum for one particle using
  // special relativity. Note that this calculation is true due to the
  // conservation of momentum in the rest frame of the parent particle plus energy conservation.
  double kaon1_momentum = sqrt((M*M - (m1 + m2)*(m1 + m2)) * (M*M - (m1 - m2)*(m1 - m2))) / (2.*M);
  // calculate the four vector of the first daughter particle
  double kaon1_px =
      kaon1_momentum * std::sin(kaon1_polar) * std::cos(kaon1_azimuth);

  double kaon1_py =
      kaon1_momentum * std::sin(kaon1_polar) * std::sin(kaon1_azimuth);

  double kaon1_pz = kaon1_momentum * std::cos(kaon1_polar);

  double kaon1_E =
      std::sqrt(m1 * m1 + kaon1_momentum * kaon1_momentum);

  double kaon2_E =
      std::sqrt(m2 * m2 + kaon1_momentum * kaon1_momentum);

  daughter_one->SetPxPyPzE(kaon1_px, kaon1_py, kaon1_pz, kaon1_E);
  // due to momentum conservation, the second daughter particle will have the
  // exact opposite signed momentum compared with the first daughter particle.
  daughter_two->SetPxPyPzE(-kaon1_px, -kaon1_py, -kaon1_pz, kaon2_E);

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
