#include "Random_routines.h"
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <cmath>
#include <iostream>
#include <string>

double pi = M_PI;
double KAON_MASS_RND = 0.493;
double PION_MASS_RND = 0.139;
double ELECTRON_MASS_RND = 0.000511;
// hardcoded starlight histogram files
TFile *kaon_file_rnd =
    new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/kaon_co.root");
TFile *kaon_inc_file =
    new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/kaon_inco.root");
TFile *pion_file =
    new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/pion_co.root");
TFile *pion_inc_file =
    new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/pion_inco.root");
TFile *electron_file =
    new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/electron.root");

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
TLorentzVector *
Random_routines::get_slight_lorentz_vector(TFile *root_file,
                                           std::string particle_type) {
  TLorentzVector *rc_vector = new TLorentzVector();
  if (particle_type == "kaon") {
    TH1F *phi_distr = (TH1F *)root_file->Get("AzimuthKa");
    TH1F *inv_mass_distr = (TH1F *)root_file->Get("InvMassKa");
    TH1F *pt_distr = (TH1F *)root_file->Get("PtKa");
    TH1F *eta_distr = (TH1F *)root_file->Get("EtaKa");
    rc_vector->SetPtEtaPhiM(pt_distr->GetRandom(), eta_distr->GetRandom(),
                            phi_distr->GetRandom(),
                            inv_mass_distr->GetRandom());
  } else if (particle_type == "pion") {
    TH1F *phi_distr = (TH1F *)root_file->Get("AzimuthPi");
    TH1F *inv_mass_distr = (TH1F *)root_file->Get("InvMassPi");
    TH1F *pt_distr = (TH1F *)root_file->Get("PtPi");
    TH1F *eta_distr = (TH1F *)root_file->Get("EtaPi");
    rc_vector->SetPtEtaPhiM(pt_distr->GetRandom(), eta_distr->GetRandom(),
                            phi_distr->GetRandom(),
                            inv_mass_distr->GetRandom());
  } else if (particle_type == "electron") {
    TH1F *phi_distr = (TH1F *)root_file->Get("AzimuthEl");
    TH1F *inv_mass_distr = (TH1F *)root_file->Get("InvMassEl");
    TH1F *pt_distr = (TH1F *)root_file->Get("PtEl");
    TH1F *eta_distr = (TH1F *)root_file->Get("EtaEl");
    rc_vector->SetPtEtaPhiM(pt_distr->GetRandom(), eta_distr->GetRandom(),
                            phi_distr->GetRandom(),
                            inv_mass_distr->GetRandom());

  } else {
    std::cout << "this particle type is not supported!" << std::endl;
  }
  return rc_vector;
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

void Random_routines::add_pt_percent_loss(TLorentzVector *target_vector,
                                          double energy_lost_percent) {
  double factor = 1. - energy_lost_percent / 100.;
  target_vector->SetPtEtaPhiM(target_vector->Pt() * factor,
                              target_vector->Eta(), target_vector->Phi(),
                              target_vector->M());
}

void Random_routines::trigger_two_body_decay(
    std::string type, int event_num, std::vector<TLorentzVector *> &mc_parent,
    std::vector<TLorentzVector *> &mc_d1, std::vector<TLorentzVector *> &mc_d2,
    std::vector<TLorentzVector *> &rc_parent,
    std::vector<TLorentzVector *> &rc_d1,
    std::vector<TLorentzVector *> &rc_d2) {

  int count = 0;
  if (type == "kaon") {
    while (count < event_num) {
      TLorentzVector *mc_parent_ptr =
          Random_routines::get_slight_lorentz_vector(kaon_file_rnd, "kaon");

      mc_parent.push_back(mc_parent_ptr);

      std::vector<TLorentzVector *> daughter_list =
          Random_routines::two_body_decay(mc_parent_ptr, KAON_MASS_RND, KAON_MASS_RND);

      mc_d1.push_back(daughter_list[0]);
      mc_d2.push_back(daughter_list[1]);

      TLorentzVector *d1_rc_ptr = new TLorentzVector;
      d1_rc_ptr->SetPtEtaPhiM(daughter_list[0]->Pt(), daughter_list[0]->Eta(),
                              daughter_list[0]->Phi(), daughter_list[0]->M());

      TLorentzVector *d2_rc_ptr = new TLorentzVector;
      d2_rc_ptr->SetPtEtaPhiM(daughter_list[1]->Pt(), daughter_list[1]->Eta(),
                              daughter_list[1]->Phi(), daughter_list[1]->M());

      Random_routines::add_pt_percent_loss(d1_rc_ptr, 2.5);
      Random_routines::add_pt_percent_loss(d2_rc_ptr, 2.5);
      Random_routines::add_gaussian_pt_error(d1_rc_ptr, 0.04 * d1_rc_ptr->Pt());
      Random_routines::add_gaussian_pt_error(d2_rc_ptr, 0.04 * d2_rc_ptr->Pt());

      // // DEBUG ONLY
      // double rand_pt = rng.Uniform(0.06,0.4);
      // d1_rc_ptr->SetPtEtaPhiM(rand_pt, d1_rc_ptr->Eta(), d1_rc_ptr->Phi(), d1_rc_ptr->M());
      // double rand_pt_2 = rng.Uniform(0.06,0.4);
      // d2_rc_ptr->SetPtEtaPhiM(rand_pt_2, d2_rc_ptr->Eta(), d2_rc_ptr->Phi(), d2_rc_ptr->M());

      rc_d1.push_back(d1_rc_ptr);
      rc_d2.push_back(d2_rc_ptr);

      TLorentzVector parent_rc = *d1_rc_ptr + *d2_rc_ptr;
      TLorentzVector *parent_rc_ptr = new TLorentzVector;
      parent_rc_ptr->SetPtEtaPhiM(parent_rc.Pt(), parent_rc.Eta(),
                                  parent_rc.Eta(), parent_rc.M());
      rc_parent.push_back(parent_rc_ptr);
      count++;
    }
  } else if (type == "pion") {
    while (count < event_num) {
      TLorentzVector *mc_parent_ptr =
          Random_routines::get_slight_lorentz_vector(pion_file, "pion");

      mc_parent.push_back(mc_parent_ptr);

      std::vector<TLorentzVector *> daughter_list =
          Random_routines::two_body_decay(mc_parent_ptr, PION_MASS_RND, PION_MASS_RND);

      mc_d1.push_back(daughter_list[0]);
      mc_d2.push_back(daughter_list[1]);

      TLorentzVector *d1_rc_ptr = new TLorentzVector;
      d1_rc_ptr->SetPtEtaPhiM(daughter_list[0]->Pt(), daughter_list[0]->Eta(),
                              daughter_list[0]->Phi(), KAON_MASS_RND);

      TLorentzVector *d2_rc_ptr = new TLorentzVector;
      d2_rc_ptr->SetPtEtaPhiM(daughter_list[1]->Pt(), daughter_list[1]->Eta(),
                              daughter_list[1]->Phi(), KAON_MASS_RND);

      Random_routines::add_pt_percent_loss(d1_rc_ptr, 2.5);
      Random_routines::add_pt_percent_loss(d2_rc_ptr, 2.5);
      Random_routines::add_gaussian_pt_error(d1_rc_ptr, 0.04 * d1_rc_ptr->Pt());
      Random_routines::add_gaussian_pt_error(d2_rc_ptr, 0.04 * d2_rc_ptr->Pt());

      rc_d1.push_back(d1_rc_ptr);
      rc_d2.push_back(d2_rc_ptr);

      TLorentzVector parent_rc = *d1_rc_ptr + *d2_rc_ptr;
      TLorentzVector *parent_rc_ptr = new TLorentzVector;
      parent_rc_ptr->SetPtEtaPhiM(parent_rc.Pt(), parent_rc.Eta(),
                                  parent_rc.Eta(), parent_rc.M());
      rc_parent.push_back(parent_rc_ptr);
      count++;
    }
  } else if (type == "electron") {
    while (count < event_num) {
      TLorentzVector *mc_parent_ptr =
          Random_routines::get_slight_lorentz_vector(electron_file, "electron");

      mc_parent.push_back(mc_parent_ptr);

      std::vector<TLorentzVector *> daughter_list =
          Random_routines::two_body_decay(mc_parent_ptr, ELECTRON_MASS_RND, ELECTRON_MASS_RND);

      mc_d1.push_back(daughter_list[0]);
      mc_d2.push_back(daughter_list[1]);

      TLorentzVector *d1_rc_ptr = new TLorentzVector;
      d1_rc_ptr->SetPtEtaPhiM(daughter_list[0]->Pt(), daughter_list[0]->Eta(),
                              daughter_list[0]->Phi(), KAON_MASS_RND);

      TLorentzVector *d2_rc_ptr = new TLorentzVector;
      d2_rc_ptr->SetPtEtaPhiM(daughter_list[1]->Pt(), daughter_list[1]->Eta(),
                              daughter_list[1]->Phi(), KAON_MASS_RND);

      Random_routines::add_pt_percent_loss(d1_rc_ptr, 2.5);
      Random_routines::add_pt_percent_loss(d2_rc_ptr, 2.5);
      Random_routines::add_gaussian_pt_error(d1_rc_ptr, 0.04 * d1_rc_ptr->Pt());
      Random_routines::add_gaussian_pt_error(d2_rc_ptr, 0.04 * d2_rc_ptr->Pt());

      rc_d1.push_back(d1_rc_ptr);
      rc_d2.push_back(d2_rc_ptr);

      TLorentzVector parent_rc = *d1_rc_ptr + *d2_rc_ptr;
      TLorentzVector *parent_rc_ptr = new TLorentzVector;
      parent_rc_ptr->SetPtEtaPhiM(parent_rc.Pt(), parent_rc.Eta(),
                                  parent_rc.Eta(), parent_rc.M());
      rc_parent.push_back(parent_rc_ptr);
      count++;
    }
  } else if (type == "kaon_inc") {
    while (count < event_num) {
      TLorentzVector *mc_parent_ptr =
          Random_routines::get_slight_lorentz_vector(kaon_inc_file, "kaon");

      mc_parent.push_back(mc_parent_ptr);

      std::vector<TLorentzVector *> daughter_list =
          Random_routines::two_body_decay(mc_parent_ptr, KAON_MASS_RND, KAON_MASS_RND);

      mc_d1.push_back(daughter_list[0]);
      mc_d2.push_back(daughter_list[1]);

      TLorentzVector *d1_rc_ptr = new TLorentzVector;
      d1_rc_ptr->SetPtEtaPhiM(daughter_list[0]->Pt(), daughter_list[0]->Eta(),
                              daughter_list[0]->Phi(), daughter_list[0]->M());

      TLorentzVector *d2_rc_ptr = new TLorentzVector;
      d2_rc_ptr->SetPtEtaPhiM(daughter_list[1]->Pt(), daughter_list[1]->Eta(),
                              daughter_list[1]->Phi(), daughter_list[1]->M());

      Random_routines::add_pt_percent_loss(d1_rc_ptr, 2.5);
      Random_routines::add_pt_percent_loss(d2_rc_ptr, 2.5);
      Random_routines::add_gaussian_pt_error(d1_rc_ptr, 0.04 * d1_rc_ptr->Pt());
      Random_routines::add_gaussian_pt_error(d2_rc_ptr, 0.04 * d2_rc_ptr->Pt());

      rc_d1.push_back(d1_rc_ptr);
      rc_d2.push_back(d2_rc_ptr);

      TLorentzVector parent_rc = *d1_rc_ptr + *d2_rc_ptr;
      TLorentzVector *parent_rc_ptr = new TLorentzVector;
      parent_rc_ptr->SetPtEtaPhiM(parent_rc.Pt(), parent_rc.Eta(),
                                  parent_rc.Eta(), parent_rc.M());
      rc_parent.push_back(parent_rc_ptr);
      count++;
    }

  } else if (type == "pion_inc") {
    while (count < event_num) {
      TLorentzVector *mc_parent_ptr =
          Random_routines::get_slight_lorentz_vector(pion_inc_file, "pion");

      mc_parent.push_back(mc_parent_ptr);

      std::vector<TLorentzVector *> daughter_list =
          Random_routines::two_body_decay(mc_parent_ptr, PION_MASS_RND, PION_MASS_RND);

      mc_d1.push_back(daughter_list[0]);
      mc_d2.push_back(daughter_list[1]);

      TLorentzVector *d1_rc_ptr = new TLorentzVector;
      d1_rc_ptr->SetPtEtaPhiM(daughter_list[0]->Pt(), daughter_list[0]->Eta(),
                              daughter_list[0]->Phi(), KAON_MASS_RND);

      TLorentzVector *d2_rc_ptr = new TLorentzVector;
      d2_rc_ptr->SetPtEtaPhiM(daughter_list[1]->Pt(), daughter_list[1]->Eta(),
                              daughter_list[1]->Phi(), KAON_MASS_RND);

      Random_routines::add_pt_percent_loss(d1_rc_ptr, 2.5);
      Random_routines::add_pt_percent_loss(d2_rc_ptr, 2.5);
      Random_routines::add_gaussian_pt_error(d1_rc_ptr, 0.04 * d1_rc_ptr->Pt());
      Random_routines::add_gaussian_pt_error(d2_rc_ptr, 0.04 * d2_rc_ptr->Pt());

      rc_d1.push_back(d1_rc_ptr);
      rc_d2.push_back(d2_rc_ptr);

      TLorentzVector parent_rc = *d1_rc_ptr + *d2_rc_ptr;
      TLorentzVector *parent_rc_ptr = new TLorentzVector;
      parent_rc_ptr->SetPtEtaPhiM(parent_rc.Pt(), parent_rc.Eta(),
                                  parent_rc.Eta(), parent_rc.M());
      rc_parent.push_back(parent_rc_ptr);
      count++;
    }
  } else {
    std::cout << "invalid decay type" << std::endl;
  }
}
