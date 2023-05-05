#include <fstream>
#include <iostream>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <cmath>
#include <TRandom3.h>

double PHI_MASS = 1.019;
double KAON_MASS = 0.493;
double PT_MIN = 0.;
double PT_MAX = 2.;
double ETA_MIN = 0.;
double ETA_MAX = 0.;
double PHI_MIN = 0.; 
double PHI_MAX = 2 * M_PI;
double SAMPLE_SIZE = 2e6;
double pi = M_PI;
thread_local TRandom3 rng(0);

int ican = 0;
void makeCan() {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican++ ), "", 900, 600 );
    can->SetTopMargin(0.06);
    can->SetRightMargin(0.01);
}


// return a randomly parepared particle lorentz vector ptr with known mass
// and bounded pt-> transverse momentum, eta-> pseudorapidity, phi-> transverse
// azimuthal angle
TLorentzVector * get_random_lorentz_vector(
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
symmetrical_two_body_decay(TLorentzVector *parent_particle,
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


void add_uniform_pt_error(TLorentzVector *target_vector,
                                           double percent_error) {
  double max_error = percent_error / 100. * target_vector->Pt();
  double random_error = rng.Uniform(-max_error, max_error);
  target_vector->SetPtEtaPhiE(target_vector->Pt() + random_error,
                              target_vector->Eta(), target_vector->Phi(),
                              target_vector->E());
}

void add_gaussian_pt_error(TLorentzVector *target_vector,
                                           double percent_error) {
  double stdev = percent_error / 100. * target_vector->Pt();
  double random_error = rng.Gaus(0., stdev);
  target_vector->SetPtEtaPhiE(target_vector->Pt() + random_error,
                              target_vector->Eta(), target_vector->Phi(),
                              target_vector->E());
}

void analysis()
{
    std::vector<TLorentzVector*> parent_vector;
    std::vector<TLorentzVector*> daughter1_vector;
    std::vector<TLorentzVector*> daughter2_vector;

    for (int i = 0; i < SAMPLE_SIZE; i++)
    {
        TLorentzVector* parent_particle_ptr =  get_random_lorentz_vector(PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, PHI_MIN, PHI_MAX, PHI_MASS);
        parent_vector.push_back(parent_particle_ptr);
        std::vector<TLorentzVector*> daughter_ptr_pair =  symmetrical_two_body_decay(parent_particle_ptr, KAON_MASS);
        add_gaussian_pt_error(daughter_ptr_pair[0], 10);
        add_gaussian_pt_error(daughter_ptr_pair[1], 10);
        daughter1_vector.push_back(daughter_ptr_pair[0]);
        daughter2_vector.push_back(daughter_ptr_pair[1]);

    }
                        
    TH1F* combined_masses = new TH1F("Combined Masses", "Toy Model Combined Masses;m_{K^+ K^-}(GeV);count", 500, 0.8, 1.3);
    for (int i = 0; i < SAMPLE_SIZE; i++)
    {
       TLorentzVector reconstructed_parent_vector = *daughter1_vector[i] + *daughter2_vector[i]; 
        
       if(daughter1_vector[i]->Pt() > 0.06 && daughter2_vector[i]->Pt() > 0.06){
            combined_masses->Fill(reconstructed_parent_vector.M());
        }
    }


    makeCan();
    combined_masses->Draw();
}
