#include "Random_routines.h"
#include <TLorentzVector.h>
#include <math.h>
#include <vector>
#include <cmath>

double PI = M_PI;
double MAX_ITREATION = 10000;
double PT_MIN = 0.;
double PT_MAX = 3.;
double ETA_MIN = 0.;
double ETA_MAX = 4.;
double M_KAON = 0.493;
double M_PHI = 1.019;

std::vector<TLorentzVector*> * phi_decay(TLorentzVector* phi_ptr);

int main(){
  // test run for phi--> K+ K- decay process
  // THIS IS TO BE FANILIZED INTO A FUNCTION! TESTING ONLY
    
  TLorentzVector* phi_ptr = (Random_routines::get_uniform_randomized_lorentz_ptr(M_PHI, PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, 0., 2. * PI, 0));
  TLorentzVector* kaon1_ptr = new TLorentzVector;
  TLorentzVector* kaon2_ptr = new TLorentzVector;

  // calculate kaon momentum in phi meson REST FRAME
  // in Phi rest frame, we have Ek1 = Ek2. This is because the rest mass
  // of Kaons are the same, and in Phi rest frame, we must have momentum conservation
  // so that the two daughter KAON has the same magnitude of momentum 

  // daughter particle angle generation
  double kaon1_polar = Random_routines::get_random_flat_std(0., 2. * PI, 0);
  double kaon1_azimuth = Random_routines::get_random_flat_std(0., 2. * PI, 0);

  // daughter particle momentum generaion
  double kaon1_momentum = sqrt((M_PHI * M_PHI - 4. * M_KAON * M_KAON) / 4.);
  double kaon1_pz = kaon1_momentum * cos(kaon1_polar);
  double kaon1_px = kaon1_momentum * sin(kaon1_polar) * cos(kaon1_azimuth);
  double kaon1_py = kaon1_momentum * sin(kaon1_polar) * sin(kaon1_azimuth);

  kaon1_ptr->SetPxPyPzE(kaon1_px, kaon1_py, kaon1_pz, sqrt(M_KAON * M_KAON + kaon1_px * kaon1_px + kaon1_py * kaon1_py + kaon1_pz * kaon1_pz));

  // trivial step to get kaon2 four vector once we determine phi and kaon1
  kaon2_ptr->SetPxPyPzE(-kaon1_px, -kaon1_py, -kaon1_pz, sqrt(M_KAON * M_KAON + kaon1_px * kaon1_px + kaon1_py * kaon1_py + kaon1_pz * kaon1_pz));
  // define the boost vector to be the four vector of the phi meson 
  // and lab_frame is staionary frame with respect to phi meson
  TVector3 boost_lab_frame = phi_ptr->BoostVector();

  // boost the daughter particles from phi meson frame to lab frame
  kaon1_ptr->Boost(boost_lab_frame);
  kaon2_ptr->Boost(boost_lab_frame);

  
  // Output kaon pT, Eta, Phi, and M
  std::cout << "Kaon 1: pT = " << kaon1_ptr->Pt() << ", Eta = " << kaon1_ptr->Eta() << ", Phi = " << kaon1_ptr->Phi() << ", M = " << kaon1_ptr->M() << std::endl;
  std::cout << "Kaon 2: pT = " << kaon2_ptr->Pt() << ", Eta = " << kaon2_ptr->Eta() << ", Phi = " << kaon2_ptr->Phi() << ", M = " << kaon2_ptr->M() << std::endl;
  std::cout << "Phi: pT = " << phi_ptr->Pt() << ", Eta = " << phi_ptr->Eta() << ", Phi = " << phi_ptr->Phi() << ", M = " << phi_ptr->M() << std::endl;
  std::cout << "Delta pT = " << fabs(phi_ptr->Pt() - (kaon1_ptr->Pt() + kaon2_ptr->Pt())) << std::endl ;

  return 0;




}

