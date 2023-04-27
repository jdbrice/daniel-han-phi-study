#include "Random_routines.h"
#include <TLorentzVector.h>
#include <iomanip>
#include <math.h>
#include <ostream>
#include <vector>
#include <cmath>
#include <fstream>

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

  int test_int = 1;
  std::ofstream fail_data;
  fail_data.open("fail_data.csv");
  fail_data << "K1 Pt" << std::setw(26) << " " 
            << "K2 Pt" << std::setw(27) << " " 
            << "Phi Pt" << std::setw(25) << " " 
            << "Delta Pt" << std::endl; 
  int fail_counter = 0;
  do{
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


    double delta_pt = fabs(phi_ptr->Pt() - (kaon1_ptr->Pt() + kaon2_ptr->Pt()));
    // Output kaon pT, Eta, Phi, and M for high delta_pt
    if(delta_pt/phi_ptr->Pt() > 0.01){

      fail_data << std::scientific << std::setprecision(10) 
                << kaon1_ptr->Pt() << std::setw(15) << " "
                << kaon2_ptr->Pt() << std::setw(15) << " "
                << phi_ptr->Pt() << std::setw(15) << " "
                << delta_pt << std::endl;
      fail_counter++;
    }

    test_int++;
  }
  while (test_int < MAX_ITREATION);
  std::cout << "There are " << fail_counter<< " cases out of " << MAX_ITREATION << " that filed the one percent error test :-(" << std::endl;
  std::cout << "Failed data are stored in \"fail_data.csv \"" << std::endl;
  fail_data.close();
  return 0;




}

