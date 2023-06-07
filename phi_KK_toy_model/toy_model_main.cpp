#include "Random_routines.h"
#include "Selection_routines.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <cmath>
#include <fstream>
#include <iostream>

thread_local TRandom3 rng_toy(0);

// add pion uniform distibution between 0.5 GeV to 2* RHo mass 
#define RHO_MASS  rng_toy.Uniform(0.3, 2. * 0.77)
// double RHO_MASS = 0.77;
double PHI_MASS = 1.019;
double KAON_MASS = 0.493;
double PION_MASS = 0.139;
double PT_MIN = 0.1;
double PT_MAX = 10.;
double ETA_MIN = 0.;
double ETA_MAX = 4.;
double PHI_MIN = 0.;
double PHI_MAX = 2 * M_PI;
int PHI_SAMPLE_SIZE = 4000;
int RHO_SAMPLE_SIZE = 10 * PHI_SAMPLE_SIZE;

int main(int argc, char **argv) {
  // create vectors for parent vector and daughter vectors
  std::vector<TLorentzVector *> parent_vector;
  std::vector<TLorentzVector *> daughter1_vector;
  std::vector<TLorentzVector *> daughter2_vector;

  // simulate the decay process for rho -> pi+ pi-
  for (int i = 0; i < RHO_SAMPLE_SIZE; i++) {
    // generate parent vector that is a 4 vector uniformly distributed between
    // pt, eta, phi bounds with rho mass
    TLorentzVector *rho_parent_particle_ptr =
        Random_routines::get_random_lorentz_vector(
            PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, PHI_MIN, PHI_MAX, RHO_MASS);
    parent_vector.push_back(rho_parent_particle_ptr);

    // simulate the decay process with the daughter mass being pion
    std::vector<TLorentzVector *> daughter_ptr_pair =
        Random_routines::symmetrical_two_body_decay(rho_parent_particle_ptr,
                                                    PION_MASS);

    // blur the daughter particles by 2 percent to simulate actual particle
    // detector accuracy
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                           0.04 * daughter_ptr_pair[0]->Pt());
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                           0.04 * daughter_ptr_pair[1]->Pt());

    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);

    // now change the daughter mass to be Kaon- the assumed particle when doing
    // Kaon selections
    daughter1_vector[i]->SetPtEtaPhiM(daughter1_vector[i]->Pt(),
                                      daughter1_vector[i]->Eta(),
                                      daughter1_vector[i]->Phi(), KAON_MASS);
    daughter2_vector[i]->SetPtEtaPhiM(daughter2_vector[i]->Pt(),
                                      daughter2_vector[i]->Eta(),
                                      daughter2_vector[i]->Phi(), KAON_MASS);
  }

  // simualte decay process for phi -> K+ K-
  for (int i = 0; i < PHI_SAMPLE_SIZE; i++) {
    // generate parent vector that is a 4 vector uniformly distributed between
    // pt, eta, phi bounds with phi mass
    TLorentzVector *phi_parent_particle_ptr =
        Random_routines::get_random_lorentz_vector(
            PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, PHI_MIN, PHI_MAX, PHI_MASS);
    parent_vector.push_back(phi_parent_particle_ptr);

    // simulate the decay process with the daughter mass being kaon
    std::vector<TLorentzVector *> daughter_ptr_pair =
        Random_routines::symmetrical_two_body_decay(phi_parent_particle_ptr,
                                                    KAON_MASS);

    // blur the daughter particles by 2 percent to simulate actual particle
    // detector accuracy
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                           0.04 * daughter_ptr_pair[0]->Pt());
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                           0.04 * daughter_ptr_pair[1]->Pt());
    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);
  }

  // create histogram to draw rc mass
  TH1F *parent_rc_mass = new TH1F(
      "Combined Masses", "Toy Model Combined Masses;m_{K^+ K^-}(GeV);count",
      100, 0.8, 2.5);

  // create an instance of particle selector
  Selector pid = Selector();

  // select daughter particle to be Kaon
  for (int i = 0; i < parent_vector.size(); i++) {
    if (std::abs(pid.get_NSigmaKaon(daughter1_vector[i])) < 5. &&
        std::abs(pid.get_NSigmaKaon(daughter2_vector[i])) < 5. &&
        std::abs(pid.get_NSigmaPion(daughter1_vector[i])) > 5. &&
        std::abs(pid.get_NSigmaPion(daughter2_vector[i])) > 5. &&
        daughter1_vector[i]->Pt() > 0.06 && daughter2_vector[i]->Pt() > 0.06) {
      // reconstruct the dauther particles only if they are kaons.
      // This effectively selects phi
      TLorentzVector reconstructed_parent =
          *(daughter1_vector[i]) + *(daughter2_vector[i]);
      // fill the reconstructed parent mass
      parent_rc_mass->Fill(reconstructed_parent.M());
    }
  }

  // drawing the result
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  parent_rc_mass->Draw();
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  app.Run();
  return 0;
}
