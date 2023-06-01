#include "Random_routines.h"
#include "Selection_routines.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <cmath>
#include <fstream>
#include <iostream>

double RHO_MASS = 0.77;
double PHI_MASS = 1.019;
double KAON_MASS = 0.493;
double PION_MASS = 0.139;
double PT_MIN = 0.1;
double PT_MAX = 10.;
double ETA_MIN = 0.;
double ETA_MAX = 4.;
double PHI_MIN = 0.;
double PHI_MAX = 2 * M_PI;
int PHI_SAMPLE_SIZE = 1e6;

int main(int argc, char **argv) {
  // create vectors for parent vector and daughter vectors
  std::vector<TLorentzVector *> parent_vector;
  std::vector<TLorentzVector *> daughter1_vector;
  std::vector<TLorentzVector *> daughter2_vector;

  TH1F *mc_pt = new TH1F("parent mc phi", "MC P_T; P_T; Count", 500, 0., 1.5);

  // create histogram to draw rc mass
  TH2F *rc_mass_pt_selection = new TH2F(
      "pt vs mass",
      "Toy RC Pt v.s. Mass with 60MeV Pt Cut-off;m_{K^+ K^-}(GeV); P_T", 100,
      0.98, 1.1, 500, 0., 1.5);

  // simualte decay process for phi -> K+ K-
  for (int i = 0; i < PHI_SAMPLE_SIZE; i++) {
    // generate parent vector that is a 4 vector uniformly distributed between
    // pt, eta, phi bounds with phi mass
    TLorentzVector *phi_parent_particle_ptr =
        Random_routines::get_random_lorentz_vector(
            PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, PHI_MIN, PHI_MAX, PHI_MASS);
    parent_vector.push_back(phi_parent_particle_ptr);
    mc_pt->Fill(phi_parent_particle_ptr->Pt());
    // simulate the decay process with the daughter mass being kaon
    std::vector<TLorentzVector *> daughter_ptr_pair =
        Random_routines::symmetrical_two_body_decay(phi_parent_particle_ptr,
                                                    KAON_MASS);

    // blur the daughter particles by 2 percent to simulate actual particle
    // detector accuracy
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                           0.07 * daughter_ptr_pair[0]->Pt());
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                           0.07 * daughter_ptr_pair[1]->Pt());
    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);
  }



  // create an instance of particle selector
  Selector pid = Selector();

  // select daughter particle to be Kaon
  for (int i = 0; i < parent_vector.size(); i++) {

    if (std::abs(pid.get_NSigmaKaon(daughter1_vector[i])) < 5. &&
        std::abs(pid.get_NSigmaKaon(daughter2_vector[i])) < 5. &&
        std::abs(pid.get_NSigmaPion(daughter1_vector[i])) > 5. &&
        std::abs(pid.get_NSigmaPion(daughter2_vector[i])) > 5. &&
        daughter1_vector[i]->Pt() > 0.0 && daughter1_vector[i]->Pt() > 0.0) {

      TLorentzVector reconstructed_parent = (*daughter1_vector[i]) + (*daughter2_vector[i]);
      // fill the reconstructed parent mass
      rc_mass_pt_selection->Fill(reconstructed_parent.M(),
                                 reconstructed_parent.Pt());
    }
  }

  // drawing the result
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 1280, 720);

  TH1F *acceptance = (TH1F*) rc_mass_pt_selection->ProjectionY()->Clone();
  acceptance->SetXTitle("P_T");
  acceptance->SetYTitle("Reco Percent");
  acceptance->SetName("acceptance_rate");
  acceptance->SetTitle("Reco Acceptance Rate");
  acceptance->Divide(mc_pt);
  acceptance->Draw("C");
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  app.Run();
  return 0;
}
