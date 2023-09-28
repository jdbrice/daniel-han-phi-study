#include "Random_routines.h"
#include "Selection_routines.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <cmath>
#include <fstream>
#include <iostream>

thread_local TRandom3 rng_toy(0);


double PHI_MASS = 1.019;
double KAON_MASS = 0.493;
double PT_MIN = 0.;
double PT_MAX = 3.;
double ETA_MIN = -1.;
double ETA_MAX = 1.;
double PHI_MIN = 0.;
double PHI_MAX = 2 * M_PI;
int PHI_SAMPLE_SIZE = 4500000;

  TFile *kaon_file =
      new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/kaon_co.root");

int main(int argc, char **argv) {
  // create vectors for parent vector and daughter vectors
  std::vector<TLorentzVector *> parent_vector;
  std::vector<TLorentzVector *> daughter1_vector;
  std::vector<TLorentzVector *> daughter2_vector;


  TH1F* phi_pt_rc_hist = new TH1F("Reco Effiency", "Reco Effiency;Phi P_T(GeV);Rate", 100, 0, 0.3 );

  TH1F* phi_pt_mc_hist = new TH1F("rc_phi_pt", "RC Phi Transverse Momenntum; P_T(GeV);Counts", 100, 0, 0.3);
  // simulate the decay process for rho -> pi+ pi-

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
        Random_routines::two_body_decay(phi_parent_particle_ptr,
                                                    KAON_MASS,KAON_MASS);

    // blur the daughter particles by 2 percent to simulate actual particle
    // detector accuracy
    // Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
    //                                        0.04 * daughter_ptr_pair[0]->Pt());
    // Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
    //                                        0.04 * daughter_ptr_pair[1]->Pt());
    //
    // Random_routines::add_pt_percent_loss(daughter_ptr_pair[0],3.);
    // Random_routines::add_pt_percent_loss(daughter_ptr_pair[1],3.);
    //
    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);

    phi_pt_mc_hist->Fill(phi_parent_particle_ptr->Pt());
  }



  // select daughter particle to be Kaon
  for (int i = 0; i < parent_vector.size(); i++) {
    if (daughter1_vector[i]->Pt() > 0.06 && daughter2_vector[i]->Pt() > 0.06)  {
      // reconstruct the dauther particles only if they are kaons.
      // This effectively selects phi
      TLorentzVector reconstructed_parent =
          *(daughter1_vector[i]) + *(daughter2_vector[i]);
      // fill the reconstructed parent mass
      phi_pt_rc_hist->Fill(reconstructed_parent.Pt());
    }
  }
  phi_pt_rc_hist->Divide(phi_pt_mc_hist);
  // drawing the result
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  phi_pt_rc_hist->Draw();
  TFile *file = new TFile("itpc_eff.root", "RECREATE");
  phi_pt_rc_hist->Write();
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  app.Run();
  return 0;
}
