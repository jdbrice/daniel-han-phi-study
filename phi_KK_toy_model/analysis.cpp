#include "Random_routines.h"
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

double PHI_MASS = 1.019;
double KAON_MASS = 0.493;
double PT_MIN = 0.;
double PT_MAX = 2.;
double ETA_MIN = 0.;
double ETA_MAX = 0.;
double PHI_MIN = 0.;
double PHI_MAX = 2 * M_PI;
double SAMPLE_SIZE = 2e6;

int main(int argc, char **argv) {
  std::vector<TLorentzVector *> parent_vector;
  std::vector<TLorentzVector *> daughter1_vector;
  std::vector<TLorentzVector *> daughter2_vector;
  TApplication app("app", &argc, argv);

  for (int i = 0; i < SAMPLE_SIZE; i++) {
    TLorentzVector *parent_particle_ptr =
        Random_routines::get_random_lorentz_vector(
            PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, PHI_MIN, PHI_MAX, PHI_MASS);
    parent_vector.push_back(parent_particle_ptr);
    std::vector<TLorentzVector *> daughter_ptr_pair =
        Random_routines::symmetrical_two_body_decay(parent_particle_ptr,
                                                    KAON_MASS);
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0]);
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1]);
    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);
  }

  TH1F *combined_masses = new TH1F(
      "Combined Masses", "Toy Model Combined Masses;m_{K^+ K^-}(GeV);count",
      500, 0.8, 1.3);
  for (int i = 0; i < SAMPLE_SIZE; i++) {
    TLorentzVector reconstructed_parent_vector =
        *daughter1_vector[i] + *daughter2_vector[i];

    if (daughter1_vector[i]->Pt() > 0.06 && daughter2_vector[i]->Pt() > 0.06) {
      combined_masses->Fill(reconstructed_parent_vector.M());
    }
  }
  gROOT->SetBatch(0);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  combined_masses->Draw();
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  app.Run();
  return 0;
}
