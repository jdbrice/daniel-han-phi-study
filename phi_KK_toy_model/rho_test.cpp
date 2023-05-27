#include "Random_routines.h"
#include "selection.h"
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

double RHO_MASS = 0.77;
double PHI_MASS = 1.019;
double KAON_MASS = 0.493;
double PION_MASS = 0.139;
double PT_MIN = 0.;
double PT_MAX = 1.;
double ETA_MIN = 0.;
double ETA_MAX = 4.;
double PHI_MIN = 0.;
double PHI_MAX = 2 * M_PI;
double PHI_SAMPLE_SIZE = 100;
double RHO_SAMPLE_SIZE = 1000;

int main(int argc, char **argv)
{
  std::vector<TLorentzVector *> parent_vector;
  std::vector<TLorentzVector *> daughter1_vector;
  std::vector<TLorentzVector *> daughter2_vector;
  TApplication app("app", &argc, argv);

  for (int i = 0; i < RHO_SAMPLE_SIZE; i++)
  {
    TLorentzVector *rho_parent_particle_ptr =
        Random_routines::get_random_lorentz_vector(
            PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, PHI_MIN, PHI_MAX, RHO_MASS);
    parent_vector.push_back(rho_parent_particle_ptr);
    std::vector<TLorentzVector *> daughter_ptr_pair =
        Random_routines::symmetrical_two_body_decay(rho_parent_particle_ptr,
                                                    PION_MASS);
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0], 0.01 * daughter_ptr_pair[0]->Pt());
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1], 0.01 * daughter_ptr_pair[1]->Pt());
    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);
    daughter1_vector[i]->SetPtEtaPhiM(daughter1_vector[i]->Pt(),
                                      daughter1_vector[i]->Eta(),
                                      daughter1_vector[i]->Phi(), KAON_MASS);
    daughter2_vector[i]->SetPtEtaPhiM(daughter2_vector[i]->Pt(),
                                      daughter2_vector[i]->Eta(),
                                      daughter2_vector[i]->Phi(), KAON_MASS);
  }

  for (int i = 0; i < PHI_SAMPLE_SIZE; i++)
  {
    TLorentzVector *phi_parent_particle_ptr =
        Random_routines::get_random_lorentz_vector(
            PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, PHI_MIN, PHI_MAX, PHI_MASS);
    parent_vector.push_back(phi_parent_particle_ptr);
    std::vector<TLorentzVector *> daughter_ptr_pair =
        Random_routines::symmetrical_two_body_decay(phi_parent_particle_ptr,
                                                    KAON_MASS);
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0], 0.01 * daughter_ptr_pair[0]->Pt());
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1], 0.01 * daughter_ptr_pair[1]->Pt());
    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);
  }

  TH1F *combined_masses = new TH1F(
      "Combined Masses", "Toy Model Combined Masses;m_{K^+ K^-}(GeV);count",
      100, 0., 4.);
  TH1F *daughter_masses = new TH1F(
      "Combined Masses", "Toy Model Combined Masses;m_{K^+ K^-}(GeV);count",
      100, 0., 4.);
  
  Selector PID = Selector();
  
  for (int i = 0; i < parent_vector.size(); i++){
    if (PID.get_NSigmaKaon(daughter1_vector[i]) <5 && PID.get_NSigmaPion(daughter1_vector[i])>5){
      TLorentzVector reconstructed_parent = *(daughter1_vector[i]) + *(daughter2_vector[i]);
      combined_masses->Fill(reconstructed_parent.M());
    }
  }

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
