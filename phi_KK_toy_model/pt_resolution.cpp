#include "Random_routines.h"
#include "Selection_routines.h"
#include <RtypesCore.h>
#include <TSystem.h>
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

double PHI_MASS = 1.019;
double KAON_MASS = 0.493;
double PT_MIN = 0.1;
double PT_MAX = 10.;
double ETA_MIN = 0.;
double ETA_MAX = 4.;
double PHI_MIN = 0.;
double PHI_MAX = 2 * M_PI;
double SAMPLE_SIZE = 3e4;

// create a instance of simulation with fixed percentage of pt blur and sample size
TH1F *get_simulation_pdf(double pt_blur_percent, int sample_size);

// return the p value from a chi2 calculation between a simulation and experiment 
double calculate_pt_blur_p_value(double pt_blur_percent, TH1F *experimental_pdf);

int main(int argc, char **argv) {
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  
  // load the processed phi experimental data. This comes from ../phi_selection_pt_resolution/pt_resolution_processor.C
  TFile* processed_experimental_data  = new TFile("processed_phi_exp_data.root");

  // load histogram from file 
  TH1F* experimental_pdf =(TH1F *) processed_experimental_data->Get("Phi Selected");

  // initial pt percent error and its corresponding p value from a chi2 test 
  double final_p_value = 0. ;
  double final_pt_percent_error = 0. ;

  // run simulation with pt percent error from 0 percent 1o 10 percent with step size 0.1  
  // this upper limit is testd experimentally 
  for (double pt_percent_error = 0.; pt_percent_error <= 10.; pt_percent_error += 0.1){

  
    double p_value_max = calculate_pt_blur_p_value(pt_percent_error, experimental_pdf)  ;
    // record highest p value and the correponding pt percent error 
    //
    if (p_value_max > final_p_value){
      final_p_value = p_value_max;
      final_pt_percent_error = pt_percent_error;
    }
  }

  std::cout<<final_p_value << "    " << final_pt_percent_error<<std::endl;

  // get an intance of the resulting histogram from simulation with optimal pt error calculated above  
  TH1F* simulated_pdf =  get_simulation_pdf(final_pt_percent_error, SAMPLE_SIZE);

  gROOT->ForceStyle();
  experimental_pdf->Draw("PLC");
  simulated_pdf->Draw("PLC SAME");
  gPad->BuildLegend();
  app.Run();

  return 0;
}

TH1F *get_simulation_pdf(double pt_blur_percent,
                         int sample_size = SAMPLE_SIZE) {

  std::vector<TLorentzVector *> parent_vector;
  std::vector<TLorentzVector *> daughter1_vector;
  std::vector<TLorentzVector *> daughter2_vector;
  
  // sample from a uniformly distributed parent particle four momentum
  // with given bounds 
  for (int i = 0; i < SAMPLE_SIZE; i++) {
    TLorentzVector *parent_particle_ptr =
        Random_routines::get_random_lorentz_vector(
            PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, PHI_MIN, PHI_MAX, PHI_MASS);
    parent_vector.push_back(parent_particle_ptr);

    // decay assuming kaon mass
    std::vector<TLorentzVector *> daughter_ptr_pair =
        Random_routines::symmetrical_two_body_decay(parent_particle_ptr,
                                                    KAON_MASS);
    // add fixed percent pt error to daughter particle
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                           pt_blur_percent / 100. *
                                               daughter_ptr_pair[0]->Pt());
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                           pt_blur_percent / 100. *
                                               daughter_ptr_pair[1]->Pt());

    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);

    // manual garbe colleciton
    daughter_ptr_pair.clear();
  }
  // histogram to return
  TH1F *combined_masses =
      new TH1F("Selected Phi Mesons", "Simulated Phi Meson Mass PDF;m_{K^+ K^-}(GeV);probability",
               100, 1., 1.04);

  Selector pid = Selector();

  for (int i = 0; i < SAMPLE_SIZE; i++) {
    TLorentzVector reconstructed_parent_vector =
        *daughter1_vector[i] + *daughter2_vector[i];
    // selection of kaon particle with conditions exactly the same as the experimental selection
    if (daughter1_vector[i]->Pt() > 0.06 && daughter2_vector[i]->Pt() > 0.06 &&
        std::abs(pid.get_NSigmaKaon(daughter1_vector[i])) < 5. &&
        std::abs(pid.get_NSigmaPion(daughter1_vector[i])) > 5. &&
        std::abs(pid.get_NSigmaPion(daughter2_vector[i])) > 5. &&
        std::abs(pid.get_NSigmaKaon(daughter2_vector[i])) < 5.) {
      combined_masses->Fill(reconstructed_parent_vector.M());
    }
  }

  daughter1_vector.clear();
  daughter2_vector.clear();
  parent_vector.clear();
  // convert histogram to pdf
  combined_masses->Scale(1. / combined_masses->Integral());

  return combined_masses;
}

double calculate_pt_blur_p_value(double pt_blur_percent, TH1F *experimental_pdf) {
  double p_value = 0;
  TH1F *rc_pdf = get_simulation_pdf(pt_blur_percent);
  p_value = rc_pdf->Chi2Test(experimental_pdf, "UU;NORM");
  return p_value;
}
