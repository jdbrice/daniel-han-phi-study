#include "Random_routines.h"
#include "Selection_routines.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

double PHI_MASS = 1.019;
double KAON_MASS = 0.493;
double PT_MIN = 0.1;
double PT_MAX = 10.;
double ETA_MIN = 0.;
double ETA_MAX = 4.;
double PHI_MIN = 0.;
double PHI_MAX = 2 * M_PI;
double SAMPLE_SIZE = 3e5;

// create a instance of simulation with fixed percentage of pt blur and sample
// size
TH1F *get_simulation_pdf(double pt_blur_percent, double energy_lost_percent,
                         int sample_size);

// return the p value from a chi2 calculation between a simulation and
// experiment
double calculate_pt_blur_p_value(double pt_blur_percent,
                                 double energy_lost_percent,
                                 TH1F *experimental_pdf);

double calculate_fitted_pt_blur_chisq(double pt_blur_percent,
                                      double energy_lost_percent,
                                      TF1 *target_function);

// funtion definition for fitting
Double_t fit_function(Double_t *x, Double_t *par) {
  // par[0] = background
  // par[1] = amplitude
  // par[2] = mean
  // par[3] = sigma

  double background = par[0];
  double gaussian = par[1] * exp(-0.5 * ((x[0] - par[2]) / (par[3] * x[0])) *
                                 ((x[0] - par[2]) / (par[3] * x[0])));
  double xmin = 1.;
  double xmax = 1.04;
  if (x[0] < xmin || x[0] > xmax) {
    gaussian = 0.0;
  }
  return background + gaussian;
}

int main(int argc, char **argv) {
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");

  double fitted_background = 0.;
  double fitted_mean = 1.01794e+00;
  double fitted_stdev = 3.31181e-03;
  double normalized_amplitude = 1. / (std::sqrt(2. * M_PI) * fitted_stdev);

  TF1 *fit_tf1 = new TF1("Experimental Fit", fit_function, 1., 1.04, 4);
  fit_tf1->SetParameter(0, fitted_background);
  fit_tf1->SetParameter(1, normalized_amplitude);
  fit_tf1->SetParameter(2, fitted_mean);
  fit_tf1->SetParameter(3, fitted_stdev);

  double final_chisq = 1e9;

  double final_pt_percent_error = 0.;
  double final_pt_lost_percent = 0.;
  // run simulation with pt percent error from 0 percent 1o 10 percent with step
  // size 0.1
  // this upper limit is testd experimentally
  for (double pt_percent_error = 1.; pt_percent_error <= 6.;
       pt_percent_error += 0.1) {
    for (double pt_lost_percent = 0.8; pt_lost_percent <= 2.;
         pt_lost_percent += 0.1) {

      double chi_sq_min = calculate_fitted_pt_blur_chisq(
          pt_percent_error, pt_lost_percent, fit_tf1);
      // record highest p value and the correponding pt percent error
      //
      if (chi_sq_min < final_chisq) {
        final_chisq = chi_sq_min;
        final_pt_percent_error = pt_percent_error;
        final_pt_lost_percent = pt_lost_percent;
      }
    }
  }
  std::cout << final_chisq << "    " << final_pt_percent_error << "    " << final_pt_lost_percent << std::endl;

  TH1F *sample_hist = get_simulation_pdf(final_pt_percent_error,
                                         final_pt_lost_percent, SAMPLE_SIZE);
  sample_hist->Draw("hist");
  fit_tf1->Draw("same");

  std::cout << sample_hist->GetMean() << std::endl;
  gPad->BuildLegend();
  app.Run();

  return 0;
}

TH1F *get_simulation_pdf(double pt_blur_percent, double energy_lost_percent,
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
        Random_routines::two_body_decay(parent_particle_ptr, KAON_MASS,
                                        KAON_MASS);
    // add fixed percent pt error to daughter particle
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                           pt_blur_percent / 100. *
                                               daughter_ptr_pair[0]->Pt());
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                           pt_blur_percent / 100. *
                                               daughter_ptr_pair[1]->Pt());

    Random_routines::add_pt_percent_loss(daughter_ptr_pair[0],
                                         energy_lost_percent);

    Random_routines::add_pt_percent_loss(daughter_ptr_pair[1],
                                         energy_lost_percent);

    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);

    // manual garbe colleciton
    daughter_ptr_pair.clear();
  }
  // histogram to return
  TH1F *combined_masses = new TH1F(
      "Selected Phi Mesons",
      "Simulated Phi Meson Mass PDF;m_{K^+ K^-}(GeV);Normalized Amplitude", 100,
      1., 1.04);

  Selector pid = Selector();

  for (int i = 0; i < SAMPLE_SIZE; i++) {
    TLorentzVector reconstructed_parent_vector =
        *daughter1_vector[i] + *daughter2_vector[i];
    // selection of kaon particle with conditions exactly the same as the
    // experimental selection
    if (daughter1_vector[i]->Pt() > 0.06 && daughter2_vector[i]->Pt() > 0.06 &&
        std::abs(pid.get_NSigmaKaon(daughter1_vector[i])) < 50. &&
        std::abs(pid.get_NSigmaPion(daughter1_vector[i])) > 50. &&
        std::abs(pid.get_NSigmaPion(daughter2_vector[i])) > 50. &&
        std::abs(pid.get_NSigmaKaon(daughter2_vector[i])) < 50.) {
      combined_masses->Fill(reconstructed_parent_vector.M());
    }
    delete parent_vector[i];
    delete daughter1_vector[i];
    delete daughter2_vector[i];
  }

  daughter1_vector.clear();
  daughter2_vector.clear();
  parent_vector.clear();
  // convert histogram to pdf
  combined_masses->Scale(1. / combined_masses->Integral(), "width");

  return combined_masses;
}

double calculate_pt_blur_p_value(double pt_blur_percent,
                                 double energy_lost_percent,
                                 TH1F *experimental_pdf) {
  double p_value = 0;
  TH1F *rc_pdf = get_simulation_pdf(pt_blur_percent, energy_lost_percent);
  p_value = rc_pdf->Chi2Test(experimental_pdf, "UU;NORM");
  delete rc_pdf;
  return p_value;
}

double calculate_fitted_pt_blur_chisq(double pt_blur_percent,
                                      double energy_lost_percent,
                                      TF1 *target_function) {
  double chisq = 0;
  TH1F *rc_pdf = get_simulation_pdf(pt_blur_percent, energy_lost_percent);
  // rc_pdf->Sumw2();
  chisq = rc_pdf->Chisquare(target_function, "R, 1., 1.04");
  delete rc_pdf;
  return chisq;
}
