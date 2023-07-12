#include "Random_routines.h"
#include "TLegend.h"
#include "Selection_routines.h"
#include <Rtypes.h>
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
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

// global parameters
double PHI_MASS = 1.019;
double KAON_MASS = 0.493;
double SAMPLE_SIZE = 15000;
double SMEARING_MIN = 1.;
double SMEARING_MAX = 7.;
double LOSS_MIN = 1.;
double LOSS_MAX = 2.5;
double STEP = 0.1;

// method to fill vectors
std::vector<double> arange(double minVal, double maxVal, double step) {
  std::vector<double> result;

  for (double i = minVal; i <= maxVal; i += step) {
    result.push_back(i);
  }

  return result;
}

TH1F *get_simulation_pdf(double pt_blur_percent, double momentum_lost_percent,
                         int sample_size);

double calculate_fitted_pt_blur_chisq(double pt_blur_percent,
                                      double momentum_lost_percent,
                                      TH1 *target_histogram);
// starlight kaon distribution
TFile *kaon_file =
    new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/kaon.root");

// data histogram with constant background subtracted
TFile *data_file = new TFile("/home/xihe/daniel-han-phi-study/phi_KK_toy_model/processed_phi_exp_data.root");

TH1F *data_rc_phi_mass = (TH1F *)data_file->Get("RC Phi Selected");

Selector pid = Selector();

int main(int argc, char** argv) {
  std::vector<double> momentum_smearing_mesh =
      arange(SMEARING_MIN, SMEARING_MAX, STEP);

  std::vector<double> pt_loss_mesh = arange(LOSS_MIN, LOSS_MAX, STEP);

  data_rc_phi_mass->Scale(1. / data_rc_phi_mass->Integral() , "width");

  double final_chisq = DBL_MAX;
  double final_pt_smearing;
  double final_pt_loss;

  for(double pt_smearing:momentum_smearing_mesh){
    for (double pt_loss:pt_loss_mesh){
      double current_chi2 = calculate_fitted_pt_blur_chisq(pt_smearing, pt_loss, data_rc_phi_mass);
      if (current_chi2 < final_chisq){
        final_pt_smearing = pt_smearing;
        final_pt_loss = pt_loss;
        final_chisq = current_chi2;
      }
    }
  }

  std::cout<< "final chisq/ndf = " << final_chisq << "     pt smearing percent = " << final_pt_smearing << "     pt loss = "<< final_pt_loss << std::endl;

  TH1F* sample_pdf = get_simulation_pdf(final_pt_smearing, final_pt_loss, SAMPLE_SIZE);
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 0, 0, 800, 600);
  data_rc_phi_mass->SetLineColor(kBlue);
  sample_pdf->SetLineColor(kRed);
  sample_pdf->Smooth();
  sample_pdf->Draw("PE");
  data_rc_phi_mass->Draw("HIST;SAME");

  // Create a legend
  TLegend *leg = new TLegend(0.1,0.7,0.48,0.9); // You can adjust these parameters to change the position of the legend
  leg->SetHeader("Legend"); // You can set a title for your legend
  leg->AddEntry(data_rc_phi_mass, "Run 19 Au + Au", "l");
  leg->AddEntry(sample_pdf, "MC Starlight with Smearing", "l");
  leg->Draw();

  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  app.Run();

  return 0;
}

TH1F *get_simulation_pdf(double pt_blur_percent, double momentum_lost_percent,
                         int sample_size = SAMPLE_SIZE) {

  std::vector<TLorentzVector *> parent_vector;
  std::vector<TLorentzVector *> daughter1_vector;
  std::vector<TLorentzVector *> daughter2_vector;

  // sample from a uniformly distributed parent particle four momentum
  // with given bounds
  for (int i = 0; i < SAMPLE_SIZE; i++) {
    TLorentzVector *parent_particle_ptr =
        Random_routines::get_slight_lorentz_vector(kaon_file, "kaon");
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
                                         momentum_lost_percent);

    Random_routines::add_pt_percent_loss(daughter_ptr_pair[1],
                                         momentum_lost_percent);

    daughter1_vector.push_back(daughter_ptr_pair[0]);
    daughter2_vector.push_back(daughter_ptr_pair[1]);

    // manual garbe colleciton
    daughter_ptr_pair.clear();
  }
  // histogram to return
  TH1F *combined_masses = new TH1F(
      "Kaon Candidates",
      "m_{K + K} PDF;m_{K^{+}} + m_{K^{-}}(GeV/c^{2});Normalized Amplitude", 100,
      1., 1.04);


  for (int i = 0; i < SAMPLE_SIZE; i++) {
    TLorentzVector reconstructed_parent_vector =
        *daughter1_vector[i] + *daughter2_vector[i];
    // selection of kaon particle with conditions exactly the same as the
    // experimental selection
    if (daughter1_vector[i]->Pt() > 0.06 && daughter2_vector[i]->Pt() > 0.06 &&
        std::abs(pid.compute_NSigmaKaon(daughter1_vector[i]->Pt(), pid.dEdxKaon,
                                        pid.sigma_meson)) < 5. &&
        std::abs(pid.compute_NSigmaPion(daughter1_vector[i]->Pt(), pid.dEdxKaon,
                                        pid.sigma_meson)) > 5. &&
        std::abs(pid.compute_NSigmaPion(daughter2_vector[i]->Pt(), pid.dEdxKaon,
                                        pid.sigma_meson)) > 5. &&
        std::abs(pid.compute_NSigmaKaon(daughter2_vector[i]->Pt(), pid.dEdxKaon,
                                        pid.sigma_meson)) < 5.) {

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


double calculate_fitted_pt_blur_chisq(double pt_blur_percent,
                                      double momentum_lost_percent,
                                      TH1* target_histogram) {
  double chisq;
  TH1F *rc_pdf = get_simulation_pdf(pt_blur_percent, momentum_lost_percent);
  chisq = target_histogram->Chi2Test(rc_pdf, "WW;CHI2/NDF");
  delete rc_pdf;
  return chisq;
}
