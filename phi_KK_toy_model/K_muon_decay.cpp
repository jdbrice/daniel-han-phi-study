#include "Random_routines.h"
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
#include <vector>

thread_local TRandom3 rng_toy(0);

double KAON_MASS = 0.493;
double MUON_MASS = 0.105;
double MUON_NUTRINO_MASS = KAON_MASS / 1000.;
double PT_MIN = 0.0;
double PT_MAX = 10.;
double ETA_MIN = -1.;
double ETA_MAX = 1.;
double PHI_MIN = 0.;
double PHI_MAX = 2 * M_PI;
int K_SAMPLE_SIZE = 4500000;

TRandom3 decay_rng;

int main(int argc, char **argv)
{
  // generate parent Kaon particles
  std::vector<TLorentzVector *> kaon_list;
  std::vector<TLorentzVector *> muon_list;
  std::vector<TLorentzVector *> nutrino_list;

  // create histogram to draw rc mass
  TH1F *TOF_trigger_rate = new TH1F(
      "Muon P_T > 200 MeV", "Muon TOF Trigger Rate; Kaon P_T (GeV);Rate", 100,
      PT_MIN, PT_MAX);

  TH1F *kaon_pt_hist =
      new TH1F("MC Kaon", "Kaon PT distibution; Kaon P_T (GeV);count", 100,
               PT_MIN, PT_MAX);

  TH1F *delta_phi =
      new TH1F("MC Kaon", "Kaon Muon DeltaPhi Distibution; DeltaPhi(rad);count",
               100, -M_PI, M_PI);

  TH1F *delta_phi_tof = new TH1F(
      "Muon P_T > 200 MeV",
      "Kaon Muon DeltaPhi Distibution; DeltaPhi(rad);count", 100, -M_PI, M_PI);

  TH2F *delta_phi_muon_pt =
      new TH2F("K-> mu + v_{mu}",
               "Kaon Muon DeltaPhi Distibution; Muon_P_T(GeV); DeltaPhi(Rad)",
               100, 0., 0.35, 100, -M_PI, M_PI);

  TH2F *delta_phi_muon_pt_tof =
      new TH2F("Muon > 200 MeV",
               "Kaon Muon DeltaPhi Distibution; Muon P_T (GeV); DeltaPhi(Rad)",
               100, 0., 0.35, 100, -M_PI, M_PI);

  // list of momentum value
  std::vector<double> branched_kaon_pt_list;
  std::vector<double> branched_kaon_energy_list;

  std::vector<double> test_energy_difference;

  // create mc kaons
  for (int i = 0; i < K_SAMPLE_SIZE; i++)
  {
    kaon_list.push_back(Random_routines::get_random_lorentz_vector(
        PT_MIN, PT_MAX, ETA_MIN, ETA_MAX, PHI_MIN, PHI_MAX, KAON_MASS));
  }

  // simulate the decay process of K->mu + v_mu
  for (TLorentzVector *kaon_ptr : kaon_list)
  {
    double decay_number = decay_rng.Uniform(); // branching ratio number

    kaon_pt_hist->Fill(
        kaon_ptr
            ->Pt()); // this histogram contains the kaon P_T for all mc Kaons.

    // simulate the branching probability
    if (decay_number < 0.6356)
    {

      std::vector<TLorentzVector *> daughter_list =
          Random_routines::two_body_decay(kaon_ptr, MUON_MASS,
                                          MUON_NUTRINO_MASS);

      muon_list.push_back(daughter_list[0]);

      nutrino_list.push_back(daughter_list[1]);

      delta_phi->Fill(kaon_ptr->DeltaPhi(*daughter_list[0]));

      delta_phi_muon_pt->Fill(daughter_list[0]->Pt(),
                              kaon_ptr->DeltaPhi(*daughter_list[0]));

      branched_kaon_pt_list.push_back(kaon_ptr->Pt()); // this list is only for kaons that branch into muon
      branched_kaon_energy_list.push_back(kaon_ptr->E());

      if (daughter_list[0]->Pt() > 0.2)
      {

        delta_phi_tof->Fill(kaon_ptr->DeltaPhi(*daughter_list[0]));

        delta_phi_muon_pt_tof->Fill(daughter_list[0]->Pt(),
                                    kaon_ptr->DeltaPhi(*daughter_list[0]));
      }
    }
  }

  // select from the branched kaon pt such that the decayed muon has pt > 0.2
  // Note that there is a index corresponance between muon_list and branched_kaon_pt_list
  for (int i = 0; i < muon_list.size(); i++)
  {
    if (muon_list[i]->Pt() > 0.2)
    {
      TOF_trigger_rate->Fill(branched_kaon_pt_list[i]);
      test_energy_difference.push_back(branched_kaon_energy_list[i] - muon_list[i]->E() - nutrino_list[i]->E());
    }
  }

  TOF_trigger_rate->Divide(kaon_pt_hist);
  // drawing the result
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 0, 0, 800, 600);
  TOF_trigger_rate->Draw();
  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 0, 0, 800, 600);
  delta_phi->Draw();
  TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 0, 0, 800, 600);
  delta_phi_tof->Draw();
  TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 0, 0, 800, 600);
  delta_phi_muon_pt->Draw("colz");
  TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 0, 0, 800, 600);
  delta_phi_muon_pt_tof->Draw("colz");
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  app.Run();
  return 0;
}
