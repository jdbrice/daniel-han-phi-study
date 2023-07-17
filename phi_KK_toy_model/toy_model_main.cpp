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
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <cmath>
#include <fstream>
#include <iostream>

thread_local TRandom3 rng_toy(0);

int PHI_SAMPLE_SIZE = 1500;
int PHI_INC_SAMPLE_SIZE =  PHI_SAMPLE_SIZE;
int RHO_SAMPLE_SIZE = 10 * PHI_SAMPLE_SIZE;
int RHO_INC_SAMPLE_SIZE = 10 * PHI_SAMPLE_SIZE;
int ELECTRON_SAMPLE_SIZE = 0;

int main(int argc, char **argv) {
  std::vector<TLorentzVector *> mc_phi, mc_k1, mc_k2, rc_phi, rc_k1, rc_k2;
  std::vector<TLorentzVector *> mc_rho, mc_pi1, mc_pi2, rc_rho, rc_pi1, rc_pi2;
  std::vector<TLorentzVector *> mc_gamma, mc_e1, mc_e2, rc_gamma, rc_e1, rc_e2;
  std::vector<TLorentzVector *> mc_inc_phi, mc_inc_k1, mc_inc_k2, rc_inc_phi,
      rc_inc_k1, rc_inc_k2;
  std::vector<TLorentzVector *> mc_inc_rho, mc_inc_pi1, mc_inc_pi2, rc_inc_rho,
      rc_inc_pi1, rc_inc_pi2;

  std::vector<TLorentzVector *> rc_daughter1, rc_daughter2, rc_parent;

  Random_routines::trigger_two_body_decay("kaon", PHI_SAMPLE_SIZE, mc_phi,
                                          mc_k1, mc_k2, rc_phi, rc_k1, rc_k2);
  Random_routines::trigger_two_body_decay(
      "pion", RHO_SAMPLE_SIZE, mc_rho, mc_pi1, mc_pi2, rc_rho, rc_pi1, rc_pi2);
  Random_routines::trigger_two_body_decay("electron", ELECTRON_SAMPLE_SIZE,
                                          mc_gamma, mc_e1, mc_e2, rc_gamma,
                                          rc_e1, rc_e2);
  Random_routines::trigger_two_body_decay("kaon_inc", PHI_INC_SAMPLE_SIZE,
                                          mc_inc_phi, mc_inc_k1, mc_inc_k2,
                                          rc_inc_phi, rc_inc_k1, rc_inc_k2);
  Random_routines::trigger_two_body_decay("pion_inc", RHO_INC_SAMPLE_SIZE,
                                          mc_inc_rho, mc_inc_pi1, mc_inc_pi2,
                                          rc_inc_rho, rc_inc_pi1, rc_inc_pi2);

  rc_daughter1.insert(rc_daughter1.end(), rc_k1.begin(), rc_k1.end());
  rc_daughter2.insert(rc_daughter2.end(), rc_k2.begin(), rc_k2.end());
  rc_daughter1.insert(rc_daughter1.end(), rc_inc_k1.begin(), rc_inc_k1.end());
  rc_daughter2.insert(rc_daughter2.end(), rc_inc_k2.begin(), rc_inc_k2.end());
  rc_daughter1.insert(rc_daughter1.end(), rc_pi1.begin(), rc_pi1.end());
  rc_daughter2.insert(rc_daughter2.end(), rc_pi2.begin(), rc_pi2.end());
  rc_daughter1.insert(rc_daughter1.end(), rc_inc_pi1.begin(), rc_inc_pi1.end());
  rc_daughter2.insert(rc_daughter2.end(), rc_inc_pi2.begin(), rc_inc_pi2.end());
  rc_daughter1.insert(rc_daughter1.end(), rc_e1.begin(), rc_e1.end());
  rc_daughter2.insert(rc_daughter2.end(), rc_e2.begin(), rc_e2.end());

  rc_parent.insert(rc_parent.end(), rc_phi.begin(), rc_phi.end());
  rc_parent.insert(rc_parent.end(), rc_inc_phi.begin(), rc_inc_phi.end());
  rc_parent.insert(rc_parent.end(), rc_rho.begin(), rc_rho.end());
  rc_parent.insert(rc_parent.end(), rc_inc_rho.begin(), rc_inc_rho.end());
  rc_parent.insert(rc_parent.end(), rc_gamma.begin(), rc_gamma.end());

  // create histogram to draw rc mass
  TH1F *parent_rc_mass_total = new TH1F(
      "Combined Masses", "Toy Model Reco Mass;m_{K^+} + m_{K^-}(GeV);count",
      100, 0.9, 1.1);
  TH1F *parent_rc_mass_mc_phi = new TH1F(
      "Combined Masses", "Toy Model Reco Mass;m_{K^+} + m_{K^-}(GeV);count",
      100, 0.9, 1.1);
  TH1F *parent_rc_mass_mc_rho = new TH1F(
      "Combined Masses", "Toy Model Reco Mass;m_{K^+} + m_{K^-}(GeV);count",
      100, 0.9, 1.1);
  TH1F *parent_rc_mass_mc_electron = new TH1F(
      "Combined Masses", "Toy Model Reco Mass;m_{K^+} + m_{K^-}(GeV);count",
      100, 0.9, 1.1);

  TH1F *nsigmapion =
      new TH1F("All Tracks", "NsigmaPion;NsigmaPion;count", 100, -4, 190);

  TH1F *nsigmakaon =
      new TH1F("All Tracks", "NsigmaKaon;NsigmaKaon;count", 100, -200, 4);

  TH1F *all_pt =
      new TH1F("All Tracks", "Toy Model All Tracks P_{T};P_{T}(GeV);count", 100,
               0., 1.8);

  TH1F *all_eta = new TH1F("All Tracks", "Toy Model All Tracks #eta;#eta;count",
                           100, -4, 4);

  TH1F *all_phi = new TH1F("All Tracks", "Toy Model All Tracks #phi;#phi;count",
                           100, -3.15, 3.15);

  TH1F *kaon_pt =
      new TH1F("kaon Tracks", "Toy Model Kaon Candidate P_{T};P_{T}(GeV);count",
               100, 0., 1.);

  TH1F *kaon_eta = new TH1F(
      "kaon Tracks", "Toy Model Kaon Candidate #eta;#eta;count", 100, -4, 4);

  TH1F *kaon_phi =
      new TH1F("kaon Tracks", "Toy Model Kaon Candidate #phi;#phi;count", 100,
               -3.15, 3.15);

  // create histogram to draw rc mass
  TH2F *parent_rc_mass_pt_total = new TH2F(
      "Combined Masses", "Toy Model Combined Masses;m_{K^+ K^-}(GeV);RC P_T",
      100, 0.98, 1.04, 100, 0, 1.);

  // create an instance of particle selector
  Selector pid = Selector();
  for (int i = 0; i < rc_k1.size(); i++) {
    double pt1 = rc_k1[i]->Pt();
    double pt2 = rc_k2[i]->Pt();
    double p1 = rc_k1[i]->P();
    double p2 = rc_k2[i]->P();
    double nsigmakaon1 =
        pid.compute_NSigmaKaon(p1, pid.dEdxKaon, pid.sigma_meson);
    double nsigmapion1 =
        pid.compute_NSigmaPion(p1, pid.dEdxKaon, pid.sigma_meson);
    double nsigmaelectron1 =
        pid.compute_NSigmaElectron(p1, pid.dEdxKaon, pid.sigma_meson);
    double nsigmakaon2 =
        pid.compute_NSigmaKaon(p2, pid.dEdxKaon, pid.sigma_meson);
    double nsigmapion2 =
        pid.compute_NSigmaPion(p2, pid.dEdxKaon, pid.sigma_meson);
    double nsigmaelectron2 =
        pid.compute_NSigmaElectron(p2, pid.dEdxKaon, pid.sigma_meson);

    if (abs(nsigmapion1) > 5 && abs(nsigmapion2) > 5 && abs(nsigmakaon1) < 5 &&
        abs(nsigmakaon2) < 5 && pt1 > 0.06 && pt2 > 0.06) {
      kaon_pt->Fill(pt1);
      kaon_pt->Fill(pt2);
      parent_rc_mass_total->Fill(rc_phi[i]->M());
      parent_rc_mass_mc_phi->Fill(rc_phi[i]->M());
    }
  }
  for (int i = 0; i < rc_inc_k1.size(); i++) {
    double pt1 = rc_inc_k1[i]->Pt();
    double pt2 = rc_inc_k2[i]->Pt();
    double p1 = rc_inc_k1[i]->P();
    double p2 = rc_inc_k2[i]->P();
    double nsigmakaon1 =
        pid.compute_NSigmaKaon(p1, pid.dEdxKaon, pid.sigma_meson);
    double nsigmapion1 =
        pid.compute_NSigmaPion(p1, pid.dEdxKaon, pid.sigma_meson);
    double nsigmaelectron1 =
        pid.compute_NSigmaElectron(p1, pid.dEdxKaon, pid.sigma_meson);
    double nsigmakaon2 =
        pid.compute_NSigmaKaon(p2, pid.dEdxKaon, pid.sigma_meson);
    double nsigmapion2 =
        pid.compute_NSigmaPion(p2, pid.dEdxKaon, pid.sigma_meson);
    double nsigmaelectron2 =
        pid.compute_NSigmaElectron(p2, pid.dEdxKaon, pid.sigma_meson);

    if (abs(nsigmapion1) > 5 && abs(nsigmapion2) > 5 && abs(nsigmakaon1) < 5 &&
        abs(nsigmakaon2) < 5 && pt1 > 0.06 && pt2 > 0.06) {
      kaon_pt->Fill(pt1);
      kaon_pt->Fill(pt2);
      parent_rc_mass_total->Fill(rc_inc_phi[i]->M());
      parent_rc_mass_mc_phi->Fill(rc_inc_phi[i]->M());
    }
  }
  for (int i = 0; i < rc_pi1.size(); i++) {
    double pt1 = rc_pi1[i]->Pt();
    double pt2 = rc_pi1[i]->Pt();
    double p1 = rc_pi1[i]->P();
    double p2 = rc_pi2[i]->P();
    double nsigmakaon1 =
        pid.compute_NSigmaKaon(p1, pid.dEdxPion, pid.sigma_meson);
    double nsigmapion1 =
        pid.compute_NSigmaPion(p1, pid.dEdxPion, pid.sigma_meson);
    double nsigmaelectron1 =
        pid.compute_NSigmaElectron(p1, pid.dEdxPion, pid.sigma_meson);
    double nsigmakaon2 =
        pid.compute_NSigmaKaon(p2, pid.dEdxPion, pid.sigma_meson);
    double nsigmapion2 =
        pid.compute_NSigmaPion(p2, pid.dEdxPion, pid.sigma_meson);
    double nsigmaelectron2 =
        pid.compute_NSigmaElectron(p2, pid.dEdxPion, pid.sigma_meson);

    if (abs(nsigmapion1) > 5 && abs(nsigmapion2) > 5 && abs(nsigmakaon1) < 5 &&
        abs(nsigmakaon2) < 5 && pt1 > 0.06 && pt2 > 0.06) {
      kaon_pt->Fill(pt1);
      kaon_pt->Fill(pt2);
      parent_rc_mass_total->Fill(rc_rho[i]->M());
      parent_rc_mass_mc_rho->Fill(rc_rho[i]->M());
    }
  }
  for (int i = 0; i < rc_inc_pi1.size(); i++) {
    double pt1 = rc_inc_pi1[i]->Pt();
    double pt2 = rc_inc_pi1[i]->Pt();
    double p1 = rc_inc_pi1[i]->P();
    double p2 = rc_inc_pi2[i]->P();
    double nsigmakaon1 =
        pid.compute_NSigmaKaon(p1, pid.dEdxPion, pid.sigma_meson);
    double nsigmapion1 =
        pid.compute_NSigmaPion(p1, pid.dEdxPion, pid.sigma_meson);
    double nsigmaelectron1 =
        pid.compute_NSigmaElectron(p1, pid.dEdxPion, pid.sigma_meson);
    double nsigmakaon2 =
        pid.compute_NSigmaKaon(p2, pid.dEdxPion, pid.sigma_meson);
    double nsigmapion2 =
        pid.compute_NSigmaPion(p2, pid.dEdxPion, pid.sigma_meson);
    double nsigmaelectron2 =
        pid.compute_NSigmaElectron(p2, pid.dEdxPion, pid.sigma_meson);

    if (abs(nsigmapion1) > 5 && abs(nsigmapion2) > 5 && abs(nsigmakaon1) < 5 &&
        abs(nsigmakaon2) < 5 && pt1 > 0.06 && pt2 > 0.06) {
      kaon_pt->Fill(pt1);
      kaon_pt->Fill(pt2);
      parent_rc_mass_total->Fill(rc_inc_rho[i]->M());
      parent_rc_mass_mc_rho->Fill(rc_rho[i]->M());
    }
  }

  // drawing the result
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 0, 0, 800, 600);
  nsigmakaon->Draw();
  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 0, 0, 800, 600);
  nsigmapion->Draw();
  TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 0, 0, 800, 600);
  all_phi->Draw();
  TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 0, 0, 800, 600);
  // nsigmakaon->Draw();
  kaon_pt->Draw();
  TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 0, 0, 800, 600);
  kaon_eta->Draw();
  TCanvas *canvas6 = new TCanvas("canvas6", "canvas6", 0, 0, 800, 600);
  kaon_phi->Draw();
  parent_rc_mass_total->Draw("pfc");
  parent_rc_mass_mc_phi->Draw("same;pfc");
  parent_rc_mass_mc_rho->Draw("same;pfc");
  parent_rc_mass_mc_electron->Draw("same;pfc");
  auto legend = new TLegend(0.7, 0.55, 0.98, 0.75);
  legend->SetHeader("Reco Parent Mass");
  legend->AddEntry(parent_rc_mass_total, "Total Entries", "f");
  legend->AddEntry(parent_rc_mass_mc_phi, "MC #phi", "f");
  legend->AddEntry(parent_rc_mass_mc_rho, "MC #rho", "f");
  legend->AddEntry(parent_rc_mass_mc_electron, "MC electron pair", "f");
  legend->Draw();
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  app.Run();
  return 0;
}

void add_electron_to_histogram(TH1F *hist, int peak_count, double mass_min,
                               double mass_max) {}
