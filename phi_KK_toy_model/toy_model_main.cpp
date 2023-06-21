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

double KAON_MASS = 0.493;
double PION_MASS = 0.139;
double MUON_MASS = 0.105658;
double BRANCHING_RATIO = 0.6356;
double NEUTRINO_MASS = MUON_MASS / 1000.;
int PHI_SAMPLE_SIZE = 45000;
int RHO_SAMPLE_SIZE = 10 * PHI_SAMPLE_SIZE;

int main(int argc, char **argv) {
  // create vectors for parent vector and daughter vectors
  std::vector<TLorentzVector *> phi_ptr_vector;
  std::vector<TLorentzVector *> kp_ptr_vector;
  std::vector<TLorentzVector *> km_ptr_vector;

  // hardcoded starlight histogram files
  TFile *kaon_file =
      new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/kaon.root");
  TFile *pion_file =
      new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/pion.root");

  // simulate the decay process for rho -> pi+ pi-
  for (int i = 0; i < RHO_SAMPLE_SIZE; i++) {

    TLorentzVector *rho_parent_particle_ptr =
        Random_routines::get_slight_lorentz_vector(pion_file, "pion");
    phi_ptr_vector.push_back(rho_parent_particle_ptr);

    // simulate the decay process with the daughter mass being pion
    std::vector<TLorentzVector *> daughter_ptr_pair =
        Random_routines::two_body_decay(rho_parent_particle_ptr, PION_MASS,
                                        PION_MASS);

    // blur the daughter particles by 2 percent to simulate actual particle
    // detector accuracy
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                           0.04 * daughter_ptr_pair[0]->Pt());
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                           0.04 * daughter_ptr_pair[1]->Pt());

    kp_ptr_vector.push_back(daughter_ptr_pair[0]);
    km_ptr_vector.push_back(daughter_ptr_pair[1]);

    // now change the daughter mass to be Kaon- the assumed particle when doing
    // Kaon selections
    kp_ptr_vector[i]->SetPtEtaPhiM(kp_ptr_vector[i]->Pt(),
                                   kp_ptr_vector[i]->Eta(),
                                   kp_ptr_vector[i]->Phi(), KAON_MASS);
    km_ptr_vector[i]->SetPtEtaPhiM(km_ptr_vector[i]->Pt(),
                                   km_ptr_vector[i]->Eta(),
                                   km_ptr_vector[i]->Phi(), KAON_MASS);
  }

  // simualte decay process for phi -> K+ K-
  for (int i = 0; i < PHI_SAMPLE_SIZE; i++) {
    // generate parent vector that is a 4 vector uniformly distributed between
    // pt, eta, phi bounds with phi mass
    TLorentzVector *phi_parent_particle_ptr =
        Random_routines::get_slight_lorentz_vector(kaon_file, "kaon");
    phi_ptr_vector.push_back(phi_parent_particle_ptr);

    // simulate the decay process with the daughter mass being kaon
    std::vector<TLorentzVector *> daughter_ptr_pair =
        Random_routines::two_body_decay(phi_parent_particle_ptr, KAON_MASS,
                                        KAON_MASS);

    // blur the daughter particles by 2 percent to simulate actual particle
    // detector accuracy
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                           0.04 * daughter_ptr_pair[0]->Pt());
    Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                           0.04 * daughter_ptr_pair[1]->Pt());
    kp_ptr_vector.push_back(daughter_ptr_pair[0]);
    km_ptr_vector.push_back(daughter_ptr_pair[1]);
  }

  // create histogram to draw rc mass
  TH1F *parent_rc_mass_total = new TH1F(
      "Combined Masses", "Toy Model Combined Masses;m_{K^+ K^-}(GeV);count",
      100, 1., 1.1);

  // create histogram to draw rc mass
  TH1F *parent_rc_mass_tof = new TH1F(
      "Combined Masses Triggered by TOF",
      "Toy Model Combined Masses;m_{K^+ K^-}(GeV);count", 100, 1., 1.1);

  // create histogram to draw rc mass
  TH2F *parent_rc_mass_pt_total = new TH2F(
      "Combined Masses", "Toy Model Combined Masses;m_{K^+ K^-}(GeV);RC P_T",
      100, 1., 1.1, 100, 0, 0.5);

  // create histogram to draw rc mass
  TH2F *parent_rc_mass_pt_tof =
      new TH2F("Combined Masses Triggered by TOF",
               "Toy Model Combined Masses;m_{K^+ K^-}(GeV);RC P_T", 100, 1.,
               1.1, 100, 0, 0.5);

  // create an instance of particle selector
  Selector pid = Selector();

<<<<<<< HEAD
  // select daughter particle to be Kaon
  for (int i = 0; i < parent_vector.size(); i++)
  {
    if (std::abs(pid.get_NSigmaKaon(daughter1_vector[i])) < 5. &&
        std::abs(pid.get_NSigmaKaon(daughter2_vector[i])) < 5. &&
        std::abs(pid.get_NSigmaPion(daughter1_vector[i])) > 50. &&
        std::abs(pid.get_NSigmaPion(daughter2_vector[i])) > 50. &&
        daughter1_vector[i]->Pt() > 0.06 && daughter2_vector[i]->Pt() > 0.06)
    {
=======
  // select daughter particle to be Kaon. This is for case where daughter
  // particle has momentum > 60 MeV
  for (int i = 0; i < phi_ptr_vector.size(); i++) {
    // rng for branching ratio of kaon -> muon decay
    double muon_rng = rng_toy.Uniform();
    if (std::abs(pid.get_NSigmaKaon(kp_ptr_vector[i])) < 50. &&
        std::abs(pid.get_NSigmaKaon(km_ptr_vector[i])) < 50. &&
        std::abs(pid.get_NSigmaPion(kp_ptr_vector[i])) > 50. &&
        std::abs(pid.get_NSigmaPion(km_ptr_vector[i])) > 50. &&
        kp_ptr_vector[i]->Pt() > 0.06 && km_ptr_vector[i]->Pt() > 0.06) {
>>>>>>> muon-decay
      // reconstruct the dauther particles only if they are kaons.
      // This effectively selects phi
      TLorentzVector reconstructed_parent =
          *(kp_ptr_vector[i]) + *(km_ptr_vector[i]);
      // fill the reconstructed parent mass
      parent_rc_mass_total->Fill(reconstructed_parent.M());
      parent_rc_mass_pt_total->Fill(reconstructed_parent.M(),
                                    reconstructed_parent.Pt());
    }

    // case where both kaon momentum < 60 MeV, both kaons decayed into muons
    else if (std::abs(pid.get_NSigmaKaon(kp_ptr_vector[i])) < 50. &&
             std::abs(pid.get_NSigmaKaon(km_ptr_vector[i])) < 50. &&
             std::abs(pid.get_NSigmaPion(kp_ptr_vector[i])) > 50. &&
             std::abs(pid.get_NSigmaPion(km_ptr_vector[i])) > 50. &&
             kp_ptr_vector[i]->Pt() < 0.06 && km_ptr_vector[i]->Pt() < 0.06 &&
             muon_rng < BRANCHING_RATIO * BRANCHING_RATIO) {

      // case for positively charged muon
      std::vector<TLorentzVector *> mup_ptr_vector =
          Random_routines::two_body_decay(kp_ptr_vector[i], MUON_MASS,
                                          NEUTRINO_MASS);
      Random_routines::add_gaussian_pt_error(mup_ptr_vector[0],
                                             0.04 * mup_ptr_vector[0]->Pt());
      Random_routines::add_gaussian_pt_error(mup_ptr_vector[1],
                                             0.04 * mup_ptr_vector[1]->Pt());

      TLorentzVector rc_kp = *mup_ptr_vector[0] + *mup_ptr_vector[1];

      // case for negatively charged muon
      std::vector<TLorentzVector *> mum_ptr_vector =
          Random_routines::two_body_decay(km_ptr_vector[i], MUON_MASS,
                                          NEUTRINO_MASS);
      Random_routines::add_gaussian_pt_error(mum_ptr_vector[0],
                                             0.04 * mum_ptr_vector[0]->Pt());
      Random_routines::add_gaussian_pt_error(mum_ptr_vector[1],
                                             0.04 * mum_ptr_vector[1]->Pt());
      TLorentzVector rc_km = *mum_ptr_vector[0] + *mum_ptr_vector[1];

      // reconstructed two rc kaon
      TLorentzVector rc_phi = rc_km + rc_kp;

      // case where at least one track is detected by tof
      if (mup_ptr_vector[0]->Pt() > 0.2 || mum_ptr_vector[0]->Pt() > 0.2) {
        parent_rc_mass_total->Fill(rc_phi.M());
      parent_rc_mass_pt_total->Fill(rc_phi.M(),
                                    rc_phi.Pt());

        parent_rc_mass_tof->Fill(rc_phi.M());
      parent_rc_mass_pt_tof->Fill(rc_phi.M(),
                                    rc_phi.Pt());
      }
    }
  }

  std::cout << parent_rc_mass_tof->GetEntries() /
                   parent_rc_mass_total->GetEntries() * 100.
            << std::endl;
  // drawing the result
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  parent_rc_mass_pt_total->Draw("colz");
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  app.Run();
  return 0;
}
