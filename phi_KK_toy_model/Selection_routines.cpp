#include "Selection_routines.h"
#include "Random_routines.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <TTreeReader.h>

// global rng
thread_local TRandom3 rng_selection(0);

// constructor and variable init
Selector::Selector(double sigma) { this->sigma = sigma; }
// constructor and variable init
Selector::Selector() {}

double Selector::sample_dedx(double p, TH1 *particle_dedx_distr) {
  double sampled_dedx =
      rng_selection.Gaus(particle_dedx_distr->Interpolate(p), sigma);
  return sampled_dedx;
}

double Selector::compute_NSigmaKaon(double p, TH1 *particle_dedx_distr) {
  if (p < 0.11){
    p = 0.11;
  }
  double sampled_dedx = Selector::sample_dedx(p, particle_dedx_distr);
  double kaon_mc = dEdxKaon->Interpolate(p);
  return (sampled_dedx - kaon_mc) / sigma;
}

double Selector::compute_NSigmaPion(double p, TH1 *particle_dedx_distr) {
  if (p < 0.11){
    p = 0.11;
  }
  double sampled_dedx = Selector::sample_dedx(p, particle_dedx_distr);
  double pion_mc = dEdxPion->Interpolate(p);
  return (sampled_dedx - pion_mc) / sigma;
}

double Selector::compute_NSigmaElectron(double p, TH1 *particle_dedx_distr) {
  double sampled_dedx = Selector::sample_dedx(p, particle_dedx_distr);
  double electron_mc = dEdxElectron->Interpolate(p);
  return (sampled_dedx - electron_mc) / sigma;
}

void Selector::draw_dEdx() {
  TH2D *dEdx = new TH2D("dEdx Pion", "dEdX v.s. Pt; Pt(GeV); dEdx", 500, 0.1,
                        10., 500, 0.0025, 0.008);
  // set up momentum mesh
  for (double p = 0.01; p <= 10.; p += 1e-4) {
    double kaon_mc = dEdxKaon->Interpolate(p);
    double pion_mc = dEdxPion->Interpolate(p);
    double electron_mc = dEdxElectron->Interpolate(p);
    dEdx->Fill(p, kaon_mc);
    dEdx->Fill(p, pion_mc);
    dEdx->Fill(p, electron_mc);
  }
  dEdx->Draw();
}

void Selector::draw_dEdx_blurred() {
  TH2D *dEdx_blurred =
      new TH2D("dEdx Kaon + Pion", "dEdX v.s. Pt; Pt(GeV); dEdx", 500, 0.1, 2.,
               500, 0., 0.014);
  // set up momentum mesh
  for (double p = 0.01; p <= 2.; p += 1e-7) {
    double kaon_rc = sample_dedx(p, dEdxKaon);
    double pion_rc = sample_dedx(p, dEdxPion);
    double electron_rc = sample_dedx(p, dEdxElectron);
    dEdx_blurred->Fill(p, kaon_rc);
    dEdx_blurred->Fill(p, pion_rc);
    dEdx_blurred->Fill(p, electron_rc);
  }
  dEdx_blurred->Draw("colz");
}

// int main(int argc, char **argv) {
//   double KAON_MASS = 0.493;
//   double PION_MASS = 0.139;
//   int PHI_SAMPLE_SIZE = 4500;
//   int RHO_SAMPLE_SIZE = PHI_SAMPLE_SIZE * 10;
//   // create vectors for parent vector and daughter vectors
//   std::vector<TLorentzVector *> parent_vector;
//   std::vector<TLorentzVector *> daughter1_vector;
//   std::vector<TLorentzVector *> daughter2_vector;
//
//   double NSigmaPion = 0;
//   double NSigmaKaon = 0;
//
//   TH1F *mc_pion_nsigmapion =
//       new TH1F("mc pion", "NSigmaPion mc Pion", 100, -10, 10);
//   TH1F *mc_pion_nsigmakaon =
//       new TH1F("mc pion", "NSigmaKaon mcPion", 100, -200, 10);
//   TH1F *mc_kaon_nsigmapion =
//       new TH1F(" mc kaon ", "NSigmaPion mc Kaon", 100, -10, 200);
//   TH1F *mc_kaon_nsigmakaon =
//       new TH1F("mc kaon", "NSigmaKaon mc Kaon", 100, -10, 10);
//   TH1F *all_p =
//       new TH1F("all tracks", "Momentum", 100, 0, 2.);
//
//   TH1F *debug =
//       new TH1F("MC KAON NSIGMAPION > 160", "DEBUG KAON MOMENTUM;MOMENTUM(GEV);", 100, 0, 0.5);
//
//   Selector PID = Selector();
//
//   // hardcoded starlight histogram files
//   TFile *kaon_file =
//       new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/kaon.root");
//   TFile *pion_file =
//       new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/pion.root");
//
//   // simulate the decay process for rho -> pi+ pi-
//   for (int i = 0; i < RHO_SAMPLE_SIZE; i++) {
//
//     TLorentzVector *rho_parent_particle_ptr =
//         Random_routines::get_slight_lorentz_vector(pion_file, "pion");
//     parent_vector.push_back(rho_parent_particle_ptr);
//
//     // simulate the decay process with the daughter mass being pion
//     std::vector<TLorentzVector *> daughter_ptr_pair =
//         Random_routines::two_body_decay(rho_parent_particle_ptr, PION_MASS,
//                                         PION_MASS);
//
//     daughter_ptr_pair[0]->SetPtEtaPhiM(daughter_ptr_pair[0]->Pt(),
//                                        daughter_ptr_pair[0]->Eta(),
//                                        daughter_ptr_pair[0]->Phi(), KAON_MASS);
//     daughter_ptr_pair[1]->SetPtEtaPhiM(daughter_ptr_pair[1]->Pt(),
//                                        daughter_ptr_pair[1]->Eta(),
//                                        daughter_ptr_pair[1]->Phi(), KAON_MASS);
//
//     daughter1_vector.push_back(daughter_ptr_pair[0]);
//     daughter2_vector.push_back(daughter_ptr_pair[1]);
//
//     mc_pion_nsigmakaon->Fill(PID.compute_NSigmaKaon(daughter_ptr_pair[0]->P(), PID.dEdxPion));
//     mc_pion_nsigmakaon->Fill(PID.compute_NSigmaKaon(daughter_ptr_pair[1]->P(), PID.dEdxPion));
//
//     mc_pion_nsigmapion->Fill(PID.compute_NSigmaPion(daughter_ptr_pair[0]->P(), PID.dEdxPion));
//     mc_pion_nsigmapion->Fill(PID.compute_NSigmaPion(daughter_ptr_pair[1]->P(), PID.dEdxPion));
//
//     all_p->Fill(daughter_ptr_pair[0]->P());
//     all_p->Fill(daughter_ptr_pair[1]->P());
//   }
//
//   // simualte decay process for phi -> K+ K-
//   for (int i = 0; i < PHI_SAMPLE_SIZE; i++) {
//     // generate parent vector that is a 4 vector uniformly distributed
//     // between
//     // pt, eta, phi bounds with phi mass
//     TLorentzVector *phi_parent_particle_ptr =
//         Random_routines::get_slight_lorentz_vector(kaon_file, "kaon");
//     parent_vector.push_back(phi_parent_particle_ptr);
//
//     // simulate the decay process with the daughter mass being kaon
//     std::vector<TLorentzVector *> daughter_ptr_pair =
//         Random_routines::two_body_decay(phi_parent_particle_ptr, KAON_MASS,
//                                         KAON_MASS);
//
//     daughter1_vector.push_back(daughter_ptr_pair[0]);
//     daughter2_vector.push_back(daughter_ptr_pair[1]);
//
//     mc_kaon_nsigmakaon->Fill(PID.compute_NSigmaKaon(daughter_ptr_pair[0]->P(), PID.dEdxKaon));
//     mc_kaon_nsigmakaon->Fill(PID.compute_NSigmaKaon(daughter_ptr_pair[1]->P(), PID.dEdxKaon));
//
//     mc_kaon_nsigmapion->Fill(PID.compute_NSigmaPion(daughter_ptr_pair[0]->P(), PID.dEdxKaon));
//     mc_kaon_nsigmapion->Fill(PID.compute_NSigmaPion(daughter_ptr_pair[1]->P(), PID.dEdxKaon));
//
//     if (PID.compute_NSigmaPion(daughter_ptr_pair[0]->P(), PID.dEdxKaon) > 160){
//       debug->Fill(daughter_ptr_pair[0]->P());
//     }
//     if (PID.compute_NSigmaPion(daughter_ptr_pair[1]->P(), PID.dEdxKaon) > 160){
//       debug->Fill(daughter_ptr_pair[1]->P());
//     }
//     all_p->Fill(daughter_ptr_pair[0]->P());
//     all_p->Fill(daughter_ptr_pair[1]->P());
//   }
//   TApplication app("app", &argc, argv);
//   TCanvas *canvas = new TCanvas("canvas", "canvas", 0, 0, 800, 600);
//   mc_pion_nsigmakaon->Draw();
//   TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 0, 0, 800, 600);
//   mc_pion_nsigmapion->Draw();
//   TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 0, 0, 800, 600);
//   mc_kaon_nsigmakaon->Draw();
//   TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 0, 0, 800, 600);
//   mc_kaon_nsigmapion->Draw();
//   TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 0, 0, 800, 600);
//   debug->Draw();
//   canvas->Modified();
//   canvas->Update();
//   TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
//   root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
//                        "Terminate()");
//   app.Run();
//   return 0;
// }
