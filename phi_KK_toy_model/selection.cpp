#include "selection.h"
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

// return the index of an element in an vector<double> such that the element of
// that index in the array is the closest (smallest 1-d eucledian norm) to the
// give double
int find_closest_index(const std::vector<double> *vec, double target) {
  double min_diff = std::abs((*vec)[0] - target);
  int best_index = 0;

  for (int i = 0; i < vec->size(); i++) {
    double diff = std::abs((*vec)[i] - target);
    if (diff < min_diff) {
      min_diff = diff;
      best_index = i;
    }
  }
  return best_index;
}

// constructor and variable init
Selector::Selector(double sigma) { this->sigma = sigma; }
// constructor and variable init
Selector::Selector() {}

void Selector::draw_dEdx() {
  std::vector<double> dEdx_kaon;
  std::vector<double> dEdx_pion;
  std::vector<double> pt;
  // set up momentum mesh
  for (double p = 0.1; p <= 2.; p += 1e-7) {
    pt.push_back(p);
  }
  // fill dEdx in the momentum mesh for Kaon
  for (double p : pt) {

    dEdx_kaon.push_back(dEdxKaon->Interpolate(p) * 1000);
  }

  // fill dEdx in the momentum mesh for Pion
  for (double p : pt) {
    dEdx_pion.push_back(dEdxPion->Interpolate(p) * 1000);
  }
  TH2D *dEdx = new TH2D("dEdx Kaon + Pion", "dEdX v.s. Pt; Pt(GeV); dEdx", 500,
                        0.1, 2., 500, 0., 15.);
  for (size_t i = 0; i < pt.size(); i++) {
    double p = pt[i];
    dEdx->Fill(p, dEdx_kaon[i]);
    dEdx->Fill(p, dEdx_pion[i]);
  }
  dEdx->Draw("colz");
}

void Selector::draw_dEdx_blurred() {
  std::vector<double> dEdx_kaon_blurred;
  std::vector<double> dEdx_pion_blurred;
  std::vector<double> pt;

  // set up momentum mesh
  for (double p = 0.1; p <= 2.; p += 1e-7) {
    pt.push_back(p);
  }
  // fill dEdx in the momentum mesh for Kaon
  for (double p : pt) {

    dEdx_kaon_blurred.push_back(dEdxKaon->Interpolate(p) * 1000 +
                                rng_selection.Gaus(0, sigma));
  }

  // fill dEdx in the momentum mesh for Pion
  for (double p : pt) {
    dEdx_pion_blurred.push_back(dEdxPion->Interpolate(p) * 1000 +
                                rng_selection.Gaus(0, sigma));
  }
  TH2D *dEdx_blurred =
      new TH2D("dEdx Kaon + Pion", "dEdX v.s. Pt; Pt(GeV); dEdx", 500, 0.1, 2.,
               500, 0., 15.);
  for (size_t i = 0; i < pt.size(); i++) {
    double p = pt[i];
    dEdx_blurred->Fill(p, dEdx_kaon_blurred[i]);
    dEdx_blurred->Fill(p, dEdx_pion_blurred[i]);
  }
  dEdx_blurred->Draw("colz");
}

double Selector::get_NSigmaKaon(TLorentzVector *mc_ptr) {
  double NSigmaKaon;
  double mc_pt = mc_ptr->Pt();
  double dEdxKaon_mc = dEdxKaon->Interpolate(mc_pt) * 1000;
  double dEdxPion_mc = dEdxPion->Interpolate(mc_pt) * 1000;
  double dEdxKaon_rc =
      dEdxKaon->Interpolate(mc_pt) * 1000 + rng_selection.Gaus(0, sigma);

  if (std::abs(mc_ptr->M() - 0.493) <= 0.001) {
    NSigmaKaon = (dEdxKaon_mc - dEdxKaon_rc) / sigma;
    return NSigmaKaon;
  } else if (std::abs(mc_ptr->M() - 0.139) <= 0.001) {
    NSigmaKaon = (dEdxPion_mc - dEdxKaon_rc) / sigma;
    return NSigmaKaon;
  } else {
    std::cout << "Error, particle given has neither Kaon mass nor Pion mass"
              << std::endl;
    return 0;
  }
}

double Selector::get_NSigmaPion(TLorentzVector *mc_ptr) {
  double NSigmaPion;
  double mc_pt = mc_ptr->Pt();
  double dEdxKaon_mc = dEdxKaon->Interpolate(mc_pt) * 1000;
  double dEdxPion_mc = dEdxPion->Interpolate(mc_pt) * 1000;
  double dEdxPion_rc =
      dEdxPion->Interpolate(mc_pt) * 1000 + rng_selection.Gaus(0, sigma);

  if (std::abs(mc_ptr->M() - 0.493) <= 0.001) {
    NSigmaPion = (dEdxKaon_mc - dEdxPion_rc) / sigma;
    return NSigmaPion;
  } else if (std::abs(mc_ptr->M() - 0.139) <= 0.001) {
    NSigmaPion = (dEdxPion_mc - dEdxPion_rc) / sigma;
    return NSigmaPion;
  } else {
    std::cout << "Error, particle given has neither Kaon mass nor Pion mass"
              << std::endl;
    return 0;
  }
}

int main(int argc, char **argv) {
  std::vector<TLorentzVector *> kaons;
  Selector pid = Selector();
  // TH1D *kaon_acceptance = new TH1D(
  //     "Acceptance rate", "Accepted Kaons from Kaon mc; Pt(Gev);Acceptance Rate", 500, 0, 1);
  //
  // TH1D kaon_acceptance_normalized = TH1D(
  //     "Acceptance rate", "Accepted Kaons from Kaon mc; Pt(Gev);Acceptance Rate", 500, 0, 1);
  //
  // int mc_produced = 1e5;
  // for (double p = 0.; p < 1.; p += 1e-3) {
  //   for (int i = 0; i < mc_produced; i++) {
  //     TLorentzVector *random_mc = Random_routines::get_random_lorentz_vector(
  //         p, p, 0., 4., 0., 2. * M_PI, 0.493);
  //     Random_routines::add_gaussian_pt_error(random_mc, 0.01 * random_mc->Pt());
  //     if (pid.get_NSigmaKaon(random_mc) < 5 &&
  //         pid.get_NSigmaPion(random_mc) > 5) {
  //       kaon_acceptance->Fill(p);
  //     }
  //     delete random_mc;
  //   }
  // }
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  // kaon_acceptance_normalized = *kaon_acceptance * (1. / 2e5);
  // kaon_acceptance_normalized.Draw();
  pid.draw_dEdx_blurred();
  canvas->Modified();
  canvas->Update();
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  app.Run();
  return 0;
}
