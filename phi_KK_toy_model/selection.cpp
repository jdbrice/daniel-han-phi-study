#include "selection.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <TTreeReader.h>

thread_local TRandom3 rng(0);
// constructor realization
Selector::Selector(double sigma) {
  this->sigma = sigma;
  // load from file
  TFile *dEdx_file = new TFile(file_path.c_str());

  // load mpmK and mpmPi from root file to a 1-D histogram
  TH1D *dEdxKaon = (TH1D *)dEdx_file->Get("mpmK");
  TH1D *dEdxPion = (TH1D *)dEdx_file->Get("mpmPi");

  // set up momentum mesh
  for (double p = 0.1; p <= 2.; p += 1e-7) {
    pt.push_back(p);
  }

  // fill dEdx in the momentum mesh for Kaon
  for (double p = 0.1; p <= 2.; p+= 1e-7) {
    dEdx_kaon.push_back(dEdxKaon->Interpolate(p) * 1000);

    dEdx_kaon_blurred.push_back(dEdxKaon->Interpolate(p) * 1000 +
                                rng.Gaus(0, sigma));
  }

  // fill dEdx in the momentum mesh for Pion
  for (double p = 0.1; p<= 2.; p+= 1e-7) {
    dEdx_pion.push_back(dEdxPion->Interpolate(p) * 1000);
    dEdx_pion_blurred.push_back(dEdxPion->Interpolate(p) * 1000 +
                                rng.Gaus(0, sigma));
  }

}

void Selector::draw_dEdx(){
  TH2D *dEdx = new TH2D("dEdx Kaon + Pion", "dEdX v.s. Pt; Pt(GeV); dEdx", 500, 0.1, 2., 500, 0. , 15.);
  for (size_t i = 0; i < pt.size(); i++){
    double p = pt[i];
    dEdx->Fill(p, dEdx_kaon[i]);
    dEdx->Fill(p, dEdx_pion[i]);
  }
  dEdx->Draw("colz");
}

void Selector::draw_dEdx_blurred(){
  TH2D *dEdx_blurred = new TH2D("dEdx Kaon + Pion", "dEdX v.s. Pt; Pt(GeV); dEdx", 500, 0.1, 2., 500, 0. , 15.);
  for (size_t i = 0; i < pt.size(); i++){
    double p = pt[i];
    dEdx_blurred->Fill(p, dEdx_kaon_blurred[i]);
    dEdx_blurred->Fill(p, dEdx_pion_blurred[i]);
  }
  dEdx_blurred->Draw("colz");
}

int main(int argc, char **argv) {
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  std::string file_name = "dEdx.root";
  Selector test_selector = Selector(0.2);
  test_selector.draw_dEdx_blurred();
  app.Run();
  return 0;
}
