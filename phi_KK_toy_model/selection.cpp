#include "selection.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <TTreeReader.h>

// global rng
thread_local TRandom3 rng(0);

// return the index of an element in an vector<double> such that the element of that index
// in the array is the closest (smallest 1-d eucledian norm) to the give double
int find_closest_index(const std::vector<double>* vec, double target){
  double min_diff = std::abs((*vec)[0] -target);
  int best_index = 0;

  for (int i = 0; i < vec->size(); i++)
  {
    double diff = std::abs((*vec)[i] - target);
    if (diff < min_diff){
      min_diff = diff;
      best_index = i;
    }
  }
  return best_index;
}

// constructor and variable init
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
  for (double p : pt) {
    dEdx_kaon.push_back(dEdxKaon->Interpolate(p) * 1000);

    dEdx_kaon_blurred.push_back(dEdxKaon->Interpolate(p) * 1000 +
                                rng.Gaus(0, sigma));
  }

  // fill dEdx in the momentum mesh for Pion
  for (double p : pt) {
    dEdx_pion.push_back(dEdxPion->Interpolate(p) * 1000);
    dEdx_pion_blurred.push_back(dEdxPion->Interpolate(p) * 1000 +
                                rng.Gaus(0, sigma));
  }
}
// constructor and variable init
Selector::Selector() {
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
  for (double p : pt) {
    dEdx_kaon.push_back(dEdxKaon->Interpolate(p) * 1000);

    dEdx_kaon_blurred.push_back(dEdxKaon->Interpolate(p) * 1000 +
                                rng.Gaus(0, sigma));
  }

  // fill dEdx in the momentum mesh for Pion
  for (double p : pt) {
    dEdx_pion.push_back(dEdxPion->Interpolate(p) * 1000);
    dEdx_pion_blurred.push_back(dEdxPion->Interpolate(p) * 1000 +
                                rng.Gaus(0, sigma));
  }
}

void Selector::draw_dEdx() {
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

double Selector::get_NSigmaKaon(TLorentzVector *mc_ptr){
  double dEdx_mc;
  double dEdx_rc;

  double pt_mc = mc_ptr->Pt();

  // detect if input four vector has greater than allowed pt 
  if (pt_mc > 2.){
    std:: cout << "Error! Pt greater than 2GeV" << std::endl;
    return 0;
  }
  int pt_mc_index = find_closest_index(&pt, pt_mc);

  // case where we assume for true Kaon
  if (std::abs(mc_ptr->M() - 0.493) < 0.001){
    dEdx_mc = dEdx_kaon[pt_mc_index]; 
    dEdx_rc = dEdx_kaon_blurred[pt_mc_index]; 
    return (dEdx_mc - dEdx_rc)/sigma; 
  }

  // case where we assume for true Pion
  else if(std::abs(mc_ptr->M() - 0.139) < 0.001){
    dEdx_mc = dEdx_kaon[pt_mc_index];
    dEdx_rc = dEdx_pion_blurred[pt_mc_index];
    return (dEdx_mc - dEdx_rc)/sigma; 
  }

  // case wher the mass of the particle is not Kaon nor Pion
  else{
    std::cout << "Error, four vector provided is neither Kaon nor Pion"<< std::endl;
    return 0;
  }
}

double Selector::get_NSigmaPion(TLorentzVector *mc_ptr){
  double dEdx_mc;
  double dEdx_rc;

  double pt_mc = mc_ptr->Pt();

  // detect if input four vector has greater than allowed pt 
  if (pt_mc > 2.){
    std:: cout << "Error! Pt greater than 2GeV" << std::endl;
    return 0;
  }
  int pt_mc_index = find_closest_index(&pt, pt_mc);

  // case where we assume for true Kaon
  if (std::abs(mc_ptr->M() - 0.493) < 0.001){
    dEdx_mc = dEdx_pion[pt_mc_index]; 
    dEdx_rc = dEdx_kaon_blurred[pt_mc_index]; 
    return (dEdx_mc - dEdx_rc)/sigma; 
  }

  // case where we assume for true Pion
  else if(std::abs(mc_ptr->M() - 0.139) < 0.001){
    dEdx_mc = dEdx_pion[pt_mc_index];
    dEdx_rc = dEdx_pion_blurred[pt_mc_index];
    return (dEdx_mc - dEdx_rc)/sigma; 
  }

  // case wher the mass of the particle is not Kaon nor Pion
  else{
    std::cout << "Error, four vector provided is neither Kaon nor Pion"<< std::endl;
    return 0;
  }
}


int main(int argc, char **argv) {
  TApplication app("app", &argc, argv);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  std::string file_name = "dEdx.root";
  Selector test_selector = Selector();
  test_selector.draw_dEdx_blurred();
  app.Run();
  return 0;
}
