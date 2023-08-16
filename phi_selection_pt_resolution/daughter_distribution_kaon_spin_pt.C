#include "FemtoPairFormat.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TH2.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <math.h>

int ican = 0;
void makeCan() {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.06);
  can->SetRightMargin(0.01);
}

TRandom3 pair_rng;
TRandom3 inco_pair_rng;
TRandom3 mixed_rng;

double calc_Phi(TLorentzVector lv1, TLorentzVector lv2) {
  TLorentzVector lv_sum = lv1 + lv2;
  lv1.Boost(-lv_sum.BoostVector());
  lv2.Boost(-lv_sum.BoostVector());
  TLorentzVector lv_diff = lv1 - lv2;
  double p_dot_m = lv_sum.Px() * lv_diff.Px() + lv_sum.Py() * lv_diff.Py();
  double p_cross_m = lv_sum.Px() * lv_diff.Py() - lv_sum.Py() * lv_diff.Px();
  double cosphi = p_dot_m / (lv_sum.Perp() * lv_diff.Perp());
  double phi = acos(cosphi);
  // flip phi based on +- z direction of \vec{p} \cross \vec{m}
  phi = p_cross_m < 0 ? -phi : phi;
  return phi;
}

TH1F *normalize_histogram(TH1 *histogram) {
  TH1F *normalized_hist = (TH1F *)histogram->Clone();
  normalized_hist->GetListOfFunctions()->Delete();
  normalized_hist->Sumw2();
  normalized_hist->Scale(1. / normalized_hist->Integral());
  return normalized_hist;
}

void change_bin_structure(TH1F *target, TH1F *temp) {

  int nBins = temp->GetNbinsX();
  double x_min = temp->GetXaxis()->GetXmin();
  double x_max = temp->GetXaxis()->GetXmax();

  target->SetBins(nBins, x_min, x_max);
}

double fit_function(double *x, double *par) {
  double fit_value = par[0] * (1 + par[1] * cos(x[0]) + par[2] * cos(2 * x[0]) +
                               par[3] * cos(3 * x[0]) + par[4] * cos(4 * x[0]));
  return fit_value;
}

void daughter_distribution_kaon_spin_pt() {

  TH2F *pt_two_phi = new TH2F("Run 19 Au + Au Mixed Event",
                              "Run 19 Au + Au Mixed Event cos2#phi v.s. Pair "
                              "P_{T};P_{T} (GeV/c);cos2#phi",
                              20, 0., 0.5, 10, -20, 20);

  // Open the file containing the tree (INPUT data).
  TFile *myFile = TFile::Open("input.root");

  // This setup the reader, access the data
  TTreeReader myReader("PairDst", myFile);

  // We need this so ROOT understands the data format
  gSystem->Load("FemtoPairFormat_h.so");
  TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

  // Lorentz vectors (4-vectors) for kinematics in special relativity
  TLorentzVector lv1, lv2, lv;

  std::vector<TLorentzVector> d1_list;
  std::vector<TLorentzVector> d2_list;

  TF1 *f1 = new TF1("f1", fit_function, -TMath::Pi(), TMath::Pi(), 5);

  // Loop over all entries of the TTree or TChain.
  for (int bin = 1; bin < pt_two_phi->GetNbinsX(); bin++) {

    TH1F *phi_distribution = new TH1F(Form("phi_dist_bin_%d", bin),
                                      "Phi Distribution", 100, -3.15, 3.15);

    // Reset the reader to loop from the start for each pt bin
    myReader.Restart();

    while (myReader.Next()) {
      lv1.SetPtEtaPhiM(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.493);
      lv2.SetPtEtaPhiM(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.493);
      // reco parent vector
      lv = lv1 + lv2;
      // lorentz vector list of daughters and parent
      std::vector<TLorentzVector> lorentz_vector_list;
      lorentz_vector_list.push_back(lv1);
      lorentz_vector_list.push_back(lv2);

      // kaon candidate filling with PId selection, and P_T selection on phi
      // such that the mesons are likely tobe diffractive
      if ((pair->d1_mNSigmaPion > 5.8 || pair->d1_mNSigmaPion < -4.2) &&
          (pair->d2_mNSigmaPion > 5.8 || pair->d2_mNSigmaPion < -4.2) &&
          (pair->d1_mNSigmaElectron > 5.424112 ||
           pair->d1_mNSigmaElectron < -4.575888) &&
          (pair->d2_mNSigmaElectron > 5.424112 ||
           pair->d2_mNSigmaElectron < -4.575888) &&
          (pair->d1_mNSigmaProton > 7.36 || pair->d1_mNSigmaProton < -2.64) &&
          (pair->d2_mNSigmaProton > 7.36 || pair->d2_mNSigmaProton < -2.64) &&
          (pair->d1_mNSigmaKaon < 7.4 && pair->d1_mNSigmaKaon > -2.6) &&
          (pair->d2_mNSigmaKaon < 7.4 && pair->d2_mNSigmaKaon > -2.6) &&
          pair->mChargeSum == 0 && lv.Pt() < 0.1 && lv.M() >= 1.0 &&
          lv.M() < 1.04) {

        // shuffle the two lv for each event
        double pair_num = pair_rng.Integer(2);
        double phi = calc_Phi(lorentz_vector_list[pair_num],
                              lorentz_vector_list[1 - pair_num]);
        double pair_pt = lv.Pt();

        if (pair_pt >= pt_two_phi->GetXaxis()->GetBinLowEdge(bin) &&
            pair_pt <= pt_two_phi->GetXaxis()->GetBinUpEdge(bin)) {
          phi_distribution->Fill(phi);
        }
      }
    }
    phi_distribution->Fit(f1, "Q");
    double modulation_value = f1->GetParameter(2);

    double pt_bin_center = pt_two_phi->GetXaxis()->GetBinCenter(bin);
    pt_two_phi->Fill(pt_bin_center, modulation_value);

    delete phi_distribution;
  }
  pt_two_phi->SetMarkerStyle(20); // e.g., a simple circle.
  pt_two_phi->SetMarkerSize(2.);
  pt_two_phi->SetMarkerColor(kBlack);
  pt_two_phi->Draw("P");
}
