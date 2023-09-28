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
#include <TProfile.h>
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

double fit_function(double *x, double *par) {
  double fit_value = par[0] * (1 + par[1] * cos(x[0]) + par[2] * cos(2 * x[0]) +
                               par[3] * cos(3 * x[0]) + par[4] * cos(4 * x[0]));
  return fit_value;
}

void daughter_distribution_kaon_spin_pt() {

  TProfile *pt_2phi = new TProfile(
      "prof",
      "Mean 2cos(2#phi) vs p_{T}; p_{T}; <cos(2#phi)> All Kaon Candidates", 10,
      0.0, 1.0);
  TProfile *pt_2phi_mixed =
      new TProfile("prof",
                   "Mean 2cos(2#phi) vs p_{T} Mixed Event; p_{T}; <cos(2#phi)> "
                   "All Kaon Candidates",
                   10, 0.0, 1.0);
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
        pair->mChargeSum == 0 && lv.Pt() < 1.0 && lv.M() >= 1.0 &&
        lv.M() < 1.04) {

      // shuffle the two lv for each event
      double pair_num = pair_rng.Integer(2);
      double phi = calc_Phi(lorentz_vector_list[pair_num],
                            lorentz_vector_list[1 - pair_num]);
      double pair_pt = lv.Pt();
      pt_2phi->Fill(pair_pt, 2. * cos(2. * phi));
      d1_list.push_back(lv1);
      d2_list.push_back(lv2);
    }
  }

  for (int i = 0; i < d1_list.size(); i++) {
    std::vector<TLorentzVector> d1_buffer;

    for (int j = 0; j < d1_list.size(); j++) {
      if (j != i) {
        d1_buffer.push_back(d1_list[j]);
        d1_buffer.push_back(d2_list[j]);
      }
    }

    for (TLorentzVector buffer_lv : d1_buffer) {
      std::vector<TLorentzVector> lorentz_vector_list;
      lorentz_vector_list.push_back(buffer_lv);
      lorentz_vector_list.push_back(d1_list[i]);
      double pair_num = mixed_rng.Integer(2);
      double phi = calc_Phi(lorentz_vector_list[pair_num],
                            lorentz_vector_list[1 - pair_num]);
      if ((lorentz_vector_list[0] + lorentz_vector_list[1]).Pt() < 1.0) {
        pt_2phi_mixed->Fill(
            (lorentz_vector_list[0] + lorentz_vector_list[1]).Pt(),
            2. * cos(2. * phi));
      }
    }
  }

  // pt_two_phi->SetMarkerStyle(20); // e.g., a simple circle.
  // pt_two_phi->SetMarkerSize(2.);
  // pt_two_phi->SetMarkerColor(kBlack);
  TH1F *corrected_2phi = new TH1F(
      "corrected", "Corrected Mean 2cos(2 #phi) vs P_{T}; Pair P_{T}; Mean 2cos(2 #phi)", 10, 0.0, 1.0);
  for (int i = 1; i <= 10; ++i) {
    double binContent1 = pt_2phi->GetBinContent(i);
    double binContent2 = pt_2phi_mixed->GetBinContent(i);
    double binError1 = pt_2phi->GetBinError(i);
    double binError2 = pt_2phi_mixed->GetBinError(i);


    if (pt_2phi->GetBinEntries(i) > 0 && pt_2phi_mixed->GetBinEntries(i) > 0)
    {
      corrected_2phi->SetBinContent(i, binContent1 - binContent2);
      double error_diff = TMath::Sqrt(binError1*binError1 + binError2*binError2);
      corrected_2phi->SetBinError(i, error_diff);
    }
  }

  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // x1, y1, x2, y2
  legend->AddEntry(pt_2phi, "Run 19 Data");
  legend->AddEntry(pt_2phi_mixed, "Mixed Event");
  legend->AddEntry(corrected_2phi, "Data - Mixed Event");

  pt_2phi_mixed->SetLineColor(kRed);
  corrected_2phi->SetLineColor(kGreen);
  makeCan();
  corrected_2phi->SetMinimum(-1.5);
  corrected_2phi->Draw("pe");
  pt_2phi->Draw("same");
  pt_2phi_mixed->Draw("same");
  legend->Draw("same");
}
