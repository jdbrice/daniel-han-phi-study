// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "FemtoPairFormat.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
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
#include <cmath>
#include <math.h>

int ican = 0;
void makeCan() {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.06);
  can->SetRightMargin(0.01);
}

TRandom3 pair_rng;

double calc_Phi(TLorentzVector lv1, TLorentzVector lv2) {
  TLorentzVector lv_sum = lv1 + lv2;
  TLorentzVector lv_diff = lv1 - lv2;
  double p_dot_m = lv_sum.Px() * lv_diff.Px() + lv_sum.Py() * lv_diff.Py();
  double p_cross_m = lv_sum.Px() * lv_diff.Py() - lv_sum.Py() * lv_diff.Px();
  double cosphi = p_dot_m / (lv_sum.Perp() * lv_diff.Perp());
  double phi = acos(cosphi);
  // flip phi based on +- z direction of \vec{p} \cross \vec{m}
  phi = p_cross_m < 0 ? -phi : phi;
  return phi;
}

double fit_function(double *x, double *par) {
  return par[0] * (1 + par[1] * TMath::Cos(2. * x[0]));
}

void daughter_distribution_kaon_spin() {

  TH1F *phi_distribution =
      new TH1F("Coherent #phi Candidates",
               "#phi distribution; #phi(rad); counts", 12, -3.15, 3.15);

  // Open the file containing the tree (INPUT data).
  TFile *myFile = TFile::Open("input.root");

  // This setup the reader, access the data
  TTreeReader myReader("PairDst", myFile);

  // We need this so ROOT understands the data format
  gSystem->Load("FemtoPairFormat_h.so");
  TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

  // Lorentz vectors (4-vectors) for kinematics in special relativity
  TLorentzVector lv1, lv2, lv;

  // Loop over all entries of the TTree or TChain.
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
      phi_distribution->Fill(phi);
    }
  } // loop on events
  TF1 *f1 = new TF1("f1", fit_function, -TMath::Pi(), TMath::Pi(), 2);
  TFitResultPtr fit_result =
      phi_distribution->Fit("f1", "S"); // Fit and save the result
  TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend->SetHeader("Fit results:", "C");
  legend->AddEntry(
      (TObject *)0,
      Form("p0: %.3f #pm %.3f", fit_result->Parameter(0), fit_result->Error(0)),
      "");
  legend->AddEntry(
      (TObject *)0,
      Form("p1: %.3f #pm %.3f", fit_result->Parameter(1), fit_result->Error(1)),
      "");
  phi_distribution->Draw("pe");
  legend->Draw("same");
}
