// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "FemtoPairFormat.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TF1.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <cmath>
#include <math.h>
#include <string>

// ClassImp( FemtoPair );

int ican = 0;
void makeCan() {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.06);
  can->SetRightMargin(0.01);
}

// funtion definition for fitting
Double_t fit_function(Double_t *x, Double_t *par) {
  // par[0] = background
  // par[1] = amplitude
  // par[2] = mean
  // par[3] = sigma

  double background = par[0];
  double gaussian = par[1] * exp(-0.5 * ((x[0] - par[2]) / (par[3] * x[0])) *
                                 ((x[0] - par[2]) / (par[3] * x[0])));
  double xmin = 1.;
  double xmax = 1.04;
  if (x[0] < xmin || x[0] > xmax) {
    gaussian = 0.0;
  }
  return background + gaussian;
}

void pt_resolution_processor() {

  TFile *fo = new TFile("processed_phi_exp_data.root", "RECREATE");

  TH1F *rc_mass = new TH1F(
      "RC Phi Selected", "Run 19 Au-Au RC Phi Mass;rc_mass (GeV); Count",
      100, 1., 1.1);

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

    //      only select += pair
    if (pair->mChargeSum == 0.) {
      // Setup the Lorentz Vectors of the daugher tracks
      lv1.SetPtEtaPhiM(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi,
                       0.493); // we use Kaon mass here.
      lv2.SetPtEtaPhiM(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.493);
      // compute parent particle lorentz vector from daughters

      // use conservation of momentum and energy to set the four-momentum
      // vector for parent particle.
      lv = lv1 + lv2;
      if (fabs(pair->d1_mNSigmaKaon) < 5 && fabs(pair->d1_mNSigmaPion) > 5 &&
          fabs(pair->d2_mNSigmaKaon) < 5 && fabs(pair->d2_mNSigmaPion) > 5 &&
          pair->d1_mPt > 0.06 && pair->d2_mPt > 0.06) {
        rc_mass->Fill(lv.M());
      }
    } // selection
  }   // loop on events

  TF1 *fit_tf1 = new TF1("Gaussian + Constant Background Fit", fit_function, 1., 1.1, 4);
  fo->cd();
  fit_tf1->SetParameter(0, 4.);
  int max_bin = rc_mass->GetMaximumBin();
  double amplitude = rc_mass->GetBinContent(max_bin) - fit_tf1->GetParameter(0);
  fit_tf1->FixParameter(1, amplitude);
  fit_tf1->SetParameter(2, 1.019);
  fit_tf1->SetParameter(3, 0.01);
  rc_mass->Fit("Gaussian + Constant Background Fit", "0");

  double background = fit_tf1->GetParameter(0);
  double mean = fit_tf1->GetParameter(2);
  double stdev = fit_tf1->GetParameter(3);


  // // Subtract the background from the histogram, accounting for bin widths
  // for (int i = 1; i <= rc_mass->GetNbinsX(); ++i) {
  //   double binContent = rc_mass->GetBinContent(i);
  //   double binWidth = rc_mass->GetBinWidth(i);
  //   rc_mass->SetBinContent(i, binContent - background);
  // }
  // rc_mass->GetXaxis()->SetRangeUser(1., 1.04)  ;

  makeCan();
  rc_mass->Draw("hist");
  fit_tf1->Draw("same");
  gPad->BuildLegend();
  fo->Write();
}
