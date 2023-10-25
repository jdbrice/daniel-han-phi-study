
// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "FemtoPairFormat.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TH2.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <cmath>
#include <complex>
#include <math.h>

const double NSigmaPionShift = 0.8;
const double NSigmaElectronShift = 0.424112;
const double NSigmaProtonShift = 2.36;
const double NSigmaKaonShift = 2.4;

// ClassImp( FemtoPair );

int ican = 0;
void makeCan(bool logz = false) {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.10);
  can->SetRightMargin(0.03);

  if (logz) {
    can->SetLogz(1); // Set Z-axis to logarithmic scale
  }
}

void KaonHist() {

  // histogram for reco phi
  TH1F *reco_phi_mass =
      new TH1F("reco_phi_mass",
               "Run 19 A+A Invariant Mass Reconstructed #phi Meson ;K^{+}K^{-} "
               "Pair Invariant Mass(GeV/c^{2});counts",
               100, 1., 1.04);

  TH1F *reco_phi_pt =
      new TH1F("reco_phi_pt",
               "Run 19 A+A Transverse Momentum #phi Meson ;K^{+}K^{-} Pair "
               "Transverse Momentum(GeV/c^{2});counts",
               100, 0., 1.);

  TH1F *reco_phi_eta = new TH1F("reco_phi_eta",
                                "Run 19 A+A #eta Reconstructed #phi Meson "
                                ";K^{+}K^{-} Pair #eta(GeV/c^{2});counts",
                                100, -2, 2);

  // Open the file containing the tree (INPUT data).
  TFile *myFile = TFile::Open("~/Documents/Research Data/Run_19_Au_Au.root");

  TFile *output_hist = new TFile("kaonDistribution.root", "Recreate");

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
    lv = lv1 + lv2;

    // kaon candidate filling with PId selection
    if (abs(pair->d1_mNSigmaPion - NSigmaPionShift) > 5 &&
        abs(pair->d2_mNSigmaPion - NSigmaPionShift) > 5 &&
        abs(pair->d1_mNSigmaElectron - NSigmaElectronShift) > 5 &&
        abs(pair->d2_mNSigmaElectron - NSigmaElectronShift) > 5 &&
        abs(pair->d1_mNSigmaProton - NSigmaProtonShift) > 5 &&
        abs(pair->d2_mNSigmaProton - NSigmaProtonShift) > 5 &&
        abs(pair->d1_mNSigmaKaon - NSigmaKaonShift) < 5 &&
        abs(pair->d2_mNSigmaKaon - NSigmaKaonShift) < 5) {

      // reco level filling
      if (pair->mChargeSum == 0 && lv.M() >= 1. && lv.M() < 1.04) {
        reco_phi_mass->Fill(lv.M());
        reco_phi_pt->Fill(lv.Pt());
        reco_phi_eta->Fill(lv.Eta());
      }

    } 
  } // loop on events
  myFile->Close();
  reco_phi_pt->Write();
  reco_phi_mass->Write();
  reco_phi_eta->Write();
  output_hist->Close();
  reco_phi_pt->Delete();
  reco_phi_mass->Delete();
  reco_phi_eta->Delete();
}
