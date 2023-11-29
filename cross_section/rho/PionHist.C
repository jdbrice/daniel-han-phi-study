
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

void PionHist() {

  // histogram for reco phi
  TH1F *reco_rho_mass =
      new TH1F("reco_rho_mass",
               "Run 19 A+A Invariant Mass Reconstructed #rho Meson ;#pi^{+}#pi^{-} "
               "Pair Invariant Mass(GeV/c^{2});counts",
               100, 0.67, 0.87);

  TH1F *reco_rho_pt =
      new TH1F("reco_rho_pt",
               "Run 19 A+A Transverse Momentum #rho Meson ;#pi^{+}#pi^{-} Pair "
               "Transverse Momentum(GeV/c^{2});counts",
               100, 0., 1.);

  TH1F *reco_rho_y = new TH1F("reco_rho_y",
                                "Run 19 A+A #y Reconstructed #rho Meson "
                                ";#pi^{+}#pi^{-} Pair #y;counts",
                                100, -1, 1);

  // Open the file containing the tree (INPUT data).
  TFile *myFile = TFile::Open("~/Documents/Research Data/Run_19_Au_Au.root");

  TFile *output_hist = new TFile("pionDistribution.root", "Recreate");

  // This setup the reader, access the data
  TTreeReader myReader("PairDst", myFile);

  // We need this so ROOT understands the data format
  gSystem->Load("FemtoPairFormat_h.so");
  TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

  // Lorentz vectors (4-vectors) for kinematics in special relativity
  TLorentzVector lv1, lv2, lv;

  // Loop over all entries of the TTree or TChain.
  while (myReader.Next()) {

    lv1.SetPtEtaPhiM(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.139);
    lv2.SetPtEtaPhiM(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.139);
    lv = lv1 + lv2;

    // kaon candidate filling with PId selection
    if (abs(pair->d1_mNSigmaPion - NSigmaPionShift) < 5 &&
        abs(pair->d2_mNSigmaPion - NSigmaPionShift) < 5 &&
        abs(pair->d1_mNSigmaElectron - NSigmaElectronShift) > 5 &&
        abs(pair->d2_mNSigmaElectron - NSigmaElectronShift) > 5 &&
        abs(pair->d1_mNSigmaProton - NSigmaProtonShift) > 5 &&
        abs(pair->d2_mNSigmaProton - NSigmaProtonShift) > 5 &&
        abs(pair->d1_mNSigmaKaon - NSigmaKaonShift) > 5 &&
        abs(pair->d2_mNSigmaKaon - NSigmaKaonShift) > 5) {

      // reco level filling
      if (pair->mChargeSum == 0 && lv.M() >= 0.67 && lv.M() < 0.87) {
        reco_rho_mass->Fill(lv.M());
        reco_rho_pt->Fill(lv.Pt());
        reco_rho_y->Fill(lv.Rapidity());
      }

    } 
  } // loop on events
  myFile->Close();
  reco_rho_pt->Write();
  reco_rho_mass->Write();
  reco_rho_y->Write();
  output_hist->Close();
  reco_rho_pt->Delete();
  reco_rho_mass->Delete();
  reco_rho_y->Delete();
}
