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

void pt_resolution_processor() {

  TFile *fo = new TFile("processed_phi_exp_data.root", "RECREATE");

  TH1F *error =
      new TH1F("Phi in Error Range", "Run19 Au+Au data;rc_mass (GeV); counts",
               100, 1.04, 1.1);

  TH1F *rc_mass = new TH1F(
      "Phi Selected", "Phi Meson Mass Distribution;rc_mass (GeV); Probability",
      100, 1., 1.04);

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

        if (lv.M() > 1.04 && lv.M() <= 1.1) {
          error->Fill(lv.M());
        }
      }
    } // selection
  }   // loop on events

  fo->cd();
  makeCan();
  error->Fit("pol0");
  TF1 *error_fit = (TF1 *)error->GetListOfFunctions()->FindObject("pol0");
  double error_mean = error_fit->GetParameter(0);
  TF1 *error_fit_result =
      new TF1("error_estimate", std::to_string(error_mean).c_str(), 0, 10.);
   rc_mass->Add(error_fit_result, -1);
   rc_mass->Scale(1. / rc_mass->Integral());
  rc_mass->Draw("HIST");
   fo->Write();
}
