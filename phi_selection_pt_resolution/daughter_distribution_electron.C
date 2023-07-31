// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "FemtoPairFormat.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TH2.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <cmath>
#include <math.h>

// ClassImp( FemtoPair );

int ican = 0;
void makeCan() {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.06);
  can->SetRightMargin(0.01);
}

void daughter_distribution_electron() {

  // histogram for pion canidataes
  TH1F *electron_DCA =
      new TH1F("Electron Candidates",
               "Run 19 A+A DCA Electron Candidates;DCA;counts", 100, 0., 3.);

  TH1F *electron_pt = new TH1F(
      "Electron Candidates",
      "Run 19 A+A P_{T} Electron Candidates;P_{T}(GeV/c);counts", 100, 0., 0.8);

  TH1F *electron_debug_pt = new TH1F(
      "Electron Candidates P_{T} < 140 MeV",
      "Run 19 A+A P_{T} Electron Candidates;P_{T}(GeV/c);counts", 100, 0.05, 0.14);

  TH1F *electron_debug_eta =
      new TH1F("Electron Candidates P_{T} < 140 MeV",
               "Run 19 A+A #eta Electron Candidates;#eta;counts", 100, -2., 2.);

  TH1F *electron_debug_phi = new TH1F(
      "Electron Candidates P_{T} < 140 MeV", "Run 19 A+A #phi Electron Candidates;#phi(rad);counts",
      100, -3.15, 3.15);

  TH1F *electron_eta =
      new TH1F("Electron Candidates",
               "Run 19 A+A #eta Electron Candidates;#eta;counts", 100, -2., 2.);
  TH1F *electron_phi = new TH1F(
      "Electron Candidates", "Run 19 A+A #phi Electron Candidates;#phi(rad);counts",
      100, -3.15, 3.15);
  TH1F *electron_Nelectron = new TH1F(
      "Electron Candidates",
      "Run 19 A+A N#sigmaElectron Electron Candidates;N#sigmaElectron;counts", 100,
      -6, 8);

  TH1F *electron_Nkaon =
      new TH1F("Electron Cnadidates",
               "Run 19 A+A N#sigmaKaon Electron Candidates;N#sigmaKaon;counts",
               100, -40, 40);

  TH1F *electron_plus_DCA =
      new TH1F("Electron_plus Candidates",
               "Run 19 A+A DCA Electron Candidates;DCA;counts", 100, 0., 3.);

  TH1F *electron_plus_pt = new TH1F(
      "Electron_plus Candidates",
      "Run 19 A+A P_{T} Electron Candidates;P_{T}(GeV/c);counts", 100, 0., 0.8);

  TH1F *electron_plus_eta =
      new TH1F("Electron_plus Candidates",
               "Run 19 A+A #eta Electron Candidates;#eta;counts", 100, -2., 2.);

  TH1F *electron_plus_phi =
      new TH1F("Electron Candidates",
               "Run 19 A+A #phi Electron_plus Candidates;#phi(rad);counts", 100,
               -3.15, 3.15);

  TH1F *electron_plus_Nelectron = new TH1F(
      "electron_plus Candidates",
      "Run 19 A+A N#sigmaelectron_plus electron Candidates;N#sigmaelectron;counts",
      100, -6, 8);

  TH1F *electron_plus_Nkaon = new TH1F(
      "Electron_plus Cnadidates",
      "Run 19 A+A N#sigmaElectron_plus Electron Candidates;N#sigmaKaon;counts", 100,
      0, 40);

  TH1F *Electron_minus_DCA =
      new TH1F("Electron_minus Candidates",
               "Run 19 A+A DCA Electron Candidates;DCA;counts", 100, 0., 3.);

  TH1F *Electron_minus_pt = new TH1F(
      "Electron_minus Candidates",
      "Run 19 A+A P_{T} Electron Candidates;P_{T}(GeV/c);counts", 100, 0., 0.8);

  TH1F *electron_minus_eta =
      new TH1F("Electron_minus Candidates",
               "Run 19 A+A #eta Electron Candidates;#eta;counts", 100, -2., 2.);

  TH1F *electron_minus_phi =
      new TH1F("Electron Candidates",
               "Run 19 A+A #phi Electron_minus Candidates;#phi(rad);counts", 100,
               -3.15, 3.15);

  TH1F *electron_minus_Nelectron = new TH1F(
      "Electron_minus Candidates",
      "Run 19 A+A N#sigmaElectron_minus Electron Candidates;N#sigmaElectron;counts",
      100, -6, 8);

  TH1F *electron_minus_Nkaon = new TH1F(
      "Electron_minus Cnadidates",
      "Run 19 A+A N#sigmaElectron_minus Electron Candidates;N#sigmaKaon;counts",
      100, 0, 40);

  TH1F *electron_pm_ratio =
      new TH1F("Electron Candidates",
               "Run 19 A+A Charge Sum Electron Candidates;Charge Sum;counts", 27,
               -2.2, 2.2);

  // histogram for same charge
  TH1F *epep_mass = new TH1F(
      "#bar{e} #bar{e}",
      "Run 19 A+A Invariant Mass #bar{e} #bar{e} ;Invariant Mass(GeV/c^2);counts",
      100, 0.9, 2.0);

  TH1F *epep_pt = new TH1F(
      "#bar{e} #bar{e}", "Run 19 A+A P_{T} #bar{e} #bar{e} ;P_{T}(GeV/c);counts", 100,
      0., 1.6);

  TH1F *epep_eta = new TH1F(
      "#bar{e} #bar{e} ", "Run 19 A+A #eta #bar{e} #bar{e} ;#eta;counts", 100, -4., 4.);

  TH1F *epep_phi =
      new TH1F("#bar{e} #bar{e}", "Run 19 A+A #phi #bar{e} #bar{e};#phi(rad);counts", 100,
               -3.15, 3.15);

  TH1F *emem_mass = new TH1F("e e",
                                 "Run 19 A+A Invariant Mass e e "
                                 ";Invariant Mass(GeV/c^2);counts",
                                 100, 0.9, 2.00);

  TH1F *emem_pt = new TH1F(
      "e e", "Run 19 A+A P_{T} e e ;P_e z{T}(GeV/c);counts", 100,
      0., 1.6);

  TH1F *emem_eta =
      new TH1F("e e ",
               "Run 19 A+A #eta e e ;#eta;counts", 100, -4., 4.);

  TH1F *emem_phi =
      new TH1F("e e", "Run 19 A+A #phi e e ;#phi(rad);counts",
               100, -3.15, 3.15);

  TH1F *like_sign_mass = new TH1F("like_sign ",
                                  "Run 19 A+A Invariant Mass Like Sign Electrons"
                                  ";Invariant Mass(GeV/c^2);counts",
                                  100, 0.9, 2.0);

  TH1F *like_sign_pt = new TH1F(
      "like_sign ", "Run 19 A+A P_{T} Like Sign Electrons ;P_{T}(GeV/c);counts",
      100, 0., 1.6);

  TH1F *like_sign_eta =
      new TH1F("like_sign ", "Run 19 A+A #eta Like Sign Electrons ;#eta;counts",
               100, -4., 4.);

  TH1F *like_sign_phi = new TH1F(
      "like_sign ", "Run 19 A+A #phi Like Sign Electrons ;#phi(rad);counts", 100,
      -3.15, 3.15);


  TH1F *reco_pair_mass = new TH1F("Reco Electron Pair",
                                  "Run 19 A+A Invariant Mass Reco Electron "
                                  "Pair;e e Invariant Mass(GeV/c^2);counts",
                                  100, 0.9, 2.0);

  TH2F *reco_pair_mass_pt = new TH2F(
      "Reco Electron Pair",
      "Run 19 A+A Invariant Mass Reco Electron Pair and P_{T} Reco Electron Pair;ee"
      "Invariant Mass(GeV/c^2);P_{T}",
      100, 0.9, 2.0, 100, 0, 0.9);

  TH1F *reco_pair_pt = new TH1F(
      "Reco Electron Pair",
      "Run 19 A+A P_{T} Reco Electron Pair;P_{T}(GeV/c);counts", 100, 0., 0.9);

  TH1F *reco_pair_eta =
      new TH1F("Reco Electron Pair",
               "Run 19 A+A #eta Reco Electron Pair;#eta;counts", 100, -4., 4.);

  TH1F *reco_pair_phi = new TH1F(
      "Reco Electron Pair", "Run 19 A+A #phi Reco Electron Pair;#phi(rad);counts",
      100, -3.15, 3.15);
  // Open the file containing the tree (INPUT data).
  TFile *myFile = TFile::Open("input.root");

  // This setup the reader, access the data
  TTreeReader myReader("PairDst", myFile);

  // We need this so ROOT understands the data format
  gSystem->Load("FemtoPairFormat_h.so");
  TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

  TLorentzVector lv1, lv2, lv;
  // Loop over all entries of the TTree or TChain.
  while (myReader.Next()) {

    // reconstrauction of parent lorentz vector
    // assuming pion mass
    lv1.SetPtEtaPhiM(pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.493);
    lv2.SetPtEtaPhiM(pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.493);
    lv = lv1 + lv2;
    // kaon candidate filling with PId selection
    if ((pair->d1_mNSigmaPion > 5.8 || pair->d1_mNSigmaPion < -4.2) &&
        (pair->d2_mNSigmaPion > 5.8 || pair->d2_mNSigmaPion < -4.2) &&
        (pair->d1_mNSigmaElectron < 5.424112 &&
         pair->d1_mNSigmaElectron > -4.575888) &&
        (pair->d2_mNSigmaElectron < 5.424112 &&
         pair->d2_mNSigmaElectron > -4.575888) &&
        (pair->d1_mNSigmaProton > 7.36 || pair->d1_mNSigmaProton < -2.64) &&
        (pair->d2_mNSigmaProton > 7.36 || pair->d2_mNSigmaProton < -2.64) &&
        (pair->d1_mNSigmaKaon > 7.4 || pair->d1_mNSigmaKaon < -2.6) &&
        (pair->d2_mNSigmaKaon > 7.4 || pair->d2_mNSigmaKaon < -2.6)) {

      electron_Nkaon->Fill(pair->d1_mNSigmaKaon);
      electron_Nkaon->Fill(pair->d2_mNSigmaKaon);

      electron_Nelectron->Fill(pair->d1_mNSigmaElectron);
      electron_Nelectron->Fill(pair->d2_mNSigmaElectron);

      // reco level filling
      if (pair->mChargeSum == 0) {
        electron_pm_ratio->Fill(pair->mChargeSum);
        reco_pair_pt->Fill(lv.Pt());
        reco_pair_eta->Fill(lv.Eta());
        reco_pair_phi->Fill(lv.Phi());
        reco_pair_mass->Fill(lv.M());
        reco_pair_mass_pt->Fill(lv.M(), lv.Pt());

        electron_pt->Fill(pair->d1_mPt);
        electron_pt->Fill(pair->d2_mPt);

        electron_DCA->Fill(pair->d1_mDCA);
        electron_DCA->Fill(pair->d2_mDCA);

        electron_eta->Fill(pair->d1_mEta);
        electron_eta->Fill(pair->d2_mEta);

        electron_phi->Fill(pair->d1_mPhi);
        electron_phi->Fill(pair->d2_mPhi);

        if(pair->d1_mPt < 0.14){
          electron_debug_pt->Fill(pair->d1_mPt);
          electron_debug_eta->Fill(pair->d1_mEta);
          electron_debug_phi->Fill(pair->d1_mPhi);
        }

        if(pair->d2_mPt < 0.14){
          electron_debug_pt->Fill(pair->d2_mPt);
          electron_debug_eta->Fill(pair->d2_mEta);
          electron_debug_phi->Fill(pair->d2_mPhi);
        }


      } else if (pair->mChargeSum == 2) {
        electron_pm_ratio->Fill(pair->mChargeSum);
        epep_pt->Fill(lv.Pt());
        epep_mass->Fill(lv.M());
        epep_eta->Fill(lv.Eta());
        epep_phi->Fill(lv.Phi());

        electron_plus_pt->Fill(pair->d1_mPt);
        electron_plus_pt->Fill(pair->d2_mPt);

        electron_plus_DCA->Fill(pair->d1_mDCA);
        electron_plus_DCA->Fill(pair->d2_mDCA);

        electron_plus_eta->Fill(pair->d1_mEta);
        electron_plus_eta->Fill(pair->d2_mEta);

        electron_plus_phi->Fill(pair->d1_mPhi);
        electron_plus_phi->Fill(pair->d2_mPhi);

        electron_plus_Nkaon->Fill(pair->d1_mNSigmaKaon);
        electron_plus_Nkaon->Fill(pair->d2_mNSigmaKaon);

        electron_plus_Nelectron->Fill(pair->d1_mNSigmaElectron);
        electron_plus_Nelectron->Fill(pair->d2_mNSigmaElectron);

        like_sign_pt->Fill(lv.Pt());
        like_sign_mass->Fill(lv.M());
        like_sign_eta->Fill(lv.Eta());
        like_sign_phi->Fill(lv.Phi());

      } else if (pair->mChargeSum == -2) {
        electron_pm_ratio->Fill(pair->mChargeSum);
        emem_pt->Fill(lv.Pt());
        emem_mass->Fill(lv.M());
        emem_eta->Fill(lv.Eta());
        emem_phi->Fill(lv.Phi());

        Electron_minus_pt->Fill(pair->d1_mPt);
        Electron_minus_pt->Fill(pair->d2_mPt);

        Electron_minus_DCA->Fill(pair->d1_mDCA);
        Electron_minus_DCA->Fill(pair->d2_mDCA);

        electron_minus_eta->Fill(pair->d1_mEta);
        electron_minus_eta->Fill(pair->d2_mEta);

        electron_minus_phi->Fill(pair->d1_mPhi);
        electron_minus_phi->Fill(pair->d2_mPhi);

        electron_minus_Nkaon->Fill(pair->d1_mNSigmaKaon);
        electron_minus_Nkaon->Fill(pair->d2_mNSigmaKaon);

        electron_minus_Nelectron->Fill(pair->d1_mNSigmaElectron);
        electron_minus_Nelectron->Fill(pair->d2_mNSigmaElectron);

        like_sign_pt->Fill(lv.Pt());
        like_sign_mass->Fill(lv.M());
        like_sign_eta->Fill(lv.Eta());
        like_sign_phi->Fill(lv.Phi());
      }
    }

  } // loop on events

  makeCan();
  electron_pm_ratio->Draw();
  gPad->Print("./Plots_sigma/Electrons/proton_charge.png");

  makeCan();
  electron_DCA->Draw();
  gPad->Print("./Plots_sigma/Electrons/electron_DCA.png");
  makeCan();
  electron_pt->Draw();
  gPad->Print("./Plots_sigma/Electrons/electron_pt.png");
  makeCan();
  electron_eta->Draw();
  gPad->Print("./Plots_sigma/Electrons/electron_eta.png");
  makeCan();
  electron_phi->Draw();
  gPad->Print("./Plots_sigma/Electrons/electron_phi.png");
  makeCan();
  electron_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Electrons/electron_NSigmaKaon.png");
  makeCan();
  electron_Nelectron->Draw();
  gPad->Print("./Plots_sigma/Electrons/electron_NSigmaElectron.png");

  makeCan();
  reco_pair_mass->Draw();
  gPad->Print("./Plots_sigma/Electrons/Reco/electron_pair_mass.png");
  makeCan();
  reco_pair_pt->Draw();
  gPad->Print("./Plots_sigma/Electrons/Reco/electron_pair_pt.png");
  makeCan();
  reco_pair_eta->Draw();
  gPad->Print("./Plots_sigma/Electrons/Reco/electron_pair_eta.png");
  makeCan();
  reco_pair_phi->Draw();
  gPad->Print("./Plots_sigma/Electrons/Reco/electron_pair_phi.png");
  makeCan();
  reco_pair_mass_pt->Draw("colz");
  gPad->Print("./Plots_sigma/Electrons/Reco/electron_pair_mass_pt.png");

  makeCan();
  electron_plus_DCA->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_plus_DCA.png");
  makeCan();
  electron_plus_pt->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_plus_pt.png");
  makeCan();
  electron_plus_eta->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_plus_eta.png");
  makeCan();
  electron_plus_phi->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_plus_phi.png");
  makeCan();
  electron_plus_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_plus_NSigmaKaon.png");
  makeCan();
  electron_plus_Nelectron->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_plus_NSigmaElectron.png");

  makeCan();
  Electron_minus_DCA->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_minus_DCA.png");
  makeCan();
  Electron_minus_pt->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_minus_pt.png");
  makeCan();
  electron_minus_eta->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_minus_eta.png");
  makeCan();
  electron_minus_phi->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_minus_phi.png");
  makeCan();
  electron_minus_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_minus_NSigmaKaon.png");
  makeCan();
  electron_minus_Nelectron->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/electron_minus_NSigmaElectron.png");

  makeCan();
  epep_pt->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/e+e+_pt.png");
  makeCan();
  epep_eta->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/e+e+_eta.png");
  makeCan();
  epep_phi->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/e+e+_phi.png");
  makeCan();
  epep_mass->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/e+e+_mass.png");
  makeCan();
  emem_pt->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/e-e-_pt.png");
  makeCan();
  emem_eta->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/e-e-_eta.png");
  makeCan();
  emem_phi->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/e-e-_phi.png");
  makeCan();
  emem_mass->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/e-e-_mass.png");

  makeCan();
  like_sign_pt->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/electron_ls_pt.png");
  makeCan();
  like_sign_mass->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/electron_ls_mass.png");
  makeCan();
  like_sign_eta->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/electron_ls_eta.png");
  makeCan();
  like_sign_phi->Draw();
  gPad->Print("./Plots_sigma/Electrons/like_signs/Reco/electron_ls_phi.png");
  makeCan();
  electron_debug_pt->Draw();
  gPad->Print("./Plots_sigma/Electrons/electron_debug_pt.png");
  makeCan();
  electron_debug_eta->Draw();
  gPad->Print("./Plots_sigma/Electrons/electron_debug_eta.png");
  makeCan();
  electron_debug_phi->Draw();
  gPad->Print("./Plots_sigma/Electrons/electron_debug_phi.png");
}
