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

void daughter_distribution_proton() {

  // histogram for pion canidataes
  TH1F *proton_DCA =
      new TH1F("Proton Candidates",
               "Run 19 A+A DCA Proton Candidates;DCA;counts", 100, 0., 3.);

  TH1F *proton_pt = new TH1F(
      "Proton Candidates",
      "Run 19 A+A P_{T} Proton Candidates;P_{T}(GeV/c);counts", 100, 0., 0.8);

  TH1F *proton_eta =
      new TH1F("Proton Candidates",
               "Run 19 A+A #eta Proton Candidates;#eta;counts", 100, -2., 2.);
  TH1F *proton_phi = new TH1F(
      "Proton Candidates", "Run 19 A+A #phi Proton Candidates;#phi(rad);counts",
      100, -3.15, 3.15);
  TH1F *proton_Nproton = new TH1F(
      "Proton Candidates",
      "Run 19 A+A N#sigmaProton Proton Candidates;N#sigmaProton;counts", 100,
      -6, 8);

  TH1F *proton_Nkaon =
      new TH1F("Proton Cnadidates",
               "Run 19 A+A N#sigmaKaon Proton Candidates;N#sigmaKaon;counts",
               100, 0, 40);

  TH1F *proton_plus_DCA =
      new TH1F("Proton_plus Candidates",
               "Run 19 A+A DCA Proton Candidates;DCA;counts", 100, 0., 3.);

  TH1F *proton_plus_pt = new TH1F(
      "Proton_plus Candidates",
      "Run 19 A+A P_{T} Proton Candidates;P_{T}(GeV/c);counts", 100, 0., 0.8);

  TH1F *proton_plus_eta =
      new TH1F("Proton_plus Candidates",
               "Run 19 A+A #eta Proton Candidates;#eta;counts", 100, -2., 2.);
  TH1F *proton_plus_phi =
      new TH1F("Proton Candidates",
               "Run 19 A+A #phi proton_plus Candidates;#phi(rad);counts", 100,
               -3.15, 3.15);
  TH1F *proton_plus_Nproton = new TH1F(
      "proton_plus Candidates",
      "Run 19 A+A N#sigmaProton_plus Proton Candidates;N#sigmaProton;counts",
      100, -6, 8);

  TH1F *proton_plus_Nkaon = new TH1F(
      "Proton_plus Cnadidates",
      "Run 19 A+A N#sigmaProton_plus Proton Candidates;N#sigmaKaon;counts", 100,
      0, 40);

  TH1F *proton_minus_DCA =
      new TH1F("Proton_minus Candidates",
               "Run 19 A+A DCA Proton Candidates;DCA;counts", 100, 0., 3.);

  TH1F *proton_minus_pt = new TH1F(
      "Proton_minus Candidates",
      "Run 19 A+A P_{T} Proton Candidates;P_{T}(GeV/c);counts", 100, 0., 0.8);

  TH1F *proton_minus_eta =
      new TH1F("Proton_minus Candidates",
               "Run 19 A+A #eta Proton Candidates;#eta;counts", 100, -2., 2.);
  TH1F *proton_minus_phi =
      new TH1F("Proton Candidates",
               "Run 19 A+A #phi Proton_minus Candidates;#phi(rad);counts", 100,
               -3.15, 3.15);
  TH1F *proton_minus_Nproton = new TH1F(
      "Proton_minus Candidates",
      "Run 19 A+A N#sigmaProton_minus Proton Candidates;N#sigmaProton;counts",
      100, -6, 8);

  TH1F *proton_minus_Nkaon = new TH1F(
      "Proton_minus Cnadidates",
      "Run 19 A+A N#sigmaProton_minus Proton Candidates;N#sigmaKaon;counts",
      100, 0, 40);

  TH1F *proton_pm_ratio =
      new TH1F("Proton Candidates",
               "Run 19 A+A Charge Sum Proton Candidates;Charge Sum;counts", 27,
               -2.2, 2.2);

  // histogram for same charge
  TH1F *PropProp_mass = new TH1F(
      "p^{+} p^{+}",
      "Run 19 A+A Invariant Mass p^{+} p^{+} ;Invariant Mass(GeV/c^2);counts",
      100, 0.9, 2.0);

  TH1F *PropProp_pt = new TH1F(
      "p^{+}p^{+}", "Run 19 A+A P_{T} p^{+} p^{+} ;P_{T}(GeV/c);counts", 100,
      0., 1.6);

  TH1F *PropProp_eta = new TH1F(
      "p^{+}p^{+} ", "Run 19 A+A #eta p^{+} p^{+} ;#eta;counts", 100, -4., 4.);

  TH1F *PropProp_phi =
      new TH1F("KpKp ", "Run 19 A+A #phi p^{+}p^{+} ;#phi(rad);counts", 100,
               -3.15, 3.15);

  TH1F *PromProm_mass = new TH1F("KmKm ",
                                 "Run 19 A+A Invariant Mass #bar{p} #bar{p} "
                                 ";Invariant Mass(GeV/c^2);counts",
                                 100, 0.9, 2.00);

  TH1F *PromProm_pt = new TH1F(
      "KmKm ", "Run 19 A+A P_{T} #bar{p} #bar{p} ;P_e z{T}(GeV/c);counts", 100,
      0., 1.6);

  TH1F *PromProm_eta =
      new TH1F("#bar{p} #bar{p} ",
               "Run 19 A+A #eta #bar{p} #bar{p} ;#eta;counts", 100, -4., 4.);

  TH1F *PromProm_phi =
      new TH1F("KmKm ", "Run 19 A+A #phi #bar{p} #bar{p} ;#phi(rad);counts",
               100, -3.15, 3.15);

  TH1F *like_sign_mass = new TH1F("like_sign ",
                                  "Run 19 A+A Invariant Mass Like Sign Pions "
                                  ";Invariant Mass(GeV/c^2);counts",
                                  100, 0.9, 2.0);

  TH1F *like_sign_pt = new TH1F(
      "like_sign ", "Run 19 A+A P_{T} Like Sign Protons ;P_{T}(GeV/c);counts",
      100, 0., 1.6);

  TH1F *like_sign_eta =
      new TH1F("like_sign ", "Run 19 A+A #eta Like Sign Protons ;#eta;counts",
               100, -4., 4.);

  TH1F *like_sign_phi = new TH1F(
      "like_sign ", "Run 19 A+A #phi Like Sign Protons ;#phi(rad);counts", 100,
      -3.15, 3.15);


  TH1F *reco_pair_mass = new TH1F("Reco Proton Pair",
                                  "Run 19 A+A Invariant Mass Reco Proton "
                                  "Pair;pp Invariant Mass(GeV/c^2);counts",
                                  100, 0.9, 2.0);

  TH2F *reco_pair_mass_pt = new TH2F(
      "Reco Proton Pair",
      "Run 19 A+A Invariant Mass Reco Proton Pair and P_{T} Reco Proton Pair;pp"
      "Invariant Mass(GeV/c^2);P_{T}",
      100, 0.9, 2.0, 100, 0, 0.9);

  TH1F *reco_pair_pt = new TH1F(
      "Reco Proton Pair",
      "Run 19 A+A P_{T} Reco Proton Pair;P_{T}(GeV/c);counts", 100, 0., 0.9);

  TH1F *reco_pair_eta =
      new TH1F("Reco Proton Pair",
               "Run 19 A+A #eta Reco Proton Pair;#eta;counts", 100, -4., 4.);

  TH1F *reco_pair_phi = new TH1F(
      "Reco Proton Pair", "Run 19 A+A #phi Reco Proton Pair;#phi(rad);counts",
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
        (pair->d1_mNSigmaElectron > 5.424112 ||
         pair->d1_mNSigmaElectron < -4.575888) &&
        (pair->d2_mNSigmaElectron > 5.424112 ||
         pair->d2_mNSigmaElectron < -4.575888) &&
        (pair->d1_mNSigmaProton < 7.36 && pair->d1_mNSigmaProton > -2.64) &&
        (pair->d2_mNSigmaProton < 7.36 && pair->d2_mNSigmaProton > -2.64) &&
        (pair->d1_mNSigmaKaon > 7.4 || pair->d1_mNSigmaKaon < -2.6) &&
        (pair->d2_mNSigmaKaon > 7.4 || pair->d2_mNSigmaKaon < -2.6)) {

      proton_Nkaon->Fill(pair->d1_mNSigmaKaon);
      proton_Nkaon->Fill(pair->d2_mNSigmaKaon);

      proton_Nproton->Fill(pair->d1_mNSigmaProton);
      proton_Nproton->Fill(pair->d2_mNSigmaProton);

      // reco level filling
      if (pair->mChargeSum == 0) {
        proton_pm_ratio->Fill(pair->mChargeSum);
        reco_pair_pt->Fill(lv.Pt());
        reco_pair_eta->Fill(lv.Eta());
        reco_pair_phi->Fill(lv.Phi());
        reco_pair_mass->Fill(lv.M());
        reco_pair_mass_pt->Fill(lv.M(), lv.Pt());

        proton_pt->Fill(pair->d1_mPt);
        proton_pt->Fill(pair->d2_mPt);

        proton_DCA->Fill(pair->d1_mDCA);
        proton_DCA->Fill(pair->d2_mDCA);

        proton_eta->Fill(pair->d1_mEta);
        proton_eta->Fill(pair->d2_mEta);

        proton_phi->Fill(pair->d1_mPhi);
        proton_phi->Fill(pair->d2_mPhi);

      } else if (pair->mChargeSum == 2) {
        proton_pm_ratio->Fill(pair->mChargeSum);
        PropProp_pt->Fill(lv.Pt());
        PropProp_mass->Fill(lv.M());
        PropProp_eta->Fill(lv.Eta());
        PropProp_phi->Fill(lv.Phi());

        proton_plus_pt->Fill(pair->d1_mPt);
        proton_plus_pt->Fill(pair->d2_mPt);

        proton_plus_DCA->Fill(pair->d1_mDCA);
        proton_plus_DCA->Fill(pair->d2_mDCA);

        proton_plus_eta->Fill(pair->d1_mEta);
        proton_plus_eta->Fill(pair->d2_mEta);

        proton_plus_phi->Fill(pair->d1_mPhi);
        proton_plus_phi->Fill(pair->d2_mPhi);

        proton_plus_Nkaon->Fill(pair->d1_mNSigmaKaon);
        proton_plus_Nkaon->Fill(pair->d2_mNSigmaKaon);

        proton_plus_Nproton->Fill(pair->d1_mNSigmaPion);
        proton_plus_Nproton->Fill(pair->d2_mNSigmaPion);

        like_sign_pt->Fill(lv.Pt());
        like_sign_mass->Fill(lv.M());
        like_sign_eta->Fill(lv.Eta());
        like_sign_phi->Fill(lv.Phi());

      } else if (pair->mChargeSum == -2) {
        proton_pm_ratio->Fill(pair->mChargeSum);
        PromProm_pt->Fill(lv.Pt());
        PromProm_mass->Fill(lv.M());
        PromProm_eta->Fill(lv.Eta());
        PromProm_phi->Fill(lv.Phi());

        proton_minus_pt->Fill(pair->d1_mPt);
        proton_minus_pt->Fill(pair->d2_mPt);

        proton_minus_DCA->Fill(pair->d1_mDCA);
        proton_minus_DCA->Fill(pair->d2_mDCA);

        proton_minus_eta->Fill(pair->d1_mEta);
        proton_minus_eta->Fill(pair->d2_mEta);

        proton_minus_phi->Fill(pair->d1_mPhi);
        proton_minus_phi->Fill(pair->d2_mPhi);

        proton_minus_Nkaon->Fill(pair->d1_mNSigmaKaon);
        proton_minus_Nkaon->Fill(pair->d2_mNSigmaKaon);

        proton_minus_Nproton->Fill(pair->d1_mNSigmaPion);
        proton_minus_Nproton->Fill(pair->d2_mNSigmaPion);

        like_sign_pt->Fill(lv.Pt());
        like_sign_mass->Fill(lv.M());
        like_sign_eta->Fill(lv.Eta());
        like_sign_phi->Fill(lv.Phi());
      }
    }

  } // loop on events

  makeCan();
  proton_pm_ratio->Draw();
  gPad->Print("./Plots_sigma/Protons/proton_charge.png");

  makeCan();
  proton_DCA->Draw();
  gPad->Print("./Plots_sigma/Protons/proton_DCA.png");
  makeCan();
  proton_pt->Draw();
  gPad->Print("./Plots_sigma/Protons/proton_pt.png");
  makeCan();
  proton_eta->Draw();
  gPad->Print("./Plots_sigma/Protons/proton_eta.png");
  makeCan();
  proton_phi->Draw();
  gPad->Print("./Plots_sigma/Protons/proton_phi.png");
  makeCan();
  proton_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Protons/proton_NSigmaKaon.png");
  makeCan();
  proton_Nproton->Draw();
  gPad->Print("./Plots_sigma/Protons/proton_NSigmaPion.png");

  makeCan();
  reco_pair_mass->Draw();
  gPad->Print("./Plots_sigma/Protons/Reco/proton_pair_mass.png");
  makeCan();
  reco_pair_pt->Draw();
  gPad->Print("./Plots_sigma/Protons/Reco/proton_pair_pt.png");
  makeCan();
  reco_pair_eta->Draw();
  gPad->Print("./Plots_sigma/Protons/Reco/proton_pair_eta.png");
  makeCan();
  reco_pair_phi->Draw();
  gPad->Print("./Plots_sigma/Protons/Reco/proton_pair_phi.png");
  makeCan();
  reco_pair_mass_pt->Draw("colz");
  gPad->Print("./Plots_sigma/Protons/Reco/proton_pair_mass_pt.png");

  makeCan();
  proton_plus_DCA->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_plus_DCA.png");
  makeCan();
  proton_plus_pt->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_plus_pt.png");
  makeCan();
  proton_plus_eta->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_plus_eta.png");
  makeCan();
  proton_plus_phi->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_plus_phi.png");
  makeCan();
  proton_plus_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_plus_NSigmaKaon.png");
  makeCan();
  proton_plus_Nproton->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_plus_NSigmaPion.png");

  makeCan();
  proton_minus_DCA->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_minus_DCA.png");
  makeCan();
  proton_minus_pt->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_minus_pt.png");
  makeCan();
  proton_minus_eta->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_minus_eta.png");
  makeCan();
  proton_minus_phi->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_minus_phi.png");
  makeCan();
  proton_minus_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_minus_NSigmaKaon.png");
  makeCan();
  proton_minus_Nproton->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/proton_minus_NSigmaPion.png");

  makeCan();
  PropProp_pt->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/p+p+_pt.png");
  makeCan();
  PropProp_eta->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/p+p+_eta.png");
  makeCan();
  PropProp_phi->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/p+p+_phi.png");
  makeCan();
  PropProp_mass->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/p+p+_mass.png");
  makeCan();
  PromProm_pt->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/p-p-_pt.png");
  makeCan();
  PromProm_eta->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/p-p-_eta.png");
  makeCan();
  PromProm_phi->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/p-p-_phi.png");
  makeCan();
  PromProm_mass->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/p-p-_mass.png");

  makeCan();
  like_sign_pt->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/proton_ls_pt.png");
  makeCan();
  like_sign_mass->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/proton_ls_mass.png");
  makeCan();
  like_sign_eta->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/proton_ls_eta.png");
  makeCan();
  like_sign_phi->Draw();
  gPad->Print("./Plots_sigma/Protons/like_signs/Reco/proton_ls_phi.png");
}
