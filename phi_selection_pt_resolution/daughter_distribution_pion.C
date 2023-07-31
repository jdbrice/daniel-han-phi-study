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

void daughter_distribution_pion() {

  // histogram for pion canidataes
  TH1F *pion_DCA =
      new TH1F("Pion Candidates", "Run 19 A+A DCA Pion Candidates;DCA;counts",
               100, 0., 3.);

  TH1F *pion_pt = new TH1F(
      "Pion Candidates", "Run 19 A+A P_{T} Pion Candidates;P_{T}(GeV/c);counts",
      100, 0., 0.8);

  TH1F *pion_eta =
      new TH1F("Pion Candidates", "Run 19 A+A #eta Pion Candidates;#eta;counts",
               100, -2., 2.);
  TH1F *pion_phi = new TH1F("Pion Candidates",
                            "Run 19 A+A #phi Pion Candidates;#phi(rad);counts",
                            100, -3.15, 3.15);
  TH1F *pion_Npion = new TH1F(
      "Pion Candidates",
      "Run 19 A+A N#sigmaPion Pion Candidates;N#sigmaPion;counts", 100, -6, 8);

  TH1F *pion_Nkaon = new TH1F(
      "pion Cnadidates",
      "Run 19 A+A N#sigmapion Pion Candidates;N#sigmaKaon;counts", 100, -40, 0);

  TH1F *pion_background = new TH1F("Pion Cnadidates",
                                   "Run 19 A+A Invariant Mass "
                                   "Background;Invaraint Mass(GeV/c^2);counts",
                                   100, 1., 1.1);
  TH1F *electron_background =
      new TH1F("Electron Cnadidates",
               "Run 19 A+A Invariant Mass Electron "
               "Background;Invaraint Mass(GeV/c^2);counts",
               100, 1., 1.1);
  TH1F *proton_background =
      new TH1F("Proton Cnadidates",
               "Run 19 A+A Invariant Mass Proton "
               "Background;Invaraint Mass(GeV/c^2);counts",
               100, 1., 1.1);

  TH1F *pion_plus_DCA =
      new TH1F("pion_plus Candidates",
               "Run 19 A+A DCA Pion Candidates;DCA;counts", 100, 0., 3.);

  TH1F *pion_plus_pt = new TH1F(
      "pion_plus Candidates",
      "Run 19 A+A P_{T} Pion Candidates;P_{T}(GeV/c);counts", 100, 0., 0.8);

  TH1F *pion_plus_eta =
      new TH1F("pion_plus Candidates",
               "Run 19 A+A #eta Pion Candidates;#eta;counts", 100, -2., 2.);
  TH1F *pion_plus_phi =
      new TH1F("Pion Candidates",
               "Run 19 A+A #phi pion_plus Candidates;#phi(rad);counts", 100,
               -3.15, 3.15);
  TH1F *pion_plus_Npion =
      new TH1F("pion_plus Candidates",
               "Run 19 A+A N#sigmapion_plus Pion Candidates;N#sigmaPion;counts",
               100, -6, 8);

  TH1F *pion_plus_Nkaon =
      new TH1F("pion_plus Cnadidates",
               "Run 19 A+A N#sigmapion_plus Pion Candidates;N#sigmaKaon;counts",
               100, -40, 0);

  TH1F *pion_minus_DCA =
      new TH1F("pion_minus Candidates",
               "Run 19 A+A DCA Pion Candidates;DCA;counts", 100, 0., 3.);

  TH1F *pion_minus_pt = new TH1F(
      "pion_minus Candidates",
      "Run 19 A+A P_{T} Pion Candidates;P_{T}(GeV/c);counts", 100, 0., 0.8);

  TH1F *pion_minus_eta =
      new TH1F("pion_minus Candidates",
               "Run 19 A+A #eta Pion Candidates;#eta;counts", 100, -2., 2.);
  TH1F *pion_minus_phi =
      new TH1F("Pion Candidates",
               "Run 19 A+A #phi pion_minus Candidates;#phi(rad);counts", 100,
               -3.15, 3.15);
  TH1F *pion_minus_Npion = new TH1F(
      "pion_minus Candidates",
      "Run 19 A+A N#sigmapion_minus Pion Candidates;N#sigmaPion;counts", 100,
      -6, 8);

  TH1F *pion_minus_Nkaon = new TH1F(
      "pion_minus Cnadidates",
      "Run 19 A+A N#sigmapion_minus Pion Candidates;N#sigmaKaon;counts", 100,
      -40, 0);

  TH1F *pion_pm_ratio = new TH1F(
      "pion Candidates",
      "Run 19 A+A Charge Sum Pion Candidates;Charge Sum;counts", 27, -2.2, 2.2);

  // histogram for same charge
  TH1F *PipPip_mass = new TH1F(
      "#pi+#pi+",
      "Run 19 A+A Invariant Mass #pi+#pi+ ;Invariant Mass(GeV/c^2);counts", 100,
      0.9, 2.0);

  TH1F *PipPip_pt =
      new TH1F("#pi+#pi+ ", "Run 19 A+A P_{T} #pi+#pi+ ;P_{T}(GeV/c);counts",
               100, 0., 1.6);

  TH1F *PipPip_eta = new TH1F(
      "#pi+#pi+ ", "Run 19 A+A #eta #pi+#pi+ ;#eta;counts", 100, -4., 4.);

  TH1F *PipPip_phi = new TH1F(
      "KpKp ", "Run 19 A+A #phi #pi+#pi+ ;#phi(rad);counts", 100, -3.15, 3.15);

  TH1F *PimPim_mass = new TH1F(
      "KmKm ",
      "Run 19 A+A Invariant Mass #pi-#pi- ;Invariant Mass(GeV/c^2);counts", 100,
      0.9, 2.00);

  TH1F *PimPim_pt =
      new TH1F("KmKm ", "Run 19 A+A P_{T} #pi-#pi- ;P_e z{T}(GeV/c);counts",
               100, 0., 1.6);

  TH1F *PimPim_eta = new TH1F(
      "#pi-#pi- ", "Run 19 A+A #eta #pi-#pi- ;#eta;counts", 100, -4., 4.);

  TH1F *PimPim_phi = new TH1F(
      "KmKm ", "Run 19 A+A #phi #pi-#pi- ;#phi(rad);counts", 100, -3.15, 3.15);

  TH1F *like_sign_mass = new TH1F("like_sign ",
                                  "Run 19 A+A Invariant Mass Like Sign Pions "
                                  ";Invariant Mass(GeV/c^2);counts",
                                  100, 0.9, 2.0);

  TH1F *like_sign_pt = new TH1F(
      "like_sign ", "Run 19 A+A P_{T} Like Sign Pions ;P_{T}(GeV/c);counts",
      100, 0., 1.6);

  TH1F *like_sign_eta =
      new TH1F("like_sign ", "Run 19 A+A #eta Like Sign Pions ;#eta;counts",
               100, -4., 4.);

  TH1F *like_sign_phi = new TH1F(
      "like_sign ", "Run 19 A+A #phi Like Sign Pions ;#phi(rad);counts", 100,
      -3.15, 3.15);

  TH1F *reco_rho_mass = new TH1F(
      "reco_rho ",
      "Run 19 A+A Invariant Mass Reco rho ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.9, 2.0);

  TH2F *reco_rho_mass_pt =
      new TH2F("reco_rho",
               "Run 19 A+A Invariant Mass Reco rho and P_{T} Reco Phi;KK "
               "Invariant Mass(GeV/c^2);P_{T}",
               100, 0.9, 2.0, 100, 0, 0.9);

  TH1F *reco_rho_pt =
      new TH1F("reco_rho ", "Run 19 A+A P_{T} Reco rho ;P_{T}(GeV/c);counts",
               100, 0., 0.9);

  TH1F *reco_rho_eta = new TH1F(
      "reco_rho ", "Run 19 A+A #eta Reco rho ;#eta;counts", 100, -4., 4.);

  TH1F *reco_rho_phi =
      new TH1F("reco_rho ", "Run 19 A+A #phi Reco rho ;#phi(rad);counts", 100,
               -3.15, 3.15);
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
    if ((pair->d1_mNSigmaPion < 5.8 && pair->d1_mNSigmaPion > -4.2) &&
        (pair->d2_mNSigmaPion < 5.8 && pair->d2_mNSigmaPion > -4.2) &&
        (pair->d1_mNSigmaElectron > 5.424112 ||
         pair->d1_mNSigmaElectron < -4.575888) &&
        (pair->d2_mNSigmaElectron > 5.424112 ||
         pair->d2_mNSigmaElectron < -4.575888) &&
        (pair->d1_mNSigmaProton > 7.36 || pair->d1_mNSigmaProton < -2.64) &&
        (pair->d2_mNSigmaProton > 7.36 || pair->d2_mNSigmaProton < -2.64) &&
        (pair->d1_mNSigmaKaon > 7.4 || pair->d1_mNSigmaKaon < -2.6) &&
        (pair->d2_mNSigmaKaon > 7.4 || pair->d2_mNSigmaKaon < -2.6)) {

      pion_Nkaon->Fill(pair->d1_mNSigmaKaon);
      pion_Nkaon->Fill(pair->d2_mNSigmaKaon);

      pion_Npion->Fill(pair->d1_mNSigmaPion);
      pion_Npion->Fill(pair->d2_mNSigmaPion);

      // reco level filling
      if (pair->mChargeSum == 0) {
        pion_pm_ratio->Fill(pair->mChargeSum);
        reco_rho_pt->Fill(lv.Pt());
        reco_rho_eta->Fill(lv.Eta());
        reco_rho_phi->Fill(lv.Phi());
        reco_rho_mass->Fill(lv.M());
        reco_rho_mass_pt->Fill(lv.M(), lv.Pt());

        reco_rho_mass->Fill(lv.M());

        pion_pt->Fill(pair->d1_mPt);
        pion_pt->Fill(pair->d2_mPt);

        pion_DCA->Fill(pair->d1_mDCA);
        pion_DCA->Fill(pair->d2_mDCA);

        pion_eta->Fill(pair->d1_mEta);
        pion_eta->Fill(pair->d2_mEta);

        pion_phi->Fill(pair->d1_mPhi);
        pion_phi->Fill(pair->d2_mPhi);

      } else if (pair->mChargeSum == 2) {
        pion_pm_ratio->Fill(pair->mChargeSum);
        PipPip_pt->Fill(lv.Pt());
        PipPip_mass->Fill(lv.M());
        PipPip_eta->Fill(lv.Eta());
        PipPip_phi->Fill(lv.Phi());

        pion_plus_pt->Fill(pair->d1_mPt);
        pion_plus_pt->Fill(pair->d2_mPt);

        pion_plus_DCA->Fill(pair->d1_mDCA);
        pion_plus_DCA->Fill(pair->d2_mDCA);

        pion_plus_eta->Fill(pair->d1_mEta);
        pion_plus_eta->Fill(pair->d2_mEta);

        pion_plus_phi->Fill(pair->d1_mPhi);
        pion_plus_phi->Fill(pair->d2_mPhi);

        pion_plus_Nkaon->Fill(pair->d1_mNSigmaKaon);
        pion_plus_Nkaon->Fill(pair->d2_mNSigmaKaon);

        pion_plus_Npion->Fill(pair->d1_mNSigmaPion);
        pion_plus_Npion->Fill(pair->d2_mNSigmaPion);

        like_sign_pt->Fill(lv.Pt());
        like_sign_mass->Fill(lv.M());
        like_sign_eta->Fill(lv.Eta());
        like_sign_phi->Fill(lv.Phi());

      } else if (pair->mChargeSum == -2) {
        pion_pm_ratio->Fill(pair->mChargeSum);
        PimPim_pt->Fill(lv.Pt());
        PimPim_mass->Fill(lv.M());
        PimPim_eta->Fill(lv.Eta());
        PimPim_phi->Fill(lv.Phi());

        pion_minus_pt->Fill(pair->d1_mPt);
        pion_minus_pt->Fill(pair->d2_mPt);

        pion_minus_DCA->Fill(pair->d1_mDCA);
        pion_minus_DCA->Fill(pair->d2_mDCA);

        pion_minus_eta->Fill(pair->d1_mEta);
        pion_minus_eta->Fill(pair->d2_mEta);

        pion_minus_phi->Fill(pair->d1_mPhi);
        pion_minus_phi->Fill(pair->d2_mPhi);

        pion_minus_Nkaon->Fill(pair->d1_mNSigmaKaon);
        pion_minus_Nkaon->Fill(pair->d2_mNSigmaKaon);

        pion_minus_Npion->Fill(pair->d1_mNSigmaPion);
        pion_minus_Npion->Fill(pair->d2_mNSigmaPion);

        like_sign_pt->Fill(lv.Pt());
        like_sign_mass->Fill(lv.M());
        like_sign_eta->Fill(lv.Eta());
        like_sign_phi->Fill(lv.Phi());
      }
    }

    if (pair->mChargeSum == 0 && lv.M() >= 1.0 && lv.M() < 1.1 &&
        (pair->d1_mNSigmaPion < 5.8 && pair->d1_mNSigmaPion > -4.2) &&
        (pair->d2_mNSigmaPion < 5.8 && pair->d2_mNSigmaPion > -4.2) &&
        (pair->d1_mNSigmaPion > 3.8 || pair->d1_mNSigmaPion < -2.2) &&
        (pair->d2_mNSigmaPion > 3.8 || pair->d2_mNSigmaPion < -2.2) &&
        (pair->d1_mNSigmaElectron > 5.424112 ||
         pair->d1_mNSigmaElectron < -4.575888) &&
        (pair->d2_mNSigmaElectron > 5.424112 ||
         pair->d2_mNSigmaElectron < -4.575888) &&
        (pair->d1_mNSigmaProton > 7.36 || pair->d1_mNSigmaProton < -2.64) &&
        (pair->d2_mNSigmaProton > 7.36 || pair->d2_mNSigmaProton < -2.64) &&
        // (pair->d1_mNSigmaKaon > -7.6 && pair->d1_mNSigmaKaon < 12.4) &&
        // (pair->d2_mNSigmaKaon > -7.6 && pair->d2_mNSigmaKaon < 12.4) &&
        (pair->d1_mNSigmaKaon > 7.4 || pair->d1_mNSigmaKaon < -2.6) &&
        (pair->d2_mNSigmaKaon > 7.4 || pair->d2_mNSigmaKaon < -2.6)) {

      pion_background->Fill(lv.M());
    }

    if (pair->mChargeSum == 0 && lv.M() >= 1.0 && lv.M() < 1.1 &&
        (pair->d1_mNSigmaPion > 5.8 || pair->d1_mNSigmaPion < -4.2) &&
        (pair->d2_mNSigmaPion > 5.8 || pair->d2_mNSigmaPion < -4.2) &&
        (pair->d1_mNSigmaElectron < 5.424112 &&
         pair->d1_mNSigmaElectron > -4.575888) &&
        (pair->d2_mNSigmaElectron < 5.424112 &&
         pair->d2_mNSigmaElectron > -4.575888) &&
        (pair->d1_mNSigmaElectron > 3.424112 ||
         pair->d1_mNSigmaElectron < -2.575888) &&
        (pair->d2_mNSigmaElectron > 3.424112 ||
         pair->d2_mNSigmaElectron < -2.575888) &&
        (pair->d1_mNSigmaProton > 7.36 || pair->d1_mNSigmaProton < -2.64) &&
        (pair->d2_mNSigmaProton > 7.36 || pair->d2_mNSigmaProton < -2.64) &&
        // (pair->d1_mNSigmaKaon > -7.6 && pair->d1_mNSigmaKaon < 12.4) &&
        // (pair->d2_mNSigmaKaon > -7.6 && pair->d2_mNSigmaKaon < 12.4) &&
        (pair->d1_mNSigmaKaon > 7.4 || pair->d1_mNSigmaKaon < -2.6) &&
        (pair->d2_mNSigmaKaon > 7.4 || pair->d2_mNSigmaKaon < -2.6)) {

      electron_background->Fill(lv.M());
    }
    if (pair->mChargeSum == 0 && lv.M() >= 1.0 && lv.M() < 1.1 &&
        (pair->d1_mNSigmaPion > 5.8 || pair->d1_mNSigmaPion < -4.2) &&
        (pair->d2_mNSigmaPion > 5.8 || pair->d2_mNSigmaPion < -4.2) &&
        (pair->d1_mNSigmaElectron < 5.424112 &&
         pair->d1_mNSigmaElectron > -4.575888) &&
        (pair->d2_mNSigmaElectron < 5.424112 &&
         pair->d2_mNSigmaElectron > -4.575888) &&
        (pair->d1_mNSigmaProton < 7.36 && pair->d1_mNSigmaProton > -2.64) &&
        (pair->d2_mNSigmaProton < 7.36 && pair->d2_mNSigmaProton > -2.64) &&
        (pair->d1_mNSigmaProton > 5.36 || pair->d1_mNSigmaProton < -0.64) &&
        (pair->d2_mNSigmaProton > 5.36 || pair->d2_mNSigmaProton < -0.64) &&
        // (pair->d1_mNSigmaKaon > -7.6 && pair->d1_mNSigmaKaon < 12.4) &&
        // (pair->d2_mNSigmaKaon > -7.6 && pair->d2_mNSigmaKaon < 12.4) &&
        (pair->d1_mNSigmaKaon > 7.4 || pair->d1_mNSigmaKaon < -2.6) &&
        (pair->d2_mNSigmaKaon > 7.4 || pair->d2_mNSigmaKaon < -2.6)) {

      proton_background->Fill(lv.M());
    }

  } // loop on events

  makeCan();
  pion_pm_ratio->Draw();
  gPad->Print("./Plots_sigma/Pions/pion_charge.png");

  makeCan();
  pion_DCA->Draw();
  gPad->Print("./Plots_sigma/Pions/pion_DCA.png");
  makeCan();
  pion_pt->Draw();
  gPad->Print("./Plots_sigma/Pions/pion_pt.png");
  makeCan();
  pion_eta->Draw();
  gPad->Print("./Plots_sigma/Pions/pion_eta.png");
  makeCan();
  pion_phi->Draw();
  gPad->Print("./Plots_sigma/Pions/pion_phi.png");
  makeCan();
  pion_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Pions/pion_NSigmaKaon.png");
  makeCan();
  pion_Npion->Draw();
  gPad->Print("./Plots_sigma/Pions/pion_NSigmaPion.png");

  makeCan();
  reco_rho_mass->Draw();
  gPad->Print("./Plots_sigma/Pions/Reco/rho_mass.png");
  makeCan();
  reco_rho_pt->Draw();
  gPad->Print("./Plots_sigma/Pions/Reco/rho_pt.png");
  makeCan();
  reco_rho_eta->Draw();
  gPad->Print("./Plots_sigma/Pions/Reco/rho_eta.png");
  makeCan();
  reco_rho_phi->Draw();
  gPad->Print("./Plots_sigma/Pions/Reco/rho_phi.png");
  makeCan();
  reco_rho_mass_pt->Draw("colz");
  gPad->Print("./Plots_sigma/Pions/Reco/rho_mass_pt.png");

  makeCan();
  pion_plus_DCA->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_plus_DCA.png");
  makeCan();
  pion_plus_pt->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_plus_pt.png");
  makeCan();
  pion_plus_eta->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_plus_eta.png");
  makeCan();
  pion_plus_phi->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_plus_phi.png");
  makeCan();
  pion_plus_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_plus_NSigmaKaon.png");
  makeCan();
  pion_plus_Npion->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_plus_NSigmaPion.png");

  makeCan();
  pion_minus_DCA->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_minus_DCA.png");
  makeCan();
  pion_minus_pt->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_minus_pt.png");
  makeCan();
  pion_minus_eta->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_minus_eta.png");
  makeCan();
  pion_minus_phi->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_minus_phi.png");
  makeCan();
  pion_minus_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_minus_NSigmaKaon.png");
  makeCan();
  pion_minus_Npion->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/pion_minus_NSigmaPion.png");

  makeCan();
  PipPip_pt->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pi+pi+_pt.png");
  makeCan();
  PipPip_eta->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pi+pi+_eta.png");
  makeCan();
  PipPip_phi->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pi+pi+_phi.png");
  makeCan();
  PipPip_mass->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pi+pi+_mass.png");
  makeCan();
  PimPim_pt->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pi-pi-_pt.png");
  makeCan();
  PimPim_eta->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pi-pi-_eta.png");
  makeCan();
  PimPim_phi->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pi-pi-_phi.png");
  makeCan();
  PimPim_mass->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pi-pi-_mass.png");

  makeCan();
  like_sign_pt->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pion_ls_pt.png");
  makeCan();
  like_sign_mass->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pion_ls_mass.png");
  makeCan();
  like_sign_eta->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pion_ls_eta.png");
  makeCan();
  like_sign_phi->Draw();
  gPad->Print("./Plots_sigma/Pions/like_signs/Reco/pion_ls_phi.png");

  makeCan();
  TLegend *legend_1 = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend_1->AddEntry(pion_background, "Pion Candidates");
  legend_1->AddEntry(proton_background, "Proton Candidates");
  legend_1->AddEntry(electron_background, "Electron Candidates");
  pion_background->Draw("pfc");
  proton_background->Draw("same;pfc");
  electron_background->Draw("same;pfc");
  legend_1->Draw("same");
  gPad->Print("./Plots_sigma/All/all_backgrounds.png");
}
