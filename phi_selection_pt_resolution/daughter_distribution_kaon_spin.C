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

void daughter_distribution_kaon() {

  // histograms for all tracks
  TH1F *all_DCA = new TH1F("All Tracks", "Run 19 A+A DCA All Tracks;DCA;counts",
                           100, 0., 3.);

  TH1F *all_pt =
      new TH1F("All Tracks", "Run 19 A+A P_{T} All Tracks;P_T(GeV/c);counts",
               100, 0., 1.8);
  TH1F *all_eta = new TH1F(
      "All Tracks", "Run 19 A+A #eta All Tracks;#eta;counts", 100, -4., 4.);
  TH1F *all_phi =
      new TH1F("All Tracks", "Run 19 A+A #phi All Tracks;#phi(rad);counts", 100,
               -3.15, 3.15);
  TH1F *all_Npion = new TH1F(
      "All Tracks", "Run 19 A+A N#sigma#pi All Tracks;(N#sigmaPion);counts",
      100, -10, 80);

  TH1F *all_Nkaon = new TH1F(
      "All Tracks", "Run 19 A+A N#sigmaK All Tracks;(N#sigmaKaon);counts", 100,
      -50, 80);

  TH1F *all_chisq_pipi= new TH1F(
      "All Tracks", "Run 19 A+A #chi^{2}_{#pi#pi} All Tracks;#chi^{2}_{#pi#pi};counts",
      100, -10, 100);

  TH1F *all_chisq_kk = new TH1F(
      "All Tracks", "Run 19 A+A #chi^{2}_{KK} All Tracks;#chi^{2}_{KK};counts",
      100, -10, 4000);

  TH1F *all_chisq_ee = new TH1F(
      "All Tracks", "Run 19 A+A #chi^{2}_{ee} All Tracks;#chi^{2}_{ee};counts",
      100, -10, 800);

  TH1F *all_chisq_pp = new TH1F(
      "All Tracks", "Run 19 A+A #chi^{2}_{pp} All Tracks;#chi^{2}_{pp};counts",
      100, -10, 400);

  TH2F *all_pt_NPion = new TH2F(
      "All Tracks", "Run 19 A+A N#sigmaPion All Tracks;P_{T}(GeV/c);N#sigmaPion", 100,
      0., 1.2, 200, -20, 80);

  TH2F *all_pt_NKaon = new TH2F(
      "All Tracks", "Run 19 A+A N#sigmaKaon All Tracks;P_{T}(GeV/c);N#sigmaKaon", 100,
      0., 1.2, 200, -60, 80);

  TH2F *all_pt_NElectron = new TH2F(
      "All Tracks", "Run 19 A+A N#sigmaElectron All Tracks;P_{T}(GeV/c);N#sigmaElectron", 100,
      0., 1.2, 200, -20, 80);

  TH2F *all_pt_NProton = new TH2F(
      "All Tracks", "Run 19 A+A N#sigmaProton All Tracks;P_{T}(GeV/c);N#sigmaProton", 100,
      0., 1.2, 200, -80, 60);


  TH2F *all_nsigmapion_nsigmakaon = new TH2F(
      "All Tracks",
      "Run 19 A+A N#sigmaPion N#sigmaKaon All Tracks;N#sigmaPion;N#sigmaKaon",
      300, -20, 80, 300, -60, 60);

  TH1F *all_pm_ratio = new TH1F(
      "All Trakcs",
      "Run 19 A+A Charge Sum All Tracks;Charge Sum;counts", 27, -2.2, 2.2);

  // histogram for kaon canidataes
  TH1F *kaon_DCA =
      new TH1F("Kaon Candidates", "Run 19 A+A DCA Kaon Candidates;DCA;counts",
               100, 0., 3.);

  TH1F *kaon_pt = new TH1F(
      "Kaon Candidates", "Run 19 A+A P_{T} Kaon Candidates;P_{T}(GeV/c);counts",
      100, 0., 0.5);
  TH1F *kaon_eta =
      new TH1F("Kaon Candidates", "Run 19 A+A #eta Kaon Candidates;#eta;counts",
               100, -2., 2.);
  TH1F *kaon_phi = new TH1F("Kaon Candidates",
                            "Run 19 A+A #phi Kaon Candidates;#phi(rad);counts",
                            100, -3.15, 3.15);
  TH1F *kaon_Npion = new TH1F(
      "Kaon Candidates",
      "Run 19 A+A N#sigmaPion Kaon Candidates;N#sigmaPion;counts", 100, 5, 30);

  TH1F *kaon_Nkaon = new TH1F(
      "Kaon Cnadidates",
      "Run 19 A+A N#sigmaKaon Kaon Candidates;N#sigmaKaon;counts", 100, -4, 12);

  TH1F *kaon_pm_ratio = new TH1F(
      "Kaon Candidates",
      "Run 19 A+A Charge Sum Kaon Candidates;Charge Sum;counts", 27, -2.2, 2.2);

  // histogram for reco phi
  TH1F *reco_phi_mass = new TH1F(
      "reco_phi ",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.9, 1.4);

  TH1F *reco_coherent_phi_mass = new TH1F(
      "reco_phi ",
      "Run 19 A+A Invariant Mass Reco Phi Cohehrent;KK Invariant Mass(GeV/c^2);counts",
      100, 1.0, 1.04);

  TH1F *reco_coherent_phi_mass_tof = new TH1F(
      "reco_phi ",
      "Run 19 A+A Invariant Mass Reco Phi Cohehrent;KK Invariant Mass(GeV/c^2);counts",
      100, 1.0, 1.04);

  TH1F *reco_phi_mass_pt_big_cut = new TH1F(
      "reco_phi P_{T} daughter < 400 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.9, 1.4);

  TH1F *reco_phi_mass_pt_cut = new TH1F(
      "reco_phi P_{T} daughter < 300 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.9, 1.4);

  TH1F *reco_phi_mass_pt_small_cut = new TH1F(
      "reco_phi P_{T} daughter < 200 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.9, 1.4);

  TH1F *reco_phi_mass_phi_pt_big_cut = new TH1F(
      "reco_phi P_{T} Pair P_{T} < 500 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.9, 1.4);

  TH1F *reco_phi_mass_phi_pt_cut = new TH1F(
      "reco_phi P_{T} Pair P_{T} < 300 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.9, 1.4);

  TH1F *reco_phi_mass_phi_pt_small_cut = new TH1F(
      "reco_phi P_{T} Pair P_{T} < 200 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.9, 1.4);

  TH1F *reco_phi_near_mass = new TH1F(
      "reco_phi",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 1., 1.04);

  TH1F *reco_phi_near_mass_pt_big_cut = new TH1F(
      "reco_phi P_{T} daughter < 400 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 1., 1.04);

  TH1F *reco_phi_near_mass_pt_cut = new TH1F(
      "reco_phi P_{T} daughter < 300 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 1., 1.04);

  TH1F *reco_phi_near_mass_pt_small_cut = new TH1F(
      "reco_phi P_{T} daughter < 200 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 1., 1.04);

  TH1F *reco_phi_near_mass_phi_pt_big_cut = new TH1F(
      "reco_phi pair P_{T} < 500 MeV",
      "Run 19 A+A Reco Phi ;KK Invariant Mass(GeV/c^2);counts", 100, 1., 1.04);

  TH1F *reco_phi_near_mass_phi_pt_cut = new TH1F(
      "reco_phi pair P_{T} < 300 MeV",
      "Run 19 A+A Reco Phi ;KK Invariant Mass(GeV/c^2);counts", 100, 1., 1.04);

  TH1F *reco_phi_near_mass_phi_pt_small_cut = new TH1F(
      "reco_phi pair P_{T} < 200 MeV",
      "Run 19 A+A Reco Phi ;KK Invariant Mass(GeV/c^2);counts", 100, 1., 1.04);

  TH2F *reco_phi_mass_pt =
      new TH2F("reco_phi",
               "Run 19 A+A Invariant Mass Reco Phi and P_{T} Reco Phi;KK "
               "Invariant Mass(GeV/c^2);P_{T}",
               40, 0.95, 1.4, 40, 0, 0.9);

  TH1F *reco_phi_pt =
      new TH1F("reco_phi ", "Run 19 A+A P_{T} Reco Phi ;P_{T}(GeV/c);counts",
               100, 0., 0.9);

  TH1F *reco_phi_pt_mass_cut =
      new TH1F("reco_phi #phi mass #in [1.0, 1.04] GeV/c^{2}",
               "Run 19 A+A P_{T} Reco Phi ;P_{T}(GeV/c);counts", 100, 0, 0.9);

  TH1F *reco_phi_eta = new TH1F(
      "reco_phi ", "Run 19 A+A #eta Reco Phi ;#eta;counts", 100, -4., 4.);

  TH1F *reco_phi_phi =
      new TH1F("reco_phi ", "Run 19 A+A #phi Reco Phi ;#phi(rad);counts", 100,
               -3.15, 3.15);

  // histogram for same charge
  TH1F *KpKp_mass = new TH1F(
      "KpKp ", "Run 19 A+A Invariant Mass K+K+ ;Invariant Mass(GeV/c^2);counts",
      100, 0.9, 1.5);

  TH1F *KpKp_pt = new TH1F(
      "KpKp ", "Run 19 A+A P_{T} K+K+ ;P_{T}(GeV/c);counts", 100, 0., 1.6);

  TH1F *KpKp_eta =
      new TH1F("KpKp ", "Run 19 A+A #eta K+K+ ;#eta;counts", 100, -4., 4.);

  TH1F *KpKp_phi = new TH1F("KpKp ", "Run 19 A+A #phi K+K+ ;#phi(rad);counts",
                            100, -3.15, 3.15);

  TH1F *KmKm_mass = new TH1F(
      "KmKm ", "Run 19 A+A Invariant Mass K-K- ;Invariant Mass(GeV/c^2);counts",
      100, 0.9, 1.8);

  TH1F *KmKm_pt = new TH1F(
      "KmKm ", "Run 19 A+A P_{T} K-K- ;P_e z{T}(GeV/c);counts", 100, 0., 1.6);

  TH1F *KmKm_eta =
      new TH1F("KmKm ", "Run 19 A+A #eta K-K- ;#eta;counts", 100, -4., 4.);

  TH1F *KmKm_phi = new TH1F("KmKm ", "Run 19 A+A #phi K-K- ;#phi(rad);counts",
                            100, -3.15, 3.15);

  TH1F *like_sign_mass = new TH1F("like_sign ",
                                  "Run 19 A+A Invariant Mass Like Sign Kaons "
                                  ";Invariant Mass(GeV/c^2);counts",
                                  100, 0.9, 1.8);

  TH1F *like_sign_pt = new TH1F(
      "like_sign ", "Run 19 A+A P_{T} Like Sign Kaons ;P_{T}(GeV/c);counts",
      100, 0., 1.6);

  TH1F *like_sign_eta =
      new TH1F("like_sign ", "Run 19 A+A #eta Like Sign Kaons ;#eta;counts",
               100, -4., 4.);

  TH1F *like_sign_phi = new TH1F(
      "like_sign ", "Run 19 A+A #phi K-K- ;#phi(rad);counts", 100, -3.15, 3.15);
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
    lv = lv1 + lv2;

    // all track histogram filling
    all_pt->Fill(pair->d1_mPt);
    all_pt->Fill(pair->d2_mPt);

    all_DCA->Fill(pair->d1_mDCA);
    all_DCA->Fill(pair->d2_mDCA);

    all_eta->Fill(pair->d1_mEta);
    all_eta->Fill(pair->d2_mEta);

    all_phi->Fill(pair->d1_mPhi);
    all_phi->Fill(pair->d2_mPhi);

    all_Nkaon->Fill(pair->d1_mNSigmaKaon);
    all_Nkaon->Fill(pair->d2_mNSigmaKaon);

    all_Npion->Fill(pair->d1_mNSigmaPion);
    all_Npion->Fill(pair->d2_mNSigmaPion);

    all_chisq_pipi->Fill(pair->d1_mNSigmaPion * pair->d1_mNSigmaPion + pair->d2_mNSigmaPion * pair->d2_mNSigmaPion);
    all_chisq_kk->Fill(pair->d1_mNSigmaKaon * pair->d1_mNSigmaKaon + pair->d2_mNSigmaKaon * pair->d2_mNSigmaKaon);
    all_chisq_ee->Fill(pair->d1_mNSigmaElectron * pair->d1_mNSigmaElectron + pair->d2_mNSigmaElectron * pair->d2_mNSigmaElectron);
    all_chisq_pp->Fill(pair->d1_mNSigmaProton * pair->d1_mNSigmaProton + pair->d2_mNSigmaProton * pair->d2_mNSigmaProton);

    all_pt_NPion->Fill(pair->d1_mPt, pair->d1_mNSigmaPion);
    all_pt_NPion->Fill(pair->d2_mPt, pair->d2_mNSigmaPion);
    all_pt_NKaon->Fill(pair->d1_mPt, pair->d1_mNSigmaKaon);
    all_pt_NKaon->Fill(pair->d1_mPt, pair->d1_mNSigmaKaon);
    all_pt_NElectron->Fill(pair->d2_mPt, pair->d2_mNSigmaElectron);
    all_pt_NElectron->Fill(pair->d2_mPt, pair->d2_mNSigmaElectron);
    all_pt_NProton->Fill(pair->d2_mPt, pair->d2_mNSigmaProton);
    all_pt_NProton->Fill(pair->d2_mPt, pair->d2_mNSigmaProton);


    all_nsigmapion_nsigmakaon->Fill(pair->d1_mNSigmaPion, pair->d1_mNSigmaKaon);
    all_nsigmapion_nsigmakaon->Fill(pair->d2_mNSigmaPion, pair->d2_mNSigmaKaon);

    all_pm_ratio->Fill(pair->mChargeSum);

    // kaon candidate filling with PId selection
    if ((pair->d1_mNSigmaPion > 5.8 || pair->d1_mNSigmaPion < -4.2) &&
        (pair->d2_mNSigmaPion > 5.8 || pair->d2_mNSigmaPion < -4.2) &&
        (pair->d1_mNSigmaElectron > 5.424112 ||
         pair->d1_mNSigmaElectron < -4.575888) &&
        (pair->d2_mNSigmaElectron > 5.424112 ||
         pair->d2_mNSigmaElectron < -4.575888) &&
        (pair->d1_mNSigmaProton > 7.36 || pair->d1_mNSigmaProton < -2.64) &&
        (pair->d2_mNSigmaProton > 7.36 || pair->d2_mNSigmaProton < -2.64) &&
        (pair->d1_mNSigmaKaon < 7.4 && pair->d1_mNSigmaKaon > -2.6) &&
        (pair->d2_mNSigmaKaon < 7.4 && pair->d2_mNSigmaKaon > -2.6)) {

      kaon_Nkaon->Fill(pair->d1_mNSigmaKaon);
      kaon_Nkaon->Fill(pair->d2_mNSigmaKaon);

      kaon_Npion->Fill(pair->d1_mNSigmaPion);
      kaon_Npion->Fill(pair->d2_mNSigmaPion);

      // reconstrauction of parent lorentz vector
      // assuming kaon mass

      // reco level filling
      if (pair->mChargeSum == 0) {
        kaon_pm_ratio->Fill(pair->mChargeSum);
        reco_phi_pt->Fill(lv.Pt());
        reco_phi_mass->Fill(lv.M());

        kaon_pt->Fill(pair->d1_mPt);
        kaon_pt->Fill(pair->d2_mPt);

        kaon_DCA->Fill(pair->d1_mDCA);
        kaon_DCA->Fill(pair->d2_mDCA);

        kaon_eta->Fill(pair->d1_mEta);
        kaon_eta->Fill(pair->d2_mEta);

        kaon_phi->Fill(pair->d1_mPhi);
        kaon_phi->Fill(pair->d2_mPhi);

        if (lv.M() >= 1. && lv.M() <= 1.04) {
          reco_phi_near_mass->Fill(lv.M());
          if (pair->d1_mPt < 0.4 && pair->d2_mPt < 0.4) {
            reco_phi_near_mass_pt_big_cut->Fill(lv.M());
          }
          if (pair->d1_mPt < 0.3 && pair->d2_mPt < 0.3) {
            reco_phi_near_mass_pt_cut->Fill(lv.M());
          }
          if (pair->d1_mPt < 0.2 && pair->d2_mPt < 0.2) {
            reco_phi_near_mass_pt_small_cut->Fill(lv.M());
          }
          if (lv.Pt() < 0.5) {
            reco_phi_near_mass_phi_pt_big_cut->Fill(lv.M());
          }
          if (lv.Pt() < 0.3) {
            reco_phi_near_mass_phi_pt_cut->Fill(lv.M());
          }
          if (lv.Pt() < 0.2) {
            reco_phi_near_mass_phi_pt_small_cut->Fill(lv.M());
          }
          if (lv.Pt() < 0.1){
            reco_coherent_phi_mass->Fill(lv.M());
            if (pair->d1_mPt > 0.2 && pair->d2_mPt > 0.2){
              reco_coherent_phi_mass_tof->Fill(lv.M());
            }
          }
        }
        reco_phi_eta->Fill(lv.Eta());
        reco_phi_phi->Fill(lv.Phi());
        reco_phi_mass_pt->Fill(lv.M(), lv.Pt());
        if (pair->d1_mPt < 0.4 && pair->d2_mPt < 0.4) {
          reco_phi_mass_pt_big_cut->Fill(lv.M());
        }
        if (pair->d1_mPt < 0.3 && pair->d2_mPt < 0.3) {
          reco_phi_mass_pt_cut->Fill(lv.M());
        }
        if (pair->d1_mPt < 0.2 && pair->d2_mPt < 0.2) {
          reco_phi_mass_pt_small_cut->Fill(lv.M());
        }
        if (lv.M() < 1.04 && lv.M() > 1.00) {
          reco_phi_pt_mass_cut->Fill(lv.Pt());
        }
        if (lv.Pt() < 0.5) {
          reco_phi_mass_phi_pt_big_cut->Fill(lv.M());
        }
        if (lv.Pt() < 0.3) {
          reco_phi_mass_phi_pt_cut->Fill(lv.M());
        }
        if (lv.Pt() < 0.2) {
          reco_phi_mass_phi_pt_small_cut->Fill(lv.M());
        }
      } else if (pair->mChargeSum == 2) {
        kaon_pm_ratio->Fill(pair->mChargeSum);
        KpKp_pt->Fill(lv.Pt());
        KpKp_mass->Fill(lv.M());
        KpKp_eta->Fill(lv.Eta());
        KpKp_phi->Fill(lv.Phi());

        like_sign_pt->Fill(lv.Pt());
        like_sign_mass->Fill(lv.M());
        like_sign_eta->Fill(lv.Eta());
        like_sign_phi->Fill(lv.Phi());

      } else if (pair->mChargeSum == -2) {
        kaon_pm_ratio->Fill(pair->mChargeSum);
        KmKm_pt->Fill(lv.Pt());
        KmKm_mass->Fill(lv.M());
        KmKm_eta->Fill(lv.Eta());
        KmKm_phi->Fill(lv.Phi());

        like_sign_pt->Fill(lv.Pt());
        like_sign_mass->Fill(lv.M());
        like_sign_eta->Fill(lv.Eta());
        like_sign_phi->Fill(lv.Phi());
      }
    }

  } // loop on events

  makeCan();
  kaon_pm_ratio->Draw();
  gPad->Print("./Plots_sigma/Kaons/kaon_charge.png");

  makeCan();
  all_DCA->Draw();
  gPad->Print("./Plots_sigma/All/all_DCA.png");
  makeCan();
  all_pt->Draw();
  gPad->Print("./Plots_sigma/All/all_pt.png");
  makeCan();
  all_eta->Draw();
  gPad->Print("./Plots_sigma/All/all_eta.png");
  makeCan();
  all_phi->Draw();
  gPad->Print("./Plots_sigma/All/all_phi.png");
  makeCan();
  all_Nkaon->Draw();
  gPad->Print("./Plots_sigma/All/all_NSigmaKaon.png");
  makeCan();
  all_Npion->Draw();
  gPad->Print("./Plots_sigma/All/all_NSigmaPion.png");
  makeCan();
  all_chisq_pipi->Draw();
  gPad->Print("./Plots_sigma/All/all_chisq_pipi.png");
  makeCan();
  all_chisq_kk->Draw();
  gPad->Print("./Plots_sigma/All/all_chisq_kk.png");
  makeCan();
  all_chisq_ee->Draw();
  gPad->Print("./Plots_sigma/All/all_chisq_ee.png");
  makeCan();
  all_chisq_pp->Draw();
  gPad->Print("./Plots_sigma/All/all_pchisq_pp.png");

  all_pt_NPion->Draw("colz");
  gPad->Print("./Plots_sigma/All/all_pt_NSigmaPion.png");
  makeCan();
  all_pt_NKaon->Draw("colz");
  gPad->Print("./Plots_sigma/All/all_pt_NSigmaKaon.png");
  makeCan();
  all_pt_NProton->Draw("colz");
  gPad->Print("./Plots_sigma/All/all_pt_NSigmaProton.png");
  makeCan();
  all_pt_NElectron->Draw("colz");
  gPad->Print("./Plots_sigma/All/all_pt_NSigmaElectron.png");
  makeCan();
  all_nsigmapion_nsigmakaon->Draw("colz");
  gPad->Print("./Plots_sigma/All/all_NSigmaPion_NSigmaKaon.png");
  makeCan();
  all_pm_ratio->Draw();
  gPad->Print("./Plots_sigma/All/all_pm_ratio.png");

  makeCan();
  kaon_DCA->Draw();
  gPad->Print("./Plots_sigma/Kaons/kaon_DCA.png");
  makeCan();
  kaon_pt->Draw();
  gPad->Print("./Plots_sigma/Kaons/kaon_pt.png");
  makeCan();
  kaon_eta->Draw();
  gPad->Print("./Plots_sigma/Kaons/kaon_eta.png");
  makeCan();
  kaon_phi->Draw();
  gPad->Print("./Plots_sigma/Kaons/kaon_phi.png");
  makeCan();
  kaon_Nkaon->Draw();
  gPad->Print("./Plots_sigma/Kaons/kaon_NSigmaKaon.png");
  makeCan();
  kaon_Npion->Draw();
  gPad->Print("./Plots_sigma/Kaons/kaon_NSigmaPion.png");

  makeCan();
  reco_phi_pt->Draw();
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_phi_pt.png");
  makeCan();
  reco_phi_eta->Draw();
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_phi_eta.png");
  makeCan();
  reco_phi_phi->Draw();
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_phi_phi.png");
  makeCan();
  TLegend *legend_1 = new TLegend(0.1, 0.7, 0.3, 0.9);
  reco_phi_mass->Draw("pfc");
  reco_phi_mass_pt_big_cut->Draw("pfc;same");
  reco_phi_mass_pt_cut->Draw("pfc;same");
  reco_phi_mass_pt_small_cut->Draw("pfc;same");
  legend_1->AddEntry(reco_phi_mass, "All P_{T}");
  legend_1->AddEntry(reco_phi_mass_pt_big_cut, "Kaon P_{T} < 400 MeV");
  legend_1->AddEntry(reco_phi_mass_pt_cut, "Kaon P_{T} < 300 MeV");
  legend_1->AddEntry(reco_phi_mass_pt_small_cut, "Kaon P_{T} < 200 MeV (TOF)");
  legend_1->Draw("same");
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_phi_mass.png");
  makeCan();
  TLegend *legend_3 = new TLegend(0.1, 0.7, 0.3, 0.9);
  reco_phi_mass->Draw("pfc");
  reco_phi_mass_phi_pt_big_cut->Draw("pfc;same");
  reco_phi_mass_phi_pt_cut->Draw("pfc;same");
  reco_phi_mass_phi_pt_small_cut->Draw("pfc;same");
  legend_3->AddEntry(reco_phi_mass, "All P_{T}");
  legend_3->AddEntry(reco_phi_mass_phi_pt_big_cut, "Pair P_{T} < 500 MeV");
  legend_3->AddEntry(reco_phi_mass_phi_pt_cut, "Pair P_{T} < 300 MeV");
  legend_3->AddEntry(reco_phi_mass_phi_pt_small_cut,
                     "Pair P_{T} < 200 MeV (TOF)");
  legend_3->Draw("same");
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_phi_mass_1.png");
  makeCan();
  TLegend *legend_2 = new TLegend(0.1, 0.7, 0.3, 0.9);
  reco_phi_near_mass->Draw("pfc");
  reco_phi_near_mass_pt_big_cut->Draw("pfc;same");
  reco_phi_near_mass_pt_cut->Draw("pfc;same");
  reco_phi_near_mass_pt_small_cut->Draw("pfc;same");
  legend_2->AddEntry(reco_phi_near_mass, "All P_{T}");
  legend_2->AddEntry(reco_phi_near_mass_pt_big_cut, "Kaon P_{T} < 400 MeV");
  legend_2->AddEntry(reco_phi_near_mass_pt_cut, "Kaon P_{T} < 300 MeV");
  legend_2->AddEntry(reco_phi_near_mass_pt_small_cut,
                     "Kaon P_{T} < 200 MeV (TOF)");
  legend_2->Draw("same");
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_phi_near_mass.png");
  makeCan();
  TLegend *legend_4 = new TLegend(0.1, 0.7, 0.3, 0.9);
  reco_phi_near_mass->Draw("pfc");
  reco_phi_near_mass_phi_pt_big_cut->Draw("pfc;same");
  reco_phi_near_mass_phi_pt_cut->Draw("pfc;same");
  reco_phi_near_mass_phi_pt_small_cut->Draw("pfc;same");
  legend_4->AddEntry(reco_phi_near_mass, "All P_{T}");
  legend_4->AddEntry(reco_phi_near_mass_phi_pt_big_cut, "Pair P_{T} < 500 MeV");
  legend_4->AddEntry(reco_phi_near_mass_phi_pt_cut, "Pair P_{T} < 300 MeV");
  legend_4->AddEntry(reco_phi_near_mass_phi_pt_small_cut,
                     "Pair P_{T} < 200 MeV (TOF)");
  legend_4->Draw("same");
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_phi_near_mass_1.png");
  makeCan();
  reco_phi_mass_pt->Draw("colz");
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_phi_mass_pt.png");
  makeCan();
  reco_phi_pt_mass_cut->Draw();
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_phi_pt_mass_cut.png");

  makeCan();
  TLegend *legend_5 = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend_5->AddEntry(reco_coherent_phi_mass, "Coherent #phi Candidates");
  legend_5->AddEntry(reco_coherent_phi_mass_tof, "Coherent #phi Candidates TOF");
  reco_coherent_phi_mass->Draw("pfc");
  reco_coherent_phi_mass_tof->Draw("same;pfc");
  legend_5->Draw("same");
  gPad->Print("./Plots_sigma/Kaons/Reco/reco_coherent_phi_mass.png");


  makeCan();
  KpKp_pt->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/kpkp_pt.png");
  makeCan();
  KpKp_eta->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/kpkp_eta.png");
  makeCan();
  KpKp_phi->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/kpkp_phi.png");
  makeCan();
  KpKp_mass->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/kpkp_mass.png");
  makeCan();
  KmKm_pt->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/kmkm_pt.png");
  makeCan();
  KmKm_eta->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/kmkm_eta.png");
  makeCan();
  KmKm_phi->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/kmkm_phi.png");
  makeCan();
  KmKm_mass->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/kmkm_mass.png");

  makeCan();
  like_sign_pt->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/ls_pt.png");
  makeCan();
  like_sign_mass->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/ls_mass.png");
  makeCan();
  like_sign_eta->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/ls_eta.png");
  makeCan();
  like_sign_phi->Draw();
  gPad->Print("./Plots_sigma/Kaons/like_signs/ls_phi.png");
}
