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
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TH2.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <cmath>
#include "TPaveStats.h"
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

void daughter_distribution_kaon_debug() {

  // histograms for all tracks
  TH1F *all_DCA = new TH1F("All Tracks", "Run 19 A+A DCA All Tracks;DCA(cm);counts",
                           100, 0., 3.);

  TH1F *all_pt =
      new TH1F("All Tracks", "Run 19 A+A P_{T} All Tracks;P_{T}(GeV/c);counts",
               100, 0., 1.8);
  TH1F *all_eta = new TH1F(
      "All Tracks", "Run 19 A+A #eta All Tracks;#eta;counts", 100, -2.2, 2.2);
  TH1F *all_phi =
      new TH1F("All Tracks", "Run 19 A+A #phi All Tracks;#phi(rad);counts", 100,
               -3.15, 3.15);
  TH1F * all_Npion = new TH1F(
      "All Tracks", "Run 19 A+A N#sigma#pi All Tracks;(N#sigmaPion);counts",
      100, -10, 80);

  TH1F *all_Nkaon = new TH1F(
      "All Tracks", "Run 19 A+A N#sigmaK All Tracks;(N#sigmaKaon);counts", 100,
      -50, 80);

  TH1F *all_chisq_pipi = new TH1F(
      "All Tracks",
      "Run 19 A+A #chi^{2}_{#pi#pi} All Tracks;#chi^{2}_{#pi#pi};counts", 100,
      -10, 100);

  TH1F *all_chisq_kk = new TH1F(
      "All Tracks", "Run 19 A+A #chi^{2}_{KK} All Tracks;#chi^{2}_{KK};counts",
      100, -10, 4000);

  TH1F *all_chisq_ee = new TH1F(
      "All Tracks", "Run 19 A+A #chi^{2}_{ee} All Tracks;#chi^{2}_{ee};counts",
      100, -10, 800);

  TH1F *all_chisq_pp = new TH1F(
      "All Tracks", "Run 19 A+A #chi^{2}_{pp} All Tracks;#chi^{2}_{pp};counts",
      100, -10, 400);

  TH2F *all_pt_NPion =
      new TH2F("All Tracks",
               "Run 19 A+A N#sigmaPion All Tracks;P_{T}(GeV/c);N#sigmaPion",
               100, 0., 1.2, 200, -20, 80);

  TH2F *all_pt_NKaon =
      new TH2F("All Tracks",
               "Run 19 A+A N#sigmaKaon All Tracks;P_{T}(GeV/c);N#sigmaKaon",
               100, 0., 1.2, 200, -60, 80);

  TH2F *all_pt_NElectron = new TH2F(
      "All Tracks",
      "Run 19 A+A N#sigmaElectron All Tracks;P_{T}(GeV/c);N#sigmaElectron", 100,
      0., 1.2, 200, -20, 80);

  TH2F *all_pt_NProton =
      new TH2F("All Tracks",
               "Run 19 A+A N#sigmaProton All Tracks;P_{T}(GeV/c);N#sigmaProton",
               100, 0., 1.2, 200, -80, 60);

  TH2F *all_chisqPiPi_chisqKK =
      new TH2F("All Tracks",
               "Run 19 A + A #chi^{2}_{K^{+}K^{-}} v.s. #chi^{2}_{#pi^{+}#pi^{-}} All "
               "Tracks;#chi^{2}_{#pi^{+}#pi^{-}};#chi^{2}_{K^{+}K^{-}}",
               100, 0, 180, 100, 0, 60);

  TH2F *all_chisqEE_chisqKK =
      new TH2F("All Tracks",
               "Run 19 Au + Au #chi^{2}_{K^{+}K^{-}} v.s. #chi^{2}_{e#bar{e}} All "
               "Tracks;#chi^{2}_{e#bar{e}};#chi^{2}_{K^{+}K^{-}}",
               100, 0, 180, 100, 0, 60);

  TH2F *all_chisqPP_chisqKK =
      new TH2F("All Tracks",
               "Run 19 Au + Au #chi^{2}_{K^{+}K^{-}} v.s. #chi^{2}_{p#bar{p}} All "
               "Tracks;#chi^{2}_{p#bar{p}};#chi^{2}_{K^{+}K^{-}}",
               100, 0, 180, 100, 0, 60);

  TH1F *all_pm_ratio = new TH1F(
      "All Trakcs", "Run 19 A+A Two-Track Charge Sum;Charge Sum;counts", 27,
      -2.2, 2.2);

  // histogram for kaon canidataes
  TH1F *kaon_DCA =
      new TH1F("Kaon Candidates", "Run 19 A+A DCA Kaon Candidates;DCA(cm);counts",
               100, 0., 3.);

  TH1F *kaon_pt = new TH1F(
      "Kaon Candidates", "Run 19 A+A P_{T} Kaon Candidates;P_{T}(GeV/c);counts",
      100, 0.05, 1.0);
  TH1F *kaon_eta =
      new TH1F("Kaon Candidates", "Run 19 A+A #eta Kaon Candidates;#eta;counts",
               100, -2., 2.);
  TH1F *kaon_phi = new TH1F("Kaon Candidates",
                            "Run 19 A+A #phi Kaon Candidates;#phi(rad);counts",
                            50, -3.15, 3.15);
  TH1F *kaon_Npion = new TH1F(
      "Kaon Candidates",
      "Run 19 A+A N#sigmaPion Kaon Candidates;N#sigmaPion;counts", 100, 5, 30);

  TH1F *kaon_Nkaon = new TH1F(
      "Kaon Cnadidates",
      "Run 19 A+A Corrected N#sigmaKaon Kaon Candidates;Corrected N#sigmaKaon;counts", 100, -6, 6);

  TH1F *kaon_pm_ratio = new TH1F(
      "Kaon Candidates",
      "Run 19 A+A Two-Track Charge Sum Kaon Candidates;Charge Sum;counts", 27, -2.2, 2.2);

  TH1F *kaon_coherent_eta =
      new TH1F("Kaon Candidates for Coherent #phi",
               "Run 19 A+A #eta Kaon Candidates;#eta;counts", 100, -2., 2.);

  // histogram for reco phi
  TH1F *reco_phi_mass = new TH1F(
      "reco_phi ",
      "Run 19 A+A Invariant Mass Reconstructed #phi Meson ;K^{+}K^{-} Pair Invariant Mass(GeV/c^{2});counts",
      100, 0.95, 1.35);

  TH1F *reco_coherent_phi_mass =
      new TH1F("reco_phi ",
               "Run 19 A+A Invariant Mass Reconstructed Coherent #phi Candidates;K^{+}K^{-} Pair Invariant "
               "Mass(GeV/c^{2});counts",
               20, 1.0, 1.04);

  TH1F *reco_coherent_phi_mass_tof =
      new TH1F("reco_phi ",
               "Run 19 A+A Invariant Mass Reco Phi Cohehrent;KK Invariant "
               "Mass(GeV/c^{2});counts",
               100, 1.0, 1.04);

  TH1F *reco_phi_mass_pt_big_cut = new TH1F(
      "reco_phi P_{T} daughter < 400 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.95, 1.35);

  TH1F *reco_phi_coherent_eta = new TH1F(
      "Coherent #phi", "Run 19 A+A #eta Coherent #phi Candidates;#eta;counts",
      100, -4., 4.);

  TH1F *reco_phi_mass_pt_cut = new TH1F(
      "reco_phi P_{T} daughter < 300 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.95, 1.35);

  TH1F *reco_phi_mass_pt_small_cut = new TH1F(
      "reco_phi P_{T} daughter < 200 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.95, 1.35);

  TH1F *reco_phi_mass_phi_pt_big_cut = new TH1F(
      "reco_phi P_{T} Pair P_{T} < 500 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.95, 1.35);

  TH1F *reco_phi_mass_phi_pt_cut = new TH1F(
      "reco_phi P_{T} Pair P_{T} < 300 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.95, 1.35);

  TH1F *reco_phi_mass_phi_pt_small_cut = new TH1F(
      "reco_phi P_{T} Pair P_{T} < 200 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      100, 0.95, 1.35);

  TH1F *reco_phi_near_mass = new TH1F(
      "reco_phi",
      "Run 19 A+A Invariant Mass Reconstructed #phi Meson ;K^{+}K^{-} Pair Invariant Mass(GeV/c^{2});counts",
      20, 0.98, 1.1);

  TH1F *reco_phi_near_mass_pt_big_cut = new TH1F(
      "reco_phi P_{T} daughter < 400 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      20, 1., 1.04);

  TH1F *reco_phi_near_mass_pt_cut = new TH1F(
      "reco_phi P_{T} daughter < 300 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      20, 1., 1.04);

  TH1F *reco_phi_near_mass_pt_small_cut = new TH1F(
      "reco_phi P_{T} daughter < 100 MeV",
      "Run 19 A+A Invariant Mass Reco Phi ;KK Invariant Mass(GeV/c^2);counts",
      20, 1., 1.04);

  TH1F *reco_phi_near_mass_phi_pt_big_cut = new TH1F(
      "reco_phi pair P_{T} < 500 MeV",
      "Run 19 A+A Reco Phi ;KK Invariant Mass(GeV/c^2);counts", 20, 0.98, 1.1);

  TH1F *reco_phi_near_mass_phi_pt_cut = new TH1F(
      "reco_phi pair P_{T} < 300 MeV",
      "Run 19 A+A Reco Phi ;KK Invariant Mass(GeV/c^2);counts", 20, 0.98, 1.1);

  TH1F *reco_phi_near_mass_phi_pt_small_cut = new TH1F(
      "reco_phi pair P_{T} < 100 MeV",
      "Run 19 A+A Reco Phi ;KK Invariant Mass(GeV/c^2);counts", 20, 0.98, 1.1);

  TH2F *reco_phi_mass_pt =
      new TH2F("reco_phi",
               "Run 19 A+A K^{+}K^{-} Pair P_{T} v.s. Invariant Mass Reconstructed #phi Meson;K^{+}K^{-}"
               "Invariant Mass(GeV/c^{2});K^{+}K^{-} Pair P_{T} (GeV/c)",
               20, 0.98, 1.35, 20, 0, 0.9);

  TH1F *reco_phi_pt =
      new TH1F("reco_phi ", "Run 19 A+A P_{T} Reco Phi ;P_{T}(GeV/c);counts",
               100, 0., 0.9);

  TH1F *reco_phi_pt_mass_cut =
      new TH1F("reco_phi #phi mass #in [1.0, 1.04] GeV/c^{2}",
               "Run 19 A+A P_{T} Reconstructed #phi Meson;K^{+}K^{-} Pair P_{T}(GeV/c);counts", 100, 0, 0.9);

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

  TH1F *like_sign_near_mass =
      new TH1F("Like Sign Kaon Candidates",
               "Run 19 A+A Invariant Mass Like Sign Kaons"
               ";Invariant Mass(GeV/c^{2});counts",
               100, 1., 1.04);

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

    all_chisq_pipi->Fill(pair->d1_mNSigmaPion * pair->d1_mNSigmaPion +
                         pair->d2_mNSigmaPion * pair->d2_mNSigmaPion);
    all_chisq_kk->Fill(pair->d1_mNSigmaKaon * pair->d1_mNSigmaKaon +
                       pair->d2_mNSigmaKaon * pair->d2_mNSigmaKaon);
    all_chisq_ee->Fill(pair->d1_mNSigmaElectron * pair->d1_mNSigmaElectron +
                       pair->d2_mNSigmaElectron * pair->d2_mNSigmaElectron);
    all_chisq_pp->Fill(pair->d1_mNSigmaProton * pair->d1_mNSigmaProton +
                       pair->d2_mNSigmaProton * pair->d2_mNSigmaProton);

    all_pt_NPion->Fill(pair->d1_mPt, pair->d1_mNSigmaPion);
    all_pt_NPion->Fill(pair->d2_mPt, pair->d2_mNSigmaPion);
    all_pt_NKaon->Fill(pair->d1_mPt, pair->d1_mNSigmaKaon);
    all_pt_NKaon->Fill(pair->d1_mPt, pair->d1_mNSigmaKaon);
    all_pt_NElectron->Fill(pair->d2_mPt, pair->d2_mNSigmaElectron);
    all_pt_NElectron->Fill(pair->d2_mPt, pair->d2_mNSigmaElectron);
    all_pt_NProton->Fill(pair->d2_mPt, pair->d2_mNSigmaProton);
    all_pt_NProton->Fill(pair->d2_mPt, pair->d2_mNSigmaProton);

    double chisq_ee = (pair->d1_mNSigmaElectron - NSigmaElectronShift) *
                          (pair->d1_mNSigmaElectron - NSigmaElectronShift) +
                      (pair->d2_mNSigmaElectron - NSigmaElectronShift) *
                          (pair->d2_mNSigmaElectron - NSigmaElectronShift);

    double chisq_pipi = (pair->d1_mNSigmaPion - NSigmaPionShift) *
                            (pair->d1_mNSigmaPion - NSigmaPionShift) +
                        (pair->d2_mNSigmaPion - NSigmaPionShift) *
                            (pair->d2_mNSigmaPion - NSigmaPionShift);

    double chisq_kk = (pair->d1_mNSigmaKaon - NSigmaKaonShift) *
                          (pair->d1_mNSigmaKaon - NSigmaKaonShift) +
                      (pair->d2_mNSigmaKaon - NSigmaKaonShift) *
                          (pair->d2_mNSigmaKaon - NSigmaKaonShift);

    double chisq_pp = (pair->d1_mNSigmaProton - NSigmaProtonShift) *
                          (pair->d1_mNSigmaProton - NSigmaProtonShift) +
                      (pair->d2_mNSigmaProton - NSigmaProtonShift) *
                          (pair->d2_mNSigmaProton - NSigmaProtonShift);

    bool electron_rejection = chisq_ee > 20;
    bool pion_rejection = chisq_pipi > 50;
    bool proton_rejection = chisq_pp > 50 && chisq_pp > (1.5 * chisq_kk + 50);
    bool kaon_acceptance = chisq_kk < 30;

    //
    // USE FOR DEBUG ONLY
    // bool pion_rejection = true;
    // bool electron_rejection = true;
    // bool proton_rejection = true;


    all_chisqPiPi_chisqKK->Fill(chisq_pipi, chisq_kk);

    // all_chisqPiPi_chisqKK->Fill(
    //     (pair->d1_mNSigmaPion ) *
    //             (pair->d1_mNSigmaPion ) +
    //         (pair->d2_mNSigmaPion ) *
    //             (pair->d2_mNSigmaPion ),
    //     (pair->d1_mNSigmaKaon ) *
    //             (pair->d1_mNSigmaKaon ) +
    //         (pair->d2_mNSigmaKaon ) *
    //             (pair->d2_mNSigmaKaon ));

    all_chisqEE_chisqKK->Fill(chisq_ee, chisq_kk);

    all_chisqPP_chisqKK->Fill(chisq_pp, chisq_kk);

    all_pm_ratio->Fill(pair->mChargeSum);

    // kaon candidate filling with PId selection
    if (abs(pair->d1_mNSigmaKaon - NSigmaKaonShift) < 5 &&
        abs(pair->d2_mNSigmaKaon - NSigmaKaonShift) < 5) {
    // if (electron_rejection && pion_rejection && proton_rejection && kaon_acceptance){


        kaon_Nkaon->Fill(pair->d1_mNSigmaKaon - NSigmaKaonShift);
        kaon_Nkaon->Fill(pair->d2_mNSigmaKaon - NSigmaKaonShift);

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
          if (lv.Pt() < 0.1) {
            reco_phi_near_mass_phi_pt_small_cut->Fill(lv.M());
          }
          if (lv.Pt() < 0.1) {
            reco_coherent_phi_mass->Fill(lv.M());
            kaon_coherent_eta->Fill(lv1.Eta());
            kaon_coherent_eta->Fill(lv2.Eta());
            reco_phi_coherent_eta->Fill(lv.Eta());
            if (pair->d1_mPt > 0.2 && pair->d2_mPt > 0.2) {
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

        if (lv.M() > 1. && lv.M() < 1.04) {
          like_sign_near_mass->Fill(lv.M());
        }

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

        if (lv.M() > 1. && lv.M() < 1.04) {
          like_sign_near_mass->Fill(lv.M());
        }
      }
    }

  } // loop on events

  makeCan();
  kaon_pm_ratio->Draw();
  kaon_pm_ratio->GetXaxis()->SetTitleSize(0.05);
  kaon_pm_ratio->GetXaxis()->CenterTitle();
  kaon_pm_ratio->GetXaxis()->SetTitleOffset(0.8);
  kaon_pm_ratio->SetLineWidth(2);

  makeCan();
  all_DCA->GetXaxis()->SetTitleSize(0.05);
  all_DCA->GetXaxis()->CenterTitle();
  all_DCA->GetXaxis()->SetTitleOffset(0.8);
  all_DCA->SetLineWidth(2);
  all_DCA->Draw();
  makeCan();
  all_pt->GetXaxis()->SetTitleSize(0.05);
  all_pt->GetXaxis()->CenterTitle();
  all_pt->GetXaxis()->SetTitleOffset(0.8);
  all_pt->SetLineWidth(2);
  all_pt->Draw();
  makeCan();
  all_eta->GetXaxis()->SetTitleSize(0.05);
  all_eta->GetXaxis()->CenterTitle();
  all_eta->GetXaxis()->SetTitleOffset(0.8);
  all_eta->SetLineWidth(2);
  all_eta->Draw();
  makeCan();
  all_phi->GetXaxis()->SetTitleSize(0.05);
  all_phi->GetXaxis()->CenterTitle();
  all_phi->GetXaxis()->SetTitleOffset(0.8);
  all_phi->SetLineWidth(2);
  all_phi->Draw();
  makeCan();
  all_Nkaon->Draw();
  makeCan();
  all_Npion->Draw();
  makeCan();
  all_chisq_pipi->Draw();
  makeCan();
  all_chisq_kk->Draw();
  makeCan();
  all_chisq_ee->Draw();
  makeCan();
  all_chisq_pp->Draw();
  //
  makeCan();
  all_pt_NPion->GetXaxis()->SetTitleSize(0.05);
  all_pt_NPion->GetXaxis()->CenterTitle();
  all_pt_NPion->GetXaxis()->SetTitleOffset(0.8);
  all_pt_NPion->GetYaxis()->SetTitleSize(0.05);
  all_pt_NPion->GetYaxis()->CenterTitle();
  all_pt_NPion->GetYaxis()->SetTitleOffset(0.8);
  all_pt_NPion->Draw("colz");
  makeCan();
  all_pt_NKaon->GetXaxis()->SetTitleSize(0.05);
  all_pt_NKaon->GetXaxis()->CenterTitle();
  all_pt_NKaon->GetXaxis()->SetTitleOffset(0.8);
  all_pt_NKaon->GetYaxis()->SetTitleSize(0.05);
  all_pt_NKaon->GetYaxis()->CenterTitle();
  all_pt_NKaon->GetYaxis()->SetTitleOffset(0.8);
  all_pt_NKaon->Draw("colz");
  makeCan();
  all_pt_NProton->GetXaxis()->SetTitleSize(0.05);
  all_pt_NProton->GetXaxis()->CenterTitle();
  all_pt_NProton->GetXaxis()->SetTitleOffset(0.8);
  all_pt_NProton->GetYaxis()->SetTitleSize(0.05);
  all_pt_NProton->GetYaxis()->CenterTitle();
  all_pt_NProton->GetYaxis()->SetTitleOffset(0.8);
  all_pt_NProton->Draw("colz");
  makeCan();
  all_pt_NElectron->GetXaxis()->SetTitleSize(0.05);
  all_pt_NElectron->GetXaxis()->CenterTitle();
  all_pt_NElectron->GetXaxis()->SetTitleOffset(0.8);
  all_pt_NElectron->GetYaxis()->SetTitleSize(0.05);
  all_pt_NElectron->GetYaxis()->CenterTitle();
  all_pt_NElectron->GetYaxis()->SetTitleOffset(0.8);
  all_pt_NElectron->Draw("colz");
  makeCan(true);
  all_chisqPiPi_chisqKK->SetContour(100);
  all_chisqPiPi_chisqKK->GetXaxis()->SetTitleSize(0.05);
  all_chisqPiPi_chisqKK->GetXaxis()->CenterTitle();
  all_chisqPiPi_chisqKK->GetXaxis()->SetTitleOffset(0.8);
  all_chisqPiPi_chisqKK->GetYaxis()->SetTitleSize(0.05);
  all_chisqPiPi_chisqKK->GetYaxis()->CenterTitle();
  all_chisqPiPi_chisqKK->GetYaxis()->SetTitleOffset(0.8);
  all_chisqPiPi_chisqKK->Draw("colz");
  // TLine *line = new TLine(150, 50, 140, 50);
  // line->SetLineColor(kRed); // Setting line color to red for visibility
  // line->SetLineWidth(2);    // Setting line width
  // line->Draw("same");
  // TF1 *f1 = new TF1("f1", "(1./3.)*x", 0, 140);
  // f1->SetLineColor(
  //     kRed); // Setting the function's line color to blue for visibility
  // f1->Draw("same");

  makeCan();
  all_chisqEE_chisqKK->SetContour(100);
  all_chisqEE_chisqKK->GetXaxis()->SetTitleSize(0.05);
  all_chisqEE_chisqKK->GetXaxis()->CenterTitle();
  all_chisqEE_chisqKK->GetXaxis()->SetTitleOffset(0.8);
  all_chisqEE_chisqKK->GetYaxis()->SetTitleSize(0.05);
  all_chisqEE_chisqKK->GetYaxis()->CenterTitle();
  all_chisqEE_chisqKK->GetYaxis()->SetTitleOffset(0.8);
  all_chisqEE_chisqKK->Draw("colz");
  // TLine *line5 = new TLine(20, 0, 20, 30);
  // TLine *line6 = new TLine(20, 30, 180, 30);
  // line5->SetLineColor(kRed); // Setting line color to red for visibility
  // line5->SetLineWidth(2);    // Setting line width
  // line5->Draw("same");
  // line6->SetLineColor(kRed); // Setting line color to red for visibility
  // line6->SetLineWidth(2);    // Setting line width
  // line6->Draw("same");

  makeCan();
  TLine *line3 = new TLine(20, 0, 20, 3.5);
  // line3->SetLineColor(kRed); // Setting line color to red for visibility
  // line3->SetLineWidth(2);    // Setting line width
  all_chisqPP_chisqKK->SetContour(100);
  all_chisqPP_chisqKK->GetXaxis()->SetTitleSize(0.05);
  all_chisqPP_chisqKK->GetXaxis()->CenterTitle();
  all_chisqPP_chisqKK->GetXaxis()->SetTitleOffset(0.8);
  all_chisqPP_chisqKK->GetYaxis()->SetTitleSize(0.05);
  all_chisqPP_chisqKK->GetYaxis()->CenterTitle();
  all_chisqPP_chisqKK->GetYaxis()->SetTitleOffset(0.8);
  all_chisqPP_chisqKK->Draw("colz");
  // line3->Draw("same");
  // TLine *line4 = new TLine(90, 50, 140, 50);
  // line4->SetLineColor(kRed); // Setting line color to red for visibility
  // line4->SetLineWidth(2);    // Setting line width
  // line4->Draw("same");
  // TF1 *f2 = new TF1("f2", "(1./1.5)*(x-15)", 20, 90);
  // f2->SetLineColor(
  //     kRed); // Setting the function's line color to blue for visibility
  // f2->Draw("same");

  makeCan();
  all_pm_ratio->GetXaxis()->SetTitleSize(0.05);
  all_pm_ratio->GetXaxis()->CenterTitle();
  all_pm_ratio->GetXaxis()->SetTitleOffset(0.8);
  all_pm_ratio->SetLineWidth(2);
  all_pm_ratio->Draw();

  makeCan();
  kaon_DCA->GetXaxis()->SetTitleSize(0.05);
  kaon_DCA->GetXaxis()->CenterTitle();
  kaon_DCA->GetXaxis()->SetTitleOffset(0.8);
  kaon_DCA->SetLineWidth(2);
  kaon_DCA->SetFillColorAlpha(kBlue, 0.35);
  kaon_DCA->SetMarkerStyle(20); // This is a full circle marker
  kaon_DCA->SetMarkerColor(kRed); // Set point color to red
  kaon_DCA->SetMarkerSize(1);  // Adjust the size of the marker as needed
  kaon_DCA->Draw("e1");
  kaon_DCA->Draw("hist;same");
  makeCan();
  kaon_pt->GetXaxis()->SetTitleSize(0.05);
  kaon_pt->GetXaxis()->CenterTitle();
  kaon_pt->GetXaxis()->SetTitleOffset(0.8);
  kaon_pt->SetLineWidth(2);
  kaon_pt->SetFillColorAlpha(kBlue, 0.35);
  kaon_pt->SetMarkerStyle(20); // This is a full circle marker
  kaon_pt->SetMarkerColor(kRed); // Set point color to red
  kaon_pt->SetMarkerSize(1);  // Adjust the size of the marker as needed
  kaon_pt->Draw("e1");
  kaon_pt->Draw("hist;same");
  makeCan();
  kaon_eta->GetXaxis()->SetTitleSize(0.05);
  kaon_eta->GetXaxis()->CenterTitle();
  kaon_eta->GetXaxis()->SetTitleOffset(0.8);
  kaon_eta->SetLineWidth(2);
  kaon_eta->SetFillColorAlpha(kBlue, 0.35);
  kaon_eta->SetMarkerStyle(20); // This is a full circle marker
  kaon_eta->SetMarkerColor(kRed); // Set point color to red
  kaon_eta->SetMarkerSize(1);  // Adjust the size of the marker as needed
  kaon_eta->Draw("e1");
  kaon_eta->Draw("hist;same");
  makeCan();
  kaon_coherent_eta->Draw();
  makeCan();
  kaon_phi->GetXaxis()->SetTitleSize(0.05);
  kaon_phi->GetXaxis()->CenterTitle();
  kaon_phi->GetXaxis()->SetTitleOffset(0.8);
  kaon_phi->SetLineWidth(2);
  kaon_phi->SetFillColorAlpha(kBlue, 0.35);
  kaon_phi->SetMarkerStyle(20); // This is a full circle marker
  kaon_phi->SetMarkerColor(kRed); // Set point color to red
  kaon_phi->SetMarkerSize(1);  // Adjust the size of the marker as needed
  kaon_phi->Draw("e1");
  kaon_phi->Draw("hist;same");
  makeCan();
  kaon_Nkaon->GetXaxis()->SetTitleSize(0.05);
  kaon_Nkaon->GetXaxis()->CenterTitle();
  kaon_Nkaon->GetXaxis()->SetTitleOffset(0.8);
  kaon_Nkaon->SetLineWidth(2);
  kaon_Nkaon->SetFillColorAlpha(kBlue, 0.35);
  kaon_Nkaon->SetMarkerStyle(20); // This is a full circle marker
  kaon_Nkaon->SetMarkerColor(kRed); // Set point color to red
  kaon_Nkaon->SetMarkerSize(1);  // Adjust the size of the marker as needed
  kaon_Nkaon->Draw("e1");
  kaon_Nkaon->Draw("hist;same");
  makeCan();
  kaon_Npion->Draw();
  makeCan();
  reco_phi_pt->Draw();
  makeCan();
  reco_phi_eta->Draw();
  makeCan();
  reco_phi_phi->Draw();
  makeCan();
  TLegend *legend_1 = new TLegend(0.1, 0.7, 0.3, 0.9);
  reco_phi_mass->GetXaxis()->SetTitleSize(0.05);
  reco_phi_mass->GetXaxis()->CenterTitle();
  reco_phi_mass->GetXaxis()->SetTitleOffset(0.8);
  reco_phi_mass->Draw("pfc");
  reco_phi_mass_pt_big_cut->Draw("pfc;same");
  reco_phi_mass_pt_cut->Draw("pfc;same");
  reco_phi_mass_pt_small_cut->Draw("pfc;same");
  legend_1->AddEntry(reco_phi_mass, "All P_{T}");
  legend_1->AddEntry(reco_phi_mass_pt_big_cut, "Kaon P_{T} < 400 MeV");
  legend_1->AddEntry(reco_phi_mass_pt_cut, "Kaon P_{T} < 300 MeV");
  legend_1->AddEntry(reco_phi_mass_pt_small_cut, "Kaon P_{T} < 100 MeV"); 
  legend_1->Draw("same");
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
                     "Pair P_{T} < 100 MeV");
  legend_3->Draw("same");
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
  makeCan();
  TLegend *legend_4 = new TLegend(0.78, 0.55, 0.98, 0.75);
  reco_phi_near_mass->GetXaxis()->SetTitleSize(0.05);
  reco_phi_near_mass->GetXaxis()->CenterTitle();
  reco_phi_near_mass->GetXaxis()->SetTitleOffset(0.8);
  reco_phi_near_mass->Draw("pfc");
  reco_phi_near_mass_phi_pt_big_cut->Draw("pfc;same");
  reco_phi_near_mass_phi_pt_cut->Draw("pfc;same");
  reco_phi_near_mass_phi_pt_small_cut->Draw("pfc;same");
  legend_4->AddEntry(reco_phi_near_mass, "All P_{T}");
  legend_4->AddEntry(reco_phi_near_mass_phi_pt_big_cut, "Pair P_{T} < 500MeV"); 
  legend_4->AddEntry(reco_phi_near_mass_phi_pt_cut, "Pair P_{T} < 300MeV"); 
  legend_4->AddEntry(reco_phi_near_mass_phi_pt_small_cut,
                     "Pair P_{T} < 100 MeV");
  legend_4->Draw("same");
  makeCan();
  reco_phi_mass_pt->GetXaxis()->SetTitleSize(0.05);
  reco_phi_mass_pt->GetXaxis()->CenterTitle();
  reco_phi_mass_pt->GetXaxis()->SetTitleOffset(0.8);
  reco_phi_mass_pt->GetYaxis()->SetTitleSize(0.05);
  reco_phi_mass_pt->GetYaxis()->CenterTitle();
  reco_phi_mass_pt->GetYaxis()->SetTitleOffset(0.8);
  reco_phi_mass_pt->Draw("colz");
  makeCan();
  reco_phi_mass_pt->ProjectionX()->Draw("e1");
  reco_phi_mass_pt->ProjectionX()->Draw("hist;same");
  makeCan();
  reco_phi_pt_mass_cut->GetXaxis()->SetTitleSize(0.05);
  reco_phi_pt_mass_cut->GetXaxis()->CenterTitle();
  reco_phi_pt_mass_cut->GetXaxis()->SetTitleOffset(0.8);
  reco_phi_pt_mass_cut->SetLineWidth(2);
  reco_phi_pt_mass_cut->SetFillColorAlpha(kBlue, 0.35);
  reco_phi_pt_mass_cut->SetMarkerStyle(20); // This is a full circle marker
  reco_phi_pt_mass_cut->SetMarkerColor(kRed); // Set point color to red
  reco_phi_pt_mass_cut->SetMarkerSize(1);  // Adjust the size of the marker as needed
  reco_phi_pt_mass_cut->Draw("e1");
  reco_phi_pt_mass_cut->Draw("hist;same");
  makeCan();
  reco_phi_coherent_eta->Draw();

  makeCan();
  TLegend *legend_5 = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend_5->AddEntry(reco_coherent_phi_mass, "Coherent #phi Candidates");
  legend_5->AddEntry(reco_coherent_phi_mass_tof,
                     "Coherent #phi Candidates TOF");
  reco_coherent_phi_mass->GetXaxis()->SetTitleSize(0.05);
  reco_coherent_phi_mass->GetXaxis()->CenterTitle();
  reco_coherent_phi_mass->GetXaxis()->SetTitleOffset(0.8);
  reco_coherent_phi_mass->SetLineWidth(2);
  reco_coherent_phi_mass->SetFillColorAlpha(kBlue, 0.35);
  reco_coherent_phi_mass->SetMarkerStyle(20); // This is a full circle marker
  reco_coherent_phi_mass->SetMarkerColor(kRed); // Set point color to red
  reco_coherent_phi_mass->SetMarkerSize(1);  // Adjust the size of the marker as needed
  reco_coherent_phi_mass->Draw("e1");
  reco_coherent_phi_mass->Draw("hist;same");
  reco_coherent_phi_mass_tof->Draw("same;pfc");
  legend_5->Draw("same");
  makeCan();
  KpKp_pt->Draw();
  makeCan();
  KpKp_eta->Draw();
  makeCan();
  KpKp_phi->Draw();
  makeCan();
  KpKp_mass->Draw();
  makeCan();
  KmKm_pt->Draw();
  makeCan();
  KmKm_eta->Draw();
  makeCan();
  KmKm_phi->Draw();
  makeCan();
  KmKm_mass->Draw();
  makeCan();
  like_sign_pt->Draw();
  makeCan();
  like_sign_mass->Draw();
  makeCan();
  like_sign_eta->Draw();
  makeCan();
  like_sign_phi->Draw();
  makeCan();
  like_sign_near_mass->GetXaxis()->SetTitleSize(0.05);
  like_sign_near_mass->GetXaxis()->CenterTitle();
  like_sign_near_mass->GetXaxis()->SetTitleOffset(0.8);
  like_sign_near_mass->SetLineWidth(2);
  like_sign_near_mass->SetFillColorAlpha(kBlue, 0.35);
  like_sign_near_mass->SetMarkerStyle(20); // This is a full circle marker
  like_sign_near_mass->SetMarkerColor(kRed); // Set point color to red
  like_sign_near_mass->SetMarkerSize(1);  // Adjust the size of the marker as needed
  like_sign_near_mass->Draw("e1");
  like_sign_near_mass->Draw("same;hist");
}
