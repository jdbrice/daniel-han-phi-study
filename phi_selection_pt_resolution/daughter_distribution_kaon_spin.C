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
#include <algorithm>
#include <cmath>
#include <map>
#include <math.h>

int ican = 0;
void makeCan() {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.06);
  can->SetRightMargin(0.01);
}

TRandom3 pair_rng;
TRandom3 inco_pair_rng;
TRandom3 mixed_rng;

double calc_Phi(TLorentzVector lv1, TLorentzVector lv2) {
  TLorentzVector lv_sum = lv1 + lv2;
  lv1.Boost(-lv_sum.BoostVector());
  lv2.Boost(-lv_sum.BoostVector());
  TLorentzVector lv_diff = lv1 - lv2;
  double p_dot_m = lv_sum.Px() * lv_diff.Px() + lv_sum.Py() * lv_diff.Py();
  double p_cross_m = lv_sum.Px() * lv_diff.Py() - lv_sum.Py() * lv_diff.Px();
  double cosphi = p_dot_m / (lv_sum.Perp() * lv_diff.Perp());
  double phi = acos(cosphi);
  // flip phi based on +- z direction of \vec{p} \cross \vec{m}
  phi = p_cross_m < 0 ? -phi : phi;
  return phi;
}

TH1F *normalize_histogram(TH1 *histogram) {
  TH1F *normalized_hist = (TH1F *)histogram->Clone();
  normalized_hist->GetListOfFunctions()->Delete();
  normalized_hist->Sumw2();
  normalized_hist->Scale(1. / normalized_hist->Integral());
  return normalized_hist;
}

void change_bin_structure(TH1F *target, TH1F *temp) {

  int nBins = temp->GetNbinsX();
  double x_min = temp->GetXaxis()->GetXmin();
  double x_max = temp->GetXaxis()->GetXmax();

  target->SetBins(nBins, x_min, x_max);
}

double fit_function(double *x, double *par) {
  double fit_value = par[0] * (1 + par[1] * cos(x[0]) + par[2] * cos(2 * x[0]) +
                               par[3] * cos(3 * x[0]) + par[4] * cos(4 * x[0]));
  return fit_value;
}

void daughter_distribution_kaon_spin() {

  TH1F *phi_distribution_coherent =
      new TH1F("Coherent #phi Candidates",
               "Run 19 Au + Au #Delta#phi Coherent #phi Meson Candidates; "
               "#phi(rad); counts",
               12, -3.15, 3.15);

  TH1F *phi_distribution_incoherent =
      new TH1F("Incoherent #phi Candidates",
               "Run 19 Au + Au #Delta#phi Incoherent #phi Meson Candidates; "
               "#phi(rad); counts",
               12, -3.15, 3.15);

  TH1F *phi_distribution_mixed_fit =
      new TH1F("Run 19 Au + Au Mixed Events",
               "Run 19 Au + Au Mixed Event #Delta#Phi; #phi(rad); counts", 100,
               -3.15, 3.15);

  TH1F *phi_distribution_mixed =
      new TH1F("Run 19 Au + Au Mixed Events",
               "Run 19 Au + Au Mixed Event #Delta#Phi; #phi(rad); counts", 12,
               -3.15, 3.15);

  TH2F *pt_two_phi = new TH2F("Run 19 Au + Au Mixed Event",
                              "Run 19 Au + Au Mixed Event cos2#phi v.s. Pair "
                              "P_{T};P_{T} (GeV/c);cos2#phi",
                              10, 0., 0.5, 10, -50, 50);

  // Open the file containing the tree (INPUT data).
  TFile *myFile = TFile::Open("input.root");

  // This setup the reader, access the data
  TTreeReader myReader("PairDst", myFile);

  // We need this so ROOT understands the data format
  gSystem->Load("FemtoPairFormat_h.so");
  TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

  // Lorentz vectors (4-vectors) for kinematics in special relativity
  TLorentzVector lv1, lv2, lv;

  std::vector<TLorentzVector> d1_list;
  std::vector<TLorentzVector> d2_list;

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
      phi_distribution_coherent->Fill(phi);
      d1_list.push_back(lv1);
      d2_list.push_back(lv2);
    }
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
        pair->mChargeSum == 0 && lv.Pt() > 0.2 && lv.Pt() < 0.9 &&
        lv.M() >= 1.0 && lv.M() < 1.04) {

      // shuffle the two lv for each events
      // FOR REAL EVENTS:
      double pair_num = inco_pair_rng.Integer(2);
      double phi = calc_Phi(lorentz_vector_list[pair_num],
                            lorentz_vector_list[1 - pair_num]);
      phi_distribution_incoherent->Fill(phi);
    }

  } // loop on events
  TF1 *f1 = new TF1("f1", fit_function, -TMath::Pi(), TMath::Pi(), 5);
  TFitResultPtr fit_result =
      phi_distribution_coherent->Fit("f1", "S"); // Fit and save the result
  TLegend *legend = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend->SetHeader("Fit results:", "C");
  legend->AddEntry((TObject *)0,
                   Form("norm: %.3f #pm %.3f", fit_result->Parameter(0),
                        fit_result->Error(0)),
                   "");
  legend->AddEntry((TObject *)0,
                   Form("1#phi: %.3f #pm %.3f", fit_result->Parameter(1),
                        fit_result->Error(1)),
                   "");
  legend->AddEntry((TObject *)0,
                   Form("2#phi: %.3f #pm %.3f", fit_result->Parameter(2),
                        fit_result->Error(2)),
                   "");
  legend->AddEntry((TObject *)0,
                   Form("3#phi: %.3f #pm %.3f", fit_result->Parameter(3),
                        fit_result->Error(3)),
                   "");
  legend->AddEntry((TObject *)0,
                   Form("4#phi: %.3f #pm %.3f", fit_result->Parameter(4),
                        fit_result->Error(4)),
                   "");
  phi_distribution_coherent->Draw("pe");
  legend->Draw("same");
  gPad->Print("./Plots_sigma/phi_delta_phi.png");

  for (int i = 0; i < d1_list.size(); i++) {
    std::vector<TLorentzVector> d1_buffer;

    for (int j = 0; j < d1_list.size(); j++) {
      if (j != i) {
        d1_buffer.push_back(d1_list[j]);
        d1_buffer.push_back(d2_list[j]);
      }
    }

    for (TLorentzVector buffer_lv : d1_buffer) {
      std::vector<TLorentzVector> lorentz_vector_list;
      lorentz_vector_list.push_back(buffer_lv);
      lorentz_vector_list.push_back(d1_list[i]);
      double pair_num = mixed_rng.Integer(2);
      double phi = calc_Phi(lorentz_vector_list[pair_num],
                            lorentz_vector_list[1 - pair_num]);
      if ((lorentz_vector_list[0] + lorentz_vector_list[1]).Pt() < 0.1) {
        phi_distribution_mixed->Fill(phi);
        phi_distribution_mixed_fit->Fill(phi);
      }
    }
  }

  makeCan();
  TFitResultPtr fit_result2 =
      phi_distribution_mixed_fit->Fit("f1", "S"); // Fit and save the result
  TLegend *legend2 = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend2->SetHeader("Fit results:", "C");
  legend2->AddEntry((TObject *)0,
                    Form("norm: %.3f #pm %.3f", fit_result2->Parameter(0),
                         fit_result2->Error(0)),
                    "");
  legend2->AddEntry((TObject *)0,
                    Form("1#phi: %.3f #pm %.3f", fit_result2->Parameter(1),
                         fit_result2->Error(1)),
                    "");
  legend2->AddEntry((TObject *)0,
                    Form("2#phi: %.3f #pm %.3f", fit_result2->Parameter(2),
                         fit_result2->Error(2)),
                    "");
  legend2->AddEntry((TObject *)0,
                    Form("3#phi: %.3f #pm %.3f", fit_result2->Parameter(3),
                         fit_result2->Error(3)),
                    "");
  legend2->AddEntry((TObject *)0,
                    Form("4#phi: %.3f #pm %.3f", fit_result2->Parameter(4),
                         fit_result2->Error(4)),
                    "");

  phi_distribution_mixed_fit->Draw("pe");
  legend2->Draw("same");
  gPad->Print("./Plots_sigma/phi_delta_phi_mixed.png");

  makeCan();
  TFitResultPtr fit_result3 =
      phi_distribution_incoherent->Fit("f1", "S"); // Fit and save the result
  TLegend *legend3 = new TLegend(0.1, 0.7, 0.3, 0.9);
  legend3->SetHeader("Fit results:", "C");
  legend3->AddEntry((TObject *)0,
                    Form("norm: %.3f #pm %.3f", fit_result3->Parameter(0),
                         fit_result3->Error(0)),
                    "");
  legend3->AddEntry((TObject *)0,
                    Form("1#phi: %.3f #pm %.3f", fit_result3->Parameter(1),
                         fit_result3->Error(1)),
                    "");
  legend3->AddEntry((TObject *)0,
                    Form("2#phi: %.3f #pm %.3f", fit_result3->Parameter(2),
                         fit_result3->Error(2)),
                    "");
  legend3->AddEntry((TObject *)0,
                    Form("3#phi: %.3f #pm %.3f", fit_result3->Parameter(3),
                         fit_result3->Error(3)),
                    "");
  legend3->AddEntry((TObject *)0,
                    Form("4#phi: %.3f #pm %.3f", fit_result3->Parameter(4),
                         fit_result3->Error(4)),
                    "");

  phi_distribution_incoherent->Draw("pe");
  legend3->Draw("same");
  gPad->Print("./Plots_sigma/phi_delta_phi_inco.png");

  TH1F *normalized_mixed = normalize_histogram(phi_distribution_mixed);

  TH1F *normalized_coherent = normalize_histogram(phi_distribution_coherent);

  TH1F *normalized_coherent_corrected = (TH1F *)normalized_coherent->Clone();
  normalized_coherent_corrected->Add(normalized_mixed, -1);

  makeCan();
  normalized_coherent_corrected->Draw();
  TFitResultPtr fit_result4 =
      normalized_coherent_corrected->Fit("f1", "S"); // Fit and save the result
  TLegend *legend4 = new TLegend(0.1, 0.7, 0.4, 0.9);
  legend4->SetHeader("Fit results:", "C");
  legend4->AddEntry((TObject *)0,
                    Form("norm: %.3f #pm %.3f", fit_result4->Parameter(0),
                         fit_result4->Error(0)),
                    "");
  legend4->AddEntry((TObject *)0,
                    Form("1#phi: %.3f #pm %.3f", fit_result4->Parameter(1),
                         fit_result4->Error(1)),
                    "");
  legend4->AddEntry((TObject *)0,
                    Form("2#phi: %.3f #pm %.3f", fit_result4->Parameter(2),
                         fit_result4->Error(2)),
                    "");
  legend4->AddEntry((TObject *)0,
                    Form("3#phi: %.3f #pm %.3f", fit_result4->Parameter(3),
                         fit_result4->Error(3)),
                    "");
  legend4->AddEntry((TObject *)0,
                    Form("4#phi: %.3f #pm %.3f", fit_result4->Parameter(4),
                         fit_result4->Error(4)),
                    "");
  legend4->Draw("same");
  gPad->Print("./Plots_sigma/phi_delta_phi_corrected.png");
}
