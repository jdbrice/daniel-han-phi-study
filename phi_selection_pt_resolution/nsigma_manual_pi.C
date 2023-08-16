// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)
#include "FemtoPairFormat.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <TDirectory.h>
#include <TFitResultPtr.h>
#include <TH2.h>
#include <TSystem.h>
#include <cmath>
#include <math.h>

// ClassImp( FemtoPair );

int ican = 0;
void makeCan() {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.06);
  can->SetRightMargin(0.01);
}

void nsigma_manual_pi() {

  // histograms for all tracks
  TH1F *all_DCA = new TH1F("Pion Candidates",
                           "Run 19 A+A DCA All Tracks;DCA;counts", 100, 0., 3.);

  TH1F *pion_Npion =
      new TH1F("Pion Candidates",
               "Run 19 A+A N#sigma#pi Pion Candidates;(N#sigmaPion);counts",
               100, -10, 10);

  TH2F *pion_pt_NPion =
      new TH2F("Pion Candidates",
               "Run 19 A+A N#sigmaPion Pion Candidates;P_{T};N#sigmaPion", 100,
               0.25, 0.5, 100, -5, 7);
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

    if (pair->mChargeSum == 0 &&
        (pair->d1_mNSigmaKaon > 7.4 || pair->d1_mNSigmaKaon < -2.6) &&
        (pair->d2_mNSigmaKaon > 7.4 || pair->d2_mNSigmaKaon < -2.6) &&
        (pair->d1_mNSigmaElectron > 5.424112 || pair->d1_mNSigmaElectron < -4.575888) &&
        (pair->d2_mNSigmaElectron > 5.424112 || pair->d2_mNSigmaElectron < -4.575888) &&
        (pair->d1_mNSigmaProton > 7.36 || pair->d1_mNSigmaProton < -2.64) &&
        (pair->d2_mNSigmaProton > 7.36 || pair->d2_mNSigmaProton < -2.64) &&
        abs(pair->d2_mNSigmaPion) < 30 &&
        abs(pair->d1_mNSigmaPion) < 30) {
      pion_Npion->Fill(pair->d1_mNSigmaPion);
      pion_Npion->Fill(pair->d2_mNSigmaPion);

      pion_pt_NPion->Fill(pair->d1_mPt * cosh(pair->d1_mEta), pair->d1_mNSigmaPion);
      pion_pt_NPion->Fill(pair->d2_mPt * cosh(pair->d2_mEta), pair->d2_mNSigmaPion);
    }
  } // loop on events

  std::vector<double> mean_values;
  std::vector<double> pt_values;
  std::vector<double> mean_errors;
  // go through all the bins
  for (int i = 1; i <= pion_pt_NPion->GetNbinsX(); ++i) {
    // Projection Y for each bin in X
    TH1D *y_proj = pion_pt_NPion->ProjectionY("_py", i, i);
    // fit gaussian quitely for each bin
    TFitResultPtr fit_result = y_proj->Fit("gaus", "QS");
    // store the mean_value
    if ((int)fit_result == 0) { // case for successful fit
      mean_values.push_back(fit_result->Parameter(1));
      pt_values.push_back(pion_pt_NPion->GetXaxis()->GetBinCenter(i));
      mean_errors.push_back(fit_result->Error(1));
    } else {
      mean_values.push_back(999);
      pt_values.push_back(pion_pt_NPion->GetXaxis()->GetBinCenter(i));
      mean_errors.push_back(999);
    }
    delete y_proj;
  }
  // get the all momentum nsigmapion result
  TH1D *all_momentum = pion_pt_NPion->ProjectionY();
  TFitResultPtr fitResult =
      all_momentum->Fit("gaus", "QS"); // Fit the projection
  double nsigmapion_mean = fitResult->Parameter(1);
  double nsigmapion_error = fitResult->Error(1);
  TF1 *all_momentum_fit = all_momentum->GetFunction("gaus");

  // graphing the momentum-dependent result
  std::cout << nsigmapion_mean << "     " << nsigmapion_error << std::endl;
  // Create TGraph
  makeCan();
  TGraphErrors *graph = new TGraphErrors(mean_values.size(), &pt_values[0],
                                         &mean_values[0], 0, &mean_errors[0]);
  graph->SetTitle(
      "NSigmaPion Mean vs. Momentum; p(GeV/c); NSigmaPion Mean");
  graph->SetMarkerSize(0.8);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.8);
  graph->SetMarkerStyle(20);
  graph->GetXaxis()->SetTitleSize(0.05);
  graph->GetXaxis()->CenterTitle();
  graph->GetXaxis()->SetTitleOffset(0.8);
  graph->GetYaxis()->SetTitleSize(0.05);
  graph->GetYaxis()->CenterTitle();
  graph->GetYaxis()->SetTitleOffset(0.8);
  graph->Draw("AP"); // Draw as points with axes

  // graphing the all momentum result
  TBox *meanBox =
      new TBox(pt_values.front(), nsigmapion_mean - nsigmapion_error,
               pt_values.back(), nsigmapion_mean + nsigmapion_error);

  meanBox->SetFillColor(kBlue);
  meanBox->SetFillStyle(3004);
  meanBox->Draw("same");

  TLine *meanLine = new TLine(pt_values.front(), nsigmapion_mean,
                              pt_values.back(), nsigmapion_mean);
  meanLine->SetLineColor(kBlue);
  meanLine->Draw("same");
  // Create a legend
  TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
  legend->SetHeader("Legend",
                    "C"); // option "C" allows to center the header
  legend->AddEntry(graph, "Indivisual Momentum Bins", "lep");
  legend->AddEntry(meanBox, "All Momentum", "f");
  legend->Draw("same");
  makeCan();
  pion_pt_NPion->Draw("colz");
  makeCan();
  all_momentum->Draw("ep");
  all_momentum_fit->SetLineColor(kRed);
  all_momentum_fit->Draw("same");

  TLine *mean_line = new TLine(nsigmapion_mean, 0, nsigmapion_mean, 24000);
  mean_line->SetLineColor(kBlue);
  mean_line->Draw("same");

  TBox *box = new TBox(nsigmapion_mean - nsigmapion_error, 0,
                       nsigmapion_mean + nsigmapion_error, 24000);
  box->SetFillColor(kBlue);
  box->SetFillStyle(3004); // semi-transparent
  //
  box->Draw("same");
  // Create a legend
  TLegend *legend2 = new TLegend(0.1, 0.7, 0.48, 0.9);
  legend2->SetHeader("Legend",
                     "C"); // option "C" allows to center the header
  legend2->AddEntry(all_momentum, "All Momentum Data", "lep");
  legend2->AddEntry(all_momentum_fit, "All Momentum Fit", "l");
  legend2->Draw("same");
}
