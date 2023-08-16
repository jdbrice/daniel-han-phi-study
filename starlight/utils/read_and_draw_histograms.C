// Include necessary ROOT headers
#include "THStack.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>

int ican = 0;
void makeCan(bool logz = false) {
  TCanvas *can = new TCanvas(TString::Format("can%d", ican++), "", 900, 600);
  can->SetTopMargin(0.10);
  can->SetRightMargin(0.03);

  if (logz) {
    can->SetLogz(1); // Set Z-axis to logarithmic scale
  }
}

void read_and_draw_histograms() {
  // 1. Read histogram `fPt1` from file "kaon_inco.root"
  TFile *fileInco = TFile::Open("./kaon_inco.root");
  TH1F *histInco = (TH1F *)fileInco->Get("fPt1");

  // 2. Read histogram `fPt1` from file "kaon_co.root"
  TFile *fileCo = TFile::Open("./kaon_co.root");
  TH1F *histCo = (TH1F *)fileCo->Get("fPt1");

  THStack *kaon_stack = new THStack(
      "STARlight Generated Kaon",
      "STARlight Generated Kaon P_{T};K^{#pm} P_{T}(GeV/c);counts");
  // Check if histograms are loaded properly
  if (!histInco || !histCo) {
    printf("Error: Could not load histograms\n");
    return;
  }

  // Create a canvas for drawing
  makeCan();
  gStyle->SetOptStat(0);

  // Draw histograms
  histInco->SetTitle("Kaon Daughter Transverse Momentuim by STARlight");
  histInco->GetXaxis()->SetTitleSize(0.05);
  histInco->GetXaxis()->CenterTitle();
  histInco->GetXaxis()->SetTitleOffset(0.8);

  histCo->SetTitle("Kaon Daughter Transverse Momentuim by STARlight");
  //
  histCo->SetFillColor(kRed);
  histInco->SetFillColor(kBlue);
  kaon_stack->Add(histInco);
  kaon_stack->Add(histCo);

  kaon_stack->Draw("e1");
  kaon_stack->Draw("same;hist");

  kaon_stack->GetHistogram()->GetXaxis()->SetTitleSize(0.05);
  kaon_stack->GetHistogram()->GetXaxis()->CenterTitle();
  kaon_stack->GetHistogram()->GetXaxis()->SetTitleOffset(0.8);


  gPad->Modified();
  gPad->Update();
  // Create a legend
  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9); // x1, y1, x2, y2
  leg->AddEntry(
      histInco, "Incoherent Production",
      "f"); // "f" denotes using the fill color/style for the legend icon
  leg->AddEntry(histCo, "Coherent Production", "f");
  leg->Draw();

  // Save the canvas to an image, if you wish
  // canvas->SaveAs("output.png");
}
