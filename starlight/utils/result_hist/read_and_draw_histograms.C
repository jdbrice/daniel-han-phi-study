// Include necessary ROOT headers
#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
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
  TFile *fileInco = TFile::Open("./kaon_incoherent.root");
  TH1F *histInco = (TH1F *)fileInco->Get("PtKa_tpc_2");
  // TH1F *histInco_2 = (TH1F *)fileInco->Get("PtKa_tpc_2");

  // 2. Read histogram `fPt1` from file "kaon_co.root"
  TFile *fileCo = TFile::Open("./kaon_coherent.root");
  TH1F *histCo = (TH1F *)fileCo->Get("PtKa_tpc_2");
  // TH1F *histCo_2 = (TH1F *)fileCo->Get("PtKa_tpc_2");

  // Check if histograms are loaded properly
  if (!histInco || !histCo) {
    printf("Error: Could not load histograms\n");
    return;
  }
  THStack *phi_stack = new THStack("STARlight Generated #phi Mesons","STARlight Generated #phi Meson P_{T};K^{+}K^{-} Pair P_{T};counts");

  TH1F *sampled_coherent = new TH1F("STARlight Coherent Production",
                                    "STARlight P_{T} #phi Meson", 100, 0., 0.9);

  TH1F *sampled_incoherent = new TH1F("STARlight Incoherent Production",
                                    "STARlight P_{T} #phi Meson", 100, 0., 0.9);

  TH1F *sampled_all = new TH1F("STARlight Coherent Production",
                                    "STARlight P_{T} #phi Meson Production; K^{+}K^{-} Pair P_{T}(GeV/c);counts", 100, 0., 0.9);

  for (int i = 0; i < 6700; i ++){
    double coherent_pt = histCo->GetRandom();
    sampled_coherent->Fill(coherent_pt);
    sampled_all->Fill(coherent_pt);
  }

  for (int i = 0; i < 20600; i++){
    double inco_pt = histInco->GetRandom();
    sampled_incoherent->Fill(inco_pt);
    sampled_all->Fill(inco_pt);
  }
  // Create a canvas for drawing
  makeCan();
  gStyle->SetOptStat(0);

  // Draw histograms

  
  sampled_all->GetXaxis()->SetTitleSize(0.05);
  sampled_all->GetXaxis()->CenterTitle();
  sampled_all->GetXaxis()->SetTitleOffset(0.8);
  sampled_all->SetLineWidth(1);

  sampled_incoherent->SetFillColorAlpha(kBlue, 0.35);
  sampled_coherent->SetFillColorAlpha(kRed, 0.35);

  sampled_coherent->Scale(100./6700);
  sampled_incoherent->Scale(100./6700);
  sampled_all->Scale(100./6700);

  phi_stack->Add(sampled_incoherent);
  phi_stack->Add(sampled_coherent);

  sampled_all->SetMarkerStyle(20);
  sampled_all->SetMarkerColor(kRed);
  sampled_all->Draw("e1");
  phi_stack->Draw("hist;same");

  // Create a legend
  TLegend *leg = new TLegend(0.7, 0.7, 0.96, 0.9); // x1, y1, x2, y2
  leg->AddEntry(
      sampled_incoherent, "Incoherent STARlight Production",
      "f"); // "f" denotes using the fill color/style for the legend icon
  leg->AddEntry(sampled_coherent, "Coherent STARlight Production", "f");
  leg->AddEntry(sampled_all, "All STARlight Production", "pe");
  leg->Draw();

  // Save the canvas to an image, if you wish
  // canvas->SaveAs("output.png");
}
