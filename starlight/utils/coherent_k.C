// Include necessary ROOT headers
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

void coherent_k() {

  // 2. Read histogram `fPt1` from file "kaon_co.root"
  TFile *fileCo = TFile::Open("./kaon_co_narrow.root");
  TH1F *kaon_co = (TH1F *)fileCo->Get("fPt1");
  TH1F *phi_co = (TH1F *)fileCo->Get("PtKa");

  // Check if histograms are loaded properly
  if (!kaon_co || !phi_co) {
    printf("Error: Could not load histograms\n");
    return;
  }

  // Create a canvas for drawing
  makeCan();
  // Draw histograms
  kaon_co->SetTitle("Kaon Daughter Transverse Momentuim by STARlight");
  kaon_co->GetXaxis()->SetTitleSize(0.05);
  kaon_co->GetXaxis()->CenterTitle();
  kaon_co->GetXaxis()->SetTitleOffset(0.8);
  kaon_co->SetLineWidth(2);
  kaon_co->Draw();
  kaon_co->Draw("pe;same");
  gPad->Print("slight_phi_coherent_pt.png");

  makeCan();
  // Draw histograms
  phi_co->SetTitle("Kaon Pair Transverse Momentuim by STARlight");
  phi_co->GetXaxis()->SetTitleSize(0.05);
  phi_co->GetXaxis()->CenterTitle();
  phi_co->GetXaxis()->SetTitleOffset(0.8);
  phi_co->SetLineWidth(2);
  phi_co->Draw();
  phi_co->Draw("pe;same");
  gPad->Print("slight_kaon_coherent_pt.png");
}
