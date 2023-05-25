#include "selection.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRootCanvas.h>
#include <TTreeReader.h>
#include <TH1F.h>

// constructor realization
Selector::Selector(std::string file_path, double sigma) {
  this->file_path = file_path;
  this->sigma = sigma;
}

void Selector::blur_dEdx_file() {
  TFile *dEdx_file = new TFile(file_path.c_str());

  dEdx_file->ls();

  TH1D *test = new TH1D("test","test", 400, 0, 5);

  test = (TH1D*)dEdx_file->Get("mpmK");

  for (int i=0; i <= test->GetNbinsX()+1; ++i){
    double binContent = test->GetBinContent(i);
    double binCenter = test->GetBinCenter(i);
    std::cout << "Bin " << i << ": Content = " << binContent << ", Center = " << binCenter << std::endl;
  }

  test->Draw();
}

int main(int argc, char **argv) {
  TApplication app("app", &argc, argv);
  std::string file_name = "dEdx.root";
  Selector test_selector = Selector(file_name, 1);
  TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
  TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
  root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                       "Terminate()");
  test_selector.blur_dEdx_file();
  app.Run();
  return 0;
}
