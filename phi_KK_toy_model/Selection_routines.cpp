#include "Selection_routines.h"
#include "Random_routines.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <TTreeReader.h>

// global rng
thread_local TRandom3 rng_selection(0);

// constructor and variable init
Selector::Selector(double sigma) { this->sigma = sigma; }
// constructor and variable init
Selector::Selector() {}

void Selector::draw_dEdx()
{
  std::vector<double> dEdx_kaon_list;
  std::vector<double> dEdx_pion_list;
  std::vector<double> pt;
  // set up momentum mesh
  for (double p = 0.01; p <= 10.; p += 1e-4)
  {
    pt.push_back(p);
  }
  // fill dEdx in the momentum mesh for Kaon
  for (double p : pt)
  {
    dEdx_kaon_list.push_back(dEdxKaon->Interpolate(p));
    dEdx_pion_list.push_back(dEdxPion->Interpolate(p));
  }

  TH2D *dEdx = new TH2D("dEdx Pion", "dEdX v.s. Pt; Pt(GeV); dEdx", 500,
                        0.1, 10., 500, 0.0025, 0.008);

  for (size_t i = 0; i < pt.size(); i++)
  {
    dEdx->Fill(pt[i], dEdx_kaon_list[i]);
    dEdx->Fill(pt[i], dEdx_pion_list[i]);
  }
  dEdx->Draw();
}

void Selector::draw_dEdx_blurred()
{
  std::vector<double> dEdx_kaon_blurred;
  std::vector<double> dEdx_pion_blurred;
  std::vector<double> pt;

  // set up momentum mesh
  for (double p = 0.01; p <= 2.; p += 1e-7)
  {
    pt.push_back(p);
  }
  // fill dEdx in the momentum mesh for Kaon
  for (double p : pt)
  {
    dEdx_kaon_blurred.push_back(rng_selection.Gaus(dEdxKaon->Interpolate(p), sigma));
  }

  // fill dEdx in the momentum mesh for Pion
  for (double p : pt)
  {
    dEdx_pion_blurred.push_back(rng_selection.Gaus(dEdxPion->Interpolate(p), sigma));
  }

  TH2D *dEdx_blurred =
      new TH2D("dEdx Kaon + Pion", "dEdX v.s. Pt; Pt(GeV); dEdx", 500, 0.1, 2.,
               500, 0., 0.014);

  for (size_t i = 0; i < pt.size(); i++)
  {
    double p = pt[i];
    dEdx_blurred->Fill(p, dEdx_kaon_blurred[i]);
    dEdx_blurred->Fill(p, dEdx_pion_blurred[i]);
  }
  dEdx_blurred->Draw("colz");
}

double Selector::get_NSigmaKaon(TLorentzVector *lv_ptr)
{
  double NSigmaKaon;
  double pt = lv_ptr->Pt();
  // retrive exact dEdxKaon from root file
  double dEdxKaon_mc = dEdxKaon->Interpolate(pt);
  // blur the dEdxKaon by an amount that eyeball to fit the
  // experimental data
  double dEdxKaon_rc = rng_selection.Gaus(dEdxKaon_mc, sigma);
  // An assumption is made that the measured particle is kaon
  if (std::abs(lv_ptr->M() - 0.493) <= 0.001)
  {
    NSigmaKaon = (dEdxKaon_rc - dEdxKaon_mc) / sigma;
    return NSigmaKaon;
  }
  else
  {
    std::cout << "Error, particle given dose not have kaon mass"
              << std::endl;
    return 0;
  }
}
double Selector::get_NSigmaElectron(TLorentzVector *lv_ptr)
{
  double NSigmaElectron;
  double pt = lv_ptr->Pt();
  // retrive exact dEdxKaon from root file
  double dEdxElectron_mc = dEdxElectron->Interpolate(pt);
  // blur the dEdxKaon by an amount that eyeball to fit the
  // experimental data
  double dEdxKaon_rc = rng_selection.Gaus(dEdxElectron_mc, sigma);
  // An assumption is made that the measured particle is kaon
  if (std::abs(lv_ptr->M() - 0.493) <= 0.001)
  {
    NSigmaElectron = (dEdxKaon_rc - dEdxElectron_mc) / sigma;
    return NSigmaElectron;
  }
  else
  {
    std::cout << "Error, particle given dose not have kaon mass"
              << std::endl;
    return 0;
  }
}

double Selector::get_NSigmaPion(TLorentzVector *lv_ptr)
{
  double NSigmaPion;
  double pt = lv_ptr->Pt();
  // retrive exact dEdxPion from root file
  double dEdxPion_mc = dEdxPion->Interpolate(pt);
  double dEdxKaon_mc = dEdxKaon->Interpolate(pt);
  // blur the dEdxKaon by an amount that eyeball to fit the
  // experimental data
  double dEdxKaon_rc = rng_selection.Gaus(dEdxKaon_mc, sigma);
  // An assumption is made that the measured particle is kaon
  if (std::abs(lv_ptr->M() - 0.493) <= 0.001)
  {
    NSigmaPion = (dEdxKaon_rc - dEdxPion_mc) / sigma;
    return NSigmaPion;
  }
  else
  {
    std::cout << "Error, particle given dose not have kaon mass"
              << std::endl;
    return 0;
  }
}

// int main(int argc, char **argv)
// {
//   double KAON_MASS = 0.493;
//   double PION_MASS = 0.139;
//   int PHI_SAMPLE_SIZE = 4500;
//   int RHO_SAMPLE_SIZE = PHI_SAMPLE_SIZE * 10;
//   // create vectors for parent vector and daughter vectors
//   std::vector<TLorentzVector *> parent_vector;
//   std::vector<TLorentzVector *> daughter1_vector;
//   std::vector<TLorentzVector *> daughter2_vector;
//
//   double NSigmaPion = 0;
//   double NSigmaKaon = 0;
//
//   Selector PID = Selector();
//
//   // hardcoded starlight histogram files
//   TFile *kaon_file = new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/kaon.root");
//   TFile *pion_file = new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/pion.root");
//
//   // simulate the decay process for rho -> pi+ pi-
//   for (int i = 0; i < RHO_SAMPLE_SIZE; i++)
//   {
//
//     TLorentzVector *rho_parent_particle_ptr =
//         Random_routines::get_slight_lorentz_vector(pion_file, "pion");
//     parent_vector.push_back(rho_parent_particle_ptr);
//
//     // simulate the decay process with the daughter mass being pion
//     std::vector<TLorentzVector *> daughter_ptr_pair =
//         Random_routines::two_body_decay(rho_parent_particle_ptr,
//                                         PION_MASS, PION_MASS);
//
//     daughter_ptr_pair[0]->SetPtEtaPhiM(daughter_ptr_pair[0]->Pt(), daughter_ptr_pair[0]->Eta(), daughter_ptr_pair[0]->Phi(), KAON_MASS);
//     daughter_ptr_pair[1]->SetPtEtaPhiM(daughter_ptr_pair[1]->Pt(), daughter_ptr_pair[1]->Eta(), daughter_ptr_pair[1]->Phi(), KAON_MASS);
//
//     daughter1_vector.push_back(daughter_ptr_pair[0]);
//     daughter2_vector.push_back(daughter_ptr_pair[1]);
//
//     NSigmaKaon += abs(PID.get_NSigmaKaon(daughter_ptr_pair[0]));
//     NSigmaPion += abs(PID.get_NSigmaPion(daughter_ptr_pair[0]));
//   }
//   std::cout << NSigmaKaon / RHO_SAMPLE_SIZE << "    " << NSigmaPion / RHO_SAMPLE_SIZE << std::endl;
//
//   NSigmaKaon = 0;
//   NSigmaPion = 0;
//
//   // simualte decay process for phi -> K+ K-
//   for (int i = 0; i < PHI_SAMPLE_SIZE; i++)
//   {
//     // generate parent vector that is a 4 vector uniformly distributed between
//     // pt, eta, phi bounds with phi mass
//     TLorentzVector *phi_parent_particle_ptr = Random_routines::get_slight_lorentz_vector(kaon_file, "kaon");
//     parent_vector.push_back(phi_parent_particle_ptr);
//
//     // simulate the decay process with the daughter mass being kaon
//     std::vector<TLorentzVector *> daughter_ptr_pair =
//         Random_routines::two_body_decay(phi_parent_particle_ptr,
//                                         KAON_MASS, KAON_MASS);
//
//     daughter1_vector.push_back(daughter_ptr_pair[0]);
//     daughter2_vector.push_back(daughter_ptr_pair[1]);
//     NSigmaKaon += abs(PID.get_NSigmaKaon(daughter_ptr_pair[0]));
//     NSigmaPion += abs(PID.get_NSigmaPion(daughter_ptr_pair[0]));
//   }
//   std::cout << NSigmaKaon / PHI_SAMPLE_SIZE << "    " << NSigmaPion / PHI_SAMPLE_SIZE << std::endl;
//
//   TApplication app("app", &argc, argv);
//   TCanvas *canvas = new TCanvas("canvas", "canvas2", 0, 0, 800, 600);
//   PID.draw_dEdx_blurred();
//   canvas->Modified();
//   canvas->Update();
//   TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
//   root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
//                        "Terminate()");
//   app.Run();
//   return 0;
// }
