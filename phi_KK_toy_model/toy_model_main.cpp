#include "Random_routines.h"
#include "Selection_routines.h"
#include <RtypesCore.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TRootCanvas.h>
#include <cmath>
#include <fstream>
#include <iostream>

thread_local TRandom3 rng_toy(0);

double KAON_MASS = 0.493;
double PION_MASS = 0.139;
double ELECTRON_MASS = 0.000511;
int PHI_SAMPLE_SIZE = 1500;
int RHO_SAMPLE_SIZE = PHI_SAMPLE_SIZE * 15;
int ElECTRON_SAMPLE_SIZE = PHI_SAMPLE_SIZE;

int main(int argc, char **argv)
{
    // create vectors for parent vector and daughter vectors
    std::vector<TLorentzVector *> parent_ptr_vector;
    std::vector<TLorentzVector *> pdaughter_ptr_vector;
    std::vector<TLorentzVector *> mdaughter_ptr_vector;

    // hardcoded starlight histogram files
    TFile *kaon_file =
        new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/kaon.root");
    TFile *pion_file =
        new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/pion.root");
    TFile *electron_file =
        new TFile("/home/xihe/daniel-han-phi-study/starlight_hist/electron.root");

    // simulate the decay process for rho -> pi+ pi-
    for (int i = 0; i < RHO_SAMPLE_SIZE; i++)
    {

        TLorentzVector *rho_parent_particle_ptr =
            Random_routines::get_slight_lorentz_vector(pion_file, "pion");
        parent_ptr_vector.push_back(rho_parent_particle_ptr);

        // simulate the decay process with the daughter mass being pion
        std::vector<TLorentzVector *> daughter_ptr_pair =
            Random_routines::two_body_decay(rho_parent_particle_ptr, PION_MASS,
                                            PION_MASS);

        // blur the daughter particles by 2 percent to simulate actual particle
        // detector accuracy
        Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                               0.04 * daughter_ptr_pair[0]->Pt());
        Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                               0.04 * daughter_ptr_pair[1]->Pt());

        pdaughter_ptr_vector.push_back(daughter_ptr_pair[0]);
        mdaughter_ptr_vector.push_back(daughter_ptr_pair[1]);

        // now change the daughter mass to be Kaon- the assumed particle when doing
        // Kaon selections
        pdaughter_ptr_vector[i]->SetPtEtaPhiM(
            pdaughter_ptr_vector[i]->Pt(), pdaughter_ptr_vector[i]->Eta(),
            pdaughter_ptr_vector[i]->Phi(), KAON_MASS);
        mdaughter_ptr_vector[i]->SetPtEtaPhiM(
            mdaughter_ptr_vector[i]->Pt(), mdaughter_ptr_vector[i]->Eta(),
            mdaughter_ptr_vector[i]->Phi(), KAON_MASS);
    }

    // simualte decay process for phi -> K+ K-
    for (int i = 0; i < PHI_SAMPLE_SIZE; i++)
    {
        // generate parent vector that is a 4 vector uniformly distributed between
        // pt, eta, phi bounds with phi mass
        TLorentzVector *phi_parent_particle_ptr =
            Random_routines::get_slight_lorentz_vector(kaon_file, "kaon");
        parent_ptr_vector.push_back(phi_parent_particle_ptr);

        // simulate the decay process with the daughter mass being kaon
        std::vector<TLorentzVector *> daughter_ptr_pair =
            Random_routines::two_body_decay(phi_parent_particle_ptr, KAON_MASS,
                                            KAON_MASS);

        // blur the daughter particles by 4 percent to simulate actual particle
        // detector accuracy

        // only the real phi daughter with pt > 0.06 is blurred because they will
        // not be selected in dedx code, but shouldn't be blurred twice in TOF
        // process.
        if (daughter_ptr_pair[0]->Pt() > 0.06)
        {
            Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                                   0.04 * daughter_ptr_pair[0]->Pt());
        }

        if (daughter_ptr_pair[1]->Pt() > 0.06)
        {
            Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                                   0.04 * daughter_ptr_pair[1]->Pt());
        }
        pdaughter_ptr_vector.push_back(daughter_ptr_pair[0]);
        mdaughter_ptr_vector.push_back(daughter_ptr_pair[1]);
    }

    for (int i = 0; i < ElECTRON_SAMPLE_SIZE; i++)
    {
        // generate parent vector that is a 4 vector uniformly distributed between
        // pt, eta, phi bounds with phi mass
        TLorentzVector *combined_electron_ptr =
            Random_routines::get_slight_lorentz_vector(electron_file, "electron");
        parent_ptr_vector.push_back(combined_electron_ptr);

        // simulate the decay process with the daughter mass being pion
        std::vector<TLorentzVector *> daughter_ptr_pair =
            Random_routines::two_body_decay(combined_electron_ptr, ELECTRON_MASS, ELECTRON_MASS);

        // blur the daughter particles by 2 percent to simulate actual particle
        // detector accuracy
        Random_routines::add_gaussian_pt_error(daughter_ptr_pair[0],
                                               0.04 * daughter_ptr_pair[0]->Pt());
        Random_routines::add_gaussian_pt_error(daughter_ptr_pair[1],
                                               0.04 * daughter_ptr_pair[1]->Pt());

        pdaughter_ptr_vector.push_back(daughter_ptr_pair[0]);
        mdaughter_ptr_vector.push_back(daughter_ptr_pair[1]);

        int electron_starting_index = PHI_SAMPLE_SIZE + RHO_SAMPLE_SIZE;
        // now change the daughter mass to be Kaon- the assumed particle when doing
        // Kaon selections
        pdaughter_ptr_vector[i + electron_starting_index]->SetPtEtaPhiM(
            pdaughter_ptr_vector[i + electron_starting_index]->Pt(), pdaughter_ptr_vector[i + electron_starting_index]->Eta(),
            pdaughter_ptr_vector[i + electron_starting_index]->Phi(), KAON_MASS);
        mdaughter_ptr_vector[i + electron_starting_index]->SetPtEtaPhiM(
            mdaughter_ptr_vector[i + electron_starting_index]->Pt(), mdaughter_ptr_vector[i + electron_starting_index]->Eta(),
            mdaughter_ptr_vector[i + electron_starting_index]->Phi(), KAON_MASS);
    }
    // create histogram to draw rc mass
    TH1F *parent_rc_mass_total = new TH1F(
        "Combined Masses", "Toy Model Reco Mass;m_{K^+} + m_{K^-}(GeV);count",
        100, 0.96, 1.8);
    TH1F *parent_rc_mass_mc_phi = new TH1F(
        "Combined Masses", "Toy Model Reco Mass;m_{K^+} + m_{K^-}(GeV);count",
        100, 0.96, 1.8);
    TH1F *parent_rc_mass_mc_rho = new TH1F(
        "Combined Masses", "Toy Model Reco Mass;m_{K^+} + m_{K^-}(GeV);count",
        100, 0.96, 1.8);
    TH1F *parent_rc_mass_mc_electron = new TH1F(
        "Combined Masses", "Toy Model Reco Mass;m_{K^+} + m_{K^-}(GeV);count",
        100, 0.96, 1.8);

    // create histogram to draw rc mass
    TH2F *parent_rc_mass_pt_total = new TH2F(
        "Combined Masses", "Toy Model Combined Masses;m_{K^+ K^-}(GeV);RC P_T",
        100, 0.98, 1.1, 100, 0, 0.5);

    std::vector<double> daughter_mass_vector;
    // create an instance of particle selector
    Selector pid = Selector();

    // select daughter particle to be Kaon. This is for case where daughter
    // particle has momentum > 60 MeV
    for (int i = 0; i < RHO_SAMPLE_SIZE; i++)
    {
        if (std::abs(pid.get_NSigmaKaon(pdaughter_ptr_vector[i])) < 50. &&
            std::abs(pid.get_NSigmaKaon(mdaughter_ptr_vector[i])) < 50. &&
            std::abs(pid.get_NSigmaPion(pdaughter_ptr_vector[i])) > 50. &&
            std::abs(pid.get_NSigmaPion(mdaughter_ptr_vector[i])) > 50. &&
            pdaughter_ptr_vector[i]->Pt() > 0.06 &&
            mdaughter_ptr_vector[i]->Pt() > 0.06)
        {
            // reconstruct the dauther particles only if they are kaons.
            // This effectively selects phi
            TLorentzVector reconstructed_parent =
                *(pdaughter_ptr_vector[i]) + *(mdaughter_ptr_vector[i]);
            // fill the reconstructed parent mass
            parent_rc_mass_total->Fill(reconstructed_parent.M());
            parent_rc_mass_mc_rho->Fill(reconstructed_parent.M());
            parent_rc_mass_pt_total->Fill(reconstructed_parent.M(),
                                          reconstructed_parent.Pt());
        }
    }

    for (int i = RHO_SAMPLE_SIZE; i < RHO_SAMPLE_SIZE + PHI_SAMPLE_SIZE; i++)
    {
        if (std::abs(pid.get_NSigmaKaon(pdaughter_ptr_vector[i])) < 50. &&
            std::abs(pid.get_NSigmaKaon(mdaughter_ptr_vector[i])) < 50. &&
            std::abs(pid.get_NSigmaPion(pdaughter_ptr_vector[i])) > 50. &&
            std::abs(pid.get_NSigmaPion(mdaughter_ptr_vector[i])) > 50. &&
            pdaughter_ptr_vector[i]->Pt() > 0.06 &&
            mdaughter_ptr_vector[i]->Pt() > 0.06)
        {
            // reconstruct the dauther particles only if they are kaons.
            // This effectively selects phi
            TLorentzVector reconstructed_parent =
                *(pdaughter_ptr_vector[i]) + *(mdaughter_ptr_vector[i]);
            // fill the reconstructed parent mass
            parent_rc_mass_total->Fill(reconstructed_parent.M());
            parent_rc_mass_mc_phi->Fill(reconstructed_parent.M());
            parent_rc_mass_pt_total->Fill(reconstructed_parent.M(),
                                          reconstructed_parent.Pt());
        }
    }
    for (int i = RHO_SAMPLE_SIZE + PHI_SAMPLE_SIZE; i < parent_ptr_vector.size(); i++)
    {
        if (std::abs(pid.get_NSigmaKaon(pdaughter_ptr_vector[i])) < 50. &&
            std::abs(pid.get_NSigmaKaon(mdaughter_ptr_vector[i])) < 50. &&
            std::abs(pid.get_NSigmaPion(pdaughter_ptr_vector[i])) > 50. &&
            std::abs(pid.get_NSigmaPion(mdaughter_ptr_vector[i])) > 50. &&
            pdaughter_ptr_vector[i]->Pt() > 0.06 &&
            mdaughter_ptr_vector[i]->Pt() > 0.06)
        {
            // reconstruct the daughter particles only if they are kaons.
            // This effectively selects phi
            TLorentzVector reconstructed_parent =
                *(pdaughter_ptr_vector[i]) + *(mdaughter_ptr_vector[i]);
            // fill the reconstructed parent mass
            parent_rc_mass_total->Fill(reconstructed_parent.M());
            parent_rc_mass_mc_electron->Fill(reconstructed_parent.M());
            parent_rc_mass_pt_total->Fill(reconstructed_parent.M(),
                                          reconstructed_parent.Pt());
        }
    }

    // drawing the result
    TApplication app("app", &argc, argv);
    TCanvas *canvas = new TCanvas("canvas", "masses", 0, 0, 1366, 768);
    parent_rc_mass_total->Draw("pfc");
    parent_rc_mass_mc_phi->Draw("same;pfc");
    parent_rc_mass_mc_rho->Draw("same;pfc");
    parent_rc_mass_mc_electron->Draw("same;pfc");

    auto legend = new TLegend(0.7, 0.55, 0.98, 0.75);
    legend->SetHeader("Reco Parent Mass");
    legend->AddEntry(parent_rc_mass_total, "Total Entries", "f");
    legend->AddEntry(parent_rc_mass_mc_phi, "MC #phi", "f");
    legend->AddEntry(parent_rc_mass_mc_rho, "MC #rho", "f");
    legend->AddEntry(parent_rc_mass_mc_electron, "MC electron pair", "f");
    legend->Draw();
    canvas->Modified();
    canvas->Update();

    TRootCanvas *root_canvas = (TRootCanvas *)canvas->GetCanvasImp();
    root_canvas->Connect("CloseWindow()", "TApplication", gApplication,
                         "Terminate()");
    app.Run();
    return 0;
}

void add_electron_to_histogram(TH1F *hist, int peak_count, double mass_min,
                               double mass_max) {}
