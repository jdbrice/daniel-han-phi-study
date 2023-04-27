#include "Random_routines.h"
#include <TLorentzVector.h>
#include <vector>
#include <cmath>

double PI = std::atan(1.) * 4.;
double MAX_ITREATION = 10000;
double PT_MIN = 0.;
double PT_MAX = 3.;
double PETA_MIN = 0.;
double PETA_MAX = 3.;
double M_KAON = 0.493;
double M_PHI = 1.019;

int main(){
  std::vector<TLorentzVector*> lv1_list;
  std::vector<TLorentzVector*> lv2_list;
  std::vector<TLorentzVector*> lv_list;

  for (int i = 0; i < MAX_ITREATION; i++){
    TLorentzVector* lv1 = Random_routines::get_uniform_randomized_lorentz_ptr(M_KAON, PT_MIN, PT_MAX, PETA_MIN, PETA_MAX, 0., 2. * PI, 0);
    lv1_list.push_back(lv1);
    TLorentzVector* lv2 = Random_routines::get_uniform_randomized_lorentz_ptr(M_KAON, PT_MIN, PT_MAX, PETA_MIN, PETA_MAX, 0., 2. * PI, 0);
    lv2_list.push_back(lv2);
    TLorentzVector* lv = new TLorentzVector(*lv1 + *lv2);
    lv_list.push_back(lv);
    delete lv;
  }

  lv1_list[100]->Print();
  lv2_list[100]->Print();
  lv_list[100]->Print();


}
