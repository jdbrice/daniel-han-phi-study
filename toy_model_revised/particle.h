#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <Math/Vector4D.h>
#include <vector>
#include <memory>

class Particle{
  static int counter;
  std::vector<std::shared_ptr<Particle>> daughters;
  int pdgCode;
  int UID;
  bool reconstructed = 0;

public:
  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiMVector> fourVector;

  Particle(int pdgCode, ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiMVector> fourVector);

  void addDaughter(std::shared_ptr<Particle>& daughter);

  std::vector<std::shared_ptr<Particle>>& getDaughters();

  int getUID();

  int getPdgCode();

  bool getRecoStatus();

  void setReconstructed(bool recoStatus);
};

