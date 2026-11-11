#ifndef SRothman_JetSetup_JetTools_h
#define SRothman_JetSetup_JetTools_h

#include "SRothman/DataFormats/interface/EECJet.h"

#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

// Framework-free jet and constituent selection and systematic variations on srothman::EECJet.
namespace jetsetup {

  // ---- constituent categories ----
  enum Category : std::uint32_t {
    kChargedHadron = 1u << 0,  // charged, |pdgId| not 11 or 13
    kNeutralHadron = 1u << 1,  // neutral, not a photon
    kPhoton = 1u << 2,         // pdgId 22
    kElectron = 1u << 3,       // |pdgId| 11
    kMuon = 1u << 4,           // |pdgId| 13
    kAll = (1u << 5) - 1
  };
  Category category(const srothman::EECParticle& p);
  // Accepts chargedHadron, neutralHadron, photon, electron, muon, all; throws std::invalid_argument otherwise.
  std::uint32_t categoryMask(const std::vector<std::string>& names);

  double pt(const srothman::EECParticle& p);

  // ---- jet selection ----
  struct JetCuts {
    double ptMin = 0;
    double absEtaMax = 1e9;
    std::uint32_t idMask = 0;  // every bit set here must be set in jetIdBits
    std::uint32_t minConstituents = 0;
  };
  bool passes(const srothman::EECJet& j, const JetCuts& c);
  // True if any (eta, phi) is within dR < dRmax of the jet axis.
  bool overlaps(const srothman::EECJet& j, const std::vector<std::pair<double, double>>& etaPhi, double dRmax);
  double deltaR(double eta1, double phi1, double eta2, double phi2);

  // ---- constituent selection ----
  struct ConstituentCuts {
    double ptMin = 0;
    bool chargedOnly = false;
    std::vector<int> keepAbsPdgIds;  // empty: keep all
    std::vector<int> dropAbsPdgIds;
    float minPuppiWeight = 0;
  };
  bool passes(const srothman::EECParticle& p, const ConstituentCuts& c);
  void selectConstituents(srothman::EECJet& j, const ConstituentCuts& c);

  // ---- variations ----
  // Seed from the event id, the source-jet index and a salt; independent of the module, so every
  // variation chain draws the same random numbers for the same jet.
  std::uint64_t seed(std::uint64_t run, std::uint64_t lumi, std::uint64_t event, std::uint64_t jetSrcIdx, std::uint64_t salt);

  // Scales the jet four-momentum (pt and mass) by f; eta and phi are unchanged; f must be >= 0.
  void scaleJet(srothman::EECJet& j, double f);
  // Scales the four-momenta of the constituents in `mask` by f.
  void scaleConstituents(srothman::EECJet& j, std::uint32_t mask, double f);
  // Drops each constituent in `mask` with probability p.
  void dropConstituents(srothman::EECJet& j, std::uint32_t mask, double p, std::mt19937_64& rng);

  // JER hybrid smearing factor (JME recipe): with a matched gen jet, 1 + (sf - 1)(pt - genPt)/pt; otherwise
  // 1 + gauss * resolution * sqrt(max(sf^2 - 1, 0)). Negative factors are clamped to 0.
  // genPt < 0 means unmatched.
  double jerFactor(double pt, double genPt, double sf, double resolution, double gauss);
  // Gen-jet match for JER: dR < R/2 and |pt - genPt| < 3 * resolution * pt.
  bool jerMatch(double dR, double pt, double genPt, double R, double resolution);

}  // namespace jetsetup

#endif
