#include "SRothman/JetSetup/interface/JetTools.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace jetsetup {

  Category category(const srothman::EECParticle& p) {
    const int a = std::abs(p.pdgId);
    if (a == 11)
      return kElectron;
    if (a == 13)
      return kMuon;
    if (a == 22)
      return kPhoton;
    return p.charge != 0 ? kChargedHadron : kNeutralHadron;
  }

  std::uint32_t categoryMask(const std::vector<std::string>& names) {
    std::uint32_t m = 0;
    for (const auto& n : names) {
      if (n == "chargedHadron")
        m |= kChargedHadron;
      else if (n == "neutralHadron")
        m |= kNeutralHadron;
      else if (n == "photon")
        m |= kPhoton;
      else if (n == "electron")
        m |= kElectron;
      else if (n == "muon")
        m |= kMuon;
      else if (n == "all")
        m |= kAll;
      else
        throw std::invalid_argument("jetsetup: unknown constituent category '" + n + "'");
    }
    return m;
  }

  double pt(const srothman::EECParticle& p) { return std::hypot(p.px, p.py); }

  bool passes(const srothman::EECJet& j, const JetCuts& c) {
    return j.pt >= c.ptMin && std::abs(j.eta) <= c.absEtaMax && (j.jetIdBits & c.idMask) == c.idMask &&
           j.constituents.size() >= c.minConstituents;
  }

  double deltaR(double eta1, double phi1, double eta2, double phi2) {
    const double dphi = std::remainder(phi1 - phi2, 2 * M_PI);
    return std::hypot(eta1 - eta2, dphi);
  }

  bool overlaps(const srothman::EECJet& j, const std::vector<std::pair<double, double>>& etaPhi, double dRmax) {
    for (const auto& [eta, phi] : etaPhi)
      if (deltaR(j.eta, j.phi, eta, phi) < dRmax)
        return true;
    return false;
  }

  bool passes(const srothman::EECParticle& p, const ConstituentCuts& c) {
    if (pt(p) < c.ptMin || p.puppiWeight < c.minPuppiWeight)
      return false;
    if (c.chargedOnly && p.charge == 0)
      return false;
    const int a = std::abs(p.pdgId);
    if (!c.keepAbsPdgIds.empty() && std::find(c.keepAbsPdgIds.begin(), c.keepAbsPdgIds.end(), a) == c.keepAbsPdgIds.end())
      return false;
    return std::find(c.dropAbsPdgIds.begin(), c.dropAbsPdgIds.end(), a) == c.dropAbsPdgIds.end();
  }

  void selectConstituents(srothman::EECJet& j, const ConstituentCuts& c) {
    auto& v = j.constituents;
    v.erase(std::remove_if(v.begin(), v.end(), [&](const auto& p) { return !passes(p, c); }), v.end());
  }

  std::uint64_t seed(std::uint64_t run, std::uint64_t lumi, std::uint64_t event, std::uint64_t jetSrcIdx, std::uint64_t salt) {
    // splitmix64 finaliser applied to each word in turn.
    auto mix = [](std::uint64_t x) {
      x += 0x9e3779b97f4a7c15ull;
      x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
      x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
      return x ^ (x >> 31);
    };
    std::uint64_t h = mix(salt);
    for (auto w : {run, lumi, event, jetSrcIdx})
      h = mix(h ^ w);
    return h;
  }

  void scaleJet(srothman::EECJet& j, double f) {
    if (!(f >= 0))
      throw std::invalid_argument("jetsetup::scaleJet: negative or NaN factor");
    j.pt *= f;
    j.mass *= f;
  }

  void scaleConstituents(srothman::EECJet& j, std::uint32_t mask, double f) {
    for (auto& p : j.constituents)
      if (category(p) & mask) {
        p.px *= f;
        p.py *= f;
        p.pz *= f;
        p.E *= f;
      }
  }

  void dropConstituents(srothman::EECJet& j, std::uint32_t mask, double prob, std::mt19937_64& rng) {
    // One draw per in-mask constituent, in order, so the outcome depends only on the seed.
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<srothman::EECParticle> kept;
    kept.reserve(j.constituents.size());
    for (const auto& p : j.constituents)
      if (!(category(p) & mask) || u(rng) >= prob)
        kept.push_back(p);
    j.constituents = std::move(kept);
  }

  double jerFactor(double pt, double genPt, double sf, double resolution, double gauss) {
    const double f = genPt >= 0 ? 1.0 + (sf - 1.0) * (pt - genPt) / pt
                                : 1.0 + gauss * resolution * std::sqrt(std::max(sf * sf - 1.0, 0.0));
    return std::max(f, 0.0);
  }

  bool jerMatch(double dR, double pt, double genPt, double R, double resolution) {
    return dR < R / 2 && std::abs(pt - genPt) < 3 * resolution * pt;
  }

}  // namespace jetsetup
