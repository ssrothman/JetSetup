#include "catch.hpp"

#include "SRothman/JetSetup/interface/JetTools.h"

#include <cmath>
#include <set>

using namespace jetsetup;
using srothman::EECJet;
using srothman::EECParticle;

namespace {
  EECParticle part(double px, double py, double pz, int q, int pdg, float w = 1.0f) {
    return {px, py, pz, std::sqrt(px * px + py * py + pz * pz), static_cast<std::int8_t>(q), pdg, w, 0};
  }
  EECJet mixedJet() {
    EECJet j;
    j.pt = 300;
    j.eta = 0.5;
    j.phi = 1.0;
    j.mass = 40;
    j.constituents = {part(10, 0, 1, 1, 211),    // charged hadron
                      part(0, 5, 2, 0, 130),     // neutral hadron
                      part(3, 3, 0, 0, 22),      // photon
                      part(1, 1, 1, -1, 11),     // electron
                      part(2, 0, 0, 1, -13),     // muon
                      part(0.5, 0, 0, -1, -211, 0.2f)};
    return j;
  }
}  // namespace

TEST_CASE("category and categoryMask", "[jettools]") {
  const auto j = mixedJet();
  CHECK(category(j.constituents[0]) == kChargedHadron);
  CHECK(category(j.constituents[1]) == kNeutralHadron);
  CHECK(category(j.constituents[2]) == kPhoton);
  CHECK(category(j.constituents[3]) == kElectron);
  CHECK(category(j.constituents[4]) == kMuon);
  CHECK(categoryMask({"photon", "muon"}) == (kPhoton | kMuon));
  CHECK(categoryMask({"all"}) == kAll);
  CHECK(categoryMask({}) == 0);
  CHECK_THROWS_AS(categoryMask({"photons"}), std::invalid_argument);
}

TEST_CASE("jet cuts and overlap", "[jettools]") {
  auto j = mixedJet();
  j.jetIdBits = 1;
  CHECK(passes(j, JetCuts{299, 0.6, 1, 6}));
  CHECK_FALSE(passes(j, JetCuts{301, 0.6, 0, 0}));
  CHECK_FALSE(passes(j, JetCuts{0, 0.4, 0, 0}));
  CHECK_FALSE(passes(j, JetCuts{0, 1e9, 2, 0}));
  CHECK_FALSE(passes(j, JetCuts{0, 1e9, 0, 7}));
  CHECK(overlaps(j, {{0.5, 1.0 + 2 * M_PI - 2 * M_PI}}, 0.1));
  CHECK(overlaps(j, {{3.0, 0.0}, {0.9, 1.3}}, 0.8));   // dR = 0.5
  CHECK_FALSE(overlaps(j, {{0.9, 1.3}}, 0.4));
  CHECK(deltaR(0.0, 3.1, 0.0, -3.1) == Approx(2 * M_PI - 6.2));  // wraps through pi
}

TEST_CASE("constituent selection", "[jettools]") {
  auto j = mixedJet();
  selectConstituents(j, ConstituentCuts{0, true, {}, {}, 0});
  CHECK(j.constituents.size() == 4);  // charged: 211, 11, -13, -211
  j = mixedJet();
  selectConstituents(j, ConstituentCuts{2.0, false, {}, {}, 0});
  CHECK(j.constituents.size() == 4);  // drops pt 1.414 electron and pt 0.5 hadron
  j = mixedJet();
  selectConstituents(j, ConstituentCuts{0, false, {211}, {}, 0.5f});
  REQUIRE(j.constituents.size() == 1);
  CHECK(j.constituents[0].pdgId == 211);
  j = mixedJet();
  selectConstituents(j, ConstituentCuts{0, false, {}, {11, 13}, 0});
  CHECK(j.constituents.size() == 4);
}

TEST_CASE("scaling jets and constituents", "[jettools]") {
  auto j = mixedJet();
  scaleJet(j, 1.1);
  CHECK(j.pt == Approx(330));
  CHECK(j.mass == Approx(44));
  CHECK(j.eta == 0.5);
  CHECK_THROWS_AS(scaleJet(j, -0.1), std::invalid_argument);
  j = mixedJet();
  scaleConstituents(j, kPhoton | kNeutralHadron, 0.5);
  CHECK(j.constituents[0].px == 10);
  CHECK(j.constituents[1].py == 2.5);
  CHECK(j.constituents[1].E == Approx(0.5 * std::sqrt(29.0)));
  CHECK(j.constituents[2].px == 1.5);
  CHECK(j.constituents[3].px == 1);
}

TEST_CASE("constituent drop is reproducible and has the right rate", "[jettools]") {
  EECJet big;
  for (int i = 0; i < 4000; ++i)
    big.constituents.push_back(part(1 + i, 0, 0, i % 2, i % 2 ? 211 : 22));
  auto a = big, b = big;
  std::mt19937_64 r1(seed(1, 2, 3, 4, 0)), r2(seed(1, 2, 3, 4, 0));
  dropConstituents(a, kChargedHadron, 0.3, r1);
  dropConstituents(b, kChargedHadron, 0.3, r2);
  REQUIRE(a.constituents.size() == b.constituents.size());
  for (std::size_t i = 0; i < a.constituents.size(); ++i)
    CHECK(a.constituents[i].px == b.constituents[i].px);
  std::size_t photons = 0, charged = 0;
  for (const auto& p : a.constituents)
    (p.pdgId == 22 ? photons : charged)++;
  CHECK(photons == 2000);
  CHECK(charged > 1300);  // expectation 1400, sigma ~ 20
  CHECK(charged < 1500);
}

TEST_CASE("seed depends on every input", "[jettools]") {
  std::set<std::uint64_t> s{seed(1, 1, 1, 1, 0), seed(2, 1, 1, 1, 0), seed(1, 2, 1, 1, 0), seed(1, 1, 2, 1, 0),
                            seed(1, 1, 1, 2, 0), seed(1, 1, 1, 1, 1)};
  CHECK(s.size() == 6);
  CHECK(seed(5, 6, 7, 8, 9) == seed(5, 6, 7, 8, 9));
}

TEST_CASE("JER hybrid factor and matching", "[jettools]") {
  // Matched: 1 + (sf - 1)(pt - genPt)/pt = 1 + 0.2 * 10 / 100.
  CHECK(jerFactor(100, 90, 1.2, 0.1, 5.0) == Approx(1.02));
  // Unmatched: 1 + g * res * sqrt(sf^2 - 1).
  CHECK(jerFactor(100, -1, 1.25, 0.1, 2.0) == Approx(1 + 2.0 * 0.1 * 0.75));
  CHECK(jerFactor(100, -1, 0.9, 0.1, 2.0) == 1.0);  // sf < 1 unmatched: no smearing
  CHECK(jerFactor(100, -1, 2.0, 1.0, -5.0) == 0.0);  // clamped
  CHECK(jerMatch(0.3, 100, 80, 0.8, 0.1));
  CHECK_FALSE(jerMatch(0.41, 100, 99, 0.8, 0.1));
  CHECK_FALSE(jerMatch(0.1, 100, 69, 0.8, 0.1));
}
