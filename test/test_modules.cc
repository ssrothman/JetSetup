#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SRothman/DataFormats/interface/EECJet.h"
#include "SRothman/JetSetup/interface/JetTools.h"

#include "correction.h"

#include <cmath>
#include <string>

using srothman::EECJet;
using srothman::EECJetCollection;

namespace {

  const std::string kHeader = "from FWCore.TestProcessor.TestProcess import *\nprocess = TestProcess()\n";
  const std::string kJME =
      "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer22EE/fatJet_jerc.json.gz";

  srothman::EECParticle part(double px, double py, double pz, int q, int pdg) {
    return {px, py, pz, std::sqrt(px * px + py * py + pz * pz), static_cast<std::int8_t>(q), pdg, 1.0f, 0};
  }

  EECJetCollection jets() {
    EECJetCollection v(3);
    v[0].pt = 450;
    v[0].eta = 0.3;
    v[0].phi = 0.2;
    v[0].mass = 90;
    v[0].srcIdx = 0;
    v[0].jetIdBits = 3;
    v[0].constituents = {part(100, 10, 30, 1, 211), part(50, 5, 20, 0, 22), part(30, 0, 5, 0, 130)};
    v[1].pt = 250;
    v[1].eta = -1.9;
    v[1].phi = 2.8;
    v[1].mass = 30;
    v[1].srcIdx = 1;
    v[1].jetIdBits = 1;
    v[1].constituents = {part(-40, 10, -100, -1, -211)};
    v[2].pt = 120;
    v[2].eta = 2.2;
    v[2].phi = -1.0;
    v[2].srcIdx = 4;
    return v;
  }

  double col(const nanoaod::FlatTable& t, const std::string& name, unsigned row) {
    const int c = t.columnIndex(name);
    REQUIRE(c >= 0);
    return t.getAnyValue(row, c);
  }

}  // namespace

TEST_CASE("EECJetSelector applies cuts and lepton cleaning", "[modules]") {
  edm::test::TestProcessor::Config config{
      kHeader +
      "process.toTest = cms.EDProducer('EECJetSelector', src = cms.InputTag('jets'), ptMin = cms.double(200), "
      "absEtaMax = cms.double(2.0), idMask = cms.uint32(1), minConstituents = cms.uint32(1), vetoSrcs = "
      "cms.VInputTag('leptons'), vetoDR = cms.double(0.8))\nprocess.moduleToTest(process.toTest)\n"};
  auto putJ = config.produces<EECJetCollection>("jets");
  auto putL = config.produces<std::vector<reco::LeafCandidate>>("leptons");
  edm::test::TestProcessor tester(config);
  {
    auto ev = tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(jets())),
                          std::make_pair(putL, std::make_unique<std::vector<reco::LeafCandidate>>()));
    const auto& out = *ev.get<EECJetCollection>();
    REQUIRE(out.size() == 2);  // jet 2 fails pt
    CHECK(out[0].srcIdx == 0);
    CHECK(out[1].srcIdx == 1);
  }
  // A lepton at dR ~ 0.5 from jet 0 removes it.
  std::vector<reco::LeafCandidate> lep(1);
  lep[0].setP4(reco::LeafCandidate::PolarLorentzVector(50, 0.6, 0.6, 0.1));
  auto ev = tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(jets())),
                        std::make_pair(putL, std::make_unique<std::vector<reco::LeafCandidate>>(lep)));
  const auto& out = *ev.get<EECJetCollection>();
  REQUIRE(out.size() == 1);
  CHECK(out[0].srcIdx == 1);
}

TEST_CASE("EECConstituentSelector keeps jets and filters constituents", "[modules]") {
  edm::test::TestProcessor::Config config{
      kHeader +
      "process.toTest = cms.EDProducer('EECConstituentSelector', src = cms.InputTag('jets'), ptMin = "
      "cms.double(1.0), chargedOnly = cms.bool(True), keepAbsPdgIds = cms.vint32(), dropAbsPdgIds = cms.vint32(), "
      "minPuppiWeight = cms.double(0))\nprocess.moduleToTest(process.toTest)\n"};
  auto put = config.produces<EECJetCollection>("jets");
  edm::test::TestProcessor tester(config);
  auto ev = tester.test(std::make_pair(put, std::make_unique<EECJetCollection>(jets())));
  const auto& out = *ev.get<EECJetCollection>();
  REQUIRE(out.size() == 3);
  REQUIRE(out[0].constituents.size() == 1);
  CHECK(out[0].constituents[0].pdgId == 211);
  CHECK(out[0].pt == 450);
  CHECK(out[1].constituents.size() == 1);
  CHECK(out[2].constituents.empty());
}

TEST_CASE("EECJetVariationProducer jes matches correctionlib", "[modules][variation]") {
  auto cset = correction::CorrectionSet::from_file(kJME);
  auto unc = cset->at("Summer22EE_22Sep2023_V2_MC_Total_AK8PFPuppi");
  for (std::string dir : {"up", "down"}) {
    edm::test::TestProcessor::Config config{kHeader +
                                            "process.toTest = cms.EDProducer('EECJetVariationProducer', src = "
                                            "cms.InputTag('jets'), type = cms.string('jes'), direction = cms.string('" +
                                            dir + "'))\nprocess.moduleToTest(process.toTest)\n"};
    auto put = config.produces<EECJetCollection>("jets");
    edm::test::TestProcessor tester(config);
    const auto in = jets();
    auto ev = tester.test(std::make_pair(put, std::make_unique<EECJetCollection>(in)));
    const auto& out = *ev.get<EECJetCollection>();
    REQUIRE(out.size() == in.size());
    for (std::size_t i = 0; i < in.size(); ++i) {
      const double u = unc->evaluate({in[i].eta, in[i].pt});
      CHECK(u > 0.001);
      const double f = dir == "up" ? 1 + u : 1 - u;
      CHECK(out[i].pt == Approx(in[i].pt * f).epsilon(1e-12));
      CHECK(out[i].mass == Approx(in[i].mass * f).epsilon(1e-12));
      CHECK(out[i].constituents.size() == in[i].constituents.size());
      if (!out[i].constituents.empty())
        CHECK(out[i].constituents[0].px == in[i].constituents[0].px);
    }
  }
}

TEST_CASE("EECJetVariationProducer jer: matched jets use the hybrid formula", "[modules][variation]") {
  auto cset = correction::CorrectionSet::from_file(kJME);
  auto sfC = cset->at("Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK8PFPuppi");
  auto resC = cset->at("Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK8PFPuppi");
  const double rho = 20.0;
  for (std::string dir : {"nominal", "up"}) {
    edm::test::TestProcessor::Config config{kHeader +
                                            "process.toTest = cms.EDProducer('EECJetVariationProducer', src = "
                                            "cms.InputTag('jets'), type = cms.string('jer'), direction = cms.string('" +
                                            dir +
                                            "'), rho = cms.InputTag('rho'), genJets = cms.InputTag('gen'), jetR = "
                                            "cms.double(0.8))\nprocess.moduleToTest(process.toTest)\n"};
    auto putJ = config.produces<EECJetCollection>("jets");
    auto putR = config.produces<double>("rho");
    auto putG = config.produces<std::vector<reco::LeafCandidate>>("gen");
    edm::test::TestProcessor tester(config);
    const auto in = jets();
    // Gen jet close to jet 0 (dR 0.05, pt 440); none near jets 1 and 2.
    std::vector<reco::LeafCandidate> gen(1);
    gen[0].setP4(reco::LeafCandidate::PolarLorentzVector(440, 0.33, 0.24, 85));
    auto ev = tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(in)),
                          std::make_pair(putR, std::make_unique<double>(rho)),
                          std::make_pair(putG, std::make_unique<std::vector<reco::LeafCandidate>>(gen)));
    const auto& out = *ev.get<EECJetCollection>();
    const std::string sys = dir == "nominal" ? "nom" : dir;
    const double sf = sfC->evaluate({in[0].eta, in[0].pt, sys});
    REQUIRE(resC->evaluate({in[0].eta, in[0].pt, rho}) * 3 * 450 > 10);  // the match criterion holds
    CHECK(out[0].pt == Approx(in[0].pt * (1 + (sf - 1) * (450 - 440.0) / 450)).epsilon(1e-9));
    CHECK(sf > 1.0);  // otherwise the matched formula would not move the jet
  }
}

TEST_CASE("EECJetVariationProducer constScale and constDrop", "[modules][variation]") {
  {
    edm::test::TestProcessor::Config config{
        kHeader +
        "process.toTest = cms.EDProducer('EECJetVariationProducer', src = cms.InputTag('jets'), type = "
        "cms.string('constScale'), categories = cms.vstring('photon'), scale = cms.double(1.05))\n"
        "process.moduleToTest(process.toTest)\n"};
    auto put = config.produces<EECJetCollection>("jets");
    edm::test::TestProcessor tester(config);
    auto ev = tester.test(std::make_pair(put, std::make_unique<EECJetCollection>(jets())));
    const auto& out = *ev.get<EECJetCollection>();
    CHECK(out[0].constituents[0].px == 100);
    CHECK(out[0].constituents[1].px == Approx(52.5));
    CHECK(out[0].pt == 450);
  }
  {
    edm::test::TestProcessor::Config config{
        kHeader +
        "process.toTest = cms.EDProducer('EECJetVariationProducer', src = cms.InputTag('jets'), type = "
        "cms.string('constDrop'), categories = cms.vstring('all'), dropProbability = cms.double(1.0))\n"
        "process.moduleToTest(process.toTest)\n"};
    auto put = config.produces<EECJetCollection>("jets");
    edm::test::TestProcessor tester(config);
    auto ev = tester.test(std::make_pair(put, std::make_unique<EECJetCollection>(jets())));
    for (const auto& j : *ev.get<EECJetCollection>())
      CHECK(j.constituents.empty());
  }
  edm::test::TestProcessor::Config bad{kHeader +
                                       "process.toTest = cms.EDProducer('EECJetVariationProducer', src = "
                                       "cms.InputTag('jets'), type = cms.string('constDrop'), categories = "
                                       "cms.vstring('gluon'))\nprocess.moduleToTest(process.toTest)\n"};
  bad.produces<EECJetCollection>("jets");
  CHECK_THROWS_AS(edm::test::TestProcessor(bad), cms::Exception);
}

namespace {
  // Unmatched-jet JER smearing of jets() for one seedSalt.
  EECJetCollection smearUnmatched(unsigned salt) {
    edm::test::TestProcessor::Config config{
        kHeader +
        "process.toTest = cms.EDProducer('EECJetVariationProducer', src = cms.InputTag('jets'), type = "
        "cms.string('jer'), direction = cms.string('nominal'), rho = cms.InputTag('rho'), genJets = "
        "cms.InputTag('gen'), seedSalt = cms.uint32(" +
        std::to_string(salt) + "))\nprocess.moduleToTest(process.toTest)\n"};
    auto putJ = config.produces<EECJetCollection>("jets");
    auto putR = config.produces<double>("rho");
    auto putG = config.produces<std::vector<reco::LeafCandidate>>("gen");
    edm::test::TestProcessor tester(config);
    auto ev = tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(jets())),
                          std::make_pair(putR, std::make_unique<double>(15.0)),
                          std::make_pair(putG, std::make_unique<std::vector<reco::LeafCandidate>>()));
    return *ev.get<EECJetCollection>();
  }
}  // namespace

TEST_CASE("EECJetVariationProducer jer smearing depends only on event, jet and seedSalt", "[modules][variation]") {
  const auto a = smearUnmatched(0), b = smearUnmatched(0), c = smearUnmatched(7);
  const auto in = jets();
  for (std::size_t i = 0; i < in.size(); ++i) {
    CHECK(a[i].pt == b[i].pt);
    CHECK(a[i].pt != c[i].pt);
  }
  // With SF > 1 an unmatched jet is smeared.
  CHECK(a[1].pt != in[1].pt);
}

TEST_CASE("EECJetTableProducer", "[modules][table]") {
  edm::test::TestProcessor::Config config{kHeader +
                                          "process.toTest = cms.EDProducer('EECJetTableProducer', src = "
                                          "cms.InputTag('jets'), name = cms.string('FatJetEEC'), writeConstituents = "
                                          "cms.bool(True))\nprocess.moduleToTest(process.toTest)\n"};
  auto put = config.produces<EECJetCollection>("jets");
  edm::test::TestProcessor tester(config);
  auto ev = tester.test(std::make_pair(put, std::make_unique<EECJetCollection>(jets())));
  const auto& jt = *ev.get<nanoaod::FlatTable>();
  const auto& ct = *ev.get<nanoaod::FlatTable>("constituents");
  CHECK(jt.name() == "FatJetEEC");
  REQUIRE(jt.nRows() == 3);
  CHECK(col(jt, "pt", 0) == 450);
  CHECK(col(jt, "srcIdx", 2) == 4);
  CHECK(col(jt, "jetId", 0) == 3);
  CHECK(col(jt, "nConstituents", 0) == 3);
  CHECK(ct.name() == "FatJetEECConstituents");
  REQUIRE(ct.nRows() == 4);
  CHECK(col(ct, "iJet", 3) == 1);
  CHECK(col(ct, "pt", 0) == Approx(std::hypot(100, 10)).epsilon(1e-6));
  CHECK(col(ct, "eta", 3) == Approx(std::asinh(-100 / std::hypot(40, 10))).epsilon(1e-6));
  CHECK(col(ct, "phi", 3) == Approx(std::atan2(10, -40)).epsilon(1e-6));
  CHECK(col(ct, "charge", 3) == -1);
  CHECK(col(ct, "pdgId", 1) == 22);
}
