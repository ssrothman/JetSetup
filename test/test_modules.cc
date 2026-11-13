#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
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
      "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/fatJet_jerc.json.gz";

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
  auto unc = cset->at("Summer22EE_22Sep2023_V4_MC_Total_AK8PFPuppi");
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
  auto sfC = cset->at("Summer22EE_22Sep2023_JRV2_MC_ScaleFactor_AK8PFPuppi");
  auto uncC = cset->at("Summer22EE_22Sep2023_JRV2_MC_SFUncertainty_AK8PFPuppi");
  auto resC = cset->at("Summer22EE_22Sep2023_JRV2_MC_PtResolution_AK8PFPuppi");
  const double rho = 20.0;
  for (std::string dir : {"nominal", "up", "down"}) {
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
    const double u = uncC->evaluate({in[0].eta, in[0].pt});
    REQUIRE(u > 0.001);
    const double sf = sfC->evaluate({in[0].eta, in[0].pt}) * (dir == "up" ? 1 + u : dir == "down" ? 1 - u : 1.0);
    REQUIRE(resC->evaluate({in[0].eta, in[0].pt, rho}) * 3 * 450 > 10);  // the match criterion holds
    CHECK(out[0].pt == Approx(in[0].pt * (1 + (sf - 1) * (450 - 440.0) / 450)).epsilon(1e-9));
    CHECK(sf > 1.0);  // otherwise the matched formula would not move the jet
    for (std::size_t i = 0; i < in.size(); ++i)
      CHECK(out[i].jerFactor == Approx(out[i].pt / in[i].pt).epsilon(1e-12));
    CHECK(out[0].jerFactor != 1.0);
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
  CHECK(a[1].jerFactor == Approx(a[1].pt / in[1].pt).epsilon(1e-12));
}

TEST_CASE("EECJetVariationProducer: jerFactor accumulates and other types leave it", "[modules][variation]") {
  auto in = jets();
  in[0].jerFactor = 1.5;  // as after an earlier JER step
  edm::test::TestProcessor::Config config{
      kHeader +
      "process.toTest = cms.EDProducer('EECJetVariationProducer', src = cms.InputTag('jets'), type = "
      "cms.string('jer'), direction = cms.string('nominal'), rho = cms.InputTag('rho'), genJets = "
      "cms.InputTag('gen'))\nprocess.moduleToTest(process.toTest)\n"};
  auto putJ = config.produces<EECJetCollection>("jets");
  auto putR = config.produces<double>("rho");
  auto putG = config.produces<std::vector<reco::LeafCandidate>>("gen");
  edm::test::TestProcessor tester(config);
  auto ev = tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(in)),
                        std::make_pair(putR, std::make_unique<double>(15.0)),
                        std::make_pair(putG, std::make_unique<std::vector<reco::LeafCandidate>>()));
  const auto& out = *ev.get<EECJetCollection>();
  CHECK(out[0].jerFactor == Approx(1.5 * out[0].pt / in[0].pt).epsilon(1e-12));
  edm::test::TestProcessor::Config jes{
      kHeader +
      "process.toTest = cms.EDProducer('EECJetVariationProducer', src = cms.InputTag('jets'), type = "
      "cms.string('jes'), direction = cms.string('up'))\nprocess.moduleToTest(process.toTest)\n"};
  auto putJ2 = jes.produces<EECJetCollection>("jets");
  edm::test::TestProcessor tester2(jes);
  auto ev2 = tester2.test(std::make_pair(putJ2, std::make_unique<EECJetCollection>(in)));
  CHECK(ev2.get<EECJetCollection>()->at(0).jerFactor == 1.5);
  CHECK(ev2.get<EECJetCollection>()->at(0).pt != in[0].pt);
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

TEST_CASE("EECJetExtensionTableProducer: columns from the source pat::Jet at srcIdx", "[modules][table]") {
  edm::test::TestProcessor::Config config{
      kHeader +
      "process.toTest = cms.EDProducer('EECJetExtensionTableProducer', src = cms.InputTag('jets'), jets = "
      "cms.InputTag('pat'), name = cms.string('RecoJet'), variables = cms.VPSet("
      "cms.PSet(name = cms.string('jecPt'), expr = cms.string('pt'), type = cms.string('float'), doc = "
      "cms.string(''), precision = cms.int32(-1)), "
      "cms.PSet(name = cms.string('area'), expr = cms.string('jetArea()'), type = cms.string('float'), doc = "
      "cms.string(''), precision = cms.int32(-1)), "
      "cms.PSet(name = cms.string('isCentral'), expr = cms.string('abs(eta) < 1'), type = cms.string('bool'), doc = "
      "cms.string(''), precision = cms.int32(-1)), "
      "cms.PSet(name = cms.string('idx'), expr = cms.string('charge + 7'), type = cms.string('int16'), doc = "
      "cms.string(''), precision = cms.int32(-1))))\nprocess.moduleToTest(process.toTest)\n"};
  auto putJ = config.produces<EECJetCollection>("jets");
  auto putP = config.produces<std::vector<pat::Jet>>("pat");
  edm::test::TestProcessor tester(config);
  std::vector<pat::Jet> pats(5);
  for (std::size_t i = 0; i < pats.size(); ++i) {
    pats[i].setP4(reco::LeafCandidate::PolarLorentzVector(100 + 10 * i, i < 2 ? 0.5 : 1.5, 0, 10));
    pats[i].setJetArea(0.4 + 0.01 * i);
  }
  auto in = jets();  // srcIdx 0, 1, 4
  in[1].jerFactor = 1.07;
  auto ev = tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(in)),
                        std::make_pair(putP, std::make_unique<std::vector<pat::Jet>>(pats)));
  const auto& t = *ev.get<nanoaod::FlatTable>();
  CHECK(t.name() == "RecoJet");
  CHECK(t.extension());
  REQUIRE(t.nRows() == 3);
  for (unsigned r = 0; r < 3; ++r) {
    const auto& p = pats[in[r].srcIdx];
    CHECK(col(t, "jecPt", r) == Approx(p.pt()).epsilon(1e-6));
    CHECK(col(t, "area", r) == Approx(p.jetArea()).epsilon(1e-6));
    CHECK(col(t, "isCentral", r) == (std::abs(p.eta()) < 1 ? 1 : 0));
    CHECK(col(t, "idx", r) == 7);
  }
  CHECK(col(t, "jecPt", 2) == Approx(140).epsilon(1e-6));  // srcIdx 4, not row 2
  CHECK(col(t, "jerFactor", 0) == 1.0);
  CHECK(col(t, "jerFactor", 1) == Approx(1.07).epsilon(1e-6));
  // srcIdx outside the source collection is an error.
  std::vector<pat::Jet> few(pats.begin(), pats.begin() + 2);
  CHECK_THROWS_AS(tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(jets())),
                              std::make_pair(putP, std::make_unique<std::vector<pat::Jet>>(few))),
                  cms::Exception);
}

TEST_CASE("EECJetExtensionTableProducer: configuration errors", "[modules][table]") {
  for (const std::string& v :
       {"name = cms.string('a'), expr = cms.string('pt'), type = cms.string('double')",
        "name = cms.string('jerFactor'), expr = cms.string('pt'), type = cms.string('float')",
        "name = cms.string('a'), expr = cms.string('noSuchMethod()'), type = cms.string('float')"}) {
    edm::test::TestProcessor::Config config{
        kHeader +
        "process.toTest = cms.EDProducer('EECJetExtensionTableProducer', src = cms.InputTag('jets'), jets = "
        "cms.InputTag('pat'), name = cms.string('RecoJet'), variables = cms.VPSet(cms.PSet(" +
        v + ", doc = cms.string(''), precision = cms.int32(-1))))\nprocess.moduleToTest(process.toTest)\n"};
    config.produces<EECJetCollection>("jets");
    config.produces<std::vector<pat::Jet>>("pat");
    CHECK_THROWS(edm::test::TestProcessor(config));
  }
}

TEST_CASE("PATJetJECProducer: re-applied compound JEC, JES uncertainty, pt order", "[modules][jec]") {
  const std::string file =
      "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jet_jerc.json.gz";
  const std::string key = "Summer22EE_22Sep2023_V4_MC_L1L2L3Res_AK4PFPuppi",
                    unc = "Summer22EE_22Sep2023_V4_MC_Total_AK4PFPuppi";
  edm::test::TestProcessor::Config config{
      kHeader +
      "process.toTest = cms.EDProducer('PATJetJECProducer', src = cms.InputTag('pat'), rho = cms.InputTag('rho'), "
      "correctionFile = cms.string('" +
      file + "'), key = cms.string('" + key + "'), uncKey = cms.string('" + unc +
      "'), jetIdVersion = cms.string(''))\nprocess.moduleToTest(process.toTest)\n"};
  auto putP = config.produces<std::vector<pat::Jet>>("pat");
  auto putR = config.produces<double>("rho");
  edm::test::TestProcessor tester(config);
  auto cs = correction::CorrectionSet::from_file(file);
  // Raw jets (no JEC levels stored): 30 GeV at eta 0.5 and 60 GeV at eta -2.0; the corrected order is by pt.
  std::vector<pat::Jet> raw(2);
  raw[0].setP4(reco::LeafCandidate::PolarLorentzVector(30, 0.5, 0.1, 5));
  raw[0].setJetArea(0.5);
  raw[1].setP4(reco::LeafCandidate::PolarLorentzVector(60, -2.0, 1.0, 8));
  raw[1].setJetArea(0.45);
  const double rho = 25;
  auto ev = tester.test(std::make_pair(putP, std::make_unique<std::vector<pat::Jet>>(raw)),
                        std::make_pair(putR, std::make_unique<double>(rho)));
  const auto& out = *ev.get<std::vector<pat::Jet>>();
  REQUIRE(out.size() == 2);
  CHECK(out[0].userInt("srcIdx") == 1);  // 60 GeV raw jet first
  for (const auto& j : out) {
    const auto& r = raw[j.userInt("srcIdx")];
    const double c = cs->compound().at(key)->evaluate({r.jetArea(), r.eta(), r.pt(), rho});
    REQUIRE(std::abs(c - 1) > 1e-3);  // the correction is not trivial
    CHECK(j.pt() == Approx(r.pt() * c).epsilon(1e-12));
    CHECK(j.mass() == Approx(r.mass() * c).epsilon(1e-9));
    CHECK(j.eta() == Approx(r.eta()).epsilon(1e-12));
    CHECK(j.userFloat("rawPt") == Approx(r.pt()).epsilon(1e-6));
    CHECK(j.userFloat("jecFactor") == Approx(c).epsilon(1e-6));
    CHECK(j.userFloat("jesUncTotal") == Approx(cs->at(unc)->evaluate({j.eta(), j.pt()})).epsilon(1e-6));
    CHECK(j.userInt("jetId") == 0);
  }
  CHECK(out[0].pt() > out[1].pt());
}

TEST_CASE("PATJetJECProducer: data keys take the run number", "[modules][jec]") {
  const std::string file =
      "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jet_jerc.json.gz";
  const std::string key = "Summer22EE_22Sep2023_V4_DATA_L1L2L3Res_AK4PFPuppi";
  edm::test::TestProcessor::Config config{
      kHeader +
      "process.toTest = cms.EDProducer('PATJetJECProducer', src = cms.InputTag('pat'), rho = cms.InputTag('rho'), "
      "correctionFile = cms.string('" +
      file + "'), key = cms.string('" + key +
      "'), uncKey = cms.string(''), jetIdVersion = cms.string(''))\nprocess.moduleToTest(process.toTest)\n"};
  auto putP = config.produces<std::vector<pat::Jet>>("pat");
  auto putR = config.produces<double>("rho");
  edm::test::TestProcessor tester(config);
  auto cs = correction::CorrectionSet::from_file(file);
  std::vector<pat::Jet> raw(1);
  raw[0].setP4(reco::LeafCandidate::PolarLorentzVector(40, 1.0, 0.1, 5));
  raw[0].setJetArea(0.5);
  for (unsigned run : {360000u, 362500u}) {
    tester.setRunNumber(run);
    auto ev = tester.test(std::make_pair(putP, std::make_unique<std::vector<pat::Jet>>(raw)),
                          std::make_pair(putR, std::make_unique<double>(20.0)));
    const auto& j = ev.get<std::vector<pat::Jet>>()->at(0);
    CHECK(j.pt() == Approx(40 * cs->compound().at(key)->evaluate({0.5, 1.0, 40.0, 20.0, double(run)})).epsilon(1e-12));
    CHECK(!j.hasUserFloat("jesUncTotal"));
  }
}

TEST_CASE("EECJetExtensionTableProducer: taggedVariables from the re-tagged copy at srcIdx", "[modules][table]") {
  edm::test::TestProcessor::Config config{
      kHeader +
      "process.toTest = cms.EDProducer('EECJetExtensionTableProducer', src = cms.InputTag('jets'), jets = "
      "cms.InputTag('pat'), name = cms.string('RecoJet'), variables = cms.VPSet(cms.PSet(name = cms.string('a'), expr "
      "= "
      "cms.string('jetArea()'), type = cms.string('float'), doc = cms.string(''), precision = cms.int32(-1))), "
      "taggedJets = cms.InputTag('tagged'), taggedVariables = cms.VPSet(cms.PSet(name = cms.string('ta'), expr = "
      "cms.string('jetArea()'), type = cms.string('float'), doc = cms.string(''), precision = "
      "cms.int32(-1))))\nprocess.moduleToTest(process.toTest)\n"};
  auto putJ = config.produces<EECJetCollection>("jets");
  auto putP = config.produces<std::vector<pat::Jet>>("pat");
  auto putT = config.produces<std::vector<pat::Jet>>("tagged");
  edm::test::TestProcessor tester(config);
  std::vector<pat::Jet> pats(5), tagged(5);
  for (std::size_t i = 0; i < pats.size(); ++i) {
    pats[i].setP4(reco::LeafCandidate::PolarLorentzVector(100 + 10 * i, 0.1 * i, 0.2 * i, 10));
    pats[i].setJetArea(0.4 + 0.01 * i);
    tagged[i] = pats[i];
    tagged[i].setP4(pats[i].p4() * 1.05);  // other JECs, same direction
    tagged[i].setJetArea(1.0 + i);
  }
  auto ev = tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(jets())),
                        std::make_pair(putP, std::make_unique<std::vector<pat::Jet>>(pats)),
                        std::make_pair(putT, std::make_unique<std::vector<pat::Jet>>(tagged)));
  const auto& t = *ev.get<nanoaod::FlatTable>();
  const auto in = jets();  // srcIdx 0, 1, 4
  for (unsigned r = 0; r < 3; ++r) {
    CHECK(col(t, "a", r) == Approx(pats[in[r].srcIdx].jetArea()).epsilon(1e-6));
    CHECK(col(t, "ta", r) == Approx(1.0 + in[r].srcIdx).epsilon(1e-6));
  }
  // A re-tagged collection in another order, or of another size, is an error.
  std::vector<pat::Jet> swapped = tagged;
  std::swap(swapped[0], swapped[1]);
  CHECK_THROWS_AS(tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(jets())),
                              std::make_pair(putP, std::make_unique<std::vector<pat::Jet>>(pats)),
                              std::make_pair(putT, std::make_unique<std::vector<pat::Jet>>(swapped))),
                  cms::Exception);
  CHECK_THROWS_AS(
      tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(jets())),
                  std::make_pair(putP, std::make_unique<std::vector<pat::Jet>>(pats)),
                  std::make_pair(putT, std::make_unique<std::vector<pat::Jet>>(tagged.begin(), tagged.begin() + 4))),
      cms::Exception);
}

TEST_CASE("EECJetOverlapProducer: per-channel flags and veto columns, jets kept", "[modules][veto]") {
  edm::test::TestProcessor::Config config{
      kHeader +
      "process.toTest = cms.EDProducer('EECJetOverlapProducer', src = cms.InputTag('jets'), channels = "
      "cms.vstring('Zmumu', 'GammaJets'), vetoSrcs = cms.VInputTag('mu', 'pho'), vetoDR = cms.double(0.8), name = "
      "cms.string('RecoJet'))\nprocess.moduleToTest(process.toTest)\n"};
  auto putJ = config.produces<EECJetCollection>("jets");
  auto putM = config.produces<std::vector<reco::LeafCandidate>>("mu");
  auto putP = config.produces<std::vector<reco::LeafCandidate>>("pho");
  edm::test::TestProcessor tester(config);
  // jets(): (0.3, 0.2), (-1.9, 2.8), (2.2, -1.0). A muon at dR 0.5 from jet 0; a photon at dR 0.79 from jet 1 (wrapping
  // in phi) and one at dR 0.81 from jet 2.
  std::vector<reco::LeafCandidate> mu(1), pho(2);
  mu[0].setP4(reco::LeafCandidate::PolarLorentzVector(50, 0.6, 0.6, 0.1));
  pho[0].setP4(reco::LeafCandidate::PolarLorentzVector(50, -1.9, 2.8 + 0.79 - 2 * M_PI, 0));
  pho[1].setP4(reco::LeafCandidate::PolarLorentzVector(50, 2.2 + 0.81, -1.0, 0));
  auto ev = tester.test(std::make_pair(putJ, std::make_unique<EECJetCollection>(jets())),
                        std::make_pair(putM, std::make_unique<std::vector<reco::LeafCandidate>>(mu)),
                        std::make_pair(putP, std::make_unique<std::vector<reco::LeafCandidate>>(pho)));
  CHECK(*ev.get<std::vector<bool>>("Zmumu") == std::vector<bool>{true, false, false});
  CHECK(*ev.get<std::vector<bool>>("GammaJets") == std::vector<bool>{false, true, false});
  const auto& t = *ev.get<nanoaod::FlatTable>();
  CHECK(t.name() == "RecoJet");
  CHECK(t.extension());
  REQUIRE(t.nRows() == 3);
  for (unsigned r = 0; r < 3; ++r) {
    CHECK(col(t, "vetoZmumu", r) == (r == 0 ? 1 : 0));
    CHECK(col(t, "vetoGammaJets", r) == (r == 1 ? 1 : 0));
  }
}
