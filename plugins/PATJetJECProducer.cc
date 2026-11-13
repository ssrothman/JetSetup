#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "correction.h"

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

// Copies of the pat::Jets with the JEC re-applied from the raw momentum: p4 = raw p4 * C, with C the correctionlib
// compound correction `key` (inputs JetA, JetEta, JetPt of the raw jet, Rho, and run for data keys), sorted by corrected pt
// if sortByPt (otherwise in source order, so that srcIdx is the position). User data:
// rawPt, jecFactor (C), jesUncTotal (relative uncertainty from uncKey at the corrected pt; absent if uncKey is
// empty), srcIdx (index in src) and jetId (bit 0 TIGHT, bit 1 TIGHTLEPVETO of PFJetIDSelectionFunctor; 0 if
// jetIdVersion is empty).
// A stream module: PFJetIDSelectionFunctor::operator() is not const.
class PATJetJECProducer : public edm::stream::EDProducer<> {
public:
  explicit PATJetJECProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        rho_(consumes(cfg.getParameter<edm::InputTag>("rho"))),
        put_(produces()),
        sortByPt_(cfg.getParameter<bool>("sortByPt")) {
    auto cset = correction::CorrectionSet::from_file(cfg.getParameter<std::string>("correctionFile"));
    jec_ = cset->compound().at(cfg.getParameter<std::string>("key"));
    for (const auto& in : jec_->inputs())
      if (in.name() != "JetA" && in.name() != "JetEta" && in.name() != "JetPt" && in.name() != "Rho" &&
          in.name() != "run")
        throw cms::Exception("Configuration") << "PATJetJECProducer: unsupported JEC input " << in.name();
    const auto unc = cfg.getParameter<std::string>("uncKey");
    if (!unc.empty())
      unc_ = cset->at(unc);
    const auto version = cfg.getParameter<std::string>("jetIdVersion");
    if (!version.empty()) {
      edm::ParameterSet p;
      p.addParameter<std::string>("version", version);
      p.addParameter<std::string>("quality", "TIGHT");
      tight_.emplace(p);
      p.addParameter<std::string>("quality", "TIGHTLEPVETO");
      tightLepVeto_.emplace(p);
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("slimmedJetsPuppi"));
    desc.add<edm::InputTag>("rho", edm::InputTag("fixedGridRhoFastjetAll"));
    desc.add<std::string>("correctionFile");
    desc.add<std::string>("key")->setComment("compound correction, e.g. ..._MC_L1L2L3Res_AK4PFPuppi");
    desc.add<std::string>("uncKey")->setComment("JES uncertainty, e.g. ..._MC_Total_AK4PFPuppi; empty: none");
    desc.add<std::string>("jetIdVersion")->setComment("PFJetIDSelectionFunctor version; empty: jetId 0");
    desc.add<bool>("sortByPt")->setComment("sort by corrected pt; false: source order");
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto& jets = event.get(src_);
    const double rho = event.get(rho_);
    auto out = std::make_unique<std::vector<pat::Jet>>();
    for (std::size_t i = 0; i < jets.size(); ++i) {
      const auto& j = jets[i];
      const auto raw = j.jecSetsAvailable() ? j.correctedP4("Uncorrected") : j.p4();
      std::vector<correction::Variable::Type> args;
      for (const auto& in : jec_->inputs()) {
        const auto& n = in.name();
        if (n == "run")
          args.emplace_back(double(event.id().run()));
        else
          args.emplace_back(n == "JetA" ? j.jetArea() : n == "JetEta" ? raw.eta() : n == "JetPt" ? raw.pt() : rho);
      }
      const double c = jec_->evaluate(args);
      // Not setP4: scaleEnergy also inserts a first level "Uncorrected" = 1/c, so correctedP4("Uncorrected") and the
      // energy fractions (jecFactor(0) * energy) still give the raw jet. The other stored levels stay MiniAOD's.
      pat::Jet cj(j.jecSetsAvailable() ? j.correctedJet("Uncorrected") : j);
      cj.scaleEnergy(c, "Uncorrected");
      cj.addUserFloat("rawPt", raw.pt());
      cj.addUserFloat("jecFactor", c);
      if (unc_)
        cj.addUserFloat("jesUncTotal", (*unc_)->evaluate({cj.eta(), cj.pt()}));
      cj.addUserInt("srcIdx", i);
      int id = 0;
      if (tight_) {
        pat::strbitset r1 = tight_->getBitTemplate(), r2 = tightLepVeto_->getBitTemplate();
        id = ((*tight_)(j, r1) ? 1 : 0) | ((*tightLepVeto_)(j, r2) ? 2 : 0);
      }
      cj.addUserInt("jetId", id);
      out->push_back(std::move(cj));
    }
    if (sortByPt_)
      std::stable_sort(out->begin(), out->end(), [](const auto& a, const auto& b) { return a.pt() > b.pt(); });
    event.put(put_, std::move(out));
  }

private:
  const edm::EDGetTokenT<edm::View<pat::Jet>> src_;
  const edm::EDGetTokenT<double> rho_;
  const edm::EDPutTokenT<std::vector<pat::Jet>> put_;
  const bool sortByPt_;
  correction::CompoundCorrection::Ref jec_;
  std::optional<correction::Correction::Ref> unc_;
  std::optional<PFJetIDSelectionFunctor> tight_, tightLepVeto_;
};

DEFINE_FWK_MODULE(PATJetJECProducer);
