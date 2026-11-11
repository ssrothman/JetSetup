#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "SRothman/JetSetup/plugins/ConstituentFill.h"

#include <memory>
#include <optional>

// pat::Jet -> EECJet. Constituents are the jet daughters, optionally scaled by their PUPPI weight.
// Throws for a jet above ptMin without daughters.
// jetIdBits: bit 0 = TIGHT, bit 1 = TIGHTLEPVETO of PFJetIDSelectionFunctor(jetIdVersion); 0 if jetIdVersion is empty.
class PATJetToEECJetProducer : public edm::stream::EDProducer<> {
public:
  explicit PATJetToEECJetProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        put_(produces()),
        applyPuppi_(cfg.getParameter<bool>("applyPuppiWeights")),
        ptMin_(cfg.getParameter<double>("ptMin")) {
    const auto vm = cfg.getParameter<edm::InputTag>("puppiValueMap");
    if (!vm.label().empty())
      vmToken_ = consumes(vm);
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
    desc.add<edm::InputTag>("src", edm::InputTag("slimmedJetsAK8"));
    desc.add<bool>("applyPuppiWeights", true);
    desc.add<edm::InputTag>("puppiValueMap", edm::InputTag(""))
        ->setComment("empty: use pat::PackedCandidate::puppiWeight()");
    desc.add<std::string>("jetIdVersion", "RUN3PUPPIruns2022FGruns2023CD")
        ->setComment("PFJetIDSelectionFunctor version; empty: no jet ID bits");
    desc.add<double>("ptMin", 0.0)->setComment("corrected pt threshold applied before conversion");
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto& jets = event.get(src_);
    const edm::ValueMap<float>* vm = vmToken_ ? &event.get(*vmToken_) : nullptr;
    srothman::EECJetCollection out;
    for (std::size_t i = 0; i < jets.size(); ++i) {
      const auto& jet = jets[i];
      if (jet.pt() < ptMin_)
        continue;
      srothman::EECJet j;
      j.pt = jet.pt();
      j.eta = jet.eta();
      j.phi = jet.phi();
      j.mass = jet.mass();
      j.rawPt = jet.jecSetsAvailable() ? jet.correctedP4("Uncorrected").pt() : jet.pt();
      j.srcIdx = i;
      if (tight_) {
        pat::strbitset r1 = tight_->getBitTemplate(), r2 = tightLepVeto_->getBitTemplate();
        j.jetIdBits = ((*tight_)(jet, r1) ? 1u : 0u) | ((*tightLepVeto_)(jet, r2) ? 2u : 0u);
      }
      // MiniAOD drops the daughters of some jets (slimmedJetsAK8: pt < 170); an empty EEC would bias silently.
      if (jet.numberOfDaughters() == 0)
        throw cms::Exception("JetSetup") << "PATJetToEECJetProducer: jet " << i << " (pt " << jet.pt()
                                         << ") has no daughters in the input; raise ptMin above the MiniAOD "
                                            "dropDaughters threshold (slimmedJetsAK8: pt < 170) or recluster";
      for (const auto& c : jet.daughterPtrVector()) {
        const float w = applyPuppi_ ? jetsetup::puppiWeight(c, vm) : 1.0f;
        jetsetup::appendConstituent(j, *c, w, c.key());
      }
      out.push_back(std::move(j));
    }
    event.emplace(put_, std::move(out));
  }

private:
  const edm::EDGetTokenT<edm::View<pat::Jet>> src_;
  const edm::EDPutTokenT<srothman::EECJetCollection> put_;
  std::optional<edm::EDGetTokenT<edm::ValueMap<float>>> vmToken_;
  const bool applyPuppi_;
  const double ptMin_;
  std::optional<PFJetIDSelectionFunctor> tight_, tightLepVeto_;
};

DEFINE_FWK_MODULE(PATJetToEECJetProducer);
