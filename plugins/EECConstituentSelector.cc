#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "SRothman/JetSetup/interface/JetTools.h"

// Removes constituents failing ConstituentCuts; jets keep their four-momentum and stay in the collection.
class EECConstituentSelector : public edm::global::EDProducer<> {
public:
  explicit EECConstituentSelector(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))), put_(produces()) {
    cuts_.ptMin = cfg.getParameter<double>("ptMin");
    cuts_.chargedOnly = cfg.getParameter<bool>("chargedOnly");
    cuts_.keepAbsPdgIds = cfg.getParameter<std::vector<int>>("keepAbsPdgIds");
    cuts_.dropAbsPdgIds = cfg.getParameter<std::vector<int>>("dropAbsPdgIds");
    cuts_.minPuppiWeight = cfg.getParameter<double>("minPuppiWeight");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("eecJets"));
    desc.add<double>("ptMin", 0.0)->setComment("on the (PUPPI-weighted) constituent pt");
    desc.add<bool>("chargedOnly", false);
    desc.add<std::vector<int>>("keepAbsPdgIds", {})->setComment("empty: keep all");
    desc.add<std::vector<int>>("dropAbsPdgIds", {});
    desc.add<double>("minPuppiWeight", 0.0);
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    auto out = event.get(src_);
    for (auto& j : out)
      jetsetup::selectConstituents(j, cuts_);
    event.emplace(put_, std::move(out));
  }

private:
  const edm::EDGetTokenT<srothman::EECJetCollection> src_;
  const edm::EDPutTokenT<srothman::EECJetCollection> put_;
  jetsetup::ConstituentCuts cuts_;
};

DEFINE_FWK_MODULE(EECConstituentSelector);
