#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "SRothman/JetSetup/plugins/ConstituentFill.h"

// reco::GenJet -> EECJet with the gen-jet daughters as constituents; rawPt = pt, no jet ID bits.
class GenJetToEECJetProducer : public edm::global::EDProducer<> {
public:
  explicit GenJetToEECJetProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        put_(produces()),
        ptMin_(cfg.getParameter<double>("ptMin")) {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("slimmedGenJetsAK8"));
    desc.add<double>("ptMin", 0.0);
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    const auto& jets = event.get(src_);
    srothman::EECJetCollection out;
    for (std::size_t i = 0; i < jets.size(); ++i) {
      const auto& jet = jets[i];
      if (jet.pt() < ptMin_)
        continue;
      srothman::EECJet j;
      j.pt = j.rawPt = jet.pt();
      j.eta = jet.eta();
      j.phi = jet.phi();
      j.mass = jet.mass();
      j.srcIdx = i;
      for (const auto& c : jet.daughterPtrVector())
        jetsetup::appendConstituent(j, *c, 1.0f, c.key());
      out.push_back(std::move(j));
    }
    event.emplace(put_, std::move(out));
  }

private:
  const edm::EDGetTokenT<edm::View<reco::GenJet>> src_;
  const edm::EDPutTokenT<srothman::EECJetCollection> put_;
  const double ptMin_;
};

DEFINE_FWK_MODULE(GenJetToEECJetProducer);
