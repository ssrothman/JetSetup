#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/transform.h"
#include "SRothman/JetSetup/interface/JetTools.h"

// Keeps the jets passing JetCuts and not within vetoDR of any candidate in vetoSrcs.
class EECJetSelector : public edm::global::EDProducer<> {
public:
  explicit EECJetSelector(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        put_(produces()),
        vetoes_(edm::vector_transform(cfg.getParameter<std::vector<edm::InputTag>>("vetoSrcs"),
                                      [this](const edm::InputTag& t) { return consumes<edm::View<reco::Candidate>>(t); })),
        vetoDR_(cfg.getParameter<double>("vetoDR")) {
    cuts_.ptMin = cfg.getParameter<double>("ptMin");
    cuts_.absEtaMax = cfg.getParameter<double>("absEtaMax");
    cuts_.idMask = cfg.getParameter<unsigned>("idMask");
    cuts_.minConstituents = cfg.getParameter<unsigned>("minConstituents");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("eecJets"));
    desc.add<double>("ptMin", 0.0);
    desc.add<double>("absEtaMax", 1e9);
    desc.add<unsigned>("idMask", 0)->setComment("required jetIdBits: 1 = tight, 2 = tightLepVeto");
    desc.add<unsigned>("minConstituents", 1);
    desc.add<std::vector<edm::InputTag>>("vetoSrcs", {});
    desc.add<double>("vetoDR", 0.8);
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    std::vector<std::pair<double, double>> veto;
    for (const auto& t : vetoes_)
      for (const auto& c : event.get(t))
        veto.emplace_back(c.eta(), c.phi());
    srothman::EECJetCollection out;
    for (const auto& j : event.get(src_))
      if (jetsetup::passes(j, cuts_) && !jetsetup::overlaps(j, veto, vetoDR_))
        out.push_back(j);
    event.emplace(put_, std::move(out));
  }

private:
  const edm::EDGetTokenT<srothman::EECJetCollection> src_;
  const edm::EDPutTokenT<srothman::EECJetCollection> put_;
  const std::vector<edm::EDGetTokenT<edm::View<reco::Candidate>>> vetoes_;
  const double vetoDR_;
  jetsetup::JetCuts cuts_;
};

DEFINE_FWK_MODULE(EECJetSelector);
