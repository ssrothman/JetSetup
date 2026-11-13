#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SRothman/JetSetup/interface/JetTools.h"

#include <string>
#include <vector>

// Per channel (channels[k], candidates vetoSrcs[k]): a std::vector<bool> (instance = the channel name), aligned with
// the jets in src, true for a jet within vetoDR of any candidate; and an extension of the jet table `name` with one
// bool column veto<channel> per channel.
class EECJetOverlapProducer : public edm::global::EDProducer<> {
public:
  explicit EECJetOverlapProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        channels_(cfg.getParameter<std::vector<std::string>>("channels")),
        name_(cfg.getParameter<std::string>("name")),
        vetoDR_(cfg.getParameter<double>("vetoDR")),
        tablePut_(produces()) {
    const auto srcs = cfg.getParameter<std::vector<edm::InputTag>>("vetoSrcs");
    if (srcs.size() != channels_.size())
      throw cms::Exception("Configuration")
          << "EECJetOverlapProducer: " << channels_.size() << " channels but " << srcs.size() << " vetoSrcs";
    for (std::size_t k = 0; k < srcs.size(); ++k) {
      vetoes_.push_back(consumes<edm::View<reco::Candidate>>(srcs[k]));
      flagPuts_.push_back(produces<std::vector<bool>>(channels_[k]));
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("jetsReco"));
    desc.add<std::vector<std::string>>("channels", {});
    desc.add<std::vector<edm::InputTag>>("vetoSrcs", {})->setComment("selected objects of each channel");
    desc.add<double>("vetoDR", 0.8);
    desc.add<std::string>("name", "RecoJet")->setComment("name of the jet table to extend");
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    const auto& jets = event.get(src_);
    auto table = std::make_unique<nanoaod::FlatTable>(jets.size(), name_, false, true);
    for (std::size_t k = 0; k < channels_.size(); ++k) {
      std::vector<std::pair<double, double>> veto;
      for (const auto& c : event.get(vetoes_[k]))
        veto.emplace_back(c.eta(), c.phi());
      std::vector<bool> flags;
      for (const auto& j : jets)
        flags.push_back(jetsetup::overlaps(j, veto, vetoDR_));
      table->addColumn<bool>("veto" + channels_[k],
                             flags,
                             "within " + std::to_string(vetoDR_) + " of a selected " + channels_[k] + " object");
      event.emplace(flagPuts_[k], std::move(flags));
    }
    event.put(tablePut_, std::move(table));
  }

private:
  const edm::EDGetTokenT<srothman::EECJetCollection> src_;
  const std::vector<std::string> channels_;
  const std::string name_;
  const double vetoDR_;
  const edm::EDPutTokenT<nanoaod::FlatTable> tablePut_;
  std::vector<edm::EDGetTokenT<edm::View<reco::Candidate>>> vetoes_;
  std::vector<edm::EDPutTokenT<std::vector<bool>>> flagPuts_;
};

DEFINE_FWK_MODULE(EECJetOverlapProducer);
