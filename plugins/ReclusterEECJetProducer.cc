#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "SRothman/JetSetup/plugins/ConstituentFill.h"

#include <fastjet/ClusterSequence.hh>

#include <algorithm>
#include <cstdlib>
#include <optional>
#include <vector>

// Clusters candidates with fastjet; jets are uncorrected (rawPt = pt), with no jet ID bits.
// Candidates whose |pdgId| is in dropAbsPdgIds (e.g. neutrinos at gen level) are not clustered.
class ReclusterEECJetProducer : public edm::global::EDProducer<> {
public:
  explicit ReclusterEECJetProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        put_(produces()),
        applyPuppi_(cfg.getParameter<bool>("applyPuppiWeights")),
        jetPtMin_(cfg.getParameter<double>("jetPtMin")),
        drop_(cfg.getParameter<std::vector<int>>("dropAbsPdgIds")),
        def_(algorithm(cfg.getParameter<std::string>("algorithm")), cfg.getParameter<double>("R")) {
    const auto vm = cfg.getParameter<edm::InputTag>("puppiValueMap");
    if (!vm.label().empty())
      vmToken_ = consumes(vm);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("packedPFCandidates"));
    desc.add<std::string>("algorithm", "antikt")->setComment("antikt, kt or cambridge");
    desc.add<double>("R", 0.8);
    desc.add<double>("jetPtMin", 100.0);
    desc.add<bool>("applyPuppiWeights", true);
    desc.add<edm::InputTag>("puppiValueMap", edm::InputTag(""));
    desc.add<std::vector<int>>("dropAbsPdgIds", {});
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    const auto cands = event.getHandle(src_);
    const edm::ValueMap<float>* vm = vmToken_ ? &event.get(*vmToken_) : nullptr;
    std::vector<fastjet::PseudoJet> in;
    std::vector<float> weights;
    for (std::size_t i = 0; i < cands->size(); ++i) {
      const auto& c = (*cands)[i];
      if (std::find(drop_.begin(), drop_.end(), std::abs(c.pdgId())) != drop_.end())
        continue;
      const float w = applyPuppi_ ? jetsetup::puppiWeight(cands->ptrAt(i), vm) : 1.0f;
      if (!(w > 0))
        continue;
      in.emplace_back(w * c.px(), w * c.py(), w * c.pz(), w * c.energy());
      in.back().set_user_index(i);
      weights.push_back(w);
    }
    fastjet::ClusterSequence cs(in, def_);
    srothman::EECJetCollection out;
    const auto jets = fastjet::sorted_by_pt(cs.inclusive_jets(jetPtMin_));
    for (std::size_t k = 0; k < jets.size(); ++k) {
      const auto& fj = jets[k];
      srothman::EECJet j;
      j.pt = j.rawPt = fj.pt();
      j.eta = fj.pseudorapidity();
      j.phi = fj.phi_std();
      j.mass = std::max(fj.m(), 0.0);
      j.srcIdx = k;
      for (const auto& p : fj.constituents()) {
        const auto idx = p.user_index();
        jetsetup::appendConstituent(j, (*cands)[idx], weightOf(idx, in, weights), idx);
      }
      out.push_back(std::move(j));
    }
    event.emplace(put_, std::move(out));
  }

private:
  static fastjet::JetAlgorithm algorithm(const std::string& a) {
    if (a == "antikt")
      return fastjet::antikt_algorithm;
    if (a == "kt")
      return fastjet::kt_algorithm;
    if (a == "cambridge")
      return fastjet::cambridge_algorithm;
    throw cms::Exception("Configuration") << "ReclusterEECJetProducer: unknown algorithm " << a;
  }
  // Weight of candidate idx, found through the input list (user_index is the candidate index).
  static float weightOf(int idx, const std::vector<fastjet::PseudoJet>& in, const std::vector<float>& w) {
    const auto it = std::lower_bound(
        in.begin(), in.end(), idx, [](const fastjet::PseudoJet& p, int i) { return p.user_index() < i; });
    return w[it - in.begin()];
  }

  const edm::EDGetTokenT<edm::View<reco::Candidate>> src_;
  const edm::EDPutTokenT<srothman::EECJetCollection> put_;
  std::optional<edm::EDGetTokenT<edm::ValueMap<float>>> vmToken_;
  const bool applyPuppi_;
  const double jetPtMin_;
  const std::vector<int> drop_;
  const fastjet::JetDefinition def_;
};

DEFINE_FWK_MODULE(ReclusterEECJetProducer);
