#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SRothman/DataFormats/interface/EECJet.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Extension of the jet table `name`: for each EECJet in src, `variables` (PAT string expressions; a cut for bool
// columns) evaluated on the source pat::Jet jets[srcIdx], `taggedVariables` on taggedJets[srcIdx] (a re-tagged copy of
// jets in the same order; a jet that differs in eta/phi is an error), and jerFactor (the JER smearing factor applied).
class EECJetExtensionTableProducer : public edm::global::EDProducer<> {
  struct Column {
    std::string name, doc, type;
    int precision;
    std::shared_ptr<StringObjectFunction<pat::Jet>> f;
    std::shared_ptr<StringCutObjectSelector<pat::Jet>> cut;
  };

public:
  explicit EECJetExtensionTableProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        jets_(consumes(cfg.getParameter<edm::InputTag>("jets"))),
        put_(produces()),
        name_(cfg.getParameter<std::string>("name")),
        columns_(parse(cfg.getParameter<std::vector<edm::ParameterSet>>("variables"))),
        tagged_(parse(cfg.getParameter<std::vector<edm::ParameterSet>>("taggedVariables"))) {
    if (!tagged_.empty())
      taggedJets_ = consumes(cfg.getParameter<edm::InputTag>("taggedJets"));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("jetsReco"))->setComment("the table's EECJetCollection");
    desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK8"))->setComment("the pat::Jet source of srcIdx");
    desc.add<std::string>("name", "RecoJet")->setComment("name of the jet table to extend");
    edm::ParameterSetDescription var;
    var.add<std::string>("name");
    var.add<std::string>("expr");
    var.add<std::string>("type")->setComment("float | int | int16 | uint8 | bool (expr is then a cut)");
    var.add<std::string>("doc");
    var.add<int>("precision")->setComment("float mantissa bits; -1 = full float");
    desc.addVPSet("variables", var, {});
    desc.add<edm::InputTag>("taggedJets", edm::InputTag(""))->setComment("re-tagged copy of jets, for taggedVariables");
    desc.addVPSet("taggedVariables", var, {});
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    const auto& eec = event.get(src_);
    const auto& jets = event.get(jets_);
    std::vector<const pat::Jet*> src, tagged;
    std::vector<float> jer;
    for (const auto& j : eec) {
      if (j.srcIdx >= jets.size())
        throw cms::Exception("EECJetExtensionTableProducer")
            << name_ << ": srcIdx " << j.srcIdx << " outside the " << jets.size() << " source jets";
      src.push_back(&jets[j.srcIdx]);
      jer.push_back(j.jerFactor);
    }
    if (!tagged_.empty()) {
      const auto& tj = event.get(taggedJets_);
      if (tj.size() != jets.size())
        throw cms::Exception("EECJetExtensionTableProducer")
            << name_ << ": " << tj.size() << " re-tagged jets for " << jets.size() << " source jets";
      for (const auto& j : eec) {
        const auto& t = tj[j.srcIdx];
        if (reco::deltaR2(jets[j.srcIdx], t) > 1e-6)
          throw cms::Exception("EECJetExtensionTableProducer")
              << name_ << ": re-tagged jets are not in the source order";
        tagged.push_back(&t);
      }
    }
    auto table = std::make_unique<nanoaod::FlatTable>(eec.size(), name_, false, true);
    for (const auto& c : columns_)
      fillColumn(*table, c, src);
    for (const auto& c : tagged_)
      fillColumn(*table, c, tagged);
    table->addColumn<float>("jerFactor", jer, "JER smearing factor applied to the jet four-momentum (1: none)");
    event.put(put_, std::move(table));
  }

private:
  static std::vector<Column> parse(const std::vector<edm::ParameterSet>& vars) {
    std::vector<Column> columns;
    for (const auto& v : vars) {
      const auto expr = v.getParameter<std::string>("expr");
      Column c{v.getParameter<std::string>("name"),
               v.getParameter<std::string>("doc"),
               v.getParameter<std::string>("type"),
               v.getParameter<int>("precision"),
               nullptr,
               nullptr};
      if (c.type == "bool")
        c.cut = std::make_shared<StringCutObjectSelector<pat::Jet>>(expr);
      else if (c.type == "float" || c.type == "int" || c.type == "int16" || c.type == "uint8")
        c.f = std::make_shared<StringObjectFunction<pat::Jet>>(expr);
      else
        throw cms::Exception("Configuration") << "EECJetExtensionTableProducer: column " << c.name << ": type "
                                              << c.type << " is not float, int, int16, uint8 or bool";
      if (c.name == "jerFactor")
        throw cms::Exception("Configuration") << "EECJetExtensionTableProducer: jerFactor is a built-in column";
      columns.push_back(std::move(c));
    }
    return columns;
  }

  static void fillColumn(nanoaod::FlatTable& t, const Column& c, const std::vector<const pat::Jet*>& jets) {
    if (c.type == "float")
      fill<float>(t, c, jets);
    else if (c.type == "int")
      fill<int>(t, c, jets);
    else if (c.type == "int16")
      fill<std::int16_t>(t, c, jets);
    else if (c.type == "uint8")
      fill<std::uint8_t>(t, c, jets);
    else
      fill<bool>(t, c, jets);
  }

  template <typename T>
  static void fill(nanoaod::FlatTable& t, const Column& c, const std::vector<const pat::Jet*>& jets) {
    std::vector<T> v;
    for (const auto* j : jets)
      v.push_back(c.cut ? T((*c.cut)(*j)) : static_cast<T>((*c.f)(*j)));
    if constexpr (std::is_same_v<T, float>)
      t.addColumn<T>(c.name, v, c.doc, c.precision);
    else
      t.addColumn<T>(c.name, v, c.doc);
  }

  const edm::EDGetTokenT<srothman::EECJetCollection> src_;
  const edm::EDGetTokenT<edm::View<pat::Jet>> jets_;
  const edm::EDPutTokenT<nanoaod::FlatTable> put_;
  const std::string name_;
  const std::vector<Column> columns_, tagged_;
  edm::EDGetTokenT<edm::View<pat::Jet>> taggedJets_;
};

DEFINE_FWK_MODULE(EECJetExtensionTableProducer);
