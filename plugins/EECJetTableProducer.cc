#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "SRothman/JetSetup/interface/JetTools.h"

#include <cmath>

// NanoAOD tables: "<name>" (one row per jet) and, if writeConstituents, "<name>Constituents" (iJet, pt, eta,
// phi, E, charge, pdgId, puppiWeight, srcIdx; momenta after PUPPI weighting and variations).
class EECJetTableProducer : public edm::global::EDProducer<> {
public:
  explicit EECJetTableProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        name_(cfg.getParameter<std::string>("name")),
        doc_(cfg.getParameter<std::string>("doc")),
        writeConstituents_(cfg.getParameter<bool>("writeConstituents")),
        precision_(cfg.getParameter<int>("precision")),
        jetPut_(produces()) {
    if (writeConstituents_)
      constPut_ = produces("constituents");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("eecJets"));
    desc.add<std::string>("name", "EECJet");
    desc.add<std::string>("doc", "");
    desc.add<bool>("writeConstituents", false);
    desc.add<int>("precision", -1)->setComment("mantissa bits of the float columns; -1 = full float");
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    const auto& jets = event.get(src_);
    std::vector<float> pt, eta, phi, mass, rawPt;
    std::vector<std::uint32_t> srcIdx, nConst;
    std::vector<std::uint8_t> id;
    std::vector<std::uint16_t> iJet;
    std::vector<float> cpt, ceta, cphi, cE, cw;
    std::vector<std::int16_t> cq;
    std::vector<std::int32_t> cpdg;
    std::vector<std::uint32_t> csrc;
    for (std::size_t k = 0; k < jets.size(); ++k) {
      const auto& j = jets[k];
      pt.push_back(j.pt);
      eta.push_back(j.eta);
      phi.push_back(j.phi);
      mass.push_back(j.mass);
      rawPt.push_back(j.rawPt);
      srcIdx.push_back(j.srcIdx);
      id.push_back(j.jetIdBits);
      nConst.push_back(j.constituents.size());
      if (!writeConstituents_)
        continue;
      for (const auto& c : j.constituents) {
        const double t = jetsetup::pt(c);
        iJet.push_back(k);
        cpt.push_back(t);
        ceta.push_back(t > 0 ? std::asinh(c.pz / t) : (c.pz >= 0 ? 1e10 : -1e10));
        cphi.push_back(std::atan2(c.py, c.px));
        cE.push_back(c.E);
        cw.push_back(c.puppiWeight);
        cq.push_back(c.charge);
        cpdg.push_back(c.pdgId);
        csrc.push_back(c.srcIdx);
      }
    }
    auto jt = std::make_unique<nanoaod::FlatTable>(jets.size(), name_, false);
    jt->setDoc(doc_);
    jt->addColumn<float>("pt", pt, "pt after corrections and variations", precision_);
    jt->addColumn<float>("eta", eta, "eta", precision_);
    jt->addColumn<float>("phi", phi, "phi", precision_);
    jt->addColumn<float>("mass", mass, "mass after corrections and variations", precision_);
    jt->addColumn<float>("rawPt", rawPt, "uncorrected pt", precision_);
    jt->addColumn<std::uint32_t>("srcIdx", srcIdx, "index in the source jet collection");
    jt->addColumn<std::uint8_t>("jetId", id, "bit 0: tight, bit 1: tightLepVeto");
    jt->addColumn<std::uint32_t>("nConstituents", nConst, "number of constituents (after selection)");
    event.put(jetPut_, std::move(jt));
    if (writeConstituents_) {
      auto ct = std::make_unique<nanoaod::FlatTable>(iJet.size(), name_ + "Constituents", false);
      ct->setDoc(doc_);
      ct->addColumn<std::uint16_t>("iJet", iJet, "index of the jet in " + name_);
      ct->addColumn<float>("pt", cpt, "pt", precision_);
      ct->addColumn<float>("eta", ceta, "eta", precision_);
      ct->addColumn<float>("phi", cphi, "phi", precision_);
      ct->addColumn<float>("E", cE, "energy", precision_);
      ct->addColumn<std::int16_t>("charge", cq, "charge");
      ct->addColumn<std::int32_t>("pdgId", cpdg, "pdgId");
      ct->addColumn<float>("puppiWeight", cw, "PUPPI weight applied to the momentum");
      ct->addColumn<std::uint32_t>("srcIdx", csrc, "index in the source candidate collection");
      event.put(constPut_, std::move(ct));
    }
  }

private:
  const edm::EDGetTokenT<srothman::EECJetCollection> src_;
  const std::string name_, doc_;
  const bool writeConstituents_;
  const int precision_;
  const edm::EDPutTokenT<nanoaod::FlatTable> jetPut_;
  edm::EDPutTokenT<nanoaod::FlatTable> constPut_;
};

DEFINE_FWK_MODULE(EECJetTableProducer);
