#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SRothman/JetSetup/interface/JetTools.h"

#include "correction.h"

#include <optional>
#include <random>

// Applies one variation to every jet:
//  none        copy
//  jes         jet p4 *= 1 +- Total uncertainty (direction up | down)
//  jer         hybrid JER smearing of the jet p4 (direction nominal | up | down), MC only
//  constScale  constituents in `categories` scaled by `scale`
//  constDrop   constituents in `categories` dropped with probability `dropProbability`
// JES/JER change only the jet four-momentum, not the constituents. Random numbers are seeded per
// (event, jet srcIdx, seedSalt), so all chains with the same salt draw the same numbers for a jet.
class EECJetVariationProducer : public edm::global::EDProducer<> {
  enum class Type { None, JES, JER, ConstScale, ConstDrop };

public:
  explicit EECJetVariationProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        put_(produces()),
        type_(parseType(cfg.getParameter<std::string>("type"))),
        direction_(cfg.getParameter<std::string>("direction")),
        jetR_(cfg.getParameter<double>("jetR")),
        scale_(cfg.getParameter<double>("scale")),
        dropP_(cfg.getParameter<double>("dropProbability")),
        salt_(cfg.getParameter<unsigned>("seedSalt")) {
    try {
      mask_ = jetsetup::categoryMask(cfg.getParameter<std::vector<std::string>>("categories"));
    } catch (const std::exception& ex) {
      throw cms::Exception("Configuration") << "EECJetVariationProducer: " << ex.what();
    }
    if (type_ == Type::JES || type_ == Type::JER) {
      const auto cset = correction::CorrectionSet::from_file(cfg.getParameter<std::string>("correctionFile"));
      if (type_ == Type::JES) {
        if (direction_ != "up" && direction_ != "down")
          throw cms::Exception("Configuration") << "EECJetVariationProducer: jes direction must be up or down";
        jes_ = cset->at(cfg.getParameter<std::string>("jesKey"));
      } else {
        if (direction_ != "nominal" && direction_ != "up" && direction_ != "down")
          throw cms::Exception("Configuration") << "EECJetVariationProducer: jer direction must be nominal, up or down";
        jerSF_ = cset->at(cfg.getParameter<std::string>("jerScaleFactorKey"));
        jerRes_ = cset->at(cfg.getParameter<std::string>("jerResolutionKey"));
        rho_ = consumes(cfg.getParameter<edm::InputTag>("rho"));
        genJets_ = consumes(cfg.getParameter<edm::InputTag>("genJets"));
      }
    }
    if ((type_ == Type::ConstScale || type_ == Type::ConstDrop) && mask_ == 0)
      throw cms::Exception("Configuration") << "EECJetVariationProducer: no constituent categories given";
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("eecJets"));
    desc.add<std::string>("type", "none")->setComment("none, jes, jer, constScale or constDrop");
    desc.add<std::string>("direction", "nominal")->setComment("jes: up | down; jer: nominal | up | down");
    desc.add<std::string>("correctionFile",
                          "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer22EE/"
                          "fatJet_jerc.json.gz");
    desc.add<std::string>("jesKey", "Summer22EE_22Sep2023_V2_MC_Total_AK8PFPuppi");
    desc.add<std::string>("jerScaleFactorKey", "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK8PFPuppi");
    desc.add<std::string>("jerResolutionKey", "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK8PFPuppi");
    desc.add<edm::InputTag>("rho", edm::InputTag("fixedGridRhoFastjetAll"));
    desc.add<edm::InputTag>("genJets", edm::InputTag("slimmedGenJetsAK8"));
    desc.add<double>("jetR", 0.8)->setComment("jet radius for the JER gen-jet match");
    desc.add<std::vector<std::string>>("categories", {})
        ->setComment("chargedHadron, neutralHadron, photon, electron, muon, all");
    desc.add<double>("scale", 1.0);
    desc.add<double>("dropProbability", 0.0);
    desc.add<unsigned>("seedSalt", 0);
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    auto out = event.get(src_);
    const auto& id = event.id();
    double rho = 0;
    std::vector<std::pair<double, std::array<double, 2>>> gen;  // (pt, (eta, phi))
    if (type_ == Type::JER) {
      rho = event.get(rho_);
      for (const auto& g : event.get(genJets_))
        gen.push_back({g.pt(), {{g.eta(), g.phi()}}});
    }
    for (auto& j : out) {
      std::mt19937_64 rng(jetsetup::seed(id.run(), id.luminosityBlock(), id.event(), j.srcIdx, salt_));
      switch (type_) {
        case Type::None:
          break;
        case Type::JES: {
          const double u = jes_->evaluate({j.eta, j.pt});
          jetsetup::scaleJet(j, direction_ == "up" ? 1 + u : 1 - u);
          break;
        }
        case Type::JER: {
          const double sf = jerSF_->evaluate({j.eta, j.pt, direction_ == "nominal" ? std::string("nom") : direction_});
          const double res = jerRes_->evaluate({j.eta, j.pt, rho});
          double genPt = -1, best = 1e9;
          for (const auto& [gpt, ep] : gen) {
            const double dR = jetsetup::deltaR(j.eta, j.phi, ep[0], ep[1]);
            if (dR < best && jetsetup::jerMatch(dR, j.pt, gpt, jetR_, res)) {
              best = dR;
              genPt = gpt;
            }
          }
          std::normal_distribution<double> gauss(0.0, 1.0);
          jetsetup::scaleJet(j, jetsetup::jerFactor(j.pt, genPt, sf, res, gauss(rng)));
          break;
        }
        case Type::ConstScale:
          jetsetup::scaleConstituents(j, mask_, scale_);
          break;
        case Type::ConstDrop:
          jetsetup::dropConstituents(j, mask_, dropP_, rng);
          break;
      }
    }
    event.emplace(put_, std::move(out));
  }

private:
  static Type parseType(const std::string& t) {
    if (t == "none")
      return Type::None;
    if (t == "jes")
      return Type::JES;
    if (t == "jer")
      return Type::JER;
    if (t == "constScale")
      return Type::ConstScale;
    if (t == "constDrop")
      return Type::ConstDrop;
    throw cms::Exception("Configuration") << "EECJetVariationProducer: unknown type " << t;
  }

  const edm::EDGetTokenT<srothman::EECJetCollection> src_;
  const edm::EDPutTokenT<srothman::EECJetCollection> put_;
  const Type type_;
  const std::string direction_;
  const double jetR_, scale_, dropP_;
  const unsigned salt_;
  std::uint32_t mask_ = 0;
  correction::Correction::Ref jes_, jerSF_, jerRes_;
  edm::EDGetTokenT<double> rho_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> genJets_;
};

DEFINE_FWK_MODULE(EECJetVariationProducer);
