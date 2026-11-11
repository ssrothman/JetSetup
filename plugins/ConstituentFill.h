#ifndef SRothman_JetSetup_plugins_ConstituentFill_h
#define SRothman_JetSetup_plugins_ConstituentFill_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SRothman/DataFormats/interface/EECJet.h"

#include <cmath>

namespace jetsetup {

  // PUPPI weight of a constituent: from the ValueMap if given, else pat::PackedCandidate::puppiWeight().
  inline float puppiWeight(const edm::Ptr<reco::Candidate>& c, const edm::ValueMap<float>* vm) {
    if (vm)
      return vm->get(c.id(), c.key());
    const auto* pc = dynamic_cast<const pat::PackedCandidate*>(c.get());
    if (!pc)
      throw cms::Exception("JetSetup") << "PUPPI weight requested for a constituent that is not a pat::PackedCandidate"
                                          " and no puppiValueMap was given";
    return pc->puppiWeight();
  }

  // Appends the candidate scaled by w; a candidate whose scaled three-momentum vanishes is skipped,
  // because its direction (and so every EEC angle) is undefined.
  inline void appendConstituent(srothman::EECJet& j, const reco::Candidate& c, float w, unsigned srcIdx) {
    const double px = w * c.px(), py = w * c.py(), pz = w * c.pz();
    if (!(px * px + py * py + pz * pz > 0))
      return;
    j.constituents.push_back(
        {px, py, pz, w * c.energy(), static_cast<std::int8_t>(c.charge()), c.pdgId(), w, srcIdx});
  }

}  // namespace jetsetup

#endif
