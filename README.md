# SRothman/JetSetup

Builds `srothman::EECJetCollection`s (see `SRothman/DataFormats`) from MiniAOD for the EEC plugins in `SRothman/EECs`:
translation, jet and constituent selection, and systematic variations, plus NanoAOD tables.

## Modules
| Module | Purpose |
|---|---|
| `PATJetToEECJetProducer` | `pat::Jet` → EECJet. Constituents are the daughters, scaled by their PUPPI weight (packed `puppiWeight()` or a ValueMap). `jetIdBits`: bit 0 TIGHT, bit 1 TIGHTLEPVETO (`PFJetIDSelectionFunctor`, default `RUN3PUPPIruns2022FGruns2023CD`). **Throws for a jet above `ptMin` without daughters.** MiniAOD drops the daughters of `slimmedJetsAK8` with pt < 170 (`PhysicsTools/PatAlgos/python/slimming/applySubstructure_cff.py`), so `ptMin` must be ≥ 170 for that collection. |
| `GenJetToEECJetProducer` | `reco::GenJet` → EECJet with its daughters (packed gen particles). |
| `ReclusterEECJetProducer` | fastjet reclustering of any candidate collection (anti-kt/kt/CA, any R). Optional PUPPI weights; `dropAbsPdgIds` (for example neutrinos at gen level). Jets are uncorrected, with no ID. |
| `EECJetSelector` | pt, \|η\|, required `jetIdBits`, minimum constituents, optional ΔR removal against `View<reco::Candidate>` inputs (unused by `setupJets`, which flags instead). |
| `EECJetOverlapProducer` | per channel: `std::vector<bool>` (instance = channel) of jets within `vetoDR` of the channel's objects, and `veto<channel>` columns extending the jet table. |
| `EECConstituentSelector` | constituent pt (after PUPPI weighting), charged only, \|pdgId\| keep/drop lists, minimum PUPPI weight. |
| `EECJetVariationProducer` | one variation per instance: `jes` (up/down, Total uncertainty), `jer` (hybrid smearing: nominal/up/down), `constScale`, `constDrop` (by constituent category). JES/JER change only the jet four-momentum. Corrections are read from the JME correctionlib JSONs. |
| `EECJetTableProducer` | NanoAOD tables `<name>` (jets) and optionally `<name>Constituents`. |
| `EECJetExtensionTableProducer` | extension of a jet table: YAML string-expression columns on the source `pat::Jet` (`src[srcIdx]`), optionally on a re-tagged copy in the same order (`taggedJets`, `taggedVariables`; a jet that differs in η/φ is an error), and `jerFactor` (the JER smearing factor applied, accumulated in `EECJet::jerFactor` by `EECJetVariationProducer`). Nominal table only. |
| `PATJetJECProducer` | pat::Jet copies with the JEC re-applied from the raw momentum (correctionlib compound key), sorted by pt if `sortByPt` (else in source order); `correctedP4("Uncorrected")` and the energy fractions still refer to the raw jet; userFloats `rawPt`, `jecFactor`, `jesUncTotal`; userInts `srcIdx`, `jetId`. A stream module (the jet-ID functor is not const). Used for the AK4 table and, with a `jec` block, before the AK8 translator. |

Random numbers (JER stochastic smearing, constituent drop) are seeded by (run, lumi, event, source-jet index, `seedSalt`)
and not by the module label. Every variation chain with the same salt therefore draws the same numbers for a given jet.

Constituent categories: `chargedHadron` (charged, not e/μ), `neutralHadron`, `photon`, `electron`, `muon`, `all`.

## Python
`SRothman.JetSetup.eecJets_cff.setupJets(process, label, spec, isMC, vetoSrcs, writeConstituents, task)` (vetoSrcs: {channel: InputTag of its selected objects}; jets are flagged per channel, not removed) builds one jet
collection from a spec dict. With a `retag` block, the listed taggers (`GlobalParT`) run on a copy of `src` through `updateJetCollection` (GlobalTag JECs of `retag.jecPayload`, as NanoAOD; geometry, magnetic field and TransientTrackBuilder are loaded); the copy only feeds the `retag.variables` columns.
The spec is one entry of `jets/collections` in the boostedV analysis YAML
(`simulation/configs/default.yaml`, schema `analysis.schema.yaml`), and the helper has no physics defaults of its own.
It returns `{variationName: final label}`, where `''` is the nominal chain:

- nominal (reco, `pat`, MC): `<label>Raw → <label>JER (nominal) → <label>Const → <label>`
- `jes`: `Raw → JES up/down → JER nominal → Const → <label><Name>`
- `jer`: `Raw → JER up/down → Const → …`
- `constScale` / `constDrop`: `<label>JER → variation → Const → …`
- `cuts`: the nominal chain with overridden `jetCuts` / `constituentCuts`.

Variations, and JER smearing, are built in MC only. The JME correction file and keys and the variation sizes come
from the spec; the constituent-variation sizes in `default.yaml` are placeholders.

## Tests
- `scram b runtests`: `testSRothmanJetSetup` (framework-free tools) and `testSRothmanJetSetupModules` (TestProcessor; the JES/JER checks use the cvmfs JME JSONs).
- Real-data validation (manual, needs a MiniAOD file):
  `cmsRun test/validation/validateJets_cfg.py config=<analysis yaml> inputFiles=file:... outputFile=jets.root maxEvents=2000` then
  `python3 -I test/validation/checkJets.py jets_numEvent2000.root`.
