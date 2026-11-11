# SRothman/JetSetup

Builds `srothman::EECJetCollection`s (see `SRothman/DataFormats`) from MiniAOD for the EEC plugins in `SRothman/EECs`:
translation, jet and constituent selection, and systematic variations, plus NanoAOD tables.

## Modules
| Module | Purpose |
|---|---|
| `PATJetToEECJetProducer` | `pat::Jet` → EECJet. Constituents are the daughters, scaled by their PUPPI weight (packed `puppiWeight()` or a ValueMap). `jetIdBits`: bit 0 TIGHT, bit 1 TIGHTLEPVETO (`PFJetIDSelectionFunctor`, default `RUN3PUPPIruns2022FGruns2023CD`). **Throws for a jet above `ptMin` without daughters.** MiniAOD drops the daughters of `slimmedJetsAK8` with pt < 170 (`PhysicsTools/PatAlgos/python/slimming/applySubstructure_cff.py`), so `ptMin` must be ≥ 170 for that collection. |
| `GenJetToEECJetProducer` | `reco::GenJet` → EECJet with its daughters (packed gen particles). |
| `ReclusterEECJetProducer` | fastjet reclustering of any candidate collection (anti-kt/kt/CA, any R). Optional PUPPI weights; `dropAbsPdgIds` (for example neutrinos at gen level). Jets are uncorrected, with no ID. |
| `EECJetSelector` | pt, \|η\|, required `jetIdBits`, minimum constituents, ΔR cleaning against any `View<reco::Candidate>` inputs. |
| `EECConstituentSelector` | constituent pt (after PUPPI weighting), charged only, \|pdgId\| keep/drop lists, minimum PUPPI weight. |
| `EECJetVariationProducer` | one variation per instance: `jes` (up/down, Total uncertainty), `jer` (hybrid smearing: nominal/up/down), `constScale`, `constDrop` (by constituent category). JES/JER change only the jet four-momentum. Corrections are read from the JME correctionlib JSONs. |
| `EECJetTableProducer` | NanoAOD tables `<name>` (jets) and optionally `<name>Constituents`. |

Random numbers (JER stochastic smearing, constituent drop) are seeded by (run, lumi, event, source-jet index, `seedSalt`)
and not by the module label. Every variation chain with the same salt therefore draws the same numbers for a given jet.

Constituent categories: `chargedHadron` (charged, not e/μ), `neutralHadron`, `photon`, `electron`, `muon`, `all`.

## Python
`SRothman.JetSetup.eecJets_cff.setupJets(process, label, level='reco'|'gen', mode='pat'|'recluster', R, ...)` adds the
chains below and returns `{variationName: final label}`, where `''` is the nominal chain:

- nominal (reco, `pat`, MC): `<label>Raw → <label>JER (nominal) → <label>Const → <label>`
- `jes`: `Raw → JES up/down → JER nominal → Const → <label><Name>`
- `jer`: `Raw → JER up/down → Const → …`
- `constScale` / `constDrop`: `<label>JER → variation → Const → …`
- `cuts`: the nominal chain with overridden `jetCuts` / `constituentCuts`.

`jmeConfig(R)` gives the correctionlib file and keys for Summer22EE AK4/AK8 PUPPI jets.

`CONSTITUENT_VARIATIONS` holds **placeholder** constituent variations with invented sizes, which are not CMS measurements: `TrkDrop` drops 2% of charged hadrons (`seedSalt` 1, independent of JER smearing); `ChargedScale`, `NeutralScale` and `PhotonScale` Up/Down shift by ±1%, ±5% and ±2%.

## Tests
- `scram b runtests`: `testSRothmanJetSetup` (framework-free tools) and `testSRothmanJetSetupModules` (TestProcessor; the JES/JER checks use the cvmfs JME JSONs).
- Real-data validation (manual, needs a MiniAOD file):
  `cmsRun test/validation/validateJets_cfg.py inputFiles=file:... outputFile=jets.root maxEvents=2000` then
  `python3 -I test/validation/checkJets.py jets_numEvent2000.root`.
