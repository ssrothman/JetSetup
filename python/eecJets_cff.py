"""Builds EECJetCollection chains (translate -> [JES/JER] -> constituent selection -> jet selection) with
systematic variations, and their NanoAOD tables.

Nominal reco chain (mode 'pat', MC):  <label>Raw -> <label>JER (nominal) -> <label>Const -> <label>
Variation V gets its own chain ending in <label><V>; see setupJets.
"""
import FWCore.ParameterSet.Config as cms

JME_DIR = '/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/'


def jmeConfig(R=0.8, era='2022_Summer22EE', tag='Summer22EE_22Sep2023', jecVersion='V2', jerVersion='JRV1'):
    """correctionlib file and keys for PUPPI jets of radius 0.4 or 0.8."""
    algo, fn = {0.4: ('AK4PFPuppi', 'jet_jerc.json.gz'), 0.8: ('AK8PFPuppi', 'fatJet_jerc.json.gz')}[R]
    return dict(correctionFile=cms.string(JME_DIR + era + '/' + fn),
                jesKey=cms.string('%s_%s_MC_Total_%s' % (tag, jecVersion, algo)),
                jerScaleFactorKey=cms.string('%s_%s_MC_ScaleFactor_%s' % (tag, jerVersion, algo)),
                jerResolutionKey=cms.string('%s_%s_MC_PtResolution_%s' % (tag, jerVersion, algo)))


# Placeholder constituent variations: the sizes are invented for development, not derived from CMS measurements.
# constDrop uses its own seedSalt so its draws are independent of the JER stochastic smearing (salt 0).
CONSTITUENT_VARIATIONS = {
    'TrkDrop': dict(type='constDrop', categories=['chargedHadron'], dropProbability=0.02, seedSalt=1),
    'ChargedScaleUp': dict(type='constScale', categories=['chargedHadron'], scale=1.01),
    'ChargedScaleDown': dict(type='constScale', categories=['chargedHadron'], scale=0.99),
    'NeutralScaleUp': dict(type='constScale', categories=['neutralHadron'], scale=1.05),
    'NeutralScaleDown': dict(type='constScale', categories=['neutralHadron'], scale=0.95),
    'PhotonScaleUp': dict(type='constScale', categories=['photon'], scale=1.02),
    'PhotonScaleDown': dict(type='constScale', categories=['photon'], scale=0.98),
}

_DEFAULT_SRC = {('reco', 'pat', 0.8): 'slimmedJetsAK8', ('reco', 'pat', 0.4): 'slimmedJetsPuppi',
                ('gen', 'pat', 0.8): 'slimmedGenJetsAK8', ('gen', 'pat', 0.4): 'slimmedGenJets',
                ('reco', 'recluster'): 'packedPFCandidates', ('gen', 'recluster'): 'packedGenParticles'}
_GEN_JETS = {0.8: 'slimmedGenJetsAK8', 0.4: 'slimmedGenJets'}
_NEUTRINOS = [12, 14, 16]


def _add(process, task, label, module):
    setattr(process, label, module)
    if task is not None:
        task.add(module)
    return label


def _jetSelector(src, jetCuts, vetoSrcs, vetoDR):
    p = cms.EDProducer('EECJetSelector', src=cms.InputTag(src),
                       vetoSrcs=cms.VInputTag(*vetoSrcs), vetoDR=cms.double(vetoDR))
    for k, v in jetCuts.items():
        setattr(p, k, cms.double(v) if k in ('ptMin', 'absEtaMax') else cms.uint32(v))
    return p


def _constSelector(src, constituentCuts):
    p = cms.EDProducer('EECConstituentSelector', src=cms.InputTag(src))
    for k, v in constituentCuts.items():
        if k in ('ptMin', 'minPuppiWeight'):
            setattr(p, k, cms.double(v))
        elif k == 'chargedOnly':
            setattr(p, k, cms.bool(v))
        else:
            setattr(p, k, cms.vint32(*v))
    return p


def _variation(src, R, jme, **spec):
    p = cms.EDProducer('EECJetVariationProducer', src=cms.InputTag(src), jetR=cms.double(R),
                       genJets=cms.InputTag(_GEN_JETS[R]), **jme)
    for k, v in spec.items():
        if k in ('type', 'direction'):
            setattr(p, k, cms.string(v))
        elif k == 'categories':
            setattr(p, k, cms.vstring(*v))
        elif k == 'seedSalt':
            setattr(p, k, cms.uint32(v))
        else:
            setattr(p, k, cms.double(v))
    return p


def setupJets(process, label, level='reco', mode='pat', R=0.8, src=None, ptMin=100.0,
              jetIdVersion='RUN3PUPPIruns2022FGruns2023CD', isMC=True, jme=None,
              jetCuts=None, constituentCuts=None, vetoSrcs=(), vetoDR=None,
              variations=None, tableName=None, writeConstituents=False, task=None):
    """Adds the nominal chain and one chain per variation; returns {variationName: final label}, '' = nominal.

    jetCuts: EECJetSelector parameters (ptMin, absEtaMax, idMask, minConstituents).
    constituentCuts: EECConstituentSelector parameters.
    variations: {name: spec}. A spec is either
      dict(type='jes', direction='up'|'down')                         (reco 'pat' MC only)
      dict(type='jer', direction='up'|'down')                         (replaces nominal JER; reco 'pat' MC only)
      dict(type='constScale'|'constDrop', categories=[...], scale=..., dropProbability=..., seedSalt=...)
      dict(type='cuts', jetCuts={...}, constituentCuts={...})         (threshold variation: overrides merged into nominal)
    tableName: if set, adds EECJetTableProducer <final label>Table named tableName + variation name.
    """
    if level not in ('reco', 'gen') or mode not in ('pat', 'recluster'):
        raise ValueError('level must be reco|gen and mode pat|recluster')
    jetCuts = dict(jetCuts or {})
    constituentCuts = dict(constituentCuts or {})
    variations = dict(variations or {})
    vetoDR = R if vetoDR is None else vetoDR
    jme = jme if jme is not None else jmeConfig(R)
    canJERC = level == 'reco' and mode == 'pat' and isMC
    for name, spec in variations.items():
        if spec['type'] in ('jes', 'jer') and not canJERC:
            raise ValueError('variation %s: JES/JER need reco MiniAOD jets in MC' % name)
        if spec['type'] not in ('jes', 'jer', 'constScale', 'constDrop', 'cuts'):
            raise ValueError('variation %s: unknown type %s' % (name, spec['type']))

    if src is None:
        src = _DEFAULT_SRC[(level, mode, R)] if mode == 'pat' else _DEFAULT_SRC[(level, mode)]
    if mode == 'pat' and level == 'reco':
        raw = cms.EDProducer('PATJetToEECJetProducer', src=cms.InputTag(src), applyPuppiWeights=cms.bool(True),
                             puppiValueMap=cms.InputTag(''), jetIdVersion=cms.string(jetIdVersion),
                             ptMin=cms.double(ptMin))
    elif mode == 'pat':
        raw = cms.EDProducer('GenJetToEECJetProducer', src=cms.InputTag(src), ptMin=cms.double(ptMin))
    else:
        raw = cms.EDProducer('ReclusterEECJetProducer', src=cms.InputTag(src), algorithm=cms.string('antikt'),
                             R=cms.double(R), jetPtMin=cms.double(ptMin),
                             applyPuppiWeights=cms.bool(level == 'reco'), puppiValueMap=cms.InputTag(''),
                             dropAbsPdgIds=cms.vint32(*(_NEUTRINOS if level == 'gen' else [])))
    rawLabel = _add(process, task, label + 'Raw', raw)

    def tail(name, start, jc, cc):
        c = _add(process, task, label + name + 'Const', _constSelector(start, cc))
        out = _add(process, task, label + name, _jetSelector(c, jc, vetoSrcs, vetoDR))
        if tableName:
            _add(process, task, out + 'Table',
                 cms.EDProducer('EECJetTableProducer', src=cms.InputTag(out), name=cms.string(tableName + name),
                                doc=cms.string(''), writeConstituents=cms.bool(writeConstituents),
                                precision=cms.int32(-1)))
        return out

    nominalStart = rawLabel
    if canJERC:
        nominalStart = _add(process, task, label + 'JER',
                            _variation(rawLabel, R, jme, type='jer', direction='nominal'))
    result = {'': tail('', nominalStart, jetCuts, constituentCuts)}

    for name, spec in variations.items():
        spec = dict(spec)
        t = spec.pop('type')
        if t == 'jes':
            j = _add(process, task, label + name + 'JES', _variation(rawLabel, R, jme, type='jes', **spec))
            start = _add(process, task, label + name + 'JER', _variation(j, R, jme, type='jer', direction='nominal'))
        elif t == 'jer':
            start = _add(process, task, label + name + 'JER', _variation(rawLabel, R, jme, type='jer', **spec))
        elif t in ('constScale', 'constDrop'):
            start = _add(process, task, label + name + 'Var', _variation(nominalStart, R, jme, type=t, **spec))
        else:
            jc = dict(jetCuts, **spec.get('jetCuts', {}))
            cc = dict(constituentCuts, **spec.get('constituentCuts', {}))
            result[name] = tail(name, nominalStart, jc, cc)
            continue
        result[name] = tail(name, start, jetCuts, constituentCuts)
    return result
