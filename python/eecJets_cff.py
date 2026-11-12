"""Builds EECJetCollection chains (translate -> [JES/JER] -> constituent selection -> jet selection) with
systematic variations, and their NanoAOD tables.

Every value comes from the caller's collection spec (see simulation/configs/default.yaml in boostedV); there are
no physics defaults here.

Nominal reco chain (mode 'pat', MC):  <label>Raw -> <label>JER (nominal) -> <label>Const -> <label>
Variation V gets its own chain ending in <label><V>:
  jes        Raw -> JES up/down -> JER nominal -> Const -> <label><V>
  jer        Raw -> JER up/down -> Const -> ...
  constScale / constDrop   <label>JER -> variation -> Const -> ...
  cuts       nominal chain with overridden jetCuts / constituentCuts
"""
import FWCore.ParameterSet.Config as cms


def _add(process, task, label, module):
    setattr(process, label, module)
    if task is not None:
        task.add(module)
    return label


def _jetSelector(src, c, vetoSrcs, vetoDR):
    return cms.EDProducer('EECJetSelector', src=cms.InputTag(src), ptMin=cms.double(c['ptMin']),
                          absEtaMax=cms.double(c['absEtaMax']), idMask=cms.uint32(c['idMask']),
                          minConstituents=cms.uint32(c['minConstituents']),
                          vetoSrcs=cms.VInputTag(*vetoSrcs), vetoDR=cms.double(vetoDR))


def _constSelector(src, c):
    return cms.EDProducer('EECConstituentSelector', src=cms.InputTag(src), ptMin=cms.double(c['ptMin']),
                          chargedOnly=cms.bool(c['chargedOnly']), keepAbsPdgIds=cms.vint32(*c['keepAbsPdgIds']),
                          dropAbsPdgIds=cms.vint32(*c['dropAbsPdgIds']), minPuppiWeight=cms.double(c['minPuppiWeight']))


def _variation(src, R, jerc, vtype, direction='nominal', categories=(), scale=1.0, dropProbability=0.0, seedSalt=0):
    # Parameters a variation type does not use are set to inert values so that every parameter is explicit.
    j = jerc or dict(correctionFile='', jesKey='', jerScaleFactorKey='', jerScaleFactorUncertaintyKey='',
                          jerResolutionKey='', rho='', genJets='')
    return cms.EDProducer('EECJetVariationProducer', src=cms.InputTag(src), type=cms.string(vtype),
                          direction=cms.string(direction), correctionFile=cms.string(j['correctionFile']),
                          jesKey=cms.string(j['jesKey']), jerScaleFactorKey=cms.string(j['jerScaleFactorKey']),
                          jerScaleFactorUncertaintyKey=cms.string(j['jerScaleFactorUncertaintyKey']),
                          jerResolutionKey=cms.string(j['jerResolutionKey']), rho=cms.InputTag(j['rho']),
                          genJets=cms.InputTag(j['genJets']), jetR=cms.double(R), categories=cms.vstring(*categories),
                          scale=cms.double(scale), dropProbability=cms.double(dropProbability),
                          seedSalt=cms.uint32(seedSalt))


def _translator(spec):
    t, level, mode = spec['translator'], spec['level'], spec['mode']
    src = cms.InputTag(spec['src'])
    if mode == 'pat' and level == 'reco':
        return cms.EDProducer('PATJetToEECJetProducer', src=src, applyPuppiWeights=cms.bool(t['applyPuppiWeights']),
                              puppiValueMap=cms.InputTag(t['puppiValueMap']),
                              jetIdVersion=cms.string(t['jetIdVersion']), ptMin=cms.double(t['ptMin']))
    if mode == 'pat':
        return cms.EDProducer('GenJetToEECJetProducer', src=src, ptMin=cms.double(t['ptMin']))
    return cms.EDProducer('ReclusterEECJetProducer', src=src, algorithm=cms.string(t['algorithm']),
                          R=cms.double(spec['R']), jetPtMin=cms.double(t['ptMin']),
                          applyPuppiWeights=cms.bool(t['applyPuppiWeights']),
                          puppiValueMap=cms.InputTag(t['puppiValueMap']), dropAbsPdgIds=cms.vint32(*t['dropAbsPdgIds']))


def setupJets(process, label, spec, isMC, vetoSrcs=(), writeConstituents=False, task=None):
    """Adds the nominal chain and, in MC, one chain per spec['variations'] entry; returns
    {variationName: final label}, '' = nominal.

    spec keys: level (reco|gen), mode (pat|recluster), src, R, translator, jetCuts, constituentCuts, cleaningDR,
    table {name, precision}, and for reco 'pat' jets jerc {correctionFile, jesKey, jerScaleFactorKey,
    jerScaleFactorUncertaintyKey, jerResolutionKey, rho, genJets, jerSeedSalt} and variations {name: {type, ...}}.
    JER smearing (nominal) is applied to reco 'pat' jets in MC only.
    """
    if spec['level'] not in ('reco', 'gen') or spec['mode'] not in ('pat', 'recluster'):
        raise ValueError('%s: level must be reco|gen and mode pat|recluster' % label)
    R = spec['R']
    canJERC = spec['level'] == 'reco' and spec['mode'] == 'pat' and isMC
    jerc = spec.get('jerc')
    variations = spec.get('variations', {}) if isMC else {}
    for name, v in variations.items():
        if v['type'] in ('jes', 'jer') and not canJERC:
            raise ValueError('%s variation %s: JES/JER need reco MiniAOD jets in MC' % (label, name))
    if canJERC and jerc is None:
        raise ValueError('%s: reco MiniAOD jets in MC need a jerc block' % label)

    rawLabel = _add(process, task, label + 'Raw', _translator(spec))

    def tail(name, start, jc, cc):
        c = _add(process, task, label + name + 'Const', _constSelector(start, cc))
        out = _add(process, task, label + name, _jetSelector(c, jc, vetoSrcs, spec['cleaningDR']))
        _add(process, task, out + 'Table',
             cms.EDProducer('EECJetTableProducer', src=cms.InputTag(out), name=cms.string(spec['table']['name'] + name),
                            doc=cms.string(''), writeConstituents=cms.bool(writeConstituents),
                            precision=cms.int32(spec['table']['precision'])))
        return out

    def jer(src, direction):
        return _variation(src, R, jerc, 'jer', direction=direction, seedSalt=jerc['jerSeedSalt'])

    nominalStart = _add(process, task, label + 'JER', jer(rawLabel, 'nominal')) if canJERC else rawLabel
    result = {'': tail('', nominalStart, spec['jetCuts'], spec['constituentCuts'])}

    for name, v in variations.items():
        v = dict(v)
        t = v.pop('type')
        if t == 'jes':
            j = _add(process, task, label + name + 'JES', _variation(rawLabel, R, jerc, 'jes', direction=v['direction']))
            start = _add(process, task, label + name + 'JER', jer(j, 'nominal'))
        elif t == 'jer':
            start = _add(process, task, label + name + 'JER', jer(rawLabel, v['direction']))
        elif t == 'constScale':
            start = _add(process, task, label + name + 'Var',
                         _variation(nominalStart, R, jerc, t, categories=v['categories'], scale=v['scale']))
        elif t == 'constDrop':
            start = _add(process, task, label + name + 'Var',
                         _variation(nominalStart, R, jerc, t, categories=v['categories'],
                                    dropProbability=v['dropProbability'], seedSalt=v['seedSalt']))
        elif t == 'cuts':
            jc = dict(spec['jetCuts'], **v.get('jetCuts', {}))
            cc = dict(spec['constituentCuts'], **v.get('constituentCuts', {}))
            result[name] = tail(name, nominalStart, jc, cc)
            continue
        else:
            raise ValueError('%s variation %s: unknown type %s' % (label, name, t))
        result[name] = tail(name, start, spec['jetCuts'], spec['constituentCuts'])
    return result
