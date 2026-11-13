"""Builds EECJetCollection chains (translate -> [JES/JER] -> constituent selection -> jet selection) with
systematic variations, and their NanoAOD tables.

Every value comes from the caller's collection spec (see simulation/configs/default.yaml in boostedV); there are
no physics defaults here.

Nominal reco chain (mode 'pat', MC):  [<label>JEC ->] <label>Raw -> <label>JER (nominal) -> <label>Const -> <label>
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


RETAG_TAGGERS = dict(GlobalParT=('RecoBTag.ONNXRuntime.pfGlobalParticleTransformerAK8_cff',
                                  '_pfGlobalParticleTransformerAK8JetTagsAll'))


def _retag(process, label, spec, task):
    """Runs the taggers of spec['retag'] on the MiniAOD jets as NanoAOD does (updateJetCollection; the tagger
    re-applies the GlobalTag JECs of retag['jecPayload']) and returns the label of the re-tagged copy, which only
    provides the retag['variables'] columns; None without taggers."""
    r = spec.get('retag', {'taggers': []})
    if not r['taggers']:
        return None
    import importlib
    from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    discs = []
    for t in r['taggers']:
        mod, name = RETAG_TAGGERS[t]
        discs += list(getattr(importlib.import_module(mod), name))
    # The tag-info producers build transient tracks, which need the geometry and the magnetic field.
    for cff in ('Configuration.StandardSequences.GeometryRecoDB_cff', 'Configuration.StandardSequences.MagneticField_cff',
                'TrackingTools.TransientTrack.TransientTrackBuilder_cfi'):
        process.load(cff)
    postfix = label[0].upper() + label[1:] + 'Retag'
    updateJetCollection(process, jetSource=cms.InputTag(spec['src']), pvSource=cms.InputTag(r['vertices']),
                        svSource=cms.InputTag(r['secondaryVertices']), rParam=spec['R'],
                        jetCorrections=(r['jecPayload'], cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
                        btagDiscriminators=discs, postfix=postfix, printWarning=False)
    if task is not None:
        task.add(getPatAlgosToolsTask(process))
    return 'selectedUpdatedPatJets' + postfix


def _jec(spec, isMC):
    # Source order (sortByPt False) keeps srcIdx valid for the re-tagged copy of spec['src'].
    j = spec['jec']
    return cms.EDProducer('PATJetJECProducer', src=cms.InputTag(spec['src']), rho=cms.InputTag(j['rho']),
                          correctionFile=cms.string(j['correctionFile']),
                          key=cms.string(j['mcKey'] if isMC else j['dataKey']), uncKey=cms.string(''),
                          jetIdVersion=cms.string(''), sortByPt=cms.bool(False))


def _translator(spec, src):
    t, level, mode = spec['translator'], spec['level'], spec['mode']
    src = cms.InputTag(src)
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


def setupJets(process, label, spec, isMC, vetoSrcs=None, writeConstituents=False, task=None):
    """Adds the nominal chain and, in MC, one chain per spec['variations'] entry; returns
    {variationName: final label}, '' = nominal.

    spec keys: level (reco|gen), mode (pat|recluster), src, R, translator, jetCuts, constituentCuts, cleaningDR,
    table {name, precision, extraVariables (pat-mode reco jets only; nominal table only)}, retag {taggers, vertices,
    secondaryVertices, jecPayload, variables} (pat-mode reco jets only; the taggers run on a copy of src that only feeds
    the retag variables), jec {correctionFile, mcKey, dataKey, rho} (pat-mode reco jets only: the JEC re-applied from the
    raw momentum before the translator; the table extraVariables then see the re-corrected jets), and for reco 'pat'
    jets jerc {correctionFile, jesKey, jerScaleFactorKey,
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

    if spec['table'].get('extraVariables') and not (spec['level'] == 'reco' and spec['mode'] == 'pat'):
        raise ValueError('%s: table extraVariables need reco MiniAOD (pat) jets' % label)
    if spec.get('retag', {'taggers': []})['taggers'] and not (spec['level'] == 'reco' and spec['mode'] == 'pat'):
        raise ValueError('%s: retag needs reco MiniAOD (pat) jets' % label)
    if 'jec' in spec and not (spec['level'] == 'reco' and spec['mode'] == 'pat'):
        raise ValueError('%s: jec needs reco MiniAOD (pat) jets' % label)
    tagged = _retag(process, label, spec, task)
    jetSrc = _add(process, task, label + 'JEC', _jec(spec, isMC)) if 'jec' in spec else spec['src']
    rawLabel = _add(process, task, label + 'Raw', _translator(spec, jetSrc))

    def tail(name, start, jc, cc):
        c = _add(process, task, label + name + 'Const', _constSelector(start, cc))
        out = _add(process, task, label + name, _jetSelector(c, jc, (), spec['cleaningDR']))
        vetoes = vetoSrcs or {}
        if vetoes:
            # Jets are kept; each channel's overlap flags feed veto<ch> columns and the channel's aggregators.
            _add(process, task, out + 'VetoTable',
                 cms.EDProducer('EECJetOverlapProducer', src=cms.InputTag(out), channels=cms.vstring(*vetoes),
                                vetoSrcs=cms.VInputTag(*vetoes.values()), vetoDR=cms.double(spec['cleaningDR']),
                                name=cms.string(spec['table']['name'] + name)))
        _add(process, task, out + 'Table',
             cms.EDProducer('EECJetTableProducer', src=cms.InputTag(out), name=cms.string(spec['table']['name'] + name),
                            doc=cms.string(''), writeConstituents=cms.bool(writeConstituents),
                            precision=cms.int32(spec['table']['precision'])))
        pset = lambda vs: cms.VPSet(*[cms.PSet(name=cms.string(v['name']), expr=cms.string(v['expr']),
                                               type=cms.string(v['type']), doc=cms.string(v['doc']),
                                               precision=cms.int32(v['precision'])) for v in vs if isMC or not v['mcOnly']])
        extra = spec['table'].get('extraVariables', [])
        taggedVars = spec['retag']['variables'] if tagged else []
        if name == '' and (extra or taggedVars):
            _add(process, task, out + 'ExtTable',
                 cms.EDProducer('EECJetExtensionTableProducer', src=cms.InputTag(out), jets=cms.InputTag(jetSrc),
                                name=cms.string(spec['table']['name']), variables=pset(extra),
                                taggedJets=cms.InputTag(tagged or ''), taggedVariables=pset(taggedVars)))
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
