"""Runs every JetSetup translator mode and a set of variations on MiniAOD and writes NANO tables.
Usage: cmsRun validateJets_cfg.py config=<analysis yaml> inputFiles=file:... outputFile=jets.root maxEvents=2000
The jet collections are derived from the config's Reco and Gen collections; checkJets.py needs the variations below."""
import copy

import FWCore.ParameterSet.Config as cms
import yaml
from FWCore.ParameterSet.VarParsing import VarParsing
from SRothman.JetSetup.eecJets_cff import setupJets

opts = VarParsing('analysis')
opts.register('config', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'analysis configuration YAML')
opts.parseArguments()
with open(opts.config) as f:
    colls = yaml.safe_load(f)['jets']['collections']

process = cms.Process('JETVAL')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring(opts.inputFiles))
process.maxEvents.input = opts.maxEvents
task = cms.Task()

reco, gen = copy.deepcopy(colls['Reco']), copy.deepcopy(colls['Gen'])
for spec in (reco, gen):
    spec['jetCuts'] = dict(spec['jetCuts'], ptMin=0.0, absEtaMax=1e9, idMask=0, minConstituents=0)
gen['translator']['ptMin'] = 150.0
# Exaggerated variations so that checkJets.py can see their effect.
reco['variations'] = {'JESUp': dict(type='jes', direction='up'), 'JERUp': dict(type='jer', direction='up'),
                      'DropCharged': dict(type='constDrop', categories=['chargedHadron'], dropProbability=0.5, seedSalt=1),
                      'ScalePhoton': dict(type='constScale', categories=['photon'], scale=1.1),
                      'TightId': dict(type='cuts', jetCuts=dict(idMask=2))}
recluster = dict(mode='recluster', mcOnly=False, src='packedPFCandidates', jerc=None, variations={},
                 translator=dict(ptMin=150.0, algorithm='antikt', applyPuppiWeights=True, puppiValueMap='', dropAbsPdgIds=[]))
recoRc = dict(copy.deepcopy(reco), **recluster, table=dict(name='RecoAK8Rc', precision=-1))
genRc = dict(copy.deepcopy(gen), **dict(recluster, src='packedGenParticles',
                                        translator=dict(recluster['translator'], applyPuppiWeights=False,
                                                        dropAbsPdgIds=[12, 14, 16])), table=dict(name='GenAK8Rc', precision=-1))
reco['table']['name'], gen['table']['name'] = 'RecoAK8', 'GenAK8'
for label, spec in (('recoAK8', reco), ('genAK8', gen), ('recoAK8Rc', recoRc), ('genAK8Rc', genRc)):
    setupJets(process, label, spec, isMC=True, writeConstituents=True, task=task)

process.out = cms.OutputModule('NanoAODOutputModule', fileName=cms.untracked.string(opts.outputFile),
                               outputCommands=cms.untracked.vstring('drop *', 'keep nanoaodFlatTable_*Table_*_*'))
process.e = cms.EndPath(process.out, task)
