"""Runs every JetSetup translator mode and a set of variations on MiniAOD and writes NANO tables.
Usage: cmsRun validateJets_cfg.py inputFiles=file:... outputFile=jets_nano.root maxEvents=500"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from SRothman.JetSetup.eecJets_cff import setupJets

opts = VarParsing('analysis')
opts.parseArguments()

process = cms.Process('JETVAL')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring(opts.inputFiles))
process.maxEvents.input = opts.maxEvents
task = cms.Task()

setupJets(process, 'recoAK8', level='reco', mode='pat', ptMin=180, tableName='RecoAK8', writeConstituents=True,
          jetCuts=dict(minConstituents=0),
          variations={'JESUp': dict(type='jes', direction='up'), 'JERUp': dict(type='jer', direction='up'),
                      'DropCharged': dict(type='constDrop', categories=['chargedHadron'], dropProbability=0.5),
                      'ScalePhoton': dict(type='constScale', categories=['photon'], scale=1.1),
                      'TightId': dict(type='cuts', jetCuts=dict(idMask=2))},
          task=task)
setupJets(process, 'genAK8', level='gen', mode='pat', ptMin=150, tableName='GenAK8', writeConstituents=True,
          jetCuts=dict(minConstituents=0), task=task)
setupJets(process, 'recoAK8Rc', level='reco', mode='recluster', ptMin=150, tableName='RecoAK8Rc',
          writeConstituents=True, jetCuts=dict(minConstituents=0), task=task)
setupJets(process, 'genAK8Rc', level='gen', mode='recluster', ptMin=150, tableName='GenAK8Rc',
          writeConstituents=True, jetCuts=dict(minConstituents=0), task=task)

process.out = cms.OutputModule('NanoAODOutputModule', fileName=cms.untracked.string(opts.outputFile),
                               outputCommands=cms.untracked.vstring('drop *', 'keep nanoaodFlatTable_*Table_*_*'))
process.e = cms.EndPath(process.out, task)
