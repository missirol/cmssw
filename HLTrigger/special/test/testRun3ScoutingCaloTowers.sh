#!/bin/bash

hltGetConfiguration /dev/CMSSW_15_0_0/GRun \
  --globaltag 150X_dataRun3_HLT_v1 \
  --data \
  --no-prescale \
  --output all \
  --max-events 10 \
  --paths "*ScoutingPF*","*PFScouting*","-MC*" \
  --input /store/data/Run2025C/HLTPhysics/RAW/v1/000/393/461/00000/836f8873-8791-4c83-8a4f-d475d676c7a9.root \
  > hlt.py

cat <<@EOF >> hlt.py

process.hltScoutingCaloTowerPacker = cms.EDProducer("HLTScoutingCaloTowerProducer",
  src = cms.InputTag('hltTowerMakerForAll'),
  mantissaPrecision = cms.int32(10),
)

process.HLTPFScoutingPackingSequence.insert(0, process.hltScoutingCaloTowerPacker)

process.hltOutputScoutingPF.outputCommands += [
    'keep *_hltScoutingCaloTowerPacker_*_*',
]

process.options.wantSummary = False
process.options.numberOfThreads = 1
process.options.numberOfStreams = 0

del process.MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

streamPaths = [foo for foo in process.endpaths_() if foo.endswith('Output') and foo != 'ScoutingPFOutput']
for foo in streamPaths:
    process.__delattr__(foo)

del process.dqmOutput
@EOF

cmsRun hlt.py 2>&1 | tee hlt.log
