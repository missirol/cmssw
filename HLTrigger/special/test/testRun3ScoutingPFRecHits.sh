#!/bin/bash

inputFiles=($(ls /eos/cms/store/data/Run2025C/HLTPhysics/RAW/v1/000/393/461/*/*.root))
printf -v joined '%s,' "${inputFiles[@]:0:1}"
inputFilesStr="${joined%,}"
inputFilesStr=${inputFilesStr//\/eos\/cms/}

hltGetConfiguration /dev/CMSSW_15_0_0/GRun \
  --globaltag 150X_dataRun3_HLT_v1 \
  --data \
  --no-prescale \
  --output all \
  --max-events -1 \
  --paths "*ScoutingPF*","*PFScouting*","-MC*" \
  --input "${inputFilesStr}" \
  > hlt1.py

cat <<@EOF >> hlt1.py

process.hltOutputScoutingPF.fileName = 'hlt1.root'

process.hltOutputScoutingPF.compressionAlgorithm = 'LZMA'
process.hltOutputScoutingPF.compressionLevel = 4

process.options.wantSummary = False
process.options.numberOfThreads = 1
process.options.numberOfStreams = 0

del process.MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

streamPaths = [foo for foo in process.endpaths_() if foo.endswith('Output') and foo != 'ScoutingPFOutput']
for foo in streamPaths:
    process.__delattr__(foo)
@EOF

cp hlt1.py hlt2.py
cat <<@EOF >> hlt2.py

process.hltOutputScoutingPF.fileName = 'hlt2.root'

process.hltScoutingPFRecHitPackerECAL = cms.EDProducer("HLTScoutingPFRecHitProducer",
  src = cms.InputTag('hltParticleFlowRecHitECALUnseeded'),
  minEnergy = cms.double(-1),
  mantissaPrecision = cms.int32(10),
)

process.hltScoutingPFRecHitPackerHBHE = cms.EDProducer("HLTScoutingPFRecHitProducer",
  src = cms.InputTag('hltParticleFlowRecHitHBHE'),
  minEnergy = cms.double(-1),
  mantissaPrecision = cms.int32(10),
)

process.hltScoutingPFRecHitPackerHF = cms.EDProducer("HLTScoutingPFRecHitProducer",
  src = cms.InputTag('hltParticleFlowRecHitHF'),
  minEnergy = cms.double(-1),
  mantissaPrecision = cms.int32(10),
)

process.HLTPFScoutingPackingSequence.insert(0, process.hltScoutingPFRecHitPackerECAL)
process.HLTPFScoutingPackingSequence.insert(1, process.hltScoutingPFRecHitPackerHBHE)
process.HLTPFScoutingPackingSequence.insert(2, process.hltScoutingPFRecHitPackerHF)

process.hltOutputScoutingPF.outputCommands += [
    'keep *_hltScoutingPFRecHitPackerECAL_*_*',
    'keep *_hltScoutingPFRecHitPackerHBHE_*_*',
    'keep *_hltScoutingPFRecHitPackerHF_*_*',
]
@EOF

echo "=================================="
echo " hlt1 (baseline)"
echo "=================================="
#cmsRun hlt1.py 2>&1 | tee hlt1.log

echo "=================================="
echo " hlt2 (target)"
echo "=================================="
cmsRun hlt2.py 2>&1 | tee hlt2.log
