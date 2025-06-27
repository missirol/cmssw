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

process.hltScoutingEcalBarrelRecHitPacker = cms.EDProducer("HLTScoutingEcalRecHitProducer",
  src = cms.InputTag('hltEcalRecHit:EcalRecHitsEB'),
  minEnergy = cms.double(0.5),
  mantissaPrecision = cms.int32(10),
)

process.hltScoutingEcalEndcapRecHitPacker = cms.EDProducer("HLTScoutingEcalRecHitProducer",
  src = cms.InputTag('hltEcalRecHit:EcalRecHitsEE'),
  minEnergy = cms.double(0.5),
  mantissaPrecision = cms.int32(10),
)

process.hltScoutingHBHERecHitPacker = cms.EDProducer("HLTScoutingHBHERecHitProducer",
  src = cms.InputTag('hltHbhereco'),
  minEnergy = cms.double(0.5),
  mantissaPrecision = cms.int32(10),
)

process.hltScoutingHORecHitPacker = cms.EDProducer("HLTScoutingHORecHitProducer",
  src = cms.InputTag('hltHoreco'),
  minEnergy = cms.double(0.5),
  mantissaPrecision = cms.int32(10),
)

process.hltScoutingHFRecHitPacker = cms.EDProducer("HLTScoutingHFRecHitProducer",
  src = cms.InputTag('hltHfreco'),
  minEnergy = cms.double(0.5),
  mantissaPrecision = cms.int32(10),
)

process.HLTPFScoutingPackingSequence.insert(0, process.hltScoutingEcalBarrelRecHitPacker)
process.HLTPFScoutingPackingSequence.insert(1, process.hltScoutingEcalEndcapRecHitPacker)
process.HLTPFScoutingPackingSequence.insert(2, process.hltScoutingHBHERecHitPacker)
process.HLTPFScoutingPackingSequence.insert(3, process.hltScoutingHORecHitPacker)
process.HLTPFScoutingPackingSequence.insert(4, process.hltScoutingHFRecHitPacker)

process.hltOutputScoutingPF.outputCommands += [
    'keep *_hltScoutingEcalBarrelRecHitPacker_*_*',
    'keep *_hltScoutingEcalEndcapRecHitPacker_*_*',
    'keep *_hltScoutingHBHERecHitPacker_*_*',
    'keep *_hltScoutingHORecHitPacker_*_*',
    'keep *_hltScoutingHFRecHitPacker_*_*',
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
