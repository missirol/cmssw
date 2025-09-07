import FWCore.ParameterSet.Config as cms

import argparse
import sys

parser = argparse.ArgumentParser(prog=sys.argv[0],
    description='Test of the L1TZDCEtSumsProducer and L1TZDCEtSumsPrinter plugins')

parser.add_argument('-f', '--fileNames', dest='fileNames', nargs='+',
    default=['file:/eos/cms/store/user/cmsbuild/store/hidata/HIRun2024B/HIForward0/RAW/v1/000/388/784/00000/a277c6d8-c445-4d2b-a45c-0e74e4ed8ce8.root'],
    help='Input EDM file(s)'
)

parser.add_argument('-l', '--rawDataLabel', dest='rawDataLabel', type=str, default='rawDataRepacker',
    help="Label of the FEDRawDataCollection product to be used as input")

parser.add_argument('-g', '--globalTag', dest='globalTag', type=str, default='150X_dataRun3_HLT_v1',
    help="Name of the GlobalTag")

parser.add_argument('-n', '--maxEvents', dest='maxEvents', type=int, default=100,
    help="Max number of events to be processed")

parser.add_argument('--skipEvents', dest='skipEvents', type=int, default=0,
    help="Value of process.source.skipEvents")

parser.add_argument('-t', '--numberOfThreads', dest='numberOfThreads', type=int, default=1,
    help="Value of process.options.numberOfThreads")

parser.add_argument('-s', '--numberOfStreams', dest='numberOfStreams', type=int, default=0,
    help="Value of process.options.numberOfStreams")

args = parser.parse_args()

process = cms.Process('TEST')

process.maxEvents.input = args.maxEvents

process.options.numberOfThreads = args.numberOfThreads
process.options.numberOfStreams = args.numberOfStreams
process.options.wantSummary = False

# MessageLogger
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Global Tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag = GlobalTag(process.GlobalTag, args.globalTag, '')

# Input source
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(args.fileNames),
    skipEvents = cms.untracked.uint32(args.skipEvents)
)

# EventSetup modules
process.HcalTopologyIdealEP = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string( "" ),
    MergePosition = cms.untracked.bool( True ),
    appendToDataLabel = cms.string( "" )
)

process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string( "" )
)

process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string( "" )
)

process.zdcTopologyEP = cms.ESProducer("ZdcTopologyEP",
    appendToDataLabel = cms.string( "" )
)

# EventData modules
from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis as _gtStage2Digis
process.gtStage2Digis = _gtStage2Digis.clone(InputLabel = args.rawDataLabel)

from EventFilter.HcalRawToDigi.HcalRawToDigi_cfi import hcalDigis as _hcalDigis
process.hcalDigis = _hcalDigis.clone(InputLabel = args.rawDataLabel)

from L1Trigger.L1TZDC.l1tZDCEtSums_cfi import l1tZDCEtSums as _l1tZDCEtSums
process.l1tZDCEtSums = _l1tZDCEtSums.clone(hcalTPDigis = 'hcalDigis')

from L1Trigger.L1TZDC.l1tEtSumsPrinter_cfi import l1tEtSumsPrinter as _l1tEtSumsPrinter
process.l1tEtSumsPrinter1 = _l1tEtSumsPrinter.clone(src = 'gtStage2Digis:EtSum')
process.l1tEtSumsPrinter2 = _l1tEtSumsPrinter.clone(src = 'simCaloStage2Digis')

process.simCaloStage2Digis = cms.EDProducer("L1TStage2Layer2Producer",
    firmware = cms.int32(1),
    towerToken = cms.InputTag("simCaloStage2Layer1Digis"),
    useStaticConfig = cms.bool(False)
)

process.simCaloStage2Layer1Digis = cms.EDProducer("L1TCaloLayer1",
    ecalToken = cms.InputTag("unpackEcal","EcalTriggerPrimitives"),
    firmwareVersion = cms.int32(3),
    hcalToken = cms.InputTag("hcalDigis"),
    unpackEcalMask = cms.bool(False),
    unpackHcalMask = cms.bool(False),
    useCalib = cms.bool(True),
    useECALLUT = cms.bool(True),
    useHCALFBLUT = cms.bool(False),
    useHCALLUT = cms.bool(True),
    useHFLUT = cms.bool(True),
    useLSB = cms.bool(True),
    verbose = cms.untracked.bool(False)
)

process.unpackEcal = cms.EDProducer("EcalRawToDigi",
    DoRegional = cms.bool(False),
    FEDs = cms.vint32(
        601, 602, 603, 604, 605,
        606, 607, 608, 609, 610,
        611, 612, 613, 614, 615,
        616, 617, 618, 619, 620,
        621, 622, 623, 624, 625,
        626, 627, 628, 629, 630,
        631, 632, 633, 634, 635,
        636, 637, 638, 639, 640,
        641, 642, 643, 644, 645,
        646, 647, 648, 649, 650,
        651, 652, 653, 654
    ),
    FedLabel = cms.InputTag("listfeds"),
    InputLabel = cms.InputTag("rawDataRepacker","","@skipCurrentProcess"),
    eventPut = cms.bool(True),
    feIdCheck = cms.bool(True),
    feUnpacking = cms.bool(True),
    forceToKeepFRData = cms.bool(False),
    headerUnpacking = cms.bool(True),
    memUnpacking = cms.bool(True),
    mightGet = cms.optional.untracked.vstring,
    numbTriggerTSamples = cms.int32(1),
    numbXtalTSamples = cms.int32(10),
    orderedDCCIdList = cms.vint32(
        1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        11, 12, 13, 14, 15,
        16, 17, 18, 19, 20,
        21, 22, 23, 24, 25,
        26, 27, 28, 29, 30,
        31, 32, 33, 34, 35,
        36, 37, 38, 39, 40,
        41, 42, 43, 44, 45,
        46, 47, 48, 49, 50,
        51, 52, 53, 54
    ),
    orderedFedList = cms.vint32(
        601, 602, 603, 604, 605,
        606, 607, 608, 609, 610,
        611, 612, 613, 614, 615,
        616, 617, 618, 619, 620,
        621, 622, 623, 624, 625,
        626, 627, 628, 629, 630,
        631, 632, 633, 634, 635,
        636, 637, 638, 639, 640,
        641, 642, 643, 644, 645,
        646, 647, 648, 649, 650,
        651, 652, 653, 654
    ),
    silentMode = cms.untracked.bool(True),
    srpUnpacking = cms.bool(True),
    syncCheck = cms.bool(True),
    tccUnpacking = cms.bool(True)
)

# Path
process.ThePath = cms.Path(
    process.gtStage2Digis
  + process.hcalDigis
  + process.unpackEcal
  + process.simCaloStage2Layer1Digis
  + process.simCaloStage2Digis
  + process.l1tEtSumsPrinter1
  + process.l1tEtSumsPrinter2
)

process.CaloGeometryBuilder = cms.ESProducer( "CaloGeometryBuilder",
  SelectedCalos = cms.vstring( 'HCAL',
    'ZDC',
    'EcalBarrel',
    'EcalEndcap',
    'EcalPreshower',
    'TOWER' )
)

#process.HcalGeometryFromDBEP = cms.ESProducer( "HcalGeometryFromDBEP",
#  applyAlignment = cms.bool( False )
#)

process.HcalTrigTowerGeometryESProducer = cms.ESProducer("HcalTrigTowerGeometryESProducer")

process.CaloTPGTranscoder = cms.ESProducer("CaloTPGTranscoderULUTs",
    LUTfactor = cms.vint32(1, 2, 5, 0),
    RCTLSB = cms.double(0.25),
    ZS = cms.vint32(4, 2, 1, 0),
    hcalLUT1 = cms.FileInPath('CalibCalorimetry/CaloTPG/data/outputLUTtranscoder_physics.dat'),
    hcalLUT2 = cms.FileInPath('CalibCalorimetry/CaloTPG/data/TPGcalcDecompress2.txt'),
    ietaLowerBound = cms.vint32(1, 18, 27, 29),
    ietaUpperBound = cms.vint32(17, 26, 28, 32),
    linearLUTs = cms.bool(True),
    nominal_gain = cms.double(0.177),
    read_Ascii_Compression_LUTs = cms.bool(False),
    read_Ascii_RCT_LUTs = cms.bool(False),
    tpScales = cms.PSet(
        HBHE = cms.PSet(
            LSBQIE11 = cms.double(0.0625),
            LSBQIE11Overlap = cms.double(0.0625),
            LSBQIE8 = cms.double(0.125)
        ),
        HF = cms.PSet(
            NCTShift = cms.int32(2),
            RCTShift = cms.int32(3)
        )
    )
)

process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)

process.EcalElectronicsMappingBuilder = cms.ESProducer( "EcalElectronicsMappingBuilder" )
