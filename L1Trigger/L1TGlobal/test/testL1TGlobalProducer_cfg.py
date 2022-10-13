import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')

process.options.numberOfThreads = 1
process.options.numberOfStreams = 0
process.options.wantSummary = False
process.maxEvents.input = 1000

# Global Tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

# Input source
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_12_6_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/125X_mcRun3_2022_realistic_v3-v1/2580000/2d96539c-b321-401f-b7b2-51884a5d421f.root',
    )
)

# EventSetup modules
process.GlobalParametersRcdSource = cms.ESSource('EmptyESSource',
    recordName = cms.string('L1TGlobalParametersRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.GlobalParameters = cms.ESProducer('StableParametersTrivialProducer',
    # trigger decision
    NumberPhysTriggers = cms.uint32(512), # number of physics trigger algorithms
    # trigger objects
    NumberL1Muon = cms.uint32(8),    # muons
    NumberL1EGamma = cms.uint32(12), # e/gamma and isolated e/gamma objects
    NumberL1Jet = cms.uint32(12),    # jets
    NumberL1Tau = cms.uint32(12),    # taus
    # hardware
    NumberChips = cms.uint32(1),  # number of maximum chips defined in the xml file
    PinsOnChip = cms.uint32(512), # number of pins on the GTL condition chips
    # correspondence 'condition chip - GTL algorithm word' in the hardware
    # e.g.: chip 2: 0 - 95;  chip 1: 96 - 128 (191)
    OrderOfChip = cms.vint32(1),
)

process.L1TUtmTriggerMenuRcdSource = cms.ESSource('EmptyESSource',
    recordName = cms.string('L1TUtmTriggerMenuRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.L1TriggerMenu = cms.ESProducer('L1TUtmTriggerMenuESProducer',
    L1TriggerMenuFile = cms.string('test/Collisions2022_FracPrescaleTest_L1Menu.xml'),
)

process.L1TGlobalPrescalesVetosFractRcdSource = cms.ESSource('EmptyESSource',
    recordName = cms.string('L1TGlobalPrescalesVetosFractRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)

process.L1TGlobalPrescalesVetosFract = cms.ESProducer('L1TGlobalPrescalesVetosFractESProducer',
    TriggerMenuLuminosity = cms.string('startup'),
    Verbosity = cms.int32(0),
    AlgoBxMaskDefault = cms.int32(1),
    PrescaleXMLFile   = cms.string('test/Collisions2022_FracPrescaleTest_UGT_BASE_RS_Prescales.xml'),
    AlgoBxMaskXMLFile = cms.string('test/Collisions2022_FracPrescaleTest_UGT_BASE_RS_AlgoBxMask.xml'),
    FinOrMaskXMLFile  = cms.string('test/Collisions2022_FracPrescaleTest_UGT_BASE_RS_FinOrMask.xml'),
    VetoMaskXMLFile   = cms.string('test/Collisions2022_FracPrescaleTest_UGT_BASE_RS_VetoMask.xml'),
)

# EventData modules
process.simGtExtFakeStage2Digis = cms.EDProducer('L1TExtCondProducer',
    bxFirst = cms.int32(-2),
    bxLast = cms.int32(2),
    setBptxAND = cms.bool(True),
    setBptxMinus = cms.bool(True),
    setBptxOR = cms.bool(True),
    setBptxPlus = cms.bool(True),
    tcdsRecordLabel = cms.InputTag('')
)

process.simGtStage2Digis = cms.EDProducer('L1TGlobalProducer',
    AlgoBlkInputTag = cms.InputTag('gtStage2Digis'),
    AlgorithmTriggersUnmasked = cms.bool(False),
    AlgorithmTriggersUnprescaled = cms.bool(False),
    EGammaInputTag = cms.InputTag('simCaloStage2Digis'),
    EtSumInputTag = cms.InputTag('simCaloStage2Digis'),
    ExtInputTag = cms.InputTag('simGtExtFakeStage2Digis'),
    GetPrescaleColumnFromData = cms.bool(False),
    JetInputTag = cms.InputTag('simCaloStage2Digis'),
    MuonInputTag = cms.InputTag('simGmtStage2Digis'),
    MuonShowerInputTag = cms.InputTag('simGmtShowerDigis'),
    PrescaleSet = cms.uint32(2),
    RequireMenuToMatchAlgoBlkInput = cms.bool(False),
    TauInputTag = cms.InputTag('simCaloStage2Digis'),
    useMuonShowers = cms.bool(True)
)

# Task definition
process.l1tTask = cms.Task( process.simGtExtFakeStage2Digis, process.simGtStage2Digis )

# Path definition
process.l1tPath = cms.Path( process.l1tTask )

# Analyser of L1T-menu results
process.l1tGlobalSummary = cms.EDAnalyzer( 'L1TGlobalSummary',
    AlgInputTag = cms.InputTag( 'simGtStage2Digis' ),
    ExtInputTag = cms.InputTag( 'simGtStage2Digis' ),
    MinBx = cms.int32( 0 ),
    MaxBx = cms.int32( 0 ),
    DumpTrigResults = cms.bool( False ),
    DumpRecord = cms.bool( False ),
    DumpTrigSummary = cms.bool( True ),
    ReadPrescalesFromFile = cms.bool( False ),
    psFileName = cms.string( '' ),
    psColumn = cms.int32( 0 )
)

# EndPath definition
process.l1tEndPath = cms.EndPath( process.l1tGlobalSummary )

# MessageLogger
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100 # only report every 100th event start
process.MessageLogger.L1TGlobalSummary = cms.untracked.PSet()
