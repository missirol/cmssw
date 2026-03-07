import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.l1scoutingrun3_cff import *

from Configuration.Eras.Modifier_run3_l1scouting_2026_cff import run3_l1scouting_2026

l1scoutingNanoTask = cms.Task(
    l1scoutingMuonPhysicalValueMap,
    l1scoutingEGammaPhysicalValueMap,
    l1scoutingTauPhysicalValueMap,
    l1scoutingJetPhysicalValueMap,
    l1scoutingMuonTable,
    l1scoutingEGammaTable,
    l1scoutingTauTable,
    l1scoutingJetTable,
    l1scoutingEtSumTable,
    l1scoutingBMTFStubTable,
)

_l1scoutingNanoTask = l1scoutingNanoTask.copy()
_l1scoutingNanoTask.add(l1scoutingCaloTowerPhysicalValueMap)
_l1scoutingNanoTask.add(l1scoutingCaloTowerTable)
run3_l1scouting_2026.toReplaceWith(l1scoutingNanoTask, _l1scoutingNanoTask)

l1scoutingNanoSequence = cms.Sequence(l1scoutingNanoTask)

def _getOutputModuleLabels(process, outputModuleType):
    return [outModLabel for outModLabel in process.outputModules_() \
        if process.outputModules_()[outModLabel].type_() == outputModuleType]

def _getPoolOutputModuleLabels(process):
    return _getOutputModuleLabels(process, 'PoolOutputModule')

def _getNanoAODOutputModuleLabels(process):
    return _getOutputModuleLabels(process, 'NanoAODOutputModule')

def _getOrbitNanoAODOutputModuleLabels(process):
    return _getOutputModuleLabels(process, 'OrbitNanoAODOutputModule')

###
### Customisation to run on the "L1Scouting" primary dataset
###
def customiseL1ScoutingNanoAOD(process):
    # NANO: convert instances of NanoAODOutputModule to instances of OrbitNanoAODOutputModule
    nanoAODOutputModuleLabels = _getNanoAODOutputModuleLabels(process)
    for outModLabel in nanoAODOutputModuleLabels:
        outMod = getattr(process, outModLabel)
        setattr(process, outModLabel, cms.OutputModule("OrbitNanoAODOutputModule",
            **outMod.parameters_(),
            skipEmptyBXs = cms.bool(True), # drop empty BXs
            selectedBx = cms.InputTag('')
        ))

    # NANO and NANOEDM: customise the event content
    poolOutputModuleLabels = _getPoolOutputModuleLabels(process)
    for outModLabel in (nanoAODOutputModuleLabels + poolOutputModuleLabels):
        outMod = getattr(process, outModLabel)
        outMod.outputCommands = [
            "drop *",
            "keep l1ScoutingRun3OrbitFlatTable_*_*_*",
        ]

    return process

###
### Customisation to run on the "L1ScoutingSelection" primary dataset
###
def customiseL1ScoutingNanoAODSelection(process):
    process = customiseL1ScoutingNanoAOD(process)

    # change input collections from the L1SCOUT data tier
    process.l1scoutingMuonPhysicalValueMap.src = "FinalBxSelectorMuon:Muon"
    process.l1scoutingEGammaPhysicalValueMap.src = "FinalBxSelectorEGamma:EGamma"
    process.l1scoutingJetPhysicalValueMap.src = "FinalBxSelectorJet:Jet"
    process.l1scoutingCaloTowerPhysicalValueMap.src = "FinalBxSelectorCaloTower:CaloTower"

    process.l1scoutingMuonTable.src = "FinalBxSelectorMuon:Muon"
    process.l1scoutingEGammaTable.src = "FinalBxSelectorEGamma:EGamma"
    process.l1scoutingJetTable.src = "FinalBxSelectorJet:Jet"
    process.l1scoutingEtSumTable.src = "FinalBxSelectorBxSums:EtSum"
    process.l1scoutingBMTFStubTable.src = "FinalBxSelectorBMTFStub:BMTFStub"
    process.l1scoutingCaloTowerTable.src = "FinalBxSelectorCaloTower:CaloTower"

    # do not throw an exception if CaloTowers are not present in the L1ScoutingSelection dataset
    process.l1scoutingCaloTowerTable.skipNonExistingSrc = True

    # drop L1Tau
    process.l1scoutingNanoTask.remove(process.l1scoutingTauTable)

    # NANO: customise instances of OrbitNanoAODOutputModule
    for outModLabel in _getOrbitNanoAODOutputModuleLabels(process):
        outMod = getattr(process, outModLabel)
        outMod.outputCommands += ["keep uints_*_SelBx_*"] # keep SelBx
        outMod.selectedBx = "FinalBxSelector:SelBx" # use to select products

    # NANOEDM: modify outputCommands of PoolOutputModule instances
    for outModLabel in _getPoolOutputModuleLabels(process):
        outMod = getattr(process, outModLabel)
        outMod.outputCommands += ["keep uints_*_SelBx_*"] # keep SelBx

    return process

###
### Additional customisations
###
###  - These functions are designed to be used with the --customise flag of cmsDriver.py,
###    e.g. "--customise PhysicsTools/NanoAOD/python/custom_l1scoutingrun3_cff.dropStub".
###
def addHardwareValues(process):
    # add hardware values to variables
    process.l1scoutingMuonTable.variables = cms.PSet(
        process.l1scoutingMuonTable.variables,
        l1scoutingMuonUnconvertedVariables
    )
    process.l1scoutingEGammaTable.variables = cms.PSet(
        process.l1scoutingEGammaTable.variables,
        l1scoutingCaloObjectUnconvertedVariables
    )
    process.l1scoutingTauTable.variables = cms.PSet(
        process.l1scoutingTauTable.variables,
        l1scoutingCaloObjectUnconvertedVariables
    )
    process.l1scoutingJetTable.variables = cms.PSet(
        process.l1scoutingJetTable.variables,
        l1scoutingCaloObjectUnconvertedVariables
    )
    process.l1scoutingCaloTowerTable.variables = cms.PSet(
        process.l1scoutingCaloTowerTable.variables,
        l1scoutingCaloTowerUnconvertedVariables
    )

    # EtSum uses dedicated EDProducer and can add hardware values by setting a boolean
    process.l1scoutingEtSumTable.writeHardwareValues = True

    return process

def keepHardwareValuesOnly(process):
    # first, add hardware values
    process = addHardwareValues(process)

    # remove physical values
    # currently external values are all physical values, so we can simple remove them
    process.l1scoutingMuonTable.externalVariables = cms.PSet()
    process.l1scoutingEGammaTable.externalVariables = cms.PSet()
    process.l1scoutingTauTable.externalVariables = cms.PSet()
    process.l1scoutingJetTable.externalVariables = cms.PSet()
    process.l1scoutingCaloTowerTable.externalVariables = cms.PSet()

    # EtSum uses dedicated EDProducer and can remove physical values by setting a boolean
    process.l1scoutingEtSumTable.writePhysicalValues = False

    return process

def outputMultipleEtSums(process):
    process.l1scoutingEtSumTable.singleton = False
    return process

def dropEmptyBXs(process):
    for outModLabel in _getOrbitNanoAODOutputModuleLabels(process):
        getattr(process, outModLabel).skipEmptyBXs = True
    return process

def keepEmptyBXs(process):
    for outModLabel in _getOrbitNanoAODOutputModuleLabels(process):
        getattr(process, outModLabel).skipEmptyBXs = False
    return process

def dropBMTFStub(process):
    process.l1scoutingNanoTask.remove(process.l1scoutingBMTFStubTable)
    return process
