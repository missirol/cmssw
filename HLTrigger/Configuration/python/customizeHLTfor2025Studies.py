import FWCore.ParameterSet.Config as cms

from HLTrigger.Configuration.common import *

def customizeHLTfor2025HCALPFCuts(process):
    try:
        process.GlobalTag.toGet += [
            cms.PSet(
                record = cms.string('HcalPFCutsRcd'),
                tag = cms.string('HcalPFCuts_2025_mc'),
            ),
        ]
    except:
        raise RuntimeError("customizeHLTfor2025HCALPFCuts -- GlobalTag ESSource could not be customized !")

    return process

def customizeHLTfor2025PFHadronCalibrations(process):
    try:
        process.GlobalTag.toGet += [
            cms.PSet(
                record = cms.string('PFCalibrationRcd'),
                tag = cms.string('PFCalibration_Run3Winter25_MC_hlt_v1'),
                label = cms.untracked.string('HLT'),
            ),
        ]
    except:
        raise RuntimeError("customizeHLTfor2025PFHadronCalibrations -- GlobalTag ESSource could not be customized !")

    return process

def customizeHLTforCMSHLT3469(process):
    for prod in producers_by_type(process, 'CaloTowersCreator'):
        prod.EcalRecHitThresh = True
    return process

def customizeHLTfor2025JECs(process):
    jecTagsDict = {
        'AK4CaloHLT': '', # not available yet
        'AK8CaloHLT': '', # not available yet
        'AK4PFHLT': '', # not available yet
        'AK8PFHLT': '', # not available yet
    }

    try:
        for (labelName, tagName) in jecTagsDict.items():
            process.GlobalTag.toGet += [
                cms.PSet(
                    record = cms.string("JetCorrectionsRecord"),
                    label = cms.untracked.string(labelName),
                    tag = cms.string(tagName),
                ),
            ]
    except:
        raise RuntimeError("customizeHLTfor2025JECs -- GlobalTag ESSource could not be customized !")

    return process

def customizeHLTfor2025Studies(process):
    process = customizeHLTfor2025HCALPFCuts(process)
    process = customizeHLTfor2025PFHadronCalibrations(process)
#    process = customizeHLTforCMSHLT3469(process)
#    process = customizeHLTfor2025JECs(process)
    return process

def customizeHLTfor2024L1TMenu(process):
    seed_replacements = {

        'L1_SingleMu5_BMTF' : 'L1_AlwaysTrue',
        'L1_SingleMu13_SQ14_BMTF': 'L1_AlwaysTrue',

        'L1_AXO_Medium' : 'L1_AXO_Nominal',
        'L1_AXO_VVTight': 'L1_AlwaysTrue',
        'L1_AXO_VVVTight': 'L1_AlwaysTrue',

        'L1_CICADA_VVTight': 'L1_AlwaysTrue',
        'L1_CICADA_VVVTight': 'L1_AlwaysTrue',
        'L1_CICADA_VVVVTight': 'L1_AlwaysTrue',

        'L1_DoubleTau_Iso34_Iso26_er2p1_Jet55_RmOvlp_dR0p5': 'L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5 OR L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5',
        'L1_DoubleTau_Iso38_Iso26_er2p1_Jet55_RmOvlp_dR0p5': 'L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5 OR L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5',
        'L1_DoubleTau_Iso40_Iso26_er2p1_Jet55_RmOvlp_dR0p5': 'L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5 OR L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5',

        'L1_DoubleEG15_11_er1p2_dR_Max0p6': 'L1_DoubleEG11_er1p2_dR_Max0p6',
        'L1_DoubleEG16_11_er1p2_dR_Max0p6': 'L1_DoubleEG11_er1p2_dR_Max0p6',
        'L1_DoubleEG17_11_er1p2_dR_Max0p6': 'L1_DoubleEG11_er1p2_dR_Max0p6',

        'L1_DoubleJet_110_35_DoubleJet35_Mass_Min1000': 'L1_AlwaysTrue',
        'L1_DoubleJet_110_35_DoubleJet35_Mass_Min1100': 'L1_AlwaysTrue',
        'L1_DoubleJet_110_35_DoubleJet35_Mass_Min1200': 'L1_AlwaysTrue',
        'L1_DoubleJet45_Mass_Min700_IsoTau45er2p1_RmOvlp_dR0p5': 'L1_AlwaysTrue',
        'L1_DoubleJet45_Mass_Min800_IsoTau45er2p1_RmOvlp_dR0p5': 'L1_AlwaysTrue',
        'L1_DoubleJet_65_35_DoubleJet35_Mass_Min750_DoubleJetCentral50': 'L1_AlwaysTrue',
        'L1_DoubleJet_65_35_DoubleJet35_Mass_Min850_DoubleJetCentral50': 'L1_AlwaysTrue',
        'L1_DoubleJet_65_35_DoubleJet35_Mass_Min950_DoubleJetCentral50': 'L1_AlwaysTrue',
        'L1_DoubleJet45_Mass_Min700_LooseIsoEG20er2p1_RmOvlp_dR0p2': 'L1_AlwaysTrue',
        'L1_DoubleJet45_Mass_Min800_LooseIsoEG20er2p1_RmOvlp_dR0p2': 'L1_AlwaysTrue',
        'L1_DoubleJet_85_35_DoubleJet35_Mass_Min700_Mu3OQ': 'L1_AlwaysTrue',
        'L1_DoubleJet_85_35_DoubleJet35_Mass_Min800_Mu3OQ': 'L1_AlwaysTrue',
        'L1_DoubleJet_85_35_DoubleJet35_Mass_Min900_Mu3OQ': 'L1_AlwaysTrue',
        'L1_DoubleJet_70_35_DoubleJet35_Mass_Min600_ETMHF65': 'L1_AlwaysTrue',
        'L1_DoubleJet_70_35_DoubleJet35_Mass_Min700_ETMHF65': 'L1_AlwaysTrue',
        'L1_DoubleJet_70_35_DoubleJet35_Mass_Min800_ETMHF65': 'L1_AlwaysTrue',
    }

    for module in filters_by_type(process, 'HLTL1TSeed'):
        l1Seed = module.L1SeedsLogicalExpression.value()
        if any(old_seed in l1Seed for old_seed in seed_replacements):
            for old_seed, new_seed in seed_replacements.items():
                l1Seed = l1Seed.replace(old_seed, new_seed)
            module.L1SeedsLogicalExpression = cms.string(l1Seed)

    return process
