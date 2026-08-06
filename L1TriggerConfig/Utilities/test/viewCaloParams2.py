import FWCore.ParameterSet.Config as cms

import argparse
import sys

parser = argparse.ArgumentParser(prog=sys.argv[0],
    description='Print to stderr the content of a l1t::CaloParams object from the EventSetup (record: L1TCaloParamsRcd)',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument('-r', '--runNumber', type=int, default=1,
    help='Run number to identify the IOV if the l1t::CaloParams object is taken from a GlobalTag')

group = parser.add_mutually_exclusive_group(required = True)
group.add_argument('-g', '--globaltag', dest='globaltag', type=str, default=None,
                   help = 'Name of the GlobalTag')

group.add_argument('-t', '--tag', dest='tag', type=str, default=None,
                   help = 'Name of the L1TCaloParams tag to be read from the conditions database')

group.add_argument('-c', '--caloParams-cfi', dest='caloParams_cfi', type=str, default=None,
                   help = 'Argument of process.load to load an instance of the "L1TCaloStage2ParamsESProducer" plugin')

args = parser.parse_args()

process = cms.Process("TEST")

process.maxEvents.input = 1

from FWCore.Modules.modules import EmptySource
process.source = EmptySource(firstRun = args.runNumber)

def add_globaltag(globaltag):
    from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = globaltag)

if args.globaltag:
    add_globaltag(args.globaltag)

elif args.tag:
    add_globaltag('161X_dataRun3_HLT_v1')
    process.GlobalTag.toGet += [
        cms.PSet(
            record = cms.string("L1TCaloParamsRcd"),
            tag = cms.string(args.tag),
        )
    ]

else:
    process.load(args.caloParams_cfi)

process.load("L1TriggerConfig.Utilities.l1tCaloParamsViewer2_cfi")
process.Path = cms.Path(process.l1tCaloParamsViewer2)
