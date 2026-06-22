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

group.add_argument('-c', '--caloParams-cfi', dest='caloParams_cfi', type=str, default=None,
                   help = 'Argument of process.load to load an instance of the "L1TCaloStage2ParamsESProducer" plugin')

args = parser.parse_args()

process = cms.Process("TEST")

process.maxEvents.input = 1

from FWCore.Modules.modules import EmptySource
process.source = EmptySource(firstRun = args.runNumber)

if args.globaltag:
    from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = args.globaltag)
else:
    process.load(args.caloParams_cfi)

process.load("L1TriggerConfig.Utilities.l1tCaloParamsViewer2_cfi")
process.Path = cms.Path(process.l1tCaloParamsViewer2)
