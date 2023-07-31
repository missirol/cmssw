import FWCore.ParameterSet.Config as cms

def createProcess(processName):
    from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
    return cms.Process(processName, Run3_2023)

def load_cff(process, cff):
    process.load(cff)
    return process

nDiff = 0

for hltMenu in ['GRun', 'HIon', 'PRef', 'PIon', 'Fake', 'Fake1', 'Fake2']:
    p1 = createProcess('TMP1')
    hlt_cff = f'HLTrigger.Configuration.HLT_{hltMenu}_cff'
    p1.load(hlt_cff)

    print(f'Process #1 loaded "{hlt_cff}"')

    for other_cff in [
        'Configuration.StandardSequences.GeometryRecoDB_cff',
        'Configuration.StandardSequences.MagneticField_cff',
        'Configuration.StandardSequences.Digi_cff',
        'Configuration.StandardSequences.DigiToRaw_cff',
        'Configuration.StandardSequences.RawToDigi_cff',
        'Configuration.StandardSequences.SimL1Emulator_cff',
        'Configuration.StandardSequences.L1Reco_cff',
        'Configuration.StandardSequences.Reconstruction_cff',
        'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff',
    ]:
        p2 = createProcess('TMP2')
        p2.load(other_cff)
        print(f'Process #2 loaded "{other_cff}"')

        # ESSources
        for ess1, mod1 in p1.es_sources_().items():
            if ess1 in p2.es_sources_():
                mod2 = p2.es_sources_()[ess1]
                print(f'  ESSource: "{ess1}"')
                if mod1.dumpPython() != mod2.dumpPython():
                    log_str = f'the ESSource "{ess1}" loaded via "{hlt_cff}" (A)'
                    log_str += f' is different from the one loaded via "{other_cff}" (B)'
                    print(f'\n>> ERROR: {log_str}')
                    print(f'A = {mod1.dumpPython()}')
                    print(f'B = {mod2.dumpPython()}\n')
                    nDiff += 1

        # ESProducers
        for ess1, mod1 in p1.es_producers_().items():
            if ess1 in p2.es_producers_():
                mod2 = p2.es_producers_()[ess1]
                print(f'  ESProducer: "{ess1}"')
                if mod1.dumpPython() != mod2.dumpPython():
                    log_str = f'the ESProducer "{ess1}" loaded via "{hlt_cff}" (A)'
                    log_str += f' is different from the one loaded via "{other_cff}" (B)'
                    print(f'\n>> ERROR: {log_str}')
                    print(f'A = {mod1.dumpPython()}')
                    print(f'B = {mod2.dumpPython()}\n')
                    nDiff += 1

    print('-'*50)

if nDiff > 0:
    raise SystemExit(1)
