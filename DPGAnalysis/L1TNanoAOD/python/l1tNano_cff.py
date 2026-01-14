import FWCore.ParameterSet.Config as cms

from DPGAnalysis.L1TNanoAOD.l1tnanotables_cff import *
from PhysicsTools.NanoAOD.l1trig_cff import l1TablesTask
from PhysicsTools.NanoAOD.nano_cff import *

l1tNanoTask = cms.Task(nanoMetadata, l1TablesTask)

l1tNanoSequence = cms.Sequence(l1tNanoTask)

def addEmulObjects(process):

    process.l1tNanoTask.add(l1EmulObjTablesTask)
    
    return process

def addUnpackedCaloTPs(process):

    process.l1tNanoTask.add(process.l1CaloTPsNanoTask)
    
    return process

def addEmulCaloTPs(process):

    process.l1tNanoTask.add(process.l1EmulCaloTPsNanoTask)

    return process

def addUnpackedCaloLayer1(process):

    process.l1tNanoTask.add(process.l1CaloLayer1NanoTask)

    return process

def addEmulCaloLayer1(process):

    process.l1tNanoTask.add(process.l1EmulCaloLayer1NanoTask)
         
    return process

def addUnpackedCaloTPsandLayer1(process):

    addUnpackedCaloTPs(process)
    addUnpackedCaloLayer1(process)

    return process

def addEmulCaloTPsandLayer1(process):

    addEmulCaloTPs(process)
    addEmulCaloLayer1(process)

    return process

def addCaloFull(process):

    addEmulCaloTPsandLayer1(process)
    addUnpackedCaloTPsandLayer1(process)
    addEmulObjects(process)

    return process
