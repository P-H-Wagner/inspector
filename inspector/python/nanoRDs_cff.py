from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *   #why do we need all these nanoaod functions?
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *

# gen matching 
from rds.inspector.inspector_cff import *
from rds.inspector.inspectorGEN_cff import *

#G: nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectBParkTables + l1bits)  #purpose?

# from PhysiscsTools.NanoAOD
nanoSequence = cms.Sequence(nanoMetadata ) #+ globalTables)

def nanoAOD_customizeGenMatching(process):
    process.nanoGenMatchingSequence = cms.Sequence( process.nanoSequence + inspectorSequence)
    return process
def nanoAOD_customizeGENMatching(process):
    process.nanoGENMatchingSequence = cms.Sequence( process.nanoSequence + inspectorGENSequence)
    return process

