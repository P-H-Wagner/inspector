from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms
#from BcToJPsiMuMu_2021_files import files
#from HbToJPsiMuMu_2021_files import files
import os

channel = 'HOOK_CHANNEL'
print('========> Processing the miniAOD from: HOOK_INPUT')

#globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
globaltag = '106X_upgrade2018_realistic_v11_L1v1'

#what's the purpose of this 
annotation = '%s nevts:%d' % ('file:/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanotest.root', 100)

from Configuration.StandardSequences.Eras import eras
process = cms.Process('RDsNANO',eras.Run2_2018)

# import of standard configurations (do we need all of them?)
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('rds.inspector.nanoRDs_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


#prints the time report
process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True),
                             useJobReport = cms.untracked.bool(True)
)

#load all the chosen options
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(HOOK_N_EVENTS)
)

# Input source

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(HOOK_FILE_IN),
    secondaryFileNames = cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

#purpose?
process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:HOOK_FILE_OUT' ),
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')


from rds.inspector.nanoRDs_cff import *


if HOOK_GENSIM:
  process = nanoAOD_customizeGENSIMMatching(process)
  process.nanoAOD_Bs_step= cms.Path(process.nanoGENSIMMatchingSequence)

else: 

  process = nanoAOD_customizeGenMatching(process)
  process.nanoAOD_Bs_step= cms.Path(process.nanoGenMatchingSequence  )


process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
    process.nanoAOD_Bs_step,
    process.endjob_step, 
    process.NANOAODoutput_step
                               )

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('nanoAOD_Bs_step')
)

## multicore, I hope this does not screw everything up
## must be consistent with the cpu number of the batch submission!
#process.options.numberOfThreads=cms.untracked.uint32(8)
#process.options.numberOfStreams=cms.untracked.uint32(0)


# ?? 
### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
