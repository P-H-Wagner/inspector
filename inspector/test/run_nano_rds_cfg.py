from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

# TODO: put different samples into parser (flag from command line)
channel = 'dsmu'
insp    = 'gen' # reco hammer gen

import os

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
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

def filesFromFolder(directory):
  filenames = os.listdir(directory)
  for name in filenames: print(name)
  return ['file:' + directory + filename for filename in filenames ]

def filesFromTxt(txtFile):
  with open(txtFile) as dataFiles: 
    filenames = [line for line in dataFiles]
  return ['file:' + 'root://cms-xrd-global.cern.ch//' + filename for filename in filenames ]

# Input source

if channel == 'sig':

  if insp == 'hammer':

    directory = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/sig_fragment_19_07_2024_12_53_00/"
    inputfiles = filesFromFolder(directory)
    naming = 'signals_hammer'

  if insp == 'reco':
    #directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1/' #old signals MA Thesis
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/signals/all_signals_request_21_11_23.txt' # new signals!!
    #inputfiles = filesFromFolder(directory)
    inputfiles = filesFromTxt(directory)

if channel == 'hb':

  if insp == "gen":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/hb_fragment_11_06_2024_18_56_33/' #higher stats
  else:
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1/' #hb (riccardos mini)

  inputfiles = filesFromFolder(directory)



if channel == 'bplus':

  if insp == 'gen':
    directory = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bplus_fragment_03_06_2024_19_15_19/" # old filter
    directory = "/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bplus_fragment_11_06_2024_09_03_41/" # old filter includes all decays for crosscheck
    inputfiles = filesFromFolder(directory)
    
  if insp == 'reco':
    txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bplus/bplus.txt' #data
    inputfiles = filesFromTxt(txtFile)[0:1]

if channel == 'lambdab':

  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/lambdab/lambdab.txt' #bs 
  inputfiles = filesFromTxt(directory)
  naming = 'lambdab'

if channel == 'bs':

  if insp == 'gen':
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_03_06_2024_19_15_47/'  # old filter
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_26_08_2024_16_13_34/'  # old filte but different channel
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_26_08_2024_19_27_38/'  # old filte but different channel and more stats!
    #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_06_06_2024_09_52_40/'  # new filter
    #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_10_06_2024_17_46_51/'#
    inputfiles = filesFromFolder(directory)
    naming = 'bs'

  if insp == 'reco':
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bs/bs.txt' #bs 
    inputfiles = filesFromTxt(directory)
    naming = 'bs'

if channel == 'b0':

  if insp == 'gen':
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/b0_fragment_06_06_2024_20_38_57/'  # new filter
    inputfiles = filesFromFolder(directory)[0:20]
    naming = 'b0'

  if insp == 'reco':
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bs/bs.txt' #bs 
    inputfiles = filesFromTxt(directory)
    naming = 'b0'

if channel == "dstau":
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dstau_fragment_27_09_2024_11_26_41/'  # old filter
  inputfiles = filesFromFolder(directory)
  naming = 'dstau'

if channel == "dsmu":
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsmu_fragment_15_10_2024_13_04_12/'  # old filter
  inputfiles = filesFromFolder(directory)
  naming = 'dsmu'

print(inputfiles)
inputfiles = inputfiles[0:1]


if channel == 'data':
  txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPark_2018_D/BPark_2018D.txt' #data
  inputfiles = filesFromTxt(txtFile)

process.source = cms.Source(
    "PoolSource",
    #fileNames = cms.untracked.vstring('file:/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1/all_signals_HbToDsPhiKKPiMuNu_MT_97.root'),
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/50003/199766A2-70D0-674D-A9CE-D5EB255BD87A.root'), # data file which is always empty, investigate this 
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/40000/56D888ED-EB2C-B24F-A5A0-8D162DAFFA25.root'), #data file to compare with riccs MA code
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/50002/D0AE1369-0D7B-554C-BBB9-7B324AACCABD.root'), #data file to compare with riccs MA code
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/data/Run2018D/ParkingBPH1/MINIAOD/UL2018_MiniAODv2-v1/50002/36CD4F31-A249-DF49-A3FF-32DCA7223D09.root'), #10
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18MiniAODv2/BsToDoubleCharm_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v1/2820000/FD3D3434-1FAA-D743-914D-D3F4A122FAAD.root'), #7
    fileNames = cms.untracked.vstring(inputfiles),# all_signals_HbToDsPhiKKPiMuNu_MT_0.root'), #automized case
    #fileNames = cms.untracked.vstring("file:/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_03_06_2024_19_15_47/bs_fragment_chunk3_03_06_2024_19_15_47.root"),
    #fileNames = cms.untracked.vstring('file:root://cms-xrd-global.cern.ch////store/mc/RunIISummer20UL18MiniAODv2/BsToDoubleCharm_DsFilter_PhiFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/Custom_RDStarPU_BParking_106X_upgrade2018_realistic_v16_L1v1-v1/2820000/CEF9599E-F695-554B-97D2-276663715216.root'),
    #fileNames = cms.untracked.vstring("file:/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_10_06_2024_17_46_51/bs_fragment_chunk9_10_06_2024_17_46_51.root"),

    #eventsToProcess  = cms.untracked.VEventRange('1:22955'),
    secondaryFileNames = cms.untracked.vstring(),
    #skipEvents=cms.untracked.uint32(0) # skip first n events   
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
    #fileName = cms.untracked.string('file:/scratch/pahwagne/nanoAOD/test.root' ),
    fileName = cms.untracked.string('file:inspectortest.root'), #used for local tests
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


# add all sequences as addtributes
from rds.inspector.nanoRDs_cff import *

#can only gen match on mc
if insp =="gen":
  process = nanoAOD_customizeGENMatching(process) 
  process.nanoAOD_Bs_step= cms.Path(process.nanoGENMatchingSequence)
elif insp =="reco":
  process = nanoAOD_customizeGenMatching(process) 
  # Path and EndPath definitions
  process.nanoAOD_Bs_step= cms.Path(process.nanoGenMatchingSequence)
else:
  print("Not a valid inspector type! There are 'gen' or 'reco'")


#process.nanoAOD_Bs_step= cms.Path(process.triggerSequence) ## to run only Trigger.cc for debugging
process.endjob_step = cms.EndPath(process.endOfProcess)

process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
    process.nanoAOD_Bs_step,
    process.endjob_step, 
    process.NANOAODoutput_step # commented out ? not saving !!
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
