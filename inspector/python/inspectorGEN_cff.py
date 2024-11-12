import FWCore.ParameterSet.Config as cms
#from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from rds.inspector.common_cff import RDsCandVars, ufloat, uint, ubool
#from rds.inspector.rds_common_cff import TableDefaultVariables, TableDefault
from rds.inspector.variables_cff import inspectorVariables, hammerVariables

#from PhysicsTools.RDsNano.primaryVertices_cff import *

#BsToDsPhiKKPiMuCfg = BuilderDefaultCfg.clone()
#BTo3MuCfg.dileptons             = cms.InputTag('JpsiMuonPairs')
#BTo3MuCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc

inspectorGEN = cms.EDProducer(
    'inspectorGEN',
    #bs           = cms.InputTag('BsToDsPhiKKPiMu', 'bs'),
    genCand   = cms.InputTag("genParticles"),
minMuPt     = cms.double(7.0),
maxMuEta    = cms.double(1.5),
maxdRHadMuon = cms.double(1.2),       # max dR between hadron and muon
mindRHadMuon = cms.double(0.005),     # min dR "
maxdzDiffHadMuon = cms.double(0.6),   # difference in dz between muon/pv and had/pv
phiMassAllowance = cms.double(0.030), # allow 30 MeV when collecting candidates for phi 
dsMassAllowance = cms.double(0.150),  # allow 150 MeV when collecting candidates for ds
drMatchGen = cms.double(0.1),         # allow 0.1, (0.05 would also be reasonable) in dR when gen matching
maxBsMass = cms.double(8.0),

piMass = cms.double(0.13957039),      # pi mass
kMass = cms.double(0.493677),         # kaon mass
phiMass = cms.double(1.019461),       # phi mass
dsMass = cms.double(1.96834),         # ds mass
dsStarMass = cms.double(2.112204),    # ds star mass
muMass = cms.double(0.105658),        # mu mass
tauMass = cms.double(1.77686),        # tau mass
bsMass = cms.double(5.36688),         # bs mass
isoCone = cms.double(0.5)             # cut on dR for the mu isolation cone
)

print( " ========> Parameters used:")
print(inspectorGEN.dumpPython)

#BsToDsPhiKKPiMuVariables.extend(vertexVariables)
combined_variablesGEN = cms.PSet(
  inspectorVariables,
  hammerVariables,
)
  
inspectorGENTable = cms.EDProducer('SimpleCompositeCandidateFlatTableProducer',

    src = cms.InputTag("inspectorGEN:gen"),
    cut = cms.string(""),
    name = cms.string("gen"),
    doc  = cms.string("inspectorGEN"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), 
    variables = combined_variablesGEN

)

#tables = cms.Sequence(BsToDsPhiKKPiMuTable) # * vertexTable)
inspectorGENSequence = cms.Sequence(inspectorGEN + inspectorGENTable)

