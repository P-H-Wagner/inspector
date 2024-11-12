import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from rds.inspector.common_cff import RDsCandVars, ufloat, uint, ubool


#BuilderDefaultCfg = cms.PSet(
#    dimuons             = cms.InputTag('JpsiMuonPairs','muonPairsForB'),
#    #pvSelected = cms.InputTag('pvSelector', 'bestVertex'),
#    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    #muons            = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
#    muonsTransientTracks = JpsiMuonPairs.transientTracksSrc,
#    #kaons                 = cms.InputTag('tracksBPark', 'SelectedTracks'),
#    #kaonsTransientTracks  = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
#    beamSpot              = cms.InputTag("offlineBeamSpot"),
#    tracks                = cms.InputTag("packedPFCandidates"),
#    lostTracks            = cms.InputTag("lostTracks"),
#    #kaonSelection         = cms.string(''),
#    isoTracksSelection    = cms.string('pt > 0.5 && abs(eta)<2.5'),
#    preVtxSelection       = cms.string(''),
#    postVtxSelection      = cms.string(' && '.join([
#        'userFloat("fitted_cos_theta_2D") >= 0',
#        'mass < 10.',
#        'userInt("sv_OK") == 1',
#        'userFloat("sv_prob") > 1e-8',
#        ])
#    ),
#    bits                  = cms.InputTag("TriggerResults","","HLT"),               
#    objects               = cms.InputTag("slimmedPatTrigger"), 
#)

TableDefault = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src       = cms.InputTag("your_cand_name"),
    cut       = cms.string(""),
    name      = cms.string("your_cand_name"),
    doc       = cms.string("your_cand_name Variable"),
    singleton = cms.bool(False),
    extension = cms.bool(False),
    variables = cms.PSet(),
)

#here we can add more!
TableDefaultVariables = cms.PSet(
    RDsCandVars,
    #sumpt = ufloat("sumpt")
)

#builder for final states with 3 particles
Final3PartTableVariables = TableDefaultVariables.clone(
    kIdx     = uint('k_idx'),
    bodies3_fit_k_pt    = ufloat('fitted_k_pt'),
    bodies3_fit_k_eta   = ufloat('fitted_k_eta'),
    bodies3_fit_k_phi   = ufloat('fitted_k_phi'),
    k_iso03     = ufloat('k_iso03'),
    k_iso04     = ufloat('k_iso04'),
    E_mu_star   = ufloat('E_mu_star'),
    E_mu_canc   = ufloat('E_mu_#'),
    n_k_used = uint('n_k_used'),
    ip3D_pvjpsi = ufloat('ip3D_pvjpsi'),
    ip3D_pvjpsi_e = ufloat('ip3D_pvjpsi_e'),


    #dz and dxy  for muon particle w.r.t. best pv.

    k_pvjpsi_dxy = ufloat('k_pvjpsi_dxy'),
    k_pvjpsi_dz = ufloat('k_pvjpsi_dz'),
    k_pvjpsi_dxyErr = ufloat('k_pvjpsi_dxyErr'),
    k_pvjpsi_dzErr = ufloat('k_pvjpsi_dzErr'),
)

#builder for final states with 3 muons
Final3MuonsTableVariables = Final3PartTableVariables.clone(
    ip3D_pvb = ufloat('ip3D_pvb'),
    ip3D_pvb_e = ufloat('ip3D_pvb_e'),

    k_pvb_dxy = ufloat('k_pvb_dxy'),
    k_pvb_dz = ufloat('k_pvb_dz'),
    k_pvb_dxyErr = ufloat('k_pvb_dxyErr'),
    k_pvb_dzErr = ufloat('k_pvb_dzErr'),

    pvb_idx = uint('pvb_idx'),
    pvb_x = ufloat('pvb_x'),
    pvb_y = ufloat('pvb_y'),
    pvb_z = ufloat('pvb_z'),
    pvb_ex = ufloat('pvb_ex'),
    pvb_ey = ufloat('pvb_ey'),
    pvb_ez = ufloat('pvb_ez'),
    pvb_exy = ufloat('pvb_exz'),
    pvb_eyz = ufloat('pvb_eyz'),
    pvb_exz = ufloat('pvb_exz'),
    pvb_chi2 = ufloat('pvb_chi2'),

    mu1_pvb_dxy = ufloat('mu1_pvb_dxy'),
    mu1_pvb_dz = ufloat('mu1_pvb_dz'),
    mu2_pvb_dxy = ufloat('mu2_pvb_dxy'),
    mu2_pvb_dz = ufloat('mu2_pvb_dz'),

    mu1_pvb_dxyErr = ufloat('mu1_pvb_dxyErr'),
    mu1_pvb_dzErr = ufloat('mu1_pvb_dzErr'),
    mu2_pvb_dxyErr = ufloat('mu2_pvb_dxyErr'),
    mu2_pvb_dzErr = ufloat('mu2_pvb_dzErr'),

    ip3D_pvfirst = ufloat('ip3D_pvfirst'),
    ip3D_pvfirst_e = ufloat('ip3D_pvfirst_e'),

    k_pvfirst_dxy = ufloat('k_pvfirst_dxy'),
    k_pvfirst_dz = ufloat('k_pvfirst_dz'),
    k_pvfirst_dxyErr = ufloat('k_pvfirst_dxyErr'),
    k_pvfirst_dzErr = ufloat('k_pvfirst_dzErr'),

    pvfirst_idx = uint('pvfirst_idx'),
    pvfirst_x = ufloat('pvfirst_x'),
    pvfirst_y = ufloat('pvfirst_y'),
    pvfirst_z = ufloat('pvfirst_z'),
    pvfirst_ex = ufloat('pvfirst_ex'),
    pvfirst_ey = ufloat('pvfirst_ey'),
    pvfirst_ez = ufloat('pvfirst_ez'),
    pvfirst_exy = ufloat('pvfirst_exz'),
    pvfirst_eyz = ufloat('pvfirst_eyz'),
    pvfirst_exz = ufloat('pvfirst_exz'),
    pvfirst_chi2 = ufloat('pvfirst_chi2'),

    mu1_pvfirst_dxy = ufloat('mu1_pvfirst_dxy'),
    mu1_pvfirst_dz = ufloat('mu1_pvfirst_dz'),
    mu2_pvfirst_dxy = ufloat('mu2_pvfirst_dxy'),
    mu2_pvfirst_dz = ufloat('mu2_pvfirst_dz'),

    mu1_pvfirst_dxyErr = ufloat('mu1_pvfirst_dxyErr'),
    mu1_pvfirst_dzErr = ufloat('mu1_pvfirst_dzErr'),
    mu2_pvfirst_dxyErr = ufloat('mu2_pvfirst_dxyErr'),
    mu2_pvfirst_dzErr = ufloat('mu2_pvfirst_dzErr'),
)
