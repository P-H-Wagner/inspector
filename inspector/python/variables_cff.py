import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from rds.inspector.common_cff import RDsCandVars, ufloat, uint, ubool

###################################################
## Here we define all variables we want to store ##
## in the output file. Remark that in python, we ##
## we can not give more than 255 arguments to a  ##
## function (here PSet (..)). Therefore, we need ##
## several PSets and combine them later in the   ##
## BsToDsPhiKKPiMu_cff.py                        ##
###################################################

inspectorVariables = cms.PSet(
        ## gen particle information
        mu_px     = Var("userFloat('mu_gen_px')",float),
        mu_py     = Var("userFloat('mu_gen_py')",float),
        mu_pz     = Var("userFloat('mu_gen_pz')",float),
        mu_pt     = Var("userFloat('mu_gen_pt')",float),
        mu_eta    = Var("userFloat('mu_gen_eta')",float),
        mu_phi    = Var("userFloat('mu_gen_phi')",float),
        mu_m      = Var("userFloat('mu_gen_m')",float),
        mu_charge = Var("userFloat('mu_gen_charge')",float),
        mu_pdgid  = Var("userInt('mu_gen_pdgid')", int),

        k1_px     = Var("userFloat('k1_gen_px')",float),
        k1_py     = Var("userFloat('k1_gen_py')",float),
        k1_pz     = Var("userFloat('k1_gen_pz')",float),
        k1_pt     = Var("userFloat('k1_gen_pt')",float),
        k1_eta    = Var("userFloat('k1_gen_eta')",float),
        k1_phi    = Var("userFloat('k1_gen_phi')",float),
        k1_m   = Var("userFloat('k1_gen_m')",float),
        k1_charge = Var("userFloat('k1_gen_charge')",float),
        k1_pdgid  = Var("userInt('k1_gen_pdgid')",int),

        k2_px     = Var("userFloat('k2_gen_px')",float),
        k2_py     = Var("userFloat('k2_gen_py')",float),
        k2_pz     = Var("userFloat('k2_gen_pz')",float),
        k2_pt     = Var("userFloat('k2_gen_pt')",float),
        k2_eta    = Var("userFloat('k2_gen_eta')",float),
        k2_phi    = Var("userFloat('k2_gen_phi')",float),
        k2_m   = Var("userFloat('k2_gen_m')",float),
        k2_charge = Var("userFloat('k2_gen_charge')",float),
        k2_pdgid  = Var("userInt('k2_gen_pdgid')",int),

        pi_px     = Var("userFloat('pi_gen_px')",float),
        pi_py     = Var("userFloat('pi_gen_py')",float),
        pi_pz     = Var("userFloat('pi_gen_pz')",float),
        pi_pt     = Var("userFloat('pi_gen_pt')",float),
        pi_eta    = Var("userFloat('pi_gen_eta')",float),
        pi_phi    = Var("userFloat('pi_gen_phi')",float),
        pi_m   = Var("userFloat('pi_gen_m')",float),
        pi_charge = Var("userFloat('pi_gen_charge')",float),
        pi_pdgid  = Var("userInt('pi_gen_pdgid')",int),

        phi_px    = Var("userFloat('phi_gen_px')",float),
        phi_py    = Var("userFloat('phi_gen_py')",float),
        phi_pz    = Var("userFloat('phi_gen_pz')",float),
        phi_pt    = Var("userFloat('phi_gen_pt')",float),
        phi_eta   = Var("userFloat('phi_gen_eta')",float),
        phi_phi   = Var("userFloat('phi_gen_phi')",float),
        tv_x      = Var("userFloat('tv_x_gen')",float),
        tv_y      = Var("userFloat('tv_y_gen')",float),
        tv_z      = Var("userFloat('tv_z_gen')",float),
        phi_charge= Var("userFloat('phi_gen_charge')",float),
        phi_pdgid = Var("userInt('phi_gen_pdgid')",int),

        ds_px_gen     = Var("userFloat('ds_gen_px')",float),
        ds_py_gen     = Var("userFloat('ds_gen_py')",float),
        ds_pz_gen     = Var("userFloat('ds_gen_pz')",float),
        ds_pt_gen     = Var("userFloat('ds_gen_pt')",float),
        ds_eta_gen    = Var("userFloat('ds_gen_eta')",float),
        ds_phi_gen    = Var("userFloat('ds_gen_phi')",float),
        ds_boost  = Var("userFloat('ds_gen_boost')",float),

        sv_x      = Var("userFloat('sv_x_gen')",float),
        sv_y      = Var("userFloat('sv_y_gen')",float),
        sv_z      = Var("userFloat('sv_z_gen')",float),
        ds_charge = Var("userFloat('ds_gen_charge')",float),
        ds_pdgid  = Var("userInt('ds_gen_pdgid')",int),

        bs_px     = Var("userFloat('bs_gen_px')",float),
        bs_py     = Var("userFloat('bs_gen_py')",float),
        bs_pz     = Var("userFloat('bs_gen_pz')",float),
        bs_pt     = Var("userFloat('bs_gen_pt')",float),
        bs_eta    = Var("userFloat('bs_gen_eta')",float),
        bs_phi    = Var("userFloat('bs_gen_phi')",float),

        pv_x      = Var("userFloat('pv_x_gen')",float),
        pv_y      = Var("userFloat('pv_y_gen')",float),
        pv_z      = Var("userFloat('pv_z_gen')",float),

        scnd_pv_x      = Var("userFloat('scnd_pv_x_gen')",float),
        scnd_pv_y      = Var("userFloat('scnd_pv_y_gen')",float),
        scnd_pv_z      = Var("userFloat('scnd_pv_z_gen')",float),
        scnd_pv_idx    = Var("userInt('scnd_pv_idx_gen')",int),

        bs_charge = Var("userFloat('bs_gen_charge')",float),
        bs_pdgid  = Var("userInt('bs_gen_pdgid')",int),
        bs_boost  = Var("userFloat('b_boost_gen')",float),
        bs_boost_pt   = Var("userFloat('b_boost_gen_pt')",float),
        bs_boost_eta  = Var("userFloat('b_boost_gen_eta')",float),
        bs_boost_phi  = Var("userFloat('b_boost_gen_phi')",float),

        
        disc_is_negative = Var("userInt('disc_is_negative_gen')",int),
        disc_negativity  = Var("userFloat('disc_negativity_gen')",float),

        bs_lhcb_pt    = Var("userFloat('bs_gen_lhcb_pt')",float),
        bs_lhcb_eta   = Var("userFloat('bs_gen_lhcb_eta')",float),
        bs_lhcb_phi   = Var("userFloat('bs_gen_lhcb_phi')",float),

        fv_x          = Var("userFloat('fv_x_gen')",float),
        fv_y          = Var("userFloat('fv_y_gen')",float),
        fv_z          = Var("userFloat('fv_z_gen')",float),

        m2_miss       = ufloat('m2_miss_gen'),
        pt_miss       = ufloat('pt_miss_gen'),
        q2            = ufloat('q2_gen'),
        e_star            = ufloat('e_star_gen'),

        angMuW        = Var("userFloat('angMuWGen')",float),
        cosMuW        = Var("userFloat('cosMuWGen')",float),
        cosMuWLhcb    = Var("userFloat('cosMuWGenLhcb')",float),
        cosMuWReco1   = Var("userFloat('cosMuWGenReco1')",float),
        cosMuWReco2   = Var("userFloat('cosMuWGenReco2')",float),

        angPiK1       = Var("userFloat('angPiK1Gen')",float),
        cosPiK1       = Var("userFloat('cosPiK1Gen')",float),

        angPiK2       = Var("userFloat('angPiK2Gen')",float),
        cosPiK2       = Var("userFloat('cosPiK2Gen')",float),

        angPiDs       = Var("userFloat('angPiDsGen')",float),
        cosPiDs       = Var("userFloat('cosPiDsGen')",float),
        cosPiDsLhcb   = Var("userFloat('cosPiDsGenLhcb')",float),

        angPhiDs      = Var("userFloat('angPhiDsGen')",float),
        cosPhiDs      = Var("userFloat('cosPhiDsGen')",float),

        angPlaneBs    = Var("userFloat('angPlaneBsGen')",float),
        cosPlaneBs    = Var("userFloat('cosPlaneBsGen')",float),
 
        angPlaneDs    = Var("userFloat('angPlaneDsGen')",float),
        cosPlaneDs    = Var("userFloat('cosPlaneDsGen')",float),

        match_success = uint('gen_match_success'),

        #mu_iso_03     = Var("userFloat('mu_iso_03_gen')", float),
        #mu_iso_04     = Var("userFloat('mu_iso_04_gen')", float),
        #mu_rel_iso_03 = Var("userFloat('mu_rel_iso_03_gen')", float),
        #mu_rel_iso_04 = Var("userFloat('mu_rel_iso_04_gen')", float),

        e_gamma       = Var("userFloat('e_gamma_gen')",float),

        ## signal id
        sig               = Var("userInt('sig')",int),
        b_mother_id       = Var("userInt('b_mother_id')",int),
)



##################################################

## this is for the ED Filter
#empty = cms.PSet(arrived = Var("userInt('arrived')",int))

