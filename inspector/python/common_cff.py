import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars


def ufloat(expr, precision = -1, doc = ''):
  return Var('userFloat("%s")' % expr, 
             float, precision = precision, doc = doc)

def uint(expr, doc = ''):
  return Var('userInt("%s")' % expr, int, doc = doc)

def ubool(expr, precision = -1, doc = ''):
  return Var('userInt("%s") == 1' % expr, bool, doc = doc)


#check CandVars() at:
#https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/common_cff.py
#by default, it contains: pt, eta, phi pdgId mass and charge
#means, these variables we do not have to define, they are already defined


RDsCandVars = CandVars.clone()
RDsCandVars.charge.precision = cms.int32(-1)
RDsCandVars.eta   .precision = cms.int32(-1)
RDsCandVars.mass  .precision = cms.int32(-1)
RDsCandVars.pdgId .precision = cms.int32(-1)
RDsCandVars.phi   .precision = cms.int32(-1)
RDsCandVars.pt    .precision = cms.int32(-1)
