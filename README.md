#Inspector

This inspector has the purpose to solely inspect the gen collections of a miniAOD. It does not apply any filter or have any cuts, it simply inspects and assigns a signal ID. There are two different plugins:

- inspectorGEN.cc
- inspector.cc

They only differ in the collections they access. For RECO, gen particles are stored in:
```
prunedCand = cms.InputTag("prunedGenParticles"),
packedCand = cms.InputTag("packedGenParticles") 
```
While for pure gen simulation, the gen particles are stored in:

```
genCand = cms.InputTag("genParticles")
```

Remark that in pruned/packed irrelevant gen particles are removed and thus they are smaller, the ntuplizer thuns runs much faster on RECO.

