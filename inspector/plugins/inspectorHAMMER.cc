#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
// for vertex fitting (both global and sequential)
#include "TrackingTools/Records/interface/TransientTrackRecord.h" 
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h" 
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h" 
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h" 
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h" 
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h" 
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicConstraint.h" 
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFitPrimitives/interface/Matrices.h" 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h" 

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // for the tracks!
#include "FWCore/Framework/interface/MakerMacros.h"
#include <limits>
#include <algorithm>
#include <cmath>
#include "helper.h" // helper functions
// include "PxPyPzMVector.h" // to new :(
#include "TLorentzVector.h" // use this instead 
#include "TVector3.h" // for boost vector
// for gen matching
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h" 

// B field
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
// for 3DPoint
#include "DataFormats/Math/interface/Point3D.h"
// for the cov matrix correction
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h" 


// counters
int nEventsHAMMER = 0;
int nMuonsHAMMER  = 0;
int nTracksHAMMER = 0;
int nGenMatchedHAMMER = 0;

int beforePhiMassHAMMER = 0;
int afterPhiMassHAMMER = 0;

int nKKPiMuHAMMER = 0;
int nFoundPhiHAMMER  = 0;
int nFoundDsHAMMER   = 0;
int nFoundBHAMMER    = 0;
int nBMassCutHAMMER   = 0;

class inspectorHAMMER : public edm::global::EDProducer<> {

public:

  //define collections which dont exist by default  
  typedef std::vector<reco::GenParticle> GenParticleCollection;
  typedef std::vector<pat::PackedGenParticle> PackedGenParticleCollection;
  //constructor
  explicit inspectorHAMMER(const edm::ParameterSet&);
  //destructor
  ~inspectorHAMMER() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  virtual void endJob() override;

  int  getPVIdx(const reco::VertexCollection*,const reco::TransientTrack&) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}


  //must be constant as it takes constant arguments!! otherwise compiler rises an error
  //because it thinks it may modify the input!!
  reco::TransientTrack getTransientTrack(const reco::Track track) const {    
      reco::TransientTrack transientTrack(track, paramField);

      return transientTrack;
    }


private:
  
  //Bfield
  OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");

  //cuts 

  const double minMuPt_;
  const double maxMuEta_;
  const double maxdRHadMuon_;
  const double mindRHadMuon_;
  const double maxdzDiffHadMuon_; 
  const double phiMassAllowance_;
  const double dsMassAllowance_;
  const double drMatchGen_;
  const double maxBsMass_;
  const double piMass_;
  const double kMass_;
  const double phiMass_;
  const double dsMass_;
  const double dsStarMass_;
  const double muMass_;
  const double tauMass_;
  const double bsMass_;
  const double isoCone_;
  //tokens to access data later
  //edm::Input tag can not be directly initialized inside the construcor! Why did it work fro Trigger.cc??
  //anyway ... 

  //for the muons

  // vertices
  const edm::InputTag primaryVtxTag;
  const edm::EDGetTokenT<reco::VertexCollection> primaryVtx_;

  //gen for gen-matching
  const edm::InputTag genTag; //pruned is a compressed packed format
  const edm::EDGetTokenT<reco::GenParticleCollection> pureGen_;

};

//define the constructor
inspectorHAMMER::inspectorHAMMER(const edm::ParameterSet& iConfig):
    // f.e. hadSelection_ = cfg.getPatameter...

    minMuPt_(iConfig.getParameter<double>("minMuPt")),
    maxMuEta_(iConfig.getParameter<double>("maxMuEta")),
    maxdRHadMuon_(iConfig.getParameter<double>("maxdRHadMuon")),
    mindRHadMuon_(iConfig.getParameter<double>("mindRHadMuon")),
    maxdzDiffHadMuon_(iConfig.getParameter<double>("maxdzDiffHadMuon")),
    phiMassAllowance_(iConfig.getParameter<double>("phiMassAllowance")),
    dsMassAllowance_(iConfig.getParameter<double>("dsMassAllowance")),
    drMatchGen_(iConfig.getParameter<double>("drMatchGen")),
    maxBsMass_(iConfig.getParameter<double>("maxBsMass")),

    piMass_(iConfig.getParameter<double>("piMass")),
    kMass_(iConfig.getParameter<double>("kMass")),
    phiMass_(iConfig.getParameter<double>("phiMass")),
    dsMass_(iConfig.getParameter<double>("dsMass")),
    dsStarMass_(iConfig.getParameter<double>("dsStarMass")),
    muMass_(iConfig.getParameter<double>("muMass")),
    tauMass_(iConfig.getParameter<double>("tauMass")),
    bsMass_(iConfig.getParameter<double>("bsMass")),
    isoCone_(iConfig.getParameter<double>("isoCone")),

    genTag(iConfig.getParameter<edm::InputTag>("genCand")),
    pureGen_(consumes<reco::GenParticleCollection>(genTag)){
       // output collection
       produces<pat::CompositeCandidateCollection>("genHAMMER");
       //produces<pat::CompositeCandidateCollection>("genHAMMER");
       //produces<TransientTrackCollection>("kkTransientTracks");
    }

//check const keywords 

// this starts the event loop
void inspectorHAMMER::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {

  //input
  edm::Handle<reco::GenParticleCollection> pureGen;
  iEvent.getByToken(pureGen_,pureGen);

  // to save 
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  //std::unique_ptr<pat::CompositeCandidateCollection> ret_value_gen(new pat::CompositeCandidateCollection());
  //std::unique_ptr<TransientTrackCollection> kkpi_ttrack(new TransientTrackCollection);

  nEventsHAMMER++;

  //////////////////////////////////////////////////////
  // Match the bsCand (and all its final states) with //
  // a gen particle. We loop over all candidates      //
  // per event!                                       //
  //////////////////////////////////////////////////////

  int sigId = -9999;
  int bId = 0;
  int genMatchSuccess = 0;

  //count the number of gen matches we find, ideally only 1
  int nGenMatches = 0;

  std::cout << "New event!" << iEvent.id().event() <<std::endl;

  ////////////////////////////////////////////////////
  // find the gen-matched muon                      //
  ////////////////////////////////////////////////////

  for(size_t muIdxGen = 0; muIdxGen < pureGen->size(); ++muIdxGen){

    nMuonsHAMMER++;

    //define a pointer to the gen muon    
    edm::Ptr<reco::GenParticle> muPtrGen(pureGen, muIdxGen);

    //select only useful gen muons -> check this selection!
    if((fabs(muPtrGen->pdgId()) != 13)) continue; // || muPtrGen->pt() < minMuPt_ || fabs(muPtrGen->eta()) > maxMuEta_ || (muBs->charge() * muPtrGen->charge() < 0)) continue; 

    ////////////////////////////////////////////////
    // find gen matched k1                        //
    ////////////////////////////////////////////////
    int kaons = 0;
    for(size_t k1IdxGen = 0; k1IdxGen < pureGen->size(); ++k1IdxGen){
         
      nTracksHAMMER++;
      //define a pointer to the gen kaon    
      edm::Ptr<reco::GenParticle> k1PtrGen(pureGen, k1IdxGen);

      //select only useful kaons -> check this selection!
      if((fabs(k1PtrGen->pdgId()) != 321)) continue; // || !hadSelectionGen_(*k1PtrGen) || (k1Bs->charge() * k1PtrGen->charge() < 0)) continue; 
      //std::cout <<"total size: " << pureGen->size()  << std::endl;
      //std::cout <<"index now:  " << k1IdxGen << std::endl;
      //std::cout <<"pt" <<  k1PtrGen->pt() << std::endl;
      kaons++;
      ////////////////////////////////////////////////
      // find gen matched k2                        //
      ////////////////////////////////////////////////
   
      for(size_t k2IdxGen = k1IdxGen + 1; k2IdxGen < pureGen->size(); ++k2IdxGen){
   
           //avoid picking the same gen particle as for k1
           if(k2IdxGen == k1IdxGen) continue; 

           //define a pointer to the gen kaon    
           edm::Ptr<reco::GenParticle> k2PtrGen(pureGen, k2IdxGen);
   
           //select only useful kaons -> check this selection!
           if((fabs(k2PtrGen->pdgId()) != 321)) continue;  // || !hadSelectionGen_(*k2PtrGen) || (k2Bs->charge() * k2PtrGen->charge() < 0 )) continue; 

           //std::cout <<"total size: " << pureGen->size()  << std::endl;
           //std::cout <<"index now:  " << k2IdxGen << std::endl;
           //std::cout <<"pt" <<  k2PtrGen->pt() << std::endl;

           beforePhiMassHAMMER++;
           // DEBUGGING!!!
           //if(k1PtrGen->charge() * k2PtrGen->charge() > 0) continue;
           //pat::CompositeCandidate kk;
           //math::PtEtaPhiMLorentzVector k1P4(k1PtrGen->pt(), k1PtrGen->eta(), k1PtrGen->phi(), K_MASS);
           //math::PtEtaPhiMLorentzVector k2P4(k2PtrGen->pt(), k2PtrGen->eta(), k2PtrGen->phi(), K_MASS);
           //kk.setP4(k1P4 + k2P4);
           //if (fabs(kk.mass() - phiMass_) > phiMassAllowance_) continue;
           //std::cout << "survived the phi mass allowance " << std::endl;
           //afterPhiMass++;
           ////////////////////////////////////////////////
           // find gen matched pion                      //
           ////////////////////////////////////////////////

           for(size_t piIdxGen = 0; piIdxGen < pureGen->size(); ++piIdxGen){
    
             //avoid picking the same gen particle as for k1 or k2
             if((piIdxGen == k1IdxGen) || (piIdxGen == k2IdxGen)) continue; 
                     
             //define a pointer to the gen kaon    
             edm::Ptr<reco::GenParticle> piPtrGen(pureGen, piIdxGen);
   
             //select only useful kaons -> check this selection!
             if((fabs(piPtrGen->pdgId()) != 211)) continue; // || !hadSelectionGen_(*piPtrGen) || (piBs->charge() * piPtrGen->charge() < 0 )) continue; 
   

             //////////////////////////////////////////////////
             // Find resonances at gen level                 //
             //////////////////////////////////////////////////
             pat::CompositeCandidate genHAMMER; 

             nKKPiMuHAMMER++;     
             //Should we pick the best gen match (in terms of dR) only? -> No, like this is better 

             const reco::Candidate* k1Reco = k1PtrGen.get(); 
             const reco::Candidate* k2Reco = k2PtrGen.get(); 
             const reco::Candidate* piReco = piPtrGen.get(); 
             const reco::Candidate* muReco = muPtrGen.get(); 
             //printMothers(k2Reco); //-> for debugging
             // searching for phi resonance 
             auto phiFromK1 = getAncestor(k1Reco,333);
             auto phiFromK2 = getAncestor(k2Reco,333);
             if( (phiFromK1 != phiFromK2) || (phiFromK1 == nullptr) || (phiFromK2 == nullptr)) continue; 
             nFoundPhiHAMMER++;               
   
             //std::cout<< "found phi candiadte!" << std::endl;
             // searching for ds resonance 
             auto dsFromPhi = getAncestor(phiFromK1,431);
             auto dsFromPi  = getAncestor(piReco,431);
             if( (dsFromPhi != dsFromPi) || (dsFromPhi == nullptr) || (dsFromPi == nullptr)) continue; 
             nFoundDsHAMMER++;               

             //std::cout << "k pts are: " << k1PtrGen->pt() << " " << k2PtrGen->pt() << std::endl;
             //printMothers(dsFromPhi);

             // we dont know what b mother we have
             int bMotherId = 0;

             // first search for b baryon (isAncestor also checks neg. Ids)
             for(int bIdx = 5000; bIdx < 6000; bIdx++){
               if (isAncestor(dsFromPhi, bIdx)){
                 bMotherId = bIdx;
                 break;
               }
             }

             // Then search for possible B meson coming from B-baryon
             for(int bIdx = 500; bIdx < 600; bIdx++){
               if (isAncestor(dsFromPhi, bIdx)){
                 bMotherId = bIdx;
                 break;
               }
             }
              
             if (bMotherId == 0) continue; // no b mother found
             // Even if the mu is not required to come brom the b mother directly (would be signal case)
             // if it comes from another D meson (double charm background case), we still want
             // that this D meson is coming from the b mother. So the muon should share
             // the same ancestor as the Ds.
             auto bsFromDs = getAncestor(dsFromPhi,bMotherId);
             auto bsFromMu = getAncestor(muReco,   bMotherId);

             if( (bsFromDs != bsFromMu) || (bsFromDs == nullptr) || (bsFromMu == nullptr)) continue; 
 
             nFoundBHAMMER++;
             
             //if (bsFromMu->mass() > maxBsMass_) continue;
             nBMassCutHAMMER++;

             //remove oscillations
             auto bsFromMuWOOsc = removeOscillations(bsFromMu);

             nGenMatches++;
             genMatchSuccess = 1;
             nGenMatchedHAMMER++;
             
             //if(nGenMatches > 1) continue; //std::cout <<"there is more than one match!!" << std::endl;

             //get gen 4 momenta
             TLorentzVector genMuTlv; 
             TLorentzVector genK1Tlv; 
             TLorentzVector genK2Tlv; 
             TLorentzVector genPiTlv; 
             TLorentzVector genPhiTlv; 
             TLorentzVector genDsTlv; 
             TLorentzVector genBsTlv; 

             TLorentzVector genMissTlv; //for m2 miss 
             TLorentzVector genQTlv;  // for q2

             genMuTlv.SetXYZM(  muPtrGen->px(),      muPtrGen->py(),      muPtrGen->pz(),      muMass_);
             genK1Tlv.SetXYZM(  k1PtrGen->px(),      k1PtrGen->py(),      k1PtrGen->pz(),      kMass_);
             genK2Tlv.SetXYZM(  k2PtrGen->px(),      k2PtrGen->py(),      k2PtrGen->pz(),      kMass_);
             genPiTlv.SetXYZM(  piPtrGen->px(),      piPtrGen->py(),      piPtrGen->pz(),      piMass_);

             genPhiTlv.SetXYZM( phiFromK1->px(),     phiFromK1->py(),     phiFromK1->pz(),     phiMass_);
             genDsTlv.SetXYZM(  dsFromPi->px(),      dsFromPi->py(),      dsFromPi->pz(),      dsMass_);
             //genBsTlv.SetXYZM(  bsFromMuWOOsc->px(), bsFromMuWOOsc->py(), bsFromMuWOOsc->pz(), bsMass_); //changed
             genBsTlv.SetXYZM(  bsFromMu->px(), bsFromMu->py(), bsFromMu->pz(), bsMass_); //changed

             genMissTlv = genBsTlv - (genDsTlv + genMuTlv); 
             genQTlv    = genBsTlv - (genDsTlv); 

             float m2_miss_gen = genMissTlv.M2();
             float pt_miss_gen = genMissTlv.Pt();
             float q2_gen = genQTlv.M2();
             float e_star_gen   = getEStar(genBsTlv,genMuTlv);
             float e_gamma_gen  = getEGamma(genDsTlv, dsMass_, dsStarMass_);

             genHAMMER.addUserFloat("m2_miss_gen",m2_miss_gen);
             genHAMMER.addUserFloat("pt_miss_gen",pt_miss_gen);
             genHAMMER.addUserFloat("e_star_gen",e_star_gen);
             genHAMMER.addUserFloat("e_gamma_gen",e_gamma_gen);
             genHAMMER.addUserFloat("q2_gen",q2_gen);

             //vertices
             float pv_x_gen = bsFromMuWOOsc->vx();//This is the bs production vertex!
             float pv_y_gen = bsFromMuWOOsc->vy();
             float pv_z_gen = bsFromMuWOOsc->vz();

             // Do a cross check on gen: Is there another primary vtx in the primary vertex collection which
             // is closer to the gen truth?


             //std::cout << dsFromPi->vtx() << std::endl;

             float sv_x_gen = dsFromPi->vx(); //This is the ds production vertex!
             float sv_y_gen = dsFromPi->vy();
             float sv_z_gen = dsFromPi->vz();

             float tv_x_gen = phiFromK1->vx(); //This is the phi production vertex!
             float tv_y_gen = phiFromK1->vy();
             float tv_z_gen = phiFromK1->vz();

             float fv_x_gen = k1PtrGen->vx(); //This is the k1 production vertex!
             float fv_y_gen = k1PtrGen->vy();
             float fv_z_gen = k1PtrGen->vz();

             //std::cout << "on gen:" << fv_x_gen << std::endl; //This is the k1 production vertex!
             //save the gen info by adding gen candidates of final states

             //gen.addUserCand("mu_gen",muPtrGen);
             //gen.addUserCand("k1_gen",k1PtrGen);
             //gen.addUserCand("k2_gen",k2PtrGen);
             //gen.addUserCand("pi_gen",piPtrGen);
             genHAMMER.addUserFloat("mu_gen_px"      ,muPtrGen->px());
             genHAMMER.addUserFloat("mu_gen_py"      ,muPtrGen->py());
             genHAMMER.addUserFloat("mu_gen_pz"      ,muPtrGen->pz());
             genHAMMER.addUserFloat("mu_gen_pt"      ,muPtrGen->pt());
             genHAMMER.addUserFloat("mu_gen_eta"     ,muPtrGen->eta());
             genHAMMER.addUserFloat("mu_gen_phi"     ,muPtrGen->phi());
             genHAMMER.addUserFloat("mu_gen_m"    ,muPtrGen->mass());
             genHAMMER.addUserFloat("mu_gen_charge"  ,muPtrGen->charge());
             genHAMMER.addUserInt(  "mu_gen_pdgid"   ,muPtrGen->pdgId());
   
             genHAMMER.addUserFloat("k1_gen_px"      ,k1PtrGen->px());
             genHAMMER.addUserFloat("k1_gen_py"      ,k1PtrGen->py());
             genHAMMER.addUserFloat("k1_gen_pz"      ,k1PtrGen->pz());
             genHAMMER.addUserFloat("k1_gen_pt"      ,k1PtrGen->pt());
             genHAMMER.addUserFloat("k1_gen_eta"     ,k1PtrGen->eta());
             genHAMMER.addUserFloat("k1_gen_phi"     ,k1PtrGen->phi());
             genHAMMER.addUserFloat("k1_gen_m"    ,k1PtrGen->mass());
             genHAMMER.addUserFloat("k1_gen_charge"  ,k1PtrGen->charge());
             genHAMMER.addUserInt(  "k1_gen_pdgid"   ,k1PtrGen->pdgId());
   
             genHAMMER.addUserFloat("k2_gen_px"      ,k2PtrGen->px());
             genHAMMER.addUserFloat("k2_gen_py"      ,k2PtrGen->py());
             genHAMMER.addUserFloat("k2_gen_pz"      ,k2PtrGen->pz());
             genHAMMER.addUserFloat("k2_gen_pt"      ,k2PtrGen->pt());
             genHAMMER.addUserFloat("k2_gen_eta"     ,k2PtrGen->eta());
             genHAMMER.addUserFloat("k2_gen_phi"     ,k2PtrGen->phi());
             genHAMMER.addUserFloat("k2_gen_m"    ,k2PtrGen->mass());
             genHAMMER.addUserFloat("k2_gen_charge"  ,k2PtrGen->charge());
             genHAMMER.addUserInt(  "k2_gen_pdgid"   ,k2PtrGen->pdgId());
   
             genHAMMER.addUserFloat("pi_gen_px"      ,piPtrGen->px());
             genHAMMER.addUserFloat("pi_gen_py"      ,piPtrGen->py());
             genHAMMER.addUserFloat("pi_gen_pz"      ,piPtrGen->pz());
             genHAMMER.addUserFloat("pi_gen_pt"      ,piPtrGen->pt());
             genHAMMER.addUserFloat("pi_gen_eta"     ,piPtrGen->eta());
             genHAMMER.addUserFloat("pi_gen_phi"     ,piPtrGen->phi());
             genHAMMER.addUserFloat("pi_gen_m"    ,piPtrGen->mass());
             genHAMMER.addUserFloat("pi_gen_charge"  ,piPtrGen->charge());
             genHAMMER.addUserInt(  "pi_gen_pdgid"   ,piPtrGen->pdgId());

             //and gen info from the resonances
             genHAMMER.addUserFloat("phi_gen_px"     ,phiFromK1->px());
             genHAMMER.addUserFloat("phi_gen_py"     ,phiFromK1->py());
             genHAMMER.addUserFloat("phi_gen_pz"     ,phiFromK1->pz());
             genHAMMER.addUserFloat("phi_gen_pt"     ,phiFromK1->pt());
             genHAMMER.addUserFloat("phi_gen_eta"    ,phiFromK1->eta());
             genHAMMER.addUserFloat("phi_gen_phi"    ,phiFromK1->phi());
             genHAMMER.addUserFloat("tv_x_gen"       ,tv_x_gen);//This is the phi production vertex!
             genHAMMER.addUserFloat("tv_y_gen"       ,tv_y_gen);
             genHAMMER.addUserFloat("tv_z_gen"       ,tv_z_gen);
             genHAMMER.addUserFloat("phi_gen_charge" ,phiFromK1->charge());
             genHAMMER.addUserInt(  "phi_gen_pdgid"  ,phiFromK1->pdgId());

             genHAMMER.addUserFloat("ds_gen_px"     ,dsFromPi->px());
             genHAMMER.addUserFloat("ds_gen_py"     ,dsFromPi->py());
             genHAMMER.addUserFloat("ds_gen_pz"     ,dsFromPi->pz());
             genHAMMER.addUserFloat("ds_gen_pt"     ,dsFromPi->pt());
             genHAMMER.addUserFloat("ds_gen_eta"    ,dsFromPi->eta());
             genHAMMER.addUserFloat("ds_gen_phi"    ,dsFromPi->phi());
             genHAMMER.addUserFloat("ds_gen_boost"  ,genDsTlv.BoostVector().Mag());


             genHAMMER.addUserFloat("sv_x_gen"      ,sv_x_gen);//This is the ds production vertex!
             genHAMMER.addUserFloat("sv_y_gen"      ,sv_y_gen);
             genHAMMER.addUserFloat("sv_z_gen"      ,sv_z_gen);
             genHAMMER.addUserFloat("ds_gen_charge" ,dsFromPi->charge());
             genHAMMER.addUserInt(  "ds_gen_pdgid"  ,dsFromPi->pdgId());

             genHAMMER.addUserFloat("bs_gen_px"     ,bsFromMu->px());
             genHAMMER.addUserFloat("bs_gen_py"     ,bsFromMu->py());
             genHAMMER.addUserFloat("bs_gen_pz"     ,bsFromMu->pz());
             genHAMMER.addUserFloat("bs_gen_pt"     ,bsFromMu->pt());
             genHAMMER.addUserFloat("bs_gen_eta"    ,bsFromMu->eta());
             genHAMMER.addUserFloat("bs_gen_phi"    ,bsFromMu->phi());

             genHAMMER.addUserFloat("pv_x_gen"      ,pv_x_gen); //This is the bs production vertex!
             genHAMMER.addUserFloat("pv_y_gen"      ,pv_y_gen);
             genHAMMER.addUserFloat("pv_z_gen"      ,pv_z_gen);

             genHAMMER.addUserFloat("bs_gen_charge" ,bsFromMu->charge());
             genHAMMER.addUserInt(  "bs_gen_pdgid"  ,bsFromMu->pdgId());
             genHAMMER.addUserFloat("b_boost_gen"   ,genBsTlv.BoostVector().Mag());
             genHAMMER.addUserFloat("b_boost_gen_pt"   ,genBsTlv.BoostVector().Pt());
             genHAMMER.addUserFloat("b_boost_gen_eta"   ,genBsTlv.BoostVector().Eta());
             genHAMMER.addUserFloat("b_boost_gen_phi"   ,genBsTlv.BoostVector().Phi());

             //lets also store the fourth vertex ( the k production vertex)
             genHAMMER.addUserFloat("fv_x_gen"       ,fv_x_gen);//This is the kaon production vertex!
             genHAMMER.addUserFloat("fv_y_gen"       ,fv_y_gen);
             genHAMMER.addUserFloat("fv_z_gen"       ,fv_z_gen);

             //define vertex variables

             //float lxyBsGen   = std::sqrt(std::pow((pv_x_gen - sv_x_gen),2) + std::pow((pv_y_gen - sv_y_gen),2) ); 
             //float lxyzBsGen  = std::sqrt(std::pow((pv_x_gen - sv_x_gen),2) + std::pow((pv_y_gen - sv_y_gen),2) + std::pow((pv_z_gen - sv_z_gen),2) ); 
     
             //float lxyDsGen   = std::sqrt(std::pow((sv_x_gen - tv_x_gen),2) + std::pow((sv_y_gen - tv_y_gen),2) ); 
             //float lxyzDsGen  = std::sqrt(std::pow((sv_x_gen - tv_x_gen),2) + std::pow((sv_y_gen - tv_y_gen),2) + std::pow((sv_z_gen - tv_z_gen),2) ); 
     
             //float lxyPhiGen  = std::sqrt(std::pow((tv_x_gen - fv_x_gen),2) + std::pow((tv_y_gen - fv_y_gen),2) ); 
             //float lxyzPhiGen = std::sqrt(std::pow((tv_x_gen - fv_x_gen),2) + std::pow((tv_y_gen - fv_y_gen),2) + std::pow((tv_z_gen - fv_z_gen),2) ); 

             /*
             //gen.addUserFloat("lxy_bs_gen"   ,lxyBsGen);
             //gen.addUserFloat("lxyz_bs_gen"  ,lxyzBsGen);
     
             //gen.addUserFloat("lxy_ds_gen"   ,lxyDsGen);
             //gen.addUserFloat("lxyz_ds_gen"  ,lxyzDsGen);
     
             //gen.addUserFloat("lxy_phi_gen"  ,lxyPhiGen);
             //gen.addUserFloat("lxyz_phi_gen" ,lxyzPhiGen);

             math::XYZPoint pvGen(pv_x_gen,pv_y_gen,pv_z_gen); 
            
              
             float dxyMuGen    = muPtrGen->bestTrack()->dxy(pvGen);  
             //float dxyMuErrGen = muPtrGen->bestTrack()->dxyError(pvGen,pvGen.covariance());  
             //float dxyMuSigGen = dxyMuGen/dxyMuErrGen;
     
             float dzMuGen     = muPtrGen->bestTrack()->dz(pvGen);  
             //float dzMuErrGen  = muPtrGen->bestTrack()->dzError();  
             //float dzMuSigGen  = dzMuGen/dzMuErrGen ; 
     
             float dxyPiGen    = piPtrGen->bestTrack()->dxy(pvGen);  //maybe useful for Ds* vs Ds ? 
             //float dxyPiErrGen = piPtrGen->bestTrack()->dxyError(pvGen,pvGen.covariance());  
             //float dxyPiSigGen = dxyPiGen/dxyPiErrGen;
     
             float dzPiGen     = piPtrGen->bestTrack()->dz(pvGen);  
             //float dzPiErrGen  = piPtrGen->bestTrack()->dzError();  
             //float dzPiSigGen  = dzPiGen/dzPiErrGen ; 
     
             float dxyK1Gen    = k1PtrGen->bestTrack()->dxy(pvGen); //needed ? 
             //float dxyK1ErrGen = k1PtrGen->bestTrack()->dxyError(pvGen,pvGen.covariance());
             //float dxyK1SigGen = dxyK1Gen/dxyK1ErrGen;
     
             float dzK1Gen     = k1PtrGen->bestTrack()->dz(pvGen);  
             //float dzK1ErrGen  = k1PtrGen->bestTrack()->dzError();  
             //float dzK1SigGen  = dzK1Gen/dzK1ErrGen ; 
     
             float dxyK2Gen    = k2PtrGen->bestTrack()->dxy(pvGen); //needed ? 
             //float dxyK2ErrGen = k2PtrGen->bestTrack()->dxyError(pvGen,pvGen.covariance());
             //float dxyK2SigGen = dxyK2Gen/dxyK2ErrGen;
     
             float dzK2Gen     = k2PtrGen->bestTrack()->dz(pvGen);  
             //float dzK2ErrGen  = k2PtrGen->bestTrack()->dzError();  
             //float dzK2SigGen  = dzK2Gen/dzK2ErrGen ; 
     
             //std::cout << "10" << std::endl; 
             gen.addUserFloat("dxy_mu_gen",     dxyMuGen);
             gen.addUserFloat("dz_mu_gen",      dzMuGen);
             //gen.addUserFloat("dxy_mu_err_gen", dxyMuErrGen);
             //gen.addUserFloat("dz_mu_err_gen",  dzMuErrGen);
             //gen.addUserFloat("dxy_mu_sig_gen", dxyMuSigGen);
             //gen.addUserFloat("dz_mu_sig_gen",  dzMuSigGen);
     
             gen.addUserFloat("dxy_pi_gen",     dxyPiGen);
             gen.addUserFloat("dz_pi_gen",      dzPiGen);
             //gen.addUserFloat("dxy_pi_err_gen", dxyPiErrGen);
             //gen.addUserFloat("dz_pi_err_gen",  dzPiErrGen);
             //gen.addUserFloat("dxy_pi_sig_gen", dxyPiSigGen);
             //gen.addUserFloat("dz_pi_sig_gen",  dzPiSigGen);
     
             gen.addUserFloat("dxy_k1_gen",     dxyK1Gen);
             gen.addUserFloat("dz_k1_gen",      dzK1Gen);
             //gen.addUserFloat("dxy_k1_err_gen", dxyK1ErrGen);
             //gen.addUserFloat("dz_k1_err_gen",  dzK1ErrGen);
             //gen.addUserFloat("dxy_k1_sig_gen", dxyK1SigGen);
             //gen.addUserFloat("dz_k1_sig_gen",  dzK1SigGen);
     
             gen.addUserFloat("dxy_k2_gen",     dxyK2Gen);
             gen.addUserFloat("dz_k2_gen",      dzK2Gen);
             //gen.addUserFloat("dxy_k2_err_gen", dxyK2ErrGen);
             //gen.addUserFloat("dz_k2_err_gen",  dzK2ErrGen);
             //gen.addUserFloat("dxy_k2_sig_gen", dxyK2SigGen);
             //gen.addUserFloat("dz_k2_sig_gen",  dzK2SigGen);
             */

             // muon isolation (simply looping and summing track pt's is too simple!)
             /* there is no isolation on gen!
             auto muIso03Gen   = muPtrGen->pfIsolationR03(); //dR = 0.3
             auto muIso04Gen   = muPtrGen->pfIsolationR04(); //dR = 0.4

             float iso03Gen    = muIso03Gen.sumChargedHadronPt() + max(muIso03Gen.sumNeutralHadronEt() + muIso03Gen.sumPhotonEt() - 0.5 * muIso03Gen.sumPUPt(), 0.0)           ,
             float iso04Gen    = muIso04Gen.sumChargedHadronPt() + max(muIso04Gen.sumNeutralHadronEt() + muIso04Gen.sumPhotonEt() - 0.5 * muIso04Gen.sumPUPt(), 0.0)           ,
             float relIso03Gen = iso03Gen / muPtrGen->pt();
             float relIso04Gen = iso04Gen / muPtrGen->pt();

             gen.addUserFloat("mu_iso_03_gen", iso03Gen);
             gen.addUserFloat("mu_iso_04_gen", iso04Gen);
             gen.addUserFloat("mu_rel_iso_03_gen", relIso03Gen);
             gen.addUserFloat("mu_rel_iso_04_gen", relIso04Gen);
             */


             ///////////////////////////////////////////////////////// 
             // now find the channel ID, we have the following scheme:
            
             // SIGNAL:
             // Bs -> Ds mu   0
             // Bs -> Ds tau  1
             //
             // Bs -> Ds* mu  5
             // Bs -> Ds* tau 6
             //
             // HB: (the mother base is defined in step1 for each mother)
             //
             // B mother -> Ds  +  D                              mother_base + 0
             // B mother -> Ds  +  D*                             mother_base + 1
             // B mother -> Ds  +  (something else producing mu)  mother_base + 2

             // B mother -> Ds* +  D                              mother_base + 5
             // B mother -> Ds* +  D*                             mother_base + 6
             // B mother -> Ds* +  (something else producing mu)  mother_base + 7

             // ground state charmed mesons:  
             std::vector<int> dMesons{
               411, // "           : yes 1869
               421, // "           : no  1864
               431, // "           : yes 1968
               4122 // "           : yes 4122
             };

             // excited charmed mesons:  
             std::vector<int> dMesonsExc{
               10411, // charged/mass: yes 2400 
               10421, // "           : no  2400
               413,   // "           : yes 2010
               423,   // "           : no  2007
               10413, // "           : yes 2420
               10423, // "           : no  2420
               20413, // "           : yes ????
               20423, // "           : no  2430
               415,   // "           : yes 2460
               425,   // "           : no  2460
               10431, // "           : yes 2317
               433,   // "           : yes todo
               10433, // "           : yes 2536
               20433, // "           : yes 2460 (resp. 2457)
               435    // "           : yes 2573
             };


             bId = bMotherId;

             // Step1: the b Mother fixes the 10s
             switch(abs(bMotherId)){
               case 521:  sigId = 100;  break;  // B+
               case 511:  sigId = 200;  break;  // B0
               case 531:  sigId = 300;  break;  // Bs
               case 5122: sigId = 400;  break;  // LambdaB
               default:   sigId = 500;          // anything else
             }

             //bool isNotDoubleCharm = false;
              
             //std::cout << "candidate nr: " << nRecoCandidates << std::endl;
             int dsID = getDsID(piPtrGen); // get charmed strange ID
             //std::cout << "ds Id is: " << dsID << std::endl;
             int dID  = getSecondCharmID(muPtrGen); // get charmed ID
             //std::cout << "d Id is: " << dID << std::endl;
             bool isTau = isAncestor(muPtrGen, 15); 

             switch(dsID){
               case 431:   sigId += 0;  break; // Ds+
               case 433:   sigId += 10; break; // Ds+*
               case 10431: sigId += 20; break; // Ds+(2317)*
               case 20433: sigId += 30; break; // Ds+(2457)*
               default:    sigId += 40; break; // anything else
             }

             // Signal candidates enter here
             if ((dID == 0) && ((dsID == 431) ||(dsID == 433))&& (abs(bMotherId) == 531)) {
               //std::cout << "i enter the signal tag!" << std::endl;
               sigId -= 300; // signal should live in 0
               if (isTau) sigId += 1; 
             }

             // special case of B+ -> K nu mu / B+ -> K tau nu 
             else if((dID == 0) && (dsID == 0) && (abs(bMotherId) == 521)){
             
               if (!isTau) sigId += 7;
               if (isTau)  sigId  += 8;

             }

             else {
               switch(dID){
                 case 411:   sigId += 0;  break; // D+
                 case 421:   sigId += 1; break;  // D0
                 case 431:   sigId += 2; break;  // Ds
                 case 413:   sigId += 3; break;  // D+*
                 case 423:   sigId += 4; break;  // D0*
                 case 433:   sigId += 5; break;  // Ds+*
                 case 4122:   sigId += 6; break;  // Lambda c

                 default:    sigId += 9; break; // anything else
               }
             }


             ////////////////////////////////////
             // SPECIAL FOR HAMMER:            //
             // Save also tau and Ds* info     //
             ////////////////////////////////////

             float tau_gen_pt;
             float tau_gen_eta;
             float tau_gen_phi;
             float tau_gen_m;

             float dsStar_gen_pt;
             float dsStar_gen_eta;
             float dsStar_gen_phi;
             float dsStar_gen_m;

             if (sigId == 0){

               //  Bs -> Ds + mu + nu 
               tau_gen_pt = std::nan("nan");
               tau_gen_eta = std::nan("nan");
               tau_gen_phi = std::nan("nan");
               tau_gen_m = std::nan("nan");

               dsStar_gen_pt = std::nan("nan");
               dsStar_gen_eta = std::nan("nan");
               dsStar_gen_phi = std::nan("nan");
               dsStar_gen_m = std::nan("nan");

             }


             else if (sigId == 1){

               //  Bs -> Ds + tau + nu 

               // get the tau (we know its there)
               auto tauFromMu = getAncestor(muReco,15);

               tau_gen_pt     = tauFromMu->pt();
               tau_gen_eta    = tauFromMu->eta();
               tau_gen_phi    = tauFromMu->phi();
               tau_gen_m      = tauMass_; 

               dsStar_gen_pt  = std::nan("nan");
               dsStar_gen_eta = std::nan("nan");
               dsStar_gen_phi = std::nan("nan");
               dsStar_gen_m   = std::nan("nan");

             }

             else if (sigId == 10){

               //  Bs -> Ds* + mu + nu 

               tau_gen_pt  = std::nan("nan");
               tau_gen_eta = std::nan("nan");
               tau_gen_phi = std::nan("nan");
               tau_gen_m   = std::nan("nan");

               // get the Ds* (we know its there)
               auto dsStarFromDs = getAncestor(dsFromPi,433);

               dsStar_gen_pt     = dsStarFromDs->pt();
               dsStar_gen_eta    = dsStarFromDs->eta();
               dsStar_gen_phi    = dsStarFromDs->phi();
               dsStar_gen_m      = dsStarMass_; 

             }

             else if (sigId == 11){

               //  Bs -> Ds* + tau + nu 

               // get the tau (we know its there)
               auto tauFromMu = getAncestor(muReco,15);

               tau_gen_pt     = tauFromMu->pt();
               tau_gen_eta    = tauFromMu->eta();
               tau_gen_phi    = tauFromMu->phi();
               tau_gen_m      = tauMass_; 

               // get the Ds* (we know its there)
               auto dsStarFromDs = getAncestor(dsFromPi,433);

               dsStar_gen_pt     = dsStarFromDs->pt();
               dsStar_gen_eta    = dsStarFromDs->eta();
               dsStar_gen_phi    = dsStarFromDs->phi();
               dsStar_gen_m      = dsStarMass_; 

             }

             else{

               tau_gen_pt  = std::nan("nan");
               tau_gen_eta = std::nan("nan");
               tau_gen_phi = std::nan("nan");
               tau_gen_m   = std::nan("nan");

               dsStar_gen_pt  = std::nan("nan");
               dsStar_gen_eta = std::nan("nan");
               dsStar_gen_phi = std::nan("nan");
               dsStar_gen_m   = std::nan("nan");

             }


             genHAMMER.addUserFloat("tau_gen_pt",  tau_gen_pt);
             genHAMMER.addUserFloat("tau_gen_eta", tau_gen_eta);
             genHAMMER.addUserFloat("tau_gen_phi", tau_gen_phi);
             genHAMMER.addUserFloat("tau_gen_m",   tau_gen_m);

             genHAMMER.addUserFloat("dsStar_gen_pt",  dsStar_gen_pt);
             genHAMMER.addUserFloat("dsStar_gen_eta", dsStar_gen_eta);
             genHAMMER.addUserFloat("dsStar_gen_phi", dsStar_gen_phi);
             genHAMMER.addUserFloat("dsStar_gen_m",   dsStar_gen_m);

             /*
             bool isDsStar    = false;
             bool isDMeson    = false;
             bool isDMesonExc = false;
             bool isSignal    = true;
             bool isTauSignal = false;

             // Step2: distinguish between Ds + charm / Ds* + charm
             if(isAncestor(piPtrGen, 433)) isDsStar = true; 

             // Step3: distinguish if the mu is coming from a second charmed meson
             int dMesonIdx = -1;

             for(size_t dIdx = 0; dIdx < dMesons.size(); dIdx++){
               if(isAncestor(muPtrGen, dMesons.at(dIdx))){
                 dMesonIdx = dIdx; 
                 isDMeson = true;
                 isSignal = false;
                 break;
               }
             } 

             if (dMesonIdx == 0) {
               //only enter here when upper loop didnt find d meson
               for(size_t dIdx = 0; dIdx < dMesonsExc.size(); dIdx++){
                 if(isAncestor(muPtrGen, dMesonsExc.at(dIdx))){
                   dMesonIdx = dIdx;
                   isDMesonExc = true;
                   isSignal = false;
                   break;
                 }
               }  
             }
                            
             // step4: check if its a true signal or B-meson -> Ds + something leptonaically decaying 

             //if its not bs and no double charm, its something else leptonically decaying
             if ((dMesonIdx == -1) && (abs(bMotherId) != 531)) isSignal = false; 
            
             if ((dMesonIdx == -1) && (abs(bMotherId) == 531)) {
               // if its Bs and no double charm event, really be sure that we have the signal!

               std::vector<int> daus;
               for( size_t dauIdx = 0; dauIdx < bsFromMu->numberOfDaughters(); dauIdx++){

                 int dauId = bsFromMu->daughter(dauIdx)->pdgId();
                 //                  Ds*                   Ds                nu_tau                   tau                   nu_mu                   mu                photon
                 if ( (abs(dauId) != 433) && (abs(dauId) != 431) &&  (abs(dauId) != 16 ) &&  (abs(dauId) != 15 ) &&  (abs(dauId) != 14 ) &&  (abs(dauId) != 13 ) && (abs(dauId) != 22)){
                   isSignal = false; // daughter is none of the involved signal particles!
                 }
                 daus.push_back(abs(dauId)); //append abs of daughter id
                 if ( abs(dauId) == 15) isTauSignal = true; // found a tau!
               }
             }

             // now we assign the signal Id based on our boolean flags
             //printDaughters(bsFromMu); //-> for debugging

             if (isSignal)                                sigId  = 0; // signal
             if (isDsStar)                                sigId += 5;
             if (isDMesonExc)                             sigId += 1; 
             if (!isSignal && !isDMeson && !isDMesonExc)  sigId += 2; // if this is true, isDMesonEx should be false!
             if (isSignal && isTauSignal)                 sigId += 1; // only for signals

             */

             printDaughters(bsFromMu); //-> for debugging
             std::cout << "mom:" << bMotherId << std::endl; //for debugging
             std::cout << "ID is:" << sigId     << std::endl;  //for debugging;

             /*
             //define helicity angles

             //test lhcb method on gen
             TLorentzVector lhcbBsTlvGen = lhcbMethod(genDsTlv + genMuTlv, pv_x_gen, pv_y_gen, pv_z_gen, sv_x_gen, sv_y_gen, sv_z_gen, bsMass_);      
             //std::cout << "this is the gen matched case: " << std::endl;
             std::tuple<std::vector<TLorentzVector>,float> recoResultGen = recoMethod(genDsTlv + genMuTlv, pv_x_gen, pv_y_gen, pv_z_gen, sv_x_gen, sv_y_gen, sv_z_gen, bsMass_);      

             std::vector<TLorentzVector> recosGen = std::get<0>(recoResultGen);
             float discNegativityGen                 = std::get<1>(recoResultGen);

             int discIsNegativeGen = 0;
             if (discNegativityGen > 0) discIsNegativeGen = 1;

             TLorentzVector reco1BsTlvGen = recosGen.at(0);
             TLorentzVector reco2BsTlvGen = recosGen.at(1);

             genHAMMER.addUserInt("disc_is_negative_gen", discIsNegativeGen); 
             genHAMMER.addUserFloat("disc_negativity_gen", discNegativityGen); 

             genHAMMER.addUserFloat("bs_gen_lhcb_pt",   lhcbBsTlvGen.Pt());
             genHAMMER.addUserFloat("bs_gen_lhcb_eta",  lhcbBsTlvGen.Eta());
             genHAMMER.addUserFloat("bs_gen_lhcb_phi",  lhcbBsTlvGen.Phi());

             //angle between Mu and W


             float angMuWGen = angMuW(genDsTlv,genBsTlv,genMuTlv); 

             float angMuWGenLhcb  = angMuW(genDsTlv,lhcbBsTlvGen,genMuTlv); 
             float angMuWGenReco1 = angMuW(genDsTlv,reco1BsTlvGen,genMuTlv); 
             float angMuWGenReco2 = angMuW(genDsTlv,reco2BsTlvGen,genMuTlv); 

             genHAMMER.addUserFloat("angMuWGen",angMuWGen);
             genHAMMER.addUserFloat("cosMuWGen",cos(angMuWGen));
             genHAMMER.addUserFloat("cosMuWGenLhcb",cos(angMuWGenLhcb));
             genHAMMER.addUserFloat("cosMuWGenReco1",cos(angMuWGenReco1));
             genHAMMER.addUserFloat("cosMuWGenReco2",cos(angMuWGenReco2));


             //angle between k1(k2) and pion in phi rest frame
             float angPiK1Gen  = angDoubleDecay(genPhiTlv, genK1Tlv,  genPiTlv);
             float angPiK2Gen  = angDoubleDecay(genPhiTlv, genK2Tlv,  genPiTlv);
             genHAMMER.addUserFloat("angPiK1Gen", angPiK1Gen);
             genHAMMER.addUserFloat("angPiK2Gen", angPiK2Gen);
             genHAMMER.addUserFloat("cosPiK1Gen", cos(angPiK1Gen));
             genHAMMER.addUserFloat("cosPiK2Gen", cos(angPiK2Gen));

             // equivalently, angle between phi(pi) and bs in ds rest frame
             float angPhiDsGen = angDsPi(genDsTlv,  genPhiTlv, genBsTlv);
             float angPiDsGen  = angDsPi(genDsTlv,  genPiTlv,  genBsTlv);
             float angPiDsGenLhcb  = angDsPi(genDsTlv,  genPiTlv,  lhcbBsTlvGen);

             genHAMMER.addUserFloat("angPhiDsGen", angPhiDsGen);
             genHAMMER.addUserFloat("angPiDsGen",  angPiDsGen);
             genHAMMER.addUserFloat("cosPhiDsGen", cos(angPhiDsGen));
             genHAMMER.addUserFloat("cosPiDsGen",  cos(angPiDsGen));
             genHAMMER.addUserFloat("cosPiDsGenLhcb",  cos(angPiDsGenLhcb));

             //plane angle
             float angPlaneBsGen = angPlane(genDsTlv, genBsTlv, genMuTlv, genPiTlv);
             genHAMMER.addUserFloat("angPlaneBsGen", angPlaneBsGen);
             genHAMMER.addUserFloat("cosPlaneBsGen", cos(angPlaneBsGen));

             float angPlaneDsGen = angPlane2(genDsTlv, genBsTlv, genK1Tlv, genPiTlv);
             genHAMMER.addUserFloat("angPlaneDsGen", angPlaneDsGen);
             genHAMMER.addUserFloat("cosPlaneDsGen", cos(angPlaneDsGen));
             */
             //if we reached this point we have found our gen match and we can stop the loop
             // since now we are not trying to find only 1 gen match for a reco-candidate
             // (we dont loop over the bs candidates) but rather over the events, we have to allow
             // several caniddates per event!! As we also allow several bs candidates per event!

             genHAMMER.addUserInt("sig",sigId);
             genHAMMER.addUserInt("b_mother_id",bId);
             genHAMMER.addUserInt("gen_match_success",genMatchSuccess);

             ret_value->emplace_back(genHAMMER);
             //std::cout << "i have size nr:" << ret_value->size() << std::endl;
             //std::cout << "pion pt:" << piPtrGen->pt() << std::endl;
             //std::cout << "k1 pt:" << k1PtrGen->pt() << std::endl;
             //std::cout << "k2 pt:" << k2PtrGen->pt() << std::endl;
             //std::cout << "mu pt:" << muPtrGen->pt() << std::endl;
             //////////////////////////////////////////////////

           }//close gen matching pi loop 
          //break;
          }//close gen matching k2 loop 
        //break;
        } //close gen matching k1 loop
      //std::cout << "found kaons:" << kaons << std::endl;
      //break;
      } //close gen matching mu loop

      //if (genMatchSuccess == 0) continue; 

     /////////////////////// END OF VARIABLE DEFINITION //////////////////////

iEvent.put(std::move(ret_value), "genHAMMER");
}//closing event loop


void inspectorHAMMER::endJob(){
// Printouts:
std::cout << "\n--------- GEN MATCHING MODULE ----------\n" << std::endl;
std::cout << "#Events in file                                           : " << nEventsHAMMER  << std::endl;
std::cout << "#Gen Muons in file                                        : " << nMuonsHAMMER   << std::endl;
std::cout << "#Gen Tracks in file                                       : " << nTracksHAMMER  << std::endl;
std::cout << "#Gen matched candidates                                   : " << nGenMatchedHAMMER << std::endl;
std::cout << " Before phi mass cut we have                              : " << beforePhiMassHAMMER << std::endl;
std::cout << " After phi mass cut we have                              : " << afterPhiMassHAMMER << std::endl;

std::cout << "\n#KKPiMu Gen combinations:" << nKKPiMuHAMMER << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a Phi         : " << nFoundPhiHAMMER  << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a Ds          : " << nFoundDsHAMMER   << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a B-mom       : " << nFoundBHAMMER    << std::endl;
std::cout << "#KKPiMu Gen combinations for which the B-mom < B mass cut : " << nBMassCutHAMMER << std::endl;
}


DEFINE_FWK_MODULE(inspectorHAMMER);
