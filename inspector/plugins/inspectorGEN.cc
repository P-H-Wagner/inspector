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
int nEventsGEN = 0;
int nMuonsGEN  = 0;
int nTracksGEN = 0;
int nGenMatchedGEN = 0;

int beforePhiMass = 0;
int afterPhiMass = 0;

int nKKPiMuGEN = 0;
int nFoundPhiGEN  = 0;
int nFoundDsGEN   = 0;
int nFoundBGEN    = 0;
int nBMassCutGEN   = 0;

class inspectorGEN : public edm::global::EDProducer<> {

public:

  //define collections which dont exist by default  
  typedef std::vector<reco::GenParticle> GenParticleCollection;
  typedef std::vector<pat::PackedGenParticle> PackedGenParticleCollection;
  //constructor
  explicit inspectorGEN(const edm::ParameterSet&);
  //destructor
  ~inspectorGEN() override {}
  
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
  //on pure gen level, these collection is called genParticle
  const edm::InputTag genTag; 
  const edm::EDGetTokenT<reco::GenParticleCollection> pureGen_;

};

//define the constructor
inspectorGEN::inspectorGEN(const edm::ParameterSet& iConfig):
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
       produces<pat::CompositeCandidateCollection>("gen");
       //produces<pat::CompositeCandidateCollection>("gen");
       //produces<TransientTrackCollection>("kkTransientTracks");
    }

//check const keywords 

// this starts the event loop
void inspectorGEN::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {

  //input
  edm::Handle<reco::GenParticleCollection> pureGen;
  iEvent.getByToken(pureGen_,pureGen);

  // to save 
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());

  nEventsGEN++;

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

  //std::cout << "New event!" << iEvent.id().event() <<std::endl;

  ////////////////////////////////////////////////////
  // find the gen-matched muon                      //
  ////////////////////////////////////////////////////

  for(size_t muIdxGen = 0; muIdxGen < pureGen->size(); ++muIdxGen){

    nMuonsGEN++;

    //define a pointer to the gen muon    
    edm::Ptr<reco::GenParticle> muPtrGen(pureGen, muIdxGen);

    // We dont want any filter, we only tag gen
    if((fabs(muPtrGen->pdgId()) != 13)) continue; 

    ////////////////////////////////////////////////
    // find gen matched k1                        //
    ////////////////////////////////////////////////
    int kaons = 0;
    for(size_t k1IdxGen = 0; k1IdxGen < pureGen->size(); ++k1IdxGen){
         
      nTracksGEN++;
      //define a pointer to the gen kaon    
      edm::Ptr<reco::GenParticle> k1PtrGen(pureGen, k1IdxGen);

      // We dont want any filter, we only tag gen
      if((fabs(k1PtrGen->pdgId()) != 321)) continue; 
      kaons++;

      ////////////////////////////////////////////////
      // find gen matched k2                        //
      ////////////////////////////////////////////////
   
      for(size_t k2IdxGen = k1IdxGen + 1; k2IdxGen < pureGen->size(); ++k2IdxGen){
   
           //avoid picking the same gen particle as for k1
           if(k2IdxGen == k1IdxGen) continue; 

           //define a pointer to the gen kaon    
           edm::Ptr<reco::GenParticle> k2PtrGen(pureGen, k2IdxGen);
   
           // We dont want any filter, we only tag gen
           if((fabs(k2PtrGen->pdgId()) != 321)) continue;   

           beforePhiMass++;
           ////////////////////////////////////////////////
           // find gen matched pion                      //
           ////////////////////////////////////////////////

           for(size_t piIdxGen = 0; piIdxGen < pureGen->size(); ++piIdxGen){
    
             //avoid picking the same gen particle as for k1 or k2
             if((piIdxGen == k1IdxGen) || (piIdxGen == k2IdxGen)) continue; 
                     
             //define a pointer to the gen kaon    
             edm::Ptr<reco::GenParticle> piPtrGen(pureGen, piIdxGen);
   
             // We dont want any filter, we only tag gen
             if((fabs(piPtrGen->pdgId()) != 211)) continue;  

             //////////////////////////////////////////////////
             // Find resonances at gen level                 //
             //////////////////////////////////////////////////
             pat::CompositeCandidate gen; 

             nKKPiMuGEN++;     

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
             nFoundPhiGEN++;               
   
             //std::cout<< "found phi candiadte!" << std::endl;

             // searching for ds resonance 
             auto dsFromPhi = getAncestor(phiFromK1,431);
             auto dsFromPi  = getAncestor(piReco,431);
             if( (dsFromPhi != dsFromPi) || (dsFromPhi == nullptr) || (dsFromPi == nullptr)) continue; 
             nFoundDsGEN++;               

             //std::cout<< "found ds candiadte!" << std::endl;
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
              
             //std::cout<< "mom idx is: " << bMotherId << std::endl;
             if (bMotherId == 0) continue; // no b mother found

             //std::cout<< "WE HAVE A B MOM" << bMotherId << std::endl;

             // Even if the mu is not required to come brom the b mother directly (would be signal case)
             // if it comes from another D meson (double charm background case), we still want
             // that this D meson is coming from the b mother. So the muon should share
             // the same ancestor as the Ds.

             auto bsFromDs = getAncestor(dsFromPhi,bMotherId);

             //std::cout << "bs from ds has chain: " << std::endl;
             //printDaughters(bsFromDs);

             auto bsFromMu = getAncestor(muReco,   bMotherId);

             //std::cout << "bs from mu has chain: " << std::endl;
             //printDaughters(bsFromMu);

             if( (bsFromDs != bsFromMu) || (bsFromDs == nullptr) || (bsFromMu == nullptr)) {
 
             //std::cout << "b pointer is null" << std::endl; 
             continue;
             } 
 
             nFoundBGEN++;
             
             //if (bsFromMu->mass() > maxBsMass_) continue;
             nBMassCutGEN++;

             //remove oscillations, this Bs is only needed for the vertex calculation!
             //the four momenta which interests us is the one from the oscillated one, the one that produces the daughters
             auto bsFromMuWOOsc = removeOscillations(bsFromMu);

             nGenMatches++;
             genMatchSuccess = 1;
             nGenMatchedGEN++;
             
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

             gen.addUserFloat("m2_miss_gen",m2_miss_gen);
             gen.addUserFloat("pt_miss_gen",pt_miss_gen);
             gen.addUserFloat("e_star_gen",e_star_gen);
             gen.addUserFloat("e_gamma_gen",e_gamma_gen);
             gen.addUserFloat("q2_gen",q2_gen);

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
             gen.addUserFloat("mu_gen_px"      ,muPtrGen->px());
             gen.addUserFloat("mu_gen_py"      ,muPtrGen->py());
             gen.addUserFloat("mu_gen_pz"      ,muPtrGen->pz());
             gen.addUserFloat("mu_gen_pt"      ,muPtrGen->pt());
             gen.addUserFloat("mu_gen_eta"     ,muPtrGen->eta());
             gen.addUserFloat("mu_gen_phi"     ,muPtrGen->phi());
             gen.addUserFloat("mu_gen_m"    ,muPtrGen->mass());
             gen.addUserFloat("mu_gen_charge"  ,muPtrGen->charge());
             gen.addUserInt(  "mu_gen_pdgid"   ,muPtrGen->pdgId());
   
             gen.addUserFloat("k1_gen_px"      ,k1PtrGen->px());
             gen.addUserFloat("k1_gen_py"      ,k1PtrGen->py());
             gen.addUserFloat("k1_gen_pz"      ,k1PtrGen->pz());
             gen.addUserFloat("k1_gen_pt"      ,k1PtrGen->pt());
             gen.addUserFloat("k1_gen_eta"     ,k1PtrGen->eta());
             gen.addUserFloat("k1_gen_phi"     ,k1PtrGen->phi());
             gen.addUserFloat("k1_gen_m"    ,k1PtrGen->mass());
             gen.addUserFloat("k1_gen_charge"  ,k1PtrGen->charge());
             gen.addUserInt(  "k1_gen_pdgid"   ,k1PtrGen->pdgId());
   
             gen.addUserFloat("k2_gen_px"      ,k2PtrGen->px());
             gen.addUserFloat("k2_gen_py"      ,k2PtrGen->py());
             gen.addUserFloat("k2_gen_pz"      ,k2PtrGen->pz());
             gen.addUserFloat("k2_gen_pt"      ,k2PtrGen->pt());
             gen.addUserFloat("k2_gen_eta"     ,k2PtrGen->eta());
             gen.addUserFloat("k2_gen_phi"     ,k2PtrGen->phi());
             gen.addUserFloat("k2_gen_m"    ,k2PtrGen->mass());
             gen.addUserFloat("k2_gen_charge"  ,k2PtrGen->charge());
             gen.addUserInt(  "k2_gen_pdgid"   ,k2PtrGen->pdgId());
   
             gen.addUserFloat("pi_gen_px"      ,piPtrGen->px());
             gen.addUserFloat("pi_gen_py"      ,piPtrGen->py());
             gen.addUserFloat("pi_gen_pz"      ,piPtrGen->pz());
             gen.addUserFloat("pi_gen_pt"      ,piPtrGen->pt());
             gen.addUserFloat("pi_gen_eta"     ,piPtrGen->eta());
             gen.addUserFloat("pi_gen_phi"     ,piPtrGen->phi());
             gen.addUserFloat("pi_gen_m"    ,piPtrGen->mass());
             gen.addUserFloat("pi_gen_charge"  ,piPtrGen->charge());
             gen.addUserInt(  "pi_gen_pdgid"   ,piPtrGen->pdgId());

             //and gen info from the resonances
             gen.addUserFloat("phi_gen_px"     ,phiFromK1->px());
             gen.addUserFloat("phi_gen_py"     ,phiFromK1->py());
             gen.addUserFloat("phi_gen_pz"     ,phiFromK1->pz());
             gen.addUserFloat("phi_gen_pt"     ,phiFromK1->pt());
             gen.addUserFloat("phi_gen_eta"    ,phiFromK1->eta());
             gen.addUserFloat("phi_gen_phi"    ,phiFromK1->phi());
             gen.addUserFloat("tv_x_gen"       ,tv_x_gen);//This is the phi production vertex!
             gen.addUserFloat("tv_y_gen"       ,tv_y_gen);
             gen.addUserFloat("tv_z_gen"       ,tv_z_gen);
             gen.addUserFloat("phi_gen_charge" ,phiFromK1->charge());
             gen.addUserInt(  "phi_gen_pdgid"  ,phiFromK1->pdgId());

             gen.addUserFloat("ds_gen_px"     ,dsFromPi->px());
             gen.addUserFloat("ds_gen_py"     ,dsFromPi->py());
             gen.addUserFloat("ds_gen_pz"     ,dsFromPi->pz());
             gen.addUserFloat("ds_gen_pt"     ,dsFromPi->pt());
             gen.addUserFloat("ds_gen_eta"    ,dsFromPi->eta());
             gen.addUserFloat("ds_gen_phi"    ,dsFromPi->phi());
             gen.addUserFloat("ds_gen_boost"  ,genDsTlv.BoostVector().Mag());


             gen.addUserFloat("sv_x_gen"      ,sv_x_gen);//This is the ds production vertex!
             gen.addUserFloat("sv_y_gen"      ,sv_y_gen);
             gen.addUserFloat("sv_z_gen"      ,sv_z_gen);
             gen.addUserFloat("ds_gen_m"      ,dsFromPi->mass());
             gen.addUserFloat("ds_gen_charge" ,dsFromPi->charge());
             gen.addUserInt(  "ds_gen_pdgid"  ,dsFromPi->pdgId());

             gen.addUserFloat("bs_gen_px"     ,bsFromMu->px());
             gen.addUserFloat("bs_gen_py"     ,bsFromMu->py());
             gen.addUserFloat("bs_gen_pz"     ,bsFromMu->pz());
             gen.addUserFloat("bs_gen_pt"     ,bsFromMu->pt());
             gen.addUserFloat("bs_gen_eta"    ,bsFromMu->eta());
             gen.addUserFloat("bs_gen_phi"    ,bsFromMu->phi());

             gen.addUserFloat("pv_x_gen"      ,pv_x_gen); //This is the bs production vertex!
             gen.addUserFloat("pv_y_gen"      ,pv_y_gen);
             gen.addUserFloat("pv_z_gen"      ,pv_z_gen);

             gen.addUserFloat("bs_gen_m"      ,bsFromMu->mass());
             gen.addUserFloat("bs_gen_charge" ,bsFromMu->charge());
             gen.addUserInt(  "bs_gen_pdgid"  ,bsFromMu->pdgId());
             gen.addUserFloat("b_boost_gen"   ,genBsTlv.BoostVector().Mag());
             gen.addUserFloat("b_boost_gen_pt"   ,genBsTlv.BoostVector().Pt());
             gen.addUserFloat("b_boost_gen_eta"   ,genBsTlv.BoostVector().Eta());
             gen.addUserFloat("b_boost_gen_phi"   ,genBsTlv.BoostVector().Phi());

             //lets also store the fourth vertex ( the k production vertex)
             gen.addUserFloat("fv_x_gen"       ,fv_x_gen);//This is the kaon production vertex!
             gen.addUserFloat("fv_y_gen"       ,fv_y_gen);
             gen.addUserFloat("fv_z_gen"       ,fv_z_gen);


             ///////////////////////////// 
             // now find the channel ID //
             /////////////////////////////

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
             auto photonPtr = printDirectDaughters(bsFromMu, false);             
 
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
               default:    sigId = 500; break; // anything else
             }

             int checkSignal = -1;
             int checkKNuMu  = -1;

             if (abs(bMotherId) == 531) checkSignal = isSignal(bsFromMu); 
             if (abs(bMotherId) == 521) checkKNuMu  = isKNuMu(bsFromMu); 

             if (checkSignal==-1){
               std::cout << "this is not tagged as signal" << std::endl;
               //printDaughters(bsFromMu);
             }

             // Signal candidates enter here
             if (checkSignal != -1)    sigId = checkSignal;
             else if(checkKNuMu != -1) sigId = checkKNuMu;

             else {
               switch(dID){
                 case 411:   sigId += 0;  break; // D+
                 case 421:   sigId += 1; break;  // D0
                 case 431:   sigId += 2; break;  // Ds
                 case 413:   sigId += 3; break;  // D+*
                 case 423:   sigId += 4; break;  // D0*
                 case 433:   sigId += 5; break;  // Ds+*
                 case 4122:  sigId += 6; break;  // Lambda c

                 default:    sigId = 500; break; // anything else
               }
             }

             // we want to be sure that we dont tag a Hb channel which was not simulated!
             // since the other b meson can decay freely, this can happen! F.e. we can have stuff like
             // Bs -> Double charm + additional pions/kaons/gammas
             std::set<int> bs_channels      = {302, 300, 303, 312, 305, 315, 0, 1, 10 ,11}; 
             std::set<int> b0_channels      = {200, 203, 210, 213, 212, 205, 215, 220, 223, 230, 233}; 
             std::set<int> bplus_channels   = { 107, 108, 101, 104, 117, 118, 121, 124, 131, 134, 111, 114}; 
             std::set<int> lambdab_channels = {406, 416}; 

             if      ((abs(bMotherId) == 531) && (bs_channels.find(sigId)       == bs_channels.end() ))      sigId = 500;
             else if ((abs(bMotherId) == 521) && (bplus_channels.find(sigId)    == bplus_channels.end() ))   sigId = 500;
             else if ((abs(bMotherId) == 511) && (b0_channels.find(sigId)       == b0_channels.end() ))      sigId = 500;
             else if ((abs(bMotherId) == 5122) && (lambdab_channels.find(sigId) == lambdab_channels.end() )) sigId = 500;

             double photonEnergy;
             if (photonPtr != nullptr) photonEnergy =photonPtr->energy();
             else photonEnergy = std::nan("nan");

             gen.addUserFloat("photon_energy", photonEnergy);

             /////////////////////////////////////////////



             ////////////////////////////////////
             // SPECIAL FOR HAMMER:            //
             // Save also tau and Ds* info     //
             ////////////////////////////////////

             float tau_gen_pt;
             float tau_gen_eta;
             float tau_gen_phi;
             float tau_gen_m;
             int   tau_gen_pdgid; 

             float dsStar_gen_pt;
             float dsStar_gen_eta;
             float dsStar_gen_phi;
             float dsStar_gen_m;
             int   dsStar_gen_pdgid; 

             if (sigId == 0){

               //  Bs -> Ds + mu + nu 
               tau_gen_pt       = std::nan("nan");
               tau_gen_eta      = std::nan("nan");
               tau_gen_phi      = std::nan("nan");
               tau_gen_m        = std::nan("nan");
               tau_gen_pdgid    = -9999;

               dsStar_gen_pt    = std::nan("nan");
               dsStar_gen_eta   = std::nan("nan");
               dsStar_gen_phi   = std::nan("nan");
               dsStar_gen_m     = std::nan("nan");
               dsStar_gen_pdgid = -9999;

             }


             else if (sigId == 1){

               //  Bs -> Ds + tau + nu 

               // get the tau (we know its there)
               auto tauFromMu   = getAncestor(muReco,15);

               tau_gen_pt       = tauFromMu->pt();
               tau_gen_eta      = tauFromMu->eta();
               tau_gen_phi      = tauFromMu->phi();
               tau_gen_m        = tauMass_; 
               tau_gen_pdgid    = tauFromMu->pdgId();

               dsStar_gen_pt    = std::nan("nan");
               dsStar_gen_eta   = std::nan("nan");
               dsStar_gen_phi   = std::nan("nan");
               dsStar_gen_m     = std::nan("nan");
               dsStar_gen_pdgid = -9999;

             }

             else if (sigId == 10){

               //  Bs -> Ds* + mu + nu 

               tau_gen_pt       = std::nan("nan");
               tau_gen_eta      = std::nan("nan");
               tau_gen_phi      = std::nan("nan");
               tau_gen_m        = std::nan("nan");
               tau_gen_pdgid    = -9999;

               // get the Ds* (we know its there)
               auto dsStarFromDs = getAncestor(dsFromPi,433);

               dsStar_gen_pt    = dsStarFromDs->pt();
               dsStar_gen_eta   = dsStarFromDs->eta();
               dsStar_gen_phi   = dsStarFromDs->phi();
               dsStar_gen_m     = dsStarMass_; 
               dsStar_gen_pdgid = dsStarFromDs->pdgId();

             }

             else if (sigId == 11){

               //  Bs -> Ds* + tau + nu 

               // get the tau (we know its there)
               auto tauFromMu = getAncestor(muReco,15);

               tau_gen_pt       = tauFromMu->pt();
               tau_gen_eta      = tauFromMu->eta();
               tau_gen_phi      = tauFromMu->phi();
               tau_gen_m        = tauMass_; 
               tau_gen_pdgid    = tauFromMu->pdgId();

               // get the Ds* (we know its there)
               auto dsStarFromDs = getAncestor(dsFromPi,433);

               dsStar_gen_pt     = dsStarFromDs->pt();
               dsStar_gen_eta    = dsStarFromDs->eta();
               dsStar_gen_phi    = dsStarFromDs->phi();
               dsStar_gen_m      = dsStarMass_; 
               dsStar_gen_pdgid  = dsStarFromDs->pdgId();

             }

             else{

               tau_gen_pt      = std::nan("nan");
               tau_gen_eta     = std::nan("nan");
               tau_gen_phi     = std::nan("nan");
               tau_gen_m       = std::nan("nan");
               tau_gen_pdgid   = -9999;

               dsStar_gen_pt     = std::nan("nan");
               dsStar_gen_eta    = std::nan("nan");
               dsStar_gen_phi    = std::nan("nan");
               dsStar_gen_m      = std::nan("nan");
               dsStar_gen_pdgid  = -9999;

             }


             gen.addUserFloat("tau_gen_pt",         tau_gen_pt);
             gen.addUserFloat("tau_gen_eta",        tau_gen_eta);
             gen.addUserFloat("tau_gen_phi",        tau_gen_phi);
             gen.addUserFloat("tau_gen_m",          tau_gen_m);
             gen.addUserInt("tau_gen_pdgid",      tau_gen_pdgid);

             gen.addUserFloat("dsStar_gen_pt",      dsStar_gen_pt);
             gen.addUserFloat("dsStar_gen_eta",     dsStar_gen_eta);
             gen.addUserFloat("dsStar_gen_phi",     dsStar_gen_phi);
             gen.addUserFloat("dsStar_gen_m",       dsStar_gen_m);
             gen.addUserInt("dsStar_gen_pdgid",   dsStar_gen_pdgid);



             //auto dummy = printDirectDaughters(bsFromMu, false); //-> for debugging
             //printDaughters(bsFromMu); //-> for debugging
             //std::cout << "mom:" << bMotherId << std::endl; //for debugging
             //std::cout << "ID is:" << sigId     << std::endl;  //for debugging;

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

             gen.addUserInt("disc_is_negative_gen", discIsNegativeGen); 
             gen.addUserFloat("disc_negativity_gen", discNegativityGen); 

             gen.addUserFloat("bs_gen_lhcb_pt",   lhcbBsTlvGen.Pt());
             gen.addUserFloat("bs_gen_lhcb_eta",  lhcbBsTlvGen.Eta());
             gen.addUserFloat("bs_gen_lhcb_phi",  lhcbBsTlvGen.Phi());

             //angle between Mu and W


             float angMuWGen = angMuW(genDsTlv,genBsTlv,genMuTlv); 

             float angMuWGenLhcb  = angMuW(genDsTlv,lhcbBsTlvGen,genMuTlv); 
             float angMuWGenReco1 = angMuW(genDsTlv,reco1BsTlvGen,genMuTlv); 
             float angMuWGenReco2 = angMuW(genDsTlv,reco2BsTlvGen,genMuTlv); 

             gen.addUserFloat("angMuWGen",angMuWGen);
             gen.addUserFloat("cosMuWGen",cos(angMuWGen));
             gen.addUserFloat("cosMuWGenLhcb",cos(angMuWGenLhcb));
             gen.addUserFloat("cosMuWGenReco1",cos(angMuWGenReco1));
             gen.addUserFloat("cosMuWGenReco2",cos(angMuWGenReco2));


             //angle between k1(k2) and pion in phi rest frame
             float angPiK1Gen  = angDoubleDecay(genPhiTlv, genK1Tlv,  genPiTlv);
             float angPiK2Gen  = angDoubleDecay(genPhiTlv, genK2Tlv,  genPiTlv);
             gen.addUserFloat("angPiK1Gen", angPiK1Gen);
             gen.addUserFloat("angPiK2Gen", angPiK2Gen);
             gen.addUserFloat("cosPiK1Gen", cos(angPiK1Gen));
             gen.addUserFloat("cosPiK2Gen", cos(angPiK2Gen));

             // equivalently, angle between phi(pi) and bs in ds rest frame
             float angPhiDsGen = angDsPi(genDsTlv,  genPhiTlv, genBsTlv);
             float angPiDsGen  = angDsPi(genDsTlv,  genPiTlv,  genBsTlv);
             float angPiDsGenLhcb  = angDsPi(genDsTlv,  genPiTlv,  lhcbBsTlvGen);

             gen.addUserFloat("angPhiDsGen", angPhiDsGen);
             gen.addUserFloat("angPiDsGen",  angPiDsGen);
             gen.addUserFloat("cosPhiDsGen", cos(angPhiDsGen));
             gen.addUserFloat("cosPiDsGen",  cos(angPiDsGen));
             gen.addUserFloat("cosPiDsGenLhcb",  cos(angPiDsGenLhcb));

             //plane angle
             float angPlaneBsGen = angPlane(genDsTlv, genBsTlv, genMuTlv, genPiTlv);
             gen.addUserFloat("angPlaneBsGen", angPlaneBsGen);
             gen.addUserFloat("cosPlaneBsGen", cos(angPlaneBsGen));

             float angPlaneDsGen = angPlane2(genDsTlv, genBsTlv, genK1Tlv, genPiTlv);
             gen.addUserFloat("angPlaneDsGen", angPlaneDsGen);
             gen.addUserFloat("cosPlaneDsGen", cos(angPlaneDsGen));

             //if we reached this point we have found our gen match and we can stop the loop
             // since now we are not trying to find only 1 gen match for a reco-candidate
             // (we dont loop over the bs candidates) but rather over the events, we have to allow
             // several caniddates per event!! As we also allow several bs candidates per event!

             gen.addUserInt("sig",sigId);
             gen.addUserInt("b_mother_id",bId);
             gen.addUserInt("gen_match_success",genMatchSuccess);

             ret_value->emplace_back(gen);
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

iEvent.put(std::move(ret_value), "gen");
}//closing event loop


void inspectorGEN::endJob(){
// Printouts:
std::cout << "\n--------- GEN MATCHING MODULE ----------\n" << std::endl;
std::cout << "#Events in file                                           : " << nEventsGEN  << std::endl;
std::cout << "#Gen Muons in file                                        : " << nMuonsGEN   << std::endl;
std::cout << "#Gen Tracks in file                                       : " << nTracksGEN  << std::endl;
std::cout << "#Gen matched candidates                                   : " << nGenMatchedGEN << std::endl;
std::cout << " Before phi mass cut we have                              : " << beforePhiMass << std::endl;
std::cout << " After phi mass cut we have                              : " << afterPhiMass << std::endl;

std::cout << "\n#KKPiMu Gen combinations:" << nKKPiMuGEN << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a Phi         : " << nFoundPhiGEN  << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a Ds          : " << nFoundDsGEN   << std::endl;
std::cout << "#KKPiMu Gen combinations for which we found a B-mom       : " << nFoundBGEN    << std::endl;
std::cout << "#KKPiMu Gen combinations for which the B-mom < B mass cut : " << nBMassCutGEN << std::endl;
}


DEFINE_FWK_MODULE(inspectorGEN);
