#ifndef PhysicsTools_BParkingNano_helpers
#define PhysicsTools_BParkingNano_helpers

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryVector/interface/PV3DBase.h"
#include "Math/LorentzVector.h"

// for the fit
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"

// 4 vectors
#include "TLorentzVector.h"
#include "TVector3.h" // for boost vector

// for the cov matrix correction
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

// for the reco function
#include <tuple>

#include <vector>
#include <algorithm>
#include <limits>
#include <memory>

typedef std::vector<reco::TransientTrack> TransientTrackCollection;

constexpr float K_MASS = 0.493677;
constexpr float PI_MASS = 0.139571;



inline auto printDirectDaughters(const auto mom, bool print){

  for(size_t dauIdx = 0; dauIdx < mom->numberOfDaughters(); ++dauIdx){
    if (print){
    std::cout << "Mom is: "<< mom->pdgId() << std::endl;
    std::cout << "and has direct daughter: " << mom->daughter(dauIdx)->pdgId() << std::endl;
    }

    if((mom->daughter(dauIdx)->pdgId() == 22)){

    auto photon = mom->daughter(dauIdx);
    //photonEnergy = mom->daughter(dauIdx)->energy();
    return photon;
    }

  }
  const reco::Candidate* empty = nullptr;
  return empty;
}



///////////////////////////////////////////////////////////////////////////////////
// checks if B+ -> Ds K Nu Mu

inline int isKNuMu(const auto mom){

  int foundSignal = -1;

  std::set<int> k_dsmu         = {321, 431, 13, 14};
  std::set<int> k_dsstarmu     = {321, 433, 13, 14};
  std::set<int> k_dstau        = {321, 431, 15, 16};
  std::set<int> k_dsstartau    = {321, 433, 15, 16};

  std::set<int> dau_list;
 
  int nDaus = mom->numberOfDaughters();

  for(size_t dauIdx = 0; dauIdx < mom->numberOfDaughters(); ++dauIdx){
    dau_list.insert( abs(mom->daughter(dauIdx)->pdgId()) );
  }

  if      ((dau_list == k_dsmu)      && (nDaus == 4)) foundSignal = 107;
  else if ((dau_list == k_dsstarmu)  && (nDaus == 4)) foundSignal = 117;
  else if ((dau_list == k_dstau)     && (nDaus == 4)) foundSignal = 108;
  else if ((dau_list == k_dsstartau) && (nDaus == 4)) foundSignal = 118;

  return foundSignal;
}

///////////////////////////////////////////////////////////////////////////////////
// checks if signal 

inline int isSignal(const auto mom){

  int foundSignal = -1;

  std::set<int> signal_dsmu         = {431, 13, 14};
  std::set<int> signal_dsstarmu     = {433, 13, 14};
  std::set<int> signal_dstau        = {431, 15, 16};
  std::set<int> signal_dsstartau    = {433, 15, 16};

  std::set<int> dau_list;
 
  size_t nDaus = 0;

  for(size_t dauIdx = 0; dauIdx < mom->numberOfDaughters(); ++dauIdx){
 
    unsigned int dauId = abs(mom->daughter(dauIdx)->pdgId());
  
    if(dauId == 22){ 
      //std::cout << "found photon!!" << std::endl;
      //auto dummy = printDirectDaughters(mom, true);
      // These are soft-photons, f.e. FSR, which however are assigned 'promptly'
      // to the Bs. I.e. the decays often look like: Bs -> Ds + mu + nu + gamma + gamma
      // But nevertheless, this is still a signal, we checked that they are low energy.
      
      continue; 
    }
   
    else{
      dau_list.insert( dauId );
      ++nDaus;
    }

  }

  // now compare the daughters with the signal daughters
  // the comparison of the number of daughters is important to avoid decays where
  // we would have the same daughters but in a different multiplicity (even if I can
  // not think of any... but lets be sure! )

  if      ((dau_list == signal_dsmu)      && (nDaus == 3)) foundSignal = 0;
  else if ((dau_list == signal_dsstarmu)  && (nDaus == 3)) foundSignal = 10;
  else if ((dau_list == signal_dstau)     && (nDaus == 3)) foundSignal = 1;
  else if ((dau_list == signal_dsstartau) && (nDaus == 3)) foundSignal = 11;

  return foundSignal;
}

///////////////////////////////////////////////////////////////////////////////////
// function which prints all daughters 

inline void printDaughters(const auto mom){

  for(size_t dauIdx = 0; dauIdx < mom->numberOfDaughters(); ++dauIdx){
    std::cout << " Now mom is: "<< mom->pdgId() << std::endl;
    std::cout << "With daughter: " << mom->daughter(dauIdx)->pdgId() << std::endl;
    printDaughters(mom->daughter(dauIdx)); 
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////
// function which prints all moms 

inline void printMothers(const auto dau){

  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){
    std::cout << " Now dau is: "<< dau->pdgId() << std::endl;
    std::cout << " And has Nr of moms: "<< dau->numberOfMothers() << std::endl;
    std::cout << " Accessing mom: " << dau->mother(momIdx)->pdgId() << std::endl;
    printMothers(dau->mother(momIdx)); 
  }
  std::cout << "done... returning back" << std::endl;
  return;
}
///////////////////////////////////////////////////////////////////////////////////
// function which checks if a genParticle has a certain ancestor 

inline bool isAncestor(const auto dau, const int id){
  //std::cout << "pdgId = "<< dau->pdgId() << std::endl;
  if (fabs(dau->pdgId()) == id){ 
    return true;
    }
  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){
    if (isAncestor(dau->mother(momIdx), id)) return true;  
  }
  return false;
}
///////////////////////////////////////////////////////////////////////////////////
// function which checks if mom is ancestor of dau

inline bool hasAncestor(const auto dau, const auto mom){
  //std::cout << "pdgId = "<< dau->pdgId() << std::endl;
  if (dau == mom){ 
    return true;
    }
  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){
    if (hasAncestor(dau->mother(momIdx), mom)) return true;  
  }
  return false;
}
///////////////////////////////////////////////////////////////////////////////////
// function which returns pt eta phi of the ancestor in order to compare ancestors.

inline std::vector<double> infoAncestor(const auto dau, const int id){

  if (fabs(dau->pdgId()) == id){
    return {dau->pt(),dau->eta(),dau->phi(),dau->vx(),dau->vy(),dau->vz()};
  }

  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){
    if (isAncestor(dau->mother(momIdx), id)) {
      std::vector<double> ptEtaPhiVxVyVz = infoAncestor(dau->mother(momIdx),id);
      return ptEtaPhiVxVyVz;
    }
  }
  std::vector<double> failedVector(6, std::numeric_limits<double>::quiet_NaN());
  return failedVector;
}

///////////////////////////////////////////////////////////////////////////////////
// function which returns pointer to ancestor with pdgid <id> such that we can match by pointer! :)

inline auto getAncestor(const auto dau, const int id){

  //the pointer type changes when accessing moms, VERY ANNOYING IN A RECURSIVE FUNCTION
  //std::cout << "I am at pdg Id = " << dau->pdgId() << " and vertex vx = " << dau->vx() << std::endl; 
  if ((fabs(dau->pdgId()) == id)){
    //std::cout << "sucess!" << std::endl;
    return dau;
  }

  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){
    //std::cout << "Now I access mom Nr " << momIdx << std::endl;
    if (isAncestor(dau->mother(momIdx), id)) {
      auto dau2 = getAncestor(dau->mother(momIdx),id);
      return dau2;
    }
  }
  const reco::Candidate* empty = nullptr; 
  return empty;
}

inline int getDsID(auto pi){
  int dsID = 0; 
  if (isAncestor(pi, 431))   dsID = 431;   // Ds+ 
  if (isAncestor(pi, 433))   dsID = 433;   // Ds+* 
  if (isAncestor(pi, 10431)) dsID = 10431; // Ds*(2317)+ 
  if (isAncestor(pi, 20433)) dsID = 20433; // Ds*(2457)+
  return dsID;
}

inline int getSecondCharmID(auto mu){
  int dID = 0; 
  if (isAncestor(mu, 411))   dID = 411;   // D+
  if (isAncestor(mu, 421))   dID = 421;   // D0
  if (isAncestor(mu, 431))   dID = 431;   // Ds+
  
  if (isAncestor(mu, 413))   dID = 413;  // D+* 
  if (isAncestor(mu, 423))   dID = 423;  // D0* 
  if (isAncestor(mu, 433))   dID = 433;  // Ds*+ 
  if (isAncestor(mu, 4122))  dID = 4122; // Lambdac+ 
  return dID;
}
///////////////////////////////////////////////////////////////////////////////////
// function which removes the un-oscillated ancestor of dau
// f.e. dau is -531 and comes from 531 via oscillation, then this function removes oscillation
// and returns a pointer to the 531 particle, which has the correct vertex!

inline const reco::Candidate* removeOscillations(const auto dau){

  //std::cout << "I am at pdg Id = " << dau->pdgId() << " and vertex vx = " << dau->vx() << std::endl; 

  for(size_t momIdx = 0; momIdx < dau->numberOfMothers(); ++momIdx){

    //check if dau has a mother with the same pdg Id but opposite sign
    if (dau->mother(momIdx)->pdgId() == (-1 * dau->pdgId())) {
      //std::cout << "oscillation!" << std::endl;

      auto dau2 = removeOscillations(dau->mother(momIdx));
      return dau2;
    }
  }
  return dau;
}
///////////////////////////////////////////////////////////////////////////////////
//function which gets the hel angle between the mu and W
inline float angMuW (TLorentzVector d, TLorentzVector b, TLorentzVector mu){       

  //get q2
  TLorentzVector q = b - d;
  double q2 = q.M2();

  //boost Ds into Bs rest frame
  TVector3 bBoost = b.BoostVector();
  d.Boost(-bBoost);

  //get W via Ds 
  TLorentzVector w;
  w.SetVectM(-d.Vect(),std::sqrt(q2));           
  // boost it back into lab frame
  w.Boost(bBoost);
  // now take the boost vector of w
  TVector3 wBoost = w.BoostVector();
  //boost the muon into the w rest frame
  mu.Boost(-wBoost); 
  //boost the W back into the bs rest frame
  w.Boost(-bBoost);
 
  //now take the angle
  return w.Angle(mu.Vect());

}
///////////////////////////////////////////////////////////////////////////////////
//function which gets the hel angle between dau1 and dau2 in the rest frame of the restFrame particle

inline float angDoubleDecay (TLorentzVector restFrame, TLorentzVector dau1, TLorentzVector dau2) {
  TVector3 restBoost = restFrame.BoostVector();
  dau1.Boost(-restBoost);
  dau2.Boost(-restBoost);
  return dau1.Angle(dau2.Vect());
}

///////////////////////////////////////////////////////////////////////////////////
//function which gets the hel angle between dau1 and dau2 in the rest frame of the restFrame particle

inline float angDsPi (TLorentzVector d, TLorentzVector dau, TLorentzVector b) {

  TVector3 bBoost = b.BoostVector();
  TVector3 dBoost = d.BoostVector();
  dau.Boost(-dBoost);
  d.Boost(-bBoost);

  return d.Angle(dau.Vect());
}
///////////////////////////////////////////////////////////////////////////////////
//function which gets angle between decay planes

inline float angPlane (TLorentzVector d, TLorentzVector b, TLorentzVector mu, TLorentzVector pi) {

  //get q2
  TLorentzVector q = b - d;
  double q2 = q.M2();

  //get Ds boost vector in the lab frame (to boost pi)
  TVector3 dBoost = d.BoostVector();

  //boost Ds into Bs rest frame
  TVector3 bBoost = b.BoostVector();
  d.Boost(-bBoost);

  //get W via Ds 
  TLorentzVector w;
  w.SetVectM(-d.Vect(),std::sqrt(q2));           
  // boost it back into lab frame
  w.Boost(bBoost);
  // now take the boost vector of w
  TVector3 wBoost = w.BoostVector();
  //boost the muon into the w rest frame
  mu.Boost(-wBoost);
  //boost the W back into the bs rest frame
  w.Boost(-bBoost);

  //now boost the pi into the ds rest frame, ds is already in Bs rest frame
  pi.Boost(-dBoost);

  // normal vector on lepton-W plane
  TVector3 n1 = mu.Vect().Cross(w.Vect());
  // normal vector on ds-pi plane
  TVector3 n2 = d.Vect().Cross(pi.Vect()); 

  return n1.Angle(n2);
}

///////////////////////////////////////////////////////////////////////////////////
//function which gets angle of kaon in DsPi plane 

inline float angPlane2 (TLorentzVector d, TLorentzVector b, TLorentzVector k, TLorentzVector pi) {

  //get Ds boost vector in the lab frame (to boost pi)
  TVector3 dBoost = d.BoostVector();

  //now boost the pi and k1 into the ds rest frame
  pi.Boost(-dBoost);
  k.Boost(-dBoost);

  //boost Ds into Bs rest frame
  TVector3 bBoost = b.BoostVector();
  d.Boost(-bBoost);

  // normal vector on Ds - pi plane
  TVector3 n1 = d.Vect().Cross(pi.Vect());
  // normal vector on k1 - pi plane
  TVector3 n2 = k.Vect().Cross(pi.Vect()); 

  return n1.Angle(n2);
}

///////////////////////////////////////////////////////////////////////////////////
//function which returns E*, the momentum of mu in b rest frame

inline float getEStar(TLorentzVector b, TLorentzVector mu){

  //get b boost vector
  TVector3 bBoost = b.BoostVector();
  mu.Boost(-bBoost);
  
  return mu.E();

}

///////////////////////////////////////////////////////////////////////////////////
//function which returns E*, the momentum of mu in b rest frame

inline float getEGamma(TLorentzVector ds, const double dsMass_, const double dsStarMass_ ){

  // collienar approx to get Ds star 
  TLorentzVector dsStar = ds;
  dsStar *= dsStarMass_ / dsMass_ ; 

  float eGamma = std::sqrt( std::pow(dsStarMass_,2) * (1 + std::pow(ds.P() / dsMass_ ,2) )) - ds.E();
 
  return eGamma;

}



///////////////////////////////////////////////////////////////////////////////////
//function which returns the phi difference of two phi variables (in cms coordinate system)
inline double phiDiff(double phi1, double phi2){

  double dPhi = phi1 - phi2;
  double pi = 3.14159265358979323846;
  while (fabs(dPhi) > pi) {
    int sgn = dPhi > 0? 1 : -1;
    dPhi -= sgn*2*pi;
  }
  return dPhi;
}

///////////////////////////////////////////////////////////////////////////////////
//function which returns the phi difference of two phi variables (in cms coordinate system)
inline TLorentzVector collMethod(TLorentzVector dMu, const double bMass_){

        TLorentzVector b = dMu;

        double dMuMass = dMu.M();
        b *= bMass_ / dMuMass; //scale it

        return b;

}

///////////////////////////////////////////////////////////////////////////////////
// LHCb method for bs reconstruction

inline TLorentzVector lhcbMethod(TLorentzVector dMu, float v1_x, float v1_y, float v1_z, float v2_x, float v2_y ,float v2_z, const double bMass_){

  TVector3 bsFlightDir;
  TVector3 beamAxis;
  TVector3 radialAxis;
           
  bsFlightDir.SetXYZ(v2_x - v1_x, v2_y - v1_y , v2_z - v1_z);
  
  beamAxis.SetXYZ(0.0,0.0,1.0); // in z direction
  radialAxis.SetXYZ(1.0,0.0,0.0); //in x direction
   
  float theta = bsFlightDir.Theta(); //angle between beam axis and bs flight dir
  float lhcbPz = dMu.Pz() *bMass_ / dMu.M();
  float lhcbPt = lhcbPz * std::tan(theta); //angle is in radians! std::tan also takes radians :)
  
  bsFlightDir.SetZ(0.0); //project on xy plane for phi calculation
  float lhcbPhi;

  // give attention that the phi component is correct
  if (bsFlightDir.Py() > 0) lhcbPhi = bsFlightDir.Angle(radialAxis); //get the phi angle
  else lhcbPhi = - bsFlightDir.Angle(radialAxis); //get the phi angle
 
  TLorentzVector b; 
  float eta = - std::log(std::tan(theta/2));
  b.SetPtEtaPhiM(lhcbPt,eta,lhcbPhi,bMass_); 
     
  return b;

} 

///////////////////////////////////////////////////////////////////////////////////
// Alternative LHCb method for bs reconstruction

inline TLorentzVector lhcbAltMethod(TLorentzVector dMu, float v1_x, float v1_y, float v1_z, float v2_x, float v2_y ,float v2_z, const double bMass_){

  TVector3 bsFlightDir;
  bsFlightDir.SetXYZ(v2_x - v1_x, v2_y - v1_y , v2_z - v1_z);
  
  TVector3 lhcbAltBs;
  TLorentzVector lhcbAltBsTlv;
  
  lhcbAltBs = bsFlightDir.Unit();  
  lhcbAltBs *= dMu.Vect().Mag() * bMass_ / dMu.M(); 
  lhcbAltBsTlv.SetXYZM(lhcbAltBs.Px(),lhcbAltBs.Py(),lhcbAltBs.Pz(), bMass_);      
  
  return lhcbAltBsTlv;

}

inline std::tuple<std::vector<TLorentzVector>, float> recoMethod(TLorentzVector dMu, float v1_x, float v1_y, float v1_z, float v2_x, float v2_y ,float v2_z, const double bMass_){

  TLorentzVector recoBsTlv1;
  TLorentzVector recoBsTlv2;
  
  double recoNeutPll_1; // neutrino momentum parallel to bs direction
  double recoNeutPll_2; // "
  
  double recoBsAbs_1; // absolute 3 momentum of bs
  double recoBsAbs_2; // "
  
  // bs flight direction
  TVector3 bsFlightDir;
  bsFlightDir.SetXYZ(v2_x - v1_x, v2_y - v1_y , v2_z - v1_z);
   
  // angle between the bs flight direction and the Dsmu system (visible)
  double recoAngle = dMu.Angle(bsFlightDir); 
                      
  //define parameters 

  // momentum of DsMu system parallel to the bs
  double recoDsMuPll = std::cos(recoAngle) * dMu.Vect().Mag();
  // momentum of DsMu system orthogonal to the bs
  double recoDsMuT = std::sin(recoAngle) * dMu.Vect().Mag();
  // energy of Dsmu system
  double recoDsMuE = dMu.E(); 
  // cocktail, drops out of equation
  double recoMix = std::pow(bMass_,2) + std::pow(recoDsMuPll,2) - std::pow(recoDsMuT,2) - std::pow(recoDsMuE,2);
  
  // define a,b,c, to give to mitternachtsformel
  double a = 4*(std::pow(recoDsMuPll,2) - std::pow(recoDsMuE,2));
  double b = 4*recoDsMuPll*recoMix;
  double c = std::pow(recoMix,2) - 4*std::pow(recoDsMuE,2)*std::pow(recoDsMuT,2);

  // discriminant
  double disc = std::pow(b,2) - 4*a*c;      
  float discNegativity;
  /*
  std::cout << "v1_x   = " << v1_x << std::endl;
  std::cout << "v1_y   = " << v1_y << std::endl;
  std::cout << "v1_z   = " << v1_z << std::endl;

  std::cout << "v2_x   = " << v2_x << std::endl;
  std::cout << "v2_y   = " << v2_y << std::endl;
  std::cout << "v2_z   = " << v2_z << std::endl;
  std::cout << "dMu Pt = " << dMu.Pt() << std::endl;
  */
 
  if( disc >= 0) {
  
    //non complex root -> nice! 
    recoNeutPll_1 = (-b + std::sqrt(disc)) / (2*a);
    recoNeutPll_2 = (-b - std::sqrt(disc)) / (2*a);
    //std::cout << "---- disc is positive!!" << std::endl;
    //std::cout << "disc = " << disc << std::endl;
    discNegativity = 0;
 
  }
  else{
    // complex root, only save -b / 2a (i.e. set disc to zero!)

    recoNeutPll_1 = -b / (2*a);
    recoNeutPll_2 = -b / (2*a);
    //std::cout << "---- disc is negative!!" << std::endl;
    //std::cout << "disc = " << disc << std::endl;
    discNegativity = std::sqrt(abs(disc)) / abs(b);
  }
 
  recoBsAbs_1 = recoDsMuPll + recoNeutPll_1; 
  recoBsAbs_2 = recoDsMuPll + recoNeutPll_2; 
  
  TVector3 recoBs_1 = bsFlightDir.Unit();
  TVector3 recoBs_2 = bsFlightDir.Unit();
  
  recoBs_1 *= recoBsAbs_1; 
  recoBs_2 *= recoBsAbs_2; 
  
  recoBsTlv1.SetXYZM(recoBs_1.Px(),recoBs_1.Py(),recoBs_1.Pz(),bMass_);
  recoBsTlv2.SetXYZM(recoBs_2.Px(),recoBs_2.Py(),recoBs_2.Pz(),bMass_);

  std::vector<TLorentzVector> recos;
  recos.push_back(recoBsTlv1);
  recos.push_back(recoBsTlv2);
  
  return std::make_tuple(recos,discNegativity); 
  
}


///////////////////////////////////////////////////////////////////////////////////
// Fix the track covariance matrix to be pos. def.


inline reco::Track correctCovMat(const reco::Track *tk, double delta){

// Parameters associated to the 5D curvilinear covariance matrix: 
// (qoverp, lambda, phi, dxy, dsz) 
// Defined as:
//   qoverp = q / abs(p) = signed inverse of momentum [1/GeV] 
//   lambda = pi/2 - polar angle at the given point 
//   phi = azimuth angle at the given point 
//   dxy = -vx*sin(phi) + vy*cos(phi) [cm] 
//   dsz = vz*cos(lambda) - (vx*cos(phi)+vy*sin(phi))*sin(lambda) [cm] 

    unsigned int i;
    unsigned int j;
    double min_eig = 1;

    // Get the original covariance matrix
    reco::TrackBase::CovarianceMatrix cov = tk->covariance();

    // Define a TMatrixDSym of the same shape as the old cov matrix.
    // Sym -> you only have to give one dimension, it will be a symm matrix
    // of shape (cov.kRows, cov.kRows)
    TMatrixDSym newCov(cov.kRows);

    // loop over old cov matrix 
    for (i = 0; i < cov.kRows; i++) {
        for (j = 0; j < cov.kRows; j++) {
            // change nan or inf values to 1e-6 
            if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
                cov(i,j) = 1e-6;
            // fill new covariacne matrix
            newCov(i,j) = cov(i,j);
        }
    }

    // Define a vector of size cov.kRows
    TVectorD eig(cov.kRows);
    // Fill it with the egienvalues of the newCov
    newCov.EigenVectors(eig);

    // loop over eigenvalues and find the minimal eigenvalue :)
    for (i = 0; i < cov.kRows; i++)
        if (eig(i) < min_eig)
            min_eig = eig(i);

    // If the minimum eigenvalue is less than zero, then subtract it from the diagonal and add `delta`.
    if (min_eig < 0) {
        for (i = 0; i < cov.kRows; i++)
            cov(i,i) -= min_eig - delta;
    }

    return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
    }

///////////////////////////////////////////////////////////////////////////////////
// Fix the track covariance matrix to be pos. def.
inline reco::Track fixTrack(const reco::TrackRef& tk)
{
    reco::Track t = reco::Track(*tk);
    return correctCovMat(&t, 1e-8);
}


///////////////////////////////////////////////////////////////////////////////////
// Vertex Fit 
inline RefCountedKinematicTree vertexFit(std::vector<RefCountedKinematicParticle> toFit, ParticleMass massConstr, bool applyConstr)
{

  //define fitter
  KinematicConstrainedVertexFitter fitter;

  //define constraint
  MultiTrackKinematicConstraint* constr = new TwoTrackMassKinematicConstraint(massConstr);

  RefCountedKinematicTree fitTree;

  if(applyConstr) {
  // perform the fit
  fitTree = fitter.fit(toFit, constr);
  }
  else{
  fitTree = fitter.fit(toFit);
  }

  return fitTree;

}
#endif
