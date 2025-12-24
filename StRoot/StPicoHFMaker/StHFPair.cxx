#include <limits>
#include <cmath>
#include <iostream>

#include "StHFPair.h"

#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"

using namespace std;

ClassImp(StHFPair)

// _________________________________________________________
StHFPair::StHFPair(): mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(std::numeric_limits<unsigned short>::max()), mParticle2Idx(std::numeric_limits<unsigned short>::max()),
  mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()),
  mThetaProdPlane(std::numeric_limits<float>::quiet_NaN()), mProdPlane(TVector3()) {
}

// _________________________________________________________
StHFPair::StHFPair(StHFPair const * t) : mLorentzVector(t->mLorentzVector), mDecayVertex(t->mDecayVertex),
   mPointingAngle(t->mPointingAngle), mDecayLength(t->mDecayLength),
   mParticle1Dca(t->mParticle1Dca), mParticle2Dca(t->mParticle2Dca),
   mParticle1Idx(t->mParticle1Idx), mParticle2Idx(t->mParticle2Idx),
   mDcaDaughters(t->mDcaDaughters), mCosThetaStar(t->mCosThetaStar),
   mThetaProdPlane(t->mThetaProdPlane), mProdPlane(t->mProdPlane) {
}

// _________________________________________________________
StHFPair::StHFPair(StPicoTrack const * const particle1, StPicoTrack const * const particle2,
		   float p1MassHypo, float p2MassHypo, unsigned short const p1Idx, unsigned short const p2Idx,
		   TVector3 const & vtx, float const bField, bool const useStraightLine) : 
  mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),
  mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()) {
  // -- Create pair out of 2 tracks
  //     prefixes code:
  //      p1 means particle 1
  //      p2 means particle 2
  //      pair means particle1-particle2  pair

  if ((!particle1 || !particle2) || (particle1->id() == particle2->id())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    return;
  }

	StPicoPhysicalHelix p1Helix = particle1->helix(bField); 
  StPicoPhysicalHelix p2Helix = particle2->helix(bField); //bFiled not in kilogauss - is properly computed inside helix(double B) function in StPicoTrack.h

  // -- move origins of helices to the primary vertex origin
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));

  // -- use straight lines approximation to get point of DCA of particle1-particle2 pair
  TVector3 const p1Mom = p1Helix.momentum(bField * kilogauss);
  TVector3 const p2Mom = p2Helix.momentum(bField * kilogauss);
  StPicoPhysicalHelix const p1StraightLine(p1Mom, p1Helix.origin(), 0, particle1->charge());
  StPicoPhysicalHelix const p2StraightLine(p2Mom, p2Helix.origin(), 0, particle2->charge());

  pair<double, double> const ss = (useStraightLine) ? p1StraightLine.pathLengths(p2StraightLine) : p1Helix.pathLengths(p2Helix);
  
  if(useStraightLine)
  {
    TVector3 const p1AtDcaToP2 = p1StraightLine.at(ss.first);
    TVector3 const p2AtDcaToP1 = p2StraightLine.at(ss.second);  
    
    // -- calculate DCA of particle1 to particl2 at their DCA
    mDcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag();
    
    // -- calculate decay vertex (secondary) 
    mDecayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
  }
  else
  {
    TVector3 const p1AtDcaToP2 = p1Helix.at(ss.first);
    TVector3 const p2AtDcaToP1 = p2Helix.at(ss.second);  
    
    // -- calculate DCA of particle1 to particl2 at their DCA
    mDcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag();
    
    // -- calculate decay vertex (secondary) 
    mDecayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
  }

  // -- calculate Lorentz vector of particle1-particle2 pair
  TVector3 const p1MomAtDca = p1Helix.momentumAt(ss.first,  bField * kilogauss);
  TVector3 const p2MomAtDca = p2Helix.momentumAt(ss.second, bField * kilogauss);

  TLorentzVector const p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2() + p1MassHypo*p1MassHypo));
  TLorentzVector const p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2() + p2MassHypo*p2MassHypo));

  mLorentzVector = p1FourMom + p2FourMom;

  // -- calculate cosThetaStar
  TLorentzVector const pairFourMomReverse(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz(), mLorentzVector.E());
  TLorentzVector p1FourMomStar = p1FourMom;
  //p1FourMomStar.Boost(pairFourMomReverse.Vect());
  p1FourMomStar.Boost(pairFourMomReverse.BoostVector());  
  mCosThetaStar = std::cos(p1FourMomStar.Vect().Angle(mLorentzVector.Vect()));

  TVector3 beamVector(0.,0.,1.); //unity vector along the beam axis
  TVector3 mProdPlane_work = beamVector.Cross(mLorentzVector.Vect());
  mProdPlane = ( mProdPlane_work )*(1./mProdPlane_work.Mag() );

  mThetaProdPlane = mProdPlane.Angle(p1FourMomStar.Vect());

  //cout<<mThetaProdPlane<<endl; 
 
  // -- calculate pointing angle and decay length with respect to primary vertex 
  //    if decay vertex is a tertiary vertex
  //    -> only rough estimate -> needs to be updated after secondary vertex is found
  TVector3 const vtxToV0 = mDecayVertex - vtx;
  mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());
  mDecayLength = vtxToV0.Mag();
  mDcaToPrimaryVertex = mDecayLength*sin(mPointingAngle); // sine law: DcaToPrimaryVertex/sin(pointingAngle) = decayLength/sin(90°)

  // -- calculate DCA of tracks to primary vertex
  //    if decay vertex is a tertiary vertex
  //    -> only rough estimate -> needs to be updated after secondary vertex is found
  mParticle1Dca = (p1Helix.origin() - vtx).Mag();
  mParticle2Dca = (p2Helix.origin() - vtx).Mag();
}

// _________________________________________________________
StHFPair::StHFPair(StPicoTrack const * const particle1, StHFPair const * const particle2,
		   float p1MassHypo, float p2MassHypo, unsigned short const p1Idx, unsigned short const p2Idx,
		   TVector3 const & vtx, float const bField, bool const useStraightLine) :
  mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),
  mDcaDaughters(std::numeric_limits<float>::max()), mCosThetaStar(std::numeric_limits<float>::quiet_NaN()) {
  // -- Create pair out of a particle and and a pair
  //     prefixes code:
  //      p1 means particle 1
  //      p2 means particle 2 - which is a pair of tertiaryP1 and tertiaryP2
  //      pair means particle1-particle2  pair
  
  // -- assume incoming pair always to be neutral 
  float p2Charge = 0.;

  // -- checking that particle1 and tertiary particle (=pair) exists, and particle1 id different from both used to reconstruct tertiary pair
  if ((!particle1 || !particle2) || (particle1->id() == particle2->particle1Idx()) || ( particle1->id() == particle2->particle2Idx())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    return;
  }

  StPicoPhysicalHelix p1Helix = particle1->helix(bField);
  TVector3 const p2Mom(particle2->px(),particle2->py(),particle2->pz());
  TVector3 const p2Origin(particle2->v0x(),particle2->v0y(),particle2->v0z());
  StPicoPhysicalHelix p2Helix(p2Mom,p2Origin,bField*kilogauss, p2Charge);

  // --move origins of helices to the primary vertex origin
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));

  // --use straight lines approximation to get point of DCA of particle1-particle2 pair
  TVector3 const p1Mom = p1Helix.momentum(bField * kilogauss);
  StPicoPhysicalHelix const p1StraightLine(p1Mom, p1Helix.origin(), 0, particle1->charge());
  StPicoPhysicalHelix const p2StraightLine(p2Mom, p2Helix.origin(), 0, p2Charge);
  
  //pair<double, double> const ss = (useStraightLine) ? p1StraightLine.pathLengths(p2StraightLine) : p1Helix.pathLengths(p2Helix);
  
  //first, uses straight line approximation to find approximate position of SV
  pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine); 
  
  // -- momenta at DCA_pair
  TVector3 p1MomAtDca;
  TVector3 p2MomAtDca;
  
    
  if(useStraightLine)
  {
    //use SV from straight line approximation
    
    TVector3 const p1AtDcaToP2 = p1StraightLine.at(ss.first);
    TVector3 const p2AtDcaToP1 = p2StraightLine.at(ss.second);  
    
    // -- calculate DCA of particle1 to particl2 at their DCA
    mDcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag();
    
    // -- calculate decay vertex (secondary) 
    mDecayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
    
    // -- get momenta at DCA_pair
    p1MomAtDca = p1Helix.momentumAt(ss.first,  bField * kilogauss);
    p2MomAtDca = p2Helix.momentumAt(ss.second, bField * kilogauss);
  }
  else
  {
    //find exact position of SV using helixes
    //use SV from straight line approxiamtion as a staring point
    
    //make a copy of particle helixes
    StPicoPhysicalHelix p1Helix_SV = p1Helix;
    StPicoPhysicalHelix p2Helix_SV = p2Helix;
    
    //move their origin to SV from straight line approxiamtion
    p1Helix_SV.moveOrigin(ss.first);
    p2Helix_SV.moveOrigin(ss.second);
    
    //get new DCA_pair using helixes
    pair<double, double> const ss_helix = p1Helix_SV.pathLengths(p2Helix_SV);
    
    TVector3 const p1AtDcaToP2 = p1Helix_SV.at(ss_helix.first);
    TVector3 const p2AtDcaToP1 = p2Helix_SV.at(ss_helix.second);  
    
    // -- calculate DCA of particle1 to particl2 at their DCA
    mDcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).Mag();
    
    // -- calculate decay vertex (secondary) 
    mDecayVertex = (p1AtDcaToP2 + p2AtDcaToP1) * 0.5 ;
    
    // -- get momenta at DCA_pair    
    p1MomAtDca = p1Helix.momentumAt(ss_helix.first,  bField * kilogauss);
    p2MomAtDca = p2Helix.momentumAt(ss_helix.second, bField * kilogauss);
  }
    
  // -- calculate Lorentz vector of particle1-particle2 pair
  
  TLorentzVector const p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2() + p1MassHypo*p1MassHypo));
  TLorentzVector const p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2() + p2MassHypo*p2MassHypo));
  
  mLorentzVector = p1FourMom + p2FourMom;

  // -- calculate cosThetaStar
  TLorentzVector const pairFourMomReverse(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz(), mLorentzVector.E());
  TLorentzVector p1FourMomStar = p1FourMom;
  p1FourMomStar.Boost(pairFourMomReverse.BoostVector()); //particle1 momentum in mother rest frame
  mCosThetaStar = std::cos(p1FourMomStar.Vect().Angle(mLorentzVector.Vect()));


  TVector3 beamVector(0.,0.,1.); //unity vector along the beam axis
  TVector3 mProdPlane_work = beamVector.Cross(mLorentzVector.Vect());
  mProdPlane = ( mProdPlane_work )*(1./mProdPlane_work.Mag() );

  mThetaProdPlane = mProdPlane.Angle(p1FourMomStar.Vect());
  

  // -- calculate pointing angle and decay length with respect to primary vertex
  TVector3 const vtxToV0 = mDecayVertex - vtx;
  mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());
  mDecayLength = vtxToV0.Mag();
   
  // -- calculate DCA of tracks to primary vertex
  mParticle1Dca = (p1Helix.origin() - vtx).Mag();
  mParticle2Dca = (p2Helix.origin() - vtx).Mag();
}
// _________________________________________________________
float StHFPair::pointingAngle(TVector3 const & vtx2) const{
  // -- Overloaded function recalculates pointing angle given secondary vertex
  TVector3 const vtx2ToTertiary(mDecayVertex - vtx2);
  float const nPointingAngle = vtx2ToTertiary.Angle(mLorentzVector.Vect());
  return nPointingAngle;
}
// _________________________________________________________
float StHFPair::decayLength(TVector3 const & vtx2) const{
  // -- Overloaded function recalculates decayLength given secondary vertex
  TVector3 const vtx2ToTertiary(mDecayVertex - vtx2); 
  float const nDecayLength = vtx2ToTertiary.Mag();  
  return nDecayLength;
}
// _________________________________________________________
float StHFPair::particle1Dca(StPicoTrack const * p1track, TVector3 const & vtx2, float bField) const{
  // -- Overloaded function recalculates daughter dca 2 updated vertex
  StPicoPhysicalHelix p1Helix = p1track->helix(bField);
  // -- move origins of helices to the primary vertex origin
  p1Helix.moveOrigin(p1Helix.pathLength(vtx2));

  float const nParticle1Dca = (p1Helix.origin() - vtx2).Mag();
  return nParticle1Dca;
}
// _________________________________________________________
float StHFPair::particle2Dca(StPicoTrack const * p2track, TVector3 const & vtx2, float bField) const{
  // -- Overloaded function recalculates daughter dca 2 updated vertex
  StPicoPhysicalHelix p2Helix = p2track->helix(bField);
  // -- move origins of helices to the primary vertex origin
  p2Helix.moveOrigin(p2Helix.pathLength(vtx2));

  float const nParticle2Dca = (p2Helix.origin()  - vtx2).Mag();
  return nParticle2Dca;
}


