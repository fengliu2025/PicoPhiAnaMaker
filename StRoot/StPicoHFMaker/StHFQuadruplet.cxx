#include <limits>
#include <cmath>

#include "StHFQuadruplet.h"

#include "StPicoEvent/StPicoPhysicalHelix.h"
//#include "StarClassLibrary/StPicoPhysicalHelix.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"

ClassImp(StHFQuadruplet)

// _________________________________________________________
StHFQuadruplet::StHFQuadruplet(): mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle3Dca(std::numeric_limits<float>::quiet_NaN()),  mParticle4Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(std::numeric_limits<unsigned short>::max()), mParticle2Idx(std::numeric_limits<unsigned short>::max()),
  mParticle3Idx(std::numeric_limits<unsigned short>::max()),mParticle4Idx(std::numeric_limits<unsigned short>::max()),
  mDcaDaughters12(std::numeric_limits<float>::max()),  mDcaDaughters13(std::numeric_limits<float>::max()),
  mDcaDaughters14(std::numeric_limits<float>::max()),  mDcaDaughters23(std::numeric_limits<float>::max()),
  mDcaDaughters24(std::numeric_limits<float>::max()),  mDcaDaughters34(std::numeric_limits<float>::max()),
  mCosThetaStar(std::numeric_limits<float>::quiet_NaN())
{
}

// _________________________________________________________
StHFQuadruplet::StHFQuadruplet(StHFQuadruplet const * t) : 
  mLorentzVector(t->mLorentzVector), mDecayVertex(t->mDecayVertex),
  mPointingAngle(t->mPointingAngle), mDecayLength(t->mDecayLength), 
  mParticle1Dca(t->mParticle1Dca), mParticle2Dca(t->mParticle2Dca), mParticle3Dca(t->mParticle3Dca),mParticle4Dca(t->mParticle4Dca),
  mParticle1Idx(t->mParticle1Idx), mParticle2Idx(t->mParticle2Idx), mParticle3Idx(t->mParticle3Idx),mParticle4Idx(t->mParticle4Idx),
  mDcaDaughters12(t->mDcaDaughters12), mDcaDaughters13(t->mDcaDaughters13), 
  mDcaDaughters14(t->mDcaDaughters14), mDcaDaughters23(t->mDcaDaughters23), 
  mDcaDaughters24(t->mDcaDaughters24), mDcaDaughters34(t->mDcaDaughters34), 
  mCosThetaStar(t->mCosThetaStar)
{
}
// _________________________________________________________
StHFQuadruplet::StHFQuadruplet(StPicoTrack const * const particle1, StPicoTrack const * const particle2, StPicoTrack const * const particle3, StPicoTrack const * particle4,
			       float p1MassHypo, float p2MassHypo, float p3MassHypo,float p4MassHypo,
			       unsigned short const p1Idx, unsigned short const p2Idx, unsigned short const p3Idx, unsigned short const p4Idx,
			 TVector3 const & vtx, float const bField)  : 
  mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle3Dca(std::numeric_limits<float>::quiet_NaN()),  mParticle4Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),  mParticle3Idx(p3Idx),mParticle4Idx(p4Idx),
  mDcaDaughters12(std::numeric_limits<float>::max()), mDcaDaughters13(std::numeric_limits<float>::max()),  
  mDcaDaughters14(std::numeric_limits<float>::max()),  mDcaDaughters23(std::numeric_limits<float>::max()),
  mDcaDaughters24(std::numeric_limits<float>::max()),  mDcaDaughters34(std::numeric_limits<float>::max()),
  mCosThetaStar(std::numeric_limits<float>::min())
{
  // -- Create quadruplet out of 4 tracks
  //     prefixes code:
  //      p1 means particle 1
  //      p2 means particle 2
  //      p3 means particle 3
  //      p4 means particle 4
  //      pair means particle1-particle2 pair||  particle2-particle3 pair ||  particle1-particle3 pair
  //      quadruplet particle1-particle2-particle3-particle4

  if ((!particle1 || !particle2 || !particle3 || !particle4) || 
      (particle1->id() == particle2->id() || particle1->id() == particle3->id() || particle1->id() == particle4->id() || particle2->id() == particle3->id() || particle2->id() == particle4->id() || particle3->id() == particle4->id())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    mParticle3Idx = std::numeric_limits<unsigned short>::max();
    mParticle4Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
	
	// -- use straight lines approximation to get point of DCA of particle1-particle2 pair
	//SL16j, Vanek
	StPicoPhysicalHelix p1Helix = particle1->helix(bField); //used at line 159
  StPicoPhysicalHelix p2Helix = particle2->helix(bField);
  StPicoPhysicalHelix p3Helix = particle3->helix(bField);
  StPicoPhysicalHelix p4Helix = particle4->helix(bField);

  // -- move origins of helices to the primary vertex origin - neede later in the code
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));
  p3Helix.moveOrigin(p3Helix.pathLength(vtx));
  p4Helix.moveOrigin(p4Helix.pathLength(vtx));
    
  TVector3 const p1Mom = particle1->gMom(); //momentum at pvtx,  computed inside StPicoTrack
  TVector3 const p2Mom = particle2->gMom();
  TVector3 const p3Mom = particle3->gMom();
  TVector3 const p4Mom = particle4->gMom();
  
  StPicoPhysicalHelix const p1StraightLine(p1Mom, particle1->origin(), 0, particle1->charge()); //origin computed iside StPicoTrack
  StPicoPhysicalHelix const p2StraightLine(p2Mom, particle2->origin(), 0, particle2->charge());
  StPicoPhysicalHelix const p3StraightLine(p3Mom, particle3->origin(), 0, particle3->charge());
  StPicoPhysicalHelix const p4StraightLine(p4Mom, particle4->origin(), 0, particle4->charge());

//-----------------------------------------------------------------------------------------

  pair<double, double> const ss12 = p1StraightLine.pathLengths(p2StraightLine);
  TVector3 const p1AtDcaToP2 = p1StraightLine.at(ss12.first);
  TVector3 const p2AtDcaToP1 = p2StraightLine.at(ss12.second);

  pair<double, double> const ss13 = p1StraightLine.pathLengths(p3StraightLine);
  TVector3 const p1AtDcaToP3 = p1StraightLine.at(ss13.first);
  TVector3 const p3AtDcaToP1 = p3StraightLine.at(ss13.second);

  pair<double, double> const ss14 = p1StraightLine.pathLengths(p4StraightLine);
  TVector3 const p1AtDcaToP4 = p1StraightLine.at(ss14.first);
  TVector3 const p4AtDcaToP1 = p4StraightLine.at(ss14.second);

  pair<double, double> const ss23 = p2StraightLine.pathLengths(p3StraightLine);
  TVector3 const p2AtDcaToP3 = p2StraightLine.at(ss23.first);
  TVector3 const p3AtDcaToP2 = p3StraightLine.at(ss23.second);

  pair<double, double> const ss24 = p2StraightLine.pathLengths(p4StraightLine);
  TVector3 const p2AtDcaToP4 = p2StraightLine.at(ss24.first);
  TVector3 const p4AtDcaToP2 = p4StraightLine.at(ss24.second);
  
  pair<double, double> const ss34 = p3StraightLine.pathLengths(p4StraightLine);
  TVector3 const p3AtDcaToP4 = p3StraightLine.at(ss34.first);
  TVector3 const p4AtDcaToP3 = p4StraightLine.at(ss34.second);
  
  // -- calculate DCA of particle1 to particl2 at their DCA
  mDcaDaughters12 = (p1AtDcaToP2 - p2AtDcaToP1).Mag();

  // -- calculate DCA of particle1 to particl3 at their DCA
  mDcaDaughters13 = (p1AtDcaToP3 - p3AtDcaToP1).Mag();
  
  // -- calculate DCA of particle1 to particl4 at their DCA
  mDcaDaughters14 = (p1AtDcaToP4 - p4AtDcaToP1).Mag();

  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters23 = (p2AtDcaToP3 - p3AtDcaToP2).Mag();

  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters24 = (p2AtDcaToP4 - p4AtDcaToP2).Mag();
  
  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters34 = (p3AtDcaToP4 - p4AtDcaToP3).Mag();

  
  // -- calculate decay vertex (secondary)
  TVector3 mDecayVertex = ( p1AtDcaToP2 + p2AtDcaToP1 + p1AtDcaToP3 + p3AtDcaToP1 + p1AtDcaToP4 + p4AtDcaToP1 + p2AtDcaToP3 + p3AtDcaToP2 + p2AtDcaToP4 + p4AtDcaToP2 +p3AtDcaToP4 + p4AtDcaToP3)*(1./12.0);
  
  // -- constructing mother daughter four momentum. Need helix (not straight line) for each daughter
  double const p1AtV0 = p1Helix.pathLength( mDecayVertex );
  TVector3 const p1MomAtDca = p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss);

  double const p2AtV0 = p2Helix.pathLength( mDecayVertex );
  TVector3 const p2MomAtDca = p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss);
  
  double const p3AtV0 = p3Helix.pathLength( mDecayVertex );
  TVector3 const p3MomAtDca = p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss);

  double const p4AtV0 = p4Helix.pathLength( mDecayVertex );
  TVector3 const p4MomAtDca = p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss);
  
  TLorentzVector const p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2() + p1MassHypo*p1MassHypo));
  TLorentzVector const p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2() + p2MassHypo*p2MassHypo));
  TLorentzVector const p3FourMom(p3MomAtDca, sqrt(p3MomAtDca.Mag2() + p3MassHypo*p3MassHypo));
  TLorentzVector const p4FourMom(p4MomAtDca, sqrt(p4MomAtDca.Mag2() + p4MassHypo*p4MassHypo));
  
  mLorentzVector = p1FourMom + p2FourMom + p3FourMom + p4FourMom;
   
  // -- calculate cosThetaStar
  //    ->> Lomnitz: Need to rethink theta star
  TLorentzVector const quadrupletFourMomReverse(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz(), mLorentzVector.E());
  TLorentzVector p1FourMomStar = p1FourMom;
  p1FourMomStar.Boost(quadrupletFourMomReverse.Vect());
  mCosThetaStar = std::cos(p1FourMomStar.Vect().Angle(mLorentzVector.Vect()));
  
  // -- calculate pointing angle and decay length
  TVector3 const vtxToV0 = mDecayVertex - vtx;
  mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());
  mDecayLength = vtxToV0.Mag();

  // --- calculate DCA of tracks to primary vertex
  mParticle1Dca = (p1Helix.origin() - vtx).Mag();
  mParticle2Dca = (p2Helix.origin() - vtx).Mag();
  mParticle3Dca = (p3Helix.origin() - vtx).Mag();
  mParticle4Dca = (p4Helix.origin() - vtx).Mag();

}
// _________________________________________________________
StHFQuadruplet::StHFQuadruplet(StPicoTrack const * const particle1, StPicoTrack const * const particle2, StPicoTrack const * const particle3, StHFPair const * particle4,
			       float p1MassHypo, float p2MassHypo, float p3MassHypo,float p4MassHypo,
			       unsigned short const p1Idx, unsigned short const p2Idx, unsigned short const p3Idx, unsigned short const p4Idx,
			 TVector3 const & vtx, float const bField)  : 
  mLorentzVector(TLorentzVector()), mDecayVertex(TVector3()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle3Dca(std::numeric_limits<float>::quiet_NaN()),  mParticle4Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),  mParticle3Idx(p3Idx),mParticle4Idx(p4Idx),
  mDcaDaughters12(std::numeric_limits<float>::max()), mDcaDaughters13(std::numeric_limits<float>::max()),  
  mDcaDaughters14(std::numeric_limits<float>::max()),  mDcaDaughters23(std::numeric_limits<float>::max()),
  mDcaDaughters24(std::numeric_limits<float>::max()),  mDcaDaughters34(std::numeric_limits<float>::max()),
  mCosThetaStar(std::numeric_limits<float>::min())
{
  // -- Create quadruplet out of 4 tracks
  //     prefixes code:
  //      p1 means particle 1
  //      p2 means particle 2
  //      p3 means particle 3
  //      p4 means particle 4
  //      pair means particle1-particle2 pair||  particle2-particle3 pair ||  particle1-particle3 pair
  //      triplet particle1-particle2-particle3

  if ((!particle1 || !particle2 || !particle3 || !particle4) || 
      (particle1->id() == particle2->id() || particle1->id() == particle3->id() || particle1->id() == (particle4->particle1Idx() || particle4->particle2Idx()) || particle2->id() == particle3->id() || particle2->id() == (particle4->particle1Idx() || particle4->particle2Idx()) || particle3->id() == (particle4->particle1Idx() || particle4->particle2Idx()))) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    mParticle3Idx = std::numeric_limits<unsigned short>::max();
    mParticle4Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
	//-----------------------------------------------------------------------------------------	

	// -- use straight lines approximation to get point of DCA of particle1-particle2 pair
	//SL16j, Vanek
	StPicoPhysicalHelix p1Helix = particle1->helix(bField); //used at line 326
  StPicoPhysicalHelix p2Helix = particle2->helix(bField); 
  StPicoPhysicalHelix p3Helix = particle3->helix(bField);

  // -- move origins of helices to the primary vertex origin - neede later in the code
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));
  p3Helix.moveOrigin(p3Helix.pathLength(vtx));
  // No need to move p4 to origin already defined as such

  TVector3 const p1Mom = particle1->gMom(); //momentum at pVtx, comuted inside StPicoTrack
  TVector3 const p2Mom = particle2->gMom(); //bFiled NOT in kilogauss - properly computed inside helix(double) function in StPicoTrack
  TVector3 const p3Mom = particle3->gMom();
  TVector3 const p4Mom(particle4->px(), particle4->py(), particle4->pz());

  // Build p4 helix from pair: assuming pair will always be neutral charge and origin is set to the given vtx 
  StPicoPhysicalHelix p4Helix(p4Mom, vtx, bField * kilogauss, 0);
  
  StPicoPhysicalHelix const p1StraightLine(p1Mom, particle1->origin(), 0, particle1->charge());
  StPicoPhysicalHelix const p2StraightLine(p2Mom, particle2->origin(), 0, particle2->charge());
  StPicoPhysicalHelix const p3StraightLine(p3Mom, particle3->origin(), 0, particle3->charge());
//-----------------------------------------------------------------------------------------
  
  pair<double, double> const ss12 = p1StraightLine.pathLengths(p2StraightLine);
  TVector3 const p1AtDcaToP2 = p1StraightLine.at(ss12.first);
  TVector3 const p2AtDcaToP1 = p2StraightLine.at(ss12.second);

  pair<double, double> const ss13 = p1StraightLine.pathLengths(p3StraightLine);
  TVector3 const p1AtDcaToP3 = p1StraightLine.at(ss13.first);
  TVector3 const p3AtDcaToP1 = p3StraightLine.at(ss13.second);

  pair<double, double> const ss14 = p1StraightLine.pathLengths(p4Helix);
  TVector3 const p1AtDcaToP4 = p1StraightLine.at(ss14.first);
  TVector3 const p4AtDcaToP1 = p4Helix.at(ss14.second);

  pair<double, double> const ss23 = p2StraightLine.pathLengths(p3StraightLine);
  TVector3 const p2AtDcaToP3 = p2StraightLine.at(ss23.first);
  TVector3 const p3AtDcaToP2 = p3StraightLine.at(ss23.second);

  pair<double, double> const ss24 = p2StraightLine.pathLengths(p4Helix);
  TVector3 const p2AtDcaToP4 = p2StraightLine.at(ss24.first);
  TVector3 const p4AtDcaToP2 = p4Helix.at(ss24.second);
  
  pair<double, double> const ss34 = p3StraightLine.pathLengths(p4Helix);
  TVector3 const p3AtDcaToP4 = p3StraightLine.at(ss34.first);
  TVector3 const p4AtDcaToP3 = p4Helix.at(ss34.second);
  
  // -- calculate DCA of particle1 to particl2 at their DCA
  mDcaDaughters12 = (p1AtDcaToP2 - p2AtDcaToP1).Mag();

  // -- calculate DCA of particle1 to particl3 at their DCA
  mDcaDaughters13 = (p1AtDcaToP3 - p3AtDcaToP1).Mag();
  
  // -- calculate DCA of particle1 to particl4 at their DCA
  mDcaDaughters14 = (p1AtDcaToP4 - p4AtDcaToP1).Mag();

  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters23 = (p2AtDcaToP3 - p3AtDcaToP2).Mag();

  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters24 = (p2AtDcaToP4 - p4AtDcaToP2).Mag();
  
  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters34 = (p3AtDcaToP4 - p4AtDcaToP3).Mag();

  
  // -- calculate decay vertex (secondary)
  TVector3 mDecayVertex = ( p1AtDcaToP2 + p2AtDcaToP1 + p1AtDcaToP3 + p3AtDcaToP1 + p1AtDcaToP4 + p4AtDcaToP1 + p2AtDcaToP3 + p3AtDcaToP2 + p2AtDcaToP4 + p4AtDcaToP2 +p3AtDcaToP4 + p4AtDcaToP3)*(1./12.0);
  
  // -- constructing mother daughter four momentum. Need helix (not straight line) for each daughter
  double const p1AtV0 = p1Helix.pathLength( mDecayVertex );
  TVector3 const p1MomAtDca = p1Helix.momentumAt(p1AtV0 ,  bField * kilogauss);

  double const p2AtV0 = p2Helix.pathLength( mDecayVertex );
  TVector3 const p2MomAtDca = p2Helix.momentumAt(p2AtV0 ,  bField * kilogauss);
  
  double const p3AtV0 = p3Helix.pathLength( mDecayVertex );
  TVector3 const p3MomAtDca = p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss);

  double const p4AtV0 = p4Helix.pathLength( mDecayVertex );
  TVector3 const p4MomAtDca = p4Helix.momentumAt(p4AtV0 ,  bField * kilogauss);
  
  TLorentzVector const p1FourMom(p1MomAtDca, sqrt(p1MomAtDca.Mag2() + p1MassHypo*p1MassHypo));
  TLorentzVector const p2FourMom(p2MomAtDca, sqrt(p2MomAtDca.Mag2() + p2MassHypo*p2MassHypo));
  TLorentzVector const p3FourMom(p3MomAtDca, sqrt(p3MomAtDca.Mag2() + p3MassHypo*p3MassHypo));
  TLorentzVector const p4FourMom = particle4->lorentzVector();
  
  mLorentzVector = p1FourMom + p2FourMom + p3FourMom + p4FourMom;
   
  // -- calculate cosThetaStar
  //    ->> Lomnitz: Need to rethink theta star
  TLorentzVector const quadrupletFourMomReverse(-mLorentzVector.Px(), -mLorentzVector.Py(), -mLorentzVector.Pz(), mLorentzVector.E());
  TLorentzVector p1FourMomStar = p1FourMom;
  p1FourMomStar.Boost(quadrupletFourMomReverse.Vect());
  mCosThetaStar = std::cos(p1FourMomStar.Vect().Angle(mLorentzVector.Vect()));
  
  // -- calculate pointing angle and decay length
  TVector3 const vtxToV0 = mDecayVertex - vtx;
  mPointingAngle = vtxToV0.Angle(mLorentzVector.Vect());
  mDecayLength = vtxToV0.Mag();
 
  // --- calculate DCA of tracks to primary vertex
  mParticle1Dca = (p1Helix.origin() - vtx).Mag();
  mParticle2Dca = (p2Helix.origin() - vtx).Mag();
  mParticle3Dca = (p3Helix.origin() - vtx).Mag();
  mParticle4Dca = (p4Helix.origin() - vtx).Mag();


}

