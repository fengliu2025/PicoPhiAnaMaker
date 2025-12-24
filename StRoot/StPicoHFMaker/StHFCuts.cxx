#include <limits>

#include "StHFCuts.h"

#include "StHFPair.h"
#include "StHFTriplet.h"

#include "StPicoEvent/StPicoTrack.h" //Vanek

ClassImp(StHFCuts)

// _________________________________________________________
StHFCuts::StHFCuts() : StPicoCutsBase("HFCutsBase"), 
  mSecondaryPairDcaDaughtersMax(std::numeric_limits<float>::max()), 
  mSecondaryPairDecayLengthMin(std::numeric_limits<float>::lowest()), mSecondaryPairDecayLengthMax(std::numeric_limits<float>::max()), 
  mSecondaryPairCosThetaMin(std::numeric_limits<float>::lowest()), 
  mSecondaryPairMassMin(std::numeric_limits<float>::lowest()), mSecondaryPairMassMax(std::numeric_limits<float>::max()), 
  mSecondaryPairDcaToPvMax(std::numeric_limits<float>::max()),
  mSecondaryPairDcaToPvMax_2(std::numeric_limits<float>::max()),

  mTertiaryPairDcaDaughtersMax(std::numeric_limits<float>::max()), 
  mTertiaryPairDecayLengthMin(std::numeric_limits<float>::lowest()), mTertiaryPairDecayLengthMax(std::numeric_limits<float>::max()), 
  mTertiaryPairCosThetaMin(std::numeric_limits<float>::lowest()), 
  mTertiaryPairMassMin(std::numeric_limits<float>::lowest()), mTertiaryPairMassMax(std::numeric_limits<float>::max()), 
  mTertiaryPairDcaToPvMax(std::numeric_limits<float>::max()),

  mSecondaryTripletDcaDaughters12Max(std::numeric_limits<float>::max()), mSecondaryTripletDcaDaughters23Max(std::numeric_limits<float>::max()), 
  mSecondaryTripletDcaDaughters31Max(std::numeric_limits<float>::max()), 
  mSecondaryTripletDecayLengthMin(std::numeric_limits<float>::lowest()), mSecondaryTripletDecayLengthMax(std::numeric_limits<float>::max()), 
  mSecondaryTripletCosThetaMin(std::numeric_limits<float>::lowest()), 
  mSecondaryTripletMassMin(std::numeric_limits<float>::lowest()), mSecondaryTripletMassMax(std::numeric_limits<float>::max()) {
  // -- default constructor
}

// _________________________________________________________
StHFCuts::StHFCuts(const Char_t *name) : StPicoCutsBase(name), 
  mSecondaryPairDcaDaughtersMax(std::numeric_limits<float>::max()), 
  mSecondaryPairDecayLengthMin(std::numeric_limits<float>::lowest()), mSecondaryPairDecayLengthMax(std::numeric_limits<float>::max()), 
  mSecondaryPairCosThetaMin(std::numeric_limits<float>::lowest()), 
  mSecondaryPairMassMin(std::numeric_limits<float>::lowest()), mSecondaryPairMassMax(std::numeric_limits<float>::max()), 
  mSecondaryPairDcaToPvMax(std::numeric_limits<float>::max()),
  mSecondaryPairDcaToPvMax_2(std::numeric_limits<float>::max()),

  mTertiaryPairDcaDaughtersMax(std::numeric_limits<float>::max()), 
  mTertiaryPairDecayLengthMin(std::numeric_limits<float>::lowest()), mTertiaryPairDecayLengthMax(std::numeric_limits<float>::max()), 
  mTertiaryPairCosThetaMin(std::numeric_limits<float>::lowest()), 
  mTertiaryPairMassMin(std::numeric_limits<float>::lowest()), mTertiaryPairMassMax(std::numeric_limits<float>::max()), 
  mTertiaryPairDcaToPvMax(std::numeric_limits<float>::max()),

  mSecondaryTripletDcaDaughters12Max(std::numeric_limits<float>::max()), mSecondaryTripletDcaDaughters23Max(std::numeric_limits<float>::max()), 
  mSecondaryTripletDcaDaughters31Max(std::numeric_limits<float>::max()), 
  mSecondaryTripletDecayLengthMin(std::numeric_limits<float>::lowest()), mSecondaryTripletDecayLengthMax(std::numeric_limits<float>::max()), 
  mSecondaryTripletCosThetaMin(std::numeric_limits<float>::lowest()), 
  mSecondaryTripletMassMin(std::numeric_limits<float>::lowest()), mSecondaryTripletMassMax(std::numeric_limits<float>::max()) {
  // -- constructor
}

// _________________________________________________________
StHFCuts::~StHFCuts() { 
  // destructor
  
}

// =======================================================================

// _________________________________________________________
bool StHFCuts::isClosePair(StHFPair const & pair) const {
  // -- check for a pair which is close in dca w/o mass constraint,
  //    using secondary vertex cuts
  return ( std::cos(pair.pointingAngle()) > mSecondaryPairCosThetaMin &&
	   pair.decayLength() > mSecondaryPairDecayLengthMin && pair.decayLength() < mSecondaryPairDecayLengthMax &&
	   pair.dcaDaughters() < mSecondaryPairDcaDaughtersMax);
}


// _________________________________________________________
bool StHFCuts::isGoodSecondaryVertexPair(StHFPair const & pair) const {
  // -- check for good secondary vertex pair

  return ( pair.m() > mSecondaryPairMassMin && pair.m() < mSecondaryPairMassMax &&
	   std::cos(pair.pointingAngle()) > mSecondaryPairCosThetaMin &&
	   pair.decayLength() > mSecondaryPairDecayLengthMin && pair.decayLength() < mSecondaryPairDecayLengthMax &&
	   pair.dcaDaughters() < mSecondaryPairDcaDaughtersMax &&
	   pair.DcaToPrimaryVertex() < mSecondaryPairDcaToPvMax);
}

// _________________________________________________________
bool StHFCuts::isGoodSecondaryVertexPair_2(StHFPair const & pair) const {
  // -- check for good secondary vertex pair

  return ( pair.m() > mSecondaryPairMassMin_2 && pair.m() < mSecondaryPairMassMax_2 &&
	   std::cos(pair.pointingAngle()) > mSecondaryPairCosThetaMin_2 &&
	   pair.decayLength() > mSecondaryPairDecayLengthMin_2 && pair.decayLength() < mSecondaryPairDecayLengthMax_2 &&
	   pair.dcaDaughters() < mSecondaryPairDcaDaughtersMax_2 &&
	   pair.DcaToPrimaryVertex() < mSecondaryPairDcaToPvMax_2);
}

// _________________________________________________________
bool StHFCuts::isGoodTertiaryVertexPair(StHFPair const & pair) const {
  // -- check for good tertiary vertex pair

  return ( pair.m() > mTertiaryPairMassMin && pair.m() < mTertiaryPairMassMax &&
	   std::cos(pair.pointingAngle()) > mTertiaryPairCosThetaMin &&
	   pair.decayLength() > mTertiaryPairDecayLengthMin && pair.decayLength() < mTertiaryPairDecayLengthMax &&
	   pair.dcaDaughters() < mTertiaryPairDcaDaughtersMax &&
	   pair.DcaToPrimaryVertex() < mTertiaryPairDcaToPvMax);
}

// _________________________________________________________
bool StHFCuts::isGoodSecondaryVertexTriplet(StHFTriplet const & triplet) const {
  // -- check for good secondary vertex triplet

  if(triplet.pt() < mHighPtThreshold || mHighPtThreshold == 0)
  {
    return ( triplet.m() > mSecondaryTripletMassMin && triplet.m() < mSecondaryTripletMassMax &&
	   std::cos(triplet.pointingAngle()) > mSecondaryTripletCosThetaMin &&
	   triplet.decayLength() > mSecondaryTripletDecayLengthMin && triplet.decayLength() < mSecondaryTripletDecayLengthMax &&
	   triplet.dcaDaughters12() < mSecondaryTripletDcaDaughters12Max &&
	   triplet.dcaDaughters23() < mSecondaryTripletDcaDaughters23Max &&
	   triplet.dcaDaughters31() < mSecondaryTripletDcaDaughters31Max);
  }
  else
  {
    return ( triplet.m() > mSecondaryTripletMassMin_02 && triplet.m() < mSecondaryTripletMassMax_02 &&
	   std::cos(triplet.pointingAngle()) > mSecondaryTripletCosThetaMin_02 &&
	   triplet.decayLength() > mSecondaryTripletDecayLengthMin_02 && triplet.decayLength() < mSecondaryTripletDecayLengthMax_02 &&
	   triplet.dcaDaughters12() < mSecondaryTripletDcaDaughters12Max_02 &&
	   triplet.dcaDaughters23() < mSecondaryTripletDcaDaughters23Max_02 &&
	   triplet.dcaDaughters31() < mSecondaryTripletDcaDaughters31Max_02);
  }  
  
}

// _________________________________________________________
bool StHFCuts::hasGoodNHitsFitMinHist(StPicoTrack const *track) const {
	return ( track->nHitsFit()>=mNHitsFitMinHist );
}

// _________________________________________________________
/* //moved to StPicoCutsBase
bool StHFCuts::hasGoodEta(TVector3 const & trkMom) const {
	return ( fabs(trkMom.PseudoRapidity())<= mEta );
}
*/
// _________________________________________________________
void StHFCuts::setCutTPCNSigmaHadronHist(float nSigHadr, int hadrFlag) {

		switch(hadrFlag)
		{
			case 1: mNSigPionHist = nSigHadr;
			case 2: mNSigKaonHist = nSigHadr;
			case 3: mNSigProtonHist = nSigHadr;
		}


}

// _________________________________________________________
bool StHFCuts::hasGoodNSigmaHist(StPicoTrack const *track, int hadrFlag) const {
	switch(hadrFlag)
		{
			case 1: return (fabs(track->nSigmaPion()) < mNSigPionHist );
			case 2: return (fabs(track->nSigmaKaon()) < mNSigKaonHist );
			case 3: return (fabs(track->nSigmaProton()) < mNSigProtonHist);
		}
	
}

// _________________________________________________________
bool StHFCuts::hasGoodTripletdV0Max(StHFTriplet const &triplet) const {
	
  if(triplet.pt() < mHighPtThreshold || mHighPtThreshold == 0)
  {
    return( triplet.dV0Max() < mdV0MaxCut );
  }
  else
  {
    return( triplet.dV0Max() < mdV0MaxCut_02 );
  } 

}

// _________________________________________________________
bool StHFCuts::hasGoodPtQA(StPicoTrack const *track) const {

	return ( track->gPt() > mPtQA  );
}

// _________________________________________________________
bool StHFCuts::hasGoodTripletDaughtersDCAtoPV(StHFTriplet const &triplet) const
{
  return ( triplet.particle1Dca() > mDcaMinDaughter && triplet.particle2Dca() > mDcaMinDaughter && triplet.particle3Dca() > mDcaMinDaughter );

}

//___________________________________________________________________________________________________________
bool StHFCuts::requireStrictTOF() const 
{ 
  return mStrictTOF; 
}

