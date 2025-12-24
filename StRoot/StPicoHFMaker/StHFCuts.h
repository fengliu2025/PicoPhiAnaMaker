#ifndef StHFCUTS_H
#define StHFCUTS_H

/* **************************************************
 *  Cut class for HF analysis
 *  - Based on PicoCuts class 
 *
 *  Initial Authors:  
 *            Xin Dong        (xdong@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *          **Jochen Thaeder  (jmthader@lbl.gov)   
 *
 *  Contributing Authors
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Guannan Xie     (guannanxie@lbl.gov)
 *            Miroslav Simko  (msimko@bnl.gov)
 *            Jan Vanek       (jvanek@bnl.gov)
 *
 *  ** Code Maintainer
 *
 * **************************************************
 */

#include "StPicoCutsBase/StPicoCutsBase.h"


class StHFPair;
class StHFTriplet;
class StPicoTrack; //Vanek

class StHFCuts : public StPicoCutsBase
{
 public:
  
  StHFCuts();
  StHFCuts(const Char_t *name);
  ~StHFCuts();
  
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  virtual void init() { initBase(); }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   

  bool isClosePair(StHFPair const & pair) const;

  bool isGoodSecondaryVertexPair(StHFPair const & pair) const;
  bool isGoodSecondaryVertexPair_2(StHFPair const & pair) const; //to check different set of pair cuts

  bool isGoodTertiaryVertexPair(StHFPair const & pair) const;
  bool isGoodSecondaryVertexTriplet(StHFTriplet const & triplet) const;

//---------MY CUTS----------------------------------------------------------

	bool hasGoodNHitsFitMinHist(StPicoTrack const *track)	const;
	//bool hasGoodEta(TVector3 const & trkMom) const; //moved to StPicoCutsBase
	bool hasGoodNSigmaHist(StPicoTrack const *track, int hadrFlag) const;
	bool hasGoodTripletdV0Max(StHFTriplet const &triplet) const;
	bool hasGoodPtQA(StPicoTrack const *track) const;
  bool hasGoodTripletDaughtersDCAtoPV(StHFTriplet const &triplet) const;

  bool requireStrictTOF() const; //


//--------------------------------------------------------------------------

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- SETTER for CUTS
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
  void setCutSecondaryPair(float dcaDaughtersMax, float decayLengthMin, float decayLengthMax, 
  float cosThetaMin, float massMin, float massMax); 

  void setCutSecondaryPair_2(float dcaDaughtersMax, float decayLengthMin, float decayLengthMax, 
  float cosThetaMin, float massMin, float massMax); 
  
  void setCutTertiaryPair(float dcaDaughtersMax, float decayLengthMin, float decayLengthMax, 
	float cosThetaMin, float massMin, float massMax); 
  
  void setCutSecondaryTriplet(float dcaDaughters12Max, float dcaDaughters23Max, float dcaDaughters31Max, 
  float decayLengthMin, float decayLengthMax, 
	float cosThetaMin, float massMin, float massMax);

  void setCutSecondaryTripletHighPt(float dcaDaughters12Max, float dcaDaughters23Max, float dcaDaughters31Max, 
  float decayLengthMin, float decayLengthMax, 
	float cosThetaMin, float massMin, float massMax,
  float PtThreshold);

  void setCutSecondaryPairDcaToPvMax(float dcaToPvMax){ mSecondaryPairDcaToPvMax = dcaToPvMax; }

  void setCutSecondaryPairDcaToPvMax_2(float dcaToPvMax){ mSecondaryPairDcaToPvMax_2 = dcaToPvMax; }

  void setCutTertiaryPairDcaToPvMax(float dcaToPvMax) { mTertiaryPairDcaToPvMax = dcaToPvMax; }

	void setCutPtQA(float PtQA) { mPtQA = PtQA; }

//----MY SETTERS------------------------------------------------------------------------------------------

	void setCutNHitsFitMinHist(int i) {mNHitsFitMinHist = i;}
	//void setCutEta(float eta) { mEta = eta; } //moved to StPicoCutsBase
	void setCutTPCNSigmaHadronHist(float nSigHadr, int hadrFlag);
	void setCutTripletdV0Max(float dV0MaxSetCut) {mdV0MaxCut = dV0MaxSetCut;}
  void setCutTripletdV0MaxHighPt(float dV0MaxSetCut) {mdV0MaxCut_02 = dV0MaxSetCut;}
  void setCutSecondaryDaughtersDCAtoPVmin(float dcaMin) { mDcaMinDaughter = dcaMin; }
  void setCutRequireStrictTOF(bool strictTOF) { mStrictTOF = strictTOF; }

  void setStream(float stream) { mStream = stream; }
//--------------------------------------------------------------------------------------------------------

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --   
  // -- GETTER for single CUTS
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

  const float&    cutSecondaryPairDcaDaughtersMax()       const;
  const float&    cutSecondaryPairDecayLengthMin()        const;
  const float&    cutSecondaryPairDecayLengthMax()        const;
  const float&    cutSecondaryPairCosThetaMin()           const;
  const float&    cutSecondaryPairMassMin()               const;
  const float&    cutSecondaryPairMassMax()               const;
  const float&    cutSecondaryPairDcaToPvMax()		  const;

  const float&    cutSecondaryPairDcaDaughtersMax_2()       const;
  const float&    cutSecondaryPairDecayLengthMin_2()        const;
  const float&    cutSecondaryPairDecayLengthMax_2()        const;
  const float&    cutSecondaryPairCosThetaMin_2()           const;
  const float&    cutSecondaryPairMassMin_2()               const;
  const float&    cutSecondaryPairMassMax_2()               const;
  const float&    cutSecondaryPairDcaToPvMax_2()		  const;

  const float&    cutTertiaryPairDcaDaughtersMax()        const;
  const float&    cutTertiaryPairDecayLengthMin()         const;
  const float&    cutTertiaryPairDecayLengthMax()         const;
  const float&    cutTertiaryPairCosThetaMin()            const;
  const float&    cutTertiaryPairMassMin()                const;
  const float&    cutTertiaryPairMassMax()                const;
  const float&    cutTertiaryPairDcaToPvMax()		  const;

  const float&    cutSecondaryTripletDcaDaughters12Max()  const;
  const float&    cutSecondaryTripletDcaDaughters23Max()  const;
  const float&    cutSecondaryTripletDcaDaughters31Max()  const;
  const float&    cutSecondaryTripletDecayLengthMin()     const;
  const float&    cutSecondaryTripletDecayLengthMax()     const;
  const float&    cutSecondaryTripletCosThetaMin()        const;
  const float&    cutSecondaryTripletMassMin()            const;
  const float&    cutSecondaryTripletMassMax()            const;

  const float&    cutHighPtThreshold()                    const;

//----MY GETTERS------------------------------------------------------------------------------------------

	const float&    getStream()                             const;
//--------------------------------------------------------------------------------------------------------

 private:
  
  StHFCuts(StHFCuts const &);       
  StHFCuts& operator=(StHFCuts const &);

	//---MY VARIABLES------------------------
	
	int mNHitsFitMinHist;
	//float mEta; //moved to StPicoCutsBase

	float mNSigPionHist;
  float mNSigKaonHist;
	float mNSigProtonHist;

	float mdV0MaxCut;

  float mDcaMinDaughter;

	float mPtQA;

  float mStream;

  bool mStrictTOF;

  //---------------------------------------

  // ------------------------------------------
  // -- Pair cuts for secondary pair
  // ------------------------------------------
  float mSecondaryPairDcaDaughtersMax;
  float mSecondaryPairDecayLengthMin; 
  float mSecondaryPairDecayLengthMax; 
  float mSecondaryPairCosThetaMin;
  float mSecondaryPairMassMin;
  float mSecondaryPairMassMax;
  float mSecondaryPairDcaToPvMax;

  // ------------------------------------------
  // -- Pair cuts for secondary pair - alternative cuts
  // ------------------------------------------
  float mSecondaryPairDcaDaughtersMax_2;
  float mSecondaryPairDecayLengthMin_2; 
  float mSecondaryPairDecayLengthMax_2; 
  float mSecondaryPairCosThetaMin_2;
  float mSecondaryPairMassMin_2;
  float mSecondaryPairMassMax_2;
  float mSecondaryPairDcaToPvMax_2;

  // ------------------------------------------
  // -- Pair cuts tertiary pair
  // ------------------------------------------
  float mTertiaryPairDcaDaughtersMax;
  float mTertiaryPairDecayLengthMin; 
  float mTertiaryPairDecayLengthMax; 
  float mTertiaryPairCosThetaMin;
  float mTertiaryPairMassMin;
  float mTertiaryPairMassMax;
  float mTertiaryPairDcaToPvMax;

  // ------------------------------------------
  // -- Cuts of secondary triplet
  // ------------------------------------------
  float mSecondaryTripletDcaDaughters12Max;
  float mSecondaryTripletDcaDaughters23Max;
  float mSecondaryTripletDcaDaughters31Max;
  float mSecondaryTripletDecayLengthMin; 
  float mSecondaryTripletDecayLengthMax; 
  float mSecondaryTripletCosThetaMin;
  float mSecondaryTripletMassMin;
  float mSecondaryTripletMassMax;

  // ------------------------------------------
  // -- Cuts of secondary triplet (high pT) //added by Vanek
  // ------------------------------------------
  float mSecondaryTripletDcaDaughters12Max_02;
  float mSecondaryTripletDcaDaughters23Max_02;
  float mSecondaryTripletDcaDaughters31Max_02;
  float mSecondaryTripletDecayLengthMin_02; 
  float mSecondaryTripletDecayLengthMax_02; 
  float mSecondaryTripletCosThetaMin_02;
  float mSecondaryTripletMassMin_02;
  float mSecondaryTripletMassMax_02;

  float mdV0MaxCut_02;

  float mHighPtThreshold;

  ClassDef(StHFCuts,1)
};

inline void StHFCuts::setCutSecondaryPair(float dcaDaughtersMax, float decayLengthMin, float decayLengthMax, 
					  float cosThetaMin, float massMin, float massMax)  {
  mSecondaryPairDcaDaughtersMax = dcaDaughtersMax;
  mSecondaryPairDecayLengthMin = decayLengthMin; mSecondaryPairDecayLengthMax = decayLengthMax;
  mSecondaryPairCosThetaMin = cosThetaMin;
  mSecondaryPairMassMin = massMin; mSecondaryPairMassMax = massMax; 
}

inline void StHFCuts::setCutSecondaryPair_2(float dcaDaughtersMax, float decayLengthMin, float decayLengthMax, 
					  float cosThetaMin, float massMin, float massMax)  {
  mSecondaryPairDcaDaughtersMax_2 = dcaDaughtersMax;
  mSecondaryPairDecayLengthMin_2 = decayLengthMin; mSecondaryPairDecayLengthMax_2 = decayLengthMax;
  mSecondaryPairCosThetaMin_2 = cosThetaMin;
  mSecondaryPairMassMin_2 = massMin; mSecondaryPairMassMax_2 = massMax; 
}

inline void StHFCuts::setCutTertiaryPair(float dcaDaughtersMax, float decayLengthMin, float decayLengthMax, 
					 float cosThetaMin, float massMin, float massMax)  {
  mTertiaryPairDcaDaughtersMax = dcaDaughtersMax;
  mTertiaryPairDecayLengthMin = decayLengthMin; mTertiaryPairDecayLengthMax = decayLengthMax;
  mTertiaryPairCosThetaMin = cosThetaMin;
  mTertiaryPairMassMin = massMin; mTertiaryPairMassMax = massMax; 
}
  
inline void StHFCuts::setCutSecondaryTriplet(float dcaDaughters12Max, float dcaDaughters23Max, float dcaDaughters31Max, 
					     float decayLengthMin, float decayLengthMax, 
					     float cosThetaMin, float massMin, float massMax)  {
  // setting up the triplet
  mSecondaryTripletDcaDaughters12Max = dcaDaughters12Max; mSecondaryTripletDcaDaughters23Max = dcaDaughters23Max; 
  mSecondaryTripletDcaDaughters31Max = dcaDaughters31Max; 
  mSecondaryTripletDecayLengthMin = decayLengthMin; mSecondaryTripletDecayLengthMax = decayLengthMax; 
  mSecondaryTripletCosThetaMin = cosThetaMin;
  mSecondaryTripletMassMin = massMin; mSecondaryTripletMassMax = massMax;

  mHighPtThreshold = 0; //Set default value. If want to use high-pT, call StHFCuts::setCutSecondaryTripletHighPt() AFTER StHFCuts::setCutSecondaryTriplet() in you run macro.

  // setting up the doublet
  mSecondaryPairDcaDaughtersMax = mSecondaryTripletDcaDaughters12Max;
  mSecondaryPairDecayLengthMin = mSecondaryTripletDecayLengthMin;
  mSecondaryPairDcaDaughtersMax = mSecondaryTripletDecayLengthMax;
}

inline void StHFCuts::setCutSecondaryTripletHighPt(float dcaDaughters12Max, float dcaDaughters23Max, float dcaDaughters31Max, 
					     float decayLengthMin, float decayLengthMax, 
					     float cosThetaMin, float massMin, float massMax,
               float PtThreshold )  {
  // setting up the triplet
  mSecondaryTripletDcaDaughters12Max_02 = dcaDaughters12Max; mSecondaryTripletDcaDaughters23Max_02 = dcaDaughters23Max; 
  mSecondaryTripletDcaDaughters31Max_02 = dcaDaughters31Max; 
  mSecondaryTripletDecayLengthMin_02 = decayLengthMin; mSecondaryTripletDecayLengthMax_02 = decayLengthMax; 
  mSecondaryTripletCosThetaMin_02 = cosThetaMin;
  mSecondaryTripletMassMin_02 = massMin; mSecondaryTripletMassMax_02 = massMax;

  mHighPtThreshold = PtThreshold;
}

inline const float&    StHFCuts::cutSecondaryPairDcaDaughtersMax()       const { return mSecondaryPairDcaDaughtersMax; }
inline const float&    StHFCuts::cutSecondaryPairDecayLengthMin()        const { return mSecondaryPairDecayLengthMin; }
inline const float&    StHFCuts::cutSecondaryPairDecayLengthMax()        const { return mSecondaryPairDecayLengthMax; }
inline const float&    StHFCuts::cutSecondaryPairCosThetaMin()           const { return mSecondaryPairCosThetaMin; }
inline const float&    StHFCuts::cutSecondaryPairMassMin()               const { return mSecondaryPairMassMin; }
inline const float&    StHFCuts::cutSecondaryPairMassMax()               const { return mSecondaryPairMassMax; }
inline const float&    StHFCuts::cutSecondaryPairDcaToPvMax()            const { return mSecondaryPairDcaToPvMax; }

inline const float&    StHFCuts::cutSecondaryPairDcaDaughtersMax_2()       const { return mSecondaryPairDcaDaughtersMax_2; }
inline const float&    StHFCuts::cutSecondaryPairDecayLengthMin_2()        const { return mSecondaryPairDecayLengthMin_2; }
inline const float&    StHFCuts::cutSecondaryPairDecayLengthMax_2()        const { return mSecondaryPairDecayLengthMax_2; }
inline const float&    StHFCuts::cutSecondaryPairCosThetaMin_2()           const { return mSecondaryPairCosThetaMin_2; }
inline const float&    StHFCuts::cutSecondaryPairMassMin_2()               const { return mSecondaryPairMassMin_2; }
inline const float&    StHFCuts::cutSecondaryPairMassMax_2()               const { return mSecondaryPairMassMax_2; }
inline const float&    StHFCuts::cutSecondaryPairDcaToPvMax_2()            const { return mSecondaryPairDcaToPvMax_2; }

inline const float&    StHFCuts::cutTertiaryPairDcaDaughtersMax()        const { return mTertiaryPairDcaDaughtersMax; }
inline const float&    StHFCuts::cutTertiaryPairDecayLengthMin()         const { return mTertiaryPairDecayLengthMin; }
inline const float&    StHFCuts::cutTertiaryPairDecayLengthMax()         const { return mTertiaryPairDecayLengthMax; }
inline const float&    StHFCuts::cutTertiaryPairCosThetaMin()            const { return mTertiaryPairCosThetaMin; }
inline const float&    StHFCuts::cutTertiaryPairMassMin()                const { return mTertiaryPairMassMin; }
inline const float&    StHFCuts::cutTertiaryPairMassMax()                const { return mTertiaryPairMassMax; }
inline const float&    StHFCuts::cutTertiaryPairDcaToPvMax()             const { return mTertiaryPairDcaToPvMax; }

inline const float&    StHFCuts::cutSecondaryTripletDcaDaughters12Max()  const { return mSecondaryTripletDcaDaughters12Max; }
inline const float&    StHFCuts::cutSecondaryTripletDcaDaughters23Max()  const { return mSecondaryTripletDcaDaughters23Max; }
inline const float&    StHFCuts::cutSecondaryTripletDcaDaughters31Max()  const { return mSecondaryTripletDcaDaughters31Max; }
inline const float&    StHFCuts::cutSecondaryTripletDecayLengthMin()     const { return mSecondaryTripletDecayLengthMin; }
inline const float&    StHFCuts::cutSecondaryTripletDecayLengthMax()     const { return mSecondaryTripletDecayLengthMax; }
inline const float&    StHFCuts::cutSecondaryTripletCosThetaMin()        const { return mSecondaryTripletCosThetaMin; }
inline const float&    StHFCuts::cutSecondaryTripletMassMin()            const { return mSecondaryTripletMassMin; }
inline const float&    StHFCuts::cutSecondaryTripletMassMax()            const { return mSecondaryTripletMassMax; }

inline const float&    StHFCuts::cutHighPtThreshold()                    const { return mHighPtThreshold; }

inline const float&    StHFCuts::getStream()                             const { return mStream; }

#endif
