//Modified local version to make compatible with PicoDst.
//Changed STAR vectors and helixes for ROOT versions.
//Jan Vanek, 2022/08/25

// $Id: StV0TofCorrection.h,v 1.1 2024/11/04 17:59:46 vanekjan Exp $
/*! \file StV0TofCorrection.h
    \brief header file of StV0TofCorrection class created by Patrick Huck (phuck@lbl.gov) [12/13/2010]

    A detailed class description and example analysis results can be found at 
    http://drupal.star.bnl.gov/STAR/system/files/stv0tofcorrection_v3.pdf
*/
// ---
// $Log: StV0TofCorrection.h,v $
// Revision 1.1  2024/11/04 17:59:46  vanekjan
// added Production
//
// Revision 1.2  2014/01/20 16:55:50  geurts
// Major update by Patrick Huck:
// - fixed memory leaks, clean-up
// - switch from a custom closest-points-of-approach-finder to the StHelix::pathLengths()
// - added a public helper function calcM2() to calculate mass^2 after TOF correction
//
// Revision 1.1  2011/01/24 21:01:49  geurts
// *** empty log message ***
//
// ---
#ifndef STV0TOFCORRECTION_H
#define STV0TOFCORRECTION_H

#include <stdio.h>
#include <stdarg.h>

#include <TROOT.h>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"

#include <iostream> 
#include "TVirtualStreamerInfo.h"

#include <vector>
#include <map>

//! template class defining an inline container 
/*! used to store 3d-vectors and mother tracks */ 
template < class T, class container = std::vector<T> >
class StInlineContainer : public container {
 private:
  Int_t NrArgs; //!< number of arguments
  StInlineContainer<T>* pContainer; //!< pointer to inline container
  
 public:
  StInlineContainer() { NrArgs = 0; pContainer = NULL; }
  virtual ~StInlineContainer() {}
  //! explicit constructor
  explicit StInlineContainer(const T vec) : container(1,vec) { } //!<
  //! operater () overloaded for filling container
  StInlineContainer<T> operator()(const T vec) {
    NrArgs++;
    pContainer->setNrArgs(NrArgs);
    pContainer->push_back(vec);
    return *pContainer;
  }
  //! setter for number of arguments
  void setNrArgs(const Int_t& c) { NrArgs = c; }
  //! getter for number of arguments
  Int_t getNrArgs() { return NrArgs; }
  //! set pointer to container
  void setPointer2Container(StInlineContainer<T>* p) { pContainer = p; }
};

typedef std::map< std::string, std::pair<float, float> > M2CutType;

//! A class providing tools to correct the time of flight of V0 particles
/*! detailed class description see
    http://drupal.star.bnl.gov/STAR/system/files/stv0tofcorrection_v3.pdf */
class StV0TofCorrection {
 private:
  //! inline function for time-of-flight calculation 
  inline Float_t calcTof(const TLorentzVector&, const Float_t&);
  //! inline function for input check
  inline Bool_t inputOk();
  //! inline function for m2 cuts
  inline Bool_t cutOnMass2(const Float_t&, const std::string&);

  Float_t clight; //!< speed of light in cm/ns
  Bool_t cleared, cleared2; //!< status of containers (cleared?) 
  Int_t NrDecays; //!< number of decays
  StInlineContainer<TVector3>* Vectors3D; //!< container for 3d-vectors
  StInlineContainer<TLorentzVector>* tracks; //!< container for mother tracks

  //! hardcoded default m2 limits for pi,K,p
  static M2CutType createM2CutMap() {
    M2CutType m;
    m["pi"] = std::make_pair(0.015, 0.025);
    m["K"] = std::make_pair(0.2, 0.3);
    m["p"] = std::make_pair(0.8085, 1.15);
    return m;
  }

 public:
  StV0TofCorrection() ; //!< constructor
  virtual ~StV0TofCorrection() ; //!< destructor

  //! function for initialization of 3d-vectors container
  StInlineContainer<TVector3> setVectors3D(const TVector3 vec) {
    if ( cleared ) {
      Vectors3D = new StInlineContainer<TVector3>;
    }
    else {
      const char* err_msg =
	"Make sure to call clearContainers() function"
	"every time a correction of a particle is done!"
	"Otherwise you're using the same container over and over again!";
      Error("StV0TofCorrection::setVectors3D", err_msg);
    }
    Vectors3D->setNrArgs(0);
    Vectors3D->setPointer2Container(Vectors3D);
    cleared = kFALSE;
    return (*Vectors3D)(vec);
  }

  //! function for initialization of mother tracks container
  StInlineContainer<TLorentzVector> setMotherTracks(const TLorentzVector tr) {
    if ( cleared2 ) {
      tracks  = new StInlineContainer<TLorentzVector>;
      tracks->setNrArgs(0);
      tracks->setPointer2Container(tracks);
      cleared2 = kFALSE;
      return (*tracks)(tr);
    }
    else {
      Error(
	  "StV0TofCorrection::setMotherTracks",
	  "Again: Don't forget clearContainers()!"
	  );
      return (*tracks);
    }
  }

  //! main function for beta correction
  Bool_t correctBeta(
      const StPicoPhysicalHelix&, const Float_t&, Float_t&,
      const Float_t& MomentumA = -999.,
      const std::string& particle = ""
      );
 
  //! function for finalization by destroying containers
  //  (free allocated memory and set pointers to NULL)
  void clearContainers() {
    delete Vectors3D; delete tracks;
    Vectors3D = NULL; tracks = NULL;
    cleared = kTRUE; cleared2 = kTRUE;
  }

  //! helper function to calculate m2 from momentum & beta
  inline Float_t calcM2(const Float_t& mom, const Float_t& beta) {
    Float_t f2 = 1./(beta*beta) - 1.;
    return mom*mom * f2;
  }

  //! default m2 cuts
  static M2CutType mM2CutMap;

  ClassDef(StV0TofCorrection,1)
};
#endif
