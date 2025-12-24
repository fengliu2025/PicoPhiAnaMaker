//Modified local version to make compatible with PicoDst.
//Added function compatible with ROOT versions of vectors and helixes to be used with new PicoDst format.
//Jan Vanek, 2022/08/25  

/***************************************************************************
 *
 * $Id: tofPathLength.cc,v 1.1 2024/11/04 17:59:46 vanekjan Exp $
 *
 * Author: Frank Geurts
 ***************************************************************************
 *
 * Description: Calculate helix path length between to points.
 *              requires begin and end point StThreeVectors and helix curvature
 *
 ***************************************************************************
 *
 **************************************************/
#include "StThreeVectorD.hh"
#include "StThreeVectorF.hh"
#include "TVector3.h"
#include "tofPathLength.hh"

// ----------------------------------------------------------------------------------------------
#if 0
double tofPathLength(const StThreeVectorD*  beginPoint, const StThreeVectorD* endPoint, const double curvature){
  // used in StTofGeometry::tofHelixToArray
  double x,y,z;
  x = (double)beginPoint->x();
  y = (double)beginPoint->y();
  z = (double)beginPoint->z();
  StThreeVector<double> bp (x,y,z);
  x = (double)endPoint->x();
  y = (double)endPoint->y();
  z = (double)endPoint->z();
  StThreeVector<double> ep(x,y,z);
  return tofPathLength(&bp,&ep,curvature);
}
#endif

double tofPathLength(const StThreeVectorD*  beginPoint, const StThreeVectorF* endPoint, const double curvature){
  double x,y,z;
  x = (double)beginPoint->x();
  y = (double)beginPoint->y();
  z = (double)beginPoint->z();
  StThreeVector<double> bp (x,y,z);
  x = (double)endPoint->x();
  y = (double)endPoint->y();
  z = (double)endPoint->z();
  StThreeVector<double> ep(x,y,z);
  return tofPathLength(&bp,&ep,curvature);
}


double tofPathLength(const StThreeVectorF* beginPoint, const StThreeVectorD* endPoint, const double curvature){
  double x,y,z;
  x = (double)beginPoint->x();
  y = (double)beginPoint->y();
  z = (double)beginPoint->z();
  StThreeVector<double> bp (x,y,z);
  x = (double)endPoint->x();
  y = (double)endPoint->y();
  z = (double)endPoint->z();
  StThreeVector<double> ep(x,y,z);
  return tofPathLength(&bp,&ep,curvature);
}

double tofPathLength(const StThreeVectorF* beginPoint, const StThreeVectorF* endPoint, const double curvature){
  double x,y,z;
  x = (double)beginPoint->x();
  y = (double)beginPoint->y();
  z = (double)beginPoint->z();
  StThreeVector<double> bp (x,y,z);
  x = (double)endPoint->x();
  y = (double)endPoint->y();
  z = (double)endPoint->z();
  StThreeVector<double> ep(x,y,z);
  return tofPathLength(&bp,&ep,curvature);
}

double tofPathLength(const TVector3* beginPoint, const TVector3* endPoint, const double curvature){
  double x,y,z;
  x = (double)beginPoint->x();
  y = (double)beginPoint->y();
  z = (double)beginPoint->z();
  StThreeVector<double> bp (x,y,z);
  x = (double)endPoint->x();
  y = (double)endPoint->y();
  z = (double)endPoint->z();
  StThreeVector<double> ep(x,y,z);
  return tofPathLength(&bp,&ep,curvature);
}

// ----------------------------------------------------------------------------------------------
double tofPathLength(const StThreeVector<double>*  beginPoint, const StThreeVector<double>* endPoint, const double curvature){

  double xdif =  endPoint->x() - beginPoint->x();
  double ydif =  endPoint->y() - beginPoint->y();
  
  double C = ::sqrt(xdif*xdif + ydif*ydif);
  double s_perp = C;
  if (curvature){
    double R = 1/curvature;
    s_perp = 2*R * asin(C/(2*R));
  }

  double s_z = fabs(endPoint->z() - beginPoint->z());
  double value = ::sqrt(s_perp*s_perp + s_z*s_z);

  return(value);
}

