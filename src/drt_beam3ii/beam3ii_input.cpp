/*!----------------------------------------------------------------------
\file beam3ii_input.cpp
\brief

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "beam3ii.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/largerotations.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3ii::ReadElement(const std::string&          eletype,
                                         const std::string&          distype,
                                         DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  linedef->ExtractDouble("CROSS",crosssec_);

  double shear_correction = 0.0;
  linedef->ExtractDouble("SHEARCORR",shear_correction);
  crosssecshear_ = crosssec_ * shear_correction;

  /*read beam moments of inertia of area; currently the beam3ii element works only with rotationally symmetric
   * crosssection so that the moment of inertia of area around both principal can be expressed by one input
   * number I_; however, the implementation itself is a general one and works also for other cases; the only
   * point which has to be made sure is that the nodal triad T_ is initialized in the registration process
   * (->beam3ii.cpp) in such a way that t1 is the unit vector along the beam axis and t2 and t3 are the principal
   * axes with moment of inertia of area Iyy_ and Izz_, respectively; so a modification to more general kinds of
   * cross sections can be done easily by allowing for more complex input right here and by calculating an approxipate
   * initial nodal triad in the frame of the registration; */

  linedef->ExtractDouble("IYY",Iyy_);
  linedef->ExtractDouble("IZZ",Izz_);
  linedef->ExtractDouble("IRR",Irr_);
  
  //set nodal tridas according to input file
  Qnew_.resize(NumNode());
  Qold_.resize(NumNode());
  Qconv_.resize(NumNode());
  
  //extract triads at element nodes in reference configuration as rotation vectors and save them as quaternions at each node, respectively
  std::vector<double> triads;
  linedef->ExtractDoubleVector("TRIADS",triads);
  LINALG::Matrix<3,1> nodeangle;
    for(int i=0; i<NumNode(); i++)
    {
      for(int j=0; j<3; j++)
        nodeangle(j) = triads[3*i+j];
        
      LARGEROTATIONS::angletoquaternion(nodeangle,Qnew_[i]);
    }
    
  Qold_  = Qnew_;
  Qconv_ = Qnew_;

  return true;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                            (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::SetIyy(const double& Iyy)
{
  Iyy_ = Iyy;
  return;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                            (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::SetIzz(const double& Izz)
{
  Izz_ = Izz;
  return;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                            (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::SetIrr(const double& Irr)
{
  Irr_ = Irr;
  return;
}

/*------------------------------------------------------------------------*
 | Set cross section area                           (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::SetCrossSec(const double& crosssec)
{
  crosssec_ = crosssec;
  return;
}

/*------------------------------------------------------------------------*
 | Set cross section area with shear correction     (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3ii::SetCrossSecShear(const double& crosssecshear)
{
  crosssecshear_ = crosssecshear;
  return;
}

