/*!----------------------------------------------------------------------
\file beam3cl_input.cpp
\brief


\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262

 *-----------------------------------------------------------------------*/


#include "beam3cl.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/largerotations.H"

//namespace with utility functions for operations with large rotations used
using namespace LARGEROTATIONS;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::BeamCL::ReadElement(const std::string& eletype,
                                        const std::string& distype,
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

  /*read beam moments of inertia of area; currently the BeamCL element works only with rotationally symmetric
   * crosssection so that the moment of inertia of area around both principal can be expressed by one input
   * number I_; however, the implementation itself is a general one and works also for other cases; the only
   * point which has to be made sure is that the nodal triad T_ is initialized in the registration process
   * (->beamCL.cpp) in such a way that t1 is the unit vector along the beam axis and t2 and t3 are the principal
   * axes with moment of inertia of area Iyy_ and Izz_, respectively; so a modification to more general kinds of
   * cross sections can be done easily by allowing for more complex input right here and by calculating an approxipate
   * initial nodal triad in the frame of the registration; */
  linedef->ExtractDouble("IYY",Iyy_);
  linedef->ExtractDouble("IZZ",Izz_);
  linedef->ExtractDouble("IRR",Irr_);
  linedef->ExtractDoubleVector("BPOS",mybindingposition_);

  //set nodal tridas according to input file
  rQnew_.resize(4);
  rQold_.resize(4);
  rQconv_.resize(4);
  Qnew_.resize(2);
  Qold_.resize(2);
  Qconv_.resize(2);

  return true;
}


/*------------------------------------------------------------------------*
 | Set moment of inertia                            (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::SetIyy(const double& Iyy)
{
  Iyy_ = Iyy;
  return;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                            (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::SetIzz(const double& Izz)
{
  Izz_ = Izz;
  return;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                            (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::SetIrr(const double& Irr)
{
  Irr_ = Irr;
  return;
}

/*------------------------------------------------------------------------*
 | Set cross section area                           (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::SetCrossSec(const double& crosssec)
{
  crosssec_ = crosssec;
  return;
}

/*------------------------------------------------------------------------*
 | Set cross section area with shear correction     (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::SetCrossSecShear(const double& crosssecshear)
{
  crosssecshear_ = crosssecshear;
  return;
}
/*------------------------------------------------------------------------*
 | Set internodal Binding Position                  (public) mueller 11/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamCL::SetBindingPosition(const double& pos1,const double& pos2)
{
  mybindingposition_.resize(2);
  mybindingposition_[0]= pos1;
  mybindingposition_[1]= pos2;

  return;
}

