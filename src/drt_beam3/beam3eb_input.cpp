/*!----------------------------------------------------------------------
\file beam3eb_input.cpp
\brief

\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262


*-----------------------------------------------------------------------------------------------------------*/

#include "beam3eb.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/largerotations.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3eb::ReadElement(const std::string& eletype,
                                       const std::string& distype,
                                       DRT::INPUT::LineDefinition* linedef)
{

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  linedef->ExtractDouble("CROSS",crosssec_);

  /*read beam moments of inertia of area; currently the beam3eb element works only with rotationally symmetric
   * crosssection so that the moment of inertia of area around both principal axes can be expressed by one input
   * number I_; however, the implementation itself is a general one and works also for other cases;*/

  //currently only rotationally symmetric profiles for beam --> Iyy = Izz
  linedef->ExtractDouble("MOMIN",Iyy_);
  linedef->ExtractDouble("MOMIN",Izz_);

  //torsional moment of inertia
  linedef->ExtractDouble("MOMINPOL",Irr_);

  return true;
}
/*------------------------------------------------------------------------*
 | Set moment of inertia                          (public) mukherjee 11/13|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::SetIyy(const double& Iyy)
{
  Iyy_ = Iyy;
  return;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                          (public) mukherjee 11/13|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::SetIzz(const double& Izz)
{
  Izz_ = Izz;
  return;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                          (public) mukherjee 11/13|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::SetIrr(const double& Irr)
{
  Irr_ = Irr;
  return;
}

/*------------------------------------------------------------------------*
 | Set cross section area                         (public) mukherjee 11/13|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3eb::SetCrossSec(const double& crosssec)
{
  crosssec_ = crosssec;
  return;
}


