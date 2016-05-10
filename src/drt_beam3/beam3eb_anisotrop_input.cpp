/*!----------------------------------------------------------------------
\file beam3eb_anisotrop_input.cpp

\brief three dimensional nonlinear rod based on a C1 curve


\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262


*-----------------------------------------------------------------------------------------------------------*/

#include "beam3eb_anisotrop.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/largerotations.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3ebanisotrop::ReadElement(const std::string& eletype,
                                       const std::string& distype,
                                       DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);
  linedef->ExtractDouble("CROSS",crosssec_);
  /*read beam moments of inertia of area; currently the Beam3ebanisotrop element works only with rotationally symmetric
   * crosssection so that the moment of inertia of area around both principal axes can be expressed by one input
   * number I_; however, the implementation itself is a general one and works also for other cases;*/
  //currently only rotationally symmetric profiles for beam --> Iyy = Izz
  linedef->ExtractDouble("MOMINY",Iyy_);
  linedef->ExtractDouble("MOMINZ",Izz_);
  //torsional moment of inertia
  linedef->ExtractDouble("MOMINPOL",Irr_);
  //Extract the tangents on the nodes of the element (currently only line2-elements are used, but the implementation
  //provides functionality for line n - elements)
  std::vector<double> tangents;
  linedef->ExtractDoubleVector("TANGENTS",tangents);
  Tref_.resize(2);
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<3; j++)
      Tref_[i](j)=tangents[3*i+j];

    //make sure that Tref_ is a unit vector, otherwise later calculations would be wrong (those who make use of the fact that the nodal values at the
    //boundaries can be directly taken from the displacement vector
    Tref_[i].Scale(1/Tref_[i].Norm2());
  }
  //Extract the curvature of the nodes of the element (NOTE: this is always done, even if it is not necessary for SR and CP system, as the mesh generation
  //is done without regard to the triad used - this choice was made as one only needs to generate a mesh once for all three triads
  std::vector<double> deriv_tangents;
  linedef->ExtractDoubleVector("CURVATURE",deriv_tangents);
  G2ref_.resize(2);
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<3; j++)
      G2ref_[i](j)=deriv_tangents[3*i+j];

    G2ref_[i].Scale(1/G2ref_[i].Norm2());
  }
  linedef->ExtractDouble("OPT",copt_);
  return true;
}
