/*----------------------------------------------------------------------------*/
/*!
\file beam3_input.cpp

\brief three dimensional nonlinear corotational Reissner beam element

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "beam3.H"

#include "../drt_lib/drt_linedefinition.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_parameter.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  if (Material()->Parameter()->Name() != "MAT_BeamReissnerElastHyper" and
      Material()->Parameter()->Name() != "MAT_BeamReissnerElastHyper_ByModes")
  {
    dserror(
        "The material parameter definition '%s' is not supported by Beam3 element! "
        "Choose MAT_BeamReissnerElastHyper or MAT_BeamReissnerElastHyper_ByModes!",
        Material()->Parameter()->Name().c_str());
  }

  /*read beam moments of inertia of area; currently the beam3 element works only with rotationally
   * symmetric crosssection so that the moment of inertia of area around both principal can be
   * expressed by one input number I_; however, the implementation itself is a general one and works
   * also for other cases; the only point which has to be made sure is that the nodal triad T_ is
   * initialized in the registration process
   * (->beam3.cpp) in such a way that t1 is the unit vector along the beam axis and t2 and t3 are
   * the principal axes with moment of inertia of area Iyy_ and Izz_, respectively; so a
   * modification to more general kinds of cross sections can be done easily by allowing for more
   * complex input right here and by calculating an approxipate
   * initial nodal triad in the frame of the registration; */

  // safety check:
  LINALG::Matrix<3, 3> Cm, Cb;
  GetConstitutiveMatrices(Cm, Cb);

  if (Cb(1, 1) != Cb(2, 2))
    dserror(
        "currently the beam3 element works only with rotationally symmetric crosssection, "
        "i.e. the area moment of inertia w.r.t. first and second principal axis of inertia must "
        "have the same value! check material parameters MOMIN2/MOMIN3 in your input file !");

  return true;
}
