/*!----------------------------------------------------------------------
\file beam3r_input.cpp
\brief

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "beam3r.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/largerotations.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3r::ReadElement(const std::string&          eletype,
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

  linedef->ExtractDouble("IYY",Iyy_);
  linedef->ExtractDouble("IZZ",Izz_);
  linedef->ExtractDouble("IRR",Irr_);

  linedef->ExtractDouble("IT",inertscaletrans_);
  linedef->ExtractDouble("IR1",inertscalerot1_);
  linedef->ExtractDouble("IR2",inertscalerot2_);

  //set nodal tridas according to input file
  Qnew_.resize(NumNode());
  Qold_.resize(NumNode());
  Qconv_.resize(NumNode());
  theta0_.resize(NumNode());

  // Attention! expression "TRIADS" in input file is misleading.
  // The 3 specified values per node define a rotational pseudovector, which
  // parameterizes the orientation of the triad at this node
  // (relative to the global reference coordinate system)
  // extract rotational pseudovectors at element nodes in reference configuration and save them as quaternions at each node, respectively
  std::vector<double> triads;
  linedef->ExtractDoubleVector("TRIADS",triads);
    for(int i=0; i<NumNode(); i++)
    {
      for(int j=0; j<3; j++)
        theta0_[i](j) = triads[3*i+j];

      LARGEROTATIONS::angletoquaternion(theta0_[i],Qnew_[i]);
    }

  Qold_  = Qnew_;
  Qconv_ = Qnew_;

  return true;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                            (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::SetIyy(const double& Iyy)
{
  Iyy_ = Iyy;
  return;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                            (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::SetIzz(const double& Izz)
{
  Izz_ = Izz;
  return;
}

/*------------------------------------------------------------------------*
 | Set moment of inertia                            (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::SetIrr(const double& Irr)
{
  Irr_ = Irr;
  return;
}

/*------------------------------------------------------------------------*
 | Set cross section area                           (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::SetCrossSec(const double& crosssec)
{
  crosssec_ = crosssec;
  return;
}

/*------------------------------------------------------------------------*
 | Set cross section area with shear correction     (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::SetCrossSecShear(const double& crosssecshear)
{
  crosssecshear_ = crosssecshear;
  return;
}
