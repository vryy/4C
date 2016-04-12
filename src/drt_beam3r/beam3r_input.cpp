/*!----------------------------------------------------------------------
\file beam3r_input.cpp

\maintainer Christoph Meier

\brief input related methods of 3D nonlinear Reissner beam element
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
  /* the triad field is discretized with Lagrange polynomials of order NumNode()-1;
   * the centerline is either discretized in the same way or with 3rd order Hermite polynomials;
   * we thus make a difference between nnodetriad and nnodecl;
   * assumptions: nnodecl<=nnodetriad
   * first nodes with local ID 0...nnodecl-1 are used for interpolation of centerline AND triad field*/
  const int nnodetriad=NumNode();

  if(linedef->HaveNamed("HERM2LIN2") or linedef->HaveNamed("HERM2LINE2") or
      linedef->HaveNamed("HERM2LIN3") or linedef->HaveNamed("HERM2LINE3") or
      linedef->HaveNamed("HERM2LIN4") or linedef->HaveNamed("HERM2LINE4") or
      linedef->HaveNamed("HERM2LIN5") or linedef->HaveNamed("HERM2LINE5"))
    centerline_hermite_=true;
  else
    centerline_hermite_=false;

  // read number of material model and cross-section specs
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

  if(linedef->HaveNamed("IT"))
    linedef->ExtractDouble("IT",inertscaletrans_);
  else
    inertscaletrans_=1.0;

  if(linedef->HaveNamed("IR1"))
    linedef->ExtractDouble("IR1",inertscalerot1_);
  else
    inertscalerot1_=1.0;

  if(linedef->HaveNamed("IR2"))
    linedef->ExtractDouble("IR2",inertscalerot2_);
  else
    inertscalerot2_=1.0;

  // set nodal triads according to input file
  Qnewnode_.resize(nnodetriad);
  Qoldnode_.resize(nnodetriad);
  Qconvnode_.resize(nnodetriad);
  theta0node_.resize(nnodetriad);

  /* Attention! expression "TRIADS" in input file is misleading.
   * The 3 specified values per node define a rotational pseudovector, which
   * parameterizes the orientation of the triad at this node
   * (relative to the global reference coordinate system)*/
  /* extract rotational pseudovectors at element nodes in reference configuration
   *  and save them as quaternions at each node, respectively*/
  std::vector<double> nodal_rotvecs;
  linedef->ExtractDoubleVector("TRIADS",nodal_rotvecs);
  for(int i=0; i<nnodetriad; i++)
  {
    for(int j=0; j<3; j++)
      theta0node_[i](j) = nodal_rotvecs[3*i+j];

    LARGEROTATIONS::angletoquaternion(theta0node_[i],Qnewnode_[i]);
  }

  Qoldnode_ = Qnewnode_;
  Qconvnode_ = Qnewnode_;

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
