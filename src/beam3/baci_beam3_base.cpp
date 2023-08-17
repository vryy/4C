/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief base class for all beam elements

\level 2

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_beam3_base.H"

#include "baci_beaminteraction_calc_utils.H"
#include "baci_beaminteraction_periodic_boundingbox.H"
#include "baci_discretization_geometric_search_bounding_volume.H"
#include "baci_discretization_geometric_search_params.H"
#include "baci_inpar_browniandyn.H"  // enums
#include "baci_lib_globalproblem.H"
#include "baci_mat_beam_templated_material_generic.H"
#include "baci_structure_new_elements_paramsinterface.H"
#include "baci_utils_fad.H"

#include <Sacado.hpp>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Base::Beam3Base(int id, int owner)
    : DRT::Element(id, owner),
      Tref_(0),
      centerline_hermite_(true),
      filamenttype_(INPAR::BEAMINTERACTION::filtype_arbitrary),
      interface_ptr_(Teuchos::null),
      browndyn_interface_ptr_(Teuchos::null)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam3Base::Beam3Base(const DRT::ELEMENTS::Beam3Base& old)
    : DRT::Element(old),
      Tref_(old.Tref_),
      centerline_hermite_(old.centerline_hermite_),
      bspotposxi_(old.bspotposxi_),
      filamenttype_(old.filamenttype_)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);

  // bspotposxi_
  AddtoPack(data, bspotposxi_);
  // filamenttype_
  AddtoPack(data, filamenttype_);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  // bspotposxi_
  ExtractfromPack(position, data, bspotposxi_);
  // filamenttype_
  filamenttype_ = static_cast<INPAR::BEAMINTERACTION::FilamentType>(ExtractInt(position, data));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::SetBrownianDynParamsInterfacePtr()
{
  browndyn_interface_ptr_ = interface_ptr_->GetBrownianDynParamInterface();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::Beam3Base::ParamsInterfacePtr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<BROWNIANDYN::ParamsInterface> DRT::ELEMENTS::Beam3Base::BrownianDynParamsInterfacePtr()
    const
{
  return browndyn_interface_ptr_;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
std::vector<int> DRT::ELEMENTS::Beam3Base::GetAdditiveDofGIDs(
    const DRT::Discretization& discret, const DRT::Node& node) const
{
  std::vector<int> dofgids;
  std::vector<int> dofindices;

  // first collect all DoF indices
  this->PositionDofIndices(dofindices, node);
  this->TangentDofIndices(dofindices, node);
  this->Rotation1DDofIndices(dofindices, node);
  this->TangentLengthDofIndices(dofindices, node);

  // now ask for the GIDs of the DoFs with collected local indices
  dofgids.reserve(dofindices.size());
  for (unsigned int i = 0; i < dofindices.size(); ++i)
    dofgids.push_back(discret.Dof(0, &node, dofindices[i]));

  return dofgids;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<int> DRT::ELEMENTS::Beam3Base::GetRotVecDofGIDs(
    const DRT::Discretization& discret, const DRT::Node& node) const
{
  std::vector<int> dofgids;
  std::vector<int> dofindices;

  // first collect all DoF indices
  this->RotationVecDofIndices(dofindices, node);

  // now ask for the GIDs of the DoFs with collected local indices
  dofgids.reserve(dofindices.size());
  for (unsigned int i = 0; i < dofindices.size(); ++i)
    dofgids.push_back(discret.Dof(0, &node, dofindices[i]));

  return dofgids;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double DRT::ELEMENTS::Beam3Base::GetCircularCrossSectionRadiusForInteractions() const
{
  return GetBeamMaterial().GetInteractionRadius();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetRefPosAtXi(
    CORE::LINALG::Matrix<3, 1>& refpos, const double& xi) const
{
  const int numclnodes = this->NumCenterlineNodes();
  const int numnodalvalues = this->HermiteCenterlineInterpolation() ? 2 : 1;

  std::vector<double> zerovec;
  zerovec.resize(3 * numnodalvalues * numclnodes);

  this->GetPosAtXi(refpos, xi, zerovec);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::BeamMaterial& DRT::ELEMENTS::Beam3Base::GetBeamMaterial() const
{
  // get the material law
  Teuchos::RCP<MAT::Material> material_ptr = Material();

  if (material_ptr->MaterialType() != INPAR::MAT::m_beam_elast_hyper_generic)
    dserror("unknown or improper type of material law! expected beam material law!");

  return *static_cast<MAT::BeamMaterial*>(material_ptr.get());
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/

template <typename T>
MAT::BeamMaterialTemplated<T>& DRT::ELEMENTS::Beam3Base::GetTemplatedBeamMaterial() const
{
  return *Teuchos::rcp_dynamic_cast<MAT::BeamMaterialTemplated<T>>(Material(), true);
};


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void DRT::ELEMENTS::Beam3Base::GetConstitutiveMatrices(
    CORE::LINALG::Matrix<3, 3, T>& CN, CORE::LINALG::Matrix<3, 3, T>& CM) const
{
  GetTemplatedBeamMaterial<T>().GetConstitutiveMatrixOfForcesMaterialFrame(CN);
  GetTemplatedBeamMaterial<T>().GetConstitutiveMatrixOfMomentsMaterialFrame(CM);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void DRT::ELEMENTS::Beam3Base::GetTranslationalAndRotationalMassInertiaTensor(
    double& mass_inertia_translational, CORE::LINALG::Matrix<3, 3, T>& J) const
{
  GetTranslationalMassInertiaFactor(mass_inertia_translational);
  GetBeamMaterial().GetMassMomentOfInertiaTensorMaterialFrame(J);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetTranslationalMassInertiaFactor(
    double& mass_inertia_translational) const
{
  mass_inertia_translational = GetBeamMaterial().GetTranslationalMassInertiaFactor();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetDampingCoefficients(CORE::LINALG::Matrix<3, 1>& gamma) const
{
  switch (BrownianDynParamsInterface().HowBeamDampingCoefficientsAreSpecified())
  {
    case INPAR::BROWNIANDYN::cylinder_geometry_approx:
    {
      /* These are coefficients for a straight cylindrical rod taken from
       * Howard, p. 107, table 6.2. The order is as follows:
       * (0) damping of translation parallel to axis,
       * (1) damping of translation orthogonal to axis,
       * (2) damping of rotation around its own axis */

      gamma(0) = 2.0 * M_PI * BrownianDynParamsInterface().GetViscosity();
      gamma(1) = 4.0 * M_PI * BrownianDynParamsInterface().GetViscosity();
      gamma(2) = 4.0 * M_PI * BrownianDynParamsInterface().GetViscosity() *
                 GetCircularCrossSectionRadiusForInteractions() *
                 GetCircularCrossSectionRadiusForInteractions();

      // huge improvement in convergence of non-linear solver in case of artificial factor 4000
      //      gamma(2) *= 4000.0;

      break;
    }

    case INPAR::BROWNIANDYN::input_file:
    {
      gamma(0) =
          BrownianDynParamsInterface().GetBeamDampingCoefficientPrefactorsFromInputFile()[0] *
          BrownianDynParamsInterface().GetViscosity();
      gamma(1) =
          BrownianDynParamsInterface().GetBeamDampingCoefficientPrefactorsFromInputFile()[1] *
          BrownianDynParamsInterface().GetViscosity();
      gamma(2) =
          BrownianDynParamsInterface().GetBeamDampingCoefficientPrefactorsFromInputFile()[2] *
          BrownianDynParamsInterface().GetViscosity();

      break;
    }

    default:
    {
      dserror("Invalid choice of how damping coefficient values for beams are specified!");

      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <unsigned int ndim, typename T>
void DRT::ELEMENTS::Beam3Base::GetBackgroundVelocity(
    Teuchos::ParameterList& params,  //!< parameter list
    const CORE::LINALG::Matrix<ndim, 1, T>&
        evaluationpoint,                              //!< point at which background velocity and
                                                      //!< its gradient has to be computed
    CORE::LINALG::Matrix<ndim, 1, T>& velbackground,  //!< velocity of background fluid
    CORE::LINALG::Matrix<ndim, ndim, T>& velbackgroundgrad)
    const  //!< gradient of velocity of background fluid
{
  /*note: this function is not yet a general one, but always assumes a shear flow, where the
   * velocity of the background fluid is always directed in direction
   * params.get<int>("DBCDISPDIR",0) and orthogonal to z-axis. In 3D the velocity increases linearly
   * in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */

  // default values for background velocity and its gradient
  velbackground.PutScalar(0.0);
  velbackgroundgrad.PutScalar(0.0);
}

/*-----------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of          |
 | periodic boundary conditions; if two nodes within one element are se-       |
 | parated by a periodic boundary, one of them is shifted such that the final  |
 | distance in R^3 is the same as the initial distance in the periodic         |
 | space; the shift affects computation on element level within that           |
 | iteration step, only (no change in global variables performed)              |
 *-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::UnShiftNodePosition(
    std::vector<double>& disp, CORE::GEO::MESHFREE::BoundingBox const& periodic_boundingbox) const
{
  /* get number of degrees of freedom per node; note:
   * the following function assumes the same number of degrees
   * of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));

  // get number of nodes that are used for centerline interpolation
  unsigned int nnodecl = NumCenterlineNodes();

  // loop through all nodes except for the first node which remains
  // fixed as reference node
  static CORE::LINALG::Matrix<3, 1> d(true), ref(true), X(true);
  d.Clear();
  ref.Clear();
  X.Clear();
  for (unsigned int i = 1; i < nnodecl; ++i)
  {
    for (int dim = 0; dim < 3; ++dim)
    {
      d(dim) = disp[numdof * i + dim];
      ref(dim) = Nodes()[0]->X()[dim] + disp[numdof * 0 + dim];
      X(dim) = Nodes()[i]->X()[dim];
    }

    periodic_boundingbox.UnShift3D(d, ref, X);

    for (unsigned int dim = 0; dim < 3; ++dim)
    {
      disp[numdof * i + dim] = d(dim);
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetDirectionsOfShifts(std::vector<double>& disp,
    CORE::GEO::MESHFREE::BoundingBox const& periodic_boundingbox,
    std::vector<bool>& shift_in_dim) const
{
  /* get number of degrees of freedom per node; note:
   * the following function assumes the same number of degrees
   * of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));
  // get number of nodes that are used for centerline interpolation
  unsigned int nnodecl = NumCenterlineNodes();

  shift_in_dim.clear();
  shift_in_dim.resize(3);

  // loop through all nodes except for the first node which remains
  // fixed as reference node
  static CORE::LINALG::Matrix<3, 1> d(true), ref(true), X(true);
  d.Clear();
  ref.Clear();
  X.Clear();
  for (unsigned int i = 1; i < nnodecl; ++i)
  {
    for (int dim = 0; dim < 3; ++dim)
    {
      d(dim) = disp[numdof * i + dim];
      ref(dim) = Nodes()[0]->X()[dim] + disp[numdof * 0 + dim];
      X(dim) = Nodes()[i]->X()[dim];
    }

    periodic_boundingbox.CheckIfShiftBetweenPoints(d, ref, shift_in_dim, X);

    for (unsigned int dim = 0; dim < 3; ++dim)
    {
      disp[numdof * i + dim] = d(dim);
    }
  }
}

/*--------------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetPosOfBindingSpot(CORE::LINALG::Matrix<3, 1>& pos,
    std::vector<double>& disp, INPAR::BEAMINTERACTION::CrosslinkerType linkertype, int bspotlocn,
    CORE::GEO::MESHFREE::BoundingBox const& periodic_boundingbox) const
{
  const double xi = bspotposxi_.at(linkertype)[bspotlocn];
  // get position
  GetPosAtXi(pos, xi, disp);

  // check if pos at xi lies outside the periodic box, if it does, shift it back in
  periodic_boundingbox.Shift3D(pos);
}

/*--------------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3Base::GetTriadOfBindingSpot(CORE::LINALG::Matrix<3, 3>& triad,
    std::vector<double>& disp, INPAR::BEAMINTERACTION::CrosslinkerType linkertype,
    int bspotlocn) const
{
  const double xi = bspotposxi_.at(linkertype)[bspotlocn];
  // get position
  GetTriadAtXi(triad, xi, disp);
}

/*--------------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------------*/
CORE::GEOMETRICSEARCH::BoundingVolume DRT::ELEMENTS::Beam3Base::GetBoundingVolume(
    const DRT::Discretization& discret,
    const Teuchos::RCP<const Epetra_Vector>& result_data_dofbased,
    const Teuchos::RCP<const CORE::GEOMETRICSEARCH::GeometricSearchParams>& params) const
{
  // Get the centerline dof values of the beam.
  std::vector<double> element_posdofvec;
  BEAMINTERACTION::UTILS::ExtractPosDofVecValues(
      discret, this, result_data_dofbased, element_posdofvec);
  CORE::GEOMETRICSEARCH::BoundingVolume bounding_volume;

  CORE::LINALG::Matrix<3, 1, double> point;

  // TODO: replace this with convex hull from bezier curve (small student project?)
  // Add a certain number of points along the beam.
  const unsigned int n_points = 5;
  for (unsigned int i_point = 0; i_point < n_points; ++i_point)
  {
    const double xi = -1.0 + 2.0 / (n_points - 1) * i_point;
    this->GetPosAtXi(point, xi, element_posdofvec);
    bounding_volume.AddPoint(point);
  }

  // Add the radius times a safety factor.
  const double safety_factor = params->GetBeamBoundingVolumeScaling();
  const double radius = GetCircularCrossSectionRadiusForInteractions();
  bounding_volume.ExtendBoundaries(radius * safety_factor);

  return bounding_volume;
}

/*--------------------------------------------------------------------------------------------*
 | explicit template instantiations                                                           |
 *--------------------------------------------------------------------------------------------*/
template void DRT::ELEMENTS::Beam3Base::GetConstitutiveMatrices<double>(
    CORE::LINALG::Matrix<3, 3, double>& CN, CORE::LINALG::Matrix<3, 3, double>& CM) const;
template void DRT::ELEMENTS::Beam3Base::GetConstitutiveMatrices<Sacado::Fad::DFad<double>>(
    CORE::LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>& CN,
    CORE::LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>& CM) const;

template void DRT::ELEMENTS::Beam3Base::GetTranslationalAndRotationalMassInertiaTensor<double>(
    double&, CORE::LINALG::Matrix<3, 3, double>&) const;
template void
DRT::ELEMENTS::Beam3Base::GetTranslationalAndRotationalMassInertiaTensor<Sacado::Fad::DFad<double>>(
    double&, CORE::LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>&) const;

template void DRT::ELEMENTS::Beam3Base::GetBackgroundVelocity<3, double>(Teuchos::ParameterList&,
    const CORE::LINALG::Matrix<3, 1, double>&, CORE::LINALG::Matrix<3, 1, double>&,
    CORE::LINALG::Matrix<3, 3, double>&) const;
template void DRT::ELEMENTS::Beam3Base::GetBackgroundVelocity<3, Sacado::Fad::DFad<double>>(
    Teuchos::ParameterList&, const CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>>&,
    CORE::LINALG::Matrix<3, 1, Sacado::Fad::DFad<double>>&,
    CORE::LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>&) const;

template MAT::BeamMaterialTemplated<double>&
DRT::ELEMENTS::Beam3Base::GetTemplatedBeamMaterial<double>() const;
template MAT::BeamMaterialTemplated<Sacado::Fad::DFad<double>>&
DRT::ELEMENTS::Beam3Base::GetTemplatedBeamMaterial<Sacado::Fad::DFad<double>>() const;
