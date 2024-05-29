/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of HDG transport element

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_scatra_ele_calc_hdg.hpp"

#include "4C_discretization_fem_general_elementtype.hpp"
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_boundary_integration.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_discretization_fem_general_utils_polynomial.hpp"
#include "4C_discretization_geometry_position_array.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_discret_hdg.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_calc.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::ScaTraEleCalcHDG(
    const int numdofpernode, const int numscal, const std::string& disname)
    : numdofpernode_(numdofpernode),
      numscal_(numscal),
      usescompletepoly_(),
      scatrapara_(DRT::ELEMENTS::ScaTraEleParameterStd::Instance(disname))

{
}


/*----------------------------------------------------------------------*
 | singleton access method                               hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname, bool create)
{
  static std::map<std::string, ScaTraEleCalcHDG<distype, probdim>*> instances;

  if (create)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcHDG<distype, probdim>(numdofpernode, numscal, disname);
  }

  else if (instances.find(disname) != instances.end())
  {
    for (typename std::map<std::string, ScaTraEleCalcHDG<distype, probdim>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
    {
      delete i->second;
      i->second = nullptr;
    }

    instances.clear();
    return nullptr;
  }

  return instances[disname];
}

/*----------------------------------------------------------------------*
 | Initialize Shapes                                     hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::InitializeShapes(
    const CORE::Elements::Element* ele, const std::string& disname)
{
  //  DRT::ELEMENTS::ScaTraHDG * hdgele =
  //  dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<CORE::Elements::Element*>(ele));
  // Check if this is an HDG element, if yes, can initialize...
  if (DRT::ELEMENTS::ScaTraHDG* hdgele =
          dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<CORE::Elements::Element*>(ele)))
  {
    usescompletepoly_ = hdgele->uses_complete_polynomial_space();

    if (shapes_ == Teuchos::null)
      shapes_ = Teuchos::rcp(new CORE::FE::ShapeValues<distype>(
          hdgele->Degree(), usescompletepoly_, 2 * hdgele->Degree()));
    else if (shapes_->degree_ != unsigned(hdgele->Degree()) ||
             shapes_->usescompletepoly_ != usescompletepoly_)
      shapes_ = Teuchos::rcp(new CORE::FE::ShapeValues<distype>(
          hdgele->Degree(), usescompletepoly_, 2 * hdgele->Degree()));

    int onfdofs = 0;
    for (unsigned int i = 0; i < nfaces_; ++i)
    {
      CORE::FE::ShapeValuesFaceParams svfparams(
          ele->Faces()[i]->Degree(), shapes_->usescompletepoly_, 2 * ele->Faces()[i]->Degree());

      shapesface_ = CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      onfdofs += shapesface_->nfdofs_;
    }
    hdgele->SetDofs(shapes_->ndofs_);
    hdgele->SetOnfDofs(onfdofs);

    // check if only one scalar is defined
    if (numscal_ > 1) FOUR_C_THROW("Not implemented for multiple scalars");

    if (local_solver_ == Teuchos::null)
      local_solver_ =
          Teuchos::rcp(new LocalSolver(ele, *shapes_, *shapesface_, usescompletepoly_, disname, 1));
  }
  else
    FOUR_C_THROW("Only works for HDG transport elements");
}

/*----------------------------------------------------------------------*
 | Evaluate                                              hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::Evaluate(CORE::Elements::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    CORE::Elements::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
    CORE::LINALG::SerialDenseMatrix&, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseVector&, CORE::LINALG::SerialDenseVector&)
{
  // check if this is an hdg element
  const DRT::ELEMENTS::ScaTraHDG* hdgele = dynamic_cast<const DRT::ELEMENTS::ScaTraHDG*>(ele);
  if (!hdgele) FOUR_C_THROW("Cannot cast element to scatra hdg element");

  InitializeShapes(ele, discretization.Name());

  shapes_->Evaluate(*ele);

  read_global_vectors(ele, discretization, la);
  get_material_params(ele);

  elevec1.putScalar(0.0);
  if (!local_solver_->scatrapara_->SemiImplicit())
  {
    elemat1.putScalar(0.0);
    local_solver_->AddDiffMat(elemat1, hdgele);
    local_solver_->AddReacMat(elemat1, hdgele);
  }
  local_solver_->ComputeResidual(
      params, elevec1, elemat1, interiorPhin_, tracenm_, tracen_, hdgele);

  return 0;
}


/*----------------------------------------------------------------------*
 * Evaluate Service                                      hoermann  09/15|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::EvaluateService(CORE::Elements::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    CORE::Elements::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // check if this is an hdg element
  DRT::ELEMENTS::ScaTraHDG* hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(ele);
  if (!hdgele) FOUR_C_THROW("cannot cast element to scatrahdg element");

  // get the action required

  SCATRA::Action act;
  if (params.get<bool>("hdg_action", false))
  {
    switch (Teuchos::getIntegralValue<DRT::HDGAction>(params, "action"))
    {
      case DRT::HDGAction::project_dirich_field:
        act = SCATRA::Action::project_dirich_field;
        break;
      default:
        FOUR_C_THROW("HDG Action type not supported");
    }
  }
  else
  {
    act = Teuchos::getIntegralValue<SCATRA::Action>(params, "action");
  }

  InitializeShapes(ele, discretization.Name());

  switch (act)
  {
    case SCATRA::Action::update_interior_variables:
    {
      shapes_->Evaluate(*ele);
      read_global_vectors(ele, discretization, la);

      return update_interior_variables(hdgele, params, elevec1_epetra);
      break;
    }

    case SCATRA::Action::interpolate_hdg_to_node:
    {
      shapes_->Evaluate(*ele);
      read_global_vectors(ele, discretization, la);
      return NodeBasedValues(ele, discretization, elevec1_epetra);
      break;
    }

    case SCATRA::Action::set_initial_field:
    {
      element_init(ele);
      prepare_material_params(ele);
      // set initial field
      return SetInitialField(ele, params, elevec1_epetra, elevec2_epetra);

      break;
    }
    case SCATRA::Action::calc_mat_initial:
    {
      if (hdgele->PadaptEle() || !hdgele->MatInit())
      {
        shapes_->Evaluate(*ele);
        element_init(ele);
        read_global_vectors(ele, discretization, la);
        prepare_material_params(ele);
        local_solver_->ComputeMatrices(ele);
        local_solver_->CondenseLocalPart(hdgele);
      }
      elemat1_epetra.putScalar(0.0);
      local_solver_->AddDiffMat(elemat1_epetra, hdgele);

      break;
    }
    case SCATRA::Action::project_material_field:
    {
      project_material_field(ele);
      break;
    }
    case SCATRA::Action::project_field:
    {
      shapes_->Evaluate(*ele);
      return ProjectField(ele, discretization, params, elevec1_epetra, elevec2_epetra, la);
      break;
    }
    case SCATRA::Action::time_update_material:
    {
      time_update_material(ele);
      break;
    }
    case SCATRA::Action::get_material_internal_state:
    {
      get_material_internal_state(ele, params, discretization);
      break;
    }
    case SCATRA::Action::set_material_internal_state:
    {
      set_material_internal_state(ele, params, discretization);
      break;
    }
    case SCATRA::Action::project_dirich_field:
    {
      if (params.isParameter("faceconsider"))
      {
        return ProjectDirichField(ele, params, discretization, la, elevec1_epetra);
      }
      break;
    }
    case SCATRA::Action::project_neumann_field:
    {
      int face = params.get<int>("face");
      int sumindex = 0;
      for (int i = 0; i < face; ++i)
      {
        CORE::FE::PolynomialSpaceParams parameter(CORE::FE::DisTypeToFaceShapeType<distype>::shape,
            ele->Faces()[i]->Degree(), shapes_->usescompletepoly_);
        int nfdofs = CORE::FE::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(parameter)->Size();
        sumindex += nfdofs;
      }
      local_solver_->ComputeNeumannBC(ele, params, face, elevec1_epetra, sumindex);
      break;
    }
    case SCATRA::Action::calc_padaptivity:
    {
      shapes_->Evaluate(*ele);
      read_global_vectors(ele, discretization, la);
      return CalcPAdaptivity(ele, discretization, params);
      break;
    }
    case SCATRA::Action::calc_error:
    {
      shapes_->Evaluate(*ele);
      read_global_vectors(ele, discretization, la);
      return CalcError(ele, params, elevec1_epetra);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of action for ScaTraHDG");
      break;
    }
  }  // end of switch(act)

  return 0;
}

/*----------------------------------------------------------------------*
 | Calculate node based values                           hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::NodeBasedValues(CORE::Elements::Element* ele,
    DRT::Discretization& discretization, CORE::LINALG::SerialDenseVector& elevec1)
{
  FOUR_C_ASSERT(elevec1.numRows() == (int)nen_ * (2 + nsd_), "Vector does not have correct size");
  CORE::LINALG::SerialDenseMatrix locations =
      CORE::FE::getEleNodeNumbering_nodes_paramspace(distype);
  CORE::LINALG::SerialDenseVector values(shapes_->ndofs_);

  DRT::ELEMENTS::ScaTraHDG* hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(ele);

  for (unsigned int i = 0; i < nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim = 0; idim < nsd_; idim++) shapes_->xsi(idim) = locations(idim, i);
    shapes_->polySpace_->Evaluate(shapes_->xsi, values);

    // compute values for concentrations (and gradients) by summing over all basis functions
    double sum = 0;
    std::vector<double> sumgrad(nsd_, 0.0);
    for (unsigned int k = 0; k < hdgele->ndofs_; ++k)
    {
      sum += values(k) * interiorPhinp_(k);
      for (unsigned int d = 0; d < nsd_; ++d)
        sumgrad[d] += values(k) * interiorPhinp_(k + (d + 1) * hdgele->ndofs_);
    }
    // store node value for concentrations and gradient in element vector
    elevec1(i) = sum;
    for (unsigned int d = 0; d < nsd_; ++d) elevec1(i + (2 + d) * nen_) = sumgrad[d];
  }

  // get trace solution values
  locations = CORE::FE::getEleNodeNumbering_nodes_paramspace(
      CORE::FE::DisTypeToFaceShapeType<distype>::shape);


  CORE::LINALG::SerialDenseVector touchcount(nen_);
  CORE::LINALG::SerialDenseVector fvalues(1);
  int sumindex = 0;
  for (unsigned int face = 0; face < nfaces_; ++face)
  {
    CORE::FE::ShapeValuesFaceParams svfparams(
        ele->Faces()[face]->Degree(), shapes_->usescompletepoly_, 2 * ele->Faces()[face]->Degree());
    shapesface_ = CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, face);

    fvalues.resize(shapesface_->nfdofs_);

    for (int i = 0; i < CORE::FE::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
    {
      // evaluate shape polynomials in node
      for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
        shapesface_->xsi(idim) = locations(idim, i);

      shapesface_->polySpace_->Evaluate(shapesface_->xsi, fvalues);

      double sum = 0.0;
      for (unsigned int k = 0; k < shapesface_->nfdofs_; ++k)
        sum += fvalues(k) * tracen_(sumindex + k);

      elevec1(nen_ + shapesface_->faceNodeOrder[face][i]) += sum;
      touchcount(shapesface_->faceNodeOrder[face][i])++;
    }
    sumindex += shapesface_->nfdofs_;
  }

  for (unsigned int i = 0; i < nen_; ++i) elevec1(nen_ + i) /= touchcount(i);

  return 0;
}  // NodeBasedValues

/*----------------------------------------------------------------------*
 * ProjectDirichField                                  berardocco 05/20 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::ProjectDirichField(
    CORE::Elements::Element* ele, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
    CORE::LINALG::SerialDenseVector& elevec1)
{
  // get actual time
  const double time = params.get<double>("time");

  Teuchos::Array<int>* func = params.getPtr<Teuchos::Array<int>>("funct");

  const int face = params.get<unsigned int>("faceconsider");
  CORE::FE::ShapeValuesFaceParams svfparams(
      ele->Faces()[face]->Degree(), shapes_->usescompletepoly_, 2 * ele->Faces()[face]->Degree());

  shapesface_ = CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);

  shapesface_->EvaluateFace(*ele, face);

  CORE::LINALG::SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
  CORE::LINALG::SerialDenseVector trVec(shapesface_->nfdofs_);

  // integration loop
  for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
  {
    // global coordinates of current Gauss point
    double coordgp[3];  // we always need three coordinates for function evaluation!
    for (int i = 0; i < 3; ++i) coordgp[i] = shapesface_->xyzreal(i, q);

    const double fac = shapesface_->jfac(q);
    // evaluate function at current Gauss point (provide always 3D coordinates!)
    const double functfac = GLOBAL::Problem::Instance()
                                ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>((*func)[0] - 1)
                                .Evaluate(coordgp, time, 0);

    // Creating the mass matrix and the RHS vector
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    {
      // Mass matrix
      for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
        mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;

      // RHS
      trVec(i) += shapesface_->shfunct(i, q) * functfac * fac;
    }

  }  // loop over integration points

  using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
  using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
  inverseMass.setMatrix(Teuchos::rcpFromRef(mass));
  inverseMass.setVectors(Teuchos::rcpFromRef(trVec), Teuchos::rcpFromRef(trVec));
  inverseMass.solve();

  for (unsigned int node = 0; node < shapesface_->nfdofs_; node++) elevec1[node] = trVec(node);

  return 0;
}  // ProjectDirichField


/*----------------------------------------------------------------------*
 * read_global_vectors                                     hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::read_global_vectors(
    CORE::Elements::Element* ele, DRT::Discretization& discretization,
    CORE::Elements::Element::LocationArray& la)
{
  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<CORE::Elements::Element*>(ele));

  // read the HDG solution vector (for traces)

  tracen_.size(hdgele->onfdofs_);
  interiorPhin_.size(shapes_->ndofs_ * (nsd_ + 1));
  interiorPhinp_.size(shapes_->ndofs_ * (nsd_ + 1));
  tracenm_.size(hdgele->onfdofs_);
  Teuchos::RCP<const Epetra_Vector> phiaf = discretization.GetState("phiaf");
  if (phiaf == Teuchos::null) FOUR_C_THROW("Cannot get state vector phiaf");
  CORE::FE::ExtractMyValues(*phiaf, tracen_, la[0].lm_);


  if (discretization.HasState("phin"))
  {
    Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
    CORE::FE::ExtractMyValues(*phin, tracenm_, la[0].lm_);
  }

  Teuchos::RCP<const Epetra_Vector> intphinp = discretization.GetState(2, "intphinp");
  if (intphinp == Teuchos::null) FOUR_C_THROW("Cannot get state vector intphinp");
  std::vector<int> localDofs = discretization.Dof(2, ele);
  CORE::FE::ExtractMyValues(*intphinp, interiorPhinp_, localDofs);

  if (discretization.HasState(2, "intphin"))
  {
    Teuchos::RCP<const Epetra_Vector> intphin = discretization.GetState(2, "intphin");
    CORE::FE::ExtractMyValues(*intphin, interiorPhin_, localDofs);
  }

  return;
}  // read_global_vectors


/*----------------------------------------------------------------------*
 * LocalSolver                                           hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::LocalSolver(
    const CORE::Elements::Element* ele, CORE::FE::ShapeValues<distype>& shapeValues,
    CORE::FE::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly,
    const std::string& disname,
    int numscal)
    :  // ndofs_ (shapeValues.ndofs_),
       //    shapes_(&shapeValues),
       //    shapesface_(&shapeValuesFace),
       //    Amat(ndofs_,ndofs_),
       //    Bmat(ndofs_,nsd_*ndofs_),
       //    Cmat(),
       //    Dmat(nsd_*ndofs_,nsd_*ndofs_),
       //    Emat(),
       //    Gmat(),
       //    Hmat(),
       //    Mmat(ndofs_,ndofs_),
       //    EmatT(),
       //    BmatMT(nsd_*ndofs_,ndofs_),
       //    Kmat(),
       //    Ivecn_(hdgele->ndofs_),
       //    Ivecnp_(hdgele->ndofs_),
       //    Imatnpderiv_(hdgele->ndofs_,hdgele->ndofs_),
       //    invAmat(ndofs_,ndofs_),
       //    invAMmat(ndofs_,ndofs_),
       //    massPart(ndofs_,shapeValues.nqpoints_),
       //    massPartW(ndofs_,shapeValues.nqpoints_),
       //    BTAMmat(ndofs_*nsd_, ndofs_),
       //    invCondmat(ndofs_*nsd_, ndofs_*nsd_),
       //    Xmat(),
      onfdofs_(0),
      scatrapara_(DRT::ELEMENTS::ScaTraEleParameterStd::Instance(disname)),
      scatraparatimint_(DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(disname), false)
//    diff_(nsd_,nsd_),                        // diffusion coefficient
//    invdiff_(nsd_,nsd_),                     // inverse diffusion coefficient
//    reacoeff_(numscal),                     // reaction coefficient
//    tau_()

{
}  // LocalSolver


/*----------------------------------------------------------------------*
 * Compute internal and face matrices
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::ComputeMatrices(
    CORE::Elements::Element* ele)
{
  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<CORE::Elements::Element*>(ele));

  // init matrices
  hdgele->Amat_.putScalar(0.0);
  hdgele->Bmat_.putScalar(0.0);
  hdgele->Cmat_.putScalar(0.0);
  hdgele->Dmat_.putScalar(0.0);
  hdgele->Emat_.putScalar(0.0);
  hdgele->Gmat_.putScalar(0.0);
  hdgele->Hmat_.putScalar(0.0);
  hdgele->Mmat_.putScalar(0.0);
  hdgele->EmatT_.putScalar(0.0);
  hdgele->BmatMT_.putScalar(0.0);
  hdgele->Kmat_.putScalar(0.0);
  hdgele->invAMmat_.putScalar(0.0);
  hdgele->BTAMmat_.putScalar(0.0);
  hdgele->invCondmat_.putScalar(0.0);
  hdgele->Xmat_.putScalar(0.0);


  bool usescompletepoly = hdgele->uses_complete_polynomial_space();

  shapes_ = Teuchos::rcp(
      new CORE::FE::ShapeValues<distype>(hdgele->Degree(), usescompletepoly, 2 * ele->Degree()));
  shapes_->Evaluate(*ele);
  compute_interior_matrices(hdgele);

  int sumindex = 0;
  for (unsigned int nface = 0; nface < nfaces_; ++nface)
  {
    CORE::FE::ShapeValuesFaceParams svfparams(ele->Faces()[nface]->Degree(),
        shapes_->usescompletepoly_, 2 * ele->Faces()[nface]->Degree());

    shapesface_ = CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);

    shapesface_->EvaluateFace(*ele, nface);

    ComputeFaceMatrices(nface, sumindex, hdgele);
    sumindex += shapesface_->nfdofs_;
  }

  // calculate AMmat = A + (1/(dt*theta))*M
  double dt = scatraparatimint_->Dt();
  double theta = scatraparatimint_->TimeFac() * (1 / dt);

  hdgele->invAMmat_ = hdgele->Mmat_;
  hdgele->invAMmat_.scale(1.0 / (dt * theta));
  hdgele->invAMmat_ += hdgele->Amat_;
  using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
  using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseAMmat;
  inverseAMmat.setMatrix(Teuchos::rcpFromRef(hdgele->invAMmat_));
  int err = inverseAMmat.invert();
  if (err != 0)
  {
    if (scatraparatimint_->IsStationary())
      FOUR_C_THROW(
          "Inversion for AMmat failed with errorcode %d. This might be due to the fact that in "
          "stationary problems Mmat_ is a zero matrix and AMat_ (if there is no convection) only "
          "has boundary integrals. Therefore, if you are using elements with internal degrees of "
          "freedom (high degree?), invAMmat_ matrix will be singular. If none of this is the case, "
          "you'll need to find the problem yourself.",
          err);
    else
      FOUR_C_THROW("Inversion for AMmat failed with errorcode %d", err);
  }
}

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices                                   hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::ComputeFaceMatrices(
    const int face, int indexstart, DRT::ELEMENTS::ScaTraHDG* hdgele)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ScaTraEleCalcHDG::ComputeFaceMatrices");

  // Compute the matrices C, E and H
  // Here, we don't consider material properties! Don't forget this during condensation

  // stabilization parameter tau is used as defined in the inputfile, i.e. it is NOT multiplied by
  // diffusion, timestep etc. Thus in the input file one has to define it and already multiply with
  // the diffusion, timstep etc. Two possibilities are to use tau=diff/dx or tau=diff/dx+dx/dt.
  double tau = 0.0;

  // set the correct stabilization parameter tau depending on the stabilization method
  switch (scatrapara_->StabType())
  {
    case (INPAR::SCATRA::stabtype_hdg_centered):
      tau = scatrapara_->TauValue();
      break;
    case (INPAR::SCATRA::stabtype_hdg_upwind):
      if (shapesface_->normal(0, 0) + shapesface_->normal(1, 0) < 0.0)
        tau = 0.0;
      else
        tau = scatrapara_->TauValue();
      break;
    case (INPAR::SCATRA::stabtype_no_stabilization):
      tau = 0.0;
      break;
    default:
      FOUR_C_THROW("Unknown definition for stabilization parameter for HDG");
      break;
  }

  // Add convection term (velocity at quadrature points on face)
  // at the moment it is set to zero
  CORE::LINALG::SerialDenseMatrix velface(nsd_, shapesface_->nqpoints_);

  // loop over number of shape functions (scalar dofs per element)
  for (unsigned int q = 0; q < hdgele->ndofs_; ++q)
  {
    if (shapesface_->shfunctI.NonzeroOnFace(q))
    {
      // loop over number of shape functions (scalar dofs per face)
      for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
      {
        // C and E
        // numerical integration: sum over quadrature points
        double tempE = 0.0;
        double tempC = 0.0;
        std::array<double, nsd_> temp_d = {0.0};
        // loop over number of quadrature points on face
        for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
        {
          double temp =
              shapesface_->jfac(i) * shapesface_->shfunct(p, i) * shapesface_->shfunctI(q, i);
          tempE += temp;
          for (unsigned int j = 0; j < nsd_; ++j)
          {
            velface(j, i) = 0.0;
            temp_d[j] += temp * shapesface_->normals(j, i);
            tempC += temp * velface(j, i) * shapesface_->normals(j, i);
          }
          for (unsigned int j = 0; j < nsd_; ++j)
          {
            hdgele->Emat_(j * hdgele->ndofs_ + q, indexstart + p) = -temp_d[j];
            hdgele->EmatT_(indexstart + p, j * hdgele->ndofs_ + q) = -temp_d[j];
          }
        }
        hdgele->Cmat_(q, indexstart + p) = tempC - tau * tempE;
        hdgele->Gmat_(indexstart + p, q) = tau * tempE;
      }
    }  // for (unsigned int q=0; q<hdgele->ndofs_; ++q)
  }    // for (unsigned int p=0; p<nfdofs_; ++p)

  // H
  // loop over number of shape functions (scalar dofs per face)
  for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
  {
    // loop over number of shape functions (scalar dofs per face)
    for (unsigned int q = 0; q < shapesface_->nfdofs_; ++q)
    {
      double tempG = 0.0;
      double tempH = 0.0;
      // loop over number of quadrature points on face
      for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
      {
        double temp =
            shapesface_->jfac(i) * shapesface_->shfunct(p, i) * shapesface_->shfunct(q, i);
        tempG += temp;
        for (unsigned int j = 0; j < nsd_; ++j)
        {
          velface(j, i) = 0.0;  // convection
          tempH += temp * velface(j, i) * shapesface_->normals(j, i);
        }
      }
      hdgele->Hmat_(indexstart + p, indexstart + q) = tempH - tau * tempG;
    }  // for (unsigned int q=0; q<nfdofs_; ++q)
  }    // for (unsigned int p=0; p<nfdofs_; ++p)


  // one term is still missing in A
  // loop over number of shape functions (scalar dofs per element)
  for (unsigned int p = 0; p < hdgele->ndofs_; ++p)
  {
    for (unsigned int q = 0; q <= p; ++q)
    {
      double tempA = 0.0;
      if (shapesface_->shfunctI.NonzeroOnFace(p) && shapesface_->shfunctI.NonzeroOnFace(q))
      {
        // loop over number of quadrature points on face
        for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
          tempA += shapesface_->jfac(i) * shapesface_->shfunctI(p, i) * shapesface_->shfunctI(q, i);
        hdgele->Amat_(p, q) += tau * tempA;
        if (p != q) hdgele->Amat_(q, p) += tau * tempA;
      }
    }
  }

  return;
}  // ComputeFaceMatrices

/*----------------------------------------------------------------------*
 * compute_interior_matrices_tet                             hoermann 01/17|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::compute_interior_matrices_tet(
    DRT::ELEMENTS::ScaTraHDG* hdgele)
{
  CORE::LINALG::SerialDenseMatrix vel(nsd_, shapes_->nqpoints_);
  CORE::LINALG::SerialDenseMatrix gradPart(hdgele->ndofs_ * nsd_, shapes_->nqpoints_);
  CORE::LINALG::SerialDenseMatrix gradPartVel(hdgele->ndofs_, shapes_->nqpoints_);

  CORE::LINALG::SerialDenseMatrix massPart(hdgele->ndofs_, shapes_->nqpoints_);
  CORE::LINALG::SerialDenseMatrix massPartW(hdgele->ndofs_, shapes_->nqpoints_);
  std::vector<CORE::LINALG::SerialDenseMatrix> DW(
      nsd_ * nsd_, CORE::LINALG::SerialDenseMatrix(hdgele->ndofs_, hdgele->ndofs_));


  // polynomial space to get the value of the shape function at the material gauss points
  CORE::FE::PolynomialSpaceParams params(distype, hdgele->Degree(), shapes_->usescompletepoly_);
  Teuchos::RCP<CORE::FE::PolynomialSpace<probdim>> polySpace =
      CORE::FE::PolynomialSpaceCache<probdim>::Instance().Create(params);

  const CORE::FE::IntPointsAndWeights<CORE::FE::dim<distype>> intpoints(
      SCATRA::DisTypeToMatGaussRule<distype>::get_gauss_rule(2 * hdgele->Degree()));

  std::vector<CORE::LINALG::SerialDenseVector> shape_gp(intpoints.IP().nquad);
  std::vector<CORE::LINALG::SerialDenseMatrix> massPartDW(
      nsd_ * nsd_, CORE::LINALG::SerialDenseMatrix(hdgele->ndofs_, intpoints.IP().nquad));


  // coordinate of gauss points
  CORE::LINALG::Matrix<probdim, 1> gp_coord(true);

  for (int q = 0; q < intpoints.IP().nquad; ++q)
  {
    shape_gp[q].size(polySpace->Size());

    // gaussian points coordinates
    for (int idim = 0; idim < CORE::FE::dim<distype>; ++idim)
      gp_coord(idim) = intpoints.IP().qxg[q][idim];
    polySpace->Evaluate(gp_coord, shape_gp[q]);
  }

  double jacdet = shapes_->xjm.Determinant();

  CORE::LINALG::SerialDenseMatrix massPartD(hdgele->ndofs_, shape_gp.size());

  // loop over quadrature points
  for (unsigned int q = 0; q < shape_gp.size(); ++q)
  {
    // loop over shape functions
    for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
    {
      massPartD(i, q) = shape_gp[q](i);

      if (hdgele->invdiff_.size() == 1)
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            massPartDW[d * nsd_ + e](i, q) =
                shape_gp[q](i) * jacdet * intpoints.IP().qwgt[q] * hdgele->invdiff_[0](d, e);
      else if (hdgele->invdiff_.size() == shape_gp.size())
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            massPartDW[d * nsd_ + e](i, q) =
                shape_gp[q](i) * jacdet * intpoints.IP().qwgt[q] * hdgele->invdiff_[q](d, e);
      else
        FOUR_C_THROW("Diffusion tensor not defined properly");
    }
  }

  // loop over quadrature points
  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    // loop over shape functions
    for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
    {
      massPart(i, q) = shapes_->shfunct(i, q);
      massPartW(i, q) = shapes_->shfunct(i, q) * shapes_->jfac(q);

      for (unsigned int d = 0; d < nsd_; ++d)
      {
        vel(d, q) = 0.0;
        gradPart(d * hdgele->ndofs_ + i, q) = shapes_->shderxy(i * nsd_ + d, q);

        gradPartVel(i, q) += shapes_->shderxy(i * nsd_ + d, q) * vel(d, q);
      }
    }
  }

  for (unsigned int d = 0; d < nsd_; ++d)
    for (unsigned int e = 0; e < nsd_; ++e)
      CORE::LINALG::multiplyNT(DW[d * nsd_ + e], massPartD, massPartDW[d * nsd_ + e]);


  // multiply matrices to perform summation over quadrature points
  if (not scatraparatimint_->IsStationary())
  {
    CORE::LINALG::multiplyNT(hdgele->Mmat_, massPart, massPartW);
  }
  CORE::LINALG::multiplyNT(0.0, hdgele->Amat_, -1.0, gradPartVel,
      massPartW);  // first part of A matrix (only if velocity field not zero)
  CORE::LINALG::multiplyNT(0.0, hdgele->Bmat_, -1.0, massPartW, gradPart);

  for (unsigned int j = 0; j < hdgele->ndofs_; ++j)
    for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
    {
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        for (unsigned int e = 0; e < nsd_; e++)
          hdgele->Dmat_(d * hdgele->ndofs_ + i, e * hdgele->ndofs_ + j) = DW[d * nsd_ + e](i, j);
        hdgele->BmatMT_(d * hdgele->ndofs_ + i, j) =
            -1.0 * hdgele->Bmat_(j, d * hdgele->ndofs_ + i);
      }
    }

  return;
}  // compute_interior_matrices_tet

/*----------------------------------------------------------------------*
 * compute_interior_matrices                                hoermann 01/17|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::compute_interior_matrices(
    DRT::ELEMENTS::ScaTraHDG* hdgele)
{
  if (distype == CORE::FE::CellType::tet4 or distype == CORE::FE::CellType::tet10)
    compute_interior_matrices_tet(hdgele);
  else
    compute_interior_matrices_all(hdgele);

  return;
}

/*----------------------------------------------------------------------*
 * compute_interior_matrices                                hoermann 09/15|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::compute_interior_matrices_all(
    DRT::ELEMENTS::ScaTraHDG* hdgele)
{
  CORE::LINALG::SerialDenseMatrix vel(nsd_, shapes_->nqpoints_);
  CORE::LINALG::SerialDenseMatrix gradPart(hdgele->ndofs_ * nsd_, shapes_->nqpoints_);
  CORE::LINALG::SerialDenseMatrix gradPartVel(hdgele->ndofs_, shapes_->nqpoints_);

  CORE::LINALG::SerialDenseMatrix massPart(hdgele->ndofs_, shapes_->nqpoints_);
  CORE::LINALG::SerialDenseMatrix massPartW(hdgele->ndofs_, shapes_->nqpoints_);
  std::vector<CORE::LINALG::SerialDenseMatrix> massPartDW(
      nsd_ * nsd_, CORE::LINALG::SerialDenseMatrix(hdgele->ndofs_, shapes_->nqpoints_));
  std::vector<CORE::LINALG::SerialDenseMatrix> DW(
      nsd_ * nsd_, CORE::LINALG::SerialDenseMatrix(hdgele->ndofs_, hdgele->ndofs_));

  // loop over quadrature points
  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    // loop over shape functions
    for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
    {
      massPart(i, q) = shapes_->shfunct(i, q);
      massPartW(i, q) = shapes_->shfunct(i, q) * shapes_->jfac(q);

      if (hdgele->invdiff_.size() == 1)
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            massPartDW[d * nsd_ + e](i, q) =
                shapes_->shfunct(i, q) * shapes_->jfac(q) * hdgele->invdiff_[0](d, e);
      else if (hdgele->invdiff_.size() == shapes_->nqpoints_)
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            massPartDW[d * nsd_ + e](i, q) =
                shapes_->shfunct(i, q) * shapes_->jfac(q) * hdgele->invdiff_[q](d, e);
      else
        FOUR_C_THROW("Diffusion tensor not defined properly");


      for (unsigned int d = 0; d < nsd_; ++d)
      {
        vel(d, q) = 0.0;
        gradPart(d * hdgele->ndofs_ + i, q) = shapes_->shderxy(i * nsd_ + d, q);

        gradPartVel(i, q) += shapes_->shderxy(i * nsd_ + d, q) * vel(d, q);
      }
    }
  }



  for (unsigned int d = 0; d < nsd_; ++d)
    for (unsigned int e = 0; e < nsd_; ++e)
      CORE::LINALG::multiplyNT(DW[d * nsd_ + e], massPart, massPartDW[d * nsd_ + e]);


  // multiply matrices to perform summation over quadrature points
  if (not scatraparatimint_->IsStationary())
  {
    CORE::LINALG::multiplyNT(hdgele->Mmat_, massPart, massPartW);
  }
  CORE::LINALG::multiplyNT(0.0, hdgele->Amat_, -1.0, gradPartVel,
      massPartW);  // first part of A matrix (only if velocity field not zero)
  CORE::LINALG::multiplyNT(0.0, hdgele->Bmat_, -1.0, massPartW, gradPart);

  for (unsigned int j = 0; j < hdgele->ndofs_; ++j)
    for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
    {
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        for (unsigned int e = 0; e < nsd_; e++)
          hdgele->Dmat_(d * hdgele->ndofs_ + i, e * hdgele->ndofs_ + j) = DW[d * nsd_ + e](i, j);
        hdgele->BmatMT_(d * hdgele->ndofs_ + i, j) =
            -1.0 * hdgele->Bmat_(j, d * hdgele->ndofs_ + i);
      }
    }

  return;
}  // compute_interior_matrices

/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::ComputeResidual(
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& elevec,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseVector& interiorPhin,
    CORE::LINALG::SerialDenseVector& tracen, CORE::LINALG::SerialDenseVector& tracenp,
    const DRT::ELEMENTS::ScaTraHDG* hdgele)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ScaTraEleCalcHDG::ComputeResidual");

  /*
   for trapezoidal rule
                                                   -1
                        +--------------------------+
   +----------------------------------------------------------------+ |                          |
   |      n                          n       n            n         |  }=:s n+1                 |
   (1/dt*theta)M + A    B  |   | M Phi   + dt (1 - theta) ( A Phi   + B Q   + C Lambda   )      | n
   n            n R     = - [ G  E^T ] |                          |   | n       n            n     |
   + dt (1 - theta) ( G Phi    + E^T Q  + H Lambda   )  - |          -B^T         D  |   |   0     +
   dt (1 - theta) ( - B^T Phi   + D Q   + E Lambda   )  |  }=:t
                        +--------------------------+
   +----------------------------------------------------------------+ -1
                                   +-------------------------+
   +---------------------------------------------+ |                         |    | n n+1   | |
   (1/dt*theta)M + A    B  |    |  - dt (1 - theta) I    -   dt theta  I      |
                      - [ G  E^T ] |                         |    | | |      -B^T            D  | |
   0                           |
                                   +-------------------------+
   +---------------------------------------------+ I..reactive term


    With

     s =   -M * U^n + dt*(1-theta) * (  A U^n + B Q^n + C L^n )   - dt*theta I^n+1 -dt*(1-theta) I^n

     t =             dt*(1-theta) * (-B^T U^n  + D Q^n + E L^n )



                              +--------------------------+^-1   +-----+
    n+1                       |  1/(dt*theta)M + A    B  |      |  s  |                           n
   n            n R     =  - [ G  E^T ]      |                          |      |     |     +
   dt(1-theta) ( G Phi    + E^T Q  + H Lambda )  = |        -B^T           D  |      |  t  |
                              +--------------------------+      +-----+

                                                                    n          n      n
                         =  - (G x + E^T y ) +   dt(1-theta) (  G  U    + E^T Q  + H L  )

               with
                                            -1    -1                                    -1
       y = ( D - (-B^T)   (1/(dt*theta)M + A)    B)     (t - (-B^T)   (1/(dt*theta)M + A)    s)

       x = (1/(dt*theta)M + A)^-1 ( s - B y)


       AM = (1/(dt*theta)M + A)
   */

  CORE::LINALG::SerialDenseVector tempinteriorphin(hdgele->ndofs_);
  for (unsigned int i = 0; i < hdgele->ndofs_; i++) tempinteriorphin(i) = interiorPhin(i);

  CORE::LINALG::SerialDenseVector tempinteriorgradphin(hdgele->ndofs_ * nsd_);
  for (unsigned int i = 0; i < hdgele->ndofs_ * nsd_; i++)
    tempinteriorgradphin(i) = interiorPhin(hdgele->ndofs_ + i);

  double dt = scatraparatimint_->Dt();
  double theta = scatraparatimint_->TimeFac() * (1 / dt);
  const double time = scatraparatimint_->Time();
  bool source = scatrapara_->IsEMD();

  CORE::LINALG::SerialDenseVector tempVec1(hdgele->ndofs_);
  CORE::LINALG::SerialDenseVector tempVec2(hdgele->ndofs_ * nsd_);

  if (theta != 1.0)
  {
    CORE::LINALG::multiply(tempVec1, hdgele->Amat_, tempinteriorphin);
    CORE::LINALG::multiply(1.0, tempVec1, 1.0, hdgele->Bmat_, tempinteriorgradphin);
    CORE::LINALG::multiply(
        1.0, tempVec1, 1.0, hdgele->Cmat_, tracen);  // = (  A U^n + B Q^n + C L^n )
    tempVec1.scale(dt * (1.0 - theta));
  }
  CORE::LINALG::multiply(1.0, tempVec1, -1.0, hdgele->Mmat_,
      tempinteriorphin);  // = s = -M * U^n + dt*(1-theta) * (  A U^n + B Q^n + C L^n )

  CORE::LINALG::SerialDenseVector tempVecI(hdgele->ndofs_);
  if (!scatrapara_->SemiImplicit())
  {
    // Reaction term
    tempVecI = hdgele->Ivecnp_;
    if (source)
    {
      ComputeSource(hdgele, tempVecI, time + dt);
    }
    tempVecI.scale(dt * theta);
    tempVec1 += tempVecI;
    if (theta != 1.0)
    {
      tempVecI = hdgele->Ivecn_;
      if (source)
      {
        ComputeSource(hdgele, tempVecI, time);
      }
      tempVecI.scale(dt * (1.0 - theta));
      tempVec1 += tempVecI;
    }
  }
  //= M * U^n + dt*(1-theta) * (  A U^n + B Q^n + C L^n )   - dt*theta I^n+1 -dt*(1-theta) I^n
  else
  {
    tempVecI = hdgele->Ivecn_;
    if (source)
    {
      ComputeSource(hdgele, tempVecI, time);
    }
    tempVecI.scale(dt);
    tempVec1 += tempVecI;
  }

  if (theta != 1.0)
  {
    CORE::LINALG::multiply(tempVec2, hdgele->BmatMT_, tempinteriorphin);
    CORE::LINALG::multiply(1.0, tempVec2, 1.0, hdgele->Dmat_, tempinteriorgradphin);
    CORE::LINALG::multiply(
        1.0, tempVec2, 1.0, hdgele->Emat_, tracen);  // = (-B^T U^n  + D Q^n + E L^n )
    tempVec2.scale(dt * (1.0 - theta));  //= t = dt*(1-theta) * (-B^T U^n  + D Q^n + E L^n )
  }

  CORE::LINALG::multiply(
      1.0, tempVec2, -1.0, hdgele->BTAMmat_, tempVec1);  // = t - (-B^T) AM^{-1} s

  CORE::LINALG::SerialDenseVector tempVec3(hdgele->ndofs_ * nsd_);
  CORE::LINALG::multiply(tempVec3, hdgele->invCondmat_,
      tempVec2);  // = y= ( D - (-B^T)   (AM)^-1    B)^-1     (t - (-B^T)   (AM^-1)    s)

  CORE::LINALG::multiply(1.0, tempVec1, -1.0, hdgele->Bmat_, tempVec3);  // = ( s - B y)
  CORE::LINALG::SerialDenseVector tempVec4(hdgele->ndofs_);
  CORE::LINALG::multiply(
      tempVec4, hdgele->invAMmat_, tempVec1);  // = x= (1/(dt*theta)M + A)^-1 ( s - B y)


  if (theta != 1.0)
  {
    CORE::LINALG::multiply(elevec, hdgele->Gmat_, tempinteriorphin);
    CORE::LINALG::multiply(1.0, elevec, 1.0, hdgele->EmatT_, tempinteriorgradphin);
    CORE::LINALG::multiply(
        1.0, elevec, 1.0, hdgele->Hmat_, tracen);  // = (  G  U^n    + E^T Q^n  + H L^n  )
    elevec.scale(dt * (1.0 - theta));  // = dt*(1-theta) * (  G  U^n    + E^T Q^n  + H L^n  )
  }

  CORE::LINALG::multiply(1.0, elevec, -1.0, hdgele->Gmat_, tempVec4);
  CORE::LINALG::multiply(1.0, elevec, -1.0, hdgele->EmatT_,
      tempVec3);  // =  - (G x + E^T y ) +   dt(1-theta) (  G  U    + E^T Q  + H L  )

  CORE::LINALG::multiply(1.0, elevec, 1.0, hdgele->Kmat_, tracenp);

  return;
}  // ComputeResidual

/*----------------------------------------------------------------------*
 * ComputeSource
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::ComputeSource(
    const CORE::Elements::Element* ele, CORE::LINALG::SerialDenseVector& elevec1, const double time)
{
  const int funcno = scatrapara_->EMDSource();

  shapes_->Evaluate(*ele);

  // CORE::LINALG::SerialDenseVector source(nsd_);
  if (nsd_ != GLOBAL::Problem::Instance()
                  ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(funcno - 1)
                  .NumberComponents())
    FOUR_C_THROW(
        "The source does not have the correct number of components.\n The correct number of "
        "components should be equal to the number of spatial dimensions.\n Fix the source "
        "function.");

  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    CORE::LINALG::Matrix<nsd_, 1> xyz;
    // add it all up
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
      {
        double source = 0;
        for (unsigned int d = 0; d < nsd_; ++d) xyz(d) = shapes_->nodexyzreal(d, j);
        for (unsigned int d = 0; d < nsd_; ++d)
          source += shapes_->shderxy(j * nsd_ + d, q) *
                    GLOBAL::Problem::Instance()
                        ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(funcno - 1)
                        .Evaluate(xyz.A(), time, d);
        elevec1(i) += shapes_->shfunct(i, q) * source * shapes_->jfac(q);
      }
  }

  return;
}

/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::CondenseLocalPart(
    DRT::ELEMENTS::ScaTraHDG* hdgele)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ScaTraEleCalcHDG::CondenseLocalPart");

  /*
   THE MATRIX

      +-------------------------------------+  +----+
      | 1/(dt*theta)M + A      B         C  |  | U  |
      |                                     |  |    |
      |         -B^T           D         E  |  | Q  |
      |                                     |  |    |
      |           G           E^T        H  |  | L  |
      +-------------------------------------+  +----+


                     (               +----------+^-1  +-----+      )
                     (               |  AM   B  |     |  C  |      )
    K = (dt*theta)   ( - [ G  E^T ]  |          |     |     |  + H )   = (dt*theta) ( - G x - E^T y
   + H )  , (               | -B^T  D  |     |  E  |      ) (               +----------+     +-----+
   )

       with
                            -1   -1                    -1
       y = ( D - (-B^T)   AM    B)     (E - (-B^T)   AM    C)

       x = AM^-1 ( C - B y)


       AM = (1/(dt*theta)M + A)

   */

  int onfdofs = hdgele->Kmat_.numRows();

  double dt = scatraparatimint_->Dt();
  double theta = scatraparatimint_->TimeFac() * (1 / dt);

  CORE::LINALG::SerialDenseMatrix tempMat1(hdgele->ndofs_ * nsd_, hdgele->ndofs_);
  CORE::LINALG::multiply(tempMat1, hdgele->BmatMT_, hdgele->invAMmat_);  // =  (-B^T) AM^{-1}

  hdgele->BTAMmat_ = tempMat1;

  CORE::LINALG::SerialDenseMatrix tempMat2(hdgele->ndofs_ * nsd_, hdgele->ndofs_ * nsd_);
  tempMat2 = hdgele->Dmat_;

  CORE::LINALG::multiply(1.0, tempMat2, -1.0, tempMat1, hdgele->Bmat_);  // = D - (-B^T) AM^{-1} B
  CORE::LINALG::SerialDenseMatrix tempMat3(hdgele->ndofs_ * nsd_, onfdofs);
  tempMat3 = hdgele->Emat_;
  CORE::LINALG::multiply(1.0, tempMat3, -1.0, tempMat1, hdgele->Cmat_);  // = E - (-B^T) AM^{-1} C

  using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
  using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseinW;
  inverseinW.setMatrix(Teuchos::rcpFromRef(tempMat2));
  int err = inverseinW.invert();
  if (err != 0)
    FOUR_C_THROW(
        "Inversion of temporary matrix for Schur complement failed with errorcode %d", err);
  // tempMat2 = (  D - H A^{-1} B )^{-1}

  hdgele->invCondmat_ = tempMat2;

  hdgele->Kmat_ = hdgele->Hmat_;  // = H

  CORE::LINALG::SerialDenseMatrix tempMat4(hdgele->ndofs_ * nsd_, onfdofs);
  CORE::LINALG::multiply(tempMat4, tempMat2, tempMat3);                        // = y
  CORE::LINALG::multiply(1.0, hdgele->Kmat_, -1.0, hdgele->EmatT_, tempMat4);  // = - E^T y + H

  CORE::LINALG::SerialDenseMatrix tempMat5(hdgele->ndofs_, onfdofs);
  tempMat5 = hdgele->Cmat_;
  CORE::LINALG::multiply(1.0, tempMat5, -1.0, hdgele->Bmat_, tempMat4);  // = C -B y

  CORE::LINALG::SerialDenseMatrix tempMat6(hdgele->ndofs_, onfdofs);
  CORE::LINALG::multiply(tempMat6, hdgele->invAMmat_, tempMat5);  // = x = AM^{-1} ( C - B y )

  // save for later use
  hdgele->Xmat_ = tempMat6;

  CORE::LINALG::multiply(
      1.0, hdgele->Kmat_, -1.0, hdgele->Gmat_, tempMat6);  // = K = H - G x - E^T y

  hdgele->Kmat_.scale(dt * theta);

  return;
}  // CondenseLocalPart

/*----------------------------------------------------------------------*
 * Add Diffusive Part to Matrix
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::AddDiffMat(
    CORE::LINALG::SerialDenseMatrix& eleMat, const DRT::ELEMENTS::ScaTraHDG* hdgele)
{
  eleMat = hdgele->Kmat_;
  eleMat.scale(-1.0);

  return;
}  // AddDiffMat


/*----------------------------------------------------------------------*
 * Add Reactive Part to Matrix
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::AddReacMat(
    CORE::LINALG::SerialDenseMatrix& eleMat, const DRT::ELEMENTS::ScaTraHDG* hdgele)
{
  double dt = scatraparatimint_->Dt();
  double theta = scatraparatimint_->TimeFac() * (1 / dt);

  // add derivative of reaction term
  CORE::LINALG::SerialDenseMatrix tempMat1(hdgele->ndofs_, hdgele->ndofs_);
  tempMat1 = hdgele->Imatnpderiv_;  // = I'

  CORE::LINALG::SerialDenseMatrix tempMat2(hdgele->ndofs_, hdgele->onfdofs_);
  CORE::LINALG::multiply(0.0, tempMat2, -1.0, tempMat1, hdgele->Xmat_);  // = I' * (-x1)

  CORE::LINALG::SerialDenseMatrix tempMat3(hdgele->ndofs_ * nsd_, hdgele->onfdofs_);
  CORE::LINALG::multiply(
      0.0, tempMat3, -1.0, hdgele->BTAMmat_, tempMat2);  // = nullptr*y1 - (-B^T) AM^{-1} I'* (-x1)
  CORE::LINALG::SerialDenseMatrix tempMat4(hdgele->ndofs_ * nsd_, hdgele->onfdofs_);
  CORE::LINALG::multiply(tempMat4, hdgele->invCondmat_,
      tempMat3);  // = y2 = ( D - (-B^T) AM^{-1} B)^-1  (nullptr*y1 - (-B^T) AM^{-1} I'* (-x1))

  CORE::LINALG::multiply(1.0, tempMat2, -1.0, hdgele->Bmat_, tempMat4);  // = I'*(-x1) -B y2

  CORE::LINALG::SerialDenseMatrix tempMat5(hdgele->ndofs_, hdgele->onfdofs_);
  CORE::LINALG::multiply(
      tempMat5, hdgele->invAMmat_, tempMat2);  // = x2 = AM^{-1} ( I'*(-x1) - B y2 )

  CORE::LINALG::multiply(1.0, eleMat, dt * theta, hdgele->EmatT_, tempMat4);  // = K - E^T y2
  CORE::LINALG::multiply(1.0, eleMat, dt * theta, hdgele->Gmat_, tempMat5);   // = K- G x2 - E^T y2

}  // AddReacMat


/*----------------------------------------------------------------------*
 |  Compute Neumann BC                                    hoermann 09/15|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::ComputeNeumannBC(
    CORE::Elements::Element* ele, Teuchos::ParameterList& params, int face,
    CORE::LINALG::SerialDenseVector& elevec, int indexstart)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ScaTraHDGEleCalc::ComputeNeumannBC");

  CORE::Conditions::Condition* condition = params.get<CORE::Conditions::Condition*>("condition");
  if (condition == nullptr) FOUR_C_THROW("Cannot access Neumann boundary condition!");

  // get actual time
  const double time = scatraparatimint_->Time();

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const auto* onoff = &condition->parameters().Get<std::vector<int>>("onoff");
  const auto* val = &condition->parameters().Get<std::vector<double>>("val");
  const auto* func = &condition->parameters().Get<std::vector<int>>("funct");


  CORE::FE::ShapeValuesFaceParams svfparams(
      ele->Faces()[face]->Degree(), shapes_->usescompletepoly_, 2 * ele->Faces()[face]->Degree());

  shapesface_ = CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);

  shapesface_->EvaluateFace(*ele, face);
  shapes_->Evaluate(*ele);

  // integration loop
  for (unsigned int iquad = 0; iquad < shapesface_->nqpoints_; ++iquad)
  {
    // factor given by spatial function
    double functfac = 1.0;

    // global coordinates of current Gauss point
    double coordgp[3];  // we always need three coordinates for function evaluation!
    for (int i = 0; i < 3; ++i) coordgp[i] = shapesface_->xyzreal(i, iquad);

    int functnum = -1;
    const double* coordgpref = coordgp;  // needed for function evaluation

    if ((*onoff)[0])  // is this dof activated?
    {
      // factor given by spatial function
      if (func) functnum = (*func)[0];

      if (functnum > 0)
      {
        // evaluate function at current Gauss point (provide always 3D coordinates!)
        functfac = GLOBAL::Problem::Instance()
                       ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(functnum - 1)
                       .Evaluate(coordgpref, time, 0);
      }
      else
        functfac = 1.;

      const double val_fac_funct_fac = (*val)[0] * shapesface_->jfac(iquad) * functfac;

      for (unsigned int node = 0; node < shapesface_->nfdofs_; node++)
        elevec[indexstart + node] += shapesface_->shfunct(node, iquad) * val_fac_funct_fac;
    }  // if ((*onoff)[dof])
  }    // loop over integration points


  return;
}  // ComputeNeumannBC

/*----------------------------------------------------------------------*
 |  prepare material parameter                            hoermann 11/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::prepare_material_params(
    CORE::Elements::Element* ele  //!< the element we are dealing with
)
{
  Teuchos::RCP<std::vector<CORE::LINALG::SerialDenseMatrix>> difftensor =
      Teuchos::rcp(new std::vector<CORE::LINALG::SerialDenseMatrix>);

  // get the material
  Teuchos::RCP<CORE::MAT::Material> material = ele->Material();

  if (material->MaterialType() == CORE::Materials::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList>& actmat =
        Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() < numscal_) FOUR_C_THROW("Not enough materials in MatList.");

    for (int k = 0; k < numscal_; ++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP<CORE::MAT::Material> singlemat = actmat->MaterialById(matid);

      for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
        prepare_materials(ele, singlemat, k, difftensor);
    }
  }
  else
    prepare_materials(ele, material, 0, difftensor);

  DRT::ELEMENTS::ScaTraHDG* hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(ele);
  for (unsigned int i = 0; i < (*difftensor).size(); ++i)
    local_solver_->prepare_material_parameter(hdgele, (*difftensor)[i]);



  return;
}  // ScaTraEleCalcHDG::get_material_params


/*----------------------------------------------------------------------*
 |  get the material parameter                            hoermann 09/15|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::get_material_params(
    CORE::Elements::Element* ele  //!< the element we are dealing with
)
{
  CORE::LINALG::SerialDenseMatrix difftensor(nsd_, nsd_);
  CORE::LINALG::SerialDenseVector ivecn(shapes_->ndofs_);
  CORE::LINALG::SerialDenseVector ivecnp(shapes_->ndofs_);
  CORE::LINALG::SerialDenseMatrix ivecnpderiv(shapes_->ndofs_, shapes_->ndofs_);

  // get the material
  Teuchos::RCP<CORE::MAT::Material> material = ele->Material();

  if (material->MaterialType() == CORE::Materials::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList>& actmat =
        Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() < numscal_) FOUR_C_THROW("Not enough materials in MatList.");

    for (int k = 0; k < numscal_; ++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP<CORE::MAT::Material> singlemat = actmat->MaterialById(matid);

      materials(singlemat, k, difftensor, ivecn, ivecnp, ivecnpderiv);
    }
  }
  else
    materials(material, 0, difftensor, ivecn, ivecnp, ivecnpderiv);

  DRT::ELEMENTS::ScaTraHDG* hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(ele);
  local_solver_->set_material_parameter(hdgele, ivecn, ivecnp, ivecnpderiv);



  return;
}  // ScaTraEleCalcHDG::get_material_params


/*----------------------------------------------------------------------*
 * update_interior_variables
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::update_interior_variables(
    DRT::ELEMENTS::ScaTraHDG* hdgele, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector& elevec
    //    double dt
)
{
  /*
   THE MATRIX

      +-------------------------------------+  +----+
      | 1/(dt*theta)M + A      B         C  |  | U  |
      |                                     |  |    |
      |         -B^T           D         E  |  | Q  |
      |                                     |  |    |
      |           G           E^T        H  |  | L  |
      +-------------------------------------+  +----+


    +---------+                       +--------------------------+^-1
  +--------------------------------------------------------------------------------------------+ |
  U^n+1  |                       |  1/(dt*theta)M + A    B  |     |   M * U^n - dt*(1-theta) * (  A
  U^n + B Q^n + C L^n )   - dt*theta I^n+1 -dt*(1-theta) I^n | |         |    = 1/(dt * theta)   |
  |     | | - |  Q^n+1  |                       |        -B^T           D  |     |  - dt*(1-theta) *
  (   -B^T U^n  + D Q^n + E L^n )                                          |
    +---------+                       +--------------------------+
  +--------------------------------------------------------------------------------------------+

                              +--------------------------+^-1    +------+
                              |  1/(dt*theta)M + A    B  |       |   C  |
                        -     |                          |       |      |  L^n+1
                              |        -B^T           D  |       |   E  |
                              +--------------------------+       +------+


  s =   M * U^n - dt*(1-theta) * (  A U^n + B Q^n + C L^n )   - dt*theta I^n+1 -dt*(1-theta) I^n -
  dt* theta C L^n+1

  t =            - dt*(1-theta) * (-B^T U^n  + D Q^n + E L^n ) -     dt*theta E L^n+1



   +---------+                        +--------------------------+^-1   +-----+ +-----+ |  U^n+1  |
  |  1/(dt*theta)M + A    B  |      |  s  |                          |  x  | |         |   =   1/(dt
  * theta)   |                          |      |     |      =   1/(dt * theta)  |     |    , | Q^n+1
  |                        |        -B^T           D  |      |  t  |                          |  y |
   +---------+                        +--------------------------+      +-----+ +-----+

       with
                                             -1    -1                                    -1
       y = ( D - (-B^T)   (1/(dt*theta)M + A)    B)     (t - (-B^T)   (1/(dt*theta)M + A)    s)

       x = (1/(dt*theta)M + A)^-1 ( s - B y)

       AM = (1/(dt*theta)M + A)

   */

  CORE::LINALG::SerialDenseVector tempinteriorphin(hdgele->ndofs_);
  for (unsigned int i = 0; i < hdgele->ndofs_; ++i) tempinteriorphin(i) = interiorPhin_(i);

  CORE::LINALG::SerialDenseVector tempinteriorgradphin(hdgele->ndofs_ * nsd_);
  for (unsigned int i = 0; i < hdgele->ndofs_ * nsd_; ++i)
    tempinteriorgradphin(i) = interiorPhin_(hdgele->ndofs_ + i);

  double dt = local_solver_->scatraparatimint_->Dt();
  double theta = local_solver_->scatraparatimint_->TimeFac() * (1 / dt);
  const double time = local_solver_->scatraparatimint_->Time();
  bool source = local_solver_->scatrapara_->IsEMD();

  CORE::LINALG::SerialDenseVector tempVec1(hdgele->ndofs_);
  if (theta != 1.0)
  {
    CORE::LINALG::multiply(tempVec1, hdgele->Amat_, tempinteriorphin);
    CORE::LINALG::multiply(1.0, tempVec1, 1.0, hdgele->Bmat_, tempinteriorgradphin);
    CORE::LINALG::multiply(
        1.0, tempVec1, 1.0, hdgele->Cmat_, tracenm_);  // = ( A U^n + B Q^n + C L^n )
    tempVec1.scale(-dt * (1. - theta));
  }
  CORE::LINALG::multiply(1.0, tempVec1, 1.0, hdgele->Mmat_,
      tempinteriorphin);  // = -M * U^n + dt*(1-theta) * (  A U^n + B Q^n + C L^n )


  // Reaction term
  CORE::LINALG::SerialDenseVector tempVecI(hdgele->ndofs_);
  if (!local_solver_->scatrapara_->SemiImplicit())
  {
    tempVecI = hdgele->Ivecnp_;
    if (source)
    {
      local_solver_->ComputeSource(hdgele, tempVecI, time + dt);
    }
    tempVecI.scale(-dt * theta);
    tempVec1 += tempVecI;
    if (theta != 1.0)
    {
      tempVecI = hdgele->Ivecn_;
      if (source)
      {
        local_solver_->ComputeSource(hdgele, tempVecI, time);
      }
      tempVecI.scale(-dt * (1.0 - theta));
      tempVec1 += tempVecI;
    }
  }
  else
  {
    tempVecI = hdgele->Ivecn_;
    if (source)
    {
      local_solver_->ComputeSource(hdgele, tempVecI, time);
    }
    tempVecI.scale(-dt);
    tempVec1 += tempVecI;
  }


  CORE::LINALG::multiply(1.0, tempVec1, -dt * theta, hdgele->Cmat_,
      tracen_);  //= s = -M * U^n + dt*(1-theta) * (A U^n + B Q^n + C L^n) - dt*theta I^n+1
                 //-dt*(1-theta) I^n - dt* theta C L^n+1

  CORE::LINALG::SerialDenseVector tempVec2(hdgele->ndofs_ * nsd_);
  if (theta != 1.0)
  {
    CORE::LINALG::multiply(tempVec2, hdgele->BmatMT_, tempinteriorphin);
    CORE::LINALG::multiply(1.0, tempVec2, 1.0, hdgele->Dmat_, tempinteriorgradphin);
    CORE::LINALG::multiply(
        1.0, tempVec2, 1.0, hdgele->Emat_, tracenm_);  // = (-B^T U^n  + D Q^n + E L^n )
    tempVec2.scale(-dt * (1. - theta));
  }
  CORE::LINALG::multiply(1.0, tempVec2, -dt * theta, hdgele->Emat_,
      tracen_);  //= t = -dt*(1-theta) * (-B^T U^n  + D Q^n + E L^n )  - dt*theta E L^n+1

  // y = ( D - (-B^T)   (AM)^-1    B)^-1     (t - (-B^T)   (AM^-1)    s)
  // x = (1/(dt*theta)M + A)^-1 ( s - B y)


  CORE::LINALG::multiply(
      1.0, tempVec2, -1.0, hdgele->BTAMmat_, tempVec1);  // = t - (-B^T) AM^{-1} s

  CORE::LINALG::SerialDenseVector tempVec3(hdgele->ndofs_ * nsd_);
  CORE::LINALG::multiply(tempVec3, hdgele->invCondmat_,
      tempVec2);  // = y= ( D - (-B^T)   (AM)^-1    B)^-1     (t - (-B^T)   (AM^-1)    s)

  CORE::LINALG::multiply(1.0, tempVec1, -1.0, hdgele->Bmat_, tempVec3);  // = ( s - B y)
  CORE::LINALG::SerialDenseVector tempVec4(hdgele->ndofs_);
  CORE::LINALG::multiply(
      tempVec4, hdgele->invAMmat_, tempVec1);  // = x= (1/(dt*theta)M + A)^-1 ( s - B y)

  tempVec3.scale(1. / (dt * theta));
  tempVec4.scale(1. / (dt * theta));

  for (unsigned int i = 0; i < hdgele->ndofs_; ++i) elevec(i) = tempVec4(i);
  for (unsigned int i = 0; i < nsd_ * hdgele->ndofs_; ++i) elevec(hdgele->ndofs_ + i) = tempVec3(i);


  return 0;
}  // update_interior_variables


/*----------------------------------------------------------------------*
 * SetInitialField
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::SetInitialField(
    const CORE::Elements::Element* ele, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2)
{
  shapes_->Evaluate(*ele);

  CORE::LINALG::SerialDenseMatrix Mmat(shapes_->ndofs_, shapes_->ndofs_);
  CORE::LINALG::SerialDenseMatrix massPart(shapes_->ndofs_, shapes_->nqpoints_);
  CORE::LINALG::SerialDenseMatrix massPartW(shapes_->ndofs_, shapes_->nqpoints_);

  // reshape elevec2 as matrix
  FOUR_C_ASSERT(
      elevec2.numRows() == 0 || unsigned(elevec2.numRows()) == shapes_->ndofs_ * (nsd_ + 1),
      "Wrong size in project vector 2");

  // get function
  const int* start_func = params.getPtr<int>("funct");

  // internal variables
  if (elevec2.numRows() > 0)
  {
    CORE::LINALG::SerialDenseMatrix localMat(
        Teuchos::View, elevec2.values(), shapes_->ndofs_, shapes_->ndofs_, nsd_ + 1);
    localMat.putScalar(0.0);

    // create mass matrix for interior by looping over quadrature points
    for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
    {
      const double fac = shapes_->jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_->xyzreal(d, q);  // coordinates of quadrature point in real coordinates
      double phi;
      double gradphi[nsd_];

      FOUR_C_ASSERT(start_func != nullptr, "funct not set for initial value");
      if (GLOBAL::Problem::Instance()
                  ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(*start_func - 1)
                  .NumberComponents() != 1 &&
          GLOBAL::Problem::Instance()
                  ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(*start_func - 1)
                  .NumberComponents() != nsd_ + 2)
        FOUR_C_THROW(
            "Impossible to initialize the field with the given number of components of the initial "
            "field. Set the number of components to either 1 or nsd_ + 2.\nThe fields are ordered "
            "as:\n- phi\n- gradphi\n- tracephi");

      phi = GLOBAL::Problem::Instance()
                ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(*start_func - 1)
                .Evaluate(xyz, 0, 0);
      for (unsigned int i = 0; i < nsd_; ++i)
        gradphi[i] = GLOBAL::Problem::Instance()
                         ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(*start_func - 1)
                         .Evaluate(xyz, 0, 1 + i);

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      {
        massPart(i, q) = shapes_->shfunct(i, q);
        massPartW(i, q) = shapes_->shfunct(i, q) * fac;
        localMat(i, 0) += shapes_->shfunct(i, q) * phi * fac;
        for (unsigned int j = 0; j < nsd_; ++j)
          localMat(i, 1 + j) += shapes_->shfunct(i, q) * gradphi[j] * fac;
      }
    }

    CORE::LINALG::multiplyNT(Mmat, massPart, massPartW);
    {
      using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
      using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
      Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
      inverseMass.setMatrix(Teuchos::rcpFromRef(Mmat));
      inverseMass.setVectors(Teuchos::rcpFromRef(localMat), Teuchos::rcpFromRef(localMat));
      inverseMass.factorWithEquilibration(true);
      int err2 = inverseMass.factor();
      int err = inverseMass.solve();
      if (err != 0 || err2 != 0) FOUR_C_THROW("Inversion of matrix failed with errorcode %d", err);
    }
  }

  // We have to set the initial trace field and gradient field,
  // because the calculation of the residual in the first time step should have a correct trace
  // field and also a correct gradient!

  // trace variable
  int nfdofs = 0;
  for (unsigned int face = 0; face < nfaces_; ++face)
  {
    CORE::FE::ShapeValuesFaceParams svfparams(
        ele->Faces()[face]->Degree(), shapes_->usescompletepoly_, 2 * ele->Faces()[face]->Degree());
    shapesface_ = CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, face);

    CORE::LINALG::SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
    CORE::LINALG::SerialDenseMatrix trVec(shapesface_->nfdofs_, 1);

    // loop over quadrature points
    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      const double fac = shapesface_->jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d) xyz[d] = shapesface_->xyzreal(d, q);

      double trphi;
      trphi = GLOBAL::Problem::Instance()
                  ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(*start_func - 1)
                  .Evaluate(xyz, 0, nsd_ + 1);

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
          mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;
        trVec(i, 0) += shapesface_->shfunct(i, q) * trphi * fac;
      }
    }

    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
    inverseMass.setMatrix(Teuchos::rcpFromRef(mass));
    inverseMass.setVectors(Teuchos::rcpFromRef(trVec), Teuchos::rcpFromRef(trVec));
    inverseMass.factorWithEquilibration(true);
    int err2 = inverseMass.factor();
    int err = inverseMass.solve();
    if (err != 0 || err2 != 0) FOUR_C_THROW("Inversion of matrix failed with errorcode %d", err);
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i) elevec1(nfdofs + i) = trVec(i, 0);

    nfdofs += shapesface_->nfdofs_;
  }

  return 0;

}  // SetInitialField


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::prepare_materials(
    CORE::Elements::Element* ele,                            //!< the element we are dealing with
    const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
    const int k,                                             //!< id of current scalar
    Teuchos::RCP<std::vector<CORE::LINALG::SerialDenseMatrix>> difftensor  //!< diffusion tensor
)
{
  const Teuchos::RCP<const MAT::ScatraMat>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

  double diffscalar;

  // get constant diffusivity
  diffscalar = actmat->Diffusivity();

  CORE::LINALG::SerialDenseMatrix difftensortmp(nsd_, nsd_);

  for (unsigned int i = 0; i < nsd_; ++i) difftensortmp(i, i) = diffscalar;

  (*difftensor).push_back(difftensortmp);

  return;
}  // ScaTraEleCalcHDG::Materials


/*----------------------------------------------------------------------*
 |  Set material parameter                               hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::set_material_parameter(
    DRT::ELEMENTS::ScaTraHDG* hdgele,             //!< hdg element
    CORE::LINALG::SerialDenseVector& ivecn,       //!< reaction term at time n
    CORE::LINALG::SerialDenseVector& ivecnp,      //!< reaction term at time n+1
    CORE::LINALG::SerialDenseMatrix& ivecnpderiv  //!< reaction term derivaitve
)
{
  // Initialize reaction and diffusion matrices
  hdgele->Ivecn_.size(hdgele->ndofs_);
  hdgele->Ivecnp_.size(hdgele->ndofs_);
  hdgele->Imatnpderiv_.shape(hdgele->ndofs_, hdgele->ndofs_);


  hdgele->Ivecn_ = ivecn;
  hdgele->Ivecnp_ = ivecnp;
  hdgele->Imatnpderiv_ = ivecnpderiv;
}

/*----------------------------------------------------------------------*
 |  Prepare material parameter                           hoermann 11/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::LocalSolver::prepare_material_parameter(
    DRT::ELEMENTS::ScaTraHDG* hdgele,            //!< hdg element
    CORE::LINALG::SerialDenseMatrix& difftensor  //!< diffusion tensor
)
{
  using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
  using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseindifftensor;
  inverseindifftensor.setMatrix(Teuchos::rcpFromRef(difftensor));
  int err = inverseindifftensor.invert();
  if (err != 0) FOUR_C_THROW("Inversion of diffusion tensor failed with errorcode %d", err);

  hdgele->invdiff_.push_back(difftensor);
}

/*----------------------------------------------------------------------*
 * Element Init
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::element_init(CORE::Elements::Element* ele)
{
  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<CORE::Elements::Element*>(ele));

  hdgele->Amat_.shape(hdgele->ndofs_, hdgele->ndofs_);
  hdgele->Bmat_.shape(hdgele->ndofs_, nsd_ * hdgele->ndofs_);
  hdgele->Cmat_.shape(hdgele->ndofs_, hdgele->onfdofs_);
  hdgele->Dmat_.shape(nsd_ * hdgele->ndofs_, nsd_ * hdgele->ndofs_);
  hdgele->Emat_.shape(hdgele->ndofs_ * nsd_, hdgele->onfdofs_);
  hdgele->Gmat_.shape(hdgele->onfdofs_, hdgele->ndofs_);
  hdgele->EmatT_.shape(hdgele->onfdofs_, nsd_ * hdgele->ndofs_);
  hdgele->Hmat_.shape(hdgele->onfdofs_, hdgele->onfdofs_);
  hdgele->Mmat_.shape(hdgele->ndofs_, hdgele->ndofs_);
  hdgele->Kmat_.shape(hdgele->onfdofs_, hdgele->onfdofs_);
  hdgele->Xmat_.shape(hdgele->ndofs_, hdgele->onfdofs_);
  hdgele->BmatMT_.shape(nsd_ * hdgele->ndofs_, hdgele->ndofs_);
  hdgele->invAMmat_.shape(hdgele->ndofs_, hdgele->ndofs_);
  hdgele->BTAMmat_.shape(hdgele->ndofs_ * nsd_, hdgele->ndofs_);
  hdgele->invCondmat_.shape(hdgele->ndofs_ * nsd_, hdgele->ndofs_ * nsd_);
  hdgele->diff_.shape(nsd_, nsd_);
  hdgele->invdiff_.clear();
  hdgele->Ivecn_.size(hdgele->ndofs_);
  hdgele->Ivecnp_.size(hdgele->ndofs_);
  hdgele->Imatnpderiv_.shape(hdgele->ndofs_, hdgele->ndofs_);

  hdgele->SetMatInit(true);

  return;
}  // ElementInit

/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::ProjectField(
    const CORE::Elements::Element* ele, DRT::Discretization& discretization,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseVector& elevec2, CORE::Elements::Element::LocationArray& la)
{
  int nds_var_old = params.get<int>("nds_var_old");
  int nds_intvar_old = params.get<int>("nds_intvar_old");

  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<CORE::Elements::Element*>(ele));

  // set change of element degree to false
  hdgele->SetPadaptEle(false);

  Teuchos::RCP<CORE::FE::ShapeValues<distype>> shapes_old =
      Teuchos::rcp(new CORE::FE::ShapeValues<distype>(
          hdgele->DegreeOld(), usescompletepoly_, 2 * hdgele->DegreeOld()));

  FOUR_C_ASSERT(
      elevec2.numRows() == 0 || unsigned(elevec2.numRows()) == shapes_->ndofs_ * (nsd_ + 1),
      "Wrong size in project vector 2");

  // polynomial space to get the value of the shape function at the new points
  CORE::FE::PolynomialSpaceParams params_old(distype, shapes_old->degree_, usescompletepoly_);
  Teuchos::RCP<CORE::FE::PolynomialSpace<probdim>> polySpace_old =
      CORE::FE::PolynomialSpaceCache<probdim>::Instance().Create(params_old);

  CORE::LINALG::SerialDenseVector interiorPhi_old(shapes_old->ndofs_ * (nsd_ + 1));

  // get node based values!
  Teuchos::RCP<const Epetra_Vector> matrix_state = params.get<Teuchos::RCP<Epetra_Vector>>("phi");

  std::vector<double> tracephi;
  CORE::FE::ExtractMyValues(*matrix_state, tracephi, la[nds_var_old].lm_);

  // get node based values!
  matrix_state = params.get<Teuchos::RCP<Epetra_Vector>>("intphi");
  std::vector<double> intphi;
  CORE::FE::ExtractMyValues(*matrix_state, intphi, la[nds_intvar_old].lm_);
  if (intphi.size() != shapes_old->ndofs_ * (nsd_ + 1))
    FOUR_C_THROW(
        "node number not matching: %d vs. %d", intphi.size(), shapes_old->ndofs_ * (nsd_ + 1));

  for (unsigned int i = 0; i < shapes_old->ndofs_ * (nsd_ + 1); ++i) interiorPhi_old(i) = intphi[i];

  if (!usescompletepoly_)
  {
    // copy values if degree stays the same instead of projecting
    if (hdgele->DegreeOld() == hdgele->Degree())
      for (unsigned int i = 0; i < shapes_old->ndofs_ * (nsd_ + 1); ++i) elevec2(i) = intphi[i];
    else
    {
      // set change of element degree to true
      hdgele->SetPadaptEle(true);
      CORE::LINALG::SerialDenseMatrix tempMat(
          shapes_->ndofs_ * (nsd_ + 1), shapes_old->ndofs_ * (nsd_ + 1));

      for (unsigned int i = 0; i < shapes_->ndofs_; i++)
      {
        CORE::LINALG::SerialDenseVector tempVec(shapes_old->ndofs_);
        CORE::LINALG::Matrix<nsd_, 1> point(shapes_->nodexyzunit[i]);
        polySpace_old->Evaluate(point, tempVec);
        for (unsigned int j = 0; j < nsd_ + 1; j++)
          for (unsigned int k = 0; k < shapes_old->ndofs_; k++)
            tempMat(j * shapes_->ndofs_ + i, j * shapes_old->ndofs_ + k) = tempVec(k);
      }
      CORE::LINALG::multiply(elevec2, tempMat, interiorPhi_old);
    }

    // project trace field
    int nfdofs = 0;
    int nfdofs_old = 0;
    for (unsigned int face = 0; face < nfaces_; face++)
    {
      // shape values of new element degree
      CORE::FE::ShapeValuesFaceParams svfparams(
          ele->Faces()[face]->Degree(), usescompletepoly_, 2 * ele->Faces()[face]->Degree());
      shapesface_ = CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      shapesface_->EvaluateFace(*ele, face);


      // shape values of old element degree
      DRT::ELEMENTS::ScaTraHDGIntFace* hdgeleface = dynamic_cast<DRT::ELEMENTS::ScaTraHDGIntFace*>(
          const_cast<CORE::Elements::FaceElement*>(ele->Faces()[face].getRawPtr()));
      CORE::FE::ShapeValuesFaceParams svfparams_old(
          hdgeleface->DegreeOld(), usescompletepoly_, 2 * hdgeleface->DegreeOld());

      Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>> shapesface_old =
          CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams_old);


      CORE::FE::PolynomialSpaceParams polyparams(CORE::FE::DisTypeToFaceShapeType<distype>::shape,
          hdgeleface->DegreeOld(), usescompletepoly_);
      Teuchos::RCP<CORE::FE::PolynomialSpace<nsd_ - 1>> polySpaceFace_old =
          CORE::FE::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(polyparams);

      CORE::LINALG::SerialDenseVector tracePhi_face_old(shapesface_old->nfdofs_);

      for (unsigned int i = 0; i < shapesface_old->nfdofs_; i++)
        tracePhi_face_old(i) = tracephi[nfdofs_old + i];

      if (ele->Faces()[face]->Degree() == hdgeleface->DegreeOld())
        for (unsigned int i = 0; i < shapesface_old->nfdofs_; i++)
          elevec1(nfdofs + i) = tracephi[nfdofs_old + i];
      else
      {
        // set change of element degree to true
        hdgele->SetPadaptEle(true);

        CORE::LINALG::SerialDenseMatrix tempMat1(shapesface_->nfdofs_, shapesface_old->nfdofs_);

        CORE::LINALG::SerialDenseVector tempVec2(shapesface_->nfdofs_);

        for (unsigned int i = 0; i < shapesface_->nfdofs_; i++)
        {
          CORE::LINALG::SerialDenseVector tempVec(shapesface_old->nfdofs_);
          CORE::LINALG::Matrix<nsd_ - 1, 1> point(shapesface_->nodexyzunit[i]);

          polySpaceFace_old->Evaluate(point, tempVec);
          for (unsigned int k = 0; k < shapesface_old->nfdofs_; k++) tempMat1(i, k) = tempVec(k);
        }

        CORE::LINALG::multiply(tempVec2, tempMat1, tracePhi_face_old);

        for (unsigned int i = 0; i < shapesface_->nfdofs_; i++) elevec1(nfdofs + i) = tempVec2(i);
      }

      nfdofs += shapesface_->nfdofs_;
      nfdofs_old += shapesface_old->nfdofs_;
    }
  }
  else
  {
    if (hdgele->DegreeOld() != hdgele->Degree())
      hdgele->SetPadaptEle(true);  // set change of element degree to true

    unsigned int size_ndofs = std::min(shapes_old->ndofs_, shapes_->ndofs_);
    for (unsigned int i = 0; i < nsd_ + 1; i++)
      for (unsigned int j = 0; j < size_ndofs; j++)
        elevec2(i * shapes_->ndofs_ + j) = interiorPhi_old(i * shapes_old->ndofs_ + j);

    int nfdofs = 0;
    int nfdofs_old = 0;
    for (unsigned int face = 0; face < nfaces_; face++)
    {
      // shape values of new element degree
      CORE::FE::ShapeValuesFaceParams svfparams(
          ele->Faces()[face]->Degree(), usescompletepoly_, 2 * ele->Faces()[face]->Degree());
      shapesface_ = CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      shapesface_->EvaluateFace(*ele, face);

      // shape values of old element degree
      DRT::ELEMENTS::ScaTraHDGIntFace* hdgeleface = dynamic_cast<DRT::ELEMENTS::ScaTraHDGIntFace*>(
          const_cast<CORE::Elements::FaceElement*>(ele->Faces()[face].getRawPtr()));
      CORE::FE::ShapeValuesFaceParams svfparams_old(
          hdgeleface->DegreeOld(), usescompletepoly_, 2 * hdgeleface->DegreeOld());

      Teuchos::RCP<CORE::FE::ShapeValuesFace<distype>> shapesface_old =
          CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams_old);


      CORE::FE::PolynomialSpaceParams polyparams(CORE::FE::DisTypeToFaceShapeType<distype>::shape,
          hdgeleface->DegreeOld(), usescompletepoly_);
      Teuchos::RCP<CORE::FE::PolynomialSpace<nsd_ - 1>> polySpaceFace_old =
          CORE::FE::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(polyparams);

      //      CORE::LINALG::SerialDenseVector tracePhi_face_old(shapesface_old->nfdofs_);

      if (ele->Faces()[face]->Degree() != hdgeleface->DegreeOld())
        hdgele->SetPadaptEle(true);  // set change of element degree to true

      unsigned int size_nfdofs = std::min(shapesface_->nfdofs_, shapesface_old->nfdofs_);
      for (unsigned int i = 0; i < size_nfdofs; i++) elevec1(nfdofs + i) = tracephi[nfdofs_old + i];

      nfdofs += shapesface_->nfdofs_;
      nfdofs_old += shapesface_old->nfdofs_;
    }
  }

  return 0;

}  // ProjectField


/*----------------------------------------------------------------------*
 * Calc P-Adaptivity
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::CalcPAdaptivity(
    const CORE::Elements::Element* ele, DRT::Discretization& discretization,
    Teuchos::ParameterList& params)
{
  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<CORE::Elements::Element*>(ele));

  CORE::LINALG::SerialDenseVector tempinteriorgradphinp(hdgele->ndofs_ * nsd_);
  for (unsigned int i = 0; i < hdgele->ndofs_ * nsd_; ++i)
    tempinteriorgradphinp(i) = interiorPhinp_(hdgele->ndofs_ + i);

  double error = 0;
  int sumindex = 0;
  for (unsigned int nface = 0; nface < nfaces_; ++nface)
  {
    CORE::FE::ShapeValuesFaceParams svfparams(ele->Faces()[nface]->Degree(),
        shapes_->usescompletepoly_, 2 * ele->Faces()[nface]->Degree());

    shapesface_ = CORE::FE::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, nface);

    CORE::LINALG::SerialDenseMatrix QMat(shapesface_->nqpoints_, hdgele->ndofs_ * nsd_);
    CORE::LINALG::SerialDenseMatrix QMatW(shapesface_->nqpoints_, hdgele->ndofs_ * nsd_);
    CORE::LINALG::SerialDenseMatrix UMat(
        shapesface_->nqpoints_, hdgele->ndofs_ + shapesface_->nfdofs_);
    CORE::LINALG::SerialDenseMatrix UMatW(
        shapesface_->nqpoints_, hdgele->ndofs_ + shapesface_->nfdofs_);

    CORE::LINALG::SerialDenseVector tempinteriorphinp(hdgele->ndofs_ + shapesface_->nfdofs_);
    for (unsigned int i = 0; i < hdgele->ndofs_; i++) tempinteriorphinp(i) = interiorPhinp_(i);
    for (unsigned int i = 0; i < shapesface_->nfdofs_; i++)
      tempinteriorphinp(hdgele->ndofs_ + i) = tracen_(sumindex + i);

    // loop over quadrature points
    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      // loop over shape functions
      for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
      {
        UMat(q, i) = shapesface_->shfunctI(i, q);
        UMatW(q, i) = shapesface_->shfunctI(i, q) * shapesface_->jfac(q);
        for (unsigned int k = 0; k < nsd_; ++k)
        {
          QMat(q, hdgele->ndofs_ * k + i) = shapesface_->shfunctI(i, q) * shapesface_->normal(k);
          QMatW(q, hdgele->ndofs_ * k + i) =
              shapesface_->shfunctI(i, q) * shapesface_->jfac(q) * shapesface_->normal(k);
        }
      }
      // loop over face shape functions
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        UMat(q, i + hdgele->ndofs_) = -shapesface_->shfunct(i, q);
        UMatW(q, i + hdgele->ndofs_) = -shapesface_->shfunct(i, q) * shapesface_->jfac(q);
      }
    }
    sumindex += shapesface_->nfdofs_;

    CORE::LINALG::SerialDenseVector tempVec1(shapesface_->nqpoints_);
    CORE::LINALG::SerialDenseVector tempVec2(shapesface_->nqpoints_);
    CORE::LINALG::SerialDenseVector tempVec3(shapesface_->nqpoints_);
    CORE::LINALG::SerialDenseVector tempVec4(shapesface_->nqpoints_);


    CORE::LINALG::multiply(tempVec1, QMatW, tempinteriorgradphinp);
    CORE::LINALG::multiply(tempVec2, QMat, tempinteriorgradphinp);
    CORE::LINALG::multiply(tempVec3, UMatW, tempinteriorphinp);
    CORE::LINALG::multiply(tempVec4, UMat, tempinteriorphinp);

    double errorface = 0;
    double facearea = 0;

    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      errorface +=
          tempVec1(q) * tempVec2(q) + tempVec3(q) * tempVec4(q) - 2 * tempVec1(q) * tempVec4(q);
      facearea += shapesface_->jfac(q);
    }

    // normalize error with surface area of face
    error += errorface / facearea;
  }

  params.set<double>("error", error);


  return 0;
}  // CalcPAdaptivity

/*----------------------------------------------------------------------*
 * Calc Error
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype, probdim>::CalcError(const CORE::Elements::Element* ele,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& elevec)
{
  DRT::ELEMENTS::ScaTraHDG* hdgele =
      dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<CORE::Elements::Element*>(ele));

  // For the calculation of the error we use a higher integration rule
  CORE::FE::ShapeValues<distype> highshapes(
      ele->Degree(), shapes_->usescompletepoly_, (ele->Degree() + 2) * 2);
  highshapes.Evaluate(*ele);

  double error_phi = 0.0, error_grad_phi = 0.0;
  double exact_phi = 0.0, exact_grad_phi = 0.0;

  // get function
  const int func = params.get<int>("error function number");
  const double time = params.get<double>("time");

  if (GLOBAL::Problem::Instance()
          ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(func - 1)
          .NumberComponents() != 1)
    FOUR_C_THROW(
        "The number of component must be one. The grandient is computed with forward auomatic "
        "differentiation.");

  CORE::LINALG::Matrix<nsd_, 1> xsi;
  double phi(nsd_);
  CORE::LINALG::SerialDenseVector gradPhi(nsd_);

  for (unsigned int q = 0; q < highshapes.nqpoints_; ++q)
  {
    phi = 0;
    gradPhi.putScalar(0.0);
    if (hdgele->invdiff_.size() == 1)
      for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      {
        phi += highshapes.shfunct(i, q) * interiorPhinp_(i);
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            gradPhi(d) += highshapes.shfunct(i, q) * interiorPhinp_(i + (e + 1) * shapes_->ndofs_) *
                          hdgele->invdiff_[0](d, e);
      }
    else if (hdgele->invdiff_.size() == highshapes.nqpoints_)
      for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      {
        phi += highshapes.shfunct(i, q) * interiorPhinp_(i);
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            gradPhi(d) += highshapes.shfunct(i, q) * interiorPhinp_(i + (e + 1) * shapes_->ndofs_) *
                          hdgele->invdiff_[q](d, e);
      }
    else
      FOUR_C_THROW("Diffusion tensor not defined properly. Impossible to compute error.");


    // Analytical function evaluation
    // Evaluate error function and its derivatives in the integration point (real) coordinates
    for (unsigned int idim = 0; idim < nsd_; idim++) xsi(idim) = highshapes.xyzreal(idim, q);
    double funct = GLOBAL::Problem::Instance()
                       ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(func - 1)
                       .Evaluate(xsi.A(), time, 0);
    std::vector<double> deriv = GLOBAL::Problem::Instance()
                                    ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(func - 1)
                                    .evaluate_spatial_derivative(xsi.A(), time, 0);

    error_phi += std::pow((funct - phi), 2) * highshapes.jfac(q);
    exact_phi += std::pow(funct, 2) * highshapes.jfac(q);
    for (unsigned int d = 0; d < nsd_; ++d)
    {
      error_grad_phi += std::pow((deriv[d] - gradPhi(d)), 2) * highshapes.jfac(q);
      exact_grad_phi += std::pow(deriv[d], 2) * highshapes.jfac(q);
    }
  }

  elevec[0] = error_phi;
  elevec[1] = exact_phi;
  elevec[2] = error_grad_phi;
  elevec[3] = exact_grad_phi;

  return 0;
}

// template classes
// 1D elements
// template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::line2,1>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::line2,2>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::line2,3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::line3,1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::tri3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcHDG<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
