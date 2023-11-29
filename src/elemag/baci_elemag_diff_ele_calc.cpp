/*----------------------------------------------------------------------*/
/*! \file
\brief Electromagnetic diffusion equation element implementation

<pre>
\level 2

</pre>
*/
/*--------------------------------------------------------------------------*/

#include "baci_elemag_diff_ele_calc.H"

#include "baci_discretization_fem_general_utils_boundary_integration.H"
#include "baci_discretization_geometry_position_array.H"
#include "baci_elemag_ele_action.H"
#include "baci_lib_discret.H"
#include "baci_lib_discret_hdg.H"
#include "baci_lib_elementtype.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_utils_densematrix_multiply.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_mat_electromagnetic.H"
#include "baci_utils_function.H"

#include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_TimeMonitor.hpp>


/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ElemagDiffEleCalc<distype>::ElemagDiffEleCalc()
{
}

/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ElemagDiffEleCalc<distype>::Evaluate(DRT::ELEMENTS::Elemag* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra, const CORE::DRT::UTILS::GaussIntegration&,
    bool offdiag)
{
  return this->Evaluate(ele, discretization, lm, params, mat, elemat1_epetra, elemat2_epetra,
      elevec1_epetra, elevec2_epetra, elevec3_epetra, offdiag);
}

/*----------------------------------------------------------------------*
 * Evaluate
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ElemagDiffEleCalc<distype>::Evaluate(DRT::ELEMENTS::Elemag* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1,
    CORE::LINALG::SerialDenseMatrix&, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector&, bool offdiag)
{
  // check if this is an hdg element and init completepoly
  DRT::ELEMENTS::ElemagDiff* hdgele = static_cast<DRT::ELEMENTS::ElemagDiff*>(ele);
  usescompletepoly_ = hdgele->UsesCompletePolynomialSpace();
  // else dserror("cannot cast element to elemag element");

  ELEMAG::Action action;
  if (params.get<bool>("hdg_action", false))
  {
    switch (Teuchos::getIntegralValue<DRT::HDGAction>(params, "action"))
    {
      case DRT::HDGAction::project_dirich_field:
        action = ELEMAG::Action::project_dirich_field;
        break;
      default:
        dserror("HDG Action type not supported");
    }
  }
  else
  {
    action = DRT::INPUT::get<ELEMAG::Action>(params, "action");
  }

  InitializeShapes(hdgele);

  bool updateonly = false;
  shapes_->Evaluate(*ele);
  switch (action)
  {
    case ELEMAG::project_field:
    {
      localSolver_->ProjectField(hdgele, params, elevec1, elevec2);
      break;
    }
    case ELEMAG::project_electric_from_scatra_field:
    {
      ElementInit(hdgele, params);
      localSolver_->ProjectElectricFieldFromScatra(hdgele, params, mat, elevec1);
      break;
    }
    case ELEMAG::compute_error:
    {
      // Postprocess the solution only if required
      if (params.get<bool>("postprocess")) localSolver_->PostProcessing(*hdgele);

      localSolver_->ComputeError(hdgele, params, elevec1);
      break;
    }
    case ELEMAG::project_field_test:
    {
      localSolver_->ProjectFieldTest(hdgele, params, elevec1, elevec2);
      break;
    }
    case ELEMAG::project_field_test_trace:
    {
      localSolver_->ProjectFieldTestTrace(hdgele, params, elevec1);

      break;
    }
    case ELEMAG::project_dirich_field:
    {
      // if (mat->MaterialType() != INPAR::MAT::m_electromagneticmat)
      //  dserror("for physical type 'lossless' please supply MAT_Electromagnetic");
      if (params.isParameter("faceconsider"))
      {
        // ElementInit(hdgele, params);
        localSolver_->ProjectDirichField(hdgele, params, elevec1);
      }
      break;
    }
    case ELEMAG::ele_init:
    {
      ElementInit(hdgele, params);
      break;
    }
    case ELEMAG::fill_restart_vecs:
    {
      // bool padapty = params.get<bool>("padaptivity");
      ReadGlobalVectors(hdgele, discretization, lm);
      FillRestartVectors(hdgele, discretization);
      break;
    }
    case ELEMAG::ele_init_from_restart:
    {
      ElementInitFromRestart(hdgele, discretization);
      break;
    }
    case ELEMAG::interpolate_hdg_to_node:
    {
      ReadGlobalVectors(hdgele, discretization, lm);
      // The post processing is only used to compute the error so far as it is computationally
      // extremely expensive.
      // TODO: Improve post processing efficiency
      // localSolver_->PostProcessing(*hdgele);
      InterpolateSolutionToNodes(hdgele, discretization, elevec1);
      break;
    }
    case ELEMAG::calc_abc:
    {
      int face = params.get<int>("face");
      int sumindex = 0;
      for (int i = 0; i < face; ++i)
      {
        CORE::DRT::UTILS::PolynomialSpaceParams params(
            CORE::DRT::UTILS::DisTypeToFaceShapeType<distype>::shape, hdgele->Faces()[i]->Degree(),
            usescompletepoly_);
        int nfdofs =
            CORE::DRT::UTILS::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(params)->Size();
        sumindex += nfdofs;
      }
      ReadGlobalVectors(hdgele, discretization, lm);
      if (!params.isParameter("nodeindices"))
        localSolver_->ComputeAbsorbingBC(
            discretization, hdgele, params, mat, face, elemat1, sumindex, elevec1);
      else
        dserror("why would you set an absorbing LINE in THREE dimensions?");

      break;
    }
    /*
    case ELEMAG::bd_integrate:
    {
      int face = params.get<int>("face");
      localSolver_->ComputeBoundaryIntegral(hdgele, params, face);

      break;
    }
    */
    case ELEMAG::calc_systemmat_and_residual:
    {
      // const bool resonly = params.get<bool>("resonly");
      // const bool padapty = params.get<bool>("padaptivity");
      double dt = params.get<double>("dt");
      const double tau = params.get<double>("tau");
      dyna_ = params.get<INPAR::ELEMAG::DynamicType>("dynamic type");

      // Compute RHS factor
      const MAT::ElectromagneticMat* elemagmat =
          static_cast<const MAT::ElectromagneticMat*>(mat.get());
      const double mu = elemagmat->mu(hdgele->Id());
      if (mu < 0.1)
        params.set<double>("mod_mu", std::pow(mu, 0.5 * (1 + std::log(dt) / std::log(mu))));
      else
        params.set<double>("mod_mu", std::pow(mu, 0.0));

      ReadGlobalVectors(hdgele, discretization, lm);
      elevec1.putScalar(0.0);
      localSolver_->ComputeMatrices(discretization, mat, *hdgele, dt, dyna_, tau);

      // if (!resonly)
      localSolver_->CondenseLocalPart(elemat1);

      // Make the matrix symmetric if the element contains dirichlet faces
      localSolver_->Symmetrify(*hdgele, elemat1);

      localSolver_->ComputeResidual(params, elevec1, dt, *hdgele);

      break;
    }
    case ELEMAG::update_secondary_solution:
      updateonly = true;  // no break here!!!
      [[fallthrough]];
    case ELEMAG::update_secondary_solution_and_calc_residual:
    {
      // bool errormaps = params.get<bool>("errormaps");
      bool errormaps = false;
      // const bool allelesequal = params.get<bool>("allelesequal");

      const double dt = params.get<double>("dt");
      const double tau = params.get<double>("tau");
      dyna_ = params.get<INPAR::ELEMAG::DynamicType>("dynamic type");

      // Compute RHS factor
      const MAT::ElectromagneticMat* elemagmat =
          static_cast<const MAT::ElectromagneticMat*>(mat.get());
      const double mu = elemagmat->mu(hdgele->Id());
      if (mu < 0.1)
        params.set<double>("mod_mu", std::pow(mu, 0.5 * (1 + std::log(dt) / std::log(mu))));
      else
        params.set<double>("mod_mu", std::pow(mu, 0.0));

      ReadGlobalVectors(hdgele, discretization, lm);

      elevec1.putScalar(0.0);
      localSolver_->ComputeMatrices(discretization, mat, *hdgele, dt, dyna_, tau);
      /* Could be useful for optimization purposes
      if(!allelesequal)
        localSolver_->ComputeMatrices(discretization, mat, *hdgele, dt, dyna_);
      */

      UpdateInteriorVariablesAndComputeResidual(
          params, *hdgele, mat, elevec1, dt, errormaps, updateonly);

      break;
    }
    case ELEMAG::get_gauss_points:
    {
      int rows = shapes_->xyzreal.numRows();
      int cols = shapes_->xyzreal.numCols();
      elemat1.shape(rows, cols);

      for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c) elemat1(r, c) = shapes_->xyzreal(r, c);

      break;
    }
    default:
    {
      std::cout << "Action: " << action << std::endl;
      dserror("unknown action supplied");
      break;
    }
  }  // switch(action)

  return 0;
}

/*----------------------------------------------------------------------*
 * Print trace
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::PrintTrace(DRT::Element* ele)
{
  std::cout << "Local trace of element: " << ele->LID() << std::endl;
  std::cout << "Number of entries: " << localtrace_.size() << std::endl;
  std::cout << "Number of spatial dimensions: " << nsd_ << std::endl;
  std::cout << "Numer of faces: " << nfaces_ << std::endl;
  std::cout << "Numer of DOF per face: " << ele->NumDofPerFace(0) << std::endl;
  unsigned int index = 0;
  unsigned int second_index = 0;
  for (std::vector<double>::iterator iter = localtrace_.begin(); iter != localtrace_.end();
       iter++, index++, second_index++)
  {
    if (index % ele->NumDofPerFace(0) == 0)
    {
      std::cout << "Face number: " << index / ele->NumDofPerFace(0) << std::endl;
      second_index = 0;
    }
    if (second_index % shapesface_->nfdofs_ == 0)
      std::cout << "\tField component: " << second_index / shapesface_->nfdofs_ << std::endl;
    std::cout << "\t\t" << *iter << std::endl;
  }
  return;
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::InitializeShapes(
    const DRT::ELEMENTS::ElemagDiff* ele)
{
  if (shapes_ == Teuchos::null || shapes_->degree_ != unsigned(ele->Degree()) ||
      shapes_->usescompletepoly_ != usescompletepoly_)
  {
    shapes_ = Teuchos::rcp(new CORE::DRT::UTILS::ShapeValues<distype>(
        ele->Degree(), usescompletepoly_, 2 * ele->Degree()));

    // TODO: Check wheter 2 * (ele->Degree()+1) is required as exact integration degree
    postproc_shapes_ = Teuchos::rcp(new CORE::DRT::UTILS::ShapeValues<distype>(
        ele->Degree() + 1, usescompletepoly_, 2 * (ele->Degree() + 1)));
  }

  if (shapesface_ == Teuchos::null)
  {
    CORE::DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele->Degree(), usescompletepoly_, 2 * ele->Degree());
    shapesface_ = Teuchos::rcp(new CORE::DRT::UTILS::ShapeValuesFace<distype>(svfparams));
  }

  if (localSolver_ == Teuchos::null || localSolver_->ndofs_ != shapes_->ndofs_)
  {
    localSolver_ =
        Teuchos::rcp(new LocalSolver(ele, *shapes_, shapesface_, dyna_, *postproc_shapes_));
  }
}

/*----------------------------------------------------------------------*
 * ReadGlobalVectors
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::ReadGlobalVectors(
    DRT::Element* ele, DRT::Discretization& discretization, const std::vector<int>& lm)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagDiffEleCalc::ReadGlobalVectors");
  DRT::ELEMENTS::ElemagDiff* elemagele = static_cast<DRT::ELEMENTS::ElemagDiff*>(ele);

  // read vectors from element storage
  interiorElectricnm_.size(elemagele->eleinteriorElectricnm1_.numRows());
  interiorElectricnp_.size(elemagele->eleinteriorElectric_.numRows());
  interiorMagneticnp_.size(elemagele->eleinteriorMagnetic_.numRows());

  interiorElectricnm_ = elemagele->eleinteriorElectricnm1_;
  interiorElectricnp_ = elemagele->eleinteriorElectric_;
  interiorMagneticnp_ = elemagele->eleinteriorMagnetic_;

  // read vectors from time integrator
  if (discretization.HasState("trace"))  // in case of "update interior variables"
  {
    elemagele->elenodeTrace2d_.size(lm.size());
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState("trace");
    DRT::UTILS::ExtractMyValues(*matrix_state, elemagele->elenodeTrace2d_, lm);
  }

  return;
}  // ReadGlobalVectors

/*----------------------------------------------------------------------*
 * FillRestartVectors
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::FillRestartVectors(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  // sort this back to the interior values vector
  int size = shapes_->ndofs_ * nsd_;

  DRT::ELEMENTS::Elemag* hdgele = static_cast<DRT::ELEMENTS::Elemag*>(ele);

  std::vector<double> interiorVar(size * 2);
  for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
  {
    interiorVar[i] = hdgele->eleinteriorMagnetic_(i);
    interiorVar[shapes_->ndofs_ * nsd_ + i] = hdgele->eleinteriorElectric_(i);
  }

  // tell this change in the interior variables the discretization
  std::vector<int> localDofs = discretization.Dof(1, ele);
  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "intVar");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for (unsigned int i = 0; i < localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorVar[i];
    }
  }

  std::vector<double> interiorVarnm(size);
  for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
  {
    interiorVarnm[i] = hdgele->eleinteriorElectricnm1_(i);
  }

  // Here the magnetic field is not stored because there is no need for the time integration
  Teuchos::RCP<const Epetra_Vector> intVarnm = discretization.GetState(1, "intVarnm");
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*intVarnm);
  for (unsigned int i = size; i < localDofs.size(); ++i)
  {
    const int lid = intdofcolmap->LID(localDofs[i]);
    secondary[lid] = interiorVarnm[i - size];
  }

  return;
}

/*----------------------------------------------------------------------*
 * ElementInitFromRestart
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::ElementInitFromRestart(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  DRT::ELEMENTS::ElemagDiff* elemagele = dynamic_cast<DRT::ELEMENTS::ElemagDiff*>(ele);
  unsigned int size = shapes_->ndofs_ * nsd_;

  std::vector<double> interiorVar(size * 2);

  Teuchos::RCP<const Epetra_Vector> intVar = discretization.GetState(1, "intVar");
  std::vector<int> localDofs1 = discretization.Dof(1, ele);
  DRT::UTILS::ExtractMyValues(*intVar, interiorVar, localDofs1);
  // now write this in corresponding eleinteriorElectric_ and eleinteriorMagnetic_
  for (unsigned int i = 0; i < size; ++i)
  {
    elemagele->eleinteriorMagnetic_(i) = interiorVar[i];
    elemagele->eleinteriorElectric_(i) = interiorVar[size + i];
  }

  std::vector<double> interiorVarnm(size * 2);

  Teuchos::RCP<const Epetra_Vector> intVarnm = discretization.GetState(1, "intVarnm");
  DRT::UTILS::ExtractMyValues(*intVarnm, interiorVarnm, localDofs1);
  for (unsigned int i = 0; i < size; ++i)
  {
    elemagele->eleinteriorElectricnm1_(i) = interiorVarnm[size + i];
  }

  return;
}

/*----------------------------------------------------------------------*
 * Element init
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::ElementInit(
    DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params)
{
  // Save the position of element's dirichlet dofs
  if (params.isParameter("dirichdof"))
    ele->lm_.lmdirich_ = *params.get<std::vector<int>*>("dirichdof");

  dyna_ = params.get<INPAR::ELEMAG::DynamicType>("dyna");

  // each element has to store the interior vectors by itseld, p-adaptivity or not
  // so, shape it, as you need it
  DRT::ELEMENTS::ElemagDiff* diff_ele = static_cast<DRT::ELEMENTS::ElemagDiff*>(ele);
  diff_ele->eleinteriorElectricnm3_.size(shapes_->ndofs_ * nsd_);
  diff_ele->eleinteriorElectricnm2_.size(shapes_->ndofs_ * nsd_);
  diff_ele->eleinteriorElectricnm1_.size(shapes_->ndofs_ * nsd_);
  diff_ele->eleinteriorElectric_.size(shapes_->ndofs_ * nsd_);
  diff_ele->eleinteriorMagnetic_.size(shapes_->ndofs_ * nsd_);

  // ele->elenodeTrace_.Size(ele->NumFace() * shapesface_->nfdofs_ * nsd_);
  ele->elenodeTrace2d_.size(ele->NumFace() * shapesface_->nfdofs_ * (nsd_ - 1));

  // Postproc
  ele->eleinteriorElectricPost_.size(postproc_shapes_->ndofs_ * nsd_);

  return;
}

/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ProjectField(
    DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2)
{
  shapes_.Evaluate(*ele);

  // get function
  const int* start_func = params.getPtr<int>("startfuncno");
  const double time = params.get<double>("time");

  // the RHS matrix has to have the row dimension equal to the number of shape
  // functions(so we have one coefficient for each) and a number of column
  // equal to the overall number of component that we want to solve for.
  // The number is nsd_ because we have one field.
  CORE::LINALG::SerialDenseMatrix localMat(shapes_.ndofs_, nsd_);
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    // Storing the values of the coordinates for the current quadrature point
    // and of the jacobian computed in that point
    const double fac = shapes_.jfac(q);
    CORE::LINALG::Matrix<nsd_, 1> xyz;
    for (unsigned int d = 0; d < nsd_; ++d)
      xyz(d) = shapes_.xyzreal(d, q);  // coordinates of quadrature point in real coordinates
    // Creating the temporary electric and magnetic field vector intVal
    // The vector is going to contain first the electric and then the magnetic
    // field such that the field will be initialized as first tree component
    // of the specified function as electric field, last three components as
    // magnetic field. If there is only one component all the components will
    // be initialized to the same value.
    CORE::LINALG::SerialDenseVector intVal(nsd_);
    dsassert(start_func != nullptr, "funct not set for initial value");
    EvaluateAll(*start_func, time, xyz, intVal);
    // now fill the components in the one-sided mass matrix and the right hand side
    for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
    {
      // Mass matrix
      massPart(i, q) = shapes_.shfunct(i, q);
      massPartW(i, q) = shapes_.shfunct(i, q) * fac;

      // RHS for the electric and magnetic field
      for (int j = 0; j < intVal.numRows(); ++j)
        localMat(i, j) += shapes_.shfunct(i, q) * intVal(j) * fac;
    }
  }
  // The integration is made by computing the matrix product
  CORE::LINALG::multiplyNT(massMat, massPart, massPartW);
  {
    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
    inverseMass.setMatrix(Teuchos::rcpFromRef(massMat));
    inverseMass.setVectors(Teuchos::rcpFromRef(localMat), Teuchos::rcpFromRef(localMat));
    inverseMass.solve();
  }

  // Here we move the values from the temporary variable to the variable
  // contained in the element
  for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
  {
    // Now we are storing the variables by component, meaning that we save for
    // each component the value for each dof and then we move to the next component.
    for (unsigned int d = 0; d < nsd_; ++d)
    {
      ele->eleinteriorElectric_(d * shapes_.ndofs_ + r) = localMat(r, d);  // Electric field
    }
  }

  if (dyna_ == INPAR::ELEMAG::elemag_bdf4)
    for (int s = 1; s < 4; s++)
    {
      localMat.putScalar(0.0);
      const double dt = params.get<double>("dt");
      for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
      {
        const double fac = shapes_.jfac(q);
        CORE::LINALG::Matrix<nsd_, 1> xyz;
        for (unsigned int d = 0; d < nsd_; ++d) xyz(d) = shapes_.xyzreal(d, q);

        CORE::LINALG::SerialDenseVector intVal(nsd_);

        EvaluateAll(*start_func, time - s * dt, xyz, intVal);
        for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
        {
          massPart(i, q) = shapes_.shfunct(i, q);
          massPartW(i, q) = shapes_.shfunct(i, q) * fac;
          for (int j = 0; j < intVal.numRows(); ++j)
            localMat(i, j) += shapes_.shfunct(i, q) * intVal(j) * fac;
        }
      }


      // The integration is made by computing the matrix product
      CORE::LINALG::multiplyNT(massMat, massPart, massPartW);
      {
        using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
        using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
        Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
        inverseMass.setMatrix(Teuchos::rcpFromRef(massMat));
        inverseMass.setVectors(Teuchos::rcpFromRef(localMat), Teuchos::rcpFromRef(localMat));
        inverseMass.solve();
      }

      for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
        for (unsigned int d = 0; d < nsd_; ++d) switch (s)
          {
            case 1:
              ele->eleinteriorElectricnm1_(d * shapes_.ndofs_ + r) = localMat(r, d);
              break;
            case 2:
              ele->eleinteriorElectricnm2_(d * shapes_.ndofs_ + r) = localMat(r, d);
              break;
            case 3:
              ele->eleinteriorElectricnm3_(d * shapes_.ndofs_ + r) = localMat(r, d);
              break;
          }
    }

  return 0;
}

/*----------------------------------------------------------------------*
 * ProjectElectricFieldFromScatra
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ProjectElectricFieldFromScatra(
    DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
    const Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseVector& elevec)
{
  shapes_.Evaluate(*ele);

  Teuchos::RCP<CORE::LINALG::SerialDenseVector> nodevals_phi =
      params.get<Teuchos::RCP<CORE::LINALG::SerialDenseVector>>("nodevals_phi");

  if (params.isParameter("ishdg"))
  {
    unsigned int scatra_ndofs = params.get<unsigned int>("ndofs");
    if (scatra_ndofs != shapes_.ndofs_)
      dserror("Scatra ndofs does not match elemag ndofs. This is not yet implemented.");

    // This is necessary as the scatra code computes the gradient of the scalar field multiplied by
    // the inverse of the conductivity. Therefore to obtain the actual value of the electric field
    // it is necessary to rescale it by the value of conductivity.
    const MAT::ElectromagneticMat* elemagmat =
        static_cast<const MAT::ElectromagneticMat*>(mat.get());
    const double sigma = elemagmat->sigma(ele->Id());

    for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        ele->eleinteriorElectric_(i + d * shapes_.ndofs_) =
            -(*nodevals_phi)(i + (d + 1) * shapes_.ndofs_) / sigma;

    return 0;
  }

  if ((*nodevals_phi).length() != nen_) dserror("node number not matching");

  CORE::LINALG::SerialDenseMatrix localMat(shapes_.ndofs_, nsd_);
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    // Storing the values of the coordinates for the current quadrature point
    // and of the jacobian computed in that point
    const double fac = shapes_.jfac(q);
    CORE::LINALG::SerialDenseVector intVal(nsd_);
    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int n = 0; n < nen_; ++n)
        intVal(d) += shapes_.derxy(n * nsd_ + d, q) * (*nodevals_phi)[n];

    // now fill the components in the one-sided mass matrix and the right hand side
    for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
    {
      // Mass matrix
      massPart(i, q) = shapes_.shfunct(i, q);
      massPartW(i, q) = shapes_.shfunct(i, q) * fac;

      // RHS for the electric and magnetic field
      for (int j = 0; j < intVal.numRows(); ++j)
        localMat(i, j) += shapes_.shfunct(i, q) * intVal(j) * fac;
    }
  }
  // The integration is made by computing the matrix product
  CORE::LINALG::multiplyNT(massMat, massPart, massPartW);
  {
    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
    inverseMass.setMatrix(Teuchos::rcpFromRef(massMat));
    inverseMass.setVectors(Teuchos::rcpFromRef(localMat), Teuchos::rcpFromRef(localMat));
    inverseMass.solve();
  }

  for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
    for (unsigned int d = 0; d < nsd_; ++d)
      ele->eleinteriorElectric_(i + d * shapes_.ndofs_) = -localMat(i, d);

  return 0;
}

/*----------------------------------------------------------------------*
 * ComputeError
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ComputeError(
    DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector& elevec1)
{
  double error_ele = 0.0, error_ele_post = 0.0, error_mag = 0.0;
  double exact_ele = 0.0, exact_mag = 0.0;
  CORE::LINALG::SerialDenseVector error_ele_grad(2), error_mag_grad(2), error_ele_post_grad(2);
  shapes_.Evaluate(*ele);
  postproc_shapes_.Evaluate(*ele);

  CORE::DRT::UTILS::ShapeValues<distype> highshapes(
      ele->Degree(), shapes_.usescompletepoly_, (ele->Degree() + 2) * 2);
  highshapes.Evaluate(*ele);

  CORE::DRT::UTILS::ShapeValues<distype> highshapes_post(
      ele->Degree() + 1, shapes_.usescompletepoly_, (ele->Degree() + 3) * 2);
  highshapes_post.Evaluate(*ele);
  // get function
  const int func = params.get<int>("funcno");
  const double time = params.get<double>("time");
  // for the calculation of the error, we use a higher integration rule
  CORE::LINALG::Matrix<nsd_, 1> xsi;
  CORE::LINALG::SerialDenseVector electric(nsd_);
  CORE::LINALG::SerialDenseVector electric_post(nsd_);
  CORE::LINALG::SerialDenseMatrix electric_grad(nsd_, nsd_);
  CORE::LINALG::SerialDenseMatrix electric_post_grad(nsd_, nsd_);
  CORE::LINALG::SerialDenseVector magnetic(nsd_);
  CORE::LINALG::SerialDenseMatrix magnetic_grad(nsd_, nsd_);
  CORE::LINALG::SerialDenseVector analytical(2 * nsd_);
  CORE::LINALG::SerialDenseMatrix analytical_grad(2 * nsd_, nsd_);

  for (unsigned int q = 0; q < highshapes.nqpoints_; ++q)
  {
    // Zero all temp vectors
    electric.putScalar(0.0), magnetic.putScalar(0.0);
    electric_grad.putScalar(0.0), magnetic_grad.putScalar(0.0);
    analytical.putScalar(0.0), analytical_grad.putScalar(0.0);

    for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        electric(d) += highshapes.shfunct(i, q) * ele->eleinteriorElectric_(d * shapes_.ndofs_ + i);
        magnetic(d) += highshapes.shfunct(i, q) * ele->eleinteriorMagnetic_(d * shapes_.ndofs_ + i);
        for (unsigned int d_grad = 0; d_grad < nsd_; ++d_grad)
        {
          electric_grad(d, d_grad) += highshapes.shderxy(i * nsd_ + d_grad, q) *
                                      ele->eleinteriorElectric_(d * shapes_.ndofs_ + i);
          magnetic_grad(d, d_grad) += highshapes.shderxy(i * nsd_ + d_grad, q) *
                                      ele->eleinteriorMagnetic_(d * shapes_.ndofs_ + i);
        }
      }

    // Evaluate error function and its derivatives in the integration point (real) coordinates
    for (unsigned int idim = 0; idim < nsd_; idim++) xsi(idim) = highshapes.xyzreal(idim, q);
    EvaluateAll(func, time, xsi, analytical);
    ComputeFunctionGradient(func, time, xsi, analytical_grad);

    for (unsigned int d = 0; d < nsd_; ++d)
    {
      // Electric error
      error_ele += std::pow((analytical(d) - electric(d)), 2) * highshapes.jfac(q);
      exact_ele += std::pow(analytical(d), 2) * highshapes.jfac(q);
      // Magnetic error
      error_mag += std::pow((analytical(d + nsd_) - magnetic(d)), 2) * highshapes.jfac(q);
      exact_mag += std::pow(analytical(d + nsd_), 2) * highshapes.jfac(q);
      // Divergence
      error_ele_grad(0) +=
          std::pow(analytical_grad(d, d) - electric_grad(d, d), 2) * highshapes.jfac(q);
      error_mag_grad(0) +=
          std::pow(analytical_grad(d + nsd_, d) - magnetic_grad(d, d), 2) * highshapes.jfac(q);
      // Rotor
      // Compute rotor components
      const double analytical_electric_rot = analytical_grad((d + 2) % nsd_, (d + 1) % nsd_) -
                                             analytical_grad((d + 1) % nsd_, (d + 2) % nsd_);
      const double analytical_magnetic_rot =
          analytical_grad((d + 2) % nsd_ + nsd_, (d + 1) % nsd_) -
          analytical_grad((d + 1) % nsd_ + nsd_, (d + 2) % nsd_);
      const double electric_rot = electric_grad((d + 2) % nsd_, (d + 1) % nsd_) -
                                  electric_grad((d + 1) % nsd_, (d + 2) % nsd_);
      const double magnetic_rot = magnetic_grad((d + 2) % nsd_, (d + 1) % nsd_) -
                                  magnetic_grad((d + 1) % nsd_, (d + 2) % nsd_);
      // Use rotor components
      error_ele_grad(1) += std::pow(analytical_electric_rot - electric_rot, 2) * highshapes.jfac(q);
      error_mag_grad(1) += std::pow(analytical_magnetic_rot - magnetic_rot, 2) * highshapes.jfac(q);
    }
  }

  if (params.get<bool>("postprocess"))
  {
    // Post-processed quantities
    for (unsigned int q = 0; q < highshapes_post.nqpoints_; ++q)
    {
      // Zero all temp vectors
      electric_post.putScalar(0.0), electric_post_grad.putScalar(0.0);
      analytical.putScalar(0.0), analytical_grad.putScalar(0.0);

      for (unsigned int i = 0; i < highshapes_post.ndofs_; ++i)
        for (unsigned int d = 0; d < nsd_; ++d)
        {
          electric_post(d) += highshapes_post.shfunct(i, q) *
                              ele->eleinteriorElectricPost_(d * highshapes_post.ndofs_ + i);
          for (unsigned int d_grad = 0; d_grad < nsd_; ++d_grad)
          {
            electric_post_grad(d, d_grad) +=
                highshapes_post.shderxy(i * nsd_ + d_grad, q) *
                ele->eleinteriorElectricPost_(d * highshapes_post.ndofs_ + i);
          }
        }

      // Evaluate error function and its derivatives in the integration point (real) coordinates
      for (unsigned int idim = 0; idim < nsd_; idim++) xsi(idim) = highshapes_post.xyzreal(idim, q);
      EvaluateAll(func, time, xsi, analytical);
      ComputeFunctionGradient(func, time, xsi, analytical_grad);

      for (unsigned int d = 0; d < nsd_; ++d)
      {
        // Electric error
        error_ele_post += std::pow((analytical(d) - electric_post(d)), 2) * highshapes_post.jfac(q);
        exact_ele += std::pow(analytical(d), 2) * highshapes_post.jfac(q);
        // Divergence
        error_ele_post_grad(0) +=
            std::pow(analytical_grad(d, d) - electric_post_grad(d, d), 2) * highshapes_post.jfac(q);
        // Rotor
        // Compute rotor components
        const double analytical_electric_rot = analytical_grad((d + 2) % nsd_, (d + 1) % nsd_) -
                                               analytical_grad((d + 1) % nsd_, (d + 2) % nsd_);
        const double electric_post_rot = electric_post_grad((d + 2) % nsd_, (d + 1) % nsd_) -
                                         electric_post_grad((d + 1) % nsd_, (d + 2) % nsd_);
        // Use rotor components
        error_ele_post_grad(1) +=
            std::pow(analytical_electric_rot - electric_post_rot, 2) * highshapes_post.jfac(q);
      }
    }
  }

  // Electric error
  elevec1[0] = error_ele;
  // Electric analytical reference
  elevec1[1] = exact_ele;
  // Magnetic error
  elevec1[2] = error_mag;
  // Magnetic analytical reference
  elevec1[3] = exact_mag;
  // Electric Hdiv error
  elevec1[4] = error_ele + error_ele_grad(0);
  // Electric Hrot error
  elevec1[5] = error_ele + error_ele_grad(1);
  // Magnetic Hdiv error
  elevec1[6] = error_mag + error_mag_grad(0);
  //// Magnetic Hrot error
  elevec1[7] = error_mag + error_mag_grad(1);
  // Electric error (post)
  elevec1[8] = error_ele_post;
  // Electric Hdiv error (post)
  elevec1[9] = error_ele_post + error_ele_post_grad(0);
  // Electric Hrot error (post)
  elevec1[10] = error_ele_post + error_ele_post_grad(1);

  return;
}

/*----------------------------------------------------------------------*
 * PostProcessing
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::PostProcessing(
    DRT::ELEMENTS::ElemagDiff& ele)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagDiffEleCalc::PostProcessing");
  shapes_.Evaluate(ele);
  postproc_shapes_.Evaluate(ele);

  CORE::DRT::UTILS::ShapeValues<distype> k2_shapes(ele.Degree() + 2, false, 2 * (ele.Degree() + 2));
  k2_shapes.Evaluate(ele);

  // Build matrix
  CORE::LINALG::SerialDenseMatrix postproc_mat(
      postproc_shapes_.ndofs_ * nsd_ + k2_shapes.ndofs_ * nsd_, postproc_shapes_.ndofs_ * nsd_);
  CORE::LINALG::SerialDenseVector postproc_rhs(
      postproc_shapes_.ndofs_ * nsd_ + k2_shapes.ndofs_ * nsd_);
  for (unsigned int q = 0; q < postproc_shapes_.nqpoints_; ++q)
    for (unsigned int i = 0; i < postproc_shapes_.ndofs_; ++i)
      for (unsigned int j = 0; j < postproc_shapes_.ndofs_; ++j)
        for (unsigned int d = 0; d < nsd_; ++d)
        {
          postproc_mat(i + d * postproc_shapes_.ndofs_, j + d * postproc_shapes_.ndofs_) +=
              (postproc_shapes_.shderxy(i * nsd_ + ((d + 2) % nsd_), q) *
                      postproc_shapes_.shderxy(j * nsd_ + ((d + 2) % nsd_), q) +
                  postproc_shapes_.shderxy(i * nsd_ + ((d + 1) % nsd_), q) *
                      postproc_shapes_.shderxy(j * nsd_ + ((d + 1) % nsd_), q)) *
              postproc_shapes_.jfac(q);
          postproc_mat(
              i + d * postproc_shapes_.ndofs_, j + ((d + 1) % nsd_) * postproc_shapes_.ndofs_) -=
              postproc_shapes_.shderxy(i * nsd_ + ((d + 1) % nsd_), q) *
              postproc_shapes_.shderxy(j * nsd_ + d, q) * postproc_shapes_.jfac(q);
          postproc_mat(
              i + d * postproc_shapes_.ndofs_, j + ((d + 2) % nsd_) * postproc_shapes_.ndofs_) -=
              postproc_shapes_.shderxy(i * nsd_ + ((d + 2) % nsd_), q) *
              postproc_shapes_.shderxy(j * nsd_ + d, q) * postproc_shapes_.jfac(q);
        }

  for (unsigned int q = 0; q < k2_shapes.nqpoints_; ++q)
  {
    CORE::LINALG::Matrix<nsd_, 1> xsi;
    const double* gpcoord = k2_shapes.quadrature_->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++) xsi(idim) = gpcoord[idim];
    CORE::LINALG::SerialDenseVector values(postproc_shapes_.ndofs_);
    postproc_shapes_.polySpace_->Evaluate(xsi, values);
    for (unsigned int i = 0; i < k2_shapes.ndofs_; ++i)
      for (unsigned int j = 0; j < postproc_shapes_.ndofs_; ++j)
        for (unsigned int d = 0; d < nsd_; ++d)
          postproc_mat(i + d * k2_shapes.ndofs_ + nsd_ * postproc_shapes_.ndofs_,
              j + d * postproc_shapes_.ndofs_) +=
              k2_shapes.shderxy(i * nsd_ + d, q) * values(j) * k2_shapes.jfac(q);
  }

  // Build RHS
  for (unsigned int q = 0; q < postproc_shapes_.nqpoints_; ++q)
  {
    CORE::LINALG::Matrix<nsd_, 1> xsi;
    const double* gpcoord = postproc_shapes_.quadrature_->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++) xsi(idim) = gpcoord[idim];
    CORE::LINALG::SerialDenseVector values(shapes_.ndofs_);
    shapes_.polySpace_->Evaluate(xsi, values);
    for (unsigned int i = 0; i < postproc_shapes_.ndofs_; ++i)
      for (unsigned int j = 0; j < shapes_.ndofs_; ++j)
        for (unsigned int d = 0; d < nsd_; ++d)
        {
          postproc_rhs(i + d * postproc_shapes_.ndofs_) +=
              postproc_shapes_.shderxy(i * nsd_ + ((d + 2) % nsd_), q) * values(j) *
              ele.eleinteriorMagnetic_(j + ((d + 1) % nsd_) * shapes_.ndofs_) *
              postproc_shapes_.jfac(q);
          postproc_rhs(i + d * postproc_shapes_.ndofs_) -=
              postproc_shapes_.shderxy(i * nsd_ + ((d + 1) % nsd_), q) * values(j) *
              ele.eleinteriorMagnetic_(j + ((d + 2) % nsd_) * shapes_.ndofs_) *
              postproc_shapes_.jfac(q);
        }
  }

  for (unsigned int q = 0; q < k2_shapes.nqpoints_; ++q)
  {
    CORE::LINALG::Matrix<nsd_, 1> xsi;
    const double* gpcoord = k2_shapes.quadrature_->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++) xsi(idim) = gpcoord[idim];
    CORE::LINALG::SerialDenseVector values(shapes_.ndofs_);
    shapes_.polySpace_->Evaluate(xsi, values);
    for (unsigned int i = 0; i < k2_shapes.ndofs_; ++i)
      for (unsigned int j = 0; j < shapes_.ndofs_; ++j)
        for (unsigned int d = 0; d < nsd_; ++d)
          postproc_rhs(i + d * k2_shapes.ndofs_ + nsd_ * postproc_shapes_.ndofs_) +=
              k2_shapes.shderxy(i * nsd_ + d, q) * values(j) *
              ele.eleinteriorElectric_(j + d * shapes_.ndofs_) * k2_shapes.jfac(q);
  }

  {
    CORE::LINALG::SerialDenseVector test(postproc_rhs.length() * 2);
    Teuchos::LAPACK<int, double> solve;
    int err;
    solve.GELS('N', postproc_mat.numRows(), postproc_mat.numCols(), 1, postproc_mat.values(),
        postproc_mat.stride(), postproc_rhs.values(), postproc_rhs.length(), test.values(),
        test.length(), &err);
    if (err != 0)
      dserror("Least-square approximation for the Postprocessing failed with error %d", err);

    for (int i = 0; i < ele.eleinteriorElectricPost_.length(); ++i)
      ele.eleinteriorElectricPost_(i) = postproc_rhs(i);
  }

  return;
}  // PostProcessing

/*----------------------------------------------------------------------*
 * ProjectFieldTest
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ProjectFieldTest(
    DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2)
{
  shapes_.Evaluate(*ele);
  postproc_shapes_.Evaluate(*ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.numRows() == 0 || unsigned(elevec2.numRows()) == nsd_ * shapes_.ndofs_,
      "Wrong size in project vector 2");

  // get function
  const int* start_func = params.getPtr<int>("startfuncno");
  const double time = params.get<double>("time");

  bool do_internal = false;
  bool do_postprocess = true;

  // internal variables
  if (do_internal)
  {
    // the RHS matrix has to have the row dimension equal to the number of shape
    // functions(so we have one coefficient for each) and a number of column
    // equal to the overall number of component that we want to solve for.
    // The number is nsd_*2 because we have two fields..
    CORE::LINALG::SerialDenseMatrix localMat(shapes_.ndofs_, nsd_ * 2);
    for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
    {
      // Storing the values of the coordinates for the current quadrature point
      // and of the jacobian computed in that point
      const double fac = shapes_.jfac(q);
      CORE::LINALG::Matrix<nsd_, 1> xyz;
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz(d) = shapes_.xyzreal(d, q);  // coordinates of quadrature point in real coordinates
      // Creating the temporary electric and magnetic field vector intVal
      // The vector is going to contain first the electric and then the magnetic
      // field such that the field will be initialized as first tree component
      // of the specified function as electric field, last three components as
      // magnetic field. If there is only one component all the components will
      // be initialized to the same value.
      CORE::LINALG::SerialDenseVector intVal(2 * nsd_);
      dsassert(start_func != nullptr, "funct not set for initial value");
      EvaluateAll(*start_func, time, xyz, intVal);
      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      {
        // Mass matrix
        massPart(i, q) = shapes_.shfunct(i, q);
        massPartW(i, q) = shapes_.shfunct(i, q) * fac;

        // RHS for the electric and magnetic field
        for (int j = 0; j < intVal.numRows(); ++j)
          localMat(i, j) += shapes_.shfunct(i, q) * intVal(j) * fac;
      }
    }
    // The integration is made by computing the matrix product
    CORE::LINALG::multiplyNT(massMat, massPart, massPartW);
    {
      using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
      using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
      Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
      inverseMass.setMatrix(Teuchos::rcpFromRef(massMat));
      inverseMass.setVectors(Teuchos::rcpFromRef(localMat), Teuchos::rcpFromRef(localMat));
      inverseMass.solve();
    }

    // Here we move the values from the temporary variable to the variable
    // contained in the element
    for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
    {
      // Now we are storing the variables by component, meaning that we save for
      // each component the value for each dof and then we move to the next component.
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        ele->eleinteriorElectric_(d * shapes_.ndofs_ + r) = localMat(r, d);  // Electric field
        ele->eleinteriorMagnetic_(d * shapes_.ndofs_ + r) = localMat(r, d + nsd_);  // magnetic
      }
    }
  }

  // postproc variables
  if (do_postprocess)
  {
    CORE::LINALG::SerialDenseMatrix localMat(postproc_shapes_.ndofs_, nsd_ * 2);
    for (unsigned int q = 0; q < postproc_shapes_.nqpoints_; ++q)
    {
      const double fac = postproc_shapes_.jfac(q);
      CORE::LINALG::Matrix<nsd_, 1> xyz;
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz(d) =
            postproc_shapes_.xyzreal(d, q);  // coordinates of quadrature point in real coordinates

      CORE::LINALG::SerialDenseVector intVal(nsd_);
      EvaluateAll(*start_func, time, xyz, intVal);

      for (unsigned int i = 0; i < postproc_shapes_.ndofs_; ++i)
      {
        // Mass matrix
        massPart(i, q) = postproc_shapes_.shfunct(i, q);
        massPartW(i, q) = postproc_shapes_.shfunct(i, q) * fac;

        // RHS for the electric and magnetic field
        for (int j = 0; j < intVal.numRows(); ++j)
          localMat(i, j) += postproc_shapes_.shfunct(i, q) * intVal(j) * fac;
      }
    }
    // The integration is made by computing the matrix product
    CORE::LINALG::multiplyNT(massMat, massPart, massPartW);
    {
      using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
      using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
      Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
      inverseMass.setMatrix(Teuchos::rcpFromRef(massMat));
      inverseMass.setVectors(Teuchos::rcpFromRef(localMat), Teuchos::rcpFromRef(localMat));
      inverseMass.solve();
    }

    // Here we move the values from the temporary variable to the variable
    // contained in the element
    for (unsigned int r = 0; r < postproc_shapes_.ndofs_; ++r)
    {
      // Now we are storing the variables by component, meaning that we save for
      // each component the value for each dof and then we move to the next component.
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        ele->eleinteriorElectricPost_(d * postproc_shapes_.ndofs_ + r) = localMat(r, d);
      }
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * ProjectFieldTestTrace
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ProjectFieldTestTrace(
    DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector& elevec1)
{
  // Here we have the projection of the field on the trace
  // mass is the mass matrix for the system to be solved
  // the dimension of the mass matrix is given by the number of shape functions
  CORE::LINALG::SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
  // TRaceVEC is the vector of the trace values
  // instead of being a vector it is a matrix so that we use the same matrix
  // to solve the projection problem on every component of the field
  CORE::LINALG::SerialDenseMatrix trVec(shapesface_->nfdofs_, nsd_);

  const int* start_func = params.getPtr<int>("startfuncno");
  const double time = params.get<double>("time");

  // Cycling through faces
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    // Updating face data
    shapesface_->EvaluateFace(*ele, f);

    // Initializing the matrices
    // It is necessary to create a matrix and a trVec for each face because the
    // dimensions of each face can differ from the previous one and the jacobian
    // contains the dimension of the face in it.
    mass.putScalar(0.0);
    trVec.putScalar(0.0);

    // Cycling through the quadrature points
    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      // For each quadrature point we have a vector containing the field
      // components and a vector containing the spatial coordinates of that point
      CORE::LINALG::SerialDenseVector trace(nsd_);
      CORE::LINALG::Matrix<nsd_, 1> xyz;

      // Temporary variable to store the jacobian of the face (contains the weigth)
      const double fac = shapesface_->jfac(q);
      // Coordinates of quadrature point in real coordinates from the face to
      // the temporary variable. It is just to make the code easier to handle
      for (unsigned int d = 0; d < nsd_; ++d) xyz(d) = shapesface_->xyzreal(d, q);

      // Evaluation of the function in the quadrature point being considered
      EvaluateAll(*start_func, time, xyz, trace);

      // Creating the mass matrix and the RHS vector
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // Mass matrix
        for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
          mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;

        // RHS
        for (unsigned int d = 0; d < nsd_; ++d)
          trVec(i, d) += shapesface_->shfunct(i, q) * trace(d) * fac;
      }
    }

    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
    inverseMass.setMatrix(Teuchos::rcpFromRef(mass));
    inverseMass.setVectors(Teuchos::rcpFromRef(trVec), Teuchos::rcpFromRef(trVec));
    inverseMass.solve();

    CORE::LINALG::SerialDenseVector tempVec(shapesface_->nfdofs_ * (nsd_));
    CORE::LINALG::SerialDenseVector faceVec(shapesface_->nfdofs_ * (nsd_ - 1));
    // Filling the vector of trace values
    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // remember that "f" is an iterator index and therefore we are
        // cycling through all the faces and all the entries of elevec1
        // except for the first one where we will put the pressure average
        tempVec(d * shapesface_->nfdofs_ + i) = trVec(i, d);
      }

    CORE::LINALG::SerialDenseMatrix transformatrix(
        (nsd_ - 1) * shapesface_->nfdofs_, nsd_ * shapesface_->nfdofs_);
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        for (unsigned int q = 0; q < nsd_ - 1; ++q)
          transformatrix(shapesface_->nfdofs_ * q + i, shapesface_->nfdofs_ * d + i) =
              shapesface_->tangent(d, q);

    CORE::LINALG::multiply(faceVec, transformatrix, tempVec);

    // Filling the vector of trace values
    for (unsigned int d = 0; d < nsd_ - 1; ++d)
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // remember that "f" is an iterator index and therefore we are
        // cycling through all the faces and all the entries of elevec1
        // except for the first one where we will put the pressure average
        elevec1(f * shapesface_->nfdofs_ * (nsd_ - 1) + d * shapesface_->nfdofs_ + i) =
            faceVec(d * shapesface_->nfdofs_ + i);
        ele->elenodeTrace2d_(f * shapesface_->nfdofs_ * (nsd_ - 1) + d * shapesface_->nfdofs_ + i) =
            faceVec(d * shapesface_->nfdofs_ + i);
      }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * ProjectDirichField
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ProjectDirichField(
    DRT::ELEMENTS::ElemagDiff* ele, Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseVector& elevec1)
{
  // Updating face data
  const int face = params.get<unsigned int>("faceconsider");
  // std::cout << face << std::endl;
  shapesface_->EvaluateFace(*ele, face);

  Teuchos::Array<int>* functno = params.getPtr<Teuchos::Array<int>>("funct");
  const double time = params.get<double>("time");

  // Here we have the projection of the field on the trace
  // mass is the mass matrix for the system to be solved
  // the dimension of the mass matrix is given by the number of shape functions
  CORE::LINALG::SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
  // TRaceVEC is the vector of the trace values
  // instead of being a vector it is a matrix so that we use the same matrix
  // to solve the projection problem on every component of the field
  CORE::LINALG::SerialDenseMatrix trVec(shapesface_->nfdofs_, nsd_);

  // Cycling through the quadrature points
  for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
  {
    // For each quadrature point we have a vector containing the field
    // components and a vector containing the spatial coordinates of that point
    CORE::LINALG::SerialDenseVector trace(nsd_);
    CORE::LINALG::Matrix<nsd_, 1> xyz;

    // Temporary variable to store the jacobian of the face (contains the weigth)
    const double fac = shapesface_->jfac(q);
    // Coordinates of quadrature point in real coordinates from the face to
    // the temporary variable. It is just to make the code easier to handle
    for (unsigned int d = 0; d < nsd_; ++d) xyz(d) = shapesface_->xyzreal(d, q);

    // Evaluation of the function in the quadrature point being considered
    EvaluateAll((*functno)[0], time, xyz, trace);

    // Creating the mass matrix and the RHS vector
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    {
      // Mass matrix
      for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
        mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;

      // RHS
      for (unsigned int d = 0; d < nsd_; ++d)
        trVec(i, d) += shapesface_->shfunct(i, q) * trace(d) * fac;
    }
  }

  using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
  using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
  inverseMass.setMatrix(Teuchos::rcpFromRef(mass));
  inverseMass.setVectors(Teuchos::rcpFromRef(trVec), Teuchos::rcpFromRef(trVec));
  inverseMass.solve();

  CORE::LINALG::SerialDenseVector tempVec(shapesface_->nfdofs_ * (nsd_));
  // Filling the vector of trace values
  for (unsigned int d = 0; d < nsd_; ++d)
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      tempVec(d * shapesface_->nfdofs_ + i) = trVec(i, d);

  CORE::LINALG::SerialDenseMatrix transformatrix(
      (nsd_ - 1) * shapesface_->nfdofs_, nsd_ * shapesface_->nfdofs_);
  for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int q = 0; q < nsd_ - 1; ++q)
        transformatrix(shapesface_->nfdofs_ * q + i, shapesface_->nfdofs_ * d + i) =
            shapesface_->tangent(d, q);

  CORE::LINALG::multiply(elevec1, transformatrix, tempVec);

  return 0;
}

/*----------------------------------------------------------------------*
 * EvaluateAll
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::EvaluateAll(const int start_func,
    const double t, const CORE::LINALG::Matrix<nsd_, 1>& xyz,
    CORE::LINALG::SerialDenseVector& v) const
{
  int numComp = DRT::Problem::Instance()
                    ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(start_func - 1)
                    .NumberComponents();

  // If the number is not recognised throw an error
  if (not(numComp == v.numRows() || numComp == 2 * v.numRows() || numComp == v.numRows() / 2 ||
          numComp == 1))
    dserror(
        "Supply ONE component for your function or NUMDIM, not anything else! With NUMDIM "
        "components the field will be initialized componentwise, if only one component is "
        "provided, every component of the field will be initialized with the same values.");

  // If there is on component for each entry of the vector use une for each
  // If the vector is half the number of the component only use the firt half
  // If the number of component is half of the vector, repeat the first half twice
  // If there is only one component always use it
  for (int d = 0; d < v.numRows(); ++d)
    v[d] = DRT::Problem::Instance()
               ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(start_func - 1)
               .Evaluate(xyz.A(), t, d % numComp);

  return;
}

/*----------------------------------------------------------------------*
 * ComputeFunctionGradient
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ComputeFunctionGradient(
    const int start_func, const double t, const CORE::LINALG::Matrix<nsd_, 1>& xyz,
    CORE::LINALG::SerialDenseMatrix& v) const
{
  int numComp = DRT::Problem::Instance()
                    ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(start_func - 1)
                    .NumberComponents();

  // If the number is not recognised throw an error
  if (not(numComp == v.numRows() || numComp == 2 * v.numRows() || numComp == v.numRows() / 2 ||
          numComp == 1))
    dserror(
        "Supply ONE component for your function or NUMDIM, not anything else! With NUMDIM "
        "components the field will be initialized componentwise, if only one component is "
        "provided, every component of the field will be initialized with the same values.");
  // If there is on component for each entry of the vector use une for each
  // If the vector is half the number of the component only use the firt half
  // If the number of component is half of the vector, repeat the first half twice
  // If there is only one component always use it
  for (int d = 0; d < v.numRows(); ++d)
  {
    std::vector<double> deriv = DRT::Problem::Instance()
                                    ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(start_func - 1)
                                    .EvaluateSpatialDerivative(xyz.A(), t, d % numComp);
    for (unsigned int d_der = 0; d_der < nsd_; ++d_der) v(d, d_der) = deriv[d_der];
  }

  return;
}

/*----------------------------------------------------------------------*
 * ComputeFunctionTimeDerivative
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ComputeFunctionTimeDerivative(
    const int start_func, const double t, const double dt, const CORE::LINALG::Matrix<nsd_, 1>& xyz,
    CORE::LINALG::SerialDenseVector& v) const
{
  int numComp = DRT::Problem::Instance()
                    ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(start_func - 1)
                    .NumberComponents();

  // If the number is not recognised throw an error
  if (not(numComp == v.numRows() || numComp == 2 * v.numRows() || numComp == v.numRows() / 2 ||
          numComp == 1))
    dserror(
        "Supply ONE component for your start function or NUMDIM, not anything else! With NUMDIM "
        "components the field will be initialized componentwise, if only one component is "
        "provided, every component of the field will be initialized with the same values.");

  // If there is on component for each entry of the vector use one for each
  // If the vector is half the number of the component only use the firt half
  // If the number of component is half of the vector, repeat the first half twice
  // If there is only one component always use it
  for (int d = 0; d < v.numRows(); ++d)
    v[d] = (DRT::Problem::Instance()
                   ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(start_func - 1)
                   .Evaluate(xyz.A(), t + (0.5 * dt), d % numComp) -
               DRT::Problem::Instance()
                   ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(start_func - 1)
                   .Evaluate(xyz.A(), t - (0.5 * dt), d % numComp)) /
           dt;

  return;
}

/*----------------------------------------------------------------------*
 | InterpolateSolutionToNodes                          berardocco 04/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ElemagDiffEleCalc<distype>::InterpolateSolutionToNodes(
    DRT::ELEMENTS::ElemagDiff* ele, DRT::Discretization& discretization,
    CORE::LINALG::SerialDenseVector& elevec1)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagDiffEleCalc::InterpolateSolutionToNodes");
  InitializeShapes(ele);

  // Check if the vector has the correct size
  dsassert(elevec1.numRows() == (int)nen_ * (4 * nsd_), "Vector does not have correct size");

  // Getting the connectivity matrix
  // Contains the (local) coordinates of the nodes belonging to the element
  CORE::LINALG::SerialDenseMatrix locations =
      CORE::DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  // This vector will contain the values of the shape functions computed in a
  // certain coordinate. In fact the lenght of the vector is given by the number
  // of shape functions, that is the same of the number of degrees of freedom of
  // an element.
  CORE::LINALG::SerialDenseVector values(shapes_->ndofs_);
  CORE::LINALG::SerialDenseVector post_values(postproc_shapes_->ndofs_);

  // EVALUATE SHAPE POLYNOMIALS IN NODE
  // In hdg we can have several more points inside the element than in the
  //"real" discretization and therefore it is necessary to compute the value
  // that the internal solution takes in the node of the discretization.

  // Cycling through all the "real" nodes of the element to get the coordinates
  // Remember that the coordinates are the local ones.
  for (unsigned int i = 0; i < nen_; ++i)
  {
    // Cycling through the spatial dimensions to get the coordinates
    for (unsigned int idim = 0; idim < nsd_; idim++) shapes_->xsi(idim) = locations(idim, i);
    // Evaluating the polinomials in the point given by "shapes_->xsi".
    // The polynomials are the internal ones.
    // The result of the evaluation is given in "values".
    shapes_->polySpace_->Evaluate(shapes_->xsi, values);
    postproc_shapes_->polySpace_->Evaluate(shapes_->xsi, post_values);

    // compute values for interior unknown by summing over all basis functions
    for (unsigned int d = 0; d < nsd_; ++d)
    {
      double sum_electric = 0.0;
      double sum_electric_post = 0.0;
      double sum_magnetic = 0.0;
      // Cycling through all the shape functions
      for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
      {
        // The overall value in the chosen point is given by the sum of the
        // values of the shape functions multiplied by their coefficients.
        sum_electric += values(k) * ele->eleinteriorElectric_[d * shapes_->ndofs_ + k];
        sum_magnetic += values(k) * ele->eleinteriorMagnetic_[d * shapes_->ndofs_ + k];
      }
      for (unsigned int k = 0; k < postproc_shapes_->ndofs_; ++k)
        sum_electric_post +=
            post_values(k) * ele->eleinteriorElectricPost_[d * postproc_shapes_->ndofs_ + k];
      // sum contains the linear combination of the shape functions times the
      // coefficients and its values are reordered in elevec1 grouped by
      // component: the first component for every node, then the following
      // component for the same nodes and so on for every component.
      elevec1(d * nen_ + i) = sum_electric;
      elevec1(nen_ * nsd_ + d * nen_ + i) = sum_magnetic;
      elevec1(nen_ * 2 * nsd_ + d * nen_ + i) = sum_electric_post;
    }
  }

  // get trace solution values
  // Same as before bu this time the dimension is nsd_-1 because we went from
  // the interior to the faces. We have to be careful because we are using a
  // part of the previous vector. The coordinates are still in the local frame.
  locations = CORE::DRT::UTILS::getEleNodeNumbering_nodes_paramspace(
      CORE::DRT::UTILS::DisTypeToFaceShapeType<distype>::shape);

  // Storing the number of nodes for each face of the element as vector
  // NumberCornerNodes
  std::vector<int> ncn = CORE::DRT::UTILS::getNumberOfFaceElementCornerNodes(distype);
  // NumberInternalNodes
  std::vector<int> nin = CORE::DRT::UTILS::getNumberOfFaceElementInternalNodes(distype);

  // Cycling the faces of the element
  CORE::LINALG::SerialDenseVector fvalues(shapesface_->nfdofs_);
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    // Checking how many nodes the face has
    const int nfn = CORE::DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace;

    shapesface_->EvaluateFace(*ele, f);

    CORE::LINALG::SerialDenseVector facetrace((nsd_ - 1) * shapesface_->nfdofs_);
    CORE::LINALG::SerialDenseVector temptrace(nsd_ * shapesface_->nfdofs_);

    // The dimension of the coordinate matrix is now nsd_ times the number of nodes in the face.
    CORE::LINALG::Matrix<nsd_ - 1, nfn> xsishuffle(true);

    // Cycling throught the nodes of the face to store the node positions in the
    // correct order using xsishuffle as a temporary vector
    for (int i = 0; i < nfn; ++i)
    {
      // cycling through the spatial dimensions
      for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
      {
        // If the face belongs to the element being considered
        if (ele->Faces()[f]->ParentMasterElement() == ele)
          xsishuffle(idim, i) = locations(idim, i);
        else
          // If the face does not belong to the element being considered it is
          // necessary to change the ordering
          xsishuffle(idim, ele->Faces()[f]->GetLocalTrafoMap()[i]) = locations(idim, i);
      }
    }

    // Transformation for the face reference system
    CORE::LINALG::SerialDenseMatrix transformatrix(
        (nsd_ - 1) * shapesface_->nfdofs_, nsd_ * shapesface_->nfdofs_);
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        for (unsigned int q = 0; q < nsd_ - 1; ++q)
          transformatrix(shapesface_->nfdofs_ * q + i, shapesface_->nfdofs_ * d + i) =
              shapesface_->tangent(d, q);

    // Storing the face part of the trace vector
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      for (unsigned int d = 0; d < nsd_ - 1; ++d)
        facetrace(shapesface_->nfdofs_ * d + i) =
            ele->elenodeTrace2d_[f * (nsd_ - 1) * shapesface_->nfdofs_ + shapesface_->nfdofs_ * d +
                                 i];

    CORE::LINALG::multiplyTN(temptrace, transformatrix, facetrace);

    // EVALUATE SHAPE POLYNOMIALS IN NODE
    // Now that we have an ordered coordinates vector we can easily compute the
    // values of the shape functions in the nodes.
    for (int i = 0; i < nfn; ++i)
    {
      // Storing the actual coordinates of the current node
      for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
        shapesface_->xsi(idim) = xsishuffle(idim, i);

      // Actually evaluating shape polynomials in node
      shapesface_->polySpace_->Evaluate(shapesface_->xsi, fvalues);

      // compute values for trace vector by summing over the shape functions
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        double sum = 0;
        // Linear combination of the values of the shape functions and
        // relative weighting coefficients. The weighting coefficients are
        // given by the value of the unknowns in the nodes.
        for (unsigned int k = 0; k < shapesface_->nfdofs_; ++k)
          sum += fvalues(k) * temptrace[d * shapesface_->nfdofs_ + k];
        // Ordering the results of the interpolation in the vector being careful
        // about the ordering of the nodes in the faces.
        if (i < ncn[f])
        {
          elevec1((nsd_ * 3 + d) * nen_ + shapesface_->faceNodeOrder[f][i]) += sum / nsd_;
        }
        else if (i < nfn - nin[f])
        {
          elevec1((nsd_ * 3 + d) * nen_ + shapesface_->faceNodeOrder[f][i]) += sum / (nsd_ - 1);
        }
        else
        {
          elevec1((nsd_ * 3 + d) * nen_ + shapesface_->faceNodeOrder[f][i]) += sum;
        }
      }
    }
  }
  return 0;
}

template <CORE::FE::CellType distype>
DRT::ELEMENTS::ElemagDiffEleCalc<distype>* DRT::ELEMENTS::ElemagDiffEleCalc<distype>::Instance(
    CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<DRT::ELEMENTS::ElemagDiffEleCalc<distype>>(
            new DRT::ELEMENTS::ElemagDiffEleCalc<distype>());
      });

  return singleton_owner.Instance(action);
}

/*----------------------------------------------------------------------*
 * Constructor LocalSolver
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::LocalSolver(
    const DRT::ELEMENTS::ElemagDiff* ele, CORE::DRT::UTILS::ShapeValues<distype>& shapeValues,
    Teuchos::RCP<CORE::DRT::UTILS::ShapeValuesFace<distype>>& shapeValuesFace,
    INPAR::ELEMAG::DynamicType& dyna, CORE::DRT::UTILS::ShapeValues<distype>& postproc_shapeValues)
    : ndofs_(shapeValues.ndofs_),
      shapes_(shapeValues),
      shapesface_(shapeValuesFace),
      postproc_shapes_(postproc_shapeValues),
      dyna_(dyna)
{
  // shape all matrices
  // Each one of these matrices is related to one equation of the formulation,
  // therefore ndofs equations in FEM terms) and one variable.
  // The number of entries is then given by ndofs time sthe dimension of the
  // space where the unknown lies. For vectorial field nsd_ gives the dimension.
  Amat.shape(nsd_ * ndofs_, nsd_ * ndofs_);
  invAmat.shape(nsd_ * ndofs_, nsd_ * ndofs_);
  Bmat.shape(nsd_ * ndofs_, nsd_ * ndofs_);
  Dmat.shape(nsd_ * ndofs_, nsd_ * ndofs_);
  Emat.shape(nsd_ * ndofs_, nsd_ * ndofs_);
  Gmat.shape(nsd_ * ndofs_, nsd_ * ndofs_);
  // These matrices have a "strange" shape because to merge them there will be
  // applied a matrix multiplication between the first one and the transposed
  // second one. The shape of the resulting matrix will therefore be ndofs x ndofs.
  massMat.shape(ndofs_, ndofs_);
  massPart.shape(ndofs_, shapeValues.nqpoints_);
  massPartW.shape(ndofs_, shapeValues.nqpoints_);

  // Matrix compriending the hybrid variable or the continuity condition
  // It is necessary to compute the overall number of degrees
  // of freedom by summing the number of degrees of freedom on every face
  // surrounding the volume of the element.
  // ONFaceDegreesOfFreedomS
  int onfdofs = 0;
  for (unsigned int i = 0; i < nfaces_; ++i)
  {
    // Evaluating the dofs number on each face of the element
    shapesface_->EvaluateFace(*ele, i);
    // Computing the dimension of the approximation space for the hybrid variable
    onfdofs += shapesface_->nfdofs_;
  }
  // The hybrid variable is vectorial and therefore the dimension of the space
  // has to be multiplied by nsd_.
  onfdofs = onfdofs * (nsd_ - 1);

  // This part is specially dependent on the formulation being used, in fact,
  // when the matrices relative to the surface integrals have to be created
  // those will have different dimensions depending on the variable that appears
  // in the integral itself. The hybrid variable is defined in the trace space
  // and therefore its shape functions belong to the same space,
  // with the consequence of being onfdof shape functions.
  // Cmat and Hmat are the matrix that belongs to the equation for u
  // and electric field but multiply the hybrid variable, therefore their dimensions are:
  // o) nsd_*ndofs_ x onfdofs
  Cmat.shape(nsd_ * ndofs_, onfdofs);
  Hmat.shape(nsd_ * ndofs_, onfdofs);
  // Finally Lmat is the matrix that belongs to the continuity condition and
  // multiplies the hybrid variable and therefore its dimensions are:
  // o) ondofs x ondofs
  Lmat.shape(onfdofs, onfdofs);
}

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidual
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::UpdateInteriorVariablesAndComputeResidual(
    Teuchos::ParameterList& params, DRT::ELEMENTS::ElemagDiff& ele,
    const Teuchos::RCP<MAT::Material>& mat, CORE::LINALG::SerialDenseVector& elevec, double dt,
    bool errormaps, bool updateonly)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "DRT::ELEMENTS::ElemagDiffEleCalc::UpdateInteriorVariablesAndComputeResidual");

  // *****************************************************
  // update interior variables first
  // *****************************************************

  CORE::LINALG::SerialDenseVector tempVec1(shapes_->ndofs_ * nsd_);
  CORE::LINALG::SerialDenseVector tempVec2(shapes_->ndofs_ * nsd_);
  CORE::LINALG::SerialDenseVector tempVec3(shapes_->ndofs_ * nsd_);
  CORE::LINALG::SerialDenseVector xVec(shapes_->ndofs_ * nsd_);
  CORE::LINALG::SerialDenseVector yVec(shapes_->ndofs_ * nsd_);
  CORE::LINALG::SerialDenseMatrix tempMat(shapes_->ndofs_ * nsd_, shapes_->ndofs_ * nsd_);
  CORE::LINALG::SerialDenseMatrix tempMat2(shapes_->ndofs_ * nsd_, shapes_->ndofs_ * nsd_);

  // The source has to be checked for bdf
  localSolver_->ComputeSource(params, tempVec2, xVec);

  CORE::LINALG::multiplyTN(tempMat, localSolver_->Bmat, localSolver_->invAmat);  // FA^{-1}

  tempMat2 += localSolver_->Emat;
  tempMat2 += localSolver_->Gmat;
  // Only if the D matrix is not zero <-> epsilon != 0
  tempMat2 += localSolver_->Dmat;
  CORE::LINALG::multiply(1.0, tempMat2, -1.0, tempMat, localSolver_->Bmat);  //(E + G) - FA^{-1}B
  {
    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> invert;
    invert.setMatrix(Teuchos::rcpFromRef(tempMat2));
    invert.invert();  //  [(E + G) - FA^{-1}B]^{-1}
  }

  if (dyna_ == INPAR::ELEMAG::elemag_bdf2)
  {
    CORE::LINALG::multiply(-1.0, tempVec2, -1.0 / 3.0, localSolver_->Emat,
        ele.eleinteriorElectricnm1_);  // (1/3)EE^{n}
    CORE::LINALG::multiply(1.0, tempVec2, 4.0 / 3.0, localSolver_->Emat,
        ele.eleinteriorElectric_);  // ^E = (4/3)EE^{n+1} - (1/3)EE^{n}
                                    // Only if the D matrix is not zero <-> epsilon != 0
    CORE::LINALG::multiply(
        1.0, tempVec2, 1.0 / 2.0, localSolver_->Dmat, ele.eleinteriorElectricnm2_);
    CORE::LINALG::multiply(1.0, tempVec2, -2.0, localSolver_->Dmat, ele.eleinteriorElectricnm1_);
    CORE::LINALG::multiply(1.0, tempVec2, 5.0 / 2.0, localSolver_->Dmat, ele.eleinteriorElectric_);
  }
  else if (dyna_ == INPAR::ELEMAG::elemag_bdf4)
  {
    CORE::LINALG::multiply(-1.0, tempVec2, -3.0 / 25.0, localSolver_->Emat,
        ele.eleinteriorElectricnm3_);  // (1/3)E E^{n} + I_s
    CORE::LINALG::multiply(1.0, tempVec2, 16.0 / 25.0, localSolver_->Emat,
        ele.eleinteriorElectricnm2_);  // ^E = (4/3)EE^{n+1} - (1/3)EE^{n} - I_s
    CORE::LINALG::multiply(
        1.0, tempVec2, -36.0 / 25.0, localSolver_->Emat, ele.eleinteriorElectricnm1_);
    CORE::LINALG::multiply(
        1.0, tempVec2, 48.0 / 25.0, localSolver_->Emat, ele.eleinteriorElectric_);
  }
  else
  {
    CORE::LINALG::multiply(-1.0, tempVec2, 1.0, localSolver_->Emat,
        ele.eleinteriorElectric_);  // EE - I_s Implicit euler
                                    // Only if the D matrix is not zero <-> epsilon != 0
    CORE::LINALG::multiply(1.0, tempVec2, -1.0, localSolver_->Dmat, ele.eleinteriorElectricnm1_);
    CORE::LINALG::multiply(1.0, tempVec2, 2.0, localSolver_->Dmat, ele.eleinteriorElectric_);
  }
  ele.eleinteriorElectricnm3_ = ele.eleinteriorElectricnm2_;
  ele.eleinteriorElectricnm2_ = ele.eleinteriorElectricnm1_;
  ele.eleinteriorElectricnm1_ = ele.eleinteriorElectric_;

  // C\lambda^{n+2}
  CORE::LINALG::multiply(tempVec1, localSolver_->Cmat, ele.elenodeTrace2d_);
  CORE::LINALG::multiply(
      1.0, tempVec2, -1.0, localSolver_->Hmat, ele.elenodeTrace2d_);  // ^E - I_s - H\lambda^{n+2}
  CORE::LINALG::multiply(1.0, tempVec2, 1.0, tempMat, tempVec1);      //  ^E + FA^{-1}C\lambda^{n+2}

  //  E^{n+2} = [(E + G) - FA^{-1}B]^{-1} (^E + (FA^{-1}C - H)\lambda^{n+2})
  CORE::LINALG::multiply(ele.eleinteriorElectric_, tempMat2, tempVec2);

  // C\lambda^{n+2} + BE^{n}
  CORE::LINALG::multiply(1.0, tempVec1, 1.0, localSolver_->Bmat, ele.eleinteriorElectric_);

  //  = -A^{-1}(C\lambda^{n+2} + BE^{n})
  CORE::LINALG::multiply(0.0, ele.eleinteriorMagnetic_, -1.0, localSolver_->invAmat, tempVec1);

  // Updateresidual

  if (dyna_ == INPAR::ELEMAG::elemag_bdf2)
  {
    //  = -1/3EE^{n+2} - I_s
    CORE::LINALG::multiply(-1.0, xVec, -1.0 / 3.0, localSolver_->Emat, ele.eleinteriorElectricnm1_);
    ////  = ^E - I_s = 4/3EE^{n} - 1/3EE^{n+2} - I_s
    CORE::LINALG::multiply(1.0, xVec, 4.0 / 3.0, localSolver_->Emat, ele.eleinteriorElectric_);
    // Only if the D matrix is not zero <-> epsilon != 0
    CORE::LINALG::multiply(1.0, xVec, 1.0 / 2.0, localSolver_->Dmat, ele.eleinteriorElectricnm2_);
    CORE::LINALG::multiply(1.0, xVec, -2.0, localSolver_->Dmat, ele.eleinteriorElectricnm1_);
    CORE::LINALG::multiply(1.0, xVec, 5.0 / 2.0, localSolver_->Dmat, ele.eleinteriorElectric_);
  }
  else if (dyna_ == INPAR::ELEMAG::elemag_bdf4)
  {
    CORE::LINALG::multiply(-1.0, xVec, -3.0 / 25.0, localSolver_->Emat,
        ele.eleinteriorElectricnm3_);  // (1/3)E E^{n} + I_s
    CORE::LINALG::multiply(1.0, xVec, 16.0 / 25.0, localSolver_->Emat,
        ele.eleinteriorElectricnm2_);  // ^E = (4/3)EE^{n+1} - (1/3)EE^{n} - I_s
    CORE::LINALG::multiply(
        1.0, xVec, -36.0 / 25.0, localSolver_->Emat, ele.eleinteriorElectricnm1_);
    CORE::LINALG::multiply(1.0, xVec, 48.0 / 25.0, localSolver_->Emat, ele.eleinteriorElectric_);
  }
  else
  {
    CORE::LINALG::multiply(-1.0, xVec, 1.0, localSolver_->Emat,
        ele.eleinteriorElectric_);  // Implicit euler
                                    // Only if the D matrix is not zero <-> epsilon != 0
    CORE::LINALG::multiply(
        1.0, xVec, -1.0, localSolver_->Dmat, ele.eleinteriorElectricnm1_);  // Implicit euler
    CORE::LINALG::multiply(
        1.0, xVec, 2.0, localSolver_->Dmat, ele.eleinteriorElectric_);  // Implicit euler}
  }
  //  y = [(E + G) - FA^{-1}B]^{-1}^(E - I_s)
  CORE::LINALG::multiply(yVec, tempMat2, xVec);

  CORE::LINALG::multiplyTN(0.0, elevec, -1.0, localSolver_->Hmat, yVec);  //  = -Jy

  CORE::LINALG::multiply(xVec, localSolver_->Bmat, yVec);     //  = By
  CORE::LINALG::multiply(yVec, localSolver_->invAmat, xVec);  //  = A^{-1} By

  CORE::LINALG::multiplyTN(1.0, elevec, 1.0, localSolver_->Cmat, yVec);  //  = Ix - Jy

  return;
}  // UpdateInteriorVariablesAndComputeResidual

/*----------------------------------------------------------------------*
 * ComputeAbsorbingBC
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ComputeAbsorbingBC(
    DRT::Discretization& discretization, DRT::ELEMENTS::ElemagDiff* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat, int face,
    CORE::LINALG::SerialDenseMatrix& elemat, int indexstart,
    CORE::LINALG::SerialDenseVector& elevec1)
{
  dserror("ComputeAbsorbingBC() not yet implemented.");
  /*
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagDiffEleCalc::ComputeAbsorbingBC");

  shapesface_->EvaluateFace(*ele, face);

  // Get the user defined functions
  Teuchos::RCP<DRT::Condition>* cond = params.getPtr<Teuchos::RCP<DRT::Condition>>("condition");
  const std::vector<int>* funct = (*cond)->Get<std::vector<int>>("funct");
  const double time = params.get<double>("time");

  CORE::LINALG::SerialDenseVector tempVec1(shapesface_->nfdofs_ * nsd_);
  CORE::LINALG::SerialDenseVector tempVec2(shapesface_->nfdofs_ * (nsd_ - 1));
  // the RHS matrix has to have the row dimension equal to the number of shape
  // functions(so we have one coefficient for each) and a number of column
  // equal to the overall number of component that we want to solve for.
  // The number is nsd_*2 because we have two fields..
  CORE::LINALG::SerialDenseMatrix localMat(shapesface_->nfdofs_, nsd_ * 2);
  {
    CORE::LINALG::SerialDenseMatrix tempMassMat(shapesface_->nfdofs_, shapesface_->nfdofs_);
    CORE::LINALG::SerialDenseMatrix tempMat(shapesface_->nfdofs_, shapesface_->nqpoints_);
    CORE::LINALG::SerialDenseMatrix tempMatW(shapesface_->nfdofs_, shapesface_->nqpoints_);
    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      // Storing the values of the coordinates for the current quadrature point
      // and of the jacobian computed in that point
      const double fac = shapesface_->jfac(q);
      CORE::LINALG::Matrix<nsd_, 1> xyz;
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz(d) = shapesface_->xyzreal(d, q);  // coordinates of quadrature point in real
  coordinates
      // Creating the temporary electric and magnetic field vector intVal
      // The vector is going to contain first the electric and then the magnetic
      // field such that the field will be initialized as first tree component
      // of the specified function as electric field, last three components as
      // magnetic field. If there is only one component all the components will
      // be initialized to the same value.
      CORE::LINALG::SerialDenseVector intVal(2 * nsd_);
      dsassert(funct != nullptr, "funct not set for initial value");
      EvaluateAll((*funct)[0], time, xyz, intVal);
      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // Mass matrix
        tempMat(i, q) = shapesface_->shfunct(i, q);
        tempMatW(i, q) = shapesface_->shfunct(i, q) * fac;

        // RHS for the electric and magnetic field
        for (int j = 0; j < intVal.numRows(); ++j)
          localMat(i, j) += shapesface_->shfunct(i, q) * intVal(j) * fac;
      }
    }
    // The integration is made by computing the matrix product
CORE::LINALG::multiplyNT(    tempMassMat,  tempMat,  tempMatW);
    {
      using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseMass;
      inverseMass.setMatrix(Teuchos::rcpFromRef(tempMassMat));
      inverseMass.setVectors(Teuchos::rcpFromRef(localMat), Teuchos::rcpFromRef(localMat));
      inverseMass.solve();
    }
  }

  for (unsigned int r = 0; r < shapesface_->nfdofs_; ++r)
    for (unsigned int d = 0; d < nsd_; ++d)
      tempVec1(d * shapesface_->nfdofs_ + r) = localMat(r, d);  // Electric field

  // Creating the matrix
  CORE::LINALG::SerialDenseMatrix transformatrix(
      (nsd_ - 1) * shapesface_->nfdofs_, nsd_ * shapesface_->nfdofs_);
  CORE::LINALG::SerialDenseMatrix inv_transformatrix(
      nsd_ * shapesface_->nfdofs_, (nsd_ - 1) * shapesface_->nfdofs_);
  for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
      for (unsigned int d = 0; d < nsd_ ; ++d)
        for (unsigned int q = 0; q < nsd_-1; ++q)
          if (i == j)
          {
            // I need tangents because I'm translating real coordinates to face ones
            transformatrix(shapesface_->nfdofs_ * q + i, shapesface_->nfdofs_ * d + j) =
                shapesface_->tangent(q, d);
          }

  const MAT::ElectromagneticMat* actmat = static_cast<const
  MAT::ElectromagneticMat*>(mat.get()); double impedance = sqrt(actmat->epsilon(ele->Id()) /
  actmat->mu(ele->Id()));

  // MIXED SHAPE FUNCTIONS
  // The matrix that are going to be build here are D,I and J
  // loop over number of internal shape functions
  // Here we need to create only the first part of tghe D and H matrix to be multiplied by the
  // transformation matrices and then put in the real D and H matrices
  CORE::LINALG::SerialDenseMatrix tempI(shapesface_->nfdofs_ * nsd_, shapesface_->nfdofs_ * nsd_);
  CORE::LINALG::SerialDenseMatrix tempJ(shapesface_->nfdofs_ * nsd_, shapesface_->nfdofs_ * nsd_);
  for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
  {
    // If the shape function is zero on the face we can just skip it. Remember
    // that the matrix have already been set to zero and therefore if nothing
    // is done the value ramains zero
    // loop over number of face shape functions
    for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
    {
      // Now that the integration has been carried on it is necessary to place
      // the value in the right position inside the matrices
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        // i internal shape functions
        // j boundary shape functions
        for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
        {
          // Storing the value of the integral without the normal components
          const double temp =
              shapesface_->jfac(q) * shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q);
        }
      }  // for (unsigned int d = 0; d < nsd_; ++d)
    }    // for (unsigned int j=0; j<shapesface_->nfdofs_; ++j)
  }      // for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)

  // Fill face values into the matrices
  CORE::LINALG::SerialDenseMatrix magneticMat(
      shapesface_->nfdofs_ * (nsd_ - 1), shapesface_->nfdofs_ * nsd_);
  CORE::LINALG::SerialDenseMatrix electricMat(
      shapesface_->nfdofs_ * (nsd_ - 1), shapesface_->nfdofs_ * nsd_);
CORE::LINALG::multiply(  magneticMat,  transformatrix,  tempI);
CORE::LINALG::multiply(  electricMat,  transformatrix,  tempJ);

CORE::LINALG::multiply(0.0,   tempVec2, impedance,  electricMat,  tempVec1);

  for (unsigned int r = 0; r < shapesface_->nfdofs_; ++r)
    for (unsigned int d = 0; d < nsd_; ++d)
      tempVec1(d * shapesface_->nfdofs_ + r) = localMat(r, d + nsd_);  // magnetic

CORE::LINALG::multiply(1.0,   tempVec2, 1.0,  magneticMat,  tempVec1);

  unsigned int newindex = shapesface_->nfdofs_ * (nsd_ - 1) * face;

  for (int i = 0; i < tempVec2.numRows(); ++i) elevec1(newindex + i) = tempVec2(i);

  bool resonly = params.get<bool>("resonly");
  if (!resonly)
  {
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
        for (unsigned int d = 0; d < nsd_ - 1; ++d)
          for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
          {
            const double temp = impedance * shapesface_->jfac(q) * shapesface_->shfunct(i, q) *
                                shapesface_->shfunct(j, q);
            elemat(newindex + shapesface_->nfdofs_ * d + i,
                newindex + shapesface_->nfdofs_ * d + j) += temp;
          }
  }
  */
  return;
}

/*----------------------------------------------------------------------*
 * ComputeSource
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ComputeSource(
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& interiorSourcen,
    CORE::LINALG::SerialDenseVector& interiorSourcenp)
{
  int funcno = params.get<int>("sourcefuncno");
  if (funcno <= 0) return;  // there is no such thing as a volume force

  const double factor = params.get<double>("mod_mu");

  // the vector to be filled
  std::vector<CORE::LINALG::SerialDenseVector> source(2, CORE::LINALG::SerialDenseVector(nsd_));

  // what time is it?
  double tn = params.get<double>("time");
  double tp = params.get<double>("timep");
  double dt = (tp - tn);

  // Cycle Gauss points
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    // Gauss point location
    CORE::LINALG::Matrix<nsd_, 1> xyz;
    for (unsigned int d = 0; d < nsd_; ++d) xyz(d) = shapes_.xyzreal(d, q);

    // Evaluate time derivative of the source term
    // We evaluate at tn and tp as they are already the next time step
    ComputeFunctionTimeDerivative(funcno, tn, dt, xyz, source[0]);
    ComputeFunctionTimeDerivative(funcno, tp, dt, xyz, source[1]);

    // add it all up
    for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        interiorSourcen(i + d * shapes_.ndofs_) +=
            shapes_.shfunct(i, q) * source[0](d) * shapes_.jfac(q) * factor;
        interiorSourcenp(i + d * shapes_.ndofs_) +=
            shapes_.shfunct(i, q) * source[1](d) * shapes_.jfac(q) * factor;
      }
  }

  return;
}

/*----------------------------------------------------------------------*
 * ComputeInteriorMatrices
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ComputeInteriorMatrices(
    double dt, double sigma, double mu, double epsilon)
{
  // The definitions of the matrices created here can be found in the internal
  // paper from Berardocco "A hybridizable discontinous Galerkin method for
  // electromagnetics in subsurface applications".
  // The explicit form of these matrices is reported for convenience?
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagDiffEleCalc::ComputeInteriorMatrices");
  // Why is this made in this order? Is it faster in this order? Or is it better
  // to have it shape_functions->quadrature_points?
  // loop quadrature points
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    // loop shape functions
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      massPart(i, q) = shapes_.shfunct(i, q);
      const double valf = shapes_.shfunct(i, q) * shapes_.jfac(q);
      massPartW(i, q) = valf;
    }
  }

  CORE::LINALG::SerialDenseMatrix tmpMat(ndofs_, ndofs_);
  // this temorary matrix is used to compute the numerical integration and the
  // values are then copied in the right places. Probably it is also possible
  // to have the matrix multiplication to obtain directly the correct matrices
  // but it would mean to compute three time sthe same value for each shape
  // function instead of computing it only omnce and then directly copying it.
  CORE::LINALG::multiplyNT(tmpMat, massPart, massPartW);
  double alpha;
  if (mu < 0.1)
    alpha = 0.5 * (1 + std::log(dt) / std::log(mu));
  else
    alpha = 0.0;
  // A, E and part of G
  for (unsigned int j = 0; j < ndofs_; ++j)
    for (unsigned int i = 0; i < ndofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        Amat(d * ndofs_ + i, d * ndofs_ + j) = -std::pow(mu, 1.0 - alpha) * tmpMat(i, j);
        Dmat(d * ndofs_ + i, d * ndofs_ + j) = epsilon * tmpMat(i, j);
        Emat(d * ndofs_ + i, d * ndofs_ + j) = sigma * std::pow(mu, alpha) * tmpMat(i, j);
      }

  if (dyna_ == INPAR::ELEMAG::elemag_bdf2)
  {
    Emat.scale(3.0 / (2.0 * dt));  // BDF2
    Dmat.scale(2.0 / (dt * dt));   // BDF2
  }
  else if (dyna_ == INPAR::ELEMAG::elemag_bdf4)
  {
    Emat.scale(25.0 / (12.0 * dt));  // BDF4
    if (epsilon != 0) dserror("Not implemented.");
  }
  else
  {
    Emat.scale(1.0 / dt);         // Implicit euler
    Dmat.scale(1.0 / (dt * dt));  // Implicit euler
  }

  {  // We are creating this scope to destroy everything related to the matrix inversion
    // We are going to need both A and its inverse and therefore we are storing both
    invAmat += Amat;
    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> invA;
    invA.setMatrix(Teuchos::rcpFromRef(invAmat));
    int err = invA.invert();
    if (err != 0) dserror("Inversion for Amat failed with errorcode %d", err);
  }

  for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
    for (unsigned int j = 0; j < shapes_.ndofs_; ++j)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
        {
          Bmat(i + d * ndofs_, j + ((d + 1) % nsd_) * ndofs_) +=
              shapes_.shderxy(i * nsd_ + ((d + 2) % nsd_), q) * shapes_.shfunct(j, q) *
              shapes_.jfac(q);
          Bmat(i + d * ndofs_, j + ((d + 2) % nsd_) * ndofs_) -=
              shapes_.shderxy(i * nsd_ + ((d + 1) % nsd_), q) * shapes_.shfunct(j, q) *
              shapes_.jfac(q);
        }
      }

  return;
}

/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ComputeResidual(
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& elevec, double dt,
    DRT::ELEMENTS::ElemagDiff& ele)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagDiffEleCalc::ComputeResidual");

  // for BDF2
  //                                -1
  //                     +---------+    +------------+
  //                     | A    B  |    |     0      |
  // R^{n}  = - [ I  J ] |         |    |            |  =  Ix - Jy
  //                     | F   E+G |    |^E-I_s^{n+2}|
  //                     +---------+    +------------+
  //
  //  x = A^{-1} (By)
  //
  //  y = ((E + G) - F A^{-1} B)^{-1} (^E - I_s^{n+2})

  const unsigned int intdofs = ndofs_ * nsd_;
  // All the vectors are initilized to zero
  CORE::LINALG::SerialDenseVector tempVec1(intdofs);
  CORE::LINALG::SerialDenseVector tempVec2(intdofs);
  // Once the compute source is ready we will need to delete these
  // The ComputeSource is necesessary to include the forcing terms
  ComputeSource(params, tempVec2, tempVec1);

  if (dyna_ == INPAR::ELEMAG::elemag_bdf2)
  {
    // The last -1.0 in the following function has not been removed such that once
    // the ComputeSource() has been created there will be no need to change it
    CORE::LINALG::multiply(
        -1.0, tempVec1, -1.0 / 3.0, Emat, ele.eleinteriorElectricnm1_);  // (1/3)E E^{n} + I_s
    CORE::LINALG::multiply(1.0, tempVec1, 4.0 / 3.0, Emat,
        ele.eleinteriorElectric_);  // ^E = (4/3)EE^{n+1} - (1/3)EE^{n} - I_s
                                    // Only if the D matrix is not zero <-> epsilon != 0
    CORE::LINALG::multiply(1.0, tempVec1, 1.0 / 2.0, Dmat, ele.eleinteriorElectricnm2_);
    CORE::LINALG::multiply(1.0, tempVec1, -2.0, Dmat, ele.eleinteriorElectricnm1_);
    CORE::LINALG::multiply(1.0, tempVec1, 5.0 / 2.0, Dmat, ele.eleinteriorElectric_);
  }
  else if (dyna_ == INPAR::ELEMAG::elemag_bdf4)
  {
    // The last -1.0 in the following function has not been removed such that once
    // the ComputeSource() has been created there will be no need to change it
    CORE::LINALG::multiply(
        -1.0, tempVec1, -3.0 / 25.0, Emat, ele.eleinteriorElectricnm3_);  // (1/3)E E^{n} + I_s
    CORE::LINALG::multiply(1.0, tempVec1, 16.0 / 25.0, Emat,
        ele.eleinteriorElectricnm2_);  // ^E = (4/3)EE^{n+1} - (1/3)EE^{n} - I_s
    CORE::LINALG::multiply(1.0, tempVec1, -36.0 / 25.0, Emat, ele.eleinteriorElectricnm1_);
    CORE::LINALG::multiply(1.0, tempVec1, 48.0 / 25.0, Emat, ele.eleinteriorElectric_);
  }
  else
  {
    // Implicit euler
    CORE::LINALG::multiply(-1.0, tempVec1, 1.0, Emat,
        ele.eleinteriorElectric_);  // E E^{n} -\dot{I}_s
                                    // Only if the D matrix is not zero <-> epsilon != 0
    CORE::LINALG::multiply(
        1.0, tempVec1, -1.0, Dmat, ele.eleinteriorElectricnm1_);  // E E^{n} -\dot{I}_s
    CORE::LINALG::multiply(
        1.0, tempVec1, 2.0, Dmat, ele.eleinteriorElectric_);  // E E^{n} -\dot{I}_s
  }

  CORE::LINALG::SerialDenseMatrix tempMat1(intdofs, intdofs);
  CORE::LINALG::multiplyTN(tempMat1, Bmat, invAmat);  // F A^{-1}

  CORE::LINALG::SerialDenseMatrix tempMat2(intdofs, intdofs);

  tempMat2 += Emat;
  tempMat2 += Gmat;
  tempMat2 += Dmat;
  CORE::LINALG::multiply(1.0, tempMat2, -1.0, tempMat1, Bmat);  // = (E + G) - F A^{-1} B
  {
    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseinW;
    inverseinW.setMatrix(Teuchos::rcpFromRef(tempMat2));
    int err = inverseinW.invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d", err);
  }
  // tempMat2 = ((E + G) - F A^{-1} B)^{-1}

  CORE::LINALG::multiply(tempVec2, tempMat2, tempVec1);         // y
  CORE::LINALG::multiplyTN(0.0, elevec, -1.0, Hmat, tempVec2);  //  -Jy

  CORE::LINALG::multiply(tempVec1, Bmat, tempVec2);            // By
  CORE::LINALG::multiply(tempVec2, invAmat, tempVec1);         //  x = A^{-1} By
  CORE::LINALG::multiplyTN(1.0, elevec, 1.0, Cmat, tempVec2);  //  Ix - Jy

  return;
}  // ComputeResidual

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ComputeFaceMatrices(const int face,
    double dt, int indexstart, int newindex, double sigma, double mu, const double tau)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagDiffEleCalc::ComputeFaceMatrices");

  // Tau is defined as (\frac{|sigma|}{\mu t_c})^0.5 where t_c is a
  // characteristic time scale
  // const double tau = tau ? tau : sqrt(sigma / mu / dt);
  // This routine seems complex but it's not (well, it is just as the others)
  // It is divided in three parts:
  //  o   Mixed shape functions integration
  //  o   Interior shaoe function integration
  //  o   Boundary shape function integration
  // The difference lays on the number of dofs per unknown (it depends on the
  // space where we are looking for solutions) and therefore there will be three
  // big groups of nested for loops

  // Be carefull about the fact that this routin is calld once per each face of
  // the element and the convention of grouping the shape functions per spatial
  // dimension (first all those for x, then those for y and so on) is respected
  // on a face basis. Therefore expect to have submatrices of
  // shapesface_->nfdofs_*shapesface_->nfdofs_ or shapesface_->nfdofs_*ndofs_
  //(or vice-versa), divided by nsd_ * "matrix dimension" submatrices with all
  // entries set to zero.
  CORE::LINALG::SerialDenseMatrix transformatrix(
      (nsd_ - 1) * shapesface_->nfdofs_, nsd_ * shapesface_->nfdofs_);
  for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int q = 0; q < nsd_ - 1; ++q)
        transformatrix(shapesface_->nfdofs_ * q + i, shapesface_->nfdofs_ * d + i) =
            shapesface_->tangent(d, q);



  // MIXED SHAPE FUNCTIONS
  // The matrix that are going to be build here are D,I and J
  // loop over number of internal shape functions
  // Here we need to create only the first part of tghe D and H matrix to be multiplied by the
  // transformation matrices and then put in the real D and H matrices
  CORE::LINALG::SerialDenseMatrix tempC(ndofs_ * nsd_, shapesface_->nfdofs_ * nsd_);
  CORE::LINALG::SerialDenseMatrix tempH(ndofs_ * nsd_, shapesface_->nfdofs_ * nsd_);
  for (unsigned int i = 0; i < ndofs_; ++i)
  {
    // If the shape function is zero on the face we can just skip it. Remember
    // that the matrix have already been set to zero and therefore if nothing
    // is done the value ramains zero
    if (shapesface_->shfunctI.NonzeroOnFace(i))
    {
      // loop over number of face shape functions
      for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
      {
        // Now that the integration has been carried on it is necessary to place
        // the value in the right position inside the matrices
        for (unsigned int d = 0; d < nsd_; ++d)
        {
          // i internal shape functions
          // j boundary shape functions
          for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
          {
            // Storing the value of the integral without the normal components
            const double temp = tau * shapesface_->jfac(q) * shapesface_->shfunctI(i, q) *
                                shapesface_->shfunct(j, q);
            const double temp2 =
                shapesface_->jfac(q) * shapesface_->shfunctI(i, q) * shapesface_->shfunct(j, q);
            // Filling the matrices
            // Cmat
            //+1 coordinate
            tempC(d * ndofs_ + i, shapesface_->nfdofs_ * ((d + 1) % nsd_) + j) -=
                temp2 * shapesface_->normals(((d + 2) % nsd_), q);
            //+2 coordinate
            tempC(d * ndofs_ + i, shapesface_->nfdofs_ * ((d + 2) % nsd_) + j) +=
                temp2 * shapesface_->normals(((d + 1) % nsd_), q);
            // Hmat
            // 0 coordinate
            tempH(d * ndofs_ + i, shapesface_->nfdofs_ * d + j) -= temp;
          }
        }  // for (unsigned int d = 0; d < nsd_; ++d)
      }    // for (unsigned int j=0; j<ndofs_; ++j)
    }      // if( shapesface_->shfunctI.NonzeroOnFace(i) )
  }        // for (unsigned int i = 0; i < ndofs_; ++i)

  // Fill face values into the matrices
  {
    CORE::LINALG::SerialDenseMatrix tempMat1(ndofs_ * nsd_, shapesface_->nfdofs_ * (nsd_ - 1));
    CORE::LINALG::SerialDenseMatrix tempMat2(ndofs_ * nsd_, shapesface_->nfdofs_ * (nsd_ - 1));
    CORE::LINALG::multiplyNT(tempMat1, tempC, transformatrix);
    CORE::LINALG::multiplyNT(tempMat2, tempH, transformatrix);

    for (unsigned int i = 0; i < ndofs_ * nsd_; ++i)
      for (unsigned int j = 0; j < shapesface_->nfdofs_ * (nsd_ - 1); ++j)
      {
        Cmat(i, newindex + j) = tempMat1(i, j);
        Hmat(i, newindex + j) = tempMat2(i, j);
      }
  }


  // BOUNDARY SHAPE FUNCTIONS
  // loop over number of shape functions
  for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
  {
    // loop over number of shape functions
    for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
    {
      for (unsigned int d = 0; d < nsd_ - 1; ++d)
      {
        // If the face is perpendicular to the d direction it is necessary to
        // enforce the component of the hybrid variable to be zero because we
        // know that the hybrid variable is defined as the perpendicular
        // component of the elctric field.
        for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
        {
          const double temp =
              tau * shapesface_->jfac(q) * shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q);
          Lmat(newindex + shapesface_->nfdofs_ * d + i, newindex + shapesface_->nfdofs_ * d + j) +=
              temp;
        }
      }
    }  // for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
  }    // for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)

  // INTERIOR SHAPE FUNCTIONS
  // Some terms are still missing in G!!
  for (unsigned int i = 0; i < ndofs_; ++i)
    for (unsigned int j = 0; j < ndofs_; ++j)
      if (shapesface_->shfunctI.NonzeroOnFace(i) && shapesface_->shfunctI.NonzeroOnFace(j))
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
          {
            const double temp = tau * shapesface_->jfac(q) * shapesface_->shfunctI(i, q) *
                                shapesface_->shfunctI(j, q);
            // Gmat
            // 0 coordinate
            Gmat(d * ndofs_ + i, d * ndofs_ + j) +=
                temp * (std::pow(shapesface_->normals(((d + 1) % nsd_), q), 2) +
                           std::pow(shapesface_->normals(((d + 2) % nsd_), q), 2));
            //+1 coordinate
            Gmat(d * ndofs_ + i, ((d + 1) % nsd_) * ndofs_ + j) -=
                temp * shapesface_->normals(d, q) * shapesface_->normals(((d + 1) % nsd_), q);
            //+2 coordinate
            Gmat(d * ndofs_ + i, ((d + 2) % nsd_) * ndofs_ + j) -=
                temp * shapesface_->normals(d, q) * shapesface_->normals(((d + 2) % nsd_), q);
          }

  return;
}  // ComputeFaceMatrices


/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::CondenseLocalPart(
    CORE::LINALG::SerialDenseMatrix& eleMat)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagDiffEleCalc::CondenseLocalPart");

  // THE MATRIX
  //                             -1
  //                   +--------+    +-----+
  //                   |        |    |     |
  //                   | A   B  |    |  C  |
  //  K = L - [ I  J ] |        |    |     |  = L - I X - J Y
  //                   | F  E+G |    |  H  |
  //                   +--------+    +-----+

  //   Y = [ (E+G) - F A^{-1} B ]^{-1} [ H - F A^{-1} C]

  //   X = A^{-1} [ C - B Y ]

  const unsigned int onfdofs = eleMat.numRows();
  const unsigned int intdofs = ndofs_ * nsd_;

  // Thi can be useful to remember when coding
  // int 	Multiply (char TransA, char TransB, double ScalarAB, Matrix &A, Matrix &B, double
  // ScalarThis) this = ScalarThis*this + ScalarAB*A*B
  CORE::LINALG::SerialDenseMatrix tempMat1(intdofs, intdofs);
  CORE::LINALG::multiplyTN(tempMat1, Bmat, invAmat);  // =  F A^{-1}

  CORE::LINALG::SerialDenseMatrix tempMat2(intdofs, intdofs);

  // This is E+G
  tempMat2 += Emat;  // = E
  tempMat2 += Gmat;  // = E + G
  tempMat2 += Dmat;

  CORE::LINALG::multiply(1.0, tempMat2, -1.0, tempMat1, Bmat);  // = (E+G) - F A^{-1} B

  CORE::LINALG::SerialDenseMatrix tempMat3(intdofs, onfdofs);
  tempMat3 += Hmat;  // = H

  CORE::LINALG::multiply(1.0, tempMat3, -1.0, tempMat1, Cmat);  // = H - F A^{-1} C

  // Inverting the first part of the Y matrix
  {
    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> inverseinW;
    inverseinW.setMatrix(Teuchos::rcpFromRef(tempMat2));
    int err = inverseinW.invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d", err);
  }
  // tempMat2 = [(E+G) - F A^{-1} B]^{-1}

  eleMat = Lmat;  // = L
  // reusing matrix that are not needed
  tempMat1.shape(intdofs, onfdofs);
  CORE::LINALG::multiply(
      tempMat1, tempMat2, tempMat3);  //  Y = [(E+G) - F A^{-1} B]^{-1}(H - F A^{-1} C)
  CORE::LINALG::multiplyTN(1.0, eleMat, -1.0, Hmat, tempMat1);  // = L - J Y

  tempMat2.shape(intdofs, onfdofs);
  tempMat2 = Cmat;
  CORE::LINALG::multiply(1.0, tempMat2, -1.0, Bmat, tempMat1);  // = C - B Y

  tempMat3.shape(intdofs, onfdofs);
  CORE::LINALG::multiply(tempMat3, invAmat, tempMat2);  // = X = A^{-1} ( C - B Y )

  CORE::LINALG::multiplyTN(1.0, eleMat, -1.0, Cmat, tempMat3);  // = K = L - I X - J y

  return;
}  // CondenseLocalPart

/*----------------------------------------------------------------------*
 * Symmetrify
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::Symmetrify(
    DRT::ELEMENTS::ElemagDiff& ele, CORE::LINALG::SerialDenseMatrix& eleMat, bool dodirich)
{
  if (ele.lm_.lmdirich_.size())
    for (int i = 0; i < eleMat.numRows(); ++i)
      if (ele.lm_.lmdirich_[i])
        for (int j = 0; j < eleMat.numCols(); ++j)
        {
          eleMat(i, j) = 0.0;
          eleMat(j, i) = 0.0;
        }

  for (int i = 0; i < eleMat.numRows(); ++i)
    for (int j = i; j < eleMat.numCols(); ++j) eleMat(j, i) = eleMat(i, j);

  return;
}  // Symmetrify

/*----------------------------------------------------------------------*
 * Compute internal and face matrices
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ElemagDiffEleCalc<distype>::LocalSolver::ComputeMatrices(
    DRT::Discretization& discretization, const Teuchos::RCP<MAT::Material>& mat,
    DRT::ELEMENTS::ElemagDiff& ele, double dt, INPAR::ELEMAG::DynamicType dyna, const double tau)
{
  // The material properties change elementwise or can also be computed pointwise?
  // Check current_informations, \chapter{Elements and materials for electromagnetics},
  // \section{Remarks}
  const MAT::ElectromagneticMat* elemagmat = static_cast<const MAT::ElectromagneticMat*>(mat.get());
  double sigma = elemagmat->sigma(ele.Id());
  double epsilon = elemagmat->epsilon(ele.Id());
  double mu = elemagmat->mu(ele.Id());

  // Why this? Why do we need to make these matrices zero here? Why not all of them?
  // init face matrices
  invAmat.putScalar(0.0);
  Amat.putScalar(0.0);
  Bmat.putScalar(0.0);
  Cmat.putScalar(0.0);
  Dmat.putScalar(0.0);
  Emat.putScalar(0.0);
  Gmat.putScalar(0.0);
  Hmat.putScalar(0.0);
  Lmat.putScalar(0.0);

  // Here is the computation for the matrices of volume integrals
  ComputeInteriorMatrices(dt, sigma, mu, epsilon);

  // sumindex is going to be used to decide where we are inside the face matrix
  // because for every face we move to different dofs
  int sumindex = 0;
  int newindex = 0;
  for (unsigned int face = 0; face < nfaces_; ++face)
  {
    /* This part is to be used for efficiency reasons, at the beginning the
    //standard procedure is used
    CORE::DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele.Faces()[face]->Degree(),
        shapes_.usescompletepoly_, 2 * ele.Faces()[face]->Degree());
    shapesface_ = CORE::DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    */

    // Updating face data
    shapesface_->EvaluateFace(ele, face);

    // Here are the matrices for the boundary integrals
    ComputeFaceMatrices(face, dt, sumindex, newindex, sigma, mu, tau);
    sumindex += nsd_ * shapesface_->nfdofs_;
    newindex += (nsd_ - 1) * shapesface_->nfdofs_;
  }

  return;
}  // ComputeMatrices


// template classes
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::pyramid5>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::nurbs9>;
template class DRT::ELEMENTS::ElemagDiffEleCalc<CORE::FE::CellType::nurbs27>;
