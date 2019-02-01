/*--------------------------------------------------------------------------*/
/*!
\file elemag_ele_calc.cpp
\brief All functionality for electromagnetic element evaluations

<pre>
\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244
</pre>
 */
/*--------------------------------------------------------------------------*/

#include "elemag_ele_calc.H"
#include "elemag_ele_action.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_geometry/position_array.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_mat/electromagnetic.H"

#include <Epetra_SerialDenseSolver.h>
#include <Teuchos_TimeMonitor.hpp>


namespace
{
  void zeroMatrix(Epetra_SerialDenseMatrix& mat)
  {
    std::memset(mat.A(), 0, sizeof(double) * mat.M() * mat.N());
  }

  void reshapeMatrixIfNecessary(Epetra_SerialDenseMatrix& matrix, const int nrows, const int ncols)
  {
    if (nrows != matrix.M() || ncols != matrix.N()) matrix.Shape(nrows, ncols);
  }
}  // namespace

/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ElemagEleCalc<distype>::ElemagEleCalc()
{
}

/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ElemagEleCalc<distype>::Evaluate(DRT::ELEMENTS::Elemag* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
    const DRT::UTILS::GaussIntegration&, bool offdiag)
{
  return this->Evaluate(ele, discretization, lm, params, mat, elemat1_epetra, elemat2_epetra,
      elevec1_epetra, elevec2_epetra, elevec3_epetra, offdiag);
}

/*----------------------------------------------------------------------*
 * Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ElemagEleCalc<distype>::Evaluate(DRT::ELEMENTS::Elemag* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix&,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector&,
    bool offdiag)
{
  // check if this is an hdg element and init completepoly
  if (const DRT::ELEMENTS::Elemag* hdgele = dynamic_cast<const DRT::ELEMENTS::Elemag*>(ele))
    usescompletepoly_ = hdgele->UsesCompletePolynomialSpace();
  else
    dserror("cannot cast element to elemag element");
  const ELEMAG::Action action = DRT::INPUT::get<ELEMAG::Action>(params, "action");

  InitializeShapes(ele);

  bool updateonly = false;
  shapes_->Evaluate(*ele);
  switch (action)
  {
    case ELEMAG::project_field:
    {
      ElementInit(ele, params);

      localSolver_->ProjectField(ele, params, elevec1, elevec2);
      break;
    }
    case ELEMAG::compute_error:
    {
      localSolver_->ComputeError(ele, params, elevec1);
      break;
    }
    case ELEMAG::project_field_test:
    {
      localSolver_->ProjectFieldTest(ele, params, elevec1, elevec2);
      break;
    }
    case ELEMAG::project_field_test_trace:
    {
      // ElementInit(ele, params);
      localSolver_->ProjectFieldTestTrace(ele, params, elevec1);

      break;
    }
    case ELEMAG::project_dirich_field:
    {
      // if (mat->MaterialType() != INPAR::MAT::m_electromagneticmat)
      //  dserror("for physical type 'lossless' please supply MAT_Electromagnetic");
      if (params.isParameter("faceconsider"))
      {
        // ElementInit(ele, params);
        localSolver_->ProjectDirichField(ele, params, elevec1);
      }
      break;
    }
    case ELEMAG::ele_init:
    {
      ElementInit(ele, params);
      break;
    }
    case ELEMAG::fill_restart_vecs:
    {
      // bool padapty = params.get<bool>("padaptivity");
      ReadGlobalVectors(ele, discretization, lm);
      FillRestartVectors(ele, discretization);
      break;
    }
    case ELEMAG::ele_init_from_restart:
    {
      ElementInit(ele, params);
      ElementInitFromRestart(ele, discretization);
      break;
    }
    case ELEMAG::interpolate_hdg_to_node:
    {
      ReadGlobalVectors(ele, discretization, lm);
      InterpolateSolutionToNodes(ele, discretization, elevec1);
      break;
    }
    case ELEMAG::calc_abc:
    {
      int face = params.get<int>("face");
      int sumindex = 0;
      for (int i = 0; i < face; ++i)
      {
        DRT::UTILS::PolynomialSpaceParams params(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,
            ele->Faces()[i]->Degree(), usescompletepoly_);
        int nfdofs = DRT::UTILS::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(params)->Size();
        sumindex += nfdofs;
      }
      ReadGlobalVectors(ele, discretization, lm);
      if (!params.isParameter("nodeindices"))
        localSolver_->ComputeAbsorbingBC(
            discretization, ele, params, mat, face, elemat1, sumindex, elevec1);
      else
        dserror("why would you set an absorbing LINE in THREE dimensions?");

      break;
    }
    /*
    case ELEMAG::bd_integrate:
    {
      int face = params.get<int>("face");
      localSolver_->ComputeBoundaryIntegral(ele, params, face);

      break;
    }
    */
    case ELEMAG::calc_systemmat_and_residual:
    {
      // const bool resonly = params.get<bool>("resonly");
      // const bool padapty = params.get<bool>("padaptivity");
      double dt = params.get<double>("dt");
      dyna_ = params.get<INPAR::ELEMAG::DynamicType>("dynamic type");

      ReadGlobalVectors(ele, discretization, lm);
      zeroMatrix(elevec1);
      localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);

      // if (!resonly)
      localSolver_->CondenseLocalPart(elemat1);

      localSolver_->ComputeResidual(params, elevec1, interiorMagneticnp_, interiorElectricnp_);

      break;
    }
    case ELEMAG::update_secondary_solution:
      updateonly = true;  // no break here!!!
    case ELEMAG::update_secondary_solution_and_calc_residual:
    {
      // bool errormaps = params.get<bool>("errormaps");
      bool errormaps = false;
      // const bool allelesequal = params.get<bool>("allelesequal");

      double dt = params.get<double>("dt");
      dyna_ = params.get<INPAR::ELEMAG::DynamicType>("dynamic type");

      ReadGlobalVectors(ele, discretization, lm);

      zeroMatrix(elevec1);
      localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);
      /* Could be useful for optimization purposes
      if(!allelesequal)
        localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);
      */

      UpdateInteriorVariablesAndComputeResidual(
          params, *ele, mat, elevec1, dt, errormaps, updateonly);

      break;
    }
    case ELEMAG::get_gauss_points:
    {
      int rows = shapes_->xyzreal.M();
      int cols = shapes_->xyzreal.N();
      elemat1.Shape(rows, cols);

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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::PrintTrace(DRT::Element* ele)
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

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::InitializeShapes(const DRT::ELEMENTS::Elemag* ele)
{
  if (shapes_ == Teuchos::null)
    shapes_ = Teuchos::rcp(
        new DRT::UTILS::ShapeValues<distype>(ele->Degree(), usescompletepoly_, 2 * ele->Degree()));
  else if (shapes_->degree_ != unsigned(ele->Degree()) ||
           shapes_->usescompletepoly_ != usescompletepoly_)
    shapes_ = Teuchos::rcp(
        new DRT::UTILS::ShapeValues<distype>(ele->Degree(), usescompletepoly_, 2 * ele->Degree()));

  if (shapesface_ == Teuchos::null)
  {
    DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele->Degree(), usescompletepoly_, 2 * ele->Degree());
    shapesface_ = Teuchos::rcp(new DRT::UTILS::ShapeValuesFace<distype>(svfparams));
  }

  if (localSolver_ == Teuchos::null)
    localSolver_ = Teuchos::rcp(new LocalSolver(ele, *shapes_, shapesface_, dyna_));
  else if (localSolver_->ndofs_ != shapes_->ndofs_)
    localSolver_ = Teuchos::rcp(new LocalSolver(ele, *shapes_, shapesface_, dyna_));
}

/*----------------------------------------------------------------------*
 * ReadGlobalVectors
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::ReadGlobalVectors(
    DRT::Element* ele, DRT::Discretization& discretization, const std::vector<int>& lm)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagEleCalc::ReadGlobalVectors");
  DRT::ELEMENTS::Elemag* elemagele = dynamic_cast<DRT::ELEMENTS::Elemag*>(ele);

  // read vectors from element storage
  reshapeMatrixIfNecessary(interiorElectricnp_, elemagele->eleinteriorElectric_.M(), 1);
  reshapeMatrixIfNecessary(interiorMagneticnp_, elemagele->eleinteriorMagnetic_.M(), 1);

  interiorElectricnp_ = elemagele->eleinteriorElectric_;
  interiorMagneticnp_ = elemagele->eleinteriorMagnetic_;

  // read vectors from time integrator
  if (discretization.HasState("trace"))  // in case of "update interior variables"
  {
    reshapeMatrixIfNecessary(elemagele->elenodeTrace2d_, lm.size(), 1);
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState("trace");
    DRT::UTILS::ExtractMyValues(*matrix_state, elemagele->elenodeTrace2d_, lm);
  }

  return;
}  // ReadGlobalVectors

/*----------------------------------------------------------------------*
 * FillRestartVectors
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::FillRestartVectors(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  // sort this back to the interior values vector
  int size = shapes_->ndofs_ * nsd_ * 2;

  std::vector<double> interiorVar(size);
  for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
  {
    interiorVar[i] = interiorMagneticnp_(i);
    interiorVar[shapes_->ndofs_ * nsd_ + i] = interiorElectricnp_(i);
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

  return;
}

/*----------------------------------------------------------------------*
 * ElementInitFromRestart
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::ElementInitFromRestart(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  DRT::ELEMENTS::Elemag* elemagele = dynamic_cast<DRT::ELEMENTS::Elemag*>(ele);
  int size = shapes_->ndofs_ * nsd_ * 2;

  std::vector<double> interiorVar(size);

  Teuchos::RCP<const Epetra_Vector> intVar = discretization.GetState(1, "intVar");
  std::vector<int> localDofs1 = discretization.Dof(1, ele);
  DRT::UTILS::ExtractMyValues(*intVar, interiorVar, localDofs1);

  // now write this in corresponding eleinteriorElectric_ and eleinteriorMagnetic_
  for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
  {
    elemagele->eleinteriorMagnetic_(i) = interiorVar[i];
    elemagele->eleinteriorElectric_(i) = interiorVar[shapes_->ndofs_ * nsd_ + i];
  }

  return;
}

/*----------------------------------------------------------------------*
 * Element init
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::ElementInit(
    DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params)
{
  // each element has to store the interior vectors by itseld, p-adaptivity or not
  // so, shape it, as you need it
  ele->eleinteriorElectric_.Shape(shapes_->ndofs_ * nsd_, 1);
  ele->eleinteriorMagnetic_.Shape(shapes_->ndofs_ * nsd_, 1);

  ele->elenodeTrace_.Shape(ele->NumFace() * shapesface_->nfdofs_ * nsd_, 1);
  ele->elenodeTrace2d_.Shape(ele->NumFace() * shapesface_->nfdofs_ * (nsd_ - 1), 1);

  return;
}

/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ProjectField(DRT::ELEMENTS::Elemag* ele,
    Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2)
{
  shapes_.Evaluate(*ele);

  // get function
  const int* start_func = params.getPtr<int>("startfuncno");
  const double time = params.get<double>("time");

  // the RHS matrix has to have the row dimension equal to the number of shape
  // functions(so we have one coefficient for each) and a number of column
  // equal to the overall number of component that we want to solve for.
  // The number is nsd_*2 because we have two fields..
  Epetra_SerialDenseMatrix localMat(shapes_.ndofs_, nsd_ * 2);
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    // Storing the values of the coordinates for the current quadrature point
    // and of the jacobian computed in that point
    const double fac = shapes_.jfac(q);
    LINALG::Matrix<nsd_, 1> xyz;
    for (unsigned int d = 0; d < nsd_; ++d)
      xyz(d) = shapes_.xyzreal(d, q);  // coordinates of quadrature point in real coordinates
    // Creating the temporary electric and magnetic field vector intVal
    // The vector is going to contain first the electric and then the magnetic
    // field such that the field will be initialized as first tree component
    // of the specified function as electric field, last three components as
    // magnetic field. If there is only one component all the components will
    // be initialized to the same value.
    Epetra_SerialDenseVector intVal(2 * nsd_);
    dsassert(start_func != NULL, "funct not set for initial value");
    EvaluateAll(*start_func, time, xyz, intVal);
    // now fill the components in the one-sided mass matrix and the right hand side
    for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
    {
      // Mass matrix
      massPart(i, q) = shapes_.shfunct(i, q);
      massPartW(i, q) = shapes_.shfunct(i, q) * fac;

      // RHS for the electric and magnetic field
      for (int j = 0; j < intVal.M(); ++j)
        localMat(i, j) += shapes_.shfunct(i, q) * intVal(j) * fac;
    }
  }
  // The integration is made by computing the matrix product
  massMat.Multiply('N', 'T', 1., massPart, massPartW, 0.);
  {
    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(massMat);
    inverseMass.SetVectors(localMat, localMat);
    inverseMass.Solve();
  }

  // Here we move the values from the temporary variable to the variable
  // contained in the element
  for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
  {
    // Now we are storing the variables by component, meaning that we save for
    // each component the value for each dof and then we move to the next component.
    for (unsigned int d = 0; d < nsd_; ++d)
    {
      ele->eleinteriorElectric_(d * shapes_.ndofs_ + r) = localMat(r, d);         // Electric field
      ele->eleinteriorMagnetic_(d * shapes_.ndofs_ + r) = localMat(r, d + nsd_);  // magnetic
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * ComputeError
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ComputeError(
    DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1)
{
  double error_ele = 0.0, error_mag = 0.0;
  double exact_ele = 0.0, exact_mag = 0.0;
  shapes_.Evaluate(*ele);

  // get function
  const int func = params.get<int>("funcno");
  const double time = params.get<double>("time");
  // for the calculation of the error, we use a higher integration rule
  Teuchos::RCP<DRT::UTILS::GaussPoints> highquad =
      DRT::UTILS::GaussPointCache::Instance().Create(distype, (ele->Degree() + 2) * 2);
  LINALG::Matrix<nsd_, 1> xsi;
  Epetra_SerialDenseVector values(shapes_.ndofs_);
  LINALG::Matrix<nsd_, nen_> deriv;
  LINALG::Matrix<nsd_, nsd_> xjm;
  Epetra_SerialDenseVector electric(nsd_);
  Epetra_SerialDenseVector magnetic(nsd_);
  Epetra_SerialDenseVector analytical(2 * nsd_);

  for (int q = 0; q < highquad->NumPoints(); ++q)
  {
    const double* gpcoord = highquad->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++) xsi(idim) = gpcoord[idim];
    shapes_.polySpace_->Evaluate(xsi, values);

    DRT::UTILS::shape_function_deriv1<distype>(xsi, deriv);
    xjm.MultiplyNT(deriv, shapes_.xyze);
    double highjfac = xjm.Determinant() * highquad->Weight(q);

    electric.Scale(0.0);
    magnetic.Scale(0.0);
    for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        electric(d) += values(i) * ele->eleinteriorElectric_(d * shapes_.ndofs_ + i);
        magnetic(d) += values(i) * ele->eleinteriorMagnetic_(d * shapes_.ndofs_ + i);
      }

    LINALG::Matrix<nen_, 1> myfunct;
    DRT::UTILS::shape_function<distype>(xsi, myfunct);
    LINALG::Matrix<nsd_, 1> xyzmat;
    xyzmat.MultiplyNN(shapes_.xyze, myfunct);

    // Creating the temporary electric and magnetic field vector intVal
    // The vector is going to contain first the electric and then the magnetic
    // field such that the field will be initialized as first tree component
    // of the specified function as electric field, last three components as
    // magnetic field. If there is only one component all the components will
    // be initialized to the same value.
    analytical.Scale(0.0);
    EvaluateAll(func, time, xyzmat, analytical);

    for (unsigned int d = 0; d < nsd_; ++d)
    {
      error_ele += pow((analytical(d) - electric(d)), 2) * highjfac;
      exact_ele += pow(analytical(d), 2) * highjfac;
      error_mag += pow((analytical(d + nsd_) - magnetic(d)), 2) * highjfac;
      exact_mag += pow(analytical(d + nsd_), 2) * highjfac;
    }
  }

  elevec1[0] = error_ele;
  elevec1[1] = exact_ele;
  elevec1[2] = error_mag;
  elevec1[3] = exact_mag;

  return;
}

/*----------------------------------------------------------------------*
 * ProjectFieldTest
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ProjectFieldTest(DRT::ELEMENTS::Elemag* ele,
    Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2)
{
  shapes_.Evaluate(*ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 || unsigned(elevec2.M()) == nsd_ * shapes_.ndofs_,
      "Wrong size in project vector 2");

  // get function
  const int* start_func = params.getPtr<int>("startfuncno");
  const double time = params.get<double>("time");

  // internal variables
  if (elevec2.M() > 0)
  {
    // the RHS matrix has to have the row dimension equal to the number of shape
    // functions(so we have one coefficient for each) and a number of column
    // equal to the overall number of component that we want to solve for.
    // The number is nsd_*2 because we have two fields..
    Epetra_SerialDenseMatrix localMat(shapes_.ndofs_, nsd_ * 2);
    for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
    {
      // Storing the values of the coordinates for the current quadrature point
      // and of the jacobian computed in that point
      const double fac = shapes_.jfac(q);
      LINALG::Matrix<nsd_, 1> xyz;
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz(d) = shapes_.xyzreal(d, q);  // coordinates of quadrature point in real coordinates
      // Creating the temporary electric and magnetic field vector intVal
      // The vector is going to contain first the electric and then the magnetic
      // field such that the field will be initialized as first tree component
      // of the specified function as electric field, last three components as
      // magnetic field. If there is only one component all the components will
      // be initialized to the same value.
      Epetra_SerialDenseVector intVal(2 * nsd_);
      dsassert(start_func != NULL, "funct not set for initial value");
      EvaluateAll(*start_func, time, xyz, intVal);
      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      {
        // Mass matrix
        massPart(i, q) = shapes_.shfunct(i, q);
        massPartW(i, q) = shapes_.shfunct(i, q) * fac;

        // RHS for the electric and magnetic field
        for (int j = 0; j < intVal.M(); ++j)
          localMat(i, j) += shapes_.shfunct(i, q) * intVal(j) * fac;
      }
    }
    // The integration is made by computing the matrix product
    massMat.Multiply('N', 'T', 1., massPart, massPartW, 0.);
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(massMat);
      inverseMass.SetVectors(localMat, localMat);
      inverseMass.Solve();
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
  return 0;
}

/*----------------------------------------------------------------------*
 * ProjectFieldTestTrace
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ProjectFieldTestTrace(
    DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1)
{
  // Here we have the projection of the field on the trace
  // mass is the mass matrix for the system to be solved
  // the dimension of the mass matrix is given by the number of shape functions
  Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
  // TRaceVEC is the vector of the trace values
  // instead of being a vector it is a matrix so that we use the same matrix
  // to solve the projection problem on every component of the field
  Epetra_SerialDenseMatrix trVec(shapesface_->nfdofs_, nsd_);

  dsassert(elevec1.M() == static_cast<int>(nsd_ * shapesface_->nfdofs_) ||
               elevec1.M() == static_cast<int>(nfaces_ * nsd_ * shapesface_->nfdofs_),
      "Wrong size in project vector 1");

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
    zeroMatrix(mass);
    zeroMatrix(trVec);

    // Cycling through the quadrature points
    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      // For each quadrature point we have a vector containing the field
      // components and a vector containing the spatial coordinates of that point
      Epetra_SerialDenseVector trace(nsd_);
      LINALG::Matrix<nsd_, 1> xyz;

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

    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec, trVec);
    inverseMass.Solve();

    Epetra_SerialDenseVector tempVec(shapesface_->nfdofs_ * (nsd_));
    Epetra_SerialDenseVector faceVec(shapesface_->nfdofs_ * (nsd_ - 1));
    // Filling the vector of trace values
    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // remember that "f" is an iterator index and therefore we are
        // cycling through all the faces and all the entries of elevec1
        // except for the first one where we will put the pressure average
        tempVec(d * shapesface_->nfdofs_ + i) = trVec(i, d);
      }

    Epetra_SerialDenseMatrix transformatrix(
        (nsd_ - 1) * shapesface_->nfdofs_, nsd_ * shapesface_->nfdofs_);
    Epetra_SerialDenseMatrix inv_transformatrix(
        nsd_ * shapesface_->nfdofs_, (nsd_ - 1) * shapesface_->nfdofs_);
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    {
      for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
      {
        for (unsigned int d = 0; d < nsd_ - 1; ++d)
        {
          for (unsigned int q = 0; q < nsd_; ++q)
          {
            // I need tangents because I'm translating real coordinates to face ones
            if (i == j)
            {
              transformatrix(shapesface_->nfdofs_ * d + i, shapesface_->nfdofs_ * q + j) =
                  shapesface_->tangents(d, q);
              inv_transformatrix(shapesface_->nfdofs_ * q + i, shapesface_->nfdofs_ * d + j) =
                  shapesface_->tangent(q, d);
            }
          }
        }
      }
    }

    faceVec.Multiply('N', 'N', 1.0, transformatrix, tempVec, 0.0);

    // Filling the vector of trace values
    for (unsigned int d = 0; d < nsd_ - 1; ++d)
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // remember that "f" is an iterator index and therefore we are
        // cycling through all the faces and all the entries of elevec1
        // except for the first one where we will put the pressure average
        elevec1(f * shapesface_->nfdofs_ * (nsd_ - 1) + d * shapesface_->nfdofs_ + i) =
            faceVec(d * shapesface_->nfdofs_ + i);
      }
  }

  return 0;
}
/*----------------------------------------------------------------------*
 * ProjectDirichField
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ProjectDirichField(
    DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1)
{
  // This routine is needed to project a field as boundary condition while for now it only gives
  // a zero value.
  /*
  shapes_.Evaluate(*ele);

  // in case this paramter "faceconsider" is set, we are applying Dirichlet values
  // and have to evaluate the trace field for the given face!
  Teuchos::Array<int> *functno = params.getPtr<Teuchos::Array<int> >("funct");
  const unsigned int *faceConsider = params.getPtr<unsigned int>("faceconsider");
  const double time = params.get<double>("time");

  //DRT::UTILS::ShapeValuesFaceParams
  svfparams(ele->Faces()[*faceConsider]->Degree(),shapes_.usescompletepoly_,2 *
  ele->Faces()[*faceConsider]->Degree());
  //shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
  shapesface_->EvaluateFace(*ele, *faceConsider);

  //We ahave to check if the face is perpendicular to some axis, if it is, we
  //have to reduce the projection system to avoid encountering a floating point
  //exception. To do that, we check if any of the normal vectors has module 1,
  //and if it does, it means that the face is perpendicular to that axis.
  //At this point we reduce the 3D projection problem to a 2D one.

  unsigned int direction = UINT_MAX;
  for (unsigned int d = 0; d < nsd_; d++)
    if (fabs(shapesface_->normals(d,0)) > 1.0 - 1e-1)
    {
      direction = d;
      break;
    }
  //Initialize the matrix to contain the value of the field given by the
  //boundary condition
  Epetra_SerialDenseVector trace(nsd_);

  Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_ * nsd_, shapesface_->nfdofs_ * nsd_);
  Epetra_SerialDenseVector trVec(shapesface_->nfdofs_ * nsd_);

  //For every quadrature point
  for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
  {
    //Saving the jacobian in the point and the coordinate in real space
    const double fac = shapesface_->jfac(q);
    LINALG::Matrix<nsd_,1> xyz;
    for (unsigned int d = 0; d < nsd_; ++d)
      xyz(d) = shapesface_->xyzreal(d, q);

    //It is necessary to write a function that takes the time as argument and
    //also the same function should accept the component of the initial
    //function to be changed.
    EvaluateAll((*functno)[0], time, xyz, trace);

    // now fill the components in the mass matrix and the right hand side
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    {
      // mass matrix
      //Be careful because this is meant to act on the tangential component
      //of the electric field in such a wahy that the tangential components
      //will assume the value specified by the boundary conditions
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
          if (d != direction)
          {
            //This is the normal projection of a field
            mass(shapesface_->nfdofs_ * d + i,shapesface_->nfdofs_ * d + j) +=
  shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;
            //mass(shapesface_->nfdofs_ * d + i,shapesface_->nfdofs_ * ((d+1)%nsd_) + j) -=
  shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac * shapesface_->normals(((d+2)%nsd_),
  q);
            //mass(shapesface_->nfdofs_ * d + i,shapesface_->nfdofs_ * ((d+2)%nsd_) + j) +=
  shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac * shapesface_->normals(((d+1)%nsd_),
  q);
          }
          else if (d == direction && i == j)
            mass(shapesface_->nfdofs_ * d + i,shapesface_->nfdofs_ * d + j) = 1.0;
        if(d!= direction)
          trVec(shapesface_->nfdofs_ * d + i) += shapesface_->shfunct(i, q) * trace(d) * fac;
        else
        trVec(shapesface_->nfdofs_ * d + i) = 0;
      }
    }
  }// for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)

  Epetra_SerialDenseSolver inverseMass;
  inverseMass.SetMatrix(mass);
  inverseMass.SetVectors(trVec, trVec);
  inverseMass.Solve();
  */
  for (unsigned int i = 0; i < shapesface_->nfdofs_ * (nsd_ - 1); ++i) elevec1(i) = 0.0;
  // elevec1(i) = trVec(i);

  return 0;
}

/*----------------------------------------------------------------------*
 * EvaluateAll
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::EvaluateAll(const int start_func,
    const double t, const LINALG::Matrix<nsd_, 1>& xyz, Epetra_SerialDenseVector& v) const
{
  int numComp = DRT::Problem::Instance()->Funct(start_func - 1).NumberComponents();

  // If there is on component for each entry of the vector use une for each
  if (numComp == v.M())
  {
    for (int d = 0; d < v.M(); ++d)
      v[d] = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(d, xyz.A(), t);
  }
  // If the vector is half the number of the component only use the firt half
  else if (numComp == 2 * v.M())
  {
    for (int d = 0; d < v.M(); ++d)
      v[d] = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(d, xyz.A(), t);
  }
  // If the number of component is half of the vector, repeat the first half twice
  else if (numComp == v.M() / 2)
  {
    for (int d = 0; d < v.M(); ++d)
      v[d] = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(d % numComp, xyz.A(), t);
  }
  // If there is only one component always use it
  else if (numComp == 1)
  {
    for (int d = 0; d < v.M(); ++d)
      v[d] = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(0, xyz.A(), t);
  }
  // If the number is not recognised throw an error
  else
    dserror(
        "Supply ONE component for your start function or NUMDIM, not anything else! With NUMDIM "
        "components the field will be initialized componentwise, if only one component is "
        "provided, every component of the field will be initialized with the same values.");
  return;
}

/*----------------------------------------------------------------------*
 | InterpolateSolutionToNodes                          berardocco 04/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ElemagEleCalc<distype>::InterpolateSolutionToNodes(DRT::ELEMENTS::Elemag* ele,
    DRT::Discretization& discretization, Epetra_SerialDenseVector& elevec1)
{
  InitializeShapes(ele);

  // Check if the vector has the correct size
  dsassert(elevec1.M() == (int)nen_ * (3 * nsd_), "Vector does not have correct size");

  // Getting the connectivity matrix
  // Contains the (local) coordinates of the nodes belonging to the element
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  // This vector will contain the values of the shape functions computed in a
  // certain coordinate. In fact the lenght of the vector is given by the number
  // of shape functions, that is the same of the number of degrees of freedom of
  // an element.
  Epetra_SerialDenseVector values(shapes_->ndofs_);

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

    // compute values for interior unknown by summing over all basis functions
    for (unsigned int d = 0; d < nsd_; ++d)
    {
      double sum_electric = 0.0;
      double sum_magnetic = 0.0;
      // Cycling through all the shape functions
      for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
      {
        // The overall value in the chosen point is given by the sum of the
        // values of the shape functions multiplied by their coefficients.
        sum_electric += values(k) * interiorElectricnp_[d * shapes_->ndofs_ + k];
        sum_magnetic += values(k) * interiorMagneticnp_[d * shapes_->ndofs_ + k];
      }
      // sum contains the linear combination of the shape functions times the
      // coefficients and its values are reordered in elevec1 grouped by
      // component: the first component for every node, then the following
      // component for the same nodes and so on for every component.
      elevec1(d * nen_ + i) = sum_electric;
      elevec1(nen_ * nsd_ + d * nen_ + i) = sum_magnetic;
    }
  }

  // get trace solution values
  // Same as before bu this time the dimension is nsd_-1 because we went from
  // the interior to the faces. We have to be careful because we are using a
  // part of the previous vector. The coordinates are still in the local frame.
  locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(
      DRT::UTILS::DisTypeToFaceShapeType<distype>::shape);

  // Storing the number of nodes for each face of the element as vector
  // NumberCornerNodes
  std::vector<int> ncn = DRT::UTILS::getNumberOfFaceElementCornerNodes(distype);
  // NumberInternalNodes
  std::vector<int> nin = DRT::UTILS::getNumberOfFaceElementInternalNodes(distype);

  // Cycling the faces of the element
  Epetra_SerialDenseVector fvalues(shapesface_->nfdofs_);
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    // Checking how many nodes the face has
    const int nfn = DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace;

    shapesface_->EvaluateFace(*ele, f);

    Epetra_SerialDenseVector facetrace((nsd_ - 1) * shapesface_->nfdofs_);
    Epetra_SerialDenseVector temptrace(nsd_ * shapesface_->nfdofs_);

    // As already said, the dimension of the coordinate matrix is now nsd_-1
    // times the number of nodes in the face.
    LINALG::Matrix<nsd_ - 1, nfn> xsishuffle(true);

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
    Epetra_SerialDenseMatrix transformatrix(
        (nsd_ - 1) * shapesface_->nfdofs_, nsd_ * shapesface_->nfdofs_);
    Epetra_SerialDenseMatrix inv_transformatrix(
        nsd_ * shapesface_->nfdofs_, (nsd_ - 1) * shapesface_->nfdofs_);
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    {
      for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
      {
        for (unsigned int d = 0; d < nsd_ - 1; ++d)
        {
          for (unsigned int q = 0; q < nsd_; ++q)
          {
            // I need tangents because I'm translating real coordinates to face ones
            if (i == j)
            {
              transformatrix(shapesface_->nfdofs_ * d + i, shapesface_->nfdofs_ * q + j) =
                  shapesface_->tangents(d, q);
              inv_transformatrix(shapesface_->nfdofs_ * q + i, shapesface_->nfdofs_ * d + j) =
                  shapesface_->tangent(q, d);
            }
          }
        }
      }
    }
    // Storing the face part of the trace vector
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    {
      for (unsigned int d = 0; d < nsd_ - 1; ++d)
      {
        facetrace(shapesface_->nfdofs_ * d + i) =
            ele->elenodeTrace2d_[f * (nsd_ - 1) * shapesface_->nfdofs_ + shapesface_->nfdofs_ * d +
                                 i];
      }
    }

    temptrace.Multiply('N', 'N', 1.0, inv_transformatrix, facetrace, 0);

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
          elevec1((nsd_ * 2 + d) * nen_ + shapesface_->faceNodeOrder[f][i]) += sum / nsd_;
        }
        else if (i < nfn - nin[f])
        {
          elevec1((nsd_ * 2 + d) * nen_ + shapesface_->faceNodeOrder[f][i]) += sum / (nsd_ - 1);
        }
        else
        {
          elevec1((nsd_ * 2 + d) * nen_ + shapesface_->faceNodeOrder[f][i]) += sum;
        }
      }
    }
  }
  return 0;
}

/*----------------------------------------------------------------------*
 * Instance
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ElemagEleCalc<distype>* DRT::ELEMENTS::ElemagEleCalc<distype>::Instance(bool create)
{
  static ElemagEleCalc<distype>* instance;
  if (create)
  {
    if (instance == NULL) instance = new ElemagEleCalc<distype>();
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 * Done
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}

/*----------------------------------------------------------------------*
 * Constructor LocalSolver
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::LocalSolver(const DRT::ELEMENTS::Elemag* ele,
    DRT::UTILS::ShapeValues<distype>& shapeValues,
    Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype>>& shapeValuesFace,
    INPAR::ELEMAG::DynamicType& dyna)
    : ndofs_(shapeValues.ndofs_), shapes_(shapeValues), shapesface_(shapeValuesFace), dyna_(dyna)
{
  // shape all matrices
  // Each one of these matrices is related to one equation of the formulation,
  // therefore ndofs equations in FEM terms) and one variable.
  // The number of entries is then given by ndofs time sthe dimension of the
  // space where the unknown lies. For vectorial field nsd_ gives the dimension.
  reshapeMatrixIfNecessary(Amat, nsd_ * ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(invAmat, nsd_ * ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(Cmat, nsd_ * ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(Emat, nsd_ * ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(Fmat, nsd_ * ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(Gmat, nsd_ * ndofs_, nsd_ * ndofs_);
  // These matrices have a "strange" shape because to merge them there will be
  // applied a matrix multiplication between the first one and the transposed
  // second one. The shape of the resulting matrix will therefore be ndofs x ndofs.
  reshapeMatrixIfNecessary(massMat, ndofs_, ndofs_);
  reshapeMatrixIfNecessary(massPart, ndofs_, shapeValues.nqpoints_);
  reshapeMatrixIfNecessary(massPartW, ndofs_, shapeValues.nqpoints_);

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
  // Dmat and Hmat are the matrix that belongs to the equation for the magnetic
  // and electric field but multiply the hybrid variable, therefore their dimensions are:
  // o) nsd_*ndofs_ x onfdofs
  reshapeMatrixIfNecessary(Dmat, nsd_ * ndofs_, onfdofs);
  reshapeMatrixIfNecessary(Hmat, nsd_ * ndofs_, onfdofs);
  // Matrices Imat and Jmat describe the part of the continuity condition that
  // multiply the electric and magnetic fields and therefore their dimensions are:
  // o) ondofs x nsd_*ndofs_
  reshapeMatrixIfNecessary(Imat, onfdofs, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(Jmat, onfdofs, nsd_ * ndofs_);
  // Finally Jmat is the matrix that belongs to the continuity condition and
  // multiplies the hybrid variable and therefore its dimensions are:
  // o) ondofs x ondofs
  reshapeMatrixIfNecessary(Lmat, onfdofs, onfdofs);
}

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::UpdateInteriorVariablesAndComputeResidual(
    Teuchos::ParameterList& params, DRT::ELEMENTS::Elemag& ele,
    const Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseVector& elevec, double dt,
    bool errormaps, bool updateonly)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "DRT::ELEMENTS::ElemagEleCalc::UpdateInteriorVariablesAndComputeResidual");

  // *****************************************************
  // update interior variables first
  // *****************************************************

  Epetra_SerialDenseVector tempVec1(shapes_->ndofs_ * nsd_);
  Epetra_SerialDenseVector tempVec2(shapes_->ndofs_ * nsd_);
  Epetra_SerialDenseVector tempVec3(shapes_->ndofs_ * nsd_);
  Epetra_SerialDenseVector xVec(shapes_->ndofs_ * nsd_);
  Epetra_SerialDenseVector yVec(shapes_->ndofs_ * nsd_);
  Epetra_SerialDenseMatrix tempMat(shapes_->ndofs_ * nsd_, shapes_->ndofs_ * nsd_);
  Epetra_SerialDenseMatrix tempMat2(shapes_->ndofs_ * nsd_, shapes_->ndofs_ * nsd_);

  localSolver_->ComputeSource(params, tempVec2, xVec);

  tempVec1.Multiply('N', 'N', 1.0, localSolver_->Amat, ele.eleinteriorMagnetic_, 0.0);  //  AH^{n−1}
  //  h = AH^{n−1} - D\lambda^{n-1}
  tempVec1.Multiply('N', 'N', -1.0, localSolver_->Dmat, ele.elenodeTrace2d_, 1.0);
  tempMat.Multiply('N', 'N', 1.0, localSolver_->Fmat, localSolver_->invAmat, 0.0);  // FA^{-1}

  tempMat2 += localSolver_->Emat;
  tempMat2 += localSolver_->Gmat;
  tempMat2.Multiply('N', 'N', -1.0, tempMat, localSolver_->Cmat, 1.0);  //(E + G) - FA^{-1}C
  {
    Epetra_SerialDenseSolver invert;
    invert.SetMatrix(tempMat2);
    invert.Invert();  //  [(E + G) - FA^{-1}C]^{-1}
  }
  // EE^{n-1} - I_s
  tempVec2.Multiply('N', 'N', 1.0, localSolver_->Emat, ele.eleinteriorElectric_, -1.0);
  // e = EE^{n-1} - I_s - H\lambda^{n-1}
  tempVec2.Multiply('N', 'N', -1.0, localSolver_->Hmat, ele.elenodeTrace2d_, 1.0);

  tempVec2.Multiply('N', 'N', -1.0, tempMat, tempVec1, 1.0);  //  e - FA^{-1}h

  //  E^{n} = [(E + G) - FA^{-1}C]^{-1} (e - FA^{-1}h)
  interiorElectricnp_.Multiply('N', 'N', 1.0, tempMat2, tempVec2, 0.0);


  tempVec1.Multiply('N', 'N', -1.0, localSolver_->Cmat, interiorElectricnp_, 1.0);  //  h - CE^{n}

  //  = A^{-1}(h - CE^{n})
  interiorMagneticnp_.Multiply('N', 'N', 1.0, localSolver_->invAmat, tempVec1, 0.0);

  ele.eleinteriorMagnetic_ = interiorMagneticnp_;
  ele.eleinteriorElectric_ = interiorElectricnp_;

  // Updateresidual

  //  = -FH^{n} - I_s
  xVec.Multiply('N', 'N', -1.0, localSolver_->Fmat, ele.eleinteriorMagnetic_, -1.0);
  //  = EE^{n} -I_s - FH^{n}
  xVec.Multiply('N', 'N', 1.0, localSolver_->Emat, ele.eleinteriorElectric_, 1.0);
  //  = [(E + G) - FA^{-1}C]^{-1}(EE^{n} - I_s^{n} - FH^{n})
  yVec.Multiply('N', 'N', 1.0, tempMat2, xVec, 0.0);

  elevec.Multiply('N', 'N', 1.0, localSolver_->Jmat, yVec, 0.0);  //  = Jy

  xVec.Multiply('N', 'N', 1.0, localSolver_->Amat, ele.eleinteriorMagnetic_, 0.0);  //  = AH^{n}
  xVec.Multiply('N', 'N', -1.0, localSolver_->Cmat, yVec, 1.0);    //  = AH^{n} - Cy
  yVec.Multiply('N', 'N', 1.0, localSolver_->invAmat, xVec, 0.0);  //  = A^{-1}(AH^{n} - Cy)

  elevec.Multiply('N', 'N', -1.0, localSolver_->Imat, yVec, -1.0);  //  = -Ix -Jy

  return;
}  // UpdateInteriorVariablesAndComputeResidual

/*----------------------------------------------------------------------*
 * ComputeAbsorbingBC
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ComputeAbsorbingBC(
    DRT::Discretization& discretization, DRT::ELEMENTS::Elemag* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, int face, Epetra_SerialDenseMatrix& elemat, int indexstart,
    Epetra_SerialDenseVector& elevec1)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagEleCalc::ComputeAbsorbingBC");

  shapesface_->EvaluateFace(*ele, face);

  // Get the user defined functions
  Teuchos::RCP<DRT::Condition>* cond = params.getPtr<Teuchos::RCP<DRT::Condition>>("condition");
  const std::vector<int>* funct = (*cond)->Get<std::vector<int>>("funct");
  const double time = params.get<double>("time");

  Epetra_SerialDenseVector tempVec1(ndofs_ * nsd_);
  Epetra_SerialDenseVector tempVec2(shapesface_->nfdofs_ * (nsd_ - 1));
  // the RHS matrix has to have the row dimension equal to the number of shape
  // functions(so we have one coefficient for each) and a number of column
  // equal to the overall number of component that we want to solve for.
  // The number is nsd_*2 because we have two fields..
  Epetra_SerialDenseMatrix localMat(ndofs_, nsd_ * 2);
  {
    for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
    {
      // Storing the values of the coordinates for the current quadrature point
      // and of the jacobian computed in that point
      const double fac = shapes_.jfac(q);
      LINALG::Matrix<nsd_, 1> xyz;
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz(d) = shapes_.xyzreal(d, q);  // coordinates of quadrature point in real coordinates
      // Creating the temporary electric and magnetic field vector intVal
      // The vector is going to contain first the electric and then the magnetic
      // field such that the field will be initialized as first tree component
      // of the specified function as electric field, last three components as
      // magnetic field. If there is only one component all the components will
      // be initialized to the same value.
      Epetra_SerialDenseVector intVal(2 * nsd_);
      dsassert(funct != NULL, "funct not set for initial value");
      EvaluateAll((*funct)[0], time, xyz, intVal);
      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < ndofs_; ++i)
      {
        // Mass matrix
        massPart(i, q) = shapes_.shfunct(i, q);
        massPartW(i, q) = shapes_.shfunct(i, q) * fac;

        // RHS for the electric and magnetic field
        for (int j = 0; j < intVal.M(); ++j)
          localMat(i, j) += shapes_.shfunct(i, q) * intVal(j) * fac;
      }
    }
    // The integration is made by computing the matrix product
    massMat.Multiply('N', 'T', 1., massPart, massPartW, 0.);
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(massMat);
      inverseMass.SetVectors(localMat, localMat);
      inverseMass.Solve();
    }
  }

  for (unsigned int r = 0; r < ndofs_; ++r)
    for (unsigned int d = 0; d < nsd_; ++d)
      tempVec1(d * ndofs_ + r) = localMat(r, d);  // Electric field

  // Creating the matrix
  Epetra_SerialDenseMatrix transformatrix(
      (nsd_ - 1) * shapesface_->nfdofs_, nsd_ * shapesface_->nfdofs_);
  Epetra_SerialDenseMatrix inv_transformatrix(
      nsd_ * shapesface_->nfdofs_, (nsd_ - 1) * shapesface_->nfdofs_);
  for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
    for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
      for (unsigned int d = 0; d < nsd_ - 1; ++d)
        for (unsigned int q = 0; q < nsd_; ++q)
          if (i == j)
          {
            // I need tangents because I'm translating real coordinates to face ones
            transformatrix(shapesface_->nfdofs_ * d + i, shapesface_->nfdofs_ * q + j) =
                shapesface_->tangents(d, q);
            inv_transformatrix(shapesface_->nfdofs_ * q + i, shapesface_->nfdofs_ * d + j) =
                shapesface_->tangent(q, d);
          }

  const MAT::ElectromagneticMat* actmat = static_cast<const MAT::ElectromagneticMat*>(mat.get());
  double impedance = sqrt(actmat->epsilon(ele->Id()) / actmat->mu(ele->Id()));

  // MIXED SHAPE FUNCTIONS
  // The matrix that are going to be build here are D,I and J
  // loop over number of internal shape functions
  // Here we need to create only the first part of tghe D and H matrix to be multiplied by the
  // transformation matrices and then put in the real D and H matrices
  Epetra_SerialDenseMatrix tempI(shapesface_->nfdofs_ * nsd_, ndofs_ * nsd_);
  Epetra_SerialDenseMatrix tempJ(shapesface_->nfdofs_ * nsd_, ndofs_ * nsd_);
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
            const double temp =
                shapesface_->jfac(q) * shapesface_->shfunctI(i, q) * shapesface_->shfunct(j, q);

            // Filling the matrices
            // Imat = -[H x n]
            //+1 coordinate
            tempI(shapesface_->nfdofs_ * d + j, ((d + 1) % nsd_) * ndofs_ + i) -=
                temp * shapesface_->normals((d + 2) % nsd_, q);
            //+2 coordinate
            tempI(shapesface_->nfdofs_ * d + j, ((d + 2) % nsd_) * ndofs_ + i) +=
                temp * shapesface_->normals((d + 1) % nsd_, q);
            // Jmat
            tempJ(shapesface_->nfdofs_ * d + j, d * shapesface_->nfdofs_ + i) += temp;
          }
        }  // for (unsigned int d = 0; d < nsd_; ++d)
      }    // for (unsigned int j=0; j<ndofs_; ++j)
    }      // if( shapesface_->shfunctI.NonzeroOnFace(i) )
  }        // for (unsigned int i = 0; i < ndofs_; ++i)

  // Fill face values into the matrices
  Epetra_SerialDenseMatrix magneticMat(shapesface_->nfdofs_ * (nsd_ - 1), ndofs_ * nsd_);
  Epetra_SerialDenseMatrix electricMat(shapesface_->nfdofs_ * (nsd_ - 1), ndofs_ * nsd_);
  magneticMat.Multiply('N', 'N', 1.0, transformatrix, tempI, 0.0);
  electricMat.Multiply('N', 'N', 1.0, transformatrix, tempJ, 0.0);

  tempVec2.Multiply('N', 'N', impedance, electricMat, tempVec1, 0.0);

  for (unsigned int r = 0; r < ndofs_; ++r)
    for (unsigned int d = 0; d < nsd_; ++d)
      tempVec1(d * ndofs_ + r) = localMat(r, d + nsd_);  // magnetic

  tempVec2.Multiply('N', 'N', 1.0, magneticMat, tempVec1, 1.0);

  unsigned int newindex = shapesface_->nfdofs_ * (nsd_ - 1) * face;

  for (int i = 0; i < tempVec2.M(); ++i) elevec1(newindex + i) = tempVec2(i);

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

  return;
}

/*----------------------------------------------------------------------*
 * ComputeSource
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ComputeSource(
    Teuchos::ParameterList& params, Epetra_SerialDenseVector& interiorSourcen,
    Epetra_SerialDenseVector& interiorSourcenp)
{
  int funcno = params.get<int>("sourcefuncno");
  if (funcno <= 0) return;  // there is no such thing as a volume force
  // if(DRT::Problem::Instance()->Funct(funcno).NumberComponents()>1) dserror("for standard
  // elemag, the source term has to be scalar, i.e. only 1 component");

  // the vector to be filled
  Epetra_SerialDenseVector sourcen(nsd_);
  Epetra_SerialDenseVector sourcenp(nsd_);

  // what time is it?
  double tn = params.get<double>("time");
  double tp = params.get<double>("timep");

  // double f_value = 0.0;
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    LINALG::Matrix<nsd_, 1> xyz;
    for (unsigned int d = 0; d < nsd_; ++d) xyz(d) = shapes_.xyzreal(d, q);

    // calculate right hand side contribution for dp/dt
    EvaluateAll(funcno, tn, xyz, sourcen);
    EvaluateAll(funcno, tp, xyz, sourcenp);

    // add it all up
    for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        interiorSourcen(i + d * shapes_.ndofs_) +=
            shapes_.shfunct(i, q) * sourcen(d) * shapes_.jfac(q);
        interiorSourcenp(i + d * shapes_.ndofs_) +=
            shapes_.shfunct(i, q) * sourcenp(d) * shapes_.jfac(q);
      }
  }

  return;
}

/*----------------------------------------------------------------------*
 * ComputeInteriorMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ComputeInteriorMatrices(
    double dt, double sigma, double mu, double epsilon)
{
  // The definitions of the matrices created here can be found in the internal
  // paper from Gravemeier "A hybridizable discontinous Galerkin method for
  // electromagnetics in subsurface applications".
  // The explicit form of these matrices is reported for convenience?
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagEleCalc::ComputeInteriorMatrices");
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

  Epetra_SerialDenseMatrix tmpMat(ndofs_, ndofs_);
  // this temorary matrix is used to compute the numerical integration and the
  // values are then copied in the right places. Probably it is also possible
  // to have the matrix multiplication to obtain directly the correct matrices
  // but it would mean to compute three time sthe same value for each shape
  // function instead of computing it only omnce and then directly copying it.
  tmpMat.Multiply('N', 'T', 1.0, massPart, massPartW, 0.0);
  // A, E and part of G
  for (unsigned int j = 0; j < ndofs_; ++j)
    for (unsigned int i = 0; i < ndofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        Amat(d * ndofs_ + i, d * ndofs_ + j) = mu * tmpMat(i, j);
        Emat(d * ndofs_ + i, d * ndofs_ + j) = epsilon * tmpMat(i, j);
        // Carefull because G is not yet complete, it is necessary to add the boundary
        Gmat(d * ndofs_ + i, d * ndofs_ + j) = sigma * tmpMat(i, j);
      }

  Amat.Scale(1 / dt);
  Emat.Scale(1 / dt);

  {  // We are creating this scope to destroy everything related to the matrix inversion
    // We are going to need both A and its inverse and therefore we are storing both
    invAmat += Amat;
    Epetra_SerialDenseSolver invA;
    invA.SetMatrix(invAmat);
    int err = invA.Invert();
    if (err != 0) dserror("Inversion for Amat failed with errorcode %d", err);
  }

  for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
    for (unsigned int j = 0; j < shapes_.ndofs_; ++j)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
        {
          Cmat(i + d * ndofs_, j + ((d + 1) % nsd_) * ndofs_) +=
              shapes_.shderxy(i * nsd_ + ((d + 2) % nsd_), q) * shapes_.shfunct(j, q) *
              shapes_.jfac(q);
          Cmat(i + d * ndofs_, j + ((d + 2) % nsd_) * ndofs_) -=
              shapes_.shderxy(i * nsd_ + ((d + 1) % nsd_), q) * shapes_.shfunct(j, q) *
              shapes_.jfac(q);
        }
      }

  // The summation over quadrature points is manually made
  for (unsigned int i = 0; i < ndofs_; ++i)
    for (unsigned int j = 0; j < ndofs_; ++j)
      for (unsigned int d = 0; d < nsd_; ++d)
        for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
        {
          const double temp = shapes_.shfunct(i, q) * shapes_.jfac(q);
          // this can be avoided but the optimization of the code comes later
          Fmat(d * ndofs_ + i, ((d + 1) % nsd_) * ndofs_ + j) +=
              temp * shapes_.shderxy(j * nsd_ + ((d + 2) % nsd_), q);
          Fmat(d * ndofs_ + i, ((d + 2) % nsd_) * ndofs_ + j) -=
              temp * shapes_.shderxy(j * nsd_ + ((d + 1) % nsd_), q);
        }

  return;
}

/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ComputeResidual(
    Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec,
    Epetra_SerialDenseVector& interiorMagneticn, Epetra_SerialDenseVector& interiorElectricn)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagEleCalc::ComputeResidual");

  // for implicit Euler
  //                                -1
  //                     +---------+    +-----------+
  //                     | A    C  |    | A H^{n-1} |
  // R^{n}  = - [ I  J ] |         |    |           |  =  -Ix - Jy
  //                     | F   E+G |    | E E^{n-1} |
  //                     +---------+    +-----------+
  //
  //  x = A^{-1} (AH^{n-1} - Cy)
  //
  //  y = ((E + G) - F A^{-1} C)^{-1} ((E E^{n-1} - I_s^{n-1}) - F A^{-1} A H^{n-1})

  const unsigned int intdofs = ndofs_ * nsd_;
  // All the vectors are initilized to zero
  Epetra_SerialDenseVector tempVec1(intdofs);
  Epetra_SerialDenseVector tempVec2(intdofs);
  Epetra_SerialDenseVector tempVec3(intdofs);
  // Once the compute source is ready we will need to delete these
  // The ComputeSource is necesessary to include the forcing terms
  ComputeSource(params, tempVec2, tempVec3);

  // The last -1.0 in the following function has not been removed such that once
  // the ComputeSource() has been created there will be no need to change it
  tempVec1.Multiply('N', 'N', 1.0, Amat, interiorMagneticn, 0.0);   // AH^{n-1}
  tempVec2.Multiply('N', 'N', 1.0, Emat, interiorElectricn, -1.0);  // E E^{n-1} - I_s^{n-1}

  Epetra_SerialDenseMatrix tempMat1(intdofs, intdofs);
  tempMat1.Multiply('N', 'N', 1.0, Fmat, invAmat, 0.0);  // F A^{-1}
  tempVec2.Multiply(
      'N', 'N', -1.0, tempMat1, tempVec1, 1.0);  // ((E E^{n-1} - I_s^{n-1}) - F A^{-1} A H^{n-1})

  Epetra_SerialDenseMatrix tempMat2(intdofs, intdofs);
  // Gmat already contains Emat in it
  tempMat2 += Emat;
  tempMat2 += Gmat;
  tempMat2.Multiply('N', 'N', -1.0, tempMat1, Cmat, 1.0);  // = (E + G) - F A^{-1} C
  {
    Epetra_SerialDenseSolver inverseinW;
    inverseinW.SetMatrix(tempMat2);
    int err = inverseinW.Invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d", err);
  }
  // tempMat2 = ((E + G) - F A^{-1} C)^{-1}

  tempVec3.Multiply('N', 'N', 1.0, tempMat2, tempVec2, 0.0);  // y
  elevec.Multiply('N', 'N', -1.0, Jmat, tempVec3, 0.0);       //  -Jy

  tempVec1.Multiply('N', 'N', -1.0, Cmat, tempVec3, 1.0);    // AH^{n-1} - Cy
  tempVec3.Multiply('N', 'N', 1.0, invAmat, tempVec1, 0.0);  //  x = A^{-1} (AH^{n-1} - Cy)
  elevec.Multiply('N', 'N', -1.0, Imat, tempVec3, 1.0);      //  -Ix - Jy

  return;
}  // ComputeResidual

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ComputeFaceMatrices(
    const int face, double dt, int indexstart, int newindex, double sigma, double mu)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagEleCalc::ComputeFaceMatrices");

  // Tau is defined as (\frac{|sigma|}{\mu t_c})^0.5 where t_c is a
  // characteristic time scale
  const double tau = 100000.0;  // sqrt(sigma/mu/dt);

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
  Epetra_SerialDenseMatrix transformatrix(
      (nsd_ - 1) * shapesface_->nfdofs_, nsd_ * shapesface_->nfdofs_);
  Epetra_SerialDenseMatrix inv_transformatrix(
      nsd_ * shapesface_->nfdofs_, (nsd_ - 1) * shapesface_->nfdofs_);
  for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
  {
    for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
    {
      for (unsigned int d = 0; d < nsd_ - 1; ++d)
      {
        for (unsigned int q = 0; q < nsd_; ++q)
        {
          // I need tangents because I'm translating real coordinates to face ones
          if (i == j)
          {
            transformatrix(shapesface_->nfdofs_ * d + i, shapesface_->nfdofs_ * q + j) =
                shapesface_->tangents(d, q);
            inv_transformatrix(shapesface_->nfdofs_ * q + i, shapesface_->nfdofs_ * d + j) =
                shapesface_->tangent(q, d);
          }
        }
      }
    }
  }

  // MIXED SHAPE FUNCTIONS
  // The matrix that are going to be build here are D,I and J
  // loop over number of internal shape functions
  // Here we need to create only the first part of tghe D and H matrix to be multiplied by the
  // transformation matrices and then put in the real D and H matrices
  Epetra_SerialDenseMatrix tempD(ndofs_ * nsd_, shapesface_->nfdofs_ * nsd_);
  Epetra_SerialDenseMatrix tempH(ndofs_ * nsd_, shapesface_->nfdofs_ * nsd_);
  Epetra_SerialDenseMatrix tempI(shapesface_->nfdofs_ * nsd_, ndofs_ * nsd_);
  Epetra_SerialDenseMatrix tempJ(shapesface_->nfdofs_ * nsd_, ndofs_ * nsd_);
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
            // Dmat
            //+1 coordinate
            tempD(d * ndofs_ + i, shapesface_->nfdofs_ * ((d + 2) % nsd_) + j) +=
                temp2 * shapesface_->normals(((d + 1) % nsd_), q);
            //+2 coordinate
            tempD(d * ndofs_ + i, shapesface_->nfdofs_ * ((d + 1) % nsd_) + j) -=
                temp2 * shapesface_->normals(((d + 2) % nsd_), q);
            // Hmat
            // 0 coordinate
            tempH(d * ndofs_ + i, shapesface_->nfdofs_ * d + j) +=
                temp * (pow(shapesface_->normals(((d + 1) % nsd_), q), 2) +
                           pow(shapesface_->normals(((d + 2) % nsd_), q), 2));
            //+1 coordinate
            tempH(d * ndofs_ + i, shapesface_->nfdofs_ * ((d + 1) % nsd_) + j) -=
                temp * shapesface_->normals(d, q) * shapesface_->normals(((d + 1) % nsd_), q);
            //+2 coordinate
            tempH(d * ndofs_ + i, shapesface_->nfdofs_ * ((d + 2) % nsd_) + j) -=
                temp * shapesface_->normals(d, q) * shapesface_->normals(((d + 2) % nsd_), q);
            // Imat
            //+1 coordinate
            tempI(shapesface_->nfdofs_ * d + j, ((d + 1) % nsd_) * ndofs_ + i) -=
                temp2 * shapesface_->normals((d + 2) % nsd_, q);
            //+2 coordinate
            tempI(shapesface_->nfdofs_ * d + j, ((d + 2) % nsd_) * ndofs_ + i) +=
                temp2 * shapesface_->normals((d + 1) % nsd_, q);
            // Jmat
            // Own coordinate
            // Jmat(shapesface_->nfdofs_*d + j, d * ndofs_ + i) += temp2;
            tempJ(shapesface_->nfdofs_ * d + j, d * ndofs_ + i) +=
                temp * (pow(shapesface_->normals((d + 1) % nsd_, q), 2) +
                           pow(shapesface_->normals((d + 2) % nsd_, q), 2));
            //+1 coordinate
            tempJ(shapesface_->nfdofs_ * d + j, ((d + 1) % nsd_) * ndofs_ + i) -=
                temp * shapesface_->normals(d, q) * shapesface_->normals((d + 1) % nsd_, q);
            //+2 coordinate
            tempJ(shapesface_->nfdofs_ * d + j, ((d + 2) % nsd_) * ndofs_ + i) -=
                temp * shapesface_->normals(d, q) * shapesface_->normals((d + 2) % nsd_, q);
          }
        }  // for (unsigned int d = 0; d < nsd_; ++d)
      }    // for (unsigned int j=0; j<ndofs_; ++j)
    }      // if( shapesface_->shfunctI.NonzeroOnFace(i) )
  }        // for (unsigned int i = 0; i < ndofs_; ++i)

  // Fill face values into the matrices
  {
    Epetra_SerialDenseMatrix tempMat1(ndofs_ * nsd_, shapesface_->nfdofs_ * (nsd_ - 1));
    Epetra_SerialDenseMatrix tempMat2(ndofs_ * nsd_, shapesface_->nfdofs_ * (nsd_ - 1));
    Epetra_SerialDenseMatrix tempMat3(shapesface_->nfdofs_ * (nsd_ - 1), ndofs_ * nsd_);
    Epetra_SerialDenseMatrix tempMat4(shapesface_->nfdofs_ * (nsd_ - 1), ndofs_ * nsd_);
    tempMat1.Multiply('N', 'N', 1.0, tempD, inv_transformatrix, 0.0);
    tempMat2.Multiply('N', 'N', 1.0, tempH, inv_transformatrix, 0.0);
    tempMat3.Multiply('N', 'N', 1.0, transformatrix, tempI, 0.0);
    tempMat4.Multiply('N', 'N', 1.0, transformatrix, tempJ, 0.0);

    for (unsigned int i = 0; i < ndofs_ * nsd_; ++i)
      for (unsigned int j = 0; j < shapesface_->nfdofs_ * (nsd_ - 1); ++j)
      {
        Dmat(i, newindex + j) = tempMat1(i, j);
        Hmat(i, newindex + j) = tempMat2(i, j);
        Imat(newindex + j, i) = tempMat3(j, i);
        Jmat(newindex + j, i) = tempMat4(j, i);
      }
  }


  // BOUNDARY SHAPE FUNCTIONS
  // Epetra_SerialDenseMatrix tempL(shapesface_->nfdofs_ * (nsd_-1), shapesface_->nfdofs_ *
  // (nsd_-1));
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
          Lmat(newindex + shapesface_->nfdofs_ * d + i, newindex + shapesface_->nfdofs_ * d + j) -=
              temp;  //* (pow(shapesface_->normals(((d+1)%nsd_), q),2) +
                     // pow(shapesface_->normals(((d+2)%nsd_), q),2));
        }
      }
    }  // for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
  }    // for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)

  // INTERIOR SHAPE FUNCTIONS
  // Some terms are still missing in G!!
  for (unsigned int i = 0; i < ndofs_; ++i)
  {
    for (unsigned int j = 0; j < ndofs_; ++j)
    {
      if (shapesface_->shfunctI.NonzeroOnFace(i) && shapesface_->shfunctI.NonzeroOnFace(j))
      {
        for (unsigned int d = 0; d < nsd_; ++d)
        {
          for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
          {
            const double temp = tau * shapesface_->jfac(q) * shapesface_->shfunctI(i, q) *
                                shapesface_->shfunctI(j, q);
            // Gmat
            // 0 coordinate
            Gmat(d * ndofs_ + i, d * ndofs_ + j) -=
                temp * (pow(shapesface_->normals(((d + 1) % nsd_), q), 2) +
                           pow(shapesface_->normals(((d + 2) % nsd_), q), 2));
            //+1 coordinate
            Gmat(d * ndofs_ + i, ((d + 1) % nsd_) * ndofs_ + j) +=
                temp * shapesface_->normals(d, q) * shapesface_->normals(((d + 1) % nsd_), q);
            //+2 coordinate
            Gmat(d * ndofs_ + i, ((d + 2) % nsd_) * ndofs_ + j) +=
                temp * shapesface_->normals(d, q) * shapesface_->normals(((d + 2) % nsd_), q);
          }
        }
      }
    }
  }

  return;
}  // ComputeFaceMatrices


/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::CondenseLocalPart(
    Epetra_SerialDenseMatrix& eleMat)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ElemagEleCalc::CondenseLocalPart");

  // THE MATRIX
  //                             -1
  //                   +--------+    +-----+
  //                   |        |    |     |
  //                   | A   C  |    |  D  |
  //  K = L - [ I  J ] |        |    |     |  = L - I X - J Y
  //                   | F  E+G |    |  H  |
  //                   +--------+    +-----+

  //   Y = [ (E+G) - F A^{-1} C ]^{-1} [ H - F A^{-1} D]

  //   X = A^{-1} [ D - C Y ]

  const unsigned int onfdofs = eleMat.M();
  const unsigned int intdofs = ndofs_ * nsd_;

  // Thi can be useful to remember when coding
  // int 	Multiply (char TransA, char TransB, double ScalarAB, Matrix &A, Matrix &B, double
  // ScalarThis) this = ScalarThis*this + ScalarAB*A*B
  Epetra_SerialDenseMatrix tempMat1(intdofs, intdofs);
  tempMat1.Multiply('N', 'N', 1.0, Fmat, invAmat, 0.0);  // =  F A^{-1}

  Epetra_SerialDenseMatrix tempMat2(intdofs, intdofs);

  // This is E+G
  tempMat2 += Emat;  // = E
  tempMat2 += Gmat;  // = E + G

  tempMat2.Multiply('N', 'N', -1.0, tempMat1, Cmat, 1.0);  // = (E+G) - F A^{-1} C

  Epetra_SerialDenseMatrix tempMat3(intdofs, onfdofs);
  tempMat3 += Hmat;                                        // = H
  tempMat3.Multiply('N', 'N', -1.0, tempMat1, Dmat, 1.0);  // = H - F A^{-1} D

  // Inverting the first part of the Y matrix
  {
    Epetra_SerialDenseSolver inverseinW;
    inverseinW.SetMatrix(tempMat2);
    int err = inverseinW.Invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d", err);
  }
  // tempMat2 = [(E+G) - F A^{-1} C]^{-1}

  eleMat = Lmat;  // = L
  // reusing matrix that are not needed
  tempMat1.Shape(intdofs, onfdofs);
  tempMat1.Multiply(
      'N', 'N', 1.0, tempMat2, tempMat3, 0.0);  //  Y = [(E+G) - F A^{-1} C]^{-1}(H - F A^{-1} D)
  eleMat.Multiply('N', 'N', -1.0, Jmat, tempMat1, 1.0);  // = L - J Y

  tempMat2.Shape(intdofs, onfdofs);
  tempMat2 = Dmat;
  tempMat2.Multiply('N', 'N', -1.0, Cmat, tempMat1, 1.0);  // = D - C Y

  tempMat3.Shape(intdofs, onfdofs);
  tempMat3.Multiply('N', 'N', 1.0, invAmat, tempMat2, 0.0);  // = X = A^{-1} ( D - C Y )

  eleMat.Multiply('N', 'N', -1.0, Imat, tempMat3, 1.0);  // = K = L - I X - J y

  return;
}  // CondenseLocalPart

/*----------------------------------------------------------------------*
 * Compute internal and face matrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ElemagEleCalc<distype>::LocalSolver::ComputeMatrices(
    DRT::Discretization& discretization, const Teuchos::RCP<MAT::Material>& mat,
    DRT::ELEMENTS::Elemag& ele, double dt, INPAR::ELEMAG::DynamicType dyna)
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
  zeroMatrix(invAmat);
  zeroMatrix(Amat);
  zeroMatrix(Cmat);
  zeroMatrix(Dmat);
  zeroMatrix(Emat);
  zeroMatrix(Fmat);
  zeroMatrix(Gmat);
  zeroMatrix(Hmat);
  zeroMatrix(Imat);
  zeroMatrix(Jmat);
  zeroMatrix(Lmat);

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
    DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele.Faces()[face]->Degree(),
        shapes_.usescompletepoly_, 2 * ele.Faces()[face]->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    */

    // Updating face data
    shapesface_->EvaluateFace(ele, face);

    // Here are the matrices for the boundary integrals
    ComputeFaceMatrices(face, dt, sumindex, newindex, sigma, mu);
    sumindex += nsd_ * shapesface_->nfdofs_;
    newindex += (nsd_ - 1) * shapesface_->nfdofs_;
  }

  return;
}

// template classes
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::tri3>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::ElemagEleCalc<DRT::Element::nurbs27>;