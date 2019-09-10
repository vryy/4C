/*----------------------------------------------------------------------------*/
/*! \file
\brief Routines for calculation of HDG weakly compressible fluid element

\level 2

\maintainer Andrea La Spina
*/
/*----------------------------------------------------------------------------*/

#include "fluid_ele_calc_hdg_weak_comp.H"
#include "fluid_ele_calc.H"
#include "fluid_ele_parameter.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_parameter_timint.H"
#include "fluid_ele_action.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_geometry/position_array.H"

#include "../drt_fluid/fluid_functions.H"
#include "../drt_mat/fluid_weakly_compressible.H"

#include "../drt_fem_general/drt_utils_polynomial.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include <Epetra_SerialDenseSolver.h>



namespace
{
  void zeroMatrix(Epetra_SerialDenseMatrix& mat)
  {
    std::memset(mat.A(), 0, sizeof(double) * mat.M() * mat.N());
  }
}  // namespace



template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::FluidEleCalcHDGWeakComp() : usescompletepoly_(true)
{
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::InitializeShapes(
    const DRT::ELEMENTS::Fluid* ele)
{
  // Check if this is an HDG weakly compressible element
  if (const DRT::ELEMENTS::FluidHDGWeakComp* hdgwkele =
          dynamic_cast<const DRT::ELEMENTS::FluidHDGWeakComp*>(ele))
  {
    // use complete polynomial space
    usescompletepoly_ = hdgwkele->UsesCompletePolynomialSpace();

    // initialize shapes
    if (shapes_ == Teuchos::null)
      shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(
          hdgwkele->Degree(), usescompletepoly_, 2 * ele->Degree()));
    else if (shapes_->degree_ != unsigned(ele->Degree()) ||
             shapes_->usescompletepoly_ != usescompletepoly_)
      shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(
          hdgwkele->Degree(), usescompletepoly_, 2 * ele->Degree()));

    // initialize shapes on faces
    if (shapesface_ == Teuchos::null)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(
          ele->Degree(), usescompletepoly_, 2 * ele->Degree());
      shapesface_ = Teuchos::rcp(new DRT::UTILS::ShapeValuesFace<distype>(svfparams));
    }

    // initialize local solver
    if (localSolver_ == Teuchos::null)
      localSolver_ = Teuchos::rcp(new LocalSolver(ele, *shapes_, *shapesface_, usescompletepoly_));
  }
  else
    dserror("Only works for HDG weakly compressible fluid elements");
}



template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
    const DRT::UTILS::GaussIntegration&, bool offdiag)
{
  return this->Evaluate(ele, discretization, lm, params, mat, elemat1_epetra, elemat2_epetra,
      elevec1_epetra, elevec2_epetra, elevec3_epetra, offdiag);
}



template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix&,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector&, Epetra_SerialDenseVector&,
    bool offdiag)
{
  // read global vectors
  ReadGlobalVectors(*ele, discretization, lm);

  // initialize shapes
  InitializeShapes(ele);
  shapes_->Evaluate(*ele, aleDis_);

  // initialize all
  localSolver_->InitializeAll();

  // compute interior residual and matrices
  localSolver_->ComputeInteriorResidual(mat, interiorVal_, interiorAcc_, aleVel_);
  localSolver_->ComputeInteriorMatrices(mat);

  // compute face residual and face matrices
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    shapesface_->EvaluateFace(*ele, f, aleDis_);
    localSolver_->ComputeFaceResidual(f, mat, interiorVal_, traceVal_, aleVel_);
    localSolver_->ComputeFaceMatrices(f, mat);
  }

  // group residuals
  localSolver_->ComputeLocalResidual();
  localSolver_->ComputeGlobalResidual(*ele);

  // group matrices
  localSolver_->ComputeLocalLocalMatrix();
  localSolver_->ComputeLocalGlobalMatrix(*ele);
  localSolver_->ComputeGlobalLocalMatrix(*ele);
  localSolver_->ComputeGlobalGlobalMatrix(*ele);

  // condense local part
  localSolver_->InvertLocalLocalMatrix();
  localSolver_->CondenseLocalResidual(elevec1);
  localSolver_->CondenseLocalMatrix(elemat1);

  // divide rhs by alpha_f
  elevec1.Scale(1.0 / localSolver_->fldparatimint_->AlphaF());

  return 0;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::ReadGlobalVectors(
    const DRT::Element& ele, DRT::Discretization& discretization, const std::vector<int>& lm)
{
  // initialize the vectors
  traceVal_.resize(nfaces_ * (1 + nsd_) * shapesface_->nfdofs_);
  interiorVal_.resize((msd_ + 1 + nsd_) * shapes_->ndofs_);
  interiorAcc_.resize((msd_ + 1 + nsd_) * shapes_->ndofs_);

  // read the trace values
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(0, "velaf");
  DRT::UTILS::ExtractMyValues(*matrix_state, traceVal_, lm);

  // get local dofs
  std::vector<int> localDofs = discretization.Dof(1, &ele);

  // read the interior values
  matrix_state = discretization.GetState(1, "intvelaf");
  DRT::UTILS::ExtractMyValues(*matrix_state, interiorVal_, localDofs);

  // read the interior time derivatives
  matrix_state = discretization.GetState(1, "intaccam");
  DRT::UTILS::ExtractMyValues(*matrix_state, interiorAcc_, localDofs);

  // read ale vectors
  ReadAleVectors(ele, discretization);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::ReadAleVectors(
    const DRT::Element& ele, DRT::Discretization& discretization)
{
  // initialize ale vectors
  aleDis_.resize(nsd_ * nen_);
  aleVel_.resize(nsd_ * nen_);

  // check presence of ale state
  if (discretization.HasState(2, "dispnp"))
  {
    const DRT::ELEMENTS::FluidHDGWeakComp* hdgwkele =
        dynamic_cast<const DRT::ELEMENTS::FluidHDGWeakComp*>(&ele);
    if (hdgwkele->IsAle())
    {
      // get ale dofs
      const int* nodeids = ele.NodeIds();
      std::vector<int> aleDofs;
      for (unsigned int n = 0; n < nen_; ++n)
      {
        std::vector<int> nodeDofs = discretization.Dof(2, discretization.gNode(nodeids[n]));
        for (unsigned int d = 0; d < nsd_; ++d) aleDofs.push_back(nodeDofs[d]);
      }

      // initialize matrix state
      Teuchos::RCP<const Epetra_Vector> matrix_state;

      // read the ale displacement
      matrix_state = discretization.GetState(2, "dispnp");
      DRT::UTILS::ExtractMyValues(*matrix_state, aleDis_, aleDofs);

      // read the ale velocity
      matrix_state = discretization.GetState(2, "gridv");
      DRT::UTILS::ExtractMyValues(*matrix_state, aleVel_, aleDofs);
    }
  }
}



template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::EvaluateService(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params, "action");

  switch (act)
  {
    case FLD::project_fluid_field:
    {
      return ProjectField(ele, params, mat, discretization, lm, elevec1, elevec2);
      break;
    }
    case FLD::update_local_solution:
    {
      return UpdateLocalSolution(ele, params, mat, discretization, lm, elevec1);
      break;
    }
    case FLD::interpolate_hdg_to_node:
    {
      return InterpolateSolutionToNodes(ele, discretization, elevec1);
      break;
    }
    case FLD::calc_fluid_error:
    {
      return ComputeError(ele, params, mat, discretization, lm, elevec1);
      break;
    }
    default:
      dserror("Unknown type of action for FluidHDG");
      break;
  }

  return 0;
}



template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::UpdateLocalSolution(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseVector& interiorinc)
{
  // read global vectors
  ReadGlobalVectors(*ele, discretization, lm);

  // initialize shapes
  InitializeShapes(ele);
  shapes_->Evaluate(*ele, aleDis_);

  // initialize all
  localSolver_->InitializeAll();

  // compute interior residual and matrices
  localSolver_->ComputeInteriorResidual(mat, interiorVal_, interiorAcc_, aleVel_);
  localSolver_->ComputeInteriorMatrices(mat);

  // compute face residual and face matrices
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    shapesface_->EvaluateFace(*ele, f, aleDis_);
    localSolver_->ComputeFaceResidual(f, mat, interiorVal_, traceVal_, aleVel_);
    localSolver_->ComputeFaceMatrices(f, mat);
  }

  // group residuals
  localSolver_->ComputeLocalResidual();

  // group matrices
  localSolver_->ComputeLocalLocalMatrix();
  localSolver_->ComputeLocalGlobalMatrix(*ele);

  // invert local-local matrix
  localSolver_->InvertLocalLocalMatrix();

  // extract local trace increments
  std::vector<double> localtraceinc_vec;
  localtraceinc_vec.resize(nfaces_ * (1 + nsd_) * shapesface_->nfdofs_);
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(0, "globaltraceinc");
  DRT::UTILS::ExtractMyValues(*matrix_state, localtraceinc_vec, lm);

  // convert local trace increments to Epetra_SerialDenseVector
  Epetra_SerialDenseVector localtraceinc(nfaces_ * (1 + nsd_) * shapesface_->nfdofs_);
  for (unsigned int i = 0; i < nfaces_ * (1 + nsd_) * shapesface_->nfdofs_; ++i)
    localtraceinc(i) = localtraceinc_vec[i];

  // compute local solver vector
  Epetra_SerialDenseVector LocalSolverVec((msd_ + 1 + nsd_) * shapes_->ndofs_);
  LocalSolverVec.Multiply('N', 'N', 1.0, localSolver_->KlocallocalInv, localSolver_->Rlocal, 0.0);

  // compute local solver matrix
  Epetra_SerialDenseMatrix LocalSolverMat(
      (msd_ + 1 + nsd_) * shapes_->ndofs_, nfaces_ * (1 + nsd_) * shapesface_->nfdofs_);
  LocalSolverMat.Multiply(
      'N', 'N', -1.0, localSolver_->KlocallocalInv, localSolver_->Klocalglobal, 0.0);

  // compute local increments
  interiorinc.Shape((msd_ + 1 + nsd_) * shapes_->ndofs_, 1);
  interiorinc = LocalSolverVec;
  interiorinc.Multiply('N', 'N', 1.0, LocalSolverMat, localtraceinc, 1.0);

  return 0;
}



template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::ComputeError(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseVector& elevec)
{
  // read ale vectors
  ReadAleVectors(*ele, discretization);

  // initialize shapes
  InitializeShapes(ele);
  shapes_->Evaluate(*ele, aleDis_);

  // get time
  const double time = localSolver_->fldparatimint_->Time();

  // get interior values
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "intvelnp");
  std::vector<int> localDofs = discretization.Dof(1, ele);
  std::vector<double> vecValues(localDofs.size());
  for (unsigned int i = 0; i < localDofs.size(); ++i)
  {
    const int lid = matrix_state->Map().LID(localDofs[i]);
    vecValues[i] = (*matrix_state)[lid];
  }

  // initialize exact solution
  LINALG::Matrix<msd_, 1> L_ex(true);
  double r_ex = 0.0;
  LINALG::Matrix<nsd_, 1> w_ex(true);

  // initialize spatial coordinates
  LINALG::Matrix<nsd_, 1> xyz(true);

  // get function number
  const int calcerrfunctno =
      DRT::INPUT::get<INPAR::FLUID::CalcError>(params, "error function number");

  // initialize errors
  double err_L = 0.0;
  double err_r = 0.0;
  double err_w = 0.0;
  double norm_L = 0.0;
  double norm_r = 0.0;
  double norm_w = 0.0;

  // initialize numerical values
  Epetra_SerialDenseVector Leg;
  Epetra_SerialDenseVector reg;
  Epetra_SerialDenseVector weg;

  // ease notation
  LINALG::SerialDenseMatrix N = shapes_->shfunct;
  LINALG::SerialDenseVector fac = shapes_->jfac;
  unsigned int nqpoints = shapes_->nqpoints_;
  unsigned int ndofs = shapes_->ndofs_;

  // loop over quadrature points
  for (unsigned int q = 0; q < nqpoints; ++q)
  {
    // clear vectors
    Leg.Shape(msd_, 1);
    reg.Shape(1, 1);
    weg.Shape(nsd_, 1);

    // interpolate values on gauss points
    for (unsigned int i = 0; i < ndofs; ++i)
    {
      for (unsigned int m = 0; m < msd_; ++m) Leg(m) += N(i, q) * vecValues[m * ndofs + i];

      reg(0) += N(i, q) * vecValues[msd_ * ndofs + i];

      for (unsigned int d = 0; d <= nsd_; ++d)
        weg(d) += N(i, q) * vecValues[(msd_ + 1 + d) * ndofs + i];
    }

    for (unsigned int d = 0; d < nsd_; ++d) xyz(d) = shapes_->xyzreal(d, q);

    // evaluate analytical solution
    EvaluateAll(calcerrfunctno, xyz, time, L_ex, r_ex, w_ex);

    // evaluate errors
    for (unsigned int m = 0; m < msd_; ++m)
      err_L += (L_ex(m) - Leg(m)) * (L_ex(m) - Leg(m)) * fac(q);

    err_r += (r_ex - reg(0)) * (r_ex - reg(0)) * fac(q);

    for (unsigned int d = 0; d < nsd_; ++d)
      err_w += (w_ex(d) - weg(d)) * (w_ex(d) - weg(d)) * fac(q);

    // evaluate norms
    for (unsigned int m = 0; m < msd_; ++m) norm_L += L_ex(m) * L_ex(m) * fac(q);

    norm_r += r_ex * r_ex * fac(q);

    for (unsigned int d = 0; d < nsd_; ++d) norm_w += w_ex(d) * w_ex(d) * fac(q);
  }

  // fill vector
  elevec[0] += err_L;
  elevec[1] += err_r;
  elevec[2] += err_w;
  elevec[3] += norm_L;
  elevec[4] += norm_r;
  elevec[5] += norm_w;

  return 0;
};



template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::ProjectField(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2)
{
  // read ale vectors
  ReadAleVectors(*ele, discretization);

  // Create the necessary objects to the solution of the problem as the solver
  // and the shape functions for both the interior, shapes_, and the trace, shapesface_.
  InitializeShapes(ele);

  // Evaluate the element at the gauss points
  shapes_->Evaluate(*ele, aleDis_);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 || elevec2.M() == static_cast<int>((msd_ + 1 + nsd_) * shapes_->ndofs_),
      "Wrong size in project vector 2");

  // get initial function and current time
  const int* initfield = params.getPtr<int>("initfield");
  const int* startfunc = params.getPtr<int>("startfuncno");
  double* time = params.getPtr<double>("time");

  if (elevec2.M() > 0)
  {
    // Create the local matrix from starting at the addres where elevec2 is with the right shape
    Epetra_SerialDenseMatrix localMat(
        View, elevec2.A(), shapes_->ndofs_, shapes_->ndofs_, msd_ + 1 + nsd_, false);
    // Initialize matrix to zeros
    zeroMatrix(localMat);

    // create mass matrix for interior by looping over quadrature points
    for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
    {
      // jfac is a vector containing the jacobian times the weight of the quadrature points
      const double fac = shapes_->jfac(q);
      // xyz is a vector containing the coordiantes of the quadrature points in real coordinates
      LINALG::Matrix<nsd_, 1> xyz(false);
      // Filling xyz with the values take from the element xyzreal matrix
      for (unsigned int d = 0; d < nsd_; ++d) xyz(d) = shapes_->xyzreal(d, q);
      // Declaring vectors for interior variables
      LINALG::Matrix<msd_, 1> L(true);
      double r = 0.0;
      LINALG::Matrix<nsd_, 1> w(true);

      dsassert(initfield != NULL && startfunc != NULL,
          "initfield or startfuncno not set for initial value");

      // This function returns th value of the interior variables from the
      // given initial field that can be a know field or a user-defined one
      EvaluateDensityMomentum(*startfunc, xyz, 0.0, r, w);

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      {
        localSolver_->massPart(i, q) = shapes_->shfunct(i, q);
        localSolver_->massPartW(i, q) = shapes_->shfunct(i, q) * fac;

        // RHS for the mixed variable
        for (unsigned int m = 0; m < msd_; ++m)
          localMat(i, m) += shapes_->shfunct(i, q) * L(m) * fac;

        // RHS for the density
        localMat(i, msd_) += shapes_->shfunct(i, q) * r * fac;

        // RHS for the momenutm
        for (unsigned int d = 0; d < nsd_; ++d)
          localMat(i, msd_ + 1 + d) += shapes_->shfunct(i, q) * w(d) * fac;
      }
    }
    // The integration is made by computing the matrix product
    localSolver_->massMat.Multiply(
        'N', 'T', 1., localSolver_->massPart, localSolver_->massPartW, 0.);
    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(localSolver_->massMat);
    inverseMass.SetVectors(localMat, localMat);
    inverseMass.Solve();
  }

  // Here we have the projection of the field on the trace
  Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
  // trVec is the vector of the trace values
  Epetra_SerialDenseMatrix trVec(shapesface_->nfdofs_, 1 + nsd_);
  dsassert(elevec1.M() == static_cast<int>((1 + nsd_) * shapesface_->nfdofs_) ||
               elevec1.M() == static_cast<int>(nfaces_ * (1 + nsd_) * shapesface_->nfdofs_),
      "Wrong size in project vector 1");

  const unsigned int* faceConsider = params.getPtr<unsigned int>("faceconsider");
  Teuchos::Array<int>* functno = params.getPtr<Teuchos::Array<int>>("funct");
  Teuchos::Array<int>* onoff = params.getPtr<Teuchos::Array<int>>("onoff");

  // Project the field for all the faces of the element
  for (unsigned int face = 0; face < nfaces_; ++face)
  {
    // check whether we are in the project phase for all faces or for boundary values
    if (initfield == NULL)
    {
      // We get here only if it is not an initial value but it is a time
      // dependant boundary value. If we are here we only want the function to run
      // for boundary faces specified in the faceConsider variable
      dsassert(faceConsider != NULL, "Unsupported operation");
      if (*faceConsider != face) continue;
    }

    // the same function as before but for the trace elements
    shapesface_->EvaluateFace(*ele, face, aleDis_);

    // Initializing the matrices
    zeroMatrix(mass);
    zeroMatrix(trVec);

    // For each quadrature point we evaluate the trace values and the shape functions
    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      // shapesface_->jfac contains the jacobian evaluated in the quadrature points
      const double fac = shapesface_->jfac(q);
      // xyz is the vector containing the coordinates of the quadrature points
      LINALG::Matrix<nsd_, 1> xyz(false);

      // Taking the real coordinates of quadrature points of the current face
      for (unsigned int d = 0; d < nsd_; ++d) xyz(d) = shapesface_->xyzreal(d, q);

      // Creating the vector of traces variables
      double r;
      LINALG::Matrix<nsd_, 1> w(false);

      // Create dummy variables
      double dummy_r;
      LINALG::Matrix<nsd_, 1> dummy_w;

      // Deciding if we are initializing a field or if it is a time dependant
      // boundary condition
      if (initfield != NULL)  // Initial function
        EvaluateDensityMomentum(*startfunc, xyz, 0.0, r, w);
      else
      {
        // This is used to project a function only on the boundary during the
        // temporal evolution of the simulation
        dsassert(functno != NULL && time != NULL && onoff != NULL, "No array with functions given");

        // Deciding if to use the function or not for the density
        if ((*onoff)[0] == 1)
        {
          const int funct_num = (*functno)[0];
          if (funct_num > 0) EvaluateDensityMomentum(funct_num, xyz, *time, r, dummy_w);
        }

        for (unsigned int d = 0; d < nsd_; ++d)
        {
          // Deciding if to use the function or not for the current component of the momentum
          if ((*onoff)[1 + d] == 0) continue;

          // If we are using the function, evaluate it
          const int funct_num = (*functno)[1 + d];
          if (funct_num > 0) EvaluateDensityMomentum(funct_num, xyz, *time, dummy_r, w);
        }
      }
      // now fill the components in the mass matrix and the right hand side

      // This is a more usual way to compute the mass matrix (double for loop)
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
          mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;

        // RHS
        trVec(i, 0) += shapesface_->shfunct(i, q) * r * fac;
        for (unsigned int d = 0; d < nsd_; ++d)
          trVec(i, 1 + d) += shapesface_->shfunct(i, q) * w(d) * fac;
      }
    }

    // Solving step, nothing fancy
    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(mass);
    // In this cas trVec is a proper vector and not a matrix used as multiple
    // RHS vectors
    inverseMass.SetVectors(trVec, trVec);
    inverseMass.Solve();

    // In this case we fill elevec1 with the values of trVec because we have not
    // defined trVec as a matrix beginning where elevec1 begins
    if (initfield != NULL)  // This is for initial functions
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        elevec1(face * shapesface_->nfdofs_ * (1 + nsd_) + i) = trVec(i, 0);
        for (unsigned int d = 0; d < nsd_; ++d)
          elevec1(face * shapesface_->nfdofs_ * (1 + nsd_) + (1 + d) * shapesface_->nfdofs_ + i) =
              trVec(i, 1 + d);
      }
    else  // This is only for boundary faces during time evolution
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        elevec1(i) = trVec(i, 0);
        for (unsigned int d = 0; d < nsd_; ++d)
          elevec1((1 + d) * shapesface_->nfdofs_ + i) = trVec(i, 1 + d);
      }
  }  // for over the faces

  return 0;
}



template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::InterpolateSolutionToNodes(
    DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
    Epetra_SerialDenseVector& elevec1)
{
  // read ale vectors
  ReadAleVectors(*ele, discretization);

  InitializeShapes(ele);
  // Check if the vector has the correct size
  dsassert(
      elevec1.M() == (int)nen_ * (msd_ + 1 + nsd_ + 1 + nsd_), "Vector does not have correct size");

  // Getting the connectivity matrix
  // Contains the (local) coordinates of the nodes belonging to the element
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  // This vector will contain the values of the shape functions computed in a
  // certain coordinate. In fact the lenght of the vector is given by the number
  // of shape functions, that is the same of the number of degrees of freedom of
  // an element.
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  // get local solution values
  // The vector "matrix_state" contains the interior velocity values following
  // the local id numbers
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "intvelnp");
  // Vector of the ids of the DOF for the element
  std::vector<int> localDofs = discretization.Dof(1, ele);
  // SOLution VALUES
  std::vector<double> solvalues(localDofs.size());

  // Filling every entry of the solvalue vector obtaining the values from the
  //"matrix_state" vector.
  for (unsigned int i = 0; i < solvalues.size(); ++i)
  {
    // Finding the local id of the current "localDofs"
    const int lid = matrix_state->Map().LID(localDofs[i]);
    // Saving the value of the "localDofs[i]" in the "solvalues" vector
    solvalues[i] = (*matrix_state)[lid];
  }
  elevec1.Scale(0.);

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

    // compute values for internal unknowns by summing over all basis functions
    for (unsigned int m = 0; m <= msd_; ++m)
    {
      double sum_mixedvar = 0.0;
      for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
        sum_mixedvar += values(k) * solvalues[m * shapes_->ndofs_ + k];
      elevec1(m * nen_ + i) = sum_mixedvar;
    }
    {
      double sum_density = 0.0;
      for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
        sum_density += values(k) * solvalues[msd_ * shapes_->ndofs_ + k];
      elevec1(msd_ * nen_ + i) = sum_density;
    }
    for (unsigned int d = 0; d <= nsd_; ++d)
    {
      double sum_momentum = 0.0;
      for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
        sum_momentum += values(k) * solvalues[(msd_ + 1 + d) * shapes_->ndofs_ + k];
      elevec1((msd_ + 1 + d) * nen_ + i) = sum_momentum;
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

  // Now the vector "matrix_state" contains the trace velocity values following
  // the local id numbers
  matrix_state = discretization.GetState(0, "velnp");

  // we have always two dofsets
  Element::LocationArray la(2);
  ele->LocationVector(discretization, la, false);
  localDofs = la[0].lm_;
  solvalues.resize(localDofs.size());

  for (unsigned int i = 0; i < solvalues.size(); ++i)
  {
    const int lid = matrix_state->Map().LID(localDofs[i]);
    solvalues[i] = (*matrix_state)[lid];
  }
  for (int i = (msd_ + 1 + nsd_) * nen_; i < elevec1.M(); ++i) elevec1(i) = 0.0;

  Epetra_SerialDenseVector fvalues(shapesface_->nfdofs_);
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    // Checking how many nodes the face has
    const int nfn = DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace;

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

      // compute values for the trace by summing over all basis functions
      {
        double sum_trace_density = 0.0;
        for (unsigned int k = 0; k < shapesface_->nfdofs_; ++k)
          sum_trace_density += fvalues(k) * solvalues[f * shapesface_->nfdofs_ * (1 + nsd_) + k];
        if (i < ncn[f])
          elevec1((msd_ + 1 + nsd_) * nen_ + shapesface_->faceNodeOrder[f][i]) +=
              sum_trace_density / nsd_;
        else if (i < nfn - nin[f])
          elevec1((msd_ + 1 + nsd_) * nen_ + shapesface_->faceNodeOrder[f][i]) +=
              sum_trace_density / (nsd_ - 1);
        else
          elevec1((msd_ + 1 + nsd_) * nen_ + shapesface_->faceNodeOrder[f][i]) += sum_trace_density;
      }
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        double sum_trace_momentum = 0.0;
        for (unsigned int k = 0; k < shapesface_->nfdofs_; ++k)
        {
          sum_trace_momentum +=
              fvalues(k) *
              solvalues[f * shapesface_->nfdofs_ * (1 + nsd_) + (1 + d) * shapesface_->nfdofs_ + k];
        }
        if (i < ncn[f])
          elevec1((msd_ + 1 + nsd_ + 1 + d) * nen_ + shapesface_->faceNodeOrder[f][i]) +=
              sum_trace_momentum / nsd_;
        else if (i < nfn - nin[f])
          elevec1((msd_ + 1 + nsd_ + 1 + d) * nen_ + shapesface_->faceNodeOrder[f][i]) +=
              sum_trace_momentum / (nsd_ - 1);
        else
          elevec1((msd_ + 1 + nsd_ + 1 + d) * nen_ + shapesface_->faceNodeOrder[f][i]) +=
              sum_trace_momentum;
      }
    }
  }

  return 0;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::EvaluateAll(const int funcnum,
    const LINALG::Matrix<nsd_, 1>& xyz, const double t, LINALG::Matrix<msd_, 1>& L, double& r,
    LINALG::Matrix<nsd_, 1>& w) const
{
  r = DRT::Problem::Instance()->Funct(funcnum - 1).Evaluate(0, xyz.A(), t);

  for (unsigned int d = 0; d < nsd_; ++d)
    w(d) = DRT::Problem::Instance()->Funct(funcnum - 1).Evaluate(1 + d, xyz.A(), t);

  for (unsigned int m = 0; m < msd_; ++m)
    L(m) = DRT::Problem::Instance()->Funct(funcnum - 1).Evaluate(1 + nsd_ + m, xyz.A(), t);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::EvaluateDensityMomentum(const int funcnum,
    const LINALG::Matrix<nsd_, 1>& xyz, const double t, double& r, LINALG::Matrix<nsd_, 1>& w) const
{
  r = DRT::Problem::Instance()->Funct(funcnum - 1).Evaluate(0, xyz.A(), t);

  for (unsigned int d = 0; d < nsd_; ++d)
    w(d) = DRT::Problem::Instance()->Funct(funcnum - 1).Evaluate(1 + d, xyz.A(), t);
}



template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>*
DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::Instance(bool create)
{
  static FluidEleCalcHDGWeakComp<distype>* instance;
  if (create)
  {
    if (instance == NULL)
    {
      instance = new FluidEleCalcHDGWeakComp<distype>();
    }
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}



template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::LocalSolver(
    const DRT::ELEMENTS::Fluid* ele, const DRT::UTILS::ShapeValues<distype>& shapeValues,
    DRT::UTILS::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly)
    : ndofs_(shapeValues.ndofs_),
      convective(true),
      unsteady(true),
      ale(ele->IsAle()),
      shapes_(shapeValues),
      shapesface_(shapeValuesFace)
{
  // count ndofs on faces
  ndofsfaces_ = 0;
  for (unsigned int i = 0; i < nfaces_; ++i)
  {
    shapesface_.EvaluateFace(*ele, i);
    ndofsfaces_ += shapesface_.nfdofs_;
  }

  // initialize stabilization parameters
  tau_r = 0.0;
  tau_w = 0.0;

  // initialize auxiliary matrices
  massPart.Shape(ndofs_, shapes_.nqpoints_);
  massPartW.Shape(ndofs_, shapes_.nqpoints_);
  massMat.Shape(ndofs_, ndofs_);

  // initialize unknowns
  Leg.Shape(msd_, shapes_.nqpoints_);
  reg.Shape(shapes_.nqpoints_, 1);
  weg.Shape(nsd_, shapes_.nqpoints_);

  // initialize ALE variables
  aeg.Shape(nsd_, shapes_.nqpoints_);
  dadxyzeg.Shape(nsd_ * nsd_, shapes_.nqpoints_);

  // initialize matrices
  ALL.Shape(msd_ * ndofs_, msd_ * ndofs_);
  ALr.Shape(msd_ * ndofs_, 1 * ndofs_);
  ALw.Shape(msd_ * ndofs_, nsd_ * ndofs_);
  ALR.Shape(msd_ * ndofs_, 1 * ndofsfaces_);
  ALW.Shape(msd_ * ndofs_, nsd_ * ndofsfaces_);
  Arr.Shape(1 * ndofs_, 1 * ndofs_);
  Arw.Shape(1 * ndofs_, nsd_ * ndofs_);
  ArR.Shape(1 * ndofs_, 1 * ndofsfaces_);
  ArW.Shape(1 * ndofs_, nsd_ * ndofsfaces_);
  AwL.Shape(nsd_ * ndofs_, msd_ * ndofs_);
  Awr.Shape(nsd_ * ndofs_, 1 * ndofs_);
  Aww.Shape(nsd_ * ndofs_, nsd_ * ndofs_);
  AwR.Shape(nsd_ * ndofs_, 1 * ndofsfaces_);
  AwW.Shape(nsd_ * ndofs_, nsd_ * ndofsfaces_);
  ARr.Shape(1 * ndofsfaces_, 1 * ndofs_);
  ARR.Shape(1 * ndofsfaces_, 1 * ndofsfaces_);
  AWL.Shape(nsd_ * ndofsfaces_, msd_ * ndofs_);
  AWw.Shape(nsd_ * ndofsfaces_, nsd_ * ndofs_);
  AWR.Shape(nsd_ * ndofsfaces_, 1 * ndofsfaces_);
  AWW.Shape(nsd_ * ndofsfaces_, nsd_ * ndofsfaces_);

  // initialize residuals
  RL.Shape(msd_ * ndofs_, 1);
  Rr.Shape(1 * ndofs_, 1);
  Rw.Shape(nsd_ * ndofs_, 1);
  RR.Shape(1 * ndofsfaces_, 1);
  RW.Shape(nsd_ * ndofsfaces_, 1);

  // initialize local/global matrices/vectors
  Klocallocal.Shape((msd_ + 1 + nsd_) * ndofs_, (msd_ + 1 + nsd_) * ndofs_);
  Klocalglobal.Shape((msd_ + 1 + nsd_) * ndofs_, (1 + nsd_) * ndofsfaces_);
  Kgloballocal.Shape((1 + nsd_) * ndofsfaces_, (msd_ + 1 + nsd_) * ndofs_);
  Kglobalglobal.Shape((1 + nsd_) * ndofsfaces_, (1 + nsd_) * ndofsfaces_);
  Rlocal.Shape((msd_ + 1 + nsd_) * ndofs_, 1);
  Rglobal.Shape((1 + nsd_) * ndofsfaces_, 1);
  KlocallocalInv.Shape((msd_ + 1 + nsd_) * ndofs_, (msd_ + 1 + nsd_) * ndofs_);

  // pair of indices in Voigt notation
  int s = 0;
  for (unsigned int d = 0; d < nsd_; ++d)
    for (unsigned int e = d + 1; e < nsd_; ++e)
    {
      VoigtP[s][0] = d;
      VoigtP[s][1] = e;
      s = s + 1;
    }

  // pointer to class FluidEleParameter (access to the general parameter)
  fldparatimint_ = Teuchos::rcp(DRT::ELEMENTS::FluidEleParameterTimInt::Instance(), false);

  // initialize also general parameter list, also it will be overwritten in derived subclasses
  fldpara_ = Teuchos::rcp(DRT::ELEMENTS::FluidEleParameterStd::Instance(), false);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::InitializeAll()
{
  // initialize unknowns
  zeroMatrix(Leg);
  zeroMatrix(reg);
  zeroMatrix(weg);

  // initialize ALE variables
  zeroMatrix(aeg);
  zeroMatrix(dadxyzeg);

  // initialize matrices
  zeroMatrix(ALL);
  zeroMatrix(ALr);
  zeroMatrix(ALw);
  zeroMatrix(ALR);
  zeroMatrix(ALW);
  zeroMatrix(Arr);
  zeroMatrix(Arw);
  zeroMatrix(ArR);
  zeroMatrix(ArW);
  zeroMatrix(AwL);
  zeroMatrix(Awr);
  zeroMatrix(Aww);
  zeroMatrix(AwR);
  zeroMatrix(AwW);
  zeroMatrix(ARr);
  zeroMatrix(ARR);
  zeroMatrix(AWL);
  zeroMatrix(AWw);
  zeroMatrix(AWR);
  zeroMatrix(AWW);

  // initialize residuals
  zeroMatrix(RL);
  zeroMatrix(Rr);
  zeroMatrix(Rw);
  zeroMatrix(RR);
  zeroMatrix(RW);

  // initialize local/global matrices/vectors
  zeroMatrix(Klocallocal);
  zeroMatrix(Klocalglobal);
  zeroMatrix(Kgloballocal);
  zeroMatrix(Kglobalglobal);
  zeroMatrix(Rlocal);
  zeroMatrix(Rglobal);
  zeroMatrix(KlocallocalInv);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeMaterialMatrix(
    const Teuchos::RCP<MAT::Material>& mat, const LINALG::Matrix<nsd_, 1>& xyz,
    Epetra_SerialDenseMatrix& DL, Epetra_SerialDenseMatrix& Dw)
{
  // initialize DL and Dw
  DL.Shape(msd_, msd_);
  Dw.Shape(msd_, msd_);

  // evaluate D_fac
  Epetra_SerialDenseMatrix D_fac(msd_, msd_);
  for (unsigned int m = 0; m < msd_; ++m) D_fac(m, m) = 1.0;
  for (unsigned int d = 0; d < nsd_; ++d)
    for (unsigned int e = 0; e < nsd_; ++e) D_fac(d, e) -= 1.0 / 3.0;

  // check for variable viscosity
  const Teuchos::ParameterList& fluidparams = DRT::Problem::Instance()->FluidDynamicParams();
  int varviscfuncnum = fluidparams.get<int>("VARVISCFUNCNO");

  // compute material matrix
  if (varviscfuncnum < 0)  // constant viscosity
  {
    // get material
    const MAT::WeaklyCompressibleFluid* actmat =
        static_cast<const MAT::WeaklyCompressibleFluid*>(mat.get());

    // get viscosity
    double mu = actmat->Viscosity();

    // evaluate Dw
    for (unsigned int d = 0; d < nsd_; ++d) Dw(d, d) = std::pow(2.0 * mu, 1.0 / 2.0);
    for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
      Dw(nsd_ + s, nsd_ + s) = std::pow(mu, 1.0 / 2.0);

    // evaluate DL
    DL.Multiply('N', 'N', 1., Dw, D_fac, 0.);
  }
  else  // variable viscosity
  {
    // get time
    const double time = fldparatimint_->Time();

    // get viscosity
    double mu = DRT::Problem::Instance()->Funct(varviscfuncnum - 1).Evaluate(0, xyz.A(), time);

    // evaluate Dw
    for (unsigned int d = 0; d < nsd_; ++d) Dw(d, d) = 2.0 * mu;
    for (unsigned int s = 0; s < (msd_ - nsd_); ++s) Dw(nsd_ + s, nsd_ + s) = mu;

    // evaluate DL
    DL = D_fac;
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeInteriorResidual(
    const Teuchos::RCP<MAT::Material>& mat, const std::vector<double>& val,
    const std::vector<double>& accel, const std::vector<double>& alevel)
{
  // set convective flag
  INPAR::FLUID::PhysicalType physicaltype = fldpara_->PhysicalType();
  convective = (physicaltype != INPAR::FLUID::weakly_compressible_stokes_dens_mom);

  // set unsteady flag
  const Teuchos::ParameterList& fluidparams = DRT::Problem::Instance()->FluidDynamicParams();
  std::string timeintegr = fluidparams.get<std::string>("TIMEINTEGR");
  unsteady = (timeintegr.compare("Stationary") != 0);

  // get material properties
  const MAT::WeaklyCompressibleFluid* actmat =
      static_cast<const MAT::WeaklyCompressibleFluid*>(mat.get());

  // get force function
  int forcefuncnum = fluidparams.get<int>("BODYFORCEFUNCNO");

  // get time
  const double time = fldparatimint_->Time();

  // ease notation
  LINALG::SerialDenseMatrix N = shapes_.shfunct;
  LINALG::SerialDenseMatrix Nn = shapes_.funct;
  LINALG::SerialDenseMatrix Nxyz = shapes_.shderxy;
  LINALG::SerialDenseMatrix Nnxyz = shapes_.derxy;
  LINALG::SerialDenseVector fac = shapes_.jfac;
  unsigned int nqpoints = shapes_.nqpoints_;

  // extract interior nodal values
  LINALG::Matrix<nsd_, 1> xyze;
  Epetra_SerialDenseMatrix Le(msd_, ndofs_);
  Epetra_SerialDenseVector re(ndofs_);
  Epetra_SerialDenseMatrix we(nsd_, ndofs_);
  Epetra_SerialDenseVector pe(ndofs_);
  Epetra_SerialDenseVector drdte(ndofs_);
  Epetra_SerialDenseMatrix dwdte(nsd_, ndofs_);
  Epetra_SerialDenseMatrix DLe(msd_, msd_);
  Epetra_SerialDenseMatrix Dwe(msd_, msd_);
  Epetra_SerialDenseMatrix DwLe(msd_, ndofs_);

  for (unsigned int i = 0; i < ndofs_; ++i)
  {
    for (unsigned int d = 0; d < nsd_; ++d) xyze(d) = shapes_.nodexyzreal(d, i);

    for (unsigned int m = 0; m < msd_; ++m) Le(m, i) = val[m * ndofs_ + i];

    re(i) = val[msd_ * ndofs_ + i];

    for (unsigned int d = 0; d < nsd_; ++d) we(d, i) = val[(msd_ + 1 + d) * ndofs_ + i];

    pe(i) = actmat->ComputePressure(re[i]);

    drdte(i) = accel[msd_ * ndofs_ + i];

    for (unsigned int d = 0; d < nsd_; ++d) dwdte(d, i) = accel[(msd_ + 1 + d) * ndofs_ + i];

    // compute material matrix
    ComputeMaterialMatrix(mat, xyze, DLe, Dwe);

    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int e = 0; e < nsd_; ++e) DwLe(d, i) += Dwe(d, e) * Le(e, i);
    for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
      DwLe(nsd_ + s, i) += Dwe(nsd_ + s, nsd_ + s) * Le(nsd_ + s, i);
  }

  // extract ale nodal values
  Epetra_SerialDenseMatrix ae(nsd_, ndofs_);

  if (ale)
    for (unsigned int n = 0; n < nen_; ++n)
      for (unsigned int d = 0; d < nsd_; ++d) ae(d, n) = alevel[n * nsd_ + d];

  // initialize values interpolated on gauss points
  LINALG::Matrix<nsd_, 1> xyzeg;
  Epetra_SerialDenseMatrix feg(1 + nsd_, nqpoints);
  Epetra_SerialDenseVector drdteg(nqpoints);
  Epetra_SerialDenseMatrix dwdteg(nsd_, nqpoints);
  Epetra_SerialDenseMatrix dDwLdxyzeg(msd_ * nsd_, nqpoints);
  Epetra_SerialDenseMatrix dpdxyzeg(1 * nsd_, nqpoints);
  Epetra_SerialDenseMatrix DLeg(msd_, msd_);
  Epetra_SerialDenseMatrix Dweg(msd_, msd_);

  // loop over quadrature points
  for (unsigned int q = 0; q < nqpoints; ++q)
  {
    // interpolate values on gauss points
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int d = 0; d < nsd_; ++d) xyzeg(d) = shapes_.xyzreal(d, q);

      for (unsigned int m = 0; m < msd_; ++m) Leg(m, q) += N(i, q) * Le(m, i);

      reg(q) += N(i, q) * re(i);

      for (unsigned int d = 0; d < nsd_; ++d) weg(d, q) += N(i, q) * we(d, i);

      if (forcefuncnum > 0)
        for (unsigned int dmod = 0; dmod < (1 + nsd_); ++dmod)
          feg(dmod, q) =
              DRT::Problem::Instance()->Funct(forcefuncnum - 1).Evaluate(dmod, xyzeg.A(), time);

      drdteg(q) += N(i, q) * drdte(i);

      for (unsigned int d = 0; d < nsd_; ++d) dwdteg(d, q) += N(i, q) * dwdte(d, i);

      for (unsigned int m = 0; m < msd_; ++m)
        for (unsigned int d = 0; d < nsd_; ++d)
          dDwLdxyzeg(m * nsd_ + d, q) += Nxyz(i * nsd_ + d, q) * DwLe(m, i);

      for (unsigned int d = 0; d < nsd_; ++d) dpdxyzeg(d, q) += Nxyz(i * nsd_ + d, q) * pe(i);
    }

    // compute material matrix
    ComputeMaterialMatrix(mat, xyzeg, DLeg, Dweg);

    // interpolate ale values on gauss points
    if (ale)
      for (unsigned int n = 0; n < nen_; ++n)
      {
        for (unsigned int d = 0; d < nsd_; ++d) aeg(d, q) += Nn(n, q) * ae(d, n);

        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            dadxyzeg(d * nsd_ + e, q) += Nnxyz(n * nsd_ + e, q) * ae(d, n);
      }

    // compute residuals
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int m = 0; m < msd_; ++m) RL(m * ndofs_ + i) += +N(i, q) * Leg(m, q) * fac(q);

      for (unsigned int d = 0; d < nsd_; ++d)
        for (unsigned int e = 0; e < nsd_; ++e)
          RL(d * ndofs_ + i) += -Nxyz(i * nsd_ + e, q) * DLeg(e, d) * weg(e, q) / reg(q) * fac(q);
      for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
        RL((nsd_ + s) * ndofs_ + i) += -Nxyz(i * nsd_ + VoigtP[s][0], q) *
                                           DLeg(nsd_ + s, nsd_ + s) * weg(VoigtP[s][1], q) /
                                           reg(q) * fac(q) -
                                       Nxyz(i * nsd_ + VoigtP[s][1], q) * DLeg(nsd_ + s, nsd_ + s) *
                                           weg(VoigtP[s][0], q) / reg(q) * fac(q);

      if (unsteady) Rr(i) += -N(i, q) * drdteg(q) * fac(q);

      if (ale)
        for (unsigned int d = 0; d < nsd_; ++d)
          Rr(i) += -N(i, q) * reg(q) * dadxyzeg(d * nsd_ + d, q) * fac(q);

      for (unsigned int d = 0; d < nsd_; ++d) Rr(i) += +Nxyz(i * nsd_ + d, q) * weg(d, q) * fac(q);

      if (ale)
        for (unsigned int d = 0; d < nsd_; ++d)
          Rr(i) += -Nxyz(i * nsd_ + d, q) * reg(q) * aeg(d, q) * fac(q);

      Rr(i) += +N(i, q) * feg(0, q) * fac(q);

      if (unsteady)
        for (unsigned int d = 0; d < nsd_; ++d)
          Rw(d * ndofs_ + i) += -N(i, q) * dwdteg(d, q) * fac(q);

      if (ale)
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            Rw(d * ndofs_ + i) += -N(i, q) * weg(d, q) * dadxyzeg(e * nsd_ + e, q) * fac(q);

      if (convective)
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            Rw(d * ndofs_ + i) += +Nxyz(i * nsd_ + e, q) * weg(d, q) * weg(e, q) / reg(q) * fac(q);

      if (ale)
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            Rw(d * ndofs_ + i) += -Nxyz(i * nsd_ + e, q) * weg(d, q) * aeg(e, q) * fac(q);

      for (unsigned int d = 0; d < nsd_; ++d)
        Rw(d * ndofs_ + i) += -N(i, q) * dDwLdxyzeg(d * nsd_ + d, q) * fac(q);
      for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
      {
        Rw(VoigtP[s][0] * ndofs_ + i) +=
            -N(i, q) * dDwLdxyzeg((nsd_ + s) * nsd_ + VoigtP[s][1], q) * fac(q);
        Rw(VoigtP[s][1] * ndofs_ + i) +=
            -N(i, q) * dDwLdxyzeg((nsd_ + s) * nsd_ + VoigtP[s][0], q) * fac(q);
      }

      for (unsigned int d = 0; d < nsd_; ++d)
        Rw(d * ndofs_ + i) += -N(i, q) * dpdxyzeg(d, q) * fac(q);

      for (unsigned int d = 0; d < nsd_; ++d)
        Rw(d * ndofs_ + i) += +N(i, q) * feg(1 + d, q) * fac(q);
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeInteriorMatrices(
    const Teuchos::RCP<MAT::Material>& mat)
{
  // get material properties
  const MAT::WeaklyCompressibleFluid* actmat =
      static_cast<const MAT::WeaklyCompressibleFluid*>(mat.get());
  double eps = actmat->ComprCoeff();

  // inverse of time factor
  const double invtimefac = 1.0 / (fldparatimint_->TimeFac());

  // ease notation
  LINALG::SerialDenseMatrix N = shapes_.shfunct;
  LINALG::SerialDenseMatrix Nxyz = shapes_.shderxy;
  LINALG::SerialDenseVector fac = shapes_.jfac;
  unsigned int nqpoints = shapes_.nqpoints_;

  // extract interior nodal values
  LINALG::Matrix<nsd_, 1> xyze;
  Epetra_SerialDenseMatrix DLe(msd_, msd_);
  Epetra_SerialDenseMatrix Dwe(msd_, msd_);
  Epetra_SerialDenseMatrix Dwemod(nsd_ * nsd_ + (msd_ - nsd_), ndofs_);

  for (unsigned int i = 0; i < ndofs_; ++i)
  {
    for (unsigned int d = 0; d < nsd_; ++d) xyze(d) = shapes_.nodexyzreal(d, i);

    // compute material matrix
    ComputeMaterialMatrix(mat, xyze, DLe, Dwe);

    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int e = 0; e < nsd_; ++e) Dwemod(d * nsd_ + e, i) = Dwe(d, e);
    for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
      Dwemod(nsd_ * nsd_ + s, i) = Dwe(nsd_ + s, nsd_ + s);
  }

  // initialize values interpolated on gauss points
  LINALG::Matrix<nsd_, 1> xyzeg;
  Epetra_SerialDenseMatrix DLeg(msd_, msd_);
  Epetra_SerialDenseMatrix Dweg(msd_, msd_);
  Epetra_SerialDenseMatrix dDwdxyzeg((nsd_ * nsd_ + (msd_ - nsd_)) * nsd_, nqpoints);

  // loop over quadrature points
  for (unsigned int q = 0; q < nqpoints; ++q)
  {
    // interpolate values on gauss points
    for (unsigned int d = 0; d < nsd_; ++d) xyzeg(d) = shapes_.xyzreal(d, q);

    // compute material matrix
    ComputeMaterialMatrix(mat, xyzeg, DLeg, Dweg);

    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int mmod = 0; mmod < (nsd_ * nsd_ + (msd_ - nsd_)); ++mmod)
        for (unsigned int d = 0; d < nsd_; ++d)
          dDwdxyzeg(mmod * nsd_ + d, q) += Nxyz(i * nsd_ + d, q) * Dwemod(mmod, i);
    }

    // compute matrices
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int j = 0; j < ndofs_; ++j)
      {
        for (unsigned int m = 0; m < msd_; ++m)
          ALL(m * ndofs_ + i, m * ndofs_ + j) += -N(i, q) * N(j, q) * fac(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            ALr(d * ndofs_ + i, j) += -Nxyz(i * nsd_ + e, q) * DLeg(e, d) * weg(e, q) /
                                      std::pow(reg(q), 2.0) * N(j, q) * fac(q);
        for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
          ALr((nsd_ + s) * ndofs_ + i, j) +=
              -Nxyz(i * nsd_ + VoigtP[s][0], q) * DLeg(nsd_ + s, nsd_ + s) * weg(VoigtP[s][1], q) /
                  std::pow(reg(q), 2.0) * N(j, q) * fac(q) -
              Nxyz(i * nsd_ + VoigtP[s][1], q) * DLeg(nsd_ + s, nsd_ + s) * weg(VoigtP[s][0], q) /
                  std::pow(reg(q), 2.0) * N(j, q) * fac(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            ALw(d * ndofs_ + i, e * ndofs_ + j) +=
                +Nxyz(i * nsd_ + e, q) * DLeg(e, d) / reg(q) * N(j, q) * fac(q);
        for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
        {
          ALw((nsd_ + s) * ndofs_ + i, VoigtP[s][0] * ndofs_ + j) +=
              +Nxyz(i * nsd_ + VoigtP[s][1], q) * DLeg(nsd_ + s, nsd_ + s) / reg(q) * N(j, q) *
              fac(q);
          ALw((nsd_ + s) * ndofs_ + i, VoigtP[s][1] * ndofs_ + j) +=
              +Nxyz(i * nsd_ + VoigtP[s][0], q) * DLeg(nsd_ + s, nsd_ + s) / reg(q) * N(j, q) *
              fac(q);
        }

        if (unsteady) Arr(i, j) += +invtimefac * N(i, q) * N(j, q) * fac(q);

        if (ale)
          for (unsigned int d = 0; d < nsd_; ++d)
            Arr(i, j) += +N(i, q) * dadxyzeg(d * nsd_ + d, q) * N(j, q) * fac(q);

        if (ale)
          for (unsigned int d = 0; d < nsd_; ++d)
            Arr(i, j) += +Nxyz(i * nsd_ + d, q) * aeg(d, q) * N(j, q) * fac(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          Arw(i, d * ndofs_ + j) += -Nxyz(i * nsd_ + d, q) * N(j, q) * fac(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            AwL(d * ndofs_ + i, e * ndofs_ + j) +=
                +N(i, q) * Dweg(d, e) * Nxyz(j * nsd_ + d, q) * fac(q);
        for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
        {
          AwL(VoigtP[s][0] * ndofs_ + i, (nsd_ + s) * ndofs_ + j) +=
              +N(i, q) * Dweg(nsd_ + s, nsd_ + s) * Nxyz(j * nsd_ + VoigtP[s][1], q) * fac(q);
          AwL(VoigtP[s][1] * ndofs_ + i, (nsd_ + s) * ndofs_ + j) +=
              +N(i, q) * Dweg(nsd_ + s, nsd_ + s) * Nxyz(j * nsd_ + VoigtP[s][0], q) * fac(q);
        }

        // TODO to test this term. It arises from the linearization of the variable viscosity term
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            AwL(d * ndofs_ + i, e * ndofs_ + j) +=
                +N(i, q) * dDwdxyzeg((d * nsd_ + e) * nsd_ + d, q) * N(j, q) * fac(q);
        for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
        {
          AwL(VoigtP[s][0] * ndofs_ + i, (nsd_ + s) * ndofs_ + j) +=
              +N(i, q) * dDwdxyzeg((nsd_ * nsd_ + s) * nsd_ + VoigtP[s][1], q) * N(j, q) * fac(q);
          AwL(VoigtP[s][1] * ndofs_ + i, (nsd_ + s) * ndofs_ + j) +=
              +N(i, q) * dDwdxyzeg((nsd_ * nsd_ + s) * nsd_ + VoigtP[s][0], q) * N(j, q) * fac(q);
        }

        if (convective)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              Awr(d * ndofs_ + i, j) += +Nxyz(i * nsd_ + e, q) * weg(d, q) * weg(e, q) /
                                        std::pow(reg(q), 2.0) * N(j, q) * fac(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          Awr(d * ndofs_ + i, j) += +1.0 / eps * N(i, q) * Nxyz(j * nsd_ + d, q) * fac(q);

        if (unsteady)
          for (unsigned int d = 0; d < nsd_; ++d)
            Aww(d * ndofs_ + i, d * ndofs_ + j) += +invtimefac * N(i, q) * N(j, q) * fac(q);

        if (ale)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              Aww(d * ndofs_ + i, d * ndofs_ + j) +=
                  +N(i, q) * dadxyzeg(e * nsd_ + e, q) * N(j, q) * fac(q);

        if (convective)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              Aww(d * ndofs_ + i, d * ndofs_ + j) +=
                  -Nxyz(i * nsd_ + e, q) * weg(e, q) / reg(q) * N(j, q) * fac(q);

        if (ale)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              Aww(d * ndofs_ + i, d * ndofs_ + j) +=
                  +Nxyz(i * nsd_ + e, q) * aeg(e, q) * N(j, q) * fac(q);

        if (convective)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              Aww(d * ndofs_ + i, e * ndofs_ + j) +=
                  -Nxyz(i * nsd_ + e, q) * weg(d, q) / reg(q) * N(j, q) * fac(q);
      }
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeFaceResidual(const int f,
    const Teuchos::RCP<MAT::Material>& mat, const std::vector<double>& val,
    const std::vector<double>& traceval, const std::vector<double>& alevel)
{
  // get material properties
  const MAT::WeaklyCompressibleFluid* actmat =
      static_cast<const MAT::WeaklyCompressibleFluid*>(mat.get());
  double eps = actmat->ComprCoeff();

  // get stabilization parameters
  const Teuchos::ParameterList& fluidparams = DRT::Problem::Instance()->FluidDynamicParams();
  double tau_r_ref = fluidparams.get<double>("STAB_DEN_REF");
  double tau_w_ref = fluidparams.get<double>("STAB_MOM_REF");
  tau_r = tau_r_ref / eps;
  tau_w = tau_w_ref;

  // ease notation
  DRT::UTILS::ShapeValuesInteriorOnFace Ni = shapesface_.shfunctI;
  LINALG::SerialDenseMatrix Nf = shapesface_.shfunct;
  LINALG::SerialDenseMatrix Nn = shapesface_.funct;
  LINALG::SerialDenseMatrix nxyz = shapesface_.normals;
  LINALG::SerialDenseVector facf = shapesface_.jfac;
  unsigned int nfqpoints = shapesface_.nqpoints_;
  unsigned int nfdofs = shapesface_.nfdofs_;
  unsigned int nfn = shapesface_.nfn_;
  std::vector<std::vector<int>> faceNodeOrder = shapesface_.faceNodeOrder;

  // extract interior nodal values
  Epetra_SerialDenseMatrix Le(msd_, ndofs_);
  Epetra_SerialDenseVector re(ndofs_);
  Epetra_SerialDenseMatrix we(nsd_, ndofs_);
  Epetra_SerialDenseVector pe(ndofs_);

  for (unsigned int i = 0; i < ndofs_; ++i)
  {
    for (unsigned int m = 0; m < msd_; ++m) Le(m, i) = val[m * ndofs_ + i];

    re(i) = val[msd_ * ndofs_ + i];

    for (unsigned int d = 0; d < nsd_; ++d) we(d, i) = val[(msd_ + 1 + d) * ndofs_ + i];

    pe(i) = actmat->ComputePressure(re[i]);
  }

  // extract ale nodal values
  Epetra_SerialDenseMatrix aef(nsd_, nfn);

  if (ale)
    for (unsigned int n = 0; n < nfn; ++n)
      for (unsigned int d = 0; d < nsd_; ++d) aef(d, n) = alevel[faceNodeOrder[f][n] * nsd_ + d];

  // extract trace nodal value
  Epetra_SerialDenseVector rhatef(nfdofs);
  Epetra_SerialDenseMatrix whatef(nsd_, nfdofs);
  for (unsigned int i = 0; i < nfdofs; ++i)
  {
    rhatef(i) = traceval[f * nfdofs * (1 + nsd_) + i];

    for (unsigned int d = 0; d < nsd_; ++d)
      whatef(d, i) = traceval[f * nfdofs * (1 + nsd_) + (1 + d) * nfdofs + i];
  }

  // initialize values interpolated on gauss points
  LINALG::Matrix<nsd_, 1> xyzefg;
  Epetra_SerialDenseMatrix Lefg(msd_, nfqpoints);
  Epetra_SerialDenseVector refg(nfqpoints);
  Epetra_SerialDenseMatrix wefg(nsd_, nfqpoints);
  Epetra_SerialDenseVector phatefg(nfqpoints);
  Epetra_SerialDenseMatrix DLefg(msd_, msd_);
  Epetra_SerialDenseMatrix Dwefg(msd_, msd_);
  rhatefg.Shape(nfqpoints, 1);
  whatefg.Shape(nsd_, nfqpoints);
  aefg.Shape(nsd_, nfqpoints);

  // loop over quadrature points
  for (unsigned int q = 0; q < nfqpoints; ++q)
  {
    // interpolate interior values on gauss points
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int m = 0; m < msd_; ++m) Lefg(m, q) += Ni(i, q) * Le(m, i);

      refg(q) += Ni(i, q) * re(i);

      for (unsigned int d = 0; d < nsd_; ++d) wefg(d, q) += Ni(i, q) * we(d, i);
    }

    // interpolate ale values on gauss points
    if (ale)
      for (unsigned int n = 0; n < nfn; ++n)
        for (unsigned int d = 0; d < nsd_; ++d) aefg(d, q) += Nn(n, q) * aef(d, n);

    // interpolate trace values on gauss points
    for (unsigned int i = 0; i < nfdofs; ++i)
    {
      rhatefg(q) += Nf(i, q) * rhatef(i);

      for (unsigned int d = 0; d < nsd_; ++d) whatefg(d, q) += Nf(i, q) * whatef(d, i);

      phatefg(q) = actmat->ComputePressure(rhatefg(q));
    }

    // interpolate values on gauss points
    for (unsigned int d = 0; d < nsd_; ++d) xyzefg(d) = shapesface_.xyzreal(d, q);

    // compute material matrix
    ComputeMaterialMatrix(mat, xyzefg, DLefg, Dwefg);

    // compute residuals (interior contributions)
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int d = 0; d < nsd_; ++d)
        for (unsigned int e = 0; e < nsd_; ++e)
          RL(d * ndofs_ + i) +=
              +Ni(i, q) * nxyz(e, q) * DLefg(e, d) * whatefg(e, q) / rhatefg(q) * facf(q);
      for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
        RL((nsd_ + s) * ndofs_ + i) +=
            +Ni(i, q) * nxyz(VoigtP[s][0], q) * DLefg(nsd_ + s, nsd_ + s) *
                whatefg(VoigtP[s][1], q) / rhatefg(q) * facf(q) +
            Ni(i, q) * nxyz(VoigtP[s][1], q) * DLefg(nsd_ + s, nsd_ + s) *
                whatefg(VoigtP[s][0], q) / rhatefg(q) * facf(q);

      for (unsigned int d = 0; d < nsd_; ++d)
        Rr(i) += -Ni(i, q) * whatefg(d, q) * nxyz(d, q) * facf(q);

      if (ale)
        for (unsigned int d = 0; d < nsd_; ++d)
          Rr(i) += +Ni(i, q) * rhatefg(q) * aefg(d, q) * nxyz(d, q) * facf(q);

      Rr(i) += -Ni(i, q) * tau_r * refg(q) * facf(q);

      Rr(i) += +Ni(i, q) * tau_r * rhatefg(q) * facf(q);

      if (convective)
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            Rw(d * ndofs_ + i) +=
                -Ni(i, q) * whatefg(d, q) * whatefg(e, q) / rhatefg(q) * nxyz(e, q) * facf(q);

      if (ale)
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            Rw(d * ndofs_ + i) += +Ni(i, q) * whatefg(d, q) * aefg(e, q) * nxyz(e, q) * facf(q);

      for (unsigned int d = 0; d < nsd_; ++d)
        Rw(d * ndofs_ + i) += -Ni(i, q) * tau_w * wefg(d, q) * facf(q);

      for (unsigned int d = 0; d < nsd_; ++d)
        Rw(d * ndofs_ + i) += +Ni(i, q) * tau_w * whatefg(d, q) * facf(q);
    }

    // compute residuals (face contributions)
    for (unsigned int i = 0; i < nfdofs; ++i)
    {
      RR(f * nfdofs + i) += -Nf(i, q) * tau_r * refg(q) * facf(q);

      RR(f * nfdofs + i) += +Nf(i, q) * tau_r * rhatefg(q) * facf(q);

      for (unsigned int d = 0; d < nsd_; ++d)
        for (unsigned int e = 0; e < nsd_; ++e)
          RW((f * nsd_ + d) * nfdofs + i) +=
              -Nf(i, q) * Dwefg(d, e) * Lefg(e, q) * nxyz(d, q) * facf(q);
      for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
      {
        RW((f * nsd_ + VoigtP[s][0]) * nfdofs + i) += -Nf(i, q) * Dwefg(nsd_ + s, nsd_ + s) *
                                                      Lefg(nsd_ + s, q) * nxyz(VoigtP[s][1], q) *
                                                      facf(q);
        RW((f * nsd_ + VoigtP[s][1]) * nfdofs + i) += -Nf(i, q) * Dwefg(nsd_ + s, nsd_ + s) *
                                                      Lefg(nsd_ + s, q) * nxyz(VoigtP[s][0], q) *
                                                      facf(q);
      }

      for (unsigned int d = 0; d < nsd_; ++d)
        RW((f * nsd_ + d) * nfdofs + i) += -Nf(i, q) * phatefg(q) * nxyz(d, q) * facf(q);

      for (unsigned int d = 0; d < nsd_; ++d)
        RW((f * nsd_ + d) * nfdofs + i) += -Nf(i, q) * tau_w * wefg(d, q) * facf(q);

      for (unsigned int d = 0; d < nsd_; ++d)
        RW((f * nsd_ + d) * nfdofs + i) += +Nf(i, q) * tau_w * whatefg(d, q) * facf(q);

      //  TODO Neumann term
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeFaceMatrices(
    const int f, const Teuchos::RCP<MAT::Material>& mat)
{
  // get material properties
  const MAT::WeaklyCompressibleFluid* actmat =
      static_cast<const MAT::WeaklyCompressibleFluid*>(mat.get());
  double eps = actmat->ComprCoeff();

  // ease notation
  DRT::UTILS::ShapeValuesInteriorOnFace Ni = shapesface_.shfunctI;
  LINALG::SerialDenseMatrix Nf = shapesface_.shfunct;
  LINALG::SerialDenseMatrix nxyz = shapesface_.normals;
  LINALG::SerialDenseVector facf = shapesface_.jfac;
  unsigned int nfqpoints = shapesface_.nqpoints_;
  unsigned int nfdofs = shapesface_.nfdofs_;

  // initialize values interpolated on gauss points
  LINALG::Matrix<nsd_, 1> xyzefg;
  Epetra_SerialDenseMatrix DLefg(msd_, msd_);
  Epetra_SerialDenseMatrix Dwefg(msd_, msd_);

  // loop over quadrature points
  for (unsigned int q = 0; q < nfqpoints; ++q)
  {
    // interpolate values on gauss points
    for (unsigned int d = 0; d < nsd_; ++d) xyzefg(d) = shapesface_.xyzreal(d, q);

    // compute material matrix
    ComputeMaterialMatrix(mat, xyzefg, DLefg, Dwefg);

    // compute matrices (interior - interior contributions)
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int j = 0; j < ndofs_; ++j)
      {
        Arr(i, j) += +Ni(i, q) * tau_r * Ni(j, q) * facf(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          Aww(d * ndofs_ + i, d * ndofs_ + j) += +Ni(i, q) * tau_w * Ni(j, q) * facf(q);
      }
    }

    // compute matrices (interior - face contributions)
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int j = 0; j < nfdofs; ++j)
      {
        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            ALR(d * ndofs_ + i, f * nfdofs + j) += +Ni(i, q) * nxyz(e, q) * DLefg(e, d) *
                                                   whatefg(e, q) / std::pow(rhatefg(q), 2.0) *
                                                   Nf(j, q) * facf(q);
        for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
          ALR((nsd_ + s) * ndofs_ + i, f * nfdofs + j) +=
              +Ni(i, q) * nxyz(VoigtP[s][0], q) * DLefg(nsd_ + s, nsd_ + s) *
                  whatefg(VoigtP[s][1], q) / std::pow(rhatefg(q), 2.0) * Nf(j, q) * facf(q) +
              Ni(i, q) * nxyz(VoigtP[s][1], q) * DLefg(nsd_ + s, nsd_ + s) *
                  whatefg(VoigtP[s][0], q) / std::pow(rhatefg(q), 2.0) * Nf(j, q) * facf(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            ALW(d * ndofs_ + i, (f * nsd_ + e) * nfdofs + j) +=
                -Ni(i, q) * nxyz(e, q) * DLefg(e, d) / rhatefg(q) * Nf(j, q) * facf(q);
        for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
        {
          ALW((nsd_ + s) * ndofs_ + i, (f * nsd_ + VoigtP[s][0]) * nfdofs + j) +=
              -Ni(i, q) * nxyz(VoigtP[s][1], q) * DLefg(nsd_ + s, nsd_ + s) / rhatefg(q) *
              Nf(j, q) * facf(q);
          ALW((nsd_ + s) * ndofs_ + i, (f * nsd_ + VoigtP[s][1]) * nfdofs + j) +=
              -Ni(i, q) * nxyz(VoigtP[s][0], q) * DLefg(nsd_ + s, nsd_ + s) / rhatefg(q) *
              Nf(j, q) * facf(q);
        }

        if (ale)
          for (unsigned int d = 0; d < nsd_; ++d)
            ArR(i, f * nfdofs + j) += -Ni(i, q) * aefg(d, q) * nxyz(d, q) * Nf(j, q) * facf(q);

        ArR(i, f * nfdofs + j) += -Ni(i, q) * tau_r * Nf(j, q) * facf(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          ArW(i, (f * nsd_ + d) * nfdofs + j) += +Ni(i, q) * nxyz(d, q) * Nf(j, q) * facf(q);

        if (convective)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              AwR(d * ndofs_ + i, f * nfdofs + j) += -Ni(i, q) * whatefg(d, q) * whatefg(e, q) /
                                                     std::pow(rhatefg(q), 2.0) * nxyz(e, q) *
                                                     Nf(j, q) * facf(q);

        if (convective)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              AwW(d * ndofs_ + i, (f * nsd_ + d) * nfdofs + j) +=
                  +Ni(i, q) * whatefg(e, q) / rhatefg(q) * nxyz(e, q) * Nf(j, q) * facf(q);

        if (ale)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              AwW(d * ndofs_ + i, (f * nsd_ + d) * nfdofs + j) +=
                  -Ni(i, q) * aefg(e, q) * nxyz(e, q) * Nf(j, q) * facf(q);

        if (convective)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              AwW(d * ndofs_ + i, (f * nsd_ + e) * nfdofs + j) +=
                  +Ni(i, q) * whatefg(d, q) / rhatefg(q) * nxyz(e, q) * Nf(j, q) * facf(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          AwW(d * ndofs_ + i, (f * nsd_ + d) * nfdofs + j) +=
              -Ni(i, q) * tau_w * Nf(j, q) * facf(q);
      }
    }

    // compute matrices (face - interior contributions)
    for (unsigned int i = 0; i < nfdofs; ++i)
    {
      for (unsigned int j = 0; j < ndofs_; ++j)
      {
        ARr(f * nfdofs + i, j) += +Nf(i, q) * tau_r * Ni(j, q) * facf(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          for (unsigned int e = 0; e < nsd_; ++e)
            AWL((f * nsd_ + d) * nfdofs + i, e * ndofs_ + j) +=
                +Nf(i, q) * nxyz(d, q) * Dwefg(d, e) * Ni(j, q) * facf(q);
        for (unsigned int s = 0; s < (msd_ - nsd_); ++s)
        {
          AWL((f * nsd_ + VoigtP[s][0]) * nfdofs + i, (nsd_ + s) * ndofs_ + j) +=
              +Nf(i, q) * nxyz(VoigtP[s][1], q) * Dwefg(nsd_ + s, nsd_ + s) * Ni(j, q) * facf(q);
          AWL((f * nsd_ + VoigtP[s][1]) * nfdofs + i, (nsd_ + s) * ndofs_ + j) +=
              +Nf(i, q) * nxyz(VoigtP[s][0], q) * Dwefg(nsd_ + s, nsd_ + s) * Ni(j, q) * facf(q);
        }

        for (unsigned int d = 0; d < nsd_; ++d)
          AWw((f * nsd_ + d) * nfdofs + i, d * ndofs_ + j) +=
              +Nf(i, q) * tau_w * Ni(j, q) * facf(q);
      }
    }

    // compute matrices (face - face contributions)
    for (unsigned int i = 0; i < nfdofs; ++i)
    {
      for (unsigned int j = 0; j < nfdofs; ++j)
      {
        ARR(f * nfdofs + i, f * nfdofs + j) += -Nf(i, q) * tau_r * Nf(j, q) * facf(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          AWR((f * nsd_ + d) * nfdofs + i, f * nfdofs + j) +=
              +Nf(i, q) * nxyz(d, q) / eps * Nf(j, q) * facf(q);

        for (unsigned int d = 0; d < nsd_; ++d)
          AWW((f * nsd_ + d) * nfdofs + i, (f * nsd_ + d) * nfdofs + j) +=
              -Nf(i, q) * tau_w * Nf(j, q) * facf(q);
      }
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeLocalResidual()
{
  // fill vector
  for (unsigned int i = 0; i < ndofs_; ++i)
  {
    for (unsigned int m = 0; m < msd_; ++m) Rlocal(m * ndofs_ + i) = RL(m * ndofs_ + i);

    Rlocal(msd_ * ndofs_ + i) = Rr(i);

    for (unsigned int d = 0; d < nsd_; ++d)
      Rlocal((msd_ + 1 + d) * ndofs_ + i) = Rw(d * ndofs_ + i);
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeGlobalResidual(
    DRT::ELEMENTS::Fluid& ele)
{
  // loop in faces
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    // evaluate face
    shapesface_.EvaluateFace(ele, f);

    // ease notation
    unsigned int nfdofs = shapesface_.nfdofs_;

    // fill vector
    for (unsigned int i = 0; i < nfdofs; ++i)
    {
      Rglobal(f * nfdofs * (1 + nsd_) + i) = RR(f * nfdofs + i);

      for (unsigned int d = 0; d < nsd_; ++d)
        Rglobal(f * nfdofs * (1 + nsd_) + (1 + d) * nfdofs + i) = RW((f * nsd_ + d) * nfdofs + i);
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeLocalLocalMatrix()
{
  // fill matrix
  for (unsigned int i = 0; i < ndofs_; ++i)
  {
    for (unsigned int j = 0; j < ndofs_; ++j)
    {
      for (unsigned int m = 0; m < msd_; ++m)
      {
        for (unsigned int n = 0; n < msd_; ++n)
          Klocallocal(m * ndofs_ + i, n * ndofs_ + j) = ALL(m * ndofs_ + i, n * ndofs_ + j);

        Klocallocal(m * ndofs_ + i, msd_ * ndofs_ + j) = ALr(m * ndofs_ + i, j);

        for (unsigned int d = 0; d < nsd_; ++d)
          Klocallocal(m * ndofs_ + i, (msd_ + 1 + d) * ndofs_ + j) =
              ALw(m * ndofs_ + i, d * ndofs_ + j);
      }

      Klocallocal(msd_ * ndofs_ + i, msd_ * ndofs_ + j) = Arr(i, j);

      for (unsigned int d = 0; d < nsd_; ++d)
        Klocallocal(msd_ * ndofs_ + i, (msd_ + 1 + d) * ndofs_ + j) = Arw(i, d * ndofs_ + j);

      for (unsigned int d = 0; d < nsd_; ++d)
      {
        for (unsigned int m = 0; m < msd_; ++m)
          Klocallocal((msd_ + 1 + d) * ndofs_ + i, m * ndofs_ + j) =
              AwL(d * ndofs_ + i, m * ndofs_ + j);

        Klocallocal((msd_ + 1 + d) * ndofs_ + i, msd_ * ndofs_ + j) = Awr(d * ndofs_ + i, j);

        for (unsigned int e = 0; e < nsd_; ++e)
          Klocallocal((msd_ + 1 + d) * ndofs_ + i, (msd_ + 1 + e) * ndofs_ + j) =
              Aww(d * ndofs_ + i, e * ndofs_ + j);
      }
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeLocalGlobalMatrix(
    DRT::ELEMENTS::Fluid& ele)
{
  // loop in faces
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    // evaluate face
    shapesface_.EvaluateFace(ele, f);

    // ease notation
    unsigned int nfdofs = shapesface_.nfdofs_;

    // fill matrix
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int j = 0; j < nfdofs; ++j)
      {
        for (unsigned int m = 0; m < msd_; ++m)
        {
          Klocalglobal(m * ndofs_ + i, f * nfdofs * (1 + nsd_) + j) =
              ALR(m * ndofs_ + i, f * nfdofs + j);

          for (unsigned int d = 0; d < nsd_; ++d)
            Klocalglobal(m * ndofs_ + i, f * nfdofs * (1 + nsd_) + (1 + d) * nfdofs + j) =
                ALW(m * ndofs_ + i, (f * nsd_ + d) * nfdofs + j);
        }

        Klocalglobal(msd_ * ndofs_ + i, f * nfdofs * (1 + nsd_) + j) = ArR(i, f * nfdofs + j);

        for (unsigned int d = 0; d < nsd_; ++d)
          Klocalglobal(msd_ * ndofs_ + i, f * nfdofs * (1 + nsd_) + (1 + d) * nfdofs + j) =
              ArW(i, (f * nsd_ + d) * nfdofs + j);

        for (unsigned int d = 0; d < nsd_; ++d)
        {
          Klocalglobal((msd_ + 1 + d) * ndofs_ + i, f * nfdofs * (1 + nsd_) + j) =
              AwR(d * ndofs_ + i, f * nfdofs + j);

          for (unsigned int e = 0; e < nsd_; ++e)
            Klocalglobal(
                (msd_ + 1 + d) * ndofs_ + i, f * nfdofs * (1 + nsd_) + (1 + e) * nfdofs + j) =
                AwW(d * ndofs_ + i, (f * nsd_ + e) * nfdofs + j);
        }
      }
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeGlobalLocalMatrix(
    DRT::ELEMENTS::Fluid& ele)
{
  // loop in faces
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    // evaluate face
    shapesface_.EvaluateFace(ele, f);

    // ease notation
    unsigned int nfdofs = shapesface_.nfdofs_;

    // fill matrix
    for (unsigned int i = 0; i < nfdofs; ++i)
    {
      for (unsigned int j = 0; j < ndofs_; ++j)
      {
        Kgloballocal(f * nfdofs * (1 + nsd_) + i, msd_ * ndofs_ + j) = ARr(f * nfdofs + i, j);

        for (unsigned int d = 0; d < nsd_; ++d)
        {
          for (unsigned int m = 0; m < msd_; ++m)
            Kgloballocal(f * nfdofs * (1 + nsd_) + (1 + d) * nfdofs + i, m * ndofs_ + j) =
                AWL((f * nsd_ + d) * nfdofs + i, m * ndofs_ + j);

          for (unsigned int e = 0; e < nsd_; ++e)
            Kgloballocal(f * nfdofs * (1 + nsd_) + (1 + d) * nfdofs + i,
                (msd_ + 1 + e) * ndofs_ + j) = AWw((f * nsd_ + d) * nfdofs + i, e * ndofs_ + j);
        }
      }
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::ComputeGlobalGlobalMatrix(
    DRT::ELEMENTS::Fluid& ele)
{
  // loop in faces
  for (unsigned int f = 0; f < nfaces_; ++f)
  {
    // evaluate face
    shapesface_.EvaluateFace(ele, f);

    // ease notation
    unsigned int nfdofs = shapesface_.nfdofs_;

    // fill matrix
    for (unsigned int i = 0; i < nfdofs; ++i)
    {
      for (unsigned int j = 0; j < nfdofs; ++j)
      {
        Kglobalglobal(f * nfdofs * (1 + nsd_) + i, f * nfdofs * (1 + nsd_) + j) =
            ARR(f * nfdofs + i, f * nfdofs + j);

        for (unsigned int d = 0; d < nsd_; ++d)
        {
          Kglobalglobal(f * nfdofs * (1 + nsd_) + (1 + d) * nfdofs + i,
              f * nfdofs * (1 + nsd_) + j) = AWR((f * nsd_ + d) * nfdofs + i, f * nfdofs + j);

          for (unsigned int e = 0; e < nsd_; ++e)
            Kglobalglobal(f * nfdofs * (1 + nsd_) + (1 + d) * nfdofs + i,
                f * nfdofs * (1 + nsd_) + (1 + e) * nfdofs + j) =
                AWW((f * nsd_ + d) * nfdofs + i, (f * nsd_ + e) * nfdofs + j);
        }
      }
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::InvertLocalLocalMatrix()
{
  KlocallocalInv = Klocallocal;
  KlocallocalInvSolver.SetMatrix(KlocallocalInv);
  int err = KlocallocalInvSolver.Invert();
  if (err != 0) dserror("Inversion of local-local matrix failed with errorcode %d", err);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::CondenseLocalResidual(
    Epetra_SerialDenseVector& eleVec)
{
  // initialize element vector
  eleVec.Shape((1 + nsd_) * ndofsfaces_, 1);

  // create auxiliary vector
  Epetra_SerialDenseVector eleVecAux;
  eleVecAux.Shape((msd_ + 1 + nsd_) * ndofs_, 1);

  // compute element vector
  eleVecAux.Multiply('N', 'N', 1.0, KlocallocalInv, Rlocal, 0.0);
  eleVec = Rglobal;
  eleVec.Multiply('N', 'N', -1.0, Kgloballocal, eleVecAux, 1.0);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::CondenseLocalMatrix(
    Epetra_SerialDenseMatrix& eleMat)
{
  // initialize element matrix
  eleMat.Shape((1 + nsd_) * ndofsfaces_, (1 + nsd_) * ndofsfaces_);

  // create auxiliary matrix
  Epetra_SerialDenseMatrix eleMatAux;
  eleMatAux.Shape((msd_ + 1 + nsd_) * ndofs_, (1 + nsd_) * ndofsfaces_);

  // compute element matrix
  eleMatAux.Multiply('N', 'N', 1.0, KlocallocalInv, Klocalglobal, 0.0);
  eleMat = Kglobalglobal;
  eleMat.Multiply('N', 'N', -1.0, Kgloballocal, eleMatAux, 1.0);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDGWeakComp<distype>::LocalSolver::PrintMatricesAndResiduals(
    DRT::ELEMENTS::Fluid& ele, Epetra_SerialDenseVector& eleVec, Epetra_SerialDenseMatrix& eleMat)
{
  // element
  std::cout << "\n\n Element number = " << ele.Id() + 1;

  // matrices
  std::cout << "\n\n ALL = \n\n";
  ALL.Print(std::cout);
  std::cout << "\n\n ALr = \n\n";
  ALr.Print(std::cout);
  std::cout << "\n\n ALw = \n\n";
  ALw.Print(std::cout);
  std::cout << "\n\n ALR = \n\n";
  ALR.Print(std::cout);
  std::cout << "\n\n ALW = \n\n";
  ALW.Print(std::cout);
  std::cout << "\n\n Arr = \n\n";
  Arr.Print(std::cout);
  std::cout << "\n\n Arw = \n\n";
  Arw.Print(std::cout);
  std::cout << "\n\n ArR = \n\n";
  ArR.Print(std::cout);
  std::cout << "\n\n ArW = \n\n";
  ArW.Print(std::cout);
  std::cout << "\n\n AwL = \n\n";
  AwL.Print(std::cout);
  std::cout << "\n\n Awr = \n\n";
  Awr.Print(std::cout);
  std::cout << "\n\n Aww = \n\n";
  Aww.Print(std::cout);
  std::cout << "\n\n AwR = \n\n";
  AwR.Print(std::cout);
  std::cout << "\n\n AwW = \n\n";
  AwW.Print(std::cout);
  std::cout << "\n\n ARr = \n\n";
  ARr.Print(std::cout);
  std::cout << "\n\n ARR = \n\n";
  ARR.Print(std::cout);
  std::cout << "\n\n AWL = \n\n";
  AWL.Print(std::cout);
  std::cout << "\n\n AWw = \n\n";
  AWw.Print(std::cout);
  std::cout << "\n\n AWR = \n\n";
  AWR.Print(std::cout);
  std::cout << "\n\n AWW = \n\n";
  AWW.Print(std::cout);

  // residuals
  std::cout << "\n\n RL = \n\n";
  RL.Print(std::cout);
  std::cout << "\n\n Rr = \n\n";
  Rr.Print(std::cout);
  std::cout << "\n\n Rw = \n\n";
  Rw.Print(std::cout);
  std::cout << "\n\n RR = \n\n";
  RR.Print(std::cout);
  std::cout << "\n\n RW = \n\n";
  RW.Print(std::cout);

  // local/global matrices/vectors
  std::cout << "\n\n Klocallocal = \n\n";
  Klocallocal.Print(std::cout);
  std::cout << "\n\n Klocalglobal = \n\n";
  Klocalglobal.Print(std::cout);
  std::cout << "\n\n Kgloballocal = \n\n";
  Kgloballocal.Print(std::cout);
  std::cout << "\n\n Kglobalglobal = \n\n";
  Kglobalglobal.Print(std::cout);
  std::cout << "\n\n Rlocal = \n\n";
  Rlocal.Print(std::cout);
  std::cout << "\n\n Rglobal = \n\n";
  Rglobal.Print(std::cout);
  std::cout << "\n\n KlocallocalInv = \n\n";
  KlocallocalInv.Print(std::cout);

  // element vector and matrix
  std::cout << "\n\n eleVec = \n\n";
  eleVec.Print(std::cout);
  std::cout << "\n\n eleMat = \n\n";
  eleMat.Print(std::cout);
}



// explicit instantiation of template classes
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::wedge15>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcHDGWeakComp<DRT::Element::nurbs27>;
