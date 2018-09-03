/*--------------------------------------------------------------------------*/
/*!
\file acou_sol_ele_calc.cpp
\brief

<pre>
\level 2

\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "acou_sol_ele_calc.H"
#include "acou_ele_calc.H"
#include "acou_ele_action.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_geometry/position_array.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_mat/acoustic_sol.H"

#include <Epetra_SerialDenseSolver.h>

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
DRT::ELEMENTS::AcouSolEleCalc<distype>::AcouSolEleCalc()
{
}

/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouSolEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou* ele,
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
int DRT::ELEMENTS::AcouSolEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix&,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector&,
    bool offdiag)
{
  // check if this is an hdg element and init completepoly
  if (const DRT::ELEMENTS::AcouSol* hdgele = dynamic_cast<const DRT::ELEMENTS::AcouSol*>(ele))
    usescompletepoly_ = hdgele->UsesCompletePolynomialSpace();
  else
    dserror("cannot cast element to acousol element");
  InitializeShapes(ele);

  const ACOU::Action action = DRT::INPUT::get<ACOU::Action>(params, "action");
  bool updateonly = false;
  shapes_->Evaluate(*ele);

  switch (action)
  {
    case ACOU::project_field:
    {
      if (mat->MaterialType() != INPAR::MAT::m_acousticsolmat)
        dserror("for physical type 'solid' please supply MAT_AcousticSol");
      ElementInit(ele, params);
      ProjectField(ele, params, elevec1, elevec2);
      break;
    }
    case ACOU::project_dirich_field:
    {
      if (mat->MaterialType() != INPAR::MAT::m_acousticsolmat)
        dserror("for physical type 'solid' please supply MAT_AcousticSol");
      ElementInit(ele, params);
      ProjectDirichField(ele, params, elevec1);
      break;
    }
    case ACOU::project_optical_field:
    {
      if (mat->MaterialType() != INPAR::MAT::m_acousticsolmat)
        dserror("for physical type 'solid' please supply MAT_AcousticSol");
      ProjectOpticalField(ele, params, elevec2);
      break;
    }
    case ACOU::ele_init:
    {
      ElementInit(ele, params);
      break;
    }
    case ACOU::fill_restart_vecs:
    {
      ReadGlobalVectors(ele, discretization, lm);
      FillRestartVectors(ele, discretization);
      break;
    }
    case ACOU::ele_init_from_restart:
    {
      ElementInit(ele, params);
      ElementInitFromRestart(ele, discretization);
      break;
    }
    case ACOU::interpolate_hdg_to_node:
    {
      bool writestress = params.get<bool>("writestress");
      ReadGlobalVectors(ele, discretization, lm);
      NodeBasedValues(mat, ele, elevec1, writestress);
      break;
    }
    case ACOU::interpolate_psi_to_node:
    {
      double dt = params.get<double>("dt");
      ReadGlobalVectors(ele, discretization, lm);
      NodeBasedPsi(mat, ele, elevec1, dt);
      break;
    }
    case ACOU::calc_acou_error:
    {
      ReadGlobalVectors(ele, discretization, lm);
      ComputeError(ele, params, elevec1);
      break;
    }
    case ACOU::calc_abc:
    {
      int face = params.get<int>("face");
      shapesface_->EvaluateFace(*ele, face);
      // note: absorbing bcs are treated fully implicit!
      localSolver_->ComputeAbsorbingBC(ele, params, mat, face, elemat1, elevec1);
      break;
    }
    case ACOU::calc_pressuremon:
    {
      // this action is not needed for solid since the source is calculated in the update routine
      dserror("solids do not need this action, why are you here?");
      break;
    }
    case ACOU::calc_systemmat_and_residual:
    {
      const bool resonly = params.get<bool>("resonly");
      double dt = params.get<double>("dt");
      dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

      ReadGlobalVectors(ele, discretization, lm);
      zeroMatrix(elevec1);
      ComputeMatrices(mat, *ele, dt, dyna_);

      if (!resonly)
      {
        localSolver_->CondenseLocalPart(elemat1, dyna_);
      }

      localSolver_->ComputeResidual(
          params, elevec1, interiorGradVelnp_, interiorVelnp_, interiorPressnp_, traceVal_, dyna_);

      break;
    }
    case ACOU::update_secondary_solution:
      updateonly = true;  // no break here!!!
    case ACOU::update_secondary_solution_and_calc_residual:
    {
      bool errormaps = params.get<bool>("errormaps");
      const bool allelesequal = params.get<bool>("allelesequal");
      double dt = params.get<double>("dt");
      dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

      ReadGlobalVectors(ele, discretization, lm);

      zeroMatrix(elevec1);
      if (!allelesequal) ComputeMatrices(mat, *ele, dt, dyna_);

      UpdateInteriorVariablesAndComputeResidual(params, ele, elevec1, dt, errormaps, updateonly);

      break;
    }
    case ACOU::get_gauss_points:
    {
      int rows = shapes_->xyzreal.M();
      int cols = shapes_->xyzreal.N();
      elemat1.Shape(rows, cols);

      for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c) elemat1(r, c) = shapes_->xyzreal(r, c);

      break;
    }
    default:
      dserror("unknown action supplied %d", action);
      break;
  }  // switch(action)

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::InitializeShapes(const DRT::ELEMENTS::Acou* ele)
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
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
  }
  // TODO: check distype

  if (localSolver_ == Teuchos::null)
    localSolver_ = Teuchos::rcp(new LocalSolver(ele, *shapes_, *shapesface_, usescompletepoly_));
  else if (localSolver_->ndofs_ != shapes_->ndofs_)
    localSolver_ = Teuchos::rcp(new LocalSolver(ele, *shapes_, *shapesface_, usescompletepoly_));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouSolEleCalc<distype>* DRT::ELEMENTS::AcouSolEleCalc<distype>::Instance(
    bool create)
{
  static AcouSolEleCalc<distype>* instance;
  if (create)
  {
    if (instance == NULL)
    {
      instance = new AcouSolEleCalc<distype>();
    }
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}

/*----------------------------------------------------------------------*
 * Constructor LocalSolver
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::LocalSolver(const DRT::ELEMENTS::Acou* ele,
    const DRT::UTILS::ShapeValues<distype>& shapeValues,
    DRT::UTILS::ShapeValuesFace<distype>& shapeValuesFace, bool completepoly)
    : ndofs_(shapeValues.ndofs_), shapes_(shapeValues), shapesface_(shapeValuesFace)
{
  int onfdofs = 0;
  for (unsigned int i = 0; i < nfaces_; ++i)
  {
    shapesface_.EvaluateFace(*ele, i);
    onfdofs += shapesface_.nfdofs_;
  }
  onfdofs *= nsd_;

  amat.Shape(nsd_ * nsd_ * ndofs_, nsd_ * nsd_ * ndofs_);
  invamat.Shape(ndofs_, ndofs_);
  bmat.Shape(nsd_ * nsd_ * ndofs_, nsd_ * ndofs_);
  cmat.Shape(nsd_ * nsd_ * ndofs_, onfdofs);
  dmat.Shape(nsd_ * ndofs_, nsd_ * nsd_ * ndofs_);
  emat.Shape(nsd_ * ndofs_, nsd_ * ndofs_);
  ehatmat.Shape(nsd_ * ndofs_, nsd_ * ndofs_);
  fmat.Shape(nsd_ * ndofs_, ndofs_);
  gmat.Shape(nsd_ * ndofs_, onfdofs);
  hmat.Shape(ndofs_, nsd_ * ndofs_);
  imat.Shape(ndofs_, ndofs_);
  jmat.Shape(ndofs_, onfdofs);
  kmat.Shape(onfdofs, nsd_ * nsd_ * ndofs_);
  lmat.Shape(onfdofs, nsd_ * ndofs_);
  mmat.Shape(onfdofs, ndofs_);
  nmat.Shape(onfdofs, onfdofs);
}


/*----------------------------------------------------------------------*
 * Element init // TODO: implement in AcouEleInterface
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::ElementInit(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params)
{
  DRT::ELEMENTS::AcouSol* solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);

  // each element has to store the interior vectors by itseld, p-adaptivity or not
  // so, shape it, as you need it
  if (params.isParameter("dynamic type"))
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");
  else
    dyna_ = INPAR::ACOU::acou_impleuler;
  solele->eleinteriorGradVelnp_.Shape(shapes_->ndofs_ * nsd_ * nsd_, 1);
  solele->eleinteriorVelnp_.Shape(shapes_->ndofs_ * nsd_, 1);
  solele->eleinteriorPressnp_.Shape(shapes_->ndofs_, 1);

  bool padaptivity = false;
  if (params.isParameter("padaptivity")) padaptivity = params.get<bool>("padaptivity");
  if (padaptivity) solele->elenodeTrace_.Shape(nen_, 1);

  return;
}  // ElementInit

/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouSolEleCalc<distype>::ProjectField(DRT::ELEMENTS::Acou* ele,
    Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2)
{
  DRT::ELEMENTS::AcouSol* solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 || unsigned(elevec2.M()) == (nsd_ * nsd_ + nsd_ + 1) * shapes_->ndofs_,
      "Wrong size in project vector 2");

  // get function
  const int* start_func = params.getPtr<int>("funct");

  // internal variables
  if (elevec2.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(shapes_->nqpoints_, nsd_ * nsd_ + nsd_ + 1);
    Epetra_SerialDenseMatrix massPart(shapes_->ndofs_, shapes_->nqpoints_);

    for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
    {
      const double fac = shapes_->jfac(q);
      const double sqrtfac = std::sqrt(fac);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_->xyzreal(d, q);  // coordinates of quadrature point in real coordinates
      double p;
      double gradient[nsd_];
      double velgrad[nsd_ * nsd_];

      dsassert(start_func != NULL, "startfuncno not set for initial value");
      EvaluateAll(*start_func, xyz, p, gradient, velgrad, 0.0);  // p at quadrature point

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      {
        massPart(i, q) = shapes_->shfunct(i, q) * sqrtfac;
        localMat(i, nsd_) += shapes_->shfunct(i, q) * p * fac;
        for (unsigned int j = 0; j < nsd_; ++j)
          localMat(i, j) += shapes_->shfunct(i, q) * gradient[j] * fac;
        for (unsigned int j = 0; j < nsd_ * nsd_; ++j)
          localMat(i, nsd_ + 1 + j) += shapes_->shfunct(i, q) * velgrad[j] * fac;
      }
    }

    Epetra_SerialDenseMatrix massMat(shapes_->ndofs_, shapes_->ndofs_);
    massMat.Multiply('N', 'T', 1., massPart, massPart, 0.);
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(massMat);
      inverseMass.SetVectors(localMat, localMat);
      inverseMass.Solve();
    }

    for (unsigned int r = 0; r < shapes_->ndofs_; ++r)
    {
      solele->eleinteriorPressnp_(r) += localMat(r, nsd_);  // pressure
      for (unsigned int i = 0; i < nsd_; ++i)
        solele->eleinteriorVelnp_(i * shapes_->ndofs_ + r) += localMat(r, i);  // velocity
      for (unsigned int j = 0; j < nsd_ * nsd_; ++j)
        solele->eleinteriorGradVelnp_(j * shapes_->ndofs_ + r) +=
            localMat(r, j + nsd_ + 1);  // velocity gradient
    }
  }  // if (elevec2.M() > 0)

  // if(dyna_==INPAR::ACOU::acou_trapezoidal) // then we have to set the initial trace field,
  // because the calculation of the residual in the first time step needs a correct trace field!
  {
    // trace variable
    int nfdofs = 0;
    for (unsigned int face = 0; face < nfaces_; ++face)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Faces()[face]->Degree(),
          shapes_->usescompletepoly_, 2 * ele->Faces()[face]->Degree());
      shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      shapesface_->EvaluateFace(*ele, face);

      Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
      Epetra_SerialDenseMatrix trVec(shapesface_->nfdofs_, nsd_);

      mass.Scale(0.);
      trVec.Scale(0.);

      for (unsigned int q = 0; q < shapesface_->nfdofs_; ++q)
      {
        const double fac = shapesface_->jfac(q);
        double xyz[nsd_];
        for (unsigned int d = 0; d < nsd_; ++d) xyz[d] = shapesface_->xyzreal(d, q);
        double p;
        double gradient[nsd_];
        double velgrad[nsd_ * nsd_];

        EvaluateAll(*start_func, xyz, p, gradient, velgrad, 0.0);  // u and p at quadrature point

        // now fill the components in the mass matrix and the right hand side
        for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
        {
          // mass matrix
          for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
            mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;
          for (unsigned int d = 0; d < nsd_; ++d)
            trVec(i, d) += shapesface_->shfunct(i, q) * gradient[d] * fac;
        }
      }

      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(mass);
      inverseMass.SetVectors(trVec, trVec);
      inverseMass.Solve();

      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
        for (unsigned int d = 0; d < nsd_; ++d)
          elevec1(nfdofs + i + d * shapesface_->nfdofs_) = trVec(i, d);

      nfdofs += shapesface_->nfdofs_ * nsd_;
    }
  }

  return 0;
}  // ProjectField

/*----------------------------------------------------------------------*
 * ProjectDirichField
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouSolEleCalc<distype>::ProjectDirichField(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1)
{
  // trace variable (zero, because no memory, no time derivative)
  dsassert(
      elevec1.M() == int(nfaces_ * shapesface_->nfdofs_ * nsd_), "Wrong size in project vector 1");
  elevec1.Scale(0.0);

  // this is it:
  if (params.isParameter("faceconsider"))
  {
    Teuchos::Array<int>* functno = params.getPtr<Teuchos::Array<int>>("funct");
    const unsigned int* faceConsider = params.getPtr<unsigned int>("faceconsider");
    double* time = params.getPtr<double>("time");

    DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Faces()[*faceConsider]->Degree(),
        shapes_->usescompletepoly_, 2 * ele->Faces()[*faceConsider]->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, *faceConsider);

    Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
    Epetra_SerialDenseMatrix trVec(shapesface_->nfdofs_, nsd_ + 1);

    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      const double fac = shapesface_->jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d) xyz[d] = shapesface_->xyzreal(d, q);
      double p;
      double gradient[nsd_];
      double velgrad[nsd_ * nsd_];

      EvaluateAll((*functno)[0], xyz, p, gradient, velgrad, *time);

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
          mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;
        trVec(i, nsd_) += shapesface_->shfunct(i, q) * p * fac;
        for (unsigned int j = 0; j < nsd_; ++j)
          trVec(i, j) += shapes_->shfunct(i, q) * gradient[j] * fac;
      }
    }

    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec, trVec);
    inverseMass.Solve();

    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      for (unsigned int j = 0; j < nsd_; ++j) elevec1(i + j * shapesface_->nfdofs_) = trVec(i, j);
  }

  return 0;
}  // ProjectDirichField

/*----------------------------------------------------------------------*
 * EvaluateAll
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::EvaluateAll(const int start_func,
    const double (&xyz)[nsd_], double& p, double (&v)[nsd_], double (&gv)[nsd_ * nsd_],
    double t) const

{
  p = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(0, xyz, 0.0);
  int numcomp = DRT::Problem::Instance()->Funct(start_func - 1).NumberComponents();
  if (numcomp != 1 && numcomp != nsd_ + 1 && numcomp != nsd_ * nsd_ + nsd_ + 1)
    dserror("supply 1, DIM+1, or DIM*DIM+DIM+1 components (pressure, velocity, velocity gradient)");

  if (numcomp == nsd_ + 1 || numcomp == nsd_ * nsd_ + nsd_ + 1)
  {
    for (unsigned int d = 0; d < nsd_; ++d)
      v[d] = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(d + 1, xyz, t);
  }
  else
  {
    for (unsigned int d = 0; d < nsd_; ++d) v[d] = 0.0;
  }
  if (numcomp == nsd_ * nsd_ + nsd_ + 1)
  {
    for (unsigned int d = 0; d < nsd_ * nsd_; ++d)
      gv[d] = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(nsd_ + 1 + d, xyz, t);
  }
  else
    for (unsigned int d = 0; d < nsd_ * nsd_; ++d) gv[d] = 0.0;

  return;
}

/*----------------------------------------------------------------------*
 * ProjectOpticalField
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouSolEleCalc<distype>::ProjectOpticalField(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec2)
{
  DRT::ELEMENTS::AcouSol* solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);
  shapes_->Evaluate(*ele);

  double* energy = params.get<double*>("gpvalues");

  Epetra_SerialDenseVector localMat(shapes_->ndofs_);
  Epetra_SerialDenseMatrix massPart(shapes_->ndofs_, shapes_->nqpoints_);

  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    const double fac = shapes_->jfac(q);
    const double sqrtfac = std::sqrt(fac);

    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
    {
      massPart(i, q) = shapes_->shfunct(i, q) * sqrtfac;
      localMat(i) -=
          shapes_->shfunct(i, q) * energy[q] *
          fac;  // negative sign for convention p_0=-Gamma mu_a phi but energy is only "mu_a phi"
    }
  }

  Epetra_SerialDenseMatrix massMat(shapes_->ndofs_, shapes_->ndofs_);
  massMat.Multiply('N', 'T', 1., massPart, massPart, 0.);
  {
    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(massMat);
    inverseMass.SetVectors(localMat, localMat);
    inverseMass.Solve();
  }

  for (unsigned int r = 0; r < shapes_->ndofs_; ++r)
    solele->eleinteriorPressnp_(r) += localMat(r, 0);  // pressure
  solele->eleinteriorVelnp_.Scale(0.0);                // velocity
  solele->eleinteriorGradVelnp_.Scale(0.0);            // velocity gradient

  return 0;
}  // ProjectOpticalField

/*----------------------------------------------------------------------*
 * ReadGlobalVectors
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::ReadGlobalVectors(
    DRT::Element* ele, DRT::Discretization& discretization, const std::vector<int>& lm)
{
  DRT::ELEMENTS::AcouSol* solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);

  reshapeMatrixIfNecessary(interiorGradVelnp_, shapes_->ndofs_ * nsd_ * nsd_, 1);
  reshapeMatrixIfNecessary(interiorVelnp_, shapes_->ndofs_ * nsd_, 1);
  reshapeMatrixIfNecessary(interiorPressnp_, shapes_->ndofs_, 1);

  interiorGradVelnp_ = solele->eleinteriorGradVelnp_;
  interiorVelnp_ = solele->eleinteriorVelnp_;
  interiorPressnp_ = solele->eleinteriorPressnp_;

  // read the HDG solution vector (for traces)
  if (discretization.HasState("trace"))
  {
    traceVal_.resize(nfaces_ * shapesface_->nfdofs_ * nsd_);
    dsassert(lm.size() == traceVal_.size(), "Internal error");
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState("trace");
    DRT::UTILS::ExtractMyValues(*matrix_state, traceVal_, lm);
  }

  return;
}  // ReadGlobalVectors

/*----------------------------------------------------------------------*
 * FillRestartVectors
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::FillRestartVectors(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  std::vector<double> interiorValnp(shapes_->ndofs_ * (nsd_ * nsd_ + nsd_ + 1));

  for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
  {
    interiorValnp[i] = interiorPressnp_(i);
    for (unsigned int d = 0; d < nsd_; ++d)
      interiorValnp[shapes_->ndofs_ + i + d * shapes_->ndofs_] =
          interiorVelnp_(i + d * shapes_->ndofs_);
    for (unsigned int d = 0; d < nsd_ * nsd_; ++d)
      interiorValnp[shapes_->ndofs_ * (1 + nsd_) + i + d * shapes_->ndofs_] =
          interiorGradVelnp_(i + d * shapes_->ndofs_);
  }

  std::vector<int> localDofs = discretization.Dof(1, ele);
  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1, "intvelnp");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for (unsigned int i = 0; i < localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  return;
}  // FillRestartVectors


/*----------------------------------------------------------------------*
 * ElementInitFromRestart
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::ElementInitFromRestart(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  DRT::ELEMENTS::AcouSol* solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);
  std::vector<double> interiorValnp(shapes_->ndofs_ * (nsd_ * nsd_ + nsd_ + 2));
  std::vector<int> localDofs1 = discretization.Dof(1, solele);

  Teuchos::RCP<const Epetra_Vector> intvel = discretization.GetState(1, "intvelnp");
  DRT::UTILS::ExtractMyValues(*intvel, interiorValnp, localDofs1);
  for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
  {
    solele->eleinteriorPressnp_(i) = interiorValnp[i];
    for (unsigned int d = 0; d < nsd_; ++d)
      solele->eleinteriorVelnp_(i + d * shapes_->ndofs_) =
          interiorValnp[shapes_->ndofs_ + i + d * shapes_->ndofs_];
    for (unsigned int d = 0; d < nsd_ * nsd_; ++d)
      solele->eleinteriorGradVelnp_(i + d * shapes_->ndofs_) =
          interiorValnp[shapes_->ndofs_ * (1 + nsd_) + i + d * shapes_->ndofs_];
  }

  return;
}  // ElementInitFromRestart


/*----------------------------------------------------------------------*
 * ComputeError
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::ComputeError(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec)
{
  // for the calculation of the error, we use a higher integration rule
  Teuchos::RCP<DRT::UTILS::GaussPoints> highquad =
      DRT::UTILS::GaussPointCache::Instance().Create(distype, (ele->Degree() + 2) * 2);
  LINALG::Matrix<nsd_, 1> xsi;
  Epetra_SerialDenseVector values(shapes_->ndofs_);
  LINALG::Matrix<nsd_, nen_> deriv;
  LINALG::Matrix<nsd_, nsd_> xjm;

  double time = params.get<double>("time");

  // get function
  int funcno = params.get<int>("funct");
  if (funcno < 0)
    dserror(
        "please provide an analytic solution for the error calculation, set CALCERRORFUNCNO and "
        "the corresponding FUNCT");
  funcno--;
  if (DRT::Problem::Instance()->Funct(funcno).NumberComponents() !=
      int(nsd_ * nsd_) + int(nsd_) + 1)
    dserror(
        "please provide numdim*numdim+numdim+1 components in the function for the analytic "
        "solution, first for pressure, following for velocity, following for time derivative of "
        "velocity gradient");

  double err_p = 0.0, norm_p = 0.0;
  double numerical_p = 0.0;
  double exact_p = 0.0;

  double err_v = 0.0, norm_v = 0.0;
  double numerical_v[nsd_] = {0.0};
  double exact_v[nsd_] = {0.0};

  double err_h = 0.0, norm_h = 0.0;
  double numerical_h[nsd_ * nsd_] = {0.0};
  double exact_h[nsd_ * nsd_] = {0.0};


  for (int q = 0; q < highquad->NumPoints(); ++q)
  {
    // initialization
    numerical_p = 0.0;
    exact_p = 0.0;
    const double* gpcoord = highquad->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++)
    {
      xsi(idim) = gpcoord[idim];
      numerical_v[idim] = 0.0;
      exact_v[idim] = 0.0;
    }
    for (unsigned int idim = 0; idim < nsd_ * nsd_; idim++)
    {
      numerical_h[idim] = 0.0;
      exact_h[idim] = 0.0;
    }
    shapes_->polySpace_->Evaluate(xsi, values);

    DRT::UTILS::shape_function_deriv1<distype>(xsi, deriv);
    xjm.MultiplyNT(deriv, shapes_->xyze);
    double highjfac = xjm.Determinant() * highquad->Weight(q);

    // evaluation of numerical values
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
    {
      numerical_p += values(i) * interiorPressnp_(i);
      for (unsigned int idim = 0; idim < nsd_; ++idim)
        numerical_v[idim] += values(i) * interiorVelnp_(i + idim * shapes_->ndofs_);
      for (unsigned int idim = 0; idim < nsd_ * nsd_; idim++)
        numerical_h[idim] += values(i) * interiorGradVelnp_(i + idim * shapes_->ndofs_);
    }

    // evaluation of analytical values
    LINALG::Matrix<nen_, 1> myfunct;
    DRT::UTILS::shape_function<distype>(xsi, myfunct);
    LINALG::Matrix<nsd_, 1> xyzmat;
    xyzmat.MultiplyNN(shapes_->xyze, myfunct);

    double xyz[nsd_];
    for (unsigned int d = 0; d < nsd_; ++d) xyz[d] = xyzmat(d, 0);

    exact_p = DRT::Problem::Instance()->Funct(funcno).Evaluate(0, xyz, time);
    for (unsigned int idim = 0; idim < nsd_; ++idim)
      exact_v[idim] = DRT::Problem::Instance()->Funct(funcno).Evaluate(int(idim) + 1, xyz, time);
    for (unsigned int idim = 0; idim < nsd_ * nsd_; idim++)
      exact_h[idim] =
          DRT::Problem::Instance()->Funct(funcno).Evaluate(int(idim) + 1 + nsd_, xyz, time);

    // error calculation
    err_p += (exact_p - numerical_p) * (exact_p - numerical_p) * highjfac;
    norm_p += exact_p * exact_p * highjfac;

    for (unsigned int idim = 0; idim < nsd_; ++idim)
    {
      err_v += (exact_v[idim] - numerical_v[idim]) * (exact_v[idim] - numerical_v[idim]) * highjfac;
      norm_v += exact_v[idim] * exact_v[idim] * highjfac;
    }
    for (unsigned int idim = 0; idim < nsd_ * nsd_; ++idim)
    {
      err_h += (exact_h[idim] - numerical_h[idim]) * (exact_h[idim] - numerical_h[idim]) * highjfac;
      norm_h += exact_h[idim] * exact_h[idim] * highjfac;
    }
  }

  elevec[0] += err_p;
  elevec[1] += norm_p;
  elevec[2] += err_v;
  elevec[3] += norm_v;
  elevec[4] += err_h;
  elevec[5] += norm_h;

  return;
}  // ComputeError

/*----------------------------------------------------------------------*
 * NodeBasedValues
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::NodeBasedValues(const Teuchos::RCP<MAT::Material>& mat,
    DRT::ELEMENTS::Acou* ele, Epetra_SerialDenseVector& elevec1, bool writestress)
{
  dsassert(elevec1.M() == (int)nen_ * (2 * nsd_ + 2 + nsd_ * nsd_) + 2,
      "Vector does not have correct size");
  elevec1.Scale(0.0);

  const MAT::AcousticSolMat* actmat = static_cast<const MAT::AcousticSolMat*>(mat.get());
  double c = actmat->SpeedofSound();
  double visc = actmat->Viscosity();
  double lambda = c * c - visc;

  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  double norm_p = 0.0, area = 0.0;
  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    double numerical_p = 0.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
    {
      numerical_p += shapes_->shfunct(i, q) * interiorPressnp_(i);
    }
    norm_p += numerical_p * shapes_->jfac(q);
    area += shapes_->jfac(q);
  }

  elevec1((2 * nsd_ + 2) * nen_) = norm_p / area;

  for (unsigned int i = 0; i < nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim = 0; idim < nsd_; idim++) shapes_->xsi(idim) = locations(idim, i);
    shapes_->polySpace_->Evaluate(shapes_->xsi, values);

    // compute values for velocity and pressure by summing over all basis functions

    // pressure
    double sump = 0.0;
    for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
    {
      sump += values(k) * interiorPressnp_(k);
    }
    elevec1(2 * nsd_ * nen_ + i) = sump;

    // velocity
    for (unsigned int d = 0; d < nsd_; ++d)
    {
      double sumv = 0.0;
      for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
        sumv += values(k) * interiorVelnp_(d * shapes_->ndofs_ + k);
      elevec1(d * nen_ + i) = sumv;
    }

    // stress or velocitygradient
    if (!writestress)
    {
      for (unsigned int d = 0; d < nsd_; ++d)
        for (unsigned int e = 0; e < nsd_; ++e)
        {
          double sumstress = 0.0;
          for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
            sumstress += values(k) *
                         interiorGradVelnp_(k + d * nsd_ * shapes_->ndofs_ + e * shapes_->ndofs_);
          elevec1(nen_ * (2 * nsd_ + 2) + 2 + nen_ * (d * nsd_ + e) + i) = sumstress;
        }
    }
    else
    {
      // visualization of the stresses sigma = mu ( H + H^T ) - lamda / ( lambda + mu ) p I
      for (unsigned int d = 0; d < nsd_; ++d)
        for (unsigned int e = 0; e <= d; ++e)
        {
          double sumstress = 0.0;
          if (d != e)
          {
            for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
              sumstress +=
                  values(k) * visc *
                  (interiorGradVelnp_(k + d * nsd_ * shapes_->ndofs_ + e * shapes_->ndofs_) +
                      interiorGradVelnp_(k + e * nsd_ * shapes_->ndofs_ + d * shapes_->ndofs_));
            if (d == 1 && e == 0)
              elevec1(nen_ * (2 * nsd_ + 2) + 2 + nen_ * (3) + i) = sumstress;
            else if (d == 2 && e == 0)
              elevec1(nen_ * (2 * nsd_ + 2) + 2 + nen_ * (5) + i) = sumstress;
            else  // (if (d==2&&e==1))
              elevec1(nen_ * (2 * nsd_ + 2) + 2 + nen_ * (4) + i) = sumstress;
          }
          else  // (if(d==e))
          {
            for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
              sumstress +=
                  values(k) *
                  (visc * interiorGradVelnp_(k + d * nsd_ * shapes_->ndofs_ + e * shapes_->ndofs_) +
                      visc *
                          interiorGradVelnp_(k + e * nsd_ * shapes_->ndofs_ + d * shapes_->ndofs_) -
                      lambda / (lambda + visc) * interiorPressnp_(k));
            elevec1(nen_ * (2 * nsd_ + 2) + 2 + nen_ * (d) + i) = sumstress;
          }
        }
    }
  }

  // trace velocity
  locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(
      DRT::UTILS::DisTypeToFaceShapeType<distype>::shape);
  Epetra_SerialDenseVector fvalues(shapesface_->nfdofs_);
  Epetra_SerialDenseVector temp(elevec1.M());
  for (unsigned int face = 0; face < nfaces_; ++face)
  {
    shapesface_->EvaluateFace(*ele, face);
    for (int i = 0; i < DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
    {
      // evaluate shape polynomials in node
      for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
        shapesface_->xsi(idim) = locations(idim, i);
      shapesface_->polySpace_->Evaluate(
          shapesface_->xsi, fvalues);  // TODO: fix face orientation here

      // compute values for velocity and pressure by summing over all basis functions
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        double sum = 0.0;
        for (unsigned int k = 0; k < shapesface_->nfdofs_; ++k)
          sum += fvalues(k) *
                 traceVal_[face * shapesface_->nfdofs_ * nsd_ + d * shapesface_->nfdofs_ + k];
        elevec1((nsd_ + d) * nen_ + shapesface_->faceNodeOrder[face][i]) += sum;
        temp((nsd_ + d) * nen_ + shapesface_->faceNodeOrder[face][i])++;
      }
    }
  }
  for (unsigned int i = 0; i < nen_ * nsd_; ++i) elevec1(nsd_ * nen_ + i) /= temp(nsd_ * nen_ + i);

  return;
}  // NodeBasedValues

/*----------------------------------------------------------------------*
 * NodeBasedPsi
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::NodeBasedPsi(const Teuchos::RCP<MAT::Material>& mat,
    DRT::ELEMENTS::Acou* ele, Epetra_SerialDenseVector& elevec1, double dt)
{
  dsassert(elevec1.M() == int(nen_), "Vector does not have correct size");
  elevec1.Scale(0.0);
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);
  const MAT::AcousticSolMat* actmat = static_cast<const MAT::AcousticSolMat*>(mat.get());
  double c = actmat->SpeedofSound();

  for (unsigned int i = 0; i < nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim = 0; idim < nsd_; idim++) shapes_->xsi(idim) = locations(idim, i);
    shapes_->polySpace_->Evaluate(shapes_->xsi, values);

    // compute values for velocity and pressure by summing over all basis functions
    double sump = 0.0;
    for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
    {
      sump += values(k) * interiorPressnp_(k);
    }
    elevec1(i) = sump / dt / c / c;
  }

  return;
}  // NodeBasedPsi

/*----------------------------------------------------------------------*
 * Compute internal and face matrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::ComputeMatrices(const Teuchos::RCP<MAT::Material>& mat,
    DRT::ELEMENTS::Acou& ele, double dt, INPAR::ACOU::DynamicType dyna)
{
  const MAT::AcousticSolMat* actmat = static_cast<const MAT::AcousticSolMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();
  double tau = 1.0;  // or rho * c (because of the units)
  double visc = actmat->Viscosity();

  zeroMatrix(localSolver_->amat);
  zeroMatrix(localSolver_->invamat);
  zeroMatrix(localSolver_->bmat);
  zeroMatrix(localSolver_->cmat);
  zeroMatrix(localSolver_->dmat);
  zeroMatrix(localSolver_->emat);
  zeroMatrix(localSolver_->ehatmat);
  zeroMatrix(localSolver_->fmat);
  zeroMatrix(localSolver_->gmat);
  zeroMatrix(localSolver_->hmat);
  zeroMatrix(localSolver_->imat);
  zeroMatrix(localSolver_->jmat);
  zeroMatrix(localSolver_->kmat);
  zeroMatrix(localSolver_->lmat);
  zeroMatrix(localSolver_->mmat);
  zeroMatrix(localSolver_->nmat);

  localSolver_->ComputeInteriorMatrices();
  for (unsigned int face = 0; face < nfaces_; ++face)
  {
    shapesface_->EvaluateFace(ele, face);
    localSolver_->ComputeFaceMatrices(face);
  }

  // scale the matrices
  localSolver_->amat.Scale(1.0 / dt);
  localSolver_->invamat.Scale(dt);
  localSolver_->emat.Scale(rho / dt);
  localSolver_->ehatmat.Scale(tau);
  localSolver_->imat.Scale(1.0 / c / c / dt);
  localSolver_->nmat.Scale(tau);

  // localSolver_->bmat.Scale(          );
  localSolver_->cmat.Scale(-1.0);
  localSolver_->dmat.Scale(-visc);
  // localSolver_->fmat.Scale(          );
  localSolver_->gmat.Scale(-tau);
  localSolver_->hmat.Scale(-rho);
  localSolver_->jmat.Scale(rho);
  localSolver_->kmat.Scale(visc);
  localSolver_->lmat.Scale(-tau);
  localSolver_->mmat.Scale(-1.0);


  return;
}  // ComputeMatrices

/*----------------------------------------------------------------------*
 * ComputeInteriorMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::ComputeInteriorMatrices()
{
  // interior matrices are a, b, d, e, ehat, f, h, i, k, l and lhat

  // standard mass matrices are for example i, k, l
  // we calculate i first, and then fill all related matrices

  Epetra_SerialDenseMatrix massPart(ndofs_, shapes_.nqpoints_);
  Epetra_SerialDenseMatrix gradPart(ndofs_ * nsd_, shapes_.nqpoints_);

  // loop quadrature points
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    const double sqrtfac = std::sqrt(shapes_.jfac(q));
    // loop shape functions
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      const double valf = shapes_.shfunct(i, q) * sqrtfac;
      massPart(i, q) = valf;
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        const double vald = shapes_.shderxy(i * nsd_ + d, q) * sqrtfac;
        gradPart(d * ndofs_ + i, q) = vald;
      }
    }
  }

  // multiply matrices to perform summation over quadrature points
  for (unsigned int i = 0; i < ndofs_; ++i)
  {
    for (unsigned int j = 0; j < ndofs_; ++j)
    {
      double sum = 0.0;
      double sums[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d) sums[d] = 0.0;
      for (unsigned int k = 0; k < shapes_.nqpoints_; ++k)
      {
        sum += massPart(i, k) * massPart(j, k);
        for (unsigned int d = 0; d < nsd_; ++d)
          sums[d] += gradPart(d * ndofs_ + i, k) * massPart(j, k);
      }
      imat(i, j) = sum;
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        hmat(i, j + d * ndofs_) = sums[d];
        fmat(j + d * ndofs_, i) = sums[d];
      }
    }
  }

  Epetra_SerialDenseMatrix tempmassgrad(nsd_ * ndofs_, ndofs_);
  tempmassgrad.Multiply('N', 'T', 1., gradPart, massPart, 0.);

  bmat.Scale(0.0);
  dmat.Scale(0.0);
  for (unsigned d = 0; d < nsd_; ++d)
    for (unsigned i = 0; i < ndofs_; ++i)
      for (unsigned j = 0; j < ndofs_ * nsd_; ++j)
      {
        bmat(j + d * nsd_ * ndofs_, i + d * ndofs_) = tempmassgrad(j, i);
        dmat(i + d * ndofs_, j + d * nsd_ * ndofs_) = tempmassgrad(j, i);
      }

  // e is also a mass matrix, but for every dimension. here, we have to consider the ordering
  // of the velocity components!
  for (unsigned int i = 0; i < ndofs_; ++i)
    for (unsigned int j = 0; j < ndofs_; ++j)
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        emat(d * ndofs_ + i, d * ndofs_ + j) = imat(i, j);
        for (unsigned int e = 0; e < nsd_; ++e)
          amat(i + (d * nsd_ + e) * ndofs_, j + (d * nsd_ + e) * ndofs_) = imat(i, j);
      }

  // we will need the inverse of amat. but we invert imat, since amat is a blockmutiple of imat
  invamat = imat;
  {
    Epetra_SerialDenseSolver inverseamat;
    inverseamat.SetMatrix(invamat);
    inverseamat.Invert();
  }

  invamat.Reshape(nsd_ * nsd_ * ndofs_, nsd_ * nsd_ * ndofs_);
  for (unsigned int d = 1; d < nsd_ * nsd_; ++d)
    for (unsigned int i = 0; i < ndofs_; ++i)
      for (unsigned int j = 0; j < ndofs_; ++j)
        invamat(d * ndofs_ + i, d * ndofs_ + j) = invamat(i, j);

  return;
}  // ComputeInteriorMatrices

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::ComputeFaceMatrices(const int face)
{
  // missing term in ehat
  for (unsigned int p = 0; p < ndofs_; ++p)
  {
    for (unsigned int q = 0; q < ndofs_; ++q)
    {
      double tempehat = 0.0;
      for (unsigned int i = 0; i < shapesface_.nqpoints_; ++i)
        tempehat += shapesface_.jfac(i) * shapesface_.shfunctI(p, i) * shapesface_.shfunctI(q, i);
      for (unsigned int d = 0; d < nsd_; ++d) ehatmat(d * ndofs_ + q, d * ndofs_ + p) += tempehat;
    }
  }

  // n
  for (unsigned int p = 0; p < shapesface_.nfdofs_; ++p)
  {
    for (unsigned int q = 0; q <= p; ++q)
    {
      double tempn = 0.0;
      for (unsigned int i = 0; i < shapesface_.nqpoints_; ++i)
        tempn += shapesface_.jfac(i) * shapesface_.shfunct(p, i) * shapesface_.shfunct(q, i);
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        nmat(face * shapesface_.nfdofs_ * nsd_ + p + d * shapesface_.nfdofs_,
            face * shapesface_.nfdofs_ * nsd_ + q + d * shapesface_.nfdofs_) = tempn;
        nmat(face * shapesface_.nfdofs_ * nsd_ + q + d * shapesface_.nfdofs_,
            face * shapesface_.nfdofs_ * nsd_ + p + d * shapesface_.nfdofs_) = tempn;
      }
    }
  }

  // c, g, j, m, n, o
  for (unsigned int p = 0; p < shapesface_.nfdofs_; ++p)
  {
    for (unsigned int q = 0; q < ndofs_; ++q)
    {
      double tempmat = 0.0;
      for (unsigned int i = 0; i < shapesface_.nqpoints_; ++i)
      {
        double temp = shapesface_.jfac(i) * shapesface_.shfunct(p, i) * shapesface_.shfunctI(q, i);
        tempmat += temp;
        for (unsigned int j = 0; j < nsd_; ++j)
        {
          double temp_d = temp * shapesface_.normals(j, i);
          jmat(q, p + j * shapesface_.nfdofs_ + face * shapesface_.nfdofs_ * nsd_) += temp_d;
          mmat(p + j * shapesface_.nfdofs_ + face * shapesface_.nfdofs_ * nsd_, q) += temp_d;

          for (unsigned int e = 0; e < nsd_; ++e)
          {
            kmat(p + e * shapesface_.nfdofs_ + face * shapesface_.nfdofs_ * nsd_,
                q + (e * nsd_ + j) * ndofs_) += temp_d;
            cmat(q + (e * nsd_ + j) * ndofs_,
                p + e * shapesface_.nfdofs_ + face * shapesface_.nfdofs_ * nsd_) += temp_d;
          }
        }
      }

      for (unsigned int j = 0; j < nsd_; ++j)
      {
        gmat(q + j * ndofs_, p + j * shapesface_.nfdofs_ + face * shapesface_.nfdofs_ * nsd_) =
            tempmat;
        lmat(p + j * shapesface_.nfdofs_ + face * shapesface_.nfdofs_ * nsd_, q + j * ndofs_) =
            tempmat;
      }
    }
  }

  return;
}  // ComputeFaceMatrices

/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::CondenseLocalPart(
    Epetra_SerialDenseMatrix& eleMat, INPAR::ACOU::DynamicType dyna)
{
  Epetra_SerialDenseMatrix dinvamat(nsd_ * ndofs_, nsd_ * nsd_ * ndofs_);
  dinvamat.Multiply('N', 'N', 1.0, dmat, invamat, 0.0);

  Epetra_SerialDenseMatrix ol(nsd_ * ndofs_, nsd_ * ndofs_);
  ol = ehatmat;
  ol += emat;
  ol.Multiply('N', 'N', -1.0, dinvamat, bmat, 1.0);

  Epetra_SerialDenseMatrix toinv((nsd_ + 1) * ndofs_, (nsd_ + 1) * ndofs_);
  for (unsigned int i = 0; i < ndofs_; ++i)
  {
    for (unsigned int j = 0; j < ndofs_; ++j)
    {
      toinv(nsd_ * ndofs_ + i, nsd_ * ndofs_ + j) = imat(i, j);
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        toinv(d * ndofs_ + i, d * ndofs_ + j) = ol(d * ndofs_ + i, d * ndofs_ + j);
        toinv(nsd_ * ndofs_ + i, d * ndofs_ + j) = hmat(i, d * ndofs_ + j);
        toinv(d * ndofs_ + i, nsd_ * ndofs_ + j) = fmat(d * ndofs_ + i, j);
      }
    }
  }
  Epetra_SerialDenseMatrix rhsinv((nsd_ + 1) * ndofs_, nsd_ * shapesface_.nfdofs_ * nfaces_);

  Epetra_SerialDenseMatrix tempmat1(nsd_ * ndofs_, nsd_ * shapesface_.nfdofs_ * nfaces_);
  tempmat1 = gmat;
  tempmat1.Multiply('N', 'N', -1.0, dinvamat, cmat, 1.0);

  for (unsigned int i = 0; i < nsd_ * shapesface_.nfdofs_ * nfaces_; ++i)
    for (unsigned int j = 0; j < nsd_ * ndofs_; ++j) rhsinv(j, i) = tempmat1(j, i);

  for (unsigned int i = 0; i < nsd_ * shapesface_.nfdofs_ * nfaces_; ++i)
    for (unsigned int j = 0; j < ndofs_; ++j) rhsinv(nsd_ * ndofs_ + j, i) = jmat(j, i);

  // invert
  {
    Epetra_SerialDenseSolver inverse;
    inverse.SetMatrix(toinv);
    inverse.SetVectors(rhsinv, rhsinv);
    inverse.Solve();
  }

  tempmat1.Shape(ndofs_, nsd_ * shapesface_.nfdofs_ * nfaces_);
  for (unsigned int i = 0; i < nsd_ * shapesface_.nfdofs_ * nfaces_; ++i)
    for (unsigned int j = 0; j < ndofs_; ++j) tempmat1(j, i) = rhsinv(nsd_ * ndofs_ + j, i);

  eleMat = nmat;
  eleMat.Multiply('N', 'N', -1.0, mmat, tempmat1, 1.0);

  tempmat1.Shape(nsd_ * ndofs_, nsd_ * shapesface_.nfdofs_ * nfaces_);
  for (unsigned int i = 0; i < nsd_ * shapesface_.nfdofs_ * nfaces_; ++i)
    for (unsigned int j = 0; j < nsd_ * ndofs_; ++j) tempmat1(j, i) = rhsinv(j, i);

  eleMat.Multiply('N', 'N', -1.0, lmat, tempmat1, 1.0);

  ol.Shape(nsd_ * nsd_ * ndofs_, nsd_ * shapesface_.nfdofs_ * nfaces_);
  ol = cmat;
  ol.Multiply('N', 'N', -1.0, bmat, tempmat1, 1.0);
  tempmat1.Shape(nsd_ * nsd_ * ndofs_, nsd_ * shapesface_.nfdofs_ * nfaces_);
  tempmat1.Multiply('N', 'N', 1.0, invamat, ol, 0.0);
  eleMat.Multiply('N', 'N', -1.0, kmat, tempmat1, 1.0);

  return;
}  // CondenseLocalPart

/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::ComputeResidual(
    Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec,
    Epetra_SerialDenseVector& interiorGradVeln, Epetra_SerialDenseVector& interiorVeln,
    Epetra_SerialDenseVector& interiorPressn, std::vector<double> traceVal,
    INPAR::ACOU::DynamicType dyna)
{
  Epetra_SerialDenseVector traceVal_SDV(nfaces_ * shapesface_.nfdofs_ * nsd_);
  for (unsigned i = 0; i < nfaces_ * shapesface_.nfdofs_ * nsd_; ++i) traceVal_SDV(i) = traceVal[i];

  Epetra_SerialDenseMatrix dinvamat(nsd_ * ndofs_, nsd_ * nsd_ * ndofs_);
  dinvamat.Multiply('N', 'N', 1.0, dmat, invamat, 0.0);

  Epetra_SerialDenseMatrix toinv((nsd_ + 1) * ndofs_, (nsd_ + 1) * ndofs_);
  {
    Epetra_SerialDenseMatrix ol(nsd_ * ndofs_, nsd_ * ndofs_);
    ol = ehatmat;
    ol += emat;
    ol.Multiply('N', 'N', -1.0, dinvamat, bmat, 1.0);

    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int j = 0; j < ndofs_; ++j)
      {
        toinv(nsd_ * ndofs_ + i, nsd_ * ndofs_ + j) = imat(i, j);
        for (unsigned int d = 0; d < nsd_; ++d)
        {
          toinv(d * ndofs_ + i, d * ndofs_ + j) = ol(d * ndofs_ + i, d * ndofs_ + j);
          toinv(nsd_ * ndofs_ + i, d * ndofs_ + j) = hmat(i, d * ndofs_ + j);
          toinv(d * ndofs_ + i, nsd_ * ndofs_ + j) = fmat(d * ndofs_ + i, j);
        }
      }
    }
  }
  Epetra_SerialDenseVector rhsinv((nsd_ + 1) * ndofs_);
  Epetra_SerialDenseVector ol(nsd_ * ndofs_);

  // source term
  Epetra_SerialDenseVector dummy(nsd_ * ndofs_);
  ComputeSource(params, dummy, ol);
  ol.Multiply('N', 'N', 1.0, emat, interiorVeln, 1.0);

  Epetra_SerialDenseVector tempvec1(nsd_ * nsd_ * ndofs_);
  tempvec1.Multiply('N', 'N', 1.0, amat, interiorGradVeln, 0.0);
  ol.Multiply('N', 'N', -1.0, dinvamat, tempvec1, 1.0);
  for (unsigned int i = 0; i < nsd_ * ndofs_; ++i) rhsinv(i) = ol(i);

  ol.Resize(ndofs_);
  ol.Multiply('N', 'N', 1.0, imat, interiorPressn, 0.0);
  for (unsigned int i = 0; i < ndofs_; ++i) rhsinv(nsd_ * ndofs_ + i) = ol(i);

  // invert
  {
    Epetra_SerialDenseSolver inverse;
    inverse.SetMatrix(toinv);
    inverse.SetVectors(rhsinv, rhsinv);
    inverse.Solve();
  }

  for (unsigned int i = 0; i < ndofs_; ++i) ol(i) = rhsinv(nsd_ * ndofs_ + i);
  elevec.Multiply('N', 'N', -1.0, mmat, ol, 1.0);

  ol.Shape(nsd_ * ndofs_, 1);
  for (unsigned int i = 0; i < nsd_ * ndofs_; ++i) ol(i, 0) = rhsinv(i);
  elevec.Multiply('N', 'N', -1.0, lmat, ol, 1.0);

  tempvec1.Multiply('N', 'N', -1.0, bmat, ol, 1.0);
  ol.Shape(nsd_ * nsd_ * ndofs_, 1);
  ol.Multiply('N', 'N', 1.0, invamat, tempvec1, 0.0);
  elevec.Multiply('N', 'N', -1.0, kmat, ol, 1.0);

  return;
}  // ComputeResidual


/*----------------------------------------------------------------------*
 * ComputeAbsorbingBC
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::ComputeAbsorbingBC(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat,
    int face, Epetra_SerialDenseMatrix& elemat, Epetra_SerialDenseVector& elevec)
{
  const MAT::AcousticSolMat* actmat = static_cast<const MAT::AcousticSolMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();
  bool resonly = params.get<bool>("resonly");

  // warning! this is not the correct first order absorbing boundary condition
  // check the hdg_acou paper, section 3.5 and implement it!

  // symmetric
  if (!resonly)
  {
    for (unsigned int p = 0; p < shapesface_.nfdofs_; ++p)
    {
      for (unsigned int q = 0; q < shapesface_.nfdofs_; ++q)
      {
        double temp[nsd_ * nsd_];
        for (unsigned int i = 0; i < nsd_ * nsd_; ++i) temp[i] = 0.0;
        for (unsigned int i = 0; i < shapesface_.nqpoints_; ++i)
          for (unsigned int d = 0; d < nsd_; ++d)
            for (unsigned int e = 0; e < nsd_; ++e)
              temp[d + e * nsd_] += shapesface_.jfac(i) * shapesface_.shfunct(p, i) *
                                    shapesface_.shfunct(q, i) * shapesface_.normals(d, i) *
                                    shapesface_.normals(e, i);

        for (unsigned int d = 0; d < nsd_; ++d)
        {
          for (unsigned int e = 0; e < nsd_; ++e)
          {
            elemat(face * shapesface_.nfdofs_ * nsd_ + p + d * shapesface_.nfdofs_,
                face * shapesface_.nfdofs_ * nsd_ + q + e * shapesface_.nfdofs_) +=
                rho * c * temp[d + e * nsd_];
          }
        }
      }
    }
  }
  return;
}  // ComputeAbsorbingBC

/*----------------------------------------------------------------------*
 * ComputeSource
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::ComputeSource(
    Teuchos::ParameterList& params, Epetra_SerialDenseVector& interiorSourcen,
    Epetra_SerialDenseVector& interiorSourcenp)
{
  int funcno = params.get<int>("sourcefuncno");
  if (funcno < 0) return;  // there is no such thing as a volume force

  if (DRT::Problem::Instance()->Funct(funcno).NumberComponents() != int(nsd_))
    dserror(
        "for acoustic sol elements, the body load is on the velocity and hence has to have numdim "
        "components");

  // the vector to be filled
  Epetra_SerialDenseVector source(shapes_.ndofs_ * nsd_);

  // what time is it?
  double tn = params.get<double>("time");
  double tp = params.get<double>("timep");

  double f_value = 0.0;
  double f_value_p = 0.0;
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    double xyz[nsd_];
    for (unsigned int d = 0; d < nsd_; ++d) xyz[d] = shapes_.xyzreal(d, q);

    for (unsigned int dim = 0; dim < nsd_; ++dim)
    {
      // calculate right hand side contribution for dp/dt
      f_value = DRT::Problem::Instance()->Funct(funcno).Evaluate(int(dim), xyz, tn);
      f_value_p = DRT::Problem::Instance()->Funct(funcno).Evaluate(int(dim), xyz, tp);

      // add it all up
      for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      {
        interiorSourcen(i + dim * shapes_.ndofs_) +=
            shapes_.shfunct(i, q) * f_value * shapes_.jfac(q);
        interiorSourcenp(i + dim * shapes_.ndofs_) +=
            shapes_.shfunct(i, q) * f_value_p * shapes_.jfac(q);
      }
    }
  }

  return;
}  // ComputeSource


/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::UpdateInteriorVariablesAndComputeResidual(
    Teuchos::ParameterList& params, DRT::ELEMENTS::Acou* ele, Epetra_SerialDenseVector& elevec,
    double dt, bool errormaps, bool updateonly)
{
  DRT::ELEMENTS::AcouSol* solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);

  Epetra_SerialDenseVector tempGradVelnp;
  Epetra_SerialDenseVector tempVelnp;
  Epetra_SerialDenseVector tempPressnp;

  Epetra_SerialDenseVector traceVal_SDV(nfaces_ * shapesface_->nfdofs_ * nsd_);
  for (unsigned i = 0; i < nfaces_ * shapesface_->nfdofs_ * nsd_; ++i)
    traceVal_SDV(i) = traceVal_[i];

  // *****************************************************
  // update interior variables first
  // *****************************************************

  Epetra_SerialDenseMatrix toinv((nsd_ + 1) * shapes_->ndofs_, (nsd_ + 1) * shapes_->ndofs_);
  Epetra_SerialDenseVector tempsource(nsd_ * shapes_->ndofs_);

  Epetra_SerialDenseMatrix dinvamat(nsd_ * shapes_->ndofs_, nsd_ * nsd_ * shapes_->ndofs_);
  dinvamat.Multiply('N', 'N', 1.0, localSolver_->dmat, localSolver_->invamat, 0.0);

  {
    Epetra_SerialDenseMatrix ol(nsd_ * shapes_->ndofs_, nsd_ * shapes_->ndofs_);
    ol = localSolver_->ehatmat;
    ol += localSolver_->emat;
    ol.Multiply('N', 'N', -1.0, dinvamat, localSolver_->bmat, 1.0);

    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
    {
      for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
      {
        toinv(nsd_ * shapes_->ndofs_ + i, nsd_ * shapes_->ndofs_ + j) = localSolver_->imat(i, j);
        for (unsigned int d = 0; d < nsd_; ++d)
        {
          toinv(d * shapes_->ndofs_ + i, d * shapes_->ndofs_ + j) =
              ol(d * shapes_->ndofs_ + i, d * shapes_->ndofs_ + j);
          toinv(nsd_ * shapes_->ndofs_ + i, d * shapes_->ndofs_ + j) =
              localSolver_->hmat(i, d * shapes_->ndofs_ + j);
          toinv(d * shapes_->ndofs_ + i, nsd_ * shapes_->ndofs_ + j) =
              localSolver_->fmat(d * shapes_->ndofs_ + i, j);
        }
      }
    }
  }

  // calculate rhs parts
  Epetra_SerialDenseVector ol(nsd_ * shapes_->ndofs_);
  localSolver_->ComputeSource(params, ol, tempsource);
  ol.Multiply('N', 'N', 1.0, localSolver_->emat, interiorVelnp_, 1.0);
  ol.Multiply('N', 'N', -1.0, localSolver_->gmat, traceVal_SDV, 1.0);


  Epetra_SerialDenseVector tempvec1(nsd_ * nsd_ * shapes_->ndofs_);
  tempvec1.Multiply('N', 'N', 1.0, localSolver_->amat, interiorGradVelnp_, 0.0);
  tempvec1.Multiply('N', 'N', -1.0, localSolver_->cmat, traceVal_SDV, 1.0);
  ol.Multiply('N', 'N', -1.0, dinvamat, tempvec1, 1.0);
  Epetra_SerialDenseVector rhsinv((nsd_ + 1) * shapes_->ndofs_);
  for (unsigned int i = 0; i < nsd_ * shapes_->ndofs_; ++i) rhsinv(i) = ol(i);

  ol.Shape(shapes_->ndofs_, 1);
  ol.Multiply('N', 'N', 1.0, localSolver_->imat, interiorPressnp_, 0.0);
  ol.Multiply('N', 'N', -1.0, localSolver_->jmat, traceVal_SDV, 1.0);

  // invert
  {
    Epetra_SerialDenseSolver inverse;
    inverse.SetMatrix(toinv);
    inverse.Invert();
  }
  Epetra_SerialDenseVector sol((nsd_ + 1) * shapes_->ndofs_);
  sol.Multiply('N', 'N', 1.0, toinv, rhsinv, 0.0);

  for (unsigned int i = 0; i < nsd_ * shapes_->ndofs_; ++i) interiorVelnp_(i) = sol(i);
  for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
    interiorPressnp_(i) = sol(nsd_ * shapes_->ndofs_ + i);

  tempvec1.Multiply('N', 'N', -1.0, localSolver_->bmat, interiorVelnp_, 1.0);
  interiorGradVelnp_.Multiply('N', 'N', 1.0, localSolver_->invamat, tempvec1, 0.0);

  // tell this change in the interior variables the discretization
  solele->eleinteriorPressnp_ = interiorPressnp_;
  solele->eleinteriorVelnp_ = interiorVelnp_;
  solele->eleinteriorGradVelnp_ = interiorGradVelnp_;

  // *****************************************************
  // local postprocessing to calculate error maps
  // *****************************************************

  if (errormaps)
  {
    // first step: calculate the displacement gradient (bold L in paper by Nguyen)
    Epetra_SerialDenseVector temp(shapes_->ndofs_ * nsd_ * nsd_);
    temp.Multiply('N', 'N', 1.0, localSolver_->bmat, interiorVelnp_, 0.0);
    temp.Multiply('N', 'N', 1.0, localSolver_->cmat, traceVal_SDV, 1.0);
    Epetra_SerialDenseVector L(shapes_->ndofs_ * nsd_ * nsd_);
    L.Multiply('N', 'N', 1.0 / dt, localSolver_->invamat, temp, 0.0);

    // second step: postprocess the velocity field:
    double err_v = EstimateError(*solele, L);

    Teuchos::RCP<std::vector<double>> values =
        params.get<Teuchos::RCP<std::vector<double>>>("elevals");
    bool padaptivity = params.get<bool>("padaptivity");
    if (padaptivity)  // time integrators or last stage of DIRK
    {
      double padaptivitytol = params.get<double>("padaptivitytol");
      double delta_k = 0.0;
      if (err_v != 0.0) delta_k = std::ceil(std::log10(err_v / padaptivitytol));

      int newdeg = (solele->Degree() + delta_k);
      if (newdeg < 1) newdeg = 1;
      if (newdeg > 6) newdeg = 6;
      (*values)[solele->Id()] = newdeg;

      // projection of the internal field (cannot yet project trace field, since the
      // face degree depends on all neighbors! (also not needed for DIRK)
      ProjectPadapField(solele, newdeg);
      updateonly = true;
      return;
    }
    else
      (*values)[solele->Id()] = err_v;
  }

  if (updateonly) return;

  // *****************************************************
  // compute residual second (reuse intermediate matrices)
  // *****************************************************

  // calculate rhs parts
  ol.Shape(nsd_ * shapes_->ndofs_, 1);
  ol += tempsource;
  ol.Multiply('N', 'N', 1.0, localSolver_->emat, interiorVelnp_, 1.0);
  tempvec1.Scale(0.0);
  tempvec1.Multiply('N', 'N', 1.0, localSolver_->amat, interiorGradVelnp_, 0.0);
  ol.Multiply('N', 'N', -1.0, dinvamat, tempvec1, 1.0);

  for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i) rhsinv(i) = ol(i);

  ol.Resize(shapes_->ndofs_);
  ol.Multiply('N', 'N', 1.0, localSolver_->imat, interiorPressnp_, 0.0);

  sol.Multiply('N', 'N', 1.0, toinv, rhsinv, 0.0);

  for (unsigned int i = 0; i < shapes_->ndofs_; ++i) ol(i) = sol(nsd_ * shapes_->ndofs_ + i);
  elevec.Multiply('N', 'N', -1.0, localSolver_->mmat, ol, 1.0);

  ol.Shape(shapes_->ndofs_ * nsd_, 1);
  for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i) ol(i) = sol(i);
  elevec.Multiply('N', 'N', -1.0, localSolver_->lmat, ol, 1.0);

  tempvec1.Multiply('N', 'N', -1.0, localSolver_->bmat, ol, 1.0);
  ol.Shape(nsd_ * nsd_ * shapes_->ndofs_, 1);
  ol.Multiply('N', 'N', 1.0, localSolver_->invamat, tempvec1, 0.0);
  elevec.Multiply('N', 'N', -1.0, localSolver_->kmat, ol, 1.0);

  return;
}  // UpdateInteriorVariablesAndComputeResidual


/*----------------------------------------------------------------------*
 * ProjectPadapField
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::ProjectPadapField(DRT::ELEMENTS::Acou* ele, int newdeg)
{
  dserror("TODO");
}

/*----------------------------------------------------------------------*
 * EstimateError
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::AcouSolEleCalc<distype>::EstimateError(
    DRT::ELEMENTS::AcouSol& ele, Epetra_SerialDenseVector& L)
{
  DRT::UTILS::PolynomialSpace<nsd_> postpoly(
      distype, ele.Degree() + 1, ele.UsesCompletePolynomialSpace());

  LINALG::Matrix<nsd_, 1> xsi;
  int ndofspost = postpoly.Size();

  Epetra_SerialDenseMatrix h(ndofspost * nsd_, ndofspost * nsd_);
  Epetra_SerialDenseVector rhs(ndofspost * nsd_);

  Epetra_SerialDenseMatrix derivs(nsd_, ndofspost);
  Epetra_SerialDenseVector myvalues(shapes_->ndofs_);

  LINALG::Matrix<nsd_, nen_> deriv;
  LINALG::Matrix<nsd_, nsd_> xjm, xji;

  Epetra_SerialDenseVector values(ndofspost);

  Teuchos::RCP<DRT::UTILS::GaussPoints> postquad =
      DRT::UTILS::GaussPointCache::Instance().Create(distype, (ele.Degree() + 1) * 2);

  LINALG::Matrix<nsd_, nen_> xyze;
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(&ele, xyze);

  for (int q = 0; q < postquad->NumPoints(); ++q)
  {
    const double* gpcoord = postquad->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++) xsi(idim) = gpcoord[idim];

    postpoly.Evaluate(xsi, values);
    postpoly.Evaluate_deriv1(xsi, derivs);

    DRT::UTILS::shape_function_deriv1<distype>(xsi, deriv);
    xjm.MultiplyNT(deriv, xyze);
    const double jfac = xji.Invert(xjm) * postquad->Weight(q);

    // transform shape functions derivatives
    for (int i = 0; i < ndofspost; ++i)
    {
      double res[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        res[d] = xji(d, 0) * derivs(0, i);
        for (unsigned int e = 1; e < nsd_; ++e) res[d] += xji(d, e) * derivs(e, i);
      }
      for (unsigned int d = 0; d < nsd_; ++d) derivs(d, i) = res[d];
    }

    shapes_->polySpace_->Evaluate(xsi, myvalues);

    for (unsigned int d = 0; d < nsd_; ++d)
    {
      for (int j = 0; j < ndofspost; ++j) h(d * ndofspost, j + d * ndofspost) += values(j) * jfac;
      for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
        rhs(d * ndofspost) += myvalues(j) * jfac * interiorVelnp_[d * shapes_->ndofs_ + j];
    }

    for (int i = 1; i < ndofspost; ++i)
    {
      for (int j = 0; j < ndofspost; ++j)
      {
        double t = 0;
        for (unsigned int d = 0; d < nsd_; ++d) t += derivs(d, i) * derivs(d, j);
        for (unsigned int d = 0; d < nsd_; ++d) h(i + d * ndofspost, j + d * ndofspost) += t * jfac;
      }
    }

    double vgrad[nsd_ * nsd_] = {0.0};

    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int e = 0; e < nsd_; ++e)
        for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
          vgrad[d + e * nsd_] += myvalues(j) * L(j + (d + e * nsd_) * shapes_->ndofs_);
    for (int i = 1; i < ndofspost; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        for (unsigned int e = 0; e < nsd_; ++e)
          rhs(i + e * ndofspost) += vgrad[d + e * nsd_] * derivs(d, i) * jfac;
  }

  {
    Epetra_SerialDenseSolver inverseH;
    inverseH.SetMatrix(h);
    inverseH.SetVectors(rhs, rhs);
    inverseH.Solve();
  }

  double err_v = 0.0;
  double numerical_post[nsd_] = {0.0};
  double numerical[nsd_] = {0.0};
  double area = 0.0;

  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    const double* gpcoord = shapes_->quadrature_->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++)
    {
      xsi(idim) = gpcoord[idim];
      numerical_post[idim] = 0.0;
      numerical[idim] = 0.0;
    }
    postpoly.Evaluate(xsi, values);
    for (int i = 0; i < ndofspost; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        numerical_post[d] += values(i) * rhs(i + d * ndofspost);
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        numerical[d] += shapes_->shfunct(i, q) * interiorVelnp_(i + d * shapes_->ndofs_);
    area += shapes_->jfac(q);

    for (unsigned int idim = 0; idim < nsd_; ++idim)
      err_v += (numerical_post[idim] - numerical[idim]) * (numerical_post[idim] - numerical[idim]) *
               shapes_->jfac(q);
  }

  err_v /= area;

  return err_v;
}

// template classes
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::tri3>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::AcouSolEleCalc<DRT::Element::nurbs27>;
