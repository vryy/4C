/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_hdg.cpp

\brief main file containing routines for calculation of HDG fluid element

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_hdg.H"
#include "fluid_ele_calc.H"
#include "fluid_ele_parameter.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_parameter_timint.H"
#include "fluid_ele_action.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_geometry/position_array.H"

#include "../drt_mat/newtonianfluid.H"

#include "../drt_fem_general/drt_utils_polynomial.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include <Epetra_SerialDenseSolver.h>



namespace
{
  void zeroMatrix (Epetra_SerialDenseMatrix &mat)
  {
    std::memset(mat.A(), 0, sizeof(double)*mat.M()*mat.N());
  }
}


/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcHDG<distype>::FluidEleCalcHDG()
{}



/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDG<distype>::Evaluate(DRT::ELEMENTS::Fluid*    ele,
                                                      DRT::Discretization & discretization,
                                                      const std::vector<int> & lm,
                                                      Teuchos::ParameterList&    params,
                                                      Teuchos::RCP<MAT::Material> & mat,
                                                      Epetra_SerialDenseMatrix&  elemat1_epetra,
                                                      Epetra_SerialDenseMatrix&  elemat2_epetra,
                                                      Epetra_SerialDenseVector&  elevec1_epetra,
                                                      Epetra_SerialDenseVector&  elevec2_epetra,
                                                      Epetra_SerialDenseVector&  elevec3_epetra,
                                                      const DRT::UTILS::GaussIntegration & ,
                                                      bool                       offdiag)
{
  return this->Evaluate( ele, discretization, lm, params, mat,
                         elemat1_epetra, elemat2_epetra,
                         elevec1_epetra, elevec2_epetra, elevec3_epetra,
                         offdiag );
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::InitializeShapes(const DRT::ELEMENTS::Fluid* ele)
{
  // Check if this is an HDG element, if yes, can initialize...
  if (const DRT::ELEMENTS::FluidHDG * hdgele = dynamic_cast<const DRT::ELEMENTS::FluidHDG*>(ele))
  {
    usescompletepoly_ = hdgele->UsesCompletePolynomialSpace();
    if (shapes_ == Teuchos::null )
      shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(hdgele->Degree(),
                                                                  usescompletepoly_,
                                                                  2*ele->Degree()));
    else if (shapes_->degree_ != unsigned(ele->Degree()) || shapes_->usescompletepoly_ != usescompletepoly_)
      shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(hdgele->Degree(),
                                                                  usescompletepoly_,
                                                                  2*ele->Degree()));

    if (shapesface_ == Teuchos::null)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Degree(), usescompletepoly_, 2 * ele->Degree());
      shapesface_ = Teuchos::rcp(new DRT::UTILS::ShapeValuesFace<distype>(svfparams));
    }
    // TODO: check distype

    if (localSolver_ == Teuchos::null)
      localSolver_ = Teuchos::rcp(new LocalSolver(ele,*shapes_,*shapesface_,usescompletepoly_));
  }
  else
    dserror("Only works for HDG fluid elements");
}



template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDG<distype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
                                                      DRT::Discretization    &   discretization,
                                                      const std::vector<int> &   lm,
                                                      Teuchos::ParameterList&    params,
                                                      Teuchos::RCP<MAT::Material> & mat,
                                                      Epetra_SerialDenseMatrix&  elemat1,
                                                      Epetra_SerialDenseMatrix&  ,
                                                      Epetra_SerialDenseVector&  elevec1,
                                                      Epetra_SerialDenseVector&  ,
                                                      Epetra_SerialDenseVector&  ,
                                                      bool                       offdiag)
{

  InitializeShapes(ele);

  const bool updateLocally = params.get<bool>("needslocalupdate");

  shapes_->Evaluate(*ele);

  ebofoaf_.PutScalar(0);
  eprescpgaf_.PutScalar(0);
  escabofoaf_.PutScalar(0);
  FluidEleCalc<distype>::BodyForce(ele, localSolver_->fldparatimint_->Time(),
      localSolver_->fldpara_->PhysicalType(), ebofoaf_, eprescpgaf_, escabofoaf_);

  //interior body force vector if applicable
  interiorebofoaf_.resize(((nsd_+1)*nsd_+1)*shapes_->ndofs_,0.0);
  if(params.get<bool>("forcing",false))
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"forcing");
    std::vector<int> localDofs = discretization.Dof(1, ele);
    DRT::UTILS::ExtractMyValues(*matrix_state,interiorebofoaf_,localDofs);
  }

  ReadGlobalVectors(*ele, discretization, lm, updateLocally);

  // solves the local problem of the nonlinear iteration before
  if (updateLocally) {
    localSolver_->ComputeInteriorResidual(mat, interiorVal_, interiorAcc_, traceVal_[0], ebofoaf_, interiorebofoaf_);
    localSolver_->ComputeInteriorMatrices(mat, false);

    dsassert(nfaces_ == static_cast<unsigned int>(ele->NumFace()), "Internal error");

    // loop over faces
    for (unsigned int f=0; f<nfaces_; ++f) {
      shapesface_->EvaluateFace(*ele,f);
      localSolver_->ComputeFaceResidual(f, mat, interiorVal_, traceVal_, elevec1);
      localSolver_->ComputeFaceMatrices(f, mat, false, elemat1);
    }

    localSolver_->EliminateVelocityGradient(elemat1);
    localSolver_->SolveResidual();
    UpdateSecondarySolution(*ele, discretization, localSolver_->gUpd, localSolver_->upUpd);
  }

  zeroMatrix(elemat1);
  zeroMatrix(elevec1);
  localSolver_->ComputeInteriorResidual(mat, interiorVal_, interiorAcc_, traceVal_[0], ebofoaf_, interiorebofoaf_);
  localSolver_->ComputeInteriorMatrices(mat, updateLocally);
  for (unsigned int f=0; f<nfaces_; ++f) {
    shapesface_->EvaluateFace(*ele,f);
    localSolver_->ComputeFaceResidual(f, mat, interiorVal_, traceVal_, elevec1);
    localSolver_->ComputeFaceMatrices(f, mat, updateLocally, elemat1);
  }

  if (!updateLocally)
    localSolver_->EliminateVelocityGradient(elemat1);

  localSolver_->CondenseLocalPart(elemat1,elevec1);

  elevec1.Scale(1./localSolver_->fldparatimint_->AlphaF());

  return 0;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::
ReadGlobalVectors(const DRT::Element     & ele,
                  DRT::Discretization    & discretization,
                  const std::vector<int> & lm,
                  const bool               updateLocally)
{
  // read the HDG solution vector (for traces)
  traceVal_.resize(1+nfaces_*nsd_*shapesface_->nfdofs_);
  interiorVal_.resize(((nsd_+1)*nsd_+1)*shapes_->ndofs_);
  interiorAcc_.resize(((nsd_+1)*nsd_+1)*shapes_->ndofs_);
  dsassert(lm.size() == traceVal_.size(), "Internal error");
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState("velaf");
  DRT::UTILS::ExtractMyValues(*matrix_state,traceVal_,lm);

  // read the interior values from solution vector
  matrix_state = discretization.GetState(1,"intvelaf");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  DRT::UTILS::ExtractMyValues(*matrix_state,interiorVal_,localDofs);

  matrix_state = discretization.GetState(1,"intaccam");
  DRT::UTILS::ExtractMyValues(*matrix_state,interiorAcc_,localDofs);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::
UpdateSecondarySolution(const DRT::Element     & ele,
                        DRT::Discretization    & discretization,
                        const Epetra_SerialDenseVector &updateG,
                        const Epetra_SerialDenseVector &updateUp)
{
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnp");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  dsassert(localDofs.size() == static_cast<std::size_t>(updateG.Length()+updateUp.Length()),
           "Internal error");

  // update vector content by making the vector writeable (need to adjust in calling site before
  // clearing the state when used in parallel)
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);

  for (unsigned int i=0; i<localDofs.size(); ++i) {
    const int lid = intdofcolmap->LID(localDofs[i]);
    double update = i<nsd_*nsd_*shapes_->ndofs_ ? updateG(i) : updateUp(i-nsd_*nsd_*shapes_->ndofs_);

    const double valfac = 1./localSolver_->fldparatimint_->AlphaF();
    secondary[lid] += update * valfac;

    // write the update back into the local vectors (when doing local update,
    // we do not re-read from the global vectors)
    interiorVal_[i] += update;

    const double accfac = localSolver_->fldparatimint_->AlphaM() * valfac /
        (localSolver_->fldparatimint_->Dt() * localSolver_->fldparatimint_->Gamma());
    interiorAcc_[i] += update * accfac;
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDG<distype>::EvaluateService(
    DRT::ELEMENTS::Fluid*     ele,
    Teuchos::ParameterList&   params,
    Teuchos::RCP<MAT::Material> & mat,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params,"action");

  switch(act)
  {
  case FLD::calc_fluid_error:
  {
    // compute error for a known analytical solution
    return ComputeError(ele, params, mat, discretization, lm, elevec1);
  }
  break;
  case FLD::interpolate_hdg_to_node:
  {
    return InterpolateSolutionToNodes(
        ele,
        discretization,
        elevec1);
    break;
  }
  case FLD::interpolate_hdg_for_hit:
  {
    InterpolateSolutionForHIT(
        ele,
        discretization,
        elevec1);
    break;
  }
  case FLD::project_hdg_force_on_dof_vec_for_hit:
  {
    ProjectForceOnDofVecForHIT(
        ele,
        discretization,
        elevec1,
        elevec2);
    break;
  }
  case FLD::project_hdg_initial_field_for_hit:
  {
    ProjectInitialFieldForHIT(
        ele,
        discretization,
        elevec1,
        elevec2,
        elevec3);
    break;
  }
  case FLD::project_fluid_field:
  {
    return ProjectField(
        ele,
        params,
        mat,
        discretization,
        lm,
        elevec1,
        elevec2);
    break;
  }
  default:
    dserror("Unknown type of action for FluidHDG");
    break;
  } // end of switch(act)

  return 0;
}



template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDG<distype>::ComputeError(
    DRT::ELEMENTS::Fluid*                ele,
    Teuchos::ParameterList&              params,
    Teuchos::RCP<MAT::Material>&         mat,
    DRT::Discretization&                 discretization,
    std::vector<int>&                    lm,
    Epetra_SerialDenseVector&            elevec)
{
  InitializeShapes(ele);

  shapes_->Evaluate(*ele);
  const double time = localSolver_->fldparatimint_->Time();

  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnp");
  std::vector<int> localDofs = discretization.Dof(1, ele);
  std::vector<double> vecValues(localDofs.size());

  for (unsigned int i=0; i<localDofs.size(); ++i) {
    const int lid = matrix_state->Map().LID(localDofs[i]);
    vecValues[i] = (*matrix_state)[lid];
  }

  // analytic solution
  LINALG::Matrix<nsd_,1> u(true);
  double p = 0.0;
  LINALG::Matrix<nsd_,nsd_> dervel(true);
  LINALG::Matrix<nsd_,1> xyz(true);

  const INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params,"calculate error");
  const int calcerrfunctno = DRT::INPUT::get<INPAR::FLUID::CalcError>(params,"error function number");

  double err_u = 0., err_p = 0., err_h = 0., norm_u = 0., norm_p = 0., norm_h = 0.;
  for (unsigned int q=0; q<shapes_->nqpoints_; ++q)
  {
    const double jfac = shapes_->jfac(q);
    double numericalGrad[nsd_][nsd_];
    double numerical[nsd_+1];
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e) {
        numericalGrad[d][e] = 0;
        for (unsigned int i=0; i<shapes_->ndofs_; ++i)
          numericalGrad[d][e] += shapes_->shfunct(i,q) * vecValues[(d*nsd_+e)*shapes_->ndofs_+i];
      }
    for (unsigned int d=0; d<=nsd_; ++d) {
      numerical[d] = 0.;
      for (unsigned int i=0; i<shapes_->ndofs_; ++i)
        numerical[d] += shapes_->shfunct(i,q) * vecValues[(nsd_*nsd_+d)*shapes_->ndofs_+i];
    }
    for (unsigned int d=0; d<nsd_; ++d)
      xyz(d) = shapes_->xyzreal(d,q);

    FluidEleCalc<distype>::EvaluateAnalyticSolutionPoint(xyz, time, calcerr, calcerrfunctno, mat, u, p, dervel);

    for (unsigned int d=0; d<nsd_; ++d)
      err_u += (u(d) - numerical[d]) * (u(d) - numerical[d]) * jfac;
    err_p += (p - numerical[nsd_]) * (p - numerical[nsd_]) * jfac;
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e)
        err_h += (dervel(d,e) - numericalGrad[d][e]) * (dervel(d,e) - numericalGrad[d][e]) * jfac;
    for (unsigned int d=0; d<nsd_; ++d)
      norm_u += u(d) * u(d) * jfac;
    norm_p += p * p * jfac;
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e)
        norm_h += dervel(e,d) * dervel(e,d) * jfac;
  }

  elevec[0] += err_u;
  elevec[1] += err_p;
  elevec[2] += err_h;
  elevec[3] += norm_u;
  elevec[4] += norm_p;
  elevec[5] += norm_h;

  return 0;
};



/// projection of function field
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDG<distype>::ProjectField(
    DRT::ELEMENTS::Fluid*                ele,
    Teuchos::ParameterList&              params,
    Teuchos::RCP<MAT::Material>&         mat,
    DRT::Discretization&                 discretization,
    std::vector<int>&                    lm,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2)
{
  InitializeShapes(ele);

  shapes_->Evaluate(*ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 ||
           elevec2.M() == static_cast<int>((nsd_*nsd_+nsd_+1)*shapes_->ndofs_), "Wrong size in project vector 2");

  // get function
  const int *initfield = params.getPtr<int>("initfield");
  const int *startfunc = params.getPtr<int>("startfuncno");
  double *time = params.getPtr<double>("time");

  double avgpre = 0., vol = 0.;
  if (elevec2.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(View, elevec2.A(), shapes_->ndofs_, shapes_->ndofs_, nsd_*nsd_+nsd_+1, false);
    zeroMatrix(localMat);

    // create mass matrix for interior by looping over quadrature points
    for (unsigned int q=0; q<shapes_->nqpoints_; ++q )
    {
      const double fac = shapes_->jfac(q);
      LINALG::Matrix<nsd_,1> xyz(false);
      for (unsigned int d=0; d<nsd_; ++d)
        xyz(d) = shapes_->xyzreal(d,q);
      LINALG::Matrix<nsd_,1>    u(false);
      LINALG::Matrix<nsd_,nsd_> grad(true); //is not necessarily set in EvaluateAll
      double p;
      dsassert(initfield != NULL && startfunc != NULL,
               "initfield or startfuncno not set for initial value");
      EvaluateAll(*startfunc, INPAR::FLUID::InitialField(*initfield), xyz, u, grad, p);

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i=0; i<shapes_->ndofs_; ++i) {
        // mass matrix part
        localSolver_->massPart(i,q)  = shapes_->shfunct(i,q);
        localSolver_->massPartW(i,q) = shapes_->shfunct(i,q) * fac;

        for (unsigned int d=0; d<nsd_; ++d)
          for (unsigned int e=0; e<nsd_; ++e)
            localMat(i,d*nsd_+e) += shapes_->shfunct(i,q) * grad(d,e) * fac;
        for (unsigned int d=0; d<nsd_; ++d)
          localMat(i,nsd_*nsd_+d) += shapes_->shfunct(i,q) * u(d) * fac;
        localMat(i,nsd_*nsd_+nsd_) += shapes_->shfunct(i,q) * p * fac;
      }

      avgpre += p * fac;
      vol += fac;
    }
    localSolver_->massMat.Multiply('N', 'T', 1., localSolver_->massPart,localSolver_->massPartW, 0.);

    // solve mass matrix system, return values in localMat = elevec2 correctly ordered
    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(localSolver_->massMat);
    inverseMass.SetVectors(localMat,localMat);
    inverseMass.Solve();
  }

  Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
  Epetra_SerialDenseMatrix trVec(shapesface_->nfdofs_, nsd_);
  dsassert(elevec1.M() == static_cast<int>(nsd_*shapesface_->nfdofs_) ||
           elevec1.M() == 1+static_cast<int>(nfaces_*nsd_*shapesface_->nfdofs_), "Wrong size in project vector 1");

  const unsigned int *faceConsider = params.getPtr<unsigned int>("faceconsider");
  Teuchos::Array<int> *functno = params.getPtr<Teuchos::Array<int> >("funct");
  Teuchos::Array<int> *onoff = params.getPtr<Teuchos::Array<int> >("onoff");

  for (unsigned int face=0; face<nfaces_; ++face)
  {
    // check whether we are in the project phase for all faces or for boundary values
    if (initfield == NULL) {
      dsassert(faceConsider != NULL, "Unsupported operation");
      if (*faceConsider != face)
        continue;
    }
    shapesface_->EvaluateFace(*ele,face);
    zeroMatrix(mass);
    zeroMatrix(trVec);


    for (unsigned int q=0; q<shapesface_->nqpoints_; ++q) {
      const double fac = shapesface_->jfac(q);
      LINALG::Matrix<nsd_,1> xyz(false);
      for (unsigned int d=0; d<nsd_; ++d)
        xyz(d) = shapesface_->xyzreal(d,q);
      LINALG::Matrix<nsd_,1> u(false);
      if (initfield != NULL)
        EvaluateVelocity(*startfunc, INPAR::FLUID::InitialField(*initfield), xyz, u);
      else
      {
        dsassert(functno != NULL && time != NULL && onoff != NULL,
                 "No array with functions given");
        for (unsigned int d=0; d<nsd_; ++d) {
          if ((*onoff)[d] == 0)
            continue;
          const int funct_num = (*functno)[d];
          if (funct_num > 0)
            u(d) = DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(d, xyz.A(), *time);
        }
      }

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i=0; i<shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j=0; j<shapesface_->nfdofs_; ++j)
          mass(i,j) += shapesface_->shfunct(i,q) * shapesface_->shfunct(j,q) * fac;

        for (unsigned int d=0; d<nsd_; ++d)
          trVec(i,d) += shapesface_->shfunct(i,q) * u(d) * fac;
      }
    }

    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec,trVec);
    inverseMass.Solve();

    if (initfield != NULL)
      for (unsigned int d=0; d<nsd_; ++d)
        for (unsigned int i=0; i<shapesface_->nfdofs_; ++i)
          elevec1(1+face*shapesface_->nfdofs_*nsd_+d*shapesface_->nfdofs_+i) = trVec(i,d);
    else
      for (unsigned int d=0; d<nsd_; ++d)
        for (unsigned int i=0; i<shapesface_->nfdofs_; ++i)
          elevec1(d*shapesface_->nfdofs_+i) = trVec(i,d);
  }
  if (initfield != NULL)
    elevec1(0) = avgpre/vol;

  return 0;
}



template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDG<distype>::InterpolateSolutionToNodes(
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization&                 discretization,
    Epetra_SerialDenseVector&            elevec1)
{
  InitializeShapes(ele);
  dsassert(elevec1.M() == (int)nen_*(2*nsd_+1)+1, "Vector does not have correct size");
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  // get local solution values
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnp");
  std::vector<int> localDofs = discretization.Dof(1, ele);
  std::vector<double> solvalues (localDofs.size());

  for (unsigned int i=0; i<solvalues.size(); ++i) {
    const int lid = matrix_state->Map().LID(localDofs[i]);
    solvalues[i] = (*matrix_state)[lid];
  }

  for (unsigned int i=0; i<nen_; ++i) {
    // evaluate shape polynomials in node
    for (unsigned int idim=0;idim<nsd_;idim++)
      shapes_->xsi(idim) = locations(idim,i);
    shapes_->polySpace_->Evaluate(shapes_->xsi,values);

    // compute values for velocity and pressure by summing over all basis functions
    for (unsigned int d=0; d<=nsd_; ++d) {
      double sum = 0;
      for (unsigned int k=0; k<shapes_->ndofs_; ++k)
        sum += values(k) * solvalues[(nsd_*nsd_+d)*shapes_->ndofs_+k];
      elevec1(d*nen_+i) = sum;
    }
  }

  // get trace solution values
  locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace
                (DRT::UTILS::DisTypeToFaceShapeType<distype>::shape);
  matrix_state = discretization.GetState(0,"velnp");

  //we have always two dofsets
  Element::LocationArray la(2);
  ele->LocationVector(discretization,la,false);
  localDofs = la[0].lm_;
  solvalues.resize (localDofs.size());

  for (unsigned int i=0; i<solvalues.size(); ++i) {
    const int lid = matrix_state->Map().LID(localDofs[i]);
    solvalues[i] = (*matrix_state)[lid];
  }

  Epetra_SerialDenseVector fvalues(shapesface_->nfdofs_);
  for (unsigned int f=0; f<nfaces_; ++f)
  {
    LINALG::Matrix<nsd_-1,DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace>   xsishuffle(true);

    for (int i=0; i<DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
    {
           // evaluate shape polynomials in node
      for (unsigned int idim=0;idim<nsd_-1;idim++)
      {
        if(ele->Faces()[f]->ParentMasterElement() == ele)
          xsishuffle(idim,i) = locations(idim,i);
        else
          xsishuffle(idim,ele->Faces()[f]->GetLocalTrafoMap()[i]) = locations(idim,i);
      }
    }

    for (int i=0; i<DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
    {
      // evaluate shape polynomials in node
      for (unsigned int idim=0;idim<nsd_-1;idim++)
        shapesface_->xsi(idim) = xsishuffle(idim,i);
      shapesface_->polySpace_->Evaluate(shapesface_->xsi,fvalues);

      // compute values for velocity and pressure by summing over all basis functions
      for (unsigned int d=0; d<nsd_; ++d) {
        double sum = 0;
        for (unsigned int k=0; k<shapesface_->nfdofs_; ++k)
          sum += fvalues(k) * solvalues[1+f*nsd_*shapesface_->nfdofs_+d*shapesface_->nfdofs_+k];
        elevec1((nsd_+1+d)*nen_+shapesface_->faceNodeOrder[f][i]) = sum;
      }
    }
  }

  elevec1((2*nsd_+1)*nen_) = solvalues[0];


  return 0;
}

/*----------------------------------------------------------------------*
 * interpolate solution for postprocessing of hit              bk 03/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDG<distype>::InterpolateSolutionForHIT(
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization&                 discretization,
    Epetra_SerialDenseVector&            elevec1)
{
  InitializeShapes(ele);
  //get coordinates of hex 8
  LINALG::Matrix<nsd_,nen_> xyze(true);

  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze);

  const int numsamppoints = 5;
  dsassert(elevec1.M() == numsamppoints*numsamppoints*numsamppoints*6, "Vector does not have correct size");
  //sampling locations in 1D in parent domain
  double loc1D[numsamppoints] = {-0.8, -0.4, 0.0, 0.4, 0.8};
  Epetra_SerialDenseMatrix locations(3,125);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  int l=0;
  for (int i=0;i<numsamppoints;i++)
    for (int j=0;j<numsamppoints;j++)
      for (int k=0;k<numsamppoints;k++)
      {
        locations(0,l)=loc1D[k];
        locations(1,l)=loc1D[j];
        locations(2,l)=loc1D[i];
        l++;
      }
  // get local solution values
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnp");
  std::vector<int> localDofs = discretization.Dof(1, ele);
  std::vector<double> solvalues (localDofs.size());

  for (unsigned int i=0; i<solvalues.size(); ++i)
  {
    const int lid = matrix_state->Map().LID(localDofs[i]);
    solvalues[i] = (*matrix_state)[lid];
  }

  for (unsigned int i=0; i<numsamppoints*numsamppoints*numsamppoints; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim=0;idim<nsd_;idim++)
      shapes_->xsi(idim) = locations(idim,i);
    shapes_->polySpace_->Evaluate(shapes_->xsi,values);

    // compute values for velocity and pressure by summing over all basis functions
    for (unsigned int d=0; d<nsd_; ++d)
    {
      double sum = 0;
      for (unsigned int k=0; k<shapes_->ndofs_; ++k)
        sum += values(k) * solvalues[(nsd_*nsd_+d)*shapes_->ndofs_+k];
      elevec1(6*i+d) = sum;
    }

    //also save coordinates
    LINALG::Matrix<nen_,1> myfunct;
    DRT::UTILS::shape_function<distype>(shapes_->xsi,myfunct);

    LINALG::Matrix<nsd_,1> mypoint(true);
    mypoint.MultiplyNN(xyze,myfunct);

    for (unsigned int d=0; d<nsd_; ++d)
      elevec1(6*i+d+3) = mypoint(d);
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * project force for hit                                       bk 03/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDG<distype>::ProjectForceOnDofVecForHIT(
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization&                 discretization,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2)
{
  const int numsamppoints = 5;
//  dsassert(elevec1.M() == numsamppoints*numsamppoints*numsamppoints*6, "Vector does not have correct size");
  //sampling locations in 1D in parent domain
  double loc1D[numsamppoints] = {-0.8, -0.4, 0.0, 0.4, 0.8};

  Epetra_SerialDenseMatrix locations;
#ifdef DEBUG
  locations.Shape(3,125);
  int l=0;
  for (int i=0;i<numsamppoints;i++)
    for (int j=0;j<numsamppoints;j++)
      for (int k=0;k<numsamppoints;k++)
      {
        locations(0,l)=loc1D[k];
        locations(1,l)=loc1D[j];
        locations(2,l)=loc1D[i];
        l++;
      }
#endif

  std::vector<DRT::UTILS::LagrangePolynomial> poly1d;
  const unsigned int degree = 4;
  std::vector<double> points(degree);
  for (unsigned int i=0; i<=degree; ++i)
  {
    for (unsigned int j=0, c=0; j<=degree; ++j)
      if (i!=j)
      {
        points[c] = loc1D[j];
        ++c;
      }
    poly1d.push_back(DRT::UTILS::LagrangePolynomial(points, loc1D[i]));
  }

  DRT::UTILS::PolynomialSpaceTensor<nsd_,DRT::UTILS::LagrangePolynomial>
  poly(poly1d);

#ifdef DEBUG
  //check if we have the right number of polynomials
  if(poly.Size()!=125)
    dserror("wrong number of polynomials");
#endif

  InitializeShapes(ele);
  shapes_->Evaluate(*ele);


  if (elevec1.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(View, elevec1.A(), shapes_->ndofs_, shapes_->ndofs_, nsd_*nsd_+nsd_+1, false);
    zeroMatrix(localMat);

    // create mass matrix for interior by looping over quadrature points
    for (unsigned int q=0; q<shapes_->nqpoints_; ++q )
    {
      LINALG::Matrix<nsd_,1>    f(false);
      const double fac = shapes_->jfac(q);
      Epetra_SerialDenseVector values(numsamppoints*numsamppoints*numsamppoints);
      LINALG::Matrix<nsd_,1>        xsi(false);
      for(unsigned int sdm=0;sdm<nsd_;sdm++)
        xsi(sdm) = shapes_->quadrature_->Point(q)[sdm];

      poly.Evaluate(xsi,values);
      // compute values for force and coordinates by summing over all basis functions
      for (unsigned int d=0; d<nsd_; ++d)
      {
        double sum = 0.0;
        for (unsigned int k=0; k<numsamppoints*numsamppoints*numsamppoints; ++k)
          sum += values(k) * elevec2(6*k+d);
        f(d) = sum;

#ifdef DEBUG
        //check plausibility via comparison of quadrature coordinates
        sum = 0.0;
        for (unsigned int k=0; k<numsamppoints*numsamppoints*numsamppoints; ++k)
          sum += values(k) * locations(d,k);
        if(not(sum+1e-9 > xsi(d) and sum-1e-9 < xsi(d)))
        {
          std::cout << "Gauss point:  " << xsi(d) << "  coordinate:  " << sum << std::endl;
          dserror("Plausibility check failed! Problem might be sequence of polynomials");
        }
#endif
      }

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i=0; i<shapes_->ndofs_; ++i) {
        // mass matrix part
        localSolver_->massPart(i,q)  = shapes_->shfunct(i,q);
        localSolver_->massPartW(i,q) = shapes_->shfunct(i,q) * fac;

        for (unsigned int d=0; d<nsd_; ++d)
          localMat(i,nsd_*nsd_+d) += shapes_->shfunct(i,q) * f(d) * fac;
      }
    }
    localSolver_->massMat.Multiply('N', 'T', 1., localSolver_->massPart,localSolver_->massPartW, 0.);

    // solve mass matrix system, return values in localMat = elevec2 correctly ordered
    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(localSolver_->massMat);
    inverseMass.SetVectors(localMat,localMat);
    inverseMass.Solve();
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * project force for hit                                       bk 03/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcHDG<distype>::ProjectInitialFieldForHIT(
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization&                 discretization,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2,
    Epetra_SerialDenseVector&            elevec3)
{
  const int numsamppoints = 5;
//  dsassert(elevec1.M() == numsamppoints*numsamppoints*numsamppoints*6, "Vector does not have correct size");
  //sampling locations in 1D in parent domain
  double loc1D[numsamppoints] = {-0.8, -0.4, 0.0, 0.4, 0.8};

  Epetra_SerialDenseMatrix locations;
#ifdef DEBUG
  locations.Shape(3,125);
  int l=0;
  for (int i=0;i<numsamppoints;i++)
    for (int j=0;j<numsamppoints;j++)
      for (int k=0;k<numsamppoints;k++)
      {
        locations(0,l)=loc1D[k];
        locations(1,l)=loc1D[j];
        locations(2,l)=loc1D[i];
        l++;
      }
#endif

  std::vector<DRT::UTILS::LagrangePolynomial> poly1d;
  const unsigned int degree = 4;
  std::vector<double> points(degree);
  for (unsigned int i=0; i<=degree; ++i)
  {
    for (unsigned int j=0, c=0; j<=degree; ++j)
      if (i!=j)
      {
        points[c] = loc1D[j];
        ++c;
      }
    poly1d.push_back(DRT::UTILS::LagrangePolynomial(points, loc1D[i]));
  }

  DRT::UTILS::PolynomialSpaceTensor<nsd_,DRT::UTILS::LagrangePolynomial>
  poly(poly1d);

  InitializeShapes(ele);
  shapes_->Evaluate(*ele);


  if (elevec1.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(View, elevec1.A(), shapes_->ndofs_, shapes_->ndofs_, nsd_*nsd_+nsd_+1, false);
    zeroMatrix(localMat);

    // create mass matrix for interior by looping over quadrature points
    for (unsigned int q=0; q<shapes_->nqpoints_; ++q )
    {
      LINALG::Matrix<nsd_,1>    f(false);
      const double fac = shapes_->jfac(q);
      Epetra_SerialDenseVector values(numsamppoints*numsamppoints*numsamppoints);
      LINALG::Matrix<nsd_,1>        xsi(false);
      for(unsigned int sdm=0;sdm<nsd_;sdm++)
        xsi(sdm) = shapes_->quadrature_->Point(q)[sdm];

      poly.Evaluate(xsi,values);
      // compute values for force and coordinates by summing over all basis functions
      for (unsigned int d=0; d<nsd_; ++d)
      {
        double sum = 0.0;
        for (unsigned int k=0; k<numsamppoints*numsamppoints*numsamppoints; ++k)
          sum += values(k) * elevec2(6*k+d);
        f(d) = sum;

#ifdef DEBUG
        //check plausibility via comparison of quadrature coordinates
        sum = 0.0;
        for (unsigned int k=0; k<numsamppoints*numsamppoints*numsamppoints; ++k)
          sum += values(k) * locations(d,k);
        if(not(sum+1e-9 > xsi(d) and sum-1e-9 < xsi(d)))
        {
          std::cout << "Gauss point:  " << xsi(d) << "  coordinate:  " << sum << std::endl;
          dserror("Plausibility check failed! Problem might be sequence of polynomials");
        }
#endif
      }

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i=0; i<shapes_->ndofs_; ++i)
      {
        // mass matrix part
        localSolver_->massPart(i,q)  = shapes_->shfunct(i,q);
        localSolver_->massPartW(i,q) = shapes_->shfunct(i,q) * fac;

        for (unsigned int d=0; d<nsd_; ++d)
          localMat(i,nsd_*nsd_+d) += shapes_->shfunct(i,q) * f(d) * fac;
      }
    }
    localSolver_->massMat.Multiply('N', 'T', 1., localSolver_->massPart,localSolver_->massPartW, 0.);

    // solve mass matrix system, return values in localMat = elevec2 correctly ordered
    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(localSolver_->massMat);
    inverseMass.SetVectors(localMat,localMat);
    inverseMass.Solve();
  }

  //traces
  Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
  Epetra_SerialDenseMatrix trVec(shapesface_->nfdofs_, nsd_);
  dsassert(elevec3.M() == static_cast<int>(nsd_*shapesface_->nfdofs_) ||
           elevec3.M() == 1+static_cast<int>(nfaces_*nsd_*shapesface_->nfdofs_), "Wrong size in project vector 1");

  for (unsigned int face=0; face<nfaces_; ++face)
  {
    shapesface_->EvaluateFace(*ele,face);
    zeroMatrix(mass);
    zeroMatrix(trVec);

    LINALG::Matrix<nsd_,nsd_> trafo;
    LINALG::SerialDenseMatrix faceQPoints;
    DRT::UTILS::BoundaryGPToParentGP<nsd_>(faceQPoints,trafo,*shapesface_->quadrature_,distype,
        DRT::UTILS::getEleFaceShapeType(distype,face), face);

    for (unsigned int q=0; q<shapesface_->nqpoints_; ++q)
    {
      const double fac = shapesface_->jfac(q);
      LINALG::Matrix<nsd_,1>        xsi(false);

      //use the location of the quadrature point in the parent element to evaluate the polynomial
      for (unsigned int d=0; d<nsd_; ++d)
        xsi(d) = faceQPoints(q,d);

      LINALG::Matrix<nsd_,1> u(false);

      Epetra_SerialDenseVector values(numsamppoints*numsamppoints*numsamppoints);

      poly.Evaluate(xsi,values);
      // compute values for force and coordinates by summing over all basis functions
      for (unsigned int d=0; d<nsd_; ++d)
      {
        double sum = 0.0;
        for (unsigned int k=0; k<numsamppoints*numsamppoints*numsamppoints; ++k)
          sum += values(k) * elevec2(6*k+d);
        u(d) = sum;
      }

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i=0; i<shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j=0; j<shapesface_->nfdofs_; ++j)
          mass(i,j) += shapesface_->shfunct(i,q) * shapesface_->shfunct(j,q) * fac;

        for (unsigned int d=0; d<nsd_; ++d)
          trVec(i,d) += shapesface_->shfunct(i,q) * u(d) * fac;
      }
    }

    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec,trVec);
    inverseMass.Solve();


    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int i=0; i<shapesface_->nfdofs_; ++i)
        elevec3(1+face*shapesface_->nfdofs_*nsd_+d*shapesface_->nfdofs_+i) = trVec(i,d);

  }

  elevec3(0) = 0.0;

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::EvaluateVelocity(const int startfunc,
    const INPAR::FLUID::InitialField initfield,
    const LINALG::Matrix<nsd_,1>  &xyz,
    LINALG::Matrix<nsd_,1>        &u) const
{
  // pass on dummy entries (costs a little but will not be significant)
  LINALG::Matrix<nsd_,nsd_> grad(true);
  double p;
  EvaluateAll(startfunc, initfield, xyz, u, grad, p);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::EvaluateAll(const int startfunc,
    const INPAR::FLUID::InitialField initfield,
    const LINALG::Matrix<nsd_,1>    &xyz,
    LINALG::Matrix<nsd_,1>          &u,
    LINALG::Matrix<nsd_,nsd_>       &grad,
    double                          &p) const
{
  switch (initfield)
  {
  case INPAR::FLUID::initfield_beltrami_flow:
    // check whether present flow is indeed three-dimensional
  {
    if (nsd_!=3) dserror("Beltrami flow is a three-dimensional flow!");

    // set constants for analytical solution
    const double a = M_PI/4.0;
    const double d = M_PI/2.0;
    u(0) = -a * ( std::exp(a*xyz(0)) * std::sin(a*xyz(1) + d*xyz(2)) +
                  std::exp(a*xyz(2)) * std::cos(a*xyz(0) + d*xyz(1)) );
    u(1) = -a * ( std::exp(a*xyz(1)) * std::sin(a*xyz(2) + d*xyz(0)) +
                  std::exp(a*xyz(0)) * std::cos(a*xyz(1) + d*xyz(2)) );
    u(2) = -a * ( std::exp(a*xyz(2)) * std::sin(a*xyz(0) + d*xyz(1)) +
                  std::exp(a*xyz(1)) * std::cos(a*xyz(2) + d*xyz(0)) );

    grad(0,0) = -a * ( a * std::exp(a*xyz(0)) * std::sin(a*xyz(1) + d*xyz(2)) -
                        a * std::exp(a*xyz(2)) * std::sin(a*xyz(0) + d*xyz(1)) );
    grad(0,1) = -a * ( a * std::exp(a*xyz(0)) * std::cos(a*xyz(1) + d*xyz(2)) -
                        d * std::exp(a*xyz(2)) * std::sin(a*xyz(0) + d*xyz(1)) );
    grad(0,2) = -a * ( d * std::exp(a*xyz(0)) * std::cos(a*xyz(1) + d*xyz(2)) +
                        a * std::exp(a*xyz(2)) * std::cos(a*xyz(0) + d*xyz(1)) );
    grad(1,0) = -a * ( d * std::exp(a*xyz(1)) * std::cos(a*xyz(2) + d*xyz(0)) +
                        a * std::exp(a*xyz(0)) * std::cos(a*xyz(1) + d*xyz(2)) );
    grad(1,1) = -a * ( a * std::exp(a*xyz(1)) * std::sin(a*xyz(2) + d*xyz(0)) -
                        a * std::exp(a*xyz(0)) * std::sin(a*xyz(1) + d*xyz(2)) );
    grad(1,2) = -a * ( a * std::exp(a*xyz(1)) * std::cos(a*xyz(2) + d*xyz(0)) -
                        d * std::exp(a*xyz(0)) * std::sin(a*xyz(1) + d*xyz(2)) );
    grad(2,0) = -a * ( a * std::exp(a*xyz(2)) * std::cos(a*xyz(0) + d*xyz(1)) -
                        d * std::exp(a*xyz(1)) * std::sin(a*xyz(2) + d*xyz(0)) );
    grad(2,1) = -a * ( d * std::exp(a*xyz(2)) * std::cos(a*xyz(0) + d*xyz(1)) +
                        a * std::exp(a*xyz(1)) * std::cos(a*xyz(2) + d*xyz(0)) );
    grad(2,2) = -a * ( a * std::exp(a*xyz(2)) * std::sin(a*xyz(0) + d*xyz(1)) -
                        a * std::exp(a*xyz(1)) * std::sin(a*xyz(2) + d*xyz(0)) );

    p = -a*a/2.0 *
      ( std::exp(2.0*a*xyz(0))
        + std::exp(2.0*a*xyz(1))
        + std::exp(2.0*a*xyz(2))
        + 2.0 * std::sin(a*xyz(0) + d*xyz(1)) * std::cos(a*xyz(2) + d*xyz(0)) * std::exp(a*(xyz(1)+xyz(2)))
        + 2.0 * std::sin(a*xyz(1) + d*xyz(2)) * std::cos(a*xyz(0) + d*xyz(1)) * std::exp(a*(xyz(2)+xyz(0)))
        + 2.0 * std::sin(a*xyz(2) + d*xyz(0)) * std::cos(a*xyz(1) + d*xyz(2)) * std::exp(a*(xyz(0)+xyz(1)))
        );
  }
  break;

  case INPAR::FLUID::initfield_field_by_function:
  case INPAR::FLUID::initfield_disturbed_field_from_function:
  {
    for(unsigned int index=0; index<nsd_; ++index)
      u(index) = DRT::Problem::Instance()->Funct(startfunc-1).Evaluate(index,xyz.A(),0);
    p = DRT::Problem::Instance()->Funct(startfunc-1).Evaluate(nsd_,xyz.A(),0);
  }
  break;

  default:
    dserror("Given field %i not yet implemented.", initfield);
    break;
  }

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcHDG<distype> *
DRT::ELEMENTS::FluidEleCalcHDG<distype>::Instance( bool create )
{
  static FluidEleCalcHDG<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcHDG<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}




/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}






template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcHDG<distype>::LocalSolver::
LocalSolver(const DRT::ELEMENTS::Fluid* ele,
            const DRT::UTILS::ShapeValues<distype> &shapeValues,
            DRT::UTILS::ShapeValuesFace<distype> &shapeValuesFace,
            bool completepoly)
:
ndofs_ (shapeValues.ndofs_),
stokes (false),
shapes_(shapeValues),
shapesface_(shapeValuesFace)
{
  uuMat.Shape((nsd_+1)*ndofs_,(nsd_+1)*ndofs_);
  uuMatFinal.Shape((nsd_+1)*ndofs_,(nsd_+1)*ndofs_);
  guMat.Shape(nsd_*ndofs_,ndofs_);
  ugMat.Shape(nsd_*ndofs_,ndofs_);

  int onfdofs = 0;
  for(unsigned int i=0; i<nfaces_; ++i)
  {
    shapesface_.EvaluateFace(*ele,i);
    onfdofs += shapesface_.nfdofs_;
  }
  onfdofs *= nsd_;

  gfMat.Shape(nsd_*nsd_*ndofs_,1+onfdofs);
  fgMat.Shape(gfMat.N(), gfMat.M());
  ufMat.Shape((nsd_+1)*ndofs_,1+onfdofs);
  fuMat.Shape(ufMat.N(), ufMat.M());

  massPart.Shape(ndofs_,shapes_.nqpoints_);
  massPartW.Shape(ndofs_,shapes_.nqpoints_);
  gradPart.Shape(nsd_*ndofs_,shapes_.nqpoints_);
  uPart.Shape(ndofs_*nsd_,shapes_.nqpoints_);

  massMat.Shape(ndofs_,ndofs_);
  uuconv.Shape(ndofs_*nsd_,ndofs_*nsd_);
  tmpMat.Shape(ndofs_*nsd_,ndofs_*nsd_);
  tmpMatGrad.Shape(nsd_*ndofs_,ndofs_);

  velnp.Shape(nsd_,shapes_.nqpoints_);
  gRes.Resize(nsd_*nsd_*ndofs_);
  upRes.Resize((nsd_+1)*ndofs_);
  gUpd.Resize(nsd_*nsd_*ndofs_);
  upUpd.Resize((nsd_+1)*ndofs_);

  // pointer to class FluidEleParameter (access to the general parameter)
  fldparatimint_ = Teuchos::rcp(DRT::ELEMENTS::FluidEleParameterTimInt::Instance(), false);
  // initialize also general parameter list, also it will be overwritten in derived subclasses
  fldpara_ = Teuchos::rcp(DRT::ELEMENTS::FluidEleParameterStd::Instance(), false);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::LocalSolver::
ComputeInteriorResidual(const Teuchos::RCP<MAT::Material> & mat,
                        const std::vector<double>         & val,
                        const std::vector<double>         & accel,
                        const double                        avgPressure,
                        const LINALG::Matrix<nsd_,nen_>   & ebodyforce,
                        const std::vector<double>         & intebodyforce)
{
  zeroMatrix(gRes);
  zeroMatrix(upRes);
  // get constant dynamic viscosity
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
  const double viscosity = actmat->Viscosity();
  const double density = actmat->Density();

  // interpolate the interior values onto quadrature points
  for (unsigned int q=0; q<shapes_.nqpoints_; ++q) {
    // interpolate L_np onto quadrature points
    double velgrad[nsd_][nsd_];
    double acceleration[nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e) {
        velgrad[d][e] = 0.;
        for (unsigned int i=0; i<ndofs_; ++i)
          velgrad[d][e] += shapes_.shfunct(i,q) * val[(d*nsd_+e)*ndofs_+i];
      }
    // interpolate u_np and acceleration
    for (unsigned int d=0; d<nsd_; ++d) {
      double sum = 0.;
      acceleration[d] = 0.;
      for (unsigned int i=0; i<ndofs_; ++i) {
        sum += shapes_.shfunct(i,q) * val[(nsd_*nsd_+d)*ndofs_+i];
        acceleration[d] += shapes_.shfunct(i,q) * accel[(nsd_*nsd_+d)*ndofs_+i];
      }
      velnp(d,q) = sum;
    }
    // interpolate p_np
    double pres = 0.;
    for (unsigned int i=0; i<ndofs_; ++i)
      pres += shapes_.shfunct(i,q) * val[(nsd_*nsd_+nsd_)*ndofs_+i];

    // interpolate body force (currently only ebofoaf_), values from input file
    double force[nsd_];
    for (unsigned int d=0; d<nsd_; ++d) {
      force[d] = 0.;
      for (unsigned int i=0; i<nen_; ++i)
        force[d] += shapes_.funct(i,q) * ebodyforce(d,i);
    }

    // interpolate body force (currently only ebofoaf_), values from forcing vector based on interior dofs
    for (unsigned int d=0; d<nsd_; ++d)
    {
      for (unsigned int i=0; i<ndofs_; ++i)
      {
        force[d] += shapes_.shfunct(i,q) * intebodyforce[(nsd_*nsd_+d)*ndofs_+i];
      }
    }
    // ---------------------------- compute interior residuals
    // residual for L_np
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e) {
        for (unsigned int i=0; i<ndofs_; ++i)
          gRes((d*nsd_+e)*ndofs_+i) -= (velgrad[d][e] * shapes_.shfunct(i,q) +
                                        velnp(d,q) * shapes_.shderxy(i*nsd_+e,q)) * shapes_.jfac(q);
      }
    // residual for u_np
    for (unsigned int d=0; d<nsd_; ++d) {
      double momresd[nsd_];
      if (stokes)
        for (unsigned int e=0; e<nsd_; ++e)
          momresd[e] = -viscosity*(velgrad[d][e]+velgrad[e][d]);
      else
        for (unsigned int e=0; e<nsd_; ++e)
          momresd[e] = -viscosity*(velgrad[d][e]+velgrad[e][d])+density*velnp(d,q)*velnp(e,q);
      momresd[d] += pres;
      if (!stokes)
        force[d] -= density*acceleration[d];
      for (unsigned int i=0; i<ndofs_; ++i) {
        double momder = 0.;
        for (unsigned int e=0; e<nsd_; ++e)
          momder += momresd[e] * shapes_.shderxy(i*nsd_+e,q);
        upRes(d*ndofs_+i) += (momder + force[d]* shapes_.shfunct(i,q)) * shapes_.jfac(q);
      }
    }
    // residual for p_np
    upRes(nsd_*ndofs_) += (avgPressure - pres) * shapes_.jfac(q); // mean pressure
    for (unsigned int i=1; i<ndofs_; ++i) {
      double sum = 0.;
      for (unsigned int d=0; d<nsd_; ++d)
        sum += velnp(d,q) * shapes_.shderxy(i*nsd_+d,q);
      upRes(nsd_*ndofs_+i) += sum * shapes_.jfac(q);
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::LocalSolver::
ComputeInteriorMatrices(const Teuchos::RCP<MAT::Material> &mat,
                        const bool                         evaluateOnlyNonlinear)
{
  const double invtimefac = 1.0/(fldparatimint_->TimeFac());
  if (evaluateOnlyNonlinear && stokes)
    return;

  if (stokes)
    zeroMatrix(uuconv);

  // get constant dynamic viscosity
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
  const double viscosity = actmat->Viscosity();
  const double density = actmat->Density();

  if (!evaluateOnlyNonlinear) {
    zeroMatrix(fgMat);
    zeroMatrix(gfMat);
    zeroMatrix(uuMat);
    zeroMatrix(fuMat);
    zeroMatrix(ufMat);
  }
  else {
    std::memset(fuMat.A(),0,sizeof(double)*fuMat.M()*ndofs_*nsd_); // clear only velocity part
    for (int f=0; f<ufMat.N(); ++f)
      for (unsigned int i=0; i<nsd_*ndofs_; ++i)
        ufMat(i,f) = 0.;
  }

  // loop over interior quadrature points
  for (unsigned int q=0; q<shapes_.nqpoints_; ++q )
  {
    if (!evaluateOnlyNonlinear) {
      for (unsigned int i=0; i<ndofs_; ++i)
        uuMat(nsd_*ndofs_,nsd_*ndofs_+i) += shapes_.shfunct(i,q) * shapes_.jfac(q);
      ufMat(nsd_*ndofs_,0) -= shapes_.jfac(q);
    }

    // now fill the components in the one-sided matrices
    for (unsigned int i=0; i<ndofs_; ++i) {
      // mass matrix part (velocity and velocity gradient use the same mass matrix)
      massPart(i,q) = shapes_.shfunct(i,q);
      const double valf = shapes_.shfunct(i,q) * shapes_.jfac(q);
      massPartW(i,q) = valf;

      // gradient of shape functions
      for (unsigned int d=0; d<nsd_; ++d) {
        if (!evaluateOnlyNonlinear) {
          const double vald = shapes_.shderxy(i*nsd_+d,q);
          gradPart (d*ndofs_+i,q) = vald;
        }

        if (!stokes)
          uPart(d*ndofs_+i,q) = -valf * velnp(d,q) * density;
      }
    }
  }

  // multiply matrices to perform summation over quadrature points
  if (!evaluateOnlyNonlinear) {
    massMat.Multiply('N', 'T', 1., massPart, massPartW, 0.);
    guMat.Multiply('N', 'T', 1., gradPart, massPartW, 0.);
    ugMat = guMat;
    ugMat.Scale(viscosity);
  }
  if (!stokes) {
    uuconv.Multiply('N', 'T', 1., gradPart, uPart, 0.);

    // compute convection: Need to add diagonal part and transpose off-diagonal blocks
    // (same trick as done when eliminating the velocity gradient)
    for (unsigned int i=0; i<ndofs_; ++i)
      for (unsigned int j=0; j<ndofs_; ++j) {
        double sumdiag = 0.;
        for (unsigned int d=0; d<nsd_; ++d) {
          sumdiag += uuconv(d*ndofs_+j,d*ndofs_+i);
          for (unsigned int e=0; e<d; ++e)
            std::swap(uuconv(d*ndofs_+j,e*ndofs_+i),uuconv(e*ndofs_+j,d*ndofs_+i));
        }
        for (unsigned int d=0; d<nsd_; ++d)
          uuconv(d*ndofs_+j,d*ndofs_+i) += sumdiag;
      }
  }

  // merge matrices (do not merge convection matrices into uuMat now but later)
  if (!evaluateOnlyNonlinear) {
    for (unsigned int i=0; i<ndofs_; ++i)
      for (unsigned int j=0; j<ndofs_; ++j) {

        // fill in mass matrix for the velocity
        if (!stokes)
          for (unsigned int d=0; d<nsd_; ++d) {
            uuMat(d*ndofs_+j,d*ndofs_+i) = density * invtimefac * massMat(j,i);
          }

        for (unsigned int d=0; d<nsd_; ++d) {
          // fill in -grad v * pI
          uuMat(d*ndofs_+j,nsd_*ndofs_+i) = -guMat(d*ndofs_+j,i);
          // fill in -u * grad q
          if (j > 0)
            uuMat(nsd_*ndofs_+j,d*ndofs_+i) = -guMat(d*ndofs_+j,i);
        }
      }
    // we want to multiply ugMat by guMat below for which we need to access the
    // entries in guMat in a transposed way
    for (unsigned int i=0; i<ndofs_; ++i)
      for (unsigned int d=0; d<nsd_; ++d)
        for (unsigned int j=0; j<i; ++j)
          std::swap(guMat(d*ndofs_+j,i), guMat(d*ndofs_+i,j));
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::LocalSolver::
ComputeFaceResidual(const int                           face,
                    const Teuchos::RCP<MAT::Material> & mat,
                    const std::vector<double>         & val,
                    const std::vector<double>         & traceval,
                    Epetra_SerialDenseVector          & elevec)
{
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
  const double viscosity = actmat->Viscosity();
  const double density = actmat->Density();

  // compute pressure average on element
  double presavg = 0.;
  for (unsigned int i=0; i<ndofs_; ++i)
    presavg += shapes_.shfunctAvg(i) * val[(nsd_*nsd_+nsd_)*ndofs_+i];

  double velnorm = 0., vol = 0.;
  for (unsigned int q=0; q<shapesface_.nqpoints_; ++q) {
    // interpolate u_n
    for (unsigned int d=0; d<nsd_; ++d) {
      double u_d = 0.;
      for (unsigned int i=0; i<ndofs_; ++i)
        u_d += shapesface_.shfunctI(i,q) * val[(nsd_*nsd_+d)*ndofs_+i];
      velnorm += u_d * u_d * shapesface_.jfac(q);
    }
    vol += shapesface_.jfac(q);
  }
  velnorm = std::sqrt(velnorm/vol);

  const double lengthScale = 1.;
  stabilization[face] = viscosity/lengthScale + (stokes ? 0. : velnorm*density);
  fvelnp.Shape(nsd_,shapesface_.nqpoints_);

  // interpolate the boundary values onto face quadrature points
  for (unsigned int q=0; q<shapesface_.nqpoints_; ++q) {
    // interpolate interior L_np onto face quadrature points
    double velgradnp[nsd_][nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e) {
        velgradnp[d][e] = 0.;
        for (unsigned int i=0; i<ndofs_; ++i)
          velgradnp[d][e] += shapesface_.shfunctI(i,q) *
                             val[(d*nsd_+e)*ndofs_+i];
      }
    // interpolate u_np
    double ifvelnp[nsd_];
    for (unsigned int d=0; d<nsd_; ++d) {
      ifvelnp[d] = 0.;
      for (unsigned int i=0; i<ndofs_; ++i)
        ifvelnp[d] += shapesface_.shfunctI(i,q) * val[(nsd_*nsd_+d)*ndofs_+i];
    }
    // interpolate p_np
    double presnp = 0.;
    for (unsigned int i=0; i<ndofs_; ++i)
      presnp += shapesface_.shfunctI(i,q) * val[(nsd_*nsd_+nsd_)*ndofs_+i];


    // interpolate trace value
    for (unsigned int d=0; d<nsd_; ++d) {
      double sum = 0.;
        for (unsigned int i=0; i<shapesface_.nfdofs_; ++i)
          sum += shapesface_.shfunct(i,q) * traceval[1+face*nsd_*shapesface_.nfdofs_+d*shapesface_.nfdofs_+i];
      fvelnp(d,q) = sum;
    }

    // ---------------------------- compute face residuals
    // residual for L_np
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e) {
        const double res = fvelnp(d,q) * shapesface_.normals(e,q) * shapesface_.jfac(q);
        for (unsigned int i=0; i<ndofs_; ++i)
          gRes((d*nsd_+e)*ndofs_+i) += shapesface_.shfunctI(i,q) * res;
      }

    // residual for u_np
    for (unsigned int d=0; d<nsd_; ++d) {
      double momres[nsd_];
      if (stokes)
        for (unsigned int e=0; e<nsd_; ++e)
          momres[e] = -viscosity*(velgradnp[d][e]+velgradnp[e][d]);
      else
        for (unsigned int e=0; e<nsd_; ++e)
          momres[e] = -viscosity*(velgradnp[d][e]+velgradnp[e][d])+density*fvelnp(d,q)*fvelnp(e,q);

      momres[d] += presnp;
      double res = 0;
      for (unsigned int e=0; e<nsd_; ++e)
        res += momres[e] * shapesface_.normals(e,q);
      res += stabilization[face]*(ifvelnp[d]-fvelnp(d,q));
      res *= shapesface_.jfac(q);
      for (unsigned int i=0; i<ndofs_; ++i)
        upRes(d*ndofs_+i) -= res * shapesface_.shfunctI(i,q);
      res -= (-traceval[0] + presavg) * shapesface_.jfac(q) * shapesface_.normals(d,q);
      for (unsigned int i=0; i<shapesface_.nfdofs_; ++i)
        elevec(1+face*nsd_*shapesface_.nfdofs_+d*shapesface_.nfdofs_+i) -= res * shapesface_.shfunct(i,q);
      elevec(0) -= fvelnp(d,q) * shapesface_.normals(d,q) * shapesface_.jfac(q);
    }

    // residual for p_np
    double presres = 0.;
    for (unsigned int d=0; d<nsd_; ++d)
      presres += fvelnp(d,q)*shapesface_.normals(d,q);
    presres *= shapesface_.jfac(q);
    for (unsigned int i=1; i<ndofs_; ++i) {
      upRes(nsd_*ndofs_+i) -=
          presres*(shapesface_.shfunctI(i,q)-shapes_.shfunctAvg(i));
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::LocalSolver::
ComputeFaceMatrices (const int                          face,
                     const Teuchos::RCP<MAT::Material> &mat,
                     const bool                         evaluateOnlyNonlinear,
                     Epetra_SerialDenseMatrix          &elemat)
{
  // get constant dynamic viscosity
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
  const double viscosity = actmat->Viscosity();
  const double density = actmat->Density();

  trMat.Shape(ndofs_*nsd_,shapesface_.nfdofs_);
  trMatAvg.Shape(ndofs_*nsd_,shapesface_.nfdofs_);

  // TODO write a fast version of this function where we build on the assumption that
  // face basis functions and interior basis functions coincide for certain indices

  // perform face quadrature
  for (unsigned int q=0; q<shapesface_.nqpoints_; ++q) {

    double velNormal = 0.;
    for (unsigned int d=0; d<nsd_; ++d)
      velNormal += shapesface_.normals(d,q) * fvelnp(d,q);
    velNormal *= density;

    double stabvel[nsd_][nsd_];
    for (unsigned int d=0; d<nsd_; ++d) {
      for (unsigned int e=0; e<nsd_; ++e) {
        stabvel[d][e] = 0.;
        if (!stokes)
          stabvel[d][e] += density * fvelnp(d,q) * shapesface_.normals(e,q);
      }
      if (!stokes)
        stabvel[d][d] += velNormal;
      stabvel[d][d] -= stabilization[face];
    }

    const double jac = shapesface_.jfac(q);

    for (unsigned int i=0; i<shapesface_.nfdofs_; ++i) {
      for (unsigned int j=0; j<shapesface_.nfdofs_; ++j) {
        const double shape = shapesface_.shfunct(i,q) * shapesface_.shfunct(j,q) * jac;
        for (unsigned int d=0; d<nsd_; ++d)
          for (unsigned int e=0; e<nsd_; ++e)
            elemat(1+face*nsd_*shapesface_.nfdofs_+shapesface_.nfdofs_*d+j,1+face*nsd_*shapesface_.nfdofs_+shapesface_.nfdofs_*e+i) += shape * stabvel[d][e];
      }

      if (!evaluateOnlyNonlinear)
        for (unsigned int j=0; j<ndofs_; ++j) {
          const double shape = shapesface_.shfunct(i,q) * jac * shapesface_.shfunctI(j,q);
          const double shapeAvg = shapesface_.shfunct(i,q) * jac * (shapesface_.shfunctI(j,q)-shapes_.shfunctAvg(j));
          for (unsigned int d=0; d<nsd_; ++d) {
            trMat(d*ndofs_+j,i) += shape * shapesface_.normals(d,q);
            trMatAvg(d*ndofs_+j,i) += shapeAvg * shapesface_.normals(d,q);
          }
        }

      for (unsigned int j=0; j<ndofs_; ++j) {
        const double shape = shapesface_.shfunct(i,q) * shapesface_.shfunctI(j,q) * jac;
        for (unsigned int d=0; d<nsd_; ++d) {
          for (unsigned int e=0; e<nsd_; ++e) {
            ufMat(d*ndofs_+j,1+face*nsd_*shapesface_.nfdofs_+shapesface_.nfdofs_*e+i) += shape * stabvel[d][e];
          }
          fuMat(1+face*nsd_*shapesface_.nfdofs_+shapesface_.nfdofs_*d+i,d*ndofs_+j) += shape * stabilization[face];
        }
      }

      // -<psi,\lambda * n>
      for (unsigned int d=0; d<nsd_; ++d)
        elemat(1+(face*nsd_+d)*shapesface_.nfdofs_+i,0) += shapesface_.shfunct(i,q) * jac * shapesface_.normals(d,q);
    }

    for (unsigned int i=0; i<ndofs_; ++i)
      for (unsigned int j=0; j<ndofs_; ++j) {
        const double shape = shapesface_.shfunctI(i,q) * shapesface_.shfunctI(j,q) * jac * stabilization[face];
        for (unsigned int d=0; d<nsd_; ++d)
          uuconv(d*ndofs_+i,d*ndofs_+j) += shape;
      }
    if (!evaluateOnlyNonlinear)
      for (unsigned int i=0; i<ndofs_; ++i) {
        for (unsigned int j=0; j<ndofs_; ++j) {
          const double shape = shapesface_.shfunctI(i,q) * shapesface_.shfunctI(j,q) * jac;
          for (unsigned int d=0; d<nsd_; ++d) {
            const double val = shape * shapesface_.normals(d,q);
            ugMat(d*ndofs_+j,i) -= viscosity*val;
            uuMat(d*ndofs_+j,nsd_*ndofs_+i) += val;
          }
        }
      }
  }

  // merge matrices
  if (!evaluateOnlyNonlinear) {
    for (unsigned int i=0; i<shapesface_.nfdofs_; ++i) {
      for (unsigned int j=0; j<ndofs_; ++j) {
        for (unsigned int d=0; d<nsd_; ++d) {
          fuMat(1+face*nsd_*shapesface_.nfdofs_+shapesface_.nfdofs_*d+i,nsd_*ndofs_+j) += trMatAvg(d*ndofs_+j,i);
          ufMat(nsd_*ndofs_+j,1+face*nsd_*shapesface_.nfdofs_+shapesface_.nfdofs_*d+i) += trMatAvg(d*ndofs_+j,i);
        }
        for (unsigned int d=0; d<nsd_; ++d)
          for (unsigned int e=0; e<nsd_; ++e) {
            gfMat((nsd_*d+e)*ndofs_+j,1+face*nsd_*shapesface_.nfdofs_+shapesface_.nfdofs_*d+i) = -trMat(e*ndofs_+j,i);
            fgMat(1+face*nsd_*shapesface_.nfdofs_+shapesface_.nfdofs_*d+i,(nsd_*d+e)*ndofs_+j) -= viscosity*trMat(e*ndofs_+j,i);
            fgMat(1+face*nsd_*shapesface_.nfdofs_+shapesface_.nfdofs_*e+i,(nsd_*d+e)*ndofs_+j) -= viscosity*trMat(d*ndofs_+j,i);
          }
      }
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::LocalSolver::EliminateVelocityGradient(Epetra_SerialDenseMatrix &elemat)
{
  // invert mass matrix. Inverse will be stored in massMat, too
  {
    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(massMat);
    inverseMass.Invert();
  }

  // add contribution of mass matrix to velocity/pressure part
  // create UG * diag(M^{-1}) * GU,

  // compute UG * M^{-1}, store result in tmpMatGrad
  tmpMatGrad.Multiply('N','N',1.0,ugMat,massMat,0.);

  // GU and UG are not fully generated, instead, only three different blocks are kept
  // to compute UG * M^{-1} * GU, therefore compute the product of reduced matrices
  // and fill the values in the
  // local matrix. Since we want to use the symmetric gradient and its block component
  // are exactly in the other order compared to what the big matrix-matrix product does,
  // need to transpose the blocks. Similarly, the Laplacian results in a sum of the
  // diagonal blocks.

  // compute (UG * M^{-1}) * GU
  tmpMat.Multiply('N','T',1.,tmpMatGrad,guMat,0.);
  for (unsigned int i=0; i<ndofs_; ++i)
    for (unsigned int j=0; j<ndofs_; ++j) {
      double diagSum = 0;
      for (unsigned int d=0; d<nsd_; ++d)
        diagSum += tmpMat(d*ndofs_+j,d*ndofs_+i);
      for (unsigned int d=0; d<nsd_; ++d)
        uuMat(d*ndofs_+j,d*ndofs_+i) -= tmpMat(d*ndofs_+j,d*ndofs_+i) + diagSum;
      for (unsigned int d=0; d<nsd_; ++d)
        for (unsigned int e=0; e<nsd_; ++e)
          if (d != e)
            uuMat(d*ndofs_+j,e*ndofs_+i) -= tmpMat(e*ndofs_+j,d*ndofs_+i);
    }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::LocalSolver::
SolveResidual()
{
  for (unsigned int i=0; i<(nsd_+1)*ndofs_; ++i)
    upUpd(i) = upRes(i);

  // compute UG * M^{-1} gRes. Since UG is not stored completely, need some loops.
  // Note: the data field tmpMatGrad contains UG * M^{-1} after EliminateVelocityGradient
  //
  // shape of UG in 3D:
  // [ x y z             ]   [ x     y     z     ]
  // [       x y z       ] + [   x     y     z   ]
  // [             x y z ]   [     x     y     z ]
  // whereas we store the following in tmpMatGrad:
  // [ x ]
  // [ y ]
  // [ z ]
  for (unsigned int d=0; d<nsd_; ++d)
    for (unsigned int i=0; i<ndofs_; ++i) {
      double sum[nsd_];
      for (unsigned int e=0; e<nsd_; ++e)
        sum[e] = 0.;
      for (unsigned int j=0; j<ndofs_; ++j)
        for (unsigned int e=0; e<nsd_; ++e)
          sum[e] += tmpMatGrad(d*ndofs_+i,j) * (gRes((e*nsd_+d)*ndofs_+j) + gRes((d*nsd_+e)*ndofs_+j));
      for (unsigned int e=0; e<nsd_; ++e)
        upUpd(e*ndofs_+i) -= sum[e];
    }

  // merge matrices for uu to get the real Schur complement matrix
  for (unsigned int i=0; i<nsd_*ndofs_; ++i) {
    for (unsigned int j=0; j<nsd_*ndofs_; ++j)
      uuMatFinal(j,i) = uuMat(j,i) + uuconv(j,i);
    for (unsigned int j=nsd_*ndofs_; j<(nsd_+1)*ndofs_; ++j)
      uuMatFinal(j,i) = uuMat(j,i);
  }
  for (unsigned int i=nsd_*ndofs_; i<(nsd_+1)*ndofs_; ++i)
    for (unsigned int j=0; j<(nsd_+1)*ndofs_; ++j)
      uuMatFinal(j,i) = uuMat(j,i);

  // clear first pressure column
  for (unsigned int i=0; i<(nsd_+1)*ndofs_; ++i)
    uuMatFinal(i,nsd_*ndofs_) = uuMatFinal(nsd_*ndofs_,i) = 0.;
  uuMatFinal(nsd_*ndofs_,nsd_*ndofs_) = 1.;

  // factorize uuMatFinal and solve. do not use LINALG::FixedSizeSerialDenseSolver because
  // we want to solve twice and reuse the factorization
  Epetra_LAPACK lapack;
  const int size = uuMatFinal.M();
  pivots.resize(size);
  int errnum;
  lapack.GETRF(size, size, uuMatFinal.A(), size, &(pivots[0]), &errnum);
  if (errnum > 0){
    uuMatFinal.Print(std::cout); uuMat.Print(std::cout);
  }
  dsassert(errnum == 0, "Factorization failed");
  lapack.GETRS('N', size, 1, uuMatFinal.A(), size, &(pivots[0]), upUpd.A(), size, &errnum);
  dsassert(errnum == 0, "Substitution failed");

  // compute Rg - GU * upUpd
  // shape of GU in 3D
  // [ x     ]
  // [ y     ]
  // [ z     ]
  // [   x   ]
  // [   y   ]
  // [   z   ]
  // [     x ]
  // [     y ]
  // [     z ]
  for (unsigned int d=0; d<nsd_; ++d)
    for (unsigned int i=0; i<ndofs_; ++i) {
      double sum[nsd_];
      for (unsigned int e=0; e<nsd_; ++e)
        sum[e] = 0;
      for (unsigned int j=0; j<ndofs_; ++j)
        for (unsigned int e=0; e<nsd_; ++e)
          sum[e] += guMat(d*ndofs_+j,i) * upUpd(e*ndofs_+j);
      for (unsigned int e=0; e<nsd_; ++e)
        gRes((e*nsd_+d)*ndofs_+i) -= sum[e];
    }

  // compute M^{-1} * Rg
  for (unsigned int i=0; i<ndofs_; ++i) {
    double sum[nsd_*nsd_];
    for (unsigned int e=0; e<nsd_*nsd_; ++e)
      sum[e] = 0.;
    for (unsigned int j=0; j<ndofs_; ++j)
      for (unsigned int e=0; e<nsd_*nsd_; ++e)
        sum[e] += massMat(j,i) * gRes(e*ndofs_+j); // use symmetry for faster matrix access
    for (unsigned int e=0; e<nsd_*nsd_; ++e)
      gUpd(e*ndofs_+i) = sum[e];
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::LocalSolver::
CondenseLocalPart(Epetra_SerialDenseMatrix &eleMat,
                  Epetra_SerialDenseVector &eleVec)
{
  for (unsigned int i=0; i<nfaces_*nsd_*shapesface_.nfdofs_; ++i)
    eleMat(0,1+i) = eleMat(1+i,0);

  // first get residual to obtain first part of condensed residual vector,
  // which will also compute and factorize uuMatFinal
  SolveResidual();

  // clear first column of pressure
  for (unsigned int i=1; i<1+nfaces_*nsd_*shapesface_.nfdofs_; ++i)
    fuMat(i,ndofs_*nsd_) = 0;

  // compute residual vector: need to multiply residual by fuMat and fgMat
  for (unsigned int i=1; i<1+nfaces_*nsd_*shapesface_.nfdofs_; ++i) {
    double sum = 0.;
    for (unsigned int j=0; j<ndofs_*(nsd_+1); ++j)
      sum += fuMat(i,j) * upUpd(j);
    eleVec(i) -= sum;
    sum = 0.;
    for (unsigned int j=0; j<ndofs_*nsd_*nsd_; ++j)
      sum += fgMat(i,j) * gUpd(j);
    eleVec(i) -= sum;
  }

  Epetra_BLAS blas;

  for (unsigned int f=1; f<1+nfaces_*shapesface_.nfdofs_*nsd_; ++f) {
    // gfMat is block-structured similarly to GU, so only use non-zero entries
    const unsigned cindex = ((f-1)/shapesface_.nfdofs_)%nsd_;
    // compute (UG * M^{-1}) * GF = tmpMatGrad * GF
    // shape of UG in 3D:
    // [ x y z             ]   [ x     y     z     ]
    // [       x y z       ] + [   x     y     z   ]
    // [             x y z ]   [     x     y     z ]
    const double *tmpPtr = tmpMatGrad.A();
    for (unsigned int i=0; i<ndofs_; ++i) {
      double sum1 = 0;
      for (unsigned int e=0; e<nsd_; ++e) {
        double sum2 = 0;
        for (unsigned int j=0; j<ndofs_; ++j) {
          sum1 += tmpPtr[i+e*ndofs_+j*nsd_*ndofs_] * gfMat((cindex*nsd_+e)*ndofs_+j,f);
          sum2 += tmpPtr[i+cindex*ndofs_+j*nsd_*ndofs_] * gfMat((cindex*nsd_+e)*ndofs_+j,f);
        }
        ufMat(e*ndofs_+i,f) -= sum2;
      }
      ufMat(cindex*ndofs_+i,f) -= sum1;
    }
  }

  // solve for velocity matrix
  Epetra_LAPACK lapack;
  int errnum;
  dsassert(pivots.size()==static_cast<unsigned int>(uuMatFinal.M()) && pivots[0]+pivots[1] > 0,
           "Matrix seems to not have been factorized");
  lapack.GETRS('N', uuMatFinal.M(), ufMat.N(), uuMatFinal.A(), uuMatFinal.M(), &(pivots[0]), ufMat.A(),
               ufMat.M(), &errnum);
  dsassert(errnum == 0, "Substitution failed");

  // put velocity/pressure part into element matrix
  blas.GEMM('N','N', fuMat.M(),ufMat.N(),fuMat.N(), -1., fuMat.A(),fuMat.M(), ufMat.A(),ufMat.M(),
            1., eleMat.A(), eleMat.M());

  // update gfMat and apply inverse mass matrix: GF <- M^{-1} (GF - GU * UF)
  for (unsigned int f=1; f<1+nfaces_*shapesface_.nfdofs_*nsd_; ++f) {
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int i=0; i<ndofs_; ++i) {
        double sum[nsd_];
        for (unsigned int e=0; e<nsd_; ++e)
          sum[e] = 0;
        for (unsigned int j=0; j<ndofs_; ++j)
          for (unsigned int e=0; e<nsd_; ++e)
            sum[e] += guMat(d*ndofs_+j,i) * ufMat(e*ndofs_+j,f); // note special structure of guMat (transposed)
        for (unsigned int e=0; e<nsd_; ++e)
          gfMat((e*nsd_+d)*ndofs_+i,f) -= sum[e];
      }
    // apply M^{-1}, store temporary result
    for (unsigned int i=0; i<ndofs_; ++i) {
      double sum[nsd_*nsd_];
      for (unsigned int e=0; e<nsd_*nsd_; ++e)
        sum[e] = 0.;
      for (unsigned int j=0; j<ndofs_; ++j)
        for (unsigned int e=0; e<nsd_*nsd_; ++e)
          sum[e] += massMat(j,i) * gfMat(e*ndofs_+j,f); // use symmetry for faster matrix access
      for (unsigned int e=0; e<nsd_*nsd_; ++e)
        gRes(e*ndofs_+i) = sum[e];
    }
    for (unsigned int i=0; i<ndofs_*nsd_*nsd_; ++i)
      gfMat(i,f) = gRes(i);
  }

  // compute FG * (M^{-1} GF)
  blas.GEMM('N','N', fgMat.M(),gfMat.N(),fgMat.N(), -1., fgMat.A(),fgMat.M(), gfMat.A(),gfMat.M(),
            1., eleMat.A(), eleMat.M());

}



// explicit instantiation of template classes
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::wedge15>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::nurbs27>;

