/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_hdg.cpp

\brief main file containing routines for calculation of HDG fluid element

*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc_hdg.H"
#include "fluid_ele_calc.H"
#include "fluid_ele_parameter.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_parameter_timint.H"
#include "fluid_ele_action.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_mat/newtonianfluid.H"



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
DRT::ELEMENTS::FluidEleCalcHDG<distype>::FluidEleCalcHDG():
    localSolver_(shapes_)
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
  const bool updateLocally = params.get<bool>("needslocalupdate");

  shapes_.Evaluate(*ele);

  ebofoaf_.PutScalar(0);
  eprescpgaf_.PutScalar(0);
  escabofoaf_.PutScalar(0);
  FluidEleCalc<distype>::BodyForce(ele, localSolver_.fldparatimint_->Time(),
      localSolver_.fldpara_->PhysicalType(), ebofoaf_, eprescpgaf_, escabofoaf_);

  ReadGlobalVectors(*ele, discretization, lm, updateLocally);

  // solves the local problem of the nonlinear iteration before
  if (updateLocally) {
    localSolver_.ComputeInteriorResidual(mat, interiorVal_, interiorAcc_, traceVal_[0], ebofoaf_);
    localSolver_.ComputeInteriorMatrices(mat, false);

    dsassert(nfaces_ == static_cast<unsigned int>(ele->NumFace()), "Internal error");

    // loop over faces
    for (unsigned int f=0; f<nfaces_; ++f) {
      shapes_.EvaluateFace(*ele, f);
      localSolver_.ComputeFaceResidual(f, mat, interiorVal_, traceVal_, elevec1);
      localSolver_.ComputeFaceMatrices(f, mat, false, elemat1);
    }

    localSolver_.EliminateVelocityGradient(elemat1);
    localSolver_.SolveResidual();
    UpdateSecondarySolution(*ele, discretization, localSolver_.gUpd, localSolver_.upUpd);
  }

  zeroMatrix(elemat1);
  zeroMatrix(elevec1);
  localSolver_.ComputeInteriorResidual(mat, interiorVal_, interiorAcc_, traceVal_[0], ebofoaf_);
  localSolver_.ComputeInteriorMatrices(mat, updateLocally);
  for (unsigned int f=0; f<nfaces_; ++f) {
    shapes_.EvaluateFace(*ele, f);
    localSolver_.ComputeFaceResidual(f, mat, interiorVal_, traceVal_, elevec1);
    localSolver_.ComputeFaceMatrices(f, mat, updateLocally, elemat1);
  }

  if (!updateLocally)
    localSolver_.EliminateVelocityGradient(elemat1);

  localSolver_.CondenseLocalPart(elemat1,elevec1);

  elevec1.Scale(1./localSolver_.fldparatimint_->AlphaF());

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
  traceVal_.resize(1+nfaces_*nsd_*nfdofs_);
  interiorVal_.resize(((nsd_+1)*nsd_+1)*ndofs_);
  interiorAcc_.resize(((nsd_+1)*nsd_+1)*ndofs_);
  dsassert(lm.size() == traceVal_.size(), "Internal error");
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState("velaf");
  DRT::UTILS::ExtractMyValues(*matrix_state,traceVal_,lm);

  // read the interior values from solution vector
  matrix_state = discretization.GetState(1,"intvelaf");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  DRT::UTILS::ExtractMyValues(*matrix_state,interiorVal_,localDofs);

  matrix_state = discretization.GetState(1,"intaccam");
  localDofs = discretization.Dof(1, &ele);
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
    double update = i<nsd_*nsd_*ndofs_ ? updateG(i) : updateUp(i-nsd_*nsd_*ndofs_);

    const double valfac = 1./localSolver_.fldparatimint_->AlphaF();
    secondary[lid] += update * valfac;

    // write the update back into the local vectors (when doing local update,
    // we do not re-read from the global vectors)
    interiorVal_[i] += update;

    const double accfac = localSolver_.fldparatimint_->AlphaM() * valfac /
        (localSolver_.fldparatimint_->Dt() * localSolver_.fldparatimint_->Gamma());
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
  shapes_.Evaluate(*ele);
  const double time = localSolver_.fldparatimint_->Time();

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

  double err_u = 0., err_p = 0., err_h = 0., norm_u = 0., norm_p = 0., norm_h = 0.;
  for (unsigned int q=0; q<ndofs_; ++q)
  {
    double numericalGrad[nsd_][nsd_];
    double numerical[nsd_+1];
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e) {
        numericalGrad[d][e] = 0;
        for (unsigned int i=0; i<ndofs_; ++i)
          numericalGrad[d][e] += shapes_.shfunct(i,q) * vecValues[(d*nsd_+e)*ndofs_+i];
      }
    for (unsigned int d=0; d<=nsd_; ++d) {
      numerical[d] = 0.;
      for (unsigned int i=0; i<ndofs_; ++i)
        numerical[d] += shapes_.shfunct(i,q) * vecValues[(nsd_*nsd_+d)*ndofs_+i];
    }
    for (unsigned int d=0; d<nsd_; ++d)
      xyz(d) = shapes_.xyzreal(d,q);

    FluidEleCalc<distype>::EvaluateAnalyticSolutionPoint(xyz, time, calcerr, mat, u, p, dervel);

    for (unsigned int d=0; d<nsd_; ++d)
      err_u += (u(d) - numerical[d]) * (u(d) - numerical[d]) * shapes_.jfac(q);
    err_p += (p - numerical[nsd_]) * (p - numerical[nsd_]) * shapes_.jfac(q);
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e)
        err_h += (dervel(d,e) - numericalGrad[d][e]) * (dervel(d,e) - numericalGrad[d][e]) * shapes_.jfac(q);
    for (unsigned int d=0; d<nsd_; ++d)
      norm_u += u(d) * u(d) * shapes_.jfac(q);
    norm_p += p * p * shapes_.jfac(q);
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e)
        norm_h += dervel(e,d) * dervel(e,d) * shapes_.jfac(q);
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

  shapes_.Evaluate(*ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 ||
           elevec2.M() == (nsd_*nsd_+nsd_+1)*ndofs_, "Wrong size in project vector 2");

  // get function
  const int *initfield = params.getPtr<int>("initfield");
  const int *start_func = params.getPtr<int>("startfuncno");

  double avgpre = 0., vol = 0.;
  if (elevec2.M() > 0) {
    LINALG::Matrix<ndofs_,nsd_*nsd_+nsd_+1> localMat (elevec2.A(),true);
    localMat.PutScalar(0.);

    // create mass matrix for interior by looping over quadrature points
    for (unsigned int q=0; q<ndofs_; ++q )
    {
      const double fac = shapes_.jfac(q);
      const double sqrtfac = std::sqrt(fac);
      double xyz[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        xyz[d] = shapes_.xyzreal(d,q);
      double u[nsd_];
      double grad[nsd_][nsd_];
      double p;
      dsassert(initfield != NULL && start_func != NULL,
               "initfield or startfuncno not set for initial value");
      EvaluateAll(*start_func, INPAR::FLUID::InitialField(*initfield), xyz, u, grad, p);

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i=0; i<ndofs_; ++i) {
        // mass matrix part
        localSolver_.massPart(i,q) = shapes_.shfunct(i,q) * sqrtfac;

        for (unsigned int d=0; d<nsd_; ++d)
          for (unsigned int e=0; e<nsd_; ++e)
            localMat(i,d*nsd_+e) += shapes_.shfunct(i,q) * grad[d][e] * fac;
        for (unsigned int d=0; d<nsd_; ++d)
          localMat(i,nsd_*nsd_+d) += shapes_.shfunct(i,q) * u[d] * fac;
        localMat(i,nsd_*nsd_+nsd_) += shapes_.shfunct(i,q) * p * fac;
      }

      avgpre += p * fac;
      vol += fac;
    }
    localSolver_.massMat.Multiply('N', 'T', 1., localSolver_.massPart,localSolver_.massPart, 0.);

    // solve mass matrix system, return values in localMat = elevec2 correctly ordered
    {
      LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_,nsd_*nsd_+nsd_+1> inverseMass;
      LINALG::Matrix<ndofs_,ndofs_> mass(localSolver_.massMat,true);
      inverseMass.SetMatrix(mass);
      inverseMass.SetVectors(localMat,localMat);
      inverseMass.Solve();
    }
  }

  LINALG::Matrix<nfdofs_,nfdofs_> mass;
  LINALG::Matrix<nfdofs_,nsd_> trVec;
  dsassert(elevec1.M() == nsd_*nfdofs_ ||
           elevec1.M() == 1+nfaces_*nsd_*nfdofs_, "Wrong size in project vector 1");

  const unsigned int *faceConsider = params.getPtr<unsigned int>("faceconsider");
  Teuchos::Array<int> *functno = params.getPtr<Teuchos::Array<int> >("funct");
  Teuchos::Array<int> *onoff = params.getPtr<Teuchos::Array<int> >("onoff");
  double *time = params.getPtr<double>("time");

  for (unsigned int face=0; face<nfaces_; ++face) {
    // check whether we are in the project phase for all faces or for boundary values
    if (initfield == NULL) {
      dsassert(faceConsider != NULL, "Unsupported operation");
      if (*faceConsider != face)
        continue;
    }
    shapes_.EvaluateFace(*ele, face);
    mass.PutScalar(0.);
    trVec.PutScalar(0.);

    for (unsigned int q=0; q<nfdofs_; ++q) {
      const double fac = shapes_.jfacF(q);
      double xyz[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        xyz[d] = shapes_.xyzFreal(d,q);
      double u[nsd_];
      if (initfield != NULL)
        EvaluateVelocity(*start_func, INPAR::FLUID::InitialField(*initfield), xyz, u);
      else {
        dsassert(functno != NULL && time != NULL && onoff != NULL,
                 "No array with functions given");
        for (unsigned int d=0; d<nsd_; ++d) {
          if ((*onoff)[d] == 0)
            continue;
          const int funct_num = (*functno)[d];
          if (funct_num > 0)
            u[d] = DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(d, xyz, *time,
                                                                         &discretization);
        }
      }

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i=0; i<nfdofs_; ++i) {
        // mass matrix
        for (unsigned int j=0; j<nfdofs_; ++j)
          mass(i,j) += shapes_.shfunctF(i,q) * shapes_.shfunctF(j,q) * fac;

        for (unsigned int d=0; d<nsd_; ++d)
          trVec(i,d) += shapes_.shfunctF(i,q) * u[d] * fac;
      }
    }

    LINALG::FixedSizeSerialDenseSolver<nfdofs_,nfdofs_,nsd_> inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec,trVec);
    inverseMass.Solve();

    if (initfield != NULL)
      for (unsigned int d=0; d<nsd_; ++d)
        for (unsigned int i=0; i<nfdofs_; ++i)
          elevec1(1+face*nfdofs_*nsd_+d*nfdofs_+i) = trVec(i,d);
    else
      for (unsigned int d=0; d<nsd_; ++d)
        for (unsigned int i=0; i<nfdofs_; ++i)
          elevec1(d*nfdofs_+i) = trVec(i,d);
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
  dsassert(elevec1.M() == (int)nen_*(2*nsd_+1)+1, "Vector does not have correct size");
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  LINALG::Matrix<1,ndofs_> values;

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
      shapes_.xsi(idim) = locations(idim,i);
    shapes_.polySpace_.Evaluate(shapes_.xsi,values);

    // compute values for velocity and pressure by summing over all basis functions
    for (unsigned int d=0; d<=nsd_; ++d) {
      double sum = 0;
      for (unsigned int k=0; k<ndofs_; ++k)
        sum += values(k) * solvalues[(nsd_*nsd_+d)*ndofs_+k];
      elevec1(d*nen_+i) = sum;
    }
  }

  // get trace solution values
  locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace
                (DRT::UTILS::DisTypeToFaceShapeType<distype>::shape);
  matrix_state = discretization.GetState(0,"velnp");
  localDofs = discretization.Dof(0, ele);
  solvalues.resize (localDofs.size());
  for (unsigned int i=0; i<solvalues.size(); ++i) {
    const int lid = matrix_state->Map().LID(localDofs[i]);
    solvalues[i] = (*matrix_state)[lid];
  }

  LINALG::Matrix<1,nfdofs_> fvalues;
  for (unsigned int f=0; f<nfaces_; ++f)
    for (int i=0; i<DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i) {
      // evaluate shape polynomials in node
      for (unsigned int idim=0;idim<nsd_-1;idim++)
        shapes_.xsiF(idim) = locations(idim,i);
      shapes_.polySpaceFace_.Evaluate(shapes_.xsiF,fvalues); // TODO: fix face orientation here

      // compute values for velocity and pressure by summing over all basis functions
      for (unsigned int d=0; d<nsd_; ++d) {
        double sum = 0;
        for (unsigned int k=0; k<nfdofs_; ++k)
          sum += fvalues(k) * solvalues[1+f*nsd_*nfdofs_+d*nfdofs_+k];
        elevec1((nsd_+1+d)*nen_+shapes_.faceNodeOrder[f][i]) = sum;
      }
  }
  elevec1((2*nsd_+1)*nen_) = solvalues[0];


  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::EvaluateVelocity(const int start_func,
    const INPAR::FLUID::InitialField initfield,
    const double (&xyz)[nsd_],
    double (&u)[nsd_]) const
{
  // pass on dummy entries (costs a little but will not be significant)
  double grad[nsd_][nsd_];
  double p;
  EvaluateAll(start_func, initfield, xyz, u, grad, p);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcHDG<distype>::EvaluateAll(const int start_func,
    const INPAR::FLUID::InitialField initfield,
    const double (&xyz)[nsd_],
    double (&u)[nsd_],
    double (&grad)[nsd_][nsd_],
    double  &p) const
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
    u[0] = -a * ( exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) +
                  exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) );
    u[1] = -a * ( exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) +
                  exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) );
    u[2] = -a * ( exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) +
                  exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) );

    grad[0][0] = -a * ( a * exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) -
                        a * exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) );
    grad[0][1] = -a * ( a * exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) -
                        d * exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) );
    grad[0][2] = -a * ( d * exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) +
                        a * exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) );
    grad[1][0] = -a * ( d * exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) +
                        a * exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) );
    grad[1][1] = -a * ( a * exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) -
                        a * exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) );
    grad[1][2] = -a * ( a * exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) -
                        d * exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) );
    grad[2][0] = -a * ( a * exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) -
                        d * exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) );
    grad[2][1] = -a * ( d * exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) +
                        a * exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) );
    grad[2][2] = -a * ( a * exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) -
                        a * exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) );

    p = -a*a/2.0 *
      ( exp(2.0*a*xyz[0])
        + exp(2.0*a*xyz[1])
        + exp(2.0*a*xyz[2])
        + 2.0 * sin(a*xyz[0] + d*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) * exp(a*(xyz[1]+xyz[2]))
        + 2.0 * sin(a*xyz[1] + d*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) * exp(a*(xyz[2]+xyz[0]))
        + 2.0 * sin(a*xyz[2] + d*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) * exp(a*(xyz[0]+xyz[1]))
        );
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
DRT::ELEMENTS::FluidEleCalcHDG<distype>::ShapeValues::
ShapeValues()
:
polySpace_(DRT::ELEMENTS::FluidHDG::degree),
polySpaceFace_(DRT::ELEMENTS::FluidHDG::degree)
{
  dsassert(polySpace_.Size() == ndofs_, "Wrong polynomial space constructed");
  dsassert(polySpaceFace_.Size() == nfdofs_, "Wrong polynomial space constructed");

  const int degree = DRT::ELEMENTS::FluidHDG::degree;
  LINALG::Matrix<1,ndofs_> values;
  LINALG::Matrix<nsd_,ndofs_> derivs;
  LINALG::Matrix<1,nfdofs_> faceValues;

  quadrature_ = DRT::UTILS::GaussPointCache::Instance().Create(distype, degree*2);
  const DRT::Element::DiscretizationType facedis = DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;
  fquadrature_ = DRT::UTILS::GaussPointCache::Instance().Create(facedis, degree*2);

  shfunct.Shape(ndofs_,ndofs_);
  shfunctAvg.Resize(ndofs_);
  shderiv.Shape(ndofs_*nsd_,ndofs_);
  shderxy.Shape(ndofs_*nsd_,ndofs_);
  jfac.Resize(ndofs_);

  dsassert(static_cast<unsigned int>(quadrature_->NumPoints())==ndofs_, "Internal error - not implemented");
  for (unsigned int q=0; q<ndofs_; ++q ) {
    // gauss point in real coordinates
    const double* gpcoord = quadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_;idim++)
      xsi(idim) = gpcoord[idim];

    polySpace_.Evaluate(xsi,values);
    polySpace_.Evaluate_deriv1(xsi,derivs);

    for (unsigned int i=0; i<ndofs_; ++i) {
      shfunct(i,q) = values(i);
      for (unsigned int d=0; d<nsd_; ++d)
        shderiv(i*nsd_+d,q) = derivs(d,i);
    }

    LINALG::Matrix<nen_,1> myfunct(funct.A()+q*nen_,true);
    DRT::UTILS::shape_function<distype>(xsi,myfunct);
  }

  shfunctFNoPermute.Shape(nfdofs_, nfdofs_);
  shfunctF.Shape(nfdofs_, nfdofs_);
  shfunctI.resize(nfaces_);
  for (unsigned int f=0;f<nfaces_; ++f)
    shfunctI[f].Shape(ndofs_, nfdofs_);
  normals.Shape(nsd_, nfdofs_);
  jfacF.Resize(nfdofs_);

  for (unsigned int q=0; q<nfdofs_; ++q ) {
    const double* gpcoord = fquadrature_->Point(q);

    const unsigned int codim = nsd_-1;
    for (unsigned int idim=0;idim<codim;idim++)
      xsiF(idim) = gpcoord[idim];

    polySpaceFace_.Evaluate(xsiF,faceValues);
    for (unsigned int i=0; i<nfdofs_; ++i)
      shfunctFNoPermute(i,q) = faceValues(i);

    LINALG::Matrix<nfn_,1> myfunct(functF.A()+q*nfn_,true);
    DRT::UTILS::shape_function<facedis>(xsiF,myfunct);
  }

  Epetra_SerialDenseMatrix quadrature(nfdofs_,nsd_,false);
  Epetra_SerialDenseMatrix trafo(nsd_,nsd_,false);
  for (unsigned int f=0; f<nfaces_; ++f) {
    DRT::UTILS::BoundaryGPToParentGP<nsd_>(quadrature,trafo,*fquadrature_,distype,
                                           DRT::UTILS::DisTypeToFaceShapeType<distype>::shape, f);
    for (unsigned int q=0; q<nfdofs_; ++q) {
      for (unsigned int d=0; d<nsd_; ++d)
        xsi(d) = quadrature(q,d);
      polySpace_.Evaluate(xsi,values);
      for (unsigned int i=0; i<ndofs_; ++i)
        shfunctI[f](i,q) = values(i);
    }
  }

  if (nsd_ == 2)
    faceNodeOrder = DRT::UTILS::getEleNodeNumberingLines(distype);
  else if (nsd_ == 3)
    faceNodeOrder = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
  else
    dserror("Not implemented for dim != 2, 3");
}



template <DRT::Element::DiscretizationType distype>
void
DRT::ELEMENTS::FluidEleCalcHDG<distype>::ShapeValues::Evaluate (const DRT::Element &ele)
{
  dsassert(ele.Shape() == distype, "Internal error");
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(&ele,xyze);

  for (unsigned int i=0; i<ndofs_; ++i)
    shfunctAvg(i) = 0.;
  double faceVol = 0.;

  // evaluate geometry
  for (unsigned int q=0; q<ndofs_; ++q) {
    const double* gpcoord = quadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_;idim++)
      xsi(idim) = gpcoord[idim];

    DRT::UTILS::shape_function_deriv1<distype>(xsi,deriv);
    xjm.MultiplyNT(deriv,xyze);
    jfac(q) = xji.Invert(xjm) * quadrature_->Weight(q);

    LINALG::Matrix<nen_,1> myfunct(funct.A()+q*nen_,true);
    LINALG::Matrix<nsd_,1> mypoint(xyzreal.A()+q*nsd_,true);
    mypoint.MultiplyNN(xyze,myfunct);

    // transform shape functions
    for (unsigned int i=0; i<ndofs_; ++i)
      for (unsigned int d=0; d<nsd_; ++d) {
        shderxy(i*nsd_+d,q) = xji(d,0) * shderiv(i*nsd_,q);
        for (unsigned int e=1; e<nsd_; ++e)
          shderxy(i*nsd_+d,q) += xji(d,e) * shderiv(i*nsd_+e,q);
      }

    for (unsigned int i=0; i<ndofs_; ++i)
      shfunctAvg(i) += shfunct(i,q) * jfac(q);
    faceVol += jfac(q);
  }
  faceVol = 1./faceVol;
  for (unsigned int i=0; i<ndofs_; ++i)
    shfunctAvg(i) *= faceVol;
}



template <DRT::Element::DiscretizationType distype>
void
DRT::ELEMENTS::FluidEleCalcHDG<distype>::ShapeValues::
EvaluateFace (const DRT::Element &ele,
              const unsigned int  face)
{
  const DRT::Element::DiscretizationType facedis = DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;

  // get face position array from element position array
  dsassert(faceNodeOrder[face].size() == nfn_,
           "Internal error");
  for (unsigned int i=0; i<nfn_; ++i)
    for (unsigned int d=0; d<nsd_; ++d)
      xyzeF(d,i) = xyze(d,faceNodeOrder[face][i]);

  // evaluate geometry
  for (unsigned int q=0; q<nfdofs_; ++q) {
    const double* gpcoord = fquadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_-1;idim++)
      xsiF(idim) = gpcoord[idim];

    DRT::UTILS::shape_function_deriv1<facedis>(xsiF,derivF);
    double jacdet = 0;
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<facedis>(xyzeF,derivF,metricTensor,jacdet,&normal);
    for (unsigned int d=0; d<nsd_; ++d)
      normals(d,q) = normal(d);
    jfacF(q) = jacdet * fquadrature_->Weight(q);

    LINALG::Matrix<nfn_,1> myfunct(functF.A()+q*nfn_,true);
    LINALG::Matrix<nsd_,1> mypoint(xyzFreal.A()+q*nsd_,true);
    mypoint.MultiplyNN(xyzeF,myfunct);
  }

  // figure out how to permute face indices by checking permutation of nodes
  const int * nodeIds = ele.NodeIds();
  const int * fnodeIds = ele.Faces()[face]->NodeIds();
  const int ndofs1d = DRT::ELEMENTS::FluidHDG::degree+1;
  // easy case: standard orientation
  bool standard = true;
  for (unsigned int i=0; i<nfn_; ++i)
    if (nodeIds[faceNodeOrder[face][i]] != fnodeIds[i])
      standard = false;
  if (standard) {
    //std::cout << "standard orientation" << std::endl;
    for (unsigned int q=0; q<nfdofs_; ++q)
      for (unsigned int i=0; i<nfdofs_; ++i)
        shfunctF(i,q) = shfunctFNoPermute(i,q);
  }
  // OK, the orientation is different from what I expect. see if we can find it
  else switch (nsd_)
  {
  case 2:
    // face flipped is the only case
    {
      dsassert(nodeIds[faceNodeOrder[face][1]] == fnodeIds[0] &&
               nodeIds[faceNodeOrder[face][0]] == fnodeIds[1], "Unknown face orientation in 2D");
      //std::cout << "flipped case: ";
      for (unsigned int i=0; i<nfdofs_; ++i) {
        //std::cout << nfdofs_-1-i << " ";
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(nfdofs_-1-i,q);
        }
      //std::cout << std::endl;
    }
    break;
  case 3:
    dsassert(distype == DRT::Element::hex8 ||
             distype == DRT::Element::hex20 ||
             distype == DRT::Element::hex27 ||
             distype == DRT::Element::nurbs8 ||
             distype == DRT::Element::nurbs27,
             "Not implemented for given shape");
    if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[1] &&
        nodeIds[faceNodeOrder[face][1]] == fnodeIds[0] &&
        nodeIds[faceNodeOrder[face][2]] == fnodeIds[3] &&
        nodeIds[faceNodeOrder[face][3]] == fnodeIds[2])    // x-direction mirrored
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = ndofs1d-1-ax + ay * ndofs1d;
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[0] &&
             nodeIds[faceNodeOrder[face][1]] == fnodeIds[3] &&
             nodeIds[faceNodeOrder[face][2]] == fnodeIds[2] &&
             nodeIds[faceNodeOrder[face][3]] == fnodeIds[1])    // permute x and y
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = ay + ax * ndofs1d;
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[3] &&
             nodeIds[faceNodeOrder[face][1]] == fnodeIds[2] &&
             nodeIds[faceNodeOrder[face][2]] == fnodeIds[1] &&
             nodeIds[faceNodeOrder[face][3]] == fnodeIds[0])    // y mirrored
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = ax + (ndofs1d-1-ay) * ndofs1d;
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[2] &&
             nodeIds[faceNodeOrder[face][1]] == fnodeIds[3] &&
             nodeIds[faceNodeOrder[face][2]] == fnodeIds[0] &&
             nodeIds[faceNodeOrder[face][3]] == fnodeIds[1])    // x and y mirrored
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay) * ndofs1d;
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[2] &&
             nodeIds[faceNodeOrder[face][1]] == fnodeIds[1] &&
             nodeIds[faceNodeOrder[face][2]] == fnodeIds[0] &&
             nodeIds[faceNodeOrder[face][3]] == fnodeIds[3])    // x and y mirrored and permuted
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = (ndofs1d-1-ay) + (ndofs1d-1-ax) * ndofs1d; // note that this is for lexicographic ordering
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else
    {
      for (unsigned int i=0; i<4; ++i)
        std::cout << nodeIds[faceNodeOrder[face][i]] << " " << fnodeIds[i] << "   ";
      std::cout << std::endl << std::flush;
      dserror("Unknown face orientation in 3D");
    }
    break;
  default:
    dserror("Only implemented in 2D and 3D");
    break;
  }
}



template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcHDG<distype>::LocalSolver::
LocalSolver(const ShapeValues &shapeValues)
:
stokes (false),
shapes_(shapeValues)
{
  uuMat.Shape((nsd_+1)*ndofs_,(nsd_+1)*ndofs_);
  uuMatFinal.Shape((nsd_+1)*ndofs_,(nsd_+1)*ndofs_);
  guMat.Shape(nsd_*ndofs_,ndofs_);
  ugMat.Shape(nsd_*ndofs_,ndofs_);

  gfMat.Shape(nsd_*nsd_*ndofs_,1+nfaces_*nsd_*nfdofs_);
  fgMat.Shape(gfMat.N(), gfMat.M());
  ufMat.Shape((nsd_+1)*ndofs_,1+nfaces_*nsd_*nfdofs_);
  fuMat.Shape(ufMat.N(), ufMat.M());

  massPart.Shape(ndofs_,ndofs_);
  gradPart.Shape(nsd_*ndofs_,ndofs_);
  uPart.Shape(ndofs_*nsd_,ndofs_);

  massMat.Shape(ndofs_,ndofs_);
  uuconv.Shape(ndofs_*nsd_,ndofs_*nsd_);
  tmpMat.Shape(ndofs_*nsd_,ndofs_*nsd_);
  tmpMatGrad.Shape(nsd_*ndofs_,ndofs_);

  trMat.Shape(ndofs_*nsd_,nfdofs_);
  trMatAvg.Shape(ndofs_*nsd_,nfdofs_);

  velnp.Shape(nsd_,ndofs_);
  fvelnp.Shape(nsd_,nfdofs_);
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
                        const LINALG::Matrix<nsd_,nen_>   & ebodyforce)
{
  zeroMatrix(gRes);
  zeroMatrix(upRes);
  // get constant dynamic viscosity
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
  const double viscosity = actmat->Viscosity();
  const double density = actmat->Density();

  // interpolate the interior values onto quadrature points
  for (unsigned int q=0; q<ndofs_; ++q) {
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

    // interpolate body force (currently only ebofoaf_)
    double force[nsd_];
    for (unsigned int d=0; d<nsd_; ++d) {
      force[d] = 0.;
      for (unsigned int i=0; i<nen_; ++i)
        force[d] += shapes_.funct(i,q) * ebodyforce(d,i);
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
  for (unsigned int q=0; q<ndofs_; ++q )
  {
    const double sqrtfac = std::sqrt(shapes_.jfac(q));

    // fix pressure average on element
    if (!evaluateOnlyNonlinear) {
      for (unsigned int i=0; i<ndofs_; ++i)
        uuMat(nsd_*ndofs_,nsd_*ndofs_+i) += shapes_.shfunct(i,q) * shapes_.jfac(q);
      ufMat(nsd_*ndofs_,0) -= shapes_.jfac(q);
    }

    // now fill the components in the one-sided matrices
    for (unsigned int i=0; i<ndofs_; ++i) {
      // mass matrix part (velocity and velocity gradient use the same mass matrix)
      const double valf = shapes_.shfunct(i,q) * sqrtfac;
      massPart(i,q) = valf;

      // gradient of shape functions
      for (unsigned int d=0; d<nsd_; ++d) {
        if (!evaluateOnlyNonlinear) {
          const double vald = shapes_.shderxy(i*nsd_+d,q) * sqrtfac;
          gradPart (d*ndofs_+i,q) = vald;
        }

        if (!stokes)
          uPart(d*ndofs_+i,q) = -valf * velnp(d,q) * density;
      }
    }
  }

  // multiply matrices to perform summation over quadrature points
  if (!evaluateOnlyNonlinear) {
    massMat.Multiply('N', 'T', 1., massPart, massPart, 0.);
    guMat.Multiply('N', 'T', 1., gradPart, massPart, 0.);
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
  for (unsigned int q=0; q<nfdofs_; ++q) {
    // interpolate u_n
    for (unsigned int d=0; d<nsd_; ++d) {
      double u_d = 0.;
      for (unsigned int i=0; i<ndofs_; ++i)
        u_d += shapes_.shfunctI[face](i,q) * val[(nsd_*nsd_+d)*ndofs_+i];
      velnorm += u_d * u_d * shapes_.jfacF(q);
    }
    vol += shapes_.jfacF(q);
  }
  velnorm = std::sqrt(velnorm/vol);

  const double lengthScale = 1.;
  stabilization[face] = viscosity/lengthScale + (stokes ? 0. : velnorm*density);

  // interpolate the boundary values onto face quadrature points
  for (unsigned int q=0; q<nfdofs_; ++q) {
    // interpolate interior L_np onto face quadrature points
    double velgradnp[nsd_][nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e) {
        velgradnp[d][e] = 0.;
        for (unsigned int i=0; i<ndofs_; ++i)
          velgradnp[d][e] += shapes_.shfunctI[face](i,q) *
                             val[(d*nsd_+e)*ndofs_+i];
      }
    // interpolate u_np
    double ifvelnp[nsd_];
    for (unsigned int d=0; d<nsd_; ++d) {
      ifvelnp[d] = 0.;
      for (unsigned int i=0; i<ndofs_; ++i)
        ifvelnp[d] += shapes_.shfunctI[face](i,q) * val[(nsd_*nsd_+d)*ndofs_+i];
    }
    // interpolate p_np
    double presnp = 0.;
    for (unsigned int i=0; i<ndofs_; ++i)
      presnp += shapes_.shfunctI[face](i,q) * val[(nsd_*nsd_+nsd_)*ndofs_+i];

    // interpolate trace value
    for (unsigned int d=0; d<nsd_; ++d) {
      double sum = 0.;
        for (unsigned int i=0; i<nfdofs_; ++i)
          sum += shapes_.shfunctF(i,q) * traceval[1+face*nsd_*nfdofs_+d*nfdofs_+i];
      fvelnp(d,q) = sum;
    }

    // ---------------------------- compute face residuals
    // residual for L_np
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int e=0; e<nsd_; ++e) {
        const double res = fvelnp(d,q) * shapes_.normals(e,q) * shapes_.jfacF(q);
        for (unsigned int i=0; i<ndofs_; ++i)
          gRes((d*nsd_+e)*ndofs_+i) += shapes_.shfunctI[face](i,q) * res;
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
        res += momres[e] * shapes_.normals(e,q);
      res += stabilization[face]*(ifvelnp[d]-fvelnp(d,q));
      res *= shapes_.jfacF(q);
      for (unsigned int i=0; i<ndofs_; ++i)
        upRes(d*ndofs_+i) -= res * shapes_.shfunctI[face](i,q);
      res -= (-traceval[0] + presavg) * shapes_.jfacF(q) * shapes_.normals(d,q);
      for (unsigned int i=0; i<nfdofs_; ++i)
        elevec(1+face*nsd_*nfdofs_+d*nfdofs_+i) -= res * shapes_.shfunctF(i,q);
      elevec(0) -= fvelnp(d,q) * shapes_.normals(d,q) * shapes_.jfacF(q);
    }

    // residual for p_np
    double presres = 0.;
    for (unsigned int d=0; d<nsd_; ++d)
      presres += fvelnp(d,q)*shapes_.normals(d,q);
    presres *= shapes_.jfacF(q);
    for (unsigned int i=1; i<ndofs_; ++i) {
      upRes(nsd_*ndofs_+i) -=
          presres*(shapes_.shfunctI[face](i,q)-shapes_.shfunctAvg(i));
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

  if (!evaluateOnlyNonlinear) {
    zeroMatrix(trMat);
    zeroMatrix(trMatAvg);
  }

  // TODO write a fast version of this function where we build on the assumption that
  // face basis functions and interior basis functions coincide for certain indices

  // perform face quadrature
  for (unsigned int q=0; q<nfdofs_; ++q) {

    double velNormal = 0.;
    for (unsigned int d=0; d<nsd_; ++d)
      velNormal += shapes_.normals(d,q) * fvelnp(d,q);
    velNormal *= density;

    double stabvel[nsd_][nsd_];
    for (unsigned int d=0; d<nsd_; ++d) {
      for (unsigned int e=0; e<nsd_; ++e) {
        stabvel[d][e] = 0.;
        if (!stokes)
          stabvel[d][e] += density * fvelnp(d,q) * shapes_.normals(e,q);
      }
      if (!stokes)
        stabvel[d][d] += velNormal;
      stabvel[d][d] -= stabilization[face];
    }

    const double jac = shapes_.jfacF(q);

    for (unsigned int i=0; i<nfdofs_; ++i) {
      for (unsigned int j=0; j<nfdofs_; ++j) {
        const double shape = shapes_.shfunctF(i,q) * shapes_.shfunctF(j,q) * jac;
        for (unsigned int d=0; d<nsd_; ++d)
          for (unsigned int e=0; e<nsd_; ++e)
            elemat(1+face*nsd_*nfdofs_+nfdofs_*d+j,1+face*nsd_*nfdofs_+nfdofs_*e+i) += shape * stabvel[d][e];
      }

      if (!evaluateOnlyNonlinear)
        for (unsigned int j=0; j<ndofs_; ++j) {
          const double shape = shapes_.shfunctF(i,q) * jac * shapes_.shfunctI[face](j,q);
          const double shapeAvg = shapes_.shfunctF(i,q) * jac * (shapes_.shfunctI[face](j,q)-shapes_.shfunctAvg(j));
          for (unsigned int d=0; d<nsd_; ++d) {
            trMat(d*ndofs_+j,i) += shape * shapes_.normals(d,q);
            trMatAvg(d*ndofs_+j,i) += shapeAvg * shapes_.normals(d,q);
          }
        }

      for (unsigned int j=0; j<ndofs_; ++j) {
        const double shape = shapes_.shfunctF(i,q) * shapes_.shfunctI[face](j,q) * jac;
        for (unsigned int d=0; d<nsd_; ++d) {
          for (unsigned int e=0; e<nsd_; ++e) {
            ufMat(d*ndofs_+j,1+face*nsd_*nfdofs_+nfdofs_*e+i) += shape * stabvel[d][e];
          }
          fuMat(1+face*nsd_*nfdofs_+nfdofs_*d+i,d*ndofs_+j) += shape * stabilization[face];
        }
      }

      // -<psi,\lambda * n>
      for (unsigned int d=0; d<nsd_; ++d)
        elemat(1+(face*nsd_+d)*nfdofs_+i,0) += shapes_.shfunctF(i,q) * jac * shapes_.normals(d,q);
    }

    for (unsigned int i=0; i<ndofs_; ++i)
      for (unsigned int j=0; j<ndofs_; ++j) {
        const double shape = shapes_.shfunctI[face](i,q) * shapes_.shfunctI[face](j,q) * jac * stabilization[face];
        for (unsigned int d=0; d<nsd_; ++d)
          uuconv(d*ndofs_+i,d*ndofs_+j) += shape;
      }
    if (!evaluateOnlyNonlinear)
      for (unsigned int i=0; i<ndofs_; ++i) {
        for (unsigned int j=0; j<ndofs_; ++j) {
          const double shape = shapes_.shfunctI[face](i,q) * shapes_.shfunctI[face](j,q) * jac;
          for (unsigned int d=0; d<nsd_; ++d) {
            const double val = shape * shapes_.normals(d,q);
            ugMat(d*ndofs_+j,i) -= viscosity*val;
            uuMat(d*ndofs_+j,nsd_*ndofs_+i) += val;
          }
        }
      }
  }

  // merge matrices
  if (!evaluateOnlyNonlinear) {
    for (unsigned int i=0; i<nfdofs_; ++i) {
      for (unsigned int j=0; j<ndofs_; ++j) {
        for (unsigned int d=0; d<nsd_; ++d) {
          fuMat(1+face*nsd_*nfdofs_+nfdofs_*d+i,nsd_*ndofs_+j) += trMatAvg(d*ndofs_+j,i);
          ufMat(nsd_*ndofs_+j,1+face*nsd_*nfdofs_+nfdofs_*d+i) += trMatAvg(d*ndofs_+j,i);
        }
        for (unsigned int d=0; d<nsd_; ++d)
          for (unsigned int e=0; e<nsd_; ++e) {
            gfMat((nsd_*d+e)*ndofs_+j,1+face*nsd_*nfdofs_+nfdofs_*d+i) = -trMat(e*ndofs_+j,i);
            fgMat(1+face*nsd_*nfdofs_+nfdofs_*d+i,(nsd_*d+e)*ndofs_+j) -= viscosity*trMat(e*ndofs_+j,i);
            fgMat(1+face*nsd_*nfdofs_+nfdofs_*e+i,(nsd_*d+e)*ndofs_+j) -= viscosity*trMat(d*ndofs_+j,i);
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
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseMass;
    LINALG::Matrix<ndofs_,ndofs_> mass(massMat,true);
    inverseMass.SetMatrix(mass);
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
  for (unsigned int i=0; i<nfaces_*nsd_*nfdofs_; ++i)
    eleMat(0,1+i) = eleMat(1+i,0);

  // first get residual to obtain first part of condensed residual vector,
  // which will also compute and factorize uuMatFinal
  SolveResidual();

  // clear first column of pressure
  for (unsigned int i=1; i<1+nfaces_*nsd_*nfdofs_; ++i)
    fuMat(i,ndofs_*nsd_) = 0;

  // compute residual vector: need to multiply residual by fuMat and fgMat
  for (unsigned int i=1; i<1+nfaces_*nsd_*nfdofs_; ++i) {
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

  for (unsigned int f=1; f<1+nfaces_*nfdofs_*nsd_; ++f) {
    // gfMat is block-structured similarly to GU, so only use non-zero entries
    const unsigned cindex = ((f-1)/nfdofs_)%nsd_;
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
  for (unsigned int f=1; f<1+nfaces_*nfdofs_*nsd_; ++f) {
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



// template classes
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::nurbs27>;

