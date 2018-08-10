/*--------------------------------------------------------------------------*/
/*!
\file acou_ele_calc.cpp
\brief all functionality for acoustic element evaluations

<pre>
\level 2

\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
 */
/*--------------------------------------------------------------------------*/

#include "acou_ele_calc.H"
#include "acou_ele_action.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_geometry/position_array.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_mat/acoustic.H"

#include <Epetra_SerialDenseSolver.h>
#include <Teuchos_TimeMonitor.hpp>


namespace
{
  void zeroMatrix(Epetra_SerialDenseMatrix &mat)
  {
    std::memset(mat.A(), 0, sizeof(double) * mat.M() * mat.N());
  }

  void reshapeMatrixIfNecessary(Epetra_SerialDenseMatrix &matrix, const int nrows,const int ncols)
  {
    if (nrows != matrix.M() || ncols != matrix.N())
      matrix.Shape(nrows, ncols);
  }
}

/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype>::AcouEleCalc()
{}

/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou* ele,
    DRT::Discretization & discretization, const std::vector<int> & lm,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material> & mat,
    Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra,
    const DRT::UTILS::GaussIntegration &, bool offdiag)
{
  return this->Evaluate(ele, discretization, lm, params, mat, elemat1_epetra,
                        elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra, offdiag);
}

/*----------------------------------------------------------------------*
 * Evaluate
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou* ele,
    DRT::Discretization & discretization, const std::vector<int> & lm,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material> & mat,
    Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix&,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector&, bool offdiag)
{
  // check if this is an hdg element and init completepoly
  if (const DRT::ELEMENTS::Acou * hdgele = dynamic_cast<const DRT::ELEMENTS::Acou*>(ele))
    usescompletepoly_ = hdgele->UsesCompletePolynomialSpace();
  else
    dserror("cannot cast element to acou element");
  const ACOU::Action action = DRT::INPUT::get<ACOU::Action>(params, "action");
  InitializeShapes(ele);

  bool updateonly = false;
  shapes_->Evaluate(*ele);

  switch (action)
  {
  case ACOU::project_field:
  {
    if (mat->MaterialType() != INPAR::MAT::m_acousticmat)
      dserror("for physical type 'lossless' please supply MAT_Acoustic");
    ElementInit(ele, params);

    //{ only for ADER TRI projection
    double dt = 1.0;
    localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");
    //}

    localSolver_->ProjectField(ele, params, elevec1,elevec2);
    break;
  }
  case ACOU::project_dirich_field:
  {
    if (mat->MaterialType() != INPAR::MAT::m_acousticmat)
      dserror("for physical type 'lossless' please supply MAT_Acoustic");
    ElementInit(ele, params);
    localSolver_->ProjectDirichField(ele, params, elevec1);
    break;
  }
  case ACOU::project_optical_field:
  {
    if (mat->MaterialType() != INPAR::MAT::m_acousticmat)
      dserror("for physical type 'lossless' please supply MAT_Acoustic");
    localSolver_->ProjectOpticalField(ele, params, elevec2);
    break;
  }
  case ACOU::ele_init:
  {
    ElementInit(ele, params);
    break;
  }
  case ACOU::fill_restart_vecs:
  {
    bool padapty = params.get<bool>("padaptivity");
    ReadGlobalVectors(ele, discretization, lm, padapty);
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
    const bool padapty = params.get<bool>("padaptivity");
    ReadGlobalVectors(ele, discretization, lm, padapty);
    NodeBasedValues(ele, elevec1, padapty);
    break;
  }
  case ACOU::interpolate_psi_to_node:
  {
    const bool padapty = params.get<bool>("padaptivity");
    double dt = params.get<double>("dt");
    const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
    double c = actmat->SpeedofSound(ele->Id());
    double rho = actmat->Density(ele->Id());
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");
    if(dyna_!=INPAR::ACOU::acou_impleuler) // explicit time integration does not need the scaling with c!
    {
      c = 1.0;
      rho = 1.0;
    }
    ReadGlobalVectors(ele, discretization, lm, padapty);
    NodeBasedPsi(elevec1, dt, c*c*rho);
    break;
  }
  case ACOU::calc_acou_error:
  {
    const bool padapty = params.get<bool>("padaptivity");
    double dt = params.get<double>("dt");
    const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
    double rho = actmat->Density(ele->Id());

    localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);

    ReadGlobalVectors(ele, discretization, lm, padapty);
    ComputeError(ele, params, elevec1,rho,dt);
    break;
  }
  case ACOU::calc_abc:
  {
    int face = params.get<int>("face");
    int sumindex = 0;
    for (int i = 0; i < face; ++i)
    {
      DRT::UTILS::PolynomialSpaceParams params(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,ele->Faces()[i]->Degree(), usescompletepoly_);
      int nfdofs = DRT::UTILS::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(params)->Size();
      sumindex += nfdofs;
    }

    if (!params.isParameter("nodeindices"))
      localSolver_->ComputeAbsorbingBC(discretization, ele, params, mat, face, elemat1, sumindex);
    else
      dserror("why would you set an absorbing LINE in THREE dimensions?");

    break;
  }
  case ACOU::bd_integrate:
  {
    int face = params.get<int>("face");
    localSolver_->ComputeBoundaryIntegral(ele, params, face);

    break;
  }
  case ACOU::calc_systemmat_and_residual:
  {
    const bool resonly = params.get<bool>("resonly");
    const bool padapty = params.get<bool>("padaptivity");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(ele, discretization, lm, padapty);
    zeroMatrix(elevec1);
    localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);

    if (!resonly)
      localSolver_->CondenseLocalPart(elemat1);

    localSolver_->ComputeResidual(params, elevec1, interiorVelnp_, interiorPressnp_);

    break;
  }
  case ACOU::update_secondary_solution:
    updateonly = true; // no break here!!!
  case ACOU::update_secondary_solution_and_calc_residual:
  {
    bool errormaps = params.get<bool>("errormaps");
    const bool padapty = params.get<bool>("padaptivity");
    const bool allelesequal = params.get<bool>("allelesequal");

    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(ele, discretization, lm, padapty);

    zeroMatrix(elevec1);
    if(!allelesequal)
      localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);

    // this happens for DIRK time integration and is disadvantageous, since history variables are shaped according to elevec size in UpdateInteriorVariablesAndComputeResidual
    if (unsigned(elevec1.M()) != lm.size())
      elevec1.Size(lm.size());

    UpdateInteriorVariablesAndComputeResidual(params, *ele, mat, elevec1, dt, errormaps, updateonly);

    break;
  }
  case ACOU::calc_average_pressure:
  {
    const bool padapty = params.get<bool>("padaptivity");
    ReadGlobalVectors(ele, discretization, lm, padapty);
    ComputePressureAverage(elevec1);
    break;
  }
  case ACOU::calc_ader_sysmat_and_residual:
  {
    //const bool resonly = params.get<bool>("resonly");
    const bool padapty = params.get<bool>("padaptivity");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(ele, discretization, lm, padapty);
    zeroMatrix(elevec1);
    localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);
    elemat1 = localSolver_->Gmat;

    const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
    double rho = actmat->Density(ele->Id());
    double c = actmat->SpeedofSound(ele->Id());
    localSolver_->ComputeADERResidual(ele, discretization, params, elevec1, interiorVelnp_, interiorPressnp_,c,rho);

    break;
  }
  case ACOU::update_ader_solution:
  {
    //bool errormaps = params.get<bool>("errormaps");
    const bool padapty = params.get<bool>("padaptivity");
    const bool allelesequal = params.get<bool>("allelesequal");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(ele, discretization, lm, padapty);

    zeroMatrix(elevec1);
    if(!allelesequal)
      localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);

    UpdateADERSolution(discretization,params, ele, mat, elevec1, dt);

    break;
  }
  case ACOU::prepare_ader_postprocessing:
  {
    //bool errormaps = params.get<bool>("errormaps");
    const bool padapty = params.get<bool>("padaptivity");
    const bool allelesequal = params.get<bool>("allelesequal");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(ele, discretization, lm, padapty);

    zeroMatrix(elevec1);
    if(!allelesequal)
      localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);

    elevec1.Multiply('N','N',-1.0,localSolver_->Imat,interiorVelnp_,0.0);
    elevec1.Multiply('N','N',-1.0,localSolver_->Jmat,interiorPressnp_,1.0);

    break;
  }
  case ACOU::ader_postprocessing:
  {
    //bool errormaps = params.get<bool>("errormaps");
    const bool padapty = params.get<bool>("padaptivity");
    const bool allelesequal = params.get<bool>("allelesequal");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");
    const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
    double rho = actmat->Density(ele->Id());

    ReadGlobalVectors(ele, discretization, lm, padapty);

    zeroMatrix(elevec1);
    if(!allelesequal)
      localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);

    ComputeADERPostProcessing(ele, params, elevec1,rho,dt);

    break;
  }
  case ACOU::ader_postpro_gradanddiv:
  {
    //bool errormaps = params.get<bool>("errormaps");
    const bool padapty = params.get<bool>("padaptivity");
    const bool allelesequal = params.get<bool>("allelesequal");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");
    const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
    double rho = actmat->Density(ele->Id());
    double c = actmat->SpeedofSound(ele->Id());
    DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(ele);

    ReadGlobalVectors(ele, discretization, lm, padapty);

    zeroMatrix(elevec1);
    if(!allelesequal)
      localSolver_->ComputeMatrices(discretization, mat, *ele, dt, dyna_);

    int onfdofs = elevec1.M();
    Epetra_SerialDenseVector traceVal_SDV(onfdofs);
    for (int i = 0; i < onfdofs; ++i)
      traceVal_SDV(i) = traceVal_[i];

    Epetra_SerialDenseVector temp(shapes_->ndofs_ * nsd_);
    temp.Multiply('N', 'N', 1.0, localSolver_->Bmat, interiorPressnp_, 0.0);
    temp.Multiply('N', 'N', 1.0, localSolver_->Cmat, traceVal_SDV, 1.0);

    // improved gradient
    acouele->eleADERimprovedGrad_.Multiply('N', 'N', rho / dt, localSolver_->invAmat, temp, 0.0);

    temp.Resize(shapes_->ndofs_);
    temp.Multiply('N','N',1.0,localSolver_->Hmat, interiorVelnp_,0.0);
    temp.Multiply('N','N',1.0,localSolver_->Dmat, interiorPressnp_,1.0);
    temp.Multiply('N','N',1.0,localSolver_->Emat, traceVal_SDV,1.0);

    Epetra_SerialDenseMatrix invMmat(shapes_->ndofs_,shapes_->ndofs_);
    invMmat += localSolver_->Mmat;
    {
      Epetra_SerialDenseSolver inverseMmat;
      inverseMmat.SetMatrix(invMmat);
      inverseMmat.Invert();
    }
    //improved divergence
    acouele->eleADERimprovedDiv_.Multiply('N','N',1.0/c/c/rho/dt,invMmat,temp,0.0);

    break;
  }
  case ACOU::get_gauss_points:
  {
    int rows = shapes_->xyzreal.M();
    int cols = shapes_->xyzreal.N();
    elemat1.Shape(rows,cols);

    for(int r=0; r<rows; ++r)
      for(int c=0; c<cols; ++c)
        elemat1(r,c) = shapes_->xyzreal(r,c);

    break;
  }
  default:
    dserror("unknown action supplied");
    break;
  } // switch(action)

  return 0;
}

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::InitializeShapes(
    const DRT::ELEMENTS::Acou* ele)
{
  if (shapes_ == Teuchos::null)
    shapes_ = Teuchos::rcp(
        new DRT::UTILS::ShapeValues<distype>(ele->Degree(), usescompletepoly_,2 * ele->Degree()));
  else if (shapes_->degree_ != unsigned(ele->Degree()) || shapes_->usescompletepoly_ != usescompletepoly_)
    shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(ele->Degree(), usescompletepoly_,2 * ele->Degree()));

  if (localSolver_ == Teuchos::null)
    localSolver_ = Teuchos::rcp(new LocalSolver(*shapes_,dyna_));
  else if (localSolver_->ndofs_ != shapes_->ndofs_)
    localSolver_ = Teuchos::rcp(new LocalSolver(*shapes_,dyna_));

  localSolver_->FaceSpecificConstruction(ele, usescompletepoly_);
}

/*----------------------------------------------------------------------*
 * ReadGlobalVectors
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ReadGlobalVectors(DRT::Element * ele,
    DRT::Discretization & discretization, const std::vector<int> & lm,
    const bool padaptivity)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ReadGlobalVectors");
  DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(ele);

  // read vectors from element storage
  reshapeMatrixIfNecessary(interiorVelnp_, acouele->eleinteriorVelnp_.M(), 1);
  reshapeMatrixIfNecessary(interiorPressnp_, acouele->eleinteriorPressnp_.M(),1);

  interiorVelnp_ = acouele->eleinteriorVelnp_;
  interiorPressnp_ = acouele->eleinteriorPressnp_;
  if(trac_with_pml)
    interiorauxiliaryPML_ = acouele->eleinteriorAuxiliaryPML_;

  // read vectors from time integrator
  if (discretization.HasState("trace")) // in case of "update interior variables"
  {
    traceVal_.resize(lm.size());
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState("trace");
    DRT::UTILS::ExtractMyValues(*matrix_state, traceVal_, lm);
  }

  return;
} // ReadGlobalVectors

/*----------------------------------------------------------------------*
 * FillRestartVectors
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::FillRestartVectors(DRT::Element * ele,
    DRT::Discretization & discretization)
{
  // sort this back to the interior values vector
  int size = shapes_->ndofs_ * (nsd_ + 1);
  if(trac_with_pml)
    size = shapes_->ndofs_ * (nsd_ + nsd_ + 1);

  std::vector<double> interiorValnp(size);
  for (unsigned int i = 0; i < interiorValnp.size(); ++i)
  {
    if(!trac_with_pml)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        interiorValnp[i] = interiorPressnp_((i + 1) / (nsd_ + 1) - 1);
      else
      {
        int xyz = i % (nsd_ + 1); // 0 for x, 1 for y and 2 for z (for 3D)
        interiorValnp[i] = interiorVelnp_(xyz * shapes_->ndofs_ + i / (nsd_ + 1));
      }
    }
    else //if(trac_with_pml)
    {
      if ((i + 1) % (nsd_ + nsd_ + 1) == 0)
        interiorValnp[i] = interiorPressnp_((i + 1) / (nsd_ +nsd_ + 1) - 1);
      else if ((i + 1) % (nsd_ + nsd_ + 1) <= nsd_)
      {
        int xyz = i % (nsd_ +nsd_ + 1); // 0 for x, 1 for y and 2 for z (for 3D)
        interiorValnp[i] = interiorVelnp_(xyz * shapes_->ndofs_ + i / (nsd_+nsd_ + 1));
      }
      else
      {
        int xyz = i % (nsd_ +nsd_ + 1) - nsd_; // 0 for x, 1 for y and 2 for z (for 3D)
        interiorValnp[i] = interiorauxiliaryPML_(xyz * shapes_->ndofs_ + i / (nsd_ +nsd_ + 1));
      }
    }
  }

  // tell this change in the interior variables the discretization
  std::vector<int> localDofs = discretization.Dof(1, ele);
  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnp");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for (unsigned int i = 0; i < localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 * ElementInitFromRestart
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ElementInitFromRestart(
    DRT::Element * ele, DRT::Discretization & discretization)
{
  DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(ele);
  int size = shapes_->ndofs_ * (nsd_ + 1);
  if(trac_with_pml)
    size = shapes_->ndofs_ * (nsd_ + nsd_ + 1);
  std::vector<double> interiorValnp(size);

  Teuchos::RCP<const Epetra_Vector> intvel = discretization.GetState(1,"intvelnp");
  std::vector<int> localDofs1 = discretization.Dof(1, ele);
  DRT::UTILS::ExtractMyValues(*intvel, interiorValnp, localDofs1);

  // now write this in corresponding interiorVelnp_ and interiorPressnp_
  for (unsigned int i = 0; i < interiorValnp.size(); ++i)
  {
    if(!trac_with_pml)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        acouele->eleinteriorPressnp_((i + 1) / (nsd_ + 1) - 1) = interiorValnp[i];
      else
      {
        int xyz = i % (nsd_ + 1);
        acouele->eleinteriorVelnp_(xyz * shapes_->ndofs_ + i / (nsd_ + 1)) = interiorValnp[i];
      }
    }
    else
    {
      if ((i + 1) % (nsd_ + nsd_ + 1) == 0)
        acouele->eleinteriorPressnp_((i + 1) / (nsd_ +nsd_ + 1) - 1) = interiorValnp[i];
      else if((i + 1) % (nsd_ + nsd_ + 1) <= nsd_)
      {
        int xyz = i % (nsd_+nsd_ + 1);
        acouele->eleinteriorVelnp_(xyz * shapes_->ndofs_ + i / (nsd_+nsd_ + 1)) = interiorValnp[i];
      }
      else
      {
        int xyz = i % (nsd_ +nsd_ + 1) - nsd_; // 0 for x, 1 for y and 2 for z (for 3D)
        acouele->eleinteriorAuxiliaryPML_(xyz * shapes_->ndofs_ + i / (nsd_ +nsd_ + 1)) =  interiorValnp[i];
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 * ComputeError
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ComputeError(DRT::ELEMENTS::Acou* ele,
    Teuchos::ParameterList& params,
    Epetra_SerialDenseVector& elevec, double rho, double dt)
{
  // for the calculation of the error, we use a higher integration rule
  Teuchos::RCP<DRT::UTILS::GaussPoints> highquad = DRT::UTILS::GaussPointCache::Instance().Create(distype,(ele->Degree() +2) * 2);
  LINALG::Matrix<nsd_, 1> xsi;
  Epetra_SerialDenseVector values(shapes_->ndofs_);
  LINALG::Matrix<nsd_, nen_> deriv;
  LINALG::Matrix<nsd_, nsd_> xjm;

  double time = params.get<double>("time");

  // get function
  int funcno = params.get<int>("funct");
  funcno--;

  double err_p = 0.0,norm_p = 0.0;
  double numerical = 0.0;
  double exact = 0.0;

  for(int q=0; q<highquad->NumPoints(); ++q)
  {
    const double* gpcoord = highquad->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++)
      xsi(idim) = gpcoord[idim];
    shapes_->polySpace_->Evaluate(xsi,values);

    DRT::UTILS::shape_function_deriv1<distype>(xsi,deriv);
    xjm.MultiplyNT(deriv,shapes_->xyze);
    double highjfac = xjm.Determinant() * highquad->Weight(q);

    numerical = 0.0;
    exact = 0.0;
    for (unsigned int i=0; i<shapes_->ndofs_; ++i)
      numerical += values(i) * interiorPressnp_(i);

    LINALG::Matrix<nen_,1>  myfunct;
    DRT::UTILS::shape_function<distype>(xsi,myfunct);
    LINALG::Matrix<nsd_,1> xyzmat;
    xyzmat.MultiplyNN(shapes_->xyze,myfunct);

    double xyz[nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
      xyz[d]=xyzmat(d,0);

    exact = DRT::Problem::Instance()->Funct(funcno).Evaluate(0,xyz,time);

    err_p += ( exact - numerical ) * ( exact - numerical ) * highjfac;
    norm_p += exact * exact * highjfac;
  }
  elevec[0] += err_p;
  elevec[1] += norm_p;

  return;
}

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ComputeADERPostProcessing(DRT::ELEMENTS::Acou* ele,
    Teuchos::ParameterList& params,
    Epetra_SerialDenseVector& elevec, double rho, double dt)
{

  // for the calculation of the error, we use a higher integration rule
  Teuchos::RCP<DRT::UTILS::GaussPoints> highquad = DRT::UTILS::GaussPointCache::Instance().Create(distype,(ele->Degree() +2) * 2);
  LINALG::Matrix<nsd_, 1> xsi;

  double time = params.get<double>("time");

  // get function
  int funcno = params.get<int>("funct");
  funcno--;

  // 1. step: compute improved gradient (bold p in paper by Nguyen)
  int onfdofs = elevec.M();
  Epetra_SerialDenseVector traceVal_SDV(onfdofs);
  for (int i = 0; i < onfdofs; ++i)
    traceVal_SDV(i) = traceVal_[i];

  Epetra_SerialDenseVector temp(shapes_->ndofs_ * nsd_);
  temp.Multiply('N', 'N', 1.0, localSolver_->Bmat, interiorPressnp_, 0.0);
  temp.Multiply('N', 'N', 1.0, localSolver_->Cmat, traceVal_SDV, 1.0);
  Epetra_SerialDenseVector p(shapes_->ndofs_ * nsd_);
  p.Multiply('N', 'N', rho / dt, localSolver_->invAmat, temp, 0.0);

  // 2. step: compute improved pressure field
  DRT::UTILS::PolynomialSpace<nsd_> postpoly(distype, ele->Degree() + 1, ele->UsesCompletePolynomialSpace());
//  int ndofspost = 1;
//  for (unsigned int i = 0; i < nsd_; ++i)
//    ndofspost *= (ele->Degree() + 2);
  int ndofspost = postpoly.Size();

  Epetra_SerialDenseMatrix h(ndofspost, ndofspost);
  Epetra_SerialDenseVector rhs(ndofspost);

  Epetra_SerialDenseMatrix derivs(nsd_, ndofspost);
  Epetra_SerialDenseVector values(ndofspost);

  Epetra_SerialDenseVector myvalues(shapes_->ndofs_);

  LINALG::Matrix<nsd_, nen_> deriv;
  LINALG::Matrix<nsd_, nsd_> xjm, xji;
  LINALG::Matrix<nsd_, nen_> xyze;
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_> >(ele, xyze);

  for (int q = 0; q < highquad->NumPoints(); ++q)
  {
    const double* gpcoord = highquad->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++)
      xsi(idim) = gpcoord[idim];

    postpoly.Evaluate(xsi, values);
    postpoly.Evaluate_deriv1(xsi, derivs);

    DRT::UTILS::shape_function_deriv1<distype>(xsi, deriv);
    xjm.MultiplyNT(deriv, xyze);
    const double jfac = xji.Invert(xjm) * highquad->Weight(q);

    // transform shape functions derivatives
    for (int i = 0; i < ndofspost; ++i)
    {
      double res[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        res[d] = xji(d, 0) * derivs(0, i);
        for (unsigned int e = 1; e < nsd_; ++e)
          res[d] += xji(d, e) * derivs(e, i);
      }
      for (unsigned int d = 0; d < nsd_; ++d)
        derivs(d, i) = res[d];
    }

    shapes_->polySpace_->Evaluate(xsi, myvalues);

    for (int j = 0; j < ndofspost; ++j)
      h(0, j) += values(j) * jfac;
    for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
      rhs(0) += myvalues(j) * jfac * interiorPressnp_[j];

    for (int i = 1; i < ndofspost; ++i)
    {
      for (int j = 0; j < ndofspost; ++j)
      {
        double t = 0;
        for (unsigned int d = 0; d < nsd_; ++d)
          t += derivs(d, i) * derivs(d, j);
        h(i, j) += t * jfac;
      }
    }
    double ugrad[nsd_]= {0.0};

    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
        ugrad[d] += myvalues(j) * p(j + d * shapes_->ndofs_);
    for (int i = 1; i < ndofspost; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        rhs(i) += ugrad[d] * derivs(d, i) * jfac;

  } // for (int q=0; q<postquad->NumPoints(); ++q)

  {
    Epetra_SerialDenseSolver inverseH;
    inverseH.SetMatrix(h);
    inverseH.SetVectors(rhs, rhs);
    inverseH.Solve();
  }
  // now, rhs contains the dof values of the superconvergent pressure solution

  // compute error in this field:
  double err_p_post = 0.0;
  double numerical_post = 0.0;
  double exact = 0.0;
  for(int q=0; q<highquad->NumPoints(); ++q)
  {
    const double* gpcoord = highquad->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++)
      xsi(idim) = gpcoord[idim];
    postpoly.Evaluate(xsi, values);

    DRT::UTILS::shape_function_deriv1<distype>(xsi,deriv);
    xjm.MultiplyNT(deriv,shapes_->xyze);
    double highjfac = xjm.Determinant() * highquad->Weight(q);

    numerical_post = 0.0;
    exact = 0.0;
    for (int i = 0; i < ndofspost; ++i)
      numerical_post += values(i) * rhs(i);

    LINALG::Matrix<nen_,1>  myfunct;
    DRT::UTILS::shape_function<distype>(xsi,myfunct);
    LINALG::Matrix<nsd_,1> xyzmat;
    xyzmat.MultiplyNN(shapes_->xyze,myfunct);

    double xyz[nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
      xyz[d]=xyzmat(d,0);

    exact = DRT::Problem::Instance()->Funct(funcno).Evaluate(0,xyz,time);

    err_p_post += (numerical_post - exact) * (numerical_post - exact) * highjfac;
  } // for (int q=0; q<postquad->NumPoints(); ++q)

  *(params.get<double*>("error_post")) += err_p_post;

  return;
}

/*----------------------------------------------------------------------*
 * Element init
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ElementInit(DRT::ELEMENTS::Acou* ele,
    Teuchos::ParameterList& params)
{
  // each element has to store the interior vectors by itseld, p-adaptivity or not
  // so, shape it, as you need it
  dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");
  bool withpmls = params.get<bool>("withPML");
  if(withpmls)
    ele->eleinteriorAuxiliaryPML_.Shape(shapes_->ndofs_ * nsd_,1);

  ele->eleinteriorVelnp_.Shape(shapes_->ndofs_ * nsd_, 1);
  ele->eleinteriorPressnp_.Shape(shapes_->ndofs_, 1);

  if(dyna_==INPAR::ACOU::acou_ader || dyna_==INPAR::ACOU::acou_ader_lts || dyna_==INPAR::ACOU::acou_ader_tritet)
  {
    ele->eleADERimprovedGrad_.Shape(shapes_->ndofs_ * nsd_, 1);
    ele->eleADERimprovedDiv_.Shape(shapes_->ndofs_, 1);
  }

  if (params.get<bool>("padaptivity"))
    ele->elenodeTrace_.Shape(nen_, 1);

  return;
}

/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ProjectField(
    DRT::ELEMENTS::Acou* ele,
    Teuchos::ParameterList& params,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2)
{
  shapes_.Evaluate(*ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 ||
      unsigned(elevec2.M()) == (nsd_+1)*shapes_.ndofs_, "Wrong size in project vector 2");

  // get function
  const int *start_func = params.getPtr<int>("funct");

  // internal variables
  if (elevec2.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(shapes_.ndofs_,nsd_+ 1);

    for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
    {
      const double fac = shapes_.jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_.xyzreal(d, q); // coordinates of quadrature point in real coordinates
      double p;
      double gradient[nsd_];

      dsassert(start_func != NULL,"funct not set for initial value");
      EvaluateAll(*start_func, xyz, p, gradient); // u and p at quadrature point

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      {
        massPart(i, q) = shapes_.shfunct(i, q);
        massPartW(i, q) = shapes_.shfunct(i, q) * fac;
        localMat(i,nsd_) += shapes_.shfunct(i, q) * p * fac;
        for(unsigned int j=0; j<nsd_; ++j)
          localMat(i,j) += shapes_.shfunct(i,q) * gradient[j] * fac;
      }
    }
    Mmat.Multiply('N', 'T', 1., massPart, massPartW, 0.);
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(Mmat);
      inverseMass.SetVectors(localMat, localMat);
      inverseMass.Solve();
    }

     //standard L2-projection
    for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
    {
      ele->eleinteriorPressnp_(r) += localMat(r, nsd_); // pressure
      for (unsigned int i=0;i<nsd_;++i)
        ele->eleinteriorVelnp_(i*shapes_.ndofs_+r) += localMat(r, i); // velocity
    }

    if(trac_with_pml)
      ele->eleinteriorAuxiliaryPML_.Scale(0.0);

    // projection according to Cockburn, formula 2.8 in "uniform-in-time superconvergence of the hdg methods for the acosutic wave equation" */
    /*
    // need smaller polynomial space
    DRT::UTILS::PolynomialSpace<nsd_> poly_kmin1(distype, ele->Degree()-1, ele->UsesCompletePolynomialSpace());
    unsigned int ndofskmin1 = poly_kmin1.Size();
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Faces()[0]->Degree(),shapes_.usescompletepoly_, 2 * ele->Faces()[0]->Degree());
      shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      shapesface_->EvaluateFace(*ele, 0);
    }
    int systemsize = shapes_.ndofs_*(nsd_+1);
    if(ndofskmin1*(nsd_+1)+shapes_.nfaces_*shapesface_->nfdofs_ != shapes_.ndofs_*(nsd_+1))
      dserror("this projection does not work for standard quad and hex");

    Epetra_SerialDenseMatrix h(systemsize,systemsize);
    Epetra_SerialDenseVector rhs(systemsize);
    LINALG::Matrix<nsd_, 1> xsi;
    Epetra_SerialDenseVector values(ndofskmin1);

    // fill system matrix and rhs vector with correspondent values (this quadrature should be enough)
    for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
    {
      // evaluation of shape functions k-1 in quadrature point
      const double* gpcoord = shapes_.quadrature_->Point(q);
      for (unsigned int idim = 0; idim < nsd_; idim++)
        xsi(idim) = gpcoord[idim];
      poly_kmin1.Evaluate(xsi, values);

      // evaluation of initial fields in quadrature point
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_.xyzreal(d, q); // coordinates of quadrature point in real coordinates
      double p;
      double gradient[nsd_];
      EvaluateAll(*start_func, xyz, p, gradient);

      for (unsigned int i = 0; i <ndofskmin1; ++i)
      {
        for(unsigned int d=0; d<nsd_; ++d)
          rhs(i+d*ndofskmin1) += values(i)*gradient[d]*shapes_.jfac(q);
        rhs(i+nsd_*ndofskmin1) += values(i)*p*shapes_.jfac(q);
        for(unsigned int j=0; j<ndofs_; ++j)
        {
          for(unsigned int d=0; d<nsd_+1; ++d)
            h(i+d*ndofskmin1,j+d*ndofs_) += values(i)*shapes_.shfunct(j,q)*shapes_.jfac(q);
        }
      }
    }
    // fill last rows of the matrix with Imat and Jmat
    for(unsigned int i=0; i<shapes_.nfaces_*shapesface_->nfdofs_;++i)
    {
      for(unsigned int j=0; j<nsd_*ndofs_; ++j)
        h((nsd_+1)*ndofskmin1+i,j) = Imat(i,j);
      for(unsigned int j=0; j<ndofs_; ++j)
        h((nsd_+1)*ndofskmin1+i,nsd_*ndofs_+j) = Jmat(i,j);
    }
    // last rows of rhs vector
    //for(unsigned int i=0; i<shapes_.nfaces_*shapesface_->nfdofs_;++i)
    //{
    //  for(unsigned int j=0; j<ndofs_; ++j)
    //    rhs((nsd_+1)*ndofskmin1+i) += Jmat(i,j)*localMat(j, nsd_);
    //}
    // L2 projection on face
    Epetra_SerialDenseVector rhstemp(shapes_.nfaces_*shapesface_->nfdofs_);
    for(unsigned int f=0; f<nfaces_; ++f)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Faces()[f]->Degree(),shapes_.usescompletepoly_, 2 * ele->Faces()[f]->Degree());
      shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      shapesface_->EvaluateFace(*ele, f);

      for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
      {
        double xyz[nsd_];
        for (unsigned int d = 0; d < nsd_; ++d)
          xyz[d] = shapesface_->xyzreal(d, q); // coordinates of quadrature point in real coordinates
        double p;
        double gradient[nsd_];
        EvaluateAll(*start_func, xyz, p, gradient);

        for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
          rhstemp(i+f*shapesface_->nfdofs_) += shapesface_->jfac(q) * p * shapesface_->shfunct(i,q);
      }
    }
    Epetra_SerialDenseMatrix Gtoinv(Gmat.M(),Gmat.N());
    Gtoinv += Gmat;
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(Gtoinv);
      inverseMass.SetVectors(rhstemp, rhstemp);
      inverseMass.Solve();
    }
    for(unsigned int i=0; i<shapes_.nfaces_*shapesface_->nfdofs_;++i)
    {
      for(unsigned int j=0; j<shapes_.nfaces_*shapesface_->nfdofs_; ++j)
        rhs((nsd_+1)*ndofskmin1+i) += Gmat(i,j)*rhstemp(j);
    }


    for(unsigned int f=0; f<nfaces_; ++f)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Faces()[f]->Degree(),shapes_.usescompletepoly_, 2 * ele->Faces()[f]->Degree());
      shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      shapesface_->EvaluateFace(*ele, f);

      for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
      {
        double xyz[nsd_];
        for (unsigned int d = 0; d < nsd_; ++d)
          xyz[d] = shapesface_->xyzreal(d, q); // coordinates of quadrature point in real coordinates
        double p;
        double gradient[nsd_];
        EvaluateAll(*start_func, xyz, p, gradient);

        for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
          for(unsigned int d=0; d<nsd_; ++d)
            rhs((nsd_+1)*ndofskmin1+i+f*shapesface_->nfdofs_) += shapesface_->jfac(q) * gradient[d]* shapesface_->normals(d,q) * shapesface_->shfunct(i,q);
      }
    }



    // solve the system
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(h);
      inverseMass.SetVectors(rhs, rhs);
      inverseMass.Solve();
    }


    // write values to element vectors
    for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
    {
      ele->eleinteriorPressnp_(r) += rhs(r+nsd_*ndofs_); // pressure
      for (unsigned int i=0;i<nsd_;++i)
        ele->eleinteriorVelnp_(i*shapes_.ndofs_+r) += rhs(r+i*ndofs_); // velocity
    }*/

    /* INTERPOLATION (-> no superconvergence possible)
    for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
    {
      double xyz[nsd_];
      double p;
      double gradient[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_.nodexyzreal(d,r);
      EvaluateAll(*start_func, xyz, p, gradient);
      ele->eleinteriorPressnp_(r) += p; //localMat(r, nsd_); // pressure
      for (unsigned int i=0;i<nsd_;++i)
        ele->eleinteriorVelnp_(i*shapes_.ndofs_+r) += gradient[i]; //velocity
    } */

  }

  return 0;
}

/*----------------------------------------------------------------------*
 * ProjectDirichField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ProjectDirichField(
    DRT::ELEMENTS::Acou* ele,
    Teuchos::ParameterList& params,
    Epetra_SerialDenseVector& elevec1)
{
  shapes_.Evaluate(*ele);

  // in case this paramter "faceconsider" is set, we are applying Dirichlet values
  // and have to evaluate the trace field for the given face!
  if (params.isParameter("faceconsider"))
  {
    Teuchos::Array<int> *functno = params.getPtr<Teuchos::Array<int> >("funct");
    const unsigned int *faceConsider = params.getPtr<unsigned int>("faceconsider");
    //Teuchos::Array<int> *onoff = params.getPtr<Teuchos::Array<int> >("onoff");
    //double *time = params.getPtr<double>("time");

    DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Faces()[*faceConsider]->Degree(),shapes_.usescompletepoly_,2 * ele->Faces()[*faceConsider]->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, *faceConsider);

    Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
    Epetra_SerialDenseVector trVec(shapesface_->nfdofs_);

    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      const double fac = shapesface_->jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapesface_->xyzreal(d, q);
      double p;
      double gradient[nsd_];

      EvaluateAll((*functno)[0], xyz, p, gradient);

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
          mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;
        trVec(i) += shapesface_->shfunct(i, q) * p * fac;
      }
    }

    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec, trVec);
    inverseMass.Solve();

    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      elevec1(i) = trVec(i);
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * ProjectOpticalField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ProjectOpticalField(
    DRT::ELEMENTS::Acou* ele,
    Teuchos::ParameterList& params,
    Epetra_SerialDenseVector& elevec2)
{
  shapes_.Evaluate(*ele);

  bool meshconform = params.get<bool>("mesh conform");

  // this is the easy case: the corresponding optical element has exactly the same nodes
  // and hence we can easily get the scatra dofs of the nodes and substitute the pressure
  // values
  int numlightnode = nen_;

  double lightxyz[numlightnode][nsd_];
  double values[numlightnode];
  double absorptioncoeff = -1.0;

  if (meshconform)
  {
    // get the absorption coefficient
    absorptioncoeff = params.get<double>("absorption");

    Teuchos::RCP<std::vector<double> > nodevals = params.get<Teuchos::RCP<std::vector<double> > >("nodevals");

    if (int((*nodevals).size() / (nsd_ + 1)) != numlightnode)
      dserror("node number not matching");

    // meanwhile: check if we are at least near the element, otherwise we can already go home
    for (int i = 0; i < numlightnode; ++i)
    {
      for (unsigned int j = 0; j < nsd_; ++j)
        lightxyz[i][j] = (*nodevals)[i * (nsd_ + 1) + j];
      values[i] = (*nodevals)[i * (nsd_ + 1) + nsd_];
    }

  }
  else
  {
    // get node based values!
    Teuchos::RCP<const Epetra_Vector> matrix_state = params.get<Teuchos::RCP<Epetra_Vector> >("pressurenode");
    std::vector<double> nodevals;
    DRT::UTILS::ExtractMyNodeBasedValues(ele, nodevals, *matrix_state);
    if (nodevals.size() != nen_)
      dserror("node number not matching: %d vs. %d", nodevals.size(), nen_);

    for (int i = 0; i < numlightnode; ++i)
    {
      for (unsigned int j = 0; j < nsd_; ++j)
        lightxyz[i][j] = shapes_.xyze(j, i);
      values[i] = nodevals[i];
    }
  }

  // now we have locations "lightxyz" and corresponding phis in "values"
  // what we want to do now is to evaluate the field in the locations of the
  // Gauss points in this element, to assign an initial distribution
  // so we do exactly the same as in ProjectInitalField BUT call an other
  // evaluate function which will interpolate from given lightxyz and values

  // internal variables
  if (elevec2.M() > 0)
  {
    /*Epetra_SerialDenseMatrix localMat(shapes_.ndofs_, 1);

    for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
    {
      const double fac = shapes_.jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_.xyzreal(d, q); // coordinates of quadrature point in real coordinates
      double p = 0.0;
      EvaluateLight(lightxyz, values, numlightnode, xyz, p, absorptioncoeff); // p at quadrature point

      for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      {
        massPart(i, q) = shapes_.shfunct(i, q);
        massPartW(i, q) = shapes_.shfunct(i, q)* fac;
        localMat(i, 0) += shapes_.shfunct(i, q) * p * fac;
      }
    }

    Mmat.Multiply('N', 'T', 1., massPart, massPartW, 0.);
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(Mmat);
      inverseMass.SetVectors(localMat, localMat);
      inverseMass.Solve();
    }

    for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
    {
      elevec2[r * (nsd_ + 1) + nsd_] += localMat(r, 0); // pressure
      ele->eleinteriorPressnp_(r) += localMat(r, 0); // pressure
    }
    */
    for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
    {
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_.nodexyzreal(d, r);
      double p = 0.0;
      EvaluateLight(lightxyz, values, numlightnode, xyz, p, absorptioncoeff); // p at quadrature point

      elevec2[r * (nsd_ + 1) + nsd_] += p; // pressure
      ele->eleinteriorPressnp_(r) += p; // pressure
    }


  }

  return 0;
}
/*----------------------------------------------------------------------*
 * EvaluateLight
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::EvaluateLight(
    double lightxyz[][nsd_], double values[], int numnode,
    const double(&xyz)[nsd_], double &p,
    double absorptioncoeff) const
{

  // interpolation from nodes
  if (distype == DRT::Element::quad4)
  {
    if (numnode != 4)
      dserror("wrong number of nodes given");

    //*******************************************************************
    /*LINALG::Matrix<4, 4> coeff(true);
    for(int i=0; i<4; ++i)
    {
      coeff(i,0)=1.0;
      coeff(i,1)=lightxyz[i][0];
      coeff(i,2)=lightxyz[i][1];
      coeff(i,3)=lightxyz[i][0]*lightxyz[i][1];
    }
    LINALG::Matrix<4,4> coeff_N(true);
    coeff_N(0,0)=1.0;
    coeff_N(1,1)=1.0;
    coeff_N(2,2)=1.0;
    coeff_N(3,3)=1.0;
    {
      LINALG::FixedSizeSerialDenseSolver<4,4,4> inverseCoeff;
      inverseCoeff.SetMatrix(coeff);
      inverseCoeff.SetVectors(coeff_N,coeff_N);
      int err = inverseCoeff.Solve();
      if(err != 0) dserror("Inversion of matrix in light evaluation failed with errorcode %d",err);
    }

    p = ( coeff_N(0,0) + coeff_N(1,0) * xyz[0] + coeff_N(2,0) * xyz[1] + coeff_N(3,0) * xyz[0] * xyz[1] ) * values[0]
      + ( coeff_N(0,1) + coeff_N(1,1) * xyz[0] + coeff_N(2,1) * xyz[1] + coeff_N(3,1) * xyz[0] * xyz[1] ) * values[1]
      + ( coeff_N(0,2) + coeff_N(1,2) * xyz[0] + coeff_N(2,2) * xyz[1] + coeff_N(3,2) * xyz[0] * xyz[1] ) * values[2]
      + ( coeff_N(0,3) + coeff_N(1,3) * xyz[0] + coeff_N(2,3) * xyz[1] + coeff_N(3,3) * xyz[0] * xyz[1] ) * values[3];
     */
    double average = (values[0]+values[1]+values[2]+values[3])/4.0;
    //if(p>10.0*average|| p<average/10.0)
    p = average;
    p *= -absorptioncoeff;
  }
  else if (distype == DRT::Element::hex8)
  {
    if(numnode != 8) dserror("wrong number of nodes given");

    //*******************************************************************
    LINALG::Matrix<8,8> coeff(true);
    for(int i=0; i<8; ++i)
    {
      coeff(i,0)=1.0;
      coeff(i,1)=lightxyz[i][0];
      coeff(i,2)=lightxyz[i][1];
      coeff(i,3)=lightxyz[i][2];
      coeff(i,4)=lightxyz[i][0]*lightxyz[i][1];
      coeff(i,5)=lightxyz[i][0]*lightxyz[i][2];
      coeff(i,6)=lightxyz[i][1]*lightxyz[i][2];
      coeff(i,7)=lightxyz[i][0]*lightxyz[i][1]*lightxyz[i][2];
    }
    LINALG::Matrix<8,8> coeff_N(true);
    for(int i=0; i<8; ++i)
      coeff_N(i,i)=1.0;

    {
      LINALG::FixedSizeSerialDenseSolver<8,8,8> inverseCoeff;
      inverseCoeff.SetMatrix(coeff);
      inverseCoeff.SetVectors(coeff_N,coeff_N);
      int err = inverseCoeff.Solve();
      if(err != 0) dserror("Inversion of matrix in light evaluation failed with errorcode %d",err);
    }

    p = ( coeff_N(0,0) + coeff_N(1,0) * xyz[0] + coeff_N(2,0) * xyz[1] + coeff_N(3,0) * xyz[2] + coeff_N(4,0) * xyz[0] * xyz[1] + coeff_N(5,0) * xyz[0] * xyz[2] + coeff_N(6,0) * xyz[1] * xyz[2] + coeff_N(7,0) * xyz[0] * xyz[1] * xyz[2] ) * values[0]
      + ( coeff_N(0,1) + coeff_N(1,1) * xyz[0] + coeff_N(2,1) * xyz[1] + coeff_N(3,1) * xyz[2] + coeff_N(4,1) * xyz[0] * xyz[1] + coeff_N(5,1) * xyz[0] * xyz[2] + coeff_N(6,1) * xyz[1] * xyz[2] + coeff_N(7,1) * xyz[0] * xyz[1] * xyz[2] ) * values[1]
      + ( coeff_N(0,2) + coeff_N(1,2) * xyz[0] + coeff_N(2,2) * xyz[1] + coeff_N(3,2) * xyz[2] + coeff_N(4,2) * xyz[0] * xyz[1] + coeff_N(5,2) * xyz[0] * xyz[2] + coeff_N(6,2) * xyz[1] * xyz[2] + coeff_N(7,2) * xyz[0] * xyz[1] * xyz[2] ) * values[2]
      + ( coeff_N(0,3) + coeff_N(1,3) * xyz[0] + coeff_N(2,3) * xyz[1] + coeff_N(3,3) * xyz[2] + coeff_N(4,3) * xyz[0] * xyz[1] + coeff_N(5,3) * xyz[0] * xyz[2] + coeff_N(6,3) * xyz[1] * xyz[2] + coeff_N(7,3) * xyz[0] * xyz[1] * xyz[2] ) * values[3]
      + ( coeff_N(0,4) + coeff_N(1,4) * xyz[0] + coeff_N(2,4) * xyz[1] + coeff_N(3,4) * xyz[2] + coeff_N(4,4) * xyz[0] * xyz[1] + coeff_N(5,4) * xyz[0] * xyz[2] + coeff_N(6,4) * xyz[1] * xyz[2] + coeff_N(7,4) * xyz[0] * xyz[1] * xyz[2] ) * values[4]
      + ( coeff_N(0,5) + coeff_N(1,5) * xyz[0] + coeff_N(2,5) * xyz[1] + coeff_N(3,5) * xyz[2] + coeff_N(4,5) * xyz[0] * xyz[1] + coeff_N(5,5) * xyz[0] * xyz[2] + coeff_N(6,5) * xyz[1] * xyz[2] + coeff_N(7,5) * xyz[0] * xyz[1] * xyz[2] ) * values[5]
      + ( coeff_N(0,6) + coeff_N(1,6) * xyz[0] + coeff_N(2,6) * xyz[1] + coeff_N(3,6) * xyz[2] + coeff_N(4,6) * xyz[0] * xyz[1] + coeff_N(5,6) * xyz[0] * xyz[2] + coeff_N(6,6) * xyz[1] * xyz[2] + coeff_N(7,6) * xyz[0] * xyz[1] * xyz[2] ) * values[6]
      + ( coeff_N(0,7) + coeff_N(1,7) * xyz[0] + coeff_N(2,7) * xyz[1] + coeff_N(3,7) * xyz[2] + coeff_N(4,7) * xyz[0] * xyz[1] + coeff_N(5,7) * xyz[0] * xyz[2] + coeff_N(6,7) * xyz[1] * xyz[2] + coeff_N(7,7) * xyz[0] * xyz[1] * xyz[2] ) * values[7];
    p *= -absorptioncoeff;
  }
  else
    dserror("not yet implemented");

  return;
}
//template<DRT::Element::DiscretizationType distype>
//int DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ProjectOpticalField(
//    DRT::ELEMENTS::Acou* ele,
//    Teuchos::ParameterList& params,
//    Epetra_SerialDenseVector& elevec2)
//{
//  shapes_.Evaluate(*ele);
//
//  double* energy = params.get<double*>("gpvalues");
//
//  Epetra_SerialDenseMatrix localMat(shapes_.ndofs_,1);
//
//  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
//  {
//    const double fac = shapes_.jfac(q);
//    for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
//    {
//      massPart(i, q) = shapes_.shfunct(i, q);
//      massPartW(i, q) = shapes_.shfunct(i, q) * fac;
//      localMat(i,0) -= shapes_.shfunct(i, q) * energy[q] * fac;
//    }
//  }
//
//  Mmat.Multiply('N', 'T', 1., massPart, massPartW, 0.);
//  {
//    Epetra_SerialDenseSolver inverseMass;
//    inverseMass.SetMatrix(Mmat);
//    inverseMass.SetVectors(localMat, localMat);
//    inverseMass.Solve();
//  }
//
//   //standard L2-projection
//  for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
//    ele->eleinteriorPressnp_(r) += localMat(r, 0); // pressure
//  ele->eleinteriorVelnp_.Scale(0.0);
//
//  if(trac_with_pml)
//    ele->eleinteriorAuxiliaryPML_.Scale(0.0);
//
//  return 0;
//}

/*----------------------------------------------------------------------*
 * EvaluateAll
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::EvaluateAll(
    const int start_func, const double(&xyz)[nsd_], double &p, double (&v)[nsd_]) const
{
  p = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(0, xyz, 0.0);
  if(DRT::Problem::Instance()->Funct(start_func-1).NumberComponents()==nsd_+1)
  {
    for(unsigned int d=0; d<nsd_; ++d)
      v[d] = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(d+1, xyz, 0.0);
  }
  else if(DRT::Problem::Instance()->Funct(start_func-1).NumberComponents()!=1)
    dserror("Supply ONE component for your start function or NUMDIM+1, not anything else! The first component is for the pressure, the others for the velocity.");
  else
  {
    for(unsigned int d=0; d<nsd_; ++d)
      v[d] = 0.0;
  }
  return;
}

/*----------------------------------------------------------------------*
 * ComputePressureAverage
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ComputePressureAverage(Epetra_SerialDenseVector& elevec)
{
  double norm_p = 0.0, area = 0.0;
  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    double numerical_p = 0.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      numerical_p += shapes_->shfunct(i, q) * interiorPressnp_(i);
    norm_p += numerical_p * shapes_->jfac(q);
    area += shapes_->jfac(q);
  }
  elevec[0] = norm_p / area;

  return;
}

/*----------------------------------------------------------------------*
 * NodeBasedValues
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::NodeBasedValues(
    DRT::ELEMENTS::Acou* ele, Epetra_SerialDenseVector& elevec1,
    bool padaptivity)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::NodeBasedValues");

  dsassert(elevec1.M() == (int)nen_*(nsd_+2)+1, "Vector does not have correct size");
  elevec1.Scale(0.0);
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  double norm_p = 0.0, area = 0.0;
  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    double numerical_p = 0.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      numerical_p += shapes_->shfunct(i, q) * interiorPressnp_(i);
    norm_p += numerical_p * shapes_->jfac(q);
    area += shapes_->jfac(q);
  }
  elevec1((nsd_ + 2) * nen_) = norm_p / area;

  for (unsigned int i = 0; i < nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim = 0; idim < nsd_; idim++)
      shapes_->xsi(idim) = locations(idim, i);
    shapes_->polySpace_->Evaluate(shapes_->xsi, values);

    // compute values for velocity and pressure by summing over all basis functions
    double sum = 0.0;
    for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
    {
      sum += values(k) * interiorPressnp_(k);
    }
    elevec1(nsd_ * nen_ + i) = sum;

    for (unsigned int d = 0; d < nsd_; ++d)
    {
      sum = 0.0;
      for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
        sum += values(k) * interiorVelnp_(d * shapes_->ndofs_ + k);
      elevec1(d * nen_ + i) = sum;
    }
  }
  if ( (!padaptivity && (traceVal_.size()>0)) && dyna_==INPAR::ACOU::acou_impleuler) // (trace val size is zero for explicit time integration)
  {
    Epetra_SerialDenseVector fvalues(1);
    Epetra_SerialDenseVector touchcount(nen_);
    int sumindex = 0;
    for (unsigned int face = 0; face < nfaces_; ++face)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(
          ele->Faces()[face]->Degree(),
          usescompletepoly_,
          2 * ele->Faces()[face]->Degree());
      localSolver_->shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      localSolver_->shapesface_->EvaluateFace(*ele, face);

      fvalues.Resize(localSolver_->shapesface_->nfdofs_);

      // const int * fnodeIds = ele->Faces()[face]->NodeIds();
      for (int i = 0;i < DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
      {
        // evaluate shape polynomials in node
        for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
          localSolver_->shapesface_->xsi(idim) = locations(idim, i);

        localSolver_->shapesface_->polySpace_->Evaluate(localSolver_->shapesface_->xsi, fvalues); // TODO: fix face orientation here

        // compute values for velocity and pressure by summing over all basis functions
        double sum = 0.0;
        for (unsigned int k = 0; k < localSolver_->shapesface_->nfdofs_; ++k)
          sum += fvalues(k) * traceVal_[sumindex + k];

        elevec1((nsd_ + 1) * nen_ + localSolver_->shapesface_->faceNodeOrder[face][i]) += sum;
        touchcount(localSolver_->shapesface_->faceNodeOrder[face][i])++;
      }
      sumindex += localSolver_->shapesface_->nfdofs_;
    }

    for (unsigned int i = 0; i < nen_; ++i)
      elevec1((nsd_ + 1) * nen_ + i) /= touchcount(i);

  }
  else if (padaptivity)
    for (unsigned int i = 0; i < nen_; ++i)
      elevec1((nsd_ + 1) * nen_ + i) = ele->elenodeTrace_(i);

  return;
} // NodeBasedValues

/*----------------------------------------------------------------------*
 * NodeBasedPsi
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::NodeBasedPsi(
    Epetra_SerialDenseVector& elevec1, double dt, double fac)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::NodeBasedValues");

  dsassert(elevec1.M() == int(nen_), "Vector does not have correct size");
  elevec1.Scale(0.0);
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  for (unsigned int i = 0; i < nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim = 0; idim < nsd_; idim++)
      shapes_->xsi(idim) = locations(idim, i);
    shapes_->polySpace_->Evaluate(shapes_->xsi, values);

    // compute values for velocity and pressure by summing over all basis functions
    double sum = 0;
    for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
      sum += values(k) * interiorPressnp_(k);

    elevec1(i) = sum / fac / dt;
  }

  return;
}

/*----------------------------------------------------------------------*
 * Instance
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype> *
DRT::ELEMENTS::AcouEleCalc<distype>::Instance(bool create)
{
  static AcouEleCalc<distype> * instance;
  if (create)
  {
    if (instance == NULL)
      instance = new AcouEleCalc<distype>();
  }
  else
  {
    if (instance != NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 * Done
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}

/*----------------------------------------------------------------------*
 * Constructor LocalSolver
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::LocalSolver(DRT::UTILS::ShapeValues<distype> &shapeValues, INPAR::ACOU::DynamicType &dyna) :
    ndofs_(shapeValues.ndofs_),
    shapes_(shapeValues),
    dyna_(dyna)
{
  // shape all matrices
  reshapeMatrixIfNecessary(Amat, nsd_ * ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(invAmat, nsd_ * ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(Bmat, nsd_ * ndofs_, ndofs_);
  reshapeMatrixIfNecessary(Mmat, ndofs_, ndofs_);
  reshapeMatrixIfNecessary(Dmat, ndofs_, ndofs_);
  reshapeMatrixIfNecessary(Hmat, ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(massPart, ndofs_, shapeValues.nqpoints_);
  reshapeMatrixIfNecessary(massPartW, ndofs_, shapeValues.nqpoints_);
}

/*----------------------------------------------------------------------*
 * FaceSpecificConstruction
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::FaceSpecificConstruction(
    const DRT::ELEMENTS::Acou* ele, bool completepoly)
{
  int onfdofs = 0;
  int nfdofs = 0;
  for (unsigned int i = 0; i < nfaces_; ++i)
  {
    DRT::UTILS::PolynomialSpaceParams params(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,ele->Faces()[i]->Degree(), completepoly);
    nfdofs = DRT::UTILS::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(params)->Size();
    onfdofs += nfdofs;
  }

  reshapeMatrixIfNecessary(Cmat, nsd_ * ndofs_, onfdofs);
  reshapeMatrixIfNecessary(Emat, ndofs_, onfdofs);
  reshapeMatrixIfNecessary(Gmat, onfdofs, onfdofs);
  reshapeMatrixIfNecessary(Imat, onfdofs, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(Jmat, onfdofs, ndofs_);

  return;
}

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidual
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::UpdateInteriorVariablesAndComputeResidual(
    Teuchos::ParameterList& params,
    DRT::ELEMENTS::Acou & ele, const Teuchos::RCP<MAT::Material> &mat,
    Epetra_SerialDenseVector & elevec, double dt, bool errormaps,
    bool updateonly)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::UpdateInteriorVariablesAndComputeResidual");

  int onfdofs = elevec.M();

  Epetra_SerialDenseVector tempVelnp;
  Epetra_SerialDenseVector tempPressnp;

  Epetra_SerialDenseVector traceVal_SDV(onfdofs);
  for (int i = 0; i < onfdofs; ++i)
    traceVal_SDV(i) = traceVal_[i];

  // *****************************************************
  // update interior variables first
  // *****************************************************

  Epetra_SerialDenseVector tempVec1(shapes_->ndofs_ * nsd_);
  tempVec1.Multiply('N', 'N', 1.0, localSolver_->Amat, interiorVelnp_, 0.0);
  tempVec1.Multiply('N', 'N', -1.0, localSolver_->Cmat, traceVal_SDV, 1.0);

  Epetra_SerialDenseVector tempVec2(shapes_->ndofs_);
  Epetra_SerialDenseVector tempVec3(shapes_->ndofs_);
  localSolver_->ComputeSource(params,tempVec2,tempVec3);

  tempVec2.Multiply('N', 'N', 1.0, localSolver_->Mmat, interiorPressnp_, 1.0);
  tempVec2.Multiply('N', 'N', -1.0, localSolver_->Emat, traceVal_SDV, 1.0);

  // now, we have to do the Schur complement thing, just as in "CondenseLocalPart" but without C and E
  Epetra_SerialDenseMatrix tempMat1(shapes_->ndofs_, nsd_ * shapes_->ndofs_);
  tempMat1.Multiply('N', 'N', 1.0, localSolver_->Hmat, localSolver_->invAmat,0.0);
  tempVec2.Multiply('N', 'N', -1.0, tempMat1, tempVec1, 1.0);

  Epetra_SerialDenseMatrix tempMat2(shapes_->ndofs_, shapes_->ndofs_);
  tempMat2 = localSolver_->Dmat;
  tempMat2.Multiply('N', 'N', -1.0, tempMat1, localSolver_->Bmat, 1.0);
  tempMat2 += localSolver_->Mmat;

  {
    Epetra_SerialDenseSolver inverseDmat;
    inverseDmat.SetMatrix(tempMat2);
    int err = inverseDmat.Invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d",err);
  }

  interiorPressnp_.Multiply('N', 'N', 1.0, tempMat2, tempVec2, 0.0);
  tempVec1.Multiply('N', 'N', -1.0, localSolver_->Bmat, interiorPressnp_,1.0);
  interiorVelnp_.Multiply('N', 'N', 1.0, localSolver_->invAmat, tempVec1, 0.0);

  // tell this change in the interior variables the discretization
  bool padaptivity = params.get<bool>("padaptivity");
  if (!padaptivity)
  {
    ele.eleinteriorPressnp_ = interiorPressnp_;
    ele.eleinteriorVelnp_ = interiorVelnp_;
  }

  // *****************************************************
  // local postprocessing to calculate error maps
  // *****************************************************

  if (errormaps)
  {
    const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
    double rho = actmat->Density();

    // first step: calculate the gradient we need (bold p in paper by Nguyen)
    Epetra_SerialDenseVector temp(shapes_->ndofs_ * nsd_);
    temp.Multiply('N', 'N', 1.0, localSolver_->Bmat, interiorPressnp_, 0.0);
    temp.Multiply('N', 'N', 1.0, localSolver_->Cmat, traceVal_SDV, 1.0);
    Epetra_SerialDenseVector p(shapes_->ndofs_ * nsd_);
    p.Multiply('N', 'N', rho / dt, localSolver_->invAmat, temp, 0.0);

    // second step: postprocess the pressure field
    double err_p = EstimateError(ele, p);

    Teuchos::RCP<std::vector<double> > values = params.get<Teuchos::RCP<std::vector<double> > >("elevals");

    if (padaptivity)
    {
      double padaptivitytol = params.get<double>("padaptivitytol");
      double delta_k = 0.0;
      if (err_p != 0.0)
        delta_k = std::ceil(std::log10(err_p / padaptivitytol));

      int newdeg = (ele.Degree() + delta_k);
      if (newdeg < 1)
        newdeg = 1;
      if (newdeg > 6)
        newdeg = 6;
      (*values)[ele.Id()] = newdeg;

      // projection of the internal field (cannot yet project trace field, since the
      // face degree depends on all neighbors! (also not needed for DIRK)
      ProjectPadapField(ele, newdeg);
      updateonly = true;
      return;
    }
    else
      (*values)[ele.Id()] = err_p;
  }

  if (updateonly)
    return;

  // *****************************************************
  // compute residual second (reuse intermediate matrices)
  // *****************************************************

  tempVec1.Resize(shapes_->ndofs_);
  tempVec1.Multiply('N', 'N', 1.0, localSolver_->Mmat, interiorPressnp_, 0.0);
  // ComputeSource(...)
  tempVec1+=tempVec3;

  Epetra_SerialDenseVector f(shapes_->ndofs_ * nsd_);
  f.Multiply('N', 'N', 1.0, localSolver_->Amat, interiorVelnp_, 0.0);
  tempVec1.Multiply('N', 'N', -1.0, tempMat1, f, 1.0);

  tempVec2.Resize(shapes_->ndofs_);
  tempVec2.Multiply('N', 'N', 1.0, tempMat2, tempVec1, 0.0);

  elevec.Multiply('N', 'N', -1.0, localSolver_->Jmat, tempVec2, 0.0);

  tempVec1.Resize(shapes_->ndofs_ * nsd_);
  tempVec1.Multiply('N', 'N', -1.0, localSolver_->Bmat, tempVec2, 0.0);
  tempVec1 += f;
  tempVec2.Resize(shapes_->ndofs_ * nsd_);
  tempVec2.Multiply('N', 'N', 1.0, localSolver_->invAmat, tempVec1, 0.0);

  elevec.Multiply('N', 'N', -1.0, localSolver_->Imat, tempVec2, 1.0);

  return;
} // UpdateInteriorVariablesAndComputeResidual

/*----------------------------------------------------------------------*
 * ProjectPadapField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ProjectPadapField(
    DRT::ELEMENTS::Acou & ele, int newdeg)
{
  // create shape values for projection
  Teuchos::RCP<DRT::UTILS::ShapeValues<distype> > proj_shapes = Teuchos::rcp(
      new DRT::UTILS::ShapeValues<distype>(newdeg, usescompletepoly_,2 * newdeg));
  proj_shapes->Evaluate(ele);

  // calculate mass matrix with new degree
  Epetra_SerialDenseMatrix massPart(proj_shapes->ndofs_,proj_shapes->nqpoints_);
  Epetra_SerialDenseMatrix massMat(proj_shapes->ndofs_, proj_shapes->ndofs_);
  for (unsigned int q = 0; q < proj_shapes->nqpoints_; ++q)
  {
    const double sqrtfac = std::sqrt(proj_shapes->jfac(q));
    for (unsigned int i = 0; i < proj_shapes->ndofs_; ++i)
      massPart(i, q) = proj_shapes->shfunct(i, q) * sqrtfac;
  }
  for (unsigned int i = 0; i < proj_shapes->ndofs_; ++i)
    for (unsigned int j = 0; j < proj_shapes->ndofs_; ++j)
    {
      double sum = 0.0;
      for (unsigned int k = 0; k < proj_shapes->nqpoints_; ++k)
        sum += massPart(i, k) * massPart(j, k);
      massMat(i, j) = sum;
    }
  // invert the mass matrix
  {
    Epetra_SerialDenseSolver inversemat;
    inversemat.SetMatrix(massMat);
    inversemat.Invert();
  }

  // calculate rectangular mass matrix with mixed degree
  LINALG::Matrix<nsd_, 1> xsi;
  Epetra_SerialDenseVector values(shapes_->ndofs_);
  Epetra_SerialDenseMatrix rectMat(proj_shapes->ndofs_, shapes_->ndofs_);

  for (unsigned int i = 0; i < proj_shapes->ndofs_; ++i)
    for (unsigned int k = 0; k < proj_shapes->nqpoints_; ++k)
    {
      // calculate values of original shape functions at new quadrature point
      const double* gpcoord = proj_shapes->quadrature_->Point(k);
      for (unsigned int idim = 0; idim < nsd_; idim++)
        xsi(idim) = gpcoord[idim];

      shapes_->polySpace_->Evaluate(xsi, values);
      for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
      {
        rectMat(i, j) += proj_shapes->jfac(k) * proj_shapes->shfunct(i, k) * values(j);
      }
    }
  Epetra_SerialDenseMatrix rectMat2(proj_shapes->ndofs_, shapes_->ndofs_);
  rectMat2.Multiply('N', 'N', 1.0, massMat, rectMat, 0.0);

  // build a nsd_ rectMat
  Epetra_SerialDenseMatrix rectMat3(proj_shapes->ndofs_ * nsd_,
      shapes_->ndofs_ * nsd_);
  for (unsigned d = 0; d < nsd_; ++d)
    for (unsigned int j = 0; j < proj_shapes->ndofs_; ++j)
      for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
        rectMat3(j + proj_shapes->ndofs_ * d, i + shapes_->ndofs_ * d) =
            rectMat2(j, i);

  // save the projected pressure to the element!
  reshapeMatrixIfNecessary(ele.eleinteriorPressnp_, proj_shapes->ndofs_, 1);
  ele.eleinteriorPressnp_.Multiply('N', 'N', 1.0, rectMat2, interiorPressnp_, 0.0);

  // project the velocity field as well!! (reuse matrices)
  reshapeMatrixIfNecessary(ele.eleinteriorVelnp_, proj_shapes->ndofs_ * nsd_, 1);
  ele.eleinteriorVelnp_.Multiply('N', 'N', 1.0, rectMat3, interiorVelnp_, 0.0);

  // calculation of node based trace values, since this is no longer possible afterwards!!
  {
    //ele.elenodeTrace_.Shape(nen_,1);
    reshapeMatrixIfNecessary(ele.elenodeTrace_, nen_, 1);
    Epetra_SerialDenseVector fvalues(1);
    Epetra_SerialDenseVector touchcount(nen_);
    Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

    int sumindex = 0;
    for (unsigned int face = 0; face < nfaces_; ++face)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(
          ele.Faces()[face]->Degree(),
          usescompletepoly_,
          2 * ele.Faces()[face]->Degree());
      localSolver_->shapesface_ =DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      localSolver_->shapesface_->EvaluateFace(ele, face);
      fvalues.Resize(localSolver_->shapesface_->nfdofs_);

      for (int i = 0; i < DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
      {
        // evaluate shape polynomials in node
        for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
          localSolver_->shapesface_->xsi(idim) = locations(idim, i);
        localSolver_->shapesface_->polySpace_->Evaluate(localSolver_->shapesface_->xsi, fvalues); // TODO: fix face orientation here

        // compute values for velocity and pressure by summing over all basis functions
        double sum = 0.0;
        for (unsigned int k = 0; k < localSolver_->shapesface_->nfdofs_; ++k)
          sum += fvalues(k) * traceVal_[sumindex + k];

        ele.elenodeTrace_(localSolver_->shapesface_->faceNodeOrder[face][i]) += sum;
        touchcount(localSolver_->shapesface_->faceNodeOrder[face][i])++;
      }
      sumindex += localSolver_->shapesface_->nfdofs_;
    }

    for (unsigned int i = 0; i < nen_; ++i)
      ele.elenodeTrace_(i) /= touchcount(i);
  }

  return;
}

/*----------------------------------------------------------------------*
 * CalculateError
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::AcouEleCalc<distype>::EstimateError(
    DRT::ELEMENTS::Acou & ele, Epetra_SerialDenseVector & p)
{
  DRT::UTILS::PolynomialSpace<nsd_> postpoly(distype, ele.Degree() + 1, ele.UsesCompletePolynomialSpace());
  LINALG::Matrix<nsd_, 1> xsi;
  int ndofspost = 1;
  for (unsigned int i = 0; i < nsd_; ++i)
    ndofspost *= (ele.Degree() + 2);

  Epetra_SerialDenseMatrix h(ndofspost, ndofspost);
  Epetra_SerialDenseVector rhs(ndofspost);

  Epetra_SerialDenseMatrix derivs(nsd_, ndofspost);
  Epetra_SerialDenseVector myvalues(shapes_->ndofs_);

  LINALG::Matrix<nsd_, nen_> deriv;
  LINALG::Matrix<nsd_, nsd_> xjm, xji;

  Epetra_SerialDenseVector values(ndofspost);

  Teuchos::RCP<DRT::UTILS::GaussPoints> postquad;
  postquad = DRT::UTILS::GaussPointCache::Instance().Create(distype,(ele.Degree() + 1) * 2);

  LINALG::Matrix<nsd_, nen_> xyze;
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_> >(&ele, xyze);

  for (int q = 0; q < postquad->NumPoints(); ++q)
  {
    const double* gpcoord = postquad->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++)
      xsi(idim) = gpcoord[idim];

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
        for (unsigned int e = 1; e < nsd_; ++e)
          res[d] += xji(d, e) * derivs(e, i);
      }
      for (unsigned int d = 0; d < nsd_; ++d)
        derivs(d, i) = res[d];
    }

    shapes_->polySpace_->Evaluate(xsi, myvalues);

    for (int j = 0; j < ndofspost; ++j)
      h(0, j) += values(j) * jfac;
    for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
      rhs(0) += myvalues(j) * jfac * interiorPressnp_[j];

    for (int i = 1; i < ndofspost; ++i)
    {
      for (int j = 0; j < ndofspost; ++j)
      {
        double t = 0;
        for (unsigned int d = 0; d < nsd_; ++d)
          t += derivs(d, i) * derivs(d, j);
        h(i, j) += t * jfac;
      }
    }
    double ugrad[nsd_]= {0.0};

    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
        ugrad[d] += myvalues(j) * p(j + d * shapes_->ndofs_);
    for (int i = 1; i < ndofspost; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        rhs(i) += ugrad[d] * derivs(d, i) * jfac;

  } // for (int q=0; q<postquad->NumPoints(); ++q)

  {
    Epetra_SerialDenseSolver inverseH;
    inverseH.SetMatrix(h);
    inverseH.SetVectors(rhs, rhs);
    inverseH.Solve();
  }

  double err_p = 0.0;
  double numerical_post = 0.0;
  double numerical = 0.0;
  double area = 0.0;
  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    numerical_post = 0.0;
    numerical = 0.0;
    const double* gpcoord = shapes_->quadrature_->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++)
      xsi(idim) = gpcoord[idim];

    postpoly.Evaluate(xsi, values);
    for (int i = 0; i < ndofspost; ++i)
      numerical_post += values(i) * rhs(i);
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      numerical += shapes_->shfunct(i, q) * interiorPressnp_(i);
    area += shapes_->jfac(q);
    err_p += (numerical_post - numerical) * (numerical_post - numerical)
        * shapes_->jfac(q);
  } // for (int q=0; q<postquad->NumPoints(); ++q)

  err_p /= area;

  return err_p;
} // FillMatrices

/*----------------------------------------------------------------------*
 * ComputeABCNodeVals
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ComputePMonNodeVals(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material> & mat, int face,
    Epetra_SerialDenseMatrix &elemat, Epetra_SerialDenseVector &elevec,
    int indexstart)
{

  localSolver_->shapesface_->EvaluateFace(*ele, face);

  Teuchos::RCP<std::vector<int> > indices;
  if (params.isParameter("nodeindices"))
  {
    indices = params.get<Teuchos::RCP<std::vector<int> > >("nodeindices");
    if (int(indices->size()) != elevec.M())
      dserror("Vector does not have correct size for line");
  }
  else if (elevec.M() != DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace)
    dserror("Vector does not have correct size");

  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  Epetra_SerialDenseVector fvalues(localSolver_->shapesface_->nfdofs_);
  LINALG::Matrix<nsd_-1, 1> xsiFl;

  for (int i = 0; i < DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
      xsiFl(idim) = locations(idim, i);
    localSolver_->shapesface_->polySpace_->Evaluate(xsiFl, fvalues); // TODO: fix face orientation here

    // compute values for velocity and pressure by summing over all basis functions
    double sum = 0.0;
    for (unsigned int k = 0; k < localSolver_->shapesface_->nfdofs_; ++k)
      sum += fvalues(k) * traceVal_[indexstart + k];

    if (params.isParameter("nodeindices"))
    {
      // is node i part of the line?
      for (unsigned int j = 0; j < indices->size(); ++j)
        if ((*indices)[j] == i)
        {
          elevec(j) = sum;
        }
    }
    else
      elevec(i) = sum;
  }

  return;
}

/*----------------------------------------------------------------------*
 * ComputeSource
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
ComputeSource(Teuchos::ParameterList&       params,
              Epetra_SerialDenseVector          & interiorSourcen,
              Epetra_SerialDenseVector          & interiorSourcenp)
{
  int funcno = params.get<int>("sourcefuncno");
  if(funcno<0) return; // there is no such thing as a volume force
  if(DRT::Problem::Instance()->Funct(funcno).NumberComponents()>1) dserror("for standard acou, the source term has to be scalar, i.e. only 1 component");

  // the vector to be filled
  Epetra_SerialDenseVector source(shapes_.ndofs_);

  // what time is it?
  double tn = params.get<double>("time");
  double tp = params.get<double>("timep");

  double f_value = 0.0;
  double f_value_p = 0.0;
  for (unsigned int q=0; q<shapes_.nqpoints_; ++q)
  {
    double xyz[nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
      xyz[d] = shapes_.xyzreal(d,q);

    // calculate right hand side contribution for dp/dt
    f_value = DRT::Problem::Instance()->Funct(funcno).Evaluate(0,xyz,tn);
    f_value_p = DRT::Problem::Instance()->Funct(funcno).Evaluate(0,xyz,tp);

    // add it all up
    for(unsigned int i=0; i<shapes_.ndofs_; ++i)
    {
      interiorSourcen(i) += shapes_.shfunct(i,q) * f_value * shapes_.jfac(q);
      interiorSourcenp(i) += shapes_.shfunct(i,q) * f_value_p * shapes_.jfac(q);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 * ComputeAbsorbingBC
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeAbsorbingBC(
    DRT::Discretization & discretization,
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material> & mat, int face,
    Epetra_SerialDenseMatrix &elemat,
    int indexstart)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeAbsorbingBC");

  shapesface_->EvaluateFace(*ele, face);

  bool resonly = params.get<bool>("resonly");

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double c = actmat->SpeedofSound(ele->Id());
  double rho = actmat->Density(ele->Id());

  if (!resonly)
  {
    // loop over number of shape functions
    for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
    {
      // loop over number of shape functions
      for (unsigned int q = 0; q <= p; ++q)
      {
        double tempG = 0.0;
        for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
          tempG += shapesface_->jfac(i) * shapesface_->shfunct(p, i) * shapesface_->shfunct(q, i);

        elemat(indexstart + p, indexstart + q) = elemat(indexstart + q,indexstart + p) -= tempG / c / rho;

      } // for (unsigned int q=0; q<nfdofs_; ++q)
    } // for (unsigned int p=0; p<nfdofs_; ++p)
  } // if(!resonly)

  return;
}

/*----------------------------------------------------------------------*
 * ComputeBoundaryIntegral
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeBoundaryIntegral(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params,
    int face)
{
  DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Faces()[face]->Degree(),shapes_.usescompletepoly_,2 * ele->Faces()[face]->Degree());
  shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
  shapesface_->EvaluateFace(*ele, face);

  // access boundary area variable with its actual value
  double boundaryint = params.get<double>("area");

  for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
    boundaryint += shapesface_->jfac(i);

  // add contribution to the global value
  params.set<double>("area",boundaryint);

  return;
}


/*----------------------------------------------------------------------*
 * ComputeInteriorMatrices
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeInteriorMatrices(double dt, double rho, double c)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeInteriorMatrices");

  Epetra_SerialDenseMatrix gradPart(ndofs_ * nsd_, shapes_.nqpoints_);

  // loop quadrature points
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    // loop shape functions
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      massPart(i, q) = shapes_.shfunct(i, q);
      const double valf = shapes_.shfunct(i, q) * shapes_.jfac(q);
      massPartW(i, q) = valf;
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        const double vald = shapes_.shderxy(i * nsd_ + d, q);
        gradPart(d * ndofs_ + i, q) = vald;
      }
    }
  }

  // multiply matrices to perform summation over quadrature points
  Mmat.Multiply('N', 'T', 1.0 / c / c / rho / dt, massPart, massPartW, 0.0);//Mmat.Multiply('N', 'T', 1.0 / c / c / dt, massPart, massPartW, 0.0);
  Dmat.Multiply('N', 'T', rho / dt, massPart, massPartW, 0.0); // only temporary for smaller inversion
  Bmat.Multiply('N', 'T', -1.0, gradPart, massPartW, 0.0);
  for (unsigned int j = 0; j < ndofs_; ++j)
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        Amat(d * ndofs_ + i, d * ndofs_ + j) = rho * rho* c * c * Mmat(i, j);// rho * c * c * Mmat(i, j);
        Hmat(j, d * ndofs_ + i) =  -1.0 * Bmat(d * ndofs_ + i, j);// -rho * Bmat(d * ndofs_ + i, j);
      }
    }

  // we will need the inverse of A
  {
    Epetra_SerialDenseSolver inverseAmat;
    inverseAmat.SetMatrix(Dmat);
    int err = inverseAmat.Invert();
    if (err != 0)
      dserror("Inversion for Amat failed with errorcode %d", err);
    for (unsigned int j = 0; j < ndofs_; ++j)
      for (unsigned int i = 0; i < ndofs_; ++i)
        for (unsigned int d = 0; d < nsd_; ++d)
          invAmat(d * ndofs_ + i, d * ndofs_ + j) = Dmat(i, j);
  }
  zeroMatrix(Dmat); // delete this, since we used Dmat only for the calculation of invAmat
  return;
}

/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeResidual(
    Teuchos::ParameterList&       params,
    Epetra_SerialDenseVector & elevec,
    Epetra_SerialDenseVector & interiorVeln,
    Epetra_SerialDenseVector & interiorPressn)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeResidual");

  /*
   for implicit Euler
   -1
   +---------+    +---------+
   |         |    |    n-1  |
   n              | A    B  |    | A V     |
   R  = - [ I  J ] |         |    |    n-1  |
   | H   D+M |    | M P     |
   +---------+    +---------+

   for trapezoidal rule
   -1
                   +---------+    +----------------------------------------------------------+
                   |         |    |    n-1                     n-1          n-1              |
    n              | A    B  |    | A V     - (1 - theta) ( B P   + C Lambda    )            |                    n-1     n-1           n-1
   R  = - [ I  J ] |         |    |    n-1                     n-1       n-1           n-1   | - (1 - theta) ( I V   + J P    + G Lambda    )
                   | H   D+M |    | M P     - (1 - theta) ( H V   + DmM P    + E Lambda    ) |
                   +---------+    +----------------------------------------------------------+
   */

  Epetra_SerialDenseVector tempVec1(ndofs_);
  Epetra_SerialDenseVector tempVec2(ndofs_);
  ComputeSource(params,tempVec2,tempVec1);
  tempVec1.Multiply('N', 'N', 1.0, Mmat, interiorPressn, 1.0);

  Epetra_SerialDenseVector f(ndofs_ * nsd_);
  f.Multiply('N', 'N', 1.0, Amat, interiorVeln, 0.0);

  Epetra_SerialDenseMatrix tempMat1(ndofs_, ndofs_ * nsd_);
  tempMat1.Multiply('N', 'N', 1.0, Hmat, invAmat, 0.0);
  tempVec1.Multiply('N', 'N', -1.0, tempMat1, f, 1.0); // right part of w

  Epetra_SerialDenseMatrix tempMat2(ndofs_, ndofs_);
  tempMat2 = Dmat;
  tempMat2 += Mmat;
  tempMat2.Multiply('N', 'N', -1.0, tempMat1, Bmat, 1.0); // = D + M - H A^{-1} B
  {
    Epetra_SerialDenseSolver inverseinW;
    inverseinW.SetMatrix(tempMat2);
    int err = inverseinW.Invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d",err);
  }
  // tempMat2 = ( D + M - H A^{-1} B )^{-1}

  tempVec2.Multiply('N', 'N', 1.0, tempMat2, tempVec1, 0.0); // = w = ( D + M - H A^{-1} B )^{-1} ( ... )
  elevec.Multiply('N', 'N', -1.0, Jmat, tempVec2, 0.0);

  tempVec1.Resize(ndofs_ * nsd_);
  tempVec1.Multiply('N', 'N', -1.0, Bmat, tempVec2, 0.0);
  tempVec1 += f;
  tempVec2.Shape(ndofs_ * nsd_, 1);
  tempVec2.Multiply('N', 'N', 1.0, invAmat, tempVec1, 0.0);
  elevec.Multiply('N', 'N', -1.0, Imat, tempVec2, 1.0);

  return;
} // ComputeResidual

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeFaceMatrices(
    const int face, double dt,
    int indexstart, double c)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeFaceMatrices");

  // Compute the matrices C, E and G
  // Here, we don't consider material properties! Don't forget this during condensation

  double tau = 1.0 / c ;

  // loop over number of shape functions
  for (unsigned int q = 0; q < ndofs_; ++q)
  {
    if( shapesface_->shfunctI.NonzeroOnFace(q) )
    {
      // loop over number of shape functions
      for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
      {
        // C and E

        // numerical integration: sum over quadrature points
        double tempE = 0.0;
        double temp_d[nsd_] = {0.0};
        for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
        {
          double temp = shapesface_->jfac(i) * shapesface_->shfunct(p, i)
                                             * shapesface_->shfunctI(q, i);
          tempE += temp;
          for (unsigned int j = 0; j < nsd_; ++j)
            temp_d[j] += temp * shapesface_->normals(j, i);
        }
        for (unsigned int j = 0; j < nsd_; ++j)
        {
          Cmat(j * ndofs_ + q, indexstart + p) = temp_d[j];
          Imat(indexstart + p, j * ndofs_ + q) = temp_d[j];
        }
        Emat(q, indexstart + p) = tempE;
        Jmat(indexstart + p, q) = tempE;
      }
    } // for (unsigned int q=0; q<ndofs_; ++q)
  } // for (unsigned int p=0; p<nfdofs_; ++p)

  // G

  // loop over number of shape functions
  for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
  {
    // loop over number of shape functions
    for (unsigned int q = 0; q <= p; ++q)
    {
      double tempG = 0.0;
      for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
      {
        tempG += shapesface_->jfac(i) * shapesface_->shfunct(p, i) * shapesface_->shfunct(q, i);
      }
      Gmat(indexstart + p, indexstart + q) = Gmat(indexstart + q,indexstart + p) -= tau * tempG;

    } // for (unsigned int q=0; q<nfdofs_; ++q)
  } // for (unsigned int p=0; p<nfdofs_; ++p)

  // one term is still missing in D!!
  for (unsigned int p = 0; p < ndofs_; ++p)
  {
    for (unsigned int q = 0; q <= p; ++q)
    {
      double tempD = 0.0;
      if( shapesface_->shfunctI.NonzeroOnFace(p) && shapesface_->shfunctI.NonzeroOnFace(q) )
      {
        for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
          tempD += shapesface_->jfac(i) * shapesface_->shfunctI(p, i) * shapesface_->shfunctI(q, i);
        Dmat(p, q) = Dmat(q, p) += tau * tempD;
      }
    }
  }

  return;
} // ComputeFaceMatrices

/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::CondenseLocalPart(
    Epetra_SerialDenseMatrix &eleMat)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::CondenseLocalPart");

  /*
   THE MATRIX
                             -1
                   +--------+    +-----+
                   |        |    |     |
                   | A   B  |    |  C  |
    K = - [ I  J ] |        |    |     | + G  = - I x - J y + G
                   | H  D+M |    |  E  |
                   +--------+    +-----+

                      -1     -1           -1
       y = [ D - H  A   B  ]   [ E - H  A   C  ]

            -1
       x = A   [  C - B  y ]

   */

  int onfdofs = eleMat.M();

  Epetra_SerialDenseMatrix tempMat1(ndofs_, ndofs_ * nsd_);
  tempMat1.Multiply('N', 'N', 1.0, Hmat, invAmat, 0.0); // =  H A^{-1}

  Epetra_SerialDenseMatrix tempMat2(ndofs_, ndofs_);

  // for nonlinear solves, this is not D+M but D+M(P^i)+dMdPP^i
  tempMat2 = Dmat;
  tempMat2.Scale(1.0);
  tempMat2 += Mmat;

  tempMat2.Multiply('N', 'N', -1.0, tempMat1, Bmat, 1.0); // = D - H A^{-1} B
  Epetra_SerialDenseMatrix tempMat3(ndofs_, onfdofs);
  tempMat3 = Emat;
  tempMat3.Multiply('N', 'N', -1.0, tempMat1, Cmat, 1.0); // = E - H A^{-1} C

  {
    Epetra_SerialDenseSolver inverseinW;
    inverseinW.SetMatrix(tempMat2);
    int err = inverseinW.Invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d", err);
  }
  // tempMat2 = (  D - H A^{-1} B )^{-1}

  eleMat = Gmat; // = G
  tempMat1.Shape(ndofs_, onfdofs);
  tempMat1.Multiply('N', 'N', 1.0, tempMat2, tempMat3, 0.0); // = y
  eleMat.Multiply('N', 'N', -1.0, Jmat, tempMat1, 1.0); // = - J y + G

  tempMat2.Shape(ndofs_ * nsd_, onfdofs);
  tempMat2 = Cmat;
  tempMat2.Multiply('N', 'N', -1.0, Bmat, tempMat1, 1.0); // = - C^T + B^T W

  tempMat3.Shape(ndofs_ * nsd_, onfdofs);
  tempMat3.Multiply('N', 'N', 1.0, invAmat, tempMat2, 0.0); // = x = A^{-1} ( C - B y )

  eleMat.Multiply('N', 'N', -1.0, Imat, tempMat3, 1.0); // = K = G - I x - J y

  return;
} // CondenseLocalPart

/*----------------------------------------------------------------------*
 * Compute internal and face matrices
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeMatrices(
    DRT::Discretization & discretization,
    const Teuchos::RCP<MAT::Material> &mat,
    DRT::ELEMENTS::Acou & ele,
    double dt,
    INPAR::ACOU::DynamicType dyna)
{
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density(ele.Id());
  double c = actmat->SpeedofSound(ele.Id());

  // init face matrices
  zeroMatrix(Cmat);
  zeroMatrix(Emat);
  zeroMatrix(Gmat);
  zeroMatrix(Imat);
  zeroMatrix(Jmat);

  double tau = 1.0 / c /rho;//double tau = 1.0 / c ;

  ComputeInteriorMatrices(dt, rho, c);

  int sumindex = 0;
  for (unsigned int face = 0; face < nfaces_; ++face)
  {
    DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele.Faces()[face]->Degree(),
        shapes_.usescompletepoly_, 2 * ele.Faces()[face]->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(ele, face);

    ComputeFaceMatrices(face, dt, sumindex, c*rho);
    sumindex += shapesface_->nfdofs_;
  }

  //Imat.Scale(rho);
  Emat.Scale(-tau);
  Jmat.Scale(tau);

  return;
}

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeADERResidual(
    DRT::ELEMENTS::Acou*       ele,
    DRT::Discretization &      discretization,
    Teuchos::ParameterList&    params,
    Epetra_SerialDenseVector & elevec,
    Epetra_SerialDenseVector & interiorVeln,
    Epetra_SerialDenseVector & interiorPressn,
    double                     c,
    double                     rho)
{
  // we will need the inverse of M
  Epetra_SerialDenseMatrix invMmat(shapes_.ndofs_,shapes_.ndofs_);
  invMmat += Mmat;
  {
    Epetra_SerialDenseSolver inverseMmat;
    inverseMmat.SetMatrix(invMmat);
    inverseMmat.Invert();
  }


  double dt = params.get<double>("dt");

  const unsigned int n_q_points = shapes_.nqpoints_;
  std::vector<Epetra_SerialDenseVector> vcontrib(n_q_points,Epetra_SerialDenseVector(nsd_));
  std::vector<double>                   pcontrib(n_q_points,0.0);
  std::vector<Epetra_SerialDenseVector> tempcontrib(n_q_points,Epetra_SerialDenseVector(nsd_+1));

  for(unsigned int q=0; q<n_q_points; ++q)
  {
    // get the values of pressure and velocity at this integration point
    Epetra_SerialDenseVector v_and_p(nsd_+1);
    for(unsigned int i=0; i<shapes_.ndofs_; ++i)
    {
      v_and_p[nsd_] += shapes_.shfunct(i,q)*interiorPressn(i);
      for(unsigned int d=0; d<nsd_; ++d)
        v_and_p[d]  += shapes_.shfunct(i,q)*interiorVeln(i+d*shapes_.ndofs_);
    }

    // get pressure and velocity gradient at this integration point
    Epetra_SerialDenseMatrix v_and_p_grad(nsd_+1,nsd_);
    for(unsigned int i=0; i<shapes_.ndofs_; ++i)
    {
      for(unsigned int d=0; d<nsd_; ++d)
      {
        v_and_p_grad(nsd_,d) += shapes_.shderxy(i*nsd_+d,q) * interiorPressn(i);
        for(unsigned int e=0; e<nsd_; ++e)
          v_and_p_grad(e,d) += shapes_.shderxy(i*nsd_+d,q) * interiorVeln(i+e*shapes_.ndofs_);
      }
    }

    Epetra_SerialDenseVector postpressure_val(nsd_+1);
    for(unsigned int i=0; i<shapes_.ndofs_; ++i)
    {
      postpressure_val[nsd_] += shapes_.shfunct(i,q)*ele->eleADERimprovedDiv_(i);
      for(unsigned int d=0; d<nsd_; ++d)
        postpressure_val[d]  += shapes_.shfunct(i,q)*ele->eleADERimprovedGrad_(i+d*shapes_.ndofs_);
    }


    // contribution from k=0
    for(unsigned int d=0; d<nsd_; ++d)
      vcontrib[q][d] = dt * v_and_p[d];
    pcontrib[q] = dt * v_and_p[nsd_];

    // contribution from k=1
    /*for(unsigned int d=0; d<nsd_; ++d)
    {
      vcontrib[q][d] -= dt*dt/2.0/rho*v_and_p_grad(nsd_,d);
      pcontrib[q] -= dt*dt/2.0*c*c*rho*v_and_p_grad(d,d);
    }*/
    for(unsigned int d=0; d<nsd_; ++d)
      vcontrib[q][d] -= dt*dt/2.0/rho*postpressure_val(d);
    pcontrib[q] -= dt*dt/2.0*c*c*rho*postpressure_val(nsd_);


    // evaluate phi_1
    /*for(unsigned int d=0; d<nsd_; ++d)
    {
      tempcontrib[q][d] = 1.0/rho*v_and_p_grad(nsd_,d);
      tempcontrib[q][nsd_] += c*c*rho*v_and_p_grad(d,d);
    }*/
    for(unsigned int d=0; d<nsd_+1; ++d)
      tempcontrib[q][d] = postpressure_val(d);
  }

  double fac = -dt*dt/2.0;
  for(unsigned int k=2; k<=shapes_.degree_+1; ++k)
  {
    fac *= -dt/(k+1);

    // integrate over element and apply inverse mass matrix
    Epetra_SerialDenseVector temp_v(shapes_.ndofs_*nsd_);
    Epetra_SerialDenseVector temp_p(shapes_.ndofs_);
    {
      for(unsigned int i=0; i<shapes_.ndofs_; ++i)
      {
        for(unsigned int q=0; q<n_q_points; ++q)
        {
          temp_p(i) += shapes_.shfunct(i,q) * tempcontrib[q][nsd_] * shapes_.jfac(q);
          for(unsigned int d=0; d<nsd_; ++d)
            temp_v(i+d*shapes_.ndofs_) += shapes_.shfunct(i,q) * tempcontrib[q][d] * shapes_.jfac(q);
        }
      }
      for(unsigned int s=0;s<tempcontrib.size();++s)
        for(unsigned int e=0;e<nsd_+1;++e)
          tempcontrib[s][e]=0.0;

      // apply inverse mass matrix to temp vectors
      Epetra_SerialDenseVector temp_vv(nsd_*shapes_.ndofs_);
      Epetra_SerialDenseVector temp_pp(shapes_.ndofs_);
      temp_vv+=temp_v;
      temp_pp+=temp_p;
      temp_v.Multiply('N','N',rho/dt,invAmat,temp_vv,0.0);
      temp_p.Multiply('N','N',1.0/c/c/rho/dt,invMmat,temp_pp,0.0);
    }

    for(unsigned int q=0; q<n_q_points; ++q)
    {
      // get the gauss point values
      Epetra_SerialDenseMatrix v_and_p_grad(nsd_+1,nsd_);
      for(unsigned int i=0; i<shapes_.ndofs_; ++i)
      {
        for(unsigned int d=0; d<nsd_; ++d)
        {
          v_and_p_grad(nsd_,d) += shapes_.shderxy(i*nsd_+d,q) * temp_p(i);
          for(unsigned int e=0; e<nsd_; ++e)
            v_and_p_grad(e,d) += shapes_.shderxy(i*nsd_+d,q) * temp_v(i+e*shapes_.ndofs_);
        }
      }

      // calculate contributions
      for(unsigned int d=0; d<nsd_; ++d)
      {
        vcontrib[q][d] += fac/rho*v_and_p_grad(nsd_,d);
        pcontrib[q] += fac*c*c*rho*v_and_p_grad(d,d);
      }

      // evaluate things phi_k+1 needs
      for(unsigned int d=0; d<nsd_; ++d)
      {
        tempcontrib[q][d] = 1.0/rho*v_and_p_grad(nsd_,d);
        tempcontrib[q][nsd_] += c*c*rho*v_and_p_grad(d,d);
      }
    }
  }

  // submit what we collected to the field!
  Epetra_SerialDenseVector temp_v(shapes_.ndofs_*nsd_);
  Epetra_SerialDenseVector temp_p(shapes_.ndofs_);
  for(unsigned int i=0; i<shapes_.ndofs_; ++i)
  {
    for(unsigned int q=0; q<n_q_points; ++q)
    {
      temp_p(i) += shapes_.shfunct(i,q) * pcontrib[q] * shapes_.jfac(q);
      for(unsigned int d=0; d<nsd_; ++d)
        temp_v(i+d*shapes_.ndofs_) += shapes_.shfunct(i,q) * vcontrib[q][d] * shapes_.jfac(q);
    }
  }
  Epetra_SerialDenseVector temp_vv(nsd_*shapes_.ndofs_);
  Epetra_SerialDenseVector temp_pp(shapes_.ndofs_);
  temp_vv+=temp_v;
  temp_pp+=temp_p;
  temp_v.Multiply('N','N',rho/dt,invAmat,temp_vv,0.0);
  temp_p.Multiply('N','N',1.0/c/c/rho/dt,invMmat,temp_pp,0.0);

  // store these values for the second update
  //{
    // sort this back to the interior values vector
    std::vector<double> interiorValnp(shapes_.ndofs_ * (nsd_ + 1));
    for (unsigned int i = 0; i < interiorValnp.size(); ++i)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        interiorValnp[i] = temp_p((i + 1) / (nsd_ + 1) - 1);
      else
      {
        int xyz = i % (nsd_ + 1); // 0 for x, 1 for y and 2 for z (for 3D)
        interiorValnp[i] = temp_v(xyz * shapes_.ndofs_ + i / (nsd_ + 1));
      }
    }

    // tell this change in the interior variables the discretization
    std::vector<int> localDofs = discretization.Dof(1, ele);
    const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
    {
      Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"tempsrc");
      Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
      for (unsigned int i = 0; i < localDofs.size(); ++i)
      {
        const int lid = intdofcolmap->LID(localDofs[i]);
        secondary[lid] = interiorValnp[i];
      }
    }
  //}

  // compute residual for trace system
  elevec.Multiply('N', 'N', 1.0, Imat, temp_v, 0.0);
  elevec.Multiply('N', 'N', 1.0, Jmat, temp_p, 1.0);

  return;
}

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::UpdateADERSolution(
    DRT::Discretization &             discretization,
    Teuchos::ParameterList& params,
    DRT::ELEMENTS::Acou*       ele, const Teuchos::RCP<MAT::Material> &mat,
    Epetra_SerialDenseVector & elevec, double dt)
{
  // we will need the inverse of M
  Epetra_SerialDenseMatrix invMmat(shapes_->ndofs_,shapes_->ndofs_);
  invMmat += localSolver_->Mmat;
  {
    Epetra_SerialDenseSolver inverseMmat;
    inverseMmat.SetMatrix(invMmat);
    inverseMmat.Invert();
  }
  // read tempsrc
  Epetra_SerialDenseVector tempsrc_v(shapes_->ndofs_*nsd_);
  Epetra_SerialDenseVector tempsrc_p(shapes_->ndofs_);
  {
    std::vector<double> interiorValnp(shapes_->ndofs_ * (nsd_ + 1));
    Teuchos::RCP<const Epetra_Vector> intvel = discretization.GetState(1,"tempsrc");
    std::vector<int> localDofs1 = discretization.Dof(1, ele);
    DRT::UTILS::ExtractMyValues(*intvel, interiorValnp, localDofs1);

    // now write this in corresponding interiorVelnp_ and interiorPressnp_
    for (unsigned int i = 0; i < interiorValnp.size(); ++i)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        tempsrc_p((i + 1) / (nsd_ + 1) - 1) = interiorValnp[i];
      else
      {
        int xyz = i % (nsd_ + 1);
        tempsrc_v(xyz * shapes_->ndofs_ + i / (nsd_ + 1)) = interiorValnp[i];
      }
    }
  }

  // required material information
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();

  // trace value
  int onfdofs = elevec.M();
  Epetra_SerialDenseVector traceVal_SDV(onfdofs);
  for (int i = 0; i < onfdofs; ++i)
    traceVal_SDV(i) = traceVal_[i];


  Epetra_SerialDenseVector tempvec_v(shapes_->ndofs_*nsd_);
  Epetra_SerialDenseVector tempvec_p(shapes_->ndofs_);

  tempvec_v.Multiply('N','N',-1.0/rho,localSolver_->Cmat,traceVal_SDV,0.0);
  tempvec_p.Multiply('N','N',-1.0,localSolver_->Emat,traceVal_SDV,0.0);

  tempvec_v.Multiply('N','N',1.0/rho,localSolver_->Bmat,tempsrc_p,1.0);
  tempvec_p.Multiply('N','N',rho*c*c,localSolver_->Hmat,tempsrc_v,1.0);
  tempvec_p.Multiply('N','N',1.0,localSolver_->Dmat,tempsrc_p,1.0);

  // apply inverse mass matrices
  Epetra_SerialDenseVector temp_vv(nsd_*shapes_->ndofs_);
  Epetra_SerialDenseVector temp_pp(shapes_->ndofs_);
  temp_vv+=tempvec_v;
  temp_pp+=tempvec_p;
  tempvec_v.Multiply('N','N',-rho/dt,localSolver_->invAmat,temp_vv,0.0);
  tempvec_p.Multiply('N','N',-1.0/c/c/rho/dt,invMmat,temp_pp,0.0);

  // update interior variables
  ele->eleinteriorVelnp_ += tempvec_v;
  ele->eleinteriorPressnp_ += tempvec_p;

  return;
}


// template classes
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tri3>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::nurbs27>;
