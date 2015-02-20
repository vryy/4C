/*--------------------------------------------------------------------------*/
/*!
\file acou_sol_ele_calc.cpp
\brief

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "acou_sol_ele_calc.H"
#include "acou_ele_calc.H"
#include "acou_ele_action.H"
#include "acou_utils.H"

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
  void zeroMatrix (Epetra_SerialDenseMatrix &mat)
  {
    std::memset(mat.A(), 0, sizeof(double)*mat.M()*mat.N());
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
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouSolEleCalc<distype>::AcouSolEleCalc()
{}

/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouSolEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou*    ele,
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

/*----------------------------------------------------------------------*
 * Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouSolEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou* ele,
                                                 DRT::Discretization    &   discretization,
                                                 const std::vector<int> &   lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elemat1,
                                                 Epetra_SerialDenseMatrix&  ,
                                                 Epetra_SerialDenseVector&  elevec1,
                                                 Epetra_SerialDenseVector&  elevec2,
                                                 Epetra_SerialDenseVector&  ,
                                                 bool                       offdiag)
{
  // check if this is an hdg element and init completepoly
  if (const DRT::ELEMENTS::AcouSol * hdgele = dynamic_cast<const DRT::ELEMENTS::AcouSol*>(ele))
    usescompletepoly_ = hdgele->UsesCompletePolynomialSpace();
  else
    dserror("cannot cast element to acousol element");
  InitializeShapes(ele);

  const ACOU::Action action = DRT::INPUT::get<ACOU::Action>(params,"action");
  bool updateonly = false;
  shapes_->Evaluate(*ele);

  switch(action)
  {
  case ACOU::project_field:
  {
    if(mat->MaterialType()!=INPAR::MAT::m_acousticsolmat) dserror("for physical type 'solous' please supply MAT_AcousticSol");
    ElementInit(ele,params);
    ProjectField(ele,params,elevec1,elevec2);
    break;
  }
  case ACOU::project_dirich_field:
  {
    if(mat->MaterialType()!=INPAR::MAT::m_acousticsolmat) dserror("for physical type 'solous' please supply MAT_AcousticSol");
    ElementInit(ele,params);
    ProjectDirichField(ele,params,elevec1);
    break;
  }
  case ACOU::project_optical_field:
  {
    if(mat->MaterialType()!=INPAR::MAT::m_acousticsolmat) dserror("for physical type 'solous' please supply MAT_AcousticSol");
    ProjectOpticalField(ele,params,elevec2);
    break;
  }
  case ACOU::ele_init:
  {
    ElementInit(ele,params);
    break;
  }
  case ACOU::fill_restart_vecs:
  {
    ReadGlobalVectors(ele,discretization,lm);
    FillRestartVectors(ele,discretization);
    break;
  }
  case ACOU::ele_init_from_restart:
  {
    ElementInit(ele,params);
    ElementInitFromRestart(ele,discretization);
    break;
  }
  case ACOU::interpolate_hdg_to_node:
  {
    ReadGlobalVectors(ele,discretization,lm);
    NodeBasedValues(mat,ele,elevec1);
    break;
  }
  case ACOU::interpolate_psi_to_node:
  {
    double dt = params.get<double>("dt");
    ReadGlobalVectors(ele,discretization,lm);
    NodeBasedPsi(mat,ele,elevec1,dt);
    break;
  }
  case ACOU::calc_acou_error:
  {
    ReadGlobalVectors(ele,discretization,lm);
    ComputeError(ele,params,elevec1);
    break;
  }
  case ACOU::calc_abc:
  {
    int face = params.get<int>("face");
    shapesface_->EvaluateFace(*ele,face);
    // note: absorbing bcs are treated fully implicit!
    localSolver_->ComputeAbsorbingBC(ele,params,mat,face,elemat1,elevec1);
    break;
  }
  case ACOU::calc_pressuremon:
  {
    int face = params.get<int>("face");
    double dt = params.get<double>("dt");

    ComputeMatrices(mat,*ele,dt,dyna_,true);

    shapesface_->EvaluateFace(*ele,face);

    if(!params.isParameter("nodeindices"))
      localSolver_->ComputeSourcePressureMonitor(ele,params,face,elevec1);
    else
      localSolver_->ComputeSourcePressureMonitorLine3D(ele,params,mat,face,elemat1,elevec1);

    break;
  }
  case ACOU::calc_pmon_nodevals:
  {
    int face = params.get<int>("face");
    ReadGlobalVectors(ele,discretization,lm);
    shapesface_->EvaluateFace(*ele,face);
    ComputePMonNodeVals(ele,face,elevec1);
    break;
  }
  case ACOU::calc_systemmat_and_residual:
  {
    const bool resonly = params.get<bool>("resonly");
    bool adjoint = params.get<bool>("adjoint");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(ele,discretization,lm);
    zeroMatrix(elevec1);
    ComputeMatrices(mat,*ele,dt,dyna_,adjoint);

    if(!resonly)
    {
      localSolver_->CondenseLocalPart(elemat1,dyna_);
    }

    VectorHandling(ele,params,dt);
    localSolver_->ComputeResidual(params,elevec1,interiorGradVelnp_,interiorVelnp_,interiorPressnp_,traceVal_,dyna_);

   break;
  }
  case ACOU::update_secondary_solution:
    updateonly = true; // no break here!!!
  case ACOU::update_secondary_solution_and_calc_residual:
  {
    bool adjoint = params.get<bool>("adjoint");
    bool errormaps = params.get<bool>("errormaps");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(ele,discretization,lm);

    zeroMatrix(elevec1);
    ComputeMatrices(mat,*ele,dt,dyna_,adjoint);

    UpdateInteriorVariablesAndComputeResidual(params,ele,elevec1,dt,errormaps,updateonly);

    break;
  }
  default:
    dserror("unknown action supplied %d",action);
    break;
  } // switch(action)

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::InitializeShapes(const DRT::ELEMENTS::Acou* ele)
{

  if (shapes_ == Teuchos::null )
    shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(ele->Degree(),
        usescompletepoly_,
        2*ele->Degree()));

  else if (shapes_->degree_ != unsigned(ele->Degree()) || shapes_->usescompletepoly_ != usescompletepoly_)
    shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(ele->Degree(),
        usescompletepoly_,
        2*ele->Degree()));

  if (shapesface_ == Teuchos::null)
  {
    DRT::UTILS::ShapeValuesFaceParams svfparams( ele->Degree(), usescompletepoly_, 2 * ele->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
  }
  // TODO: check distype

  localSolver_ = Teuchos::rcp(new LocalSolver(ele,*shapes_,*shapesface_,usescompletepoly_));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouSolEleCalc<distype> *
DRT::ELEMENTS::AcouSolEleCalc<distype>::Instance( bool create )
{
  static AcouSolEleCalc<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new AcouSolEleCalc<distype>();
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
void DRT::ELEMENTS::AcouSolEleCalc<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}

/*----------------------------------------------------------------------*
 * Constructor LocalSolver
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::
LocalSolver(const DRT::ELEMENTS::Acou* ele,
            const DRT::UTILS::ShapeValues<distype> &shapeValues,
            DRT::UTILS::ShapeValuesFace<distype> &shapeValuesFace,
            bool completepoly)
:
ndofs_ (shapeValues.ndofs_),
shapes_(shapeValues),
shapesface_(shapeValuesFace)
{
  int onfdofs = 0;
  for(unsigned int i=0; i<nfaces_; ++i)
  {
    shapesface_.EvaluateFace(*ele,i);
    onfdofs += shapesface_.nfdofs_;
  }
  onfdofs *= nsd_;

  amat.Shape(nsd_*nsd_*ndofs_,nsd_*nsd_*ndofs_);
  invamat.Shape(nsd_*nsd_*ndofs_,nsd_*nsd_*ndofs_);
  bmat.Shape(nsd_*nsd_*ndofs_,nsd_*ndofs_);
  cmat.Shape(nsd_*nsd_*ndofs_,onfdofs);
  dmat.Shape(nsd_*ndofs_,nsd_*nsd_*ndofs_);
  emat.Shape(nsd_*ndofs_,nsd_*ndofs_);
  ehatmat.Shape(nsd_*ndofs_,nsd_*ndofs_);
  fmat.Shape(nsd_*ndofs_,ndofs_);
  gmat.Shape(nsd_*ndofs_,onfdofs);
  hmat.Shape(ndofs_,nsd_*ndofs_);
  imat.Shape(ndofs_,ndofs_);
  jmat.Shape(ndofs_,onfdofs);
  kmat.Shape(onfdofs,nsd_*nsd_*ndofs_);
  lmat.Shape(onfdofs,nsd_*ndofs_);
  mmat.Shape(onfdofs,ndofs_);
  nmat.Shape(onfdofs,onfdofs);
}


/*----------------------------------------------------------------------*
 * Element init // TODO: implement in AcouEleInterface
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::ElementInit(DRT::ELEMENTS::Acou* ele,
                                                      Teuchos::ParameterList& params)
{
  DRT::ELEMENTS::AcouSol * solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);

  // each element has to store the interior vectors by itseld, p-adaptivity or not
  // so, shape it, as you need it
  if(params.isParameter("dynamic type")) dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");
  else dyna_ = INPAR::ACOU::acou_impleuler;
  solele->eleinteriorGradVelnp_.Shape(shapes_->ndofs_ * nsd_* nsd_, 1);
  solele->eleinteriorVelnp_.Shape(shapes_->ndofs_ * nsd_, 1);
  solele->eleinteriorPressnp_.Shape(shapes_->ndofs_, 1);
  solele->eleinteriorGradVeln_.Shape(shapes_->ndofs_ * nsd_* nsd_, 1);
  solele->eleinteriorVeln_.Shape(shapes_->ndofs_ * nsd_, 1);
  solele->eleinteriorPressn_.Shape(shapes_->ndofs_, 1);
  switch (dyna_)
  {
  case INPAR::ACOU::acou_bdf4:
    solele->eleinteriorGradVelnmmm_.Shape(shapes_->ndofs_ * nsd_* nsd_, 1);
    solele->eleinteriorVelnmmm_.Shape(shapes_->ndofs_ * nsd_, 1);
    solele->eleinteriorPressnmmm_.Shape(shapes_->ndofs_, 1); // no break here!
  case INPAR::ACOU::acou_bdf3:
    solele->eleinteriorGradVelnmm_.Shape(shapes_->ndofs_ * nsd_* nsd_, 1);
    solele->eleinteriorVelnmm_.Shape(shapes_->ndofs_ * nsd_, 1);
    solele->eleinteriorPressnmm_.Shape(shapes_->ndofs_, 1); // no break here!
  case INPAR::ACOU::acou_bdf2:
    solele->eleinteriorGradVelnm_.Shape(shapes_->ndofs_ * nsd_* nsd_, 1);
    solele->eleinteriorVelnm_.Shape(shapes_->ndofs_ * nsd_, 1);
    solele->eleinteriorPressnm_.Shape(shapes_->ndofs_, 1);
    break; // here you go!
  case INPAR::ACOU::acou_dirk23:
    solele->elesp_.resize(2);
    solele->eleyp_.resize(2);
    solele->elefp_.resize(2);
    solele->elesv_.resize(2);
    solele->eleyv_.resize(2);
    solele->elefv_.resize(2);
    solele->elesgv_.resize(2);
    solele->eleygv_.resize(2);
    solele->elefgv_.resize(2);
    break;
  case INPAR::ACOU::acou_dirk33:
  case INPAR::ACOU::acou_dirk34:
    solele->elesv_.resize(3);
    solele->eleyv_.resize(3);
    solele->elefv_.resize(3);
    solele->elesgv_.resize(3);
    solele->eleygv_.resize(3);
    solele->elefgv_.resize(3);
    solele->elesp_.resize(3);
    solele->eleyp_.resize(3);
    solele->elefp_.resize(3);
    break;
  case INPAR::ACOU::acou_dirk54:
    solele->elesv_.resize(5);
    solele->eleyv_.resize(5);
    solele->elefv_.resize(5);
    solele->elesgv_.resize(5);
    solele->eleygv_.resize(5);
    solele->elefgv_.resize(5);
    solele->elesp_.resize(5);
    solele->eleyp_.resize(5);
    solele->elefp_.resize(5);
    break;
  default:
    break; // do nothing for impl and trap
  }
  for (unsigned int i = 0; i < solele->elesp_.size(); ++i) // zero size if not dirk
  {
    solele->elesp_[i].Shape(shapes_->ndofs_, 1);
    solele->eleyp_[i].Shape(shapes_->ndofs_, 1);
    solele->elefp_[i].Shape(shapes_->ndofs_, 1);
    solele->elesv_[i].Shape(shapes_->ndofs_ * nsd_, 1);
    solele->eleyv_[i].Shape(shapes_->ndofs_ * nsd_, 1);
    solele->elefv_[i].Shape(shapes_->ndofs_ * nsd_, 1);
    solele->elesgv_[i].Shape(shapes_->ndofs_ * nsd_* nsd_, 1);
    solele->eleygv_[i].Shape(shapes_->ndofs_ * nsd_* nsd_, 1);
    solele->elefgv_[i].Shape(shapes_->ndofs_ * nsd_* nsd_, 1);
  }

  bool padaptivity = false;
  if(params.isParameter("padaptivity")) padaptivity = params.get<bool>("padaptivity");
  if (padaptivity)
    solele->elenodeTrace_.Shape(nen_, 1);

  return;
} // ElementInit


/*----------------------------------------------------------------------*
 * VectorHandling
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::
VectorHandling(DRT::ELEMENTS::Acou   * ele,
                   Teuchos::ParameterList&       params,
                   double dt)
{
  if( dyna_ == INPAR::ACOU::acou_dirk23 ||
      dyna_ == INPAR::ACOU::acou_dirk33 ||
      dyna_ == INPAR::ACOU::acou_dirk34 ||
      dyna_ == INPAR::ACOU::acou_dirk54 )
  {
    // when we are here, the time integrator is dirk, and we have to do the update things
    DRT::ELEMENTS::AcouSol * solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);

    double dirk_a[6][6];
    double dirk_b[6];
    double dirk_c[6];
    int dirk_q;
    ACOU::FillDIRKValues(dyna_, dirk_a, dirk_b, dirk_c, dirk_q);
    int stage = params.get<int>("stage");

    if(stage == 0)
    {
      solele->eleinteriorPressnp_ = solele->eleinteriorPressn_;
      solele->eleinteriorVelnp_ = solele->eleinteriorVeln_;
      solele->eleinteriorGradVelnp_ = solele->eleinteriorGradVeln_;
    }

    Epetra_SerialDenseVector tempVecp1(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecp2(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecv1(shapes_->ndofs_ * nsd_);
    Epetra_SerialDenseVector tempVecv2(shapes_->ndofs_ * nsd_);
    Epetra_SerialDenseVector tempVecgv1(shapes_->ndofs_ * nsd_* nsd_);
    Epetra_SerialDenseVector tempVecgv2(shapes_->ndofs_ * nsd_* nsd_);

    solele->elesp_[stage] = solele->eleinteriorPressn_;
    solele->elesp_[stage].Scale(1.0 / dt); // dt includes a_ii
    solele->elesv_[stage] = solele->eleinteriorVeln_;
    solele->elesv_[stage].Scale(1.0 / dt); // dt includes a_ii
    solele->elesgv_[stage] = solele->eleinteriorGradVeln_;
    solele->elesgv_[stage].Scale(1.0 / dt); // dt includes a_ii

    for(int i=0; i<stage; ++i)
    {
      tempVecp1 = solele->elesp_[i];
      tempVecp1.Scale(-1.0);
      tempVecp2 = solele->eleyp_[i];
      tempVecp2.Scale(1.0 / dt);
      tempVecp1 += tempVecp2;
      tempVecp1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      solele->elesp_[stage] += tempVecp1;
      tempVecv1 = solele->elesv_[i];
      tempVecv1.Scale(-1.0);
      tempVecv2 = solele->eleyv_[i];
      tempVecv2.Scale(1.0 / dt);
      tempVecv1 += tempVecv2;
      tempVecv1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      solele->elesv_[stage] += tempVecv1;
      tempVecgv1 = solele->elesgv_[i];
      tempVecgv1.Scale(-1.0);
      tempVecgv2 = solele->eleygv_[i];
      tempVecgv2.Scale(1.0 / dt);
      tempVecgv1 += tempVecgv2;
      tempVecgv1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      solele->elesgv_[stage] += tempVecgv1;
    }

    // these vectors s are used for the calculation of the residual -> write them to used local solver variable
    interiorPressnp_ = solele->elesp_[stage];
    interiorVelnp_ = solele->elesv_[stage];
    interiorGradVelnp_ = solele->elesgv_[stage];
    interiorPressnp_.Scale(dt);
    interiorVelnp_.Scale(dt);
    interiorGradVelnp_.Scale(dt);
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf2)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_ * nsd_; ++i)
      interiorGradVelnp_[i] = 4.0 / 3.0 * interiorGradVelnp_[i] - 1.0 / 3.0 * interiorGradVelnm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 4.0 / 3.0 * interiorVelnp_[i] - 1.0 / 3.0 * interiorVelnm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 4.0 / 3.0 * interiorPressnp_[i] - 1.0 / 3.0 * interiorPressnm_[i];
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf3)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_ * nsd_; ++i)
      interiorGradVelnp_[i] = 18.0 / 11.0 * interiorGradVelnp_[i] - 9.0 / 11.0 * interiorGradVelnm_[i] + 2.0 / 11.0 * interiorGradVelnmm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 18.0 / 11.0 * interiorVelnp_[i] - 9.0 / 11.0 * interiorVelnm_[i] + 2.0 / 11.0 * interiorVelnmm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 18.0 / 11.0 * interiorPressnp_[i] - 9.0 / 11.0 * interiorPressnm_[i] + 2.0 / 11.0 * interiorPressnmm_[i];
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf4)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_ * nsd_; ++i)
      interiorGradVelnp_[i] = 48.0 / 25.0 * interiorGradVelnp_[i] - 36.0 / 25.0 * interiorGradVelnm_[i] + 16.0 / 25.0 * interiorGradVelnmm_[i] - 3.0 / 25.0 * interiorGradVelnmmm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 48.0 / 25.0 * interiorVelnp_[i] - 36.0 / 25.0 * interiorVelnm_[i] + 16.0 / 25.0 * interiorVelnmm_[i] - 3.0 / 25.0 * interiorVelnmmm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 48.0 / 25.0 * interiorPressnp_[i] - 36.0 / 25.0 * interiorPressnm_[i] + 16.0 / 25.0 * interiorPressnmm_[i] - 3.0 / 25.0 * interiorPressnmmm_[i];
  }

  return;
} // VectorHandling

/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouSolEleCalc<distype>::ProjectField(
    DRT::ELEMENTS::Acou*                 ele,
    Teuchos::ParameterList&              params,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2)
{
  DRT::ELEMENTS::AcouSol * solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 || unsigned(elevec2.M()) == (nsd_*nsd_+nsd_+1)*shapes_->ndofs_, "Wrong size in project vector 2");

  // get function
  const int *start_func = params.getPtr<int>("funct");

  // internal variables
  if (elevec2.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(shapes_->nqpoints_,nsd_*nsd_+nsd_+1);
    Epetra_SerialDenseMatrix massPart(shapes_->ndofs_,shapes_->nqpoints_);

    for (unsigned int q=0; q<shapes_->nqpoints_; ++q )
    {
      const double fac = shapes_->jfac(q);
      const double sqrtfac = std::sqrt(fac);
      double xyz[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        xyz[d] = shapes_->xyzreal(d,q); // coordinates of quadrature point in real coordinates
      double p;
      double gradient[nsd_];
      double velgrad[nsd_*nsd_];

      dsassert(start_func != NULL,"startfuncno not set for initial value");
      EvaluateAll(*start_func, xyz,  p, gradient, velgrad, 0.0); // p at quadrature point

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i=0; i<shapes_->ndofs_; ++i)
      {
        massPart(i,q) = shapes_->shfunct(i,q) * sqrtfac;
        localMat(i,nsd_) += shapes_->shfunct(i,q) * p * fac;
        for(unsigned int j=0; j<nsd_; ++j)
          localMat(i,j) += shapes_->shfunct(i,q) * gradient[j] * fac;
        for(unsigned int j=0; j<nsd_*nsd_; ++j)
          localMat(i,nsd_+1+j) += shapes_->shfunct(i,q) * velgrad[j] * fac;
      }
    }

    Epetra_SerialDenseMatrix massMat(shapes_->ndofs_,shapes_->ndofs_);
    massMat.Multiply('N','T',1.,massPart,massPart,0.);
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(massMat);
      inverseMass.SetVectors(localMat,localMat);
      inverseMass.Solve();
    }

    for (unsigned int r = 0; r < shapes_->ndofs_; ++r)
    {
      solele->eleinteriorPressnp_(r) += localMat(r, nsd_); // pressure
      for (unsigned int i=0;i<nsd_;++i)
        solele->eleinteriorVelnp_(i*shapes_->ndofs_+r) += localMat(r, i); // velocity
      for (unsigned int j=0;j<nsd_*nsd_;++j)
        solele->eleinteriorGradVelnp_(j*shapes_->ndofs_+r) += localMat(r, j+nsd_+1); // velocity gradient
    }
  } // if (elevec2.M() > 0)

  switch (dyna_)
  {
  case INPAR::ACOU::acou_bdf4:
    solele->eleinteriorGradVelnmmm_ = solele->eleinteriorGradVelnp_;
    solele->eleinteriorVelnmmm_ = solele->eleinteriorVelnp_;
    solele->eleinteriorPressnmmm_ = solele->eleinteriorPressnp_; // no break here
  case INPAR::ACOU::acou_bdf3:
    solele->eleinteriorGradVelnmm_ = solele->eleinteriorGradVelnp_;
    solele->eleinteriorVelnmm_ = solele->eleinteriorVelnp_;
    solele->eleinteriorPressnmm_ = solele->eleinteriorPressnp_; // no break here
  case INPAR::ACOU::acou_bdf2:
    solele->eleinteriorGradVelnm_ = solele->eleinteriorGradVelnp_;
    solele->eleinteriorVelnm_ = solele->eleinteriorVelnp_;
    solele->eleinteriorPressnm_ = solele->eleinteriorPressnp_; // no break here                                          // here you go
  default:
    solele->eleinteriorGradVeln_ = solele->eleinteriorGradVelnp_;
    solele->eleinteriorVeln_ = solele->eleinteriorVelnp_;
    solele->eleinteriorPressn_ = solele->eleinteriorPressnp_;
    break; // here you go
  }

  //if(dyna_==INPAR::ACOU::acou_trapezoidal) // then we have to set the initial trace field, because the calculation of the residual in the first time step needs a correct trace field!
  {
    // trace variable
    int nfdofs = 0;
    for (unsigned int face=0; face<nfaces_; ++face)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Faces()[face]->Degree(),shapes_->usescompletepoly_,2 * ele->Faces()[face]->Degree());
      shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      shapesface_->EvaluateFace(*ele, face);

      Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_,shapesface_->nfdofs_);
      Epetra_SerialDenseMatrix trVec(shapesface_->nfdofs_,nsd_);

      mass.Scale(0.);
      trVec.Scale(0.);

      for (unsigned int q=0; q<shapesface_->nfdofs_; ++q)
      {
        const double fac = shapesface_->jfac(q);
        double xyz[nsd_];
        for (unsigned int d=0; d<nsd_; ++d)
          xyz[d] = shapesface_->xyzreal(d,q);
        double p;
        double gradient[nsd_];
        double velgrad[nsd_*nsd_];

        EvaluateAll(*start_func, xyz, p, gradient, velgrad, 0.0); // u and p at quadrature point

        // now fill the components in the mass matrix and the right hand side
        for (unsigned int i=0; i<shapesface_->nfdofs_; ++i)
        {
          // mass matrix
          for (unsigned int j=0; j<shapesface_->nfdofs_; ++j)
            mass(i,j) += shapesface_->shfunct(i,q) * shapesface_->shfunct(j,q) * fac;
          for(unsigned int d=0; d<nsd_; ++d)
            trVec(i,d) += shapesface_->shfunct(i,q) * gradient[d] * fac;
        }
      }

      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(mass);
      inverseMass.SetVectors(trVec,trVec);
      inverseMass.Solve();

      for (unsigned int i=0; i<shapesface_->nfdofs_; ++i)
        for(unsigned int d=0; d<nsd_; ++d)
          elevec1(nfdofs+i+d*shapesface_->nfdofs_) = trVec(i,d);

      nfdofs += shapesface_->nfdofs_*nsd_;
    }
  }

  return 0;
} // ProjectField

/*----------------------------------------------------------------------*
 * ProjectDirichField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouSolEleCalc<distype>::ProjectDirichField(
    DRT::ELEMENTS::Acou*                 ele,
    Teuchos::ParameterList&              params,
    Epetra_SerialDenseVector&            elevec1)
{
  // trace variable (zero, because no memory, no time derivative)
  dsassert(elevec1.M() == int(nfaces_*shapesface_->nfdofs_*nsd_), "Wrong size in project vector 1");
  elevec1.Scale(0.0);

  // this is it:
  if(params.isParameter("faceconsider"))
  {
    Teuchos::Array<int> *functno = params.getPtr<Teuchos::Array<int> >("funct");
    const unsigned int *faceConsider = params.getPtr<unsigned int>("faceconsider");
    double *time = params.getPtr<double>("time");

    DRT::UTILS::ShapeValuesFaceParams svfparams(ele->Faces()[*faceConsider]->Degree(),shapes_->usescompletepoly_,2 * ele->Faces()[*faceConsider]->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, *faceConsider);

    Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
    Epetra_SerialDenseMatrix trVec(shapesface_->nfdofs_, nsd_+1);

    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      const double fac = shapesface_->jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapesface_->xyzreal(d, q);
      double p;
      double gradient[nsd_];
      double velgrad[nsd_*nsd_];

      EvaluateAll((*functno)[0], xyz,  p, gradient, velgrad, *time);

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
          mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;
        trVec(i,nsd_) += shapesface_->shfunct(i, q) * p * fac;
        for(unsigned int j=0; j<nsd_; ++j)
          trVec(i,j) += shapes_->shfunct(i,q) * gradient[j] * fac;
      }
    }

    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec, trVec);
    inverseMass.Solve();

    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      for(unsigned int j=0; j<nsd_; ++j)
        elevec1(i+j*shapesface_->nfdofs_) = trVec(i,j);
  }

  return 0;
} // ProjectDirichField

/*----------------------------------------------------------------------*
 * EvaluateAll
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::EvaluateAll(const int start_func,
                                                          const double (&xyz)[nsd_],
                                                          double  &p,
                                                          double (&v)[nsd_],
                                                          double (&gv)[nsd_*nsd_],
                                                          double t) const

{
  p = DRT::Problem::Instance()->Funct(start_func-1).Evaluate(0,xyz,0.0,NULL);
  int numcomp = DRT::Problem::Instance()->Funct(start_func-1).NumberComponents();
  if(numcomp != 1 && numcomp != nsd_+1 && numcomp != nsd_*nsd_+nsd_+1) dserror("supply 1, DIM+1, or DIM*DIM+DIM+1 components (pressure, velocity, velocity gradient)");

  if(numcomp==nsd_+1 || numcomp==nsd_*nsd_+nsd_+1)
  {
    for(unsigned int d=0; d<nsd_; ++d)
      v[d] = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(d+1, xyz, t, NULL);
  }
  else
  {
    for(unsigned int d=0; d<nsd_; ++d)
      v[d] = 0.0;
  }
  if(numcomp==nsd_*nsd_+nsd_+1)
  {
    for(unsigned int d=0; d<nsd_*nsd_; ++d)
      gv[d] = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(nsd_+1+d, xyz, t, NULL);
  }
  else
    for(unsigned int d=0; d<nsd_*nsd_; ++d)
      gv[d] = 0.0;

  return;
}

/*----------------------------------------------------------------------*
 * ProjectOpticalField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouSolEleCalc<distype>::ProjectOpticalField(
    DRT::ELEMENTS::Acou*                 ele,
    Teuchos::ParameterList&              params,
    Epetra_SerialDenseVector&            elevec2)
{
  bool meshconform = params.get<bool>("mesh conform");
  DRT::ELEMENTS::AcouSol * solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);

  if(meshconform)
  {
    // get the absorption coefficient
    double absorptioncoeff = params.get<double>("absorption");

    Teuchos::RCP<std::vector<double> > nodevals = params.get<Teuchos::RCP<std::vector<double> > >("nodevals");
    int numlightnode = (*nodevals).size()/(nsd_+1);

    double lightxyz[numlightnode][nsd_];
    double values[numlightnode];

    for(int i=0; i<numlightnode; ++i)
    {
      for(unsigned int j=0; j<nsd_; ++j)
        lightxyz[i][j] = (*nodevals)[i*(nsd_+1)+j];
      values[i] = (*nodevals)[i*(nsd_+1)+nsd_];
    }

    // now we have locations "lightxyz" and corresponding phis in "values"
    // what we want to do now is to evaluate the field in the locations of the
    // Gauss points in this element, to assign an initial distribution
    // so we do exactly the same as in ProjectInitalField BUT call an other
    // evaluate function which will interpolate from given lightxyz and values

    shapes_->Evaluate(*ele);

    // reshape elevec2 as matrix
    dsassert(elevec2.M() == 0 ||
             unsigned(elevec2.M()) == (nsd_*nsd_+nsd_+2)*shapes_->ndofs_, "Wrong size in project vector 2");
    // internal variables
    if (elevec2.M() > 0)
    {
      Epetra_SerialDenseVector localMat(shapes_->ndofs_);
      Epetra_SerialDenseMatrix massPart(shapes_->ndofs_,shapes_->nqpoints_);

      for (unsigned int q=0; q<shapes_->nqpoints_; ++q )
      {
        const double fac = shapes_->jfac(q);
        const double sqrtfac = std::sqrt(fac);
        double xyz[nsd_];
        for (unsigned int d=0; d<nsd_; ++d)
          xyz[d] = shapes_->xyzreal(d,q); // coordinates of quadrature point in real coordinates
        double p = 0.0;
        EvaluateLight(lightxyz,values,numlightnode, xyz, p, absorptioncoeff); // p at quadrature point

        for (unsigned int i=0; i<shapes_->ndofs_; ++i)
        {
          massPart(i,q) = shapes_->shfunct(i,q) * sqrtfac;
          localMat(i) += shapes_->shfunct(i,q) * p * fac;
        }
      }
      Epetra_SerialDenseMatrix massMat(shapes_->ndofs_,shapes_->ndofs_);
      massMat.Multiply('N','T',1.,massPart,massPart,0.);
      {
        Epetra_SerialDenseSolver inverseMass;
        inverseMass.SetMatrix(massMat);
        inverseMass.SetVectors(localMat,localMat);
        inverseMass.Solve();
      }

      for (unsigned int r = 0; r < shapes_->ndofs_; ++r)
      {
        elevec2[r * (nsd_ + 1) + nsd_] += localMat(r, 0); // pressure
        solele->eleinteriorPressnp_(r) += localMat(r, 0); // pressure
      }
    }

  } // if(meshconform)

  switch(dyna_)
  {
  case INPAR::ACOU::acou_bdf4:
    solele->eleinteriorPressnmmm_ = solele->eleinteriorPressnp_; // no break here
  case INPAR::ACOU::acou_bdf3:
    solele->eleinteriorPressnmm_ = solele->eleinteriorPressnp_; // no break here
  case INPAR::ACOU::acou_bdf2:
    solele->eleinteriorPressnm_ = solele->eleinteriorPressnp_; // no break here                                          // here you go
  default:
    solele->eleinteriorPressn_ = solele->eleinteriorPressnp_;
    break;
  }

  // trace variable (zero, because no memory, no time derivative)
  // dsassert(unsigned(elevec1.M()) == nfaces_*shapes_->nfdofs_*nsd_, "Wrong size in project vector 1");
  // elevec1.Scale(0.0);

  return 0;
} // ProjectOpticalField

/*----------------------------------------------------------------------*
 * EvaluateLight
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::EvaluateLight(
                   double lightxyz[][nsd_],
                   double values[],
                   int    numnode,
                   const double (&xyz)[nsd_],
                   double  &p,
                   double absorptioncoeff
) const
{
  // interpolation from nodes
  if (distype == DRT::Element::quad4)
  {
    if(numnode!=4) dserror("wrong number of nodes given");

    //*******************************************************************
    LINALG::Matrix<4,4> coeff(true);
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
      inverseCoeff.Solve();
    }

    p = ( coeff_N(0,0) + coeff_N(1,0) * xyz[0] + coeff_N(2,0) * xyz[1] + coeff_N(3,0) * xyz[0] * xyz[1] ) * values[0]
      + ( coeff_N(0,1) + coeff_N(1,1) * xyz[0] + coeff_N(2,1) * xyz[1] + coeff_N(3,1) * xyz[0] * xyz[1] ) * values[1]
      + ( coeff_N(0,2) + coeff_N(1,2) * xyz[0] + coeff_N(2,2) * xyz[1] + coeff_N(3,2) * xyz[0] * xyz[1] ) * values[2]
      + ( coeff_N(0,3) + coeff_N(1,3) * xyz[0] + coeff_N(2,3) * xyz[1] + coeff_N(3,3) * xyz[0] * xyz[1] ) * values[3];
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
      inverseCoeff.Solve();
    }

    p = ( coeff_N(0,0) + coeff_N(1,0) * xyz[0] + coeff_N(2,0) * xyz[1]  + coeff_N(3,0) * xyz[2] + coeff_N(4,0) * xyz[0] * xyz[1] + coeff_N(5,0) * xyz[0] * xyz[2] + coeff_N(6,0) * xyz[1] * xyz[2] + coeff_N(7,0) * xyz[0] * xyz[1] * xyz[2] ) * values[0]
      + ( coeff_N(0,1) + coeff_N(1,1) * xyz[0] + coeff_N(2,1) * xyz[1]  + coeff_N(3,1) * xyz[2] + coeff_N(4,1) * xyz[0] * xyz[1] + coeff_N(5,1) * xyz[0] * xyz[2] + coeff_N(6,1) * xyz[1] * xyz[2] + coeff_N(7,1) * xyz[0] * xyz[1] * xyz[2] ) * values[1]
      + ( coeff_N(0,2) + coeff_N(1,2) * xyz[0] + coeff_N(2,2) * xyz[1]  + coeff_N(3,2) * xyz[2] + coeff_N(4,2) * xyz[0] * xyz[1] + coeff_N(5,2) * xyz[0] * xyz[2] + coeff_N(6,2) * xyz[1] * xyz[2] + coeff_N(7,2) * xyz[0] * xyz[1] * xyz[2] ) * values[2]
      + ( coeff_N(0,3) + coeff_N(1,3) * xyz[0] + coeff_N(2,3) * xyz[1]  + coeff_N(3,3) * xyz[2] + coeff_N(4,3) * xyz[0] * xyz[1] + coeff_N(5,3) * xyz[0] * xyz[2] + coeff_N(6,3) * xyz[1] * xyz[2] + coeff_N(7,3) * xyz[0] * xyz[1] * xyz[2] ) * values[3]
      + ( coeff_N(0,4) + coeff_N(1,4) * xyz[0] + coeff_N(2,4) * xyz[1]  + coeff_N(3,4) * xyz[2] + coeff_N(4,4) * xyz[0] * xyz[1] + coeff_N(5,4) * xyz[0] * xyz[2] + coeff_N(6,4) * xyz[1] * xyz[2] + coeff_N(7,4) * xyz[0] * xyz[1] * xyz[2] ) * values[4]
      + ( coeff_N(0,5) + coeff_N(1,5) * xyz[0] + coeff_N(2,5) * xyz[1]  + coeff_N(3,5) * xyz[2] + coeff_N(4,5) * xyz[0] * xyz[1] + coeff_N(5,5) * xyz[0] * xyz[2] + coeff_N(6,5) * xyz[1] * xyz[2] + coeff_N(7,5) * xyz[0] * xyz[1] * xyz[2] ) * values[5]
      + ( coeff_N(0,6) + coeff_N(1,6) * xyz[0] + coeff_N(2,6) * xyz[1]  + coeff_N(3,6) * xyz[2] + coeff_N(4,6) * xyz[0] * xyz[1] + coeff_N(5,6) * xyz[0] * xyz[2] + coeff_N(6,6) * xyz[1] * xyz[2] + coeff_N(7,6) * xyz[0] * xyz[1] * xyz[2] ) * values[6]
      + ( coeff_N(0,7) + coeff_N(1,7) * xyz[0] + coeff_N(2,7) * xyz[1]  + coeff_N(3,7) * xyz[2] + coeff_N(4,7) * xyz[0] * xyz[1] + coeff_N(5,7) * xyz[0] * xyz[2] + coeff_N(6,7) * xyz[1] * xyz[2] + coeff_N(7,7) * xyz[0] * xyz[1] * xyz[2] ) * values[7];
    p *= -absorptioncoeff;
  }
  else
    dserror("not yet implemented"); // TODO

  return;
}

/*----------------------------------------------------------------------*
 * ReadGlobalVectors
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::
ReadGlobalVectors(DRT::Element     * ele,
                  DRT::Discretization    & discretization,
                  const std::vector<int> & lm)
{
  DRT::ELEMENTS::AcouSol * solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);

  reshapeMatrixIfNecessary(interiorGradVelnp_,shapes_->ndofs_*nsd_*nsd_,1);
  reshapeMatrixIfNecessary(interiorVelnp_,shapes_->ndofs_*nsd_,1);
  reshapeMatrixIfNecessary(interiorPressnp_,shapes_->ndofs_,1);

  reshapeMatrixIfNecessary(interiorGradVeln_,shapes_->ndofs_*nsd_*nsd_,1);
  reshapeMatrixIfNecessary(interiorVeln_,shapes_->ndofs_*nsd_,1);
  reshapeMatrixIfNecessary(interiorPressn_,shapes_->ndofs_,1);

  reshapeMatrixIfNecessary(interiorGradVelnm_,shapes_->ndofs_*nsd_*nsd_,1);
  reshapeMatrixIfNecessary(interiorVelnm_,shapes_->ndofs_*nsd_,1);
  reshapeMatrixIfNecessary(interiorPressnm_,shapes_->ndofs_,1);

  reshapeMatrixIfNecessary(interiorGradVelnmm_,shapes_->ndofs_*nsd_*nsd_,1);
  reshapeMatrixIfNecessary(interiorVelnmm_,shapes_->ndofs_*nsd_,1);
  reshapeMatrixIfNecessary(interiorPressnmm_,shapes_->ndofs_,1);

  reshapeMatrixIfNecessary(interiorGradVelnmmm_,shapes_->ndofs_*nsd_*nsd_,1);
  reshapeMatrixIfNecessary(interiorVelnmmm_,shapes_->ndofs_*nsd_,1);
  reshapeMatrixIfNecessary(interiorPressnmmm_,shapes_->ndofs_,1);

  interiorGradVelnp_ = solele->eleinteriorGradVelnp_;
  interiorVelnp_ = solele->eleinteriorVelnp_;
  interiorPressnp_ = solele->eleinteriorPressnp_;
  interiorGradVeln_ = solele->eleinteriorGradVeln_;
  interiorVeln_ = solele->eleinteriorVeln_;
  interiorPressn_ = solele->eleinteriorPressn_;
  interiorGradVelnm_ = solele->eleinteriorGradVelnm_;
  interiorVelnm_ = solele->eleinteriorVelnm_;
  interiorPressnm_ = solele->eleinteriorPressnm_;
  interiorGradVelnmm_ = solele->eleinteriorGradVelnmm_;
  interiorVelnmm_ = solele->eleinteriorVelnmm_;
  interiorPressnmm_ = solele->eleinteriorPressnmm_;
  interiorGradVelnmmm_ = solele->eleinteriorGradVelnmmm_;
  interiorVelnmmm_ = solele->eleinteriorVelnmmm_;
  interiorPressnmmm_ = solele->eleinteriorPressnmmm_;

  // read the HDG solution vector (for traces)
  if(discretization.HasState("trace"))
  {
    traceVal_.resize(nfaces_*shapesface_->nfdofs_*nsd_);
    dsassert(lm.size() == traceVal_.size(), "Internal error");
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState("trace");
    DRT::UTILS::ExtractMyValues(*matrix_state,traceVal_,lm);
  }
  if(discretization.HasState("trace_m"))
  {
    traceValm_.resize(nfaces_*shapesface_->nfdofs_*nsd_);
    Teuchos::RCP<const Epetra_Vector> matrix_state_m = discretization.GetState("trace_m");
    DRT::UTILS::ExtractMyValues(*matrix_state_m,traceValm_,lm);
  }

  return;
} // ReadGlobalVectors

/*----------------------------------------------------------------------*
 * FillRestartVectors
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::
FillRestartVectors(DRT::Element           * ele,
                   DRT::Discretization    & discretization)
{
  std::vector<double> interiorValnp(shapes_->ndofs_*(nsd_*nsd_+nsd_+1));

  for(unsigned int i=0; i<shapes_->ndofs_; ++i)
  {
    interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1] = interiorPressnp_(i);
    for(unsigned int d=0; d<nsd_; ++d)
      interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d] = interiorVelnp_(i+d*shapes_->ndofs_);
    for(unsigned int d=0; d<nsd_; ++d)
      for(unsigned int e=0; e<nsd_; ++e)
        interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d] = interiorGradVelnp_(i+(e*nsd_+d)*shapes_->ndofs_);
  }

  std::vector<int> localDofs = discretization.Dof(1,ele);
  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnp");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for (unsigned int i=0; i<localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intveln");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1] = interiorPressn_(i);
      for(unsigned int d=0; d<nsd_; ++d)
        interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d] = interiorVeln_(i+d*shapes_->ndofs_);
      for(unsigned int d=0; d<nsd_; ++d)
        for(unsigned int e=0; e<nsd_; ++e)
          interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d] = interiorGradVeln_(i+(e*nsd_+d)*shapes_->ndofs_);
    }
    for (unsigned int i=0; i<localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  // only do this if desired!
  if(interiorPressnm_.M()>0 && discretization.HasState(1,"intvelnm"))
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnm");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1] = interiorPressnm_(i);
      for(unsigned int d=0; d<nsd_; ++d)
        interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d] = interiorVelnm_(i+d*shapes_->ndofs_);
      for(unsigned int d=0; d<nsd_; ++d)
        for(unsigned int e=0; e<nsd_; ++e)
          interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d] = interiorGradVelnm_(i+(e*nsd_+d)*shapes_->ndofs_);
    }
    for (unsigned int i=0; i<localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  // only do this if desired!
  if(interiorPressnmm_.M()>0 && discretization.HasState(1,"intvelnmm"))
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnmm");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1] = interiorPressnmm_(i);
      for(unsigned int d=0; d<nsd_; ++d)
        interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d] = interiorVelnmm_(i+d*shapes_->ndofs_);
      for(unsigned int d=0; d<nsd_; ++d)
        for(unsigned int e=0; e<nsd_; ++e)
          interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d] = interiorGradVelnmm_(i+(e*nsd_+d)*shapes_->ndofs_);
    }
    for (unsigned int i=0; i<localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  // only do this if desired!
  if(interiorPressnmmm_.M()>0 && discretization.HasState(1,"intvelnmmm"))
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnmmm");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1] = interiorPressnmmm_(i);
      for(unsigned int d=0; d<nsd_; ++d)
        interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d] = interiorVelnmmm_(i+d*shapes_->ndofs_);
      for(unsigned int d=0; d<nsd_; ++d)
        for(unsigned int e=0; e<nsd_; ++e)
          interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d] = interiorGradVelnmmm_(i+(e*nsd_+d)*shapes_->ndofs_);
    }
    for (unsigned int i=0; i<localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }
  return;
} // FillRestartVectors


/*----------------------------------------------------------------------*
 * ElementInitFromRestart
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::
ElementInitFromRestart(DRT::Element           * ele,
                   DRT::Discretization    & discretization)
{
  DRT::ELEMENTS::AcouSol * solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);
  std::vector<double> interiorValnp(shapes_->ndofs_*(nsd_*nsd_+nsd_+2));
  std::vector<int> localDofs1 = discretization.Dof(1,solele);

  Teuchos::RCP<const Epetra_Vector> intvel = discretization.GetState(1,"intvelnp");
  DRT::UTILS::ExtractMyValues(*intvel,interiorValnp,localDofs1);
  for(unsigned int i=0; i<shapes_->ndofs_; ++i)
  {
    solele->eleinteriorPressnp_(i) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1];
    for(unsigned int d=0; d<nsd_; ++d)
      solele->eleinteriorVelnp_(i+d*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d];
    for(unsigned int d=0; d<nsd_; ++d)
      for(unsigned int e=0; e<nsd_; ++e)
        solele->eleinteriorGradVelnp_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d];
  }

  intvel = discretization.GetState(1,"intveln");
  DRT::UTILS::ExtractMyValues(*intvel,interiorValnp,localDofs1);
  for(unsigned int i=0; i<shapes_->ndofs_; ++i)
  {
    solele->eleinteriorPressn_(i) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1];
    for(unsigned int d=0; d<nsd_; ++d)
      solele->eleinteriorVeln_(i+d*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d];
    for(unsigned int d=0; d<nsd_; ++d)
      for(unsigned int e=0; e<nsd_; ++e)
        solele->eleinteriorGradVeln_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d];
  }

  if(discretization.HasState(1,"intvelnm")&&solele->eleinteriorPressnm_.M()>0)
  {
    intvel = discretization.GetState(1,"intvelnm");
    DRT::UTILS::ExtractMyValues(*intvel,interiorValnp,localDofs1);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      solele->eleinteriorPressnm_(i) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1];
      for(unsigned int d=0; d<nsd_; ++d)
        solele->eleinteriorVelnm_(i+d*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d];
      for(unsigned int d=0; d<nsd_; ++d)
        for(unsigned int e=0; e<nsd_; ++e)
          solele->eleinteriorGradVelnm_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d];
    }
  }

  if(discretization.HasState(1,"intvelnmm")&&solele->eleinteriorPressnmm_.M()>0)
  {
    intvel = discretization.GetState(1,"intvelnmm");
    DRT::UTILS::ExtractMyValues(*intvel,interiorValnp,localDofs1);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      solele->eleinteriorPressnmm_(i) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1];
      for(unsigned int d=0; d<nsd_; ++d)
        solele->eleinteriorVelnmm_(i+d*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d];
      for(unsigned int d=0; d<nsd_; ++d)
        for(unsigned int e=0; e<nsd_; ++e)
          solele->eleinteriorGradVelnmm_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d];
    }
  }

  if(discretization.HasState(1,"intvelnmmm")&&solele->eleinteriorPressnmmm_.M()>0)
  {
    intvel = discretization.GetState(1,"intvelnmmm");
    DRT::UTILS::ExtractMyValues(*intvel,interiorValnp,localDofs1);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      solele->eleinteriorPressnmmm_(i) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+nsd_+1];
      for(unsigned int d=0; d<nsd_; ++d)
        solele->eleinteriorVelnmmm_(i+d*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+nsd_*nsd_+d];
      for(unsigned int d=0; d<nsd_; ++d)
        for(unsigned int e=0; e<nsd_; ++e)
          solele->eleinteriorGradVelnmmm_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnp[i*(nsd_*nsd_+nsd_+1)+e*nsd_+d];
    }
  }
  return;
} // ElementInitFromRestart


/*----------------------------------------------------------------------*
 * ComputeError
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::ComputeError(
    DRT::ELEMENTS::Acou*          ele,
    Teuchos::ParameterList&       params,
    Epetra_SerialDenseVector&     elevec)
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
  if(funcno<0) dserror("please provide an analytic solution for the error calculation, set CALCERRORFUNCNO and the corresponding FUNCT");
  funcno--;
  if(DRT::Problem::Instance()->Funct(funcno).NumberComponents()!=int(nsd_*nsd_)+int(nsd_)+1) dserror("please provide numdim*numdim+numdim+1 components in the function for the analytic solution, first for pressure, following for velocity, following for time derivative of velocity gradient");

  double err_p = 0.0, norm_p = 0.0;
  double numerical_p = 0.0;
  double exact_p = 0.0;

  double err_v = 0.0, norm_v = 0.0;
  double numerical_v[nsd_] = {0.0};
  double exact_v[nsd_] = {0.0};

  double err_h = 0.0, norm_h = 0.0;
  double numerical_h[nsd_*nsd_] = {0.0};
  double exact_h[nsd_*nsd_] = {0.0};


  for(int q=0; q<highquad->NumPoints(); ++q)
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
    for (unsigned int idim = 0; idim < nsd_*nsd_; idim++)
    {
      numerical_h[idim] = 0.0;
      exact_h[idim] = 0.0;
    }
    shapes_->polySpace_->Evaluate(xsi,values);

    DRT::UTILS::shape_function_deriv1<distype>(xsi,deriv);
    xjm.MultiplyNT(deriv,shapes_->xyze);
    double highjfac = xjm.Determinant() * highquad->Weight(q);

    // evaluation of numerical values
    for (unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      numerical_p += values(i) * interiorPressnp_(i);
      for (unsigned int idim = 0; idim < nsd_; ++idim)
        numerical_v[idim] += values(i) * interiorVelnp_(i+idim*shapes_->ndofs_);
      for (unsigned int idim = 0; idim < nsd_*nsd_; idim++)
        numerical_h[idim] += values(i) * interiorGradVelnp_(i+idim*shapes_->ndofs_);
    }

    // evaluation of analytical values
    LINALG::Matrix<nen_,1>  myfunct;
    DRT::UTILS::shape_function<distype>(xsi,myfunct);
    LINALG::Matrix<nsd_,1> xyzmat;
    xyzmat.MultiplyNN(shapes_->xyze,myfunct);

    double xyz[nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
      xyz[d]=xyzmat(d,0);

    exact_p = DRT::Problem::Instance()->Funct(funcno).Evaluate(0,xyz,time,NULL);
    for (unsigned int idim = 0; idim < nsd_; ++idim)
      exact_v[idim] = DRT::Problem::Instance()->Funct(funcno).Evaluate(int(idim)+1,xyz,time,NULL);
    for (unsigned int idim = 0; idim < nsd_*nsd_; idim++)
      exact_h[idim] = DRT::Problem::Instance()->Funct(funcno).Evaluate(int(idim)+1+nsd_,xyz,time,NULL);

    // error calculation
    err_p += ( exact_p - numerical_p ) * ( exact_p - numerical_p ) * highjfac;
    norm_p += exact_p * exact_p * highjfac;

    for (unsigned int idim = 0; idim < nsd_; ++idim)
    {
      err_v += (exact_v[idim] - numerical_v[idim]) * (exact_v[idim] - numerical_v[idim]) * highjfac;
      norm_v += exact_v[idim] * exact_v[idim] * highjfac;
    }
    for (unsigned int idim = 0; idim < nsd_*nsd_; ++idim)
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
} // ComputeError

/*----------------------------------------------------------------------*
 * NodeBasedValues
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::NodeBasedValues(
                          const Teuchos::RCP<MAT::Material>    &mat,
                          DRT::ELEMENTS::Acou*                 ele,
                          Epetra_SerialDenseVector&            elevec1)
{
  dsassert(elevec1.M() == (int)nen_*(2*nsd_+2+6)+2, "Vector does not have correct size");
  elevec1.Scale(0.0);

  const MAT::AcousticSolMat* actmat = static_cast<const MAT::AcousticSolMat*>(mat.get());
  double c = actmat->SpeedofSound();
  double visc = actmat->Viscosity();
  double lambda = c*c - visc;

  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  double norm_p = 0.0, area = 0.0;
  for (unsigned int q=0; q<shapes_->nqpoints_; ++q)
  {
    double numerical_p = 0.0;
    for (unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      numerical_p += shapes_->shfunct(i,q) * interiorPressnp_(i);
    }
    norm_p += numerical_p * shapes_->jfac(q);
    area += shapes_->jfac(q);
  }

  elevec1((2*nsd_+2)*nen_) = norm_p/area;

  for (unsigned int i=0; i<nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim=0;idim<nsd_;idim++)
      shapes_->xsi(idim) = locations(idim,i);
    shapes_->polySpace_->Evaluate(shapes_->xsi,values);

    // compute values for velocity and pressure by summing over all basis functions
    double sump = 0;
    for (unsigned int k=0; k<shapes_->ndofs_; ++k)
    {
      sump += values(k) * interiorPressnp_(k);
    }
    elevec1(2*nsd_*nen_+i) = sump;

    for (unsigned int d=0; d<nsd_; ++d)
    {
      double sumv = 0.0;
      for (unsigned int k=0; k<shapes_->ndofs_; ++k)
        sumv += values(k) * interiorVelnp_(d*shapes_->ndofs_+k);
      elevec1(d*nen_+i) = sumv;
    }
    /* visualization of the velocity gradient
    for (unsigned int d=0; d<nsd_*nsd_; ++d)
    {
      double sumvg = 0.0;
      for (unsigned int k=0; k<shapes_->ndofs_; ++k)
        sumvg += values(k) * interiorGradVelnp_(d*shapes_->ndofs_+k);

      if(d>6) break;
      elevec1(nen_*(2*nsd_+2)+2+d*nen_+i) = sumvg;
    }*/
    // visualization of the stresses sigma = mu ( H + H^T ) - lamda / ( lambda + mu ) p I
    for( unsigned int d=0; d<nsd_; ++d)
      for(unsigned int e=0; e<=d; ++e)
      {
        double sumstress = 0.0;
        if(d != e)
        {
          for (unsigned int k=0; k<shapes_->ndofs_; ++k)
            sumstress += values(k) * visc * ( interiorGradVelnp_(k+d*nsd_*shapes_->ndofs_+e*shapes_->ndofs_) + interiorGradVelnp_(k+e*nsd_*shapes_->ndofs_+d*shapes_->ndofs_)  );
          if(d==1&&e==0)
            elevec1(nen_*(2*nsd_+2)+2+nen_*(3)+i) = sumstress;
          else if(d==2&&e==0)
            elevec1(nen_*(2*nsd_+2)+2+nen_*(5)+i) = sumstress;
          else // (if (d==2&&e==1))
            elevec1(nen_*(2*nsd_+2)+2+nen_*(4)+i) = sumstress;
        }
        else // (if(d==e))
        {
          for (unsigned int k=0; k<shapes_->ndofs_; ++k)
            sumstress += values(k) * ( visc * interiorGradVelnp_(k+d*nsd_*shapes_->ndofs_+e*shapes_->ndofs_) + visc * interiorGradVelnp_(k+e*nsd_*shapes_->ndofs_+d*shapes_->ndofs_) - lambda / (lambda+visc) * interiorPressnp_(k) );
          elevec1(nen_*(2*nsd_+2)+2+nen_*(d)+i) = sumstress;
        }
      }
  }

  locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape);
  Epetra_SerialDenseVector fvalues(shapesface_->nfdofs_);
  Epetra_SerialDenseVector temp(elevec1.M());
  for (unsigned int face=0; face<nfaces_; ++face)
  {
    shapesface_->EvaluateFace(*ele, face);
    for (int i=0; i<DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
    {
      // evaluate shape polynomials in node
      for (unsigned int idim=0;idim<nsd_-1;idim++)
        shapesface_->xsi(idim) = locations(idim,i);
      shapesface_->polySpace_->Evaluate(shapesface_->xsi,fvalues); // TODO: fix face orientation here

      // compute values for velocity and pressure by summing over all basis functions
      for (unsigned int d=0; d<nsd_; ++d)
      {
        double sum = 0.0;
        for (unsigned int k=0; k<shapesface_->nfdofs_; ++k)
          sum += fvalues(k) * traceVal_[face*shapesface_->nfdofs_*nsd_+d*shapesface_->nfdofs_+k];
        elevec1((nsd_+d)*nen_+shapesface_->faceNodeOrder[face][i]) += sum;
        temp((nsd_+d)*nen_+shapesface_->faceNodeOrder[face][i])++;
      }
    }
  }

  for(unsigned int i=0; i<nen_*nsd_; ++i)
    elevec1(nsd_*nen_+i) /= temp(nsd_*nen_+i);

  return;
} // NodeBasedValues

/*----------------------------------------------------------------------*
 * NodeBasedPsi
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::NodeBasedPsi(
                          const Teuchos::RCP<MAT::Material>    &mat,
                          DRT::ELEMENTS::Acou*                 ele,
                          Epetra_SerialDenseVector&            elevec1,
                          double                               dt)
{
  dsassert(elevec1.M() == int(nen_), "Vector does not have correct size");
  elevec1.Scale(0.0);
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  // calculate mass matrix
  //const MAT::AcousticSolMat* actmat = static_cast<const MAT::AcousticSolMat*>(mat.get());
  //double c = actmat->SpeedofSound();

  for (unsigned int i=0; i<nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim=0;idim<nsd_;idim++)
      shapes_->xsi(idim) = locations(idim,i);
    shapes_->polySpace_->Evaluate(shapes_->xsi,values);

    // compute values for velocity and pressure by summing over all basis functions
    double sump = 0.0;
    for (unsigned int k=0; k<shapes_->ndofs_; ++k)
    {
      sump += values(k) * interiorPressnp_(k);
    }
    elevec1(i) = sump  / dt;
  }

  return;
} // NodeBasedPsi

/*----------------------------------------------------------------------*
 * Compute internal and face matrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::ComputeMatrices(const Teuchos::RCP<MAT::Material> &mat,
                                                             DRT::ELEMENTS::Acou &             ele,
                                                             double                            dt,
                                                             INPAR::ACOU::DynamicType          dyna,
                                                             bool                              adjoint)
{
  const MAT::AcousticSolMat* actmat = static_cast<const MAT::AcousticSolMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();
  double tau = 1.0;
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

  switch(dyna_)
  {
  case INPAR::ACOU::acou_bdf4: dt *= 12.0/25.0; break;
  case INPAR::ACOU::acou_bdf3: dt *= 6.0/11.0; break;
  case INPAR::ACOU::acou_bdf2: dt *= 2.0/3.0; break;
  default: break;
  }

  localSolver_->ComputeInteriorMatrices();
  for (unsigned int face=0; face<nfaces_; ++face)
  {
    shapesface_->EvaluateFace(ele, face);
    localSolver_->ComputeFaceMatrices(face);
  }

  // scale the matrices
  localSolver_->amat.Scale(     1.0/dt );
  localSolver_->invamat.Scale(      dt );
  //localSolver_->bmat.Scale(          );
  localSolver_->cmat.Scale(       -1.0 );
  localSolver_->dmat.Scale(      -visc );
  localSolver_->emat.Scale(     rho/dt );
  localSolver_->ehatmat.Scale(     tau );
  //localSolver_->fmat.Scale(          );
  localSolver_->gmat.Scale(       -tau );
  localSolver_->hmat.Scale(       -rho );
  localSolver_->imat.Scale( 1.0/c/c/dt );
  localSolver_->jmat.Scale(        rho );
  localSolver_->kmat.Scale(       visc );
  localSolver_->lmat.Scale(       -tau );
  localSolver_->mmat.Scale(       -1.0 );
  localSolver_->nmat.Scale(        tau );

  return;
} // ComputeMatrices

/*----------------------------------------------------------------------*
 * ComputeInteriorMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::
ComputeInteriorMatrices()
{
  // interior matrices are a, b, d, e, ehat, f, h, i, k, l and lhat

  // standard mass matrices are for example i, k, l
  // we calculate i first, and then fill all related matrices

  Epetra_SerialDenseMatrix  massPart(ndofs_,shapes_.nqpoints_);
  Epetra_SerialDenseMatrix  gradPart(ndofs_*nsd_,shapes_.nqpoints_);

  // loop quadrature points
  for (unsigned int q=0; q<shapes_.nqpoints_; ++q)
  {
    const double sqrtfac = std::sqrt(shapes_.jfac(q));
    // loop shape functions
    for (unsigned int i=0; i<ndofs_; ++i)
    {
      const double valf = shapes_.shfunct(i,q) * sqrtfac;
      massPart(i,q) = valf;
      for (unsigned int d=0; d<nsd_; ++d)
      {
        const double vald = shapes_.shderxy(i*nsd_+d,q) * sqrtfac;
        gradPart (d*ndofs_+i,q) = vald;
      }
    }
  }


  // multiply matrices to perform summation over quadrature points
  for (unsigned int i=0; i<ndofs_; ++i)
  {
    for (unsigned int j=0; j<ndofs_; ++j)
    {
      double sum = 0.0;
      double sums[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        sums[d] = 0.0;
      for (unsigned int k=0; k<shapes_.nqpoints_; ++k)
      {
        sum += massPart(i,k) * massPart(j,k);
        for (unsigned int d=0; d<nsd_; ++d)
          sums[d] += gradPart(d*ndofs_+i,k) * massPart(j,k);
      }
      imat(i,j) = sum;
      for (unsigned int d=0; d<nsd_; ++d)
      {
        hmat(i,j+d*ndofs_) = sums[d];
        fmat(j+d*ndofs_,i) = sums[d];
      }
    }
  }

  Epetra_SerialDenseMatrix tempmassgrad(nsd_*ndofs_,ndofs_);
  tempmassgrad.Multiply('N', 'T', 1., gradPart, massPart, 0.);

  bmat.Scale(0.0);
  dmat.Scale(0.0);
  for(unsigned d=0; d<nsd_; ++d)
    for(unsigned i=0; i<ndofs_;++i)
      for(unsigned j=0; j<ndofs_*nsd_; ++j)
      {
        bmat(j+d*nsd_*ndofs_,i+d*ndofs_)=tempmassgrad(j,i);
        dmat(i+d*ndofs_,j+d*nsd_*ndofs_)=tempmassgrad(j,i);
      }

  // e is also a mass matrix, but for every dimension. here, we have to consider the ordering
  // of the velocity components!
  for (unsigned int i=0; i<ndofs_; ++i)
    for (unsigned int j=0; j<ndofs_; ++j)
      for ( unsigned int d=0; d<nsd_; ++d)
      {
        emat(d*ndofs_+i,d*ndofs_+j) = imat(i,j);
        for (unsigned int e=0; e<nsd_; ++e)
          amat(i+(d*nsd_+e)*ndofs_,j+(d*nsd_+e)*ndofs_) = imat(i,j);
      }

  // we will need the inverse of amat. also, maybe, the element should store the
  // inverse of amat, since this is a quite expensive inversion, especially for
  // three spatial dimensions
  invamat = amat;
  {
    Epetra_SerialDenseSolver inverseamat;
    inverseamat.SetMatrix(invamat);
    inverseamat.Invert();
  }

  return;
} // ComputeInteriorMatrices

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::
ComputeFaceMatrices(const int face)
{

  // missing term in ehat
  for (unsigned int p=0; p<ndofs_; ++p)
  {
    for (unsigned int q=0; q<ndofs_; ++q)
    {
      double tempehat = 0.0;
      for (unsigned int i=0; i<shapesface_.nqpoints_; ++i)
        tempehat += shapesface_.jfac(i) * shapesface_.shfunctI(p,i) * shapesface_.shfunctI(q,i);
      for (unsigned int d=0; d<nsd_; ++d)
        ehatmat(d*ndofs_+q,d*ndofs_+p) += tempehat;
    }
  }

  // n
  for (unsigned int p=0; p<shapesface_.nfdofs_; ++p)
  {
    for (unsigned int q=0; q<=p; ++q)
    {
      double tempn = 0.0;
      for (unsigned int i=0; i<shapesface_.nqpoints_; ++i)
        tempn += shapesface_.jfac(i) * shapesface_.shfunct(p,i) * shapesface_.shfunct(q,i);
      for (unsigned int d=0; d<nsd_; ++d)
      {
        nmat(face*shapesface_.nfdofs_*nsd_+p+d*shapesface_.nfdofs_, face*shapesface_.nfdofs_*nsd_+q+d*shapesface_.nfdofs_) = tempn;
        nmat(face*shapesface_.nfdofs_*nsd_+q+d*shapesface_.nfdofs_, face*shapesface_.nfdofs_*nsd_+p+d*shapesface_.nfdofs_) = tempn;
      }
    }
  }

  // c, g, j, m, n, o
  for (unsigned int p=0; p<shapesface_.nfdofs_; ++p)
  {
    for (unsigned int q=0; q<ndofs_; ++q)
    {
      double tempmat = 0.0;
      for (unsigned int i=0; i<shapesface_.nqpoints_; ++i)
      {
        double temp = shapesface_.jfac(i) * shapesface_.shfunct(p,i) * shapesface_.shfunctI(q,i);
        tempmat += temp;
        for(unsigned int j=0; j<nsd_; ++j)
        {
          double temp_d = temp * shapesface_.normals(j,i);
          jmat(q,p+j*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_) += temp_d;
          mmat(p+j*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_,q) += temp_d;

          for(unsigned int e=0; e<nsd_; ++e)
          {
            kmat(p+e*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_,q+(e*nsd_+j)*ndofs_) += temp_d;
            cmat(q+(e*nsd_+j)*ndofs_,p+e*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_) += temp_d;
          }
        }
      }

      for(unsigned int j=0; j<nsd_; ++j)
      {
        gmat(q+j*ndofs_,p+j*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_) = tempmat;
        lmat(p+j*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_,q+j*ndofs_) = tempmat;
      }
    }
  }

  return;
} // ComputeFaceMatrices

/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::
CondenseLocalPart(Epetra_SerialDenseMatrix &eleMat,
                  INPAR::ACOU::DynamicType dyna)
{

  double theta = 1.0;
  if(dyna==INPAR::ACOU::acou_trapezoidal) theta = 0.66;

  Epetra_SerialDenseMatrix dinvamat(nsd_*ndofs_,nsd_*nsd_*ndofs_);
  dinvamat.Multiply('N','N',theta,dmat,invamat,0.0);

  Epetra_SerialDenseMatrix ol(nsd_*ndofs_,nsd_*ndofs_);
  ol = ehatmat;
  ol.Scale(theta);
  ol += emat;
  ol.Multiply('N','N',-theta,dinvamat,bmat,1.0);

  Epetra_SerialDenseMatrix toinv((nsd_+1)*ndofs_,(nsd_+1)*ndofs_);
  for(unsigned int i=0; i<ndofs_; ++i)
  {
    for(unsigned int j=0; j<ndofs_; ++j)
    {
      toinv(nsd_*ndofs_+i,nsd_*ndofs_+j) = imat(i,j);
      for(unsigned int d=0; d<nsd_; ++d)
      {
        toinv(d*ndofs_+i,d*ndofs_+j) = ol(d*ndofs_+i,d*ndofs_+j);
        toinv(nsd_*ndofs_+i,d*ndofs_+j) = theta * hmat(i,d*ndofs_+j);
        toinv(d*ndofs_+i,nsd_*ndofs_+j) = theta * fmat(d*ndofs_+i,j);
      }
    }
  }
  Epetra_SerialDenseMatrix rhsinv((nsd_+1)*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);

  Epetra_SerialDenseMatrix tempmat1(nsd_*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
  tempmat1 = gmat;
  tempmat1.Multiply('N','N',-theta,dinvamat,cmat,theta);

  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
    for(unsigned int j=0; j<nsd_*ndofs_; ++j)
      rhsinv(j,i) = tempmat1(j,i);

  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
    for(unsigned int j=0; j<ndofs_; ++j)
      rhsinv(nsd_*ndofs_+j,i) = theta * jmat(j,i);

  // invert
  {
      Epetra_SerialDenseSolver inverse;
      inverse.SetMatrix(toinv);
      inverse.SetVectors(rhsinv,rhsinv);
      inverse.Solve();
  }

  tempmat1.Shape(ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
    for(unsigned int j=0; j<ndofs_; ++j)
      tempmat1(j,i) = rhsinv(nsd_*ndofs_+j,i);

  eleMat = nmat;
  eleMat.Multiply('N','N',-theta,mmat,tempmat1,theta);

  tempmat1.Shape(nsd_*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
    for(unsigned int j=0; j<nsd_*ndofs_; ++j)
      tempmat1(j,i) = rhsinv(j,i);

  eleMat.Multiply('N','N',-theta,lmat,tempmat1,1.0);

  ol.Shape(nsd_*nsd_*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
  ol = cmat;
  ol.Multiply('N','N',-theta,bmat,tempmat1,theta);
  tempmat1.Shape(nsd_*nsd_*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
  tempmat1.Multiply('N','N',1.0,invamat,ol,0.0);
  eleMat.Multiply('N','N',-theta,kmat,tempmat1,1.0);

  return;
} // CondenseLocalPart

/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::
ComputeResidual(Teuchos::ParameterList&           params,
                Epetra_SerialDenseVector          & elevec,
                Epetra_SerialDenseVector          & interiorGradVeln,
                Epetra_SerialDenseVector          & interiorVeln,
                Epetra_SerialDenseVector          & interiorPressn,
                std::vector<double>               traceVal,
                INPAR::ACOU::DynamicType          dyna)
{
  bool adjoint = params.get<bool>("adjoint");
  double theta = 1.0;
  if(dyna==INPAR::ACOU::acou_trapezoidal) theta = 0.66;

  Epetra_SerialDenseVector traceVal_SDV(nfaces_*shapesface_.nfdofs_*nsd_);
  for(unsigned i=0; i<nfaces_*shapesface_.nfdofs_*nsd_; ++i)
    traceVal_SDV(i) = traceVal[i];

  if(!adjoint)
  {
    Epetra_SerialDenseMatrix dinvamat(nsd_*ndofs_,nsd_*nsd_*ndofs_);
    dinvamat.Multiply('N','N',theta,dmat,invamat,0.0);

    Epetra_SerialDenseMatrix toinv((nsd_+1)*ndofs_,(nsd_+1)*ndofs_);
    {
      Epetra_SerialDenseMatrix ol(nsd_*ndofs_,nsd_*ndofs_);
      ol = ehatmat;
      ol.Scale(theta);
      ol += emat;
      ol.Multiply('N','N',-theta,dinvamat,bmat,1.0);

      for(unsigned int i=0; i<ndofs_; ++i)
      {
        for(unsigned int j=0; j<ndofs_; ++j)
        {
          toinv(nsd_*ndofs_+i,nsd_*ndofs_+j) = imat(i,j);
          for(unsigned int d=0; d<nsd_; ++d)
          {
            toinv(d*ndofs_+i,d*ndofs_+j) = ol(d*ndofs_+i,d*ndofs_+j);
            toinv(nsd_*ndofs_+i,d*ndofs_+j) = theta * hmat(i,d*ndofs_+j);
            toinv(d*ndofs_+i,nsd_*ndofs_+j) = theta * fmat(d*ndofs_+i,j);
          }
        }
      }
    }
    Epetra_SerialDenseVector rhsinv((nsd_+1)*ndofs_);
    Epetra_SerialDenseVector ol(nsd_*ndofs_);

    // source term
    Epetra_SerialDenseVector dummy(nsd_*ndofs_);
    ComputeSource(params,dummy,ol);

    ol.Multiply('N','N',1.0,emat,interiorVeln,1.0);
    ol.Multiply('N','N',-(1.0-theta),dmat,interiorGradVeln,1.0);
    ol.Multiply('N','N',-(1.0-theta),ehatmat,interiorVeln,1.0);
    ol.Multiply('N','N',-(1.0-theta),fmat,interiorPressn,1.0);
    ol.Multiply('N','N',-(1.0-theta),gmat,traceVal_SDV,1.0);
    Epetra_SerialDenseVector tempvec1(nsd_*nsd_*ndofs_);
    tempvec1.Multiply('N','N',1.0,amat,interiorGradVeln,0.0);
    tempvec1.Multiply('N','N',-(1.0-theta),bmat,interiorVeln,1.0);
    tempvec1.Multiply('N','N',-(1.0-theta),cmat,traceVal_SDV,1.0);
    ol.Multiply('N','N',-1.0,dinvamat,tempvec1,1.0);
    for(unsigned int i=0; i<nsd_*ndofs_; ++i)
      rhsinv(i) = ol(i);

    ol.Resize(ndofs_);
    ol.Multiply('N','N',1.0,imat,interiorPressn,0.0);
    ol.Multiply('N','N',-(1.0-theta),hmat,interiorVeln,1.0);
    ol.Multiply('N','N',-(1.0-theta),jmat,traceVal_SDV,1.0);
    for(unsigned int i=0; i<ndofs_; ++i)
      rhsinv(nsd_*ndofs_+i) = ol(i);

    // invert
    {
      Epetra_SerialDenseSolver inverse;
      inverse.SetMatrix(toinv);
      inverse.SetVectors(rhsinv,rhsinv);
      inverse.Solve();
    }

    elevec.Multiply('N','N',-(1.0-theta),kmat,interiorGradVeln,0.0);
    elevec.Multiply('N','N',-(1.0-theta),lmat,interiorVeln,1.0);
    elevec.Multiply('N','N',-(1.0-theta),mmat,interiorPressn,1.0);
    elevec.Multiply('N','N',-(1.0-theta),nmat,traceVal_SDV,1.0);

    for(unsigned int i=0; i<ndofs_; ++i)
      ol(i) = rhsinv(nsd_*ndofs_+i);
    elevec.Multiply('N','N',-theta,mmat,ol,1.0);

    ol.Shape(nsd_*ndofs_,1);
    for(unsigned int i=0; i<nsd_*ndofs_; ++i)
      ol(i,0) = rhsinv(i);
    elevec.Multiply('N','N',-theta,lmat,ol,1.0);

    tempvec1.Multiply('N','N',-theta,bmat,ol,1.0);
    ol.Shape(nsd_*nsd_*ndofs_,1);
    ol.Multiply('N','N',1.0,invamat,tempvec1,0.0);
    elevec.Multiply('N','N',-theta,kmat,ol,1.0);
  }
  else
  {
    dserror("adjoint not yet implemented for solid acoustics");
  } // else ** if(adjoint)

  return;
} // ComputeResidual


/*----------------------------------------------------------------------*
 * ComputeAbsorbingBC
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::
ComputeAbsorbingBC(DRT::ELEMENTS::Acou*        ele,
                   Teuchos::ParameterList&     params,
                   Teuchos::RCP<MAT::Material> & mat,
                   int                         face,
                   Epetra_SerialDenseMatrix    &elemat,
                   Epetra_SerialDenseVector    &elevec)
{
  const MAT::AcousticSolMat* actmat = static_cast<const MAT::AcousticSolMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();
  bool resonly = params.get<bool>("resonly");

  // warning! this is not the correct first order absorbing boundary condition
  // check the hdg_acou paper, section 3.5 and implement it!

  if(!resonly)
  {
    for (unsigned int p=0; p<shapesface_.nfdofs_; ++p)
    {
      for (unsigned int q=0; q<shapesface_.nfdofs_; ++q)
      {
        double temp[nsd_*nsd_];
        for(unsigned int i=0; i<nsd_*nsd_; ++i)
          temp[i] = 0.0;
        for(unsigned int i=0; i<shapesface_.nqpoints_; ++i)
          for(unsigned int d=0; d<nsd_; ++d)
            for(unsigned int e=0; e<nsd_; ++e)
              temp[d+e*nsd_] += shapesface_.jfac(i) * shapesface_.shfunct(p,i) * shapesface_.shfunct(q,i) * shapesface_.normals(d,i) * shapesface_.normals(e,i);

        for (unsigned int d=0; d<nsd_; ++d)
        {
          for(unsigned int e=0; e<nsd_; ++e)
          {
            elemat(face*shapesface_.nfdofs_*nsd_+p+d*shapesface_.nfdofs_, face*shapesface_.nfdofs_*nsd_+q+e*shapesface_.nfdofs_) += rho * c * temp[d+e*nsd_];
          }
        }
      }
    }
  }

  return;
} // ComputeAbsorbingBC

/*----------------------------------------------------------------------*
 * ComputeSource
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::
ComputeSource(Teuchos::ParameterList&           params,
              Epetra_SerialDenseVector          & interiorSourcen,
              Epetra_SerialDenseVector          & interiorSourcenp)
{
  int funcno = params.get<int>("sourcefuncno");
  if(funcno<0) return; // there is no such thing as a volume force

  if(DRT::Problem::Instance()->Funct(funcno).NumberComponents()!=int(nsd_))
    dserror("for acoustic sol elements, the body load is on the velocity and hence has to have numdim components");

  // the vector to be filled
  Epetra_SerialDenseVector source(shapes_.ndofs_*nsd_);

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

    for(unsigned int dim=0; dim<nsd_; ++dim)
    {
      // calculate right hand side contribution for dp/dt
      f_value = DRT::Problem::Instance()->Funct(funcno).Evaluate(int(dim),xyz,tn,NULL);
      f_value_p = DRT::Problem::Instance()->Funct(funcno).Evaluate(int(dim),xyz,tp,NULL);

      // add it all up
      for(unsigned int i=0; i<shapes_.ndofs_; ++i)
      {
        interiorSourcen(i+dim*shapes_.ndofs_) += shapes_.shfunct(i,q) * f_value * shapes_.jfac(q);
        interiorSourcenp(i+dim*shapes_.ndofs_) += shapes_.shfunct(i,q) * f_value_p * shapes_.jfac(q);
      }
    }
  }

  return;
} // ComputeSource


/*----------------------------------------------------------------------*
 * ComputeSourcePressureMonitor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::
ComputeSourcePressureMonitor(DRT::ELEMENTS::Acou*        ele,
                   Teuchos::ParameterList&     params,
                   int                         face,
                   Epetra_SerialDenseVector    &elevec)
{
  // get the values for the source term!
  Teuchos::RCP<Epetra_MultiVector> adjointrhs = params.get<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs");
  int step = params.get<int>("step");
  int stepmax = adjointrhs->NumVectors();
  if(step<stepmax)
  {
    // here, we have to perform the Schur complement thing, because the source term from the adjoint
    // is not applied to the trace field but to the pressure field, which is an internal field
    // also, we have to evaluate the source at the given face....
    Epetra_SerialDenseVector sourceterm(ndofs_);

    bool smoothing = true;
    if(!smoothing)
    {
    }
    else
    {
      const int * fnodeIds = ele->Faces()[face]->NodeIds();
      int numfnode = ele->Faces()[face]->NumNode();
      double fnodexyz[numfnode][nsd_];
      double values[numfnode];

      for(int i=0; i<numfnode; ++i)
      {
        int localnodeid = adjointrhs->Map().LID(fnodeIds[i]);
        if(localnodeid < 0) dserror("could not find entry %d on this proc",fnodeIds[i]);
        values[i] = adjointrhs->operator ()(stepmax-step-1)->operator [](localnodeid); // in inverse order -> we're integrating backwards in time
        for(unsigned int d=0; d<nsd_; ++d)
          fnodexyz[i][d] = shapesface_.xyze(d,i);

      } // for(int i=0; i<numfnode; ++i)

      // so the source term is face based, hence, we calculate the face contribution
      // just as in AcouEleCalc::LocalSolver::ComputeSourcePressureMonitor
      // but then have to deal with the complement thing
      Epetra_SerialDenseMatrix mass(shapesface_.nfdofs_,shapesface_.nfdofs_);
      Epetra_SerialDenseVector trVec(shapesface_.nfdofs_);

      for(unsigned int q=0; q<shapesface_.nqpoints_; ++q)
      {
        const double fac = shapesface_.jfac(q);
        double xyz[nsd_];
        for (unsigned int d=0; d<nsd_; ++d)
          xyz[d] = shapesface_.xyzreal(d,q);
        double val = 0.0;

       EvaluateFaceAdjoint(fnodexyz,values,numfnode,xyz,val);

        for (unsigned int i=0; i<shapesface_.nfdofs_; ++i)
        {
          // mass matrix
          for (unsigned int j=0; j<shapesface_.nfdofs_; ++j)
            mass(i,j) += shapesface_.shfunct(i,q) * shapesface_.shfunct(j,q) * fac;
          trVec(i,0) += shapesface_.shfunct(i,q) * val * fac * double(numfnode)/double(shapesface_.nfdofs_);
        }
      } // for(unsigned int q=0; q<shapesface_.nfdofs_; ++q)

      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(mass);
      inverseMass.SetVectors(trVec,trVec);
      inverseMass.Solve();

      // NOW: sort trVec to sourceterm: therefore use shape functions of internal
      // field evaluated at face
      // TODO think of a way to transform face to parent OR to do it right at the first time!!!

      // those shape functions which are nonzero, are those which values should be nonzero
      int count = 0;
      for(unsigned int q=0; q<ndofs_; ++q)
      {
        if((shapesface_.shfunctI(q,0))>0.00001 ||(shapesface_.shfunctI(q,0))<-0.00001)
        {
          sourceterm(q) = trVec(count);
          count++;
        }
      }

    } // else ** if(!smoothing)

    // now, we have to do the condensation!
    Epetra_SerialDenseMatrix toinv((nsd_*nsd_+nsd_+1)*ndofs_,(nsd_*nsd_+nsd_+1)*ndofs_);
    for(unsigned int i=0; i<ndofs_; ++i)
    {
      for(unsigned int j=0; j<ndofs_; ++j)
      {
        toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+j) = imat(j,i);
        for(unsigned int d=0; d<nsd_; ++d)
        {
          toinv((nsd_*nsd_)*ndofs_+d*ndofs_+i,(nsd_*nsd_)*ndofs_+d*ndofs_+j) = emat(d*ndofs_+j,d*ndofs_+i) + ehatmat(d*ndofs_+j,d*ndofs_+i);
          toinv((nsd_*nsd_)*ndofs_+d*ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+j) = hmat(j,d*ndofs_+i);
          toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+i,(nsd_*nsd_)*ndofs_+d*ndofs_+j) = fmat(d*ndofs_+j,i);
        }
        for(unsigned int e=0; e<nsd_*nsd_; ++e)
        {
          for(unsigned int d=0; d<nsd_*nsd_; ++d)
            toinv(e*ndofs_+i,d*ndofs_+j) = amat(e*ndofs_+j,d*ndofs_+i);
          for(unsigned int d=0; d<nsd_; ++d)
          {
            toinv(e*ndofs_+i,(nsd_*nsd_)*ndofs_+d*ndofs_+j) = dmat(d*ndofs_+j,e*ndofs_+i);
            toinv((nsd_*nsd_)*ndofs_+d*ndofs_+i,e*ndofs_+j) = bmat(e*ndofs_+j,d*ndofs_+i);
          }
        }
      }
    }
    Epetra_SerialDenseVector rhsinv((nsd_*nsd_+nsd_+1)*ndofs_);
    for(unsigned int i=0; i<ndofs_; ++i)
      rhsinv(nsd_*nsd_*ndofs_+nsd_*ndofs_+i) = sourceterm(i,0);

    // invert
    {
      Epetra_SerialDenseSolver inverse;
      inverse.SetMatrix(toinv);
      inverse.SetVectors(rhsinv,rhsinv);
      inverse.Solve();
    }

    Epetra_SerialDenseMatrix ol(ndofs_,1);
    for(unsigned int i=0; i<ndofs_; ++i)
      ol(i,0) = rhsinv((nsd_*nsd_+nsd_)*ndofs_+i);
    elevec.Multiply('T','N',-1.0,jmat,ol,0.0);

    ol.Shape(nsd_*ndofs_,1);
    for(unsigned int i=0; i<nsd_*ndofs_; ++i)
      ol(i,0) = rhsinv(nsd_*nsd_*ndofs_+i);
    elevec.Multiply('T','N',-1.0,gmat,ol,1.0);

    ol.Shape(nsd_*nsd_*ndofs_,1);
    for(unsigned int i=0; i<nsd_*nsd_*ndofs_; ++i)
      ol(i,0) = rhsinv(i);
    elevec.Multiply('T','N',-1.0,cmat,ol,1.0);

  } // if(step<stepmax)

  return;
} // ComputeSourcePressureMonitor

/*----------------------------------------------------------------------*
 * EvaluateFaceAdjoint
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::EvaluateFaceAdjoint(
                    double fnodexyz[][nsd_],
                    double values[],
                    int numfnode,
                    const double (&xyz)[nsd_],
                    double &val) const
{
  const DRT::Element::DiscretizationType facedis = DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;
  if (facedis == DRT::Element::line2)
  {
    if(numfnode!=2) dserror("number of nodes per face should be 2 for face discretization = line2");
    val = values[1] * sqrt( (fnodexyz[0][0]-xyz[0])*(fnodexyz[0][0]-xyz[0]) + (fnodexyz[0][1]-xyz[1])*(fnodexyz[0][1]-xyz[1]) )
        + values[0] * sqrt( (fnodexyz[1][0]-xyz[0])*(fnodexyz[1][0]-xyz[0]) + (fnodexyz[1][1]-xyz[1])*(fnodexyz[1][1]-xyz[1]) );
    double dist = sqrt( (fnodexyz[0][0]-fnodexyz[1][0])*(fnodexyz[0][0]-fnodexyz[1][0]) + (fnodexyz[0][1]-fnodexyz[1][1])*(fnodexyz[0][1]-fnodexyz[1][1]) );
    val /= dist;
  }
  else if(facedis == DRT::Element::quad4)
  {
    if(numfnode!=4) dserror("number of nodes per face should be 4 for face discretization = quad4");
    dserror("this is not good, since the inversion does not work");
    LINALG::Matrix<4,4> mat(true);
    for(int i=0; i<4; ++i)
    {
      mat(i,0) = 1.0;
      mat(i,1) = fnodexyz[i][0];
      mat(i,2) = fnodexyz[i][1];
      mat(i,3) = fnodexyz[i][2];
    }
    LINALG::FixedSizeSerialDenseSolver<4,4> inversemat;
    inversemat.SetMatrix(mat);
    inversemat.Invert();

    val = 0.0;
    for(int i=0; i<4; ++i)
    {
      LINALG::Matrix<4,1> lvec(true);
      LINALG::Matrix<4,1> rvec(true);
      rvec.PutScalar(0.0);
      rvec(i) = 1.0;
      lvec.Multiply(mat,rvec);
      val += ( lvec(0) + lvec(1) * xyz[0] + lvec(2) * xyz[1] + lvec(3) * xyz[2] ) * values[i];
    }
  }
  else
    dserror("not yet implemented");

  return;
} // EvaluateFaceAdjoint

/*----------------------------------------------------------------------*
 * ComputeSourcePressureMonitor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::LocalSolver::
ComputeSourcePressureMonitorLine3D(DRT::ELEMENTS::Acou*        ele,
                   Teuchos::ParameterList&     params,
                   Teuchos::RCP<MAT::Material> & mat,
                   int                         face,
                   Epetra_SerialDenseMatrix    &elemat,
                   Epetra_SerialDenseVector    &elevec)
{
  dserror("here is something TODO for you!"); // TODO
} // ComputeSourcePressureMonitorLine3D

/*----------------------------------------------------------------------*
 * ComputePMonNodeVals
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::
ComputePMonNodeVals(DRT::ELEMENTS::Acou*        ele,
                   int                         face,
                   Epetra_SerialDenseVector    &elevec)
{
  if(elevec.M()!=DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace)  dserror("Vector does not have correct size");

  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  Epetra_SerialDenseVector values(shapes_->ndofs_);
  LINALG::Matrix<nsd_,1> xsiFl;

  const int* nodeidsface = ele->Faces()[face]->NodeIds();
  int numfacenode = ele->Faces()[face]->NumNode();
  for (unsigned int i=0; i<nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim=0;idim<nsd_;idim++)
      xsiFl(idim) = locations(idim,i);
    shapes_->polySpace_->Evaluate(xsiFl,values);

    // is node part of this face element?
    int nodeid = ele->NodeIds()[i];
    bool partof = false;
    int facenode = -1;
    for(int j=0; j<numfacenode; ++j)
      if(nodeid == nodeidsface[j])
      {
        partof = true;
        facenode = j;
        break;
      }
    if(!partof) continue;

    // compute values for velocity and pressure by summing over all basis functions
    double sum = 0.0;
    for (unsigned int k=0; k<shapes_->ndofs_; ++k)
      sum += values(k) * interiorPressnp_(k);

    elevec(facenode) = sum;
  }

  return;
} // ComputePMonNodeVals

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouSolEleCalc<distype>::
UpdateInteriorVariablesAndComputeResidual(Teuchos::ParameterList&    params,
                                          DRT::ELEMENTS::Acou *                ele,
                                          Epetra_SerialDenseVector          & elevec,
                                          double dt,
                                          bool errormaps,
                                          bool updateonly)
{
  DRT::ELEMENTS::AcouSol * solele = dynamic_cast<DRT::ELEMENTS::AcouSol*>(ele);
  bool adjoint = params.get<bool>("adjoint");
  double theta = 1.0;
  if(dyna_==INPAR::ACOU::acou_trapezoidal) theta = 0.66;
  int stage = -1;

  Epetra_SerialDenseVector tempGradVelnp;
  Epetra_SerialDenseVector tempVelnp;
  Epetra_SerialDenseVector tempPressnp;
  if (dyna_ == INPAR::ACOU::acou_bdf2)
  {
    tempGradVelnp.Shape(shapes_->ndofs_ * nsd_ * nsd_, 1);
    tempGradVelnp = interiorGradVelnp_;
    tempVelnp.Shape(shapes_->ndofs_ * nsd_, 1);
    tempVelnp = interiorVelnp_;
    tempPressnp.Shape(shapes_->ndofs_, 1);
    tempPressnp = interiorPressnp_;
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_ * nsd_; ++i)
      interiorGradVeln_[i] = interiorGradVelnp_[i] * 4.0 / 3.0 - interiorGradVelnm_[i] / 3.0;
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVeln_[i] = interiorVelnp_[i] * 4.0 / 3.0 - interiorVelnm_[i] / 3.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressn_[i] = interiorPressnp_[i] * 4.0 / 3.0 - interiorPressnm_[i] / 3.0;
    dt *= 2.0 / 3.0;
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf3)
  {
    tempGradVelnp.Shape(shapes_->ndofs_ * nsd_ * nsd_, 1);
    tempGradVelnp = interiorGradVelnp_;
    tempVelnp.Shape(shapes_->ndofs_ * nsd_, 1);
    tempVelnp = interiorVelnp_;
    tempPressnp.Shape(shapes_->ndofs_, 1);
    tempPressnp = interiorPressnp_;
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_* nsd_; ++i)
      interiorGradVeln_[i] = interiorVelnp_[i] * 18.0 / 11.0 - interiorGradVelnm_[i] * 9.0 / 11.0 + interiorGradVelnmm_[i] * 2.0 / 11.0;
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVeln_[i] = interiorVelnp_[i] * 18.0 / 11.0 - interiorVelnm_[i] * 9.0 / 11.0 + interiorVelnmm_[i] * 2.0 / 11.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressn_[i] = interiorPressnp_[i] * 18.0 / 11.0 - interiorPressnm_[i] * 9.0 / 11.0 + interiorPressnmm_[i] * 2.0 / 11.0;
    dt *= 6.0 / 11.0;
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf4)
  {
    tempGradVelnp.Shape(shapes_->ndofs_ * nsd_ * nsd_, 1);
    tempGradVelnp = interiorGradVelnp_;
    tempVelnp.Shape(shapes_->ndofs_ * nsd_, 1);
    tempVelnp = interiorVelnp_;
    tempPressnp.Shape(shapes_->ndofs_, 1);
    tempPressnp = interiorPressnp_;
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_ * nsd_; ++i)
      interiorGradVeln_[i] = interiorGradVelnp_[i] * 48.0 / 25.0 - interiorGradVelnm_[i] * 36.0 / 25.0 + interiorGradVelnmm_[i] * 16.0 / 25.0 - interiorGradVelnmmm_[i] * 3.0 / 25.0;
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVeln_[i] = interiorVelnp_[i] * 48.0 / 25.0 - interiorVelnm_[i] * 36.0 / 25.0 + interiorVelnmm_[i] * 16.0 / 25.0 - interiorVelnmmm_[i] * 3.0 / 25.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressn_[i] = interiorPressnp_[i] * 48.0 / 25.0 - interiorPressnm_[i] * 36.0 / 25.0 + interiorPressnmm_[i] * 16.0 / 25.0 - interiorPressnmmm_[i] * 3.0 / 25.0;
    dt *= 12.0 / 25.0;
  }
  else if (dyna_ == INPAR::ACOU::acou_dirk23 ||
           dyna_ == INPAR::ACOU::acou_dirk33 ||
           dyna_ == INPAR::ACOU::acou_dirk34 ||
           dyna_ == INPAR::ACOU::acou_dirk54)
  {
    // do the dirk
    stage = params.get<int>("stage");

    interiorGradVeln_ = solele->elesgv_[stage];
    interiorVeln_ = solele->elesv_[stage];
    interiorPressn_ = solele->elesp_[stage];
    interiorGradVeln_.Scale(dt);
    interiorVeln_.Scale(dt);
    interiorPressn_.Scale(dt);
  }
  else
  {
    interiorGradVeln_ = interiorGradVelnp_;
    interiorVeln_ = interiorVelnp_;
    interiorPressn_ = interiorPressnp_;

  }

  Epetra_SerialDenseVector traceVal_SDV(nfaces_*shapesface_->nfdofs_*nsd_);
  for(unsigned i=0; i<nfaces_*shapesface_->nfdofs_*nsd_; ++i)
    traceVal_SDV(i) = traceVal_[i];

  Epetra_SerialDenseVector traceVal_SDV_m(nfaces_*shapesface_->nfdofs_*nsd_);
  if(dyna_ == INPAR::ACOU::acou_trapezoidal)
  {
    for(unsigned i=0; i<nfaces_*shapesface_->nfdofs_*nsd_; ++i)
      traceVal_SDV_m(i) = traceValm_[i];
  }

  // *****************************************************
  // update interior variables first
  // *****************************************************

  Epetra_SerialDenseMatrix toinv((nsd_+1)*shapes_->ndofs_,(nsd_+1)*shapes_->ndofs_);
  Epetra_SerialDenseVector tempsource(nsd_*shapes_->ndofs_);
  if(!adjoint)
  {
    Epetra_SerialDenseMatrix dinvamat(nsd_*shapes_->ndofs_,nsd_*nsd_*shapes_->ndofs_);
    dinvamat.Multiply('N','N',theta,localSolver_->dmat,localSolver_->invamat,0.0);

    {
      Epetra_SerialDenseMatrix ol(nsd_*shapes_->ndofs_,nsd_*shapes_->ndofs_);
      ol = localSolver_->ehatmat;
      ol.Scale(theta);
      ol += localSolver_->emat;
      ol.Multiply('N','N',-theta,dinvamat,localSolver_->bmat,1.0);

      for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      {
        for(unsigned int j=0; j<shapes_->ndofs_; ++j)
        {
          toinv(nsd_*shapes_->ndofs_+i,nsd_*shapes_->ndofs_+j) = localSolver_->imat(i,j);
          for(unsigned int d=0; d<nsd_; ++d)
          {
            toinv(d*shapes_->ndofs_+i,d*shapes_->ndofs_+j) = ol(d*shapes_->ndofs_+i,d*shapes_->ndofs_+j);
            toinv(nsd_*shapes_->ndofs_+i,d*shapes_->ndofs_+j) = theta * localSolver_->hmat(i,d*shapes_->ndofs_+j);
            toinv(d*shapes_->ndofs_+i,nsd_*shapes_->ndofs_+j) = theta * localSolver_->fmat(d*shapes_->ndofs_+i,j);
          }
        }
      }
    }

    // calculate rhs parts
    Epetra_SerialDenseVector ol(nsd_*shapes_->ndofs_);
    localSolver_->ComputeSource(params,ol,tempsource);
    ol.Multiply('N','N',1.0,localSolver_->emat,interiorVeln_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->dmat,interiorGradVeln_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->ehatmat,interiorVeln_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->fmat,interiorPressn_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->gmat,traceVal_SDV_m,1.0);
    ol.Multiply('N','N',-theta,localSolver_->gmat,traceVal_SDV,1.0);


    Epetra_SerialDenseVector tempvec1(nsd_*nsd_*shapes_->ndofs_);
    tempvec1.Multiply('N','N',1.0,localSolver_->amat,interiorGradVeln_,0.0);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->bmat,interiorVeln_,1.0);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->cmat,traceVal_SDV_m,1.0);
    tempvec1.Multiply('N','N',-theta,localSolver_->cmat,traceVal_SDV,1.0);
    ol.Multiply('N','N',-1.0,dinvamat,tempvec1,1.0);
    Epetra_SerialDenseVector rhsinv((nsd_+1)*shapes_->ndofs_);
    for(unsigned int i=0; i<nsd_*shapes_->ndofs_; ++i)
      rhsinv(i) = ol(i);

    ol.Shape(shapes_->ndofs_,1);
    ol.Multiply('N','N',1.0,localSolver_->imat,interiorPressn_,0.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->hmat,interiorVeln_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->jmat,traceVal_SDV_m,1.0);
    ol.Multiply('N','N',-theta,localSolver_->jmat,traceVal_SDV,1.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv(nsd_*shapes_->ndofs_+i) = ol(i);

    // invert
    {
      Epetra_SerialDenseSolver inverse;
      inverse.SetMatrix(toinv);
      inverse.Invert();
    }
    Epetra_SerialDenseVector sol((nsd_+1)*shapes_->ndofs_);
    sol.Multiply('N','N',1.0,toinv,rhsinv,0.0);

    for(unsigned int i=0; i<nsd_*shapes_->ndofs_; ++i)
      interiorVelnp_(i) = sol(i);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      interiorPressnp_(i) = sol(nsd_*shapes_->ndofs_+i);

    tempvec1.Multiply('N','N',-theta,localSolver_->bmat,interiorVelnp_,1.0);
    interiorGradVelnp_.Multiply('N','N',1.0,localSolver_->invamat,tempvec1,0.0);
  } // if(!adjoint)
  else
  {
    dserror("adjoint not implemented for solid acoustics");
  } // else ** if(!adjoint)

  // tell this change in the interior variables the discretization
  if(dyna_ == INPAR::ACOU::acou_dirk23 ||
     dyna_ == INPAR::ACOU::acou_dirk33 ||
     dyna_ == INPAR::ACOU::acou_dirk34 ||
     dyna_ == INPAR::ACOU::acou_dirk54 )
  {
    double dirk_a[6][6];
    double dirk_b[6];
    double dirk_c[6];
    int dirk_q;
    ACOU::FillDIRKValues(dyna_, dirk_a, dirk_b, dirk_c, dirk_q);

    solele->eleyp_[stage] = interiorPressnp_;
    solele->eleyv_[stage] = interiorVelnp_;
    solele->eleygv_[stage] = interiorGradVelnp_;

    solele->elefp_[stage] = solele->eleinteriorPressn_;
    solele->elefp_[stage].Scale(-1.0);
    solele->elefp_[stage] += solele->eleyp_[stage];
    solele->elefp_[stage].Scale(1.0 / dt); // dt includes a_ii
    solele->elefv_[stage] = solele->eleinteriorVeln_;
    solele->elefv_[stage].Scale(-1.0);
    solele->elefv_[stage] += solele->eleyv_[stage];
    solele->elefv_[stage].Scale(1.0 / dt);
    solele->elefgv_[stage] = solele->eleinteriorGradVeln_;
    solele->elefgv_[stage].Scale(-1.0);
    solele->elefgv_[stage] += solele->eleygv_[stage];
    solele->elefgv_[stage].Scale(1.0 / dt);

    Epetra_SerialDenseVector tempVecp1(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecv1(shapes_->ndofs_ * nsd_);
    Epetra_SerialDenseVector tempVecgv1(shapes_->ndofs_ * nsd_ * nsd_);
    for (int i = 0; i < stage; ++i)
    {
      tempVecp1 = solele->elefp_[i];
      tempVecp1.Scale(-dirk_a[stage][i] / dirk_a[stage][stage]);
      solele->elefp_[stage] += tempVecp1;
      tempVecv1 = solele->elefv_[i];
      tempVecv1.Scale(-dirk_a[stage][i] / dirk_a[stage][stage]);
      solele->elefv_[stage] += tempVecv1;
      tempVecgv1 = solele->elefgv_[i];
      tempVecgv1.Scale(-dirk_a[stage][i] / dirk_a[stage][stage]);
      solele->elefgv_[stage] += tempVecgv1;
    }
    tempVecp1 = solele->elefp_[stage];
    tempVecp1.Scale(dt * dirk_b[stage] / dirk_a[stage][stage]);
    solele->eleinteriorPressnp_ += tempVecp1;
    tempVecv1 = solele->elefv_[stage];
    tempVecv1.Scale(dt * dirk_b[stage] / dirk_a[stage][stage]);
    solele->eleinteriorVelnp_ += tempVecv1;
    tempVecgv1 = solele->elefgv_[stage];
    tempVecgv1.Scale(dt * dirk_b[stage] / dirk_a[stage][stage]);
    solele->eleinteriorGradVelnp_ += tempVecgv1;

    if (stage == dirk_q - 1)
    {
      solele->eleinteriorPressn_ = solele->eleinteriorPressnp_;
      solele->eleinteriorVeln_ = solele->eleinteriorVelnp_;
      solele->eleinteriorGradVeln_ = solele->eleinteriorGradVelnp_;
    }
  } // if(dyna_ == INPAR::ACOU::acou_dirk?? )
  else
  {
    solele->eleinteriorPressnp_ = interiorPressnp_;
    solele->eleinteriorVelnp_ = interiorVelnp_;
    solele->eleinteriorGradVelnp_ = interiorGradVelnp_;
  }

  // *****************************************************
  // local postprocessing to calculate error maps
  // *****************************************************

  // TODO: local postprocessing

  if (updateonly) return;

  // *****************************************************
  // compute residual second (reuse intermediate matrices)
  // *****************************************************

  if (dyna_ == INPAR::ACOU::acou_dirk23 ||
      dyna_ == INPAR::ACOU::acou_dirk33 ||
      dyna_ == INPAR::ACOU::acou_dirk34 ||
      dyna_ == INPAR::ACOU::acou_dirk54)
  {
    double dirk_a[6][6];
    double dirk_b[6];
    double dirk_c[6];
    int dirk_q;
    ACOU::FillDIRKValues(dyna_, dirk_a, dirk_b, dirk_c, dirk_q);

    stage++;
    if (stage == dirk_q)
    {
      stage = 0;
      solele->eleinteriorPressnp_ = solele->eleinteriorPressn_;
      solele->eleinteriorVelnp_ = solele->eleinteriorVeln_;
      solele->eleinteriorGradVelnp_ = solele->eleinteriorGradVeln_;
    }
    Epetra_SerialDenseVector tempVecp1(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecp2(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecv1(shapes_->ndofs_ * nsd_);
    Epetra_SerialDenseVector tempVecv2(shapes_->ndofs_ * nsd_);
    Epetra_SerialDenseVector tempVecgv1(shapes_->ndofs_ * nsd_ * nsd_);
    Epetra_SerialDenseVector tempVecgv2(shapes_->ndofs_ * nsd_ * nsd_);

    solele->elesp_[stage] = solele->eleinteriorPressn_;
    solele->elesp_[stage].Scale(1.0 / dt); // dt includes a_ii
    solele->elesv_[stage] = solele->eleinteriorVeln_;
    solele->elesv_[stage].Scale(1.0 / dt); // dt includes a_ii
    solele->elesgv_[stage] = solele->eleinteriorGradVeln_;
    solele->elesgv_[stage].Scale(1.0 / dt); // dt includes a_ii

    for (int i = 0; i < stage; ++i)
    {
      tempVecp1 = solele->elesp_[i];
      tempVecp1.Scale(-1.0);
      tempVecp2 = solele->eleyp_[i];
      tempVecp2.Scale(1.0 / dt);
      tempVecp1 += tempVecp2;
      tempVecp1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      solele->elesp_[stage] += tempVecp1;
      tempVecv1 = solele->elesv_[i];
      tempVecv1.Scale(-1.0);
      tempVecv2 = solele->eleyv_[i];
      tempVecv2.Scale(1.0 / dt);
      tempVecv1 += tempVecv2;
      tempVecv1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      solele->elesv_[stage] += tempVecv1;
      tempVecgv1 = solele->elesgv_[i];
      tempVecgv1.Scale(-1.0);
      tempVecgv2 = solele->eleygv_[i];
      tempVecgv2.Scale(1.0 / dt);
      tempVecgv1 += tempVecgv2;
      tempVecgv1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      solele->elesgv_[stage] += tempVecgv1;
    }

    // these vectors s are used for the calculation of the residual -> write them to used local solver variable
    interiorPressnp_ = solele->elesp_[stage];
    interiorVelnp_ = solele->elesv_[stage];
    interiorGradVelnp_ = solele->elesgv_[stage];
    interiorPressnp_.Scale(dt);
    interiorVelnp_.Scale(dt);
    interiorGradVelnp_.Scale(dt);
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf2)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_ * nsd_; ++i)
      interiorGradVelnp_[i] = 4.0 / 3.0 * interiorGradVelnp_[i] - 1.0 / 3.0 * tempGradVelnp[i];
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 4.0 / 3.0 * interiorVelnp_[i] - 1.0 / 3.0 * tempVelnp[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 4.0 / 3.0 * interiorPressnp_[i] - 1.0 / 3.0 * tempPressnp[i];
    solele->eleinteriorPressnm_ = tempPressnp;
    solele->eleinteriorVelnm_ = tempVelnp;
    solele->eleinteriorGradVelnm_ = tempGradVelnp;
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf3)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_ * nsd_; ++i)
      interiorGradVelnp_[i] = 18.0 / 11.0 * interiorGradVelnp_[i] - 9.0 / 11.0 * tempGradVelnp[i] + 2.0 / 11.0 * interiorGradVelnm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 18.0 / 11.0 * interiorVelnp_[i] - 9.0 / 11.0 * tempVelnp[i] + 2.0 / 11.0 * interiorVelnm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 18.0 / 11.0 * interiorPressnp_[i] - 9.0 / 11.0 * tempPressnp[i] + 2.0 / 11.0 * interiorPressnm_[i];
    solele->eleinteriorPressnmm_ = solele->eleinteriorPressnm_;
    solele->eleinteriorPressnm_ = tempPressnp;
    solele->eleinteriorVelnmm_ = solele->eleinteriorVelnm_;
    solele->eleinteriorVelnm_ = tempVelnp;
    solele->eleinteriorGradVelnmm_ = solele->eleinteriorGradVelnm_;
    solele->eleinteriorGradVelnm_ = tempGradVelnp;
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf4)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_ * nsd_; ++i)
      interiorGradVelnp_[i] = 48.0 / 25.0 * interiorGradVelnp_[i] - 36.0 / 25.0 * tempGradVelnp[i] + 16.0 / 25.0 * interiorGradVelnm_[i] - 3.0 / 25.0 * interiorGradVelnmm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 48.0 / 25.0 * interiorVelnp_[i] - 36.0 / 25.0 * tempVelnp[i] + 16.0 / 25.0 * interiorVelnm_[i] - 3.0 / 25.0 * interiorVelnmm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 48.0 / 25.0 * interiorPressnp_[i] - 36.0 / 25.0 * tempPressnp[i] + 16.0 / 25.0 * interiorPressnm_[i] - 3.0 / 25.0 * interiorPressnmm_[i];
    solele->eleinteriorPressnmmm_ = solele->eleinteriorPressnmm_;
    solele->eleinteriorPressnmm_ = solele->eleinteriorPressnm_;
    solele->eleinteriorPressnm_ = tempPressnp;
    solele->eleinteriorVelnmmm_ = solele->eleinteriorVelnmm_;
    solele->eleinteriorVelnmm_ = solele->eleinteriorVelnm_;
    solele->eleinteriorVelnm_ = tempVelnp;
    solele->eleinteriorGradVelnmmm_ = solele->eleinteriorGradVelnmm_;
    solele->eleinteriorGradVelnmm_ = solele->eleinteriorGradVelnm_;
    solele->eleinteriorGradVelnm_ = tempGradVelnp;
  }

  // calculate rhs parts
  if(!adjoint)
  {
    Epetra_SerialDenseMatrix dinvamat(nsd_*shapes_->ndofs_,nsd_*nsd_*shapes_->ndofs_);
    dinvamat.Multiply('N','N',theta,localSolver_->dmat,localSolver_->invamat,0.0);

    Epetra_SerialDenseVector ol(nsd_*shapes_->ndofs_);
    ol+=tempsource;
    ol.Multiply('N','N',1.0,localSolver_->emat,interiorVelnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->dmat,interiorGradVelnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->ehatmat,interiorVelnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->fmat,interiorPressnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->gmat,traceVal_SDV,1.0);
    Epetra_SerialDenseVector tempvec1(nsd_*nsd_*shapes_->ndofs_);
    tempvec1.Multiply('N','N',1.0,localSolver_->amat,interiorGradVelnp_,0.0);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->bmat,interiorVelnp_,1.0);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->cmat,traceVal_SDV,1.0);
    ol.Multiply('N','N',-1.0,dinvamat,tempvec1,1.0);
    Epetra_SerialDenseVector rhsinv((nsd_+1)*shapes_->ndofs_);
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      rhsinv(i) = ol(i);

    ol.Resize(shapes_->ndofs_);
    ol.Multiply('N','N',1.0,localSolver_->imat,interiorPressnp_,0.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->hmat,interiorVelnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->jmat,traceVal_SDV,1.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv(shapes_->ndofs_*nsd_+i) = ol(i);

    Epetra_SerialDenseVector sol((nsd_+1)*shapes_->ndofs_);
    sol.Multiply('N','N',1.0,toinv,rhsinv,0.0);

    elevec.Multiply('N','N',-(1.0-theta),localSolver_->kmat,interiorGradVelnp_,0.0);
    elevec.Multiply('N','N',-(1.0-theta),localSolver_->lmat,interiorVelnp_,1.0);
    elevec.Multiply('N','N',-(1.0-theta),localSolver_->mmat,interiorPressnp_,1.0);
    elevec.Multiply('N','N',-(1.0-theta),localSolver_->nmat,traceVal_SDV,1.0);

    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      ol(i) = sol(nsd_*shapes_->ndofs_+i);
    elevec.Multiply('N','N',-theta,localSolver_->mmat,ol,1.0);

    ol.Shape(shapes_->ndofs_*nsd_,1);
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      ol(i) = sol(i);
    elevec.Multiply('N','N',-theta,localSolver_->lmat,ol,1.0);

    tempvec1.Multiply('N','N',-theta,localSolver_->bmat,ol,1.0);
    ol.Shape(nsd_*nsd_*shapes_->ndofs_,1);
    ol.Multiply('N','N',1.0,localSolver_->invamat,tempvec1,0.0);
    elevec.Multiply('N','N',-theta,localSolver_->kmat,ol,1.0);
  } // if(!adjoint)
  else
  {
    dserror("adjoint not yet implemented for solid acoustic");
  } // else ** if(!adjoint)
  return;
} // UpdateInteriorVariablesAndComputeResidual

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
