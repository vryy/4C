/*--------------------------------------------------------------------------*/
/*!
\file acou_visc_ele_calc.cpp
\brief

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "acou_visc_ele_calc.H"
#include "acou_ele_calc.H"
#include "acou_ele_action.H"
#include "acou_utils.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_geometry/position_array.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_mat/acoustic_visc.H"

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
DRT::ELEMENTS::AcouViscEleCalc<distype>::AcouViscEleCalc()
{}

/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouViscEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou*    ele,
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
int DRT::ELEMENTS::AcouViscEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou* ele,
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
  if (const DRT::ELEMENTS::AcouVisc * hdgele = dynamic_cast<const DRT::ELEMENTS::AcouVisc*>(ele))
    usescompletepoly_ = hdgele->UsesCompletePolynomialSpace();
  else
    dserror("cannot cast element to acouvisc element");
  InitializeShapes(ele);

  const ACOU::Action action = DRT::INPUT::get<ACOU::Action>(params,"action");
  bool updateonly = false;
  shapes_->Evaluate(*ele);

  switch(action)
  {
  case ACOU::project_field:
  {
    if(mat->MaterialType()!=INPAR::MAT::m_acousticviscmat) dserror("for physical type 'viscous' please supply MAT_AcousticVisc");
    ElementInit(ele,params);
    ProjectField(ele,params,mat,discretization,lm,elevec1,elevec2);
    break;
  }
  case ACOU::project_optical_field:
  {
    if(mat->MaterialType()!=INPAR::MAT::m_acousticviscmat) dserror("for physical type 'viscous' please supply MAT_AcousticVisc");
    ProjectOpticalField(ele,params,mat,discretization,lm,elevec1,elevec2);
    break;
  }
  case ACOU::ele_init:
  {
    ElementInit(ele,params);
    break;
  }
  case ACOU::interpolate_hdg_to_node:
  {
    ReadGlobalVectors(ele,discretization,lm);
    NodeBasedValues(ele,discretization,lm,elevec1);
    break;
  }
  case ACOU::interpolate_psi_to_node:
  {
    double dt = params.get<double>("dt");
    ReadGlobalVectors(ele,discretization,lm);
    NodeBasedPsi(mat,ele,discretization,lm,elevec1,dt);
    break;
  }
  case ACOU::calc_acou_error:
  {
    ReadGlobalVectors(ele,discretization,lm);
    ComputeError(ele,params,mat,discretization,lm,elevec1);
    break;
  }
  case ACOU::calc_abc:
  {
    int face = params.get<int>("face");
    shapesface_->EvaluateFace(ele->Degree(),usescompletepoly_,2*ele->Degree(),*ele,face,*shapes_);
    // note: absorbing bcs are treated fully implicitly!
    localSolver_->ComputeAbsorbingBC(ele,params,mat,face,elemat1,elevec1);
    break;
  }
  case ACOU::calc_pressuremon:
  {
    int face = params.get<int>("face");
    double dt = params.get<double>("dt");

    ComputeMatrices(mat,*ele,dt,dyna_,true);

    shapesface_->EvaluateFace(ele->Degree(),usescompletepoly_,2*ele->Degree(),*ele,face,*shapes_);

    if(!params.isParameter("nodeindices"))
      localSolver_->ComputeSourcePressureMonitor(ele,params,mat,face,elemat1,elevec1);
    else
      localSolver_->ComputeSourcePressureMonitorLine3D(ele,params,mat,face,elemat1,elevec1);

    break;
  }
  case ACOU::calc_pmon_nodevals:
  {
    int face = params.get<int>("face");
    ReadGlobalVectors(ele,discretization,lm);
    shapesface_->EvaluateFace(ele->Degree(),usescompletepoly_,2*ele->Degree(),*ele,face,*shapes_);
    ComputePMonNodeVals(ele,params,mat,face,elemat1,elevec1);

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
      localSolver_->CondenseLocalPart(elemat1,dyna_);

    DIRKVectorHandling(ele,params,dt);

    localSolver_->ComputeResidual(params,*ele,elevec1,interiorGradVelnp_,interiorVelnp_,interiorPressnp_,interiorDensnp_,traceVal_,dyna_);

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

    UpdateInteriorVariablesAndComputeResidual(discretization,params,ele,mat,elevec1,dt,errormaps,updateonly);

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
void DRT::ELEMENTS::AcouViscEleCalc<distype>::InitializeShapes(const DRT::ELEMENTS::Acou* ele)
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
      shapesface_ = Teuchos::rcp(new DRT::UTILS::ShapeValuesFace<distype>());
    // TODO: check distype

    localSolver_ = Teuchos::rcp(new LocalSolver(ele,*shapes_,*shapesface_,usescompletepoly_));

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouViscEleCalc<distype> *
DRT::ELEMENTS::AcouViscEleCalc<distype>::Instance( bool create )
{
  static AcouViscEleCalc<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new AcouViscEleCalc<distype>();
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
void DRT::ELEMENTS::AcouViscEleCalc<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}

/*----------------------------------------------------------------------*
 * Constructor LocalSolver
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouViscEleCalc<distype>::LocalSolver::
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
    shapesface_.EvaluateFace(ele->Degree(),
                             completepoly,
                             2*ele->Degree(),
                             *ele,
                             i,
                             shapes_);
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
  invimat.Shape(ndofs_,ndofs_);
  jmat.Shape(ndofs_,onfdofs);
  kmat.Shape(ndofs_,ndofs_);
  lmat.Shape(ndofs_,ndofs_);
  lhatmat.Shape(ndofs_,ndofs_);
  mmat.Shape(onfdofs,nsd_*nsd_*ndofs_);
  nmat.Shape(onfdofs,nsd_*ndofs_);
  omat.Shape(onfdofs,ndofs_);
  pmat.Shape(onfdofs,onfdofs);
}


/*----------------------------------------------------------------------*
 * Element init
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::ElementInit(DRT::ELEMENTS::Acou* ele,
                                                      Teuchos::ParameterList& params)
{
  shapes_->Evaluate(*ele);
  DRT::ELEMENTS::AcouVisc * hdgele = dynamic_cast<DRT::ELEMENTS::AcouVisc*>(ele);

      // each element has to store the interior vectors by itseld, p-adaptivity or not
  // so, shape it, as you need it
  dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");
  hdgele->eleinteriorValnp_.Shape(shapes_->ndofs_*(nsd_*nsd_+nsd_+2),1);
  switch(dyna_)
  {
  case INPAR::ACOU::acou_bdf4:
    hdgele->eleinteriorValnmmm_.Shape(shapes_->ndofs_*(nsd_*nsd_+nsd_+2),1); // no break here!
  case INPAR::ACOU::acou_bdf3:
    hdgele->eleinteriorValnmm_.Shape(shapes_->ndofs_*(nsd_*nsd_+nsd_+2),1); // no break here!
  case INPAR::ACOU::acou_bdf2:
    hdgele->eleinteriorValnm_.Shape(shapes_->ndofs_*(nsd_*nsd_+nsd_+2),1);
    break;                                                 // here you go!
  case INPAR::ACOU::acou_dirk23:
    hdgele->eles_.resize(2);
    hdgele->elef_.resize(2);
    hdgele->eley_.resize(2);
    break;
  case INPAR::ACOU::acou_dirk33:
  case INPAR::ACOU::acou_dirk34:
    hdgele->eles_.resize(3);
    hdgele->elef_.resize(3);
    hdgele->eley_.resize(3);
    break;
  case INPAR::ACOU::acou_dirk54:
    hdgele->eles_.resize(5);
    hdgele->elef_.resize(5);
    hdgele->eley_.resize(5);
    break;
  default:
    break; // do nothing for impl and trap
  }
  for(unsigned int i=0; i<hdgele->eles_.size(); ++i) // zero size if not dirk
  {
    hdgele->eles_[i].Shape(shapes_->ndofs_*(nsd_*nsd_+nsd_+2),1);
    hdgele->eley_[i].Shape(shapes_->ndofs_*(nsd_*nsd_+nsd_+2),1);
    hdgele->elef_[i].Shape(shapes_->ndofs_*(nsd_*nsd_+nsd_+2),1);
  }
  return;
}


/*----------------------------------------------------------------------*
 * DIRKVectorHandling
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::
DIRKVectorHandling(DRT::ELEMENTS::Acou   * ele,
                   Teuchos::ParameterList&       params,
                   double dt)
{
  if( dyna_ != INPAR::ACOU::acou_dirk23 &&
      dyna_ != INPAR::ACOU::acou_dirk33 &&
      dyna_ != INPAR::ACOU::acou_dirk34 &&
      dyna_ != INPAR::ACOU::acou_dirk54 )
    return; // do nothing for none dirk time integrators

  // when we are here, the time integrator is dirk, and we have to do the update things
  DRT::ELEMENTS::AcouVisc * acouele = dynamic_cast<DRT::ELEMENTS::AcouVisc*>(ele);

  double dirk_a[6][6];
  double dirk_b[6];
  int dirk_q;
  ACOU::FillDIRKValues(dyna_,dirk_a,dirk_b,dirk_q);
  int stage = params.get<int>("stage");

  if(stage == 0)
    acouele->eleinteriorValnp_ = acouele->eleinteriorValn_;


  Epetra_SerialDenseVector tempVec1(shapes_->ndofs_*(nsd_*nsd_+nsd_+2));
  Epetra_SerialDenseVector tempVec2(shapes_->ndofs_*(nsd_*nsd_+nsd_+2));

  acouele->eles_[stage] = acouele->eleinteriorValn_;
  acouele->eles_[stage].Scale(1.0/dt); // dt includes a_ii

  for(int i=0; i<stage; ++i)
  {
    tempVec1 = acouele->eles_[i];
    tempVec1.Scale(-1.0);
    tempVec2 = acouele->eley_[i];
    tempVec2.Scale(1.0/dt);
    tempVec1 += tempVec2;
    tempVec1.Scale(dirk_a[stage][i]/dirk_a[stage][stage]);
    acouele->eles_[stage] += tempVec1;

  }

  // these vectors s are used for the calculation of the residual -> write them to used local solver variable
  for(unsigned i=0; i<interiorValnp_.size(); ++i)
    interiorValnp_[i] = acouele->eles_[stage](i) * dt;
  for(unsigned int i=0; i<shapes_->ndofs_; ++i)
  {
    interiorPressnp_(i) = interiorValnp_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1];
    interiorDensnp_(i)  = interiorValnp_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_];
    for(unsigned int d=0; d<nsd_; ++d)
      interiorVelnp_(i+d*shapes_->ndofs_) = interiorValnp_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+d];
    for(unsigned int d=0; d<nsd_; ++d)
      for(unsigned int e=0; e<nsd_; ++e)
        interiorGradVelnp_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnp_[i*(nsd_*nsd_+nsd_+2)+e*nsd_+d];
  }

  return;
}

/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouViscEleCalc<distype>::ProjectField(
    DRT::ELEMENTS::Acou*                 ele,
    Teuchos::ParameterList&              params,
    Teuchos::RCP<MAT::Material>&         mat,
    DRT::Discretization&                 discretization,
    const std::vector<int>&              lm,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2)
{
  shapes_->Evaluate(*ele);
  DRT::ELEMENTS::AcouVisc * hdgele = dynamic_cast<DRT::ELEMENTS::AcouVisc*>(ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 ||
           unsigned(elevec2.M()) == (nsd_*nsd_+nsd_+2)*shapes_->ndofs_, "Wrong size in project vector 2");

  const MAT::AcousticViscMat* actmat = dynamic_cast<const MAT::AcousticViscMat*>(mat.get());
  double c = actmat->SpeedofSound();

  // get function
  const int *start_func = params.getPtr<int>("funct");

  // internal variables
  if (elevec2.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(shapes_->nqpoints_,1);
    Epetra_SerialDenseMatrix massPart(shapes_->ndofs_,shapes_->nqpoints_);

    for (unsigned int q=0; q<shapes_->nqpoints_; ++q )
    {
      const double fac = shapes_->jfac(q);
      const double sqrtfac = std::sqrt(fac);
      double xyz[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        xyz[d] = shapes_->xyzreal(d,q); // coordinates of quadrature point in real coordinates
      double p;
      dsassert(start_func != NULL,"startfuncno not set for initial value");
      EvaluateAll(*start_func, xyz,  p); // p at quadrature point

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i=0; i<shapes_->ndofs_; ++i)
      {
        massPart(i,q) = shapes_->shfunct(i,q) * sqrtfac;
        localMat(i,0) += shapes_->shfunct(i,q) * p * fac;
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

    for (unsigned int r=0; r<shapes_->ndofs_; ++r )
    {
      //elevec2[r*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1] += localMat(r,0); // pressure
      hdgele->eleinteriorValnp_(r*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1) +=localMat(r,0); // pressure
      //elevec2[r*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_] += localMat(r,0) / c / c; // density
      hdgele->eleinteriorValnp_(r*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_) += localMat(r,0) / c / c; // density
    }
  } // if (elevec2.M() > 0)

  switch(dyna_)
  {
  case INPAR::ACOU::acou_bdf4:
    hdgele->eleinteriorValnmmm_ = hdgele->eleinteriorValnp_; // no break here
  case INPAR::ACOU::acou_bdf3:
    hdgele->eleinteriorValnmm_ = hdgele->eleinteriorValnp_;  // no break here
  case INPAR::ACOU::acou_bdf2:
    hdgele->eleinteriorValnm_ = hdgele->eleinteriorValnp_;
    break;                                                   // here you go
  default:
    hdgele->eleinteriorValn_ = hdgele->eleinteriorValnp_;
    break;
  }

  // trace variable (zero, because no memory, no time derivative)
  // dsassert(elevec1.M() == nfaces_*nfdofs_*nsd_, "Wrong size in project vector 1");
  // elevec1.Scale(0.0);
  if(params.isParameter("faceconsider"))
    dserror("Evaluation of heterogeneous Dirichlet boundary conditions not yet implemented for viscous acoustics");

  return 0;
}

/*----------------------------------------------------------------------*
 * EvaluateAll
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::EvaluateAll(const int start_func,
                                                          const double (&xyz)[nsd_],
                                                          double  &p) const
{
  p = DRT::Problem::Instance()->Funct(start_func-1).Evaluate(0,xyz,0.0,NULL);
  return;
}

/*----------------------------------------------------------------------*
 * ProjectOpticalField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouViscEleCalc<distype>::ProjectOpticalField(
    DRT::ELEMENTS::Acou*                 ele,
    Teuchos::ParameterList&              params,
    Teuchos::RCP<MAT::Material>&         mat,
    DRT::Discretization&                 discretization,
    const std::vector<int>&              lm,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2)
{
  bool meshconform = params.get<bool>("mesh conform");
  DRT::ELEMENTS::AcouVisc * hdgele = dynamic_cast<DRT::ELEMENTS::AcouVisc*>(ele);

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
    const MAT::AcousticViscMat* actmat = static_cast<const MAT::AcousticViscMat*>(mat.get());
    double c = actmat->SpeedofSound();

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
        double p=0.0;
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

      for (unsigned int r=0; r<shapes_->ndofs_; ++r )
      {
        // elevec2[r*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1] += localMat(r,0); // pressure
        // elevec2[r*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_] += localMat(r,0) / c / c; // density
        hdgele->eleinteriorValnp_(r*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1) += localMat(r,0); // pressure
        hdgele->eleinteriorValnp_(r*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_) += localMat(r,0) / c / c; // density
      }
    }

  } // if(meshconform)

  // trace variable (zero, because no memory, no time derivative)
  // dsassert(unsigned(elevec1.M()) == nfaces_*shapes_->nfdofs_*nsd_, "Wrong size in project vector 1");
  // elevec1.Scale(0.0);

  return 0;
}

/*----------------------------------------------------------------------*
 * EvaluateLight
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::EvaluateLight(
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
void DRT::ELEMENTS::AcouViscEleCalc<distype>::
ReadGlobalVectors(DRT::Element     * ele,
                  DRT::Discretization    & discretization,
                  const std::vector<int> & lm)
{
  DRT::ELEMENTS::AcouVisc * acouele = dynamic_cast<DRT::ELEMENTS::AcouVisc*>(ele);

  interiorGradVelnp_.Shape(shapes_->ndofs_*nsd_*nsd_,1);
  interiorDensnp_.Shape(shapes_->ndofs_,1);
  interiorVelnp_.Shape(shapes_->ndofs_*nsd_,1);
  interiorPressnp_.Shape(shapes_->ndofs_,1);

  interiorGradVeln_.Shape(shapes_->ndofs_*nsd_*nsd_,1);
  interiorDensn_.Shape(shapes_->ndofs_,1);
  interiorVeln_.Shape(shapes_->ndofs_*nsd_,1);
  interiorPressn_.Shape(shapes_->ndofs_,1);

  interiorGradVelnm_.Shape(shapes_->ndofs_*nsd_*nsd_,1);
  interiorDensnm_.Shape(shapes_->ndofs_,1);
  interiorVelnm_.Shape(shapes_->ndofs_*nsd_,1);
  interiorPressnm_.Shape(shapes_->ndofs_,1);

  interiorGradVelnmm_.Shape(shapes_->ndofs_*nsd_*nsd_,1);
  interiorDensnmm_.Shape(shapes_->ndofs_,1);
  interiorVelnmm_.Shape(shapes_->ndofs_*nsd_,1);
  interiorPressnmm_.Shape(shapes_->ndofs_,1);

  interiorGradVelnmmm_.Shape(shapes_->ndofs_*nsd_*nsd_,1);
  interiorDensnmmm_.Shape(shapes_->ndofs_,1);
  interiorVelnmmm_.Shape(shapes_->ndofs_*nsd_,1);
  interiorPressnmmm_.Shape(shapes_->ndofs_,1);

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

  // read the internal vectors
  {
    interiorValnp_.resize(shapes_->ndofs_*(nsd_*nsd_+nsd_+2));
    for(unsigned i=0; i<interiorValnp_.size(); ++i)
      interiorValnp_[i] = acouele->eleinteriorValnp_(i);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorPressnp_(i) = interiorValnp_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1];
      interiorDensnp_(i)  = interiorValnp_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_];
      for(unsigned int d=0; d<nsd_; ++d)
        interiorVelnp_(i+d*shapes_->ndofs_) = interiorValnp_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+d];
      for(unsigned int d=0; d<nsd_; ++d)
        for(unsigned int e=0; e<nsd_; ++e)
          interiorGradVelnp_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnp_[i*(nsd_*nsd_+nsd_+2)+e*nsd_+d];
    }
    if(acouele->eleinteriorValnm_.M()>0) // bdf2 and bdf3 need this
    {
      interiorValnm_.resize(shapes_->ndofs_*(nsd_*nsd_+nsd_+2));
      for(unsigned i=0; i<interiorValnp_.size(); ++i)
        interiorValnm_[i] = acouele->eleinteriorValnm_(i);
      for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      {
        interiorPressnm_(i) = interiorValnm_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1];
        interiorDensnm_(i)  = interiorValnm_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_];
        for(unsigned int d=0; d<nsd_; ++d)
          interiorVelnm_(i+d*shapes_->ndofs_) = interiorValnm_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+d];
        for(unsigned int d=0; d<nsd_; ++d)
          for(unsigned int e=0; e<nsd_; ++e)
            interiorGradVelnm_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnm_[i*(nsd_*nsd_+nsd_+2)+e*nsd_+d];
      }
    }
    if(acouele->eleinteriorValnmm_.M()>0) // bdf2 and bdf3 need this
    {
      interiorValnmm_.resize(shapes_->ndofs_*(nsd_*nsd_+nsd_+2));
      for(unsigned i=0; i<interiorValnp_.size(); ++i)
        interiorValnmm_[i] = acouele->eleinteriorValnmm_(i);
      for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      {
        interiorPressnmm_(i) = interiorValnmm_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1];
        interiorDensnmm_(i)  = interiorValnmm_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_];
        for(unsigned int d=0; d<nsd_; ++d)
          interiorVelnmm_(i+d*shapes_->ndofs_) = interiorValnmm_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+d];
        for(unsigned int d=0; d<nsd_; ++d)
          for(unsigned int e=0; e<nsd_; ++e)
            interiorGradVelnmm_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnmm_[i*(nsd_*nsd_+nsd_+2)+e*nsd_+d];
      }
    }
    if(acouele->eleinteriorValnmmm_.M()>0) // bdf2 and bdf3 need this
    {
      interiorValnmmm_.resize(shapes_->ndofs_*(nsd_*nsd_+nsd_+2));
      for(unsigned i=0; i<interiorValnp_.size(); ++i)
        interiorValnmmm_[i] = acouele->eleinteriorValnmmm_(i);
      for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      {
        interiorPressnmmm_(i) = interiorValnmmm_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1];
        interiorDensnmmm_(i)  = interiorValnmmm_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_];
        for(unsigned int d=0; d<nsd_; ++d)
          interiorVelnmmm_(i+d*shapes_->ndofs_) = interiorValnmmm_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+d];
        for(unsigned int d=0; d<nsd_; ++d)
          for(unsigned int e=0; e<nsd_; ++e)
            interiorGradVelnmmm_(i+(e*nsd_+d)*shapes_->ndofs_) = interiorValnmmm_[i*(nsd_*nsd_+nsd_+2)+e*nsd_+d];
      }
    }
  } // read the internal vectors
  return;
}

/*----------------------------------------------------------------------*
 * ComputeError
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::ComputeError(
    DRT::ELEMENTS::Acou*          ele,
    Teuchos::ParameterList&       params,
    Teuchos::RCP<MAT::Material>&  mat,
    DRT::Discretization&          discretization,
    const std::vector<int>&       lm,
    Epetra_SerialDenseVector&     elevec)
{
  shapes_->Evaluate(*ele);

  double time = params.get<double>("time");

  if(params.get<INPAR::ACOU::CalcError>("error calculation") != INPAR::ACOU::calcerror_1d) dserror("no analytical solution available");

  const MAT::AcousticViscMat* actmat = static_cast<const MAT::AcousticViscMat*>(mat.get());
  if ( actmat->Viscosity()!= 0.0 ) dserror("no analytical solution available for viscosity != 0.0");
  if ( actmat->Therm()!= 0.0 ) dserror("no analytical solution available for thermal coefficient != 0.0");
  double c = actmat->SpeedofSound();

  // get function
  const int *start_func = params.getPtr<int>("startfuncno");

  double err_p = 0.0, norm_p = 0.0;
  for (unsigned int q=0; q<shapes_->nqpoints_; ++q)
  {
    double numerical = 0.0;
    for (unsigned int i=0; i<shapes_->ndofs_; ++i)
      numerical += shapes_->shfunct(i,q) * interiorPressnp_(i);

    double exact = 0.0;
    double xyz[nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
      xyz[d] = shapes_->xyzreal(d,q);

    // careful: we assume that the x-axis is the orientation of your 1D problem
    double xyzp[nsd_];
    double xyzm[nsd_];
    xyzp[0]=xyz[0]+c*time;
    xyzm[0]=xyz[0]-c*time;
    for (unsigned int d=1; d<nsd_; ++d)
    {
      xyzp[d] = xyz[d];
      xyzm[d] = xyz[d];
    }

    // calculation of the analytical solution
    exact = 0.5 * DRT::Problem::Instance()->Funct(*start_func-1).Evaluate(0,xyzp,0.0,NULL)
          + 0.5 * DRT::Problem::Instance()->Funct(*start_func-1).Evaluate(0,xyzm,0.0,NULL);

    err_p += ( exact - numerical ) * ( exact - numerical ) * shapes_->jfac(q);
    norm_p += exact * exact * shapes_->jfac(q);
  }

  elevec[1] += err_p;
  elevec[3] += norm_p;

  return;
}

/*----------------------------------------------------------------------*
 * NodeBasedValues
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::NodeBasedValues(
                          DRT::ELEMENTS::Acou*                 ele,
                          DRT::Discretization&                 discretization,
                          const std::vector<int>&              lm,
                          Epetra_SerialDenseVector&            elevec1)
{

  dsassert(elevec1.M() == (int)nen_*(2*nsd_+2+6)+2, "Vector does not have correct size");
  elevec1.Scale(0.0);
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  double norm_d = 0.0, norm_p = 0.0, area = 0.0;
  for (unsigned int q=0; q<shapes_->nqpoints_; ++q)
  {
    double numerical_p = 0.0, numerical_d = 0.0;
    for (unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      numerical_p += shapes_->shfunct(i,q) * interiorPressnp_(i);
      numerical_d += shapes_->shfunct(i,q) * interiorDensnp_(i);
    }
    norm_d += numerical_d * shapes_->jfac(q);
    norm_p += numerical_p * shapes_->jfac(q);
    area += shapes_->jfac(q);
  }

  elevec1((2*nsd_+2)*nen_) = norm_p/area;
  elevec1((2*nsd_+2)*nen_+1) = norm_d/area;

  for (unsigned int i=0; i<nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim=0;idim<nsd_;idim++)
      shapes_->xsi(idim) = locations(idim,i);
    shapes_->polySpace_->Evaluate(shapes_->xsi,values);

    // compute values for velocity and pressure by summing over all basis functions
    double sump = 0;
    double sumd = 0;
    for (unsigned int k=0; k<shapes_->ndofs_; ++k)
    {
      sump += values(k) * interiorPressnp_(k);
      sumd += values(k) * interiorDensnp_(k);
    }
    elevec1(2*nsd_*nen_+i) = sump;
    elevec1((2*nsd_+1)*nen_+i) = sumd;

    for (unsigned int d=0; d<nsd_; ++d)
    {
      double sumv = 0.0;
      for (unsigned int k=0; k<shapes_->ndofs_; ++k)
        sumv += values(k) * interiorVelnp_(d*shapes_->ndofs_+k);
      elevec1(d*nen_+i) = sumv;
    }
    for (unsigned int d=0; d<nsd_*nsd_; ++d)
    {
      double sumvg = 0.0;
      for (unsigned int k=0; k<shapes_->ndofs_; ++k)
        sumvg += values(k) * interiorGradVelnp_(d*shapes_->ndofs_+k);

      if(d>6) break;
      elevec1(nen_*(2*nsd_+2)+2+d*nen_+i) = sumvg;
    }
  }

  locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape);
  Epetra_SerialDenseVector fvalues(shapesface_->nfdofs_);
  for (unsigned int face=0; face<nfaces_; ++face)
  {
    for (int i=0; i<DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
    {
      // evaluate shape polynomials in node
      for (unsigned int idim=0;idim<nsd_-1;idim++)
        shapesface_->xsiF(idim) = locations(idim,i);
      shapesface_->polySpaceFace_->Evaluate(shapesface_->xsiF,fvalues); // TODO: fix face orientation here

      // compute values for velocity and pressure by summing over all basis functions
      for (unsigned int d=0; d<nsd_; ++d)
      {
        double sum = 0;
        for (unsigned int k=0; k<shapesface_->nfdofs_; ++k)
          sum += fvalues(k) * traceVal_[face*shapesface_->nfdofs_*nsd_+d*shapesface_->nfdofs_+k];
        elevec1((nsd_+d)*nen_+shapesface_->faceNodeOrder[face][i]) = sum;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 * NodeBasedPsi
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::NodeBasedPsi(
                          const Teuchos::RCP<MAT::Material>    &mat,
                          DRT::ELEMENTS::Acou*                 ele,
                          DRT::Discretization&                 discretization,
                          const std::vector<int>&              lm,
                          Epetra_SerialDenseVector&            elevec1,
                          double                               dt)
{
  dsassert(elevec1.M() == int(nen_), "Vector does not have correct size");
  elevec1.Scale(0.0);
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  // calculate mass matrix
  const MAT::AcousticViscMat* actmat = static_cast<const MAT::AcousticViscMat*>(mat.get());
  double c = actmat->SpeedofSound();
  double ktherm = actmat->Therm();

  for (unsigned int i=0; i<nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim=0;idim<nsd_;idim++)
      shapes_->xsi(idim) = locations(idim,i);
    shapes_->polySpace_->Evaluate(shapes_->xsi,values);

    // compute values for velocity and pressure by summing over all basis functions
    double sump = 0.0;
    double sumr = 0.0;
    for (unsigned int k=0; k<shapes_->ndofs_; ++k)
    {
      sump += values(k) * interiorPressnp_(k);
      sumr += values(k) * interiorDensnp_(k);
    }
    elevec1(i) = sump * ktherm / dt + sumr / c / c / dt;
  }

  return;
}

/*----------------------------------------------------------------------*
 * Compute internal and face matrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::ComputeMatrices(const Teuchos::RCP<MAT::Material> &mat,
                                                                           DRT::ELEMENTS::Acou &             ele,
                                                                           double                            dt,
                                                                           INPAR::ACOU::DynamicType          dyna,
                                                                           bool                              adjoint)
{
  const MAT::AcousticViscMat* actmat = static_cast<const MAT::AcousticViscMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();
  double tau = 1.0;//1.0;///dt;//
  double visc = actmat->Viscosity();
  double therm = actmat->Therm();

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
  zeroMatrix(localSolver_->invimat);
  zeroMatrix(localSolver_->jmat);
  zeroMatrix(localSolver_->kmat);
  zeroMatrix(localSolver_->lmat);
  zeroMatrix(localSolver_->lhatmat);
  zeroMatrix(localSolver_->mmat);
  zeroMatrix(localSolver_->nmat);
  zeroMatrix(localSolver_->omat);
  zeroMatrix(localSolver_->pmat);

  switch(dyna_)
  {
  case INPAR::ACOU::acou_bdf4: dt *= 12.0/25.0; break;
  case INPAR::ACOU::acou_bdf3: dt *= 6.0/11.0; break;
  case INPAR::ACOU::acou_bdf2: dt *= 2.0/3.0; break;
  default: break;
  }

  localSolver_->ComputeInteriorMatrices(mat,dt,dyna_);
  for (unsigned int face=0; face<nfaces_; ++face)
  {
    shapesface_->EvaluateFace(ele.Degree(),usescompletepoly_,2*ele.Degree(),ele, face, *shapes_);
//    std::cout<<"ele "<<ele.Id()<<" face "<<face<<std::endl;
//    shapes_.normals.Print(std::cout);
//    const int* nodeids = ele.Faces()[face]->NodeIds();
//    std::cout<<"nodeids ";
//    for(int i=0; i<ele.Faces()[face]->NumNode(); ++i)
//      std::cout<<nodeids[i]<<" ";
//    std::cout<<std::endl;

    localSolver_->ComputeFaceMatrices(face,mat,dt);
  }
  // scale the matrices

  localSolver_->cmat.Scale(-1.0);
  localSolver_->dmat.Scale(-visc);
  localSolver_->emat.Scale(rho/dt);
  localSolver_->ehatmat.Scale(tau);
  localSolver_->gmat.Scale(-tau);
  localSolver_->hmat.Scale(-rho);
  localSolver_->imat.Scale(1.0/dt);
  localSolver_->invimat.Scale(dt);
  localSolver_->jmat.Scale(rho);
  localSolver_->lmat.Scale(therm/rho/c/c/dt);
  localSolver_->lhatmat.Scale(-1.0/c/c);
  localSolver_->mmat.Scale(-visc);
  localSolver_->nmat.Scale(tau);
  localSolver_->pmat.Scale(-tau);

  return;
}


/*----------------------------------------------------------------------*
 * ComputeInteriorMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::LocalSolver::
ComputeInteriorMatrices(const Teuchos::RCP<MAT::Material> &mat,
                        double dt,
                        INPAR::ACOU::DynamicType dyna)
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
        for (unsigned int e=0; e<nsd_; ++e)
        {
          bmat(i+(e+d*nsd_)*ndofs_,j+e*ndofs_) = sums[d];
          dmat(j+e*ndofs_,i+(e+d*nsd_)*ndofs_) = sums[d];
        }
      }
    }
  }

  // here, i is ready, and we can update related matrices!
  kmat = imat;
  lmat = imat;
  lhatmat = imat;

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

  // we will also need the inverse of imat, however, since this matrix is not taht big,
  // the element does not have to store it
  invimat = imat;
  {
    Epetra_SerialDenseSolver inverseimat;
    inverseimat.SetMatrix(invimat);
    inverseimat.Invert();
  }
  return;
}

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::LocalSolver::
ComputeFaceMatrices(const int                          face,
                    const Teuchos::RCP<MAT::Material> &mat,
                    double dt)
{

  // missing term in ehat
  for (unsigned int p=0; p<ndofs_; ++p)
  {
    for (unsigned int q=0; q<=p; ++q)
    {
      double tempehat = 0.0;
      for (unsigned int i=0; i<shapesface_.nfqpoints_; ++i)
      {
        tempehat += shapesface_.jfacF(i) * shapesface_.shfunctI[face](p,i) * shapesface_.shfunctI[face](q,i);
      }
      for (unsigned int d=0; d<nsd_; ++d)
        ehatmat(d*ndofs_+p,d*ndofs_+q) = ehatmat(d*ndofs_+q,d*ndofs_+p) += tempehat;
    }
  }

  // p
  for (unsigned int p=0; p<shapesface_.nfdofs_; ++p)
  {
    for (unsigned int q=0; q<=p; ++q)
    {
      double tempp = 0.0;
      for (unsigned int i=0; i<shapesface_.nfqpoints_; ++i)
        tempp += shapesface_.jfacF(i) * shapesface_.shfunctF(p,i) * shapesface_.shfunctF(q,i);
      for (unsigned int d=0; d<nsd_; ++d)
      {
        pmat(face*shapesface_.nfdofs_*nsd_+p+d*shapesface_.nfdofs_, face*shapesface_.nfdofs_*nsd_+q+d*shapesface_.nfdofs_) = tempp;
        pmat(face*shapesface_.nfdofs_*nsd_+q+d*shapesface_.nfdofs_, face*shapesface_.nfdofs_*nsd_+p+d*shapesface_.nfdofs_) = tempp;
      }
    }
  }

  // c, g, j, m, n, o
  for (unsigned int p=0; p<shapesface_.nfdofs_; ++p)
  {
    for (unsigned int q=0; q<ndofs_; ++q)
    {
      double tempmat = 0.0;
      for (unsigned int i=0; i<shapesface_.nfqpoints_; ++i)
      {
        double temp = shapesface_.jfacF(i) * shapesface_.shfunctF(p,i) * shapesface_.shfunctI[face](q,i);
        tempmat += temp;
        for(unsigned int j=0; j<nsd_; ++j)
        {
          double temp_d = temp*shapesface_.normals(j,i);
          jmat(q,p+j*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_) += temp_d;
          omat(p+j*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_,q) += temp_d;

          for(unsigned int e=0; e<nsd_; ++e)
          {
            // mmat(p+e*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_,q+(e*nsd_+j)*ndofs_) += temp_d;
            // cmat(q+(e*nsd_+j)*ndofs_,p+e*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_) += temp_d;
            mmat(p+e*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_,q+(j*nsd_+e)*ndofs_) += temp_d;
            cmat(q+(j*nsd_+e)*ndofs_,p+e*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_) += temp_d;
          }
        }
      }
      for(unsigned int j=0; j<nsd_; ++j)
      {
        gmat(q+j*ndofs_,p+j*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_) = tempmat;
        nmat(p+j*shapesface_.nfdofs_+face*shapesface_.nfdofs_*nsd_,q+j*ndofs_) = tempmat;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::LocalSolver::
CondenseLocalPart(Epetra_SerialDenseMatrix &eleMat,
                  INPAR::ACOU::DynamicType dyna)
{
  {
  double theta = 1.0;
  if(dyna==INPAR::ACOU::acou_trapezoidal) theta = 0.66;

  Epetra_SerialDenseMatrix dinvamat(nsd_*ndofs_,nsd_*nsd_*ndofs_);
  dinvamat.Multiply('N','N',1.0,dmat,invamat,0.0);

  Epetra_SerialDenseMatrix ol(nsd_*ndofs_,nsd_*ndofs_);
  ol = ehatmat;
  ol.Scale(theta);
  ol += emat;
  ol.Multiply('N','N',-theta,dinvamat,bmat,1.0);

  Epetra_SerialDenseMatrix toinv((nsd_+2)*ndofs_,(nsd_+2)*ndofs_);
  for(unsigned int i=0; i<ndofs_; ++i)
  {
    for(unsigned int j=0; j<ndofs_; ++j)
    {
      toinv(nsd_*ndofs_+i,nsd_*ndofs_+j) = imat(i,j);
      toinv((nsd_+1)*ndofs_+i,(nsd_+1)*ndofs_+j) = lmat(i,j) + theta * lhatmat(i,j) ;
      toinv((nsd_+1)*ndofs_+i,nsd_*ndofs_+j) = theta * kmat(i,j);
      for(unsigned int d=0; d<nsd_; ++d)
      {
        toinv(d*ndofs_+i,d*ndofs_+j) = ol(d*ndofs_+i,d*ndofs_+j);
        toinv(nsd_*ndofs_+i,d*ndofs_+j) = theta * hmat(i,d*ndofs_+j);
        toinv(d*ndofs_+i,(nsd_+1)*ndofs_+j) = theta * fmat(d*ndofs_+i,j);
      }
    }
  }
  Epetra_SerialDenseMatrix rhsinv((nsd_+2)*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);

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
      tempmat1(j,i) = rhsinv((nsd_+1)*ndofs_+j,i);

  eleMat = pmat;
  eleMat.Multiply('N','N',-theta,omat,tempmat1,theta);

  tempmat1.Shape(nsd_*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
    for(unsigned int j=0; j<nsd_*ndofs_; ++j)
      tempmat1(j,i) = rhsinv(j,i);

  eleMat.Multiply('N','N',-theta,nmat,tempmat1,1.0);

  ol.Shape(nsd_*nsd_*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
  ol = cmat;
  ol.Multiply('N','N',-theta,bmat,tempmat1,theta);
  tempmat1.Shape(nsd_*nsd_*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
  tempmat1.Multiply('N','N',1.0/theta,invamat,ol,0.0);
  eleMat.Multiply('N','N',-theta,mmat,tempmat1,1.0);
  }

//  LINALG::Matrix<(nsd_*nsd_+nsd_+2)*ndofs_,(nsd_*nsd_+nsd_+2)*ndofs_> toinv(true);
//  for(unsigned int i=0; i<ndofs_*nsd_*nsd_; ++i)
//    for(unsigned int j=0; j<ndofs_*nsd_*nsd_; ++j)
//      toinv(i,j) = amat(j,i);
//  for(unsigned int i=0; i<ndofs_*nsd_*nsd_; ++i)
//    for(unsigned int j=0; j<ndofs_*nsd_; ++j)
//      toinv(i,nsd_*nsd_*ndofs_+j) = dmat(j,i);
//  for(unsigned int i=0; i<ndofs_*nsd_*nsd_; ++i)
//    for(unsigned int j=0; j<ndofs_*nsd_; ++j)
//      toinv(nsd_*nsd_*ndofs_+j,i) = bmat(i,j);
//  for(unsigned int i=0; i<ndofs_*nsd_; ++i)
//    for(unsigned int j=0; j<ndofs_*nsd_; ++j)
//      toinv(nsd_*nsd_*ndofs_+i,nsd_*nsd_*ndofs_+j) = emat(j,i) + ehatmat(j,i);
//  for(unsigned int i=0; i<ndofs_*nsd_; ++i)
//    for(unsigned int j=0; j<ndofs_; ++j)
//      toinv(nsd_*nsd_*ndofs_+i,(nsd_*nsd_+nsd_)*ndofs_+j) = hmat(j,i);
//  for(unsigned int i=0; i<ndofs_; ++i)
//    for(unsigned int j=0; j<ndofs_; ++j)
//      toinv((nsd_*nsd_+nsd_)*ndofs_+i,(nsd_*nsd_+nsd_)*ndofs_+j) = imat(j,i);
//  for(unsigned int i=0; i<ndofs_; ++i)
//    for(unsigned int j=0; j<ndofs_; ++j)
//      toinv((nsd_*nsd_+nsd_)*ndofs_+i,(nsd_*nsd_+nsd_+1)*ndofs_+j) = kmat(j,i);
//  for(unsigned int i=0; i<ndofs_; ++i)
//    for(unsigned int j=0; j<ndofs_*nsd_; ++j)
//      toinv((nsd_*nsd_+nsd_+1)*ndofs_+i,nsd_*nsd_*ndofs_+j) = fmat(j,i);
//  for(unsigned int i=0; i<ndofs_; ++i)
//    for(unsigned int j=0; j<ndofs_; ++j)
//      toinv((nsd_*nsd_+nsd_+1)*ndofs_+i,(nsd_*nsd_+nsd_+1)*ndofs_+j) = lmat(j,i) +lhatmat(j,i);
//  toinv.Print(std::cout);
//
//  LINALG::Matrix<(nsd_*nsd_+nsd_+2)*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_> rhsinv(true);
//  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
//    for(unsigned int j=0; j<ndofs_*nsd_*nsd_; ++j)
//      rhsinv(j,i) = mmat(i,j);
//  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
//    for(unsigned int j=0; j<ndofs_*nsd_; ++j)
//      rhsinv(nsd_*nsd_*ndofs_+j,i) = nmat(i,j);
//  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
//    for(unsigned int j=0; j<ndofs_; ++j)
//      rhsinv((nsd_*nsd_+nsd_+1)*ndofs_+j,i) = omat(i,j);
//  rhsinv.Print(std::cout);
//
//  // invert
//  {
//      LINALG::FixedSizeSerialDenseSolver<(nsd_*nsd_+nsd_+2)*ndofs_,(nsd_*nsd_+nsd_+2)*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_> inverse;
//      inverse.SetMatrix(toinv);
//      inverse.SetVectors(rhsinv,rhsinv);
//      inverse.Solve();
//  }
//  rhsinv.Print(std::cout);
//
//  eleMat = pmat;
//
//  Epetra_SerialDenseMatrix tempmat1(nsd_*nsd_*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
//  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
//    for(unsigned int j=0; j<ndofs_*nsd_*nsd_; ++j)
//      tempmat1(j,i) = rhsinv(j,i);
//  tempmat1.Print(std::cout);
//  eleMat.Multiply('T','N',-1.0,cmat,tempmat1,1.0);
//
//  tempmat1.Shape(nsd_*ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
//  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
//    for(unsigned int j=0; j<ndofs_*nsd_; ++j)
//      tempmat1(j,i) = rhsinv(nsd_*nsd_*ndofs_+j,i);
//  tempmat1.Print(std::cout);
//  eleMat.Multiply('T','N',-1.0,gmat,tempmat1,1.0);
//
//  tempmat1.Shape(ndofs_,nsd_*shapesface_.nfdofs_*nfaces_);
//  for(unsigned int i=0; i<nsd_*shapesface_.nfdofs_*nfaces_; ++i)
//    for(unsigned int j=0; j<ndofs_; ++j)
//      tempmat1(j,i) = rhsinv((nsd_*nsd_+nsd_)*ndofs_+j,i);
//  tempmat1.Print(std::cout);
//  eleMat.Multiply('T','N',-1.0,jmat,tempmat1,1.0);
//
//  std::cout<<"Version B"<<std::endl;
//  eleMat.Print(std::cout);
//  std::cin.get();
  return;
}

/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::LocalSolver::
ComputeResidual(Teuchos::ParameterList&     params,
                DRT::ELEMENTS::Acou &                ele,
                Epetra_SerialDenseVector          & elevec,
                Epetra_SerialDenseVector          & interiorGradVeln,
                Epetra_SerialDenseVector          & interiorVeln,
                Epetra_SerialDenseVector          & interiorPressn,
                Epetra_SerialDenseVector          & interiorDensn,
                std::vector<double>               traceVal,
                INPAR::ACOU::DynamicType          dyna)
{
  bool adjoint = params.get<bool>("adjoint");
  double theta = 1.0;
  if(dyna==INPAR::ACOU::acou_trapezoidal) theta = 0.66;

  Epetra_SerialDenseVector traceVal_SDV(nfaces_*shapesface_.nfdofs_);
  for(unsigned i=0; i<nfaces_*shapesface_.nfdofs_; ++i)
    traceVal_SDV(i) = traceVal[i];
  if(!adjoint)
  {
    Epetra_SerialDenseMatrix dinvamat(nsd_*ndofs_,nsd_*nsd_*ndofs_);
    dinvamat.Multiply('N','N',1.0,dmat,invamat,0.0);

    Epetra_SerialDenseMatrix ol(nsd_*ndofs_,nsd_*ndofs_);
    ol = ehatmat;
    ol.Scale(theta);
    ol += emat;
    ol.Multiply('N','N',-theta,dinvamat,bmat,1.0);

    Epetra_SerialDenseMatrix toinv((nsd_+2)*ndofs_,(nsd_+2)*ndofs_);
    for(unsigned int i=0; i<ndofs_; ++i)
    {
      for(unsigned int j=0; j<ndofs_; ++j)
      {
        toinv(nsd_*ndofs_+i,nsd_*ndofs_+j) = imat(i,j);
        toinv((nsd_+1)*ndofs_+i,(nsd_+1)*ndofs_+j) = lmat(i,j) + theta * lhatmat(i,j) ;
        toinv((nsd_+1)*ndofs_+i,nsd_*ndofs_+j) = theta * kmat(i,j);
        for(unsigned int d=0; d<nsd_; ++d)
        {
          toinv(d*ndofs_+i,d*ndofs_+j) = ol(d*ndofs_+i,d*ndofs_+j);
          toinv(nsd_*ndofs_+i,d*ndofs_+j) = theta * hmat(i,d*ndofs_+j);
          toinv(d*ndofs_+i,(nsd_+1)*ndofs_+j) = theta * fmat(d*ndofs_+i,j);
        }
      }
    }

    Epetra_SerialDenseVector rhsinv((nsd_+2)*ndofs_);

    ol.Shape(nsd_*ndofs_,1);
    ol.Multiply('N','N',1.0,emat,interiorVeln,0.0);
    ol.Multiply('N','N',-(1.0-theta),dmat,interiorGradVeln,1.0);
    ol.Multiply('N','N',-(1.0-theta),ehatmat,interiorVeln,1.0);
    ol.Multiply('N','N',-(1.0-theta),fmat,interiorPressn,1.0);
    ol.Multiply('N','N',-(1.0-theta),gmat,traceVal_SDV,1.0);
    Epetra_SerialDenseVector tempvec1(nsd_*nsd_*ndofs_);
    tempvec1.Multiply('N','N',-(1.0-theta),amat,interiorGradVeln,0.0);
    tempvec1.Multiply('N','N',-(1.0-theta),bmat,interiorVeln,1.0);
    tempvec1.Multiply('N','N',-(1.0-theta),cmat,traceVal_SDV,1.0);
    ol.Multiply('N','N',-1.0,dinvamat,tempvec1,1.0);
    for(unsigned int i=0; i<nsd_*ndofs_; ++i)
      rhsinv(i) = ol(i,0);

    ol.Shape(ndofs_,1);
    ol.Multiply('N','N',1.0,imat,interiorDensn,0.0);
    ol.Multiply('N','N',-(1.0-theta),hmat,interiorVeln,1.0);
    ol.Multiply('N','N',-(1.0-theta),jmat,traceVal_SDV,1.0);
    for(unsigned int i=0; i<ndofs_; ++i)
      rhsinv(nsd_*ndofs_+i) = ol(i,0);

    ol.Multiply('N','N',1.0,lmat,interiorPressn,0.0);
    ol.Multiply('N','N',-(1.0-theta),kmat,interiorDensn,1.0);
    ol.Multiply('N','N',-(1.0-theta),lhatmat,interiorPressn,1.0);
    for(unsigned int i=0; i<ndofs_; ++i)
      rhsinv((nsd_+1)*ndofs_+i) = ol(i,0);

    // invert
    {
      Epetra_SerialDenseSolver inverse;
      inverse.SetMatrix(toinv);
      inverse.SetVectors(rhsinv,rhsinv);
      inverse.Solve();
    }

    elevec.Multiply('N','N',-(1.0-theta),mmat,interiorGradVeln,0.0);
    elevec.Multiply('N','N',-(1.0-theta),nmat,interiorVeln,1.0);
    elevec.Multiply('N','N',-(1.0-theta),omat,interiorPressn,1.0);
    elevec.Multiply('N','N',-(1.0-theta),pmat,traceVal_SDV,1.0);

    for(unsigned int i=0; i<ndofs_; ++i)
      ol(i,0) = rhsinv((nsd_+1)*ndofs_+i);
    elevec.Multiply('N','N',-theta,omat,ol,1.0);

    ol.Shape(nsd_*ndofs_,1);
    for(unsigned int i=0; i<nsd_*ndofs_; ++i)
      ol(i,0) = rhsinv(i);
    elevec.Multiply('N','N',-theta,nmat,ol,1.0);

    tempvec1.Multiply('N','N',-theta,bmat,ol,1.0);
    ol.Shape(nsd_*nsd_*ndofs_,1);
    ol.Multiply('N','N',1.0/theta,invamat,tempvec1,0.0);
    elevec.Multiply('N','N',-theta,mmat,ol,1.0);
  }
  else
  {
    Epetra_SerialDenseMatrix toinv((nsd_*nsd_+nsd_+2)*ndofs_,(nsd_*nsd_+nsd_+2)*ndofs_);

    for(unsigned int i=0; i<ndofs_; ++i)
    {
      for(unsigned int j=0; j<ndofs_; ++j)
      {
        toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+j) = imat(j,i);
        toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+ndofs_+j) = kmat(j,i);
        toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+ndofs_+j) = lmat(j,i) + lhatmat(j,i);
        for(unsigned int d=0; d<nsd_; ++d)
        {
          toinv((nsd_*nsd_)*ndofs_+d*ndofs_+i,(nsd_*nsd_)*ndofs_+d*ndofs_+j) = emat(d*ndofs_+j,d*ndofs_+i) + ehatmat(d*ndofs_+j,d*ndofs_+i);
          toinv((nsd_*nsd_)*ndofs_+d*ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+j) = hmat(j,d*ndofs_+i);
          toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+ndofs_+i,(nsd_*nsd_)*ndofs_+d*ndofs_+j) = fmat(d*ndofs_+j,i);
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
    Epetra_SerialDenseVector rhsinv((nsd_*nsd_+nsd_+2)*ndofs_);

    Epetra_SerialDenseMatrix ol(nsd_*ndofs_,1);
    ol.Multiply('N','N',1.0,emat,interiorVeln,0.0);
    for(unsigned int i=0; i<nsd_*ndofs_; ++i)
      rhsinv(nsd_*nsd_*ndofs_+i) = ol(i,0);

    ol.Shape(ndofs_,1);
    ol.Multiply('N','N',1.0,imat,interiorDensn,0.0);
    for(unsigned int i=0; i<ndofs_; ++i)
      rhsinv(nsd_*nsd_*ndofs_+nsd_*ndofs_+i) = ol(i,0);

    ol.Multiply('N','N',1.0,lmat,interiorPressn,0.0);
    for(unsigned int i=0; i<ndofs_; ++i)
      rhsinv((nsd_*nsd_+nsd_+1)*ndofs_+i) = ol(i,0);

    // invert
    {
      Epetra_SerialDenseSolver inverse;
      inverse.SetMatrix(toinv);
      inverse.SetVectors(rhsinv,rhsinv);
      inverse.Solve();
    }

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

  } // else ** if(adjoint)

  return;
}


/*----------------------------------------------------------------------*
 * ComputeAbsorbingBC
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::LocalSolver::
ComputeAbsorbingBC(DRT::ELEMENTS::Acou*        ele,
                   Teuchos::ParameterList&     params,
                   Teuchos::RCP<MAT::Material> & mat,
                   int                         face,
                   Epetra_SerialDenseMatrix    &elemat,
                   Epetra_SerialDenseVector    &elevec)
{
  const MAT::AcousticViscMat* actmat = static_cast<const MAT::AcousticViscMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();
  bool resonly = params.get<bool>("resonly");

  if(!resonly)
  {
    for (unsigned int p=0; p<shapesface_.nfdofs_; ++p)
    {
      for (unsigned int q=0; q<shapesface_.nfdofs_; ++q)
      {
        double temp[nsd_*nsd_];
        for(unsigned int i=0; i<nsd_*nsd_; ++i)
          temp[i] = 0.0;
        for(unsigned int i=0; i<shapesface_.nfqpoints_; ++i)
          for(unsigned int d=0; d<nsd_; ++d)
            for(unsigned int e=0; e<nsd_; ++e)
              temp[d+e*nsd_] += shapesface_.jfacF(i) * shapesface_.shfunctF(p,i) * shapesface_.shfunctF(q,i) * shapesface_.normals(d,i) * shapesface_.normals(e,i);

        for (unsigned int d=0; d<nsd_; ++d)
        {
          for(unsigned int e=0; e<nsd_; ++e)
          {
            elemat(face*shapesface_.nfdofs_*nsd_+p+d*shapesface_.nfdofs_, face*shapesface_.nfdofs_*nsd_+q+e*shapesface_.nfdofs_) -= rho * c * temp[d+e*nsd_];
          }
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 * ComputeSourcePressureMonitor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::LocalSolver::
ComputeSourcePressureMonitor(DRT::ELEMENTS::Acou*        ele,
                   Teuchos::ParameterList&     params,
                   Teuchos::RCP<MAT::Material> & mat,
                   int                         face,
                   Epetra_SerialDenseMatrix    &elemat,
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
          fnodexyz[i][d] = shapesface_.xyzeF(d,i);

      } // for(int i=0; i<numfnode; ++i)

      // so the source term is face based, hence, we calculate the face contribution
      // just as in AcouEleCalc::LocalSolver::ComputeSourcePressureMonitor
      // but then have to deal with the complement thing
      Epetra_SerialDenseMatrix mass(shapesface_.nfdofs_,shapesface_.nfdofs_);
      Epetra_SerialDenseVector trVec(shapesface_.nfdofs_);

      for(unsigned int q=0; q<shapesface_.nfqpoints_; ++q)
      {
        const double fac = shapesface_.jfacF(q);
        double xyz[nsd_];
        for (unsigned int d=0; d<nsd_; ++d)
          xyz[d] = shapesface_.xyzFreal(d,q);
        double val = 0.0;

       EvaluateFaceAdjoint(fnodexyz,values,numfnode,xyz,val);

        for (unsigned int i=0; i<shapesface_.nfdofs_; ++i)
        {
          // mass matrix
          for (unsigned int j=0; j<shapesface_.nfdofs_; ++j)
            mass(i,j) += shapesface_.shfunctF(i,q) * shapesface_.shfunctF(j,q) * fac;
          trVec(i,0) += shapesface_.shfunctF(i,q) * val * fac * double(numfnode)/double(shapesface_.nfdofs_);
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
        if((shapesface_.shfunctI[face](q,0))>0.00001 ||(shapesface_.shfunctI[face](q,0))<-0.00001)
        {
          sourceterm(q) = trVec(count);
          count++;
        }
      }

    } // else ** if(!smoothing)

    // now, we have to do the condensation!
    Epetra_SerialDenseMatrix toinv((nsd_*nsd_+nsd_+2)*ndofs_,(nsd_*nsd_+nsd_+2)*ndofs_);
    for(unsigned int i=0; i<ndofs_; ++i)
    {
      for(unsigned int j=0; j<ndofs_; ++j)
      {
        toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+j) = imat(j,i);
        toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+ndofs_+j) = kmat(j,i);
        toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+ndofs_+j) = lmat(j,i) + lhatmat(j,i);
        for(unsigned int d=0; d<nsd_; ++d)
        {
          toinv((nsd_*nsd_)*ndofs_+d*ndofs_+i,(nsd_*nsd_)*ndofs_+d*ndofs_+j) = emat(d*ndofs_+j,d*ndofs_+i) + ehatmat(d*ndofs_+j,d*ndofs_+i);
          toinv((nsd_*nsd_)*ndofs_+d*ndofs_+i,(nsd_*nsd_)*ndofs_+nsd_*ndofs_+j) = hmat(j,d*ndofs_+i);
          toinv((nsd_*nsd_)*ndofs_+nsd_*ndofs_+ndofs_+i,(nsd_*nsd_)*ndofs_+d*ndofs_+j) = fmat(d*ndofs_+j,i);
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
    Epetra_SerialDenseVector rhsinv((nsd_*nsd_+nsd_+2)*ndofs_);
    for(unsigned int i=0; i<ndofs_; ++i)
      rhsinv(nsd_*nsd_*ndofs_+nsd_*ndofs_+ndofs_+i) = sourceterm(i,0);

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

    // std::cout<<"sourceterm"<<std::endl;
    // sourceterm.Print(std::cout);
    // std::cout<<"elevec "<<std::endl;
    // elevec.Print(std::cout);
  } // if(step<stepmax)

  return;
}

/*----------------------------------------------------------------------*
 * EvaluateFaceAdjoint
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::LocalSolver::EvaluateFaceAdjoint(
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
}

/*----------------------------------------------------------------------*
 * ComputeSourcePressureMonitor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::LocalSolver::
ComputeSourcePressureMonitorLine3D(DRT::ELEMENTS::Acou*        ele,
                   Teuchos::ParameterList&     params,
                   Teuchos::RCP<MAT::Material> & mat,
                   int                         face,
                   Epetra_SerialDenseMatrix    &elemat,
                   Epetra_SerialDenseVector    &elevec)
{
  dserror("here is something TODO for you!"); // TODO
}



/*----------------------------------------------------------------------*
 * ComputeABCNodeVals
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::
ComputePMonNodeVals(DRT::ELEMENTS::Acou*        ele,
                   Teuchos::ParameterList&     params,
                   Teuchos::RCP<MAT::Material> & mat,
                   int                         face,
                   Epetra_SerialDenseMatrix    &elemat,
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
}

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouViscEleCalc<distype>::
UpdateInteriorVariablesAndComputeResidual(DRT::Discretization &     discretization,
                                          Teuchos::ParameterList&    params,
                                          DRT::ELEMENTS::Acou *                ele,
                                          const Teuchos::RCP<MAT::Material> &mat,
                                          Epetra_SerialDenseVector          & elevec,
                                          double dt,
                                          bool errormaps,
                                          bool updateonly)
{
  DRT::ELEMENTS::AcouVisc * hdgele = dynamic_cast<DRT::ELEMENTS::AcouVisc*>(ele);
  bool adjoint = params.get<bool>("adjoint");
  double theta = 1.0;
  if(dyna_==INPAR::ACOU::acou_trapezoidal) theta = 0.66;

  Epetra_SerialDenseVector tempVelnp;
  Epetra_SerialDenseVector tempPressnp;
  Epetra_SerialDenseVector tempDensnp;
  Epetra_SerialDenseVector tempGradVelnp;
  if(dyna_ == INPAR::ACOU::acou_bdf2)
  {
    hdgele->eleinteriorValn_ = hdgele->eleinteriorValnp_;
    tempVelnp.Shape(shapes_->ndofs_*nsd_,1);                tempVelnp = interiorVelnp_;
    tempPressnp.Shape(shapes_->ndofs_,1);                   tempPressnp  = interiorPressnp_;
    tempDensnp.Shape(shapes_->ndofs_,1);                    tempDensnp  = interiorDensnp_;
    tempGradVelnp.Shape(shapes_->ndofs_*nsd_*nsd_,1);       tempGradVelnp  = interiorGradVelnp_;

    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      interiorVeln_[i] = interiorVelnp_[i] * 4.0 / 3.0 - interiorVelnm_[i] / 3.0;
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorPressn_[i]  = interiorPressnp_[i]  * 4.0 / 3.0 - interiorPressnm_[i]  / 3.0;
      interiorDensn_[i]   = interiorDensnp_[i]   * 4.0 / 3.0 - interiorDensnm_[i]   / 3.0;
    }
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_*nsd_; ++i)
      interiorGradVeln_[i] = interiorGradVelnp_[i] * 4.0 / 3.0 - interiorGradVelnm_[i] / 3.0;
    dt *= 2.0/3.0;
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf3)
  {
    hdgele->eleinteriorValn_ = hdgele->eleinteriorValnp_;
    tempVelnp.Shape(shapes_->ndofs_*nsd_,1);    tempVelnp = interiorVelnp_;
    tempPressnp.Shape(shapes_->ndofs_,1);       tempPressnp  = interiorPressnp_;
    tempDensnp.Shape(shapes_->ndofs_,1);                    tempDensnp  = interiorDensnp_;
    tempGradVelnp.Shape(shapes_->ndofs_*nsd_*nsd_,1);       tempGradVelnp  = interiorGradVelnp_;
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      interiorVeln_[i] = interiorVelnp_[i] * 18.0 / 11.0 - interiorVelnm_[i] * 9.0 / 11.0 + interiorVelnmm_[i] * 2.0 / 11.0;
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorPressn_[i]  = interiorPressnp_[i]  * 18.0 / 11.0 - interiorPressnm_[i]  * 9.0 / 11.0 + interiorPressnmm_[i]  * 2.0 / 11.0;
      interiorDensn_[i]  = interiorDensnp_[i]  * 18.0 / 11.0 - interiorDensnm_[i]  * 9.0 / 11.0 + interiorDensnmm_[i]  * 2.0 / 11.0;
    }
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_*nsd_; ++i)
      interiorGradVeln_[i] = interiorGradVelnp_[i] * 18.0 / 11.0 - interiorGradVelnm_[i] * 9.0 / 11.0 + interiorGradVelnmm_[i] * 2.0 / 11.0;
    dt *= 6.0/11.0;
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf4)
  {
    hdgele->eleinteriorValn_ = hdgele->eleinteriorValnp_;
    tempVelnp.Shape(shapes_->ndofs_*nsd_,1);    tempVelnp = interiorVelnp_;
    tempPressnp.Shape(shapes_->ndofs_,1);       tempPressnp  = interiorPressnp_;
    tempDensnp.Shape(shapes_->ndofs_,1);                    tempDensnp  = interiorDensnp_;
    tempGradVelnp.Shape(shapes_->ndofs_*nsd_*nsd_,1);       tempGradVelnp  = interiorGradVelnp_;
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      interiorVeln_[i] = interiorVelnp_[i] * 48.0 / 25.0 - interiorVelnm_[i] * 36.0 / 25.0 + interiorVelnmm_[i] * 16.0 / 25.0 - interiorVelnmmm_[i] * 3.0 / 25.0;
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorPressn_[i]  = interiorPressnp_[i]  * 48.0 / 25.0 - interiorPressnm_[i]  * 36.0 / 25.0 + interiorPressnmm_[i]  * 16.0 / 25.0 - interiorPressnmmm_[i]  * 3.0 / 25.0;
      interiorDensn_[i]  = interiorDensnp_[i]  * 48.0 / 25.0 - interiorDensnm_[i]  * 36.0 / 25.0 + interiorDensnmm_[i]  * 16.0 / 25.0 - interiorDensnmmm_[i]  * 3.0 / 25.0;
    }
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_*nsd_; ++i)
      interiorGradVeln_[i] =  interiorGradVelnp_[i] * 48.0 / 25.0 - interiorGradVelnm_[i] * 36.0 / 25.0 + interiorGradVelnmm_[i] * 16.0 / 25.0 - interiorGradVelnmmm_[i] * 3.0 / 25.0;
    dt *= 12.0/25.0;
  }
  else if(dyna_ == INPAR::ACOU::acou_dirk23 ||dyna_ == INPAR::ACOU::acou_dirk33 ||dyna_ == INPAR::ACOU::acou_dirk34 ||dyna_ == INPAR::ACOU::acou_dirk54 )
  {
    // do the dirk
    int stage = params.get<int>("stage");
    for(unsigned i=0; i<shapes_->ndofs_; ++i)
    {
      interiorPressn_[i]  = hdgele->eles_[stage][i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1];
      interiorDensn_[i]  = hdgele->eles_[stage][i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_];
      for(unsigned int d=0; d<nsd_; ++d)
        interiorVeln_(i+d*shapes_->ndofs_) = hdgele->eles_[stage][i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+d];
      for(unsigned int d=0; d<nsd_; ++d)
        for(unsigned int e=0; e<nsd_; ++e)
          interiorGradVeln_(i+(e*nsd_+d)*shapes_->ndofs_) = hdgele->eles_[stage][i*(nsd_*nsd_+nsd_+2)+e*nsd_+d];
    }
    interiorPressn_.Scale(dt);
    interiorDensn_.Scale(dt);
    interiorVeln_.Scale(dt);
    interiorGradVeln_.Scale(dt);
  }
  else
  {
    interiorVeln_ = interiorVelnp_;
    interiorPressn_ = interiorPressnp_;
    interiorDensn_ = interiorDensnp_;
    interiorGradVeln_ = interiorGradVelnp_;
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

  Epetra_SerialDenseMatrix toinv;
  if(!adjoint)
  {
    Epetra_SerialDenseMatrix dinvamat(nsd_*shapes_->ndofs_,nsd_*nsd_*shapes_->ndofs_);
    dinvamat.Multiply('N','N',1.0,localSolver_->dmat,localSolver_->invamat,0.0);

    Epetra_SerialDenseMatrix ol(nsd_*shapes_->ndofs_,nsd_*shapes_->ndofs_);
    ol = localSolver_->ehatmat;
    ol.Scale(theta);
    ol += localSolver_->emat;
    ol.Multiply('N','N',-theta,dinvamat,localSolver_->bmat,1.0);

    toinv.Shape((nsd_+2)*shapes_->ndofs_,(nsd_+2)*shapes_->ndofs_);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      for(unsigned int j=0; j<shapes_->ndofs_; ++j)
      {
        toinv(nsd_*shapes_->ndofs_+i,nsd_*shapes_->ndofs_+j) = localSolver_->imat(i,j);
        toinv((nsd_+1)*shapes_->ndofs_+i,(nsd_+1)*shapes_->ndofs_+j) = localSolver_->lmat(i,j) + theta * localSolver_->lhatmat(i,j) ;
        toinv((nsd_+1)*shapes_->ndofs_+i,nsd_*shapes_->ndofs_+j) = theta * localSolver_->kmat(i,j);
        for(unsigned int d=0; d<nsd_; ++d)
        {
          toinv(d*shapes_->ndofs_+i,d*shapes_->ndofs_+j) = ol(d*shapes_->ndofs_+i,d*shapes_->ndofs_+j);
          toinv(nsd_*shapes_->ndofs_+i,d*shapes_->ndofs_+j) = theta * localSolver_->hmat(i,d*shapes_->ndofs_+j);
          toinv(d*shapes_->ndofs_+i,(nsd_+1)*shapes_->ndofs_+j) = theta * localSolver_->fmat(d*shapes_->ndofs_+i,j);
        }
      }
    }

    // calculate rhs parts
    ol.Shape(nsd_*shapes_->ndofs_,1);
    ol.Multiply('N','N',1.0,localSolver_->emat,interiorVeln_,0.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->dmat,interiorGradVeln_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->ehatmat,interiorVeln_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->fmat,interiorPressn_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->gmat,traceVal_SDV_m,1.0);
    ol.Multiply('N','N',-theta,localSolver_->gmat,traceVal_SDV,1.0);
    Epetra_SerialDenseVector tempvec1(nsd_*nsd_*shapes_->ndofs_);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->amat,interiorGradVeln_,0.0);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->bmat,interiorVeln_,1.0);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->cmat,traceVal_SDV_m,1.0);
    tempvec1.Multiply('N','N',-theta,localSolver_->cmat,traceVal_SDV,1.0);
    ol.Multiply('N','N',-1.0,dinvamat,tempvec1,1.0);
    Epetra_SerialDenseVector rhsinv((nsd_+2)*shapes_->ndofs_);
    for(unsigned int i=0; i<nsd_*shapes_->ndofs_; ++i)
      rhsinv(i,0) = ol(i,0);

    ol.Shape(shapes_->ndofs_,1);
    ol.Multiply('N','N',1.0,localSolver_->imat,interiorDensn_,0.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->hmat,interiorVeln_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->jmat,traceVal_SDV_m,1.0);
    ol.Multiply('N','N',-theta,localSolver_->jmat,traceVal_SDV,1.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv(nsd_*shapes_->ndofs_+i,0) = ol(i,0);

    ol.Multiply('N','N',1.0,localSolver_->lmat,interiorPressn_,0.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->kmat,interiorDensn_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->lhatmat,interiorPressn_,1.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv((nsd_+1)*shapes_->ndofs_+i,0) = ol(i,0);

    // invert
    {
      Epetra_SerialDenseSolver inverse;
      inverse.SetMatrix(toinv);
      inverse.Invert();
    }
    Epetra_SerialDenseVector sol((nsd_+2)*shapes_->ndofs_);
    sol.Multiply('N','N',1.0,toinv,rhsinv,0.0);

    for(unsigned int i=0; i<nsd_*shapes_->ndofs_; ++i)
      interiorVelnp_(i,0) = sol(i,0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorDensnp_(i,0)  = sol(nsd_*shapes_->ndofs_+i,0);
      interiorPressnp_(i,0) = sol((nsd_+1)*shapes_->ndofs_+i,0);
    }

    tempvec1.Multiply('N','N',-theta,localSolver_->bmat,interiorVelnp_,1.0);
    interiorGradVelnp_.Multiply('N','N',1.0/theta,localSolver_->invamat,tempvec1,0.0);
  } // if(!adjoint)
  else
  {
    toinv.Shape((nsd_*nsd_+nsd_+2)*shapes_->ndofs_,(nsd_*nsd_+nsd_+2)*shapes_->ndofs_);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      for(unsigned int j=0; j<shapes_->ndofs_; ++j)
      {
        toinv((nsd_*nsd_)*shapes_->ndofs_+nsd_*shapes_->ndofs_+i,(nsd_*nsd_)*shapes_->ndofs_+nsd_*shapes_->ndofs_+j) = localSolver_->imat(j,i);
        toinv((nsd_*nsd_)*shapes_->ndofs_+nsd_*shapes_->ndofs_+i,(nsd_*nsd_)*shapes_->ndofs_+nsd_*shapes_->ndofs_+shapes_->ndofs_+j) = localSolver_->kmat(j,i);
        toinv((nsd_*nsd_)*shapes_->ndofs_+nsd_*shapes_->ndofs_+shapes_->ndofs_+i,(nsd_*nsd_)*shapes_->ndofs_+nsd_*shapes_->ndofs_+shapes_->ndofs_+j) = localSolver_->lmat(j,i) + localSolver_->lhatmat(j,i);
        for(unsigned int d=0; d<nsd_; ++d)
        {
          toinv((nsd_*nsd_)*shapes_->ndofs_+d*shapes_->ndofs_+i,(nsd_*nsd_)*shapes_->ndofs_+d*shapes_->ndofs_+j) = localSolver_->emat(d*shapes_->ndofs_+j,d*shapes_->ndofs_+i) + localSolver_->ehatmat(d*shapes_->ndofs_+j,d*shapes_->ndofs_+i);
          toinv((nsd_*nsd_)*shapes_->ndofs_+d*shapes_->ndofs_+i,(nsd_*nsd_)*shapes_->ndofs_+nsd_*shapes_->ndofs_+j) = localSolver_->hmat(j,d*shapes_->ndofs_+i);
          toinv((nsd_*nsd_)*shapes_->ndofs_+nsd_*shapes_->ndofs_+shapes_->ndofs_+i,(nsd_*nsd_)*shapes_->ndofs_+d*shapes_->ndofs_+j) = localSolver_->fmat(d*shapes_->ndofs_+j,i);
        }
        for(unsigned int e=0; e<nsd_*nsd_; ++e)
        {
          for(unsigned int d=0; d<nsd_*nsd_; ++d)
            toinv(e*shapes_->ndofs_+i,d*shapes_->ndofs_+j) = localSolver_->amat(d*shapes_->ndofs_+j,e*shapes_->ndofs_+i);
          for(unsigned int d=0; d<nsd_; ++d)
          {
            toinv(e*shapes_->ndofs_+i,(nsd_*nsd_)*shapes_->ndofs_+d*shapes_->ndofs_+j) = localSolver_->dmat(d*shapes_->ndofs_+j,e*shapes_->ndofs_+i);
            toinv((nsd_*nsd_)*shapes_->ndofs_+d*shapes_->ndofs_+i,e*shapes_->ndofs_+j) = localSolver_->bmat(e*shapes_->ndofs_+j,d*shapes_->ndofs_+i);
          }
        }
      }
    }

    Epetra_SerialDenseMatrix ol(nsd_*nsd_*shapes_->ndofs_,1);
    ol.Multiply('T','N',-1.0,localSolver_->mmat,traceVal_SDV,0.0);
    Epetra_SerialDenseVector rhsinv((nsd_*nsd_+nsd_+2)*shapes_->ndofs_);
    for(unsigned int i=0; i<nsd_*nsd_*shapes_->ndofs_; ++i)
      rhsinv(i,0) = ol(i,0);

    ol.Shape(nsd_*shapes_->ndofs_,1);
    ol.Multiply('T','N',1.0,localSolver_->emat,interiorVeln_,0.0);
    ol.Multiply('T','N',-1.0,localSolver_->nmat,traceVal_SDV,1.0);
    for(unsigned int i=0; i<nsd_*shapes_->ndofs_; ++i)
      rhsinv(nsd_*nsd_*shapes_->ndofs_+i,0) = ol(i,0);

    ol.Shape(shapes_->ndofs_,1);
    ol.Multiply('T','N',1.0,localSolver_->imat,interiorDensn_,0.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv((nsd_*nsd_+nsd_)*shapes_->ndofs_+i,0) = ol(i,0);

    ol.Multiply('T','N',1.0,localSolver_->lmat,interiorPressn_,0.0);
    ol.Multiply('T','N',-1.0,localSolver_->omat,traceVal_SDV,1.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv((nsd_*nsd_+nsd_+1)*shapes_->ndofs_+i,0) = ol(i,0);

    // invert
    {
      Epetra_SerialDenseSolver inverse;
      inverse.SetMatrix(toinv);
      inverse.Invert();
      double cond = 0.0;
      inverse.ReciprocalConditionEstimate(cond);
    }
    Epetra_SerialDenseVector sol((nsd_*nsd_+nsd_+2)*shapes_->ndofs_);
    sol.Multiply('N','N',1.0,toinv,rhsinv,0.0);

    for(unsigned int i=0; i<nsd_*nsd_*shapes_->ndofs_; ++i)
      interiorGradVelnp_(i,0) = sol(i,0);
    for(unsigned int i=0; i<nsd_*shapes_->ndofs_; ++i)
      interiorVelnp_(i,0) = sol(nsd_*nsd_*shapes_->ndofs_+i,0);

    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorDensnp_(i,0)  = sol(nsd_*nsd_*shapes_->ndofs_+nsd_*shapes_->ndofs_+i,0);
      interiorPressnp_(i,0) = sol(nsd_*nsd_*shapes_->ndofs_+(nsd_+1)*shapes_->ndofs_+i,0);
    }

  } // else ** if(!adjoint)

  // sort them all back to the interior values vector
  for(unsigned int i=0; i<shapes_->ndofs_; ++i)
  {
    interiorValnp_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_+1] = interiorPressnp_(i);
    interiorValnp_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+nsd_] = interiorDensnp_(i);
    for(unsigned int d=0; d<nsd_; ++d)
      interiorValnp_[i*(nsd_*nsd_+nsd_+2)+nsd_*nsd_+d] = interiorVelnp_(i+d*shapes_->ndofs_);
    for(unsigned int d=0; d<nsd_; ++d)
      for(unsigned int e=0; e<nsd_; ++e)
        interiorValnp_[i*(nsd_*nsd_+nsd_+2)+e*nsd_+d] = interiorGradVelnp_(i+(e*nsd_+d)*shapes_->ndofs_);
  }

  if(dyna_ == INPAR::ACOU::acou_dirk23 ||dyna_ == INPAR::ACOU::acou_dirk33 ||dyna_ == INPAR::ACOU::acou_dirk34 ||dyna_ == INPAR::ACOU::acou_dirk54 )
  {
    double dirk_a[6][6];
    double dirk_b[6];
    int dirk_q;
    ACOU::FillDIRKValues(dyna_,dirk_a,dirk_b,dirk_q);

    int stage = params.get<int>("stage");
    for(unsigned i=0; i<interiorValnp_.size(); ++i)
      hdgele->eley_[stage](i) = interiorValnp_[i];
    hdgele->elef_[stage] = hdgele->eleinteriorValn_;
    hdgele->elef_[stage].Scale(-1.0);
    hdgele->elef_[stage] += hdgele->eley_[stage];
    hdgele->elef_[stage].Scale(1.0/dt);  // dt includes a_ii

    Epetra_SerialDenseVector tempVec1(shapes_->ndofs_*(nsd_*nsd_+nsd_+2));
    for(int i=0; i<stage; ++i)
    {
      tempVec1 = hdgele->elef_[i];
      tempVec1.Scale(-dirk_a[stage][i]/dirk_a[stage][stage]);
      hdgele->elef_[stage] += tempVec1;
    }
    tempVec1 = hdgele->elef_[stage];
    tempVec1.Scale(dt*dirk_b[stage]/dirk_a[stage][stage]);
    hdgele->eleinteriorValnp_ += tempVec1;

    if(unsigned(stage)==dirk_q-1)
      hdgele->eleinteriorValn_ = hdgele->eleinteriorValnp_;
  }
  else
  {
    for(unsigned i=0; i<interiorValnp_.size(); ++i)
      hdgele->eleinteriorValnp_[i] = interiorValnp_[i];
  }

  if (updateonly) return;

  // *****************************************************
  // local postprocessing to calculate error maps
  // *****************************************************

  // TODO: local postprocessing

  // *****************************************************
  // compute residual second (reuse intermediate matrices)
  // *****************************************************

  if(dyna_ == INPAR::ACOU::acou_bdf2)
  {
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      interiorVelnp_[i] = 4.0 / 3.0 * interiorVelnp_[i] - 1.0 / 3.0 * tempVelnp[i];
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorPressnp_[i]  = 4.0 / 3.0 * interiorPressnp_[i]  - 1.0 / 3.0 * tempPressnp[i];
      interiorDensnp_[i]  = 4.0 / 3.0 * interiorDensnp_[i]  - 1.0 / 3.0 * tempDensnp[i];
    }
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_*nsd_; ++i)
      interiorGradVelnp_[i] = 4.0 / 3.0 * interiorGradVelnp_[i] - 1.0 / 3.0 * tempGradVelnp[i];
    hdgele->eleinteriorValnm_ = hdgele->eleinteriorValn_;
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf3)
  {
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      interiorVelnp_[i] = 18.0 / 11.0 * interiorVelnp_[i] - 9.0 / 11.0 * tempVelnp[i] + 2.0 / 11.0 * interiorVelnm_[i];
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorPressnp_[i]  = 18.0 / 11.0 * interiorPressnp_[i]  - 9.0 / 11.0 * tempPressnp[i]  + 2.0 / 11.0 * interiorPressnm_[i];
      interiorDensnp_[i]  = 18.0 / 11.0 * interiorDensnp_[i]  - 9.0 / 11.0 * tempDensnp[i]  + 2.0 / 11.0 * interiorDensnm_[i];
    }
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_*nsd_; ++i)
      interiorGradVelnp_[i] =  18.0 / 11.0 * interiorGradVelnp_[i] - 9.0 / 11.0 * tempGradVelnp[i] + 2.0 / 11.0 * interiorGradVelnm_[i];
    hdgele->eleinteriorValnm_ = hdgele->eleinteriorValn_;
    hdgele->eleinteriorValnmm_ = hdgele->eleinteriorValnm_;
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf4)
  {
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      interiorVelnp_[i] = 48.0 / 25.0 * interiorVelnp_[i] - 36.0 / 25.0 * tempVelnp[i] + 16.0 / 25.0 * interiorVelnm_[i] - 3.0 / 25.0 * interiorVelnmm_[i];
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
    {
      interiorPressnp_[i]  = 48.0 / 25.0 * interiorPressnp_[i]  - 36.0 / 25.0 * tempPressnp[i]  + 16.0 / 25.0 * interiorPressnm_[i]  - 3.0 / 25.0 * interiorPressnmm_[i];
      interiorDensnp_[i]   = 48.0 / 25.0 * interiorDensnp_[i]  - 36.0 / 25.0 * tempDensnp[i]  + 16.0 / 25.0 * interiorDensnm_[i]  - 3.0 / 25.0 * interiorDensnmm_[i];
    }
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_*nsd_; ++i)
      interiorGradVelnp_[i] = 48.0 / 25.0 * interiorGradVelnp_[i] - 36.0 / 25.0 * tempGradVelnp[i] + 16.0 / 25.0 * interiorGradVelnm_[i] - 3.0 / 25.0 * interiorGradVelnmm_[i];
    hdgele->eleinteriorValnm_ = hdgele->eleinteriorValn_;
    hdgele->eleinteriorValnmm_ = hdgele->eleinteriorValnm_;
    hdgele->eleinteriorValnmmm_ = hdgele->eleinteriorValnmm_;
  }

  // calculate rhs parts
  if(!adjoint)
  {
    Epetra_SerialDenseMatrix dinvamat(nsd_*shapes_->ndofs_,nsd_*nsd_*shapes_->ndofs_);
    dinvamat.Multiply('N','N',1.0,localSolver_->dmat,localSolver_->invamat,0.0);

    Epetra_SerialDenseMatrix ol(nsd_*shapes_->ndofs_,1);
    ol.Multiply('N','N',1.0,localSolver_->emat,interiorVelnp_,0.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->dmat,interiorGradVelnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->ehatmat,interiorVelnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->fmat,interiorPressnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->gmat,traceVal_SDV,1.0);
    Epetra_SerialDenseVector tempvec1(nsd_*nsd_*shapes_->ndofs_);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->amat,interiorGradVelnp_,0.0);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->bmat,interiorVelnp_,1.0);
    tempvec1.Multiply('N','N',-(1.0-theta),localSolver_->cmat,traceVal_SDV,1.0);
    ol.Multiply('N','N',-1.0,dinvamat,tempvec1,1.0);
    Epetra_SerialDenseVector rhsinv((nsd_+2)*shapes_->ndofs_);
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      rhsinv(i,0) = ol(i,0);

    ol.Shape(shapes_->ndofs_,1);
    ol.Multiply('N','N',1.0,localSolver_->imat,interiorDensnp_,0.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->hmat,interiorVelnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->jmat,traceVal_SDV,1.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv(shapes_->ndofs_*nsd_+i,0) = ol(i,0);

    ol.Multiply('N','N',1.0,localSolver_->lmat,interiorPressnp_,0.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->kmat,interiorDensnp_,1.0);
    ol.Multiply('N','N',-(1.0-theta),localSolver_->lhatmat,interiorPressnp_,1.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv(shapes_->ndofs_*(nsd_+1)+i,0) = ol(i,0);

    Epetra_SerialDenseVector sol((nsd_+2)*shapes_->ndofs_);
    sol.Multiply('N','N',1.0,toinv,rhsinv,0.0);

    elevec.Multiply('N','N',-(1.0-theta),localSolver_->mmat,interiorGradVelnp_,0.0);
    elevec.Multiply('N','N',-(1.0-theta),localSolver_->nmat,interiorVelnp_,1.0);
    elevec.Multiply('N','N',-(1.0-theta),localSolver_->omat,interiorPressnp_,1.0);
    elevec.Multiply('N','N',-(1.0-theta),localSolver_->pmat,traceVal_SDV,1.0);

    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      ol(i,0) = sol(shapes_->ndofs_*(nsd_+1)+i,0);
    elevec.Multiply('N','N',-theta,localSolver_->omat,ol,1.0);

    ol.Shape(shapes_->ndofs_*nsd_,1);
    for(unsigned int i=0; i<shapes_->ndofs_*nsd_; ++i)
      ol(i,0) = sol(i,0);
    elevec.Multiply('N','N',-theta,localSolver_->nmat,ol,1.0);

    tempvec1.Multiply('N','N',-theta,localSolver_->bmat,ol,1.0);
    ol.Shape(nsd_*nsd_*shapes_->ndofs_,1);
    ol.Multiply('N','N',1.0/theta,localSolver_->invamat,tempvec1,0.0);
    elevec.Multiply('N','N',-theta,localSolver_->mmat,ol,1.0);
  } // if(!adjoint)
  else
  {
    Epetra_SerialDenseVector rhsinv((nsd_*nsd_+nsd_+2)*shapes_->ndofs_);
    Epetra_SerialDenseMatrix ol(nsd_*shapes_->ndofs_,1);
    ol.Multiply('N','N',1.0,localSolver_->emat,interiorVelnp_,0.0);
    for(unsigned int i=0; i<nsd_*shapes_->ndofs_; ++i)
      rhsinv(nsd_*nsd_*shapes_->ndofs_+i,0) = ol(i,0);

    ol.Shape(shapes_->ndofs_,1);
    ol.Multiply('N','N',1.0,localSolver_->imat,interiorDensnp_,0.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv(nsd_*nsd_*shapes_->ndofs_+nsd_*shapes_->ndofs_+i,0) = ol(i,0);

    ol.Multiply('N','N',1.0,localSolver_->lmat,interiorPressnp_,0.0);
    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      rhsinv((nsd_*nsd_+nsd_+1)*shapes_->ndofs_+i,0) = ol(i,0);

    Epetra_SerialDenseVector sol((nsd_*nsd_+nsd_+2)*shapes_->ndofs_);
    sol.Multiply('N','N',1.0,toinv,rhsinv,0.0);

    for(unsigned int i=0; i<shapes_->ndofs_; ++i)
      ol(i,0) = sol((nsd_*nsd_+nsd_)*shapes_->ndofs_+i,0);
    elevec.Multiply('T','N',-1.0,localSolver_->jmat,ol,0.0);

    ol.Shape(nsd_*shapes_->ndofs_,1);
    for(unsigned int i=0; i<nsd_*shapes_->ndofs_; ++i)
      ol(i,0) = sol(nsd_*nsd_*shapes_->ndofs_+i,0);
    elevec.Multiply('T','N',-1.0,localSolver_->gmat,ol,1.0);

    ol.Shape(nsd_*nsd_*shapes_->ndofs_,1);
    for(unsigned int i=0; i<nsd_*nsd_*shapes_->ndofs_; ++i)
      ol(i,0) = sol(i,0);
    elevec.Multiply('T','N',-1.0,localSolver_->cmat,ol,1.0);
  } // else ** if(!adjoint)
  return;
}

// template classes
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::tri3>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::AcouViscEleCalc<DRT::Element::nurbs27>;
