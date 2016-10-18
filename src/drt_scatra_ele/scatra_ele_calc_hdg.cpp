/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_hdg.cpp

\brief main file containing routines for calculation of HDG transport element

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_calc_hdg.H"
#include "scatra_ele_calc.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_action.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_geometry/position_array.H"

#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/matlist.H"

#include "../drt_fem_general/drt_utils_polynomial.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include <Epetra_SerialDenseSolver.h>
#include <Teuchos_TimeMonitor.hpp>


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
template <DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::ScaTraEleCalcHDG(const int numdofpernode,
                                                                   const int numscal,
                                                                   const std::string& disname)
  : numdofpernode_(numdofpernode),
    numscal_(numscal),
    usescompletepoly_()

{}


/*----------------------------------------------------------------------*
 | singleton access method                               hoermann 09/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>* DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    bool create
    )
{
  static std::map<std::string,ScaTraEleCalcHDG<distype,probdim>* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcHDG<distype,probdim>(numdofpernode,numscal,disname);
  }

  else if(instances.find(disname) != instances.end())
  {
    for( typename std::map<std::string,ScaTraEleCalcHDG<distype,probdim>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
     {
      delete i->second;
      i->second = NULL;
     }

    instances.clear();
    return NULL;
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                 hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0,0,"",false);

  return;
}


/*----------------------------------------------------------------------*
 | Initialize Shapes                                     hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::InitializeShapes(const DRT::Element*   ele,
                                                                        const std::string& disname
                                                                        )
{

//  DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));
  // Check if this is an HDG element, if yes, can initialize...
  if (DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele)))
  {
    usescompletepoly_ = hdgele->UsesCompletePolynomialSpace();

    if (shapes_ == Teuchos::null)
      shapes_ = Teuchos::rcp(
          new DRT::UTILS::ShapeValues<distype>(hdgele->Degree(), usescompletepoly_,2 * hdgele->Degree()));
    else if (shapes_->degree_ != unsigned(hdgele->Degree()) || shapes_->usescompletepoly_ != usescompletepoly_)
      shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(hdgele->Degree(), usescompletepoly_,2 * hdgele->Degree()));

    int onfdofs = 0;
    for(unsigned int i=0; i<nfaces_; ++i)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele->Faces()[i]->Degree(),
        shapes_->usescompletepoly_,
        2 * ele->Faces()[i]->Degree());

      shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      onfdofs += shapesface_->nfdofs_;
    }
    hdgele->SetDofs(shapes_->ndofs_);
    hdgele->SetOnfDofs(onfdofs);

    // check if only one scalar is defined
    if (numscal_ > 1)
      dserror("Not implemented for multiple scalars");

    if (localSolver_ == Teuchos::null)
      localSolver_ = Teuchos::rcp(new LocalSolver(ele,*shapes_,*shapesface_,usescompletepoly_,disname,1));
  }
  else
    dserror("Only works for HDG transport elements");

}

/*----------------------------------------------------------------------*
 | Evaluate                                              hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::Evaluate(DRT::Element*          ele,
                                                       Teuchos::ParameterList&        params,
                                                       DRT::Discretization&           discretization,
                                                       DRT::Element::LocationArray&   la,
                                                       Epetra_SerialDenseMatrix&      elemat1,
                                                       Epetra_SerialDenseMatrix&      ,
                                                       Epetra_SerialDenseVector&      elevec1,
                                                       Epetra_SerialDenseVector&      ,
                                                       Epetra_SerialDenseVector&      )
{

  // check if this is an hdg element
  const DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<const DRT::ELEMENTS::ScaTraHDG*>(ele);
  if (!hdgele)
    dserror("Cannot cast element to scatra hdg element");

  InitializeShapes(ele, discretization.Name());

  shapes_->Evaluate(*ele);

  ReadGlobalVectors(ele, discretization, la);
  GetMaterialParams(ele);

  zeroMatrix(elevec1);
  if(!localSolver_->scatrapara_->SemiImplicit())
  {
    zeroMatrix(elemat1);
    localSolver_->AddDiffMat(elemat1,hdgele);
    localSolver_->AddReacMat(elemat1,hdgele);
  }
  localSolver_->ComputeResidual(params, elevec1, elemat1, interiorPhin_, tracenm_, tracen_, hdgele);

  return 0;
}


/*----------------------------------------------------------------------*
 * Evaluate Service                                      heormann  09/15|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::EvaluateService(
                                                    DRT::Element*                  ele,
                                                    Teuchos::ParameterList&        params,
                                                    DRT::Discretization&           discretization,
                                                    DRT::Element::LocationArray&   la,
                                                    Epetra_SerialDenseMatrix&      elemat1_epetra,
                                                    Epetra_SerialDenseMatrix&      elemat2_epetra,
                                                    Epetra_SerialDenseVector&      elevec1_epetra,
                                                    Epetra_SerialDenseVector&      elevec2_epetra,
                                                    Epetra_SerialDenseVector&      elevec3_epetra
    )
{
  // check if this is an hdg element
  DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(ele);
  if (!hdgele)
    dserror("cannot cast element to scatrahdg element");

  // get the action required
  const SCATRA::Action act = DRT::INPUT::get<SCATRA::Action>(params,"action");

  InitializeShapes(ele, discretization.Name());

  switch(act)
  {
  case SCATRA::update_interior_variables:
  {
    shapes_->Evaluate(*ele);
    ReadGlobalVectors(ele, discretization, la);
    DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(ele);

    return UpdateInteriorVariables(
        hdgele,
        params,
        elevec1_epetra);
    break;
  }

  case SCATRA::interpolate_hdg_to_node:
  {
    shapes_->Evaluate(*ele);
    ReadGlobalVectors(ele, discretization, la);
    return NodeBasedValues(
        ele,
        discretization,
        elevec1_epetra);
    break;
  }

  case SCATRA::set_initial_field:
  {
    ElementInit(ele);
    // set initial field
    return SetInitialField(
        ele,
        params,
        elevec1_epetra,
        elevec2_epetra);

    break;
  }
  case SCATRA::calc_mat_initial:
  {
    if (hdgele->PadaptEle() || !hdgele->MatInit())
    {
      shapes_->Evaluate(*ele);
      ElementInit(ele);
      ReadGlobalVectors(ele, discretization, la);
      GetMaterialParams(ele);
      localSolver_->ComputeMatrices(ele);
      localSolver_->CondenseLocalPart(hdgele);
    }
    zeroMatrix(elemat1_epetra);
    localSolver_->AddDiffMat(elemat1_epetra,hdgele);

    break;
  }
  case SCATRA::project_field:
  {
    shapes_->Evaluate(*ele);
    return ProjectField(
        ele,
        discretization,
        params,
        elevec1_epetra,
        elevec2_epetra,
        la);
    break;
  }
  case SCATRA::time_update_material:
  {
    TimeUpdateMaterial(ele);
    break;
  }
  case SCATRA::get_material_internal_state:
  {
    GetMaterialInternalState(
        ele,
        params,
        discretization);
    break;
  }
  case SCATRA::set_material_internal_state:
  {
    SetMaterialInternalState(
        ele,
        params,
        discretization);
    break;
  }
  case SCATRA::project_dirich_field:
  {
    return ProjectDirichField(
        ele,
        params,
        discretization,
        la,
        elevec1_epetra);
    break;
  }
  case SCATRA::project_neumann_field:
  {
    int face = params.get<int>("face");
    int sumindex = 0;
    for (int i = 0; i < face; ++i)
    {
      DRT::UTILS::PolynomialSpaceParams parameter(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,ele->Faces()[i]->Degree(), shapes_->usescompletepoly_);
      int nfdofs = DRT::UTILS::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(parameter)->Size();
      sumindex += nfdofs;
    }
    localSolver_->ComputeNeumannBC(
      ele,
      params,
      face,
      elevec1_epetra,
      sumindex);
    break;
  }
  case SCATRA::calc_padaptivity:
  {
    shapes_->Evaluate(*ele);
    ReadGlobalVectors(ele, discretization, la);
    return CalcPAdaptivity(
        ele,
        discretization,
        params);
    break;
  }
  default:
  {
    dserror("Unknown type of action for ScaTraHDG");
    break;
  }
  } // end of switch(act)

  return 0;
}

/*----------------------------------------------------------------------*
 | Calculate node based values                           hoermann 09/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::NodeBasedValues(
    DRT::Element*                        ele,
    DRT::Discretization&                 discretization,
    Epetra_SerialDenseVector&            elevec1)
{
  dsassert(elevec1.M() == (int)nen_*(2+nsd_), "Vector does not have correct size");
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(ele);

  for (unsigned int i=0; i<nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim=0;idim<nsd_;idim++)
      shapes_->xsi(idim) = locations(idim,i);
    shapes_->polySpace_->Evaluate(shapes_->xsi,values);

    // compute values for concentrations (and gradients) by summing over all basis functions
    double sum = 0;
    std::vector<double> sumgrad(nsd_,0.0);
    for (unsigned int k=0; k<hdgele->ndofs_; ++k)
    {
      sum += values(k) * interiorPhinp_(k);
      for (unsigned int d=0; d<nsd_; ++d)
        sumgrad[d] += values(k) * interiorPhinp_(k+(d+1)*hdgele->ndofs_);
    }
    // store node value for concentrations and gradient in element vector
    elevec1(i) = sum;
    for (unsigned int d=0; d<nsd_; ++d)
      elevec1(i+(2+d)*nen_) = sumgrad[d];
  }

  // get trace solution values
  locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace
      (DRT::UTILS::DisTypeToFaceShapeType<distype>::shape);


  Epetra_SerialDenseVector touchcount(nen_);
  Epetra_SerialDenseVector fvalues(1);
  int sumindex = 0;
  for (unsigned int face=0; face<nfaces_; ++face)
  {
    DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele->Faces()[face]->Degree(),
        shapes_->usescompletepoly_,
        2 * ele->Faces()[face]->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, face);

    fvalues.Resize(shapesface_->nfdofs_);

    for (int i=0; i<DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
    {
           // evaluate shape polynomials in node
      for (unsigned int idim=0;idim<nsd_-1;idim++)
        shapesface_->xsi(idim) = locations(idim, i);

      shapesface_->polySpace_->Evaluate(shapesface_->xsi, fvalues);

      double sum = 0.0;
      for (unsigned int k = 0; k < shapesface_->nfdofs_; ++k)
        sum += fvalues(k) * tracen_(sumindex + k);

      elevec1( nen_ + shapesface_->faceNodeOrder[face][i]) += sum;
      touchcount(shapesface_->faceNodeOrder[face][i])++;
    }
    sumindex += shapesface_->nfdofs_;
  }

  for (unsigned int i = 0; i < nen_; ++i)
    elevec1(nen_ + i) /= touchcount(i);

  return 0;
}//NodeBasedValues

/*----------------------------------------------------------------------*
 * ProjectDirichField                                    hoermann 09/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::ProjectDirichField(
    DRT::Element*                  ele,
    Teuchos::ParameterList&        params,
    DRT::Discretization&           discretization,
    DRT::Element::LocationArray&   la,
    Epetra_SerialDenseVector&      elevec1
    )
{

  dserror("Dirichlet boundary conditions not implemented!");

  return 0;
} //ProjectDirichField


/*----------------------------------------------------------------------*
 * ReadGlobalVectors                                     hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::ReadGlobalVectors(
    DRT::Element *                  ele,
    DRT::Discretization &           discretization,
    DRT::Element::LocationArray &   la)
{
  DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));

  // read the HDG solution vector (for traces)
  tracen_.Shape(hdgele->onfdofs_,1);
  interiorPhin_.Shape(shapes_->ndofs_*(nsd_+1),1);
  interiorPhinp_.Shape(shapes_->ndofs_*(nsd_+1),1);
  tracenm_.Shape(hdgele->onfdofs_,1);

  Teuchos::RCP<const Epetra_Vector> phiaf = discretization.GetState("phiaf");
  if(phiaf==Teuchos::null)
    dserror("Cannot get state vector phiaf");
  DRT::UTILS::ExtractMyValues(*phiaf,tracen_,la[0].lm_);


  if (discretization.HasState("phin"))
  {
    Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
    DRT::UTILS::ExtractMyValues(*phin,tracenm_,la[0].lm_);
  }

  Teuchos::RCP<const Epetra_Vector> intphinp = discretization.GetState(2,"intphinp");
  if(intphinp==Teuchos::null)
    dserror("Cannot get state vector intphinp");
  std::vector<int> localDofs = discretization.Dof(2, ele);
  DRT::UTILS::ExtractMyValues(*intphinp,interiorPhinp_,localDofs);

  if (discretization.HasState(2,"intphin"))
  {
    Teuchos::RCP<const Epetra_Vector> intphin = discretization.GetState(2,"intphin");
    DRT::UTILS::ExtractMyValues(*intphin,interiorPhin_,localDofs);
  }

  return;
} // ReadGlobalVectors


/*----------------------------------------------------------------------*
 * LocalSolver                                           hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::
LocalSolver(
    const DRT::Element*                        ele,
    DRT::UTILS::ShapeValues<distype> &         shapeValues,
    DRT::UTILS::ShapeValuesFace<distype> &     shapeValuesFace,
    bool                                       completepoly,
    const std::string&                         disname,
    int                                        numscal
    )
 :  //ndofs_ (shapeValues.ndofs_),
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
    scatraparatimint_(DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(disname),false)
//    diff_(nsd_,nsd_),                        // diffusion coefficient
//    invdiff_(nsd_,nsd_),                     // inverse diffusion coefficient
//    reacoeff_(numscal),                     // reaction coefficient
//    tau_()

{

}//LocalSolver


/*----------------------------------------------------------------------*
 * Compute internal and face matrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::
ComputeMatrices(DRT::Element*   ele
    )
{
  DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));

  // init matrices
  zeroMatrix(hdgele->Amat_);
  zeroMatrix(hdgele->Bmat_);
  zeroMatrix(hdgele->Cmat_);
  zeroMatrix(hdgele->Dmat_);
  zeroMatrix(hdgele->Emat_);
  zeroMatrix(hdgele->Gmat_);
  zeroMatrix(hdgele->Hmat_);
  zeroMatrix(hdgele->Mmat_);
  zeroMatrix(hdgele->EmatT_);
  zeroMatrix(hdgele->BmatMT_);
  zeroMatrix(hdgele->Kmat_);
  zeroMatrix(hdgele->invAMmat_);
  zeroMatrix(hdgele->BTAMmat_);
  zeroMatrix(hdgele->invCondmat_);
  zeroMatrix(hdgele->Xmat_);


  bool usescompletepoly = hdgele->UsesCompletePolynomialSpace();

  shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(hdgele->Degree(),
                                                              usescompletepoly,
                                                              2*ele->Degree()));
  shapes_->Evaluate(*ele);
  ComputeInteriorMatrices(hdgele);

  int sumindex = 0;
  for (unsigned int nface = 0; nface < nfaces_; ++nface)
  {
    DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele->Faces()[nface]->Degree(),
        shapes_->usescompletepoly_,
        2 * ele->Faces()[nface]->Degree());

    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);

    shapesface_->EvaluateFace(*ele, nface);

    ComputeFaceMatrices(nface, sumindex,hdgele);
    sumindex += shapesface_->nfdofs_;
  }

  // calculate AMmat = A + (1/(dt*theta))*M
  double dt = scatraparatimint_->Dt();
  double theta = scatraparatimint_->TimeFac()*(1/dt);

  hdgele->invAMmat_ = hdgele->Mmat_;
  hdgele->invAMmat_.Scale(1.0/(dt*theta));
  hdgele->invAMmat_ += hdgele->Amat_;
  Epetra_SerialDenseSolver inverseAMmat;
  inverseAMmat.SetMatrix(hdgele->invAMmat_);
  int err = inverseAMmat.Invert();
  if (err != 0)
    dserror("Inversion for AMmat failed with errorcode %d", err);

}

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices                                   hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::
ComputeFaceMatrices (const int                    face,
                     int                          indexstart,
                     DRT::ELEMENTS::ScaTraHDG *   hdgele
                     )
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ScaTraEleCalcHDG::ComputeFaceMatrices");

  // Compute the matrices C, E and H
  // Here, we don't consider material properties! Don't forget this during condensation

  // stabilization parameter tau is used as defined in the inputfile, i.e. it is NOT multiplied by diffusion, timestep etc.
  // Thus in the input file one has to define it and already multiply with the diffusion, timstep etc.
  // Two possibilities are to use tau=diff/dx or tau=diff/dx+dx/dt.
  double tau=0.0;

  // set the correct stabilization parameter tau depending on the stabilization method
  switch(scatrapara_->StabType())
  {
  case(INPAR::SCATRA::stabtype_hdg_centered):
    tau = scatrapara_->TauValue();
    break;
  case(INPAR::SCATRA::stabtype_hdg_upwind):
    if (shapesface_->normal(0,0) + shapesface_->normal(1,0) < 0.0)
      tau = 0.0;
    else
      tau = scatrapara_->TauValue();
    break;
  case(INPAR::SCATRA::stabtype_no_stabilization):
    tau = 0.0;
    break;
  default:
    dserror("Unknown definition for stabilization parameter for HDG");
    break;
  }

  // Add convection term (velocity at quadrature points on face)
  // at the moment it is set to zero
  Epetra_SerialDenseMatrix velface(nsd_,shapesface_->nqpoints_);

  // loop over number of shape functions (scalar dofs per element)
  for (unsigned int q = 0; q < hdgele->ndofs_; ++q)
  {
    if( shapesface_->shfunctI.NonzeroOnFace(q) )
    {
      // loop over number of shape functions (scalar dofs per face)
      for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
      {
        // C and E
        // numerical integration: sum over quadrature points
        double tempE = 0.0;
        double tempC = 0.0;
        double temp_d[nsd_] = {0.0};
        // loop over number of quadrature points on face
        for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
        {
          double temp = shapesface_->jfac(i) * shapesface_->shfunct(p, i)
                                             * shapesface_->shfunctI(q, i);
          tempE += temp;
          for (unsigned int j = 0; j < nsd_; ++j)
          {
            velface(j,i)=0.0;
            temp_d[j] += temp * shapesface_->normals(j,i);
            tempC += temp * velface(j,i) * shapesface_->normals(j,i);
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
    } // for (unsigned int q=0; q<hdgele->ndofs_; ++q)
  } // for (unsigned int p=0; p<nfdofs_; ++p)

  // H
  // loop over number of shape functions (scalar dofs per face)
  for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
  {
    // loop over number of shape functions (scalar dofs per face)
    for (unsigned int q = 0; q <shapesface_->nfdofs_; ++q)
    {
      double tempG = 0.0;
      double tempH = 0.0;
      // loop over number of quadrature points on face
      for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
      {
        double temp = shapesface_->jfac(i) * shapesface_->shfunct(p, i) * shapesface_->shfunct(q, i);
        tempG += temp;
        for (unsigned int j = 0; j < nsd_; ++j)
        {
          velface(j,i)=0.0; //convection
          tempH += temp * velface(j,i) * shapesface_->normals(j,i);
        }
      }
      hdgele->Hmat_(indexstart + p, indexstart + q) = tempH - tau * tempG;
    } // for (unsigned int q=0; q<nfdofs_; ++q)
  } // for (unsigned int p=0; p<nfdofs_; ++p)


  // one term is still missing in A
  // loop over number of shape functions (scalar dofs per element)
  for (unsigned int p = 0; p < hdgele->ndofs_; ++p)
  {
    for (unsigned int q = 0; q <= p; ++q)
    {
      double tempA = 0.0;
      if( shapesface_->shfunctI.NonzeroOnFace(p) && shapesface_->shfunctI.NonzeroOnFace(q) )
      {
        // loop over number of quadrature points on face
        for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
          tempA += shapesface_->jfac(i) * shapesface_->shfunctI(p, i) * shapesface_->shfunctI(q, i);
        hdgele->Amat_(p, q) += tau * tempA;
        if (p != q)
          hdgele->Amat_(q, p) += tau * tempA;
      }
    }
  }

  return;
} // ComputeFaceMatrices


/*----------------------------------------------------------------------*
 * ComputeInteriorMatrices                                hoermann 09/15|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::
ComputeInteriorMatrices(DRT::ELEMENTS::ScaTraHDG * hdgele)
{

  Epetra_SerialDenseMatrix vel(nsd_,shapes_->nqpoints_);
  Epetra_SerialDenseMatrix gradPart(hdgele->ndofs_*nsd_, shapes_->nqpoints_);
  Epetra_SerialDenseMatrix gradPartVel(hdgele->ndofs_,shapes_->nqpoints_);

  Epetra_SerialDenseMatrix  massPart(hdgele->ndofs_,shapes_->nqpoints_);
  Epetra_SerialDenseMatrix  massPartW(hdgele->ndofs_,shapes_->nqpoints_);

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
        vel(d,q)=0.0;
        gradPart(d * hdgele->ndofs_ + i, q) = shapes_->shderxy(i * nsd_ + d, q);

        gradPartVel(i,q) += shapes_->shderxy(i * nsd_ + d, q) * vel(d,q);
      }
    }
  }

  // multiply matrices to perform summation over quadrature points
  hdgele->Mmat_.Multiply('N', 'T', 1.0 , massPart, massPartW, 0.0);
  hdgele->Amat_.Multiply('N', 'T', -1.0 , gradPartVel, massPartW, 0.0); // first part of A matrix (only if velocity field not zero)
  hdgele->Bmat_.Multiply('N', 'T', -1.0, massPartW, gradPart, 0.0);

  for (unsigned int j = 0; j < hdgele->ndofs_; ++j)
    for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
    {
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        for (unsigned int e = 0; e<nsd_; e++)
          hdgele->Dmat_(d * hdgele->ndofs_ +i, e * hdgele->ndofs_ + j) = hdgele->invdiff_(d,e) * hdgele->Mmat_(i,j);
        hdgele->BmatMT_(d * hdgele->ndofs_ + i, j) = -1.0 * hdgele->Bmat_(j, d * hdgele->ndofs_ + i);
      }
    }

  return;
}//ComputeInteriorMatrices

/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::ComputeResidual(
    Teuchos::ParameterList&           params,
    Epetra_SerialDenseVector &        elevec,
    Epetra_SerialDenseMatrix &        elemat1,
    Epetra_SerialDenseVector          interiorPhin,
    Epetra_SerialDenseVector          tracen,
    Epetra_SerialDenseVector          tracenp,
    const DRT::ELEMENTS::ScaTraHDG *  hdgele
    )
{

  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ScaTraEleCalcHDG::ComputeResidual");

  /*
   for trapezoidal rule
                                                   -1
                        +--------------------------+   +----------------------------------------------------------------+
                        |                          |   |      n                          n       n            n         |  }=:s
    n+1                 |  (1/dt*theta)M + A    B  |   | M Phi   + dt (1 - theta) ( A Phi   + B Q   + C Lambda   )      |                                  n          n            n
   R     = - [ G  E^T ] |                          |   |                                     n       n            n     |          + dt (1 - theta) ( G Phi    + E^T Q  + H Lambda   )  -
                        |          -B^T         D  |   |   0     + dt (1 - theta) ( - B^T Phi   + D Q   + E Lambda   )  |  }=:t
                        +--------------------------+   +----------------------------------------------------------------+
                                                              -1
                                   +-------------------------+    +---------------------------------------------+
                                   |                         |    |                    n                  n+1   |
                                   | (1/dt*theta)M + A    B  |    |  - dt (1 - theta) I    -   dt theta  I      |
                      - [ G  E^T ] |                         |    |                                             |
                                   |      -B^T            D  |    |                 0                           |
                                   +-------------------------+    +---------------------------------------------+
   I..reactive term


    With

     s =   -M * U^n + dt*(1-theta) * (  A U^n + B Q^n + C L^n )   - dt*theta I^n+1 -dt*(1-theta) I^n

     t =             dt*(1-theta) * (-B^T U^n  + D Q^n + E L^n )



                              +--------------------------+^-1   +-----+
    n+1                       |  1/(dt*theta)M + A    B  |      |  s  |                           n          n            n
   R     =  - [ G  E^T ]      |                          |      |     |     +  dt(1-theta) ( G Phi    + E^T Q  + H Lambda )  =
                              |        -B^T           D  |      |  t  |
                              +--------------------------+      +-----+

                                                                    n          n      n
                         =  - (G x + E^T y ) +   dt(1-theta) (  G  U    + E^T Q  + H L  )

               with
                                            -1    -1                                    -1
       y = ( D - (-B^T)   (1/(dt*theta)M + A)    B)     (t - (-B^T)   (1/(dt*theta)M + A)    s)

       x = (1/(dt*theta)M + A)^-1 ( s - B y)


       AM = (1/(dt*theta)M + A)
   */

  Epetra_SerialDenseVector tempinteriorphin(hdgele->ndofs_);
  for (unsigned int i = 0; i < hdgele->ndofs_; i++)
    tempinteriorphin(i) = interiorPhin(i);

  Epetra_SerialDenseVector tempinteriorgradphin(hdgele->ndofs_*nsd_);
  for (unsigned int i = 0; i < hdgele->ndofs_*nsd_; i++)
    tempinteriorgradphin(i) = interiorPhin(hdgele->ndofs_+i);

  int onfdofs = elevec.M();
  Epetra_SerialDenseVector trace_SDV(onfdofs);
  for (int i = 0; i < onfdofs; ++i)
    trace_SDV(i) = tracen(i);

  Epetra_SerialDenseVector tracenp_SDV(onfdofs);
  for (int i = 0; i < onfdofs; ++i)
    tracenp_SDV(i) = tracenp(i);

  double dt = scatraparatimint_->Dt();
  double theta = scatraparatimint_->TimeFac()*(1/dt);

  Epetra_SerialDenseVector tempVec1(hdgele->ndofs_);
  Epetra_SerialDenseVector tempVec2(hdgele->ndofs_*nsd_);

  if (theta != 1.0)
  {
    tempVec1.Multiply('N', 'N', 1.0, hdgele->Amat_, tempinteriorphin, 0.0);
    tempVec1.Multiply('N', 'N', 1.0, hdgele->Bmat_, tempinteriorgradphin, 1.0);
    tempVec1.Multiply('N', 'N', 1.0, hdgele->Cmat_, trace_SDV, 1.0); // = (  A U^n + B Q^n + C L^n )
    tempVec1.Scale(dt*(1.0-theta));
  }
  tempVec1.Multiply('N', 'N', -1.0, hdgele->Mmat_, tempinteriorphin, 1.0); // = s = -M * U^n + dt*(1-theta) * (  A U^n + B Q^n + C L^n )

  Epetra_SerialDenseVector tempVecI(hdgele->ndofs_);
  if(!scatrapara_->SemiImplicit())
  {
    // Reaction term
    tempVecI = hdgele->Ivecnp_;
    tempVecI.Scale(dt*theta);
    tempVec1 += tempVecI;
    if (theta != 1.0)
    {
      tempVecI = hdgele->Ivecn_;
      tempVecI.Scale(dt*(1.0-theta));
      tempVec1 += tempVecI;
    }
  }
  //= M * U^n + dt*(1-theta) * (  A U^n + B Q^n + C L^n )   - dt*theta I^n+1 -dt*(1-theta) I^n
  else
  {
    tempVecI = hdgele->Ivecn_;
    tempVecI.Scale(dt);
    tempVec1 += tempVecI;
  }

  if (theta != 1.0)
  {
    tempVec2.Multiply('N', 'N', 1.0, hdgele->BmatMT_, tempinteriorphin, 0.0);
    tempVec2.Multiply('N', 'N', 1.0, hdgele->Dmat_, tempinteriorgradphin, 1.0);
    tempVec2.Multiply('N', 'N', 1.0, hdgele->Emat_, trace_SDV, 1.0); // = (-B^T U^n  + D Q^n + E L^n )
    tempVec2.Scale(dt*(1.0-theta));//= t = dt*(1-theta) * (-B^T U^n  + D Q^n + E L^n )
  }

  tempVec2.Multiply('N', 'N', -1.0, hdgele->BTAMmat_, tempVec1, 1.0); // = t - (-B^T) AM^{-1} s

  Epetra_SerialDenseVector tempVec3(hdgele->ndofs_*nsd_);
  tempVec3.Multiply('N', 'N', 1.0, hdgele->invCondmat_, tempVec2, 0.0); // = y= ( D - (-B^T)   (AM)^-1    B)^-1     (t - (-B^T)   (AM^-1)    s)

  tempVec1.Multiply('N', 'N', -1.0, hdgele->Bmat_, tempVec3, 1.0); // = ( s - B y)
  Epetra_SerialDenseVector tempVec4(hdgele->ndofs_);
  tempVec4.Multiply('N', 'N', 1.0, hdgele->invAMmat_, tempVec1, 0.0); // = x= (1/(dt*theta)M + A)^-1 ( s - B y)


  if (theta != 1.0)
  {
    elevec.Multiply('N', 'N', 1.0, hdgele->Gmat_, tempinteriorphin, 0.0);
    elevec.Multiply('N', 'N', 1.0, hdgele->EmatT_, tempinteriorgradphin, 1.0);
    elevec.Multiply('N', 'N', 1.0, hdgele->Hmat_, trace_SDV, 1.0); // = (  G  U^n    + E^T Q^n  + H L^n  )
    elevec.Scale(dt*(1.0-theta));// = dt*(1-theta) * (  G  U^n    + E^T Q^n  + H L^n  )
  }

  elevec.Multiply('N', 'N', -1.0, hdgele->Gmat_, tempVec4, 1.0);
  elevec.Multiply('N', 'N', -1.0, hdgele->EmatT_, tempVec3, 1.0); // =  - (G x + E^T y ) +   dt(1-theta) (  G  U    + E^T Q  + H L  )

  elevec.Multiply('N', 'N', 1.0, hdgele->Kmat_, tracenp_SDV, 1.0);

  return;
} // ComputeResidual


/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::
CondenseLocalPart(DRT::ELEMENTS::ScaTraHDG * hdgele)
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
    K = (dt*theta)   ( - [ G  E^T ]  |          |     |     |  + H )   = (dt*theta) ( - G x - E^T y + H )  ,
                     (               | -B^T  D  |     |  E  |      )
                     (               +----------+     +-----+      )

       with
                            -1   -1                    -1
       y = ( D - (-B^T)   AM    B)     (E - (-B^T)   AM    C)

       x = AM^-1 ( C - B y)


       AM = (1/(dt*theta)M + A)

   */

  int onfdofs = hdgele->Kmat_.M();

  double dt = scatraparatimint_->Dt();
  double theta = scatraparatimint_->TimeFac()*(1/dt);

  Epetra_SerialDenseMatrix tempMat1(hdgele->ndofs_*nsd_, hdgele->ndofs_);
  tempMat1.Multiply('N', 'N', 1.0, hdgele->BmatMT_, hdgele->invAMmat_, 0.0); // =  (-B^T) AM^{-1}

  hdgele->BTAMmat_ = tempMat1;

  Epetra_SerialDenseMatrix tempMat2(hdgele->ndofs_*nsd_, hdgele->ndofs_*nsd_);
  tempMat2 = hdgele->Dmat_;

  tempMat2.Multiply('N', 'N', -1.0, tempMat1, hdgele->Bmat_, 1.0); // = D - (-B^T) AM^{-1} B
  Epetra_SerialDenseMatrix tempMat3(hdgele->ndofs_*nsd_,onfdofs);
  tempMat3 = hdgele->Emat_;
  tempMat3.Multiply('N', 'N', -1.0, tempMat1, hdgele->Cmat_, 1.0); // = E - (-B^T) AM^{-1} C

  Epetra_SerialDenseSolver inverseinW;
  inverseinW.SetMatrix(tempMat2);
  int err = inverseinW.Invert();
  if (err != 0)
    dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d", err);
  // tempMat2 = (  D - H A^{-1} B )^{-1}

  hdgele->invCondmat_ = tempMat2;

  hdgele->Kmat_ = hdgele->Hmat_; // = H

  Epetra_SerialDenseMatrix tempMat4(hdgele->ndofs_*nsd_, onfdofs);
  tempMat4.Multiply('N', 'N', 1.0, tempMat2, tempMat3, 0.0); // = y
  hdgele->Kmat_.Multiply('N', 'N', -1.0, hdgele->EmatT_, tempMat4, 1.0); // = - E^T y + H

  Epetra_SerialDenseMatrix tempMat5(hdgele->ndofs_, onfdofs);
  tempMat5 = hdgele->Cmat_;
  tempMat5.Multiply('N', 'N', -1.0, hdgele->Bmat_, tempMat4, 1.0); // = C -B y

  Epetra_SerialDenseMatrix tempMat6(hdgele->ndofs_, onfdofs);
  tempMat6.Multiply('N', 'N', 1.0, hdgele->invAMmat_, tempMat5, 0.0); // = x = AM^{-1} ( C - B y )

  // save for later use
  hdgele->Xmat_ = tempMat6;

  hdgele->Kmat_.Multiply('N', 'N', -1.0, hdgele->Gmat_, tempMat6, 1.0); // = K = H - G x - E^T y

  hdgele->Kmat_.Scale(dt*theta);

  return;
} // CondenseLocalPart

/*----------------------------------------------------------------------*
 * Add Diffusive Part to Matrix
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::
AddDiffMat(Epetra_SerialDenseMatrix &eleMat,
    const DRT::ELEMENTS::ScaTraHDG * hdgele)
{
  eleMat = hdgele->Kmat_;
  eleMat.Scale(-1.0);

  return;
}//AddDiffMat


/*----------------------------------------------------------------------*
 * Add Reactive Part to Matrix
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::
AddReacMat(Epetra_SerialDenseMatrix &eleMat,
    const DRT::ELEMENTS::ScaTraHDG * hdgele)
{
  double dt = scatraparatimint_->Dt();
  double theta = scatraparatimint_->TimeFac()*(1/dt);

  // add derivative of reaction term
  Epetra_SerialDenseMatrix tempMat1(hdgele->ndofs_,hdgele->ndofs_);
  tempMat1 = hdgele->Imatnpderiv_; // = I'

  Epetra_SerialDenseMatrix tempMat2(hdgele->ndofs_, hdgele->onfdofs_);
  tempMat2.Multiply('N', 'N', -1.0, tempMat1, hdgele->Xmat_, 0.0); // = I' * (-x1)

  Epetra_SerialDenseMatrix tempMat3(hdgele->ndofs_*nsd_,hdgele->onfdofs_);
  tempMat3.Multiply('N', 'N', -1.0, hdgele->BTAMmat_, tempMat2, 0.0); // = NULL*y1 - (-B^T) AM^{-1} I'* (-x1)
  Epetra_SerialDenseMatrix tempMat4(hdgele->ndofs_*nsd_, hdgele->onfdofs_);
  tempMat4.Multiply('N', 'N', 1.0, hdgele->invCondmat_, tempMat3, 0.0); // = y2 = ( D - (-B^T) AM^{-1} B)^-1  (NULL*y1 - (-B^T) AM^{-1} I'* (-x1))

  tempMat2.Multiply('N', 'N', -1.0, hdgele->Bmat_, tempMat4, 1.0); // = I'*(-x1) -B y2

  Epetra_SerialDenseMatrix tempMat5(hdgele->ndofs_, hdgele->onfdofs_);
  tempMat5.Multiply('N', 'N', 1.0, hdgele->invAMmat_, tempMat2, 0.0); // = x2 = AM^{-1} ( I'*(-x1) - B y2 )

  eleMat.Multiply('N', 'N', dt*theta, hdgele->EmatT_, tempMat4, 1.0); // = K - E^T y2
  eleMat.Multiply('N', 'N', dt*theta, hdgele->Gmat_, tempMat5, 1.0); // = K- G x2 - E^T y2

}//AddReacMat


/*----------------------------------------------------------------------*
 |  Compute Neumann BC                                    hoermann 09/15|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::
ComputeNeumannBC(DRT::Element*               ele,
                 Teuchos::ParameterList&     params,
                 int                         face,
                 Epetra_SerialDenseVector&   elevec,
                 int                         indexstart
                 )
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::ScaTraHDGEleCalc::ComputeNeumannBC");

  DRT::Condition* condition = params.get<DRT::Condition*>("condition");
  if(condition == NULL)
    dserror("Cannot access Neumann boundary condition!");

  // find out whether we will use a time curve
  bool usetime = true;

  // get actual time
  const double time = scatraparatimint_->Time();
  if (time<0.0) usetime = false;

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const std::vector<int>*    onoff  = condition->Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val    = condition->Get<std::vector<double> >("val"  );
  const std::vector<int>*    curve  = condition->Get<std::vector<int> >   ("curve");
  const std::vector<int>*    func   = condition->Get<std::vector<int> >   ("funct");


  DRT::UTILS::ShapeValuesFaceParams svfparams(
      ele->Faces()[face]->Degree(),
      shapes_->usescompletepoly_,
      2 * ele->Faces()[face]->Degree());

  shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);

  shapesface_->EvaluateFace(*ele, face);
  shapes_->Evaluate(*ele);

  // integration loop
  for (unsigned int iquad=0; iquad<shapesface_->nqpoints_; ++iquad)
  {
    // factor given by spatial function
    double functfac = 1.0;
    // factor given by temporal curve
    double curvefac = 1.0;

    // global coordinates of current Gauss point
    double coordgp[3];   // we always need three coordinates for function evaluation!
    for(int i=0; i<3; ++i)
      coordgp[i] = shapesface_->xyzreal(i,iquad);

    int functnum = -1;
    int curvenum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation

    if ((*onoff)[0]) // is this dof activated?
    {

      // find out whether we will use a time curve and get the factor
      if (curve)
        curvenum = (*curve)[0];

      if (curvenum>=0 && usetime)
      {
        // evaluate curve at current time
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
      }
      else
        curvefac = 1.0;

      // factor given by spatial function
      if(func)
        functnum = (*func)[0];

      if(functnum>0)
      {
        // evaluate function at current Gauss point (provide always 3D coordinates!)
        functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(0,coordgpref,time,NULL);
      }
      else
        functfac = 1.;

      const double val_fac_funct_curve_fac = (*val)[0]*shapesface_->jfac(iquad)*functfac*curvefac;

      for(unsigned int node=0; node<shapesface_->nfdofs_; node++)
        elevec[face*shapesface_->nfdofs_+node] += shapesface_->shfunct(node,iquad) * val_fac_funct_curve_fac;
    } // if ((*onoff)[dof])
  } // loop over integration points


  return;
}//ComputeNeumannBC


/*----------------------------------------------------------------------*
 |  get the material parameter                            hoermann 09/15|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::GetMaterialParams(
  DRT::Element*   ele   //!< the element we are dealing with
  )
{

  Epetra_SerialDenseMatrix difftensor(nsd_,nsd_);
  Epetra_SerialDenseVector ivecn(shapes_->ndofs_);
  Epetra_SerialDenseVector ivecnp(shapes_->ndofs_);
  Epetra_SerialDenseMatrix ivecnpderiv(shapes_->ndofs_,shapes_->ndofs_);

  double diff1 = 0;
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() < numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      Materials(singlemat,k,&difftensor,&ivecn,&ivecnp,&ivecnpderiv,&diff1);
    }
  }
  else
    Materials(material,0,&difftensor,&ivecn,&ivecnp,&ivecnpderiv,&diff1);

  DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(ele);
  localSolver_->SetMaterialParameter(hdgele,&difftensor,&ivecn,&ivecnp,&ivecnpderiv,&diff1);



  return;
} //ScaTraEleCalcHDG::GetMaterialParams


/*----------------------------------------------------------------------*
 * UpdateInteriorVariables
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::UpdateInteriorVariables(
    DRT::ELEMENTS::ScaTraHDG *   hdgele,
    Teuchos::ParameterList&      params,
    Epetra_SerialDenseVector &   elevec
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


    +---------+                       +--------------------------+^-1  +--------------------------------------------------------------------------------------------+
    |  U^n+1  |                       |  1/(dt*theta)M + A    B  |     |   M * U^n - dt*(1-theta) * (  A U^n + B Q^n + C L^n )   - dt*theta I^n+1 -dt*(1-theta) I^n |
    |         |    = 1/(dt * theta)   |                          |     |                                                                                            | -
    |  Q^n+1  |                       |        -B^T           D  |     |  - dt*(1-theta) * (   -B^T U^n  + D Q^n + E L^n )                                          |
    +---------+                       +--------------------------+     +--------------------------------------------------------------------------------------------+

                              +--------------------------+^-1    +------+
                              |  1/(dt*theta)M + A    B  |       |   C  |
                        -     |                          |       |      |  L^n+1
                              |        -B^T           D  |       |   E  |
                              +--------------------------+       +------+


  s =   M * U^n - dt*(1-theta) * (  A U^n + B Q^n + C L^n )   - dt*theta I^n+1 -dt*(1-theta) I^n     -     dt* theta C L^n+1

  t =            - dt*(1-theta) * (-B^T U^n  + D Q^n + E L^n )                                        -     dt*theta E L^n+1



   +---------+                        +--------------------------+^-1   +-----+                          +-----+
   |  U^n+1  |                        |  1/(dt*theta)M + A    B  |      |  s  |                          |  x  |
   |         |   =   1/(dt * theta)   |                          |      |     |      =   1/(dt * theta)  |     |    ,
   |  Q^n+1  |                        |        -B^T           D  |      |  t  |                          |  y  |
   +---------+                        +--------------------------+      +-----+                          +-----+

       with
                                             -1    -1                                    -1
       y = ( D - (-B^T)   (1/(dt*theta)M + A)    B)     (t - (-B^T)   (1/(dt*theta)M + A)    s)

       x = (1/(dt*theta)M + A)^-1 ( s - B y)

       AM = (1/(dt*theta)M + A)

   */

  Epetra_SerialDenseVector tempinteriorphin(hdgele->ndofs_);
  for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
    tempinteriorphin(i) = interiorPhin_(i);

  Epetra_SerialDenseVector tempinteriorgradphin(hdgele->ndofs_*nsd_);
  for (unsigned int i = 0; i < hdgele->ndofs_*nsd_; ++i)
    tempinteriorgradphin(i) = interiorPhin_(hdgele->ndofs_+i);

  double dt = localSolver_->scatraparatimint_->Dt();
  double theta = localSolver_->scatraparatimint_->TimeFac()*(1/dt);

  Epetra_SerialDenseVector tempVec1(hdgele->ndofs_);
  if (theta != 1.0)
  {
    tempVec1.Multiply('N', 'N', 1.0, hdgele->Amat_, tempinteriorphin, 0.0);
    tempVec1.Multiply('N', 'N', 1.0, hdgele->Bmat_, tempinteriorgradphin, 1.0);
    tempVec1.Multiply('N', 'N', 1.0, hdgele->Cmat_, tracenm_, 1.0); // = ( A U^n + B Q^n + C L^n )
    tempVec1.Scale(-dt*(1.-theta));
  }
  tempVec1.Multiply('N', 'N', 1.0, hdgele->Mmat_, tempinteriorphin, 1.0); // = -M * U^n + dt*(1-theta) * (  A U^n + B Q^n + C L^n )


  // Reaction term
  Epetra_SerialDenseVector tempVecI(hdgele->ndofs_);
  if(!localSolver_->scatrapara_->SemiImplicit())
  {
    tempVecI = hdgele->Ivecnp_;
    tempVecI.Scale(-dt*theta);
    tempVec1 += tempVecI;
    if (theta != 1.0)
    {
      tempVecI = hdgele->Ivecn_;
      tempVecI.Scale(-dt*(1.0-theta));
      tempVec1 += tempVecI;
    }
  }
  else
  {
    tempVecI = hdgele->Ivecn_;
    tempVecI.Scale(-dt);
    tempVec1 += tempVecI;
  }


  tempVec1.Multiply('N', 'N', -dt*theta, hdgele->Cmat_, tracen_, 1.0);//= s = -M * U^n + dt*(1-theta) * (A U^n + B Q^n + C L^n) - dt*theta I^n+1 -dt*(1-theta) I^n - dt* theta C L^n+1

  Epetra_SerialDenseVector tempVec2(hdgele->ndofs_*nsd_);
  if (theta != 1.0)
  {
    tempVec2.Multiply('N', 'N', 1.0, hdgele->BmatMT_, tempinteriorphin, 0.0);
    tempVec2.Multiply('N', 'N', 1.0, hdgele->Dmat_, tempinteriorgradphin, 1.0);
    tempVec2.Multiply('N', 'N', 1.0, hdgele->Emat_, tracenm_, 1.0); // = (-B^T U^n  + D Q^n + E L^n )
    tempVec2.Scale(-dt*(1.-theta));
  }
  tempVec2.Multiply('N', 'N', -dt*theta, hdgele->Emat_, tracen_, 1.0); //= t = -dt*(1-theta) * (-B^T U^n  + D Q^n + E L^n )  - dt*theta E L^n+1

  // y = ( D - (-B^T)   (AM)^-1    B)^-1     (t - (-B^T)   (AM^-1)    s)
  //x = (1/(dt*theta)M + A)^-1 ( s - B y)


  tempVec2.Multiply('N', 'N', -1.0, hdgele->BTAMmat_, tempVec1, 1.0); // = t - (-B^T) AM^{-1} s

  Epetra_SerialDenseVector tempVec3(hdgele->ndofs_*nsd_);
  tempVec3.Multiply('N', 'N', 1.0, hdgele->invCondmat_, tempVec2, 0.0); // = y= ( D - (-B^T)   (AM)^-1    B)^-1     (t - (-B^T)   (AM^-1)    s)

  tempVec1.Multiply('N', 'N', -1.0, hdgele->Bmat_, tempVec3, 1.0); // = ( s - B y)
  Epetra_SerialDenseVector tempVec4(hdgele->ndofs_);
  tempVec4.Multiply('N', 'N', 1.0, hdgele->invAMmat_, tempVec1, 0.0); // = x= (1/(dt*theta)M + A)^-1 ( s - B y)

  tempVec3.Scale(1./(dt*theta));
  tempVec4.Scale(1./(dt*theta));

  for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
    elevec(i) = tempVec4(i);
  for (unsigned int i = 0; i <nsd_*hdgele->ndofs_; ++i)
    elevec(hdgele->ndofs_+i) = tempVec3(i);


  return 0;
} // UpdateInteriorVariables


/*----------------------------------------------------------------------*
 * SetInitialField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::SetInitialField(
    const DRT::Element*                  ele,
    Teuchos::ParameterList&              params,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2
)
{
  shapes_->Evaluate(*ele);

  Epetra_SerialDenseMatrix Mmat(shapes_->ndofs_,shapes_->ndofs_);
  Epetra_SerialDenseMatrix massPart(shapes_->ndofs_,shapes_->nqpoints_);
  Epetra_SerialDenseMatrix  massPartW(shapes_->ndofs_,shapes_->nqpoints_);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 ||
      unsigned(elevec2.M()) == shapes_->ndofs_*(nsd_+1), "Wrong size in project vector 2");

  // get function
  const int *start_func = params.getPtr<int>("funct");

  // internal variables
  if (elevec2.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(View, elevec2.A(), shapes_->ndofs_, shapes_->ndofs_, nsd_+1, false);
    zeroMatrix(localMat);

    //create mass matrix for interior by looping over quadrature points
    for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
    {
      const double fac = shapes_->jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_->xyzreal(d, q); // coordinates of quadrature point in real coordinates
      double phi;
      double gradphi[nsd_];

      dsassert(start_func != NULL,"funct not set for initial value");
      phi = DRT::Problem::Instance()->Funct(*start_func-1).Evaluate(0,xyz,0,NULL);
      for (unsigned int i = 0; i < nsd_; ++i)
        gradphi[i] = DRT::Problem::Instance()->Funct(*start_func-1).Evaluate(1+i,xyz,0,NULL);

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      {
        massPart(i, q) = shapes_->shfunct(i, q);
        massPartW(i, q) = shapes_->shfunct(i, q) * fac;
        localMat(i,0) += shapes_->shfunct(i, q) * phi * fac;
        for(unsigned int j=0; j<nsd_; ++j)
          localMat(i,1+j) += shapes_->shfunct(i,q) * gradphi[j] * fac;
      }
    }

    Mmat.Multiply('N', 'T', 1., massPart, massPartW, 0.);
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(Mmat);
      inverseMass.SetVectors(localMat, localMat);
      inverseMass.Solve();
    }
  }

// We have to set the initial trace field and gradient field,
// because the calculation of the residual in the first time step should have a correct trace field
// and also a correct gradient!

// trace variable
  int nfdofs = 0;
  for (unsigned int face=0; face<nfaces_; ++face)
  {

    DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele->Faces()[face]->Degree(),
        shapes_->usescompletepoly_,
        2 * ele->Faces()[face]->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, face);

    Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_,shapesface_->nfdofs_);
    Epetra_SerialDenseMatrix trVec(shapesface_->nfdofs_,1);

    // loop over quadrature points
    for (unsigned int q=0; q<shapesface_->nqpoints_; ++q)
    {
      const double fac = shapesface_->jfac(q);
      double xyz[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        xyz[d] = shapesface_->xyzreal(d,q);

      double trphi;
      trphi = DRT::Problem::Instance()->Funct(*start_func-1).Evaluate(nsd_+1,xyz,0,NULL);

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i=0; i<shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j=0; j<shapesface_->nfdofs_; ++j)
          mass(i,j) += shapesface_->shfunct(i,q) * shapesface_->shfunct(j,q) * fac;
        trVec(i,0) += shapesface_->shfunct(i,q) * trphi * fac;
      }
    }

    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec,trVec);
    inverseMass.Solve();
    for (unsigned int i=0; i<shapesface_->nfdofs_; ++i)
      elevec1(nfdofs+i) = trVec(i,0);

    nfdofs += shapesface_->nfdofs_;
  }

  return 0;

}//SetInitialField


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::Materials(
    const Teuchos::RCP<const MAT::Material>   material,      //!< pointer to current material
    const int                                 k,             //!< id of current scalar
    Epetra_SerialDenseMatrix*                 difftensor,    //!< diffusion tensor
    Epetra_SerialDenseVector*                 ivecn,         //!< reaction term at time n
    Epetra_SerialDenseVector*                 ivecnp,        //!< reaction term at time n+1
    Epetra_SerialDenseMatrix*                 ivecnpderiv,   //!< reaction term derivative
    double*                                   diff1          //!< main diffusivity
    )
{
  const Teuchos::RCP<const MAT::ScatraMat>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

  double diffscalar;

  // get constant diffusivity
  diffscalar = actmat->Diffusivity();

  for (unsigned int i=0; i<nsd_;++i)
    (*difftensor)(i,i) = diffscalar;

  *diff1 = diffscalar;

  return;
} //ScaTraEleCalcHDG::Materials


/*----------------------------------------------------------------------*
 |  Set material parameter                               hoermann 09/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::LocalSolver::SetMaterialParameter(
    DRT::ELEMENTS::ScaTraHDG *   hdgele,        //!< hdg element
    Epetra_SerialDenseMatrix*    difftensor,    //!< diffusion tensor
    Epetra_SerialDenseVector*    ivecn,         //!< reaction term at time n
    Epetra_SerialDenseVector*    ivecnp,        //!< reaction term at time n+1
    Epetra_SerialDenseMatrix*    ivecnpderiv,   //!< reaction term derivaitve
    double*                      diff1          //!< main diffusivity
    )
{

  // Initialize reaction and diffusion matrices
  hdgele->diff_.Shape(nsd_,nsd_);
  hdgele->invdiff_.Shape(nsd_,nsd_);
  hdgele->Ivecn_.Shape(hdgele->ndofs_,1);
  hdgele->Ivecnp_.Shape(hdgele->ndofs_,1);
  hdgele->Imatnpderiv_.Shape(hdgele->ndofs_,hdgele->ndofs_);

  for (unsigned int i=0; i<nsd_; ++i)
    for (unsigned int j = 0; j < nsd_; ++j)
    {
      hdgele->diff_(i,j) = (*difftensor)(i,j);
      hdgele->invdiff_(i,j) = (*difftensor)(i,j);
    }
  hdgele->diff1_ = *diff1;

  Epetra_SerialDenseSolver inverseindifftensor;
  inverseindifftensor.SetMatrix(hdgele->invdiff_);
  int err = inverseindifftensor.Invert();
  if (err != 0)
    dserror("Inversion of diffusion tensor failed with errorcode %d", err);

  hdgele->Ivecn_=*ivecn;
  hdgele->Ivecnp_=*ivecnp;
  hdgele->Imatnpderiv_=*ivecnpderiv;

}

/*----------------------------------------------------------------------*
 * Element Init
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::ElementInit(
    DRT::Element *                  ele)
{
  DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));

  hdgele->Amat_.Shape(hdgele->ndofs_,hdgele->ndofs_);
  hdgele->Bmat_.Shape(hdgele->ndofs_,nsd_*hdgele->ndofs_);
  hdgele->Cmat_.Shape(hdgele->ndofs_,hdgele->onfdofs_);
  hdgele->Dmat_.Shape(nsd_*hdgele->ndofs_,nsd_*hdgele->ndofs_);
  hdgele->Emat_.Shape(hdgele->ndofs_*nsd_,hdgele->onfdofs_);
  hdgele->Gmat_.Shape(hdgele->onfdofs_,hdgele->ndofs_);
  hdgele->EmatT_.Shape(hdgele->onfdofs_,nsd_*hdgele->ndofs_);
  hdgele->Hmat_.Shape(hdgele->onfdofs_,hdgele->onfdofs_);
  hdgele->Mmat_.Shape(hdgele->ndofs_,hdgele->ndofs_);
  hdgele->Kmat_.Shape(hdgele->onfdofs_,hdgele->onfdofs_);
  hdgele->Xmat_.Shape(hdgele->ndofs_, hdgele->onfdofs_);
  hdgele->BmatMT_.Shape(nsd_*hdgele->ndofs_,hdgele->ndofs_);
  hdgele->invAMmat_.Shape(hdgele->ndofs_,hdgele->ndofs_);
  hdgele->BTAMmat_.Shape(hdgele->ndofs_*nsd_, hdgele->ndofs_);
  hdgele->invCondmat_.Shape(hdgele->ndofs_*nsd_, hdgele->ndofs_*nsd_);
  hdgele->diff_.Shape(nsd_,nsd_);
  hdgele->invdiff_.Shape(nsd_,nsd_);
  hdgele->Ivecn_.Shape(hdgele->ndofs_,1);
  hdgele->Ivecnp_.Shape(hdgele->ndofs_,1);
  hdgele->Imatnpderiv_.Shape(hdgele->ndofs_,hdgele->ndofs_);

  hdgele->SetMatInit(true);

  return;
} // ElementInit

/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::ProjectField(
    const DRT::Element*                  ele,
    DRT::Discretization &                discretization,
    Teuchos::ParameterList&              params,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2,
    DRT::Element::LocationArray&         la
)
{
  int nds_var_old = params.get<int>("nds_var_old");
  int nds_intvar_old = params.get<int>("nds_intvar_old");

  DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));

  // set change of element degree to false
  hdgele->SetPadaptEle(false);

  Teuchos::RCP<DRT::UTILS::ShapeValues<distype> > shapes_old = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(hdgele->DegreeOld(),
                                                              usescompletepoly_,
                                                              2*hdgele->DegreeOld()));

  dsassert(elevec2.M() == 0 ||
       unsigned(elevec2.M()) == shapes_->ndofs_*(nsd_+1), "Wrong size in project vector 2");

   // polynomial space to get the value of the shape function at the new points
   DRT::UTILS::PolynomialSpaceParams params_old(distype,shapes_old->degree_,usescompletepoly_);
   Teuchos::RCP<DRT::UTILS::PolynomialSpace<probdim> > polySpace_old = DRT::UTILS::PolynomialSpaceCache<probdim>::Instance().Create(params_old);

  Epetra_SerialDenseVector interiorPhi_old(shapes_old->ndofs_*(nsd_+1));

  // get node based values!
  Teuchos::RCP<const Epetra_Vector> matrix_state = params.get<Teuchos::RCP<Epetra_Vector> >("phi");

  std::vector<double> tracephi;
  DRT::UTILS::ExtractMyValues(*matrix_state,tracephi,la[nds_var_old].lm_);

  // get node based values!
  matrix_state = params.get<Teuchos::RCP<Epetra_Vector> >("intphi");
  std::vector<double> intphi;
  DRT::UTILS::ExtractMyValues(*matrix_state,intphi,la[nds_intvar_old].lm_);
  if (intphi.size() != shapes_old->ndofs_*(nsd_+1))
    dserror("node number not matching: %d vs. %d", intphi.size(), shapes_old->ndofs_*(nsd_+1));

  for (unsigned int i = 0; i < shapes_old->ndofs_*(nsd_+1); ++i)
    interiorPhi_old(i)=intphi[i];

  if (!usescompletepoly_)
  {
    // copy values if degree stays the same instead of projecting
    if(hdgele->DegreeOld() == hdgele->Degree())
      for(unsigned int i = 0; i < shapes_old->ndofs_*(nsd_+1); ++i)
        elevec2(i)=intphi[i];
    else
    {
      // set change of element degree to true
      hdgele->SetPadaptEle(true);
      Epetra_SerialDenseMatrix tempMat(shapes_->ndofs_*(nsd_+1),shapes_old->ndofs_*(nsd_+1));

      for(unsigned int i=0; i<shapes_->ndofs_;  i++)
      {
        Epetra_SerialDenseVector tempVec(shapes_old->ndofs_);
        LINALG::Matrix<nsd_,1> point(shapes_->nodexyzunit[i]);
        polySpace_old->Evaluate(point,tempVec);
        for(unsigned int j=0; j<nsd_+1; j++)
          for(unsigned int k=0; k<shapes_old->ndofs_; k++)
            tempMat(j*shapes_->ndofs_+i,j*shapes_old->ndofs_+k)=tempVec(k);
      }
      elevec2.Multiply('N','N',1.0,tempMat,interiorPhi_old,0.0);
    }

    // project trace field
    int nfdofs = 0;
    int nfdofs_old = 0;
    for (unsigned int face=0; face<nfaces_; face++)
    {
      //shape values of new element degree
      DRT::UTILS::ShapeValuesFaceParams svfparams(
          ele->Faces()[face]->Degree(),
          usescompletepoly_,
          2 * ele->Faces()[face]->Degree());
      shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      shapesface_->EvaluateFace(*ele, face);


      //shape values of old element degree
      DRT::ELEMENTS::ScaTraHDGIntFace * hdgeleface = dynamic_cast<DRT::ELEMENTS::ScaTraHDGIntFace*>(const_cast<DRT::FaceElement*>(ele->Faces()[face]));
      DRT::UTILS::ShapeValuesFaceParams svfparams_old(
          hdgeleface->DegreeOld(),
          usescompletepoly_,
          2*hdgeleface->DegreeOld());

      Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype> >shapesface_old= DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams_old);


      DRT::UTILS::PolynomialSpaceParams polyparams(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,hdgeleface->DegreeOld(),usescompletepoly_);
      Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_-1> > polySpaceFace_old = DRT::UTILS::PolynomialSpaceCache<nsd_-1>::Instance().Create(polyparams);

      Epetra_SerialDenseVector tracePhi_face_old(shapesface_old->nfdofs_);

      for(unsigned int i=0; i<shapesface_old->nfdofs_; i++)
        tracePhi_face_old(i) = tracephi[nfdofs_old+i];

      if(ele->Faces()[face]->Degree() == hdgeleface->DegreeOld())
        for(unsigned int i=0; i<shapesface_old->nfdofs_; i++)
          elevec1(nfdofs+i)=tracephi[nfdofs_old+i];
      else
      {
        // set change of element degree to true
        hdgele->SetPadaptEle(true);

        Epetra_SerialDenseMatrix tempMat1(shapesface_->nfdofs_,shapesface_old->nfdofs_);

        Epetra_SerialDenseVector tempVec2(shapesface_->nfdofs_);

        for(unsigned int i=0; i<shapesface_->nfdofs_;  i++)
        {
          Epetra_SerialDenseVector tempVec(shapesface_old->nfdofs_);
          LINALG::Matrix<nsd_-1,1> point(shapesface_->nodexyzunit[i]);

          polySpaceFace_old->Evaluate(point,tempVec);
          for(unsigned int k=0; k<shapesface_old->nfdofs_; k++)
            tempMat1(i,k)=tempVec(k);
        }

        tempVec2.Multiply('N','N',1.0,tempMat1,tracePhi_face_old,0.0);

        for (unsigned int i=0; i<shapesface_->nfdofs_; i++)
          elevec1(nfdofs+i) = tempVec2(i);
      }

      nfdofs += shapesface_->nfdofs_;
      nfdofs_old += shapesface_old->nfdofs_;
    }
  }
  else
  {
    if(hdgele->DegreeOld() != hdgele->Degree())
      hdgele->SetPadaptEle(true); // set change of element degree to true

    unsigned int size_ndofs = std::min(shapes_old->ndofs_,shapes_->ndofs_);
    for(unsigned int i=0; i<nsd_+1; i++)
      for(unsigned int j=0; j<size_ndofs; j++)
        elevec2(i*shapes_->ndofs_+j)=interiorPhi_old(i*shapes_old->ndofs_+j);

    int nfdofs = 0;
    int nfdofs_old = 0;
    for (unsigned int face=0; face<nfaces_; face++)
    {
      //shape values of new element degree
      DRT::UTILS::ShapeValuesFaceParams svfparams(
          ele->Faces()[face]->Degree(),
          usescompletepoly_,
          2 * ele->Faces()[face]->Degree());
      shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      shapesface_->EvaluateFace(*ele, face);

      //shape values of old element degree
      DRT::ELEMENTS::ScaTraHDGIntFace * hdgeleface = dynamic_cast<DRT::ELEMENTS::ScaTraHDGIntFace*>(const_cast<DRT::FaceElement*>(ele->Faces()[face]));
      DRT::UTILS::ShapeValuesFaceParams svfparams_old(
          hdgeleface->DegreeOld(),
          usescompletepoly_,
          2*hdgeleface->DegreeOld());

      Teuchos::RCP<DRT::UTILS::ShapeValuesFace<distype> >shapesface_old= DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams_old);


      DRT::UTILS::PolynomialSpaceParams polyparams(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,hdgeleface->DegreeOld(),usescompletepoly_);
      Teuchos::RCP<DRT::UTILS::PolynomialSpace<nsd_-1> > polySpaceFace_old = DRT::UTILS::PolynomialSpaceCache<nsd_-1>::Instance().Create(polyparams);

//      Epetra_SerialDenseVector tracePhi_face_old(shapesface_old->nfdofs_);

      if(ele->Faces()[face]->Degree() != hdgeleface->DegreeOld())
        hdgele->SetPadaptEle(true); // set change of element degree to true

      unsigned int size_nfdofs = std::min(shapesface_->nfdofs_,shapesface_old->nfdofs_);
      for (unsigned int i=0; i<size_nfdofs; i++)
        elevec1(nfdofs+i) = tracephi[nfdofs_old+i];

      nfdofs += shapesface_->nfdofs_;
      nfdofs_old += shapesface_old->nfdofs_;

    }

  }

  return 0;

}//ProjectField


/*----------------------------------------------------------------------*
 * Calc P-Adaptivity
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ScaTraEleCalcHDG<distype,probdim>::CalcPAdaptivity(
    const DRT::Element*                  ele,
    DRT::Discretization &                discretization,
    Teuchos::ParameterList&              params
)
{

  DRT::ELEMENTS::ScaTraHDG * hdgele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(const_cast<DRT::Element*>(ele));

  Epetra_SerialDenseVector tempinteriorgradphinp(hdgele->ndofs_*nsd_);
  for (unsigned int i = 0; i < hdgele->ndofs_*nsd_; ++i)
    tempinteriorgradphinp(i) = interiorPhinp_(hdgele->ndofs_+i);

  double error = 0;
  int sumindex = 0;
  for (unsigned int nface = 0; nface < nfaces_; ++nface)
  {

    DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele->Faces()[nface]->Degree(),
        shapes_->usescompletepoly_,
        2 * ele->Faces()[nface]->Degree());

    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, nface);

    Epetra_SerialDenseMatrix QMat(shapesface_->nqpoints_,hdgele->ndofs_*nsd_);
    Epetra_SerialDenseMatrix QMatW(shapesface_->nqpoints_,hdgele->ndofs_*nsd_);
    Epetra_SerialDenseMatrix UMat(shapesface_->nqpoints_,hdgele->ndofs_+shapesface_->nfdofs_);
    Epetra_SerialDenseMatrix UMatW(shapesface_->nqpoints_,hdgele->ndofs_+shapesface_->nfdofs_);

    Epetra_SerialDenseVector tempinteriorphinp(hdgele->ndofs_+shapesface_->nfdofs_);
    for (unsigned int i = 0; i < hdgele->ndofs_; i++)
      tempinteriorphinp(i) = interiorPhinp_(i);
    for (unsigned int i = 0; i < shapesface_->nfdofs_; i++)
      tempinteriorphinp(hdgele->ndofs_+i) = tracen_(sumindex+i);

    // loop over quadrature points
    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      // loop over shape functions
      for (unsigned int i = 0; i < hdgele->ndofs_; ++i)
      {
        UMat(q , i) = shapesface_->shfunctI(i, q);
        UMatW(q , i) = shapesface_->shfunctI(i, q) * shapesface_->jfac(q);
        for (unsigned int k = 0; k < nsd_; ++k)
        {
          QMat(q ,hdgele->ndofs_ * k + i) = shapesface_->shfunctI(i, q) * shapesface_->normal(k);
          QMatW(q ,hdgele->ndofs_ * k + i) = shapesface_->shfunctI(i, q) * shapesface_->jfac(q) * shapesface_->normal(k);
        }
      }
      // loop over face shape functions
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        UMat(q,i+hdgele->ndofs_) = -shapesface_->shfunct(i, q);
        UMatW(q,i+hdgele->ndofs_) = -shapesface_->shfunct(i, q) * shapesface_->jfac(q);
      }
    }
    sumindex += shapesface_->nfdofs_;

    Epetra_SerialDenseVector tempVec1(shapesface_->nqpoints_);
    Epetra_SerialDenseVector tempVec2(shapesface_->nqpoints_);
    Epetra_SerialDenseVector tempVec3(shapesface_->nqpoints_);
    Epetra_SerialDenseVector tempVec4(shapesface_->nqpoints_);


    tempVec1.Multiply('N','N',1.0,QMatW,tempinteriorgradphinp,0.0);
    tempVec2.Multiply('N','N',1.0,QMat,tempinteriorgradphinp,0.0);
    tempVec3.Multiply('N','N',1.0,UMatW,tempinteriorphinp,0.0);
    tempVec4.Multiply('N','N',1.0,UMat,tempinteriorphinp,0.0);

    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      error += tempVec1(q)*tempVec2(q)+tempVec3(q)*tempVec4(q)-2*tempVec1(q)*tempVec4(q);
    }

  }

  params.set<double>("error",error);


  return 0;
}//CalcPAdaptivity


// template classes
// 1D elements
//template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::line2,1>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::line2,2>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::line2,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::line3,1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::quad4,2>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::quad4,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::quad9,2>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::nurbs9,2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::hex8,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::hex27,3>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::tet4,3>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::tet10,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::pyramid5,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcHDG<DRT::Element::nurbs27>;
