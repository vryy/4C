/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc.cpp

\brief evaluation of scatra boundary terms at integration points

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/

// general Butler-Volmer is activated if the define-flag ButlerVolmer_Shifted is off
// the shifted version of Butler-Volmer is activated if the define-flag ButlerVolmer_Shifted is on
// define-flag PERCENT: how much is the curve shifted

#include <cstdlib>

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_lib/drt_globalproblem.H" // for curves and functions
#include "../drt_lib/standardtypes_cpp.H" // for EPS12 and so on

#include "../drt_mat/fourieriso.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/thermostvenantkirchhoff.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "scatra_ele.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_boundary_calc.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ScaTraEleBoundaryCalc(const int numdofpernode, const int numscal)
 : scatraparamstimint_(DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance()), // params for time integration
   scatraparams_(DRT::ELEMENTS::ScaTraEleParameterStd::Instance()),
   numdofpernode_(numdofpernode),
   numscal_(numscal),
   xyze_(true),  // initialize to zero
   weights_(true),
   myknots_(nsd_),
   mypknots_(nsd_+1),
   normalfac_(1.0),
   edispnp_(true),
   diffus_(numscal_,0),
   //valence_(numscal_,0),
   shcacp_(0.0),
   xsi_(true),
   funct_(true),
   deriv_(true),
   derxy_(true),
   normal_(true),
   velint_(true),
   metrictensor_(true)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::SetupCalc(DRT::FaceElement*              ele,
                                                          Teuchos::ParameterList&           params,
                                                          DRT::Discretization&              discretization)
{
  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

  // get additional state vector for ALE case: grid displacement
  if(scatraparams_->IsAle())
  {
    const Teuchos::RCP<Epetra_MultiVector> dispnp = params.get< Teuchos::RCP<Epetra_MultiVector> >("dispnp",Teuchos::null);
    if (dispnp==Teuchos::null) dserror("Cannot get state vector 'dispnp'");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,edispnp_,dispnp,nsd_+1);
    // add nodal displacements
    xyze_ += edispnp_;
  }
  else edispnp_.Clear();

  // Now do the nurbs specific stuff (for isogeometric elements)
  if(DRT::NURBS::IsNurbs(distype))
  {
    // for isogeometric elements --- get knotvectors for parent
    // element and boundary element, get weights
    bool zero_size = DRT::NURBS::GetKnotVectorAndWeightsForNurbsBoundary(
        ele, ele->FaceParentNumber(), ele->ParentElement()->Id(), discretization, mypknots_, myknots_, weights_, normalfac_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size) return -1;
  } // Nurbs specific stuff

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate element                                          fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::Evaluate(
    DRT::FaceElement*                   ele,
    Teuchos::ParameterList&             params,
    DRT::Discretization&                discretization,
    std::vector<int>&                   lm,
    Epetra_SerialDenseMatrix&           elemat1_epetra,
    Epetra_SerialDenseMatrix&           elemat2_epetra,
    Epetra_SerialDenseVector&           elevec1_epetra,
    Epetra_SerialDenseVector&           elevec2_epetra,
    Epetra_SerialDenseVector&           elevec3_epetra
)
{
  // check for the action parameter
  const SCATRA::BoundaryAction action = DRT::INPUT::get<SCATRA::BoundaryAction>(params,"action");

  // setup
  if(SetupCalc(ele,params,discretization) == -1)
    return 0;

  // evaluate action
  EvaluateAction(
      ele,
      params,
      discretization,
      action,
      lm,
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateAction(
    DRT::FaceElement*                 ele,
    Teuchos::ParameterList&           params,
    DRT::Discretization&              discretization,
    SCATRA::BoundaryAction            action,
    std::vector<int>&                 lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseMatrix&         elemat2_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    Epetra_SerialDenseVector&         elevec2_epetra,
    Epetra_SerialDenseVector&         elevec3_epetra
)
{
  switch (action)
  {
  case SCATRA::bd_calc_normal_vectors:
  {
    CalcNormalVectors(params,ele);
    break;
  }
  case SCATRA::bd_integrate_shape_functions:
  {
    // NOTE: add area value only for elements which are NOT ghosted!
    const bool addarea = (ele->Owner() == discretization.Comm().MyPID());
    IntegrateShapeFunctions(ele,params,elevec1_epetra,addarea);

    break;
  }
  case SCATRA::bd_calc_Neumann:
  {
    DRT::Condition* condition = params.get<DRT::Condition*>("condition");
    if(condition == NULL)
      dserror("Cannot access Neumann boundary condition!");

    EvaluateNeumann(ele,params,discretization,*condition,lm,elevec1_epetra,1.);

    break;
  }
  case SCATRA::bd_calc_Neumann_inflow:
  {
    NeumannInflow(
        ele,
        params,
        discretization,
        lm,
        elemat1_epetra,
        elevec1_epetra
        );

    break;
  }
  case SCATRA::bd_calc_convective_heat_transfer:
  {
    // get the parent element including its material
    DRT::Element* parentele = ele->ParentElement();
    Teuchos::RCP<MAT::Material> mat = parentele->Material();

    // get values of scalar
    Teuchos::RCP<const Epetra_Vector> phinp  = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // get condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null)
      dserror("Cannot access condition 'ThermoConvections'!");

    // get heat transfer coefficient and surrounding temperature
    const double heatranscoeff = cond->GetDouble("coeff");
    const double surtemp       = cond->GetDouble("surtemp");

    ConvectiveHeatTransfer(
        ele,
        mat,
        ephinp,
        elemat1_epetra,
        elevec1_epetra,
        heatranscoeff,
        surtemp
        );

    break;
  }
  case SCATRA::bd_calc_weak_Dirichlet:
  {
    // get the parent element including its material
    DRT::Element* parentele = ele->ParentElement();
    Teuchos::RCP<MAT::Material> mat = parentele->Material();

    switch (distype)
    {
    if (numscal_>1) dserror("not yet implemented for more than one scalar\n");

    // 2D:
    case DRT::Element::line2:
    {
      if(ele->ParentElement()->Shape()==DRT::Element::quad4)
      {
        WeakDirichlet<DRT::Element::line2,DRT::Element::quad4>(ele,
            params,
            discretization,
            mat,
            elemat1_epetra,
            elevec1_epetra);
      }
      else
      {
        dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
      }
      break;
    }
    // 3D:
    case DRT::Element::quad4:
    {
      if(ele->ParentElement()->Shape()==DRT::Element::hex8)
      {
        WeakDirichlet<DRT::Element::quad4,DRT::Element::hex8>(ele,
            params,
            discretization,
            mat,
            elemat1_epetra,
            elevec1_epetra);
      }
      else
      {
        dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
      }
      break;
    }
    default:
    {
      dserror("not implemented yet\n");
    }
    }

    break;
  }
  case SCATRA::bd_calc_fs3i_surface_permeability:
  {
    if(scatraparamstimint_->IsGenAlpha() or not scatraparamstimint_->IsIncremental())
      dserror("calc_surface_permeability: chosen option not available");

    // get values of scalar
    Teuchos::RCP<const Epetra_Vector> phinp  = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    //get values of wall shear stress (the euclidean norm of the stress tensor)
    Teuchos::RCP<const Epetra_Vector> wss  = discretization.GetState("WallShearStress");
    if (wss==Teuchos::null)
      dserror("Cannot get state vector 'WallShearStress'");

    // extract local values from global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);
    std::vector<double> ewss(lm.size());
    DRT::UTILS::ExtractMyValues(*wss,ewss,lm);

    // get current condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null) dserror("Cannot access condition 'SurfacePermeability'");

    const std::vector<int>* onoff = cond->Get<std::vector<int> > ("onoff");

    const double perm = cond->GetDouble("permeability coefficient");

    //get flag if concentration flux across membrane is affected by local wall shear stresses: 0->no 1->yes
    bool wss_onoff = (bool)cond->GetInt("wss onoff");
    const std::vector<double>* coeffs = cond->Get<std::vector<double> > ("wss coeffs");

    //calculate factor that accounts for WSS for each integration point
    std::vector<double> f_WSS = WSSinfluence(ewss,wss_onoff,coeffs);


    EvaluateSurfacePermeability(
        ele,
        ephinp,
        f_WSS,
        elemat1_epetra,
        elevec1_epetra,
        onoff,
        perm
    );

    break;
  }
  case SCATRA::bd_calc_fps3i_surface_permeability:
  {
    // safety checks
    if(scatraparamstimint_->IsGenAlpha() or not scatraparamstimint_->IsIncremental())
      dserror("calc_surface_permeability: chosen option not available");

    //// get velocity values at nodes
    //const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("velocity field",Teuchos::null);

    //// we deal with a (nsd_+1)-dimensional flow field
    //LINALG::Matrix<nsd_+1,nen_>  evel(true);
    //DRT::UTILS::ExtractMyNodeBasedValues(ele,evel,velocity,nsd_+1);

    // get values of scalar transport
    Teuchos::RCP<const Epetra_Vector> phinp  = discretization.GetState("phinp");
    if (phinp==Teuchos::null)
      dserror("Cannot get state vector 'phinp'");

    //get values of wall shear stress (the euclidean norm of the stress tensor)
    Teuchos::RCP<const Epetra_Vector> wss  = discretization.GetState("WallShearStress");
    if (wss==Teuchos::null)
      dserror("Cannot get state vector 'WallShearStress'");

    //get values of pressure
    Teuchos::RCP<const Epetra_Vector> pressure  = discretization.GetState("Pressure");
    if (pressure==Teuchos::null)
      dserror("Cannot get state vector 'Pressure'");

    //get mean concentration in the interface (i.e. within the membrane)
    Teuchos::RCP<const Epetra_Vector> phibar  = discretization.GetState("MeanConcentration");
    if (phibar==Teuchos::null)
      dserror("Cannot get state vector 'MeanConcentration'");

    // extract local values from global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);
    std::vector<double> ewss(lm.size());
    DRT::UTILS::ExtractMyValues(*wss,ewss,lm);
    std::vector<double> epressure(lm.size());
    DRT::UTILS::ExtractMyValues(*pressure,epressure,lm);
    std::vector<double> ephibar(lm.size());
    DRT::UTILS::ExtractMyValues(*phibar,ephibar,lm);


    // get current condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null)
      dserror("Cannot access condition 'DESIGN SCATRA COUPLING SURF CONDITIONS'");

    const std::vector<int>* onoff = cond->Get<std::vector<int> > ("onoff");

    //get the standard permeability of the interface
    double perm = cond->GetDouble("permeability coefficient");

    //get flag if concentration flux across membrane is affected by local wall shear stresses: 0->no 1->yes
    bool wss_onoff = (bool)cond->GetInt("wss onoff");
    const std::vector<double>* coeffs = cond->Get<std::vector<double> > ("wss coeffs");

    //calculate factor that accounts for WSS for each integration point
    std::vector<double> f_WSS = WSSinfluence(ewss,wss_onoff,coeffs);

    //hydraulic conductivity at interface
    const double conductivity = cond->GetDouble("hydraulic conductivity");

    //Staverman filtration coefficient at interface
    const double sigma = cond->GetDouble("filtration coefficient");

    EvaluateKedemKatchalsky(
        ele,
        ephinp,
        epressure,
        ephibar,
        f_WSS,
        elemat1_epetra,
        elevec1_epetra,
        onoff,
        perm,
        conductivity,
        sigma
    );

    break;
  }
  case SCATRA::bd_add_convective_mass_flux:
  {
    //calculate integral of convective mass/heat flux
    // NOTE: since results are added to a global vector via normal assembly
    //       it would be wrong to suppress results for a ghosted boundary!

    // get actual values of transported scalars
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // get velocity values at nodes
    const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("velocity field",Teuchos::null);

    // we deal with a (nsd_+1)-dimensional flow field
    LINALG::Matrix<nsd_+1,nen_>  evel(true);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evel,velocity,nsd_+1);

    // for the moment we ignore the return values of this method
    CalcConvectiveFlux(ele,ephinp,evel,elevec1_epetra);
    //vector<double> locfluxintegral = CalcConvectiveFlux(ele,ephinp,evel,elevec1_epetra);
    //std::cout<<"locfluxintegral[0] = "<<locfluxintegral[0]<<std::endl;

    break;
  }
  case SCATRA::bd_calc_s2icoupling:
  {
    EvaluateS2ICoupling(ele,
        params,
        discretization,
        lm,
        elemat1_epetra,
        elemat2_epetra,
        elevec1_epetra);

    break;
  }
  case SCATRA::bd_calc_Robin:
  {
    CalcRobinBoundary(
        ele,
        params,
        discretization,
        lm,
        elemat1_epetra,
        elevec1_epetra,
        1.
        );
    break;
  }
  default:
  {
    dserror("Not acting on this boundary action. Forgot implementation?");
    break;
  }
  }

  return 0;

}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition                        gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateNeumann(
    DRT::FaceElement*                   ele,
    Teuchos::ParameterList&             params,
    DRT::Discretization&                discretization,
    DRT::Condition&                     condition,
    std::vector<int>&                   lm,
    Epetra_SerialDenseVector&           elevec1,
    const double                        scalar
    )
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = scatraparamstimint_->Time();
  if (time<0.0) usetime = false;

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const int numdof = condition.GetInt("numdof");
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>* curve    = condition.Get<std::vector<int> >   ("curve");
  const std::vector<int>*    func  = condition.Get<std::vector<int> >   ("funct");

  if (numdofpernode_!=numdof)
    dserror("The NUMDOF you have entered in your TRANSPORT NEUMANN CONDITION does not equal the number of scalars.");

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    double fac = EvalShapeFuncAndIntFac(intpoints,iquad);

    // factor given by spatial function
    double functfac = 1.0;
    // factor given by temporal curve
    double curvefac = 1.0;

    // determine global coordinates of current Gauss point
    double coordgp[3];   // we always need three coordinates for function evaluation!
    for(int i=0; i<3; ++i)
      coordgp[i] = 0.;
    for(int i=0; i<nsd_; ++i)
    {
      coordgp[i] = 0.;
      for(int j=0; j<nen_; ++j)
        coordgp[i] += xyze_(i,j)*funct_(j);
    }

    int functnum = -1;
    int curvenum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation

    for(int dof=0; dof<numdofpernode_; ++dof)
    {
      if ((*onoff)[dof]) // is this dof activated?
      {
        // find out whether we will use a time curve and get the factor
        if (curve) curvenum = (*curve)[dof];

        if (curvenum>=0 && usetime)
        {
          // evaluate curve at current time
          curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
        }
        else
          curvefac = 1.0;


        // factor given by spatial function
        if(func)
          functnum = (*func)[dof];

        if(functnum>0)
        {
          // evaluate function at current Gauss point (provide always 3D coordinates!)
          functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof,coordgpref,time,NULL);
        }
        else
          functfac = 1.;

        const double val_fac_funct_curve_fac = (*val)[dof]*fac*functfac*curvefac;

        for(int node=0; node<nen_; ++node)
          //TODO: with or without eps_
          elevec1[node*numdofpernode_+dof] += scalar*funct_(node)*val_fac_funct_curve_fac;
      } // if ((*onoff)[dof])
    } // loop over dofs
  } // loop over integration points

  return 0;
}


/*----------------------------------------------------------------------*
 | calculate normals vectors                                   vg 03/09 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcNormalVectors(
    Teuchos::ParameterList&             params,
    DRT::FaceElement*                   ele
    )
{
  // access the global vector
  const Teuchos::RCP<Epetra_MultiVector> normals = params.get< Teuchos::RCP<Epetra_MultiVector> >("normal vectors",Teuchos::null);
  if (normals == Teuchos::null) dserror("Could not access vector 'normal vectors'");

  // determine constant outer normal to this element
  GetConstNormal(normal_,xyze_);

  // loop over the element nodes
  for (int j=0;j<nen_;j++)
  {
    const int nodegid = (ele->Nodes()[j])->Id();
    if (normals->Map().MyGID(nodegid) )
    {// OK, the node belongs to this processor

      // scaling to a unit vector is performed on the global level after
      // assembly of nodal contributions since we have no reliable information
      // about the number of boundary elements adjacent to a node
      for (int dim=0; dim<(nsd_+1); dim++)
      {
        normals->SumIntoGlobalValue(nodegid,dim,normal_(dim));
      }
    }
    //else: the node belongs to another processor; the ghosted
    //      element will contribute the right value on that proc
  }

}


/*----------------------------------------------------------------------*
 | calculate Neumann inflow boundary conditions                vg 03/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::NeumannInflow(
    const DRT::FaceElement*                   ele,
    Teuchos::ParameterList&                   params,
    DRT::Discretization&                      discretization,
    std::vector<int>&                         lm,
    Epetra_SerialDenseMatrix&                 emat,
    Epetra_SerialDenseVector&                 erhs
    )
{
  // get parent element
  DRT::Element* parentele = ele->ParentElement();

  // get material of parent element
  Teuchos::RCP<MAT::Material> material = parentele->Material();

  // we don't know the parent element's lm vector; so we have to build it here
  const int nenparent = parentele->NumNode();
  std::vector<int> lmparent(nenparent);
  std::vector<int> lmparentowner;
  std::vector<int> lmparentstride;
  parentele->LocationVector(discretization,lmparent,lmparentowner,lmparentstride);

  // get values of scalar
  Teuchos::RCP<const Epetra_Vector> phinp  = discretization.GetState("phinp");
  if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

  // extract local values from global vector
  std::vector<double> ephinp(lm.size());
  DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

  // get velocity values at nodes
  const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("convective velocity field",Teuchos::null);

  // we deal with a (nsd_+1)-dimensional flow field
  Epetra_SerialDenseVector evel((nsd_+1)*nenparent);
  DRT::UTILS::ExtractMyNodeBasedValues(parentele,evel,velocity,nsd_+1);

  // insert velocity field into element array
  LINALG::Matrix<nsd_+1,nen_> evelnp;
  for (int i=0;i<nen_;++i)
  {
    for (int idim=0 ; idim < nsd_+1 ; idim++)
    {
      evelnp(idim,i) = evel[idim + i*(nsd_+1)];
    }
  }

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<nen_,1> phinod(true);

  // loop over all scalars
  for(int k=0;k<numdofpernode_;++k)
  {
    // compute scalar values at nodes
    for (int inode=0; inode< nen_;++inode)
    {
      phinod(inode) = ephinp[inode*numdofpernode_+k];
    }

    // loop over all integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints,iquad,&normal_);

      // get velocity at integration point
      velint_.Multiply(evelnp,funct_);

      // normal velocity
      const double normvel = velint_.Dot(normal_);

      if (normvel<-0.0001)
      {
        // set density to 1.0
        double dens = GetDensity(material,ephinp,phinod);

        // integration factor for left-hand side
        const double lhsfac = dens*normvel*scatraparamstimint_->TimeFac()*fac;

        // integration factor for right-hand side
        double rhsfac = 0.0;
        if(scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac/scatraparamstimint_->AlphaF();
        else if(not scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac*(1.0-scatraparamstimint_->AlphaF())/scatraparamstimint_->AlphaF();
        else if(scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac;

        // matrix
        for (int vi=0; vi<nen_; ++vi)
        {
          const double vlhs = lhsfac*funct_(vi);

          const int fvi = vi*numdofpernode_+k;

          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) -= vlhs*funct_(ui);
          }
        }

        // scalar at integration point
        const double phi = funct_.Dot(phinod);

        // rhs
        const double vrhs = rhsfac*phi;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] += vrhs*funct_(vi);
        }
      }
    }
  }

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::NeumannInflow


/*----------------------------------------------------------------------*
 | get density at integration point                          fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const double DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::GetDensity(
    Teuchos::RCP<const MAT::Material>   material,
    const std::vector<double>&          ephinp,
    const LINALG::Matrix<nen_,1>&       phinod
    )
{
  // initialization
  double density(0.);

  // get density depending on material
  switch(material->MaterialType())
  {
  case INPAR::MAT::m_matlist:
  {
    const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

    const int matid = actmat->MatID(0);

    if(actmat->MaterialById(matid)->MaterialType() == INPAR::MAT::m_scatra)
      // set density to unity
      density = 1.;
    else
      dserror("type of material found in material list is not supported");

    break;
  }

  case INPAR::MAT::m_matlist_reactions:
  case INPAR::MAT::m_scatra:
  {
    // set density to unity
    density = 1.;

    break;
  }

  default:
  {
    dserror("Invalid material type!");
    break;
  }
  }

  return density;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::GetDensity


/*----------------------------------------------------------------------*
 | calculate integral of convective flux across boundary      gjb 11/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
std::vector<double> DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcConvectiveFlux(
    const DRT::Element*                 ele,
    const std::vector<double>&          ephinp,
    const LINALG::Matrix<nsd_+1,nen_>&  evelnp,
    Epetra_SerialDenseVector&           erhs
)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<nen_,1> phinod(true);

  std::vector<double> integralflux(numscal_);

  // loop over all scalars
  for(int k=0;k<numscal_;++k)
  {
    integralflux[k] = 0.0;

    // compute scalar values at nodes
    for (int inode=0; inode< nen_;++inode)
    {
      phinod(inode) = ephinp[inode*numdofpernode_+k];
    }

    // loop over all integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints,iquad,&normal_);

      // get velocity at integration point
      velint_.Multiply(evelnp,funct_);

      // normal velocity (note: normal_ is already a unit(!) normal)
      const double normvel = velint_.Dot(normal_);

      // scalar at integration point
      const double phi = funct_.Dot(phinod);

      const double val = phi*normvel*fac;
      integralflux[k] += val;
      // add contribution to provided vector (distribute over nodes using shape fct.)
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;
        erhs[fvi] += val*funct_(vi);
      }
    }
  }

  return integralflux;

} //ScaTraEleBoundaryCalc<distype>::ConvectiveFlux

/*----------------------------------------------------------------------*
 | calculate boundary cond. due to convective heat transfer    vg 10/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ConvectiveHeatTransfer(
    const DRT::Element*                 ele,
    Teuchos::RCP<const MAT::Material>   material,
    const std::vector<double>&          ephinp,
    Epetra_SerialDenseMatrix&           emat,
    Epetra_SerialDenseVector&           erhs,
    const double                        heatranscoeff,
    const double                        surtemp
    )
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<nen_,1> phinod(true);

  // loop over all scalars
  for(int k=0;k<numdofpernode_;++k)
  {
    // compute scalar values at nodes
    for (int inode=0; inode< nen_;++inode)
    {
      phinod(inode) = ephinp[inode*numdofpernode_+k];
    }

    // loop over all integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndIntFac(intpoints,iquad,&normal_);

      // get specific heat capacity at constant volume
      double shc = 0.0;
      if (material->MaterialType() == INPAR::MAT::m_th_fourier_iso)
      {
        const MAT::FourierIso* actmat = static_cast<const MAT::FourierIso*>(material.get());

        shc = actmat->Capacity();
      }
      else if (material->MaterialType() == INPAR::MAT::m_thermostvenant)
      {
        const MAT::ThermoStVenantKirchhoff* actmat = static_cast<const MAT::ThermoStVenantKirchhoff*>(material.get());

        shc = actmat->Capacity();
      }
      else dserror("Material type is not supported for convective heat transfer!");

      // integration factor for left-hand side
      const double lhsfac = heatranscoeff*scatraparamstimint_->TimeFac()*fac/shc;

      // integration factor for right-hand side
      double rhsfac = 0.0;
      if(scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac/scatraparamstimint_->AlphaF();
      else if(not scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac*(1.0-scatraparamstimint_->AlphaF())/scatraparamstimint_->AlphaF();
      else if(scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac;

      // matrix
      for (int vi=0; vi<nen_; ++vi)
      {
        const double vlhs = lhsfac*funct_(vi);

        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) -= vlhs*funct_(ui);
        }
      }

      // scalar at integration point
      const double phi = funct_.Dot(phinod);

      // rhs
      const double vrhs = rhsfac*(phi-surtemp);
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] += vrhs*funct_(vi);
      }
    }
  }

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ConvectiveHeatTransfer


/*----------------------------------------------------------------------*
 | evaluate shape functions and int. factor at int. point     gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeFuncAndIntFac(
    const DRT::UTILS::IntPointsAndWeights<nsd_>&   intpoints,   ///< integration points
    const int                                      iquad,       ///< id of current Gauss point
    LINALG::Matrix<1 + nsd_,1>*                    normalvec    ///< normal vector at Gauss point(optional)
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {xsi_(idim) = gpcoord[idim];}

  if(not DRT::NURBS::IsNurbs(distype))
  {
    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  }
  else // nurbs elements are always somewhat special...
  {
    DRT::NURBS::UTILS::nurbs_get_funct_deriv(
        funct_  ,
        deriv_  ,
        xsi_    ,
        myknots_,
        weights_,
        distype );
  }

  // the metric tensor and the area of an infinitesimal surface/line element
  // optional: get normal at integration point as well
  // Note: this is NOT yet a unit normal. Its norm corresponds to the area/length of the element
  double drs(0.0);
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor_,drs,normalvec);

  // for nurbs elements the normal vector must be scaled with a special orientation factor!!
  if(DRT::NURBS::IsNurbs(distype))
  {
    if (normalvec != NULL)
      normal_.Scale(normalfac_);
  }

  // return the integration factor
  return intpoints.IP().qwgt[iquad] * drs;
}


/*----------------------------------------------------------------------*
 | get constant normal                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::GetConstNormal(
    LINALG::Matrix<nsd_+1,1>&          normal,
    const LINALG::Matrix<nsd_+1,nen_>&  xyze
    )
{
  // determine normal to this element
  if(not DRT::NURBS::IsNurbs(distype))
  {
  switch(nsd_)
  {
  case 2:
  {
    LINALG::Matrix<3,1> dist1(true), dist2(true);
    for (int i=0; i<3; i++)
    {
      dist1(i) = xyze(i,1)-xyze(i,0);
      dist2(i) = xyze(i,2)-xyze(i,0);
    }

    normal(0) = dist1(1)*dist2(2) - dist1(2)*dist2(1);
    normal(1) = dist1(2)*dist2(0) - dist1(0)*dist2(2);
    normal(2) = dist1(0)*dist2(1) - dist1(1)*dist2(0);
  }
  break;
  case 1:
  {
    normal(0) = xyze(1,1) - xyze(1,0);
    normal(1) = (-1.0)*(xyze(0,1) - xyze(0,0));
  }
  break;
  default:
    dserror("Illegal number of space dimensions: %d",nsd_);
    break;
  } // switch(nsd)
  }
  else // NURBS case
  {
    // ToDo: this is only a temporary solution in order to have something here.
    // Current handling of node-based normal vectors not applicable in NURBS case
#if 0
    // use one integration point at element center
    const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToStabGaussRule<distype>::rule);
    // hack: ele-id = -1
    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    // this is already part of this function call
    EvalShapeFuncAndIntFac(intpoints,0,&normal);
#endif
    normal(0)= 1.0;
  }

  // length of normal to this element
  const double length = normal.Norm2();
  // outward-pointing normal of length 1.0
  if (length > EPS10)
    normal.Scale(1/length);
  else
    dserror("Zero length for element normal");

  return;
} // ScaTraEleBoundaryCalc<distype>::GetConstNormal


/*----------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition       fang 10/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateS2ICoupling(
    const DRT::Element*         ele,              ///< current boundary element
    Teuchos::ParameterList&     params,           ///< parameter list
    DRT::Discretization&        discretization,   ///< discretization
    std::vector<int>&           lm,               ///< location vector
    Epetra_SerialDenseMatrix&   eslavematrix,     ///< element matrix for slave side
    Epetra_SerialDenseMatrix&   emastermatrix,    ///< element matrix for master side
    Epetra_SerialDenseVector&   eslaveresidual    ///< element residual for slave side
    )
{
  // get global and interface state vectors
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  Teuchos::RCP<const Epetra_Vector> imasterphinp = discretization.GetState("imasterphinp");
  if (phinp == Teuchos::null or imasterphinp == Teuchos::null)
    dserror("Cannot get state vector 'phinp' or 'imasterphinp'!");

  // extract local nodal values on present and opposite sides of scatra-scatra interface
  std::vector<double> eslavephinpvec(lm.size());
  DRT::UTILS::ExtractMyValues(*phinp,eslavephinpvec,lm);
  std::vector<LINALG::Matrix<nen_,1> > eslavephinp(numscal_);
  std::vector<double> emasterphinpvec(lm.size());
  DRT::UTILS::ExtractMyValues(*imasterphinp,emasterphinpvec,lm);
  std::vector<LINALG::Matrix<nen_,1> > emasterphinp(numscal_);
  for(int inode=0; inode<nen_; ++inode)
  {
    for(int k=0; k<numscal_; ++k)
    {
      eslavephinp[k](inode,0) = eslavephinpvec[inode*numscal_+k];
      emasterphinp[k](inode,0) = emasterphinpvec[inode*numscal_+k];
    }
  }

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // access kinetic model associated with current condition
  const int kineticmodel = s2icondition->GetInt("kinetic model");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeFuncAndIntFac(intpoints,gpid);

    // evaluate overall integration factors
    const double timefacfac = scatraparamstimint_->TimeFac()*fac;
    const double timefacrhsfac = scatraparamstimint_->TimeFacRhs()*fac;
    if (timefacfac < 0. or timefacrhsfac < 0.)
      dserror("Integration factor is negative!");

    // loop over scalars
    for(int k=0; k<numscal_; ++k)
    {
      // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
      const double slavephiint = funct_.Dot(eslavephinp[k]);
      const double masterphiint = funct_.Dot(emasterphinp[k]);

      // compute matrix and vector contributions according to kinetic model for current scatra-scatra interface coupling condition
      switch(kineticmodel)
      {
        // constant permeability model
        case INPAR::S2I::kinetics_constperm:
        {
          // access real vector of constant permeabilities
          const std::vector<double>* permeabilities = s2icondition->GetMutable<std::vector<double> >("permeabilities");
          if(permeabilities == NULL)
            dserror("Cannot access vector of permeabilities for scatra-scatra interface coupling!");
          if(permeabilities->size() != (unsigned) numscal_)
            dserror("Number of permeabilities does not match number of scalars!");

          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi = vi*numscal_+k;

            for (int ui=0; ui<nen_; ++ui)
            {
              eslavematrix(fvi,ui*numscal_+k) += funct_(vi)*(*permeabilities)[k]*funct_(ui)*timefacfac;
              emastermatrix(fvi,ui*numscal_+k) -= funct_(vi)*(*permeabilities)[k]*funct_(ui)*timefacfac;
            }

            eslaveresidual[fvi] -= funct_(vi)*(*permeabilities)[k]*(slavephiint-masterphiint)*timefacrhsfac;
          }

          break;
        }

        default:
        {
          dserror("Kinetic model for scatra-scatra interface coupling not yet implemented!");
          break;
        }
      }
    } // end of loop over scalars
  } // end of loop over integration points

  return;
}

/*----------------------------------------------------------------------*
 | evaluate Robin boundary condition                    schoeder 03/15  |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::CalcRobinBoundary(
    DRT::FaceElement*                 ele,
    Teuchos::ParameterList&           params,
    DRT::Discretization&              discretization,
    std::vector<int>&                 lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    const double                      scalar
    )
{
  // get current condition
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(cond == Teuchos::null)
    dserror("Cannot access condition 'ScatraRobin'");

  double prefac = cond->GetDouble("Prefactor");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over all scalars
  for(int k=0; k<numscal_; ++k)
  {
    for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
    {
      const double fac = DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvalShapeFuncAndIntFac(intpoints,gpid);
      // evaluate overall integration factors
      const double timefac = scatraparamstimint_->TimeFac()*fac;
      const double timefacprefac = timefac * prefac;

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numscal_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          elemat1_epetra(fvi,ui*numscal_+k) += funct_(vi)*funct_(ui)* timefacprefac;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate surface/interface permeability                  Thon 11/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateSurfacePermeability(
        const DRT::Element*        ele,    ///< current boundary element
        const std::vector<double>& ephinp, ///< scalar values at element nodes
        const std::vector<double>& f_wss,  ///< factor for WSS at element nodes
        Epetra_SerialDenseMatrix&  emat,   ///< element-matrix
        Epetra_SerialDenseVector&  erhs,   ///< element-rhs
        const std::vector<int>*    onoff,  ///<flag for dofs to be considered by membrane equations of Kedem and Katchalsky
        const double               perm    ///< surface/interface permeability coefficient
    )
{

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<nen_,1> phinod(true);

  // define vector for wss concentration values at nodes
  LINALG::Matrix<nen_,1> fwssnod(true);

  // loop over all scalars
  for(int k=0;k<numdofpernode_;++k)
  {
    //flag for dofs to be considered by membrane equations of Kedem and Katchalsky
    if((*onoff)[k]==1)
    {
      for (int inode = 0; inode < nen_; ++inode)
      {
        // get scalar values at nodes
        phinod(inode) = ephinp[inode * numdofpernode_ + k];

        // get factor for WSS influence at nodes
        fwssnod(inode) = f_wss[inode * numdofpernode_ + k];
      }

      // loop over all integration points
      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {
        const double fac = EvalShapeFuncAndIntFac(intpoints,iquad,&normal_);

        // integration factor for right-hand side
        double facfac = 0.0;
        if(scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
          facfac = scatraparamstimint_->TimeFac()*fac;
        else
          dserror("EvaluateSurfacePermeability: Requested scheme not yet implemented");

        // scalar at integration point
        const double phi = funct_.Dot(phinod);

        // mean concentration at integration point
        const double facWSS = funct_.Dot(fwssnod);


        // matrix
        for (int vi=0; vi<nen_; ++vi)
        {
          const double vlhs = facfac * facWSS * perm * funct_(vi);
          const int fvi = vi*numdofpernode_+k;

          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) += vlhs*funct_(ui);
          }
        }

        // rhs
        const double vrhs = facfac * facWSS * perm * phi;
//        std::cout<<"k= "<<k<<std::endl;
//        std::cout<<"perm= "<<perm<<std::endl;
//        std::cout<<"phi= "<<phi<<std::endl;
//        std::cout<<"facWSS= "<<facWSS<<std::endl;
//        std::cout<<"vrhs= "<<vrhs<<"\n"<<std::endl;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] -= vrhs*funct_(vi);
        }
      }
    }// if((*onoff)[k]==1)
    // else //in the case of "OFF", a no flux condition is automatically applied
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Kedem-Katchalsky interface                   Thon 11/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::EvaluateKedemKatchalsky(
    const DRT::Element*        ele,         ///< the actual boundary element
    const std::vector<double>& ephinp,      ///< scalar values at element nodes
    const std::vector<double>& epressure,   ///< pressure values at element nodes
    const std::vector<double>& ephibar,     ///< mean concentration values at element nodes
    const std::vector<double>& f_wss,       ///< factor for WSS at element nodes
    Epetra_SerialDenseMatrix&  emat,        ///< element-matrix
    Epetra_SerialDenseVector&  erhs,        ///< element-rhs
    const std::vector<int>*    onoff,       ///<flag for dofs to be considered by membrane equations of Kedem and Katchalsky
    const double               perm,        ///< surface/interface permeability coefficient
    const double               conductivity,///< hydraulic conductivity at interface
    const double               sigma ///< Staverman filtration coefficient at interface
    )
{

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // define vector for scalar values at nodes
  LINALG::Matrix<nen_,1> phinod(true);

  // define vector for pressure values at nodes
  LINALG::Matrix<nen_,1> pnod(true);

  // define vector for mean concentration values at nodes
  LINALG::Matrix<nen_,1> phibarnod(true);

  // define vector for wss concentration values at nodes
  LINALG::Matrix<nen_,1> fwssnod(true);

  // loop over all scalars
  for(int k=0;k<numdofpernode_;++k)   //numdofpernode_//1
  {
    //flag for dofs to be considered by membrane equations of Kedem and Katchalsky
    if((*onoff)[k]==1)
    {
      for (int inode = 0; inode < nen_; ++inode)
      {
        // get scalar values at nodes
        phinod(inode) = ephinp[inode * numdofpernode_ + k];

        // get pressure values at nodes
        pnod(inode) = epressure[inode * numdofpernode_ + k];

        // get mean concentrations at nodes
        phibarnod(inode) = ephibar[inode * numdofpernode_ + k];

        // get factor for WSS influence at nodes
        fwssnod(inode) = f_wss[inode * numdofpernode_ + k];
      }

      // loop over all integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
      {
        const double fac = EvalShapeFuncAndIntFac(intpoints, iquad, &normal_);

        // integration factor
        double facfac = 0.0;
        if(scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
          facfac = scatraparamstimint_->TimeFac() * fac;
        else
          dserror("Kedem-Katchalsky: Requested time integration scheme not yet implemented");

        // scalar at integration point
        const double phi = funct_.Dot(phinod);

        // pressure at integration point
        const double p = funct_.Dot(pnod);

        // mean concentration at integration point
        const double phibar = funct_.Dot(phibarnod);

        // mean concentration at integration point
        const double facWSS = funct_.Dot(fwssnod);


        // matrix
        for (int vi = 0; vi < nen_; ++vi)
        {
          const double vlhs = facfac * facWSS * perm * funct_(vi);

          const int fvi = vi * numdofpernode_ + k;

          for (int ui = 0; ui < nen_; ++ui)
          {
            const int fui = ui * numdofpernode_ + k;

            emat(fvi, fui) += vlhs * funct_(ui);
          }
        }

        // rhs
        const double vrhs = facfac * facWSS * (perm * phi + (1 - sigma) * phibar * conductivity * p);
                        // J_s =f_WSS*[perm*(phi1-phi2)+(1-sigma)*phibar*conductivity*(p1-p2)]
                        // J_s := solute flux through scalar scalar interface
                        // perm:=membrane permeability
                        // sigma:=Staverman filtration coefficient
                        // phibar:= mean concentration within the membrane (for now: simply linear interpolated, but other interpolations also possible)
                        // conductivity:=local hydraulic conductivity of membrane

        for (int vi = 0; vi < nen_; ++vi)
        {
          const int fvi = vi * numdofpernode_ + k;

          erhs[fvi] -= vrhs * funct_(vi);
        }
      }
    } // if((*onoff)[k]==1)
    // else //in the case of "OFF", a no flux condition is automatically applied
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)           gjb 02/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::IntegrateShapeFunctions(
    const DRT::Element*        ele,
    Teuchos::ParameterList&    params,
    Epetra_SerialDenseVector&  elevec1,
    const bool                 addarea
)
{
  // access boundary area variable with its actual value
  double boundaryint = params.get<double>("boundaryint");

  bool outputall = false;
  if(params.isParameter("alldof"))
    outputall = true;

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid);

    // compute integral of shape functions
    for (int node=0;node<nen_;++node)
    {
      for (int k=0; k< numscal_; k++)
      {
        elevec1[node*numdofpernode_+k] += funct_(node) * fac;
      }
      if(outputall==true)
        elevec1[node*numdofpernode_+numdofpernode_-1] += funct_(node) * fac;
    }

    if (addarea)
    {
      // area calculation
      boundaryint += fac;
    }

  } //loop over integration points

  // add contribution to the global value
  params.set<double>("boundaryint",boundaryint);

  return;

} //ScaTraEleBoundaryCalc<distype>::IntegrateShapeFunction


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bdistype,
          DRT::Element::DiscretizationType pdistype>
   void  DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::WeakDirichlet(
     DRT::FaceElement*                  ele,
     Teuchos::ParameterList&            params,
     DRT::Discretization&               discretization,
     Teuchos::RCP<const MAT::Material>  material,
     Epetra_SerialDenseMatrix&          elemat_epetra,
     Epetra_SerialDenseVector&          elevec_epetra)
{
  //------------------------------------------------------------------------
  // Dirichlet boundary condition
  //------------------------------------------------------------------------
  Teuchos::RCP<DRT::Condition> dbc = params.get<Teuchos::RCP<DRT::Condition> >("condition");

  // check of total time
  bool usetime = true;
  const double time = scatraparamstimint_->Time();
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = (*dbc).Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values and spatial functions from condition
  // (assumed to be constant on element boundary)
  const std::vector<double>* val  = (*dbc).Get<std::vector<double> >("val"  );
  const std::vector<int>*    func = (*dbc).Get<std::vector<int> >   ("funct");

  // assign boundary value multiplied by time-curve factor
  double dirichval=(*val)[0]*curvefac;

  // spatial function number
  const int funcnum = (*func)[0];

  //------------------------------------------------------------------------
  // preliminary definitions for (boundary) and parent element and
  // evaluation of nodal values of velocity and scalar based on parent
  // element nodes
  //------------------------------------------------------------------------
  // get the parent element
  DRT::Element* pele = ele->ParentElement();

  // number of spatial dimensions regarding (boundary) element
  static const int bnsd = DRT::UTILS::DisTypeToDim<bdistype>::dim;

  // number of spatial dimensions regarding parent element
  static const int pnsd = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  // number of (boundary) element nodes
  static const int bnen = DRT::UTILS::DisTypeToNumNodePerEle<bdistype>::numNodePerElement;

  // number of parent element nodes
  static const int pnen = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // parent element lm vector
  std::vector<int>  plm ;
  std::vector<int>  plmowner;
  std::vector<int>  plmstride;
  pele->LocationVector(discretization,plm,plmowner,plmstride);

  // get velocity values at parent element nodes
  const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("convective velocity field",Teuchos::null);
  Epetra_SerialDenseVector evel(pnsd*pnen);
  DRT::UTILS::ExtractMyNodeBasedValues(pele,evel,velocity,pnsd);

  // get scalar values at parent element nodes
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

  // extract local values from global vectors for parent element
  std::vector<double> myphinp(plm.size());
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,plm);

  // matrix and vector definition
  LINALG::Matrix<pnsd,pnen>       evelnp;
  std::vector<LINALG::Matrix<pnen,1> > ephinp(numscal_);

  // insert into element arrays
  for (int i=0;i<pnen;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      // split for each tranported scalar, insert into element arrays
      ephinp[k](i,0) = myphinp[k+(i*numdofpernode_)];
    }

    // insert velocity field into element array
    for (int idim=0 ; idim < pnsd; idim++)
    {
      evelnp(idim,i) = evel[idim + i*pnsd];
    }
  }

  //------------------------------------------------------------------------
  // preliminary definitions for integration loop
  //------------------------------------------------------------------------
  // reshape element matrices and vectors and init to zero, construct views
  elemat_epetra.Shape(pnen,pnen);
  elevec_epetra.Size (pnen);
  LINALG::Matrix<pnen,pnen> emat(elemat_epetra.A(),true);
  LINALG::Matrix<pnen,   1> erhs(elevec_epetra.A(),true);

  // (boundary) element local node coordinates
  LINALG::Matrix<pnsd,bnen>  bxyze(true);
  GEO::fillInitialPositionArray<bdistype,pnsd,LINALG::Matrix<pnsd,bnen> >(ele,bxyze);

  // parent element local node coordinates
  LINALG::Matrix<pnsd,pnen>  pxyze(true);
  GEO::fillInitialPositionArray<pdistype,pnsd,LINALG::Matrix<pnsd,pnen> >(pele,pxyze);

  // coordinates of integration points for (boundary) and parent element
  LINALG::Matrix<bnsd,   1>  bxsi(true);
  LINALG::Matrix<pnsd,   1>  pxsi(true);

  // transposed jacobian "dx/ds" and inverse of transposed jacobian "ds/dx"
  // for parent element
  LINALG::Matrix<pnsd,pnsd>  pxjm(true);
  LINALG::Matrix<pnsd,pnsd>  pxji(true);

  // metric tensor for (boundary) element
  LINALG::Matrix<bnsd,bnsd>  bmetrictensor(true);

  // (outward-pointing) unit normal vector to (boundary) element
  LINALG::Matrix<pnsd,   1>  bnormal(true);

  // velocity vector at integration point
  LINALG::Matrix<pnsd,   1>  velint;

  // gradient of scalar value at integration point
  LINALG::Matrix<pnsd,1> gradphi;

  // (boundary) element shape functions, local and global derivatives
  LINALG::Matrix<bnen,   1>  bfunct(true);
  LINALG::Matrix<bnsd,bnen>  bderiv(true);
  LINALG::Matrix<bnsd,bnen>  bderxy(true);

  // parent element shape functions, local and global derivatives
  LINALG::Matrix<pnen,   1>  pfunct(true);
  LINALG::Matrix<pnsd,pnen>  pderiv(true);
  LINALG::Matrix<pnsd,pnen>  pderxy(true);

  //------------------------------------------------------------------------
  // additional matrices and vectors for mixed-hybrid formulation
  //------------------------------------------------------------------------
  // for volume integrals
  LINALG::Matrix<pnsd*pnen,pnsd*pnen> mat_s_q(true);
  LINALG::Matrix<pnsd*pnen,     pnen> mat_s_gradphi(true);

  LINALG::Matrix<pnsd*pnen,        1> vec_s_gradphi(true);

  // for boundary integrals
  LINALG::Matrix<     pnen,pnsd*pnen> mat_w_q_o_n(true);
  LINALG::Matrix<pnsd*pnen,     pnen> mat_s_o_n_phi(true);

  LINALG::Matrix<pnsd*pnen,        1> vec_s_o_n_phi_minus_g(true);

  // inverse matrix
  LINALG::Matrix<pnsd*pnen,pnsd*pnen> inv_s_q(true);

  //------------------------------------------------------------------------
  // check whether Nitsche (default) or mixed-hybrid formulation as well as
  // preliminary definitions and computations for Nitsche stabilization term
  //------------------------------------------------------------------------
  // default is Nitsche formulation
  bool mixhyb = false;

  // stabilization parameter for Nitsche term
  const double nitsche_stab_para = (*dbc).GetDouble("TauBscaling");

  // if stabilization parameter negative: mixed-hybrid formulation
  if (nitsche_stab_para < 0.0) mixhyb = true;

  // pre-factor for adjoint-consistency term:
  // either 1.0 (adjoint-consistent, default) or -1.0 (adjoint-inconsistent)
  double gamma = 1.0;
  const std::string* consistency = (*dbc).Get<std::string>("Choice of gamma parameter");
  if      (*consistency=="adjoint-consistent") gamma = 1.0;
  else if (*consistency=="diffusive-optimal")  gamma = -1.0;
  else dserror("unknown definition for gamma parameter: %s",(*consistency).c_str());

  // use one-point Gauss rule to do calculations at element center
  const DRT::UTILS::IntPointsAndWeights<bnsd> intpoints_tau(SCATRA::DisTypeToStabGaussRule<bdistype>::rule);

  // element surface area (1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double* gpcoord = (intpoints_tau.IP().qxg)[0];
  for (int idim=0;idim<bnsd;idim++)
  {
    bxsi(idim) = gpcoord[idim];
  }
  DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);
  double drs = 0.0;
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);
  const double area = intpoints_tau.IP().qwgt[0]*drs;

  // get number of dimensions for (boundary) element (convert from int to double)
  const double dim = (double) bnsd;

  // computation of characteristic length of (boundary) element
  // (2D: square root of element area, 1D: element length)
  const double h = std::pow(area,(1.0/dim));

  //------------------------------------------------------------------------
  // preliminary computations for integration loop
  //------------------------------------------------------------------------
  // integration points and weights for (boundary) element and parent element
  const DRT::UTILS::IntPointsAndWeights<bnsd> bintpoints(SCATRA::DisTypeToOptGaussRule<bdistype>::rule);

  const DRT::UTILS::IntPointsAndWeights<pnsd> pintpoints(SCATRA::DisTypeToOptGaussRule<pdistype>::rule);

  // transfer integration-point coordinates of (boundary) element to parent element
  Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,pnsd);
  {
    Epetra_SerialDenseMatrix gps(bintpoints.IP().nquad,bnsd);

    for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
    {
      const double* gpcoord = (bintpoints.IP().qxg)[iquad];

      for (int idim=0;idim<bnsd ;idim++)
      {
        gps(iquad,idim) = gpcoord[idim];
      }
    }
    if(pnsd==2)
    {
      DRT::UTILS::BoundaryGPToParentGP2(pqxg,gps,pdistype,bdistype,ele->FaceParentNumber());
    }
    else if (pnsd==3)
    {
      DRT::UTILS::BoundaryGPToParentGP3(pqxg,gps,pdistype,bdistype,ele->FaceParentNumber());
    }

  }

  //------------------------------------------------------------------------
  // integration loop 1: volume integrals (only for mixed-hybrid formulation)
  //------------------------------------------------------------------------
  if (mixhyb)
  {
    for (int iquad=0; iquad<pintpoints.IP().nquad; ++iquad)
    {
      // reference coordinates of integration point from (boundary) element
      const double* gpcoord = (pintpoints.IP().qxg)[iquad];
      for (int idim=0;idim<pnsd;idim++)
      {
        pxsi(idim) = gpcoord[idim];
      }

      // parent element shape functions and local derivatives
      DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

      // Jacobian matrix and determinant of parent element (including check)
      pxjm.MultiplyNT(pderiv,pxyze);
      const double det = pxji.Invert(pxjm);
      if (det < 1E-16) dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

      // compute integration factor
      const double fac = pintpoints.IP().qwgt[iquad]*det;

      // compute global derivatives
      pderxy.Multiply(pxji,pderiv);

      //--------------------------------------------------------------------
      // loop over scalars (not yet implemented for more than one scalar)
      //--------------------------------------------------------------------
      // for(int k=0;k<numdofpernode_;++k)
      int k=0;
      {
        // get viscosity
        if (material->MaterialType() == INPAR::MAT::m_scatra)
        {
          const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

          dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");

          // get constant diffusivity
          diffus_[k] = actmat->Diffusivity();
        }
        else dserror("Material type is not supported");

        // gradient of current scalar value
        gradphi.Multiply(pderxy,ephinp[k]);

        // integration factor for left-hand side
        const double lhsfac = scatraparamstimint_->TimeFac()*fac;

        // integration factor for right-hand side
        double rhsfac = 0.0;
        if(scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac/scatraparamstimint_->AlphaF();
        else if(not scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac*(1.0-scatraparamstimint_->AlphaF())/scatraparamstimint_->AlphaF();
        else if(scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
          rhsfac = lhsfac;

        //--------------------------------------------------------------------
        //  matrix and vector additions due to mixed-hybrid formulation
        //--------------------------------------------------------------------
        /*
                       /         \
                  1   |   h   h  |
              - ----- |  s , q   |
                kappa |          |
                      \          / Omega
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          //const double vlhs = lhsfac*pfunct(vi);
          const double vlhs = lhsfac*(1.0/diffus_[k])*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_s_q(fvi*pnsd+i,fui*pnsd+i) -= vlhs*pfunct(ui);
            }
          }
        }

        /*
                       /                  \
                      |  h         /   h\  |
                    + | s  , grad | phi  | |
                      |            \    /  |
                       \                  / Omega
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          //const double vlhs = lhsfac*diffus_[k]*pfunct(vi);
          const double vlhs = lhsfac*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_s_gradphi(fvi*pnsd+i,fui) += vlhs*pderxy(i,ui);
            }
          }

          //const double vrhs = rhsfac*diffus_[k]*pfunct(vi);
          const double vrhs = rhsfac*pfunct(vi);

          for(int i=0;i<pnsd;++i)
          {
            vec_s_gradphi(fvi*pnsd+i) += vrhs*gradphi(i);
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // integration loop 2: boundary integrals
  //------------------------------------------------------------------------
  for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
  {
    // reference coordinates of integration point from (boundary) element
    const double* gpcoord = (bintpoints.IP().qxg)[iquad];
    for (int idim=0;idim<bnsd;idim++)
    {
      bxsi(idim) = gpcoord[idim];
    }

    // (boundary) element shape functions
    DRT::UTILS::shape_function       <bdistype>(bxsi,bfunct);
    DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);

    // global coordinates of current integration point from (boundary) element
    LINALG::Matrix<pnsd,1> coordgp(true);
    for (int A=0;A<bnen;++A)
    {
      for(int j=0;j<pnsd;++j)
      {
        coordgp(j)+=bxyze(j,A)*bfunct(A);
      }
    }

    // reference coordinates of integration point from parent element
    for (int idim=0;idim<pnsd;idim++)
    {
      pxsi(idim) = pqxg(iquad,idim);
    }

    // parent element shape functions and local derivatives
    DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

    // Jacobian matrix and determinant of parent element (including check)
    pxjm.MultiplyNT(pderiv,pxyze);
    const double det = pxji.Invert(pxjm);
    if (det < 1E-16) dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

    // compute measure tensor for surface element, infinitesimal area element drs
    // and (outward-pointing) unit normal vector
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);

    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    if(DRT::NURBS::IsNurbs(distype))
      bnormal.Scale(normalfac_);

    // compute integration factor
    const double fac = bintpoints.IP().qwgt[iquad]*drs;

    // compute global derivatives
    pderxy.Multiply(pxji,pderiv);

#if 1
    //--------------------------------------------------------------------
    // check whether integration-point coordinates evaluated from
    // (boundary) and parent element match
    //--------------------------------------------------------------------
    LINALG::Matrix<pnsd,1> check(true);
    LINALG::Matrix<pnsd,1> diff(true);

    for (int A=0;A<pnen;++A)
    {
      for(int j=0;j<pnsd;++j)
      {
        check(j)+=pxyze(j,A)*pfunct(A);
      }
    }

    diff=check;
    diff-=coordgp;

    const double norm=diff.Norm2();

    if (norm>1e-9)
    {
      for (int j=0;j<pnsd;++j)
      {
        printf("%12.5e %12.5e\n",check(j),coordgp(j));
      }
      dserror("Gausspoint matching error %12.5e\n",norm);
    }
#endif

    //--------------------------------------------------------------------
    // factor for Dirichlet boundary condition given by spatial function
    //--------------------------------------------------------------------
    double functfac = 1.0;
    if (funcnum > 0)
    {
      // evaluate function at current integration point (important: a 3D position vector is required)
      double coordgp3D[3];
      coordgp3D[0]=0.0;
      coordgp3D[1]=0.0;
      coordgp3D[2]=0.0;
      for (int i=0; i<pnsd;i++)
        coordgp3D[i]=coordgp(i);

      functfac = DRT::Problem::Instance()->Funct(funcnum-1).Evaluate(0,&(coordgp3D[0]),time,NULL);
    }
    else functfac = 1.0;
    dirichval *= functfac;

    //--------------------------------------------------------------------
    // loop over scalars (not yet implemented for more than one scalar)
    //--------------------------------------------------------------------
    // for(int k=0;k<numdofpernode_;++k)
    int k=0;
    {
      // get viscosity
      if (material->MaterialType() == INPAR::MAT::m_scatra)
      {
        const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

        dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");

        // get constant diffusivity
        diffus_[k] = actmat->Diffusivity();
      }
      else dserror("Material type is not supported");

      // get scalar value at integration point
      const double phi = pfunct.Dot(ephinp[k]);

      // integration factor for left-hand side
      const double lhsfac = scatraparamstimint_->TimeFac()*fac;

      // integration factor for right-hand side
      double rhsfac = 0.0;
      if(scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac/scatraparamstimint_->AlphaF();
      else if(not scatraparamstimint_->IsIncremental() and scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac*(1.0-scatraparamstimint_->AlphaF())/scatraparamstimint_->AlphaF();
      else if(scatraparamstimint_->IsIncremental() and not scatraparamstimint_->IsGenAlpha())
        rhsfac = lhsfac;

      if (mixhyb)
      {
        //--------------------------------------------------------------------
        //  matrix and vector additions due to mixed-hybrid formulation
        //--------------------------------------------------------------------
        /*  consistency term
                    /           \
                   |  h   h     |
                 - | w , q  o n |
                   |            |
                   \            / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double vlhs = lhsfac*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_w_q_o_n(fvi,fui*pnsd+i) -= vlhs*pfunct(ui)*bnormal(i);
            }
          }
        }

        /*  adjoint consistency term
                    /                 \
                   |  h          h    |
                 - | s  o n , phi - g |
                   |                  |
                   \                  / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double vlhs = lhsfac*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_s_o_n_phi(fvi*pnsd+i,fui) -= vlhs*pfunct(ui)*bnormal(i);
            }
          }

          for(int i=0;i<pnsd;++i)
          {
            vec_s_o_n_phi_minus_g(fvi*pnsd+i) -= pfunct(vi)*bnormal(i)*(rhsfac*phi - scatraparamstimint_->TimeFac()*fac*dirichval);
          }
        }
      }
      else
      {
        // parameter alpha for Nitsche stabilization term
        const double alpha = nitsche_stab_para*diffus_[k]/h;

        // get velocity at integration point
        velint.Multiply(evelnp,pfunct);

        // normal velocity
        const double normvel = velint.Dot(bnormal);

        // gradient of current scalar value
        gradphi.Multiply(pderxy,ephinp[k]);

        // gradient of current scalar value in normal direction
        const double gradphi_norm = bnormal.Dot(gradphi);

        //--------------------------------------------------------------------
        //  matrix and vector additions due to Nitsche formulation
        //--------------------------------------------------------------------
        /*  consistency term
                    /                           \
                   |  h                  h      |
                 - | w , kappa * grad(phi ) o n |
                   |                            |
                   \                            / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double vlhs = lhsfac*pfunct(vi)*diffus_[k];

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              emat(fvi,fui) -= vlhs*pderxy(i,ui)*bnormal(i);
            }
          }

          const double vrhs = rhsfac*diffus_[k];

          erhs(fvi) += vrhs*pfunct(vi)*gradphi_norm;
        }

        /*  adjoint consistency term, inflow/outflow part
              / --          --                                        \
             |  |         h  |                      h           h     |
           - |  |(a o n) w  +| gamma * kappa *grad(w ) o n , phi - g  |
             |  |            |                                        |
             \  --          --                                        / Gamma_in/out
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          // compute diffusive part
          double prefac = 0.0;
          for(int i=0;i<pnsd;++i)
          {
            prefac += gamma*diffus_[k]*pderxy(i,vi)*bnormal(i);
          }

          // add convective part in case of inflow boundary
          if (normvel<-0.0001) prefac += normvel*pfunct(vi);

          const double vlhs = lhsfac*prefac;

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) -= vlhs*pfunct(ui);
          }

          erhs(fvi) += prefac*(rhsfac*phi - scatraparamstimint_->TimeFac()*fac*dirichval);
        }

        /*  stabilization term
                            /             \
                           |  h     h     |
                 + alpha * | w , phi - g  |
                           |              |
                           \              / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double prefac = alpha*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) += lhsfac*prefac*pfunct(ui);
          }

          erhs(fvi) -= prefac*(rhsfac*phi - scatraparamstimint_->TimeFac()*fac*dirichval);
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // local condensation (only for mixed-hybrid formulation)
  //------------------------------------------------------------------------
  if (mixhyb)
  {
    // matrix inversion of flux-flux block
    inv_s_q = mat_s_q;

    LINALG::FixedSizeSerialDenseSolver<pnsd*pnen,pnsd*pnen> solver;

    solver.SetMatrix(inv_s_q);
    solver.Invert();

    // computation of matrix-matrix and matrix vector products, local assembly
    for (int vi=0; vi<pnen; ++vi)
    {
      for (int ui=0; ui<pnen; ++ui)
      {
        for(int rr=0; rr<pnsd*pnen; ++rr)
        {
          for(int mm=0; mm<pnsd*pnen; ++mm)
          {
            emat(vi,ui)
              -=mat_w_q_o_n(vi,rr)*inv_s_q(rr,mm)*(mat_s_gradphi(mm,ui)+mat_s_o_n_phi(mm,ui));
          }
        }
      }
    }

    for (int vi=0; vi<pnen; ++vi)
    {
      for(int rr=0; rr<pnsd*pnen; ++rr)
      {
        for(int mm=0; mm<pnsd*pnen; ++mm)
        {
          erhs(vi)-= mat_w_q_o_n(vi,rr)*inv_s_q(rr,mm)*(-vec_s_o_n_phi_minus_g(mm)-vec_s_gradphi(mm));
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | calculate boundary conditions for                                    |
 | impl. Characteristic Galerkin (2nd order)                            |
 | time integration just for the reinitialization equation schott 04/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bdistype,
          DRT::Element::DiscretizationType pdistype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::ReinitCharacteristicGalerkinBoundary(
    DRT::FaceElement*                  ele,                  //!< transport element
    Teuchos::ParameterList&            params,               //!< parameter list
    DRT::Discretization&               discretization,       //!< discretization
    Teuchos::RCP<const MAT::Material>  material,             //!< material
    Epetra_SerialDenseMatrix&          elemat_epetra,        //!< ele sysmat
    Epetra_SerialDenseVector&          elevec_epetra         //!< ele rhs
    )
{

  //------------------------------------------------------------------------
  // preliminary definitions for (boundary) and parent element and
  // evaluation of nodal values of velocity and scalar based on parent
  // element nodes
  //------------------------------------------------------------------------
  // get the parent element
  DRT::Element* pele = ele->ParentElement();

  // number of spatial dimensions regarding (boundary) element
  static const int bnsd = DRT::UTILS::DisTypeToDim<bdistype>::dim;

  // number of spatial dimensions regarding parent element
  static const int pnsd = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  // number of (boundary) element nodes
  static const int bnen = DRT::UTILS::DisTypeToNumNodePerEle<bdistype>::numNodePerElement;

  // number of parent element nodes
  static const int pnen = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // parent element lm vector
  std::vector<int>  plm ;
  std::vector<int>  plmowner;
  std::vector<int>  plmstride;
  pele->LocationVector(discretization,plm,plmowner,plmstride);

  // get scalar values at parent element nodes
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");
  Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
  if (phinp==Teuchos::null) dserror("Cannot get state vector 'phin'");

  // extract local values from global vectors for parent element
  std::vector<double> myphinp(plm.size());
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,plm);

  std::vector<double> myphin(plm.size());
  DRT::UTILS::ExtractMyValues(*phin,myphin,plm);

  //    // matrix and vector definition
  //    LINALG::Matrix<pnsd,pnen>       evelnp;
  std::vector<LINALG::Matrix<pnen,1> > ephinp(numscal_);
  std::vector<LINALG::Matrix<pnen,1> > ephin(numscal_);

  // insert into element arrays
  for (int i=0;i<pnen;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      // split for each tranported scalar, insert into element arrays
      ephinp[k](i,0) = myphinp[k+(i*numdofpernode_)];
      ephin[k](i,0)  = myphin[k+(i*numdofpernode_)];
    }
  }

  //------------------------------------------------------------------------
  // preliminary definitions for integration loop
  //------------------------------------------------------------------------
  // reshape element matrices and vectors and init to zero, construct views
  elemat_epetra.Shape(pnen,pnen);
  elevec_epetra.Size (pnen);
  LINALG::Matrix<pnen,pnen> emat(elemat_epetra.A(),true);
  LINALG::Matrix<pnen,   1> erhs(elevec_epetra.A(),true);

  // (boundary) element local node coordinates
  LINALG::Matrix<pnsd,bnen>  bxyze(true);
  GEO::fillInitialPositionArray<bdistype,pnsd,LINALG::Matrix<pnsd,bnen> >(ele,bxyze);

  // parent element local node coordinates
  LINALG::Matrix<pnsd,pnen>  pxyze(true);
  GEO::fillInitialPositionArray<pdistype,pnsd,LINALG::Matrix<pnsd,pnen> >(pele,pxyze);

  // coordinates of integration points for (boundary) and parent element
  LINALG::Matrix<bnsd,   1>  bxsi(true);
  LINALG::Matrix<pnsd,   1>  pxsi(true);

  // transposed jacobian "dx/ds" and inverse of transposed jacobian "ds/dx"
  // for parent element
  LINALG::Matrix<pnsd,pnsd>  pxjm(true);
  LINALG::Matrix<pnsd,pnsd>  pxji(true);

  // metric tensor for (boundary) element
  LINALG::Matrix<bnsd,bnsd>  bmetrictensor(true);

  // (outward-pointing) unit normal vector to (boundary) element
  LINALG::Matrix<pnsd,   1>  bnormal(true);

  // velocity vector at integration point
  LINALG::Matrix<pnsd,   1>  velint;

  // gradient of scalar value at integration point
  LINALG::Matrix<pnsd,1> gradphi;

  // (boundary) element shape functions, local and global derivatives
  LINALG::Matrix<bnen,   1>  bfunct(true);
  LINALG::Matrix<bnsd,bnen>  bderiv(true);
  LINALG::Matrix<bnsd,bnen>  bderxy(true);

  // parent element shape functions, local and global derivatives
  LINALG::Matrix<pnen,   1>  pfunct(true);
  LINALG::Matrix<pnsd,pnen>  pderiv(true);
  LINALG::Matrix<pnsd,pnen>  pderxy(true);


  // use one-point Gauss rule to do calculations at element center
  const DRT::UTILS::IntPointsAndWeights<bnsd> intpoints_tau(SCATRA::DisTypeToStabGaussRule<bdistype>::rule);

  // element surface area (1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double* gpcoord = (intpoints_tau.IP().qxg)[0];
  for (int idim=0;idim<bnsd;idim++)
  {
    bxsi(idim) = gpcoord[idim];
  }
  DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);
  double drs = 0.0;
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);

  //------------------------------------------------------------------------
  // preliminary computations for integration loop
  //------------------------------------------------------------------------
  // integration points and weights for (boundary) element and parent element
  const DRT::UTILS::IntPointsAndWeights<bnsd> bintpoints(SCATRA::DisTypeToOptGaussRule<bdistype>::rule);

  const DRT::UTILS::IntPointsAndWeights<pnsd> pintpoints(SCATRA::DisTypeToOptGaussRule<pdistype>::rule);

  // transfer integration-point coordinates of (boundary) element to parent element
  Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,pnsd);
  {
    Epetra_SerialDenseMatrix gps(bintpoints.IP().nquad,bnsd);

    for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
    {
      const double* gpcoord = (bintpoints.IP().qxg)[iquad];

      for (int idim=0;idim<bnsd ;idim++)
      {
        gps(iquad,idim) = gpcoord[idim];
      }
    }
    if(pnsd==2)
    {
      DRT::UTILS::BoundaryGPToParentGP2(pqxg,gps,pdistype,bdistype,ele->FaceParentNumber());
    }
    else if (pnsd==3)
    {
      DRT::UTILS::BoundaryGPToParentGP3(pqxg,gps,pdistype,bdistype,ele->FaceParentNumber());
    }
  }


  const double reinit_pseudo_timestepsize_factor = params.get<double>("pseudotimestepsize_factor");

  const double meshsize = getEleDiameter<pdistype>(pxyze);

  const double pseudo_timestep_size = meshsize * reinit_pseudo_timestepsize_factor;

  //------------------------------------------------------------------------
  // integration loop: boundary integrals
  //------------------------------------------------------------------------
  for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
  {
    // reference coordinates of integration point from (boundary) element
    const double* gpcoord = (bintpoints.IP().qxg)[iquad];
    for (int idim=0;idim<bnsd;idim++)
    {
      bxsi(idim) = gpcoord[idim];
    }

    // (boundary) element shape functions
    DRT::UTILS::shape_function       <bdistype>(bxsi,bfunct);
    DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);

    // global coordinates of current integration point from (boundary) element
    LINALG::Matrix<pnsd,1> coordgp(true);
    for (int A=0;A<bnen;++A)
    {
      for(int j=0;j<pnsd;++j)
      {
        coordgp(j)+=bxyze(j,A)*bfunct(A);
      }
    }

    // reference coordinates of integration point from parent element
    for (int idim=0;idim<pnsd;idim++)
    {
      pxsi(idim) = pqxg(iquad,idim);
    }

    // parent element shape functions and local derivatives
    DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

    // Jacobian matrix and determinant of parent element (including check)
    pxjm.MultiplyNT(pderiv,pxyze);
    const double det = pxji.Invert(pxjm);
    if (det < 1E-16) dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

    // compute measure tensor for surface element, infinitesimal area element drs
    // and (outward-pointing) unit normal vector
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);

    // for nurbs elements the normal vector must be scaled with a special orientation factor!!
    if(DRT::NURBS::IsNurbs(distype))
      bnormal.Scale(normalfac_);

    // compute integration factor
    const double fac_surface = bintpoints.IP().qwgt[iquad]*drs;

    // compute global derivatives
    pderxy.Multiply(pxji,pderiv);

    //--------------------------------------------------------------------
    // loop over scalars (not yet implemented for more than one scalar)
    //--------------------------------------------------------------------
    for(int dofindex=0;dofindex<numdofpernode_;++dofindex)
    {
      //----------  --------------      |                    |
      //  mat              -1/4* dtau^2 | w, n*grad(D(psi) ) |
      //--------------------------      |                    |

      LINALG::Matrix<1,pnen> derxy_normal;
      derxy_normal.Clear();
      derxy_normal.MultiplyTN(bnormal,pderxy);

      for (int vi=0; vi<pnen; ++vi)
      {
        const int fvi = vi*numdofpernode_+dofindex;

        for (int ui=0; ui<pnen; ++ui)
        {
          const int fui = ui*numdofpernode_+dofindex;

          emat(fvi,fui) -= pfunct(vi)* (fac_surface*pseudo_timestep_size*pseudo_timestep_size/4.0) * derxy_normal(0,ui);
        }
      }

      //----------  --------------      |              m     |
      //  rhs               0.5* dtau^2 | w, n*grad(psi )    |
      //--------------------------      |                    |

      // update grad_dist_n
      LINALG::Matrix<pnsd,1> grad_dist_n(true);
      grad_dist_n.Multiply(pderxy,ephin[dofindex]);

      LINALG::Matrix<1,1> grad_dist_n_normal(true);
      grad_dist_n_normal.MultiplyTN(bnormal,grad_dist_n);

      for (int vi=0; vi<pnen; ++vi)
      {
        const int fvi = vi*numdofpernode_+dofindex;

        erhs(fvi) += pfunct(vi)*pseudo_timestep_size*pseudo_timestep_size*fac_surface/2.0 * grad_dist_n_normal(0,0);
      }


      //                    |              m+1     m  |
      //    1/4*delta_tau^2 | w, n*grad(psi   - psi ) |
      //                    |              i          |
      // update grad_dist_n
      LINALG::Matrix<pnsd,1> grad_dist_npi(true);
      grad_dist_npi.Multiply(pderxy,ephinp[dofindex]);

      LINALG::Matrix<1,1> grad_dist_npi_normal;
      grad_dist_npi_normal.Clear();
      grad_dist_npi_normal.MultiplyTN(bnormal,grad_dist_npi);

      double Grad_Dpsi_normal = grad_dist_npi_normal(0,0) - grad_dist_n_normal(0,0);


      for (int vi=0; vi<pnen; ++vi)
      {
        const int fvi = vi*numdofpernode_+dofindex;

        erhs(fvi) += pfunct(vi)*Grad_Dpsi_normal*fac_surface*pseudo_timestep_size*pseudo_timestep_size/4.0;
      }

    } // loop over scalars
  } // loop over integration points

  return;
}


/*----------------------------------------------------------------------*
 |  Factor of WSS dependent interface transport           hemmler 07/14 |
 |  Calculated as in Calvez, V., "Mathematical and numerical modeling of early atherosclerotic lesions",ESAIM: Proceedings. Vol. 30. EDP Sciences, 2010
 *----------------------------------------------------------------------*/

template <DRT::Element::DiscretizationType distype>
std::vector<double> DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype>::WSSinfluence(std::vector<double>  ewss, const bool wssonoff, const std::vector<double>* coeffs)
{
  std::vector<double> f_wss(ewss.size());

  for(int i=0; i<(int)ewss.size();i++)
  {
    if (wssonoff)
      f_wss[i] =((*coeffs)[0])*log10(1+((*coeffs)[1])/(ewss[i]+(*coeffs)[2])); //empirical function (log law) to account for influence of WSS;
    else //no WSS influence
      f_wss[i] = 1.0;
  }

  return f_wss;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalc<DRT::Element::nurbs9>;

