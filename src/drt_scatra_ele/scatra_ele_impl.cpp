/*----------------------------------------------------------------------*/
/*!
  \file scatra_ele_impl.cpp

  \brief Internal implementation of scalar transport elements

  <pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
  </pre>
*/
/*----------------------------------------------------------------------*/
#include "scatra_ele_action.H"
#include "scatra_ele_boundary_impl.H"

#include "scatra_ele_impl.H"

#include "../drt_lib/drt_globalproblem.H"  // for time curve in body force
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"
//#include <Teuchos_StandardParameterEntryValidators.hpp>  // included by inpar files
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_turbulence.H"

#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/myocard.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_spec.H"
#include "../drt_mat/arrhenius_temp.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/ion.H"
#include "../drt_mat/biofilm.H"
#include "../drt_mat/fourieriso.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/yoghurt.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/structporo.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/elchmat.H"
#include "../drt_mat/elchphase.H"

// activate debug screen output
//#define PRINT_ELCH_DEBUG
// use effective diffusion coefficient for stabilization
#define ACTIVATEBINARYELECTROLYTE


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraImplInterface* DRT::ELEMENTS::ScaTraImplInterface::Impl(
  const DRT::Element* ele,
  const enum INPAR::SCATRA::ScaTraType scatratype,
  const bool tg_or_reinit
  )
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));
  int numscal = numdofpernode;
  if (SCATRA::IsElchProblem(scatratype))
  {
    numscal -= 1;

    // get the material of the first element
    // we assume here, that the material is equal for all elements in this discretization
    Teuchos::RCP<MAT::Material> material = ele->Material();
    if (material->MaterialType() == INPAR::MAT::m_elchmat)
    {
      const MAT::ElchMat* actmat = static_cast<const MAT::ElchMat*>(material.get());

      if (actmat->Current())
        numscal -= DRT::UTILS::getDimension(ele->Shape());
    }
  }

  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    if (tg_or_reinit)
      return ReInitImpl<DRT::Element::hex8>::Instance(numdofpernode,numscal);
    else
      return ScaTraImpl<DRT::Element::hex8>::Instance(numdofpernode,numscal);
  } /*
  case DRT::Element::hex20:
  {
    return ScaTraImpl<DRT::Element::hex20>::Instance(numdofpernode,numscal);
  } */
  case DRT::Element::hex27:
  {
    return ScaTraImpl<DRT::Element::hex27>::Instance(numdofpernode,numscal);
  } /*
  case DRT::Element::nurbs8:
  {
    return ScaTraImpl<DRT::Element::nurbs8>::Instance(numdofpernode,numscal);
  } */
  case DRT::Element::nurbs27:
  {
    return ScaTraImpl<DRT::Element::nurbs27>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::tet4:
  {
    return ScaTraImpl<DRT::Element::tet4>::Instance(numdofpernode,numscal);
  } /*
  case DRT::Element::tet10:
  {
    return ScaTraImpl<DRT::Element::tet10>::Instance(numdofpernode,numscal);
  }*/
  case DRT::Element::wedge6:
  {
    return ScaTraImpl<DRT::Element::wedge6>::Instance(numdofpernode,numscal);
  } /*
  case DRT::Element::wedge15:
  {
    return ScaTraImpl<DRT::Element::wedge15>::Instance(numdofpernode,numscal);
  } */
  case DRT::Element::pyramid5:
  {
    return ScaTraImpl<DRT::Element::pyramid5>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::quad4:
  {
    return ScaTraImpl<DRT::Element::quad4>::Instance(numdofpernode,numscal);
  } /*
  case DRT::Element::quad8:
  {
    return ScaTraImpl<DRT::Element::quad8>::Instance(numdofpernode,numscal);
  } */
  case DRT::Element::quad9:
  {
    return ScaTraImpl<DRT::Element::quad9>::Instance(numdofpernode,numscal);
  } /*
  case DRT::Element::nurbs4:
  {
    return ScaTraImpl<DRT::Element::nurbs4>::Instance(numdofpernode,numscal);
  } */
  case DRT::Element::nurbs9:
  {
    return ScaTraImpl<DRT::Element::nurbs9>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::tri3:
  {
    return ScaTraImpl<DRT::Element::tri3>::Instance(numdofpernode,numscal);
  } /*
  case DRT::Element::tri6:
  {
    return ScaTraImpl<DRT::Element::tri6>::Instance(numdofpernode,numscal);
  }*/
  case DRT::Element::line2:
  {
    return ScaTraImpl<DRT::Element::line2>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::line3:
  {
    return ScaTraImpl<DRT::Element::line3>::Instance(numdofpernode,numscal);
  }
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(ele->Shape()).c_str());
    break;
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraImpl<distype> * DRT::ELEMENTS::ScaTraImpl<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create
  )
{
  static ScaTraImpl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraImpl<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0,0,false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraImpl<distype>::ScaTraImpl(const int numdofpernode, const int numscal)
  : numdofpernode_(numdofpernode),
    numscal_(numscal),
    is_elch_((numdofpernode_ - numscal_) >= 1),  // bool set implicitely
    is_ale_(false),           // bool set
    is_reactive_(false),      // bool set
    is_coupled_(false),            //bool set
    diffreastafac_(0.0),       // set double (SUPG)
    is_anisotropic_(false),   // bool set
    is_stationary_(false),    // bool set
    is_genalpha_(false),      // bool set
    is_incremental_(false),   // bool set
    is_conservative_(false),  // bool set
    sgvel_(false),            // bool set
    betterconsistency_(false), // bool set
    migrationintau_(true),     // bool set
    migrationstab_(true),      // bool set
    migrationinresidual_(true),// bool set
    update_mat_(false),        // bool set
    // whichtau_ not initialized
    turbmodel_(INPAR::FLUID::no_model), // enum initialized
    mfs_conservative_(false), // set false
    HSTCConds_(0,NULL),   //vector containing homogeneous coupling conditions
    phi_(numscal_,0.0),   // size of vector + initialized to zero
    sgphi_(numscal_,0.0),   // size of vector + initialized to zero
    mfssgphi_(numscal_,0.0),// size of vector + initialized to zero
    gradphi_(true),     // initialized to zero
    fsgradphi_(true),   // initialized to zero
    mfsggradphi_(true), // initialized to zero
    ephin_(numscal_),   // size of vector
    ephinp_(numscal_),  // size of vector
    ephiam_(numscal_),  // size of vector
    hist_(numscal_,0.0),// size of vector + initialized to zero
    ehist_(numscal_),   // size of vector
    ephi0_Reinit_Reference_(numscal_),// size of vector
    ephi0_penalty_(numscal_),         // size of vector
    fsphinp_(numscal_), // size of vector
    conint_(numscal_,0.0),  // size of vector + initialized to zero
    epotnp_(true),      // initialized to zero
    emagnetnp_(true),   // initialized to zero
    gradpot_(true),     // initialized to zero
    evelnp_(true),      // initialized to zero
    econvelnp_(true),   // initialized to zero
    efsvel_(true),      // initialized to zero
    eaccnp_(true),      // initialized to zero
    edispnp_(true),     // initialized to zero
    velint_(true),      // initialized to zero
    convelint_(true),   // initialized to zero
    sgvelint_(true),    // initialized to zero
    fsvelint_(true),    // initialized to zero
    mfsgvelint_(true),  // initialized to zero
    migvelint_(true),   // initialized to zero
    conv_(true),        // initialized to zero
    sgconv_(true),      // initialized to zero
    vdiv_(0.0),         // set double
    mfsvdiv_(0.0),      // set double
    eprenp_(true),      // initialized to zero
    thermpressnp_(0.0), // set double
    thermpressam_(0.0), // set double
    thermpressdt_(0.0), // set double
    densn_(numscal_,1.0),        // size of vector + initialized to zero
    densnp_(numscal_,1.0),       // size of vector + initialized to zero
    densam_(numscal_,1.0),       // size of vector + initialized to zero
    densgradfac_(numscal_,0.0),  // size of vector + initialized to zero
    diffus_(numscal_,0.0),       // size of vector + initialized to zero
    diffus3_(numscal_),          // size of vector
    sgdiff_(numscal_,0.0),       // size of vector + initialized to zero
    reacterm_(numscal_,0.0),     // size of vector + initialized to zero
    reacoeff_(numscal_,0.0),     // size of vector + initialized to zero
    reacoeffderiv_(numscal_,0.0),// size of vector + initialized to zero
    reacoeffderivmatrix_(numscal_,std::vector<double>(numscal_,0.0 )),// size of matrix + initialized to zero
    valence_(numscal_,0.0),      // size of vector + initialized to zero
    diffusvalence_(numscal_,0.0),// size of vector + initialized to zero
    shc_(0.0),      // set double
    visc_(0.0),     // set double
    diff_(true),    // initialized to zero
    migconv_(true), // initialized to zero
    migrea_(true),  // initialized to zero
    xsi_(true),     // initialized to zero
    xyze_(true),    // initialized to zero
    xyze3D_(true),  // initialized to zero
    lengthLine3D_(0),  // initialized to zero
    funct_(true),   // initialized to zero
    deriv_(true),   // initialized to zero
    deriv2_(true),  // initialized to zero
    derxy_(true),   // initialized to zero
    derxy2_(true),  // initialized to zero
    xjm_(true),     // initialized to zero
    xij_(true),     // initialized to zero
    xder2_(true),   // initialized to zero
    laplace_(true), // initialized to zero
    rhs_(numdofpernode_,0.0),       // size of vector + initialized to zero
    reatemprhs_(numdofpernode_,0.0),// size of vector + initialized to zero
    bodyforce_(numdofpernode_), // size of vector
    scatrares_(numscal_,0.0),  // size of vector + initialized to zero
    conv_phi_(numscal_,0.0),   // size of vector + initialized to zero
    diff_phi_(numscal_,0.0),   // size of vector + initialized to zero
    rea_phi_(numscal_,0.0),    // size of vector + initialized to zero
    tau_(numscal_,0.0),        // size of vector + initialized to zero
    tauderpot_(numscal_),      // size of vector
    efluxreconstr_(numscal_),  // size of vector
    weights_(true),      // initialized to zero
    myknots_(nsd_),       // size of vector
    eid_(0),
    diffcond_(false),
    cursolvar_(false),
    chemdiffcoupltransp_(true),
    chemdiffcouplcurr_(true),
    constparams_(true),
    newman_(false),
    diffbased_(true),
    gradphicoupling_(numscal_),
    curdiv_(0.0),
    trans_(numscal_,0.0),
    transelim_(0.0),
    transderiv_(numscal_,std::vector<double>(numscal_,0.0 )),
    cond_(1,0.0),
    condderiv_(numscal_,0.0),
    therm_(1,0.0),
    thermderiv_(numscal_,1.0),
    diffusderiv_(numscal_,0.0),
    diffuselimderiv_(numscal_,0.0),
    diffuselim_(0.0),
    eps_(1,1.0),
    tort_(1,1.0),
    epstort_(1,1.0),
    a_(0.0),
    b_(0.0),
    c_(0.0),
    ecurnp_(true),
    curint_(true),
    equpot_(INPAR::ELCH::equpot_enc),
    frt_(0.0)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraImpl<distype>::Evaluate(
  DRT::Element*              ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra
  )
{
    //TODO: SCATRA_ELE_CLEANING: CARDIO
  // --------mandatory are performed here at first ------------
  // get node coordinates (we do this for all actions!)
  if(nsd_ == 1 && nen_ == 2){ // Rotates line element in order to perform computations in 3D
    GEO::fillInitialPositionArray<distype,3,LINALG::Matrix<3,nen_> >(ele,xyze3D_);
    lengthLine3D_ = sqrt( (xyze3D_(0,0)-xyze3D_(0,1))*(xyze3D_(0,0)-xyze3D_(0,1)) + (xyze3D_(1,0)-xyze3D_(1,1))*(xyze3D_(1,0)-xyze3D_(1,1)) + (xyze3D_(2,0)-xyze3D_(2,1))*(xyze3D_(2,0)-xyze3D_(2,1)));
    xyze_(0,0)=0.0;
    xyze_(0,1)=lengthLine3D_;
  }  else{
    GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);
  }





  // set diffusion tensor to identity
  for(int k=0; k<numscal_; k++)
  {
    for(int i=0; i<nsd_; i++){ for(int j=0; j<nsd_; j++){ diffus3_[k](i,j)=0; }}
    for(int i=0; i<nsd_; i++){ diffus3_[k](i,i) = 1; }
  }


  
  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");
  switch (action)
  {
  case SCATRA::calc_mat_and_rhs:
  {
    // set flag for including reactive terms to false initially
    // flag will be set to true below when reactive material is included
    is_reactive_ = false;

    // flag will be set to true and data is read and saved below when homogeneous scatra coupling volume condition is included
    is_coupled_ = IsCoupledAndRead(discretization);



    if (scatratype == INPAR::SCATRA::scatratype_loma)
    {
      thermpressnp_ = params.get<double>("thermodynamic pressure");
      thermpressdt_ = params.get<double>("time derivative of thermodynamic pressure");
      if (is_genalpha_)
        thermpressam_ = params.get<double>("thermodynamic pressure at n+alpha_M");

      // update material with subgrid-scale scalar
      update_mat_ = params.get<bool>("update material", false);
    }

    if (scatratype == INPAR::SCATRA::scatratype_loma or
        scatratype == INPAR::SCATRA::scatratype_turbpassivesca)
    {

      // as the scalar field is constant in the turbulent inflow section
      // we do not need any turbulence model
      if (params.get<bool>("turbulent inflow",false))
      {
        if (SCATRA::InflowElement(ele))
          turbmodel_ = INPAR::FLUID::no_model;
      }
    }




    if ((scatratype == INPAR::SCATRA::scatratype_loma) and is_genalpha_)
    {
      // extract additional local values from global vector
      Teuchos::RCP<const Epetra_Vector> phiam = discretization.GetState("phiam");
      if (phiam==Teuchos::null) dserror("Cannot get state vector 'phiam'");
      std::vector<double> myphiam(lm.size());
      DRT::UTILS::ExtractMyValues(*phiam,myphiam,lm);

      // fill element array
      for (int i=0;i<nen_;++i)
      {
        for (int k = 0; k< numscal_; ++k)
        {
          // split for each transported scalar, insert into element arrays
          ephiam_[k](i,0) = myphiam[k+(i*numdofpernode_)];
        }
      } // for i
    }



    double frt(0.0);
    if (is_elch_)
    {
      // set specific parameter used in diffusion conduction formulation
      // this method need to be located inside ELCH
      DiffCondParams(ele, params);

      // safety check - only stabilization of SUPG-type available
      if ((stabinp !=INPAR::SCATRA::stabtype_no_stabilization) and (stabinp !=INPAR::SCATRA::stabtype_SUPG))
        dserror("Only SUPG-type stabilization available for ELCH.");

      // get values for el. potential at element nodes
      for (int i=0;i<nen_;++i)
      {
        epotnp_(i) = myphinp[i*numdofpernode_+numscal_];
      }
      // get parameter F/RT needed for ELCH ;-)
      frt = params.get<double>("frt");

      const INPAR::SCATRA::Consistency consistency
        = DRT::INPUT::IntegralValue<INPAR::SCATRA::Consistency>(stablist,"CONSISTENCY");
      betterconsistency_=(consistency==INPAR::SCATRA::consistency_l2_projection_lumped);

      for (int k=0; k < numscal_; k++)
      {
        if (betterconsistency_)
        {
          std::ostringstream temp;
          temp << k;
          std::string name = "flux_phi_"+temp.str();
          // try to get the pointer to the entry (and check if type is RCP<Epetra_MultiVector>)
          Teuchos::RCP<Epetra_MultiVector>* f = params.getPtr< Teuchos::RCP<Epetra_MultiVector> >(name);
          if (f!= NULL) // field has been set and is not of type Teuchos::null
          {
            DRT::UTILS::ExtractMyNodeBasedValues(ele,efluxreconstr_[k],*f,nsd_);
          }
          else
            dserror("Could not extract values of flux approximation");
        }
        else
          efluxreconstr_[k].Clear();
      }

      // get magnetic field at nodes (if available)
      // try to get the pointer to the entry (and check if type is RCP<Epetra_MultiVector>)
      Teuchos::RCP<Epetra_MultiVector>* b = params.getPtr< Teuchos::RCP<Epetra_MultiVector> >("magnetic field");
      if (b!= NULL) // magnetic field has been set and is not of type Teuchos::null
        DRT::UTILS::ExtractMyNodeBasedValues(ele,emagnetnp_,*b,nsd_);
      else
        emagnetnp_.Clear();

      if(diffcond_ == true)
      {
        // safety check - no stabilization for diffusion-conduction formulation
        if(stabinp !=INPAR::SCATRA::stabtype_no_stabilization)
          dserror("No stabilization available for the diffusion-conduction formulation.");

        if(whichtau_ != INPAR::SCATRA::tau_zero)
          dserror("No stabilization available for the diffusion-conduction formulation.");

        if((tau_gp_==false) or (mat_gp_==false))
          dserror("Evaluation of material (and stabilization parameter) available only at Gausspoints.");
      }

      if (cursolvar_ == true)
      {
        // get values for current at element nodes
        for (int i=0;i<nen_;++i)
        {
          for(int idim=0; idim<nsd_; ++idim)
          {
            //current is located after potential
            ecurnp_(idim,i) = myphinp[i*numdofpernode_+(numscal_+1)+idim];
          }
        }
      }
      else
        ecurnp_.Clear();
    }
    else
    {
      epotnp_.Clear();
      emagnetnp_.Clear();
    }


    break;
  }



  case SCATRA::get_material_parameters:
  {
    // get the material
    Teuchos::RCP<MAT::Material> material = ele->Material();

    if (material->MaterialType() == INPAR::MAT::m_sutherland)
    {
      const Teuchos::RCP<const MAT::Sutherland>& actmat
        = Teuchos::rcp_dynamic_cast<const MAT::Sutherland>(material);
      params.set("thermodynamic pressure",actmat->ThermPress());
    }
    else params.set("thermodynamic pressure",0.0);

    break;
  }



  case SCATRA::time_update_material:
  {
    std::vector<Teuchos::RCP<MAT::Myocard> > updatemat;
    updatemat.reserve(numscal_);

    // access the general material
    Teuchos::RCP<MAT::Material> material = ele->Material();

    // first, determine the materials which need a time update, i.e. myocard materials
    if (material->MaterialType() == INPAR::MAT::m_matlist)
    {
      const Teuchos::RCP<MAT::MatList> actmat
      = Teuchos::rcp_dynamic_cast<MAT::MatList>(material);
      if (actmat->NumMat() < numscal_) dserror("Not enough materials in MatList.");

      for (int k = 0;k<numscal_;++k)
      {
        const int matid = actmat->MatID(k);
        Teuchos::RCP<MAT::Material> singlemat = actmat->MaterialById(matid);

        if (singlemat->MaterialType() == INPAR::MAT::m_myocard)
        {
          // reference to Teuchos::rcp not possible here, since the material is required to be
          // not const for this application
          updatemat.push_back(Teuchos::rcp_dynamic_cast<MAT::Myocard>(singlemat));
        }
      }
    }
    if (material->MaterialType() == INPAR::MAT::m_myocard)
    {      // reference to Teuchos::rcp not possible here, since the material is required to be
      // not const for this application
      updatemat.push_back(Teuchos::rcp_dynamic_cast<MAT::Myocard>(material));
    }

    if (updatemat.size()>0) // found at least one material to be updated
    {
      // all materials in the matlist should be of the same kind
      if (updatemat.size()!= (unsigned) numscal_) dserror("Not allowed");

      // get time-step length
      const double dt   = params.get<double>("time-step length");

      // extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp==Teuchos::null)
        dserror("Cannot get state vector 'phinp'");
      std::vector<double> myphinp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

      // fill all element arrays
      for (int i=0;i<nen_;++i)
      {
        for (int k = 0; k< numscal_; ++k)
        {
          // split for each transported scalar, insert into element arrays
          ephinp_[k](i,0) = myphinp[k+(i*numdofpernode_)];
        }
      }

      // use one-point Gauss rule to do calculations at the element center
      DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

      EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0,ele->Id());

      for (unsigned i=0;i<updatemat.size();i++)
      {
        const double csnp = funct_.Dot(ephinp_[i]); // be careful, we assume k==i here
        updatemat[i]->Update(csnp, dt);
      }
    }
    break;
  }
  case SCATRA::get_material_internal_state:
  {

    // NOTE: add integral values only for elements which are NOT ghosted!
    if (ele->Owner() == discretization.Comm().MyPID())
    {

      // access the general material
      Teuchos::RCP<MAT::Material> material = ele->Material();
      Teuchos::RCP<Epetra_MultiVector> material_internal_state = params.get< Teuchos::RCP<Epetra_MultiVector> >("material_internal_state");

      if( material->MaterialType() == INPAR::MAT::m_myocard)
      {
        Teuchos::RCP<MAT::Myocard> material = Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
        for (int k = 0; k< material_internal_state->NumVectors(); ++k)
        {
          int err = material_internal_state->ReplaceGlobalValue(ele->Id(),k, material->GetInternalState(k));
          if(err != 0) dserror("%i",err);
        }
      }
      params.set< Teuchos::RCP<Epetra_MultiVector> >("material_internal_state", material_internal_state);
    }

    break;

  }

  case SCATRA::set_material_internal_state:
  {
    // NOTE: add integral values only for elements which are NOT ghosted!
    if (ele->Owner() == discretization.Comm().MyPID())
    {
      // access the general material
      Teuchos::RCP<MAT::Material> material = ele->Material();
      Teuchos::RCP<Epetra_Vector> material_internal_state_component = params.get< Teuchos::RCP<Epetra_Vector> >("material_internal_state_component");

      if( material->MaterialType() == INPAR::MAT::m_myocard)
      {
        Teuchos::RCP<MAT::Myocard> material = Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
        int k =  params.get< int >("k");
        material->SetInternalState(k,(*material_internal_state_component)[ele->Id()]);
      }

    }

    break;

  }

  case SCATRA::get_material_ionic_currents:
  {
    // NOTE: add integral values only for elements which are NOT ghosted!
    if (ele->Owner() == discretization.Comm().MyPID())
    {
      // access the general material
      Teuchos::RCP<MAT::Material> material = ele->Material();
      Teuchos::RCP<Epetra_MultiVector> material_ionic_currents = params.get< Teuchos::RCP<Epetra_MultiVector> >("material_ionic_currents");

      if( material->MaterialType() == INPAR::MAT::m_myocard)
      {
        Teuchos::RCP<MAT::Myocard> material = Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
        for (int k = 0; k< material_ionic_currents->NumVectors(); ++k)
        {
          int err = material_ionic_currents->ReplaceGlobalValue(ele->Id(),k, material->GetIonicCurrents(k));
          if(err != 0) dserror("%i",err);
        }
      }
      params.set< Teuchos::RCP<Epetra_MultiVector> >("material_ionic_currents", material_ionic_currents);
    }

    break;

  }









  case SCATRA::calc_error:
  {
    // check if length suffices
    if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

    // set specific parameter used in diffusion conduction formulation
    // this method need to be located inside ELCH
    DiffCondParams(ele, params);

    // need current solution
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    std::vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // fill element arrays
    for (int i=0;i<nen_;++i)
    {
      // split for each transported scalar, insert into element arrays
      for (int k = 0; k< numscal_; ++k)
      {
        ephinp_[k](i) = myphinp[k+(i*numdofpernode_)];
      }
      // get values for el. potential at element nodes
      epotnp_(i) = myphinp[i*numdofpernode_+numscal_];
    } // for i

    CalErrorComparedToAnalytSolution(
      ele,
      scatratype,
      params,
      elevec1_epetra);

    break;
  }
  case SCATRA::calc_elch_conductivity:
  {
    if(is_elch_)
    {
      // set specific parameter used in diffusion conduction formulation
      // this method need to be located inside ELCH
      DiffCondParams(ele, params);

      // calculate conductivity of electrolyte solution
      const double frt = params.get<double>("frt");
      // extract local values from the global vector
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      std::vector<double> myphinp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

      // fill element arrays
      for (int i=0;i<nen_;++i)
      {
        for (int k = 0; k< numscal_; ++k)
        {
          // split for each transported scalar, insert into element arrays
          ephinp_[k](i,0) = myphinp[k+(i*numdofpernode_)];
        }
      } // for i

      // global element of processor 0 is printed to the screen
      if(discretization.Comm().MyPID()==0)
        std::cout << "Electrolyte conductivity evaluated at global element " << ele->Id() << ":" << std::endl;

      CalculateConductivity(ele,frt,scatratype,elevec1_epetra);
    }
    else // conductivity = diffusivity for a electric potential field
    {
     // cout << "TEST CRISTOBAL evaluate calc_elch_conductivity" << endl;

      GetMaterialParams(ele,scatratype,0.0); // use dt=0.0 dymmy value
      elevec1_epetra(0)=diffus_[0];
      elevec1_epetra(1)=diffus_[0];
    }

    break;
  }
  case SCATRA::calc_elch_electrode_kinetics:
  {
    dserror(" ");
    break;
  }
  // calculate initial electric potential field caused by initial ion concentrations
  case SCATRA::calc_elch_initial_potential:
  {
    // set specific parameter used in diffusion conduction formulation
    // this method need to be located inside ELCH
    DiffCondParams(ele, params);

    // need initial field -> extract local values from the global vector
    Teuchos::RCP<const Epetra_Vector> phi0 = discretization.GetState("phi0");
    if (phi0==Teuchos::null) dserror("Cannot get state vector 'phi0'");
    std::vector<double> myphi0(lm.size());
    DRT::UTILS::ExtractMyValues(*phi0,myphi0,lm);

    // fill element arrays
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephinp_[k](i,0) = myphi0[k+(i*numdofpernode_)];
      }
      // get values for el. potential at element nodes
      epotnp_(i) = myphi0[i*numdofpernode_+numscal_];
    } // for i
    const double frt = params.get<double>("frt");

    if(diffcond_==false)
      CalculateElectricPotentialField(ele,frt,scatratype,elemat1_epetra,elevec1_epetra);
    else
      CalculateElectricPotentialField(ele,frt,scatratype,elemat1_epetra,elevec1_epetra,newman_);

    break;
  }

  case SCATRA::calc_integr_reaction:
  {
    // NOTE: add integral values only for elements which are NOT ghosted!
    // non reactive element actually not removed
    if (ele->Owner() == discretization.Comm().MyPID())
    {
      const double dt = params.get<double>("time-step length");
      Teuchos::RCP<std::vector<double> > myreacnp =
          params.get<Teuchos::RCP<std::vector<double> > >("local reaction integral");

      // integrations points and weights
      const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // integration loop
      for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {
        for (int k = 0; k< numscal_; ++k)
        {
          // cout << "TEST CRISTOBAL calc_integr_reaction" << endl;
          GetMaterialParams(ele,scatratype,dt);

          const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

          // scalar at integration point
          const double phi = funct_.Dot(ephinp_[k]);

          (*myreacnp)[k] += reacoeff_[k]*phi*fac;
        }
      } // loop over integration points

      params.set<Teuchos::RCP<std::vector<double> > >("local reaction integral", myreacnp);
    }
    break;
  }


  // work is done
  return 0;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 g.bau 08/08|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::Sysmat(
  DRT::Element*                         ele, ///< the element those matrix is calculated
  Epetra_SerialDenseMatrix&             emat,///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs, ///< element rhs to calculate
  Epetra_SerialDenseVector&             subgrdiff, ///< subgrid-diff.-scaling vector
  const double                          time, ///< current simulation time
  const double                          dt, ///< current time-step length
  const double                          timefac, ///< time discretization factor
  const double                          alphaF, ///< factor for generalized-alpha time integration
  const enum INPAR::SCATRA::AssgdType   whichassgd, ///< all-scale subgrid-diffusivity definition
  const enum INPAR::SCATRA::FSSUGRDIFF  whichfssgd, ///< fine-scale subgrid-diffusivity definition
  const bool                            assgd, ///< all-scale subgrid-diff. flag
  const bool                            fssgd, ///< fine-scale subgrid-diff. flag
  const double                          Cs, ///< Smagorinsky constant
  const double                          tpn, ///< turbulent Prandtl number
  const double                          Csgs_sgvel, ///< parameter of multifractal subgrid-scales
  const double                          alpha, ///< grid-filter to test-filter ratio
  const bool                            calc_N, ///< flag to activate calculation of N
  const double                          N_vel, ///< value for N if not calculated
  const enum INPAR::FLUID::RefVelocity  refvel, ///< reference velocity
  const enum INPAR::FLUID::RefLength    reflength, ///< reference length
  const double                          c_nu, ///< scaling for Re
  const bool                            nwl, ///< flag to activate near-wall limit
  const bool                            nwl_scatra, ///< flag to activate near-wall limit for scalar field
  const double                          Csgs_sgphi, ///< parameter of multifractal subgrid-scales
  const double                          c_diff, ///< scaling for Re*Pr
  const bool                            BD_gp, ///< evaluation of model coefficient at gp
//  const bool                            mfs_conservative, ///< conservative formulation of multifractal subgrid-scale modeling
  const double                          frt, ///< factor F/RT needed for ELCH calculations
  const enum INPAR::SCATRA::ScaTraType  scatratype ///< type of scalar transport problem
  )
{

//  //----------------------------------------------------------------------
//  // calculation of subgrid diffusivity and stabilization parameter(s)
//  // at element center
//  //----------------------------------------------------------------------
//
//  if (not tau_gp_)
//  {
//
//    // get velocity at element center
//    velint_.Multiply(evelnp_,funct_);
//    convelint_.Multiply(econvelnp_,funct_);
//
//    bool twoionsystem(false);
//    double resdiffus(diffus_[0]);
//    if (is_elch_) // electrochemistry problem
//    {
//      // when migration velocity is included to tau (we provide always now)
//      {
//        // compute global derivatives
//        derxy_.Multiply(xij_,deriv_);
//
//        // get "migration velocity" divided by D_k*z_k at element center
//        migvelint_.Multiply(-frt,derxy_,epotnp_);
//      }
//
//      // ELCH: special stabilization in case of binary electrolytes
//      twoionsystem= SCATRA::IsBinaryElectrolyte(valence_);
//      if (twoionsystem)
//      {
//        std::vector<int> indices_twoions = SCATRA::GetIndicesBinaryElectrolyte(valence_);
//        resdiffus = SCATRA::CalResDiffCoeff(valence_,diffus_,indices_twoions);
//#ifdef ACTIVATEBINARYELECTROLYTE
//        migrationstab_=false;
//        migrationintau_=false;
//#endif
//      }
//    }
//
//    for (int k = 0;k<numscal_;++k) // loop of each transported scalar
//    {
//      // calculation of all-scale subgrid diffusivity (artificial or due to
//      // constant-coefficient Smagorinsky model) at element center
//      if (assgd or turbmodel_ == INPAR::FLUID::smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_vreman)
//      {
//        if (assgd) dserror("Evaluation of artificial diffusion at element center not supported!");
//
//        CalcSubgrDiff(dt,timefac,whichassgd,assgd,Cs,tpn,vol,k);
//      }
//
//      // calculation of fine-scale artificial subgrid diffusivity at element center
//      if (fssgd) CalcFineScaleSubgrDiff(ele,subgrdiff,whichfssgd,Cs,tpn,vol,k);
//
//#ifdef ACTIVATEBINARYELECTROLYTE
//      if (twoionsystem && (abs(valence_[k])>EPS10))
//        CalTau(ele,resdiffus,dt,timefac,vol,k,frt,false);
//      else
//#endif
//      // calculation of stabilization parameter at element center
//      CalTau(ele,diffus_[k],dt,timefac,vol,k,frt,migrationintau_);
//    }
//
//    // compute stabilization parameter for eliminated ion species
//    if (is_elch_)
//    {
//      if(scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
//      {
//#ifdef ACTIVATEBINARYELECTROLYTE
//        if (twoionsystem && (abs(valence_[numscal_])>EPS10))
//          CalTau(ele,resdiffus,dt,timefac,vol,numscal_,frt,false);
//        else
//#endif
//        // calculation of stabilization parameter at element center
//        CalTau(ele,diffus_[numscal_],dt,timefac,vol,numscal_,frt,migrationintau_);
//      }
//    }
//  }





  // integration loop
  if (is_elch_) // electrochemistry problem
  {
    // Some safety checks. Do it here before it's too late.
    if (abs(diffreastafac_) > EPS10) dserror("Only SUPG is supported for ELCH problems");

    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

      // get concentration of transported scalar k at integration point
      // evaluation of all concentrations is necessary at this point since
      // -> homogeneous reactions of scalar k may depend on all concentrations
      // -> concentration depending material parameters for the diffusion-confection formulation
      // -> avoiding of possible errors (concentration was always defined as a vector where only on
      //    entry was filled)
      for (int k = 0;k<numscal_;++k)
      {
        //std::cout <<"ephinp_ "<< k<< ":  " <<ephinp_[k] << std::endl;
        conint_[k] = funct_.Dot(ephinp_[k]);
      }

      //----------------------------------------------------------------------
      // get material parameters (evaluation at integration point)
      //----------------------------------------------------------------------
      if (mat_gp_) GetMaterialParams(ele, scatratype,dt);

      // get velocity at integration point
      velint_.Multiply(evelnp_,funct_);
      convelint_.Multiply(econvelnp_,funct_);

      // convective part in convective form: u_x*N,x + u_y*N,y + u_z*N,z
      conv_.MultiplyTN(derxy_,convelint_);

      // momentum divergence required for conservative form
      if (is_conservative_) GetDivergence(vdiv_,evelnp_,derxy_);

      //--------------------------------------------------------------------
      // calculation of subgrid diffusivity and stabilization parameter(s)
      // at integration point
      //--------------------------------------------------------------------
      if (tau_gp_)
      {
        // compute global derivatives
        derxy_.Multiply(xij_,deriv_);

        // get "migration velocity" divided by D_k*z_k at element center
        migvelint_.Multiply(-frt,derxy_,epotnp_);

        // ELCH: special stabilization in case of binary electrolytes
        bool twoionsystem(false);
        double resdiffus(diffus_[0]);
        twoionsystem = SCATRA::IsBinaryElectrolyte(valence_);
        if (twoionsystem)
        {
          std::vector<int> indices_twoions = SCATRA::GetIndicesBinaryElectrolyte(valence_);
          resdiffus = SCATRA::CalResDiffCoeff(valence_,diffus_,indices_twoions);
#ifdef ACTIVATEBINARYELECTROLYTE
          migrationstab_=false;
          migrationintau_=false;
#endif
        }

        // loop of each transported scalar
        for (int k = 0;k<numscal_;++k)
        {
          // calculation of all-scale subgrid diffusivity (artificial or due to
          // constant-coefficient Smagorinsky model) at integration point
          if (assgd or turbmodel_ == INPAR::FLUID::smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_vreman)
          {
            if (assgd) dserror("Artificial diffusion not supported for elch!");

            CalcSubgrDiff(dt,timefac,whichassgd,assgd,Cs,tpn,vol,k);
          }

          // calculation of fine-scale artificial subgrid diffusivity
          // at integration point
          if (fssgd)
          {
            CalcFineScaleSubgrDiff(ele,subgrdiff,whichfssgd,Cs,tpn,vol,k);
            // compute gradient of fine-scale part of scalar value
            fsgradphi_.Multiply(derxy_,fsphinp_[k]);
          }

#ifdef ACTIVATEBINARYELECTROLYTE
          // use resulting diffusion coefficient for binary electrolyte solutions
          if (twoionsystem && (abs(valence_[k])>EPS10))
            CalTau(ele,resdiffus,dt,timefac,vol,k,frt,false);
          else
#endif
            // calculation of stabilization parameter at integration point
            CalTau(ele,diffus_[k],dt,timefac,vol,k,frt,migrationintau_);
        }

        if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
          dserror("Multifractal subgrid-scales not available for elch!");

        // compute stabilization parameter for eliminated ion species
        if(scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
        {
#ifdef ACTIVATEBINARYELECTROLYTE
          if (twoionsystem && (abs(valence_[numscal_])>EPS10))
            CalTau(ele,resdiffus,dt,timefac,vol,numscal_,frt,false);
          else
#endif
            // calculation of stabilization parameter at element center
            CalTau(ele,diffus_[numscal_],dt,timefac,vol,numscal_,frt,migrationintau_);
        }
      }

      for (int k = 0;k<numscal_;++k) // loop of each transported scalar
      {
        // get history data at integration point
        hist_[k] = funct_.Dot(ehist_[k]);
        // get bodyforce at integration point
        rhs_[k] = bodyforce_[k].Dot(funct_);
      }

      // safety check
      if (!is_incremental_)
        dserror("ELCH problems are always in incremental formulation");

      // compute matrix and rhs for electrochemistry problem
      if (diffcond_==false)
        CalMatElch(emat,erhs,frt,timefac,alphaF,fac,scatratype);
      // compute matrix and rhs for diffusion-conduction formulation
      else
      {
        CalMatElchBat(emat,erhs,frt,timefac,alphaF,fac,scatratype);
        //if((ele->Id()==2) and (time < 0.5 and time > 0.4))
        //  PrintEleMatToExcel(emat,erhs);
      }
    }
  } // end if (is_elch_)
  else // 'standard' scalar transport
  {
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {


      for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
      {








//        //--------------------------------------------------------------------
//        // calculation of (fine-scale) subgrid diffusivity, subgrid-scale
//        // velocity and stabilization parameter(s) at integration point
//        //--------------------------------------------------------------------
//        if (tau_gp_)
//        {
//          // calculation of all-scale subgrid diffusivity (artificial or due to
//          // constant-coefficient Smagorinsky model) at integration point
//          if (assgd or turbmodel_ == INPAR::FLUID::smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_smagorinsky or turbmodel_ == INPAR::FLUID::dynamic_vreman)
//          {
//            if (assgd) CalTau(ele,diffus_[k],dt,timefac,vol,k,0.0,false);
//
//            CalcSubgrDiff(dt,timefac,whichassgd,assgd,Cs,tpn,vol,k);
//          }
//
//          // calculation of fine-scale artificial subgrid diffusivity
//          // at integration point
//          if (fssgd)
//          {
//            CalcFineScaleSubgrDiff(ele,subgrdiff,whichfssgd,Cs,tpn,vol,k);
//          }
//
//          // calculation of subgrid-scale velocity at integration point if required
//          if (sgvel_)
//          {
//            // calculation of stabilization parameter related to fluid momentum
//            // equation at integration point
//            CalTau(ele,visc_,dt,timefac,vol,k,0.0,false);
//
//            if (scatratype != INPAR::SCATRA::scatratype_levelset)
//              CalcSubgrVelocity(ele,time,dt,timefac,k,scatratype);
//            else dserror("CalcSubgrVelocityLevelSet not available anymore");
//            //CalcSubgrVelocityLevelSet(ele,time,dt,timefac,k,ele->Id(),iquad,intpoints, iquad);
//
//            // calculation of subgrid-scale convective part
//            sgconv_.MultiplyTN(derxy_,sgvelint_);
//          }
//
//          // calculation of stabilization parameter at integration point
//          CalTau(ele,diffus_[k],dt,timefac,vol,k,0.0,false);
//        }


        // update material parameters based on inclusion of subgrid-scale
        // part of scalar (active only for mixture fraction,
        // Sutherland law and progress variable, for the time being)
        if (update_mat_)
        {
          if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
            UpdateMaterialParams(ele,mfssgphi_[k],k);
          else
            UpdateMaterialParams(ele,sgphi_[k],k);

          // recompute rhs based on updated material parameters
          rhs_[k] = bodyforce_[k].Dot(funct_)/shc_;
          rhs_[k] += thermpressdt_/shc_;
          rhs_[k] += densnp_[k]*reatemprhs_[k];
        }

        if (scatratype == INPAR::SCATRA::scatratype_poro)
        {
          Epetra_SerialDenseMatrix tempmat(emat.M(),emat.N());
          Epetra_SerialDenseVector tempvec(erhs.Length());

          // compute matrix and rhs
          CalMatAndRHS(tempmat,tempvec,fac,fssgd,timefac,dt,alphaF,k);

          //modify the elment matrix and rhs for scalar transport through porous media
          //NOTE: no stabilization terms implemented
          CalMatAndRHS_PoroScatraMod(tempmat,tempvec,fac,timefac,k,ele->Id(),iquad);

          emat += tempmat;
          erhs += tempvec;
        }
        else
          // compute matrix and rhs (standard case)
          CalMatAndRHS(emat,erhs,fac,fssgd,timefac,dt,alphaF,k);
      } // loop over each scalar
    }
  } // integration loop

  // Todo: Is there a way to implemented it nicer
  // usually, we are done here, but
  // for two certain ELCH problem formulations we have to provide
  // additional flux terms / currents across Dirichlet boundaries
  if (is_elch_)
  {
    if((scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim) or
       (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde))
    {
      double val(0.0);
      const DRT::Node* const* nodes = ele->Nodes();
      std::string condname = "Dirichlet";

      for (int vi=0; vi<nen_; ++vi)
      {
        std::vector<DRT::Condition*> dirichcond0;
        nodes[vi]->GetCondition(condname,dirichcond0);

        // there is at least one Dirichlet condition on this node
        if (dirichcond0.size() > 0)
        {
          //std::cout<<"Ele Id = "<<ele->Id()<<"  Found one Dirichlet node for vi="<<vi<<std::endl;
          const std::vector<int>*    onoff = dirichcond0[0]->Get<std::vector<int> >   ("onoff");
          for (int k=0; k<numscal_; ++k)
          {
            if ((*onoff)[k])
            {
              //std::cout<<"Dirichlet is on for k="<<k<<std::endl;
              //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<valence_[k]<<std::endl;
              const int fvi = vi*numdofpernode_+k;
              // We use the fact, that the rhs vector value for boundary nodes
              // is equivalent to the integrated negative normal flux
              // due to diffusion and migration
              val = erhs[fvi];
              erhs[vi*numdofpernode_+numscal_] += valence_[k]*(-val);
              // corresponding linearization
              for (int ui=0; ui<nen_; ++ui)
              {
                val = emat(vi*numdofpernode_+k,ui*numdofpernode_+k);
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+k)+=valence_[k]*(-val);
                val = emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_);
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+numscal_)+=valence_[k]*(-val);
              }
            }
          } // for k
        } // if Dirichlet at node vi
      } // for vi
    }  // elim

    // for concentrated solution theory (using div i as closing term for the potential)
    // additional flux terms / currents across Dirichlet boundaries
    if(diffcond_==true and newman_==true and equpot_==INPAR::ELCH::equpot_divi)
    {
      //const double faraday = INPAR::SCATRA::faraday_const;
      double val(0.0);
      const DRT::Node* const* nodes = ele->Nodes();
      std::string condname = "Dirichlet";

      for (int vi=0; vi<nen_; ++vi)
      {
        std::vector<DRT::Condition*> dirichcond0;
        nodes[vi]->GetCondition(condname,dirichcond0);

        // there is at least one Dirichlet condition on this node
        if (dirichcond0.size() > 0)
        {
          //std::cout<<"Ele Id = "<<ele->Id()<<"  Found one Dirichlet node for vi="<<vi<<std::endl;
          const std::vector<int>*    onoff = dirichcond0[0]->Get<std::vector<int> >   ("onoff");
          for (int k=0; k<numscal_; ++k)
          {
            if ((*onoff)[k])
            {
              //std::cout<<"Dirichlet is on for k="<<k<<std::endl;
              //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<valence_[k]<<std::endl;
              const int fvi = vi*numdofpernode_+k;
              // We use the fact, that the rhs vector value for boundary nodes
              // is equivalent to the integrated negative normal flux
              // due to diffusion and migration

              // scaling of div i results in a matrix with better condition number
              val = erhs[fvi];
              erhs[vi*numdofpernode_+numscal_] += valence_[k]*(-val);
              // corresponding linearization
              for (int ui=0; ui<nen_; ++ui)
              {
                val = emat(vi*numdofpernode_+k,ui*numdofpernode_+k);
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+k)+=valence_[k]*(-val);
                val = emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_);
                emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+numscal_)+=valence_[k]*(-val);
              }
            }
          } // for k
          // dirichlet condition for potential
          if ((*onoff)[numscal_])
          {
            //reacting species 0
            int k =0;

            //std::cout<<"Dirichlet is on for k="<<k<<std::endl;
            //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<valence_[k]<<std::endl;
            const int fvi = vi*numdofpernode_+numscal_;
            // We use the fact, that the rhs vector value for boundary nodes
            // is equivalent to the integrated negative normal flux
            // due to diffusion and migration

            // scaling of div i results in a matrix with better condition number
            val = erhs[fvi];
            erhs[vi*numdofpernode_+k] += 1.0/valence_[k]*(-val);
            // corresponding linearization
            for (int ui=0; ui<nen_; ++ui)
            {
              val = emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+k);
              emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += 1.0/valence_[k]*(-val);
              val = emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+numscal_);
              emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) +=1.0/valence_[k]*(-val);
            }
          }
        } // if Dirichlet at node vi
      } // for vi

      // Nernst boundary conditions have to be handled like Dirichlet conditions!!!
      std::string condname2 = "ElectrodeKinetics";
      for (int vi=0; vi<nen_; ++vi)
      {
        std::vector<DRT::Condition*> elctrodeKinetics;
        nodes[vi]->GetCondition(condname2,elctrodeKinetics);

        // there is at least one Dirichlet condition on this node
        if (elctrodeKinetics.size() == 1)
        {
          const int  kinetics = elctrodeKinetics[0]->GetInt("kinetic model");

          if (kinetics==INPAR::SCATRA::nernst)
          {
            //reacting species 0
            int k = 0;

            //std::cout<<"Dirichlet is on for k="<<k<<std::endl;
            //std::cout<<"k="<<k<<"  val="<<val<<" valence_k="<<valence_[k]<<std::endl;
            const int fvi = vi*numdofpernode_+numscal_;
            // We use the fact, that the rhs vector value for boundary nodes
            // is equivalent to the integrated negative normal flux
            // due to diffusion and migration

            // scaling of div i results in a matrix with better condition number
            val = erhs[fvi];
            erhs[vi*numdofpernode_+k] += 1.0/valence_[k]*(-val);
            // corresponding linearization
            for (int ui=0; ui<nen_; ++ui)
            {
              val = emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+k);
              emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += 1.0/valence_[k]*(-val);
              val = emat(vi*numdofpernode_+numscal_,ui*numdofpernode_+numscal_);
              emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) +=1.0/valence_[k]*(-val);
            }
          }
        } // if Dirichlet at node vi
      } // for vi
    }
    //PrintEleMatToExcel(emat,erhs);
  } // is_elch_
  return;
}

/*-----------------------------------------------------------------------*
  |  is there a reaction coupling? (private)                     mthon 09/13|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::ScaTraImpl<distype>::IsCoupledAndRead(
  const DRT::Discretization&            ScaTraDiscretization //discretisation of the ScaTra field
  )
{
    ScaTraDiscretization.GetCondition("HomoScaTraCoupling",HSTCConds_);
    int numcond = HSTCConds_.size();

    if (numcond > 0)  //if there exists condition "DESIGN HOMOGENEOUS SCATRA COUPLING VOLUME CONDITIONS"
        return true; //1;
    else
        return false; //0;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalcReacTerm(
          const std::vector<int>       stoich,                 //!<stoichometrie of current condition
          const std::string            couplingtype,           //!<type of coupling the stoichometry coefficients
          const std::string            step,                   //!<which step to evalutate the concentrations
          double*                      phis                    //!<Output
)
{
bool allpositive = true;
    for (int i=0; (unsigned)i < stoich.size(); i++)
      {
        if (stoich[i]<-EPS14)
          {
            allpositive = false;

            double ephi=0;
            if (step=="np")
                ephi = funct_.Dot(ephinp_[i]);
            else if (step=="n")
                ephi = funct_.Dot(ephin_[i]);
            else
                dserror("no valid step-input");

            if (couplingtype == "simple_multiplicative")
                *phis *=ephi;
            else //insert new Couplings here
                    dserror("couplingtype OTHER is just a dummy and hence not implemented");
          }
      }
if (allpositive)
	dserror("there must be at least one negative entry in each stoich list");
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalcReacDerivTerm(
          const std::vector<int>       stoich,                  //!<stoichometrie of current condition
          const std::string            couplingtype,            //!<type of coupling the stoichometry coefficients
          int                          toderive,                //!<concentration to be derived
          double*                      phisderiv                //!<Output
)
{
    if (stoich[toderive]<-EPS14)
      {
        for (int ii=0; (unsigned)ii < stoich.size(); ii++)
          {
            if (stoich[ii]<-EPS14)
              {
                if (couplingtype == "simple_multiplicative")
                  {
                    if (ii!=toderive) {
                        *phisderiv *= funct_.Dot(ephinp_[ii]);
                    }
                  }
                else //insert new Couplings here
                    dserror("couplingtype OTHER is just a dummy and hence not implemented");
              }
          }
      }
    else
        *phisderiv = 0;
}


///*----------------------------------------------------------------------*
//  |  get the material constants  (private)                      gjb 10/08|
//  *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void DRT::ELEMENTS::ScaTraImpl<distype>::GetMaterialParams(
//  const DRT::Element*  ele,
//  const enum INPAR::SCATRA::ScaTraType  scatratype,
//  const double dt // current time-step length
//  )
//{
//// get the material
//  Teuchos::RCP<MAT::Material> material = ele->Material();
//
//// get diffusivity / diffusivities
//  if (material->MaterialType() == INPAR::MAT::m_matlist)
//  {
//    const Teuchos::RCP<const MAT::MatList>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
//    if (actmat->NumMat() < numscal_) dserror("Not enough materials in MatList.");
//
//    for (int k = 0;k<numscal_;++k)
//    {
//      // set reaction coeff. and temperature rhs for reactive equation system to zero
//      reacoeff_[k]   = 0.0;
//      reacoeffderiv_[k]   = 0.0;
//      for (int j=0; j<numscal_; j++) {
//          (reacoeffderivmatrix_[k])[j] = 0; //reset to zero
//      }
//      reacterm_[k]   = 0.0;
//      reatemprhs_[k] = 0.0;
//
//      // set specific heat capacity at constant pressure to 1.0
//      shc_ = 1.0;
//
//      // set density at various time steps and density gradient factor to 1.0/0.0
//      densn_[k]       = 1.0;
//      densnp_[k]      = 1.0;
//      densam_[k]      = 1.0;
//      densgradfac_[k] = 0.0;
//
//      int matid = actmat->MatID(k);
//      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);
//
//      if (singlemat->MaterialType() == INPAR::MAT::m_ion)
//      {
//        const Teuchos::RCP<const MAT::Ion>& actsinglemat
//          = Teuchos::rcp_dynamic_cast<const MAT::Ion>(singlemat);
//        valence_[k] = actsinglemat->Valence();
//        diffus_[k] = actsinglemat->Diffusivity();
//        diffusvalence_[k] = valence_[k]*diffus_[k];
//
//        // Material data of eliminated ion species is read from the LAST ion material
//        // in the matlist!
//        if ((scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim) and (k==(numscal_-1)))
//        {
//          if (diffus_.size() == (unsigned) numscal_)
//          {
//            // For storing additional data, we increase the vector for
//            // diffusivity and valences by one!
//            std::cout<<"k = "<<k<<"   Did push back for diffus_ and valence_!"<<std::endl;
//            diffus_.push_back(actsinglemat->ElimDiffusivity());
//            valence_.push_back(actsinglemat->ElimValence());
//            diffusvalence_.push_back(valence_[numscal_]*diffus_[numscal_]);
//            // we also enlarge some other vectors by one
//            tau_.push_back(0.0);
//            LINALG::Matrix<nen_,1> mat(true);
//            tauderpot_.push_back(mat);
//          }
//          else if (diffus_.size() == (unsigned) (numscal_+1))
//          {
//            diffus_[numscal_]  = actsinglemat->ElimDiffusivity();
//            valence_[numscal_] = actsinglemat->ElimValence();
//            diffusvalence_[numscal_] = valence_[numscal_]*diffus_[numscal_];
//          }
//          else
//            dserror("Something is wrong with eliminated ion species data");
//          //if (ele->Id()==0)
//          //  std::cout<<"data: "<<diffus_[numscal_]<<"   "<<valence_[numscal_]<<std::endl;
//          // data check:
//          if (abs(diffus_[numscal_])< EPS13) dserror("No diffusivity for eliminated species read!");
//          if (abs(valence_[numscal_])< EPS13) dserror("No valence for eliminated species read!");
//        }
//      }
//      else if (singlemat->MaterialType() == INPAR::MAT::m_arrhenius_spec)
//      {
//        const Teuchos::RCP<const MAT::ArrheniusSpec>& actsinglemat
//          = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusSpec>(singlemat);
//
//        // compute temperature
//        const double tempnp = funct_.Dot(ephinp_[numscal_-1]);
//
//        // compute diffusivity according to Sutherland law
//        diffus_[k] = actsinglemat->ComputeDiffusivity(tempnp);
//
//        // compute reaction coefficient for species equation
//        reacoeff_[k] = actsinglemat->ComputeReactionCoeff(tempnp);
//        reacoeffderiv_[k] = reacoeff_[k];
//
//        // scalar at integration point
//        const double phi = funct_.Dot(ephinp_[k]);
//        reacterm_[k]=reacoeff_[k]*phi;
//
//        // set reaction flag to true
//        is_reactive_ = true;
//      }
//      else if (singlemat->MaterialType() == INPAR::MAT::m_arrhenius_temp)
//      {
//        if (k != numscal_-1) dserror("Temperature equation always needs to be the last variable for reactive equation system!");
//
//        const Teuchos::RCP<const MAT::ArrheniusTemp>& actsinglemat
//          = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusTemp>(singlemat);
//
//        // get specific heat capacity at constant pressure
//        shc_ = actsinglemat->Shc();
//
//        // compute species mass fraction and temperature
//        const double spmf   = funct_.Dot(ephinp_[0]);
//        const double tempnp = funct_.Dot(ephinp_[k]);
//
//        // compute diffusivity according to Sutherland law
//        diffus_[k] = actsinglemat->ComputeDiffusivity(tempnp);
//
//        // compute density based on temperature and thermodynamic pressure
//        densnp_[k] = actsinglemat->ComputeDensity(tempnp,thermpressnp_);
//
//        if (is_genalpha_)
//        {
//          // compute density at n+alpha_M
//          const double tempam = funct_.Dot(ephiam_[k]);
//          densam_[k] = actsinglemat->ComputeDensity(tempam,thermpressam_);
//
//          if (not is_incremental_)
//          {
//            // compute density at n (thermodynamic pressure approximated at n+alpha_M)
//            const double tempn = funct_.Dot(ephin_[k]);
//            densn_[k] = actsinglemat->ComputeDensity(tempn,thermpressam_);
//          }
//          else densn_[k] = 1.0;
//        }
//        else densam_[k] = densnp_[k];
//
//        // factor for density gradient
//        densgradfac_[k] = -densnp_[k]/tempnp;
//
//        // compute sum of reaction rates for temperature equation divided by specific
//        // heat capacity -> will be considered a right-hand side contribution
//        reatemprhs_[k] = actsinglemat->ComputeReactionRHS(spmf,tempnp)/shc_;
//
//        // set reaction flag to true
//        is_reactive_ = true;
//      }
//      else if (singlemat->MaterialType() == INPAR::MAT::m_scatra)
//      {
//        const Teuchos::RCP<const MAT::ScatraMat>& actsinglemat
//          = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(singlemat);
//
//        diffus_[k] = actsinglemat->Diffusivity();
//
//        // in case of reaction with constant coefficient, read coefficient and
//        // set reaction flag to true
//        reacoeff_[k] = actsinglemat->ReaCoeff();
//        if (reacoeff_[k] > EPS14) is_reactive_ = true;
//        if (reacoeff_[k] < -EPS14)
//          dserror("Reaction coefficient for species %d is not positive: %f",k, reacoeff_[k]);
//        if (is_reactive_ && is_coupled_)
//          {
//            dserror("No support of REACOEFF (in MATerials) HOMOGENEOUS SCATRA COUPLING VOLUME CONDITION at the same time");
//          }
//        else if (is_reactive_)
//        {
//            reacoeffderiv_[k] = reacoeff_[k];
//
//            // scalar at integration point
//            const double phi = funct_.Dot(ephinp_[k]);
//            reacterm_[k]=reacoeff_[k]*phi;
//        }
//        else if (is_coupled_)
//          {
//            //std::cout <<"---------------IsCoupled-Loop!\n";
//            for (int condnum = 1; (unsigned)condnum <= HSTCConds_.size(); condnum++)
//              {
//                //reading of conditions here, because of future implemantion of nonhomogeneous couplings
//                //get stoichometrie
//                const std::vector<int>   stoich = *HSTCConds_[condnum-1]->GetMutable<std::vector<int> >("stoich");
//                //get coupling type
//                const std::string  couplingtype = *HSTCConds_[condnum-1]->Get<std::string>("coupling");
//                //get reactioncoefficient
//                const double           reaccoeff = HSTCConds_[condnum-1]->GetDouble("reaccoeff");
//                if (reaccoeff<-EPS14)
//                    dserror("reaccoeff of Condition %d is negativ",condnum);
//                if (stoich[k] != 0)
//                  {
//                    double phis = 1;// scalar at integration point np
//
//                    CalcReacTerm(stoich,couplingtype,"np",&phis);
//                    reacterm_[k] += -reaccoeff*stoich[k]*phis;
//
//                    for (int j=0; j<numscal_;j++)
//                      {
//                        double phisderiv = 1; //scalarderivative at integration point np
//
//                        CalcReacDerivTerm(stoich,couplingtype,j,&phisderiv);
//                        (reacoeffderivmatrix_[k])[j] += -reaccoeff*stoich[k]*phisderiv;
//                      }
//                  } //end if(stoich[k] != 0)
//              }
//          }
//      }
//      else if (singlemat->MaterialType() == INPAR::MAT::m_biofilm)
//      {
//        const Teuchos::RCP<const MAT::Biofilm>& actsinglemat
//          = Teuchos::rcp_dynamic_cast<const MAT::Biofilm>(singlemat);
//
//        diffus_[k] = actsinglemat->Diffusivity();
//        // double rearate_k = actsinglemat->ReaRate();
//        // double satcoeff_k = actsinglemat->SatCoeff();
//
//        // set reaction flag to true
//        is_reactive_ = true;
//
//        // get substrate concentration at n+1 or n+alpha_F at integration point
//        const double csnp = funct_.Dot(ephinp_[k]);
//        //const double conp = funct_.Dot(ephinp_[1]);
//
//        // compute reaction coefficient for species equation
//        reacoeff_[k] = actsinglemat->ComputeReactionCoeff(csnp);
//        reacoeffderiv_[k] = actsinglemat->ComputeReactionCoeffDeriv(csnp);
//
//        // scalar at integration point
//        const double phi = funct_.Dot(ephinp_[k]);
//        reacterm_[k]=reacoeff_[k]*phi;
//      }
//      else if (singlemat->MaterialType() == INPAR::MAT::m_myocard)
//      {
//        Teuchos::RCP< MAT::Myocard> actsinglemat
//        = Teuchos::rcp_dynamic_cast<MAT::Myocard>(singlemat);
//
//        dsassert(numdofpernode_==1,"more than 1 dof per node for Myocard material");
//
//        // set specific heat capacity at constant pressure to 1.0
//        shc_ = 1.0;
//
//        // compute diffusivity
//        diffus_[k] = 1.0;
//        actsinglemat->ComputeDiffusivity(diffus3_[k]);
//
//        // set constant density
//        densnp_[k] = 1.0;
//        densam_[k] = 1.0;
//        densn_[k] = 1.0;
//        densgradfac_[k] = 0.0;
//
//        // set reaction and anisotropic flag to true
//        is_reactive_ = true;
//        is_anisotropic_ = true;
//
//        // get reaction coeff. and set temperature rhs for reactive equation system to zero
//        double csnp = funct_.Dot(ephinp_[k]);
//        reacoeffderiv_[k] = actsinglemat->ComputeReactionCoeffDeriv(csnp, dt);
//        reacterm_[k] = actsinglemat->ComputeReactionCoeff(csnp, dt);
//        reatemprhs_[k] = 0.0;
//      }
//      else dserror("material type not allowed");
//
//      // check whether there is negative (physical) diffusivity
//      if (diffus_[k] < -EPS15) dserror("negative (physical) diffusivity");
//    }
//  } // end if(material->MaterialType() == INPAR::MAT::m_matlist)
//  else if (material->MaterialType() == INPAR::MAT::m_scatra)
//  {
//    const Teuchos::RCP<const MAT::ScatraMat>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);
//
//    dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");
//
//    // get constant diffusivity
//    diffus_[0] = actmat->Diffusivity();
//
//    // in case of reaction with (non-zero) constant coefficient:
//    // read coefficient and set reaction flag to true
//    reacoeff_[0] = actmat->ReaCoeff();
//    if (reacoeff_[0] > EPS14) is_reactive_ = true;
//    if (reacoeff_[0] < -EPS14)
//      dserror("Reaction coefficient is not positive: %f",0, reacoeff_[0]);
//
//    if (is_reactive_ && is_coupled_)
//      {
//        dserror("No support of REACOEFF (in MATerials) HOMOGENEOUS SCATRA COUPLING VOLUME CONDITION at the same time");
//      }
//    else if (is_reactive_)
//    {
//        reacoeffderiv_[0] = reacoeff_[0];
//
//        // scalar at integration point
//        const double phi = funct_.Dot(ephinp_[0]);
//        reacterm_[0]=reacoeff_[0]*phi;
//    }
//    else if (is_coupled_)
//      {
//        //std::cout <<"---------------IsCoupled-Loop!\n";
//        for (int condnum = 1; (unsigned)condnum <= HSTCConds_.size(); condnum++)
//          {
//            //reading of conditions here, because of possible implementation of nonhomogeneous couplings in the future
//            //For higher efficiency read only once instead of for every concentration
//            //get stoichometrie
//            const std::vector<int>   stoich = *HSTCConds_[condnum-1]->GetMutable<std::vector<int> >("stoich");
//            //get coupling type
//            const std::string  couplingtype = *HSTCConds_[condnum-1]->Get<std::string>("coupling");
//            //get reactioncoefficient
//            const double           reaccoeff = HSTCConds_[condnum-1]->GetDouble("reaccoeff");
//            if (reaccoeff<-EPS14)
//                dserror("reaccoeff of Condition %d is negativ",condnum);
//            if (stoich[0] != 0)
//              {
//                double phis = 1;// scalar at integration point np
//
//                CalcReacTerm(stoich,couplingtype,"np",&phis);
//                reacterm_[0] += -reaccoeff*stoich[0]*phis;
//
//                for (int j=0; j<numscal_;j++)
//                  {
//                    double phisderiv = 1; //scalarderivative at integration point np
//
//                    CalcReacDerivTerm(stoich,couplingtype,j,&phisderiv);
//                    (reacoeffderivmatrix_[0])[j] += -reaccoeff*stoich[0]*phisderiv;
//                  }
//              } //end if(stoich[0] != 0)
//          }
//      }
//
//    // set specific heat capacity at constant pressure to 1.0
//    shc_ = 1.0;
//
//    // set temperature rhs for reactive equation system to zero
//    reatemprhs_[0] = 0.0;
//
//    // set density at various time steps and density gradient factor to 1.0/0.0
//    densn_[0]       = 1.0;
//    densnp_[0]      = 1.0;
//    densam_[0]      = 1.0;
//    densgradfac_[0] = 0.0;
//
//    // in case of multifrcatal subgrid-scales, read Schmidt number
//    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales or sgvel_
//        or turbmodel_ == INPAR::FLUID::dynamic_smagorinsky)
//    {
//      //access fluid discretization
//      Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
//      fluiddis = DRT::Problem::Instance()->GetDis("fluid");
//      //get corresponding fluid element (it has the same global ID as the scatra element)
//      DRT::Element* fluidele = fluiddis->gElement(ele->Id());
//      if (fluidele == NULL)
//        dserror("Fluid element %i not on local processor", ele->Id());
//
//      // get fluid material
//      Teuchos::RCP<MAT::Material> fluidmat = fluidele->Material();
//      if(fluidmat->MaterialType() != INPAR::MAT::m_fluid)
//        dserror("Invalid fluid material for passive scalar transport in turbulent flow!");
//
//      const Teuchos::RCP<const MAT::NewtonianFluid>& actfluidmat = Teuchos::rcp_dynamic_cast<const MAT::NewtonianFluid>(fluidmat);
//
//      // get constant dynamic viscosity
//      visc_ = actfluidmat->Viscosity();
//      densn_[0] = actfluidmat->Density();
//      densnp_[0] = actfluidmat->Density();
//      densam_[0] = actfluidmat->Density();
//
//      if (densam_[0] != 1.0 or densnp_[0] != 1.0 or densn_[0] != 1.0)
//         dserror("Check your diffusivity! Dynamic diffusivity required!");
//     }
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_ion)
//  {
//    const Teuchos::RCP<const MAT::Ion>& actsinglemat
//      = Teuchos::rcp_dynamic_cast<const MAT::Ion>(material);
//
//    dsassert(numdofpernode_==1,"more than 1 dof per node for single ion material");
//
//    // set reaction coeff. and temperature rhs for reactive equation system to zero
//    reacoeff_[0]   = 0.0;
//    reacoeffderiv_[0]   = 0.0;
//    reacterm_[0]   = 0.0;
//    reatemprhs_[0] = 0.0;
//    // set specific heat capacity at constant pressure to 1.0
//    shc_ = 1.0;
//    // set density at various time steps and density gradient factor to 1.0/0.0
//    densn_[0]       = 1.0;
//    densnp_[0]      = 1.0;
//    densam_[0]      = 1.0;
//    densgradfac_[0] = 0.0;
//
//    // get constant diffusivity
//    diffus_[0] = actsinglemat->Diffusivity();
//    valence_[0] = 0.0; // remains unused -> we only do convection-diffusion in this case!
//    diffusvalence_[0] = 0.0; // remains unused
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_mixfrac)
//  {
//    const Teuchos::RCP<const MAT::MixFrac>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::MixFrac>(material);
//
//    dsassert(numdofpernode_==1,"more than 1 dof per node for mixture-fraction material");
//
//    // compute mixture fraction at n+1 or n+alpha_F
//    const double mixfracnp = funct_.Dot(ephinp_[0]);
//
//    // compute dynamic diffusivity at n+1 or n+alpha_F based on mixture fraction
//    diffus_[0] = actmat->ComputeDiffusivity(mixfracnp);
//
//    // compute density at n+1 or n+alpha_F based on mixture fraction
//    densnp_[0] = actmat->ComputeDensity(mixfracnp);
//
//    // set specific heat capacity at constant pressure to 1.0
//    shc_ = 1.0;
//
//    if (is_genalpha_)
//    {
//      // compute density at n+alpha_M
//      const double mixfracam = funct_.Dot(ephiam_[0]);
//      densam_[0] = actmat->ComputeDensity(mixfracam);
//
//      if (not is_incremental_)
//      {
//        // compute density at n
//        const double mixfracn = funct_.Dot(ephin_[0]);
//        densn_[0] = actmat->ComputeDensity(mixfracn);
//      }
//      else densn_[0] = 1.0;
//    }
//    else densam_[0] = densnp_[0];
//
//    // factor for density gradient
//    densgradfac_[0] = -densnp_[0]*densnp_[0]*actmat->EosFacA();
//
//    // set reaction coeff. and temperature rhs for reactive equation system to zero
//    reacoeff_[0] = 0.0;
//    reacoeffderiv_[0]   = 0.0;
//    reacterm_[0]   = 0.0;
//    reatemprhs_[0] = 0.0;
//
//    // get also fluid viscosity if subgrid-scale velocity is to be included
//    // or multifractal subgrid-scales are used
//    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
//      visc_ = actmat->ComputeViscosity(mixfracnp);
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_sutherland)
//  {
//    const Teuchos::RCP<const MAT::Sutherland>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::Sutherland>(material);
//
//    dsassert(numdofpernode_==1,"more than 1 dof per node for Sutherland material");
//
//    // get specific heat capacity at constant pressure
//    shc_ = actmat->Shc();
//
//    // compute temperature at n+1 or n+alpha_F
//    const double tempnp = funct_.Dot(ephinp_[0]);
//    if (tempnp < 0.0)
//      dserror("Negative temperature occurred! Sutherland's law is defined for positive temperatures, only!");
//
//    // compute diffusivity according to Sutherland law
//    diffus_[0] = actmat->ComputeDiffusivity(tempnp);
//
//    // compute density at n+1 or n+alpha_F based on temperature
//    // and thermodynamic pressure
//    densnp_[0] = actmat->ComputeDensity(tempnp,thermpressnp_);
//
//    if (is_genalpha_)
//    {
//      // compute density at n+alpha_M
//      const double tempam = funct_.Dot(ephiam_[0]);
//      densam_[0] = actmat->ComputeDensity(tempam,thermpressam_);
//
//      if (not is_incremental_)
//      {
//        // compute density at n (thermodynamic pressure approximated at n+alpha_M)
//        const double tempn = funct_.Dot(ephin_[0]);
//        densn_[0] = actmat->ComputeDensity(tempn,thermpressam_);
//      }
//      else densn_[0] = 1.0;
//    }
//    else densam_[0] = densnp_[0];
//
//    // factor for density gradient
//    densgradfac_[0] = -densnp_[0]/tempnp;
//
//    // set reaction coeff. and temperature rhs for reactive equation system to zero
//    reacoeff_[0] = 0.0;
//    reacoeffderiv_[0] = 0.0;
//    reacterm_[0] = 0.0;
//    reatemprhs_[0] = 0.0;
//
//    // get also fluid viscosity if subgrid-scale velocity is to be included
//    // or multifractal subgrid-scales are used
//    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
//      visc_ = actmat->ComputeViscosity(tempnp);
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
//  {
//    const Teuchos::RCP<const MAT::ArrheniusPV>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusPV>(material);
//
//    dsassert(numdofpernode_==1,"more than 1 dof per node for progress-variable material");
//
//    // get progress variable at n+1 or n+alpha_F
//    const double provarnp = funct_.Dot(ephinp_[0]);
//
//    // get specific heat capacity at constant pressure and
//    // compute temperature based on progress variable
//    shc_ = actmat->ComputeShc(provarnp);
//    const double tempnp = actmat->ComputeTemperature(provarnp);
//
//    // compute density at n+1 or n+alpha_F
//    densnp_[0] = actmat->ComputeDensity(provarnp);
//
//    if (is_genalpha_)
//    {
//      // compute density at n+alpha_M
//      const double provaram = funct_.Dot(ephiam_[0]);
//      densam_[0] = actmat->ComputeDensity(provaram);
//
//      if (not is_incremental_)
//      {
//        // compute density at n
//        const double provarn = funct_.Dot(ephin_[0]);
//        densn_[0] = actmat->ComputeDensity(provarn);
//      }
//      else densn_[0] = 1.0;
//    }
//    else densam_[0] = densnp_[0];
//
//    // factor for density gradient
//    densgradfac_[0] = -densnp_[0]*actmat->ComputeFactor(provarnp);
//
//    // compute diffusivity according to Sutherland law
//    diffus_[0] = actmat->ComputeDiffusivity(tempnp);
//
//    // compute reaction coefficient for progress variable
//    reacoeff_[0] = actmat->ComputeReactionCoeff(tempnp);
//    reacoeffderiv_[0] = reacoeff_[0];
//    // compute right-hand side contribution for progress variable
//    // -> equal to reaction coefficient
//    reatemprhs_[0] = reacoeff_[0];
//
//    // scalar at integration point
//    const double phi = funct_.Dot(ephinp_[0]);
//    reacterm_[0]=reacoeff_[0]*phi;
//
//    // set reaction flag to true
//    is_reactive_ = true;
//
//    // get also fluid viscosity if subgrid-scale velocity is to be included
//    // or multifractal subgrid-scales are used
//    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
//      visc_ = actmat->ComputeViscosity(tempnp);
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
//  {
//    const Teuchos::RCP<const MAT::FerEchPV>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::FerEchPV>(material);
//
//    dsassert(numdofpernode_==1,"more than 1 dof per node for progress-variable material");
//
//    // get progress variable at n+1 or n+alpha_F
//    const double provarnp = funct_.Dot(ephinp_[0]);
//
//    // get specific heat capacity at constant pressure and
//    // compute temperature based on progress variable
//    shc_ = actmat->ComputeShc(provarnp);
//    const double tempnp = actmat->ComputeTemperature(provarnp);
//
//    // compute density at n+1 or n+alpha_F
//    densnp_[0] = actmat->ComputeDensity(provarnp);
//
//    if (is_genalpha_)
//    {
//      // compute density at n+alpha_M
//      const double provaram = funct_.Dot(ephiam_[0]);
//      densam_[0] = actmat->ComputeDensity(provaram);
//
//      if (not is_incremental_)
//      {
//        // compute density at n
//        const double provarn = funct_.Dot(ephin_[0]);
//        densn_[0] = actmat->ComputeDensity(provarn);
//      }
//      else densn_[0] = 1.0;
//    }
//    else densam_[0] = densnp_[0];
//
//    // factor for density gradient
//    densgradfac_[0] = -densnp_[0]*actmat->ComputeFactor(provarnp);
//
//    // compute diffusivity according to Sutherland law
//    diffus_[0] = actmat->ComputeDiffusivity(tempnp);
//
//    // compute reaction coefficient for progress variable
//    reacoeff_[0] = actmat->ComputeReactionCoeff(provarnp);
//    reacoeffderiv_[0] = reacoeff_[0];
//
//    // scalar at integration point
//    const double phi = funct_.Dot(ephinp_[0]);
//    reacterm_[0]=reacoeff_[0]*phi;
//
//    // compute right-hand side contribution for progress variable
//    // -> equal to reaction coefficient
//    reatemprhs_[0] = reacoeff_[0];
//
//    // set reaction flag to true
//    is_reactive_ = true;
//
//    // get also fluid viscosity if subgrid-scale velocity is to be included
//    // or multifractal subgrid-scales are used
//    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
//      visc_ = actmat->ComputeViscosity(tempnp);
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_biofilm)
//  {
//    dsassert(numdofpernode_==1,"more than 1 dof per node for BIOFILM material");
//
//    const Teuchos::RCP<const MAT::Biofilm>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::Biofilm>(material);
//
//    diffus_[0] = actmat->Diffusivity();
//    // double rearate_k = actmat->ReaRate();
//    // double satcoeff_k = actmat->SatCoeff();
//
//    // set reaction flag to true
//    is_reactive_ = true;
//
//    // get substrate concentration at n+1 or n+alpha_F at integration point
//    const double csnp = funct_.Dot(ephinp_[0]);
//    //const double conp = funct_.Dot(ephinp_[1]);
//
//    // compute reaction coefficient for species equation
//    reacoeff_[0] = actmat->ComputeReactionCoeff(csnp);
//    reacoeffderiv_[0] = actmat->ComputeReactionCoeffDeriv(csnp);
//
//    // scalar at integration point
//    const double phi = funct_.Dot(ephinp_[0]);
//    reacterm_[0]=reacoeff_[0]*phi;
//
//    // set specific heat capacity at constant pressure to 1.0
//    shc_ = 1.0;
//
//    // set temperature rhs for reactive equation system to zero
//    reatemprhs_[0] = 0.0;
//
//    // set density at various time steps and density gradient factor to 1.0/0.0
//    densn_[0]       = 1.0;
//    densnp_[0]      = 1.0;
//    densam_[0]      = 1.0;
//    densgradfac_[0] = 0.0;
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_th_fourier_iso)
//  {
//    dsassert(numdofpernode_==1,"more than 1 dof per node for isotropic Fourier material");
//
//    const Teuchos::RCP<const MAT::FourierIso>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::FourierIso>(material);
//
//    // get constant diffusivity (conductivity divided by specific heat capacity)
//    diffus_[0] = actmat->Conductivity()/actmat->Capacity();
//
//    // set density at various time steps and density gradient factor to 1.0/0.0
//    densn_[0]       = 1.0;
//    densnp_[0]      = 1.0;
//    densam_[0]      = 1.0;
//    densgradfac_[0] = 0.0;
//
//    // set specific heat capacity at constant volume
//    // (value divided by density here for its intended use on right-hand side)
//    shc_ = actmat->Capacity()/densnp_[0];
//
//    // set reaction coeff. and temperature rhs for reactive equation system to zero
//    reacterm_[0]      = 0.0;
//    reacoeff_[0]      = 0.0;
//    reacoeffderiv_[0] = 0.0;
//    reatemprhs_[0]    = 0.0;
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_thermostvenant)
//  {
//    dsassert(numdofpernode_==1,"more than 1 dof per node for thermo St. Venant-Kirchhoff material");
//
//    const Teuchos::RCP<const MAT::ThermoStVenantKirchhoff>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::ThermoStVenantKirchhoff>(material);
//
//    // get constant diffusivity (conductivity divided by specific heat capacity)
//    diffus_[0] = actmat->Conductivity()/actmat->Capacity();
//
//    // set density at various time steps and set density gradient factor to 1.0/0.0
//    densnp_[0]      = actmat->Density();
//    densam_[0]      = densnp_[0];
//    densn_[0]       = densnp_[0];
//    densgradfac_[0] = 0.0;
//
//    // set specific heat capacity at constant volume
//    // (value divided by density here for its intended use on right-hand side)
//    shc_ = actmat->Capacity()/densnp_[0];
//
//    // compute reaction coefficient
//    // (divided by density due to later multiplication by density in CalMatAndRHS)
//    Teuchos::ParameterList dummyparams;
//    const double stmodulus = actmat->STModulus(dummyparams);
//    reacoeff_[0] = -vdiv_*stmodulus/(actmat->Capacity()*densnp_[0]);
//
//    // set reaction flag to true, check whether reaction coefficient is positive
//    // and set derivative of reaction coefficient
//    if (reacoeff_[0] > EPS14 or reacoeff_[0] < -EPS14) is_reactive_ = true;
//    reacoeffderiv_[0] = reacoeff_[0];
//
//    // set temperature rhs for reactive equation system to zero
//    reatemprhs_[0] = 0.0;
//
//    // set temporal derivative of thermodynamic pressure to zero for
//    // the present structure-based scalar transport
//    thermpressdt_ = 0.0;
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_yoghurt)
//  {
//    const Teuchos::RCP<const MAT::Yoghurt>& actmat
//      = Teuchos::rcp_dynamic_cast<const MAT::Yoghurt>(material);
//
//    dsassert(numdofpernode_==1,"more than 1 dof per node for Yoghurt material");
//
//    // get specific heat capacity at constant pressure
//    shc_ = actmat->Shc();
//
//    // compute diffusivity
//    diffus_[0] = actmat->ComputeDiffusivity();
//
//    // get constant density
//    densnp_[0] = actmat->Density();
//    densam_[0] = densnp_[0];
//    densn_[0] = densnp_[0];
//
//    // set reaction coeff. and temperature rhs for reactive equation system to zero
//    reacoeff_[0] = 0.0;
//    reacoeffderiv_[0] = 0.0;
//    reacterm_[0] = 0.0;
//    reatemprhs_[0] = 0.0;
//
//    // get also fluid viscosity if subgrid-scale velocity is to be included
//    // or multifractal subgrid-scales are used
//    if (sgvel_ or turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales)
//    {
//      // compute temperature at n+1 or n+alpha_F
//      const double tempnp = funct_.Dot(ephinp_[0]);
//
//      // compute rate of strain
//      double rateofstrain = -1.0e30;
//      rateofstrain = GetStrainRate(evelnp_);
//
//      // compute viscosity for Yoghurt-like flows according to Afonso et al. (2003)
//      visc_ = actmat->ComputeViscosity(rateofstrain,tempnp);
//    }
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_myocard)
//  {
//   // cout << "TEST CRISTOBAL: ESTOY EN EL SEGUNDO myocard en GetMaterialParams" << endl;
//    //Teuchos::RCP<MAT::Myocard>& actmat
//    Teuchos::RCP<MAT::Myocard> actmat
//      = Teuchos::rcp_dynamic_cast<MAT::Myocard>(material);
//    if (scatratype != INPAR::SCATRA::scatratype_cardio_monodomain){
//      dserror("ScatraType not compatible with myocard material. Change ScatraType to some Cardio_* type.");
//    }
//
//    dsassert(numdofpernode_==1,"more than 1 dof per node for Myocard material");
//
//    // set specific heat capacity at constant pressure to 1.0
//    shc_ = 1.0;
//
//    // set diffusivity to one
//    diffus_[0] = 1.0;
//
//    //get diffusion tensor
//   actmat->ComputeDiffusivity(diffus3_[0]);
//
//    // set constant density
//    densnp_[0] = 1.0;
//    densam_[0] = 1.0;
//    densn_[0] = 1.0;
//    densgradfac_[0] = 0.0;
//
//    // set reaction and anisotropic flag to true
//    is_reactive_ = true;
//    is_anisotropic_ = true;
//
//    // get reaction coeff. and set temperature rhs for reactive equation system to zero
//    double csnp = funct_.Dot(ephinp_[0]);
//    //cout << "Calling reaction coefficients" << endl;
//    reacoeffderiv_[0] = actmat->ComputeReactionCoeffDeriv(csnp, dt);
//    reacterm_[0] = actmat->ComputeReactionCoeff(csnp, dt);
//    // cout << "TEST CRISTOBAL: reacterm_ " << reacterm_[0] << endl;
//    reatemprhs_[0] = 0.0;
//  }
//  else if (material->MaterialType() == INPAR::MAT::m_elchmat)
//    GetMaterialParamsDiffCond(material);
//  else dserror("Material type is not supported");
//
//// check whether there is negative (physical) diffusivity
//  if (diffus_[0] < -EPS15) dserror("negative (physical) diffusivity");
//  return;
//} //ScaTraImpl::GetMaterialParams

/*----------------------------------------------------------------------*
  |  evaluate element matrix and rhs (private)                   vg 02/09|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalMatAndRHS(
  Epetra_SerialDenseMatrix&             emat,
  Epetra_SerialDenseVector&             erhs,
  const double                          fac,
  const bool                            fssgd,
  const double                          timefac,
  const double                          dt,
  const double                          alphaF,
  const int                             k
  )
{










  return;
} //ScaTraImpl::CalMatAndRHS














#if 0
/*----------------------------------------------------------------------*
 |  calculate the Laplacian (weak form)             (private) ljag 11/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::GetLaplacianWeakForm(
  double& val,
  const LINALG::Matrix<nsd_,nen_>& derxy,
  const LINALG::Matrix<nsd_,nsd_>& diffus3,
  const int vi,
  const int ui)
{
  val = 0.0;
  for (int j = 0; j<nsd_; j++)
  {
  for (int i = 0; i<nsd_; i++)
      {
      val += derxy(j, vi)*diffus3(j,i)*derxy(i, ui);
      }
  }
  return;
} // ScaTraImpl<distype>::GetLaplacianWeakForm
#endif

#if 0
/*----------------------------------------------------------------------*
 |  calculate the Laplacian (weak form)             (private) gjb 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::GetLaplacianWeakForm(
  double& val,
  const LINALG::Matrix<nsd_,nen_>& derxy,
  const int vi,
  const int ui)
{
  val = 0.0;
  for (int j = 0; j<nsd_; j++)
  {
    val += derxy(j, vi)*derxy(j, ui);
  }
  return;
} // ScaTraImpl<distype>::GetLaplacianWeakForm
#endif


#if 0
/*----------------------------------------------------------------------*
 |  calculate rhs of Laplacian (weak form)          (private) gjb 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::GetLaplacianWeakFormRHS(
    double& val,
    const LINALG::Matrix<nsd_,nen_>& derxy,
    const LINALG::Matrix<nsd_,1>&   gradphi,
    const int vi)
{
  val = 0.0;
  for (int j = 0; j<nsd_; j++)
  {
    val += derxy(j,vi)*gradphi(j);
  }
  return;
} // ScaTraImpl<distype>::GetLaplacianWeakFormRHS
#endif








/*---------------------------------------------------------------------*
  |  calculate error compared to analytical solution           gjb 10/08|
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalErrorComparedToAnalytSolution(
  const DRT::Element*                   ele,
  const enum INPAR::SCATRA::ScaTraType  scatratype,
  Teuchos::ParameterList&               params,
  Epetra_SerialDenseVector&             errors
  )
{
  //at the moment, there is only one analytical test problem available!
  if (DRT::INPUT::get<SCATRA::Action>(params,"action") != SCATRA::calc_error)
    dserror("How did you get here?");

  // -------------- prepare common things first ! -----------------------
  // in the ALE case add nodal displacements
  if (is_ale_) dserror("No ALE for Kwok & Wu error calculation allowed.");

  // set constants for analytical solution
  const double t = params.get<double>("total time");
  const double frt = params.get<double>("frt");

  // get material constants
  GetMaterialParams(ele,scatratype,0.0); // use dt=0.0 dymmy value

  if(diffcond_==true)
  {
    dserror("Analytical solution for Kwok and Wu is only valid for dilute electrolyte solutions!!\n"
            "Compute corresponding transport properties on your on and activate it here");

    diffus_[0] = 2.0e-3;
    diffus_[1] = 4.0e-3;
    valence_[0] = 1.0;
    valence_[1] = -2.0;
  }

  // integrations points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  const INPAR::SCATRA::CalcError errortype = DRT::INPUT::get<INPAR::SCATRA::CalcError>(params, "calcerrorflag");
  switch(errortype)
  {
  case INPAR::SCATRA::calcerror_Kwok_Wu:
  {
    //   References:
    //   Kwok, Yue-Kuen and Wu, Charles C. K.
    //   "Fractional step algorithm for solving a multi-dimensional diffusion-migration equation"
    //   Numerical Methods for Partial Differential Equations
    //   1995, Vol 11, 389-397

    //   G. Bauer, V. Gravemeier, W.A. Wall,
    //   A 3D finite element approach for the coupled numerical simulation of
    //   electrochemical systems and fluid flow, IJNME, 86 (2011) 1339–1359.

    //if (numscal_ != 2)
    //  dserror("Numscal_ != 2 for desired error calculation.");

    // working arrays
    double                  potint(0.0);
    LINALG::Matrix<2,1>     conint(true);
    LINALG::Matrix<nsd_,1>  xint(true);
    LINALG::Matrix<2,1>     c(true);
    double                  deltapot(0.0);
    LINALG::Matrix<2,1>     deltacon(true);

    // start loop over integration points
    for (int iquad=0;iquad<intpoints.IP().nquad;iquad++)
    {
      const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

      // get values of all transported scalars at integration point
      for (int k=0; k<numscal_; ++k)
      {
        conint(k) = funct_.Dot(ephinp_[k]);
      }

      // get el. potential solution at integration point
      potint = funct_.Dot(epotnp_);

      // get global coordinate of integration point
      xint.Multiply(xyze_,funct_);

      // compute various constants
      const double d = frt*((diffus_[0]*valence_[0]) - (diffus_[1]*valence_[1]));
      if (abs(d) == 0.0) dserror("division by zero");
      const double D = frt*((valence_[0]*diffus_[0]*diffus_[1]) - (valence_[1]*diffus_[1]*diffus_[0]))/d;

      // compute analytical solution for cation and anion concentrations
      const double A0 = 2.0;
      const double m = 1.0;
      const double n = 2.0;
      const double k = 3.0;
      const double A_mnk = 1.0;
      double expterm;
      double c_0_0_0_t;

      if (nsd_==3)
      {
        expterm = exp((-D)*(m*m + n*n + k*k)*t*PI*PI);
        c(0) = A0 + (A_mnk*(cos(m*PI*xint(0))*cos(n*PI*xint(1))*cos(k*PI*xint(2)))*expterm);
        c_0_0_0_t = A0 + (A_mnk*exp((-D)*(m*m + n*n + k*k)*t*PI*PI));
      }
      else if (nsd_==2)
      {
        expterm = exp((-D)*(m*m + n*n)*t*PI*PI);
        c(0) = A0 + (A_mnk*(cos(m*PI*xint(0))*cos(n*PI*xint(1)))*expterm);
        c_0_0_0_t = A0 + (A_mnk*exp((-D)*(m*m + n*n)*t*PI*PI));
      }
      else if (nsd_==1)
      {
        expterm = exp((-D)*(m*m)*t*PI*PI);
        c(0) = A0 + (A_mnk*(cos(m*PI*xint(0)))*expterm);
        c_0_0_0_t = A0 + (A_mnk*exp((-D)*(m*m)*t*PI*PI));
      }
      else
        dserror("Illegal number of space dimensions for analyt. solution: %d",nsd_);

      // compute analytical solution for anion concentration
      c(1) = (-valence_[0]/valence_[1])* c(0);
      // compute analytical solution for el. potential
      const double pot = ((diffus_[1]-diffus_[0])/d) * log(c(0)/c_0_0_0_t);

      // compute differences between analytical solution and numerical solution
      deltapot = potint - pot;
      deltacon.Update(1.0,conint,-1.0,c);

      // add square to L2 error
      errors[0] += deltacon(0)*deltacon(0)*fac; // cation concentration
      errors[1] += deltacon(1)*deltacon(1)*fac; // anion concentration
      errors[2] += deltapot*deltapot*fac; // electric potential in electrolyte solution

    } // end of loop over integration points
  } // Kwok and Wu
  break;
  case INPAR::SCATRA::calcerror_cylinder:
  {
    // two-ion system with Butler-Volmer kinetics between two concentric cylinders
    //   G. Bauer, V. Gravemeier, W.A. Wall,
    //   A 3D finite element approach for the coupled numerical simulation of
    //   electrochemical systems and fluid flow, IJNME, 86 (2011) 1339–1359.

    if (numscal_ != 2)
      dserror("Numscal_ != 2 for desired error calculation.");

    // working arrays
    LINALG::Matrix<2,1>     conint(true);
    LINALG::Matrix<nsd_,1>  xint(true);
    LINALG::Matrix<2,1>     c(true);
    LINALG::Matrix<2,1>     deltacon(true);

    // some constants that are needed
    const double c0_inner = 0.6147737641011396;
    const double c0_outer = 1.244249192148809;
    const double r_inner = 1.0;
    const double r_outer = 2.0;
    const double pot_inner = 2.758240847314454;
    const double b = log(r_outer/r_inner);

    // start loop over integration points
    for (int iquad=0;iquad<intpoints.IP().nquad;iquad++)
    {
      const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

      // get values of all transported scalars at integration point
      for (int k=0; k<numscal_; ++k)
      {
        conint(k) = funct_.Dot(ephinp_[k]);
      }

      // get el. potential solution at integration point
      const double potint = funct_.Dot(epotnp_);

      // get global coordinate of integration point
      xint.Multiply(xyze_,funct_);

      // evaluate analytical solution for cation concentration at radial position r
      if (nsd_==3)
      {
        const double r = sqrt(xint(0)*xint(0) + xint(1)*xint(1));
        c(0) = c0_inner + ((c0_outer- c0_inner)*(log(r) - log(r_inner))/b);
      }
      else
        dserror("Illegal number of space dimensions for analyt. solution: %d",nsd_);

      // compute analytical solution for anion concentration
      c(1) = (-valence_[0]/valence_[1])* c(0);
      // compute analytical solution for el. potential
      const double d = frt*((diffus_[0]*valence_[0]) - (diffus_[1]*valence_[1]));
      if (abs(d) == 0.0) dserror("division by zero");
      // reference value + ohmic resistance + concentration potential
      const double pot = pot_inner + log(c(0)/c0_inner); // + (((diffus_[1]-diffus_[0])/d) * log(c(0)/c0_inner));

      // compute differences between analytical solution and numerical solution
      double deltapot = potint - pot;
      deltacon.Update(1.0,conint,-1.0,c);

      // add square to L2 error
      errors[0] += deltacon(0)*deltacon(0)*fac; // cation concentration
      errors[1] += deltacon(1)*deltacon(1)*fac; // anion concentration
      errors[2] += deltapot*deltapot*fac; // electric potential in electrolyte solution

    } // end of loop over integration points
  } // concentric cylinders
  break;
  case INPAR::SCATRA::calcerror_electroneutrality:
  {
    // start loop over integration points
    for (int iquad=0;iquad<intpoints.IP().nquad;iquad++)
    {
      const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

      // get values of transported scalars at integration point
      // and compute electroneutrality
      double deviation(0.0);
      for (int k=0; k<numscal_; ++k)
      {
        const double conint_k = funct_.Dot(ephinp_[k]);
        deviation += valence_[k]*conint_k;
      }

    // add square to L2 error
    errors[0] += deviation*deviation*fac;
    } // loop over integration points
  }
  break;
  default: dserror("Unknown analytical solution!"); break;
  } //switch(errortype)

  return;
} // ScaTraImpl::CalErrorComparedToAnalytSolution










/*----------------------------------------------------------------------*
 |  Do a finite difference check for a given element id. Meant for      |
 |  debugging only!                                 (private) gjb 04/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::FDcheck(
  DRT::Element*                         ele,
  Epetra_SerialDenseMatrix&             emat,
  Epetra_SerialDenseVector&             erhs,
  Epetra_SerialDenseVector&             subgrdiff,
  const double                          time,
  const double                          dt,
  const double                          timefac,
  const double                          alphaF,
  const enum INPAR::SCATRA::AssgdType   whichassgd,
  const enum INPAR::SCATRA::FSSUGRDIFF  whichfssgd,
  const bool                            assgd,
  const bool                            fssgd,
  const enum INPAR::FLUID::TurbModelAction turbmodel,
  const double                          Cs,
  const double                          tpn,
  const double                          frt,
  const enum INPAR::SCATRA::ScaTraType  scatratype
  )
{
  // magnitude of dof perturbation
  const double epsilon=1e-6; // 1.e-8 seems already too small!

  // make a copy of all input parameters potentially modified by Sysmat
  // call --- they are not intended to be modified

  // alloc the vectors that will store the original, non-perturbed values
  std::vector<LINALG::Matrix<nen_,1> > origephinp(numscal_);
  LINALG::Matrix<nsd_,nen_>       origecurnp(true);
  LINALG::Matrix<nen_,1>          origepotnp(true);
  std::vector<LINALG::Matrix<nen_,1> > origehist(numscal_);

  // copy original concentrations and potentials to these storage arrays
  for (int i=0;i<nen_;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      origephinp[k](i,0) = ephinp_[k](i,0);
      origehist[k](i,0)  = ehist_[k](i,0);
    }
    if(cursolvar_ == true)
    {
      for(int idim=0; idim<nsd_; ++idim)
        origecurnp(idim,i) = ecurnp_(idim,i);
    }

    origepotnp(i) = epotnp_(i);
  } // for i

  // allocate arrays to compute element matrices and vectors at perturbed positions
  Epetra_SerialDenseMatrix  checkmat1(emat);
  Epetra_SerialDenseVector  checkvec1(erhs);
  Epetra_SerialDenseVector  checkvec2(subgrdiff);

  // echo to screen
  printf("+-------------------------------------------+\n");
  printf("| FINITE DIFFERENCE CHECK FOR ELEMENT %5d |\n",ele->Id());
  printf("+-------------------------------------------+\n");
  printf("\n");

  // loop columns of matrix by looping nodes and then dof per nodes
  // loop nodes
  for(int nn=0;nn<nen_;++nn)
  {
    printf("-------------------------------------\n");
    printf("-------------------------------------\n");
    printf("NODE of element local id %d\n",nn);
    // loop dofs
    for(int rr=0;rr<numdofpernode_;++rr)
    {
      // number of the matrix column to check
      int dof=nn*(numdofpernode_)+rr;

      // clear element matrices and vectors to assemble
      checkmat1.Scale(0.0);
      checkvec1.Scale(0.0);
      checkvec2.Scale(0.0);

      // first put the non-perturbed values to the working arrays
      for (int i=0;i<nen_;++i)
      {
        for (int k = 0; k< numscal_; ++k)
        {
          ephinp_[k](i,0) = origephinp[k](i,0);
          ehist_[k](i,0)  = origehist[k](i,0);
        }
        if(cursolvar_ == true)
        {
          for(int idim=0; idim<nsd_; ++idim)
            ecurnp_(idim,i) = origecurnp(idim,i);
        }

        epotnp_(i) = origepotnp(i);
      } // for i

      // now perturb the respective elemental quantities
      // perturbation of the potential
      if((is_elch_) and (rr==(numscal_)))
      {
        printf("potential dof (%d). eps=%g\n",nn,epsilon);

        if (is_genalpha_)
        {
          // we want to disturb phi(n+1) with epsilon
          // => we have to disturb phi(n+alphaF) with alphaF*epsilon
          epotnp_(nn)+=(alphaF*epsilon);
        }
        else
        {
          epotnp_(nn)+=epsilon;
        }
      }
      // current only for 2D (ehrl)
      // perturbation of the current
      else if(cursolvar_==true and (rr==(numscal_+1)+0))
      {
        printf("current x dof (%d). eps=%g\n",nn,epsilon);

        if (is_genalpha_)
        {
          // we want to disturb phi(n+1) with epsilon
          // => we have to disturb phi(n+alphaF) with alphaF*epsilon
          ecurnp_(0,nn)+=(alphaF*epsilon);
        }
        else
        {
          ecurnp_(0,nn)+=epsilon;
        }
      }
      else if(cursolvar_==true and (rr==(numscal_+1)+1))
      {
        printf("current y dof (%d). eps=%g\n",nn,epsilon);

        if (is_genalpha_)
        {
          // we want to disturb phi(n+1) with epsilon
          // => we have to disturb phi(n+alphaF) with alphaF*epsilon
          ecurnp_(1,nn)+=(alphaF*epsilon);
        }
        else
        {
          ecurnp_(1,nn)+=epsilon;
        }
      }
      else
      {
        printf("concentration dof %d (%d)\n",rr,nn);

        if (is_genalpha_)
        {
          // perturbation of phi(n+1) in phi(n+alphaF) => additional factor alphaF
          ephinp_[rr](nn,0)+=(alphaF*epsilon);

          // perturbation of solution variable phi(n+1) for gen.alpha
          // leads to perturbation of phidtam (stored in ehist_)
          // with epsilon*alphaM/(gamma*dt)
          const double factor = alphaF/timefac; // = alphaM/(gamma*dt)
          ehist_[rr](nn,0)+=(factor*epsilon);

        }
        else
        {
          ephinp_[rr](nn,0)+=epsilon;
        }
      }

      // calculate the right hand side for the perturbed vector

      Sysmat(
        ele,
        checkmat1,
        checkvec1,
        checkvec2,
        time,
        dt,
        timefac,
        alphaF,
        whichassgd,
        whichfssgd,
        assgd,
        fssgd,
        turbmodel,
        Cs,
        tpn,
        0.0, // Dummy
        false, // Dummy
        0.0, // Dummy
        INPAR::FLUID::strainrate, // Dummy
        INPAR::FLUID::cube_edge, // Dummy
        0.0, // Dummy
        false,
        false, // Dummy
        0.0, // Dummy
        0.0, // Dummy
        false, // Dummy
        frt,
        scatratype);

      // compare the difference between linaer approximation and
      // (nonlinear) right hand side evaluation

      // note that it makes more sense to compare these quantities
      // than to compare the matrix entry to the difference of the
      // the right hand sides --- the latter causes numerical problems
      // do to deletion //gammi

      // however, matrix entries delivered from the element are compared
      // with the finite-difference suggestion, too. It works surprisingly well
      // for epsilon set to 1e-6 (all displayed digits nearly correct)
      // and allows a more obvious comparison!
      // when matrix entries are small, lin. and nonlin. approximation
      // look identical, although the matrix entry may be rubbish!
      // gjb

      for(int mm=0;mm<(numdofpernode_*nen_);++mm)
      {
        double val   =-erhs(mm)/epsilon;
        double lin   =-erhs(mm)/epsilon+emat(mm,dof);
        double nonlin=-checkvec1(mm)/epsilon;

        double norm=abs(lin);
        if(norm<1e-12)
        {
          norm=1e-12;
          std::cout<<"warning norm of lin is set to 10e-12"<<std::endl;
        }

        // output to screen
        {
          printf("relerr  %+12.5e   ",(lin-nonlin)/norm);
          printf("abserr  %+12.5e   ",lin-nonlin);
          printf("orig. value  %+12.5e   ",val);
          printf("lin. approx. %+12.5e   ",lin);
          printf("nonlin. funct.  %+12.5e   ",nonlin);
          printf("matrix[%d,%d]  %+12.5e   ",mm,dof,emat(mm,dof));
          // finite difference approximation (FIRST divide by epsilon and THEN subtract!)
          // ill-conditioned operation has to be done as late as possible!
          printf("FD suggestion  %+12.5e ",((erhs(mm)/epsilon)-(checkvec1(mm)/epsilon)) );
          printf("\n");
        }
      }
    }
  } // loop nodes

  // undo changes in state variables
  for (int i=0;i<nen_;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      ephinp_[k](i,0) = origephinp[k](i,0);
      ehist_[k](i,0)  = origehist[k](i,0);
    }

    if(cursolvar_== true)
    {
      for(int idim=0; idim<nsd_; ++idim)
        ecurnp_(idim,i) = origecurnp(idim,i);
    }

    epotnp_(i) = origepotnp(i);
  } // for i

  return;
}

/*----------------------------------------------------------------------*
  | calculate normalized subgrid-diffusivity matrix              vg 10/08|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalcSubgrDiffMatrix(
  const DRT::Element*           ele,
  Epetra_SerialDenseMatrix&     emat,
  const double                  timefac
  )
{
/*----------------------------------------------------------------------*/
// integration loop for one element
/*----------------------------------------------------------------------*/
// integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

// integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    for (int k=0;k<numscal_;++k)
    {
      // parameter for artificial diffusivity (scaled to one here)
      double kartfac = fac;
      if (not is_stationary_) kartfac *= timefac;

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;
          double laplawf(0.0);
          if(is_anisotropic_) dserror("Subgrid diffusivity not implemented for anisotropic materials.");
          GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
          emat(fvi,fui) += kartfac*laplawf;

          /*subtract SUPG term */
          //emat(fvi,fui) -= taufac*conv(vi)*conv(ui);
        }
      }
    }
  } // integration loop

  return;
} // ScaTraImpl::CalcSubgrDiffMatrix


/*----------------------------------------------------------------------*
  | update material parameters including s.-s. part of scalar   vg 10/11 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::UpdateMaterialParams(
  const DRT::Element*  ele,
  const double         sgphi,
  const int            k
  )
{
// get material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  if (material->MaterialType() == INPAR::MAT::m_mixfrac)
  {
    const Teuchos::RCP<const MAT::MixFrac>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::MixFrac>(material);

    // compute mixture fraction at n+1 or n+alpha_F
    double mixfracnp = funct_.Dot(ephinp_[k]);

    // add subgrid-scale part to obtain complete mixture fraction
    mixfracnp += sgphi;

    // compute dynamic diffusivity at n+1 or n+alpha_F based on mixture fraction
    diffus_[k] = actmat->ComputeDiffusivity(mixfracnp);

    // compute density at n+1 or n+alpha_F based on mixture fraction
    densnp_[k] = actmat->ComputeDensity(mixfracnp);

    if (is_genalpha_)
    {
      // compute density at n+alpha_M
      double mixfracam = funct_.Dot(ephiam_[k]);
      mixfracam += sgphi;
      densam_[k] = actmat->ComputeDensity(mixfracam);

      if (not is_incremental_)
      {
        // compute density at n
        double mixfracn = funct_.Dot(ephin_[k]);
        mixfracn += sgphi;
        densn_[k] = actmat->ComputeDensity(mixfracn);
      }
      else densn_[k] = 1.0;
    }
    else densam_[k] = densnp_[k];

    // factor for density gradient
    densgradfac_[k] = -densnp_[k]*densnp_[k]*actmat->EosFacA();
  }
  else if (material->MaterialType() == INPAR::MAT::m_sutherland)
  {
    const Teuchos::RCP<const MAT::Sutherland>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::Sutherland>(material);

    // compute temperature at n+1 or n+alpha_F
    double tempnp = funct_.Dot(ephinp_[k]);

    // add subgrid-scale part to obtain complete temperature
    tempnp += sgphi;
    if (tempnp < 0.0)
      dserror("Negative temperature occurred! Sutherland's law is defined for positive temperatures, only!");

    // compute diffusivity according to Sutherland law
    diffus_[k] = actmat->ComputeDiffusivity(tempnp);

    // compute density at n+1 or n+alpha_F based on temperature
    // and thermodynamic pressure
    densnp_[k] = actmat->ComputeDensity(tempnp,thermpressnp_);

    if (is_genalpha_)
    {
      // compute density at n+alpha_M
      double tempam = funct_.Dot(ephiam_[k]);
      tempam += sgphi;
      densam_[k] = actmat->ComputeDensity(tempam,thermpressam_);

      if (not is_incremental_)
      {
        // compute density at n (thermodynamic pressure approximated at n+alpha_M)
        double tempn = funct_.Dot(ephin_[k]);
        tempn += sgphi;
        densn_[k] = actmat->ComputeDensity(tempn,thermpressam_);
      }
      else densn_[k] = 1.0;
    }
    else densam_[k] = densnp_[k];

    // factor for density gradient
    densgradfac_[k] = -densnp_[k]/tempnp;
  }
  else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
  {
    const Teuchos::RCP<const MAT::ArrheniusPV>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::ArrheniusPV>(material);

    // get progress variable at n+1 or n+alpha_F
    double provarnp = funct_.Dot(ephinp_[k]);

    // add subgrid-scale part to obtain complete progress variable
    provarnp += sgphi;

    // get specific heat capacity at constant pressure and
    // compute temperature based on progress variable
    shc_ = actmat->ComputeShc(provarnp);
    const double tempnp = actmat->ComputeTemperature(provarnp);

    // compute density at n+1 or n+alpha_F
    densnp_[k] = actmat->ComputeDensity(provarnp);

    if (is_genalpha_)
    {
      // compute density at n+alpha_M
      double provaram = funct_.Dot(ephiam_[k]);
      provaram += sgphi;
      densam_[k] = actmat->ComputeDensity(provaram);

      if (not is_incremental_)
      {
        // compute density at n
        double provarn = funct_.Dot(ephin_[k]);
        provarn += sgphi;
        densn_[k] = actmat->ComputeDensity(provarn);
      }
      else densn_[k] = 1.0;
    }
    else densam_[k] = densnp_[k];

    // factor for density gradient
    densgradfac_[k] = -densnp_[k]*actmat->ComputeFactor(provarnp);

    // compute diffusivity according to Sutherland law
    diffus_[k] = actmat->ComputeDiffusivity(tempnp);

    // compute reaction coefficient for progress variable
    reacoeff_[k] = actmat->ComputeReactionCoeff(tempnp);
    reacoeffderiv_[k] = reacoeff_[k];
    // compute right-hand side contribution for progress variable
    // -> equal to reaction coefficient
    reatemprhs_[k] = reacoeff_[k];
  }
  else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
  {
    const Teuchos::RCP<const MAT::FerEchPV>& actmat
             = Teuchos::rcp_dynamic_cast<const MAT::FerEchPV>(material);

    // get progress variable at n+1 or n+alpha_F
    double provarnp = funct_.Dot(ephinp_[k]);

    // add subgrid-scale part to obtain complete progress variable
    provarnp += sgphi;

    // get specific heat capacity at constant pressure and
    // compute temperature based on progress variable
    shc_ = actmat->ComputeShc(provarnp);
    const double tempnp = actmat->ComputeTemperature(provarnp);

    // compute density at n+1 or n+alpha_F
    densnp_[k] = actmat->ComputeDensity(provarnp);

    if (is_genalpha_)
    {
      // compute density at n+alpha_M
      double provaram = funct_.Dot(ephiam_[k]);
      provaram += sgphi;
      densam_[k] = actmat->ComputeDensity(provaram);

      if (not is_incremental_)
      {
        // compute density at n
        double provarn = funct_.Dot(ephin_[k]);
        provarn += sgphi;
        densn_[k] = actmat->ComputeDensity(provarn);
      }
      else densn_[k] = 1.0;
    }
    else densam_[k] = densnp_[k];

    // factor for density gradient
    densgradfac_[k] = -densnp_[k]*actmat->ComputeFactor(provarnp);

    // compute diffusivity according to Sutherland law
    diffus_[k] = actmat->ComputeDiffusivity(tempnp);

    // compute reaction coefficient for progress variable
    reacoeff_[k] = actmat->ComputeReactionCoeff(provarnp);
    reacoeffderiv_[k] = reacoeff_[k];
    // compute right-hand side contribution for progress variable
    // -> equal to reaction coefficient
    reatemprhs_[k] = reacoeff_[k];
  }

  return;
} //ScaTraImpl::UpdateMaterialParams


///*----------------------------------------------------------------------*
//  |  calculate stabilization parameter  (private)              gjb 06/08 |
//  *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void DRT::ELEMENTS::ScaTraImpl<distype>::CalTau(
//  DRT::Element*                         ele,
//  double                                diffus,
//  const double                          dt,
//  const double                          timefac,
//  const double                          vol,
//  const int                             k,
//  const double                          frt,
//  const bool                            migrationintau
//  )
//{
//  // get element-type constant for tau
//  const double mk = SCATRA::MK<distype>();
//  // reset
//  tauderpot_[k].Clear();
//
//  //----------------------------------------------------------------------
//  // computation of stabilization parameters depending on definition used
//  //----------------------------------------------------------------------
//  switch (whichtau_)
//  {
//  case INPAR::SCATRA::tau_taylor_hughes_zarins:
//  case INPAR::SCATRA::tau_taylor_hughes_zarins_wo_dt:
//  {
//    /*
//
//    literature:
//    1) C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
//    of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
//    (1998) 155-196.
//    2) V. Gravemeier, W.A. Wall, An algebraic variational multiscale-
//    multigrid method for large-eddy simulation of turbulent variable-
//    density flow at low Mach number, J. Comput. Phys. 229 (2010)
//    6047-6070.
//    -> version for variable-density scalar transport equation as
//    implemented here, which corresponds to constant-density
//    version as given in the previous publication when density
//    is constant
//
//    1
//    +-                                               -+ - -
//    |        2                                        |   2
//    | c_1*rho                                  2      |
//    tau = C * | -------   +  c_2*rho*u*G*rho*u  +  c_3*mu *G:G  |
//    |     2                                           |
//    |   dt                                            |
//    +-                                               -+
//
//    with the constants and covariant metric tensor defined as follows:
//
//    C   = 1.0 (not explicitly defined here),
//    c_1 = 4.0,
//    c_2 = 1.0 (not explicitly defined here),
//    c_3 = 12.0/m_k (36.0 for linear and 144.0 for quadratic elements)
//
//    +-           -+   +-           -+   +-           -+
//    |             |   |             |   |             |
//    |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
//    G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
//    ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
//    |    i     j  |   |    i     j  |   |    i     j  |
//    +-           -+   +-           -+   +-           -+
//
//    +----
//    \
//    G : G =   +   G   * G
//    /     ij    ij
//    +----
//    i,j
//    +----
//    \
//    rho*u*G*rho*u  =   +   rho*u * G  *rho*u
//    /        i   ij      j
//    +----
//    i,j
//    */
//    // effective velocity at element center:
//    // (weighted) convective velocity + individual migration velocity
//    LINALG::Matrix<nsd_,1> veleff(convelint_,false);
//    if (is_elch_)
//    {
//      if (migrationintau) veleff.Update(diffusvalence_[k],migvelint_,1.0);
//    }
//
//    // total reaction coefficient sigma_tot: sum of "artificial" reaction
//    // due to time factor and reaction coefficient (reaction coefficient
//    // ensured to be zero in GetMaterialParams for non-reactive material)
//    double sigma_tot = reacoeff_[k];
//    if (whichtau_ == INPAR::SCATRA::tau_taylor_hughes_zarins) sigma_tot += 1.0/dt;
//
//    // computation of various values derived from covariant metric tensor
//    double G;
//    double normG(0.0);
//    double Gnormu(0.0);
//    const double dens_sqr = densnp_[k]*densnp_[k];
//    for (int nn=0;nn<nsd_;++nn)
//    {
//      for (int rr=0;rr<nsd_;++rr)
//      {
//        G = xij_(nn,0)*xij_(rr,0);
//        for(int tt=1;tt<nsd_;tt++)
//        {
//          G += xij_(nn,tt)*xij_(rr,tt);
//        }
//        normG+=G*G;
//        Gnormu+=dens_sqr*veleff(nn,0)*G*veleff(rr,0);
//        if (is_elch_) // ELCH
//        {
//          if (migrationintau)
//          {
//            // for calculation of partial derivative of tau
//            for (int jj=0;jj < nen_; jj++)
//              (tauderpot_[k])(jj,0) += dens_sqr*frt*diffusvalence_[k]*((derxy_(nn,jj)*G*veleff(rr,0))+(veleff(nn,0)*G*derxy_(rr,jj)));
//          }
//        } // ELCH
//      }
//    }
//
//    // definition of constants as described above
//    const double c1 = 4.0;
//    const double c3 = 12.0/mk;
//
//    // compute diffusive part
//    const double Gdiff = c3*diffus*diffus*normG;
//
//    // computation of stabilization parameter tau
//    tau_[k] = 1.0/(sqrt(c1*dens_sqr*DSQR(sigma_tot) + Gnormu + Gdiff));
//
//    // finalize derivative of present tau w.r.t electric potential
//    if (is_elch_)
//    {
//      if (migrationintau) tauderpot_[k].Scale(0.5*tau_[k]*tau_[k]*tau_[k]);
//    }
//  }
//  break;
//  case INPAR::SCATRA::tau_franca_valentin:
//  {
//    /*
//
//    literature:
//    L.P. Franca, F. Valentin, On an improved unusual stabilized
//    finite element method for the advective-reactive-diffusive
//    equation, Comput. Methods Appl. Mech. Engrg. 190 (2000) 1785-1800.
//
//
//    xi1,xi2 ^
//    |      /
//    |     /
//    |    /
//    1 +---+
//    |
//    |
//    |
//    +--------------> re1,re2
//    1
//
//    */
//    // get Euclidean norm of (weighted) velocity at element center
//    double vel_norm;
//    if (is_elch_ and migrationintau) migrationstab_=false;
//    // dserror("FrancaValentin with migrationintau not available at the moment");
//    /*
//    // get Euclidean norm of effective velocity at element center:
//    // (weighted) convective velocity + individual migration velocity
//    LINALG::Matrix<nsd_,1> veleff(velint_,false);
//
//    veleff.Update(diffusvalence_[k],migvelint_,1.0);
//    vel_norm = veleff.Norm2();
//
//    #ifdef VISUALIZE_ELEMENT_DATA
//    veleff.Update(diffusvalence_[k],migvelint_,0.0);
//    double vel_norm_mig = veleff.Norm2();
//    double migepe2 = mk * vel_norm_mig * h / diffus;
//
//    DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
//    if (!actele) dserror("cast to Transport* failed");
//    std::vector<double> v(1,migepe2);
//    std::ostringstream temp;
//    temp << k;
//    std::string name = "Pe_mig_"+temp.str();
//    actele->AddToData(name,v);
//    name = "hk_"+temp.str();
//    v[0] = h;
//    actele->AddToData(name,v);
//    #endif
//    }
//    else*/
//    vel_norm = convelint_.Norm2();
//
//    // total reaction coefficient sigma_tot: sum of "artificial" reaction
//    // due to time factor and reaction coefficient (reaction coefficient
//    // ensured to be zero in GetMaterialParams for non-reactive material)
//    const double sigma_tot = 1.0/timefac + reacoeff_[k];
//
//    // calculate characteristic element length
//    const double h = CalcCharEleLength(vol,vel_norm);
//
//    // various parameter computations:
//    // relating convective to viscous part
//    if (diffus < EPS14) dserror("Invalid diffusion coefficent");
//    const double epe = mk * densnp_[k] * vel_norm * h / diffus;
//    // relating viscous to reactive part
//    const double epe1 = 2.0*diffus/(mk*densnp_[k]*sigma_tot*DSQR(h));
//
//    // respective "switching" parameters
//    const double xi  = std::max(epe,1.0);
//    const double xi1 = std::max(epe1,1.0);
//
//    tau_[k] = DSQR(h)/(DSQR(h)*densnp_[k]*sigma_tot*xi1 + 2.0*diffus*xi/mk);
//
//#ifdef VISUALIZE_ELEMENT_DATA
//    // visualize resultant Pe number
//    DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
//    if (!actele) dserror("cast to Transport* failed");
//    std::vector<double> v(1,epe);
//    std::ostringstream temp;
//    temp << k;
//    std::string name = "Pe_"+temp.str();
//    actele->AddToData(name,v);
//#endif
//  }
//  break;
//  case INPAR::SCATRA::tau_franca_valentin_wo_dt:
//  {
//    /*
//
//    stabilization parameter as above without inclusion of dt-part
//
//    */
//    // get Euclidean norm of (weighted) velocity at element center
//    double vel_norm;
//    if (is_elch_ and migrationintau) migrationstab_=false;
//    // dserror("FrancaValentin with migrationintau not available at the moment");
//    /*
//    // get Euclidean norm of effective velocity at element center:
//    // (weighted) convective velocity + individual migration velocity
//    LINALG::Matrix<nsd_,1> veleff(velint_,false);
//
//    veleff.Update(diffusvalence_[k],migvelint_,1.0);
//    vel_norm = veleff.Norm2();
//
//    #ifdef VISUALIZE_ELEMENT_DATA
//    veleff.Update(diffusvalence_[k],migvelint_,0.0);
//    double vel_norm_mig = veleff.Norm2();
//    double migepe2 = mk * vel_norm_mig * h / diffus;
//
//    DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
//    if (!actele) dserror("cast to Transport* failed");
//    std::vector<double> v(1,migepe2);
//    std::ostringstream temp;
//    temp << k;
//    std::string name = "Pe_mig_"+temp.str();
//    actele->AddToData(name,v);
//    name = "hk_"+temp.str();
//    v[0] = h;
//    actele->AddToData(name,v);
//    #endif
//    }
//    else*/
//    vel_norm = convelint_.Norm2();
//
//    // calculate characteristic element length
//    const double h = CalcCharEleLength(vol,vel_norm);
//
//    // various parameter computations for case without dt:
//    // relating convective to viscous part
//    if (diffus < EPS14) dserror("Invalid diffusion coefficent");
//    const double epe = mk * densnp_[k] * vel_norm * h / diffus;
//    // relating viscous to reactive part
//    double epe1 = 0.0;
//    if (is_reactive_) epe1 = 2.0*diffus/(mk*densnp_[k]*reacoeff_[k]*DSQR(h));
//    if (is_coupled_) dserror("somthing to do here for homogeneous scatra coupling");
//
//    // respective "switching" parameters
//    const double xi  = std::max(epe,1.0);
//    const double xi1 = std::max(epe1,1.0);
//
//    tau_[k] = DSQR(h)/(DSQR(h)*densnp_[k]*reacoeff_[k]*xi1 + 2.0*diffus*xi/mk);
//
//#ifdef VISUALIZE_ELEMENT_DATA
//    // visualize resultant Pe number
//    DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
//    if (!actele) dserror("cast to Transport* failed");
//    std::vector<double> v(1,epe);
//    std::ostringstream temp;
//    temp << k;
//    std::string name = "Pe_"+temp.str();
//    actele->AddToData(name,v);
//#endif
//  }
//  break;
//  case INPAR::SCATRA::tau_shakib_hughes_codina:
//  case INPAR::SCATRA::tau_shakib_hughes_codina_wo_dt:
//  {
//    /*
//
//    literature:
//    1) F. Shakib, Finite element analysis of the compressible Euler and
//    Navier-Stokes equations, PhD thesis, Division of Applied Mechanics,
//    Stanford University, Stanford, CA, USA, 1989.
//    2) F. Shakib, T.J.R. Hughes, A new finite element formulation for
//    computational fluid dynamics: IX. Fourier analysis of space-time
//    Galerkin/least-squares algorithms, Comput. Methods Appl. Mech.
//    Engrg. 87 (1991) 35-58.
//    3) R. Codina, Stabilized finite element approximation of transient
//    incompressible flows using orthogonal subscales, Comput. Methods
//    Appl. Mech. Engrg. 191 (2002) 4295-4321.
//
//    All those proposed definitions were for non-reactive incompressible
//    flow; they are adapted to potentially reactive scalar transport
//    equations with potential density variations here.
//
//    constants defined as in Shakib (1989) / Shakib and Hughes (1991),
//    merely slightly different with respect to c_3:
//
//    c_1 = 4.0,
//    c_2 = 4.0,
//    c_3 = 4.0/(m_k*m_k) (36.0 for linear, 576.0 for quadratic ele.)
//
//    Codina (2002) proposed present version without dt and explicit
//    definition of constants
//    (condition for constants as defined here: c_2 <= sqrt(c_3)).
//
//    */
//    // get Euclidean norm of velocity
//    const double vel_norm = convelint_.Norm2();
//    if (is_elch_ and migrationintau) migrationstab_=false;
//
//    // total reaction coefficient sigma_tot: sum of "artificial" reaction
//    // due to time factor and reaction coefficient (reaction coefficient
//    // ensured to be zero in GetMaterialParams for non-reactive material)
//    double sigma_tot = reacoeff_[k];
//    if (whichtau_ == INPAR::SCATRA::tau_shakib_hughes_codina) sigma_tot += 1.0/dt;
//
//    // calculate characteristic element length
//    const double h = CalcCharEleLength(vol,vel_norm);
//
//    // definition of constants as described above
//    const double c1 = 4.0;
//    const double c2 = 4.0;
//    const double c3 = 4.0/(mk*mk);
//    // alternative value as proposed in Shakib (1989): c3 = 16.0/(mk*mk);
//
//    tau_[k] = 1.0/(sqrt(c1*DSQR(densnp_[k])*DSQR(sigma_tot)
//                        + c2*DSQR(densnp_[k])*DSQR(vel_norm)/DSQR(h)
//                        + c3*DSQR(diffus)/(DSQR(h)*DSQR(h))));
//  }
//  break;
//  case INPAR::SCATRA::tau_codina:
//  case INPAR::SCATRA::tau_codina_wo_dt:
//  {
//    /*
//
//    literature:
//    R. Codina, Comparison of some finite element methods for solving
//    the diffusion-convection-reaction equation, Comput. Methods
//    Appl. Mech. Engrg. 156 (1998) 185-210.
//
//    constants:
//    c_1 = 1.0,
//    c_2 = 2.0,
//    c_3 = 4.0/m_k (12.0 for linear, 48.0 for quadratic elements)
//
//    Codina (1998) proposed present version without dt.
//
//    */
//    // get Euclidean norm of velocity
//    const double vel_norm = convelint_.Norm2();
//
//    // total reaction coefficient sigma_tot: sum of "artificial" reaction
//    // due to time factor and reaction coefficient (reaction coefficient
//    // ensured to be zero in GetMaterialParams for non-reactive material)
//    double sigma_tot = reacoeff_[k];
//    if (whichtau_ == INPAR::SCATRA::tau_codina) sigma_tot += 1.0/dt;
//
//    // calculate characteristic element length
//    const double h = CalcCharEleLength(vol,vel_norm);
//
//    // definition of constants as described above
//    const double c1 = 1.0;
//    const double c2 = 2.0;
//    const double c3 = 4.0/mk;
//
//    tau_[k] = 1.0/(c1*densnp_[k]*sigma_tot
//                   + c2*densnp_[k]*vel_norm/h
//                   + c3*diffus/(h*h));
//  }
//  break;
//  case INPAR::SCATRA::tau_franca_madureira_valentin:
//  case INPAR::SCATRA::tau_franca_madureira_valentin_wo_dt:
//  {
//    /*
//
//    This stabilization parameter is only intended to be used for
//    reactive-diffusive problems such as structure-based scalar
//    transport problems in case of potentially dominating reaction.
//
//    literature:
//    L.P. Franca, A.L. Madureira, F. Valentin, Towards multiscale
//    functions: enriching finite element spaces with local but not
//    bubble-like functions, Comput. Methods Appl. Mech. Engrg. 194
//    (2005) 3006-3021.
//
//    */
//    // get Euclidean norm of velocity at element center
////    double vel_norm = 0.0;
////    vel_norm = convelint_.Norm2();
//
//    // total reaction coefficient sigma_tot: sum of "artificial" reaction
//    // due to time factor and reaction coefficient (reaction coefficient
//    // ensured to be zero in GetMaterialParams for non-reactive material)
//    double sigma_tot = reacoeff_[k];
//    if (whichtau_ == INPAR::SCATRA::tau_franca_madureira_valentin)
//      sigma_tot += 1.0/timefac;
//
//    // calculate characteristic element length
//    // -> currently: cubic/square root of element volume/area or
//    //    element length (3-/2-/1-D)
//    // cast dimension to a double variable -> pow()
//    const double dim = (double) nsd_;
//    const double h = std::pow(vol,1/dim);
//
//
//    // parameter relating reactive to diffusive part
//    const double epe = 2.0*diffus/(mk*densnp_[k]*sigma_tot*DSQR(h));
//
//    // respective "switching" parameter
//    const double xi = std::max(epe,1.0);
//
//    // constant c_u as suggested in Badia and Codina (2010), method A
//    // is set to be 1.0 here as in Franca et al. (2005)
//    // alternative: 4.0 as suggested in Badia and Codina (2010) for
//    // Darcy flow
//    const double c_u = 1.0;
//
//    tau_[k] = DSQR(h)/(c_u*DSQR(h)*densnp_[k]*sigma_tot*xi + (2.0*diffus/mk));
//  }
//  break;
//  case INPAR::SCATRA::tau_exact_1d:
//  {
//    // get number of dimensions (convert from int to double)
//    const double dim = (double) nsd_;
//
//    // get characteristic element length
//    double h = std::pow(vol,(1.0/dim)); // equals streamlength in 1D
//
//    // get Euclidean norm of (weighted) velocity at element center
//    double vel_norm(0.0);
//
//    if (is_elch_ and migrationintau) // ELCH
//    {
//      dserror("Migration in tau not considered in Tau_Exact_1d");
//    }
//    else
//      vel_norm = convelint_.Norm2();
//
//    if (diffus < EPS14) dserror("Invalid diffusion coefficent");
//    double epe = 0.5 * densnp_[k] * vel_norm * h / diffus;
//
//    const double pp = exp(epe);
//    const double pm = exp(-epe);
//    double xi = 0.0;
//    if (epe >= 700.0)
//      tau_[k] = 0.5*h/vel_norm;
//    else if (epe < 700.0 and epe > EPS15)
//    {
//      xi = (((pp+pm)/(pp-pm))-(1.0/epe)); // xi = coth(epe) - 1/epe
//      // compute optimal stabilization parameter
//      tau_[k] = 0.5*h*xi/vel_norm;
//
//#if 0
//      std::cout<<"epe = "<<epe<<std::endl;
//      std::cout<<"xi_opt  = "<<xi<<std::endl;
//      std::cout<<"vel_norm  = "<<vel_norm<<std::endl;
//      std::cout<<"tau_opt = "<<tau_[k]<<std::endl<<std::endl;
//#endif
//    }
//    else tau_[k] = 0.0;
//  }
//  break;
//  case INPAR::SCATRA::tau_zero:
//  {
//    // set tau's to zero (-> no stabilization effect)
//    tau_[k] = 0.0;
//  }
//  break;
//  default: dserror("unknown definition for stabilization parameter tau\n"); break;
//  } //switch (whichtau_)
//
//#if 0
//  std::cout<<"diffus  for k "<<k <<" is = "<<diffus<<std::endl;
//#endif
//#ifdef VISUALIZE_ELEMENT_DATA
//  // visualize stabilization parameter
//  DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(ele);
//  if (!actele) dserror("cast to Transport* failed");
//  std::vector<double> v(1,tau_[k]);
//  std::ostringstream temp;
//  temp << k;
//  std::string name = "tau_"+ temp.str();
//  actele->AddToData(name,v);
//#endif
//
//  return;
//} //ScaTraImpl::CalTau
//
//
///*----------------------------------------------------------------------*
//  |  calculation of characteristic element length               vg 01/11 |
//  *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//double DRT::ELEMENTS::ScaTraImpl<distype>::CalcCharEleLength(
//  const double  vol,
//  const double  vel_norm
//  )
//{
//  // define and initialize streamlength
//  double h = 0.0;
//
//  //---------------------------------------------------------------------
//  // select from various definitions for characteristic element length
//  //---------------------------------------------------------------------
//  switch (charelelength_)
//  {
//    // a) streamlength due to Tezduyar et al. (1992) -> default
//    // normed velocity vector
//    case INPAR::SCATRA::streamlength:
//    {
//      LINALG::Matrix<nsd_,1> velino(true);
//      if (vel_norm>=1e-6) velino.Update(1.0/vel_norm,convelint_);
//      else
//      {
//        velino.Clear();
//        velino(0,0) = 1.0;
//      }
//
//      // get streamlength using the normed velocity at element centre
//      LINALG::Matrix<nen_,1> tmp;
//      tmp.MultiplyTN(derxy_,velino);
//      const double val = tmp.Norm1();
//      h = 2.0/val; // h=streamlength
//    }
//    break;
//
//    // b) volume-equivalent diameter (warning: 3-D formula!)
//    case INPAR::SCATRA::volume_equivalent_diameter:
//    {
//      h = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);
//    }
//    break;
//
//    // c) cubic/square root of element volume/area or element length (3-/2-/1-D)
//    case INPAR::SCATRA::root_of_volume:
//    {
//      // cast dimension to a double varibale -> pow()
//      const double dim = double (nsd_);
//      h = std::pow(vol,1/dim);
//    }
//    break;
//
//    default: dserror("unknown characteristic element length\n");
//    break;
//  } //switch (charelelength_)
//
//  return h;
//}








/*----------------------------------------------------------------------*
  | calculate matrix and rhs for electrochemistry problem      gjb 10/08 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalMatElch(
  Epetra_SerialDenseMatrix&             emat,
  Epetra_SerialDenseVector&             erhs,
  const double                          frt,
  const double                          timefac,
  const double                          alphaF,
  const double                          fac,
  const enum INPAR::SCATRA::ScaTraType  scatratype
  )
{
  //dielectical constant
  //TODO(ehrl): good practice?
  const double epsilon = 1.e-4;
  const double faraday = INPAR::SCATRA::faraday_const;

  // get gradient of electric potential at integration point
  gradpot_.Multiply(derxy_,epotnp_);

  // migration term (convective part without z_k D_k): -F/RT\grad{\Phi}\grad
  migconv_.MultiplyTN(-frt,derxy_,gradpot_);

  // Laplacian of shape functions at integration point
  if (use2ndderiv_)
  {
    GetLaplacianStrongForm(laplace_, derxy2_);
  }

#if 0
  // DEBUG output
  std::cout<<std::endl<<"values at GP:"<<std::endl;
  std::cout<<"factor F/RT = "<<frt<<std::endl;
  for (int k=0;k<numscal_;++k)
  {std::cout<<"conint_["<<k<<"] = "<<conint_[k]<<std::endl;}
  for (int k=0;k<nsd_;++k)
  {std::cout<<"gradpot_["<<k<<"] = "<<gradpot_(k)<<std::endl;}
#endif


  for (int k = 0; k < numscal_;++k) // loop over all transported scalars
  {
    // get value of transported scalar k at integration point

    // compute gradient of scalar k at integration point
    gradphi_.Multiply(derxy_,ephinp_[k]);

    // factor D_k * z_k
    const double diffus_valence_k = diffusvalence_[k];

    double diff_ephinp_k(0.0);
    double migrea_k(0.0);
    if (use2ndderiv_) // only necessary for higher order elements
    {
      diff_.Clear();
      migrea_.Clear();

      // diffusive part:  diffus_k * ( N,xx  +  N,yy +  N,zz )
      diff_.Update(diffus_[k],laplace_);

      // get Laplacian of electric potential at integration point
      double lappot = laplace_.Dot(epotnp_);
      // reactive part of migration term
      migrea_.Update(-frt*diffus_valence_k*lappot,funct_);

      diff_ephinp_k = diff_.Dot(ephinp_[k]);   // diffusion
      migrea_k      = migrea_.Dot(ephinp_[k]); // reactive part of migration term
    }
    else
    {
      diff_.Clear();
      migrea_.Clear();
    }

    // further short cuts and definitions
    const double conv_ephinp_k = conv_.Dot(ephinp_[k]);
    const double Dkzk_mig_ephinp_k = diffus_valence_k*(migconv_.Dot(ephinp_[k]));
    const double conv_eff_k = conv_ephinp_k + Dkzk_mig_ephinp_k;

    const double taufac = tau_[k]*fac;  // corresponding stabilization parameter
    double rhsint       = rhs_[k]; // source/sink terms at int. point
    double residual     = 0.0;
    double timefacfac   = 0.0;
    double timetaufac   = 0.0;
    double rhsfac       = 0.0;
    double rhstaufac    = 0.0;

    //double residual_elim = 0.0;

    // perform time-integration specific actions
    if (is_stationary_)
    {
      // do not include any timefac for stationary calculations!
      timefacfac  = fac;
      timetaufac  = taufac;

      if (migrationinresidual_)
      {
        residual  = conv_eff_k - diff_ephinp_k + migrea_k - rhsint;
        //if (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
        //  residual_elim = (-valence_[k]/valence_[numscal_])*(conv_ephinp_k+diffusvalence_[numscal_]*(migconv_.Dot(ephinp_[k])) -((diffus_[numscal_]/diffus_[k])*diff_ephinp_k));
      }
      else
      {
        residual  = conv_ephinp_k - diff_ephinp_k - rhsint;
        //if (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
        //{
        //  residual_elim = (-valence_[k]/valence_[numscal_])*(conv_ephinp_k -((diffus_[numscal_]/diffus_[k])*diff_ephinp_k));
        //}
      }

      rhsfac      = fac;
      rhstaufac   = taufac;
    }
    else
    {
      timefacfac  = timefac * fac;
      timetaufac  = timefac * taufac;

      if (is_genalpha_)
      {
        // note: in hist_ we receive the time derivative phidtam at time t_{n+alpha_M} !!
        if (migrationinresidual_)
          residual  = hist_[k] + conv_eff_k - diff_ephinp_k + migrea_k - rhsint;
        else
          residual  = hist_[k] + conv_ephinp_k - diff_ephinp_k - rhsint;

        rhsfac    = timefacfac/alphaF;
        rhstaufac = timetaufac/alphaF;
        rhsint   *= (timefac/alphaF);  // not nice, but necessary !

        // rhs contribution due to incremental formulation (phidtam)
        // Standard Galerkin term
        const double vtrans = rhsfac*hist_[k];
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] -= vtrans*funct_(vi);
        }

        // ToDo: conservative form!!!!

      }
      else
      {
        rhsint = hist_[k] + (rhs_[k]*timefac); // contributions from t_n and \theta*dt*bodyforce(t_{n+1})

        if (migrationinresidual_)
          residual  = conint_[k] + timefac*(conv_eff_k - diff_ephinp_k + migrea_k) - rhsint;
        else
          residual  = conint_[k] + timefac*(conv_ephinp_k - diff_ephinp_k) - rhsint;

        rhsfac    = timefacfac;
        rhstaufac = taufac;

        // rhs contribution due to incremental formulation (phinp)
        // Standard Galerkin term
        const double vtrans = fac*conint_[k];
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] -= vtrans*funct_(vi);
        }


        // ToDo: conservative form!!!!

      } // if(is_genalpha_)

      //----------------------------------------------------------------
      // 1) element matrix: instationary terms
      //----------------------------------------------------------------
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;
        const double fac_funct_vi = fac*funct_(vi);

        // compute effective convective stabilization operator
        double conv_eff_vi = conv_(vi);
        if (migrationstab_)
        {
          conv_eff_vi += diffus_valence_k*migconv_(vi);
        }

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          /* Standard Galerkin term: */
          emat(fvi, fui) += fac_funct_vi*funct_(ui) ;

          /* 1) convective stabilization of transient term*/
          emat(fvi, fui) += taufac*conv_eff_vi*funct_(ui);

          /* 2) diffusive stabilization */
          // not implemented. Only stabilization of SUPG type

          /* 3) reactive stabilization (reactive part of migration term) */
          // not implemented. Only stabilization of SUPG type

        } // for ui
      } // for vi

    } // if (is_stationary_)

#ifdef PRINT_ELCH_DEBUG
    std::cout<<"tau["<<k<<"]    = "<<tau_[k]<<std::endl;
    std::cout<<"taufac["<<k<<"] = "<<taufac<<std::endl;
    if (tau_[k] != 0.0)
      std::cout<<"residual["<<k<<"] = "<< residual<<std::endl;
    std::cout<<"conv_eff_k    = "<<conv_eff_k<<std::endl;
    std::cout<<"conv_ephinp_k  = "<<conv_ephinp_k<<std::endl;
    std::cout<<"Dkzk_mig_ephinp_k = "<<Dkzk_mig_ephinp_k<<std::endl;
    std::cout<<"diff_ephinp_k = "<<diff_ephinp_k<<std::endl;
    std::cout<<"migrea_k      = "<<migrea_k <<std::endl;
    std::cout<<std::endl;
#endif

    // experimental code part
    if (betterconsistency_)
    {
      dserror("Has to be re-implemented!");
      //double fdiv(0.0); // we get the negative(!) reconstructed flux from outside!
      // compute divergence of approximated diffusive and migrative fluxes
      //GetDivergence(fdiv,efluxreconstr_[k],derxy_);
      //double taufacresidual = taufac*rhsint - timetaufac*(conv_ephinp_k + fdiv);
    } // betterconsistency

    //----------------------------------------------------------------
    // 2) element matrix: stationary terms
    //----------------------------------------------------------------
    for (int vi=0; vi<nen_; ++vi)
    {
      const int    fvi = vi*numdofpernode_+k;

      // compute effective convective stabilization operator
      double conv_eff_vi = conv_(vi);
      if (migrationstab_)
      {
        conv_eff_vi += diffus_valence_k*migconv_(vi);
      }

      const double timefacfac_funct_vi = timefacfac*funct_(vi);
      const double timefacfac_diffus_valence_k_mig_vi = timefacfac*diffus_valence_k*migconv_(vi);

      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+k;

        //----------------------------------------------------------------
        // standard Galerkin terms
        //----------------------------------------------------------------

        // matrix entries
        double matvalconc = 0.0;
        double matvalpot = 0.0;

        // convective term
        matvalconc += timefacfac_funct_vi*conv_(ui) ;

        // addition to convective term for conservative form
        if (is_conservative_)
        {
          // convective term using current scalar value
          matvalconc += timefacfac_funct_vi*vdiv_*funct_(ui);
        }

        // diffusive term
        double laplawf(0.0);
        GetLaplacianWeakForm(laplawf, derxy_,ui,vi); // compute once, reuse below!
        matvalconc += timefacfac*diffus_[k]*laplawf;

        // migration term
        // a) derivative w.r.t. concentration c_k
        matvalconc -= timefacfac_diffus_valence_k_mig_vi*funct_(ui);
        // b) derivative w.r.t. electric potential
        matvalpot += frt*timefacfac*diffus_valence_k*conint_[k]*laplawf;


        // TODO (ehrl)
        // Including stabilization yields in different results for the uncharged particle and
        // the binary electrolyte solution
        // -> Check calculation procedure of the method

        //----------------------------------------------------------------
        // Stabilization terms
        //----------------------------------------------------------------

        /* 0) transient stabilization */
        // not implemented. Only stabilization of SUPG type

        /* 1) convective stabilization */

        /* convective term */

        // I) linearization of residual part of stabilization term

        // effective convective stabilization of convective term
        // derivative of convective term in residual w.r.t. concentration c_k
        matvalconc += timetaufac*conv_eff_vi*conv_(ui);

        // migration convective stabilization of convective term
        double val_ui; GetLaplacianWeakFormRHS(val_ui, derxy_,gradphi_,ui);
        if (migrationinresidual_)
        {
          // a) derivative w.r.t. concentration_k
          matvalconc += timetaufac*conv_eff_vi*diffus_valence_k*migconv_(ui);

          // b) derivative w.r.t. electric potential
          matvalpot -= timetaufac*conv_eff_vi*diffus_valence_k*frt*val_ui;

          // note: higher-order and instationary parts of residuum part are linearized elsewhere!
        }

        // II) linearization of convective stabilization operator part of stabilization term
        if (migrationstab_)
        {
          // a) derivative w.r.t. concentration_k
          //    not necessary -> zero

          // b) derivative w.r.t. electric potential
          double laplacewf(0.0);
          GetLaplacianWeakForm(laplacewf, derxy_, ui,vi);
          matvalpot -= timetaufac*residual*diffus_valence_k*frt*laplacewf;
        }

        // III) linearization of tau part of stabilization term
        if (migrationintau_)
        {
          // derivative of tau (only effective for Taylor_Hughes_Zarins) w.r.t. electric potential
          const double tauderiv_ui = ((tauderpot_[k])(ui,0));
          matvalpot += timefacfac*tauderiv_ui*conv_eff_vi*residual;
        }

        // try to access the element matrix not too often. Can be costly
        emat(fvi,fui)                        += matvalconc;
        emat(fvi,ui*numdofpernode_+numscal_) += matvalpot;

      } // for ui

    } // for vi

    //-------------------------------------------------------------------------
    // 2b) element matrix: stationary terms (governing equation for potential)
    //-------------------------------------------------------------------------
    // what's the governing equation for the electric potential field?
    // we provide a lot of different options here:
    if (scatratype==INPAR::SCATRA::scatratype_elch_enc)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;
        const double alphaF_valence_k_fac_funct_vi = alphaF*valence_[k]*fac*funct_(vi);

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          // electroneutrality condition (only derivative w.r.t. concentration c_k)
          emat(pvi, fui) += alphaF_valence_k_fac_funct_vi*funct_(ui);
        } // for ui
      } // for vi
    }
    else if (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;
        const double timefacfac_diffus_valence_k_mig_vi = timefacfac*diffus_valence_k*migconv_(vi);

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          double laplawf(0.0);
          GetLaplacianWeakForm(laplawf, derxy_,ui,vi);

          // use 2nd order pde derived from electroneutrality condition (k=1,...,m)
          // a) derivative w.r.t. concentration c_k
          emat(pvi, fui) -= valence_[k]*(timefacfac_diffus_valence_k_mig_vi*funct_(ui));
          emat(pvi, fui) += valence_[k]*(timefacfac*diffus_[k]*laplawf);
          // b) derivative w.r.t. electric potential
          emat(pvi, ui*numdofpernode_+numscal_) += valence_[k]*(frt*timefacfac*diffus_valence_k*conint_[k]*laplawf);
        } // for ui
      } // for vi
    }
    else if (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;
        const double timefacfac_diffus_valence_k_mig_vi = timefacfac*diffus_valence_k*migconv_(vi);
        const double timefacfac_diffus_valence_m_mig_vi = timefacfac*diffus_[numscal_]*valence_[numscal_]*migconv_(vi);

        for (int ui=0; ui<nen_; ++ui)
        {
          // matrix entries
          double matvalconc = 0.0;
          double matvalpot = 0.0;

          double laplawf(0.0);
          GetLaplacianWeakForm(laplawf, derxy_,ui,vi);

          // use 2nd order pde derived from electroneutrality condition (k=1,...,m-1)
          // a) derivative w.r.t. concentration c_k
          matvalconc -= (timefacfac_diffus_valence_k_mig_vi*funct_(ui));
          matvalconc += (timefacfac*diffus_[k]*laplawf);
          // b) derivative w.r.t. electric potential
          matvalpot += (frt*timefacfac*diffus_valence_k*conint_[k]*laplawf);

          // care for eliminated species with index m
          //(diffus_ and valence_ vector were extended in GetMaterialParams()!)
          // a) derivative w.r.t. concentration c_k
          matvalconc += (timefacfac_diffus_valence_m_mig_vi*funct_(ui));
          matvalconc -= (timefacfac*diffus_[numscal_]*laplawf);
          // b) derivative w.r.t. electric potential
          matvalpot -= (frt*timefacfac*diffus_[numscal_]*valence_[numscal_]*conint_[k]*laplawf);

          // try to access the element matrix not too often. Can be costly
          const int fui = ui*numdofpernode_+k;
          emat(pvi,fui) += valence_[k]*matvalconc;
          const int pui = ui*numdofpernode_+numscal_;
          emat(pvi,pui) += valence_[k]*matvalpot;

        } // for ui
      } // for vi
    }
    else if (scatratype==INPAR::SCATRA::scatratype_elch_poisson)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;
        const double alphaF_valence_k_fac_funct_vi = alphaF*valence_[k]*fac*funct_(vi);

        for (int ui=0; ui<nen_; ++ui)
        {
          // we have a loop over k around. So prevent that the potential
          // term is added more than once!!
          if (k==0)
          {
            const int pui = ui*numdofpernode_+numscal_;
            double laplawf(0.0);
            GetLaplacianWeakForm(laplawf, derxy_,ui,vi);

            const double epsbyF = epsilon/faraday;
            emat(pvi,pui) += alphaF*fac*epsbyF*laplawf;
          }
          const int fui = ui*numdofpernode_+k;
          // electroneutrality condition (only derivative w.r.t. concentration c_k)
          emat(pvi, fui) += alphaF_valence_k_fac_funct_vi*funct_(ui);
        } // for ui
      } // for vi
    }
    else if (scatratype==INPAR::SCATRA::scatratype_elch_laplace)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;

        for (int ui=0; ui<nen_; ++ui)
        {
          // we have a loop over k around. So prevent that the potential
          // term is added more than once!!
          if (k==0)
          {
            const int pui = ui*numdofpernode_+numscal_;
            double laplawf(0.0);
            GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
            emat(pvi,pui) += alphaF*fac*laplawf;
          }
        } // for ui
      } // for vi
    }
    else
      dserror ("How did you reach this point?");


    if (use2ndderiv_)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        // compute effective convective stabilization operator
        double conv_eff_vi = conv_(vi);
        if (migrationstab_)
        {
          conv_eff_vi += diffus_valence_k*migconv_(vi);
        }

        const double timetaufac_conv_eff_vi = timetaufac*conv_eff_vi;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          // 1) convective stabilization

          // diffusive term
          // derivative w.r.t. concentration c_k
          emat(fvi,fui) -= timetaufac_conv_eff_vi*diff_(ui) ;

        } // for ui

        // reactive part of migration term
        if (migrationinresidual_)
        {
          const double timetaufac_conv_eff_vi_conint_k_frt_valence_k =timetaufac_conv_eff_vi*conint_[k]*frt*valence_[k];
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            // a) derivative w.r.t. concentration_k
            emat(fvi,fui) += timetaufac_conv_eff_vi*migrea_(ui) ;
            // note: migrea_ already contains frt*diffus_valence!!!

            // b) derivative w.r.t. electric potential
            emat(fvi, ui*numdofpernode_+numscal_) -= timetaufac_conv_eff_vi_conint_k_frt_valence_k*diff_(ui);
            // note: diff_ already includes factor D_k

          } // for ui
        }

        // 2) diffusive stabilization
        // not implemented. Only stabilization of SUPG type

        // 3) reactive stabilization (reactive part of migration term)
        // not implemented. Only stabilization of SUPG type

      } // for vi
    } // use2ndderiv


    //-----------------------------------------------------------------------
    // 3) element right hand side vector (neg. residual of nonlinear problem)
    //-----------------------------------------------------------------------
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi = vi*numdofpernode_+k;

      //----------------------------------------------------------------
      // standard Galerkin terms (ion transport equations)
      //----------------------------------------------------------------

      // RHS source term (contains old part of rhs for OST / BDF2)
      erhs[fvi] += fac*funct_(vi)*rhsint ;

      // nonlinear migration term
      erhs[fvi] += rhsfac*conint_[k]*diffus_valence_k*migconv_(vi);

      // convective term
      erhs[fvi] -= rhsfac*funct_(vi)*conv_ephinp_k;

      // addition to convective term for conservative form
      // (not included in residual)
      if (is_conservative_)
      {
        // convective term in conservative form
        erhs[fvi] -= rhsfac*funct_(vi)*conint_[k]*vdiv_;
      }

      // diffusive term
      double laplawf(0.0);
      GetLaplacianWeakFormRHS(laplawf,derxy_,gradphi_,vi);
      erhs[fvi] -= rhsfac*diffus_[k]*laplawf;


      //----------------------------------------------------------------
      // Stabilization terms
      //----------------------------------------------------------------

      // 0) transient stabilization
      //    not implemented. Only stabilization of SUPG type

      // 1) convective stabilization

      erhs[fvi] -= rhstaufac*conv_(vi)*residual;
      if (migrationstab_)
      {
        erhs[fvi] -=  rhstaufac*diffus_valence_k*migconv_(vi)*residual;
      }

      // 2) diffusive stabilization
      //    not implemented. Only stabilization of SUPG type

      // 3) reactive stabilization (reactive part of migration term)
      //    not implemented. Only stabilization of SUPG type

    } // for vi

      //----------------------------------------------------------------
      // standard Galerkin terms (equation for electric potential)
      //----------------------------------------------------------------
      // what's the governing equation for the electric potential field ?
    if (scatratype==INPAR::SCATRA::scatratype_elch_enc)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;

        // electroneutrality condition
        // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
        erhs[pvi] -= valence_[k]*fac*funct_(vi)*conint_[k];
      } // for vi
    }
    else if (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;

        double laplawf(0.0);
        GetLaplacianWeakFormRHS(laplawf,derxy_,gradphi_,vi);

        // use 2nd order pde derived from electroneutrality condition (k=1,...,m)
        erhs[pvi] += rhsfac*valence_[k]*((diffus_valence_k*conint_[k]*migconv_(vi))-(diffus_[k]*laplawf));
      } // for vi
    }
    else if (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;

        double laplawf(0.0);
        GetLaplacianWeakFormRHS(laplawf,derxy_,gradphi_,vi);

        // use 2nd order pde derived from electroneutrality condition (k=0,...,m-1)
        erhs[pvi] += rhsfac*valence_[k]*((diffus_valence_k*conint_[k]*migconv_(vi))-(diffus_[k]*laplawf));
        // care for eliminated species with index m
        //(diffus_ and valence_ vector were extended in GetMaterialParams()!)
        erhs[pvi] -= rhsfac*valence_[k]*((diffus_[numscal_]*valence_[numscal_]*conint_[k]*migconv_(vi))-(diffus_[numscal_]*laplawf));

      } // for vi
    }
    else if (scatratype==INPAR::SCATRA::scatratype_elch_poisson)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;

        // we have a loop over k around. So prevent that the potential
        // term is added more than once!!
        if (k==0)
        {
          double laplawf(0.0);
          GetLaplacianWeakFormRHS(laplawf,derxy_,gradpot_,vi);
          const double epsbyF = epsilon/faraday;
          erhs[pvi] -= fac*epsbyF*laplawf;
        }

        // electroneutrality condition
        // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
        erhs[pvi] -= valence_[k]*fac*funct_(vi)*conint_[k];

      } // for vi
    }
    else if (scatratype==INPAR::SCATRA::scatratype_elch_laplace)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        const int pvi = vi*numdofpernode_+numscal_;

        // we have a loop over k around. So prevent that the potential
        // term is added more than once!!
        if (k==0)
        {
          double laplawf(0.0);
          GetLaplacianWeakFormRHS(laplawf,derxy_,gradpot_,vi);
          erhs[pvi] -= fac*laplawf;

        }

      } // for vi
    }
    else
      dserror ("How did you reach this point?");

    // RHS vector finished


  } // loop over scalars

  return;
} // ScaTraImpl::CalMatElch


/*----------------------------------------------------------------------*
  |  Calculate conductivity (ELCH) (private)                   gjb 07/09 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalculateConductivity(
  const DRT::Element*  ele,
  const double         frt,
  const enum INPAR::SCATRA::ScaTraType  scatratype,
  Epetra_SerialDenseVector& sigma
  )
{
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_tau(SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // evaluate shape functions (and not needed derivatives) at element center
  EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0,ele->Id());

  // get concentration of transported scalar k at integration point
  for (int k = 0;k<numscal_;++k)
    conint_[k] = funct_.Dot(ephinp_[k]);

  GetMaterialParams(ele,scatratype,0.0); // use dt=0.0 dymmy value

  // compute the conductivity (1/(\Omega m) = 1 Siemens / m)
  double sigma_all(0.0);
  const double factor = frt*INPAR::SCATRA::faraday_const; // = F^2/RT

  // Dilute solution theory:
  // Conductivity is computed by
  // sigma = F^2/RT*Sum(z_k^2 D_k c_k)
  if(not diffcond_)
  {
    for(int k=0; k < numscal_; k++)
    {
      double sigma_k = factor*valence_[k]*diffusvalence_[k]*conint_[k];
      sigma[k] += sigma_k; // insert value for this ionic species
      sigma_all += sigma_k;

      // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
      if(scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
      {
        sigma_all += factor*diffusvalence_[numscal_]*valence_[k]*(-conint_[k]);
      }
    }
  }
  // Concentrated solution theory:
  // Conductivity given by a function is evaluated at bulk concentration
  else
  {
    Teuchos::RCP<MAT::Material> material = ele->Material();
    const Teuchos::RCP<const MAT::ElchMat>& actmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(material);
    // loop over single phases
    for (int iphase=0; iphase < actmat->NumPhase();++iphase)
    {
      const int phaseid = actmat->PhaseID(iphase);
      Teuchos::RCP<const MAT::Material> singlemat = actmat->PhaseById(phaseid);

      if(singlemat->MaterialType() == INPAR::MAT::m_elchphase)
      {
        const MAT::ElchPhase* actsinglemat = static_cast<const MAT::ElchPhase*>(singlemat.get());
        sigma_all = actsinglemat->ComputeConductivity(conint_[0]);
      }
      else
        dserror("Conductivity has to be defined in m_elchphase! There is no material m_elchphase");
    }
  }
  // conductivity based on ALL ionic species (even eliminated ones!)
  sigma[numscal_] += sigma_all;

  return;

} //ScaTraImpl<distype>::CalculateConductivity


/*----------------------------------------------------------------------*
  |  CalculateElectricPotentialField (ELCH) (private)          gjb 04/10 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalculateElectricPotentialField(
  const DRT::Element*         ele,
  const double                frt,
  const enum INPAR::SCATRA::ScaTraType  scatratype,
  Epetra_SerialDenseMatrix&   emat,
  Epetra_SerialDenseVector&   erhs
  )
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    // get concentration of transported scalar k at integration point
    for (int k = 0;k<numscal_;++k)
      conint_[k] = funct_.Dot(ephinp_[k]);

    // get gradient of electric potential at integration point
    gradpot_.Multiply(derxy_,epotnp_);

    // access material parameters
    GetMaterialParams(ele,scatratype,0.0); // use dt=0.0 dymmy value

    double sigmaint(0.0);
    for (int k=0; k<numscal_; ++k)
    {
      double sigma_k = frt*valence_[k]*diffusvalence_[k]*conint_[k];
      sigmaint += sigma_k;

      // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
      if(scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
        sigmaint += frt*valence_[k]*diffusvalence_[numscal_]*(-conint_[k]);

      // diffusive terms on rhs
      // gradient of current scalar value
      gradphi_.Multiply(derxy_,ephinp_[k]);
      const double vrhs = fac*diffusvalence_[k];
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = vi*numdofpernode_+numscal_;
        double laplawf(0.0);
        GetLaplacianWeakFormRHS(laplawf,derxy_,gradphi_,vi);
        erhs[fvi] -= vrhs*laplawf;
        // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
        if(scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
          erhs[fvi] -= -fac*valence_[k]*diffus_[numscal_]*laplawf;
      }

      // provide something for conc. dofs: a standard mass matrix
      for (int vi=0; vi<nen_; ++vi)
      {
        const int    fvi = vi*numdofpernode_+k;
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui = ui*numdofpernode_+k;
          emat(fvi,fui) += fac*funct_(vi)*funct_(ui);
        }
      }
    } // for k

    // ----------------------------------------matrix entries
    for (int vi=0; vi<nen_; ++vi)
    {
      const int    fvi = vi*numdofpernode_+numscal_;
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui = ui*numdofpernode_+numscal_;
        double laplawf(0.0);
        GetLaplacianWeakForm(laplawf, derxy_,ui,vi);
        emat(fvi,fui) += fac*sigmaint*laplawf;
      }

      double laplawf(0.0);
      GetLaplacianWeakFormRHS(laplawf,derxy_,gradpot_,vi);
      erhs[fvi] -= fac*sigmaint*laplawf;
    }
  } // integration loop

  return;

} //ScaTraImpl<distype>::CalculateElectricPotentialField

/*------------------------------------------------------------------------*
  |  calculate residual of scalar transport equation for the homogenized  |
  |  transport equation in poroelastic problem. (depending on respective  |
  |  stationary or time-integration scheme)                vuong 04/12  |
  *-----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalcResidual_PoroScatraMod(
  const double   dt,
  const double   timefac,
  const int      k,
  const double   porosity,
  const double   dporodt,
  LINALG::Matrix<3,1>& gradporosity
  )
{
  dserror("CalcResidual_PoroScatraMod not implemented");

  return;
} //end of CalcResidual_Poroscatra


/*---------------------------------------------------------------------------*
 |  modify element matrix and rhs for scatra in porous media (private)  vuong 06/12|
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraImpl<distype>::CalMatAndRHS_PoroScatraMod(
  Epetra_SerialDenseMatrix&             emat,
  Epetra_SerialDenseVector&             erhs,
  const double                          fac,      ///< integration factor
  const double                          timefac,  ///< time discretization factor
  const int                             k,
  const int                             eleid,
  const int                             iquad
  )
{
  //access structure discretization
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  //get corresponding structure element (it has the same global ID as the scatra element)
  DRT::Element* structele = structdis->gElement(eleid);
  if (structele == NULL)
    dserror("Structure element %i not on local processor", eleid);

  const Teuchos::RCP<const MAT::StructPoro>& structmat
            = Teuchos::rcp_dynamic_cast<const MAT::StructPoro>(structele->Material());
  if(structmat == Teuchos::null)
    dserror("invalid structure material for poroelasticity");

  const double           porosity   = structmat->GetPorosityAtGP(iquad);
  //const double           dporodt    = structmat ->GetDPoroDtAtGP(iquad);
 // LINALG::Matrix<3,1>  gradporosity = structmat->GetGradPorosityAtGP(iquad);

  //const double timefacfac = timefac * fac;
  
  //----------------------------------------------------------------
  // 1) Modification of emat due to the homogenized equation employed for
  //    the poro-scatra problem.The standard equation is multiplied by the
  //    porosity, and some other terms must be added.
  //----------------------------------------------------------------

  emat.Scale(porosity);

  /*
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = timefacfac*funct_(vi);
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui = ui*numdofpernode_+k;
      emat(fvi,fui) += v*dporodt*funct_(ui);

      double tmp=0.0;
      for(int i = 0; i<nsd_; i++)
      {
        tmp += v*funct_(ui)*convelint_(i,0)*gradporosity(i);
        tmp -= v*diffus_[k]*(derxy_(i,ui)*gradporosity(i));
      }
      emat(fvi,fui) += tmp;
    }
  }
  */

  //----------------------------------------------------------------
  // 2) Modification of the residual due to the homogenized equation employed for
  //    the poro-scatra problem.The standard equation is multiplied by the
  //    porosity, and some other terms must be added.
  //----------------------------------------------------------------

  erhs.Scale(porosity);

  /*
  // compute scalar at integration point
  const double phi = funct_.Dot(ephinp_[k]);

  double tmp = 0.0;
  for (int i=0; i<nsd_; i++)
  {
    tmp += phi*convelint_(i,0)*(gradporosity(i)) - diffus_[k]*gradphi_(i,0)*gradporosity(i);
  }
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+k;
    erhs[fvi] -= funct_(vi)* timefacfac*( phi*dporodt + tmp);
  }
  */

  return;
} //ScaTraImpl::CalMatAndRHS_Poroscatra


