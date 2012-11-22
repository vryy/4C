/*----------------------------------------------------------------------*/
/*!
\file thermo_ele_impl.cpp

\brief Internal implementation of thermo elements

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 08/09 |
 *----------------------------------------------------------------------*/
#ifdef D_THERMO

/*----------------------------------------------------------------------*
 | headers                                                   dano 08/09 |
 *----------------------------------------------------------------------*/
#include "../drt_inpar/inpar_thermo.H"
#include "../drt_inpar/inpar_structure.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

// material headers
#include "../drt_mat/fourieriso.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/thermoplasticlinelast.H"

#include "thermo_element.H" // only for visualization of element data
#include "thermo_ele_impl.H"

#include "../drt_tsi/tsi_defines.H"


/*----------------------------------------------------------------------*
 |                                                           dano 09/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TemperImplInterface* DRT::ELEMENTS::TemperImplInterface::Impl(
  DRT::Element* ele
  )
{
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    return TemperImpl<DRT::Element::hex8>::Instance();
  }
  case DRT::Element::hex20:
  {
    return TemperImpl<DRT::Element::hex20>::Instance();
  }
  case DRT::Element::hex27:
  {
    return TemperImpl<DRT::Element::hex27>::Instance();
  }
  case DRT::Element::tet4:
  {
    return TemperImpl<DRT::Element::tet4>::Instance();
  }
 /* case DRT::Element::tet10:
  {
    return TemperImpl<DRT::Element::tet10>::Instance();
  } */
  case DRT::Element::wedge6:
  {
    return TemperImpl<DRT::Element::wedge6>::Instance();
  }
/*  case DRT::Element::wedge15:
  {
    return TemperImpl<DRT::Element::wedge15>::Instance();
  } */
  case DRT::Element::pyramid5:
  {
    return TemperImpl<DRT::Element::pyramid5>::Instance();
  }
  case DRT::Element::quad4:
  {
    return TemperImpl<DRT::Element::quad4>::Instance();
  }
/*  case DRT::Element::quad8:
  {
    return TemperImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return TemperImpl<DRT::Element::quad9>::Instance();
  }*/
  case DRT::Element::tri3:
  {
    return TemperImpl<DRT::Element::tri3>::Instance();
  }
/*  case DRT::Element::tri6:
  {
    return TemperImpl<DRT::Element::tri6>::Instance();
  }*/
  case DRT::Element::line2:
  {
    return TemperImpl<DRT::Element::line2>::Instance();
  }/*
  case DRT::Element::line3:
  {
    return TemperImpl<DRT::Element::line3>::Instance();
  }*/
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
  }
  return NULL;

}  // TemperImperInterface::Impl()


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TemperImpl<distype> * DRT::ELEMENTS::TemperImpl<distype>::Instance(
  bool create
  )
{
  static TemperImpl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new TemperImpl<distype>();
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
void DRT::ELEMENTS::TemperImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 | initialisation of the data with respect to the declaration           |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TemperImpl<distype>::TemperImpl()
: etemp_(false),
  ecapa_(true),
  xyze_(true),
  radiation_(false),
  xsi_(true),
  funct_(true),
  deriv_(true),
  xjm_(true),
  xij_(true),
  derxy_(true),
  fac_(0.0),
  gradtemp_(true),
  heatflux_(false),
  cmat_(false),
  plasticmat_(false)
{
  return;
}


/*----------------------------------------------------------------------*
 | evaluate for multiple dofsets                             dano 02/10 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TemperImpl<distype>::Evaluate(
  DRT::Element* ele,
  Teuchos::ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Element::LocationArray& la,
  Epetra_SerialDenseMatrix& elemat1_epetra,  // Tangent ("stiffness")
  Epetra_SerialDenseMatrix& elemat2_epetra,  // Capacity ("mass")
  Epetra_SerialDenseVector& elevec1_epetra,  // internal force vector
  Epetra_SerialDenseVector& elevec2_epetra,  // external force vector
  Epetra_SerialDenseVector& elevec3_epetra  // capacity vector
  )
{
  // what actions are available
  // (action == "calc_thermo_fint")
  // (action == "calc_thermo_fintcapa")
  // (action == "calc_thermo_finttang")
  // (action == "proc_thermo_heatflux")
  // (action == "postproc_thermo_heatflux")
  // (action == "integrate_shape_functions")
  // (action == "calc_thermo_update_istep")
  // (action == "calc_thermo_reset_istep")
  // (action == "calc_thermo_energy")
  // (action == "calc_thermo_coupltang")

  // check length
  if (la[0].Size() != nen_*numdofpernode_)
    dserror("Location vector length does not match!");

  // disassemble temperature
  if ( discretization.HasState(0, "temperature") )
  {
    std::vector<double> mytempnp((la[0].lm_).size());
    Teuchos::RCP<const Epetra_Vector> tempnp
      = discretization.GetState(0,"temperature");
    if (tempnp == Teuchos::null)
      dserror("Cannot get state vector 'tempnp'");
    DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[0].lm_);
    // build the element temperature
    LINALG::Matrix<nen_*numdofpernode_,1> etemp(&(mytempnp[0]),true);  // view only!
    etemp_.Update(etemp);  // copy
  }
  // initialize capacity matrix
  ecapa_.Clear();
  // check for the action parameter
  const std::string action = params.get<std::string>("action","none");

  double time = 0.0;

  if(action != "calc_thermo_energy")
  {
    // extract time
    time = params.get<double>("total time");
  }

  // if ele is a thermo element --> the THR element method KinType() exists
  DRT::ELEMENTS::Thermo* therm = dynamic_cast<DRT::ELEMENTS::Thermo*>(ele);
  // kintype = 0: geo_linear or purely thermal problem
  // kintype = 1: geo_nonlinear
  int kintype = therm->KinType();

  //============================================================================
  // calculate tangent K and internal force F_int = K * Theta
  // --> for static case
  if (action == "calc_thermo_fintcond")
  {
    // set views
    LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_> etang(elemat1_epetra.A(),true);  // view only!
    LINALG::Matrix<nen_*numdofpernode_,1> efint(elevec1_epetra.A(),true);  // view only!
    // ecapa, efext, efcap not needed for this action

    // default: purely thermal / geometrically linear TSI problem
    if (kintype == 0)
    {
      CalculateFintCondCapa(
        ele,
        time,
        &etang,
        NULL,
        &efint,
        NULL,
        NULL,
        NULL
        );
    }

    // initialise the vectors
    // Evaluate() is called the first time in ThermoBaseAlgorithm: at this stage
    // the coupling field is not yet known. Pass coupling vectors filled with zeros
    // the size of the vectors is the length of the location vector*nsd_
    std::vector<double> mydisp( ( (la[0].lm_).size() )*nsd_, 0.0 );
    std::vector<double> myvel( ( (la[0].lm_).size() )*nsd_, 0.0 );

    // if it's a TSI problem with displacementcoupling_ --> go on here!
    if (la.Size()>1)
    {
      // and now get the current displacements/velocities
      if ( (discretization.HasState(1,"displacement")) ||
           (discretization.HasState(1,"velocity"))
         )
      {
        // get the displacements
        Teuchos::RCP<const Epetra_Vector> disp
          = discretization.GetState(1,"displacement");
        if (disp==Teuchos::null)
          dserror("Cannot get state vectors 'displacement'");
        // extract the displacements
        DRT::UTILS::ExtractMyValues(*disp,mydisp,la[1].lm_);

        // get the velocities
        Teuchos::RCP<const Epetra_Vector> vel
          = discretization.GetState(1,"velocity");
        if (vel==Teuchos::null)
          dserror("Cannot get state vectors 'velocity'");
        // extract the displacements
        DRT::UTILS::ExtractMyValues(*vel,myvel,la[1].lm_);

        // if there is a displacement/velocity vector
        // --> calculate coupling term

#ifndef MonTSIwithoutSTR

        // access the structure discretization, needed later for calling the solid
        // material and getting its tangent
        Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
        structdis = DRT::Problem::Instance()->GetDis("structure");

        // get GID of the first solid element (by using an homogenous material this
        // is enough)
        // ask the partner element about his Id
        const int structgid = ele->Id();

        // get the pointer to the adequate structure element based on GIDs
        DRT::Element* structele = structdis->gElement(structgid);

        // call ThermoStVenantKirchhoff material and get the temperature dependent
        // tangent ctemp
        Teuchos::RCP<MAT::Material> structmat = structele->Material();
        if (structmat->MaterialType() == INPAR::MAT::m_thermopllinelast)
          plasticmat_ = true;
      }  // disp!=0 & vel!=0
    }  // la.Size>1

    // get the time step size
    const double stepsize = params.get<double>("delta time");

    // geometrically linear TSI problem
    if ( (kintype == 0) && (la.Size()>1) )
    {
      CalculateCouplFintCondCapa(
        ele,
        time,
        mydisp,
        myvel,
        stepsize,
        &etang,
        &efint,
        NULL,
        NULL,
        NULL
        );

      if (plasticmat_)
        CalculateInternalDissipation(
          ele,
          myvel,
          stepsize,
          &etang,
          &efint
          );
    }
    // geometrically nonlinear TSI problem
    else if (kintype == 1)
    {
      CalculateNlnCouplFintCondCapa(
        ele,
        time,
        mydisp,
        myvel,
        stepsize,
        &etang,
        NULL,
        &efint,
        NULL,
        NULL,
        INPAR::THR::heatflux_none,
        INPAR::THR::tempgrad_none
        );
    }
#endif // MonTSIwithoutSTR

  }  // action == "calc_thermo_fintcond"

  //============================================================================
  // calculate only the internal force F_int, needed for restart
  else if (action == "calc_thermo_fint")
  {
    // set views
    LINALG::Matrix<nen_*numdofpernode_,1> efint(elevec1_epetra.A(),true);  // view only!
    // etang, ecapa, efext, efcap not needed for this action

    // purely thermal / geometrically linear TSI problem
    if (kintype == 0)
    {
      CalculateFintCondCapa(
        ele,
        time,
        NULL,
        NULL,
        &efint,
        NULL,
        NULL,
        NULL
        );
    }

    // initialise the vectors
    // Evaluate() is called the first time in ThermoBaseAlgorithm: at this stage
    // the coupling field is not yet known. Pass coupling vectors filled with zeros
    // the size of the vectors is the length of the location vector*nsd_
    std::vector<double> mydisp( ( (la[0].lm_).size() )*nsd_, 0.0 );
    std::vector<double> myvel( ( (la[0].lm_).size() )*nsd_, 0.0 );

    // if it's a TSI problem with displacementcoupling_ --> go on here!
    if (la.Size()>1)
    {
      // and now get the current displacements/velocities
      if ( (discretization.HasState(1,"displacement")) ||
           (discretization.HasState(1,"velocity"))
         )
      {
        // get the displacements
        Teuchos::RCP<const Epetra_Vector> disp
          = discretization.GetState(1,"displacement");
        if (disp==Teuchos::null)
          dserror("Cannot get state vectors 'displacement'");
        // extract the displacements
        DRT::UTILS::ExtractMyValues(*disp,mydisp,la[1].lm_);

        // get the velocities
        Teuchos::RCP<const Epetra_Vector> vel
          = discretization.GetState(1,"velocity");
        if (vel==Teuchos::null)
          dserror("Cannot get state vectors 'velocity'");
        // extract the velocities
        DRT::UTILS::ExtractMyValues(*vel,myvel,la[1].lm_);

        // if there is a displacement/velocity vector available go on here
        // --> calculate coupling

#ifndef MonTSIwithoutSTR

        // access the structure discretization, needed later for calling the solid
        // material and getting its tangent
        Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
        structdis = DRT::Problem::Instance()->GetDis("structure");

        // get GID of the first solid element (by using an homogenous material this
        // is enough)
        // ask the partner element about his Id
        const int structgid = ele->Id();

        // get the pointer to the adequate structure element based on GIDs
        DRT::Element* structele = structdis->gElement(structgid);

        // call ThermoStVenantKirchhoff material and get the temperature dependent
        // tangent ctemp
        Teuchos::RCP<MAT::Material> structmat = structele->Material();
        if (structmat->MaterialType() == INPAR::MAT::m_thermopllinelast)
          plasticmat_ = true;
      }  // disp!=0 & vel!=0
    }  // end la.Size>1

    // get the time step size
    const double stepsize = params.get<double>("delta time");

    // geometrically linear TSI problem
    if ( (kintype == 0) && (la.Size()>1) )
    {
      CalculateCouplFintCondCapa(
        ele,
        time,
        mydisp,
        myvel,
        stepsize,
        NULL,
        &efint,
        NULL,
        NULL,
        NULL
        );

      if (plasticmat_)
        CalculateInternalDissipation(
          ele,
          myvel,
          stepsize,
          NULL,
          &efint
          );
    }
    // geometrically nonlinear TSI problem
    else if (kintype == 1)
    {
      CalculateNlnCouplFintCondCapa(
        ele,  // current element
        time,  // current time
        mydisp,  // current displacements
        myvel,  // current velocities
        stepsize,  // time increment
        NULL,  // element conductivity matrix
        NULL,  // element capacity matrix
        &efint,  // element internal force vector
        NULL,  // heat flux at GP
        NULL,  // temperature gradients at GP
        INPAR::THR::heatflux_none,  // output option for q
        INPAR::THR::tempgrad_none  // output option for grad T
        );
    }

#endif // MonTSIwithoutSTR

  }  // action == "calc_thermo_fint"

  //============================================================================
  // calculate the capacity matrix and the internal force F_int
  // --> for dynamic case, called only once in DetermineCapaConsistTempRate()
  else if (action == "calc_thermo_fintcapa")
  {
    // set views
    LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_> ecapa(elemat2_epetra.A(),true);  // view only!
    LINALG::Matrix<nen_*numdofpernode_,1> efint(elevec1_epetra.A(),true);  // view only!
    // etang, efext, efcap not needed for this action

    // default: purely thermal / geometrically linear TSI problem
    if (kintype == 0)
    {
      CalculateFintCondCapa(
        ele,
        time,
        NULL,
        &ecapa,  // delivers capacity matrix
        &efint,
        NULL,
        NULL,
        NULL
        );
    }

    // initialise the vectors
    // Evaluate() is called the first time in ThermoBaseAlgorithm: at this stage
    // the coupling field is not yet known. Pass coupling vectors filled with zeros
    // the size of the vectors is the length of the location vector*nsd_
    std::vector<double> mydisp( ( (la[0].lm_).size() )*nsd_, 0.0 );
    std::vector<double> myvel( ( (la[0].lm_).size() )*nsd_, 0.0 );

    // if it's a TSI problem with displacementcoupling_ --> go on here!

    // and now get the current displacements/velocities
    if (la.Size()>1)
    {
      // and now get the current displacements/velocities
      if ( (discretization.HasState(1,"displacement")) ||
           (discretization.HasState(1,"velocity"))
         )
      {
        // get the displacements
        Teuchos::RCP<const Epetra_Vector> disp
          = discretization.GetState(1,"displacement");
        if (disp==Teuchos::null)
          dserror("Cannot get state vectors 'displacement'");
        // extract the displacements
        DRT::UTILS::ExtractMyValues(*disp,mydisp,la[1].lm_);

        // get the velocities
        Teuchos::RCP<const Epetra_Vector> vel
          = discretization.GetState(1,"velocity");
        if (vel==Teuchos::null)
          dserror("Cannot get state vectors 'velocity'");
        // extract the velocities
        DRT::UTILS::ExtractMyValues(*vel,myvel,la[1].lm_);

        // if there is a strucutural vector available go on here
        // --> calculate coupling

#ifndef MonTSIwithoutSTR
        // access the structure discretization, needed later for calling the solid
        // material and getting its tangent
        Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
        structdis = DRT::Problem::Instance()->GetDis("structure");

        // get GID of the first solid element (by using an homogenous material this
        // is enough)
        // ask the partner element about his Id
        const int structgid = ele->Id();

        // get the pointer to the adequate structure element based on GIDs
        DRT::Element* structele = structdis->gElement(structgid);

        // call ThermoStVenantKirchhoff material and get the temperature dependent
        // tangent ctemp
        Teuchos::RCP<MAT::Material> structmat = structele->Material();
        if (structmat->MaterialType() == INPAR::MAT::m_thermopllinelast)
          plasticmat_ = true;

#endif // MonTSIwithoutSTR

      }  // disp!=0 & vel!=0
    }  // la.Size>1

    // get the time step size
    const double stepsize = params.get<double>("delta time");

    // geometrically linear TSI problem
    if ( (kintype == 0) && (la.Size()>1) )
    {
      CalculateCouplFintCondCapa(
        ele,
        time,
        mydisp,
        myvel,
        stepsize,
        NULL,
        &efint,
        NULL,
        NULL,
        NULL
        );

      if (plasticmat_)
        CalculateInternalDissipation(
          ele,
          myvel,
          stepsize,
          NULL,
          &efint
          );
    }  // end geo_linear TSI

    // geometrically nonlinear TSI problem
    else if (kintype == 1)
    {
      CalculateNlnCouplFintCondCapa(
        ele,
        time,
        mydisp,
        myvel,
        stepsize,
        NULL,
        &ecapa,  // element capacity matrix
        &efint,  // element internal force vector
        NULL,
        NULL,
        INPAR::THR::heatflux_none,
        INPAR::THR::tempgrad_none
        );
    }

    // lumping
    if(params.get<bool>("lump capa matrix"))
    {
      const INPAR::THR::DynamicType timint
        = DRT::INPUT::get<INPAR::THR::DynamicType>(params, "time integrator",INPAR::THR::dyna_undefined);
      switch (timint)
      {
        case INPAR::THR::dyna_expleuler :
        {
          CalculateLumpMatrix(&ecapa);
          break;
        }
        case INPAR::THR::dyna_genalpha :
        case INPAR::THR::dyna_onesteptheta :
        case INPAR::THR::dyna_statics :
        {
          dserror("Lumped capacity matrix has not yet been tested");
          break;
        }
        case INPAR::THR::dyna_undefined :
        default :
        {
          dserror("Undefined time integration scheme for thermal problem!");
          break;
        }
      }  // end of switch(timint)
    }  // end of lumping
  }  // action == "calc_thermo_fintcapa"

  //============================================================================
  // called from overloaded function ApplyForceTangInternal(), exclusively for OST-timint
  // calculate tangent matrix K and consistent capacity matrix C
  // --> for dynamic case
  else if (action == "calc_thermo_finttang")
  {
    // set views
    LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_> etang(elemat1_epetra.A(),true);  // view only!
    LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_> ecapa(elemat2_epetra.A(),true);  // view only!
    LINALG::Matrix<nen_*numdofpernode_,1> efint(elevec1_epetra.A(),true);  // view only!
    LINALG::Matrix<nen_*numdofpernode_,1> efcap(elevec3_epetra.A(),true);  // view only!
    // efext not needed for this action

    // purely thermal / geometrically linear TSI problem
    if (kintype == 0)
    {
      CalculateFintCondCapa(
        ele,
        time,
        &etang,
        &ecapa_,
        &efint,
        NULL,
        NULL,
        NULL
        );
    }

    // initialise the vectors
    // Evaluate() is called the first time in ThermoBaseAlgorithm: at this stage
    // the coupling field is not yet known. Pass coupling vectors filled with zeros
    // the size of the vectors is the length of the location vector*nsd_
    std::vector<double> mydisp( ( (la[0].lm_).size() )*nsd_, 0.0 );
    std::vector<double> myvel( ( (la[0].lm_).size() )*nsd_, 0.0 );

    // if it's a TSI problem and if there are current displacements/velocities
    if (la.Size()>1)
    {
      if ( (discretization.HasState(1,"displacement")) &&
         (discretization.HasState(1,"velocity")) )
      {
        // get the displacements
        Teuchos::RCP<const Epetra_Vector> disp
          = discretization.GetState(1,"displacement");
        if (disp==Teuchos::null)
          dserror("Cannot get state vectors 'displacement'");
        // extract the displacements
        DRT::UTILS::ExtractMyValues(*disp,mydisp,la[1].lm_);

        // get the velocities
        Teuchos::RCP<const Epetra_Vector> vel
          = discretization.GetState(1,"velocity");
        if (vel==Teuchos::null)
          dserror("Cannot get state vectors 'velocity'");
        // extract the velocities
        DRT::UTILS::ExtractMyValues(*vel,myvel,la[1].lm_);

        // if there is a strucutural vector available go on here
        // --> calculate coupling

  #ifndef MonTSIwithoutSTR
        // access the structure discretization, needed later for calling the solid
        // material and getting its tangent
        Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
        structdis = DRT::Problem::Instance()->GetDis("structure");

        // get GID of the first solid element (by using an homogenous material this
        // is enough)
        // ask the partner element about his Id
        const int structgid = ele->Id();

        // get the pointer to the adequate structure element based on GIDs
        DRT::Element* structele = structdis->gElement(structgid);

        // call ThermoStVenantKirchhoff material and get the temperature dependent
        // tangent ctemp
        Teuchos::RCP<MAT::Material> structmat = structele->Material();
        if (structmat->MaterialType() == INPAR::MAT::m_thermopllinelast)
          plasticmat_ = true;

      }  // disp!=0 & vel!=0
    }  // la.Size>1

    // get the time step size
    const double stepsize = params.get<double>("delta time");

    // geometrically linear TSI problem
    if ( (kintype == 0) && (la.Size()>1) )
    {
      CalculateCouplFintCondCapa(
        ele,
        time,
        mydisp,
        myvel,
        stepsize,
        &etang,
        &efint,
        NULL,
        NULL,
        NULL
        );

      if (plasticmat_)
        CalculateInternalDissipation(
          ele,
          myvel,
          stepsize,
          &etang,
          &efint
          );
    }
    // geometrically nonlinear TSI problem
    else if (kintype == 1)
    {
      CalculateNlnCouplFintCondCapa(
        ele,  // current element
        time,  // current time
        mydisp,  // current displacements
        myvel,  // current velocities
        stepsize,  // time increment
        &etang,  // element conductivity matrix
        &ecapa_,  // element capacity matrix
        &efint,  // element internal force vector
        NULL,  // heat flux at GP
        NULL,  // temperature gradients at GP
        INPAR::THR::heatflux_none,  // output option for q
        INPAR::THR::tempgrad_none  // output option for grad T
        );
    }
#endif // MonTSIwithoutSTR

    // copy capacity matrix if available
    if (ecapa.A() != NULL) ecapa.Update(ecapa_);

    // BUILD EFFECTIVE TANGENT AND RESIDUAL ACC TO TIME INTEGRATOR
    // check the time integrator
    // K_T = 1/dt . C + theta . K
    const INPAR::THR::DynamicType timint
      = DRT::INPUT::get<INPAR::THR::DynamicType>(params, "time integrator",INPAR::THR::dyna_undefined);
    switch (timint)
    {
      case INPAR::THR::dyna_statics :
      {
        // continue
        break;
      }
      case INPAR::THR::dyna_onesteptheta :
      {
        const double theta = params.get<double>("theta");
        const double stepsize = params.get<double>("delta time");
        // combined tangent and conductivity matrix to one global matrix
        etang.Update(1.0/stepsize,ecapa_,theta);
        efcap.Multiply(ecapa_,etemp_);
        break;
      }
      case INPAR::THR::dyna_genalpha :
      {
        dserror("Genalpha not yet implemented");
        break;
      }
      case INPAR::THR::dyna_undefined :
      default :
      {
        dserror("Don't know what to do...");
        break;
      }
    }  // end of switch(timint)
  }  // action == "calc_thermo_finttang"

  //============================================================================
  // Calculate/ evaluate heatflux q and temperature gradients gradtemp at
  // gauss points
  else if (action == "proc_thermo_heatflux")
  {
    // set views
    // efext, efcap not needed for this action, elemat1+2,elevec1-3 are not used anyway

    // get storage arrays of Gauss-point-wise vectors
    Teuchos::RCP<std::vector<char> > heatfluxdata
      = params.get<Teuchos::RCP<std::vector<char> > >("heatflux");
    Teuchos::RCP<std::vector<char> > tempgraddata
      = params.get<Teuchos::RCP<std::vector<char> > >("tempgrad");
    // working arrays
    LINALG::Matrix<nquad_,nsd_> eheatflux;
    LINALG::Matrix<nquad_,nsd_> etempgrad;

    // thermal problem or geometrically linear TSI problem
    if (kintype == 0)
    {
      CalculateFintCondCapa(
        ele,
        time,
        NULL,
        NULL,
        NULL,
        NULL,
        &eheatflux,
        &etempgrad
        );
    }
    // geometrically nonlinear TSI problem
    if (kintype == 1)
    {
      // specific choice of heat flux / temperature gradient
      const INPAR::THR::HeatFluxType ioheatflux
        = DRT::INPUT::get<INPAR::THR::HeatFluxType>(params,"ioheatflux",
            INPAR::THR::heatflux_none);
      const INPAR::THR::TempGradType iotempgrad
        = DRT::INPUT::get<INPAR::THR::TempGradType>(params,"iotempgrad",
            INPAR::THR::tempgrad_none);

      // if it's a TSI problem and there are current displacements/velocities
      if (la.Size()>1)
      {
        if ( (discretization.HasState(1,"displacement")) &&
           (discretization.HasState(1,"velocity")) )
        {
          // get the displacements
          std::vector<double> mydisp((la[1].lm_).size());
          Teuchos::RCP<const Epetra_Vector> disp
            = discretization.GetState(1,"displacement");
          if (disp==Teuchos::null)
            dserror("Cannot get state vectors 'displacement'");
          DRT::UTILS::ExtractMyValues(*disp,mydisp,la[1].lm_);

          // get the velocities
          std::vector<double> myvel((la[1].lm_).size());
          Teuchos::RCP<const Epetra_Vector> vel
            = discretization.GetState(1,"velocity");
          if (vel==Teuchos::null)
            dserror("Cannot get state vectors 'velocity'");
          DRT::UTILS::ExtractMyValues(*vel,myvel,la[1].lm_);

          const double stepsize = params.get<double>("delta time");

          CalculateNlnCouplFintCondCapa(
            ele,  // current element
            time,  // current time
            mydisp,  // current displacements
            myvel,  // current velocities
            stepsize,  // time increment
            NULL,  // element conductivity matrix
            NULL,  // element capacity matrix
            NULL,  // element internal force vector
            &eheatflux,  // heat flux at GP
            &etempgrad,  // temperature gradients at GP
            ioheatflux,  // output option for q
            iotempgrad  // output option for grad T
            );
        }  // disp!=0 & vel!=0
      }  // la.Size>1
    }  // end geo_nonlinear TSI

    // scale the heatflux with (-1)
    // for the calculation the heatflux enters as positive value, but
    // q = -k * (grad T)
    eheatflux.Scale(-1);

    // store in
    DRT::PackBuffer hfdata;
    ParObject::AddtoPack(hfdata, eheatflux);
    hfdata.StartPacking();
    ParObject::AddtoPack(hfdata, eheatflux);
    std::copy(hfdata().begin(),hfdata().end(),std::back_inserter(*heatfluxdata));

    DRT::PackBuffer tgdata;
    ParObject::AddtoPack(tgdata, etempgrad);
    tgdata.StartPacking();
    ParObject::AddtoPack(tgdata, etempgrad);
    std::copy(tgdata().begin(),tgdata().end(),std::back_inserter(*tempgraddata));

  }  // action == "proc_thermo_heatflux"

  //============================================================================
  // Calculate heatflux q and temperature gradients gradtemp at gauss points
  else if (action == "postproc_thermo_heatflux")
  {
    // set views
    LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_> etang(elemat1_epetra.A(),true);  // view only!
    LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_> ecapa(elemat2_epetra.A(),true);  // view only!
    LINALG::Matrix<nen_*numdofpernode_,1> efint(elevec1_epetra.A(),true);  // view only!
    // efext, efcap not needed for this action

    const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpheatfluxmap
      = params.get<Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpheatfluxmap");
    std::string heatfluxtype = params.get<std::string>("heatfluxtype","ndxyz");
    const int gid = ele->Id();
    LINALG::Matrix<nquad_,nsd_> gpheatflux(((*gpheatfluxmap)[gid])->A(),true);  // view only!

    // set views to components
    LINALG::Matrix<nen_*numdofpernode_,1> efluxx(elevec1_epetra,true);  // view only!
    LINALG::Matrix<nen_*numdofpernode_,1> efluxy(elevec2_epetra,true);  // view only!
    LINALG::Matrix<nen_*numdofpernode_,1> efluxz(elevec3_epetra,true);  // view only!

    // catch unknown heatflux types
    bool processed = false;

    // nodally
    // extrapolate heatflux q and temperature gradient gradtemp stored at GP
    if ( (heatfluxtype == "ndxyz") or (heatfluxtype == "cxyz_ndxyz") )
    {
      processed = processed or true;
      // extrapolate heatfluxes/temperature gradients at Gauss points to nodes
      // and store results in
      ExtrapolateFromGaussPointsToNodes(ele, gpheatflux, efluxx, efluxy, efluxz);
      // method only applicable if number GP == number nodes
    }

    // centered
    if ( (heatfluxtype == "cxyz") or (heatfluxtype == "cxyz_ndxyz") )
    {
      processed = processed or true;

      Teuchos::RCP<Epetra_MultiVector> eleheatflux
        = params.get<Teuchos::RCP<Epetra_MultiVector> >("eleheatflux");
      const Epetra_BlockMap& elemap = eleheatflux->Map();
      int lid = elemap.LID(gid);
      if (lid != -1)
      {
        for (int idim=0; idim<nsd_; ++idim)
        {
          //double& s = ; // resolve pointer for faster access
          double s = 0.0;
          // nquad_: number of Gauss points
          for (int jquad=0; jquad<nquad_; ++jquad)
            s += gpheatflux(jquad,idim);
          s /= nquad_;
          ( *( (*eleheatflux)(idim) ) )[lid] = s;
        }
      }
    }
    // catch unknown heatflux types
    if (not processed)
      dserror("unknown type of heatflux/temperature gradient output on element level");

  }  // action == "postproc_thermo_heatflux"

  //============================================================================
  else if (action == "integrate_shape_functions")
  {
    // calculate integral of shape functions
    const Epetra_IntSerialDenseVector dofids = params.get<Epetra_IntSerialDenseVector>("dofids");
    IntegrateShapeFunctions(ele,elevec1_epetra,dofids);
  }

  //============================================================================
  else if (action == "calc_thermo_update_istep")
  {
    ;  // do nothing
  }

  //==================================================================================
  // allowing the predictor TangTemp in .dat --> can be decisive in compressible case!
  else if (action== "calc_thermo_reset_istep")
  {
    // do nothing
  }

  //============================================================================
  // evaluation of internal thermal energy
  else if (action == "calc_thermo_energy")
  {
    // check length of elevec1
    if (elevec1_epetra.Length() < 1) dserror("The given result vector is too short.");

    double kappa = 0.0;

    // material
    Teuchos::RCP<MAT::Material> material = ele->Material();

    // get Fourier´s law (for "ordinary" thermal problem)
    if (material->MaterialType() == INPAR::MAT::m_th_fourier_iso)
    {
      const MAT::FourierIso* actmat = static_cast<const MAT::FourierIso*>(material.get());
      kappa = actmat->Capacity();
    }
    else
      dserror("Material type is not supported");

    // get node coordinates
    GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

    // declaration of internal variables
    double intenergy = 0.0;

    // ----------------------------- integration loop for one element

    // integrations points and weights
    DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
    if (intpoints.IP().nquad != nquad_)
      dserror("Trouble with number of Gauss points");

    // --------------------------------------- loop over Gauss Points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

      // call material law => cmat_,heatflux_
      Materialize(ele);

      LINALG::Matrix<1,1> temp;
      temp.MultiplyTN(1.0,funct_,etemp_,0.0);

      // internal energy
      intenergy +=kappa*fac_*temp(0,0);

    }  // -------------------------------- end loop over Gauss Points

    elevec1_epetra(0) = intenergy;

  } // evaluation of internal energy

  //============================================================================
  // add linearistion of velocity for dynamic time integration to the stiffness term
  // calculate thermal mechanical tangent matrix K_Td
  else if (action == "calc_thermo_coupltang")
  {
    LINALG::Matrix<nen_*numdofpernode_,nen_*nsd_*numdofpernode_> etangcoupl(elemat1_epetra.A(),true);

    // if it's a TSI problem and there are the current displacements/velocities
    if (la.Size()>1)
    {
      // and now get the current displacements/velocities
      if ( (discretization.HasState(1,"displacement")) &&
         (discretization.HasState(1,"velocity")) )
      {
        std::vector<double> mydisp((la[1].lm_).size());
        // get the displacements
        Teuchos::RCP<const Epetra_Vector> disp
          = discretization.GetState(1,"displacement");
        if (disp==Teuchos::null)
          dserror("Cannot get state vectors 'displacement'");
        // extract the displacements
        DRT::UTILS::ExtractMyValues(*disp,mydisp,la[1].lm_);

        std::vector<double> myvel((la[1].lm_).size());
        // get the velocities
        Teuchos::RCP<const Epetra_Vector> vel
          = discretization.GetState(1,"velocity");
        if (vel==Teuchos::null)
          dserror("Cannot get state vectors 'velocity'");
        // extract the velocities
        DRT::UTILS::ExtractMyValues(*vel,myvel,la[1].lm_);

        // if there is a strucutural vector available go on here
        // --> calculate coupling stiffness term in case of monolithic TSI

        // BUILD EFFECTIVE TANGENT ACC TO TIME INTEGRATOR
        // check the time integrator
        const INPAR::THR::DynamicType timint
         = DRT::INPUT::get<INPAR::THR::DynamicType>(params, "time integrator",INPAR::THR::dyna_undefined);

        // get step size dt
        const double stepsize = params.get<double>("delta time");

        // geometrically linear TSI problem
        if (kintype == 0)
        {
          CalculateCouplCond(
            ele,
            mydisp,
            myvel,
            &etangcoupl
            );

          // consider linearisation of velocities due to displacements
          // major switch to different time integrators
          switch (timint)
          {
            case INPAR::THR::dyna_statics :
            {
              // Lin (v_n+1) . \Delta d_n+1 = 1/dt, cf. Diss N. Karajan (2009) for quasistatic approach
              const double fac = 1.0/stepsize;
              etangcoupl.Scale(fac);
              break;
            }
            case INPAR::THR::dyna_onesteptheta :
            {
              // K_Td = theta . k_Td^e * 1/(theta * dt)
              // K_Td = k_Td^e * 1/(dt)
              const double theta = params.get<double>("theta");
              const double str_theta = params.get<double>("str_THETA");
              etangcoupl.Scale(theta/(str_theta*stepsize));
              break;
            }
            case INPAR::THR::dyna_genalpha :
            {
              dserror("Genalpha not yet implemented");

              // TODO check scaling factor for genalpha again!
              const double str_beta = params.get<double>("str_BETA");
              const double str_gamma = params.get<double>("str_GAMMA");
              // Lin (v_n+1) . \Delta d_n+1 = (gamma) / (beta . dt)
              const double fac =  str_gamma / ( str_beta * stepsize );
              etangcoupl.Scale(fac);
              break;
            }
            case INPAR::THR::dyna_undefined :
            default :
            {
              dserror("Don't know what to do...");
              break;
            }
          }  // end of switch(timint)
        }

        // geometrically nonlinear TSI problem
        if (kintype == 1)
        {
          CalculateNlnCouplCond(
            ele,
            mydisp,
            myvel,
            params,
            &etangcoupl
            );

          switch (timint)
          {
            case INPAR::THR::dyna_statics :
            {
              // continue
              break;
            }
            case INPAR::THR::dyna_onesteptheta :
            {
              // k^e_Td = k^e_Td
              //   + theta . N_T^T . (-C_temp) . 1/2 C'_lin . N_T . T . detJ . w(gp)
              //   - theta . ( B_T^T . C_mat . C^{-1}_lin . B_T . T . detJ . w(gp) )
              const double theta = params.get<double>("theta");
              // K_Td = theta . K_Td
              etangcoupl.Scale(theta);
              break;
            }
            case INPAR::THR::dyna_genalpha :
            {
              // Lin (v_n+1) . \Delta d_n+1 = (gamma) / (beta . dt)
              const double str_beta = params.get<double>("str_BETA");
              const double str_gamma = params.get<double>("str_GAMMA");
              const double fac =  str_gamma / ( str_beta * stepsize );
              etangcoupl.Scale(fac);
              dserror("Genalpha not yet implemented");
              break;
            }
            case INPAR::THR::dyna_undefined :
            default :
            {
              dserror("Don't know what to do...");
              break;
            }
          }  // end of switch(timint)

        }
      }   // disp!=0 & vel!=0
    }  // la.Size>1

  }  // action == "calc_thermo_coupltang"

  //============================================================================
  else
  {
    dserror("Unknown type of action for Temperature Implementation: %s",action.c_str());
  }

#ifdef THRASOUTPUT
  std::cout << "etemp_ end of Evaluate thermo_ele_impl\n" << etemp_ << std::endl;
#endif // THRASOUTPUT

  return 0;

}  // Evaluate()


/*----------------------------------------------------------------------*
 | evaluate the external volume load                        bborn 09/09 |
 | condition corresponding to radiation r^ over dv with                 |
 | r^ = rho . r = scalar, i.e. even for geo nln: no difference          |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TemperImpl<distype>::EvaluateNeumann(
  DRT::Element* ele,
  Teuchos::ParameterList& params,
  DRT::Discretization& discretization,
  std::vector<int>& lm,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseMatrix* elemat1_epetra
  )
{
  // check length
  if (lm.size() != nen_*numdofpernode_)
    dserror("Location vector length does not match!");
  // set views
  LINALG::Matrix<nen_*numdofpernode_,1> efext(elevec1_epetra,true);  // view only!
  // disassemble temperature
  if (discretization.HasState(0,"temperature"))
  {
    std::vector<double> mytempnp(lm.size());
    Teuchos::RCP<const Epetra_Vector> tempnp
      = discretization.GetState("temperature");
    if (tempnp == Teuchos::null)
      dserror("Cannot get state vector 'tempnp'");
    DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,lm);
    LINALG::Matrix<nen_*numdofpernode_,1> etemp(&(mytempnp[0]),true);  // view only!
    etemp_.Update(etemp);  // copy
  }
  // check for the action parameter
  const std::string action = params.get<std::string>("action","none");
  // extract time
  const double time = params.get<double>("total time");

  // perform actions
  if (action == "calc_thermo_fext")
  {
    // so far we assume deformation INdependent external loads, i.e. NO
    // difference between geometrically (non)linear TSI

    // we prescribe a scalar value on the volume, constant for (non)linear analysis
    CalculateFintCondCapa(
      ele,
      time,
      NULL,
      NULL,
      NULL,
      &efext,
      NULL,
      NULL
      );
  }
  else
  {
    dserror("Unknown type of action for Temperature Implementation: %s",action.c_str());
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | calculate system matrix and rhs r_T(T), k_TT(T) (public) g.bau 08/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateFintCondCapa(
  DRT::Element* ele,  // the element whose matrix is calculated
  const double& time,  // current time
  LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* etang,  // conductivity matrix
  LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* ecapa,  // capacity matrix
  LINALG::Matrix<nen_*numdofpernode_,1>* efint,  // internal force
  LINALG::Matrix<nen_*numdofpernode_,1>* efext,  // external force
  LINALG::Matrix<nquad_,nsd_>* eheatflux,  // heat fluxes at Gauss points
  LINALG::Matrix<nquad_,nsd_>* etempgrad  // temperature gradients at Gauss points
  )
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // ---------------------------------------------------------------------
  // call routine for calculation of radiation in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  if (efext != NULL)
  {
    Radiation(ele,time);
  }

  // ------------------------------- integration loop for one element

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_)
    dserror("Trouble with number of Gauss points");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    // get radiation in gausspoint
    if (efext != NULL)
    {
      // fext = fext + N . r. detJ * w(gp)
      // with funct_: shape functions, fac_:detJ * w(gp)
      efext->MultiplyNN(fac_,funct_,radiation_,1.0);
    }

    // gradient of current temperature value
    // grad T = d T_j / d x_i = L . N . T = B_ij T_j
    gradtemp_.MultiplyNN(derxy_,etemp_);

    // store the temperature gradient for postprocessing
    if (etempgrad != NULL)
      for (int idim=0; idim<nsd_; ++idim)
        // (8x3)                = (3x1)
        (*etempgrad)(iquad,idim) = gradtemp_(idim);

    // call material law => cmat_,heatflux_
    // negative q is used for balance equation: -q = -(-k gradtemp)= k * gradtemp
    Materialize(ele);

    // store the heat flux for postprocessing
    if (eheatflux != NULL)
      for (int idim=0; idim<nsd_; ++idim)
        (*eheatflux)(iquad,idim) = heatflux_(idim);

#ifdef THRASOUTPUT
    std::cout << "CalculateFintCondCapa heatflux_ = " << heatflux_ << std::endl;
    std::cout << "CalculateFintCondCapa gradtemp_ = " << gradtemp_ << std::endl;
    if (etempgrad != NULL) std::cout << "CalculateFintCondCapa Nln etempgrad = " << *etempgrad << std::endl;
    if (eheatflux != NULL) std::cout << "CalculateFintCondCapa Nln eheatflux = " << *eheatflux << std::endl;
#endif  // THRASOUTPUT

    // internal force vector
    if (efint != NULL)
    {
      // fint = fint + B^T . q . detJ * w(gp)
      efint->MultiplyTN(fac_,derxy_,heatflux_,1.0);
    }

    // conductivity matrix
    if (etang != NULL)
    {
      // ke = ke + ( B^T . C_mat . B ) * detJ * w(gp)  with C_mat = k * I
      LINALG::Matrix<nsd_,nen_> aop(false); // (3x8)
      // -q = C * B
      aop.MultiplyNN(cmat_,derxy_); //(nsd_xnsd_)(nsd_xnen_)
      etang->MultiplyTN(fac_,derxy_,aop,1.0); //(nen_xnen_)=(nen_xnsd_)(nsd_xnen_)
    }

    // capacity matrix (equates the mass matrix in the structural field)
    if (ecapa != NULL)
    {
      // ce = ce + ( N^T .  (rho * C_V) . N ) * detJ * w(gp)
      // (8x8)      (8x1)               (1x8)
      // caution: funct_ implemented as (8,1)--> use transposed in code for
      // theoretic part
      ecapa->MultiplyNT((fac_*capacoeff_),funct_,funct_,1.0);
    }

   }  // --------------------------------- end loop over Gauss Points

#ifdef THRASOUTPUT
  if (efint != NULL) std::cout << "element No. = "<< ele->Id() << " efint f_Td CalculateFintCondCapa"<< *efint << std::endl;
  if (etang != NULL) std::cout << "element No. = " << ele->Id() <<  " etang nach CalculateFintCondCapa"<< *etang << std::endl;
#endif  // THRASOUTPUT
} // CalculateFintCondCapa


/*----------------------------------------------------------------------*
 | calculate coupled fraction for the system matrix          dano 05/10 |
 | and rhs: r_T(d), k_TT(d) (public)                                    |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateCouplFintCondCapa(
  DRT::Element* ele,  // the element whose matrix is calculated
  const double& time,  // current time
  std::vector<double>& disp,  // current displacements
  std::vector<double>& vel,  // current velocities
  const double& stepsize,
  LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* etang,  // conductivity matrix
  LINALG::Matrix<nen_*numdofpernode_,1>* efint,  // internal force
  LINALG::Matrix<nen_*numdofpernode_,1>* efext,  // external force
  LINALG::Matrix<nquad_,nsd_>* eheatflux,  // heat fluxes at Gauss points
  LINALG::Matrix<nquad_,nsd_>* etempgrad  // temperature gradients at Gauss points
  )
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // now get current element displacements
  LINALG::Matrix<nen_*nsd_,1> edisp;
  LINALG::Matrix<nen_*nsd_,1> evel;
  for (int i=0; i<nen_*nsd_; i++)
  {
    edisp(i,0) = disp[i+0];
    evel(i,0) = vel[i+0];
  }

#ifdef THRASOUTPUT
  std::cout << "CalculateCoupl evel\n" << evel << std::endl;
  std::cout << "edisp\n" << edisp << std::endl;
#endif // THRASOUTPUT

  // thermal material tangent
  LINALG::Matrix<6,1> ctemp(true);
  // get constant initial temperature from the material
  double thetainit = 0.0;
  Teuchos::RCP<MAT::Material> structmat = Teuchos::null;
  GetStrMaterial(ele, &ctemp, &thetainit, structmat);
  // insert the negative value of the coupling term (c.f. energy balance)
  ctemp.Scale(-1.0);

  // ------------------------------- integration loop for one element

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_)
    dserror("Trouble with number of Gauss points");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // compute inverse Jacobian matrix and derivatives at GP w.r.t material
    // coordinates
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    // calculate the linear B-operator
    LINALG::Matrix<6,nsd_*nen_*numdofpernode_> boplin;
    CalculateBoplin(&boplin,&derxy_);

    // now build the strain rates / velocities
    LINALG::Matrix<6,1> strainvel(true);
    // e' = B . d' = B . v = 0.5 * (Grad u' + Grad^T u')
    strainvel.Multiply(boplin,evel);  // (6x24)(24x1)=(6x1)

    // in case of a thermo-elasto-plastic solid material, strainvel != elastic strains
    // e' = (e^e)' + (e^p)'
    // split strainvel (=total strain) into elastic and plastic terms
    // --> thermomechanical coupling term requires elastic strain rates and
    // --> dissipation term requires the plastic strain rates
    // call the structural material
    LINALG::Matrix<6,1> plasticstrainlinrate(true);

    if (plasticmat_)
    {
      MAT::ThermoPlasticLinElast* thrpllinelast
        = static_cast <MAT::ThermoPlasticLinElast*>(structmat.get());

      // strainvel includes the elastic strain rates of the structural velocity vel
      // --> if CalcVelocity() is used in tsi
      //     --> calculate the elastic velocity with the given veln_ again
      thrpllinelast->StrainRateSplit(iquad,stepsize,strainvel);
      // overwrite the total strain rate by the elastic strain rate
      strainvel.Update(1.0, (thrpllinelast->ElasticStrainRate(iquad)), 0.0);
      plasticstrainlinrate.Update(1.0, (thrpllinelast->PlasticStrainRate(iquad)), 0.0);
    }

#ifdef CALCSTABILOFREACTTERM
    // scalar product Ctemp : (B . (d^e)')
    // in case of elastic step Ctemp : (B . (d^e)') ==  Ctemp : (B . d')
    double cbv = 0.0;
    for (int i=0; i<6; ++i)
      cbv += ctemp(i,0)*strainvel(i,0);

    // ------------------------------------ start reactive term check
    // check reactive term for stability
    // check critical parameter of reactive term
    // K = sigma / ( kappa * h^2 ) > 1 --> problems occur
    // kappa: kinematic diffusitivity
    // sigma = m I : (B . (d^e)') = Ctemp : (B . (d^e)')
    double sigma = cbv;
    std::cout << "sigma = " << sigma << std::endl;
    std::cout << "h = " << h << std::endl;
    std::cout << "h^2 = " << h*h << std::endl;
    std::cout << "kappa = " << kappa << std::endl;
    std::cout << "strainvel = " << strainvel << std::endl;
    // critical parameter for reactive dominated problem
    double K_thr = sigma / ( kappa * (h*h) );
    std::cout << "K_thr abs = " << abs(K_thr) << std::endl;
    if (abs(K_thr) > 1.0)
      std::cout << "stability problems can occur: abs(K_thr) = " << abs(K_thr) << std::endl;
    // -------------------------------------- end reactive term check
#endif  // CALCSTABILOFREACTTERM

    // N^T . (- Ctemp) : ( B .  (d^e)' )
    LINALG::Matrix<nen_,6> nctemp(true); // (8x1)(1x6)
    nctemp.MultiplyNT(funct_,ctemp);
    LINALG::Matrix<nen_,1> ncBv;
    ncBv.Multiply(nctemp,strainvel);

    // integrate internal force vector (coupling fraction towards displacements)
    if (efint != NULL)
    {
      // build the product of the shapefunctions and element temperatures
      LINALG::Matrix<1,1> nt(true);
#ifdef COUPLEINITTEMPERATURE
      // for TSI validation/verification: change nt to Theta_0 here!!!! 14.01.11
      if (ele->Id()==0)
        std::cout << "ele Id= " << ele->Id() " coupling term in thermo field with T_0" << std::endl;
      nt(0,0) = thetainit;
#else
      // default: use scalar-valued current temperature T = N . T
      nt.MultiplyTN(funct_,etemp_);
#endif

      // fintdisp = fintdisp - N^T . Ctemp : (B .  (d^e)') . N . T
      efint->Multiply(fac_,ncBv,nt,1.0);

#ifdef TSIMONOLITHASOUTPUT
      if (ele->Id()==0)
      {
        std::cout << "efint nach CalculateCoupl"<< *efint << std::endl;
        std::cout << "CouplFint\n" << std::endl;
        std::cout << "ele Id= " << ele->Id() << std::endl;
        std::cout << "boplin\n" << boplin << std::endl;
        std::cout << "etemp_ Ende CalculateCouplFintCondCapa\n" << etemp_  << std::endl;
        std::cout << "ctemp_\n" << ctemp << std::endl;
        std::cout << "ncBv\n" << ncBv << std::endl;
      }
#endif  // TSIMONOLITHASOUTPUT
    }  // if (efint!=NULL)

    // update conductivity matrix (with displacement dependent term)
    if (etang != NULL)
    {
      // k^e = k^e - ( N^T . (-m * I) . (B_d . (d^e)') . N ) * detJ * w(gp)
      // --> negative term enters the tangent (cf. L923) ctemp.Scale(-1.0);
      etang->MultiplyNT(fac_,ncBv,funct_,1.0);
    }  // if (etang!=NULL)

#ifdef THRASOUTPUT
    if (efint != NULL) std::cout << "element No. = "<< ele->Id() << " efint f_Td CalculateCouplFintCondCapa"<< *efint << std::endl;
    if (etang != NULL) std::cout << "element No. = " << ele->Id() <<  " etang nach CalculateCouplFintCondCapa"<< *etang << std::endl;
#endif  // THRASOUTPUT

  }  // ---------------------------------- end loop over Gauss Points

} // CalculateCouplFintCondCapa()


/*----------------------------------------------------------------------*
 | calculate thermal-mechanical system matrix k_Td needed    dano 03/11 |
 | in monolithic TSI (private)                                          |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateCouplCond(
  DRT::Element* ele,  // the element whose matrix is calculated
  std::vector<double>& disp,  // current displacements
  std::vector<double>& vel,  // current velocities
  LINALG::Matrix<nen_*numdofpernode_,nsd_*nen_*numdofpernode_>* etangcoupl  // conductivity matrix
  )
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // now get current element displacements
  LINALG::Matrix<nen_*nsd_,1> edisp;
  LINALG::Matrix<nen_*nsd_,1> evel;
  for (int i=0; i<nen_*nsd_; i++)
  {
    edisp(i,0) = disp[i+0];
    evel(i,0) = vel[i+0];
  }

  // in case of thermo-elasto-plastic material: elasto-plastic tangent modulus
  LINALG::Matrix<6,6> cmat(true);
  // thermal material tangent
  LINALG::Matrix<6,1> ctemp(true);
  // get constant initial temperature from the material
  double thetainit = 0.0;
  // TODO 2012-11-14 in case of different material, pass structmat here, too
  GetStrMaterial(ele, &ctemp, &thetainit, Teuchos::null);
  // insert the negative value of the coupling term (c.f. energy balance)
  ctemp.Scale(-1.0);

  // ------------------------------- integration loop for one element

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_)
    dserror("Trouble with number of Gauss points");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    // GEOMETRIC LINEAR problem the deformation gradient is equal to identity

    // calculate the linear B-operator
    LINALG::Matrix<6,nsd_*nen_*numdofpernode_> boplin;
    CalculateBoplin(&boplin,&derxy_);

    // non-symmetric stiffness matrix
    // current element temperatures
    // N_temp . temp (funct_ defined as <nen,1>
    LINALG::Matrix<1,1> nt(false);
    nt.MultiplyTN(funct_,etemp_);  // (1x8)(8x1)= (1x1)
    // N_temp^T . N_temp . temp
    LINALG::Matrix<nen_,1> nnt(false);
    nnt.Multiply(funct_,nt); // (8x1)(1x1) = (8x1)

    // N_T^T . N_T . T . Ctemp
    LINALG::Matrix<nen_,6> nntc(false); // (8x1)(1x6)
    nntc.MultiplyNT(nnt,ctemp);  // (8x6)

#ifdef TSIMONOLITHASOUTPUT
      if (ele->Id()==0)
      {
        std::cout << "Coupl Cond\n" << std::endl;
        std::cout << "ele Id= " << ele->Id() << std::endl;
        std::cout << "boplin \n" << boplin << std::endl;
        std::cout << "etemp_ Ende CalculateCouplCond\n" << etemp_  << std::endl;
        std::cout << "ctemp_\n" << ctemp << std::endl;
        std::cout << "nntc\n" << nntc << std::endl;
      }
#endif  // TSIMONOLITHASOUTPUT

    // 03.09.11 TODO in case of thermo-elasto-plastic material:
    // get elasto-plastic tangent modulus out of material C_ep
    // K_Td = K_Td - C_ep . strain_p'

    // coupling stiffness matrix
    if (etangcoupl != NULL)
    {
      // k_Td^e = k_Td^e + ( N_T^T . N_T . T . (-C_temp) . B_L ) * detJ * w(gp)
      // with C_temp = m * I
      // (8x24) = (8x6) . (6x24)
      etangcoupl->MultiplyNN(fac_,nntc,boplin,1.0);
      // TODO: 19.08.11. term for K_Td of Internal Dissipation has to be added here as well!!!
    }

   }  //---------------------------------- end loop over Gauss Points

} // CalculateCouplCond()


/*----------------------------------------------------------------------*
 | calculate coupled fraction for the system matrix          dano 11/12 |
 | and rhs: r_T(T,d), k_TT(T,d) (public)                                |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateNlnCouplFintCondCapa(
  DRT::Element* ele,  // the element whose matrix is calculated
  const double& time,  // current time
  std::vector<double>& disp,  // current displacements
  std::vector<double>& vel,  // current velocities
  const double& stepsize,  // time increment
  LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* etang,  // conductivity matrix
  LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* ecapa,  // capacity matrix
  LINALG::Matrix<nen_*numdofpernode_,1>* efint,  // internal force
  LINALG::Matrix<nquad_,nsd_>* eheatflux,  // heat fluxes at Gauss points
  LINALG::Matrix<nquad_,nsd_>* etempgrad,  // temperature gradients at Gauss points
  const INPAR::THR::HeatFluxType ioheatflux,  // heat flux output option
  const INPAR::THR::TempGradType iotempgrad  // tempgrad output option
  )
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // now get current element displacements
  LINALG::Matrix<nen_*nsd_,1> edisp;
  LINALG::Matrix<nen_*nsd_,1> evel;
  for (int i=0; i<nen_*nsd_; i++)
  {
    // (24x1) = (nen_*nsd_x1)
    edisp(i,0) = disp[i+0];
    evel(i,0) = vel[i+0];
  }

  // update element geometry
  LINALG::Matrix<nen_,nsd_> xrefe;  // material coord. of element
  LINALG::Matrix<nen_,nsd_> xcurr;  // current  coord. of element
  LINALG::Matrix<nen_,nsd_> xcurrrate;  // current  coord. of element

  DRT::Node** nodes = ele->Nodes();
  for (int i=0; i<nen_; ++i)
  {
    const double* x = nodes[i]->X();
    // (8x3) = (nen_xnsd_)
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*nsd_+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*nsd_+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*nsd_+2];

    xcurrrate(i,0) = vel[i*nsd_+0];
    xcurrrate(i,1) = vel[i*nsd_+1];
    xcurrrate(i,2) = vel[i*nsd_+2];
  }

#ifdef THRASOUTPUT
  std::cout << "CalculateNlnCoupl evel\n" << evel << std::endl;
  std::cout << "edisp\n" << edisp << std::endl;
  std::cout << "edisp" << edisp << std::endl;
  std::cout << "evel" << evel << std::endl;
  std::cout << "xrefe" << xrefe << std::endl;
  std::cout << "xcurr" << xcurr << std::endl;
  std::cout << "xcurrrate" << xcurrrate << std::endl;
  std::cout << "derxy_" << derxy_ << std::endl;
#endif // THRASOUTPUT

  // thermal material tangent
  LINALG::Matrix<6,1> ctemp(true);
  // get constant initial temperature from the material
  double thetainit = 0.0;
  // TODO 2012-11-14 in case of different material, pass structmat here, too
  GetStrMaterial(ele, &ctemp, &thetainit, Teuchos::null);
  // insert the negative value of the coupling term (c.f. energy balance)
  ctemp.Scale(-1.0);

  // build the deformation gradient w.r.t material configuration
  LINALG::Matrix<nsd_,nsd_> defgrd(false);
  // build the rate of the deformation gradient w.r.t material configuration
  LINALG::Matrix<nsd_,nsd_> defgrdrate(false);
  // inverse of deformation gradient
  LINALG::Matrix<nsd_,nsd_> invdefgrd(false);

  // ------------------------------- integration loop for one element

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_)
    dserror("Trouble with number of Gauss points");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // compute inverse Jacobian matrix and derivatives at GP w.r.t material
    // coordinates
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    // --------------------------------------------- thermal gradient
    // gradient of current temperature value
    // Grad T = d T_j / d x_i = L . N . T = B_ij T_j
    gradtemp_.MultiplyNN(derxy_,etemp_);

    // ---------------------------------------- coupling to mechanics
    // (material) deformation gradient F
    // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,derxy_);
    // rate of (material) deformation gradient F'
    // F' = d xcurr' / d xrefe = (xcurr')^T * N_XYZ^T
    defgrdrate.MultiplyTT(xcurrrate,derxy_);
    // inverse of deformation gradient
    invdefgrd.Invert(defgrd);
    // build the linear B-operator
    LINALG::Matrix<6,nen_*nsd_*numdofpernode_> boplin;
    CalculateBoplin(&boplin,&derxy_);
    // build the nonlinear B-operator
    LINALG::Matrix<6,nen_*nsd_*numdofpernode_> bop;
    CalculateBop(&bop,&defgrd,&derxy_);

    // ------- derivatives of right Cauchy-Green deformation tensor C
    // build the rate of C: C'= F^T . F' + (F')^T . F
    // OR: C' = F^T . F' if applied to symmetric tensor
    // save C' as rate vector Crate
    // C' = { C11', C22', C33', C12', C23', C31 }
    LINALG::Matrix<6,1> Cratevct(false);
    // build the inverse C: C^{-1} = F^{-1} . F^{-T}
    LINALG::Matrix<nsd_,nsd_> invC(false);
    // invCvct: C^{-1} in Voight-/vector notation
    // C^{-1} = { C11^{-1}, C22^{-1}, C33^{-1}, C12^{-1}, C23^{-1}, C31^{-1} }
    LINALG::Matrix<6,1> invCvct(false);
    CalculateCauchyGreens(Cratevct,invCvct,invC,&defgrd,&defgrdrate,&invdefgrd);

    // -------------------------- calculate strain rates / velocities
    LINALG::Matrix<6,1> strainvel(true);
    // e' = B_L. d' = B_L . v = 0.5 * (Grad u' + Grad^T u')
    strainvel.Multiply(boplin,evel);  // (6x24)(24x1)=(6x1)

    // ------------------------------------ call thermal material law
    // call material law => cmat_,heatflux_
    // negative q is used for balance equation:
    // heatflux_ = k_0 * Grad T
    Materialize(ele);
    // heatflux_ := qintermediate = k_0 * Grad T

    // initial heatflux Q = C^{-1} . qintermediate = k_0 . C^{-1} . B_T . T
    // the current heatflux q = detF . F^{-1} . q
    // store heatflux
    // (3x1)  (3x3) . (3x1)
    LINALG::Matrix<nsd_,1> initialheatflux(true);
    initialheatflux.Multiply(invC,heatflux_);
    // put the initial, material heatflux onto heatflux_
    heatflux_.Update(0.0, initialheatflux, 1.0);
    // from here on heatflux_ == -Q

    // ---------------------------------------------- post processing
    // store the temperature gradient for postprocessing
    // return gp tempgrad (only in case of tempgrad output)
    // RK: Grad T
    // AK: grad T --> grad T = Grad T . F^{-1}
    switch (iotempgrad)
    {
    case INPAR::THR::tempgrad_initial:
    {
     if (etempgrad == NULL) dserror("tempgrad data not available");
     // etempgrad = Grad T
     for (int idim=0; idim<nsd_; ++idim)
       (*etempgrad)(iquad,idim) = gradtemp_(idim);
    }
    break;
    case INPAR::THR::tempgrad_current:
    {
      if (etempgrad == NULL) dserror("tempgrad data not available");
      // etempgrad = grad T = Grad T . F^{-1} =  F^{-T} . Grad T
      // (8x3)        (3x1)   (3x1)    (3x3)     (3x3)    (3x1)
      // spatial temperature gradient
      LINALG::Matrix<nsd_,1> currentgradT(false);
      currentgradT.MultiplyTN(invdefgrd,gradtemp_);
      for (int idim=0; idim<nsd_; ++idim)
        (*etempgrad)(iquad,idim) = currentgradT(idim);
    }
    break;
    case INPAR::THR::tempgrad_none:
    {
      // no postprocessing of temperature gradients
    }
    break;
    default:
     dserror("requested tempgrad type not available");
    break;
    }  // iotempgrad

    // return gp heatfluxes (only in case of heatflux/tempgrad output)
    // RK: Q = -k_0 . invC . Grad T
    // AK: q = -k . grad T --> q =  1/(detF) . F .Q
    // with k_0 = J . k = detF . k
    switch (ioheatflux)
    {
    case INPAR::THR::heatflux_initial:
    {
      // eheatflux := Q = -k_0 . invC . Grad T
      if (eheatflux == NULL) dserror("heat flux data not available");
      for (int idim=0; idim<nsd_; ++idim)
        (*eheatflux)(iquad,idim) = heatflux_(idim);
    }
    break;
    case INPAR::THR::heatflux_current:
    {
      if (eheatflux == NULL) dserror("heat flux data not available");
      // eheatflux := q = 1/(detF) . F .Q
      // (8x3)     (3x1)            (3x3)  (3x1)
      const double detF = defgrd.Determinant();
      LINALG::Matrix<nsd_,1> spatialq;
      spatialq.MultiplyNN((1.0/detF), defgrd, heatflux_);
      for (int idim=0; idim<nsd_; ++idim)
        (*eheatflux)(iquad,idim) = spatialq(idim);
    }
    break;
    case INPAR::THR::heatflux_none:
    {
      // no postprocessing of heat fluxes, continue!
    }
    break;
    default:
      dserror("requested heat flux type not available");
    break;
    }  // ioheatflux

#ifdef THRASOUTPUT
    if (etempgrad != NULL)
      std::cout << "CalculateNlnCouplFintCondCapa Nln etempgrad = " << *etempgrad << std::endl;
    if (eheatflux != NULL)
      std::cout << "CalculateNlnCouplFintCondCapa Nln eheatflux = " << *eheatflux << std::endl;
#endif  // THRASOUTPUT

    // ----------------------------------------- terms for r_T / k_TT
    // scalar product: ctempcdot = -(m * I) : 1/2 C' = -C_temp : 1/2 C'
    double CtempCdot = 0.0;
    for (int i=0; i<6; ++i)
      CtempCdot += ctemp(i,0)*(1/2.0)*Cratevct(i,0);

    // build the product of the shapefunctions and element temperatures
    LINALG::Matrix<1,1> nt(true);
#ifdef COUPLEINITTEMPERATURE
    // for TSI validation/verification: change nt to Theta_0 here!!!!
    if (ele->Id()==0)
      std::cout << "ele Id= " << ele->Id() ": coupling term in thermo field with T_0" << std::endl;
    nt(0,0) = thetainit;
#else
    // default: use scalar-valued current temperature T = N . T
    nt.MultiplyTN(funct_,etemp_);
#endif

    // -------------------------- integrate internal force vector r_T
    // add the displacement-dependent terms to fint
    // fint = fint + fint_{Td}
    if (efint != NULL)
    {
      // fint = fint + B^T . Q . detJ * w(gp)
      //      = fint + B_T^T . (k_0) . C^{-1} . B_T . T . detJ * w(gp)
      // (8x1)        (8x3) (3x1)
      efint->MultiplyTN(fac_,derxy_,heatflux_,1.0);

      // fint_{Td} = - N^T . Ctemp : 1/2 . C' . N . T
      //              (1x8)  (6x1)       (6x1)(8x1)(8x1)
      //              (1x8)        (1x1)        (1x1)
      // fint = fint + fint_{Td}
      // with fint_{Td} = fint_{Td} - 1/2 . N^T . Ctemp : C' . N . T
      //                  + B^T . k_0 . F^{-1} . F^{-T} . B . T
      efint->Multiply((fac_*CtempCdot),funct_,nt,1.0);
    }  // if (efint!=NULL)

    // --------------------------- integrate conductivity matrix k_TT
    // update conductivity matrix k_TT (with displacement dependent term)
    if (etang != NULL)
    {
      // k^e_TT = k^e_TT + ( B_T^T . C_mat . C^{-1} . B_T ) * detJ * w(gp)
      // 3D:      (8x8)      (8x3)   (3x3)   (3x3)   (3x8)
      // with C_mat = k * I
      LINALG::Matrix<nsd_,nen_> aop(false); // (3x8)
      // -q = C_mat . C^{-1} . B
      aop.MultiplyNN(invC,derxy_); //(nsd_xnsd_)(nsd_x8)
      LINALG::Matrix<nsd_,nen_> aop1(false); // (3x8)
      aop1.MultiplyNN(cmat_,aop); //(nsd_Xnsd_)(nsd_x8)
      etang->MultiplyTN(fac_,derxy_,aop1,1.0); //(8x8)=(8x3)(3x8)
      // k^e_TT = k^e_TT + ( N^T . (-C_temp) : 1/2  C' . N ) * detJ * w(gp)
      // --> negative term enters the tangent (cf. L923) ctemp.Scale(-1.0);
      // with CtempCdot = (-C_temp) : 1/2  C'
      etang->MultiplyNT((fac_*CtempCdot),funct_,funct_,1.0);
    }  // if (etang!=NULL)

    // --------------------------------------- capacity matrix m_capa
    // capacity matrix is idependent of deformation
    // m_capa corresponds to the mass matrix of the structural field
    if (ecapa != NULL)
    {
      // m_capa = m_capa + ( N_T^T .  (rho * C_V) . N_T ) * detJ * w(gp)
      //           (8x8)     (8x1)                 (1x8)
      // caution: funct_ implemented as (8,1)--> use transposed in code for
      // theoretic part
      ecapa->MultiplyNT((fac_*capacoeff_),funct_,funct_,1.0);
    }
  }  // ---------------------------------- end loop over Gauss Points

#ifdef TSIMONOLITHASOUTPUT
    if (ele->Id()==0)
    {
      std::cout << "CouplNlnFintCondCapa\n" << std::endl;
      std::cout << "ele Id= " << ele->Id() << std::endl;
      std::cout << "ctemp_\n" << ctemp << std::endl;
      std::cout << "Cratevct\n" << Cratevct << std::endl;
      std::cout << "nccdot\n" << nccdot << std::endl;
      std::cout << "nctemp\n" << nctemp << std::endl;
      std::cout << "defgrd\n" << defgrd << std::endl;
    }
    std::cout << "heatflux_ CalculateNlnCouplFintCondCapa"<< heatflux_ << std::endl;
    std::cout << "etemp_ CalculateNlnCouplFintCondCapa\n" << etemp_  << std::endl;

    if (efint != NULL) std::cout << "element No. = " << ele->Id() << " efint f_Td CalculateNlnCouplFintCondCapa"<< *efint << std::endl;
    if (etang != NULL) std::cout << "element No. = " << ele->Id() << " etang nach CalculateNlnCouplFintCondCapa"<< *etang << std::endl;
#endif  // TSIMONOLITHASOUTPUT

}  // CalculateNlnCouplFintCondCapa()


/*----------------------------------------------------------------------*
 | calculate thermal-mechanical system matrix k_Td(d)        dano 11/12 |
 | needed in monolithic TSI (private)                                   |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateNlnCouplCond(
  DRT::Element* ele,  // the element whose matrix is calculated
  std::vector<double>& disp,  // current displacements
  std::vector<double>& vel,  // current velocities
  Teuchos::ParameterList& params,  // parameter list
  LINALG::Matrix<nen_*numdofpernode_,nsd_*nen_*numdofpernode_>* etangcoupl
  )
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // now get current element displacements
  LINALG::Matrix<nen_*nsd_,1> edisp;
  LINALG::Matrix<nen_*nsd_,1> evel;
  for (int i=0; i<nen_*nsd_; i++)
  {
    // (24x1) = (nen_*nsd_x1)
    edisp(i,0) = disp[i+0];
    evel(i,0) = vel[i+0];
  }

  // update element geometry
  LINALG::Matrix<nen_,nsd_> xrefe;  // material coord. of element
  LINALG::Matrix<nen_,nsd_> xcurr;  // current  coord. of element
  LINALG::Matrix<nen_,nsd_> xcurrrate;  // current  coord. of element

  DRT::Node** nodes = ele->Nodes();
  for (int i=0; i<nen_; ++i)
  {
    const double* x = nodes[i]->X();
    // (8x3) = (nen_xnsd_)
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*nsd_+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*nsd_+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*nsd_+2];

    xcurrrate(i,0) = vel[i*nsd_+0];
    xcurrrate(i,1) = vel[i*nsd_+1];
    xcurrrate(i,2) = vel[i*nsd_+2];
  }

#ifdef THRASOUTPUT
  std::cout << "CalculateNlnCoupl evel\n" << evel << std::endl;
  std::cout << "edisp\n" << edisp << std::endl;
  std::cout << "edisp" << edisp << std::endl;
  std::cout << "evel" << evel << std::endl;
  std::cout << "xrefe" << xrefe << std::endl;
  std::cout << "xcurr" << xcurr << std::endl;
  std::cout << "xcurrrate" << xcurrrate << std::endl;
  std::cout << "derxy_" << derxy_ << std::endl;
#endif // THRASOUTPUT

  // thermal material tangent
  LINALG::Matrix<6,1> ctemp(true);
  // get constant initial temperature from the material
  double thetainit = 0.0;
  // TODO 2012-11-14 in case of different material, pass structmat here, too
  GetStrMaterial(ele, &ctemp, &thetainit, Teuchos::null);
  // insert the negative value of the coupling term (c.f. energy balance)
  ctemp.Scale(-1.0);

  // build the deformation gradient w.r.t material configuration
  LINALG::Matrix<nsd_,nsd_> defgrd(false);
  // build the rate of the deformation gradient w.r.t material configuration
  LINALG::Matrix<nsd_,nsd_> defgrdrate(false);
  // inverse of deformation gradient
  LINALG::Matrix<nsd_,nsd_> invdefgrd(true);

  // ------------------------------- integration loop for one element

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_)
    dserror("Trouble with number of Gauss points");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // compute inverse Jacobian matrix and derivatives at GP w.r.t material
    // coordinates
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    // ------------------------------------------------ thermal terms

    // gradient of current temperature value
    // grad T = d T_j / d x_i = L . N . T = B_ij T_j
    gradtemp_.MultiplyNN(derxy_,etemp_);

    // call material law => cmat_,heatflux_
    // negative q is used for balance equation: -q = -(-k gradtemp)= k * gradtemp
    Materialize(ele);

    // put thermal material tangent in vector notation
    LINALG::Matrix<6,1> cmat_vct(true);
    if (nsd_==3)
    {
      cmat_vct(0) = cmat_(0,0);
      cmat_vct(1) = cmat_(1,1);
      cmat_vct(2) = cmat_(2,2);
    }
    else if (nsd_==2)
    {
      cmat_vct(0) = cmat_(0,0);
      cmat_vct(1) = cmat_(1,1);
    }
    else if (nsd_==1)
    {
      cmat_vct(0) = cmat_(0,0);
    }

    // current element temperatures
    // N_T . T (funct_ defined as <nen,1>
    LINALG::Matrix<1,1> nt(false);
    nt.MultiplyTN(funct_,etemp_);  // (1x8)(8x1)= (1x1)
    // N_T^T . N_T . T
    LINALG::Matrix<nen_,1> nnt(false);
    nnt.Multiply(funct_,nt); // (8x1)(1x1) = (8x1)
    // N_T^T . N_T . T . Ctemp
    LINALG::Matrix<nen_,6> nntc(false); // (8x1)(1x6)
    nntc.MultiplyNT(nnt,ctemp);  // (8x6)
    // B_T^T . B_T . T
    LINALG::Matrix<nen_,1> bgradT(false);
    bgradT.MultiplyTN(derxy_,gradtemp_); // (8x1)(1x1) = (8x1)
    // B_T^T . B_T . T . Cmat_
    LINALG::Matrix<nen_,6> bgradTcmat(false); // (8x1)(1x6)
    bgradTcmat.MultiplyNT(bgradT,cmat_vct);  // (8x6)

    // ---------------------------------------- coupling to mechanics
    // (material) deformation gradient F
    // F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,derxy_);
    // rate of (material) deformation gradient F'
    // F' = d xcurr' / d xrefe = (xcurr')^T * N_XYZ^T
    defgrdrate.MultiplyTT(xcurrrate,derxy_);
    // inverse of deformation gradient
    invdefgrd.Invert(defgrd);
    // build the linear B-operator
    LINALG::Matrix<6,nsd_*nen_*numdofpernode_> boplin;
    CalculateBoplin(&boplin,&derxy_);
    // build the nonlinear B-operator
    LINALG::Matrix<6,nen_*nsd_*numdofpernode_> bop;
    CalculateBop(&bop,&defgrd,&derxy_);

    // ------- derivatives of right Cauchy-Green deformation tensor C

    // build the rate of C: C'= F^T . F' + (F')^T . F
    // save C' as rate vector Crate
    // C' = { C11', C22', C33', C12', C23', C31 }
    LINALG::Matrix<6,1> Cratevct(false);
    // build the inverse C: C^{-1} = F^{-1} . F^{-T}
    LINALG::Matrix<nsd_,nsd_> invC(false);
    // invCvct: C^{-1} in Voight-/vector notation
    // C^{-1} = { C11^{-1}, C22^{-1}, C33^{-1}, C12^{-1}, C23^{-1}, C31^{-1} }
    LINALG::Matrix<6,1> invCvct(false);
    // calculation is done in CalculateCauchyGreens, return C', C^{-1}
    CalculateCauchyGreens(
      Cratevct,
      invCvct,
      invC,
      &defgrd,
      &defgrdrate,
      &invdefgrd
      );

    // -------------------------------- calculate linearisation of C'
    // C_T : 1/2 C'_lin --> symmetric part of C'_lin is sufficient
    // C'_lin = dCrate/dd = 1/2 . [ 1/(theta . dt) . (B^T + B) + (F')^T . B_L + B_L^T . F' ]
    //        = 1/(theta . dt) [ B^T + B ] + [ (F')^T . B_L + ( (F')^T . B_L )^T ]
    // C_T : 1/2 C'_lin = C_T : [ 1/(theta . Dt) B + B' ]
    // --> use only the symmetric part of C'_lin

    // with B' = (F')^T . B_L: calculate rate of B
    LINALG::Matrix<6,nen_*nsd_> boprate;  // (6x24)
    CalculateBop(&boprate,&defgrdrate,&derxy_);

    // ---------------------------- calculate linearisation of C^{-1}
    // calculate linearisation of C^{-1} according to so3_poro_evaluate
    // C^{-1}_lin = dCinv_dd = - F^{-1} . ( B_L . F^{-1} + F^{-T} . B_L^T ) . F^{-T}
    //                       = - F^{-1} . ( B_L . F^{-1} + (B_L . F^{-1})^T ) . F^{-T}
    LINALG::Matrix<6,nen_*nsd_> dCinv_dd (true);
    for (int n=0; n<nen_; ++n)
      for (int k=0; k<nsd_; ++k)
      {
        const int gid = n*nsd_+k;
        for (int i=0; i<nsd_; ++i)
        {
          dCinv_dd(0,gid) += -2*invC(0,i)*derxy_(i,n)*invdefgrd(0,k);
          dCinv_dd(1,gid) += -2*invC(1,i)*derxy_(i,n)*invdefgrd(1,k);
          dCinv_dd(2,gid) += -2*invC(2,i)*derxy_(i,n)*invdefgrd(2,k);
          /* ~~~ */
          dCinv_dd(3,gid) += -invC(0,i)*derxy_(i,n)*invdefgrd(1,k)-invdefgrd(0,k)*derxy_(i,n)*invC(1,i);
          dCinv_dd(4,gid) += -invC(1,i)*derxy_(i,n)*invdefgrd(2,k)-invdefgrd(1,k)*derxy_(i,n)*invC(2,i);
          dCinv_dd(5,gid) += -invC(2,i)*derxy_(i,n)*invdefgrd(0,k)-invdefgrd(2,k)*derxy_(i,n)*invC(0,i);
        }
      }

    // BUILD EFFECTIVE TANGENT ACC TO TIME INTEGRATOR
    // check the time integrator
    const INPAR::THR::DynamicType timint
     = DRT::INPUT::get<INPAR::THR::DynamicType>(params, "time integrator",INPAR::THR::dyna_undefined);
    double theta = 1.0;
    switch (timint)
    {
      case INPAR::THR::dyna_onesteptheta :
      {
        theta = params.get<double>("theta");
        break;
      }
      case INPAR::THR::dyna_undefined :
      default :
      {
        // put theta = 1.0
        break;
      }
    }  // end of switch(timint)
    // get step size dt
    const double stepsize = params.get<double>("delta time");

    // ----------------- coupling matrix k_Td only for monolithic TSI
    if (etangcoupl != NULL)
    {
      // k^e_Td = k^e_Td + theta . N_T^T . (-C_temp) . 1/2 C'_lin . N_T . T . detJ . w(gp)
      //          - theta . ( B_T^T . C_mat . C^{-1}_lin . B_T . T . detJ . w(gp) )
      // with
      // B_T: thermal gradient matrix
      // B_L: linear B-operator, gradient matrix == B_T
      // B: nonlinear B-operator, i.e. B = F^T . B_L
      // C'_lin = 1/(theta . Dt) ( B^T + B ) + F'T . B_L + B_L^T . F'
      // --> 1/2 C'_lin = sym C'_lin = 1/(theta . Dt) . B + B'
      // with boprate := B' = F'T . B_L
      // C^{-1}_lin = - F^{-1} . (B_L . F^{-1} + B_L^{T} . F^{-T}) . F^{-T}
      //
      // C_mat = k_0 * I

      // k^e_Td = theta . N_T^T . N_T . T . (-C_temp) . 1/2 C'_lin . detJ . w(gp)
      // (8x24)            (8x3) (3x8)(8x1)   (6x1)       (6x24)
      // (8x24)               (8x8)   (8x1)   (1x6)       (6x24)
      // (8x24)                  (8x1)        (1x6)       (6x24)
      // (8x24)                         (8x6)             (6x24)
      etangcoupl->Multiply(fac_,nntc,boprate,1.0);
      etangcoupl->Multiply( (fac_/(theta*stepsize)), nntc, bop, 1.0);
      // k^e_Td = k^e_Td - theta . ( B_T^T . C_mat . C^{-1}_lin . B_T . T . detJ . w(gp) )
      //        = k^e_Td + theta . ( B_T^T . C_mat . B_T . T . C^{-1}_lin . detJ . w(gp) )
      // (8x24)                      (8x3)   (3x3)  (3x8)(8x1)  (6x24)
      //                                 (8x3)        (3x1)
      //                                       (8x1) (1x24)
      //        = k^e_Td + theta . ( B_T^T . B_T . T . C_mat . C^{-1}_lin . detJ . w(gp) )
      // (8x24)                      (8x3)  (3x8)(8x1) (1x6) (6x24)
      etangcoupl->MultiplyNN(fac_,bgradTcmat,dCinv_dd,1.0);
    }  // (etangcoupl != NULL)

  }  // ---------------------------------- end loop over Gauss Points

}


/*----------------------------------------------------------------------*
 | calculate internal dissipation term, used in case of      dano 08/11 |
 | plastic material (private)                                           |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateInternalDissipation(
  DRT::Element* ele,  // the element whose matrix is calculated
  std::vector<double>& vel,  // current velocities
  const double& stepsize,
  LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* etang,  // conductivity matrix
  LINALG::Matrix<nen_*numdofpernode_,1>* efint  // internal force
  )
{

  // 2012-11-14 TODO adapt to geometrically nonlinear analysis

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // get current element displacements
  LINALG::Matrix<nen_*nsd_,1> evel;
  for (int i=0; i<nen_*nsd_; i++)
  {
    evel(i,0) = vel[i+0];
  }
#ifdef THRASOUTPUT
  std::cout << "CalculateInternalDissipation evel\n" << evel << std::endl;
#endif // THRASOUTPUT

  // thermal material tangent
  LINALG::Matrix<6,1> ctemp(true);
  // get constant initial temperature from the material
  double thetainit = 0.0;
  Teuchos::RCP<MAT::Material> structmat = Teuchos::null;
  GetStrMaterial(ele, &ctemp, &thetainit, structmat);
  // insert the negative value of the coupling term (c.f. energy balance)
  // TODO 2012-11-14 so far no scaling was used, correct??
  // ctemp.Scale(-1.0);

  // ------------------------------- integration loop for one element

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_)
    dserror("Trouble with number of Gauss points");

  // ----------------------------------------- loop over Gauss Points
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // compute inverse Jacobian matrix and derivatives at GP w.r.t material
    // coordinates
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    // GEOMETRIC LINEAR problem the deformation gradient is equal to identity

    // build the linear B-operator
    LINALG::Matrix<6,nsd_*nen_*numdofpernode_> boplin;
    CalculateBoplin(&boplin,&derxy_);

    // now build the strains and the strain velocities/rates
    LINALG::Matrix<6,1> strainvel(true);
    // e' = B . d' = B . v = 0.5 * (Grad u' + Grad^T u')
    strainvel.Multiply(boplin,evel);  // (6x24)(24x1)=(6x1)

    // call material ThermoPlasticLinElast and get the plastic strain rates
    LINALG::Matrix<6,1> plasticstrainlinrate(true);
    LINALG::Matrix<6,1> stresstemp(true);

    // plastic power
    double plasticpower = 0.0;
    // plastic temperature dependent stress power --> for K_tt
    double plastictemppower = 0.0;

    if (plasticmat_)
    {
      MAT::ThermoPlasticLinElast* thrpllinelast
        = static_cast <MAT::ThermoPlasticLinElast*>(structmat.get());
      // strainvel includes the elastic strain rates of the structural velocity vel
      // --> if CalcVelocity() is used in tsi
      //     --> calculate the elastic velocity with the given veln_ again
      //     --> build the plastic strain rate by dividing by stepsize
      thrpllinelast->StrainRateSplit(iquad,stepsize,strainvel);
      plasticstrainlinrate.Update(1.0, (thrpllinelast->PlasticStrainRate(iquad)), 0.0);
      // in the material the increment is calculated, NOT the temporal rate
      // --> divide by dt
      // plasticpower_ = (sigma^{mech}+sigma^{theta} - beta) : strain^p_{n+1}'

      // compute element temperature N_T . T
      LINALG::Matrix<1,1> Ntemp(true);
      Ntemp.MultiplyTN(funct_,etemp_);
      // in the material: 1.) Delta T = subtract ( N . T - T_0 )
      //                  2.) stresstemp = C . Delta T
      thrpllinelast->Evaluate(*(&Ntemp),*(&ctemp),*(&stresstemp));

      // plasticpower_ = (sigma^{mech} - beta) : strain^p_{n+1}'
      plasticpower = thrpllinelast->PlasticPower(iquad);
      plasticpower *= (1.0/stepsize);
      // plasticpower_ += sigma^{theta}_n+1 : strain^p_{n+1}'
      double tempstressplstrain = stresstemp(0)*plasticstrainlinrate(0) +
                                  stresstemp(1)*plasticstrainlinrate(1) +
                                  stresstemp(2)*plasticstrainlinrate(2) +
                                  stresstemp(3)*plasticstrainlinrate(3) +
                                  stresstemp(4)*plasticstrainlinrate(4) +
                                  stresstemp(5)*plasticstrainlinrate(5);
      plasticpower += tempstressplstrain;

      // C_mat^theta : epsilon_p'
      // in the material the increment is calculated, NOT the temporal rate
      // --> divide by dt
      plastictemppower = thrpllinelast->PlasticTempPower(iquad);
      // to be a power the term has to be multiplied by 1/dt
      plastictemppower *= 1.0/stepsize;
    }
    // ---------------------------------------------------- END OF COUPLING

    // update/integrate internal force vector (coupling fraction towards displacements)
    if (efint != NULL)
    {
      // update of the internal "force" vector
      efint->Update((fac_*(-plasticpower)),funct_,1.0);
    }
    
    // update/integrate coupling conductivity matrix (displacement dependent)
    if (etang != NULL)
    {
      // k = k - ( N^T . C_mat^theta . epsilon_p' ) * detJ * w(gp)
      // with C_mat^theta = m * I
      etang->MultiplyNT((fac_*(-plastictemppower)),funct_,funct_,1.0); // (8x8) = (8x1)(8x1)^T
    }

#ifdef TSIMONOLITHASOUTPUT
      if (ele->Id()==0)
      {
        std::cout << "CouplFint\n" << std::endl;
        std::cout << "boplin\n" << boplin << std::endl;
        std::cout << "etemp_ Ende InternalDiss\n" << etemp_  << std::endl;
      }
#endif  // TSIMONOLITHASOUTPUT

  }  // ---------------------------------- end loop over Gauss Points

} // CalculateInternalDissipation()


/*----------------------------------------------------------------------*
 | lump capacity matrix (private)                            dano 01/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateLumpMatrix(LINALG::Matrix<nen_*numdofpernode_,nen_*numdofpernode_>* ecapa)
{
  // lump capacity matrix
  if (ecapa != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c=0; c<(*ecapa).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r=0; r<(*ecapa).M(); ++r)  // parse rows
      {
        d += (*ecapa)(r,c);  // accumulate row entries
        (*ecapa)(r,c) = 0.0;
      }
      (*ecapa)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
}  // CalculateLumpMatrix()


/*----------------------------------------------------------------------*
 | get the radiation  (private)                              dano 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::Radiation(
  const DRT::Element* ele,
  const double time
  )
{
  std::vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  switch(nsd_)
  {
  case 3:
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);
  break;
  case 2:
    DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);
  break;
  case 1:
    DRT::UTILS::FindElementConditions(ele, "LineNeumann", myneumcond);
  break;
  default:
    dserror("Illegal number of space dimensions: %d",nsd_);
  }

  if (myneumcond.size()>1)
    dserror("more than one VolumeNeumann cond on one node");

  if (myneumcond.size()==1)
  {
    // find out whether we will use a time curve
    const std::vector<int>* curve = myneumcond[0]->Get<std::vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac(0.0);

    if (curvenum >= 0) // yes, we have a timecurve
    {
      // time factor for the intermediate step
      if(time >= 0.0)
      {
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
      }
      else
      {
        // A negative time value indicates an error.
        dserror("Negative time value in body force calculation: time = %f",time);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // get values and switches from the condition
    const std::vector<int>*    onoff = myneumcond[0]->Get<std::vector<int> >   ("onoff");
    const std::vector<double>* val   = myneumcond[0]->Get<std::vector<double> >("val"  );

    // set this condition to the radiation array
    for (int idof=0; idof<numdofpernode_; idof++) {
      radiation_(idof) = (*onoff)[idof]*(*val)[idof]*curvefac;
    }

  }
  else
  {
    // we have no dead load
    radiation_.Clear();
  }
  return;

} //TemperImpl::Radiation


/*----------------------------------------------------------------------*
 | get the material                                          dano 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::Materialize(
  const DRT::Element* ele
  )
{
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  // get Fourier´s law (for "ordinary" thermal problem)
  if (material->MaterialType() == INPAR::MAT::m_th_fourier_iso)
  {
    const MAT::FourierIso* actmat = static_cast<const MAT::FourierIso*>(material.get());
    actmat->Evaluate(gradtemp_,cmat_,heatflux_);
    capacoeff_ = actmat->Capacity();
  }
  else
  {
    dserror("Material type is not supported");
  }

  return;
} //TemperImpl::Materialize


/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at int. point     gjb 08/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
  const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  // integration points
  const int iquad,  // id of current Gauss point
  const int eleid  // the element id
  )
{
  // coordinates of the current integration point (xsi_)
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
    {xsi_(idim) = gpcoord[idim];}

  // shape functions (funct_) and their first derivatives (deriv_)
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
   */

  // derivatives at gp w.r.t material coordinates (N_XYZ in solid)
  xjm_.MultiplyNT(deriv_,xyze_);
  // xij_ = J^(T-1)
  // det = J^(T-1) *
  // J = (N_rst * X)^T (6.24 NiliFEM)
  const double det = xij_.Invert(xjm_);

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det);

  // set integration factor: fac = Gauss weight * det(J)
  fac_ = intpoints.IP().qwgt[iquad]*det;

  // compute global derivatives
  derxy_.Multiply(xij_,deriv_);

  // say goodbye
  return;

} // EvalShapeFuncAndDerivsAtIntPoint


/*----------------------------------------------------------------------*
 | integrate shape functions over domain (private)            gjb 07/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::IntegrateShapeFunctions(
  const DRT::Element* ele,
  Epetra_SerialDenseVector& elevec1,
  const Epetra_IntSerialDenseVector& dofids
  )
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,gpid,ele->Id());

    // compute integral of shape functions (only for dofid)
    for (int k=0; k<numdofpernode_; k++)
    {
      if (dofids[k] >= 0)
      {
        for (int node=0; node<nen_; node++)
        {
          elevec1[node*numdofpernode_+k] += funct_(node) * fac_;
        }
      }
    }
  } //loop over integration points

  return;
} //TemperImpl<distype>::IntegrateShapeFunction

/*----------------------------------------------------------------------*
 | extrapolateFromGaussPointsToNodes                         dano 11/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::ExtrapolateFromGaussPointsToNodes(
  DRT::Element* ele,  // the element whose matrix is calculated
  const LINALG::Matrix<nquad_,nsd_>& gpheatflux,
  LINALG::Matrix<nen_*numdofpernode_,1>& efluxx,
  LINALG::Matrix<nen_*numdofpernode_,1>& efluxy,
  LINALG::Matrix<nen_*numdofpernode_,1>& efluxz
  )
{
  // this quick'n'dirty hack functions only for hex8
  // (number of gauss points equals number of nodes)
  if ( not ( (distype == DRT::Element::hex8)
             or (distype == DRT::Element::quad4)
             or (distype == DRT::Element::line2) ) )
    dserror("Sorry, not implemented for element shape");

  // another check
  if (nen_*numdofpernode_ != nquad_)
    dserror("Works only if number of gauss points and nodes match");

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_)
    dserror("Trouble with number of Gauss points");

  // build matrix of shape functions at Gauss points
  LINALG::Matrix<nquad_,nquad_> shpfctatgps;
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // coordinates of the current integration point
    const double* gpcoord = (intpoints.IP().qxg)[iquad];
    for (int idim=0; idim<nsd_; idim++)
      xsi_(idim) = gpcoord[idim];

    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_,funct_);

    for (int inode=0; inode<nen_; ++inode)
      shpfctatgps(iquad,inode) = funct_(inode);
  }

  // extrapolation
  LINALG::Matrix<nquad_,nsd_> ndheatflux;  //  objective nodal heatflux
  LINALG::Matrix<nquad_,nsd_> gpheatflux2(gpheatflux);  // copy the heatflux at the Gauss point
  {
    LINALG::FixedSizeSerialDenseSolver<nquad_,nquad_,nsd_> solver;  // must be quadratic
    solver.SetMatrix(shpfctatgps);
    solver.SetVectors(ndheatflux,gpheatflux2);
    solver.Solve();
  }

  // copy into component vectors
  for (int idof=0; idof<nen_*numdofpernode_; ++idof)
  {
    efluxx(idof) = ndheatflux(idof,0);
    if (nsd_>1) efluxy(idof) = ndheatflux(idof,1);
    if (nsd_>2) efluxz(idof) = ndheatflux(idof,2);
  }

  // bye
  return;
} // ExtrapolateFromGaussPointsToNodes


/*----------------------------------------------------------------------*
 | calculation of characteristic element length              dano 02/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::TemperImpl<distype>::CalculateCharEleLength()
{
  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = fac_;

  // as shown in CalcCharEleLength() in ScaTraImpl
  // c) cubic/square root of element volume/area or element length (3-/2-/1-D)
  // cast dimension to a double varible -> pow()

  // get number of dimensions
  const double dim = (double) nsd_;

  // get characteristic element length as cubic root of element volume
  // (2D: square root of element area, 1D: element length)
  // h = vol^(1/dim)
  double h = std::pow(vol,(1.0/dim));

  return h;

}  // CalculateCharEleLength()


/*----------------------------------------------------------------------*
 | calculate the linear B-operator                           dano 11/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateBoplin(
  LINALG::Matrix<6,nsd_*nen_*numdofpernode_>* boplin,
  LINALG::Matrix<nsd_,nen_>* N_XYZ
  )
{
  // in thermo element derxy_ == N_XYZ in structural element (i.e. So3_Thermo)
  // lump mass matrix
  if (boplin != NULL)
  {
    // linear B-operator B_L = N_XYZ
    // disperse global derivatives to bop-lines
    // bop is arranged as usual (refer to script FE or elsewhere):
    // [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
    // [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
    // [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
    // [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
    // [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
    // [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ]
    for (int i=0; i<nen_; ++i)
    {
      (*boplin)(0,nsd_*numdofpernode_*i+0) = (*N_XYZ)(0,i);
      (*boplin)(0,nsd_*numdofpernode_*i+1) = 0.0;
      (*boplin)(0,nsd_*numdofpernode_*i+2) = 0.0;
      (*boplin)(1,nsd_*numdofpernode_*i+0) = 0.0;
      (*boplin)(1,nsd_*numdofpernode_*i+1) = (*N_XYZ)(1,i);
      (*boplin)(1,nsd_*numdofpernode_*i+2) = 0.0;
      (*boplin)(2,nsd_*numdofpernode_*i+0) = 0.0;
      (*boplin)(2,nsd_*numdofpernode_*i+1) = 0.0;
      (*boplin)(2,nsd_*numdofpernode_*i+2) = (*N_XYZ)(2,i);
      /* ~~~ */
      (*boplin)(3,nsd_*numdofpernode_*i+0) = (*N_XYZ)(1,i);
      (*boplin)(3,nsd_*numdofpernode_*i+1) = (*N_XYZ)(0,i);
      (*boplin)(3,nsd_*numdofpernode_*i+2) = 0.0;
      (*boplin)(4,nsd_*numdofpernode_*i+0) = 0.0;
      (*boplin)(4,nsd_*numdofpernode_*i+1) = (*N_XYZ)(2,i);
      (*boplin)(4,nsd_*numdofpernode_*i+2) = (*N_XYZ)(1,i);
      (*boplin)(5,nsd_*numdofpernode_*i+0) = (*N_XYZ)(2,i);
      (*boplin)(5,nsd_*numdofpernode_*i+1) = 0.0;
      (*boplin)(5,nsd_*numdofpernode_*i+2) = (*N_XYZ)(0,i);
    }
  }
}  // CalculateBoplin()


/*----------------------------------------------------------------------*
 | calculate the nonlinear B-operator                        dano 11/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateBop(
  LINALG::Matrix<6,nsd_*nen_*numdofpernode_>* bop,
  LINALG::Matrix<nsd_,nsd_>* defgrd,
  LINALG::Matrix<nsd_,nen_>* N_XYZ
  )
{
  // lump mass matrix
  if (bop != NULL)
  {
    /* non-linear B-operator (may so be called, meaning of B-operator is not so
    ** sharp in the non-linear realm) *
    ** B = F . B_L *
    ** with linear B-operator B_L =  N_XYZ (6x24) = (3x8)
    **
    **   B    =   F  . N_XYZ
    ** (6x24)   (3x3) (3x8)
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    for (int i=0; i<nen_; ++i)
    {
      (*bop)(0,nsd_*numdofpernode_*i+0) = (*defgrd)(0,0)*(*N_XYZ)(0,i);
      (*bop)(0,nsd_*numdofpernode_*i+1) = (*defgrd)(1,0)*(*N_XYZ)(0,i);
      (*bop)(0,nsd_*numdofpernode_*i+2) = (*defgrd)(2,0)*(*N_XYZ)(0,i);
      (*bop)(1,nsd_*numdofpernode_*i+0) = (*defgrd)(0,1)*(*N_XYZ)(1,i);
      (*bop)(1,nsd_*numdofpernode_*i+1) = (*defgrd)(1,1)*(*N_XYZ)(1,i);
      (*bop)(1,nsd_*numdofpernode_*i+2) = (*defgrd)(2,1)*(*N_XYZ)(1,i);
      (*bop)(2,nsd_*numdofpernode_*i+0) = (*defgrd)(0,2)*(*N_XYZ)(2,i);
      (*bop)(2,nsd_*numdofpernode_*i+1) = (*defgrd)(1,2)*(*N_XYZ)(2,i);
      (*bop)(2,nsd_*numdofpernode_*i+2) = (*defgrd)(2,2)*(*N_XYZ)(2,i);
      /* ~~~ */
      (*bop)(3,nsd_*numdofpernode_*i+0) = (*defgrd)(0,0)*(*N_XYZ)(1,i) + (*defgrd)(0,1)*(*N_XYZ)(0,i);
      (*bop)(3,nsd_*numdofpernode_*i+1) = (*defgrd)(1,0)*(*N_XYZ)(1,i) + (*defgrd)(1,1)*(*N_XYZ)(0,i);
      (*bop)(3,nsd_*numdofpernode_*i+2) = (*defgrd)(2,0)*(*N_XYZ)(1,i) + (*defgrd)(2,1)*(*N_XYZ)(0,i);
      (*bop)(4,nsd_*numdofpernode_*i+0) = (*defgrd)(0,1)*(*N_XYZ)(2,i) + (*defgrd)(0,2)*(*N_XYZ)(1,i);
      (*bop)(4,nsd_*numdofpernode_*i+1) = (*defgrd)(1,1)*(*N_XYZ)(2,i) + (*defgrd)(1,2)*(*N_XYZ)(1,i);
      (*bop)(4,nsd_*numdofpernode_*i+2) = (*defgrd)(2,1)*(*N_XYZ)(2,i) + (*defgrd)(2,2)*(*N_XYZ)(1,i);
      (*bop)(5,nsd_*numdofpernode_*i+0) = (*defgrd)(0,2)*(*N_XYZ)(0,i) + (*defgrd)(0,0)*(*N_XYZ)(2,i);
      (*bop)(5,nsd_*numdofpernode_*i+1) = (*defgrd)(1,2)*(*N_XYZ)(0,i) + (*defgrd)(1,0)*(*N_XYZ)(2,i);
      (*bop)(5,nsd_*numdofpernode_*i+2) = (*defgrd)(2,2)*(*N_XYZ)(0,i) + (*defgrd)(2,0)*(*N_XYZ)(2,i);
    }
  }
}  // CalculateBop()

/*----------------------------------------------------------------------*
 | calculate the nonlinear B-operator                        dano 11/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateCauchyGreens(
  LINALG::Matrix<6,1>& Cratevct,  // (io) C' in vector notation
  LINALG::Matrix<6,1>& invCvct,  // (io) C^{-1} in vector notation
  LINALG::Matrix<nsd_,nsd_>& invC,  // (io) C^{-1} in tensor notation
  LINALG::Matrix<nsd_,nsd_>* defgrd,  // (i) deformation gradient
  LINALG::Matrix<nsd_,nsd_>* defgrdrate,  // (i) rate of deformation gradient
  LINALG::Matrix<nsd_,nsd_>* invdefgrd  // (i) inverse of deformation gradient
  )
{
  // calculate the rate of the right Cauchy-Green deformation gradient C'
  // rate of right Cauchy-Green tensor C' = F^T . F' + (F')^T . F
  // C'= F^T . F' + (F')^T . F
  // OR: C' = F^T . F' when applied to symmetric tensor
  LINALG::Matrix<nsd_,nsd_> Crate(false);
  Crate.MultiplyTN((*defgrd), (*defgrdrate));
  Crate.MultiplyTN(1.0, (*defgrdrate), (*defgrd), 1.0);

  // copy to matrix notation
  // rate vector Crate C'
  // C' = { C11', C22', C33', C12', C23', C31' }
  Cratevct(0) = Crate(0,0);
  Cratevct(1) = Crate(1,1);
  Cratevct(2) = Crate(2,2);
  Cratevct(3) = Crate(0,1);
  Cratevct(4) = Crate(1,2);
  Cratevct(5) = Crate(2,0);

  // build the inverse of the right Cauchy-Green deformation gradient C^{-1}
  // C^{-1} = F^{-1} . F^{-T}
  invC.MultiplyNT((*invdefgrd), (*invdefgrd));
  // invCvct: C^{-1} in Voight-/vector notation
  // C^{-1} = { C11^{-1}, C22^{-1}, C33^{-1}, C12^{-1}, C23^{-1}, C31^{-1} }
  invCvct(0) = invC(0,0);
  invCvct(1) = invC(1,1);
  invCvct(2) = invC(2,2);
  invCvct(3) = invC(0,1);
  invCvct(4) = invC(1,2);
  invCvct(5) = invC(2,0);

}  // CalculateCauchyGreens()


/*----------------------------------------------------------------------*
 | get the corresponding structural material                 dano 11/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::GetStrMaterial(
  DRT::Element* ele,  // the element whose matrix is calculated
  LINALG::Matrix<6,1>* ctemp,  // temperature-dependent material tangent,
  double* thetainit,
  Teuchos::RCP<MAT::Material> structmat
  )
{
  if (DRT::Problem::Instance()->DoesExistDis("structure"))
  {
    // access the structure discretization, needed later for calling the solid
    // material and getting its tangent
    Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
    structdis = DRT::Problem::Instance()->GetDis("structure");

    // get GID of the first solid element (by using an homogenous material this
    // is enough)
    // ask the partner element about his Id
    const int structgid = ele->Id();

    // get the pointer to the adequate structure element based on GIDs
    DRT::Element* structele = structdis->gElement(structgid);

#ifdef CALCSTABILOFREACTTERM
  // check critical parameter of reactive term
  // initialise kinematic diffusivity for checking stability of reactive term
  // kappa = k/(rho C_V) = Conductivity()/Capacitity()
  double kappa = 0.0;
  // calculate element length h = (vol)^(dim)
  double h = CalculateCharEleLength();
//  std::cout << "h = " << h << std::endl;
//  double h2 = h^2;
#endif  // CALCSTABILOFREACTTERM

    // call ThermoStVenantKirchhoff material and get the temperature dependent
    // tangent ctemp
    Teuchos::RCP<MAT::Material> structmat = structele->Material();
    if (structmat->MaterialType() == INPAR::MAT::m_thermostvenant)
    {
      MAT::ThermoStVenantKirchhoff* thrstvk
        = static_cast <MAT::ThermoStVenantKirchhoff*>(structmat.get());
      thrstvk->SetupCthermo(*ctemp);
  #ifdef CALCSTABILOFREACTTERM
      // kappa = k / (rho C_V)
      kappa = thrstvk->Conductivity();
      kappa /= thrstvk->Capacity();
  #endif  // CALCSTABILOFREACTTERM

      //  for TSI validation/verification (2nd Danilovskaya problem): use COUPLEINITTEMPERATURE
  #ifdef COUPLEINITTEMPERATURE
      thetainit = thrstvk->InitTemp();
  #endif // COUPLEINITTEMPERATURE
    }  // m_thermostvenant

    else if (structmat->MaterialType() == INPAR::MAT::m_thermopllinelast)
    {
      MAT::ThermoPlasticLinElast* thrpllinelast
        = static_cast <MAT::ThermoPlasticLinElast*>(structmat.get());
      thrpllinelast->SetupCthermo(*ctemp);

      //  for TSI validation/verification (2nd Danilovskaya problem): use COUPLEINITTEMPERATURE
  #ifdef COUPLEINITTEMPERATURE
      thetainit = thrpllinelast->InitTemp();
  #endif // COUPLEINITTEMPERATURE
    }  // m_thermopllinelast

  }  // if structural discretisation exists
}  // GetSTRMaterial()


#ifdef CALCSTABILOFREACTTERM
/*----------------------------------------------------------------------*
 | get the corresponding structural material                 dano 11/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateReactiveTerm(
  LINALG::Matrix<6,1>* ctemp,  // temperature-dependent material tangent
  LINALG::Matrix<6,1>* strainvel  // strain rate
  )
{
  // scalar product Ctemp : (B . (d^e)')
  // in case of elastic step Ctemp : (B . (d^e)') ==  Ctemp : (B . d')
  double cbv = 0.0;
  for (int i=0; i<6; ++i)
    cbv += ctemp(i,0)*strainvel(i,0);

  // ------------------------------------ start reactive term check
  // check reactive term for stability
  // check critical parameter of reactive term
  // K = sigma / ( kappa * h^2 ) > 1 --> problems occur
  // kappa: kinematic diffusitivity
  // sigma = m I : (B . (d^e)') = Ctemp : (B . (d^e)')
  double sigma = cbv;
  std::cout << "sigma = " << sigma << std::endl;
  std::cout << "h = " << h << std::endl;
  std::cout << "h^2 = " << h*h << std::endl;
  std::cout << "kappa = " << kappa << std::endl;
  std::cout << "strainvel = " << strainvel << std::endl;
  // critical parameter for reactive dominated problem
  double K_thr = sigma / ( kappa * (h*h) );
  std::cout << "K_thr abs = " << abs(K_thr) << std::endl;
  if (abs(K_thr) > 1.0)
    std::cout << "stability problems can occur: abs(K_thr) = " << abs(K_thr) << std::endl;
  // -------------------------------------- end reactive term check
}
#endif  // CALCSTABILOFREACTTERM


/*----------------------------------------------------------------------*/
#endif // D_THERMO
