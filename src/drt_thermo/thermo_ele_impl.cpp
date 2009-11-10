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
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef D_THERMO

#include "../drt_inpar/inpar_thermo.H"

#include "thermo_ele_impl.H"
#include "../drt_mat/fourieriso.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_condition_utils.H"

//#define VISUALIZE_ELEMENT_DATA
#include "thermo_element.H" // only for visualization of element data

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TemperImplInterface* DRT::ELEMENTS::TemperImplInterface::Impl(
    DRT::Element* ele
    )
 {
  //! we assume here, that numdofpernode is equal for every node within
  //! the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    static TemperImpl<DRT::Element::hex8>* ch8;
    if (ch8==NULL)
      ch8 = new TemperImpl<DRT::Element::hex8>(numdofpernode);
    return ch8;
  }
/*  case DRT::Element::hex20:
  {
    static TemperImpl<DRT::Element::hex20>* ch20;
    if (ch20==NULL)
      ch20 = new TemperImpl<DRT::Element::hex20>(numdofpernode);
    return ch20;
  }
  case DRT::Element::hex27:
  {
    static TemperImpl<DRT::Element::hex27>* ch27;
    if (ch27==NULL)
      ch27 = new TemperImpl<DRT::Element::hex27>(numdofpernode);
    return ch27;
  }*/
  case DRT::Element::tet4:
  {
    static TemperImpl<DRT::Element::tet4>* ct4;
    if (ct4==NULL)
      ct4 = new TemperImpl<DRT::Element::tet4>(numdofpernode);
    return ct4;
  }
 /* case DRT::Element::tet10:
  {
    static TemperImpl<DRT::Element::tet10>* ct10;
    if (ct10==NULL)
      ct10 = new TemperImpl<DRT::Element::tet10>(numdofpernode);
    return ct10;
  } */
  case DRT::Element::wedge6:
  {
    static TemperImpl<DRT::Element::wedge6>* cw6;
    if (cw6==NULL)
      cw6 = new TemperImpl<DRT::Element::wedge6>(numdofpernode);
    return cw6;
  }
/*  case DRT::Element::wedge15:
  {
    static TemperImpl<DRT::Element::wedge15>* cw15;
    if (cw15==NULL)
      cw15 = new TemperImpl<DRT::Element::wedge15>(numdofpernode);
    return cw15;
  } */
  case DRT::Element::pyramid5:
  {
    static TemperImpl<DRT::Element::pyramid5>* cp5;
    if (cp5==NULL)
      cp5 = new TemperImpl<DRT::Element::pyramid5>(numdofpernode);
    return cp5;
  }
  case DRT::Element::quad4:
  {
    static TemperImpl<DRT::Element::quad4>* cp4;
    if (cp4==NULL)
      cp4 = new TemperImpl<DRT::Element::quad4>(numdofpernode);
    return cp4;
  }
/*  case DRT::Element::quad8:
  {
    static TemperImpl<DRT::Element::quad8>* cp8;
    if (cp8==NULL)
      cp8 = new TemperImpl<DRT::Element::quad8>(numdofpernode);
    return cp8;
  }
  case DRT::Element::quad9:
  {
    static TemperImpl<DRT::Element::quad9>* cp9;
    if (cp9==NULL)
      cp9 = new TemperImpl<DRT::Element::quad9>(numdofpernode);
    return cp9;
  }*/
  case DRT::Element::tri3:
  {
    static TemperImpl<DRT::Element::tri3>* cp3;
    if (cp3==NULL)
      cp3 = new TemperImpl<DRT::Element::tri3>(numdofpernode);
    return cp3;
  }
/*  case DRT::Element::tri6:
  {
    static TemperImpl<DRT::Element::tri6>* cp6;
    if (cp6==NULL)
      cp6 = new TemperImpl<DRT::Element::tri6>(numdofpernode);
    return cp6;
  }*/
  case DRT::Element::line2:
  {
    static TemperImpl<DRT::Element::line2>* cl2;
    if (cl2==NULL)
      cl2 = new TemperImpl<DRT::Element::line2>(numdofpernode);
    return cl2;
  }/*
  case DRT::Element::line3:
  {
    static TemperImpl<DRT::Element::line3>* cl3;
    if (cl3==NULL)
      cl3 = new TemperImpl<DRT::Element::line3>(numdofpernode);
    return cl3;
  }*/
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
  }
  return NULL;
 }

/*----------------------------------------------------------------------*
 |  Initialization of the data with respect to the declaration          |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TemperImpl<distype>::TemperImpl(int numdofpernode)
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
    cmat_(false)
 {
  return;
 }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TemperImpl<distype>::Evaluate(
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
  // check length
  if (lm.size() != iel*numdofpernode_)
    dserror("Location vector length does not match!");
  // set views
  LINALG::Matrix<iel*numdofpernode_,iel*numdofpernode_> etang(elemat1_epetra,true);  // view only!
  LINALG::Matrix<iel*numdofpernode_,iel*numdofpernode_> ecapa(elemat2_epetra,true);  // view only!
  LINALG::Matrix<iel*numdofpernode_,1> efint(elevec1_epetra,true);  // view only!
  //LINALG::Matrix<iel*numdofpernode_,1> efext(elevec2_epetra,true);  // view only!
  LINALG::Matrix<iel*numdofpernode_,1> efcap(elevec3_epetra,true);  // view only!
  // disassemble temperature
  if (discretization.HasState("temperature")) {
    std::vector<double> mytempnp(lm.size());
    Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState("temperature");
    if (tempnp==Teuchos::null)
      dserror("Cannot get state vector 'tempnp'");
    DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,lm);
    LINALG::Matrix<iel*numdofpernode_,1> etemp(&(mytempnp[0]),true);  // view only!
    etemp_.Update(etemp);  // copy
  }
  // initialize capacity matrix
  ecapa_.Clear();
  // check for the action parameter
  const std::string action = params.get<std::string>("action","none");
  // extract time
  const double time = params.get<double>("total time");
  // calculate tangent K and internal force F_int = K * Theta
  if (action=="calc_thermo_fintcond")
  {
    CalculateFintCondCapa(ele,
                          time,
                          &etang,
                          NULL,
                          &efint,
                          NULL,
                          NULL,NULL);
  }
  // calculate only the internal force F_int
  else if (action=="calc_thermo_fint")
  {
    CalculateFintCondCapa(ele,
                          time,
                          NULL,
                          NULL,
                          &efint,
                          NULL,
                          NULL,NULL);
  }
  // calculate only the internal force F_int
  else if (action=="calc_thermo_fintcapa")
  {
    CalculateFintCondCapa(ele,
                          time,
                          NULL,
                          &ecapa,  // delivers capacity matrix
                          &efint,
                          NULL,
                          NULL,NULL);
    // copy capacity matrix if available
    //if (ecapa.A() != NULL) ecapa.Update(ecapa_);
  }
  // calculate tangent matrix K and consistent capacity matrix C
  else if (action=="calc_thermo_finttang")
  {
    CalculateFintCondCapa(ele,
                          time,
                          &etang,
                          &ecapa_,
                          &efint,
                          NULL,
                          NULL,NULL);
    // copy capacity matrix if available
    if (ecapa.A() != NULL) ecapa.Update(ecapa_);
    // BUILD EFFECTIVE TANGENT AND RESIDUAL ACC TO TIME INTEGRATOR
    // check the time integrator
    const INPAR::THR::DynamicType timint
    = params.get<INPAR::THR::DynamicType>("time integrator",INPAR::THR::dyna_undefined);
    switch (timint) {
    case INPAR::THR::dyna_statics :
    {
      // continue
      break;
    }
    case INPAR::THR::dyna_onesteptheta :
    {
      const double theta = params.get<double>("theta");
      const double stepsize = params.get<double>("delta time");
      etang.Update(1.0/stepsize,ecapa_,theta);  // combined tangent and conductivity matrix to one global matrix
      efcap.Multiply(ecapa_,etemp_);
      break;
    }
    case INPAR::THR::dyna_genalpha :
    {
//      //const double theta = params.get<double>("theta");
//      const double stepsize = params.get<double>("delta time");
//      etang.Update(theta/stepsize,ecapa,1.0);
      break;
    }
    case INPAR::THR::dyna_undefined :
    default :
    {
      dserror("Don't know what to do...");
      break;
    }
    }
  }
  // Calculate heatflux q and temperature gradients gradtemp at gauss points
  else if (action=="proc_thermo_heatflux")
  {
    // get storage arrays of Gauss-point-wise vectors
    Teuchos::RCP<std::vector<char> > heatfluxdata = params.get<Teuchos::RCP<std::vector<char> > >("heatflux");
    Teuchos::RCP<std::vector<char> > tempgraddata = params.get<Teuchos::RCP<std::vector<char> > >("tempgrad");
    // working arrays
    LINALG::Matrix<nquad_,nsd_> eheatflux;
    LINALG::Matrix<nquad_,nsd_> etempgrad;
    // specific choice of heat flux / temperature gradient
    //const INPAR::THR::HeatFluxType ioheatflux = params.get<INPAR::THR::HeatFluxType>("ioheatflux");
    //const INPAR::THR::TempGradType iotempgrad = params.get<INPAR::THR::TempGradType>("iotempgrad");
    //
    CalculateFintCondCapa(ele,
                          time,
                          NULL,
                          NULL,
                          NULL,
                          NULL,
                          &eheatflux,&etempgrad);
     // store in
     ParObject::AddtoPack(*heatfluxdata, eheatflux);
     ParObject::AddtoPack(*tempgraddata, etempgrad);
  }
  // Calculate heatflux q and temperature gradients gradtemp at gauss points
  else if (action=="postproc_thermo_heatflux")
  {
    const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpheatfluxmap
       = params.get<Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpheatfluxmap");
    std::string heatfluxtype = params.get<std::string>("heatfluxtype","ndxyz");
    const int gid = ele->Id();
    LINALG::Matrix<nquad_,nsd_> gpheatflux(((*gpheatfluxmap)[gid])->A(),true);  // view only!

    // set views to components
    LINALG::Matrix<iel*numdofpernode_,1> efluxx(elevec1_epetra,true);  // view only!
    LINALG::Matrix<iel*numdofpernode_,1> efluxy(elevec2_epetra,true);  // view only!
    LINALG::Matrix<iel*numdofpernode_,1> efluxz(elevec3_epetra,true);  // view only!

    // catch unknown heatflux types
    bool processed = false;

    // nodally
    if ( (heatfluxtype=="ndxyz") or (heatfluxtype=="cxyz_ndxyz") )
    {
      processed = processed or true;
      // extrapolate stresses/strains at Gauss points to nodes
      ExtrapolateFromGaussPointsToNodes(ele, gpheatflux, efluxx, efluxy, efluxz);
    }

    // centered
    if ( (heatfluxtype=="cxyz") or (heatfluxtype=="cxyz_ndxyz") )
    {
      processed = processed or true;

      Teuchos::RCP<Epetra_MultiVector> eleheatflux = params.get<Teuchos::RCP<Epetra_MultiVector> >("eleheatflux");
      const Epetra_BlockMap& elemap = eleheatflux->Map();
      int lid = elemap.LID(gid);
      if (lid != -1)
      {
        for (int idim=0; idim<nsd_; ++idim)
        {
          //double& s = ; // resolve pointer for faster access
          double s = 0.0;
          for (int jquad=0; jquad<nquad_; ++jquad)
            s += gpheatflux(jquad,idim);
          s /= nquad_;
          (*((*eleheatflux)(idim)))[lid] = s;
        }
      }
    }

    // catch unknown heatflux types
    if (not processed)
      dserror("unknown type of stress/strain output on element level");
  }
/*  // Calculate heatflux q and temperature gradients gradtemp at gauss points
  else if (action=="calc_thermo_heatflux")
  {
   // extract local values from the global vectors
//   Teuchos::RefCountPtr<const Epetra_Vector> tempnp = discretization.GetState("tempnp");
//   if (tempnp==Teuchos::null)
//        dserror("Cannot get state vector 'tempnp'");
//   std::vector<double> mytempnp(lm.size());
//   DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,lm);

   // we always get an 3D flux vector for each node
   LINALG::Matrix<3,iel> eflux(true);

   // do a loop for systems of the temperatures
   // calculate heatflux vectors for actual temperature
   eflux.Clear();
   CalculateHeatFlux(eflux,ele);

   // assembly
   for (int k=0;k<iel;k++)
   { // form arithmetic mean of assembled nodal flux vectors
     // => factor is the number of adjacent elements for each node
     const double factor = ((ele->Nodes())[k])->NumElement();
     elevec1_epetra[k*numdofpernode_] += eflux(0,k)/factor;
     elevec2_epetra[k*numdofpernode_] += eflux(1,k)/factor;
     elevec3_epetra[k*numdofpernode_] += eflux(2,k)/factor;
   } // loop over elements

  }
  // Extrapolate heatflux q and temperature gradients gradtemp stored at gauss points
  // to element nodes
  else if (action=="postprocess_heatflux")
  {
    dserror("Action is not implemented: %s", action.c_str());
  } // calc_therm_heatflux
  else if (action =="calc_initial_time_deriv")
    // calculate time derivative for time value t_0
  {

    // get initial temperature rate values at the nodes
    const Teuchos::RCP<Epetra_MultiVector> rate = params.get< Teuchos::RCP<Epetra_MultiVector> >("temperature rate field",Teuchos::null);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,eratenp_,rate,nsd_);

    // need initial field
    Teuchos::RefCountPtr<const Epetra_Vector> temp0 = discretization.GetState("temp0");
    if (temp0==Teuchos::null)
      dserror("Cannot get state vector 'temp0'");

    // extract local values from the global vector
    std::vector<double> mytemp0(lm.size());
    DRT::UTILS::ExtractMyValues(*temp0,mytemp0,lm);

    // fill element arrays
    for (int i=0;i<iel;++i)
    {
      for (int k = 0; k< numdofpernode_; ++k)
      {
        // split for each temperature, insert into element arrays
        etempnp_[k](i,0) = mytemp0[k+(i*numdofpernode_)];
      } // for k
    } // for i

    // calculate matrix and rhs
    InitialTimeDerivative(
        ele,
        elemat1_epetra,
        elevec1_epetra
        );
  }*/
  else if (action=="integrate_shape_functions")
  {
    // calculate integral of shape functions
    const Epetra_IntSerialDenseVector dofids = params.get<Epetra_IntSerialDenseVector>("dofids");
    IntegrateShapeFunctions(ele,elevec1_epetra,dofids);
  }
  else if (action=="calc_thermo_update_istep")
  {
    ;  // do nothing
  }
  else
  {
    dserror("Unknown type of action for Temperature Implementation: %s",action.c_str());
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TemperImpl<distype>::EvaluateNeumann(
    DRT::Element*              ele,
    Teuchos::ParameterList&    params,
    DRT::Discretization&       discretization,
    std::vector<int>&          lm,
    Epetra_SerialDenseVector&  elevec1_epetra,
    Epetra_SerialDenseMatrix*  elemat1_epetra
    )
 {
  // check length
  if (lm.size() != iel*numdofpernode_)
    dserror("Location vector length does not match!");
  // set views
  LINALG::Matrix<iel*numdofpernode_,1> efext(elevec1_epetra,true);  // view only!
  // disassemble temperature
  if (discretization.HasState("temperature"))
  {
    std::vector<double> mytempnp(lm.size());
    Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState("temperature");
    if (tempnp==Teuchos::null)
      dserror("Cannot get state vector 'tempnp'");
    DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,lm);
    LINALG::Matrix<iel*numdofpernode_,1> etemp(&(mytempnp[0]),true);  // view only!
    etemp_.Update(etemp);  // copy
  }
  // check for the action parameter
  const std::string action = params.get<std::string>("action","none");
  // extract time
  const double time = params.get<double>("total time");

  // 26.10.09
  cout << "params \n" << params;

  // perform actions
  if (action=="calc_thermo_fext")
  {
    CalculateFintCondCapa(ele,
                          time,
                          NULL,
                          NULL,
                          NULL,
                          &efext,
                          NULL,NULL);
  }
  else
  {
    dserror("Unknown type of action for Temperature Implementation: %s",action.c_str());
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (public)                 g.bau 08/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateFintCondCapa(
  DRT::Element* ele,  ///< the element whose matrix is calculated
  const double& time,  ///< current time
  LINALG::Matrix<iel*numdofpernode_,iel*numdofpernode_>* etang,  ///< conductivity matrix
  LINALG::Matrix<iel*numdofpernode_,iel*numdofpernode_>* ecapa,  ///< capacity matrix
  LINALG::Matrix<iel*numdofpernode_,1>* efint,  ///< internal force
  LINALG::Matrix<iel*numdofpernode_,1>* efext,  ///< external force
  LINALG::Matrix<nquad_,nsd_>* eheatflux,  ///< heat fluxes at Gauss points
  LINALG::Matrix<nquad_,nsd_>* etempgrad  ///< temperature gradients at Gauss points
  )
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // ---------------------------------------------------------------------
  // call routine for calculation of radiation in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  if (efext != NULL) {
    Radiation(ele,time);
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_) dserror("Trouble with number of Gauss points");

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    // get radiation in gausspoint
    if (efext != NULL) {
      efext->MultiplyNN(fac_,funct_,radiation_,1.0);
    }

    // gradient of current temperature value
    gradtemp_.MultiplyNN(derxy_,etemp_);
    // store it
    if (etempgrad != NULL)
      for (int idim=0; idim<nsd_; ++idim)
        (*etempgrad)(iquad,idim) = gradtemp_(idim);

    // call material law => cmat_,heatflux_
    Materialize(ele);
    // store heatflux
    if (eheatflux != NULL)
      for (int idim=0; idim<nsd_; ++idim)
        (*eheatflux)(iquad,idim) = heatflux_(idim);

    // internal force vector
    if (efint != NULL) {
      efint->MultiplyTN(fac_,derxy_,heatflux_,1.0);
    }

    // conductivity matrix
    if (etang != NULL) {
      LINALG::Matrix<nsd_,iel> aop(false);
      aop.MultiplyNN(cmat_,derxy_);
      etang->MultiplyTN(fac_,derxy_,aop,1.0);
    }

    // capacity matrix
    if (ecapa != NULL) {
      ecapa->MultiplyNT(fac_*capacoeff_,funct_,funct_,1.0);
    }

  } // integration loop
}

/*----------------------------------------------------------------------*
 |  get the radiation  (private)                             dano 09/09 |
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
 |  get the material                                         dano 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::Materialize(
  const DRT::Element* ele
  )
{
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  // get Fourier´s law
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
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
    const int                                    iquad,      ///< id of current Gauss point
    const int                                    eleid       ///< the element id
    )
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
    {xsi_(idim) = gpcoord[idim];}

  // shape functions and their first derivatives
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

  xjm_.MultiplyNT(deriv_,xyze_);
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
 | calculate mass matrix + rhs for determ. initial time deriv. gjb 08/08|
 *----------------------------------------------------------------------*/
/*
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::InitialTimeDerivative(
    DRT::Element*                         ele,
    Epetra_SerialDenseMatrix&             emat,
    Epetra_SerialDenseVector&             erhs
    )
 {
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // dead load in element nodes at initial point in time
  Radiation(ele);

   // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);

  //----------------------------------------------------------------------
  // element integration loop
  //----------------------------------------------------------------------
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------

    //------------ get values of variables at integration point
    for (int k=0;k<numdofpernode_;++k) // deal with a system of transported scalars
    {
      // get radiation in gausspoint
      rhs_[k] = bodyforce_[k].Dot(funct_);

      // time dependend part: rho*u_x*N,x+ rho*u_y*N,y
      conv_.MultiplyTN(derxy_,rateint_);

       // get value of current temperature
      conint_[k] = funct_.Dot(etempnp_[k]);

      // gradient of current temperature value
      gradtemp_.Multiply(derxy_,etempnp_[k]);

      //----------------------------------------------------------------
      // element matrix: transient term
      //----------------------------------------------------------------
      // transient term
      for (int vi=0; vi<iel; ++vi)
      {
        const double v = fac_*funct_(vi);
        const int fvi = vi*numdofpernode_+k;

        for (int ui=0; ui<iel; ++ui)
        {
          const int fui = ui*numdofpernode_+k;

          emat(fvi,fui) += v*densfunct_(ui);
        }
      }

      //----------------------------------------------------------------
      // element right hand side: convective term in convective form
      //----------------------------------------------------------------
      double vrhs = fac_*conv_etemp0_k;
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] -= vrhs*funct_(vi);
      }

      //----------------------------------------------------------------
      // element right hand side: radiation term
      //----------------------------------------------------------------
      vrhs = fac_*rad_[k];
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi = vi*numdofpernode_+k;

        erhs[fvi] += vrhs*funct_(vi);
      }
    } // loop over each scalar k

  } // integration loop

  return;
 } // TemperImpl::InitialTimeDerivative
*/


/*----------------------------------------------------------------------*
 |  calculate heatflux                             (private) dano 11/09 |
 *----------------------------------------------------------------------*/
/*
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::CalculateFlux(
    LINALG::Matrix<3,iel>&          flux,
    const DRT::Element*             ele,
    const std::vector<double>&      etempnp,
    const Epetra_SerialDenseVector& erate
    )
  {
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // get material constants
  double capa(0.0); // capa is not a double?? 10.11.09
  double conduct(0.0);

  // get the material
  Teuchos::RefCountPtr<MAT::Material> material = ele->Material();

  // use one-point Gauss rule to do calculations at element center
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtIntPoint(intpoints,0,ele->Id());

  if (material->MaterialType() == INPAR::MAT::m_th_fourier_iso)
  {
    const MAT::FourierIso* actmat = static_cast<const MAT::FourierIso*>(material.get());

    // get material constants ???? 10.11.09
    capa = actmat->Capacity();
    conduct = actmat->conduct_;
  }
  else
    dserror("Material type is not supported");

  //----------------------------------------- declaration of variables
  LINALG::SerialDenseMatrix nodecoords;
  nodecoords = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  if ((int) nodecoords.N() != iel) dserror("number of nodes does not match");

  // loop over all nodes
  for (int iquad=0; iquad<iel; ++iquad)
  {
    // reference coordinates of the current node
    for (int idim=0;idim<nsd_;idim++)
      {xsi_(idim) = nodecoords(idim, iquad);}

    // first derivatives
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

    // compute Jacobian matrix and determinant
    // actually compute its transpose....
    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xij_.Invert(xjm_);

    if (det < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f",ele->Id(), det);

    // compute global derivatives
    derxy_.Multiply(xij_,deriv_);

//    // gradient of electric potential
//    gradpot_.Clear();
//    if (frt > 0.0) // ELCH
//    {
//      for (int k=0;k<iel;k++)
//      {
//        for (int idim=0; idim<nsd_ ;idim++)
//        {
//          gradpot_(idim) += derxy_(idim,k)*ephinp[k*numdofpernode_+numscal_];
//        }
//      }
//    }
//
//    const double ephinpatnode = ephinp[iquad*numdofpernode_+dofindex];
//    // add different flux contributions as specified by user input
//    switch (fluxtype)
//    {
//      case INPAR::SCATRA::flux_total_domain:
//        if (frt > 0.0) // ELCH
//        {
//          // migration flux terms
//          for (int idim=0; idim<nsd_ ;idim++)
//          {
//            flux(idim,iquad) += diffus_valence_frt*gradpot_(idim)*ephinpatnode;
//          }
//        }
//        // convective flux terms
//        for (int idim=0; idim<nsd_ ;idim++)
//        {
//          flux(idim,iquad) -= dens*evel[idim+iquad*nsd_]*ephinpatnode;
//        }
//        // no break statement here!
//      case INPAR::SCATRA::flux_diffusive_domain:
//        //diffusive flux terms
//        for (int k=0;k<iel;k++)
//        {
//          for (int idim=0; idim<nsd_ ;idim++)
//          {
//            flux(idim,iquad) += diffus*derxy_(idim,k)*ephinp[k*numdofpernode_+dofindex];
//          }
//        }
//        break;
//      default:
//        dserror("received illegal flag inside flux evaluation for whole domain");
//    };

    //set zeros for unused space dimenions
    for (int idim=nsd_; idim<3; idim++)
    {
      flux(idim,iquad) = 0.0;
    }
  } // loop over nodes

  return;
 } // TemperImpl::CalculateFlux
*/

/*----------------------------------------------------------------------*
 |  Integrate shape functions over domain (private)           gjb 07/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::IntegrateShapeFunctions(
    const DRT::Element*                ele,
    Epetra_SerialDenseVector&          elevec1,
    const Epetra_IntSerialDenseVector& dofids
    )
  {
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,iel> >(ele,xyze_);

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,gpid,ele->Id());

    // compute integral of shape functions (only for dofid)
    for (int k=0;k<numdofpernode_;k++)
    {
      if (dofids[k] >= 0)
      {
        for (int node=0;node<iel;node++)
        {
          elevec1[node*numdofpernode_+k] += funct_(node) * fac_;
        }
      }
    }

  } //loop over integration points

  return;
 } //TemperImpl<distype>::IntegrateShapeFunction

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TemperImpl<distype>::ExtrapolateFromGaussPointsToNodes(
  DRT::Element* ele,  ///< the element whose matrix is calculated
  const LINALG::Matrix<nquad_,nsd_>& gpheatflux,
  LINALG::Matrix<iel*numdofpernode_,1>& efluxx,
  LINALG::Matrix<iel*numdofpernode_,1>& efluxy,
  LINALG::Matrix<iel*numdofpernode_,1>& efluxz
  )
{
  // this quick'n'dirty hack functions only for hex8
  if ( not ( (distype == DRT::Element::hex8)
             or (distype == DRT::Element::quad4)
             or (distype == DRT::Element::line2) ) )
    dserror("Sorry, not implemented for element shape");

  // another check
  if (iel*numdofpernode_ != nquad_)
    dserror("Works only if number of gauss points and nodes match");

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
  if (intpoints.IP().nquad != nquad_) dserror("Trouble with number of Gauss points");

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

    for (int inode=0; inode<iel; ++inode)
      shpfctatgps(iquad,inode) = funct_(inode);
  }

  // extrapolation
  LINALG::Matrix<nquad_,nsd_> ndheatflux;
  LINALG::Matrix<nquad_,nsd_> gpheatflux2(gpheatflux);
  {
    LINALG::FixedSizeSerialDenseSolver<nquad_,nquad_,nsd_> solver;
    solver.SetMatrix(shpfctatgps);
    solver.SetVectors(ndheatflux,gpheatflux2);
    solver.Solve();
  }

  // copy into component vectors
  for (int idof=0; idof<iel*numdofpernode_; ++idof)
  {
    efluxx(idof) = ndheatflux(idof,0);
    if (nsd_>1) efluxy(idof) = ndheatflux(idof,1);
    if (nsd_>2) efluxz(idof) = ndheatflux(idof,2);
  }

  // bye
  return;
}

#endif // D_THERMO
#endif // CCADISCRET
