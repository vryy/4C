/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro_p2_evaluate.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#include "wall1_poro_p2.H"
#include "wall1_poro_p2_eletypes.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fem_general/drt_utils_gder2.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP2<distype>::ComputePorosityAndLinearization
(   Teuchos::ParameterList&                 params,
    const double&                           press,
    const double&                           J,
    const int&                              gp,
    const LINALG::Matrix<my::numnod_,1>&    shapfct,
    const LINALG::Matrix<my::numnod_,1>*    myporosity,
    const LINALG::Matrix<1,my::numdof_>&    dJ_dus,
    double &                                porosity,
    LINALG::Matrix<1,my::numdof_>&          dphi_dus)
{
  if(myporosity == NULL)
    dserror("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);

  dphi_dus.PutScalar( 0.0 );

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP2<distype>::ComputePorosityAndLinearizationOD
(   Teuchos::ParameterList&          params,
    const double&                    press,
    const double&                    J,
    const int&                       gp,
    const LINALG::Matrix<my::numnod_,1>&    shapfct,
    const LINALG::Matrix<my::numnod_,1>*    myporosity,
    double &                         porosity,
    double &                         dphi_dp)
{
  dphi_dp = 0.0;

  if(myporosity == NULL)
    dserror("no porosity values given!");
  else
    porosity = shapfct.Dot(*myporosity);

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (protected)                                       |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_PoroP2<distype>::MyEvaluate(
                                    Teuchos::ParameterList&      params,
                                    DRT::Discretization&         discretization,
                                    DRT::Element::LocationArray& la,
                                    Epetra_SerialDenseMatrix&    elemat1_epetra,
                                    Epetra_SerialDenseMatrix&    elemat2_epetra,
                                    Epetra_SerialDenseVector&    elevec1_epetra,
                                    Epetra_SerialDenseVector&    elevec2_epetra,
                                    Epetra_SerialDenseVector&    elevec3_epetra
                                    )
{
  // start with "none"
  typename my::ActionType act = my::none;

  // get the required action
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_internalforce")         act = my::calc_struct_internalforce;
  else if (action=="calc_struct_nlnstiff")              act = my::calc_struct_nlnstiff;
  else if (action=="calc_struct_nlnstiffmass")          act = my::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_multidofsetcoupling")   act = my::calc_struct_multidofsetcoupling;
  else if (action=="calc_struct_poroscatracoupling")    act = my::calc_struct_poroscatra_coupling;
  //else if (action=="calc_struct_stress")                act = calc_struct_stress;
  //else if (action=="postprocess_stress")                act = postprocess_stress;

  // what should the element do
  switch(act)
  {
  //==================================================================================
  // nonlinear stiffness, damping and internal force vector for poroelasticity
  case my::calc_struct_nlnstiff:
  case my::calc_struct_nlnstiffmass:
  {
    if(la.Size()>1)
    {
      // stiffness
      LINALG::Matrix<numdof_,numdof_> elemat1(elemat1_epetra.A(),true);
      //damping
      LINALG::Matrix<numdof_,numdof_> elemat2(elemat2_epetra.A(),true);
      // internal force vector
      LINALG::Matrix<numdof_,1> elevec1(elevec1_epetra.A(),true);
      //LINALG::Matrix<numdof_,1> elevec2(elevec2_epetra.A(),true);
      // elemat2,elevec2+3 are not used anyway

      // build the location vector only for the structure field
      std::vector<int> lm = la[0].lm_;

      LINALG::Matrix<my::numdim_,my::numnod_> mydisp(true);
      LINALG::Matrix<my::numnod_,1> myporosity(true);
      my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, NULL,"displacement");
      my::ExtractValuesFromGlobalVector(discretization,2,la[2].lm_, NULL, &myporosity,"scalar");

      LINALG::Matrix<numdof_,numdof_>* matptr = NULL;
      if (elemat1.IsInitialized())   matptr = &elemat1;

      enum INPAR::STR::DampKind damping = params.get<enum INPAR::STR::DampKind>("damping",INPAR::STR::damp_none);
      LINALG::Matrix<numdof_,numdof_>* matptr2 = NULL;
      if (elemat2.IsInitialized() and (damping==INPAR::STR::damp_material) )  matptr2 = &elemat2;

      // need current fluid state,
      // call the fluid discretization: fluid equates 2nd dofset
      // disassemble velocities and pressures

      LINALG::Matrix<my::numdim_,my::numnod_> myvel(true);

      LINALG::Matrix<my::numdim_,my::numnod_> myfluidvel(true);
      LINALG::Matrix<my::numnod_,1> myepreaf(true);

      if (discretization.HasState(0,"velocity"))
        my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &myvel, NULL,"velocity");

      if (discretization.HasState(1,"fluidvel"))
        // extract local values of the global vectors
        my::ExtractValuesFromGlobalVector(discretization,1,la[1].lm_, &myfluidvel, &myepreaf,"fluidvel");

      //calculate tangent stiffness matrix
      nlnstiff_poroelast(lm,mydisp,myvel,&myporosity,myfluidvel,myepreaf,matptr,matptr2,&elevec1,params);
    }

  }
  break;

  //==================================================================================
  // coupling terms in force-vector and stiffness matrix for poroelasticity
  case my::calc_struct_multidofsetcoupling:
  {
    // stiffness
    LINALG::Matrix<numdof_,(my::numdim_+1)*my::numnod_> elemat1(elemat1_epetra.A(),true);
    //LINALG::Matrix<numdof_,(numdim_+1)*numnod_> elemat2(elemat2_epetra.A(),true);

    // internal force vector
    //LINALG::Matrix<numdof_,1> elevec1(elevec1_epetra.A(),true);
    //LINALG::Matrix<numdof_,1> elevec2(elevec2_epetra.A(),true);

    // elemat2,elevec2+3 are not used anyway

    // build the location vector only for the structure field
    std::vector<int> lm = la[0].lm_;

    LINALG::Matrix<numdof_,(my::numdim_+1)*my::numnod_>* matptr = NULL;
    if (elemat1.IsInitialized()) matptr = &elemat1;

    // need current fluid state,
    // call the fluid discretization: fluid equates 2nd dofset
    // disassemble velocities and pressures
    if (discretization.HasState(1,"fluidvel"))
    {
      LINALG::Matrix<my::numdim_,my::numnod_> myvel(true);
      LINALG::Matrix<my::numdim_,my::numnod_> myfluidvel(true);
      LINALG::Matrix<my::numnod_,1> myepreaf(true);

      LINALG::Matrix<my::numdim_,my::numnod_> mydisp(true);
      LINALG::Matrix<my::numnod_,1> myporosity(true);
      my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, NULL,"displacement");
      my::ExtractValuesFromGlobalVector(discretization,2,la[2].lm_, NULL, &myporosity,"scalar");

      if (discretization.HasState(0,"velocity"))
        my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &myvel, NULL,"velocity");

      if (discretization.HasState(1,"fluidvel"))
        // extract local values of the global vectors
        my::ExtractValuesFromGlobalVector(discretization,1,la[1].lm_, &myfluidvel, &myepreaf,"fluidvel");

      coupling_poroelast(lm,mydisp,myvel,&myporosity,myfluidvel,myepreaf,matptr,//NULL,
          NULL,NULL,params);
    }

  }
  break;

  //==================================================================================
  // coupling terms in force-vector and stiffness matrix for poroelasticity
  case my::calc_struct_poroscatra_coupling:
  {
    // stiffness
    LINALG::Matrix<numdof_, my::numnod_> elemat1(elemat1_epetra.A(),true);
    //LINALG::Matrix<numdof_,(numdim_+1)*numnod_> elemat2(elemat2_epetra.A(),true);

    // internal force vector
    //LINALG::Matrix<numdof_,1> elevec1(elevec1_epetra.A(),true);
    //LINALG::Matrix<numdof_,1> elevec2(elevec2_epetra.A(),true);

    // elemat2,elevec2+3 are not used anyway

    // build the location vector only for the structure field
    std::vector<int> lm = la[0].lm_;

    LINALG::Matrix<numdof_,my::numnod_>* matptr = NULL;
    if (elemat1.IsInitialized()) matptr = &elemat1;

    // need current fluid state,
    // call the fluid discretization: fluid equates 2nd dofset
    // disassemble velocities and pressures
    if (discretization.HasState(1,"fluidvel"))
    {
      LINALG::Matrix<my::numdim_,my::numnod_> myvel(true);
      LINALG::Matrix<my::numdim_,my::numnod_> myfluidvel(true);
      LINALG::Matrix<my::numnod_,1> myepreaf(true);

      LINALG::Matrix<my::numdim_,my::numnod_> mydisp(true);
      LINALG::Matrix<my::numnod_,1> myporosity(true);
      my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, NULL,"displacement");
      my::ExtractValuesFromGlobalVector(discretization,2,la[2].lm_, NULL, &myporosity,"scalar");

      if (discretization.HasState(0,"velocity"))
        my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &myvel, NULL,"velocity");

      if (discretization.HasState(1,"fluidvel"))
        // extract local values of the global vectors
        my::ExtractValuesFromGlobalVector(discretization,1,la[1].lm_, &myfluidvel, &myepreaf,"fluidvel");

      coupling_poroscatra(lm,mydisp,myvel,&myporosity,myfluidvel,myepreaf,matptr,//NULL,
          NULL,NULL,params);
    }

  }
  break;

  //==================================================================================
  // nonlinear stiffness and internal force vector for poroelasticity
  case my::calc_struct_internalforce:
  {
    // internal force vector
    LINALG::Matrix<numdof_,1> elevec1(elevec1_epetra.A(),true);
    // elemat2,elevec2+3 are not used anyway

    // build the location vector only for the structure field
    std::vector<int> lm = la[0].lm_;

    LINALG::Matrix<my::numdim_,my::numnod_> mydisp(true);
    LINALG::Matrix<my::numnod_,1> myporosity(true);
    my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, NULL,"displacement");
    my::ExtractValuesFromGlobalVector(discretization,2,la[2].lm_, NULL, &myporosity,"scalar");

    LINALG::Matrix<my::numdim_,my::numnod_> myvel(true);

    LINALG::Matrix<my::numdim_,my::numnod_> myfluidvel(true);
    LINALG::Matrix<my::numnod_,1> myepreaf(true);

    // need current fluid state,
    // call the fluid discretization: fluid equates 2nd dofset
    // disassemble velocities and pressures
    if (discretization.HasState(1,"fluidvel"))
    {
      // extract local values of the global vectors
      my::ExtractValuesFromGlobalVector(discretization,1,la[1].lm_, &myfluidvel, &myepreaf,"fluidvel");

      my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &myvel, NULL,"velocity");

      nlnstiff_poroelast(lm,mydisp,myvel,&myporosity,myfluidvel,myepreaf,NULL,NULL,&elevec1,params);
    }
  }
  break;
  //==================================================================================
  default:
    //do nothing (no error because there are some actions the poro element is supposed to ignore)
    break;
  } // action
  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP2<distype>::InitElement()
{
  //initialize base element
  my::InitElement();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP2<distype>::nlnstiff_poroelast(
    std::vector<int>& lm, ///< location matrix
    LINALG::Matrix<my::numdim_, my::numnod_>&   disp,         ///< current displacements
    LINALG::Matrix<my::numdim_, my::numnod_>&   vel,          ///< current velocities
    LINALG::Matrix<my::numnod_, 1>*             porosity_dof,
    LINALG::Matrix<my::numdim_,my::numnod_> &   evelnp,       ///< fluid velocity of element
    LINALG::Matrix<my::numnod_,1> &             epreaf,       ///< fluid pressure of element
    LINALG::Matrix<numdof_,numdof_>*            stiffmatrix,  ///< element stiffness matrix
    LINALG::Matrix<numdof_,numdof_>*            reamatrix,    // element reactive matrix
    LINALG::Matrix<numdof_,1>*                  force,        ///< element internal force vector
    Teuchos::ParameterList&                     params        ///< algorithmic parameters e.g. time
    )
{

  my::GetMaterials();

  // update element geometry
  LINALG::Matrix<my::numdim_,my::numnod_> xrefe; // material coord. of element
  LINALG::Matrix<my::numdim_,my::numnod_> xcurr; // current  coord. of element

  DRT::Node** nodes = my::Nodes();
  for (int i=0; i<my::numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for(int j=0; j<my::numdim_;j++)
    {
      xrefe(j,i) = x[j];
      xcurr(j,i) = xrefe(j,i) + disp(j,i);
    }
  }

  LINALG::Matrix<numdof_,numdof_> erea_v(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  my::GaussPointLoop(  params,
                        xrefe,
                        xcurr,
                        disp,
                        vel,
                        evelnp,
                        epreaf,
                        porosity_dof,
                        erea_v,
                        stiffmatrix,
                        reamatrix,
                        force);

  if ( reamatrix != NULL )
  {
    /* additional "reactive darcy-term"
     detJ * w(gp) * ( J * reacoeff * phi^2  ) * D(v_s)
     */
    reamatrix->Update(1.0,erea_v,1.0);
  }
        //if the reaction part is not supposed to be computed separately, we add it to the stiffness
    //(this is not the best way to do it, but it only happens once during initialization)
//  else
//  {
//   const double dt = params.get<double>("delta time");
//      stiffmatrix->Update(1.0/dt,erea_v,1.0);
//  }

  return;
}  // nlnstiff_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP2<distype>::coupling_poroelast(
    std::vector<int>&                                         lm,           // location matrix
    LINALG::Matrix<my::numdim_, my::numnod_>&                 disp,         // current displacements
    LINALG::Matrix<my::numdim_, my::numnod_>&                 vel,          // current velocities
    LINALG::Matrix<my::numnod_, 1>*                           porosity,
    LINALG::Matrix<my::numdim_, my::numnod_> &                evelnp,       //current fluid velocity
    LINALG::Matrix<my::numnod_, 1> &                          epreaf,       //current fluid pressure
    LINALG::Matrix<numdof_, (my::numdim_ + 1) * my::numnod_>* stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<numdof_, (my::numdim_ + 1) * my::numnod_>* reamatrix,    // element reactive matrix
    LINALG::Matrix<numdof_, 1>*                               force,        // element internal force vector
    Teuchos::ParameterList&                                   params)       // algorithmic parameters e.g. time
{
  my::GetMaterials();

  //get structure material
  if(my::Material()==Teuchos::null) dserror("no structure material available!");
  my::structmat_ = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(my::Material());
  if(my::structmat_->MaterialType() != INPAR::MAT::m_structporo and
      my::structmat_->MaterialType() != INPAR::MAT::m_structpororeaction)
    dserror("invalid structure material for poroelasticity");

  //=======================================================================

  // update element geometry
  LINALG::Matrix<my::numdim_,my::numnod_> xrefe; // material coord. of element
  LINALG::Matrix<my::numdim_,my::numnod_> xcurr; // current  coord. of element

  DRT::Node** nodes = my::Nodes();
  for (int i=0; i<my::numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for(int j=0; j<my::numdim_;j++)
    {
      xrefe(j,i) = x[j];
      xcurr(j,i) = xrefe(j,i) + disp(j,i);
    }
  }

  LINALG::Matrix<numdof_,(my::numdim_+1)*my::numnod_> ecoupl(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  my::GaussPointLoopOD( params,
                         xrefe,
                         xcurr,
                         disp,
                         vel,
                         evelnp,
                         epreaf,
                         porosity,
                         ecoupl );

  if (stiffmatrix != NULL)
  {//todo
    // build tangent coupling matrix : effective dynamic stiffness coupling matrix
    //    K_{Teffdyn} = 1/dt C
    //                + theta K_{T}
    double theta = params.get<double>("theta");
    stiffmatrix->Update(theta,ecoupl,1.0);
  }

  return;

}  // coupling_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP2<distype>::coupling_poroscatra(
    std::vector<int>&                                         lm,           // location matrix
    LINALG::Matrix<my::numdim_, my::numnod_>&                 disp,         // current displacements
    LINALG::Matrix<my::numdim_, my::numnod_>&                 vel,          // current velocities
    LINALG::Matrix<my::numnod_, 1>*                           porosity,
    LINALG::Matrix<my::numdim_, my::numnod_> &                evelnp,       //current fluid velocity
    LINALG::Matrix<my::numnod_, 1> &                          epreaf,       //current fluid pressure
    LINALG::Matrix<numdof_,  my::numnod_>* stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<numdof_,  my::numnod_>* reamatrix,    // element reactive matrix
    LINALG::Matrix<numdof_, 1>*                               force,        // element internal force vector
    Teuchos::ParameterList&                                   params)       // algorithmic parameters e.g. time
{
  my::GetMaterials();

  //get structure material
  if(my::Material()==Teuchos::null) dserror("no structure material available!");
  my::structmat_ = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(my::Material());
  if(my::structmat_->MaterialType() != INPAR::MAT::m_structporo and
      my::structmat_->MaterialType() != INPAR::MAT::m_structpororeaction)
    dserror("invalid structure material for poroelasticity");

  //=======================================================================

  // update element geometry
  LINALG::Matrix<my::numdim_,my::numnod_> xrefe; // material coord. of element
  LINALG::Matrix<my::numdim_,my::numnod_> xcurr; // current  coord. of element

  DRT::Node** nodes = my::Nodes();
  for (int i=0; i<my::numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    for(int j=0; j<my::numdim_;j++)
    {
      xrefe(j,i) = x[j];
      xcurr(j,i) = xrefe(j,i) + disp(j,i);
    }
  }

  LINALG::Matrix<numdof_,my::numnod_> ecoupl(true);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  GaussPointLoopPoroScatraOD( params,
                         xrefe,
                         xcurr,
                         disp,
                         vel,
                         evelnp,
                         epreaf,
                         porosity,
                         ecoupl );

  if (stiffmatrix != NULL)
  {//todo
    // build tangent coupling matrix : effective dynamic stiffness coupling matrix
    //    K_{Teffdyn} = 1/dt C
    //                + theta K_{T}
    double theta = params.get<double>("theta");
    stiffmatrix->Update(theta,ecoupl,1.0);
  }

  return;

}  // coupling_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate constitutive relation for porosity and compute derivatives (protected) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP2<distype>::ConsitutiveDerivatives(
                                              Teuchos::ParameterList& params,       ///< (i) parameter list
                                              double     press,   ///< (i) fluid pressure at gauss point
                                              double     J,       ///< (i) Jacobian determinant at gauss point
                                              double     porosity,///< (i) porosity at gauss point
                                              double*    dW_dp,   ///< (o) derivative of potential w.r.t. pressure
                                              double*    dW_dphi, ///< (o) derivative of potential w.r.t. porosity
                                              double*    dW_dJ,   ///< (o) derivative of potential w.r.t. jacobian
                                              double*    W        ///< (o) inner potential
                                              )
{
  my::structmat_->ConsitutiveDerivatives(
    params,
    press,
    J,
    porosity,
    dW_dp,
    dW_dphi,
    dW_dJ,
    W
    );
  return;
};

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP2<distype>::GaussPointLoopPoroScatraOD(
    Teuchos::ParameterList&                 params,
    const LINALG::Matrix<my::numdim_,my::numnod_>&  xrefe,
    const LINALG::Matrix<my::numdim_,my::numnod_>&  xcurr,
    const LINALG::Matrix<my::numdim_,my::numnod_>&  nodaldisp,
    const LINALG::Matrix<my::numdim_,my::numnod_>&  nodalvel,
    const LINALG::Matrix<my::numdim_,my::numnod_>&  evelnp,
    const LINALG::Matrix<my::numnod_,1> &       epreaf,
    const LINALG::Matrix<my::numnod_, 1>*       porosity_dof,
    LINALG::Matrix<numdof_, my::numnod_>& stiffmatrix
                                        )
{

  LINALG::Matrix<my::numdim_,my::numnod_> N_XYZ;       //  first derivatives at gausspoint w.r.t. X, Y,Z
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<my::numdim_,my::numdim_> defgrd(true); //  deformation gradiant evaluated at gauss point
  LINALG::Matrix<my::numnod_,1> shapefct;           //  shape functions evalulated at gauss point
  LINALG::Matrix<my::numdim_,my::numnod_> deriv(true);  //  first derivatives at gausspoint w.r.t. r,s,t
  //LINALG::Matrix<numdim_,1> xsi;
  for (int gp=0; gp<my::numgpt_; ++gp)
  {
    //evaluate shape functions and derivatives at integration point
    this->ComputeShapeFunctionsAndDerivatives(gp,shapefct,deriv,N_XYZ);

    const double J = this->ComputeJacobianDeterminant(gp,xcurr,deriv);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    defgrd.MultiplyNT(xcurr,N_XYZ); //  (6.17)

    // inverse deformation gradient F^-1
    LINALG::Matrix<my::numdim_,my::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    //---------------- get pressure at integration point
    double press = shapefct.Dot(epreaf);

    //------------------ get material pressure gradient at integration point
    LINALG::Matrix<my::numdim_,1> Gradp;
    Gradp.Multiply(N_XYZ,epreaf);

    //--------------------- get fluid velocity at integration point
    LINALG::Matrix<my::numdim_,1> fvelint;
    fvelint.Multiply(evelnp,shapefct);

    //! ----------------structure displacement and velocity at integration point
    LINALG::Matrix<my::numdim_,1> velint(true);
    for(int i=0; i<my::numnod_; i++)
      for(int j=0; j<my::numdim_; j++)
        velint(j) += nodalvel(j,i) * shapefct(i);

    //F^-T * grad p
    LINALG::Matrix<my::numdim_,1> Finvgradp;
    Finvgradp.MultiplyTN(defgrd_inv, Gradp);

    //**************************************************+auxilary variables for computing the porosity and linearization
    double dphi_dp=0.0;
    double porosity=0.0;

    ComputePorosityAndLinearizationOD(params,
                                      press,
                                      J,
                                      gp,
                                      shapefct,
                                      porosity_dof,
                                      porosity,
                                      dphi_dp);

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++

    const double reacoeff = my::fluidmat_->ComputeReactionCoeff();
    const double detJ_w = my::detJ_[gp]*my::intpoints_.Weight(gp);

    for (int i=0; i<my::numnod_; i++)
    {
      const int fi = my::numdim_*i;
      const double fac = detJ_w* shapefct(i);

      for(int j=0; j<my::numdim_; j++)
      {
        for(int k=0; k<my::numnod_; k++)
        {

          /*-------structure- porosity coupling: "stress terms" + "pressure gradient terms"
           -B^T . ( -1*J*C^-1 ) * Dp
           - J * F^-T * Grad(p) * Dphi - J * F^-T * d(Grad((p))/(dp) * phi * Dp
           */
          stiffmatrix(fi+j, k) += - fac * J *  Finvgradp(j) * shapefct(k);

          /*-------structure- porosity coupling:  "dracy-terms" + "reactive darcy-terms"
           - 2 * reacoeff * J * v^f * phi *  Dphi
           + 2 * reacoeff * J * v^s * phi *  Dphi
           */
          const double tmp = fac * reacoeff * J * 2 * porosity * shapefct(k);
          stiffmatrix(fi+j, k) += -tmp * fvelint(j);

          stiffmatrix(fi+j, k) += tmp * velint(j);
        }
      }
    }

    if(my::fluidmat_->Type() == "Darcy-Brinkman")
    {
      const double visc = my::fluidmat_->Viscosity();

      // non-linear B-operator
      LINALG::Matrix<my::numstr_,numdof_> bop;
      this->ComputeBOperator(bop,defgrd,N_XYZ);

      // material fluid velocity gradient at integration point
      LINALG::Matrix<my::numdim_,my::numdim_>              fvelder;
      fvelder.MultiplyNT(evelnp,N_XYZ);

      // Right Cauchy-Green tensor = F^T * F
      LINALG::Matrix<my::numdim_,my::numdim_> cauchygreen;
      cauchygreen.MultiplyTN(defgrd,defgrd);

      // inverse Right Cauchy-Green tensor
      LINALG::Matrix<my::numdim_,my::numdim_> C_inv(false);
      C_inv.Invert(cauchygreen);

      LINALG::Matrix<my::numstr_,1> fstress;

      LINALG::Matrix<my::numdim_,my::numdim_> CinvFvel;
      LINALG::Matrix<my::numdim_,my::numdim_> tmp;
      CinvFvel.Multiply(C_inv,fvelder);
      tmp.MultiplyNT(CinvFvel,defgrd_inv);
      LINALG::Matrix<my::numdim_,my::numdim_> tmp2(tmp);
      tmp.UpdateT(1.0,tmp2,1.0);

      fstress(0) = tmp(0,0);
      fstress(1) = tmp(1,1);
      fstress(2) = tmp(0,1);

      //B^T . \sigma
      LINALG::Matrix<numdof_,1> fstressb;
      fstressb.MultiplyTN(bop,fstress);
      LINALG::Matrix<my::numdim_,my::numnod_> N_XYZ_Finv;
      N_XYZ_Finv.Multiply(defgrd_inv,N_XYZ);

      for (int i=0; i<my::numnod_; i++)
      {
        const int fi = my::numdim_*i;

        for(int j=0; j<my::numdim_; j++)
        {
          for(int k=0; k<my::numnod_; k++)
          {

            /*-------structure- porosity coupling: "darcy-brinkman stress terms"
             B^T . ( \mu*J - d(phi)/(dp) * fstress ) * Dp
             */
            stiffmatrix(fi+j, k) += detJ_w * fstressb(fi+j) * visc * J * shapefct(k);
          }
        }
      }
    }//darcy-brinkman
    /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::Wall1_PoroP2<DRT::Element::quad4>;
template class DRT::ELEMENTS::Wall1_PoroP2<DRT::Element::quad9>;


