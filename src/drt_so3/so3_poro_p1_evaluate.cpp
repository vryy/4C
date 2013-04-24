/*----------------------------------------------------------------------*/
/*!
 \file so3_poro_p1_evaluate.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include "so3_poro_p1.H"
#include "so3_poro_p1_eletypes.H"

#include "../drt_lib/drt_utils.H"

#include "../drt_mat/fluidporo.H"
#include "../drt_mat/structporo.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fem_general/drt_utils_gder2.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::ComputePorosityAndLinearization
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
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::ComputePorosityAndLinearizationOD
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
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Poro_P1< so3_ele, distype>::Evaluate(Teuchos::ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    DRT::Element::LocationArray& la,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  // start with "none"
  typename my::ActionType act = my::none;

  // get the required action
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_multidofsetcoupling")   act = my::calc_struct_multidofsetcoupling;

  // what should the element do
  switch(act)
  {
  //==================================================================================
  // off diagonal terms in stiffness matrix for monolithic coupling
  case my::calc_struct_multidofsetcoupling:
  {
    MyEvaluate(params,
                      discretization,
                      la,
                      elemat1_epetra,
                      elemat2_epetra,
                      elevec1_epetra,
                      elevec2_epetra,
                      elevec3_epetra);
  }
  break;
  //==================================================================================
  default:
  {
    //in some cases we need to write/change some data before evaluating
    my::PreEvaluate(params,
                discretization,
                la);

    Epetra_SerialDenseMatrix elemat1_sub;
    Epetra_SerialDenseMatrix elemat2_sub;
    Epetra_SerialDenseVector elevec1_sub;
    Epetra_SerialDenseVector elevec2_sub;
    Epetra_SerialDenseVector elevec3_sub;

    if(elemat1_epetra.A())
      elemat1_sub.Shape(my::numdof_,my::numdof_);
    if(elemat2_epetra.A())
      elemat2_sub.Shape(my::numdof_,my::numdof_);
    if(elevec1_epetra.A())
      elevec1_sub.Resize(my::numdof_);
    if(elevec2_epetra.A())
      elevec2_sub.Resize(my::numdof_);
    if(elevec3_epetra.A())
      elevec3_sub.Resize(my::numdof_);

    std::vector<int> lm_sub;
    for(int i=0;i<my::numnod_;i++)
      for(int j=0;j<my::numdim_;j++)
        lm_sub.push_back(la[0].lm_[i*noddof_+j]);

    //evaluate parent solid element
    so3_ele::Evaluate(params,
                      discretization,
                      lm_sub,
                      elemat1_sub,
                      elemat2_sub,
                      elevec1_sub,
                      elevec2_sub,
                      elevec3_sub);

     if(elemat1_epetra.A())
     for(int i=0;i<my::numnod_;i++)
       for(int j=0;j<my::numdim_;j++)
         for(int k=0;k<my::numnod_;k++)
           for(int l=0;l<my::numdim_;l++)
             elemat1_epetra(i*noddof_+j,k*noddof_+l) = elemat1_sub(i*my::noddof_+j,k*my::noddof_+l);

     if(elemat2_epetra.A())
     for(int i=0;i<my::numnod_;i++)
       for(int j=0;j<my::numdim_;j++)
         for(int k=0;k<my::numnod_;k++)
           for(int l=0;l<my::numdim_;l++)
             elemat2_epetra(i*noddof_+j,k*noddof_+l) = elemat2_sub(i*my::noddof_+j,k*my::noddof_+l);

     if(elevec1_epetra.A())
       for(int i=0;i<my::numnod_;i++)
         for(int j=0;j<my::numdim_;j++)
               elevec1_epetra(i*noddof_+j) = elevec1_sub(i*my::noddof_+j);

     if(elevec2_epetra.A())
       for(int i=0;i<my::numnod_;i++)
         for(int j=0;j<my::numdim_;j++)
               elevec2_epetra(i*noddof_+j) = elevec2_sub(i*my::noddof_+j);

     if(elevec3_epetra.A())
       for(int i=0;i<my::numnod_;i++)
         for(int j=0;j<my::numdim_;j++)
               elevec3_epetra(i*noddof_+j) = elevec3_sub(i*my::noddof_+j);

     //add volume coupling specific terms
      MyEvaluate(params,
               discretization,
               la,
               elemat1_epetra,
               elemat2_epetra,
               elevec1_epetra,
               elevec2_epetra,
               elevec3_epetra);
  }
  break;
  } // action

  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (protected)                                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::MyEvaluate(
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
      my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, &myporosity,"displacement");

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
      my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, &myporosity,"displacement");

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
    my::ExtractValuesFromGlobalVector(discretization,0,la[0].lm_, &mydisp, &myporosity,"displacement");

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
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::InitElement()
{
  //initialize base element
  my::InitElement();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::nlnstiff_poroelast(
    std::vector<int>& lm, ///< location matrix
    LINALG::Matrix<my::numdim_, my::numnod_>&               disp,         // current displacements
    LINALG::Matrix<my::numdim_, my::numnod_>&               vel,          // current velocities
    LINALG::Matrix<my::numnod_, 1>*               porosity_dof,
    LINALG::Matrix<my::numdim_,my::numnod_> & evelnp, //< fluid velocity of element
    LINALG::Matrix<my::numnod_,1> & epreaf, //< fluid pressure of element
    LINALG::Matrix<numdof_,numdof_>* stiffmatrix, ///< element stiffness matrix
    LINALG::Matrix<numdof_,numdof_>* reamatrix, // element reactive matrix
    LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
    Teuchos::ParameterList& params ///< algorithmic parameters e.g. time
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

  //initialize element matrizes and vectors
  LINALG::Matrix<my::numdof_,my::numdof_> estiff_stat(true);
  LINALG::Matrix<my::numdof_,my::numdof_> erea_u(true);
  LINALG::Matrix<my::numdof_,my::numdof_> erea_v(true);
  LINALG::Matrix<my::numdof_,my::numdof_> estiff_stress(true);
  LINALG::Matrix<my::numdof_,1> erea_force(true);
  LINALG::Matrix<my::numdof_,1> ecoupl_force_p(true);
  LINALG::Matrix<my::numdof_,1> ecoupl_force_v(true);
  LINALG::Matrix<my::numdof_,1> ecoupl_force_stress(true);
  LINALG::Matrix<my::numstr_,1> fstress(true);

  LINALG::Matrix<my::numdof_,my::numnod_> ecoupl_p1(true);
  LINALG::Matrix<my::numnod_,numdof_> estiff_p1(true);
  LINALG::Matrix<my::numnod_,1> ecoupl_force_p1(true);

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
                        estiff_stat,
                        erea_u,
                        erea_v,
                        estiff_stress,
                        erea_force,
                        ecoupl_force_p,
                        ecoupl_force_v,
                        ecoupl_force_stress,
                        fstress);

  GaussPointLoopP1(  params,
                          xrefe,
                          xcurr,
                          disp,
                          vel,
                          evelnp,
                          epreaf,
                          porosity_dof,
                          ecoupl_p1,
                          estiff_p1,
                          ecoupl_force_p1
                        );

  // update stiffness matrix
  if (stiffmatrix != NULL)
  {
    if ( reamatrix != NULL )
    {
      /* additional "reactive darcy-term"
       detJ * w(gp) * ( J * reacoeff * phi^2  ) * D(v_s)
       */
      for (int k=0; k<my::numnod_; k++)
        for(int l=0; l<my::numdim_; l++)
          for (int i=0; i<my::numnod_; i++)
            for(int j=0; j<my::numdim_; j++)
              (*reamatrix)(i*noddof_+j,k*noddof_+l) += erea_v(i*my::numdim_+j,k*my::numdim_+l);
    }
    else
    {
      const double dt = params.get<double>("delta time");
      //if the reaction part is not supposed to be computed separately, we add it to the stiffness
      //(this is not the best way to do it, but it only happens once during initialization)
      for (int k=0; k<my::numnod_; k++)
        for(int l=0; l<my::numdim_; l++)
          for (int i=0; i<my::numnod_; i++)
            for(int j=0; j<my::numdim_; j++)
              (*stiffmatrix)(i*noddof_+j,k*noddof_+l) += erea_v(i*my::numdim_+j,k*my::numdim_+l) / dt;
    }

    for (int k=0; k<my::numnod_; k++)
    {
      for(int l=0; l<my::numdim_; l++)
        for (int i=0; i<my::numnod_; i++)
        {
          for(int j=0; j<my::numdim_; j++)
            (*stiffmatrix)(i*noddof_+j,k*noddof_+l) +=   estiff_stat(i*my::numdim_+j,k*my::numdim_+l)
                                                       + erea_u(i*my::numdim_+j,k*my::numdim_+l)
                                                       + estiff_stress(i*my::numdim_+j,k*my::numdim_+l);
        }
      for (int i=0; i<my::numnod_; i++)
        for(int j=0; j<my::numdim_; j++)
          (*stiffmatrix)(i*noddof_+j,k*noddof_+my::numdim_) += ecoupl_p1(i*my::noddof_+j,k);
    }

    for (int i=0; i<my::numnod_; i++)
      for(int j=0; j<my::numnod_; j++)
        for (int k=0; k<noddof_; k++)
           (*stiffmatrix)(i*noddof_+my::numdim_,j*noddof_+k) += estiff_p1(i,j*noddof_+k);
  }

  // update internal force vector
  if (force != NULL )
  {
    for (int i=0; i<my::numnod_; i++)
    {
      for(int j=0; j<my::numdim_; j++)
        (*force)(i*noddof_+j) +=    ecoupl_force_stress(i*my::numdim_+j)
                                  + ecoupl_force_p(i*my::numdim_+j)
                                  + ecoupl_force_v(i*my::numdim_+j)
                                  + erea_force(i*my::numdim_+j);

      (*force)(i*noddof_+my::numdim_) += ecoupl_force_p1(i);
    }
  }  // if (force != NULL )

  return;
}  // nlnstiff_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::GaussPointLoopP1(
                                    Teuchos::ParameterList& params,
                                    const LINALG::Matrix<my::numdim_,my::numnod_>& xrefe,
                                    const LINALG::Matrix<my::numdim_,my::numnod_>& xcurr,
                                    const LINALG::Matrix<my::numdim_,my::numnod_>& nodaldisp,
                                    const LINALG::Matrix<my::numdim_,my::numnod_>& nodalvel,
                                    const LINALG::Matrix<my::numdim_,my::numnod_> & evelnp,
                                    const LINALG::Matrix<my::numnod_,1> & epreaf,
                                    const LINALG::Matrix<my::numnod_, 1>*  porosity_dof,
                                    LINALG::Matrix<my::numdof_,my::numnod_>& ecoupl_p1,
                                    LINALG::Matrix<my::numnod_,numdof_>& estiff_p1,
                                    LINALG::Matrix<my::numnod_,1>& ecoupl_force_p1
                                        )
{
  double reacoeff = my::fluidmat_->ComputeReactionCoeff();

  LINALG::Matrix<my::numdim_,my::numnod_> N_XYZ;
  LINALG::Matrix<my::numderiv2_,my::numnod_> N_XYZ2;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<my::numdim_,my::numdim_> defgrd(true);
  LINALG::Matrix<my::numnod_,1> shapefct;
  LINALG::Matrix<my::numdim_,my::numnod_> deriv ;
  LINALG::Matrix<my::numderiv2_,my::numnod_> deriv2;

  for (int gp=0; gp<my::numgpt_; ++gp)
  {
    DRT::UTILS::shape_function<distype>(my::xsi_[gp],shapefct);
    DRT::UTILS::shape_function_deriv1<distype>(my::xsi_[gp],deriv);

    /* get the inverse of the Jacobian matrix which looks like:
     **            [ X_,r  Y_,r  Z_,r ]^-1
     **     J^-1 = [ X_,s  Y_,s  Z_,s ]
     **            [ X_,t  Y_,t  Z_,t ]
     */

    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(my::invJ_[gp],deriv); // (6.21)
    double detJ = my::detJ_[gp]; // (6.22)

    if( my::ishigherorder_ )
    {
      // transposed jacobian "dX/ds"
      LINALG::Matrix<my::numdim_,my::numdim_> xjm0;
      xjm0.MultiplyNT(deriv,xrefe);

      // get the second derivatives of standard element at current GP w.r.t. rst
      DRT::UTILS::shape_function_deriv2<distype>(my::xsi_[gp],deriv2);
      // get the second derivatives of standard element at current GP w.r.t. XYZ
      DRT::UTILS::gder2<distype>(xjm0,N_XYZ,deriv2,xrefe,N_XYZ2);
    }
    else
    {
      deriv2.Clear();
      N_XYZ2.Clear();
    }

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<my::numdim_,my::numdim_> xjm;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<my::numdim_,my::numdim_> xji;
    xjm.MultiplyNT(deriv,xcurr);
    const double det = xji.Invert(xjm);

    // determinant of deformationgradient: det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
    const double J = det/detJ;

    //----------------------------------------------------
    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    // pressure gradient at integration point
    LINALG::Matrix<my::numdim_,1> Gradp;
    Gradp.Multiply(N_XYZ,epreaf);

    // fluid velocity at integration point
    LINALG::Matrix<my::numdim_,1> fvelint;
    fvelint.Multiply(evelnp,shapefct);

    // material fluid velocity gradient at integration point
    LINALG::Matrix<my::numdim_,my::numdim_>              fvelder;
    fvelder.MultiplyNT(evelnp,N_XYZ);

    // structure displacement and velocity at integration point
    LINALG::Matrix<my::numdim_,1> dispint(true);
    LINALG::Matrix<my::numdim_,1> velint(true);

    for(int i=0; i<my::numnod_; i++)
      for(int j=0; j<my::numdim_; j++)
      {
        dispint(j) += nodaldisp(j,i) * shapefct(i);
        velint(j) += nodalvel(j,i) * shapefct(i);
      }

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
    defgrd.MultiplyNT(xcurr,N_XYZ); //  (6.17)

    // inverse deformation gradient F^-1
    LINALG::Matrix<my::numdim_,my::numdim_> defgrd_inv(false);
    defgrd_inv.Invert(defgrd);

    //------------------------------------ build F^-1 as vector 9x1
    LINALG::Matrix<my::numdim_*my::numdim_,1> defgrd_inv_vec;
    defgrd_inv_vec(0)=defgrd_inv(0,0);
    defgrd_inv_vec(1)=defgrd_inv(0,1);
    defgrd_inv_vec(2)=defgrd_inv(0,2);
    defgrd_inv_vec(3)=defgrd_inv(1,0);
    defgrd_inv_vec(4)=defgrd_inv(1,1);
    defgrd_inv_vec(5)=defgrd_inv(1,2);
    defgrd_inv_vec(6)=defgrd_inv(2,0);
    defgrd_inv_vec(7)=defgrd_inv(2,1);
    defgrd_inv_vec(8)=defgrd_inv(2,2);

//    //------------------------------------ build F^-T as vector 9x1
//    LINALG::Matrix<9,1> defgrd_IT_vec;
//    defgrd_IT_vec(0)=defgrd_inv(0,0);
//    defgrd_IT_vec(1)=defgrd_inv(1,0);
//    defgrd_IT_vec(2)=defgrd_inv(2,0);
//    defgrd_IT_vec(3)=defgrd_inv(0,1);
//    defgrd_IT_vec(4)=defgrd_inv(1,1);
//    defgrd_IT_vec(5)=defgrd_inv(2,1);
//    defgrd_IT_vec(6)=defgrd_inv(0,2);
//    defgrd_IT_vec(7)=defgrd_inv(1,2);
//    defgrd_IT_vec(8)=defgrd_inv(2,2);

    //--------------------------- build N_X operator (wrt material config)
    LINALG::Matrix<9,my::numdof_> N_X(true); // set to zero
    for (int i=0; i<my::numnod_; ++i)
    {
      N_X(0,3*i+0) = N_XYZ(0,i);
      N_X(1,3*i+1) = N_XYZ(0,i);
      N_X(2,3*i+2) = N_XYZ(0,i);

      N_X(3,3*i+0) = N_XYZ(1,i);
      N_X(4,3*i+1) = N_XYZ(1,i);
      N_X(5,3*i+2) = N_XYZ(1,i);

      N_X(6,3*i+0) = N_XYZ(2,i);
      N_X(7,3*i+1) = N_XYZ(2,i);
      N_X(8,3*i+2) = N_XYZ(2,i);
    }

    //------linearization of jacobi determinant detF=J w.r.t. strucuture displacement   dJ/d(us) = dJ/dF : dF/dus = J * F^-T * N,X
    LINALG::Matrix<1,my::numdof_> dJ_dus;
    dJ_dus.MultiplyTN(J,defgrd_inv_vec,N_X);

    //F^-T * Grad p
    LINALG::Matrix<my::numdim_,1> Finvgradp;
    Finvgradp.MultiplyTN(defgrd_inv, Gradp);

    //--------------------------------------------------------------------

    //linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
    LINALG::Matrix<1,my::numdof_> dphi_dus(true);
    double porosity=0.0;

    ComputePorosityAndLinearization(params,press,J,gp,shapefct,porosity_dof,dJ_dus,porosity,dphi_dus);

    double    dW_dphi = 0.0;
    double    dW_dJ   = 0.0;
    double    dW_dp   = 0.0;
    double    W       = 0.0;
    my::structmat_->ConsitutiveDerivatives(params,
                                      press,
                                      J,
                                      porosity,
                                      &dW_dp, //dW_dp not needed
                                      &dW_dphi,
                                      &dW_dJ,
                                      &W);
    //--------------------------------------------------------

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    double detJ_w = detJ*my::intpoints_.Weight(gp);//gpweights[gp];

    //if (force != NULL or stiffmatrix != NULL or reamatrix != NULL )
    {
      for (int k=0; k<my::numnod_; k++)
      {
        //const int fk = my::numdim_*k;
        const double fac = detJ_w* shapefct(k);
        //const double v = fac * reacoeff * porosity * porosity* J;

        ecoupl_force_p1(k) += fac*W;//fac*(dW_dp*press+dW_dphi*porosity+dW_dJ*J);//fac*W;

        for (int i=0; i<my::numnod_; i++)
        {
          for(int j=0; j<my::numdim_; j++)
          {
            estiff_p1(k,i*noddof_+j) += fac * dW_dJ * dJ_dus(i*my::numdim_+j);

            ecoupl_p1(i*my::numdim_+j,k) += fac * ( 2 * J * reacoeff * porosity * ( velint(j)-fvelint(j) )
                                                  + J * Finvgradp(j) )
                                              * shapefct(i);
          }
          estiff_p1(k,i*noddof_+my::numdim_) += fac * dW_dphi * shapefct(i);
        }
      }

    }
    /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::coupling_poroelast(
    std::vector<int>& lm,                                                 // location matrix
    LINALG::Matrix<my::numdim_, my::numnod_>&               disp,         // current displacements
    LINALG::Matrix<my::numdim_, my::numnod_>&               vel,          // current velocities
    LINALG::Matrix<my::numnod_, 1>*               porosity,
    LINALG::Matrix<my::numdim_, my::numnod_> & evelnp,                       //current fluid velocity
    LINALG::Matrix<my::numnod_, 1> & epreaf,                             //current fluid pressure
    LINALG::Matrix<numdof_, (my::numdim_ + 1) * my::numnod_>* stiffmatrix,   // element stiffness matrix
    LINALG::Matrix<numdof_, (my::numdim_ + 1) * my::numnod_>* reamatrix,     // element reactive matrix
    LINALG::Matrix<numdof_, 1>* force,                               // element internal force vector
    Teuchos::ParameterList& params)                                           // algorithmic parameters e.g. time
{
  //=============================get parameters

  my::GetMaterials();

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
  //initialize element matrizes
  LINALG::Matrix<my::numdof_,(my::numdim_+1)*my::numnod_> ecoupl(true);
  LINALG::Matrix<my::numdof_,my::numnod_> ecoupl_p(true);
  LINALG::Matrix<my::numdof_,my::numdof_> ecoupl_v(true);

  LINALG::Matrix<my::numnod_,my::numnod_> ecoupl_p1_p(true);

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
                         ecoupl_p,
                         ecoupl_v  );

  GaussPointLoopP1OD( params,
                      xrefe,
                      xcurr,
                      disp,
                      vel,
                      evelnp,
                      epreaf,
                      porosity,
                      ecoupl_p1_p);

  //if (force != NULL )
  //{
    //all rhs terms are added in nlnstiff_poroelast
  //}

  if (stiffmatrix != NULL )
  {

    // add structure displacement - fluid velocity part to matrix
    for (int ui=0; ui<my::numnod_; ++ui)
    {
      const int dim_ui = my::numdim_*ui;

      for (int jdim=0; jdim < my::numdim_;++jdim)
      {
        const int dim_ui_jdim = dim_ui+jdim;

        for (int vi=0; vi<my::numnod_; ++vi)
        {
          const int numdof_vi = (my::numdim_+1)*vi;
          const int dim_vi = my::numdim_*vi;

          for (int idim=0; idim <my::numdim_; ++idim)
            ecoupl(dim_ui_jdim , numdof_vi+idim) += ecoupl_v(dim_ui_jdim , dim_vi+idim);
        }
      }
    }

    // add structure displacement - fluid pressure part to matrix
    for (int ui=0; ui<my::numnod_; ++ui)
    {
      const int dim_ui = my::numdim_*ui;

      for (int jdim=0; jdim < my::numdim_;++jdim)
      {
        const int dim_ui_jdim = dim_ui+jdim;

        for (int vi=0; vi<my::numnod_; ++vi)
          ecoupl( dim_ui_jdim , (my::numdim_+1)*vi+my::numdim_ ) += ecoupl_p( dim_ui_jdim , vi);
      }
    }

    //TODO theta should be on time integrator level
    // build tangent coupling matrix : effective dynamic stiffness coupling matrix
    //    K_{Teffdyn} = 1/dt C
    //                + theta K_{T}
    const double theta = params.get<double>("theta");

    for (int k=0; k<my::numnod_; k++)
      for(int l=0; l<(my::numdim_+1); l++)
        for (int i=0; i<my::numnod_; i++)
          for(int j=0; j<my::numdim_; j++)
            (*stiffmatrix)(i*noddof_+j,k*(my::numdim_+1)+l) += theta * ecoupl(i*my::numdim_+j,k*(my::numdim_+1)+l);

    for (int ui=0; ui<my::numnod_; ++ui)
      for (int ni=0; ni<my::numnod_; ++ni)
        (*stiffmatrix)( noddof_*ui+my::numdim_ , (my::numdim_+1)*ni+my::numdim_ ) += theta * ecoupl_p1_p( ui , ni);
  }

  //cout<<"ecoupl_p1_p"<<endl<<ecoupl_p1_p<<endl;
  return;

}  // coupling_poroelast()

/*----------------------------------------------------------------------*
 |  evaluate only the poroelasticity fraction for the element (protected) |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::GaussPointLoopP1OD(
                                    Teuchos::ParameterList& params,
                                    const LINALG::Matrix<my::numdim_,my::numnod_>& xrefe,
                                    const LINALG::Matrix<my::numdim_,my::numnod_>& xcurr,
                                    const LINALG::Matrix<my::numdim_,my::numnod_>& nodaldisp,
                                    const LINALG::Matrix<my::numdim_,my::numnod_>& nodalvel,
                                    const LINALG::Matrix<my::numdim_,my::numnod_> & evelnp,
                                    const LINALG::Matrix<my::numnod_,1> & epreaf,
                                    const LINALG::Matrix<my::numnod_, 1>*  porosity_dof,
                                    LINALG::Matrix<my::numnod_,my::numnod_>& ecoupl_p1
                                        )
{

  LINALG::Matrix<my::numnod_,1> shapefct;
  LINALG::Matrix<my::numdim_,my::numnod_> deriv ;

  for (int gp=0; gp<my::numgpt_; ++gp)
  {
    DRT::UTILS::shape_function<distype>(my::xsi_[gp],shapefct);
    DRT::UTILS::shape_function_deriv1<distype>(my::xsi_[gp],deriv);

    /* get the inverse of the Jacobian matrix which looks like:
     **            [ X_,r  Y_,r  Z_,r ]^-1
     **     J^-1 = [ X_,s  Y_,s  Z_,s ]
     **            [ X_,t  Y_,t  Z_,t ]
     */

    double detJ = my::detJ_[gp]; // (6.22)

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    //! transposed jacobian "dx/ds"
    LINALG::Matrix<my::numdim_,my::numdim_> xjm;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<my::numdim_,my::numdim_> xji;
    xjm.MultiplyNT(deriv,xcurr);
    const double det = xji.Invert(xjm);

    // determinant of deformationgradient: det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
    const double J = det/detJ;

    //----------------------------------------------------
    // pressure at integration point
    double press = shapefct.Dot(epreaf);

    //--------------------------------------------------------------------

    //linearization of porosity w.r.t structure displacement d\phi/d(us) = d\phi/dJ*dJ/d(us)
    LINALG::Matrix<1,my::numdof_> dphi_dus(true);
    double porosity=0.0;
    //dummy
    LINALG::Matrix<1,my::numdof_>    dJ_dus(true);

    ComputePorosityAndLinearization(params,press,J,gp,shapefct,porosity_dof,dJ_dus,porosity,dphi_dus);

    double    dW_dp   = 0.0;
    my::structmat_->ConsitutiveDerivatives(params,
                                      press,
                                      J,
                                      porosity,
                                      &dW_dp,
                                      NULL, // not needed
                                      NULL, // not needed
                                      NULL  // not needed
                                      );
    //--------------------------------------------------------

    // **********************evaluate stiffness matrix and force vector+++++++++++++++++++++++++
    double detJ_w = detJ*my::intpoints_.Weight(gp);//gpweights[gp];

    for (int k=0; k<my::numnod_; k++)
    {
      const double fac = detJ_w* shapefct(k);

      for (int i=0; i<my::numnod_; i++)
        ecoupl_p1(k,i) += fac * dW_dp * shapefct(i);
    }

    /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
  /* =========================================================================*/
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
#include "so3_poro_p1_fwd.hpp"
