/*!----------------------------------------------------------------------
\file so_sh8_evaluate.cpp
\brief

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#include "so_sh18.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mat/so3_material.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                           seitz 11/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh18::Evaluate(Teuchos::ParameterList&  params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH18,1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH18,1> elevec2(elevec2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH18,1> elevec3(elevec3_epetra.A(),true);

  // start with "none"
  DRT::ELEMENTS::So_sh18::ActionType act = So_sh18::none;

  // get the required action
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")                        act = So_sh18::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")                        act = So_sh18::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")                   act = So_sh18::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")                    act = So_sh18::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")                    act = So_sh18::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")                   act = So_sh18::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")                          act = So_sh18::calc_struct_stress;
  else if (action=="calc_struct_eleload")                         act = So_sh18::calc_struct_eleload;
  else if (action=="calc_struct_update_istep")                    act = So_sh18::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")                     act = So_sh18::calc_struct_reset_istep;
  else if (action=="calc_struct_reset_all")                       act = So_sh18::calc_struct_reset_all;
  else if (action=="calc_struct_energy")                          act = So_sh18::calc_struct_energy;
  else if (action=="calc_struct_errornorms")                      act = So_sh18::calc_struct_errornorms;
  else if (action=="postprocess_stress")                          act = So_sh18::postprocess_stress;
  else dserror("Unknown type of action for So_hex8: %s", action.c_str());

  // what should the element do
  switch(act)
  {

    //==================================================================================
    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

        nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
                                    INPAR::STR::stress_none,INPAR::STR::strain_none);
      break;
    }


    //==================================================================================
    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18> myemat(true);

      nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
        INPAR::STR::stress_none,INPAR::STR::strain_none);

      break;
    }

    //==================================================================================
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    case calc_struct_linstiffmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      // need current velocities and accelerations (for non constant mass matrix)
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");

      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

        nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
          INPAR::STR::stress_none,INPAR::STR::strain_none);

      if (act==calc_struct_nlnstifflmass) soh18_lumpmass(&elemat2);

      break;
    }

    //==================================================================================
    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID()==Owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
        Teuchos::RCP<std::vector<char> > stressdata = params.get<Teuchos::RCP<std::vector<char> > >("stress",Teuchos::null);
        Teuchos::RCP<std::vector<char> > straindata = params.get<Teuchos::RCP<std::vector<char> > >("strain",Teuchos::null);
        if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata==Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata==Teuchos::null) dserror("Cannot get 'strain' data");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res,myres,lm);
        LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D> stress;
        LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D> strain;
        INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

          nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,iostress,iostrain);

        {
          DRT::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(),data().end(),std::back_inserter(*stressdata));
        }

        {
          DRT::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(),data().end(),std::back_inserter(*straindata));
        }
      }
    }
    break;

    //==================================================================================
    // postprocess stresses/strains at gauss points
    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
      if (gpstressmap==Teuchos::null)
        dserror("no gp stress/strain map available for postprocessing");
      std::string stresstype = params.get<std::string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D> gpstress(((*gpstressmap)[gid])->A(),true);
      Teuchos::RCP<Epetra_MultiVector> poststress=params.get<Teuchos::RCP<Epetra_MultiVector> >("poststress",Teuchos::null);
      if (poststress==Teuchos::null)
        dserror("No element stress/strain vector available");
      if (stresstype=="ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        dserror("no node-based stress output");
      }
      else if (stresstype=="cxyz")
      {
        const Epetra_BlockMap& elemap = poststress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1)
        {
          for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
          {
            double& s = (*((*poststress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_SOH18; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOH18;
          }
        }
      }
      else
      {
        dserror("unknown type of stress/strain output on element level");
      }
    }
    break;

    //==================================================================================
    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    //==================================================================================
    case calc_struct_update_istep:
    {
      Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
      so3mat->Update();

      if (eas_)
      {
        feas_.Clear();
        Kad_.Clear();
        KaaInv_.Clear();
      }
    }
    break;

    //==================================================================================
    case calc_struct_reset_istep:
    {
      // Reset of history (if needed)
      Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
      so3mat->ResetStep();
    }
    break;

    //==================================================================================
    case calc_struct_reset_all:
    {
      // Reset of history for materials
      Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
      so3mat->ResetAll(NUMGPT_SOH18);
    }
    break;

  //==================================================================================
  default:
    dserror("Unknown type of action for So_hex18");
    break;
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)  seitz 11/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh18::EvaluateNeumann(Teuchos::ParameterList&   params,
                                             DRT::Discretization&      discretization,
                                             DRT::Condition&           condition,
                                             std::vector<int>&         lm,
                                             Epetra_SerialDenseVector& elevec1,
                                             Epetra_SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve = condition.Get<std::vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // (SPATIAL) FUNCTION BUSINESS
  const std::vector<int>* funct = condition.Get<std::vector<int> >("funct");
  LINALG::Matrix<NUMDIM_SOH18,1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim=0; dim<NUMDIM_SOH18; dim++)
      if ((*funct)[dim] > 0)
        havefunct = havefunct or true;

  /* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH18; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOH18; ++gp) {

    // shape function and derivatives
    const LINALG::Matrix<NUMNOD_SOH18,1> shapefunct = sh18_shapefcts(gp);
    const LINALG::Matrix<NUMDIM_SOH18,NUMNOD_SOH18> deriv = sh18_derivs(gp);

    // compute the Jacobian matrix
    LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18> jac;
    jac.Multiply(deriv,xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct) {
      for (int dim=0; dim<NUMDIM_SOH18; dim++) {
        xrefegp(dim) = 0.0;
        for (int nodid=0; nodid<NUMNOD_SOH18; ++nodid)
          xrefegp(dim) += shapefunct(nodid) * xrefe(nodid,dim);
      }
    }

    // integration factor
    const double fac = wgt_[gp] * curvefac * detJ;
    // distribute/add over element load vector
    for(int dim=0; dim<NUMDIM_SOH18; dim++) {
      // function evaluation
      const int functnum = (funct) ? (*funct)[dim] : -1;
      const double functfac
        = (functnum>0)
        ? DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,xrefegp.A(),time,NULL)
        : 1.0;
      const double dim_fac = (*onoff)[dim] * (*val)[dim] * fac * functfac;
      for (int nodid=0; nodid<NUMNOD_SOH18; ++nodid) {
        elevec1[nodid*NUMDIM_SOH18+dim] += shapefunct(nodid) * dim_fac;
      }
    }

  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::So_hex8::EvaluateNeumann


int DRT::ELEMENTS::So_sh18::InitJacobianMapping()
{

  LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18> xrefe;
  for (int i=0; i<NUMNOD_SOH18; ++i)
  {
    Node** nodes=Nodes();
    if(!nodes) dserror("Nodes() returned null pointer");
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
//  std::cout << "ele: " << Id() << " xrefe: " << xrefe ;
  invJ_.resize(NUMGPT_SOH18);
  detJ_.resize(NUMGPT_SOH18);


  for (int gp=0; gp<NUMGPT_SOH18; ++gp)
  {
    // reset
    invJ_[gp].Clear();
    detJ_[gp]=0.;

    // in-plane shape functions and derivatives
    const LINALG::Matrix<9,1> shapefunct_q9 = sh18_shapefcts_q9(gp);
    const LINALG::Matrix<2,9> deriv_q9      = sh18_derivs_q9(gp);

    for (int dim=0; dim<3; ++dim)
      for (int k=0; k<9; ++k)
      {
        invJ_[gp](0,dim) += .5            *deriv_q9(0,k)*(xrefe(k+9,dim)+xrefe(k,dim))
                           +.5*xsi_[gp](2)*deriv_q9(0,k)*(xrefe(k+9,dim)-xrefe(k,dim));

        invJ_[gp](1,dim) += .5            *deriv_q9(1,k)*(xrefe(k+9,dim)+xrefe(k,dim))
                           +.5*xsi_[gp](2)*deriv_q9(1,k)*(xrefe(k+9,dim)-xrefe(k,dim));

        invJ_[gp](2,dim) += .5*shapefunct_q9(k)*(xrefe(k+9,dim)-xrefe(k,dim));
      }
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp]<0.) return 1;
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::nlnstiffmass(
    std::vector<int>& lm,  ///< location matrix
    std::vector<double>& disp,  ///< current displacements
    std::vector<double>& residual,  ///< current residual displ
    LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18>* stiffmatrix,  ///< element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18>* massmatrix,  ///< element mass matrix
    LINALG::Matrix<NUMDOF_SOH18,1>* force,  ///< element internal force vector
    LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
    LINALG::Matrix<NUMGPT_SOH18,MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
    Teuchos::ParameterList& params,  ///< algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,  ///< stress output option
    const INPAR::STR::StrainType iostrain   ///< strain output option
    )
{
  // update element geometry
  LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18> xrefe;  // reference coord. of element
  LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18> xcurr;  // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH18; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH18+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH18+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH18+2];
  }

  // we need the (residual) displacement at the previous step
  LINALG::Matrix<NUMDOF_SOH18,1> res_d;
  for (int i=0; i<NUMDOF_SOH18; ++i)
    res_d(i) = residual[i];

  // EAS stuff
  std::vector<LINALG::Matrix<6,num_eas> > M_gp(num_eas);
  LINALG::Matrix<3,1> G3_0_contra;
  LINALG::Matrix<6,num_eas> M;
  if (eas_)
  {
  // recover EAS **************************************
    if (stiffmatrix)
    {
      feas_.Multiply(1.,Kad_,res_d,1.);
      alpha_eas_.Multiply(-1.,KaaInv_,feas_,1.);
    }
    // recover EAS **************************************

    // prepare EAS***************************************
    EasSetup(M_gp,G3_0_contra,xrefe);
    feas_.Clear();
    KaaInv_.Clear();
    Kad_.Clear();
    // prepare EAS***************************************
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH18; ++gp)
  {
    // in-plane shape functions and derivatives
    const LINALG::Matrix<9,1> shapefunct_q9 = sh18_shapefcts_q9(gp);
    const LINALG::Matrix<2,9> deriv_q9      = sh18_derivs_q9(gp);

    /* get the inverse of the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    // compute the Jacobian shell-style (G^T)
    LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18> jac;
    for (int dim=0; dim<3; ++dim)
      for (int k=0; k<9; ++k)
      {
        jac(0,dim) += .5            *deriv_q9(0,k)*(xrefe(k+9,dim)+xrefe(k,dim))
                     +.5*xsi_[gp](2)*deriv_q9(0,k)*(xrefe(k+9,dim)-xrefe(k,dim));

        jac(1,dim) += .5            *deriv_q9(1,k)*(xrefe(k+9,dim)+xrefe(k,dim))
                     +.5*xsi_[gp](2)*deriv_q9(1,k)*(xrefe(k+9,dim)-xrefe(k,dim));

        jac(2,dim) += .5*shapefunct_q9(k)*(xrefe(k+9,dim)-xrefe(k,dim));
      }
    double detJ = jac.Determinant();

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D> TinvT;
    EvaluateT(jac,TinvT);

    // **********************************************************************
    // set up B-Operator in local(parameter) element space including ANS
    // **********************************************************************
    LINALG::Matrix<MAT::NUM_STRESS_3D,NUMDOF_SOH18> bop_loc(true);
    CalculateBopLoc(xcurr,xrefe,shapefunct_q9,deriv_q9,gp,bop_loc);
    LINALG::Matrix<MAT::NUM_STRESS_3D,NUMDOF_SOH18> bop;
    bop.Multiply(TinvT,bop_loc);

    // **************************************************************************
    // shell-like calculation of strains
    // see Diss. Koschnik page 41
    // **************************************************************************
    LINALG::Matrix<MAT::NUM_STRESS_3D,1> lstrain(true);
    CalculateLocStrain(xcurr,xrefe,shapefunct_q9,deriv_q9,gp,lstrain);
    LINALG::Matrix<MAT::NUM_STRESS_3D,1> glstrain;
    glstrain.Multiply(TinvT,lstrain);
    // **************************************************************************
    // shell-like calculation of strains
    // **************************************************************************

      // EAS: enhance the strains ***********************************************
      if(eas_)
      {
        double t33=0.;
        for (int dim=0; dim<3; ++dim)
          t33 += jac(2,dim)*G3_0_contra(dim);

        M.Multiply(t33*t33/detJ,TinvT,M_gp[gp],0.);
        glstrain.Multiply(1.,M,alpha_eas_,1.);
      }
      // end EAS: enhance the strains *******************************************

      // calculate the deformation gradient consistent to the modified strains
      // but only if the material needs a deformation gradient (e.g. plasticity)
      LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18> defgrd;
      if (Teuchos::rcp_static_cast<MAT::So3Material>(Material())->NeedsDefgrd())
      {
        // compute the deformation gradient - shell-style
        // deformation gradient with derivatives w.r.t. local basis
        LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18> defgrd_loc(true);
        for (int k=0; k<9; ++k)
          for (int dim=0; dim<NUMDIM_SOH18; ++dim)
          {
            defgrd_loc(dim,0) += .5*deriv_q9(0,k)
                                  *(             (xcurr(k+9,dim)+xcurr(k,dim))
                                    +xsi_[gp](2)*(xcurr(k+9,dim)-xcurr(k,dim)));
            defgrd_loc(dim,1) += .5*deriv_q9(1,k)
                                  *(             (xcurr(k+9,dim)+xcurr(k,dim))
                                    +xsi_[gp](2)*(xcurr(k+9,dim)-xcurr(k,dim)));
            defgrd_loc(dim,2) += .5*shapefunct_q9(k)*(xcurr(k+9,dim)-xcurr(k,dim));
          }

        // displacement-based deformation gradient
        LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18> defgrd_disp;
        defgrd_disp.MultiplyNT(defgrd_loc,invJ_[gp]);
        if (eas_ || dsg_shear_ || dsg_membrane_ || dsg_ctl_)
          CalcConsistentDefgrd(defgrd_disp,glstrain,defgrd);
      }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix must be retrieved,
    ** all necessary data must be passed.
    */
    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D> cmat(true);
    LINALG::Matrix<MAT::NUM_STRESS_3D,1> stress(true);
    params.set<int>("gp",gp);
    Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_static_cast<MAT::So3Material>(Material());
    so3mat->Evaluate(&defgrd,&glstrain,params,&stress,&cmat,Id());
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    double detJ_w = detJ*wgt_[gp];
    // update internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }

    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<MAT::NUM_STRESS_3D,NUMDOF_SOH18> cb;
      cb.Multiply(cmat,bop);
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);  // standard hex8 evaluation
      // intergrate `geometric' stiffness matrix and add to keu *****************
      CalculateGeoStiff(shapefunct_q9,deriv_q9,TinvT,gp,detJ_w,stress,stiffmatrix);

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eas_)
      {
        LINALG::Matrix<6,num_eas> cM;
        cM.Multiply(cmat,M);
        KaaInv_.MultiplyTN(detJ_w,M,cM,1.);
        Kad_.MultiplyTN(detJ_w,M,cb,1.);
        feas_.MultiplyTN(detJ_w,M,stress,1.);
      }
      // EAS technology: integrate matrices --------------------------------- EAS
    }

    if (massmatrix != NULL) // evaluate mass matrix +++++++++++++++++++++++++
    {
      // shape function and derivatives
      const LINALG::Matrix<NUMNOD_SOH18,1> shapefunct = sh18_shapefcts(gp);

      double density = Material()->Density(gp);

      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOH18; ++inod)
      {
        ifactor = shapefunct(inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOH18; ++jnod)
        {
          massfactor = shapefunct(jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_SOH18*inod+0,NUMDIM_SOH18*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH18*inod+1,NUMDIM_SOH18*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH18*inod+2,NUMDIM_SOH18*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  if (stiffmatrix && eas_)
  {
    LINALG::FixedSizeSerialDenseSolver<num_eas,num_eas,1> solve_for_KaaInv;
    solve_for_KaaInv.SetMatrix(KaaInv_);
    int err2 = solve_for_KaaInv.Factor();
    int err = solve_for_KaaInv.Invert();
    if ((err != 0) || (err2!=0))
      dserror("Inversion of Kaa failed");

    LINALG::Matrix<NUMDOF_SOH18,num_eas> KdaKaa;
    KdaKaa.MultiplyTN(Kad_,KaaInv_);
    stiffmatrix->Multiply(-1.,KdaKaa,Kad_,1.);
    force->Multiply(-1.,KdaKaa,feas_,1.);
  }

  return;
} // DRT::ELEMENTS::So_hex8::nlnstiffmass


/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                              seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::soh18_lumpmass(LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r=0; r<(*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r,c);  // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
}


const LINALG::Matrix<NUMNOD_SOH18,1> DRT::ELEMENTS::So_sh18::sh18_shapefcts(const int gp)
{ return sh18_shapefcts(xsi_.at(gp)(0),xsi_.at(gp)(1),xsi_.at(gp)(2)); }

const LINALG::Matrix<9,1> DRT::ELEMENTS::So_sh18::sh18_shapefcts_q9(const int gp)
{ return sh18_shapefcts_q9(xsi_.at(gp)(0),xsi_.at(gp)(1)); }

const LINALG::Matrix<NUMNOD_SOH18,1> DRT::ELEMENTS::So_sh18::sh18_shapefcts(
  const double r,
  const double s,
  const double t)
{
  LINALG::Matrix<NUMNOD_SOH18,1> shapefcts;

  // fill up nodal f

  const double rp = 1.0 + r;
  const double rm = 1.0 - r;
  const double sp = 1.0 + s;
  const double sm = 1.0 - s;
  const double r2 = 1.0 - r * r;
  const double s2 = 1.0 - s * s;
  const double rh = 0.5 * r;
  const double sh = 0.5 * s;
  const double rs = rh * sh;

  shapefcts(0) =  rs * rm * sm *(1.0-t)*0.5;
  shapefcts(1) = -rs * rp * sm *(1.0-t)*0.5;
  shapefcts(2) =  rs * rp * sp *(1.0-t)*0.5;
  shapefcts(3) = -rs * rm * sp *(1.0-t)*0.5;
  shapefcts(4) = -sh * sm * r2 *(1.0-t)*0.5;
  shapefcts(5) =  rh * rp * s2 *(1.0-t)*0.5;
  shapefcts(6) =  sh * sp * r2 *(1.0-t)*0.5;
  shapefcts(7) = -rh * rm * s2 *(1.0-t)*0.5;
  shapefcts(8) =  r2 * s2      *(1.0-t)*0.5;

  shapefcts(9)  =  rs * rm * sm *(1.0+t)*0.5;
  shapefcts(10) = -rs * rp * sm *(1.0+t)*0.5;
  shapefcts(11) =  rs * rp * sp *(1.0+t)*0.5;
  shapefcts(12) = -rs * rm * sp *(1.0+t)*0.5;
  shapefcts(13) = -sh * sm * r2 *(1.0+t)*0.5;
  shapefcts(14) =  rh * rp * s2 *(1.0+t)*0.5;
  shapefcts(15) =  sh * sp * r2 *(1.0+t)*0.5;
  shapefcts(16) = -rh * rm * s2 *(1.0+t)*0.5;
  shapefcts(17) =  r2 * s2      *(1.0+t)*0.5;

  return shapefcts;
}

const LINALG::Matrix<9,1> DRT::ELEMENTS::So_sh18::sh18_shapefcts_q9(
  const double r,
  const double s)
{
  LINALG::Matrix<9,1> shapefcts;

  // fill up nodal f

  const double rp = 1.0 + r;
  const double rm = 1.0 - r;
  const double sp = 1.0 + s;
  const double sm = 1.0 - s;
  const double r2 = 1.0 - r * r;
  const double s2 = 1.0 - s * s;
  const double rh = 0.5 * r;
  const double sh = 0.5 * s;
  const double rs = rh * sh;

  shapefcts(0) =  rs * rm * sm;
  shapefcts(1) = -rs * rp * sm;
  shapefcts(2) =  rs * rp * sp;
  shapefcts(3) = -rs * rm * sp;
  shapefcts(4) = -sh * sm * r2;
  shapefcts(5) =  rh * rp * s2;
  shapefcts(6) =  sh * sp * r2;
  shapefcts(7) = -rh * rm * s2;
  shapefcts(8) =  r2 * s2     ;

  return shapefcts;
}


const LINALG::Matrix<NUMDIM_SOH18,NUMNOD_SOH18> DRT::ELEMENTS::So_sh18::sh18_derivs(const int gp)
{
  return sh18_derivs(xsi_.at(gp)(0), xsi_.at(gp)(1), xsi_.at(gp)(2));
}

const LINALG::Matrix<2,9> DRT::ELEMENTS::So_sh18::sh18_derivs_q9(const int gp)
{ return sh18_derivs_q9(xsi_.at(gp)(0), xsi_.at(gp)(1));}

const LINALG::Matrix<NUMDIM_SOH18,NUMNOD_SOH18> DRT::ELEMENTS::So_sh18::sh18_derivs(
    const double r,
    const double s,
    const double t)
{
  LINALG::Matrix<NUMDIM_SOH18,NUMNOD_SOH18> derivs;

  const double rp = 1.0 + r;
  const double rm = 1.0 - r;
  const double sp = 1.0 + s;
  const double sm = 1.0 - s;
  const double r2 = 1.0 - r * r;
  const double s2 = 1.0 - s * s;
  const double rh = 0.5 * r;
  const double sh = 0.5 * s;
  const double rs = rh * sh;
  const double rhp = r + 0.5;
  const double rhm = r - 0.5;
  const double shp = s + 0.5;
  const double shm = s - 0.5;

  derivs(0, 0) = -rhm * sh * sm*(1.0-t)*0.5;
  derivs(1, 0) = -shm * rh * rm*(1.0-t)*0.5;
  derivs(0, 1) = -rhp * sh * sm*(1.0-t)*0.5;
  derivs(1, 1) =  shm * rh * rp*(1.0-t)*0.5;
  derivs(0, 2) =  rhp * sh * sp*(1.0-t)*0.5;
  derivs(1, 2) =  shp * rh * rp*(1.0-t)*0.5;
  derivs(0, 3) =  rhm * sh * sp*(1.0-t)*0.5;
  derivs(1, 3) = -shp * rh * rm*(1.0-t)*0.5;
  derivs(0, 4) =  2.0 * r * sh * sm*(1.0-t)*0.5;
  derivs(1, 4) =  shm * r2*(1.0-t)*0.5;
  derivs(0, 5) =  rhp * s2*(1.0-t)*0.5;
  derivs(1, 5) = -2.0 * s * rh * rp*(1.0-t)*0.5;
  derivs(0, 6) = -2.0 * r * sh * sp*(1.0-t)*0.5;
  derivs(1, 6) =  shp * r2*(1.0-t)*0.5;
  derivs(0, 7) =  rhm * s2*(1.0-t)*0.5;
  derivs(1, 7) =  2.0 * s * rh * rm*(1.0-t)*0.5;
  derivs(0, 8) = -2.0 * r * s2*(1.0-t)*0.5;
  derivs(1, 8) = -2.0 * s * r2*(1.0-t)*0.5;

  derivs(0, 9) = -rhm * sh * sm*(1.0+t)*0.5;
  derivs(1, 9) = -shm * rh * rm*(1.0+t)*0.5;
  derivs(0, 10) = -rhp * sh * sm*(1.0+t)*0.5;
  derivs(1, 10) =  shm * rh * rp*(1.0+t)*0.5;
  derivs(0, 11) =  rhp * sh * sp*(1.0+t)*0.5;
  derivs(1, 11) =  shp * rh * rp*(1.0+t)*0.5;
  derivs(0, 12) =  rhm * sh * sp*(1.0+t)*0.5;
  derivs(1, 12) = -shp * rh * rm*(1.0+t)*0.5;
  derivs(0, 13) =  2.0 * r * sh * sm*(1.0+t)*0.5;
  derivs(1, 13) =  shm * r2*(1.0+t)*0.5;
  derivs(0, 14) =  rhp * s2*(1.0+t)*0.5;
  derivs(1, 14) = -2.0 * s * rh * rp*(1.0+t)*0.5;
  derivs(0, 15) = -2.0 * r * sh * sp*(1.0+t)*0.5;
  derivs(1, 15) =  shp * r2*(1.0+t)*0.5;
  derivs(0, 16) =  rhm * s2*(1.0+t)*0.5;
  derivs(1, 16) =  2.0 * s * rh * rm*(1.0+t)*0.5;
  derivs(0, 17) = -2.0 * r * s2*(1.0+t)*0.5;
  derivs(1, 17) = -2.0 * s * r2*(1.0+t)*0.5;

  derivs(2,0) =  rs * rm * sm *(-0.5);
  derivs(2,1) = -rs * rp * sm *(-0.5);
  derivs(2,2) =  rs * rp * sp *(-0.5);
  derivs(2,3) = -rs * rm * sp *(-0.5);
  derivs(2,4) = -sh * sm * r2 *(-0.5);
  derivs(2,5) =  rh * rp * s2 *(-0.5);
  derivs(2,6) =  sh * sp * r2 *(-0.5);
  derivs(2,7) = -rh * rm * s2 *(-0.5);
  derivs(2,8) =  r2 * s2      *(-0.5);

  derivs(2,9)  =  rs * rm * sm *0.5;
  derivs(2,10) = -rs * rp * sm *0.5;
  derivs(2,11) =  rs * rp * sp *0.5;
  derivs(2,12) = -rs * rm * sp *0.5;
  derivs(2,13) = -sh * sm * r2 *0.5;
  derivs(2,14) =  rh * rp * s2 *0.5;
  derivs(2,15) =  sh * sp * r2 *0.5;
  derivs(2,16) = -rh * rm * s2 *0.5;
  derivs(2,17) =  r2 * s2      *0.5;

  return derivs;
}

const LINALG::Matrix<2,9> DRT::ELEMENTS::So_sh18::sh18_derivs_q9(
    const double r,
    const double s)
{
  LINALG::Matrix<2,9> derivs;

  const double rp = 1.0 + r;
  const double rm = 1.0 - r;
  const double sp = 1.0 + s;
  const double sm = 1.0 - s;
  const double r2 = 1.0 - r * r;
  const double s2 = 1.0 - s * s;
  const double rh = 0.5 * r;
  const double sh = 0.5 * s;
  const double rhp = r + 0.5;
  const double rhm = r - 0.5;
  const double shp = s + 0.5;
  const double shm = s - 0.5;

  derivs(0, 0) = -rhm * sh * sm;
  derivs(1, 0) = -shm * rh * rm;
  derivs(0, 1) = -rhp * sh * sm;
  derivs(1, 1) =  shm * rh * rp;
  derivs(0, 2) =  rhp * sh * sp;
  derivs(1, 2) =  shp * rh * rp;
  derivs(0, 3) =  rhm * sh * sp;
  derivs(1, 3) = -shp * rh * rm;
  derivs(0, 4) =  2.0 * r * sh * sm;
  derivs(1, 4) =  shm * r2;
  derivs(0, 5) =  rhp * s2;
  derivs(1, 5) = -2.0 * s * rh * rp;
  derivs(0, 6) = -2.0 * r * sh * sp;
  derivs(1, 6) =  shp * r2;
  derivs(0, 7) =  rhm * s2;
  derivs(1, 7) =  2.0 * s * rh * rm;
  derivs(0, 8) = -2.0 * r * s2;
  derivs(1, 8) = -2.0 * s * r2;
  return derivs;
}

/*----------------------------------------------------------------------*
 |  init the element (public)                               seitz 11/14 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh18Type::Initialize(DRT::Discretization& dis)
{
  // here we order the nodes such that we have a positive definite jacobian
  //       maybe the python script generating the hex18 elements would be a better place for this.
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_sh18* actele = dynamic_cast<DRT::ELEMENTS::So_sh18*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex18* failed");
    if (actele->InitJacobianMapping()==1)
      actele->FlipT();
  }
  dis.FillComplete(false,false,false);

  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_sh18* actele = dynamic_cast<DRT::ELEMENTS::So_sh18*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex18* failed");
    if (actele->InitJacobianMapping()==1)
      dserror("why");
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  revert the 3rd parameter direction                      seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::FlipT()
{
  if (NodeIds()==NULL) dserror("couldn't get node ids");
  // reorder nodes
  int new_nodeids[NUMNOD_SOH18];
  new_nodeids[0] = NodeIds()[9];
  new_nodeids[1] = NodeIds()[10];
  new_nodeids[2] = NodeIds()[11];
  new_nodeids[3] = NodeIds()[12];
  new_nodeids[4] = NodeIds()[13];
  new_nodeids[5] = NodeIds()[14];
  new_nodeids[6] = NodeIds()[15];
  new_nodeids[7] = NodeIds()[16];
  new_nodeids[8] = NodeIds()[17];

  new_nodeids[9]  = NodeIds()[0];
  new_nodeids[10] = NodeIds()[1];
  new_nodeids[11] = NodeIds()[2];
  new_nodeids[12] = NodeIds()[3];
  new_nodeids[13] = NodeIds()[4];
  new_nodeids[14] = NodeIds()[5];
  new_nodeids[15] = NodeIds()[6];
  new_nodeids[16] = NodeIds()[7];
  new_nodeids[17] = NodeIds()[8];

  SetNodeIds(NUMNOD_SOH18,new_nodeids);
  return;
}


void DRT::ELEMENTS::So_sh18::EvaluateT(const LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18>& jac,
                                              LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>& TinvT)
{
  // build T^T transformation matrix which maps
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  TinvT(0,0) = jac(0,0) * jac(0,0);
  TinvT(1,0) = jac(1,0) * jac(1,0);
  TinvT(2,0) = jac(2,0) * jac(2,0);
  TinvT(3,0) = 2 * jac(0,0) * jac(1,0);
  TinvT(4,0) = 2 * jac(1,0) * jac(2,0);
  TinvT(5,0) = 2 * jac(0,0) * jac(2,0);

  TinvT(0,1) = jac(0,1) * jac(0,1);
  TinvT(1,1) = jac(1,1) * jac(1,1);
  TinvT(2,1) = jac(2,1) * jac(2,1);
  TinvT(3,1) = 2 * jac(0,1) * jac(1,1);
  TinvT(4,1) = 2 * jac(1,1) * jac(2,1);
  TinvT(5,1) = 2 * jac(0,1) * jac(2,1);

  TinvT(0,2) = jac(0,2) * jac(0,2);
  TinvT(1,2) = jac(1,2) * jac(1,2);
  TinvT(2,2) = jac(2,2) * jac(2,2);
  TinvT(3,2) = 2 * jac(0,2) * jac(1,2);
  TinvT(4,2) = 2 * jac(1,2) * jac(2,2);
  TinvT(5,2) = 2 * jac(0,2) * jac(2,2);

  TinvT(0,3) = jac(0,0) * jac(0,1);
  TinvT(1,3) = jac(1,0) * jac(1,1);
  TinvT(2,3) = jac(2,0) * jac(2,1);
  TinvT(3,3) = jac(0,0) * jac(1,1) + jac(1,0) * jac(0,1);
  TinvT(4,3) = jac(1,0) * jac(2,1) + jac(2,0) * jac(1,1);
  TinvT(5,3) = jac(0,0) * jac(2,1) + jac(2,0) * jac(0,1);


  TinvT(0,4) = jac(0,1) * jac(0,2);
  TinvT(1,4) = jac(1,1) * jac(1,2);
  TinvT(2,4) = jac(2,1) * jac(2,2);
  TinvT(3,4) = jac(0,1) * jac(1,2) + jac(1,1) * jac(0,2);
  TinvT(4,4) = jac(1,1) * jac(2,2) + jac(2,1) * jac(1,2);
  TinvT(5,4) = jac(0,1) * jac(2,2) + jac(2,1) * jac(0,2);

  TinvT(0,5) = jac(0,0) * jac(0,2);
  TinvT(1,5) = jac(1,0) * jac(1,2);
  TinvT(2,5) = jac(2,0) * jac(2,2);
  TinvT(3,5) = jac(0,0) * jac(1,2) + jac(1,0) * jac(0,2);
  TinvT(4,5) = jac(1,0) * jac(2,2) + jac(2,0) * jac(1,2);
  TinvT(5,5) = jac(0,0) * jac(2,2) + jac(2,0) * jac(0,2);

  // now evaluate T^{-T} with solver
  LINALG::FixedSizeSerialDenseSolver<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D,1> solve_for_inverseT;
  solve_for_inverseT.SetMatrix(TinvT);
  int err2 = solve_for_inverseT.Factor();
  int err = solve_for_inverseT.Invert();
  if ((err != 0) && (err2!=0)) dserror("Inversion of Tinv (Jacobian) failed");
  return;
}

void DRT::ELEMENTS::So_sh18::CalculateLocStrain(const LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18>& xcurr,
    const LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18>& xrefe,
    const LINALG::Matrix<9,1>& shape_q9,
    const LINALG::Matrix<2,9>& deriv_q9,
    const int gp,
    LINALG::Matrix<MAT::NUM_STRESS_3D,1>& lstrain)
{
  for (int dim=0; dim<3; ++dim)
    for (int k=0; k<9; ++k)
      for (int l=0; l<9; ++l)
      {
        if (dsg_membrane_)
        {
          // constant normal strain
          lstrain(0) += .125*dsg_membrane_r_[gp%9](k,l)
                       *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)+xcurr(l,dim))     // a_alphaalpha
                         -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)+xrefe(l,dim))     // A_alphaalpha
                        );
          lstrain(1) += .125*dsg_membrane_s_[gp%9](k,l)
                       *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)+xcurr(l,dim))     // a_alphaalpha
                         -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)+xrefe(l,dim))     // A_alphaalpha
                        );
          // constant in-plane shear strain

          lstrain(3) += .25*dsg_membrane_rs_[gp%9](k,l)
                        *(  (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)+xcurr(l,dim))           // a_01
                           -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)+xrefe(l,dim))           // A_01
                         );
        }
        else
        {
          // constant normal strain
          for (int alpha=0; alpha<2;++alpha)
            lstrain(alpha) += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l)
                              *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)+xcurr(l,dim))     // a_alphaalpha
                                -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)+xrefe(l,dim))     // A_alphaalpha
                               );

          // constant in-plane shear strain
          lstrain(3) += .25*deriv_q9(0,k)*deriv_q9(1,l)
                        *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)+xcurr(l,dim))           // a_01
                          -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)+xrefe(l,dim))           // A_01
                         );
        }

        // linear normal strain
          for (int alpha=0; alpha<2; ++alpha)
            lstrain(alpha) += xsi_[gp](2)
                              *.125*deriv_q9(alpha,k)*deriv_q9(alpha,l)
                              *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))
                                +(xcurr(k+9,dim)-xcurr(k,dim))*(xcurr(l+9,dim)+xcurr(l,dim))     // b_alphaalpha
                                -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))
                                -(xrefe(k+9,dim)-xrefe(k,dim))*(xrefe(l+9,dim)+xrefe(l,dim))     // B_alphaalpha
                               );

          // linear in-plane shear strain
          lstrain(3) += xsi_[gp](2)
                        *.25*deriv_q9(0,k)*deriv_q9(1,l)
                        *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))
                          +(xcurr(k+9,dim)-xcurr(k,dim))*(xcurr(l+9,dim)+xcurr(l,dim))     // b_01
                          -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))
                          -(xrefe(k+9,dim)-xrefe(k,dim))*(xrefe(l+9,dim)+xrefe(l,dim))     // B_01
                        );

          if (dsg_shear_)
          {
            // constant transverse shear strain
            lstrain(4) += .25*dsg_shear_s_[gp%9](k,l)
                          *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))
                            -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))
                           );
            lstrain(5) += .25*dsg_shear_r_[gp%9](k,l)
                          *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))
                            -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))
                           );
          }
          else
          {
            // constant transverse shear strain
            lstrain(4) += .25*deriv_q9(1,k)*shape_q9(l)
                          *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))     // a_12
                            -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))     // A_12
                           );
            lstrain(5) += .25*deriv_q9(0,k)*shape_q9(l)
                          *( (xcurr(k+9,dim)+xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))     // a_02
                            -(xrefe(k+9,dim)+xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))     // A_02
                          );
          }

          // linear transverse shear strain
         lstrain(4) += xsi_[gp](2)
                       *.25*deriv_q9(1,k)*shape_q9(l)
                       *( (xcurr(k+9,dim)-xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))    // a_2 dot a_2,1
                         -(xrefe(k+9,dim)-xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))    // A_2 dot A_2,1
                        );
         lstrain(5) += xsi_[gp](2)
                       *.25*deriv_q9(0,k)*shape_q9(l)
                       *( (xcurr(k+9,dim)-xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))    // a_2 dot a_2,0
                         -(xrefe(k+9,dim)-xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))    // A_2 dot A_2,0
                        );

         // transverse normal strain
         if (dsg_ctl_)
         {
           lstrain(2) += .125*dsg_transverse_t_[gp%9](k,l)
                         *( (xcurr(k+9,dim)-xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))    // a_2 dot a_2
                           -(xrefe(k+9,dim)-xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))    // A_2 dot A_2
                         );
         }
         else
         {
           lstrain(2) += .125*shape_q9(k)*shape_q9(l)
                         *( (xcurr(k+9,dim)-xcurr(k,dim))*(xcurr(l+9,dim)-xcurr(l,dim))    // a_2 dot a_2
                           -(xrefe(k+9,dim)-xrefe(k,dim))*(xrefe(l+9,dim)-xrefe(l,dim))    // A_2 dot A_2
                          );
         }
    } // k=0..8; l=0..8; dim=0..2

  return;
}

void DRT::ELEMENTS::So_sh18::CalculateBopLoc(const LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18>& xcurr,
    const LINALG::Matrix<NUMNOD_SOH18,NUMDIM_SOH18>& xrefe,
    const LINALG::Matrix<9,1>& shape_q9,
    const LINALG::Matrix<2,9>& deriv_q9,
    const int gp,
    LINALG::Matrix<MAT::NUM_STRESS_3D,NUMDOF_SOH18>& bop_loc)
{
  for (int dim = 0; dim < NUMDIM_SOH18; ++dim)
    for (int k=0; k<9; ++k)
      for (int l=0; l<9; ++l)
      {
        if (dsg_membrane_)
        {
          // constant normal strain
          bop_loc(0,(k+9)*3+dim) += .125*dsg_membrane_r_[gp%9](k,l)*(xcurr(l+9,dim)+xcurr(l,dim));
          bop_loc(0,k*3+dim)     += .125*dsg_membrane_r_[gp%9](k,l)*(xcurr(l+9,dim)+xcurr(l,dim));
          bop_loc(0,(l+9)*3+dim) += .125*dsg_membrane_r_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));
          bop_loc(0,l*3+dim)     += .125*dsg_membrane_r_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));

          bop_loc(1,(k+9)*3+dim) += .125*dsg_membrane_s_[gp%9](k,l)*(xcurr(l+9,dim)+xcurr(l,dim));
          bop_loc(1,k*3+dim)     += .125*dsg_membrane_s_[gp%9](k,l)*(xcurr(l+9,dim)+xcurr(l,dim));
          bop_loc(1,(l+9)*3+dim) += .125*dsg_membrane_s_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));
          bop_loc(1,l*3+dim)     += .125*dsg_membrane_s_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));

          // constant in-plane shear strain
          bop_loc(3,(k+9)*3+dim) += .25*dsg_membrane_rs_[gp%9](k,l)*(xcurr(l+9,dim)+xcurr(l,dim));
          bop_loc(3,k*3+dim)     += .25*dsg_membrane_rs_[gp%9](k,l)*(xcurr(l+9,dim)+xcurr(l,dim));
          bop_loc(3,(l+9)*3+dim) += .25*dsg_membrane_rs_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));
          bop_loc(3,l*3+dim)     += .25*dsg_membrane_rs_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));
        }
        else
        {
          // constant normal strain
          for (int alpha=0; alpha<2;++alpha)
          {
            bop_loc(alpha,(k+9)*3+dim) += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l)*(xcurr(l+9,dim)+xcurr(l,dim));
            bop_loc(alpha,k*3+dim)     += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l)*(xcurr(l+9,dim)+xcurr(l,dim));
            bop_loc(alpha,(l+9)*3+dim) += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l)*(xcurr(k+9,dim)+xcurr(k,dim));
            bop_loc(alpha,l*3+dim)     += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l)*(xcurr(k+9,dim)+xcurr(k,dim));
          }

          // constant in-plane shear strain
          bop_loc(3,(k+9)*3+dim) += .25*deriv_q9(0,k)*deriv_q9(1,l)*(xcurr(l+9,dim)+xcurr(l,dim));
          bop_loc(3,k*3+dim)     += .25*deriv_q9(0,k)*deriv_q9(1,l)*(xcurr(l+9,dim)+xcurr(l,dim));
          bop_loc(3,(l+9)*3+dim) += .25*deriv_q9(0,k)*deriv_q9(1,l)*(xcurr(k+9,dim)+xcurr(k,dim));
          bop_loc(3,l*3+dim)     += .25*deriv_q9(0,k)*deriv_q9(1,l)*(xcurr(k+9,dim)+xcurr(k,dim));
        }

        // linear normal strain
        for (int alpha=0; alpha<2; ++alpha)
        {
          bop_loc(alpha,(k+9)*3+dim) += xsi_[gp](2)*.25*deriv_q9(alpha,k)*deriv_q9(alpha,l)*xcurr(l+9,dim);
          bop_loc(alpha,k*3+dim)     +=-xsi_[gp](2)*.25*deriv_q9(alpha,k)*deriv_q9(alpha,l)*xcurr(l,dim);
          bop_loc(alpha,(l+9)*3+dim) += xsi_[gp](2)*.25*deriv_q9(alpha,k)*deriv_q9(alpha,l)*xcurr(k+9,dim);
          bop_loc(alpha,l*3+dim)     +=-xsi_[gp](2)*.25*deriv_q9(alpha,k)*deriv_q9(alpha,l)*xcurr(k,dim);
        }

        // linear in-plane shear strain
        bop_loc(3,(k+9)*3+dim) += xsi_[gp](2)*.5*deriv_q9(0,k)*deriv_q9(1,l)*xcurr(l+9,dim);
        bop_loc(3,k*3+dim)     +=-xsi_[gp](2)*.5*deriv_q9(0,k)*deriv_q9(1,l)*xcurr(l,dim);
        bop_loc(3,(l+9)*3+dim) += xsi_[gp](2)*.5*deriv_q9(0,k)*deriv_q9(1,l)*xcurr(k+9,dim);
        bop_loc(3,l*3+dim)     +=-xsi_[gp](2)*.5*deriv_q9(0,k)*deriv_q9(1,l)*xcurr(k,dim);

        if (dsg_shear_)
        {
          // constant transverse shear strain
          bop_loc(4,(k+9)*3+dim) += .25*dsg_shear_s_[gp%9](k,l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(4,k*3+dim)     += .25*dsg_shear_s_[gp%9](k,l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(4,(l+9)*3+dim) += .25*dsg_shear_s_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));
          bop_loc(4,l*3+dim)     +=-.25*dsg_shear_s_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));

          bop_loc(5,(k+9)*3+dim) += .25*dsg_shear_r_[gp%9](k,l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(5,k*3+dim)     += .25*dsg_shear_r_[gp%9](k,l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(5,(l+9)*3+dim) += .25*dsg_shear_r_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));
          bop_loc(5,l*3+dim)     +=-.25*dsg_shear_r_[gp%9](k,l)*(xcurr(k+9,dim)+xcurr(k,dim));
        }
        else
        {
          // constant transverse shear strain
          bop_loc(4,(k+9)*3+dim) += .25*deriv_q9(1,k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(4,k*3+dim)     += .25*deriv_q9(1,k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(4,(l+9)*3+dim) += .25*deriv_q9(1,k)*shape_q9(l)*(xcurr(k+9,dim)+xcurr(k,dim));
          bop_loc(4,l*3+dim)     +=-.25*deriv_q9(1,k)*shape_q9(l)*(xcurr(k+9,dim)+xcurr(k,dim));
          bop_loc(5,(k+9)*3+dim) += .25*deriv_q9(0,k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(5,k*3+dim)     += .25*deriv_q9(0,k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(5,(l+9)*3+dim) += .25*deriv_q9(0,k)*shape_q9(l)*(xcurr(k+9,dim)+xcurr(k,dim));
          bop_loc(5,l*3+dim)     +=-.25*deriv_q9(0,k)*shape_q9(l)*(xcurr(k+9,dim)+xcurr(k,dim));
        }

        // linear transverse shear strain
        bop_loc(4,(k+9)*3+dim) += xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
        bop_loc(4,k*3+dim)     +=-xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
        bop_loc(4,(l+9)*3+dim) += xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l)*(xcurr(k+9,dim)-xcurr(k,dim));
        bop_loc(4,l*3+dim)     +=-xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l)*(xcurr(k+9,dim)-xcurr(k,dim));
        bop_loc(5,(k+9)*3+dim) += xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
        bop_loc(5,k*3+dim)     +=-xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
        bop_loc(5,(l+9)*3+dim) += xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l)*(xcurr(k+9,dim)-xcurr(k,dim));
        bop_loc(5,l*3+dim)     +=-xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l)*(xcurr(k+9,dim)-xcurr(k,dim));

        // transverse normal strain
        if (dsg_ctl_)
        {
          bop_loc(2,(k+9)*3+dim) += .125*dsg_transverse_t_[gp%9](k,l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(2,k*3+dim)     -= .125*dsg_transverse_t_[gp%9](k,l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(2,(l+9)*3+dim) += .125*dsg_transverse_t_[gp%9](k,l)*(xcurr(k+9,dim)-xcurr(k,dim));
          bop_loc(2,l*3+dim)     -= .125*dsg_transverse_t_[gp%9](k,l)*(xcurr(k+9,dim)-xcurr(k,dim));
        }
        else
        {
          bop_loc(2,(k+9)*3+dim) += .125*shape_q9(k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(2,k*3+dim)     +=-.125*shape_q9(k)*shape_q9(l)*(xcurr(l+9,dim)-xcurr(l,dim));
          bop_loc(2,(l+9)*3+dim) += .125*shape_q9(k)*shape_q9(l)*(xcurr(k+9,dim)-xcurr(k,dim));
          bop_loc(2,l*3+dim)     +=-.125*shape_q9(k)*shape_q9(l)*(xcurr(k+9,dim)-xcurr(k,dim));
        }
    } // k=0..8; l=0..8; dim=0..2
  return;
}

void DRT::ELEMENTS::So_sh18::CalculateGeoStiff(
    const LINALG::Matrix<9,1>& shape_q9,
    const LINALG::Matrix<2,9>& deriv_q9,
    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>& TinvT,
    const int gp,
    const double detJ_w,
    const LINALG::Matrix<MAT::NUM_STRESS_3D,1>& stress,
    LINALG::Matrix<NUMDOF_SOH18,NUMDOF_SOH18>* stiffmatrix)
{
  // intergrate `geometric' stiffness matrix and add to keu *****************
  for (int k=0; k<9; ++k)
    for (int l=0; l<9; ++l)
    {
      LINALG::Matrix<6,1> G_kl,G_klp,G_kpl,G_kplp,G_lk,G_lkp,G_lpk,G_lpkp;

      // Normalverzerrungen
      if (dsg_membrane_)
      {
        // constant normal strain
        G_kl(0)    += .125*dsg_membrane_r_[gp%9](k,l);
        G_klp(0)   += .125*dsg_membrane_r_[gp%9](k,l);
        G_kpl(0)   += .125*dsg_membrane_r_[gp%9](k,l);
        G_kplp(0)  += .125*dsg_membrane_r_[gp%9](k,l);
        G_lk(0)    += .125*dsg_membrane_r_[gp%9](k,l);
        G_lkp(0)   += .125*dsg_membrane_r_[gp%9](k,l);
        G_lpk(0)   += .125*dsg_membrane_r_[gp%9](k,l);
        G_lpkp(0)  += .125*dsg_membrane_r_[gp%9](k,l);

        G_kl(1)    += .125*dsg_membrane_s_[gp%9](k,l);
        G_klp(1)   += .125*dsg_membrane_s_[gp%9](k,l);
        G_kpl(1)   += .125*dsg_membrane_s_[gp%9](k,l);
        G_kplp(1)  += .125*dsg_membrane_s_[gp%9](k,l);
        G_lk(1)    += .125*dsg_membrane_s_[gp%9](k,l);
        G_lkp(1)   += .125*dsg_membrane_s_[gp%9](k,l);
        G_lpk(1)   += .125*dsg_membrane_s_[gp%9](k,l);
        G_lpkp(1)  += .125*dsg_membrane_s_[gp%9](k,l);

        // constant in-plane shear strain
        G_kl(3)           += .25*dsg_membrane_rs_[gp%9](k,l);
        G_klp(3)          += .25*dsg_membrane_rs_[gp%9](k,l);
        G_kpl(3)          += .25*dsg_membrane_rs_[gp%9](k,l);
        G_kplp(3)         += .25*dsg_membrane_rs_[gp%9](k,l);
        G_lk(3)           += .25*dsg_membrane_rs_[gp%9](k,l);
        G_lkp(3)          += .25*dsg_membrane_rs_[gp%9](k,l);
        G_lpk(3)          += .25*dsg_membrane_rs_[gp%9](k,l);
        G_lpkp(3)         += .25*dsg_membrane_rs_[gp%9](k,l);
      }
      else
      {
        // constant normal strain
        for (int alpha=0; alpha<2; ++alpha)
        {
          G_kl(alpha)          += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l);
          G_klp(alpha)         += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l);
          G_kpl(alpha)         += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l);
          G_kplp(alpha)        += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l);
          G_lk(alpha)          += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l);
          G_lkp(alpha)         += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l);
          G_lpk(alpha)         += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l);
          G_lpkp(alpha)        += .125*deriv_q9(alpha,k)*deriv_q9(alpha,l);
        }

        // constant in-plane shear strain
        G_kl(3)           += .25*deriv_q9(0,k)*deriv_q9(1,l);
        G_klp(3)          += .25*deriv_q9(0,k)*deriv_q9(1,l);
        G_kpl(3)          += .25*deriv_q9(0,k)*deriv_q9(1,l);
        G_kplp(3)         += .25*deriv_q9(0,k)*deriv_q9(1,l);
        G_lk(3)           += .25*deriv_q9(0,k)*deriv_q9(1,l);
        G_lkp(3)          += .25*deriv_q9(0,k)*deriv_q9(1,l);
        G_lpk(3)          += .25*deriv_q9(0,k)*deriv_q9(1,l);
        G_lpkp(3)         += .25*deriv_q9(0,k)*deriv_q9(1,l);
      }

      // linear normal strain
      for (int alpha=0; alpha<2; ++alpha)
      {
        G_kl(alpha)     -= xsi_[gp](2)*.25*deriv_q9(alpha,k)*deriv_q9(alpha,l);
        G_klp(alpha)    +=0.;
        G_kpl(alpha)    +=0.;
        G_kplp(alpha)   += xsi_[gp](2)*.25*deriv_q9(alpha,k)*deriv_q9(alpha,l);
        G_lk(alpha)     -= xsi_[gp](2)*.25*deriv_q9(alpha,k)*deriv_q9(alpha,l);
        G_lkp(alpha)    +=0.;
        G_lpk(alpha)    +=0.;
        G_lpkp(alpha)   += xsi_[gp](2)*.25*deriv_q9(alpha,k)*deriv_q9(alpha,l);
      }

      // linear in-plane shear strain
      G_kl(3)       -= xsi_[gp](2)*.5*deriv_q9(0,k)*deriv_q9(1,l);
      G_klp(3)      +=0.;
      G_kpl(3)      +=0.;
      G_kplp(3)     += xsi_[gp](2)*.5*deriv_q9(0,k)*deriv_q9(1,l);
      G_lk(3)       -= xsi_[gp](2)*.5*deriv_q9(0,k)*deriv_q9(1,l);
      G_lkp(3)      +=0.;
      G_lpk(3)      +=0.;
      G_lpkp(3)     += xsi_[gp](2)*.5*deriv_q9(0,k)*deriv_q9(1,l);

      if (dsg_shear_)
      {
        // constant transverse shear strain
        G_kl(4)    -= .25*dsg_shear_s_[gp%9](k,l);
        G_klp(4)   += .25*dsg_shear_s_[gp%9](k,l);
        G_kpl(4)   -= .25*dsg_shear_s_[gp%9](k,l);
        G_kplp(4)  += .25*dsg_shear_s_[gp%9](k,l);
        G_lk(4)    -= .25*dsg_shear_s_[gp%9](k,l);
        G_lkp(4)   -= .25*dsg_shear_s_[gp%9](k,l);
        G_lpk(4)   += .25*dsg_shear_s_[gp%9](k,l);
        G_lpkp(4)  += .25*dsg_shear_s_[gp%9](k,l);

        G_kl(5)    -= .25*dsg_shear_r_[gp%9](k,l);
        G_klp(5)   += .25*dsg_shear_r_[gp%9](k,l);
        G_kpl(5)   -= .25*dsg_shear_r_[gp%9](k,l);
        G_kplp(5)  += .25*dsg_shear_r_[gp%9](k,l);
        G_lk(5)    -= .25*dsg_shear_r_[gp%9](k,l);
        G_lkp(5)   -= .25*dsg_shear_r_[gp%9](k,l);
        G_lpk(5)   += .25*dsg_shear_r_[gp%9](k,l);
        G_lpkp(5)  += .25*dsg_shear_r_[gp%9](k,l);
      }
      else
      {
        // constant transverse shear strain
        G_kl(4)      -= .25*deriv_q9(1,k)*shape_q9(l);
        G_klp(4)     += .25*deriv_q9(1,k)*shape_q9(l);
        G_kpl(4)     -= .25*deriv_q9(1,k)*shape_q9(l);
        G_kplp(4)    += .25*deriv_q9(1,k)*shape_q9(l);
        G_lk(4)      -= .25*deriv_q9(1,k)*shape_q9(l);
        G_lkp(4)     -= .25*deriv_q9(1,k)*shape_q9(l);
        G_lpk(4)     += .25*deriv_q9(1,k)*shape_q9(l);
        G_lpkp(4)    += .25*deriv_q9(1,k)*shape_q9(l);

        G_kl(5)      -= .25*deriv_q9(0,k)*shape_q9(l);
        G_klp(5)     += .25*deriv_q9(0,k)*shape_q9(l);
        G_kpl(5)     -= .25*deriv_q9(0,k)*shape_q9(l);
        G_kplp(5)    += .25*deriv_q9(0,k)*shape_q9(l);
        G_lk(5)      -= .25*deriv_q9(0,k)*shape_q9(l);
        G_lkp(5)     -= .25*deriv_q9(0,k)*shape_q9(l);
        G_lpk(5)     += .25*deriv_q9(0,k)*shape_q9(l);
        G_lpkp(5)    += .25*deriv_q9(0,k)*shape_q9(l);
      }

      // linear transverse shear strain
      G_kl(4)       += xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l);
      G_klp(4)      -= xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l);
      G_kpl(4)      -= xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l);
      G_kplp(4)     += xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l);
      G_lk(4)       += xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l);
      G_lkp(4)      -= xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l);
      G_lpk(4)      -= xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l);
      G_lpkp(4)     += xsi_[gp](2)*.25*deriv_q9(1,k)*shape_q9(l);


      G_kl(5)       += xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l);
      G_klp(5)      -= xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l);
      G_kpl(5)      -= xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l);
      G_kplp(5)     += xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l);
      G_lk(5)       += xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l);
      G_lkp(5)      -= xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l);
      G_lpk(5)      -= xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l);
      G_lpkp(5)     += xsi_[gp](2)*.25*deriv_q9(0,k)*shape_q9(l);

      // transverse normal strain
      if (dsg_ctl_)
      {
        G_kl(2)     += .125*dsg_transverse_t_[gp%9](k,l);
        G_klp(2)    -= .125*dsg_transverse_t_[gp%9](k,l);
        G_kpl(2)    -= .125*dsg_transverse_t_[gp%9](k,l);
        G_kplp(2)   += .125*dsg_transverse_t_[gp%9](k,l);
        G_lk(2)     += .125*dsg_transverse_t_[gp%9](k,l);
        G_lkp(2)    -= .125*dsg_transverse_t_[gp%9](k,l);
        G_lpk(2)    -= .125*dsg_transverse_t_[gp%9](k,l);
        G_lpkp(2)   += .125*dsg_transverse_t_[gp%9](k,l);
      }
      else
      {
        G_kl(2)      += .125*shape_q9(k)*shape_q9(l);
        G_klp(2)     -= .125*shape_q9(k)*shape_q9(l);
        G_kpl(2)     -= .125*shape_q9(k)*shape_q9(l);
        G_kplp(2)    += .125*shape_q9(k)*shape_q9(l);
        G_lk(2)      += .125*shape_q9(k)*shape_q9(l);
        G_lkp(2)     -= .125*shape_q9(k)*shape_q9(l);
        G_lpk(2)     -= .125*shape_q9(k)*shape_q9(l);
        G_lpkp(2)    += .125*shape_q9(k)*shape_q9(l);
      }

      LINALG::Matrix<6,1> G_kl_g;       G_kl_g  .Multiply(TinvT,G_kl);    const double Gkl    = detJ_w*stress.Dot(G_kl_g);
      LINALG::Matrix<6,1> G_klp_g;      G_klp_g .Multiply(TinvT,G_klp);   const double Gklp   = detJ_w*stress.Dot(G_klp_g);
      LINALG::Matrix<6,1> G_kpl_g;      G_kpl_g .Multiply(TinvT,G_kpl);   const double Gkpl   = detJ_w*stress.Dot(G_kpl_g);
      LINALG::Matrix<6,1> G_kplp_g;     G_kplp_g.Multiply(TinvT,G_kplp);  const double Gkplp  = detJ_w*stress.Dot(G_kplp_g);
      LINALG::Matrix<6,1> G_lk_g;       G_lk_g  .Multiply(TinvT,G_lk);    const double Glk    = detJ_w*stress.Dot(G_lk_g);
      LINALG::Matrix<6,1> G_lkp_g;      G_lkp_g .Multiply(TinvT,G_lkp);   const double Glkp   = detJ_w*stress.Dot(G_lkp_g);
      LINALG::Matrix<6,1> G_lpk_g;      G_lpk_g .Multiply(TinvT,G_lpk);   const double Glpk   = detJ_w*stress.Dot(G_lpk_g);
      LINALG::Matrix<6,1> G_lpkp_g;     G_lpkp_g.Multiply(TinvT,G_lpkp);  const double Glpkp  = detJ_w*stress.Dot(G_lpkp_g);
      for (int dim=0; dim<3; ++dim)
      {
        (*stiffmatrix)(k*3+dim,l*3+dim)         += Gkl;
        (*stiffmatrix)(k*3+dim,(l+9)*3+dim)     += Gklp;
        (*stiffmatrix)((k+9)*3+dim,l*3+dim)     += Gkpl;
        (*stiffmatrix)((k+9)*3+dim,(l+9)*3+dim) += Gkplp;

        (*stiffmatrix)(l*3+dim,k*3+dim)         += Glk;
        (*stiffmatrix)(l*3+dim,(k+9)*3+dim)     += Glkp;
        (*stiffmatrix)((l+9)*3+dim,k*3+dim)     += Glpk;
        (*stiffmatrix)((l+9)*3+dim,(k+9)*3+dim) += Glpkp;
      }
  } // k=0..8 l=0..8
  // end of integrate `geometric' stiffness******************************
}

void DRT::ELEMENTS::So_sh18::CalcConsistentDefgrd(LINALG::Matrix<3,3> defgrd_disp,
    LINALG::Matrix<6,1> glstrain_mod,
    LINALG::Matrix<3,3>& defgrd_mod)
{
  LINALG::Matrix<3,3> R;      // rotation tensor
  LINALG::Matrix<3,3> U_mod;  // modified right stretch tensor
  LINALG::Matrix<3,3> U_disp; // displacement-based right stretch tensor
  LINALG::Matrix<3,3> EW;     // temporarily store eigenvalues
  LINALG::Matrix<3,3> tmp;    // temporary matrix for matrix matrix matrix products
  LINALG::Matrix<3,3> tmp2;    // temporary matrix for matrix matrix matrix products

  // ******************************************************************
  // calculate modified right stretch tensor
  // ******************************************************************
  for (int i=0; i<3; i++)
    U_mod(i,i) = 2.*glstrain_mod(i) + 1.;
  U_mod(0,1) = glstrain_mod(3);
  U_mod(1,0) = glstrain_mod(3);
  U_mod(1,2) = glstrain_mod(4);
  U_mod(2,1) = glstrain_mod(4);
  U_mod(0,2) = glstrain_mod(5);
  U_mod(2,0) = glstrain_mod(5);

  LINALG::SYEV(U_mod,EW,U_mod);
  for (int i=0; i<3; ++i)
    EW(i,i) = sqrt(EW(i,i));
  tmp.Multiply(U_mod,EW);
  tmp2.MultiplyNT(tmp,U_mod);
  U_mod.Update(tmp2);

  // ******************************************************************
  // calculate displacement-based right stretch tensor
  // ******************************************************************
  U_disp.MultiplyTN(defgrd_disp,defgrd_disp);

  LINALG::SYEV(U_disp,EW,U_disp);
  for (int i=0; i<3; ++i)
    EW(i,i) = sqrt(EW(i,i));
  tmp.Multiply(U_disp,EW);
  tmp2.MultiplyNT(tmp,U_disp);
  U_disp.Update(tmp2);

  // ******************************************************************
  // compose consistent deformation gradient
  // ******************************************************************
  U_disp.Invert();
  R.Multiply(defgrd_disp,U_disp);
  defgrd_mod.Multiply(R,U_mod);

  // you're done here
  return;

}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Integrate_dsg_shear_r(const double xi,const double eta, LINALG::Matrix<9,9>& dsg_shear_r)
{
  dsg_shear_r.Clear();
  const double coord_refNode[2]={0.,0.};

  const LINALG::Matrix<2,9> deriv = sh18_derivs_q9(xi,eta);

  DRT::UTILS::IntPointsAndWeights<1> ip(DRT::UTILS::intrule_line_2point);
  for (int i=0; i<9; ++i)
  {
    const LINALG::Matrix<3,1> coord_i=NodeParamCoord(i);

    // integrations with non-empty integration domain
    if (coord_i(0)!=coord_refNode[0])
    // perform integration
      for (int gp=0; gp<ip.IP().nquad; ++gp)
      {
        const double jac = .5*(coord_i(0)-coord_refNode[0]);

        // gauss point coordinates in element parameter space
        double xi_gp[2];
        xi_gp[0] = .5*(coord_i(0)+coord_refNode[0])+jac*ip.IP().qxg[gp][0];
        xi_gp[1] = coord_i(1);

        // shape function
        const LINALG::Matrix<9,1> shape_gp = sh18_shapefcts_q9(xi_gp[0],xi_gp[1]);
        // derivative
        const LINALG::Matrix<2,9> deriv_gp = sh18_derivs_q9(xi_gp[0],xi_gp[1]);

        for (int k=0; k<9; ++k)
          for (int l=0; l<9; ++l)
            dsg_shear_r(k,l) += deriv(0,i)*deriv_gp(0,k)*shape_gp(l)*jac*ip.IP().qwgt[gp];
      }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Integrate_dsg_shear_s(const double xi,const double eta, LINALG::Matrix<9,9>& dsg_shear_s)
{
  dsg_shear_s.Clear();
  const double coord_refNode[2]={0.,0.};

  const LINALG::Matrix<2,9> deriv = sh18_derivs_q9(xi,eta);

  DRT::UTILS::IntPointsAndWeights<1> ip(DRT::UTILS::intrule_line_2point);
  for (int i=0; i<9; ++i)
  {
    const LINALG::Matrix<3,1> coord_i=NodeParamCoord(i);

    // integrations with non-empty integration domain
    if (coord_i(1)!=coord_refNode[1])
      // perform integration
        for (int gp=0; gp<ip.IP().nquad; ++gp)
        {
          // integration jacobian
          const double jac=.5*(coord_i(1)-coord_refNode[1]);

          // gauss point coordinates in element parameter space
          double xi_gp[2];
          xi_gp[0] = coord_i(0);
          xi_gp[1] = .5*(coord_i(1)+coord_refNode[1])+jac*ip.IP().qxg[gp][0];

          // shape function
          LINALG::Matrix<9,1> shape_gp = sh18_shapefcts_q9(xi_gp[0],xi_gp[1]);
          // derivative
          LINALG::Matrix<2,9> deriv_gp = sh18_derivs_q9(xi_gp[0],xi_gp[1]);

          for (int k=0; k<9; ++k)
            for (int l=0; l<9; ++l)
              dsg_shear_s(k,l) += deriv(1,i)*deriv_gp(1,k)*shape_gp(l)*jac*ip.IP().qwgt[gp];
        }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Integrate_dsg_membrane_rs(const double xi,const double eta, LINALG::Matrix<9,9>& dsg_membrane_rs)
{
  dsg_membrane_rs.Clear();

  // integration
  const double coord_refNode[2]={0.,0.};
  const LINALG::Matrix<18,3> coords = NodeParamCoord();
  DRT::UTILS::IntPointsAndWeights<1> ip(DRT::UTILS::intrule_line_2point);
  const LINALG::Matrix<2,9> deriv_xieta = sh18_derivs_q9(xi,eta);
  for (int k=0; k<9; ++k)
    for (int l=0; l<9; ++l)
    {
      for (int r=0; r<9; ++r)
        for (int s=0; s<9; ++s)

          for (int g=0; g<ip.IP().nquad; ++g)
            for (int h=0; h<ip.IP().nquad; ++h)
            {
              const double jac_g = .5*(coords(r,0)-coord_refNode[0]);
              const double jac_h = .5*(coords(s,1)-coord_refNode[1]);
              const double g_loc = .5*(coords(r,0)+coord_refNode[0])+jac_g*ip.IP().qxg[g][0];
              const double h_loc = .5*(coords(s,1)+coord_refNode[1])+jac_h*ip.IP().qxg[h][0];

              const LINALG::Matrix<2,9> deriv_g_eta = sh18_derivs_q9(g_loc,eta);
              const LINALG::Matrix<2,9> deriv_g_h   = sh18_derivs_q9(g_loc,h_loc);

              if (coords(r,0)!=coord_refNode[0] && coords(s,1)!=coord_refNode[1])
                dsg_membrane_rs(k,l) += deriv_xieta(0,r)*deriv_g_eta(1,s)*deriv_g_h(0,k)*deriv_g_h(1,l)
                        *jac_g*ip.IP().qwgt[g]
                        *jac_h*ip.IP().qwgt[h];
            }
    }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Integrate_dsg_membrane_r(const double xi,const double eta, LINALG::Matrix<9,9>& dsg_membrane_r)
{
  dsg_membrane_r.Clear();
  const double coord_refNode[2]={0.,0.};

  const LINALG::Matrix<2,9> deriv = sh18_derivs_q9(xi,eta);

  DRT::UTILS::IntPointsAndWeights<1> ip(DRT::UTILS::intrule_line_2point);
  for (int i=0; i<9; ++i)
  {
    const LINALG::Matrix<3,1> coord_i=NodeParamCoord(i);

    // integrations with non-empty integration domain
    if (coord_i(0)!=coord_refNode[0])
      // perform integration
        for (int gp=0; gp<ip.IP().nquad; ++gp)
        {
          // integration jacobian
          const double jac=.5*(coord_i(0)-coord_refNode[0]);

          // gauss point coordinates in element parameter space
          double xi_gp[2];
          xi_gp[0] = .5*(coord_i(0)+coord_refNode[0])+jac*ip.IP().qxg[gp][0];
          xi_gp[1] = coord_i(1);

          // derivative
          const LINALG::Matrix<2,9> deriv_gp = sh18_derivs_q9(xi_gp[0],xi_gp[1]);

          // fill up array
          for (int k=0; k<9; ++k)
            for (int l=0; l<9; ++l)
              dsg_membrane_r(k,l) += deriv(0,i)*deriv_gp(0,k)*deriv_gp(0,l)*jac*ip.IP().qwgt[gp];
        }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Integrate_dsg_membrane_s(const double xi,const double eta, LINALG::Matrix<9,9>& dsg_membrane_s)
{
  dsg_membrane_s.Clear();
  const double coord_refNode[2]={0.,0.};

  const LINALG::Matrix<2,9> deriv = sh18_derivs_q9(xi,eta);

  DRT::UTILS::IntPointsAndWeights<1> ip(DRT::UTILS::intrule_line_2point);
  for (int i=0; i<9; ++i)
  {
    const LINALG::Matrix<3,1> coord_i=NodeParamCoord(i);

    // integrations with non-empty integration domain
    if (coord_i(1)!=coord_refNode[1])
      // perform integration
        for (int gp=0; gp<ip.IP().nquad; ++gp)
        {
          // integration jacobian
          const double jac=.5*(coord_i(1)-coord_refNode[1]);

          // gauss point coordinates in element parameter space
          double xi_gp[2];
          xi_gp[0] = coord_i(0);
          xi_gp[1] = .5*(coord_i(1)+coord_refNode[1])+jac*ip.IP().qxg[gp][0];

          // derivative
          LINALG::Matrix<2,9> deriv_gp = sh18_derivs_q9(xi_gp[0],xi_gp[1]);

          // fill up array
          for (int k=0; k<9; ++k)
            for (int l=0; l<9; ++l)
              dsg_membrane_s(k,l) += deriv(1,i)*deriv_gp(1,k)*deriv_gp(1,l)*jac*ip.IP().qwgt[gp];
        }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Integrate_dsg_transverse_t(const double xi,const double eta, LINALG::Matrix<9,9>& dsg_transverse_t)
{
  // reset
  dsg_transverse_t.Clear();
  const LINALG::Matrix<9,1> shape=sh18_shapefcts_q9(xi,eta);

  for (int i=0; i<9; ++i)
  {
    const LINALG::Matrix<3,1> coord_i=NodeParamCoord(i);
    LINALG::Matrix<9,1> shape_gp = sh18_shapefcts_q9(coord_i(0),coord_i(1));
    for (int k=0; k<9; ++k)
      for (int l=0; l<9; ++l)
        dsg_transverse_t(k,l) += shape(i)*shape_gp(k)*shape_gp(l);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  setup all necessary DSG terms                           seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::SetupDSG()
{
  if (dsg_shear_)
  {
    dsg_shear_r_     .resize(9);
    dsg_shear_s_     .resize(9);
  }
  else
  {
    dsg_shear_r_     .resize(0);
    dsg_shear_s_     .resize(0);
  }
  if (dsg_membrane_)
  {
    dsg_membrane_r_  .resize(9);
    dsg_membrane_s_  .resize(9);
    dsg_membrane_rs_ .resize(9);
  }
  else
  {
    dsg_membrane_r_  .resize(0);
    dsg_membrane_s_  .resize(0);
    dsg_membrane_rs_ .resize(0);
  }
  if (dsg_ctl_)
    dsg_transverse_t_.resize(9);
  else
    dsg_transverse_t_.resize(0);

  // NOTE: For Lagrange polynomials the following arrays
  // are very sparse. To increase performance of Langrange
  // FE one might rather just code the non-zeros instead
  // of all those loops.
  // However, they represent what actually happens
  // for DSG elements. Moreover, other ansatz functions
  // e.g. Nurbs should hopefully be "unlocked" efficiently
  // without any further modification.
  for (int gp=0; gp<9; ++gp)
  {
    const double r=xsi_[gp](0);
    const double s=xsi_[gp](1);
    if (dsg_shear_)
    {
      Integrate_dsg_shear_r(r,s,dsg_shear_r_[gp]);
      Integrate_dsg_shear_s(r,s,dsg_shear_s_[gp]);
    }
    if (dsg_membrane_)
    {
      Integrate_dsg_membrane_r(r,s,dsg_membrane_r_[gp]);
      Integrate_dsg_membrane_s(r,s,dsg_membrane_s_[gp]);
      Integrate_dsg_membrane_rs(r,s,dsg_membrane_rs_[gp]);
    }
    if (dsg_ctl_)
    {
      Integrate_dsg_transverse_t(r,s,dsg_transverse_t_[gp]);
    }
  }

  return;
}

LINALG::Matrix<18,3> DRT::ELEMENTS::So_sh18::NodeParamCoord()
{
  LINALG::Matrix<18,3> coord;
  for (int node=0; node<NUMNOD_SOH18; ++node)
  {
    LINALG::Matrix<3,1> nodeCoord = NodeParamCoord(node);
    for (int i=0; i<3; ++i)
      coord(node,i) = nodeCoord(i);
  }
  return coord;
}

LINALG::Matrix<3,1> DRT::ELEMENTS::So_sh18::NodeParamCoord(const int node)
{
  LINALG::Matrix<3,1> coord;

  switch (node%9)
  {
  case 0:
  case 3:
  case 7:
    coord(0)=-1.;break;
  case 4:
  case 6:
  case 8:
    coord(0)=+0.;break;
  case 1:
  case 2:
  case 5:
    coord(0)=+1.;break;
  }
  switch (node%9)
  {
  case 0:
  case 1:
  case 4:
    coord(1)=-1.;break;
  case 5:
  case 7:
  case 8:
    coord(1)=+0.;break;
  case 2:
  case 3:
  case 6:
    coord(1)=+1.;break;
  }

  if (node<9) coord(2)=-1.;
  else        coord(2)=+1.;

  return coord;
}

/*----------------------------------------------------------------------*
 |  setup EAS terms                                         seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::EasSetup(
    std::vector<LINALG::Matrix<6,num_eas> >& M_gp, // M-matrix evaluated at GPs
    LINALG::Matrix<3,1>& G3_contra,
    const LINALG::Matrix<NUMNOD_SOH18,3> xrefe)    // material element coords
{
  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18> jac0inv;
  const LINALG::Matrix<2,9> deriv_q9      = sh18_derivs_q9(0,0);
    const LINALG::Matrix<9,1> shapefunct_q9 = sh18_shapefcts_q9(0,0);
  for (int dim=0; dim<3; ++dim)
    for (int k=0; k<9; ++k)
    {
      jac0inv(0,dim) += .5            *deriv_q9(0,k)*(xrefe(k+9,dim)+xrefe(k,dim));
      jac0inv(1,dim) += .5            *deriv_q9(1,k)*(xrefe(k+9,dim)+xrefe(k,dim));
      jac0inv(2,dim) += .5*shapefunct_q9(k)*(xrefe(k+9,dim)-xrefe(k,dim));
    }
  jac0inv.Invert();

  for (int dim=0;dim<3; ++dim)
    G3_contra(dim) = jac0inv(2,dim);

  // build EAS interpolation matrix M, evaluated at the GPs
  static std::vector<LINALG::Matrix<6,num_eas> > M(NUMGPT_SOH18);
  static bool eval;
  if (!eval)
  {
    for (int gp=0; gp<NUMGPT_SOH18; ++gp)
    {
      double r=xsi_[gp](0);
      double s=xsi_[gp](1);
      double t=xsi_[gp](2);
      M[gp](2,0) = 1.;
      M[gp](2,1) = r;
      M[gp](2,2) = s;
      M[gp](2,3) = r*s;
      M[gp](2,4) = 1.-3.*r*r;
      M[gp](2,5) = 1.-3.*s*s;
      M[gp](2,6) = r*r*s;
      M[gp](2,7) = r*s*s;
      M[gp](2,8) = 1.-9.*r*r*s*s;

      M[gp].Scale(t);
    }
    eval=true;
  }
  M_gp = M;

  return;
}
