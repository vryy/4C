/*!----------------------------------------------------------------------
\file so_hex8_thermo_lin_evaluate.cpp
\brief

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*----------------------------------------------------------------------*/

#include "so_hex8.H"

#include "../drt_lib/drt_utils.H"

// include the headers of temperature-dependent materials with history only
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/thermoplasticlinelast.H"
#include "../drt_mat/robinson.H"
#include "../drt_mat/plasticlinelast.H"
#include "../drt_mat/micromaterial.H"


/*----------------------------------------------------------------------*
 | evaluate the element (public)                             dano 05/10 |
 | originally by maf 04/07                                              |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::LinEvaluate(
  ParameterList& params,
  DRT::Discretization& discretization,
  DRT::Element::LocationArray& la,
  Epetra_SerialDenseMatrix& elemat1_epetra,
  Epetra_SerialDenseMatrix& elemat2_epetra,
  Epetra_SerialDenseVector& elevec1_epetra,
  Epetra_SerialDenseVector& elevec2_epetra,
  Epetra_SerialDenseVector& elevec3_epetra
  )
{

  // start with "none"
  DRT::ELEMENTS::So_hex8::ActionType act = So_hex8::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_internalforce") act = So_hex8::calc_struct_internalforce;
  else if (action=="calc_struct_nlnstiff")      act = So_hex8::calc_struct_nlnstiff;
  else if (action=="calc_struct_nlnstiffmass")  act = So_hex8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = So_hex8::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")        act = So_hex8::calc_struct_stress;
  else if (action=="calc_struct_update_istep")  act = So_hex8::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")   act = So_hex8::calc_struct_reset_istep;  // needed for TangDis predictor
  else if (action=="postprocess_stress")        act = So_hex8::postprocess_stress;
  else if (action=="calc_struct_stifftemp")     act = So_hex8::calc_struct_stifftemp;
  else dserror("Unknown type of action for So_hex8: %s",action.c_str());

  // what should the element do
  switch(act)
  {
  //==================================================================================
  // internal force vector only
  case calc_struct_internalforce:
  {
    // internal force vector
    LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
    // elemat1+2, elevec2+3 are not used anyway

    // need current displacement and residual forces
    RCP<const Epetra_Vector> disp = discretization.GetState(0,"displacement");
    RCP<const Epetra_Vector> res  = discretization.GetState(0,"residual displacement");
    if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
    // build the location vector only for the structure field
    vector<int> lm = la[0].lm_;
    vector<double> mydisp((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    vector<double> myres((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*res,myres,lm);
    // create a dummy element matrix to apply linearised EAS-stuff onto
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> myemat(true);

    // need current temperature state, call the temperature discretization
    // disassemble temperature
    if (discretization.HasState(1,"temperature"))
    {
      // check if you can get the temperature state
      Teuchos::RCP<const Epetra_Vector> tempnp
       = discretization.GetState(1,"temperature");
      if (tempnp == Teuchos::null)
        dserror("Cannot get state vector 'tempnp'");
      // call the temperature discretization in the location array
      std::vector<double> mytempnp((la[1].lm_).size());

      // the temperature field has only one dof per node, disregarded by the
      // dimension of the problem
      const int numdofpernode_ = NumDofPerNode(1,*(Nodes()[0]));
      // number of nodes per element
      const int nen_ = 8;
      if (la[1].Size() != nen_*numdofpernode_)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);

      // extract current temperatures declared as RCP<vector> 
      Teuchos::RCP<std::vector<double> >robtempnp = rcp(new std::vector<double>(la[1].lm_.size()) );
      DRT::UTILS::ExtractMyValues(*tempnp,*robtempnp,la[1].lm_);
      params.set<Teuchos::RCP<vector<double> > >("robinson_tempnp",robtempnp);

      // calculate the THERMOmechanical term for fint
      soh8_finttemp(la,mydisp,myres,mytempnp,&elevec1,NULL,NULL,params,
        INPAR::STR::stress_none,INPAR::STR::strain_none);
    }  // has temperature state

    // call the purely structural method
    linstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,NULL,params,
      INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);
  }
  break;

  //==================================================================================
  // linear stiffness and internal force vector
  case calc_struct_nlnstiff:
  {
    // stiffness
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
    // elemat2,elevec2+3 are not used anyway

    // need current displacement and residual forces
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0,"displacement");
    Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState(0,"residual displacement");
    if (disp==null || res==null)
      dserror("Cannot get state vectors 'displacement' and/or residual");
    vector<double> mydisp((la[0].lm_).size());
    // build the location vector only for the structure field
    vector<int> lm = la[0].lm_;
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);  // global, local, lm
    vector<double> myres((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*res,myres,lm);
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* matptr = NULL;
    // build a matrix dummy
    if (elemat1.IsInitialized()) matptr = &elemat1;

    // need current temperature state, call the temperature discretization
    // disassemble temperature
    if(discretization.HasState(1,"temperature"))
    {
      // check if you can get the temperature state
      Teuchos::RCP<const Epetra_Vector> tempnp
       = discretization.GetState(1,"temperature");
      if (tempnp == Teuchos::null)
        dserror("Cannot get state vector 'tempnp'");
      // call the temperature discretization in the location array
      std::vector<double> mytempnp((la[1].lm_).size());

      // the temperature field has only one dof per node, disregarded by the
      // dimension of the problem
      const int numdofpernode_ = NumDofPerNode(1,*(Nodes()[0]));
      // number of nodes per element
      const int nen_ = 8;
      if (la[1].Size() != nen_*numdofpernode_)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double> >robtempnp = rcp(new std::vector<double>(la[1].lm_.size()) );
      DRT::UTILS::ExtractMyValues(*tempnp,*robtempnp,la[1].lm_);
      params.set<Teuchos::RCP<vector<double> > >("robinson_tempnp",robtempnp);

      // calculate the THERMOmechanical term for fint
      soh8_finttemp(la,mydisp,myres,mytempnp,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    // call the purely structural method
    linstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,NULL,params,
      INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);
  }
  break;

  //==================================================================================
  // linear stiffness, internal force vector, and consistent mass matrix
  case calc_struct_nlnstiffmass:
  case calc_struct_nlnstifflmass:
  {
    // stiffness
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
    // massmatrix
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
    // elevec2+3 are not used anyway

    // need current displacement and residual forces
    Teuchos::RCP<const Epetra_Vector> disp
      = discretization.GetState(0, "displacement");
    Teuchos::RCP<const Epetra_Vector> res
      = discretization.GetState(0, "residual displacement");
    if (disp==Teuchos::null || res==Teuchos::null)
      dserror("Cannot get state vectors 'displacement' and/or residual");

    // build the location vector only for the structure field
    vector<int> lm = la[0].lm_;
    vector<double> mydisp((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm); // lm now contains only u-dofs
    vector<double> myres((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*res,myres,lm); // lm now contains only u-dofs

    // need current temperature state,
    // call the temperature discretization: thermo equates 2nd dofset
    // disassemble temperature
    if (discretization.HasState(1,"temperature"))
    {
      // call the temperature discretization in the location array
      std::vector<double> mytempnp((la[1].lm_).size());
      // check if you can get the temperature state
      Teuchos::RCP<const Epetra_Vector> tempnp
       = discretization.GetState(1,"temperature");
      if (tempnp==Teuchos::null)
        dserror("Cannot get state vector 'tempnp'");

      // the temperature field has only one dof per node, disregarded by the
      // dimension of the problem
      const int numdofpernode_ = 1;
      // number of nodes per element
      const int nen_ = 8;

      if (la[1].Size() != nen_*numdofpernode_)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double> >robtempnp = rcp(new std::vector<double>(la[1].lm_.size()) );
      DRT::UTILS::ExtractMyValues(*tempnp,*robtempnp,la[1].lm_);
      params.set<Teuchos::RCP<vector<double> > >("robinson_tempnp",robtempnp);

      // build the current temperature vector
      LINALG::Matrix<nen_*numdofpernode_,1> etemp(&(mytempnp[1]),true);  // view only!
      // calculate the THERMOmechanical term for fint
      soh8_finttemp(la,mydisp,myres,mytempnp,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    // call the purely structural method
    linstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,NULL,
      params,INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);

    if (act==calc_struct_nlnstifflmass) soh8_lumpmass(&elemat2);
  }
  break;

  //==================================================================================
  //predictor TangDis
  case calc_struct_reset_istep:
  {}
  break;

  //==================================================================================
  // evaluate stresses and strains at gauss points
  case calc_struct_stress:
  {
    // elemat1+2,elevec1-3 are not used anyway

    // nothing to do for ghost elements
    if (discretization.Comm().MyPID()==Owner())
    {
      // build the location vector only for the structure field
      vector<int> lm = la[0].lm_;
      Teuchos::RCP<const Epetra_Vector> disp
        = discretization.GetState(0,"displacement");
      Teuchos::RCP<const Epetra_Vector> res
        = discretization.GetState(0,"residual displacement");
      Teuchos::RCP<vector<char> > stressdata
        = params.get<Teuchos::RCP<vector<char> > >("stress", Teuchos::null);
      Teuchos::RCP<vector<char> > straindata
        = params.get<Teuchos::RCP<vector<char> > >("strain", Teuchos::null);
      // plastic strain data
      Teuchos::RCP<vector<char> > plstraindata
        = params.get<Teuchos::RCP<vector<char> > >("plstrain", Teuchos::null);
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      if (stressdata==Teuchos::null) dserror("Cannot get 'stress' data");
      if (straindata==Teuchos::null) dserror("Cannot get 'strain' data");
      if (plstraindata==Teuchos::null) dserror("Cannot get 'plastic strain' data");
      vector<double> mydisp((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> stress;

      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> strain;
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> plstrain;
      INPAR::STR::StressType iostress
        = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress",
            INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain
        = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain",
            INPAR::STR::strain_none);
      INPAR::STR::StrainType ioplstrain
        = DRT::INPUT::get<INPAR::STR::StrainType>(params, "ioplstrain",
            INPAR::STR::strain_none);

      // initialise the temperature-dependent stress
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> stresstemp(true);

      // need current temperature state,
      // call the temperature discretization: thermo equates 2nd dofset
      // disassemble temperature
      if (discretization.HasState(1,"temperature"))
      {
        // call the temperature discretization in the location array
        std::vector<double> mytempnp((la[1].lm_).size());
        // check if you can get the temperature state
        Teuchos::RCP<const Epetra_Vector> tempnp
          = discretization.GetState(1,"temperature");
        if (tempnp==Teuchos::null)
          dserror("Cannot get state vector 'tempnp'");

        // the temperature field has only one dof per node, disregarded by the
        // dimension of the problem
        const int numdofpernode_ = 1;
        // number of nodes per element
        const int nen_ = 8;

        if (la[1].Size() != nen_*numdofpernode_)
          dserror("Location vector length for temperature does not match!");
        // extract the current temperatures
        DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);

        // extract local values of the global vectors
        Teuchos::RCP<std::vector<double> >robtempnp = rcp(new std::vector<double>(la[1].lm_.size()) );
        DRT::UTILS::ExtractMyValues(*tempnp,*robtempnp,la[1].lm_);
        params.set<Teuchos::RCP<vector<double> > >("robinson_tempnp",robtempnp);

        // calculate the THERMOmechanical term for fint: temperature stresses
        soh8_finttemp(la,mydisp,myres,mytempnp,NULL,&stresstemp,NULL,params,
          iostress,INPAR::STR::strain_none);
      }

      // call the purely structural method
      linstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,
        &plstrain,params,iostress,iostrain,ioplstrain);

      // total stress
      // add stresstemp to the mechanical stress
      stress.Update(1.0,stresstemp,1.0);

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

      {
        DRT::PackBuffer data;
        AddtoPack(data, plstrain);
        data.StartPacking();
        AddtoPack(data, plstrain);
        std::copy(data().begin(),data().end(),std::back_inserter(*plstraindata));
      }

    }
  }
  break;

  //==================================================================================
  case calc_struct_update_istep:
  {
    // Update of history for visco material if they exist
    Teuchos::RCP<MAT::Material> mat = Material();
    if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
    {
      dserror("check if you wanna be here!!");
    }
    else if (mat->MaterialType() == INPAR::MAT::m_pllinelast)
    {
      MAT::PlasticLinElast* pllinelast = static_cast <MAT::PlasticLinElast*>(mat.get());
      pllinelast->Update();
    }
    else if (mat->MaterialType() == INPAR::MAT::m_thermopllinelast)
    {
      MAT::ThermoPlasticLinElast* thrpllinelast = static_cast <MAT::ThermoPlasticLinElast*>(mat.get());
      thrpllinelast->Update();
    }
    // incremental update of internal variables/history
    else if (mat->MaterialType() == INPAR::MAT::m_vp_robinson)
    {
      MAT::Robinson* robinson = static_cast <MAT::Robinson*>(mat.get());
      bool imrlike = false;
      robinson->Update(imrlike, 0.0);
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
    // elemat1+2,elevec1-3 are not used anyway

    const Teuchos::RCP<map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap
      = params.get<Teuchos::RCP<map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
    if (gpstressmap==Teuchos::null)
      dserror("no gp stress/strain map available for postprocessing");
    string stresstype = params.get<string>("stresstype","ndxyz");
    int gid = Id();
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> gpstress(((*gpstressmap)[gid])->A(),true);
    Teuchos::RCP<Epetra_MultiVector> poststress
      = params.get<RCP<Epetra_MultiVector> >("poststress",null);
    if (poststress==null)
      dserror("No element stress/strain vector available");

    if (stresstype=="ndxyz")
    {
      // extrapolate stresses/strains at Gauss points to nodes
      soh8_expol(gpstress, *poststress);
    }
    else if (stresstype=="cxyz")
    {
      const Epetra_BlockMap& elemap = poststress->Map();
      int lid = elemap.LID(Id());
      if (lid!=-1)
      {
        for (int i = 0; i < NUMSTR_SOH8; ++i)
        {
          double& s = (*((*poststress)(i)))[lid]; // resolve pointer for faster access
          s = 0.;
          for (int j = 0; j < NUMGPT_SOH8; ++j)
          {
            s += gpstress(j,i);
          }
          s *= 1.0/NUMGPT_SOH8;
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
  // linear stiffness, internal force vector, and consistent mass matrix
  case calc_struct_stifftemp:
  {
    // mechanical-thermal system matrix
    LINALG::Matrix<NUMDOF_SOH8,NUMNOD_SOH8> stiffmatrixcoupl(elemat1_epetra.A(),true);
    // elemat2,elevec1-3 are not used anyway

    // need current displacement and residual forces
    Teuchos::RCP<const Epetra_Vector> disp
      = discretization.GetState(0,"displacement");
    if (disp==null)
      dserror("Cannot get state vectors 'displacement'");
    vector<double> mydisp((la[0].lm_).size());
    // build the location vector only for the structure field
    DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);

    // calculate the mechanical-thermal system matrix
    soh8_stifftemp(la,mydisp,&stiffmatrixcoupl);
  }
  break;

  //==================================================================================
  default:
    dserror("Unknown type of action for So_hex8");

  } // action

  cout << "so_hex8_thermo Ende LinEvaluate\n "  << endl;

  return 0;

} // Evaluate


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                           dano 05/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::linstiffmass(
  vector<int>& lm,  // location matrix
  vector<double>& disp,  // current displacements
  vector<double>& residual,  // current residual displacements or displacement increment
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,  // element stiffness matrix
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  // element mass matrix
  LINALG::Matrix<NUMDOF_SOH8,1>* force,  // element internal force vector
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,  // stresses at GP
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,  // strains at GP
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* eleplstrain, // plastic strains at GP
  ParameterList& params,  // algorithmic parameters e.g. time
  const INPAR::STR::StressType iostress,  // stress output option
  const INPAR::STR::StrainType iostrain,  // strain output option
  const INPAR::STR::StrainType ioplstrain  // plastic strain output option
  )
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  const static vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xcurr;  // current  coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xdisp;

  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH8+2];
  }

  LINALG::Matrix<NUMDOF_SOH8,1> nodaldisp;
  // in case of Robinson's material, the (residual) displacements are required
  // residual displacements correspond to current displacement increment
  LINALG::Matrix<NUMDOF_SOH8,1> res_d;
  for (int i=0; i<NUMDOF_SOH8; ++i)
  {
    nodaldisp(i,0) = disp[i];
    res_d(i) = residual[i];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(true);
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],derivs[gp]);
    double detJ = detJ_[gp];

    // WATCH OUT: here is the difference to the non-linear method nlnstiffmass()
    // in geometrically linear analysis the deformation gradient is equal to identity
    // no difference between reference and current state
    for (int i=0; i<3; ++i) defgrd(i,i) = 1.0;

    // linear B-operator B = N_XYZ
    // disperse global derivatives to bop-lines
    // bop is arranged as usual (refer to script FE or elsewhere):
    //
    // [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
    // [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
    // [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
    // [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
    // [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
    // [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ]
    LINALG::Matrix<NUMSTR_SOH8,NUMDOF_SOH8> boplin;
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
      boplin(0,NODDOF_SOH8*i+0) = N_XYZ(0,i);
      boplin(0,NODDOF_SOH8*i+1) = 0.0;
      boplin(0,NODDOF_SOH8*i+2) = 0.0;
      boplin(1,NODDOF_SOH8*i+0) = 0.0;
      boplin(1,NODDOF_SOH8*i+1) = N_XYZ(1,i);
      boplin(1,NODDOF_SOH8*i+2) = 0.0;
      boplin(2,NODDOF_SOH8*i+0) = 0.0;
      boplin(2,NODDOF_SOH8*i+1) = 0.0;
      boplin(2,NODDOF_SOH8*i+2) = N_XYZ(2,i);
      /* ~~~ */
      boplin(3,NODDOF_SOH8*i+0) = N_XYZ(1,i);
      boplin(3,NODDOF_SOH8*i+1) = N_XYZ(0,i);
      boplin(3,NODDOF_SOH8*i+2) = 0.0;
      boplin(4,NODDOF_SOH8*i+0) = 0.0;
      boplin(4,NODDOF_SOH8*i+1) = N_XYZ(2,i);
      boplin(4,NODDOF_SOH8*i+2) = N_XYZ(1,i);
      boplin(5,NODDOF_SOH8*i+0) = N_XYZ(2,i);
      boplin(5,NODDOF_SOH8*i+1) = 0.0;
      boplin(5,NODDOF_SOH8*i+2) = N_XYZ(0,i);
    }

    // approximate linearised strain tensor using common naming of strain vector
    // glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(NUMSTR_SOH8);
    LINALG::Matrix<NUMSTR_SOH8,1> glstrain(glstrain_epetra.A(),true);
    // E = epsilon_GL == epsilon_1
    // build the linearised strain epsilon = B . d
    glstrain.Multiply(boplin,nodaldisp);

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
    case INPAR::STR::strain_gl:
    {
      if (elestrain == NULL) dserror("strain data not available");
      for (int i = 0; i < 3; ++i)
        (*elestrain)(gp,i) = glstrain(i);
      for (int i = 3; i < 6; ++i)
        (*elestrain)(gp,i) = 0.5 * glstrain(i);
    }
    break;
    case INPAR::STR::strain_ea:
    {
      if (elestrain == NULL) dserror("strain data not available");
      // rewriting Green-Lagrange strains in matrix format
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> gl;
      gl(0,0) = glstrain(0);
      gl(0,1) = 0.5*glstrain(3);
      gl(0,2) = 0.5*glstrain(5);
      gl(1,0) = gl(0,1);
      gl(1,1) = glstrain(1);
      gl(1,2) = 0.5*glstrain(4);
      gl(2,0) = gl(0,2);
      gl(2,1) = gl(1,2);
      gl(2,2) = glstrain(2);

      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> euler_almansi;
      // in geometrically linear analysis: invdefgrd == defgrd
      temp.Multiply(gl,defgrd);
      euler_almansi.MultiplyTN(defgrd,temp);

      (*elestrain)(gp,0) = euler_almansi(0,0);
      (*elestrain)(gp,1) = euler_almansi(1,1);
      (*elestrain)(gp,2) = euler_almansi(2,2);
      (*elestrain)(gp,3) = euler_almansi(0,1);
      (*elestrain)(gp,4) = euler_almansi(1,2);
      (*elestrain)(gp,5) = euler_almansi(0,2);
    }
    break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror("requested strain type not available");
    }

    // build incremental strains
    // Delta strain = B . Delta disp
    LINALG::Matrix<NUMSTR_SOH8,1> straininc(true);
    straininc.Multiply(boplin,res_d);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    double scalartemp = 0.0;
    LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(true);
    LINALG::Matrix<NUMSTR_SOH8,1> stress(true);
    LINALG::Matrix<NUMSTR_SOH8,1> plglstrain(true);

    // default: material call in structural function is purely deformation dependent
    if ( Material()->MaterialType() != INPAR::MAT::m_vp_robinson )
      soh8_mat_sel(&stress,&cmat,&density,&glstrain,&plglstrain,&defgrd,gp,params);
    // if Robinson's material --> pass the current temperature to the material
    else if ( (Material()->MaterialType() == INPAR::MAT::m_vp_robinson) )
    {
      // scalar-valued temperature: T = shapefunctions . element temperatures
      // T = N_T^(e) . T^(e)
      // get the temperature vector by extraction from parameter list
      LINALG::Matrix<NUMNOD_SOH8,1> etemp(true);
      LINALG::Matrix<1,1> Ntemp(false);
      LINALG::Matrix<NUMSTR_SOH8,1> ctemp(true);

      Teuchos::RCP<vector<double> > temperature_vector
        = params.get<Teuchos::RCP<vector<double> > >("robinson_tempnp",Teuchos::null);
      // in StructureBaseAlgorithm() temperature not yet available, i.e. ==null
      if (temperature_vector==Teuchos::null)
      {
        MAT::Robinson* robinson
          = static_cast <MAT::Robinson*>(Material().get());
        // initialise the temperature field
        scalartemp = robinson->InitTemp();
      }
      // temperature vector is available
      else  // (temperature_vector!=Teuchos::null)
      {
        for (int i=0; i<NUMNOD_SOH8; ++i)
        {
          etemp(i,0) = (*temperature_vector)[i+0];
        }
        // copy structural shape functions needed for the thermo field
        // identical shapefunctions for the displacements and the temperatures
        scalartemp  = (shapefcts[gp]).Dot(etemp);
      }
      soh8_mat_temp(&stress,&ctemp,NULL,&cmat,&defgrd,&glstrain,&plglstrain,straininc,scalartemp,&density,gp,params);

    } // end Robinson's material

    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp plastic strains (only in case of plastic strain output)
    switch (ioplstrain)
    {
    case INPAR::STR::strain_gl:
    {
     if (eleplstrain == NULL) dserror("plastic strain data not available");
     for (int i = 0; i < 3; ++i)
       (*eleplstrain)(gp,i) = plglstrain(i);
     for (int i = 3; i < 6; ++i)
       (*eleplstrain)(gp,i) = 0.5 * plglstrain(i);
    }
    break;
    case INPAR::STR::strain_ea:
    {
     if (eleplstrain == NULL) dserror("plastic strain data not available");
     // rewriting Green-Lagrange strains in matrix format
     LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> gl;
     gl(0,0) = plglstrain(0);
     gl(0,1) = 0.5*plglstrain(3);
     gl(0,2) = 0.5*plglstrain(5);
     gl(1,0) = gl(0,1);
     gl(1,1) = plglstrain(1);
     gl(1,2) = 0.5*plglstrain(4);
     gl(2,0) = gl(0,2);
     gl(2,1) = gl(1,2);
     gl(2,2) = plglstrain(2);

     LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
     LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> euler_almansi;
     // in geometrically linear analysis: invdefgrd == defgrd
     temp.Multiply(gl,defgrd);
     euler_almansi.MultiplyTN(defgrd,temp);

     (*eleplstrain)(gp,0) = euler_almansi(0,0);
     (*eleplstrain)(gp,1) = euler_almansi(1,1);
     (*eleplstrain)(gp,2) = euler_almansi(2,2);
     (*eleplstrain)(gp,3) = euler_almansi(0,1);
     (*eleplstrain)(gp,4) = euler_almansi(1,2);
     (*eleplstrain)(gp,5) = euler_almansi(0,2);
    }
    break;
    case INPAR::STR::strain_none:
     break;

    default:
     dserror("requested plastic strain type not available");
    }

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < NUMSTR_SOH8; ++i)
        (*elestress)(gp,i) = stress(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");

      LINALG::Matrix<3,3> pkstress;
      pkstress(0,0) = stress(0);
      pkstress(0,1) = stress(3);
      pkstress(0,2) = stress(5);
      pkstress(1,0) = pkstress(0,1);
      pkstress(1,1) = stress(1);
      pkstress(1,2) = stress(4);
      pkstress(2,0) = pkstress(0,2);
      pkstress(2,1) = pkstress(1,2);
      pkstress(2,2) = stress(2);

      LINALG::Matrix<3,3> temp;
      LINALG::Matrix<3,3> cauchystress;
      temp.Multiply(1.0,defgrd,pkstress);
      cauchystress.MultiplyNT(temp,defgrd);

      (*elestress)(gp,0) = cauchystress(0,0);
      (*elestress)(gp,1) = cauchystress(1,1);
      (*elestress)(gp,2) = cauchystress(2,2);
      (*elestress)(gp,3) = cauchystress(0,1);
      (*elestress)(gp,4) = cauchystress(1,2);
      (*elestress)(gp,5) = cauchystress(0,2);
    }
    break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror("requested stress type not available");
    }

    double detJ_w = detJ*gpweights[gp];

    // update/integrate internal force vector
    if (force != NULL)
    {
      // f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, boplin, stress, 1.0);
    }

    // update/integrate `elastic' and `initial-displacement' stiffness matrix
    if (stiffmatrix != NULL)
    {
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,NUMDOF_SOH8> cb;
      cb.Multiply(cmat,boplin);
      stiffmatrix->MultiplyTN(detJ_w,boplin,cb,1.0);
    }

    if (massmatrix != NULL) // evaluate mass matrix +++++++++++++++++++++++++
    {
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOH8; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += massfactor;
        }
      }

    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
} // DRT::ELEMENTS::So_hex8::linstiffmass


/*----------------------------------------------------------------------*
 | evaluate only the temperature fraction for the element    dano 05/10 |
 | originally by maf 04/07  (private)                                   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_finttemp(
  DRT::Element::LocationArray& la,  // location array
  vector<double>& disp,  // current displacements
  vector<double>& residual,  // current residual displ
  vector<double>& temp, // current temperature
  LINALG::Matrix<NUMDOF_SOH8,1>* force,  // element internal force vector
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,  // stresses at GP
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,  // strains at GP
  Teuchos::ParameterList& params,  // algorithmic parameters e.g. time
  const INPAR::STR::StressType iostress,  // stress output option
  const INPAR::STR::StrainType iostrain  // strain output option
  )
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  const static vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry (8x3)
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // X, material coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xcurr;  // x, current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH8+2];
  }

  // we need the (residual) displacement -- current increment of displacement
  LINALG::Matrix<NUMDOF_SOH8,1> res_d;
  for (int i = 0; i < NUMDOF_SOH8; ++i)
  {
    res_d(i) = residual[i];
  }

  // vector of the current element temperatures
  LINALG::Matrix<NUMNOD_SOH8,1> etemp(true);
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    etemp(i,0) = temp[i+0];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(true);

  // identical shapefunctions for the displacements and the temperatures
  LINALG::Matrix<NUMNOD_SOH8,1> shapetemp;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],derivs[gp]); // (6.21)
    double detJ = detJ_[gp]; // (6.22)

    // set to initial state (defgrd == identity)
    for (int i=0; i<3; ++i)
      defgrd(i,i) = 1.0;

    // linear B-operator B = N_XYZ
    // disperse global derivatives to bop-lines
    // bop is arranged as usual (refer to script FE or elsewhere):
    // [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
    // [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
    // [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
    // [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
    // [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
    // [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ]
    LINALG::Matrix<NUMSTR_SOH8,NUMDOF_SOH8> boplin;
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
      boplin(0,NODDOF_SOH8*i+0) = N_XYZ(0,i);
      boplin(0,NODDOF_SOH8*i+1) = 0.0;
      boplin(0,NODDOF_SOH8*i+2) = 0.0;
      boplin(1,NODDOF_SOH8*i+0) = 0.0;
      boplin(1,NODDOF_SOH8*i+1) = N_XYZ(1,i);
      boplin(1,NODDOF_SOH8*i+2) = 0.0;
      boplin(2,NODDOF_SOH8*i+0) = 0.0;
      boplin(2,NODDOF_SOH8*i+1) = 0.0;
      boplin(2,NODDOF_SOH8*i+2) = N_XYZ(2,i);
      /* ~~~ */
      boplin(3,NODDOF_SOH8*i+0) = N_XYZ(1,i);
      boplin(3,NODDOF_SOH8*i+1) = N_XYZ(0,i);
      boplin(3,NODDOF_SOH8*i+2) = 0.0;
      boplin(4,NODDOF_SOH8*i+0) = 0.0;
      boplin(4,NODDOF_SOH8*i+1) = N_XYZ(2,i);
      boplin(4,NODDOF_SOH8*i+2) = N_XYZ(1,i);
      boplin(5,NODDOF_SOH8*i+0) = N_XYZ(2,i);
      boplin(5,NODDOF_SOH8*i+1) = 0.0;
      boplin(5,NODDOF_SOH8*i+2) = N_XYZ(0,i);
    }
    // copy structural shape functions needed for the thermo field
    shapetemp.Update(shapefcts[gp]);

    // product of shapefunctions and element temperatures for stresstemp
    // N_T . T
    LINALG::Matrix<1,1> Ntemp(true);
    Ntemp.MultiplyTN(shapetemp,etemp);
    const double scalartemp  = shapetemp.Dot(etemp);

    // build iterative strains
    LINALG::Matrix<NUMSTR_SOH8,1> straininc(true);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    // calculate the stress part dependent on the temperature in the material
    LINALG::Matrix<NUMSTR_SOH8,1> ctemp(true);
    LINALG::Matrix<NUMSTR_SOH8,1> stresstemp(true);
    LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(true);
    LINALG::Matrix<NUMSTR_SOH8,1> glstrain(true);
    LINALG::Matrix<NUMSTR_SOH8,1> plglstrain(true);
    // take care: current temperature ( N . T ) is passed to the element
    //            in the material: 1.) Delta T = subtract ( N . T - T_0 )
    //                             2.) stresstemp = C . Delta T
    // do not call the material for Robinson's material
    if ( !(Material()->MaterialType() == INPAR::MAT::m_vp_robinson) )
      soh8_mat_temp(&stresstemp,&ctemp,&Ntemp,&cmat,&defgrd,&glstrain,&plglstrain,
        straininc,scalartemp,&density,gp,params);

    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < NUMSTR_SOH8; ++i)
        (*elestress)(gp,i) = stresstemp(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");

      LINALG::Matrix<3,3> pkstresstemp;
      pkstresstemp(0,0) = stresstemp(0);
      pkstresstemp(0,1) = stresstemp(3);
      pkstresstemp(0,2) = stresstemp(5);
      pkstresstemp(1,0) = pkstresstemp(0,1);
      pkstresstemp(1,1) = stresstemp(1);
      pkstresstemp(1,2) = stresstemp(4);
      pkstresstemp(2,0) = pkstresstemp(0,2);
      pkstresstemp(2,1) = pkstresstemp(1,2);
      pkstresstemp(2,2) = stresstemp(2);

      LINALG::Matrix<3,3> temp;
      LINALG::Matrix<3,3> cauchystresstemp;
      temp.Multiply(1.0,defgrd,pkstresstemp);
      cauchystresstemp.MultiplyNT(temp,defgrd);

      (*elestress)(gp,0) = cauchystresstemp(0,0);
      (*elestress)(gp,1) = cauchystresstemp(1,1);
      (*elestress)(gp,2) = cauchystresstemp(2,2);
      (*elestress)(gp,3) = cauchystresstemp(0,1);
      (*elestress)(gp,4) = cauchystresstemp(1,2);
      (*elestress)(gp,5) = cauchystresstemp(0,2);
    }
    break;
    case INPAR::STR::stress_none:
      break;

    default:
      dserror("requested stress type not available");
    }

    // integrate internal force vector
    // f = f + (B^T . sigma_temp) * detJ * w(gp)
    if (force != NULL)
    {
      double detJ_w = detJ*gpweights[gp];
      force->MultiplyTN(detJ_w, boplin, stresstemp, 1.0);
    }  // if (force != NULL)
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
} // DRT::ELEMENTS::So_hex8::soh8_finttemp


/*----------------------------------------------------------------------*
 | evaluate only the mechanical-thermal stiffness term       dano 03/11 |
 | for monolithic TSI (private)                                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_stifftemp(
  DRT::Element::LocationArray& la,
  vector<double>& disp,
  // element mechanical-thermal stiffness matrix
  LINALG::Matrix<NUMDOF_SOH8,NUMNOD_SOH8>* stiffmatrixcoupl // (nsd_*nen_ x nen_)
  )
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  const static vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry (8x3)
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // X, material coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xcurr;  // x, current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH8+2];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(true);

  // identical shapefunctions for the displacements and the temperatures
  LINALG::Matrix<NUMNOD_SOH8,1> shapetemp;

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],derivs[gp]); // (6.21)
    double detJ = detJ_[gp]; // (6.22)

    // set to initial state (defgrd == identity)
    for (int i=0; i<3; ++i) defgrd(i,i) = 1.0;

    // linear B-operator B = N_XYZ
    // disperse global derivatives to bop-lines
    // bop is arranged as usual (refer to script FE or elsewhere):
    //
    // [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
    // [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
    // [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
    // [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
    // [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
    // [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ]
    LINALG::Matrix<NUMSTR_SOH8,NUMDOF_SOH8> boplin;
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
      boplin(0,NODDOF_SOH8*i+0) = N_XYZ(0,i);
      boplin(0,NODDOF_SOH8*i+1) = 0.0;
      boplin(0,NODDOF_SOH8*i+2) = 0.0;
      boplin(1,NODDOF_SOH8*i+0) = 0.0;
      boplin(1,NODDOF_SOH8*i+1) = N_XYZ(1,i);
      boplin(1,NODDOF_SOH8*i+2) = 0.0;
      boplin(2,NODDOF_SOH8*i+0) = 0.0;
      boplin(2,NODDOF_SOH8*i+1) = 0.0;
      boplin(2,NODDOF_SOH8*i+2) = N_XYZ(2,i);
      /* ~~~ */
      boplin(3,NODDOF_SOH8*i+0) = N_XYZ(1,i);
      boplin(3,NODDOF_SOH8*i+1) = N_XYZ(0,i);
      boplin(3,NODDOF_SOH8*i+2) = 0.0;
      boplin(4,NODDOF_SOH8*i+0) = 0.0;
      boplin(4,NODDOF_SOH8*i+1) = N_XYZ(2,i);
      boplin(4,NODDOF_SOH8*i+2) = N_XYZ(1,i);
      boplin(5,NODDOF_SOH8*i+0) = N_XYZ(2,i);
      boplin(5,NODDOF_SOH8*i+1) = 0.0;
      boplin(5,NODDOF_SOH8*i+2) = N_XYZ(0,i);
    }

    // copy structural shape functions needed for the thermo field
    shapetemp.Update(1.0,shapefcts[gp],0.0); // (8x1)

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated
    */
    // get the thermal material tangent
    LINALG::Matrix<NUMSTR_SOH8,1> ctemp(true);
    Ctemp(&ctemp);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    double detJ_w = detJ*gpweights[gp];
    // update coupling matrix K_dT
    if (stiffmatrixcoupl != NULL)
    {
      // C_temp . N_temp
      LINALG::Matrix<NUMSTR_SOH8,NUMNOD_SOH8> cn(true);
      cn.MultiplyNT(ctemp,shapetemp); // (6x8)=(6x1)(1x8)
      // integrate stiffness term
      // k_st = k_st + (B^T . C_temp . N_temp) * detJ * w(gp)
      stiffmatrixcoupl->MultiplyTN(detJ_w, boplin, cn, 1.0);
    }
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
} // DRT::ELEMENTS::So_hex8::soh8_stifftemp

/*----------------------------------------------------------------------*/

