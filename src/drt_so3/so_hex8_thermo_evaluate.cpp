/*!----------------------------------------------------------------------
\file so_hex8_thermo_evaluate.cpp
\brief

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*----------------------------------------------------------------------*/

#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/thermoplasticlinelast.H"
#include "../drt_mat/robinson.H"
#include "../drt_mat/micromaterial.H"
#include <iterator>

#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | evaluate the element (public)                             dano 02/10 |
 | originally by maf 04/07                                              |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::Evaluate(
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
  // call the simple Evaluate(lm) if there is only the structure discretization
  if (la.Size()==1)
  {
    // "normal" Evaluate
    return Evaluate(
      params,
      discretization,
      la[0].lm_, // location vector is build by the first column of la
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );
  }

  // TSI: volume coupling stuff

  // type of kinematic of the problem
  // solve a geometric linear system
  if (kintype_ == DRT::ELEMENTS::So_hex8::soh8_linear)
  {
    return LinEvaluate(
      params,
      discretization,
      la,
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );
  }
  // else: geometrically non-linear with Total Lagrangean approach

  // start with "none"
  DRT::ELEMENTS::So_hex8::ActionType act = So_hex8::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_nlnstiff")              act = So_hex8::calc_struct_nlnstiff;
  else if (action=="calc_struct_nlnstiffmass")          act = So_hex8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_internalforce")         act = So_hex8::calc_struct_internalforce;
  else if (action=="calc_struct_stress")                act = So_hex8::calc_struct_stress;
  else if (action=="calc_struct_update_istep")          act = So_hex8::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")           act = So_hex8::calc_struct_reset_istep;  // needed for TangDis predictor
  else if (action=="postprocess_stress")                act = So_hex8::postprocess_stress;
  else if (action=="calc_struct_stifftemp")             act = So_hex8::calc_struct_stifftemp;
  else dserror("Unknown type of action for So_hex8: %s",action.c_str());
  // what should the element do
  switch(act)
  {
  //==================================================================================
  // nonlinear stiffness and internal force vector
  case calc_struct_nlnstiff:
  {
    // stiffness
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
    // elemat2, elevec2+3 is not used anyway

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
    if (elemat1.IsInitialized()) matptr = &elemat1;
    // call the well-known soh8_nlnstiffmass for the normal structure solution
    soh8_nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,NULL,params,
      INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);

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

      // the coupling stiffness matrix (3nx1)
      LINALG::Matrix<NUMDOF_SOH8,1> elemat1(elemat1_epetra.A(),true);

      // calculate the THERMOmechanical solution
      soh8_nlnstifftemp(la,mydisp,myres,mytempnp,&elemat1,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);

      // calculate the THERMOmechanical term for fint
      soh8_finttemp(la,mydisp,myres,mytempnp,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
  }
  break;

  //==================================================================================
  // nonlinear stiffness, internal force vector, and consistent mass matrix
  case calc_struct_nlnstiffmass:
  {
    // stiffness
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
    // massmatrix
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
    // internal force vector
    LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
    // eleve2+elevec3 is not used anyway

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
    // call the well-known soh8_nlnstiffmass for the normal structure solution
    soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,NULL,
      params,INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);

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
      // build the current temperature vector
      LINALG::Matrix<nen_*numdofpernode_,1> etemp(&(mytempnp[1]),true);  // view only!
      // the coupling stiffness matrix (3nx1)
      LINALG::Matrix<NUMDOF_SOH8,1> elemat1(elemat1_epetra.A(),true);
      // calculate the THERMOmechanical solution
      soh8_nlnstifftemp(la,mydisp,myres,mytempnp,&elemat1,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);

      // calculate the THERMOmechanical term for fint
      soh8_finttemp(la,mydisp,myres,mytempnp,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);

    }
  }
  break;

  //==================================================================================
  // internal force vector only
  case calc_struct_internalforce:
  {
    // internal force vector
    LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);

    // need current displacement and residual forces
    RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
    if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
    // build the location vector only for the structure field
    vector<int> lm = la[0].lm_;
    vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    vector<double> myres(lm.size());
    DRT::UTILS::ExtractMyValues(*res,myres,lm);
    // create a dummy element matrix to apply linearised EAS-stuff onto
    LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> myemat(true);

    // default: geometrically non-linear analysis with Total Lagrangean approach
      soh8_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,NULL,params,
                      INPAR::STR::stress_none,INPAR::STR::strain_none,INPAR::STR::strain_none);
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
  // allowing the predictor TangDis in .dat --> can be decisive in compressible case!
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
        = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain
        = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
      INPAR::STR::StrainType ioplstrain
        = DRT::INPUT::get<INPAR::STR::StrainType>(params, "ioplstrain", INPAR::STR::strain_none);
      // call the well-known soh8_nlnstiffmass for the normal structure solution
      soh8_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,&plstrain,
        params,iostress,iostrain,ioplstrain);

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
        // get the temperature dependent stress
        LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> stresstemp;
        // calculate the THERMOmechanical solution: temperature stresses
        soh8_nlnstifftemp(la,mydisp,myres,mytempnp,NULL,NULL,&stresstemp,NULL,
          params,iostress,INPAR::STR::strain_none);

        // calculate the THERMOmechanical term for fint: temperature stresses
        soh8_finttemp(la,mydisp,myres,mytempnp,NULL,&stresstemp,NULL,params,
          iostress,INPAR::STR::strain_none);

        // total stress
        // add stresstemp to the mechanical stress
        stress.Update(1.0,stresstemp,1.0);
      }

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
    RefCountPtr<MAT::Material> mat = Material();
    if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
      micro->Update();
    }
    // incremental update of internal variables/history
    if (mat->MaterialType() == INPAR::MAT::m_vp_robinson)
    {
      MAT::Robinson* robinson = static_cast<MAT::Robinson*>(mat.get());
      robinson->Update();
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

    // nothing to do for ghost elements
    if (discretization.Comm().MyPID()==Owner())
    {
      const Teuchos::RCP<map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap
        = params.get<Teuchos::RCP<map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
      if (gpstressmap==Teuchos::null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> gpstress(((*gpstressmap)[gid])->A(),true);
      Teuchos::RCP<Epetra_MultiVector> poststress
        = params.get<RCP<Epetra_MultiVector> >("poststress",Teuchos::null);
      if (poststress==Teuchos::null)
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
  }
  break;

  //==================================================================================
  default:
  dserror("Unknown type of action for So_hex8");
  break;
  } // action

  return 0;
} // Evaluate


/*----------------------------------------------------------------------*
 | evaluate only the temperature fraction for the element    dano 03/10 |
 | originally by maf 04/07  (private)                                   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_nlnstifftemp(
  DRT::Element::LocationArray& la,  // location array
  vector<double>& disp,  // current displacements
  vector<double>& residual,  // current residual displ
  vector<double>& temp, // current temperature
  LINALG::Matrix<NUMDOF_SOH8,1>* tempstiffmatrix, // coupling stiffness matrix
  LINALG::Matrix<NUMDOF_SOH8,1>* force,  // element internal force vector
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,  // stresses at GP
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,  // strains at GP
  Teuchos::ParameterList& params,  // algorithmic parameters e.g. time
  const INPAR::STR::StressType iostress,  // stress output option
  const INPAR::STR::StrainType iostrain  // strain output option
  )
{

  dserror("TSI with total Lagrangean approach not yet available!");

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
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(false);

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

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,N_XYZ); //  (6.17)

    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
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
    LINALG::Matrix<NUMSTR_SOH8,NUMDOF_SOH8> bop;
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
      bop(0,NODDOF_SOH8*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH8*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH8*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NODDOF_SOH8*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH8*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH8*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NODDOF_SOH8*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH8*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH8*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NODDOF_SOH8*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH8*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH8*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NODDOF_SOH8*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH8*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH8*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NODDOF_SOH8*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH8*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH8*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    // copy structural shape functions needed for the thermo field
    shapetemp.Update(shapefcts[gp]);

    // product of shapefunctions and element temperatures for stresstemp
    LINALG::Matrix<1,1> Ntemp(true);
    Ntemp.MultiplyTN(shapetemp,etemp);
    double scalartemp  = shapetemp.Dot(etemp);

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
    LINALG::Matrix<NUMSTR_SOH8,1> straininc(true);
    // take care: current temperature ( N . T ) is passed to the element
    //            in the material: 1.) Delta T = subtract ( N . T - T_0 )
    //                             2.) stresstemp = C . Delta T
    soh8_mat_temp(&stresstemp,&ctemp,&Ntemp,&cmat,&defgrd,&glstrain,&plglstrain,straininc,scalartemp,&density,gp,params);

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
      const double detF = defgrd.Determinant();

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
      temp.Multiply(1.0/detF,defgrd,pkstresstemp);
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
      break;
    }

    double detJ_w = detJ*gpweights[gp];
    // update internal force vector
    if (force != NULL)
    {
      // integrate internal force vector
      // f = f + (B^T . sigma_temp) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stresstemp, 1.0);
    }

    // update stiffness matrix
    if (tempstiffmatrix != NULL)
    {
      //-----------------------------------------------------------------------
      // integrate `initial-displacement-temperature' stiffness matrix  02/10
      // update the stiffness part depending on the temperature
      // k_theta += 1.0 * k_theta + detJ * w(gp) * (B^T . C_Theta . N_Theta)
      //-----------------------------------------------------------------------
      LINALG::Matrix<NUMSTR_SOH8,1> cn; // 6x1 = 6x1n
      // extract the shapefunctions on the gp_i and multiply it with cn
      for (int stresscomp=0; stresscomp<NUMSTR_SOH8; ++stresscomp)
      {
        cn.Multiply(ctemp,Ntemp); // (6x1) = (6x1)(1x1)
      }
      tempstiffmatrix->MultiplyTN(detJ_w,bop,cn,1.0); // [(3nx6)(6x1)=(3nx1)] (24x6)(6x1)=(24x1)

      //-----------------------------------------------------------------------
      // integrate `geometric' stiffness matrix and add to keu            04/10
      // kgeo += (B_L^T . sigma_Theta . B_L^T) * detJ * w(gp)  with B_L = Ni,Xj
      //-----------------------------------------------------------------------
      // \f {\mathbf S}:\Delta\delta{\mathbf E} \f
      // \f = \delta {\mathbf d}^T \,{\mathbf k}_{geo}\, \Delta{\mathbf d} \f
      // auxiliary integrated temperature stress
      LINALG::Matrix<6,1> sfac(stresstemp);

      // detJ * w(gp) * [S11,S22,S33,S12=S21,S23=S32,S13=S31]
      sfac.Scale(detJ_w);

      // intermediate Sm.B_L
      vector<double> SmB_L(3);

      // like in so_hex8_evaluate no changes for tempstress
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)
      // with B_L = Ni,Xj: LINEAR B-operator see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOH8; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
                     + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
                     + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
                     + sfac(2) * N_XYZ(2, inod);
        SmB_L[0] = sfac(3) * N_XYZ(1, inod) + sfac(0) * N_XYZ(0, inod)
                     + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
                     + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
                     + sfac(2) * N_XYZ(2, inod);

        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod)
        {
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<NUMDIM_SOH8; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          // <3n,1>
          // 3n: displacement:3(3FHG/node)*8(8nodes/element)+3(dim),
          //  1: temperature is a scalar
          (*tempstiffmatrix)(3*inod+0,0) += bopstrbop;
          (*tempstiffmatrix)(3*inod+1,0) += bopstrbop;
          (*tempstiffmatrix)(3*inod+2,0) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************
    }
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
} // DRT::ELEMENTS::So_hex8::soh8_nlnstifftemp


/*----------------------------------------------------------------------*
 | material law with temperature part for So_hex8            dano 05/10 |
 | originally by gee in so_material.cpp 10/08                           |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_mat_temp(
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stresstemp,
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>* ctemp,
  LINALG::Matrix<1,1>* Ntemp,  // temperature of element
  LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
  LINALG::Matrix<3,3>* defgrd, //
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>* plglstrain,
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>& straininc,
  const double& scalartemp,
  double* density,
  const int gp,
  Teuchos::ParameterList& params
  )
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stresstemp) dserror("No stress vector supplied");
  if (!ctemp) dserror("No material tangent matrix supplied");
  if (!Ntemp) dserror("No temperature supplied");
  if (!defgrd) dserror("No defgrd supplied");
#endif

  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here.
  // the old interface does not exist anymore...
  Teuchos::RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    // st.venant-kirchhoff-material with temperature
    case INPAR::MAT::m_thermostvenant:
    {
      MAT::ThermoStVenantKirchhoff* thrstvk
        = static_cast <MAT::ThermoStVenantKirchhoff*>(mat.get());
      thrstvk->Evaluate(*Ntemp,*ctemp,*stresstemp);
      *density = thrstvk->Density();
      return;
      break;
    }
    // small strain von Mises thermoelastoplastic material
    case INPAR::MAT::m_thermopllinelast:
    {
      MAT::ThermoPlasticLinElast* thrpllinelast
        = static_cast <MAT::ThermoPlasticLinElast*>(mat.get());
      thrpllinelast->Evaluate(*Ntemp,*ctemp,*stresstemp);
      *density = thrpllinelast->Density();
      return;
      break;
    }
    case INPAR::MAT::m_vp_robinson: /*-- visco-plastic Robinson's material */
    {
      MAT::Robinson* robinson = static_cast <MAT::Robinson*>(mat.get());
      robinson->Evaluate(*glstrain,*plglstrain,straininc,scalartemp,gp,params,*cmat,*stresstemp);
      *density = robinson->Density();
      return;
      break;
    }
    default:
      dserror("Unknown type of temperature dependent material");
    break;
  } // switch (mat->MaterialType())

  return;
} // of soh8_mat_temp()


/*----------------------------------------------------------------------*
 | get the constant temperature fraction for stresstemp      dano 05/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Ctemp(LINALG::Matrix<6,1>* ctemp)
{
  Teuchos::RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    /*-------- thermo st.venant-kirchhoff-material */
    case INPAR::MAT::m_thermostvenant:
    {
      MAT::ThermoStVenantKirchhoff* thrstvk
        = static_cast<MAT::ThermoStVenantKirchhoff*>(mat.get());
       return thrstvk->SetupCthermo(*ctemp);
       break;
    }
    // small strain von Mises thermoelastoplastic material
    case INPAR::MAT::m_thermopllinelast:
    {
      MAT::ThermoPlasticLinElast* thrpllinelast
        = static_cast <MAT::ThermoPlasticLinElast*>(mat.get());
      return thrpllinelast->SetupCthermo(*ctemp);
      break;
    }
    default:
      dserror("Cannot ask material for the temperature rhs");
      break;

  } // switch (mat->MaterialType())

}  // Ctemp()

/*----------------------------------------------------------------------*/

