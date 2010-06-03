/*!----------------------------------------------------------------------
\file so_hex8_thermo_LIN_evaluate.cpp
\brief

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include <iterator>

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

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

  // TSI: volume coupling stuff

  // stiffness
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
  // massmatrix
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
  // internal force vector
  LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH8,1> elevec2(elevec2_epetra.A(),true);
  // elevec3 is not used anyway

  // start with "none"
  DRT::ELEMENTS::So_hex8::ActionType act = So_hex8::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_nlnstiff")      act = So_hex8::calc_struct_nlnstiff;
  else if (action=="calc_struct_nlnstiffmass")  act = So_hex8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = So_hex8::calc_struct_stress;
  else if (action=="calc_struct_update_istep")  act = So_hex8::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep")   act = So_hex8::calc_struct_reset_istep;  // needed for TangDis predictor
  else if (action=="postprocess_stress")        act = So_hex8::postprocess_stress;
  else dserror("Unknown type of action for So_hex8");
  // what should the element do
  switch(act)
  {
  //==================================================================================
  // linear stiffness and internal force vector
  case calc_struct_nlnstiff:
  {
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
    soh8_linstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
                       INPAR::STR::stress_none,INPAR::STR::strain_none);

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
      const int iel = 8;
      if (la[1].Size() != iel*numdofpernode_)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);

      // the coupling stiffness matrix (3nx1)
      LINALG::Matrix<NUMDOF_SOH8,1> elemat1(elemat1_epetra.A(),true);

      // calculate the THERMOmechanical solution
      soh8_linstiffmasstemp(la,mydisp,myres,mytempnp,&elemat1,&elemat2,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
    }

  }
  break;

  //==================================================================================
  // linear stiffness, internal force vector, and consistent mass matrix
  case calc_struct_nlnstiffmass:
  {
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
    soh8_linstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
                      INPAR::STR::stress_none,INPAR::STR::strain_none);

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
      const int iel = 8;

      if (la[1].Size() != iel*numdofpernode_)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);
      // build the current temperature vector
      LINALG::Matrix<iel*numdofpernode_,1> etemp(&(mytempnp[1]),true);  // view only!
      // the coupling stiffness matrix (3nx1)
      LINALG::Matrix<NUMDOF_SOH8,1> elemat1(elemat1_epetra.A(),true);
      // calculate the THERMOmechanical solution
      soh8_linstiffmasstemp(la,mydisp,myres,mytempnp,&elemat1,&elemat2,&elevec1,
        NULL,NULL,params,INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
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
    if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
    if (stressdata==Teuchos::null) dserror("Cannot get 'stress' data");
    if (straindata==Teuchos::null) dserror("Cannot get 'strain' data");
    vector<double> mydisp((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    vector<double> myres((la[0].lm_).size());
    DRT::UTILS::ExtractMyValues(*res,myres,lm);
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> stress;
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> strain;
    INPAR::STR::StressType iostress
      = params.get<INPAR::STR::StressType>("iostress", INPAR::STR::stress_none);
    INPAR::STR::StrainType iostrain
      = params.get<INPAR::STR::StrainType>("iostrain", INPAR::STR::strain_none);
    // call the well-known soh8_nlnstiffmass for the normal structure solution
    soh8_linstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,
      iostress,iostrain);

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
      const int iel = 8;

      if (la[1].Size() != iel*numdofpernode_)
        dserror("Location vector length for temperature does not match!");
      // extract the current temperatures
      DRT::UTILS::ExtractMyValues(*tempnp,mytempnp,la[1].lm_);
      // get the temperature dependent stress
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> stresstemp;
      // calculate the THERMOmechanical solution: temperature stresses
      soh8_linstiffmasstemp(la,mydisp,myres,mytempnp,NULL,NULL,NULL,
        &stresstemp,NULL,params,iostress,INPAR::STR::strain_none);

      // total stress
      // add stresstemp to the mechanical stress
      stress.Update(1.0,stresstemp,1.0);

    }

    AddtoPack(*stressdata, stress);
    AddtoPack(*straindata, strain);

  }
  break;

  //==================================================================================
  case calc_struct_update_istep:
  {
    // Update of history for visco material if they exist
  }
  break;

  //==================================================================================
  // postprocess stresses/strains at gauss points

  // note that in the following, quantities are always referred to as
  // "stresses" etc. although they might also apply to strains
  // (depending on what this routine is called for from the post filter)
  case postprocess_stress:
  {
    const Teuchos::RCP<map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap
      = params.get<Teuchos::RCP<map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
    if (gpstressmap==Teuchos::null)
      dserror("no gp stress/strain map available for postprocessing");
    string stresstype = params.get<string>("stresstype","ndxyz");
    int gid = Id();
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> gpstress(((*gpstressmap)[gid])->A(),true);

    if (stresstype=="ndxyz")
    {
      // extrapolate stresses/strains at Gauss points to nodes
      soh8_expol(gpstress, elevec1, elevec2);

    }
    else if (stresstype=="cxyz")
    {
      Teuchos::RCP<Epetra_MultiVector> elestress =
        params.get<Teuchos::RCP<Epetra_MultiVector> >("elestress",Teuchos::null);
      if (elestress==Teuchos::null)
        dserror("No element stress/strain vector available");
      const Epetra_BlockMap& elemap = elestress->Map();
      int lid = elemap.LID(Id());
      if (lid!=-1)
      {
        for (int i = 0; i < NUMSTR_SOH8; ++i)
        {
          double& s = (*((*elestress)(i)))[lid]; // resolve pointer for faster access
          s = 0.;
          for (int j = 0; j < NUMGPT_SOH8; ++j)
          {
            s += gpstress(j,i);
          }
          s *= 1.0/NUMGPT_SOH8;
        }
      }
    }
    else if (stresstype=="cxyz_ndxyz")
    {
      // extrapolate stresses/strains at Gauss points to nodes
      soh8_expol(gpstress, elevec1, elevec2);

      Teuchos::RCP<Epetra_MultiVector> elestress =
        params.get<Teuchos::RCP<Epetra_MultiVector> >("elestress",Teuchos::null);
      if (elestress==Teuchos::null)
        dserror("No element stress/strain vector available");
      const Epetra_BlockMap elemap = elestress->Map();
      int lid = elemap.LID(Id());
      if (lid!=-1) {
        for (int i = 0; i < NUMSTR_SOH8; ++i)
        {
          double& s = (*((*elestress)(i)))[lid];
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
  default:
    dserror("Unknown type of action for So_hex8");
  } // action

  return 0;

} // Evaluate


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                           dano 05/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_linstiffmass(
  vector<int>& lm,  // location matrix
  vector<double>& disp,  // current displacements
  vector<double>& residual,  // current residual displ
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,  // element stiffness matrix
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  // element mass matrix
  LINALG::Matrix<NUMDOF_SOH8,1>* force,  // element internal force vector
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,  // stresses at GP
  LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,  // strains at GP
  ParameterList& params,  // algorithmic parameters e.g. time
  const INPAR::STR::StressType iostress,  // stress output option
  const INPAR::STR::StrainType iostrain
  )  // strain output option
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
  for (int i=0; i<NUMDOF_SOH8; ++i)
  {
    nodaldisp(i,0) = disp[i];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
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

    // set to initial state as test to receive a linear solution
    for (int i=0; i<3; ++i) defgrd(i,i) = 1;

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

    // now build the linear strain
    LINALG::Matrix<NUMSTR_SOH8,1> strainlin(true);
    strainlin.Multiply(bop,nodaldisp);

    // and rename it as glstrain to use the common methods further on

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(NUMSTR_SOH8);
    LINALG::Matrix<NUMSTR_SOH8,1> glstrain(glstrain_epetra.A(),true);
    glstrain.Update(1.0,strainlin);

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

      // inverse of deformation gradient
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invdefgrd;
      invdefgrd.Invert(defgrd);

      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> euler_almansi;
      temp.Multiply(gl,invdefgrd);
      euler_almansi.MultiplyTN(invdefgrd,temp);

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

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(true);
    LINALG::Matrix<NUMSTR_SOH8,1> stress(true);
    soh8_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

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
      const double detF = defgrd.Determinant();

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
      temp.Multiply(1.0/detF,defgrd,pkstress,0.0);
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
    if (force != NULL && stiffmatrix != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,NUMDOF_SOH8> cb;
      cb.Multiply(cmat,bop);
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);
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
} // DRT::ELEMENTS::So_hex8::soh8_linstiffmass


/*----------------------------------------------------------------------*
 | evaluate only the temperature fraction for the element    dano 05/10 |
 | originally by maf 04/07  (private)                                   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_linstiffmasstemp(
  DRT::Element::LocationArray& la,  // location array
  vector<double>& disp,  // current displacements
  vector<double>& residual,  // current residual displ
  vector<double>& temp, // current temperature
  LINALG::Matrix<NUMDOF_SOH8,1>* tempstiffmatrix, // coupling stiffness matrix
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  // element mass matrix
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

  // vector of the current element temperatures
  LINALG::Matrix<NUMNOD_SOH8,1> etemp;
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
      defgrd(i,i) = 1;

    // Linear B-operator B = N_XYZ, e.g. with defgrd == I
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
    LINALG::Matrix<1,1> Ntemp;
    Ntemp.MultiplyTN(shapetemp,etemp);

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    // calculate the stress part dependent on the temperature in the material
    LINALG::Matrix<NUMSTR_SOH8,1> ctemp(true);
    LINALG::Matrix<NUMSTR_SOH8,1> stresstemp(true);
    soh8_mat_temp(&stresstemp,&ctemp,&density,&Ntemp,&defgrd,gp,params);

    // and now add the constant temperature fraction to stresstemp, too
    LINALG::Matrix<NUMSTR_SOH8,1> stempconst(true);
    Stempconst(&ctemp,&stempconst);
    // total temperature stress
    stresstemp.Update(1.0,stempconst,1.0);

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
      temp.Multiply(1.0/detF,defgrd,pkstresstemp,0.0);
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

    double detJ_w = detJ*gpweights[gp];
    if (force != NULL && tempstiffmatrix != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma_temp) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stresstemp, 1.0);

      //-----------------------------------------------------------------------
      // integrate `initial-displacement-temperature' stiffness matrix  02/10
      // update the stiffness part depending of the temperature
      // k_theta += 1.0 * k_theta + detJ * w(gp) * (B^T . C_Theta . N_Theta)
      //-----------------------------------------------------------------------
      LINALG::Matrix<NUMSTR_SOH8,1> cn; // 6x1 = 6x1n
      // extract the shapefunctions on the gp_i and multiply it with cn
      for (int stresscomp=0; stresscomp<NUMSTR_SOH8; ++stresscomp)
      {
        cn.Multiply(ctemp,Ntemp); // (6x1) = (6x1)(1x1)
      }
      tempstiffmatrix->MultiplyTN(detJ_w,bop,cn,1.0); // [(3nx6)(6x1)=(3nx1)] (24x6)(6x1)=(24x1)
    }

   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/
  return;
} // DRT::ELEMENTS::So_hex8::soh8_linstiffmasstemp


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

