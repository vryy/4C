/*!----------------------------------------------------------------------
\file so_hex8fbar_evaluate.cpp
\brief

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_hex8fbar.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbar::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH8,1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH8,1> elevec2(elevec2_epetra.A(),true);
  // elevec3 is not used anyway

  // start with "none"
  DRT::ELEMENTS::So_hex8fbar::ActionType act = So_hex8fbar::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")                        act = So_hex8fbar::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")                        act = So_hex8fbar::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")                   act = So_hex8fbar::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")                    act = So_hex8fbar::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")                    act = So_hex8fbar::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")                   act = So_hex8fbar::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")                          act = So_hex8fbar::calc_struct_stress;
  else if (action=="calc_struct_eleload")                         act = So_hex8fbar::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")                         act = So_hex8fbar::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")                    act = So_hex8fbar::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike")                  act = So_hex8fbar::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")                     act = So_hex8fbar::calc_struct_reset_istep;
  else if (action=="calc_homog_dens")                             act = So_hex8fbar::calc_homog_dens;
  else if (action=="postprocess_stress")                          act = So_hex8fbar::postprocess_stress;
  else if (action=="multi_readrestart")                           act = So_hex8fbar::multi_readrestart;
  else dserror("Unknown type of action for So_hex8fbar");
  // what should the element do
  switch(act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (unsigned i=0; i<mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (unsigned i=0; i<myres.size(); ++i) myres[i] = 0.0;
      soh8fbar_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      soh8fbar_nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> myemat(true);
      soh8fbar_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      soh8fbar_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      if (act==calc_struct_nlnstifflmass) soh8_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID()==Owner())
      {
        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
        RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
        RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
        if (disp==null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata==null) dserror("Cannot get 'stress' data");
        if (straindata==null) dserror("Cannot get 'strain' data");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res,myres,lm);
        LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> stress;
        LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> strain;
        INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
        soh8fbar_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,iostress,iostrain);
        AddtoPack(*stressdata, stress);
        AddtoPack(*straindata, strain);
      }
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID()==Owner())
      {
        const RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
          params.get<RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
        if (gpstressmap==null)
          dserror("no gp stress/strain map available for postprocessing");
        string stresstype = params.get<string>("stresstype","ndxyz");
        int gid = Id();
        LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8> gpstress(((*gpstressmap)[gid])->A(),true);

        Teuchos::RCP<Epetra_MultiVector> poststress=params.get<Teuchos::RCP<Epetra_MultiVector> >("poststress",null);
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

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    case calc_struct_update_istep:
    {
      // Update of history for plastic material
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_plneohooke)
      {
        MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(mat.get());
        plastic->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_growth)
      {
        MAT::Growth* grow = static_cast <MAT::Growth*>(mat.get());
        grow->Update();
      }
    }
    break;

    case calc_struct_update_imrlike:
    {
      // Update of history for plastic material
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_plneohooke)
      {
        MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(mat.get());
        plastic->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_growth)
      {
        MAT::Growth* grow = static_cast <MAT::Growth*>(mat.get());
        grow->Update();
      }
    }
    break;

    case calc_struct_reset_istep:
    {
    	// Update of history for plastic material
			RCP<MAT::Material> mat = Material();
			if (mat->MaterialType() == INPAR::MAT::m_plneohooke)
			{
				MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(mat.get());
				plastic->Update();
			}
    }
    break;

    case calc_homog_dens:
    {
      soh8_homog(params);
    }
    break;


    // read restart of microscale
    case multi_readrestart:
    {
      RefCountPtr<MAT::Material> mat = Material();

      if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
        soh8_read_restart_multi(params);
    }
    break;

    default:
      dserror("Unknown type of action for So_hex8fbar");
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)               |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbar::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
  // **

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_20 with 20 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  const static vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry
   LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) {

    // compute the Jacobian matrix
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> jac;
    jac.Multiply(derivs[gp],xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    double fac = gpweights[gp] * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
      for(int dim=0; dim<NUMDIM_SOH8; dim++) {
      double dim_fac = (*onoff)[dim] * (*val)[dim] * fac;
      for (int nodid=0; nodid<NUMNOD_SOH8; ++nodid) {
        elevec1[nodid*NUMDIM_SOH8+dim] += shapefcts[gp](nodid) * dim_fac;
      }
    }

  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::So_hex8fbar::EvaluateNeumann

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8fbar::soh8fbar_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residual displ
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix, // element stiffness matrix
      LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  // element mass matrix
      LINALG::Matrix<NUMDOF_SOH8,1>* force,                 // element internal force vector
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,   // stresses at GP
      LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,   // strains at GP
      ParameterList&            params,         // algorithmic parameters e.g. time
      const INPAR::STR::StressType   iostress,  // stress output option
      const INPAR::STR::StrainType   iostrain)  // strain output option
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

  //****************************************************************************
	// deformation gradient at centroid of element
  //****************************************************************************
	//element coordinate derivatives at centroid
	LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_rst_0;
	DRT::UTILS::shape_function_3D_deriv1(N_rst_0, 0, 0, 0, hex8);
	//inverse jacobian matrix at centroid
	LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invJ_0;
	invJ_0.Multiply(N_rst_0,xrefe);
	invJ_0.Invert();
	//material derivatives at centroid
	LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ_0;
	N_XYZ_0.Multiply(invJ_0,N_rst_0);
	//deformation gradient and its determinant at centroid
	LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd_0(false);
	defgrd_0.MultiplyTT(xcurr,N_XYZ_0);
	LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invdefgrd_0;
	invdefgrd_0.Invert(defgrd_0);
	double detF_0=defgrd_0.Determinant();

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(false);
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

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,N_XYZ);
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invdefgrd;
    invdefgrd.Invert(defgrd);
    double detF=defgrd.Determinant();

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // F_bar deformation gradient =(detF_0/detF)^1/3*F
		LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd_bar(defgrd);
		double f_bar_factor=pow(detF_0/detF,1.0/3.0);
		defgrd_bar.Scale(f_bar_factor);

		// Right Cauchy-Green tensor(Fbar) = F_bar^T * F_bar
		LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> cauchygreen_bar;
		cauchygreen_bar.MultiplyTN(defgrd_bar,defgrd_bar);

    // Green-Lagrange strains(F_bar) matrix E = 0.5 * (Cauchygreen(F_bar) - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
		Epetra_SerialDenseVector glstrain_bar_epetra(NUMSTR_SOH8);
		LINALG::Matrix<NUMSTR_SOH8,1> glstrain_bar(glstrain_bar_epetra.A(),true);
		glstrain_bar(0) = 0.5 * (cauchygreen_bar(0,0) - 1.0);
		glstrain_bar(1) = 0.5 * (cauchygreen_bar(1,1) - 1.0);
		glstrain_bar(2) = 0.5 * (cauchygreen_bar(2,2) - 1.0);
		glstrain_bar(3) = cauchygreen_bar(0,1);
		glstrain_bar(4) = cauchygreen_bar(1,2);
		glstrain_bar(5) = cauchygreen_bar(2,0);

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
    case INPAR::STR::strain_gl:
    {
      if (elestrain == NULL) dserror("strain data not available");
      for (int i = 0; i < 3; ++i)
        (*elestrain)(gp,i) = glstrain_bar(i);
      for (int i = 3; i < 6; ++i)
        (*elestrain)(gp,i) = 0.5 * glstrain_bar(i);
    }
    break;
    case INPAR::STR::strain_ea:
    {
      if (elestrain == NULL) dserror("strain data not available");
      // rewriting Green-Lagrange strains in matrix format
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> gl_bar;
      gl_bar(0,0) = glstrain_bar(0);
			gl_bar(0,1) = 0.5*glstrain_bar(3);
			gl_bar(0,2) = 0.5*glstrain_bar(5);
			gl_bar(1,0) = gl_bar(0,1);
			gl_bar(1,1) = glstrain_bar(1);
			gl_bar(1,2) = 0.5*glstrain_bar(4);
			gl_bar(2,0) = gl_bar(0,2);
			gl_bar(2,1) = gl_bar(1,2);
			gl_bar(2,2) = glstrain_bar(2);

	    // inverse of fbar deformation gradient
			LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invdefgrd_bar;
			invdefgrd_bar.Invert(defgrd_bar);

			LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
			LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> euler_almansi_bar;
			temp.Multiply(gl_bar,invdefgrd_bar);
			euler_almansi_bar.MultiplyTN(invdefgrd_bar,temp);

			(*elestrain)(gp,0) = euler_almansi_bar(0,0);
			(*elestrain)(gp,1) = euler_almansi_bar(1,1);
			(*elestrain)(gp,2) = euler_almansi_bar(2,2);
			(*elestrain)(gp,3) = euler_almansi_bar(0,1);
			(*elestrain)(gp,4) = euler_almansi_bar(1,2);
			(*elestrain)(gp,5) = euler_almansi_bar(0,2);
    }
    break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror("requested strain type not available");
    }

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

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    LINALG::Matrix<NUMSTR_SOH8,NUMSTR_SOH8> cmat(true);
    LINALG::Matrix<NUMSTR_SOH8,1> stress_bar(true);
    soh8_mat_sel(&stress_bar,&cmat,&density,&glstrain_bar,&defgrd_bar,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < NUMSTR_SOH8; ++i)
        (*elestress)(gp,i) = stress_bar(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      const double detF_bar = defgrd_bar.Determinant();

      LINALG::Matrix<3,3> pkstress_bar;
			pkstress_bar(0,0) = stress_bar(0);
			pkstress_bar(0,1) = stress_bar(3);
			pkstress_bar(0,2) = stress_bar(5);
			pkstress_bar(1,0) = pkstress_bar(0,1);
			pkstress_bar(1,1) = stress_bar(1);
			pkstress_bar(1,2) = stress_bar(4);
			pkstress_bar(2,0) = pkstress_bar(0,2);
			pkstress_bar(2,1) = pkstress_bar(1,2);
			pkstress_bar(2,2) = stress_bar(2);

			LINALG::Matrix<3,3> temp;
			LINALG::Matrix<3,3> cauchystress_bar;
			temp.Multiply(1.0/detF_bar,defgrd_bar,pkstress_bar,0.0);
			cauchystress_bar.MultiplyNT(temp,defgrd_bar);

			(*elestress)(gp,0) = cauchystress_bar(0,0);
			(*elestress)(gp,1) = cauchystress_bar(1,1);
			(*elestress)(gp,2) = cauchystress_bar(2,2);
			(*elestress)(gp,3) = cauchystress_bar(0,1);
			(*elestress)(gp,4) = cauchystress_bar(1,2);
			(*elestress)(gp,5) = cauchystress_bar(0,2);
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
      force->MultiplyTN(detJ_w/f_bar_factor, bop, stress_bar, 1.0);

      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,NUMDOF_SOH8> cb;
      cb.Multiply(cmat,bop);
      stiffmatrix->MultiplyTN(detJ_w*f_bar_factor,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6,1> sfac(stress_bar); // auxiliary integrated stress
      sfac.Scale(detJ_w/f_bar_factor); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      vector<double> SmB_L(3); // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOH8; ++inod) {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
            + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
            + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
            + sfac(2) * N_XYZ(2, inod);
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod) {
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<NUMDIM_SOH8; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3*inod+0,3*jnod+0) += bopstrbop;
          (*stiffmatrix)(3*inod+1,3*jnod+1) += bopstrbop;
          (*stiffmatrix)(3*inod+2,3*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************

      // integrate additional fbar matrix
			LINALG::Matrix<NUMSTR_SOH8,1> cauchygreenvector;
			cauchygreenvector(0) = cauchygreen(0,0);
			cauchygreenvector(1) = cauchygreen(1,1);
			cauchygreenvector(2) = cauchygreen(2,2);
			cauchygreenvector(3) = 2*cauchygreen(0,1);
			cauchygreenvector(4) = 2*cauchygreen(1,2);
			cauchygreenvector(5) = 2*cauchygreen(2,0);

			LINALG::Matrix<NUMSTR_SOH8,1> ccg;
			ccg.Multiply(cmat,cauchygreenvector);

			LINALG::Matrix<NUMDOF_SOH8,1> bopccg(false); // auxiliary integrated stress
			bopccg.MultiplyTN(detJ_w*f_bar_factor/3.0,bop,ccg);

			double htensor[NUMDOF_SOH8];
			for(int n=0;n<NUMDOF_SOH8;n++)
			{
				htensor[n]=0;
				for(int i=0;i<NUMDIM_SOH8;i++)
					{
						htensor[n] += invdefgrd_0(i,n%3)*N_XYZ_0(i,n/3)-invdefgrd(i,n%3)*N_XYZ(i,n/3);
					}
			}

			LINALG::Matrix<NUMDOF_SOH8,1> bops(false); // auxiliary integrated stress
			bops.MultiplyTN(-detJ_w/f_bar_factor/3.0,bop,stress_bar);
			for(int i=0;i<NUMDOF_SOH8;i++)
			{
				for (int j=0;j<NUMDOF_SOH8;j++)
					{
						(*stiffmatrix)(i,j) += htensor[j]*(bops(i,0)+bopccg(i,0));
					}
			} // end of integrate additional `fbar' stiffness**********************
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

  }/* ==================================================== end of Loop over GP */

  return;
} // DRT::ELEMENTS::So_hex8fbar::soh8fbar_nlnstiffmass

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8fbarType::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_hex8fbar* actele = dynamic_cast<DRT::ELEMENTS::So_hex8fbar*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8fbar* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
