/*!----------------------------------------------------------------------
\file so_disp_evaluate.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/


#include "so_disp.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_contact/contact_analytical.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/micromaterial.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoDisp::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // start with "none"
  DRT::ELEMENTS::SoDisp::ActionType act = SoDisp::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = SoDisp::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = SoDisp::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = SoDisp::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = SoDisp::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = SoDisp::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = SoDisp::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = SoDisp::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = SoDisp::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = SoDisp::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = SoDisp::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")   act = SoDisp::calc_struct_reset_istep;
  else if (action=="calc_struct_errornorms")    act = SoDisp::calc_struct_errornorms;
  else if (action=="calc_init_vol")             act = SoDisp::calc_init_vol;
  else dserror("Unknown type of action for SoDisp");

  // what should the element do
  switch(act) {
    // linear stiffness
    case calc_struct_linstiff: {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      sodisp_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1, params);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff: {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      sodisp_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,params);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce: {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      Epetra_SerialDenseMatrix myemat(lm.size(),lm.size());
      sodisp_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,params);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass: {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      sodisp_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1, params);
    }
    break;

    // evaluate stresses
    case calc_struct_stress: {
      dserror("Case calc_struct_stress not yet implemented");
    }
    break;

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    case calc_struct_update_istep: {
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
      {
        MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
        micro->Update();
      }
    }
    break;

    case calc_struct_update_imrlike: {
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
      {
        MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
        micro->Update();
      }
    }
    break;

    case calc_struct_reset_istep: {
      ;// there is nothing to do here at the moment
    }
    break;

    //==================================================================================
		case calc_struct_errornorms:
		{
			// IMPORTANT NOTES (popp 11/2010):
			// - error norms are based on a small deformation assumption (linear elasticity)
			// - extension to finite deformations would be possible without difficulties,
			//   however analytical solutions are extremely rare in the nonlinear realm
			// - only implemented for SVK material (relevant for energy norm only, L2 and
			//   H1 norms are of course valid for arbitrary materials)
			// - analytical solutions are currently stored in a repository in the CONTACT
			//   namespace, however they could (should?) be moved to a more general location

			// check length of elevec1
			if (elevec1.Length() < 3) dserror("The given result vector is too short.");

			// check material law
			RCP<MAT::Material> mat = Material();

			//******************************************************************
			// only for St.Venant Kirchhoff material
			//******************************************************************
			if (mat->MaterialType() == INPAR::MAT::m_stvenant)
			{
				// declaration of variables
				double l2norm = 0.0;
				double h1norm = 0.0;
				double energynorm = 0.0;

				// shape functions, derivatives and integration weights
				vector<Epetra_SerialDenseVector> shapefcts(numgpt_disp_);
				vector<Epetra_SerialDenseMatrix> derivs(numgpt_disp_);
				vector<double> weights(numgpt_disp_);
				sodisp_shapederiv(shapefcts,derivs,weights);

				// get displacements and extract values of this element
				RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
				if (disp==null) dserror("Cannot get state displacement vector");
				vector<double> mydisp(lm.size());
				DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

				// nodal displacement vector
				Epetra_SerialDenseVector nodaldisp(numdof_disp_);
				for (int i=0; i<numdof_disp_; ++i) nodaldisp[i] = mydisp[i];

				// reference geometry (nodal positions)
				Epetra_SerialDenseMatrix xrefe(numnod_disp_,NUMDIM_DISP);
				DRT::Node** nodes = Nodes();
				for (int i=0; i<numnod_disp_; ++i)
				{
					xrefe(i,0) = nodes[i]->X()[0];
					xrefe(i,1) = nodes[i]->X()[1];
					xrefe(i,2) = nodes[i]->X()[2];
				}

				//----------------------------------------------------------------
				// loop over all Gauss points
				//----------------------------------------------------------------
				for (int gp=0; gp<numgpt_disp_; gp++)
				{
					/* compute the Jacobian matrix which looks like:
					**         [ x_,r  y_,r  z_,r ]
					**     J = [ x_,s  y_,s  z_,s ]
					**         [ x_,t  y_,t  z_,t ]
					*/
					Epetra_SerialDenseMatrix jac(NUMDIM_DISP,NUMDIM_DISP);
					jac.Multiply('N','N',1.0,derivs[gp],xrefe,1.0);

					// compute determinant of Jacobian by Sarrus' rule
					double detJ= jac(0,0) * jac(1,1) * jac(2,2)
										 + jac(0,1) * jac(1,2) * jac(2,0)
										 + jac(0,2) * jac(1,0) * jac(2,1)
										 - jac(0,0) * jac(1,2) * jac(2,1)
										 - jac(0,1) * jac(1,0) * jac(2,2)
										 - jac(0,2) * jac(1,1) * jac(2,0);
					if (abs(detJ) < 1E-16) dserror("ZERO JACOBIAN DETERMINANT");
					else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

					// Gauss weights and Jacobian determinant
					double fac = detJ * weights[gp];

					// Gauss point in reference configuration
					LINALG::Matrix<NUMDIM_DISP,1> xgp(true);
					for (int k=0;k<NUMDIM_DISP;++k)
						for (int n=0;n<numnod_disp_;++n)
							xgp(k,0) += (shapefcts[gp])(n) * xrefe(n,k);

					//**************************************************************
					// get analytical solution
					LINALG::Matrix<NUMDIM_DISP,1> uanalyt(true);
					LINALG::Matrix<NUMSTR_DISP,1> strainanalyt(true);
					LINALG::Matrix<NUMDIM_DISP,NUMDIM_DISP> derivanalyt(true);

					CONTACT::AnalyticalSolutions3D(xgp,uanalyt,strainanalyt,derivanalyt);
					//**************************************************************

					//--------------------------------------------------------------
					// (1) L2 norm
					//--------------------------------------------------------------

					// compute displacements at GP
					LINALG::Matrix<NUMDIM_DISP,1> ugp(true);
					for (int k=0;k<NUMDIM_DISP;++k)
						for (int n=0;n<numnod_disp_;++n)
							ugp(k,0) += (shapefcts[gp])(n) * nodaldisp(NODDOF_DISP*n+k,0);

					// displacement error
					LINALG::Matrix<NUMDIM_DISP,1> uerror(true);
					for (int k=0;k<NUMDIM_DISP;++k)
						uerror(k,0) = uanalyt(k,0) - ugp(k,0);

					// compute GP contribution to L2 error norm
					l2norm += fac * uerror.Dot(uerror);

					//--------------------------------------------------------------
					// (2) H1 norm
					//--------------------------------------------------------------

					/* compute derivatives N_XYZ at gp w.r.t. material coordinates
					** by solving   Jac . N_XYZ = N_rst   for N_XYZ
					** Inverse of Jacobian is therefore not explicitly computed
					*/
					Epetra_SerialDenseMatrix N_XYZ(NUMDIM_DISP,numnod_disp_);
					Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
					solve_for_inverseJac.SetMatrix(jac);            // set A=jac
					solve_for_inverseJac.SetVectors(N_XYZ,derivs.at(gp));// set X=N_XYZ, B=deriv_gp
					solve_for_inverseJac.FactorWithEquilibration(true);
					int err2 = solve_for_inverseJac.Factor();
					int err = solve_for_inverseJac.Solve();         // N_XYZ = J^-1.N_rst
					if ((err != 0) && (err2!=0)) dserror("Inversion of Jacobian failed");

					// compute partial derivatives at GP
					LINALG::Matrix<NUMDIM_DISP,NUMDIM_DISP> derivgp(true);
					for (int l=0;l<NUMDIM_DISP;++l)
						for (int m=0;m<NUMDIM_DISP;++m)
							for (int k=0;k<numnod_disp_;++k)
								derivgp(l,m) += N_XYZ(m,k) * nodaldisp(NODDOF_DISP*k+l,0);

					// derivative error
					LINALG::Matrix<NUMDIM_DISP,NUMDIM_DISP> deriverror(true);
					for (int k=0;k<NUMDIM_DISP;++k)
						for (int m=0;m<NUMDIM_DISP;++m)
							deriverror(k,m) = derivanalyt(k,m) - derivgp(k,m);

					// compute GP contribution to H1 error norm
					h1norm += fac * deriverror.Dot(deriverror);
					h1norm += fac * uerror.Dot(uerror);

					//--------------------------------------------------------------
					// (3) Energy norm
					//--------------------------------------------------------------

					// compute linear B-operator
					Epetra_SerialDenseMatrix bop(NUMSTR_DISP,numdof_disp_);
					for (int i=0; i<numnod_disp_; ++i)
					{
						bop(0,NODDOF_DISP*i+0) = N_XYZ(0,i);
						bop(0,NODDOF_DISP*i+1) = 0.0;
						bop(0,NODDOF_DISP*i+2) = 0.0;
						bop(1,NODDOF_DISP*i+0) = 0.0;
						bop(1,NODDOF_DISP*i+1) = N_XYZ(1,i);
						bop(1,NODDOF_DISP*i+2) = 0.0;
						bop(2,NODDOF_DISP*i+0) = 0.0;
						bop(2,NODDOF_DISP*i+1) = 0.0;
						bop(2,NODDOF_DISP*i+2) = N_XYZ(2,i);

						bop(3,NODDOF_DISP*i+0) = N_XYZ(1,i);
						bop(3,NODDOF_DISP*i+1) = N_XYZ(0,i);
						bop(3,NODDOF_DISP*i+2) = 0.0;
						bop(4,NODDOF_DISP*i+0) = 0.0;
						bop(4,NODDOF_DISP*i+1) = N_XYZ(2,i);
						bop(4,NODDOF_DISP*i+2) = N_XYZ(1,i);
						bop(5,NODDOF_DISP*i+0) = N_XYZ(2,i);
						bop(5,NODDOF_DISP*i+1) = 0.0;
						bop(5,NODDOF_DISP*i+2) = N_XYZ(0,i);
					}

					// compute linear strain at GP
					Epetra_SerialDenseVector straingptmp(NUMSTR_DISP);
					bop.Multiply(false,nodaldisp,straingptmp);
		      LINALG::Matrix<NUMSTR_DISP,1> straingp(true);
					for (int k=0;k<NUMSTR_DISP;++k) straingp(k,0) = straingptmp[k];

					// strain error
					LINALG::Matrix<NUMSTR_DISP,1> strainerror(true);
					for (int k=0;k<NUMSTR_DISP;++k)
						strainerror(k,0) = strainanalyt(k,0) - straingp(k,0);

					// compute stress vector and constitutive matrix
					double density = 0.0;
					LINALG::Matrix<NUMSTR_DISP,NUMSTR_DISP> cmat(true);
					LINALG::Matrix<NUMSTR_DISP,1> stress(true);
					sodisp_mat_sel(&stress,&cmat,&density,&strainerror,params);

					// compute GP contribution to energy error norm
					energynorm += fac * stress.Dot(strainerror);

					//cout << "UAnalytical:      " << ugp << endl;
					//cout << "UDiscrete:        " << uanalyt << endl;
					//cout << "StrainAnalytical: " << strainanalyt << endl;
					//cout << "StrainDiscrete:   " << straingp << endl;
					//cout << "DerivAnalytical:  " << derivanalyt << endl;
					//cout << "DerivDiscrete:    " << derivgp << endl;
				}
				//----------------------------------------------------------------

				// return results
				elevec1(0) = l2norm;
				elevec1(1) = h1norm;
				elevec1(2) = energynorm;
			}
			else
				dserror("ERROR: Error norms only implemented for SVK material");
		}
		break;

    default:
      dserror("Unknown type of action for Solid3");
  }
  return 0;
}



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoDisp::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  dserror("This element does not do body force or similar");
  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::sodisp_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residual displ
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      ParameterList&            params)         // algorithmic parameters e.g. time
{

    vector<Epetra_SerialDenseVector> shapefcts(numgpt_disp_);
    vector<Epetra_SerialDenseMatrix> derivs(numgpt_disp_);
    vector<double> weights(numgpt_disp_);
    sodisp_shapederiv(shapefcts,derivs,weights);   // call to evaluate
  /* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(numnod_disp_,NUMDIM_DISP);  // material coord. of element
  Epetra_SerialDenseMatrix xcurr(numnod_disp_,NUMDIM_DISP);  // current  coord. of element
  for (int i=0; i<numnod_disp_; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_DISP+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_DISP+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_DISP+2];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<numgpt_disp_; ++gp) {

    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    Epetra_SerialDenseMatrix jac(NUMDIM_DISP,NUMDIM_DISP);
    jac.Multiply('N','N',1.0,derivs.at(gp),xrefe,1.0);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ= jac(0,0) * jac(1,1) * jac(2,2)
               + jac(0,1) * jac(1,2) * jac(2,0)
               + jac(0,2) * jac(1,0) * jac(2,1)
               - jac(0,0) * jac(1,2) * jac(2,1)
               - jac(0,1) * jac(1,0) * jac(2,2)
               - jac(0,2) * jac(1,1) * jac(2,0);
    if (abs(detJ) < 1E-16) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    /* compute derivatives N_XYZ at gp w.r.t. material coordinates
    ** by solving   Jac . N_XYZ = N_rst   for N_XYZ
    ** Inverse of Jacobian is therefore not explicitly computed
    */
    Epetra_SerialDenseMatrix N_XYZ(NUMDIM_DISP,numnod_disp_);
    Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    solve_for_inverseJac.SetVectors(N_XYZ,derivs.at(gp));// set X=N_XYZ, B=deriv_gp
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();         // N_XYZ = J^-1.N_rst
    if ((err != 0) && (err2!=0)) dserror("Inversion of Jacobian failed");

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    Epetra_SerialDenseMatrix defgrd(NUMDIM_DISP,NUMDIM_DISP);
    defgrd.Multiply('T','T',1.0,xcurr,N_XYZ,1.0);

    // Right Cauchy-Green tensor = F^T * F
    Epetra_SerialDenseMatrix cauchygreen(NUMDIM_DISP,NUMDIM_DISP);
    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,1.0);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain(NUMSTR_DISP);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);


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
    Epetra_SerialDenseMatrix bop(NUMSTR_DISP,numdof_disp_);
    for (int i=0; i<numnod_disp_; ++i) {
      bop(0,NODDOF_DISP*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NODDOF_DISP*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NODDOF_DISP*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NODDOF_DISP*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NODDOF_DISP*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NODDOF_DISP*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NODDOF_DISP*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NODDOF_DISP*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NODDOF_DISP*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NODDOF_DISP*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NODDOF_DISP*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NODDOF_DISP*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NODDOF_DISP*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NODDOF_DISP*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NODDOF_DISP*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NODDOF_DISP*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NODDOF_DISP*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NODDOF_DISP*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    LINALG::Matrix<NUMSTR_DISP,NUMSTR_DISP> cmat_f(true);
    LINALG::Matrix<NUMSTR_DISP,1> stress_f(true);
    LINALG::Matrix<NUMSTR_DISP,1> glstrain_f(glstrain.A());
    double density = 0.0;
    sodisp_mat_sel(&stress_f,&cmat_f,&density,&glstrain_f, params);
    Epetra_SerialDenseMatrix cmat(View,cmat_f.A(),cmat_f.Rows(),cmat_f.Rows(),cmat_f.Columns());
    Epetra_SerialDenseVector stress(View,stress_f.A(),stress_f.Rows());
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
    (*force).Multiply('T','N',detJ * weights.at(gp),bop,stress,1.0);

    // integrate `elastic' and `initial-displacement' stiffness matrix
    // keu = keu + (B^T . C . B) * detJ * w(gp)
    Epetra_SerialDenseMatrix cb(NUMSTR_DISP,numdof_disp_);
    cb.Multiply('N','N',1.0,cmat,bop,1.0);          // temporary C . B
    (*stiffmatrix).Multiply('T','N',detJ * weights.at(gp),bop,cb,1.0);

    // integrate `geometric' stiffness matrix and add to keu *****************
    Epetra_SerialDenseVector sfac(stress); // auxiliary integrated stress
    sfac.Scale(detJ * weights.at(gp));     // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
    vector<double> SmB_L(NUMDIM_DISP);     // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod=0; inod<numnod_disp_; ++inod){
      SmB_L[0] = sfac(0) * N_XYZ(0,inod) + sfac(3) * N_XYZ(1,inod) + sfac(5) * N_XYZ(2,inod);
      SmB_L[1] = sfac(3) * N_XYZ(0,inod) + sfac(1) * N_XYZ(1,inod) + sfac(4) * N_XYZ(2,inod);
      SmB_L[2] = sfac(5) * N_XYZ(0,inod) + sfac(4) * N_XYZ(1,inod) + sfac(2) * N_XYZ(2,inod);
      for (int jnod=0; jnod<numnod_disp_; ++jnod){
        double bopstrbop = 0.0;            // intermediate value
        for (int idim=0; idim<NUMDIM_DISP; ++idim) bopstrbop += N_XYZ(idim,jnod) * SmB_L[idim];
        (*stiffmatrix)(NUMDIM_DISP*inod+0,NUMDIM_DISP*jnod+0) += bopstrbop;
        (*stiffmatrix)(NUMDIM_DISP*inod+1,NUMDIM_DISP*jnod+1) += bopstrbop;
        (*stiffmatrix)(NUMDIM_DISP*inod+2,NUMDIM_DISP*jnod+2) += bopstrbop;
      }
    } // end of integrate `geometric' stiffness ******************************


    if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
      // integrate consistent mass matrix
      for (int inod=0; inod<numnod_disp_; ++inod) {
        for (int jnod=0; jnod<numnod_disp_; ++jnod) {
          double massfactor = (shapefcts.at(gp))(inod) * density * (shapefcts.at(gp))(jnod)
                            * detJ * weights.at(gp);     // intermediate factor
          (*massmatrix)(NUMDIM_DISP*inod+0,NUMDIM_DISP*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_DISP*inod+1,NUMDIM_DISP*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_DISP*inod+2,NUMDIM_DISP*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
}

/*----------------------------------------------------------------------*
 |  shape functions and derivatives for SoDisp                 maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::sodisp_shapederiv(
      vector<Epetra_SerialDenseVector>& shapefcts,  // pointer to pointer of shapefct
      vector<Epetra_SerialDenseMatrix>& derivs,     // pointer to pointer of derivs
      vector<double>& weights)   // pointer to pointer of weights
{

  const DRT::Element::DiscretizationType distype = Shape();

  // (r,s,t) gp-locations of fully integrated linear 6-node Wedge
  // fill up nodal f at each gp
  // fill up df w.r.t. rst directions (NUMDIM) at each gp
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule_);
  for (int igp = 0; igp < intpoints.nquad; ++igp) {
    const double r = intpoints.qxg[igp][0];
    const double s = intpoints.qxg[igp][1];
    const double t = intpoints.qxg[igp][2];

    shapefcts[igp].Size(numnod_disp_);
    derivs[igp].Shape(NUMDIM_DISP, numnod_disp_);
    DRT::UTILS::shape_function_3D(shapefcts[igp], r, s, t, distype);
    DRT::UTILS::shape_function_3D_deriv1(derivs[igp], r, s, t, distype);
    weights[igp] = intpoints.qwgt[igp];
  }

  return;
}



