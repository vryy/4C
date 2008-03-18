/*!----------------------------------------------------------------------*
\file so_tet4_evaluate.cpp
\brief quadratic nonlinear tetrahedron 

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by : Alexander Volf
			alexander.volf@mytum.de     
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_tet4.H"
#include "so_hex8.H"
#include "so_integrator.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"

extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../drt_lib/dstrc.H"

//#define VERBOSE_OUTPUT
using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |                                                         maf 04/07    |
 | vector of material laws                                              |
 | defined in global_control.c											|
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              vlf 06/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // start with "none"
  DRT::ELEMENTS::So_tet4::ActionType act = So_tet4::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = So_tet4::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = So_tet4::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = So_tet4::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = So_tet4::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = So_tet4::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = So_tet4::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = So_tet4::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = So_tet4::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = So_tet4::calc_struct_update_istep;
  else if (action=="calc_struct_update_genalpha_imrlike") act = So_tet4::calc_struct_update_genalpha_imrlike;
  else dserror("Unknown type of action for So_tet4");

  // get the material law
  MATERIAL* actmat = &(mat[material_-1]);

  // what should the element do
  switch(act) {
    // linear stiffness
    case calc_struct_linstiff: {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      so_tet4_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,actmat);

    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff: {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      so_tet4_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,actmat);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
      dserror("Case 'calc_struct_internalforce' not yet implemented");
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass: {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      so_tet4_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
    }
    break;

    // evaluate stresses
    case calc_struct_stress: {
      dserror("Case calc_struct_stress not yet implemented");
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      Epetra_SerialDenseMatrix stresses(NUMGPT_SOTET4,NUMSTR_SOTET4);//*
      so_tet4_stress(actmat,mydisp,&stresses); //*
    }
    break;

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    case calc_struct_update_istep: {
      ;// there is nothing to do here at the moment
    }
    break;

    case calc_struct_update_genalpha_imrlike: {
      ;// there is nothing to do here at the moment
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
int DRT::ELEMENTS::So_tet4::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  dserror("DRT::ELEMENTS::So_tet4::EvaluateNeumann not implemented");
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
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
  // **
  
/* =============================================================================*
 * CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_4 with 1 GAUSS POINTS*
 * =============================================================================*/
  const static DRT::ELEMENTS::Integrator_tet4_1point tet4_int;
/* ============================================================================*/

/* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOTET4; gp++) {
    // get submatrix of deriv at actual gp
   
    /* get the matrix of the coordinates of edges needed to compute the volume,
    ** which is used here as detJ in the quadrature rule.  
    ** ("Jacobian matrix") for the quadrarture rule:
    **             [  1    1    1  	 1  ]
    ** jac_coord = [ x_1  x_2  x_3  x_4 ]
    **             [ y_1  y_2  y_3  y_4 ]
    **		       [ z_1  z_2  z_3  z_4 ]
    */
    Epetra_SerialDenseMatrix jac_coord(NUMCOORD_SOTET4,NUMCOORD_SOTET4);
    for (int i=0; i<4; i++) jac_coord(0,i)=1;
    for (int row=0;row<3;row++)
    {
    	for (int col=0;col<4;col++){
    	jac_coord(row+1,col)= Nodes()[col]->X()[row];}   	
    }
    
    // compute determinant of Jacobian with own algorithm
    // !!warning detJ is not the actual determinant of the jacobian (here needed for the quadrature rule)
    // but rather the volume of the tetrahedara
    double detJ=det_volf(jac_coord);
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    double fac = tet4_int.weights[gp] * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
    for (int nodid=0; nodid<NUMNOD_SOTET4; ++nodid) {
      for(int dim=0; dim<NUMDIM_SOTET4; dim++) {
        elevec1[nodid*NUMDIM_SOTET4+dim] += tet4_int.shapefct_gp[gp](nodid) * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }
  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::So_tet4::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            vlf 08/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      struct _MATERIAL*         material)       // element material data
{
/* =============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_4  with 1 GAUSS POINTS*
** =============================================================================*/
   const static DRT::ELEMENTS::Integrator_tet4_1point tet4_dis;
/* ============================================================================*/
  double density;
  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOTET4,NUMDIM_SOTET4);  // material coord. of element
  /* structure of xrefe:
    **             [  X_1   Y_1   Z_1  ]
    **     xrefe = [  X_2   Y_2   Z_2  ]
    **             [   |     |     |   ]
    **             [  X_4   Y_4   Z_4  ]
    */
  Epetra_SerialDenseMatrix xcurr(NUMNOD_SOTET4,NUMDIM_SOTET4);  // current  coord. of element
  /* structure of xcurr:
    **             [  x_1   y_1   z_1  ]
    **     xcurr = [  x_2   y_2   z_2  ]
    **             [   |     |     |   ]
    **             [  x_4   y_4   z_4  ]
    */   
  Epetra_SerialDenseMatrix xdisp(NUMNOD_SOTET4,NUMDIM_SOTET4); // current  displacements of element
  for (int i=0; i<NUMNOD_SOTET4; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOTET4+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOTET4+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOTET4+2];
    
    xdisp(i,0) = disp[i*NODDOF_SOTET4+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET4+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET4+2];
  }
 
 
  /* get the matrix of the coordinates of edges needed to compute the volume,
  ** which is used here as detJ in the quadrature rule.  
  ** ("Jacobian matrix") for the quadrarture rule:
  **             [  1    1    1    1  ]
  ** jac_coord = [ X_1  X_2  X_3  X_4 ]
  **             [ Y_1  Y_2  Y_3  Y_4 ]
  **		     [ Z_1  Z_2  Z_3  Z_4 ]
  */
  // compute determinant of Jacobian with own algorithm
  // !!warning detJ is not the actual determinant of the jacobian (here needed for the quadrature rule)
  // but rather the volume of 
  
  Epetra_SerialDenseMatrix jac_coord(NUMCOORD_SOTET4,NUMCOORD_SOTET4);
  for (int i=0; i<4; i++)  jac_coord(0,i)=1;
  for (int row=0;row<3;row++)
  {
   	for (int col=0;col<4;col++)
       jac_coord(row+1,col)= xrefe(col,row);	
  }
  
  double detJ=det_volf(jac_coord)/(double) 6;    //volume of a tetrahedra
  if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
  else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");
  
  
  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOTET4; gp++) {
    
    /* compute the Jacobian matrix which looks like this:
    **         [  1        1        1  	     1      ]
    **   jac = [ X_,xsi1  X_,xsi2  X_,xsi3  X_,xsi4 ]
    **         [ Y_,xsi1  Y_,xsi2  Y_,xsi3  Y_,xsi4 ]
    **		   [ Z_,xsi1  Z_,xsi2  Z_,xsi3  Z_,xsi4 ]
    */
    
    Epetra_SerialDenseMatrix jac_temp(NUMCOORD_SOTET4-1,NUMCOORD_SOTET4);
    Epetra_SerialDenseMatrix jac(NUMCOORD_SOTET4,NUMCOORD_SOTET4);
    jac_temp.Multiply('T','N',1.0,xrefe,tet4_dis.deriv_gp[gp],0.0);
   
    for (int i=0; i<4; i++) jac(0,i)=1;
    for (int row=0;row<3;row++)
    {
    	for (int col=0;col<4;col++)
    	jac(row+1,col)=jac_temp(row,col);	
    }

    #ifdef VERBOSE_OUTPUT
    cout << "jac\n" << jac;
    cout << "xrefe\n" << xrefe;
	#endif //VERBOSE_OUTPUT

    /* compute partial derivatives at gp xsi_1, xsi_2, xsi_3, xsi_4 material coordinates
    ** by solving   Jac . partials = I_aug   for partials
    ** Inverse of Jacobian is therefore not explicitly computed
    */ 
    /* structure of partials:
    **             [  dxsi_1   dxsi_1   dxsi_1  ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    **             [     |        |        |    ]
    ** partials =  [                            ]
    **             [  dxsi_4   dxsi_4   dxsi_4  ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    */   
    
    Epetra_SerialDenseMatrix I_aug(NUMCOORD_SOTET4,NUMDIM_SOTET4);
    Epetra_SerialDenseMatrix partials(NUMCOORD_SOTET4,NUMDIM_SOTET4);
    Epetra_SerialDenseMatrix N_XYZ(NUMNOD_SOTET4,NUMDIM_SOTET4);
    I_aug(1,0)=1;
	I_aug(2,1)=1;
	I_aug(3,2)=1;
	
	#ifdef VERBOSE_OUTPUT
	cout << "jac\n" << jac;
	#endif //VERBOSE_OUTPUT
	
    Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    solve_for_inverseJac.SetVectors(partials,I_aug);// set X=partials, B=I_aug
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();        
    int err = solve_for_inverseJac.Solve();         // partials = jac^-1.I_aug
    if ((err != 0) && (err2!=0)){
    	dserror("Inversion of Jacobian failed");
    }

    #ifdef VERBOSE_OUTPUT
    cout << "I_aug\n" << I_aug;
    cout << "partials\n" << partials;
    cout << "deriv_gp\n" << deriv_gp;
    #endif //VERBOSE_OUTPUT

    N_XYZ.Multiply('N','N',1.0,tet4_dis.deriv_gp[gp],partials,0.0); //N_XYZ = N_xsi_k*partials
      
    /* structure of N_XYZ:
    **             [   dN_1     dN_1     dN_1   ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    **    N_XYZ =  [     |        |        |    ]
    **             [                            ]
    **             [   dN_4     dN_4     dN_4   ]
    **             [  -------  -------  ------- ]
    **             [    dX       dY       dZ    ]
    */
    
    #ifdef VERBOSE_OUTPUT
    cout << "N_XYZ\n" << N_XYZ;
    #endif //VERBOSE_OUTPUT
    //                                      d xcurr 
    // (material) deformation gradient F = --------- = xcurr^T * N_XYZ^T
    //                                      d xrefe  
    
    /*structure of F
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dX       dX       dX    ]
    **             [                            ]
    **      F   =  [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dY       dY       dY    ]
    **             [                            ]
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dZ       dZ       dZ    ]
    */

    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOTET4,NUMDIM_SOTET4);
    defgrd.Multiply('T','N',1.0,xdisp,N_XYZ,0.0);
    defgrd(0,0)+=1;
    defgrd(1,1)+=1;
    defgrd(2,2)+=1;
    
    //defgrd.Multiply('T','N',1.0,xcurr,N_XYZ,0.0);
    
    #ifdef VERBOSE_OUTPUT
	cout << "defgr\n " << defgrd;
	#endif //VERBOSE_OUTPUT
	
    // Right Cauchy-Green tensor = F^T * F 
    Epetra_SerialDenseMatrix cauchygreen(NUMDIM_SOTET4,NUMDIM_SOTET4);
   
    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,0.0);
    
	#ifdef VERBOSE_OUTPUT
	cout << "cauchygreen\n" << cauchygreen;
	getchar();
	#endif //VERBOSE_OUTPUT
	
    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain(NUMSTR_SOTET4);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);
    
    #ifdef VERBOSE_OUTPUT
	cout << "glstrain\n" << glstrain;
	#endif //VERBOSE_OUTPUT
	
	/*----------------------------------------------------------------------*
      the B-operator used is equivalent to the one used in hex8, this needs
      to be checked if it is ok, but from the mathematics point of view, the only
      thing that needed to be changed is tho NUMDOF       vlf 07/07
      ----------------------------------------------------------------------*/
    /* non-linear B-operator (may be so called, meaning 
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
    
    Epetra_SerialDenseMatrix bop(NUMSTR_SOTET4,NUMDOF_SOTET4);
    #ifdef VERBOSE_OUTPUT
    cout << bop;
    cout << defgrd;
    cout << N_XYZ;
    #endif //VERBOSE_OUTPUT
    
    for (int i=0; i<NUMNOD_SOTET4; i++) {
      bop(0,NODDOF_SOTET4*i+0) = defgrd(0,0)*N_XYZ(i,0);
      bop(0,NODDOF_SOTET4*i+1) = defgrd(1,0)*N_XYZ(i,0);
      bop(0,NODDOF_SOTET4*i+2) = defgrd(2,0)*N_XYZ(i,0);
      bop(1,NODDOF_SOTET4*i+0) = defgrd(0,1)*N_XYZ(i,1);
      bop(1,NODDOF_SOTET4*i+1) = defgrd(1,1)*N_XYZ(i,1);
      bop(1,NODDOF_SOTET4*i+2) = defgrd(2,1)*N_XYZ(i,1);
      bop(2,NODDOF_SOTET4*i+0) = defgrd(0,2)*N_XYZ(i,2);
      bop(2,NODDOF_SOTET4*i+1) = defgrd(1,2)*N_XYZ(i,2);
      bop(2,NODDOF_SOTET4*i+2) = defgrd(2,2)*N_XYZ(i,2);
      /* ~~~ */
      bop(3,NODDOF_SOTET4*i+0) = defgrd(0,0)*N_XYZ(i,1) + defgrd(0,1)*N_XYZ(i,0);
      bop(3,NODDOF_SOTET4*i+1) = defgrd(1,0)*N_XYZ(i,1) + defgrd(1,1)*N_XYZ(i,0);
      bop(3,NODDOF_SOTET4*i+2) = defgrd(2,0)*N_XYZ(i,1) + defgrd(2,1)*N_XYZ(i,0);
      bop(4,NODDOF_SOTET4*i+0) = defgrd(0,1)*N_XYZ(i,2) + defgrd(0,2)*N_XYZ(i,1);
      bop(4,NODDOF_SOTET4*i+1) = defgrd(1,1)*N_XYZ(i,2) + defgrd(1,2)*N_XYZ(i,1);
      bop(4,NODDOF_SOTET4*i+2) = defgrd(2,1)*N_XYZ(i,2) + defgrd(2,2)*N_XYZ(i,1);
      bop(5,NODDOF_SOTET4*i+0) = defgrd(0,2)*N_XYZ(i,0) + defgrd(0,0)*N_XYZ(i,2);
      bop(5,NODDOF_SOTET4*i+1) = defgrd(1,2)*N_XYZ(i,0) + defgrd(1,0)*N_XYZ(i,2);
      bop(5,NODDOF_SOTET4*i+2) = defgrd(2,2)*N_XYZ(i,0) + defgrd(2,0)*N_XYZ(i,2);
    }
    
  	#ifdef VERBOSE_OUTPUT
	cout << "bop\n" << bop;
	#endif //VERBOSE_OUTPUT
	
    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */

    Epetra_SerialDenseMatrix cmat(NUMSTR_SOTET4,NUMSTR_SOTET4);
    Epetra_SerialDenseVector stress(NUMSTR_SOTET4);
    
    so_tet4_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp);
	#ifdef VERBOSE_OUTPUT    
    cout << "material input\n";
   	#endif //VERBOSE_OUTPUT
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
    (*force).Multiply('T','N',detJ * (tet4_dis.weights)(gp),bop,stress,1.0);

    // integrate `elastic' and `initial-displacement' stiffness matrix
    // keu = keu + (B^T . C . B) * detJ * w(gp)
    Epetra_SerialDenseMatrix cb(NUMSTR_SOTET4,NUMDOF_SOTET4);
    cb.Multiply('N','N',1.0,cmat,bop,0.0);          // temporary C . B
    stiffmatrix->Multiply('T','N',detJ * (tet4_dis.weights)(gp),bop,cb,1.0);

    // intergrate `geometric' stiffness matrix and add to keu *****************
    Epetra_SerialDenseVector sfac(stress); // auxiliary integrated stress
    sfac.Scale(detJ * (tet4_dis.weights)(gp));     // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
    vector<double> SmB_L(NUMDIM_SOTET4);     // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod=0; inod<NUMNOD_SOTET4; ++inod){
      SmB_L[0] = sfac(0) * N_XYZ(inod,0) + sfac(3) * N_XYZ(inod,1) + sfac(5) * N_XYZ(inod,2);
      SmB_L[1] = sfac(3) * N_XYZ(inod,0) + sfac(1) * N_XYZ(inod,1) + sfac(4) * N_XYZ(inod,2);
      SmB_L[2] = sfac(5) * N_XYZ(inod,0) + sfac(4) * N_XYZ(inod,1) + sfac(2) * N_XYZ(inod,2);
      for (int jnod=0; jnod<NUMNOD_SOTET4; ++jnod){
        double bopstrbop = 0.0;            // intermediate value
        for (int idim=0; idim<NUMDIM_SOTET4; ++idim) bopstrbop += N_XYZ(jnod,idim)* SmB_L[idim];
        (*stiffmatrix)(NUMDIM_SOTET4*inod+0,NUMDIM_SOTET4*jnod+0) += bopstrbop;
        (*stiffmatrix)(NUMDIM_SOTET4*inod+1,NUMDIM_SOTET4*jnod+1) += bopstrbop;
        (*stiffmatrix)(NUMDIM_SOTET4*inod+2,NUMDIM_SOTET4*jnod+2) += bopstrbop;
      }
    } // end of intergrate `geometric' stiffness *******************************
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/
   
   
  // static integrator created in any case to safe "if-case" 
  const static DRT::ELEMENTS::Integrator_tet4_4point tet4_mass;
  if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
    //consistent mass matrix evaluated using a 4-point rule
    for (int gp=0; gp<tet4_mass.num_gp; gp++) {
      for (int inod=0; inod<NUMNOD_SOTET4; ++inod){
        for (int jnod=0; jnod<NUMNOD_SOTET4; ++jnod) {
          double massfactor = (tet4_mass.shapefct_gp[gp])(inod) * density *\
                              (tet4_mass.shapefct_gp[gp])(jnod) * detJ *\
                              (tet4_mass.weights)(gp);     // intermediate factor
          (*massmatrix)(NUMDIM_SOTET4*inod+0,NUMDIM_SOTET4*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4*inod+1,NUMDIM_SOTET4*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4*inod+2,NUMDIM_SOTET4*jnod+2) += massfactor;
        }
      }
    } 
  }// end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
    
  #ifdef VERBOSE_OUTPUT    
  cout << (*stiffmatrix);
  #endif //VERBOSE_OUTPUT
  return;
} // DRT::ELEMENTS::So_tet4::so_tet4_nlnstiffmass

int DRT::ELEMENTS::Sotet4Register::Initialize(DRT::Discretization& dis)
{
  //-------------------- loop all my column elements and check rewinding
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->Type() != DRT::Element::element_so_tet4) continue;
    DRT::ELEMENTS::So_tet4* actele = dynamic_cast<DRT::ELEMENTS::So_tet4*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4* failed");

    if (!actele->donerewinding_) {
      const bool rewind = DRT::UTILS::checkRewinding3D(actele);

      if (rewind) {
        int new_nodeids[NUMNOD_SOTET4];
        const int* old_nodeids;
        old_nodeids = actele->NodeIds();
        // rewinding of nodes to arrive at mathematically positive element
        new_nodeids[0] = old_nodeids[0];
        new_nodeids[1] = old_nodeids[2];
        new_nodeids[2] = old_nodeids[1];
        new_nodeids[3] = old_nodeids[3];
        actele->SetNodeIds(NUMNOD_SOTET4, new_nodeids);
      }
      // process of rewinding done
      actele->donerewinding_ = true;
    }
  }
  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  dis.FillComplete(false,false,false);

  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
