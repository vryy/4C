/*!----------------------------------------------------------------------*
\file so_hex8_evaluate.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by : Alexander Volf
			alexander.volf@mytum.de     
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOTET10
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_tet10.H"
#include "so_hex8.H"
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
 |                                                         vlf 06/07    |
 | cuts a block out of in_matrix                                        |
 |                              										|
 *----------------------------------------------------------------------*/
void cut_volf( 
     Epetra_SerialDenseMatrix& in_matrix,
	 int A_row,int A_col,
	 int B_row,int B_col,
	 Epetra_SerialDenseMatrix& out_matrix)
{
	out_matrix.Reshape(B_row-A_row,B_col-A_col);
	
	for (int i_row=0;i_row < out_matrix.M();i_row++)
	for (int i_col; i_col < out_matrix.N();i_col++)
	{
		out_matrix(i_row,i_col)= in_matrix(i_row+A_row,i_col+A_col);			
	}	
}

/*----------------------------------------------------------------------*
 |                                                         vlf 06/07    |
 | recursive computation of determinant of a  matrix using Sarrus rule  |
 |                              										|
 *----------------------------------------------------------------------*/
double det_volf(Epetra_SerialDenseMatrix& in_matrix)
{
	//Epetra_SerialDenseMatrix temp_matrix(in_matrix.N()-1,in_matrix.N()-1);
	
	if (in_matrix.N()==1)
	{
		return in_matrix(0,0);	
	}	
	else if (in_matrix.N()==2)
	{
		return 	((in_matrix(0,0)*in_matrix(1,1))-(in_matrix(0,1)*in_matrix(1,0)));
	}
	else if (in_matrix.N()>2)
	{
		double out_det=0;
		int sign=1;
		for (int i_col=0;i_col < in_matrix.N();i_col++)
		{
			
			Epetra_SerialDenseMatrix temp_matrix(in_matrix.N()-1,in_matrix.N()-1);
			for (int c_col=0;c_col < i_col;c_col++)
			{				
				for(int row=1;row<in_matrix.N();row++)
				temp_matrix(row-1,c_col)=in_matrix(row,c_col);							
			}
			for (int c_col=i_col+1;c_col <  in_matrix.N();c_col++)
			{
			for(int row=1;row<in_matrix.N();row++)
			temp_matrix(row-1,c_col-1)=in_matrix(row,c_col);	
			}
 
			out_det=out_det+(sign* in_matrix(0,i_col)*det_volf(temp_matrix));
			sign*=-1;
		}
		return out_det;
	}	
	else return 0;
}


/*----------------------------------------------------------------------***+
 |  evaluate the element (public)                              vlf 06/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::So_tet10::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DSTraceHelper dst("So_tet10::Evaluate");
  // start with "none"
  DRT::Elements::So_tet10::ActionType act = So_tet10::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = So_tet10::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = So_tet10::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = So_tet10::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = So_tet10::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = So_tet10::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = So_tet10::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = So_tet10::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = So_tet10::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = So_tet10::calc_struct_update_istep;
  else dserror("Unknown type of action for So_tet10");

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
      so_tet10_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,actmat);//*

    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff: {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);

      so_tet10_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,actmat);
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
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);

      so_tet10_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
    }
    break;

    // evaluate stresses
    case calc_struct_stress: {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);

      Epetra_SerialDenseMatrix stresses(NUMGPT_SOTET10,NUMSTR_SOTET10);//*
      so_tet10_stress(actmat,mydisp,&stresses); //*
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

    default:
      dserror("Unknown type of action for Solid3");
  }

  return 0;
}


/*----------------------------------------------------------------------**
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::So_tet10::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  DSTraceHelper dst("So_tet10::EvaluateNeumann");

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
    curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(time);
  // **
  
/* =============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
** =============================================================================*/
/* pointer to (static) array of shape function stored in Epetra_Vectors
 * evaluated at each gp , for each node
 * it looks like this:
 * shapefct[gp] = [N1(gp)  N2(gp)  -- N10(gp)]
 */
  Epetra_SerialDenseVector* shapefct; //[NUMNOD_SOTET10][NUMGPT_SOTET10]

/* pointer to (static) array of shape function derivatives stored in Epetra_Matrices
 * evaluated at each gp, for xsi1-4
 * it looks like this:
 *                [  N1_,xsi1   N1_,xsi2    N1_,xsi3   N1_,xsi4 ]
 * deriv_gp[gp] = [  N2_,xsi1   N2_,xsi2    N2_,xsi3   N2_,xsi4 ]
 *		          [     |		   |           |          |	    ]
 *		          [ N10_,xsi1  N10_,xsi2   N10_,xsi3  N10_,xsi4 ]   
 */
  Epetra_SerialDenseMatrix* deriv;    //[NUMGPT_SOTET10*NUMDIM][NUMNOD_SOTET10]

/* pointer to (static) weight factors at each gp */
  Epetra_SerialDenseVector* weights;  //[NUMGPT_SOTET10]
  so_tet10_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOTET10,NUMDIM_SOTET10);  // material coord. of element
  for (int i=0; i<NUMNOD_SOTET10; i++){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }

  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOTET10; gp++) {
    // get submatrix of deriv at actual gp
   
    /* get the matrix of the coordinates of edges needed to compute the volume,
    ** which is used here as detJ in the quadrature rule.  
    ** ("Jacobian matrix") for the quadrarture rule:
    **             [  1    1    1  	 1  ]
    ** jac_coord = [ x_1  x_2  x_3  x_4 ]
    **             [ y_1  y_2  y_3  y_4 ]
    **		       [ z_1  z_2  z_3  z_4 ]
    */
    Epetra_SerialDenseMatrix jac_coord(NUMCOORD_SOTET10,NUMCOORD_SOTET10);
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

    double fac = (*weights)(gp) * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
    for (int nodid=0; nodid<NUMNOD_SOH8; ++nodid) {
      for(int dim=0; dim<NUMDIM_SOH8; dim++) {
        elevec1[nodid+dim] += (*shapefct)(nodid,gp) * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }
  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::Elements::So_tet10::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            vlf 06/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::So_tet10::so_tet10_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      struct _MATERIAL*         material)       // element material data
{
  DSTraceHelper dst("So_tet10::so_tet10_nlnstiffmass");

/* =============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
** =============================================================================*/
/* pointer to (static) array of shape function stored in Epetra_Vectors
 * evaluated at each gp , for each node
 * it looks like this:
 * shapefct[gp] = [N1(gp)  N2(gp)  -- N10(gp)]
 */
  Epetra_SerialDenseVector* shapefct; //[NUMNOD_SOTET10][NUMGPT_SOTET10]

/* pointer to (static) array of shape function derivatives stored in Epetra_Matrices
 * evaluated at each gp, for xsi1-4
 * it looks like this:
 *                [  N1_,xsi1   N1_,xsi2    N1_,xsi3   N1_,xsi4 ]
 * deriv_gp[gp] = [  N2_,xsi1   N2_,xsi2    N2_,xsi3   N2_,xsi4 ]
 *		          [     |		   |           |          |	    ]
 *		          [ N10_,xsi1  N10_,xsi2   N10_,xsi3  N10_,xsi4 ]   
 */
  Epetra_SerialDenseMatrix* deriv;    //[NUMGPT_SOTET10*NUMDIM][NUMNOD_SOTET10]

/* pointer to (static) weight factors at each gp */
  Epetra_SerialDenseVector* weights;  //[NUMGPT_SOTET10]
  so_tet10_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOTET10,NUMDIM_SOTET10);  // material coord. of element
  /* structure of xrefe:
    **             [  X_1   Y_1   Z_1  ]
    **     xrefe = [  X_2   Y_2   Z_2  ]
    **             [   |     |     |   ]
    **             [  X_10  Y_10  Z_10 ]
    */
  Epetra_SerialDenseMatrix xcurr(NUMNOD_SOTET10,NUMDIM_SOTET10);  // current  coord. of element
  /* structure of xcurr:
    **             [  x_1   y_1   z_1  ]
    **     xcurr = [  x_2   y_2   z_2  ]
    **             [   |     |     |   ]
    **             [  x_10  y_10  z_10 ]
    */   
  for (int i=0; i<NUMNOD_SOTET10; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOTET10+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOTET10+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOTET10+2];
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
  
  Epetra_SerialDenseMatrix jac_coord(NUMCOORD_SOTET10,NUMCOORD_SOTET10);
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
  for (int gp=0; gp<NUMGPT_SOTET10; gp++) {
    
    /* compute the Jacobian matrix which looks like this:
    **         [  1        1        1  	     1      ]
    **   jac = [ X_,xsi1  X_,xsi2  X_,xsi3  X_,xsi4 ]
    **         [ Y_,xsi1  Y_,xsi2  Y_,xsi3  Y_,xsi4 ]
    **		   [ Z_,xsi1  Z_,xsi2  Z_,xsi3  Z_,xsi4 ]
    */
    
    Epetra_SerialDenseMatrix jac_temp(NUMCOORD_SOTET10-1,NUMCOORD_SOTET10);
    Epetra_SerialDenseMatrix jac(NUMCOORD_SOTET10,NUMCOORD_SOTET10);
    jac_temp.Multiply('T','N',1.0,xrefe,deriv[gp],0.0);
   
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
    
    Epetra_SerialDenseMatrix I_aug(NUMCOORD_SOTET10,NUMDIM_SOTET10);
    Epetra_SerialDenseMatrix partials(NUMCOORD_SOTET10,NUMDIM_SOTET10);
    Epetra_SerialDenseMatrix N_XYZ(NUMNOD_SOTET10,NUMDIM_SOTET10);
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

    N_XYZ.Multiply('N','N',1.0,deriv[gp],partials,0.0); //N_XYZ = N_xsi_k*partials
      
    /* structure of N_XYZ:
    **             [   dN_1     dN_1     dN_1   ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    **    N_XYZ =  [     |        |        |    ]
    **             [                            ]
    **             [   dN_10    dN_10    dN_10  ]
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

    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOTET10,NUMDIM_SOTET10);
    defgrd.Multiply('T','N',1.0,xcurr,N_XYZ,0.0);
    
    #ifdef VERBOSE_OUTPUT
	cout << "defgr\n " << defgrd;
	#endif //VERBOSE_OUTPUT
	
    // Right Cauchy-Green tensor = F^T * F 
    Epetra_SerialDenseMatrix cauchygreen(NUMDIM_SOTET10,NUMDIM_SOTET10);
   
    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,0.0);
    
	#ifdef VERBOSE_OUTPUT
	cout << "cauchygreen\n" << cauchygreen;
	getchar();
	#endif //VERBOSE_OUTPUT
	
    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain(NUMSTR_SOTET10);
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
    
    Epetra_SerialDenseMatrix bop(NUMSTR_SOTET10,NUMDOF_SOTET10);
    #ifdef VERBOSE_OUTPUT
    cout << bop;
    cout << defgrd;
    cout << N_XYZ;
    #endif //VERBOSE_OUTPUT
    
    for (int i=0; i<NUMNOD_SOTET10; i++) {
      bop(0,NODDOF_SOTET10*i+0) = defgrd(0,0)*N_XYZ(i,0);
      bop(0,NODDOF_SOTET10*i+1) = defgrd(1,0)*N_XYZ(i,0);
      bop(0,NODDOF_SOTET10*i+2) = defgrd(2,0)*N_XYZ(i,0);
      bop(1,NODDOF_SOTET10*i+0) = defgrd(0,1)*N_XYZ(i,1);
      bop(1,NODDOF_SOTET10*i+1) = defgrd(1,1)*N_XYZ(i,1);
      bop(1,NODDOF_SOTET10*i+2) = defgrd(2,1)*N_XYZ(i,1);
      bop(2,NODDOF_SOTET10*i+0) = defgrd(0,2)*N_XYZ(i,2);
      bop(2,NODDOF_SOTET10*i+1) = defgrd(1,2)*N_XYZ(i,2);
      bop(2,NODDOF_SOTET10*i+2) = defgrd(2,2)*N_XYZ(i,2);
      /* ~~~ */
      bop(3,NODDOF_SOTET10*i+0) = defgrd(0,0)*N_XYZ(i,1) + defgrd(0,1)*N_XYZ(i,0);
      bop(3,NODDOF_SOTET10*i+1) = defgrd(1,0)*N_XYZ(i,1) + defgrd(1,1)*N_XYZ(i,0);
      bop(3,NODDOF_SOTET10*i+2) = defgrd(2,0)*N_XYZ(i,1) + defgrd(2,1)*N_XYZ(i,0);
      bop(4,NODDOF_SOTET10*i+0) = defgrd(0,1)*N_XYZ(i,2) + defgrd(0,2)*N_XYZ(i,1);
      bop(4,NODDOF_SOTET10*i+1) = defgrd(1,1)*N_XYZ(i,2) + defgrd(1,2)*N_XYZ(i,1);
      bop(4,NODDOF_SOTET10*i+2) = defgrd(2,1)*N_XYZ(i,2) + defgrd(2,2)*N_XYZ(i,1);
      bop(5,NODDOF_SOTET10*i+0) = defgrd(0,2)*N_XYZ(i,0) + defgrd(0,0)*N_XYZ(i,2);
      bop(5,NODDOF_SOTET10*i+1) = defgrd(1,2)*N_XYZ(i,0) + defgrd(1,0)*N_XYZ(i,2);
      bop(5,NODDOF_SOTET10*i+2) = defgrd(2,2)*N_XYZ(i,0) + defgrd(2,0)*N_XYZ(i,2);
    }
    
  	#ifdef VERBOSE_OUTPUT
	cout << "bop\n" << bop;
	#endif //VERBOSE_OUTPUT
	
    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */

    Epetra_SerialDenseMatrix cmat(NUMSTR_SOTET10,NUMSTR_SOTET10);
    Epetra_SerialDenseVector stress(NUMSTR_SOTET10);
    double density;
    so_tet10_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp);
	#ifdef VERBOSE_OUTPUT    
    cout << "material input\n";
   	#endif //VERBOSE_OUTPUT
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
    (*force).Multiply('T','N',detJ * (*weights)(gp),bop,stress,1.0);

    // integrate `elastic' and `initial-displacement' stiffness matrix
    // keu = keu + (B^T . C . B) * detJ * w(gp)
    Epetra_SerialDenseMatrix cb(NUMSTR_SOTET10,NUMDOF_SOTET10);
    cb.Multiply('N','N',1.0,cmat,bop,0.0);          // temporary C . B
    stiffmatrix->Multiply('T','N',detJ * (*weights)(gp),bop,cb,1.0);

    // intergrate `geometric' stiffness matrix and add to keu *****************
    Epetra_SerialDenseVector sfac(stress); // auxiliary integrated stress
    sfac.Scale(detJ * (*weights)(gp));     // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
    vector<double> SmB_L(NUMDIM_SOTET10);     // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod=0; inod<NUMNOD_SOTET10; ++inod){
      SmB_L[0] = sfac(0) * N_XYZ(inod,0) + sfac(3) * N_XYZ(inod,1) + sfac(5) * N_XYZ(inod,2);
      SmB_L[1] = sfac(3) * N_XYZ(inod,0) + sfac(1) * N_XYZ(inod,1) + sfac(4) * N_XYZ(inod,2);
      SmB_L[2] = sfac(5) * N_XYZ(inod,0) + sfac(4) * N_XYZ(inod,1) + sfac(2) * N_XYZ(inod,2);
      for (int jnod=0; jnod<NUMNOD_SOTET10; ++jnod){
        double bopstrbop = 0.0;            // intermediate value
        for (int idim=0; idim<NUMDIM_SOTET10; ++idim) bopstrbop += N_XYZ(jnod,idim)* SmB_L[idim];
        (*stiffmatrix)(NUMDIM_SOTET10*inod+0,NUMDIM_SOTET10*jnod+0) += bopstrbop;
        (*stiffmatrix)(NUMDIM_SOTET10*inod+1,NUMDIM_SOTET10*jnod+1) += bopstrbop;
        (*stiffmatrix)(NUMDIM_SOTET10*inod+2,NUMDIM_SOTET10*jnod+2) += bopstrbop;
      }
    } // end of intergrate `geometric' stiffness ******************************
	
	//MASS matrix is yet to be implemented!!!!
    if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
      //integrate concistent mass matrix
      for (int inod=0; inod<NUMNOD_SOTET10; ++inod) {
        for (int jnod=0; jnod<NUMNOD_SOTET10; ++jnod) {
          double massfactor = (*shapefct)(inod,gp) * density * (*shapefct)(jnod,gp)
                            * detJ * (*weights)(gp);     // intermediate factor
          (*massmatrix)(NUMDIM_SOTET10*inod+0,NUMDIM_SOTET10*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET10*inod+1,NUMDIM_SOTET10*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET10*inod+2,NUMDIM_SOTET10*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
    
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/  
  #ifdef VERBOSE_OUTPUT    
  cout << (*stiffmatrix);
  #endif //VERBOSE_OUTPUT
  return;
} // DRT::Elements::So_tet10::so_tet10_nlnstiffmass


/*----------------------------------------------------------------------*
 | shape functions and derivatives for So_tet10               volf 06/07|
 | shape functions of a quadratic tetrahedra using so-called            |
 | "natural coordinates" as described by Carlos A. Felippa in Adv. FEM  |
 | Aerospace Engineering Sciences - University of Colorado at Boulder   |
 *----------------------------------------------------------------------*/
void DRT::Elements::So_tet10::so_tet10_shapederiv(
      Epetra_SerialDenseVector** shapefct_ptr,  // pointer to pointer of shapefct matrix filed
      Epetra_SerialDenseMatrix** deriv_gp_ptr,  // pointer to pointer of deriv_gp matrix field
      Epetra_SerialDenseVector** weights)       // pointer to pointer of weights       			
{
  DSTraceHelper dst("So_tet10::so_tet10_shapederiv");

  // static matrix objects, kept in memory
  static Epetra_SerialDenseVector weightfactors(NUMGPT_SOTET10);   // weights for each gp
  static Epetra_SerialDenseMatrix deriv_gp[NUMGPT_SOTET10];
  static Epetra_SerialDenseVector shape_gp[NUMGPT_SOTET10];
  static bool fdf_eval;                      // flag for re-evaluate everything

  Epetra_SerialDenseMatrix df(NUMNOD_SOTET10*NUMCOORD_SOTET10,NUMGPT_SOTET10);  // derivatives

  //Quadrature rule from Carlos A. Felippa: Adv. FEM  §16.4 
  const double gploc_alpha    = (5.0 + 3.0*sqrt(5))/20.0;    // gp sampling point value for quadr. fct
  const double gploc_beta     = (5.0 - sqrt(5))/20.0;
  const double gpw      = 0.25;              // weight at every gp for linear fct
 
  if (fdf_eval==true){          // if true f,df already evaluated
    *weights  = &weightfactors; // return adress of static object to target of pointer
    *deriv_gp_ptr = deriv_gp;   // return adress of static object to target of pointer
    *shapefct_ptr = shape_gp;
    return;
  } else {
    // (xsi1, xsi2, xsi3 ,xsi4) gp-locations of fully integrated linear 10-node Tet
    const double xsi1[NUMGPT_SOTET10] = {gploc_alpha, gploc_beta , gploc_beta , gploc_beta };
    const double xsi2[NUMGPT_SOTET10] = {gploc_beta , gploc_alpha, gploc_beta , gploc_beta };
    const double xsi3[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_alpha, gploc_beta };
    const double xsi4[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_beta , gploc_alpha};
    const double w[NUMGPT_SOTET10]    = {   gpw,   gpw,   gpw,   gpw};

    // fill up nodal f at each gp 
    for (int i=0; i<NUMGPT_SOTET10; i++) {
      (shape_gp[i]).Size(NUMNOD_SOTET10);
      (shape_gp[i])[0] = xsi1[i] * (2*xsi1[i] -1);
      (shape_gp[i])[1] = xsi2[i] * (2*xsi2[i] -1);
      (shape_gp[i])[2] = xsi3[i] * (2*xsi3[i] -1);
      (shape_gp[i])[3] = xsi4[i] * (2*xsi4[i] -1);
      (shape_gp[i])[4] = 4 * xsi1[i] * xsi2[i];
      (shape_gp[i])[5] = 4 * xsi2[i] * xsi3[i];
      (shape_gp[i])[6] = 4 * xsi3[i] * xsi1[i];
      (shape_gp[i])[7] = 4 * xsi1[i] * xsi4[i];  
      (shape_gp[i])[8] = 4 * xsi2[i] * xsi4[i];
      (shape_gp[i])[9] = 4 * xsi3[i] * xsi4[i];
      weightfactors[i] = w[i]; // just for clarity how to get weight factors
    }

    // fill up df xsi1, xsi2, xsi3, xsi4 directions (NUMDIM) at each gp
    for (int i=0; i<NUMGPT_SOTET10; i++) {
        // df wrt to xsi1 "+0" for each node(0..9) at each gp [i]
        df(NUMNOD_SOTET10*0+0,i) = 4 * xsi1[i]-1;
        df(NUMNOD_SOTET10*0+1,i) = 0;
      	df(NUMNOD_SOTET10*0+2,i) = 0;
      	df(NUMNOD_SOTET10*0+3,i) = 0;
    	
     	  df(NUMNOD_SOTET10*0+4,i) = 4 * xsi2[i];
      	df(NUMNOD_SOTET10*0+5,i) = 0;
      	df(NUMNOD_SOTET10*0+6,i) = 4 * xsi3[i];
      	df(NUMNOD_SOTET10*0+7,i) = 4 * xsi4[i];  
      	df(NUMNOD_SOTET10*0+8,i) = 0;
      	df(NUMNOD_SOTET10*0+9,i) = 0;

        // df wrt to xsi2 "+1" for each node(0..9) at each gp [i]
        df(NUMNOD_SOTET10*1+0,i) = 0;
      	df(NUMNOD_SOTET10*1+1,i) = 4 * xsi2[i] - 1;
      	df(NUMNOD_SOTET10*1+2,i) = 0;
      	df(NUMNOD_SOTET10*1+3,i) = 0;
      	
      	df(NUMNOD_SOTET10*1+4,i) = 4 * xsi1[i];
      	df(NUMNOD_SOTET10*1+5,i) = 4 * xsi3[i];
      	df(NUMNOD_SOTET10*1+6,i) = 0;
      	df(NUMNOD_SOTET10*1+7,i) = 0;  
      	df(NUMNOD_SOTET10*1+8,i) = 4 * xsi4[i];
      	df(NUMNOD_SOTET10*1+9,i) = 0;

        // df wrt to xsi3 "+2" for each node(0..9) at each gp [i]
        df(NUMNOD_SOTET10*2+0,i) = 0;
      	df(NUMNOD_SOTET10*2+1,i) = 0;
      	df(NUMNOD_SOTET10*2+2,i) = 4 * xsi3[i] - 1;
      	df(NUMNOD_SOTET10*2+3,i) = 0;
      
      	df(NUMNOD_SOTET10*2+4,i) = 0;
      	df(NUMNOD_SOTET10*2+5,i) = 4 * xsi2[i];
      	df(NUMNOD_SOTET10*2+6,i) = 4 * xsi1[i];
      	df(NUMNOD_SOTET10*2+7,i) = 0;  
      	df(NUMNOD_SOTET10*2+8,i) = 0;
      	df(NUMNOD_SOTET10*2+9,i) = 4 * xsi4[i];
        
        // df wrt to xsi4 "+2" for each node(0..9) at each gp [i]
        df(NUMNOD_SOTET10*3+0,i) = 0;
      	df(NUMNOD_SOTET10*3+1,i) = 0;
      	df(NUMNOD_SOTET10*3+2,i) = 0;
      	df(NUMNOD_SOTET10*3+3,i) = 4 * xsi4[i] - 1;
      
      	df(NUMNOD_SOTET10*3+4,i) = 0;
      	df(NUMNOD_SOTET10*3+5,i) = 0;
      	df(NUMNOD_SOTET10*3+6,i) = 0;
      	df(NUMNOD_SOTET10*3+7,i) = 4 * xsi1[i];  
      	df(NUMNOD_SOTET10*3+8,i) = 4 * xsi2[i];
      	df(NUMNOD_SOTET10*3+9,i) = 4 * xsi3[i];
    }
    
    for (int gp=0; gp<NUMGPT_SOTET10; gp++) {
       // get submatrix of deriv at actual gp
    
       /* compute the deriv_gp which looks like this:
       **             [  N1_,xsi1   N1_,xsi2    N1_,xsi3   N1_,xsi4 ]
       **  deriv_gp = [  N2_,xsi1   N2_,xsi2    N2_,xsi3   N2_,xsi4 ]
       **		      [     |		   |           |          |	    ]
       **		      [ N10_,xsi1  N10_,xsi2   N10_,xsi3  N10_,xsi4 ]
       */    	
    	(deriv_gp[gp]).Shape(NUMNOD_SOTET10,NUMCOORD_SOTET10);
    	for (int shape_func=0; shape_func<NUMNOD_SOTET10; shape_func++) {
      		for (int xsi_num=0;xsi_num<NUMCOORD_SOTET10; xsi_num++) {
        		(deriv_gp[gp])(shape_func,xsi_num)= df(NUMNOD_SOTET10*xsi_num+shape_func,gp);
      		}
    	}
    }

    // return adresses of just evaluated matrices
    *weights  = &weightfactors; // return adress of static object to target of pointer
    *deriv_gp_ptr = deriv_gp;	// return adress of static object array 
    *shapefct_ptr = shape_gp;   // return adress of static object array 
    fdf_eval = true;            // now all arrays are filled statically
  }
  return;
}  // So_tet10::sotet10_shapederiv

int DRT::Elements::Sotet10Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOTET10
