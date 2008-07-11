/*!----------------------------------------------------------------------*
\file so_tet10_evaluate.cpp
\brief quadratic nonlinear tetrahedron

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by: Alexander Volf
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
#include "so_integrator.H"
#include "so_weg6.H"
#include "so_hex8.H"
#include "so_tet10.H"
#include "so_tet4.H"
#include "so_disp.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"

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
int DRT::ELEMENTS::So_tet10::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // start with "none"
  DRT::ELEMENTS::So_tet10::ActionType act = So_tet10::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = So_tet10::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = So_tet10::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = So_tet10::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = So_tet10::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = So_tet10::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = So_tet10::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")        act = So_tet10::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = So_tet10::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = So_tet10::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = So_tet10::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = So_tet10::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")   act = So_tet10::calc_struct_reset_istep;
  else if (action=="postprocess_stress")        act = So_tet10::postprocess_stress;
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
      so_tet10_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat);//*

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

      so_tet10_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce: {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      Epetra_SerialDenseMatrix myemat(lm.size(),lm.size());

      so_tet10_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,actmat);
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

      so_tet10_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat);

      if (act==calc_struct_nlnstifflmass) so_tet10_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:{
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
      RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      if (stressdata==null) dserror("Cannot get stress 'data'");
      if (straindata==null) dserror("Cannot get strain 'data'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      Epetra_SerialDenseMatrix stress(NUMGPT_SOTET10,NUMSTR_SOTET10);
      Epetra_SerialDenseMatrix strain(NUMGPT_SOTET10,NUMSTR_SOTET10);
      bool cauchy = params.get<bool>("cauchy", false);
      string iostrain = params.get<string>("iostrain", "none");
      if (iostrain!="euler_almansi") so_tet10_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,actmat,cauchy);
      else dserror("requested option not yet implemented for tet10");
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:{

      const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      RCP<Epetra_SerialDenseMatrix> gpstress = (*gpstressmap)[gid];

      if (stresstype=="ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_SOTET10,NUMSTR_SOTET10);
        so_tet10_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOTET10);

        for (int i=0;i<NUMNOD_SOTET10;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_SOTET10;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_SOTET10;++i){
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
      }
      else if (stresstype=="cxyz") {
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_SOTET10; ++i) {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_SOTET10; ++j) {
              (*((*elestress)(i)))[lid] += 1.0/NUMGPT_SOTET10 * (*gpstress)(j,i);
            }
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_SOTET10,NUMSTR_SOTET10);
        so_tet10_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOTET10);

        for (int i=0;i<NUMNOD_SOTET10;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_SOTET10;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_SOTET10;++i){
          elevec2(3*i)=nodalstresses(i,3)/numadjele[i];
          elevec2(3*i+1)=nodalstresses(i,4)/numadjele[i];
          elevec2(3*i+2)=nodalstresses(i,5)/numadjele[i];
        }
        RCP<Epetra_MultiVector> elestress=params.get<RCP<Epetra_MultiVector> >("elestress",null);
        if (elestress==null)
          dserror("No element stress/strain vector available");
        const Epetra_BlockMap elemap = elestress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1) {
          for (int i = 0; i < NUMSTR_SOTET10; ++i) {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_SOTET10; ++j) {
              (*((*elestress)(i)))[lid] += 1.0/NUMGPT_SOTET10 * (*gpstress)(j,i);
            }
          }
        }
      }
      else {
        dserror("unknown type of stress/strain output on element level");
      }
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

    case calc_struct_update_imrlike: {
      ;// there is nothing to do here at the moment
    }
    break;

    case calc_struct_reset_istep: {
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
int DRT::ELEMENTS::So_tet10::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
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
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
  // **

/* =============================================================================*
 * CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
 * =============================================================================*/
   const static DRT::ELEMENTS::Integrator_tet10_4point tet10_dis;
/* ============================================================================*/

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
    double detJ=(double) det_volf(jac_coord);
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    double fac = tet10_dis.weights(gp) * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
    for (int nodid=0; nodid<NUMNOD_SOTET10; ++nodid) {
      for(int dim=0; dim<NUMDIM_SOTET10; dim++) {
        elevec1[nodid*NUMDIM_SOTET10+dim] += (tet10_dis.shapefct_gp[gp])(nodid) * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }
  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::So_tet10::EvaluateNeumann



/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            vlf 06/07 |
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::So_tet10::so_tet10_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    // element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     // element mass matrix
      Epetra_SerialDenseVector* force,          // element internal force vector
      Epetra_SerialDenseMatrix* elestress,      // stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      // strains at GP
      struct _MATERIAL*         material,       // element material data
      const bool                cauchy)         // stress output option
{
/* =============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
** =============================================================================*/
  const static DRT::ELEMENTS::Integrator_tet10_4point tet10_dis;
/* =============================================================================*/
  double density; //forward declaration for mass matrix
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
  Epetra_SerialDenseMatrix xdisp(NUMNOD_SOTET10,NUMDIM_SOTET10);  // current  coord. of element

  for (int i=0; i<NUMNOD_SOTET10; ++i){
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOTET10+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOTET10+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOTET10+2];

    xdisp(i,0) = disp[i*NODDOF_SOTET10+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET10+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET10+2];
  }


  /* get the matrix of the coordinates of edges needed to compute the volume,
  ** which is used here as detJ in the quadrature rule.
  ** ("Jacobian matrix") for the quadrarture rule:
  **             [  1    1    1    1  ]
  ** jac_coord = [ X_1  X_2  X_3  X_4 ]
  **             [ Y_1  Y_2  Y_3  Y_4 ]
  **             [ Z_1  Z_2  Z_3  Z_4 ]
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

  double detJ=(double)( (long double)det_volf(jac_coord)/(long double) 6);    //volume of a tetrahedra
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
    **         [ Z_,xsi1  Z_,xsi2  Z_,xsi3  Z_,xsi4 ]
    */

    Epetra_SerialDenseMatrix jac_temp(NUMCOORD_SOTET10-1,NUMCOORD_SOTET10);
    Epetra_SerialDenseMatrix jac(NUMCOORD_SOTET10,NUMCOORD_SOTET10);
    jac_temp.Multiply('T','N',1.0,xrefe,tet10_dis.deriv_gp[gp],0.0);

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

    N_XYZ.Multiply('N','N',1.0,tet10_dis.deriv_gp[gp],partials,0.0); //N_XYZ = N_xsi_k*partials

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
    defgrd.Multiply('T','N',1.0,xdisp,N_XYZ,0.0);
    defgrd(0,0)+=1;
    defgrd(1,1)+=1;
    defgrd(2,2)+=1;

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

    // return gp strains (only in case of stress/strain output)
   if (elestrain != NULL){
      for (int i = 0; i < 3; ++i) {
        (*elestrain)(gp,i) = glstrain(i);
      }
      for (int i = 3; i < 6; ++i) {
        (*elestrain)(gp,i) = 0.5 * glstrain(i);
      }
    }


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

    for (int numnode=0; numnode<NUMNOD_SOTET10; numnode++) {
    	for (int numdof=0; numdof<NODDOF_SOTET10; numdof++) {
      	bop(0,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,0)*N_XYZ(numnode,0);
      	bop(1,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,1)*N_XYZ(numnode,1);
      	bop(2,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,2)*N_XYZ(numnode,2);
      	bop(3,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,0)*N_XYZ(numnode,1) + \
      			    						   defgrd(numdof,1)*N_XYZ(numnode,0);
      	bop(4,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,1)*N_XYZ(numnode,2) + \
      										   defgrd(numdof,2)*N_XYZ(numnode,1);
      	bop(5,NODDOF_SOTET10*numnode+numdof) = defgrd(numdof,2)*N_XYZ(numnode,0) + \
      										   defgrd(numdof,0)*N_XYZ(numnode,2);
    	}
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

    so_tet10_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp);
	#ifdef VERBOSE_OUTPUT
    cout << "material input\n";
   	#endif //VERBOSE_OUTPUT
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    if (elestress != NULL){                // return 2nd Piola-Kirchhoff stresses
      if (!cauchy) {
        for (int i = 0; i < NUMSTR_SOTET10; ++i) {
          (*elestress)(gp,i) = stress(i);
        }
      }
      else {                               // return Cauchy stresses
        double detF = defgrd(0,0)*defgrd(1,1)*defgrd(2,2) +
                      defgrd(0,1)*defgrd(1,2)*defgrd(2,0) +
                      defgrd(0,2)*defgrd(1,0)*defgrd(2,1) -
                      defgrd(0,2)*defgrd(1,1)*defgrd(2,0) -
                      defgrd(0,0)*defgrd(1,2)*defgrd(2,1) -
                      defgrd(0,1)*defgrd(1,0)*defgrd(2,2);

        LINALG::SerialDenseMatrix pkstress(NUMDIM_SOTET10,NUMDIM_SOTET10);
        pkstress(0,0) = stress(0);
        pkstress(0,1) = stress(3);
        pkstress(0,2) = stress(5);
        pkstress(1,0) = pkstress(0,1);
        pkstress(1,1) = stress(1);
        pkstress(1,2) = stress(4);
        pkstress(2,0) = pkstress(0,2);
        pkstress(2,1) = pkstress(1,2);
        pkstress(2,2) = stress(2);

        LINALG::SerialDenseMatrix temp(NUMDIM_SOTET10,NUMDIM_SOTET10);
        LINALG::SerialDenseMatrix cauchystress(NUMDIM_SOTET10,NUMDIM_SOTET10);
        temp.Multiply('N','N',1.0/detF,defgrd,pkstress,0.);
        cauchystress.Multiply('N','T',1.0,temp,defgrd,0.);

        (*elestress)(gp,0) = cauchystress(0,0);
        (*elestress)(gp,1) = cauchystress(1,1);
        (*elestress)(gp,2) = cauchystress(2,2);
        (*elestress)(gp,3) = cauchystress(0,1);
        (*elestress)(gp,4) = cauchystress(1,2);
        (*elestress)(gp,5) = cauchystress(0,2);
      }
    }

    if (force != NULL && stiffmatrix != NULL) {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).Multiply('T','N',detJ * tet10_dis.weights[gp],bop,stress,1.0);

      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      Epetra_SerialDenseMatrix cb(NUMSTR_SOTET10,NUMDOF_SOTET10);
      cb.Multiply('N','N',1.0,cmat,bop,0.0);          // temporary C . B
      stiffmatrix->Multiply('T','N',detJ * tet10_dis.weights(gp),bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      Epetra_SerialDenseVector sfac(stress); // auxiliary integrated stress
      sfac.Scale(detJ * tet10_dis.weights(gp));     // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
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
    }
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/
   	//Epetra_SerialDenseVector L((*stiffmatrix).N());
 //  cout << (*stiffmatrix).N() << endl;
  // LINALG::SymmetricEigen((*stiffmatrix),L, (*stiffmatrix).N());
   //cout << L;
   //getchar();
   	//my personal output routine, identifies a
	/*static DRT::Element* my_static = this;
	static int is_static=0;
	static int my_I=0;
	if (is_static!=0)
	for (int x=0;x<6;x++){
		if (Nodes()[x]->X()[1]== 1){
			 my_static = this;
			 is_static=1;
			 my_I=x;
		}
	}

    if(this ==my_static ) {
    	cout << disp[my_I*NODDOF_SOTET10+0] << " ";
    	cout << disp[my_I*NODDOF_SOTET10+1] << " ";
    	cout << disp[my_I*NODDOF_SOTET10+2] << endl << endl;
    }
  */
/* =============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 14 GAUSS POINTS*
** =============================================================================*/
  const static DRT::ELEMENTS::Integrator_tet10_14point tet10_mass;
/* ============================================================================*/
  if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
    for (int gp=0; gp<tet10_mass.num_gp; gp++) {
      //integrate concistent mass matrix
      Epetra_SerialDenseMatrix jac(NUMDIM_SOTET10,NUMDIM_SOTET10);
      jac.Multiply('N','N',1.0,tet10_mass.deriv_gp[gp],xrefe,0.0);

      // compute determinant of Jacobian by Sarrus' rule
      double detJ2= jac(0,0) * jac(1,1) * jac(2,2)
                 + jac(0,1) * jac(1,2) * jac(2,0)
                 + jac(0,2) * jac(1,0) * jac(2,1)
                 - jac(0,0) * jac(1,2) * jac(2,1)
                 - jac(0,1) * jac(1,0) * jac(2,2)
                 - jac(0,2) * jac(1,1) * jac(2,0);
      if (abs(detJ2) < 1E-16) dserror("ZERO JACOBIAN DETERMINANT");
      else if (detJ2 < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");
      for (int inod=0; inod<NUMNOD_SOTET10; ++inod) {
        for (int jnod=0; jnod<NUMNOD_SOTET10; ++jnod) {
          double massfactor = (tet10_mass.shapefct_gp[gp])(inod) * density *\
          				            (tet10_mass.shapefct_gp[gp])(jnod)  * detJ2 *\
          				            (tet10_mass.weights)(gp);    // intermediate factor
          (*massmatrix)(NUMDIM_SOTET10*inod+0,NUMDIM_SOTET10*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET10*inod+1,NUMDIM_SOTET10*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET10*inod+2,NUMDIM_SOTET10*jnod+2) += massfactor;
        }
      }
    }
//    double masssum = 0.0;
//    for (int i=0; i<NUMDOF_SOTET10; ++i){
//      for (int j=0; j<NUMDOF_SOTET10; ++j) {
//        masssum += (*massmatrix)(i,j);
//      }
//    }
//    cout << "masssum: " << masssum << endl;
  } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
  #ifdef VERBOSE_OUTPUT
  cout << (*stiffmatrix);
  #endif //VERBOSE_OUTPUT
  return;
} // DRT::ELEMENTS::So_tet10::so_tet10_nlnstiffmass

/*----------------------------------------------------------------------*
 |  lump mass matrix                                         bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_lumpmass(Epetra_SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (int c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;  
      for (int r=0; r<(*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r,c);  // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
}


int DRT::ELEMENTS::Sotet10Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
