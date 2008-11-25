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
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "so_tet10.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
//#include "Epetra_SerialDenseSolver.h"

//#define VERBOSE_OUTPUT
using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |                                                         maf 04/07    |
 | vector of material laws                                              |
 | defined in global_control.c											|
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;  ///< C-style material struct


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              vlf 06/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_SOTET10,NUMDOF_SOTET10> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOTET10,NUMDOF_SOTET10> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOTET10,1>              elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOTET10,1>              elevec2(elevec2_epetra.A(),true);
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
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (unsigned i=0; i<mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (unsigned i=0; i<myres.size(); ++i) myres[i] = 0.0;
      so_tet10_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat,
                            INPAR::STR::stress_none,INPAR::STR::strain_none);//*

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
      so_tet10_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat,
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
      LINALG::Matrix<NUMDOF_SOTET10,NUMDOF_SOTET10> myemat(true); // set to zero
      so_tet10_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,actmat,
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
      so_tet10_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat,
                            INPAR::STR::stress_none,INPAR::STR::strain_none);

      if (act==calc_struct_nlnstifflmass) so_tet10_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
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
      LINALG::Matrix<NUMGPT_SOTET10,NUMSTR_SOTET10> stress;
      LINALG::Matrix<NUMGPT_SOTET10,NUMSTR_SOTET10> strain;
      INPAR::STR::StressType iostress = params.get<INPAR::STR::StressType>("iostress", INPAR::STR::stress_none);
      INPAR::STR::StrainType iostrain = params.get<INPAR::STR::StrainType>("iostrain", INPAR::STR::strain_none);
      so_tet10_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,actmat,iostress,iostrain);
      AddtoPack(*stressdata, stress);
      AddtoPack(*straindata, strain);
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_SOTET10,NUMSTR_SOTET10> gpstress(((*gpstressmap)[gid])->A(),true);

      if (stresstype=="ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::Matrix<NUMNOD_SOTET10,NUMSTR_SOTET10> nodalstresses;
        so_tet10_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOTET10);

        DRT::Node** nodes = Nodes();
        for (int i=0;i<NUMNOD_SOTET10;++i){
          DRT::Node* node=nodes[i];
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
            double& s = (*((*elestress)(i)))[lid];
            s = 0.;
            for (int j = 0; j < NUMGPT_SOTET10; ++j) {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOTET10;
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        LINALG::Matrix<NUMNOD_SOTET10,NUMSTR_SOTET10> nodalstresses;
        so_tet10_expol(gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOTET10);

        DRT::Node** nodes = Nodes();
        for (int i=0;i<NUMNOD_SOTET10;++i){
          DRT::Node* node=nodes[i];
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
            double& s = (*((*elestress)(i)))[lid];
            s = 0.;
            for (int j = 0; j < NUMGPT_SOTET10; ++j) {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOTET10;
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
  const vector<LINALG::Matrix<NUMNOD_SOTET10,1> >& shapefcts = so_tet10_4gp_shapefcts();
  const vector<double>& gpweights = so_tet10_4gp_weights();
/* ============================================================================*/

/* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOTET10; gp++) {
    // get submatrix of deriv at actual gp

    /* get the matrix of the coordinates of edges needed to compute the volume,
    ** which is used here as detJ in the quadrature rule.
    ** ("Jacobian matrix") for the quadrarture rule:
    **             [  1    1    1    1  ]
    ** jac_coord = [ x_1  x_2  x_3  x_4 ]
    **             [ y_1  y_2  y_3  y_4 ]
    **		   [ z_1  z_2  z_3  z_4 ]
    */
    LINALG::Matrix<NUMCOORD_SOTET10,NUMCOORD_SOTET10> jac_coord;
    for (int i=0; i<4; i++) jac_coord(0,i)=1;
    for (int col=0;col<4;col++)
    {
      const double* x = Nodes()[col]->X();
      for (int row=0;row<3;row++)
          jac_coord(row+1,col)= x[row];
    }

    // compute determinant of Jacobian with own algorithm
    // !!warning detJ is not the actual determinant of the jacobian (here needed for the quadrature rule)
    // but rather the volume of the tetrahedara
    double detJ=jac_coord.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    double fac = gpweights[gp] * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
    for(int dim=0; dim<NUMDIM_SOTET10; dim++) {
      double dim_fac = (*onoff)[dim] * (*val)[dim] * fac;
      for (int nodid=0; nodid<NUMNOD_SOTET10; ++nodid) {
        elevec1(nodid*NUMDIM_SOTET10+dim) += (shapefcts[gp])(nodid) * dim_fac;
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
      vector<double>&           residual,       // current residual displ
      LINALG::Matrix<NUMDOF_SOTET10,NUMDOF_SOTET10>* stiffmatrix,    // element stiffness matrix
      LINALG::Matrix<NUMDOF_SOTET10,NUMDOF_SOTET10>* massmatrix,     // element mass matrix
      LINALG::Matrix<NUMDOF_SOTET10,1>* force,          // element internal force vector
      LINALG::Matrix<NUMGPT_SOTET10,NUMSTR_SOTET10>* elestress,      // stresses at GP
      LINALG::Matrix<NUMGPT_SOTET10,NUMSTR_SOTET10>* elestrain,      // strains at GP
      struct _MATERIAL*         material,       // element material data
      const INPAR::STR::StressType                   iostress,       // stress output option
      const INPAR::STR::StrainType                   iostrain)       // strain output option
{
/* =============================================================================*
** CONST DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
** =============================================================================*/
  const vector<LINALG::Matrix<NUMNOD_SOTET10,NUMCOORD_SOTET10> >& derivs = so_tet10_4gp_derivs();
  const vector<double>& gpweights = so_tet10_4gp_weights();
/* =============================================================================*/
  double density; //forward declaration for mass matrix
  // update element geometry
  LINALG::Matrix<NUMNOD_SOTET10,NUMDIM_SOTET10> xrefe;  // material coord. of element
  /* structure of xrefe:
    **             [  X_1   Y_1   Z_1  ]
    **     xrefe = [  X_2   Y_2   Z_2  ]
    **             [   |     |     |   ]
    **             [  X_10  Y_10  Z_10 ]
    */
  LINALG::Matrix<NUMNOD_SOTET10,NUMDIM_SOTET10> xcurr;  // current  coord. of element
  /* structure of xcurr:
    **             [  x_1   y_1   z_1  ]
    **     xcurr = [  x_2   y_2   z_2  ]
    **             [   |     |     |   ]
    **             [  x_10  y_10  z_10 ]
    */
  LINALG::Matrix<NUMNOD_SOTET10,NUMDIM_SOTET10> xdisp;  // current  coord. of element

  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOTET10; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

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

  LINALG::Matrix<NUMCOORD_SOTET10,NUMCOORD_SOTET10> jac; // this was jac_coord
  for (int i=0; i<4; i++)  jac(0,i)=1;
  for (int row=0;row<3;row++)
  {
    for (int col=0;col<4;col++)
      jac(row+1,col)= xrefe(col,row);
  }

  double detJ=jac.Determinant()/6;    //volume of a tetrahedra
  if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
  else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");



  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  LINALG::Matrix<NUMCOORD_SOTET10-1,NUMCOORD_SOTET10> jac_temp;
  //Epetra_SerialDenseMatrix jac(NUMCOORD_SOTET10,NUMCOORD_SOTET10);
  for (int gp=0; gp<NUMGPT_SOTET10; gp++) {

    /* compute the Jacobian matrix which looks like this:
    **         [  1        1        1  	     1      ]
    **   jac = [ X_,xsi1  X_,xsi2  X_,xsi3  X_,xsi4 ]
    **         [ Y_,xsi1  Y_,xsi2  Y_,xsi3  Y_,xsi4 ]
    **         [ Z_,xsi1  Z_,xsi2  Z_,xsi3  Z_,xsi4 ]
    */

    jac_temp.MultiplyTN(xrefe,derivs[gp]);
    //multiply<NUMCOORD_SOTET10-1,NUMNOD_SOTET10,NUMCOORD_SOTET10,'T','N'>
    //(jac_temp,xrefe,tet10_dis.deriv_gp[gp]);

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

    LINALG::Matrix<NUMCOORD_SOTET10,NUMDIM_SOTET10> I_aug(true);
    LINALG::Matrix<NUMCOORD_SOTET10,NUMDIM_SOTET10> partials;
    LINALG::Matrix<NUMNOD_SOTET10,NUMDIM_SOTET10> N_XYZ;
    I_aug(1,0)=1;
    I_aug(2,1)=1;
    I_aug(3,2)=1;

    #ifdef VERBOSE_OUTPUT
    cout << "jac\n" << jac;
    #endif //VERBOSE_OUTPUT

    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<NUMCOORD_SOTET10,NUMCOORD_SOTET10,NUMDIM_SOTET10> solve_for_inverseJac;
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    solve_for_inverseJac.SetVectors(partials,I_aug);// set X=partials, B=I_aug
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();         // partials = jac^-1.I_aug
    if ((err != 0) || (err2!=0)){
    	dserror("Inversion of Jacobian failed");
    }

    #ifdef VERBOSE_OUTPUT
    cout << "I_aug\n" << I_aug;
    cout << "partials\n" << partials;
    cout << "deriv_gp\n" << derivs[gp];
    #endif //VERBOSE_OUTPUT

    N_XYZ.Multiply(derivs[gp],partials); //N_XYZ = N_xsi_k*partials
    //multiply<NUMNOD_SOTET10,NUMCOORD_SOTET10,NUMDIM_SOTET10,'N','N'>(N_XYZ,tet10_dis.deriv_gp[gp],partials);

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

    LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10> defgrd(true);
    defgrd.MultiplyTN(xdisp,N_XYZ);
    //multiply<NUMDIM_SOTET10,NUMNOD_SOTET10,NUMDIM_SOTET10,'T','N'>(defgrd,xdisp,N_XYZ);
    defgrd(0,0)+=1;
    defgrd(1,1)+=1;
    defgrd(2,2)+=1;

    #ifdef VERBOSE_OUTPUT
    cout << "defgr\n " << defgrd;
    #endif //VERBOSE_OUTPUT

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10> cauchygreen;

    cauchygreen.MultiplyTN(defgrd,defgrd);
    //multiply<NUMDIM_SOTET10,NUMDIM_SOTET10,NUMDIM_SOTET10,'T','N'>(cauchygreen,defgrd,defgrd);

    #ifdef VERBOSE_OUTPUT
    cout << "cauchygreen\n" << cauchygreen;
    getchar();
    #endif //VERBOSE_OUTPUT

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain_epetra(NUMSTR_SOTET10);
    LINALG::Matrix<NUMSTR_SOTET10,1> glstrain(glstrain_epetra.A(),true);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    #ifdef VERBOSE_OUTPUT
    cout << "glstrain\n" << glstrain;
    #endif //VERBOSE_OUTPUT

    // return gp strains if necessary
    switch (iostrain)
    {
    case INPAR::STR::strain_gl:
    {
      if (elestrain == NULL) dserror("no strain data available");
      for (int i = 0; i < 3; ++i) {
        (*elestrain)(gp,i) = glstrain(i);
      }
      for (int i = 3; i < 6; ++i) {
        (*elestrain)(gp,i) = 0.5 * glstrain(i);
      }
    }
    break;
    case INPAR::STR::strain_ea:
      dserror("no Euler-Almansi strains available for tet10");
    break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror("requested strain option not available");
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

    LINALG::Matrix<NUMSTR_SOTET10,NUMDOF_SOTET10> bop;
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

    LINALG::Matrix<NUMSTR_SOTET10,NUMSTR_SOTET10> cmat(true);
    LINALG::Matrix<NUMSTR_SOTET10,1> stress(true);

    so_tet10_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp);
    #ifdef VERBOSE_OUTPUT
    cout << "material input\n";
    #endif //VERBOSE_OUTPUT
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses if necessary
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("no stress data available");
      for (int i = 0; i < NUMSTR_SOTET10; ++i) {
        (*elestress)(gp,i) = stress(i);
      }
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("no stress data available");
      double detF = defgrd.Determinant();

      LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10> pkstress;
      pkstress(0,0) = stress(0);
      pkstress(0,1) = stress(3);
      pkstress(0,2) = stress(5);
      pkstress(1,0) = pkstress(0,1);
      pkstress(1,1) = stress(1);
      pkstress(1,2) = stress(4);
      pkstress(2,0) = pkstress(0,2);
      pkstress(2,1) = pkstress(1,2);
      pkstress(2,2) = stress(2);

      LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10> temp;
      LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10> cauchystress;
      temp.Multiply(1.0/detF,defgrd,pkstress,0.);
      //multiply<NUMDIM_SOTET10,NUMDIM_SOTET10,NUMDIM_SOTET10,'N','N'>(0.0,temp,1.0/detF,defgrd,pkstress);
      cauchystress.MultiplyNT(temp,defgrd);
      //multiply<NUMDIM_SOTET10,NUMDIM_SOTET10,NUMDIM_SOTET10,'N','T'>(cauchystress,temp,defgrd);

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
      dserror("requested stress option not available");
    }

    if (force != NULL && stiffmatrix != NULL) {
      double detJ_w = detJ * gpweights[gp];
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w,bop,stress,1.0);
      //multiply<NUMDOF_SOTET10,NUMSTR_SOTET10,1,'T','N'>(1.0,*force,detJ_w,bop,stress);

      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<NUMSTR_SOTET10,NUMDOF_SOTET10> cb;
      cb.Multiply(cmat,bop);          // temporary C . B
      //multiply<NUMSTR_SOTET10,NUMSTR_SOTET10,NUMDOF_SOTET10,'N','N'>(cb,cmat,bop);
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);
      //multiply<NUMDOF_SOTET10,NUMSTR_SOTET10,NUMDOF_SOTET10,'T','N'>
      //  (1.0,*stiffmatrix,detJ_w,bop,cb);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<NUMSTR_SOTET10,1> sfac; // auxiliary integrated stress
      sfac.Update(detJ_w,stress); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
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
** CONST SHAPE FUNCTIONS and WEIGHTS for TET_10 with 10 GAUSS POINTS            *
** =============================================================================*/
  const vector<LINALG::Matrix<NUMNOD_SOTET10,1> >& shapefcts10gp = so_tet10_10gp_shapefcts();
  const vector<LINALG::Matrix<NUMDIM_SOTET10,NUMNOD_SOTET10> >& derivs10gp = so_tet10_10gp_derivs();
  const vector<double>& gpweights10gp = so_tet10_10gp_weights();
/* ============================================================================*/
  if (massmatrix != NULL){ // evaluate mass matrix +++++++++++++++++++++++++
    for (int gp=0; gp<10; gp++) {
      //integrate concistent mass matrix
      LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10> jac;
      jac.Multiply(derivs10gp[gp],xrefe);
      //multiply<NUMDIM_SOTET10,NUMNOD_SOTET10,NUMDIM_SOTET10,'N','N'>(jac,tet10_mass.deriv_gp[gp],xrefe);

      // compute determinant of Jacobian by Sarrus' rule
      double detJ2= jac.Determinant();
      if (abs(detJ2) < 1E-16) dserror("ZERO JACOBIAN DETERMINANT");
      else if (detJ2 < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");
      double factor = density * detJ2 * gpweights10gp[gp];
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOTET10; ++inod)
      {
        ifactor = (shapefcts10gp[gp])(inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOTET10; ++jnod)
        {
          massfactor = (shapefcts10gp[gp])(jnod) * ifactor;    // intermediate factor
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
 |  Evaluate Tet10 Shape fcts at 4 Gauss Points                         |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMNOD_SOTET10,1> >& DRT::ELEMENTS::So_tet10::so_tet10_4gp_shapefcts()
{
  static vector<LINALG::Matrix<NUMNOD_SOTET10,1> > shapefcts(NUMGPT_SOTET10);
  static bool shapefcts_done = false;
  if (shapefcts_done) return shapefcts;

  const double gploc_alpha    = (5.0 + 3.0*sqrt(5.0))/20.0;    // gp sampling point value for quadr. fct
  const double gploc_beta     = (5.0 - sqrt(5.0))/20.0;

  // (xsi1, xsi2, xsi3 ,xsi4) gp-locations of fully integrated linear 10-node Tet
  const double xsi1[NUMGPT_SOTET10] = {gploc_alpha, gploc_beta , gploc_beta , gploc_beta };
  const double xsi2[NUMGPT_SOTET10] = {gploc_beta , gploc_alpha, gploc_beta , gploc_beta };
  const double xsi3[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_alpha, gploc_beta };
  const double xsi4[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_beta , gploc_alpha};

  for (int gp=0; gp<NUMGPT_SOTET10; gp++) {
    (shapefcts[gp])(0) = xsi1[gp] * (2*xsi1[gp] -1);
    (shapefcts[gp])(1) = xsi2[gp] * (2*xsi2[gp] -1);
    (shapefcts[gp])(2) = xsi3[gp] * (2*xsi3[gp] -1);
    (shapefcts[gp])(3) = xsi4[gp] * (2*xsi4[gp] -1);
    (shapefcts[gp])(4) = 4 * xsi1[gp] * xsi2[gp];
    (shapefcts[gp])(5) = 4 * xsi2[gp] * xsi3[gp];
    (shapefcts[gp])(6) = 4 * xsi3[gp] * xsi1[gp];
    (shapefcts[gp])(7) = 4 * xsi1[gp] * xsi4[gp];
    (shapefcts[gp])(8) = 4 * xsi2[gp] * xsi4[gp];
    (shapefcts[gp])(9) = 4 * xsi3[gp] * xsi4[gp];
  }
  shapefcts_done = true;

  return shapefcts;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fct derivs at 4 Gauss Points                   |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMNOD_SOTET10,NUMCOORD_SOTET10> >& DRT::ELEMENTS::So_tet10::so_tet10_4gp_derivs()
{
  static vector<LINALG::Matrix<NUMNOD_SOTET10, NUMCOORD_SOTET10> > derivs(NUMGPT_SOTET10);
  static bool derivs_done = false;
  if (derivs_done) return derivs;

  const double gploc_alpha    = (5.0 + 3.0*sqrt(5.0))/20.0;    // gp sampling point value for quadr. fct
  const double gploc_beta     = (5.0 - sqrt(5.0))/20.0;

  // (xsi1, xsi2, xsi3 ,xsi4) gp-locations of fully integrated linear 10-node Tet
  const double xsi1[NUMGPT_SOTET10] = {gploc_alpha, gploc_beta , gploc_beta , gploc_beta };
  const double xsi2[NUMGPT_SOTET10] = {gploc_beta , gploc_alpha, gploc_beta , gploc_beta };
  const double xsi3[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_alpha, gploc_beta };
  const double xsi4[NUMGPT_SOTET10] = {gploc_beta , gploc_beta , gploc_beta , gploc_alpha};

  for (int gp=0; gp<NUMGPT_SOTET10; gp++) {
    (derivs[gp])(0,0) = 4 * xsi1[gp]-1;
    (derivs[gp])(1,0) = 0;
    (derivs[gp])(2,0) = 0;
    (derivs[gp])(3,0) = 0;

    (derivs[gp])(4,0) = 4 * xsi2[gp];
    (derivs[gp])(5,0) = 0;
    (derivs[gp])(6,0) = 4 * xsi3[gp];
    (derivs[gp])(7,0) = 4 * xsi4[gp];
    (derivs[gp])(8,0) = 0;
    (derivs[gp])(9,0) = 0;

    (derivs[gp])(0,1) = 0;
    (derivs[gp])(1,1) = 4 * xsi2[gp] - 1;
    (derivs[gp])(2,1) = 0;
    (derivs[gp])(3,1) = 0;

    (derivs[gp])(4,1) = 4 * xsi1[gp];
    (derivs[gp])(5,1) = 4 * xsi3[gp];
    (derivs[gp])(6,1) = 0;
    (derivs[gp])(7,1) = 0;
    (derivs[gp])(8,1) = 4 * xsi4[gp];
    (derivs[gp])(9,1) = 0;

    (derivs[gp])(0,2) = 0;
    (derivs[gp])(1,2) = 0;
    (derivs[gp])(2,2) = 4 * xsi3[gp] - 1;
    (derivs[gp])(3,2) = 0;

    (derivs[gp])(4,2) = 0;
    (derivs[gp])(5,2) = 4 * xsi2[gp];
    (derivs[gp])(6,2) = 4 * xsi1[gp];
    (derivs[gp])(7,2) = 0;
    (derivs[gp])(8,2) = 0;
    (derivs[gp])(9,2) = 4 * xsi4[gp];

    (derivs[gp])(0,3) = 0;
    (derivs[gp])(1,3) = 0;
    (derivs[gp])(2,3) = 0;
    (derivs[gp])(3,3) = 4 * xsi4[gp] - 1;

    (derivs[gp])(4,3) = 0;
    (derivs[gp])(5,3) = 0;
    (derivs[gp])(6,3) = 0;
    (derivs[gp])(7,3) = 4 * xsi1[gp];
    (derivs[gp])(8,3) = 4 * xsi2[gp];
    (derivs[gp])(9,3) = 4 * xsi3[gp];
  }
  derivs_done = true;

  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Weights at 4 Gauss Points                            |
 *----------------------------------------------------------------------*/
const vector<double>& DRT::ELEMENTS::So_tet10::so_tet10_4gp_weights()
{
  static vector<double> weights(NUMGPT_SOTET10);
  static bool weights_done = false;
  if (weights_done) return weights;
  for (int gp=0; gp<NUMGPT_SOTET10; gp++)
    weights[gp] = 0.25;
  weights_done = true;

  return weights;
}

// This is copied from the Integrator_tet10_14point, but it uses only
// 10 gauss points for calculation, so I changed the name.
/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fcts at 10 Gauss Points                        |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMNOD_SOTET10,1> >& DRT::ELEMENTS::So_tet10::so_tet10_10gp_shapefcts()
{
  const int num_gp = 10;
  static vector<LINALG::Matrix<NUMNOD_SOTET10,1> > shapefcts(num_gp);
  static bool shapefcts_done = false;
  if (shapefcts_done) return shapefcts;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_10point;
  const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule);
  for (int gp=0; gp<num_gp; gp++) {
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    const double t = intpoints.qxg[gp][2];

    DRT::UTILS::shape_function_3D(shapefcts[gp], r, s, t, DRT::Element::tet10);
  }
  shapefcts_done = true;

  return shapefcts;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fct derivs at 10 Gauss Points                  |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMDIM_SOTET10,NUMNOD_SOTET10> >& DRT::ELEMENTS::So_tet10::so_tet10_10gp_derivs()
{
  const int num_gp = 10;
  static vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> > derivs(num_gp);
  static bool derivs_done = false;
  if (derivs_done) return derivs;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_10point;
  const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule);
  for (int gp=0; gp<num_gp; gp++) {
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    const double t = intpoints.qxg[gp][2];

    DRT::UTILS::shape_function_3D_deriv1(derivs[gp], r, s, t, DRT::Element::tet10);
  }
  derivs_done = true;

  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Weights at 10 Gauss Points                           |
 *----------------------------------------------------------------------*/
const vector<double>& DRT::ELEMENTS::So_tet10::so_tet10_10gp_weights()
{
  const int num_gp = 10;
  static vector<double> weights(num_gp);
  static bool weights_done = false;
  if (weights_done) return weights;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_10point;
  const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule);
  for (int gp=0; gp<num_gp; gp++)
    weights[gp] = intpoints.qwgt[gp];
  weights_done = true;

  return weights;
}

/*----------------------------------------------------------------------*
 |  lump mass matrix                                         bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_lumpmass(LINALG::Matrix<NUMDOF_SOTET10,NUMDOF_SOTET10>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned r=0; r<(*emass).M(); ++r)  // parse rows
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
