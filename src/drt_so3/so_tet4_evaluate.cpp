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

#include "so_tet4.H"
#include "so_integrator.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"


using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |                                                         maf 04/07    |
 | vector of material laws                                              |
 | defined in global_control.c											|
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat; ///< C-style material struct


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
  else if (action=="calc_struct_nlnstifflmass") act = So_tet4::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")        act = So_tet4::calc_struct_stress;
  else if (action=="postprocess_stress")        act = So_tet4::postprocess_stress;
  else if (action=="calc_struct_eleload")       act = So_tet4::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = So_tet4::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = So_tet4::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = So_tet4::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")  act = So_tet4::calc_struct_reset_istep;
#ifdef PRESTRESS
  else if (action=="calc_struct_prestress_update") act = So_tet4::prestress_update;
#endif
  else dserror("Unknown type of action for So_tet4");

  // get the material law
  MATERIAL* actmat = &(mat[material_-1]);

  // what should the element do
  switch(act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      so_tet4_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat);

    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      so_tet4_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat);
    }
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      so_tet4_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat);

      if (act==calc_struct_nlnstifflmass) so_tet4_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
      RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      if (stressdata==null) dserror("Cannot get stress 'data'");
      if (straindata==null) dserror("Cannot get strain 'data'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      Epetra_SerialDenseMatrix stress(NUMGPT_SOTET4,NUMSTR_SOTET4);
      Epetra_SerialDenseMatrix strain(NUMGPT_SOTET4,NUMSTR_SOTET4);
      bool cauchy = params.get<bool>("cauchy", false);
      string iostrain = params.get<string>("iostrain", "none");
      bool ea = (iostrain == "euler_almansi");
      so_tet4_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,actmat,cauchy,ea);
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
      RCP<Epetra_SerialDenseMatrix> gpstress = (*gpstressmap)[gid];

      if (stresstype=="ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_SOTET4,NUMSTR_SOTET4);
        so_tet4_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOTET4);

        for (int i=0;i<NUMNOD_SOTET4;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_SOTET4;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_SOTET4;++i){
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
          for (int i = 0; i < NUMSTR_SOTET4; ++i) {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_SOTET4; ++j) {
              (*((*elestress)(i)))[lid] += 1.0/NUMGPT_SOTET4 * (*gpstress)(j,i);
            }
          }
        }
      }
      else if (stresstype=="cxyz_ndxyz") {
        // extrapolate stresses/strains at Gauss points to nodes
        Epetra_SerialDenseMatrix nodalstresses(NUMNOD_SOTET4,NUMSTR_SOTET4);
        so_tet4_expol(*gpstress,nodalstresses);

        // average nodal stresses/strains between elements
        // -> divide by number of adjacent elements
        vector<int> numadjele(NUMNOD_SOTET4);

        for (int i=0;i<NUMNOD_SOTET4;++i){
          DRT::Node* node=Nodes()[i];
          numadjele[i]=node->NumElement();
        }

        for (int i=0;i<NUMNOD_SOTET4;++i){
          elevec1(3*i)=nodalstresses(i,0)/numadjele[i];
          elevec1(3*i+1)=nodalstresses(i,1)/numadjele[i];
          elevec1(3*i+2)=nodalstresses(i,2)/numadjele[i];
        }
        for (int i=0;i<NUMNOD_SOTET4;++i){
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
          for (int i = 0; i < NUMSTR_SOTET4; ++i) {
            (*((*elestress)(i)))[lid] = 0.;
            for (int j = 0; j < NUMGPT_SOTET4; ++j) {
              (*((*elestress)(i)))[lid] += 1.0/NUMGPT_SOTET4 * (*gpstress)(j,i);
            }
          }
        }
      }
      else {
        dserror("unknown type of stress/strain output on element level");
      }
    }
    break;

#ifdef PRESTRESS
    case prestress_update:
    {
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get displacement state");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // build incremental def gradient for every gauss point
      LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOTET4,9);
      DefGradient(mydisp,gpdefgrd,*prestress_);

      // update deformation gradient and put back to storage
      LINALG::SerialDenseMatrix deltaF(3,3);
      LINALG::SerialDenseMatrix Fhist(3,3);
      LINALG::SerialDenseMatrix Fnew(3,3);
      for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
      {
        prestress_->StoragetoMatrix(gp,deltaF,gpdefgrd);
        prestress_->StoragetoMatrix(gp,Fhist,prestress_->FHistory());
        Fnew.Multiply('N','N',1.0,deltaF,Fhist,0.0);
        prestress_->MatrixtoStorage(gp,Fnew,prestress_->FHistory());
      }

      // push-forward invJ for every gaussian point
      UpdateJacobianMapping(mydisp,*prestress_);
    }
    break;
#endif

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    case calc_struct_update_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;

    case calc_struct_update_imrlike:
    {
      ;// there is nothing to do here at the moment
    }
    break;

    case calc_struct_reset_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
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

      so_tet4_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,actmat);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
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
 |  init the element jacobian mapping (protected)              gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::InitJacobianMapping()
{
  LINALG::SerialDenseMatrix xrefe(NUMNOD_SOTET4,NUMDIM_SOTET4);
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  /* get the matrix of the coordinates of nodes needed to compute the volume,
  ** which is used here as detJ in the quadrature rule.
  ** ("Jacobian matrix") for the quadrarture rule:
  **             [  1    1    1    1  ]
  **         J = [ X_1  X_2  X_3  X_4 ]
  **             [ Y_1  Y_2  Y_3  Y_4 ]
  **		 [ Z_1  Z_2  Z_3  Z_4 ]
  */
  LINALG::SerialDenseMatrix J(NUMCOORD_SOTET4,NUMCOORD_SOTET4);
  for (int i=0; i<4; i++)  J(0,i)=1;
  for (int row=0;row<3;row++)
    for (int col=0;col<4;col++)
      J(row+1,col)= xrefe(col,row);
  // volume of the element
  V_ = LINALG::DeterminantLU(J)/6.0;

  nxyz_.resize(NUMGPT_SOTET4);
  const static DRT::ELEMENTS::Integrator_tet4_1point tet4_dis;
  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    LINALG::SerialDenseMatrix jac(NUMCOORD_SOTET4,NUMCOORD_SOTET4);
    {
      LINALG::SerialDenseMatrix tmp(NUMCOORD_SOTET4-1,NUMCOORD_SOTET4);
      tmp.Multiply('T','N',1.0,xrefe,tet4_dis.deriv_gp[gp],0.0);
      for (int i=0; i<4; i++) jac(0,i)=1;
      for (int row=0;row<3;row++)
        for (int col=0;col<4;col++)
          jac(row+1,col)=tmp(row,col);
    }
    // size is 4x3
    Epetra_SerialDenseMatrix  I_aug(NUMCOORD_SOTET4,NUMDIM_SOTET4);
    // size is 4x3
    Epetra_SerialDenseMatrix partials(NUMCOORD_SOTET4,NUMDIM_SOTET4);
    I_aug(1,0)=1;
    I_aug(2,1)=1;
    I_aug(3,2)=1;

    Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    solve_for_inverseJac.SetVectors(partials,I_aug);// set X=partials, B=I_aug
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();         // partials = jac^-1.I_aug
    if ((err != 0) && (err2!=0))
    	dserror("Inversion of Jacobian failed");

    //nxyz_[gp] = N_xsi_k*partials
    // size is 4x3
    nxyz_[gp].Shape(NUMNOD_SOTET4,NUMDIM_SOTET4);
    nxyz_[gp].Multiply('N','N',1.0,tet4_dis.deriv_gp[gp],partials,0.0);
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

#ifdef PRESTRESS
    if (!(prestress_->IsInit()))
      prestress_->MatrixtoStorage(gp,nxyz_[gp],prestress_->JHistory());
#endif

  } // for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
#ifdef PRESTRESS
  prestress_->IsInit() = true;
#endif
  return;
}

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
      Epetra_SerialDenseMatrix* elestress,      // stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      // strains at GP
      struct _MATERIAL*         material,       // element material data
      const bool                cauchy,         // stress output options
      const bool                ea)             // strain output options
{
/* =============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_4  with 1 GAUSS POINTS*
** =============================================================================*/
   const static DRT::ELEMENTS::Integrator_tet4_1point tet4_dis;
/* ============================================================================*/
  double density;
  // element geometry
  /* structure of xrefe:
    **             [  X_1   Y_1   Z_1  ]
    **     xrefe = [  X_2   Y_2   Z_2  ]
    **             [   |     |     |   ]
    **             [  X_4   Y_4   Z_4  ]
    */
  /* structure of xcurr:
    **             [  x_1   y_1   z_1  ]
    **     xcurr = [  x_2   y_2   z_2  ]
    **             [   |     |     |   ]
    **             [  x_4   y_4   z_4  ]
    */
  // current  displacements of element
  LINALG::SerialDenseMatrix xdisp(NUMNOD_SOTET4,NUMDIM_SOTET4);

  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOTET4+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET4+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET4+2];
  }


  //volume of a tetrahedra
  double detJ = V_;

  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOTET4; gp++)
  {
    const Epetra_SerialDenseMatrix& nxyz = nxyz_[gp];

    //                                      d xcurr
    // (material) deformation gradient F = --------- = xcurr^T * nxyz^T
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

    // size is 3x3
    LINALG::SerialDenseMatrix defgrd(NUMDIM_SOTET4,NUMDIM_SOTET4);
#if defined(PRESTRESS) || defined(POSTSTRESS)
    {
      // get derivatives wrt to last spatial configuration
      LINALG::SerialDenseMatrix N_xyz(NUMNOD_SOTET4,NUMDIM_SOTET4);
      prestress_->StoragetoMatrix(gp,N_xyz,prestress_->JHistory());

      // build multiplicative incremental defgrd
      defgrd.Multiply('T','N',1.0,xdisp,N_xyz,0.0);
      defgrd(0,0) += 1.0;
      defgrd(1,1) += 1.0;
      defgrd(2,2) += 1.0;

      // get stored old incremental F
      LINALG::SerialDenseMatrix Fhist(3,3);
      prestress_->StoragetoMatrix(gp,Fhist,prestress_->FHistory());

      // build total defgrd = delta F * F_old
      LINALG::SerialDenseMatrix Fnew(3,3);
      Fnew.Multiply('N','N',1.0,defgrd,Fhist,0.0);
      defgrd = Fnew;
    }
#else
    defgrd.Multiply('T','N',1.0,xdisp,nxyz,0.0);
    defgrd(0,0)+=1;
    defgrd(1,1)+=1;
    defgrd(2,2)+=1;
#endif

    // Right Cauchy-Green tensor = F^T * F
    // size is 3x3
    LINALG::SerialDenseMatrix cauchygreen(NUMDIM_SOTET4,NUMDIM_SOTET4);
    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,0.0);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain(NUMSTR_SOTET4);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    // return gp strains (only in case of stress/strain output)
    if (elestrain != NULL)
    {
      if (!ea) // output Green-Lagrange strains
      {
        for (int i = 0; i < 3; ++i)
          (*elestrain)(gp,i) = glstrain(i);
        for (int i = 3; i < 6; ++i)
          (*elestrain)(gp,i) = 0.5 * glstrain(i);
      }
      else
      {
        // rewriting Green-Lagrange strains in matrix format
        LINALG::SerialDenseMatrix gl(NUMDIM_SOTET4,NUMDIM_SOTET4);
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
        Epetra_SerialDenseMatrix invdefgrd(defgrd); // make a copy here otherwise defgrd is destroyed!
        LINALG::NonsymInverse3x3(invdefgrd);

        LINALG::SerialDenseMatrix temp(NUMDIM_SOTET4,NUMDIM_SOTET4);
        LINALG::SerialDenseMatrix euler_almansi(NUMDIM_SOTET4,NUMDIM_SOTET4);
        temp.Multiply('N','N',1.0,gl,invdefgrd,0.0);
        euler_almansi.Multiply('T','N',1.0,invdefgrd,temp,0.0);

        (*elestrain)(gp,0) = euler_almansi(0,0);
        (*elestrain)(gp,1) = euler_almansi(1,1);
        (*elestrain)(gp,2) = euler_almansi(2,2);
        (*elestrain)(gp,3) = euler_almansi(0,1);
        (*elestrain)(gp,4) = euler_almansi(1,2);
        (*elestrain)(gp,5) = euler_almansi(0,2);
      }
    }

    /*----------------------------------------------------------------------*
      the B-operator used is equivalent to the one used in hex8, this needs
      to be checked if it is ok, but from the mathematics point of view, the only
      thing that needed to be changed is the NUMDOF
      ----------------------------------------------------------------------*/
    /*
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
    // size is 6x12
    LINALG::SerialDenseMatrix bop(NUMSTR_SOTET4,NUMDOF_SOTET4);
    for (int i=0; i<NUMNOD_SOTET4; i++)
    {
      bop(0,NODDOF_SOTET4*i+0) = defgrd(0,0)*nxyz(i,0);
      bop(0,NODDOF_SOTET4*i+1) = defgrd(1,0)*nxyz(i,0);
      bop(0,NODDOF_SOTET4*i+2) = defgrd(2,0)*nxyz(i,0);
      bop(1,NODDOF_SOTET4*i+0) = defgrd(0,1)*nxyz(i,1);
      bop(1,NODDOF_SOTET4*i+1) = defgrd(1,1)*nxyz(i,1);
      bop(1,NODDOF_SOTET4*i+2) = defgrd(2,1)*nxyz(i,1);
      bop(2,NODDOF_SOTET4*i+0) = defgrd(0,2)*nxyz(i,2);
      bop(2,NODDOF_SOTET4*i+1) = defgrd(1,2)*nxyz(i,2);
      bop(2,NODDOF_SOTET4*i+2) = defgrd(2,2)*nxyz(i,2);
      /* ~~~ */
      bop(3,NODDOF_SOTET4*i+0) = defgrd(0,0)*nxyz(i,1) + defgrd(0,1)*nxyz(i,0);
      bop(3,NODDOF_SOTET4*i+1) = defgrd(1,0)*nxyz(i,1) + defgrd(1,1)*nxyz(i,0);
      bop(3,NODDOF_SOTET4*i+2) = defgrd(2,0)*nxyz(i,1) + defgrd(2,1)*nxyz(i,0);
      bop(4,NODDOF_SOTET4*i+0) = defgrd(0,1)*nxyz(i,2) + defgrd(0,2)*nxyz(i,1);
      bop(4,NODDOF_SOTET4*i+1) = defgrd(1,1)*nxyz(i,2) + defgrd(1,2)*nxyz(i,1);
      bop(4,NODDOF_SOTET4*i+2) = defgrd(2,1)*nxyz(i,2) + defgrd(2,2)*nxyz(i,1);
      bop(5,NODDOF_SOTET4*i+0) = defgrd(0,2)*nxyz(i,0) + defgrd(0,0)*nxyz(i,2);
      bop(5,NODDOF_SOTET4*i+1) = defgrd(1,2)*nxyz(i,0) + defgrd(1,0)*nxyz(i,2);
      bop(5,NODDOF_SOTET4*i+2) = defgrd(2,2)*nxyz(i,0) + defgrd(2,0)*nxyz(i,2);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    // size is 6x6
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOTET4,NUMSTR_SOTET4);
    // size is 6
    Epetra_SerialDenseVector stress(NUMSTR_SOTET4);
    so_tet4_mat_sel(&stress,&cmat,&density,&glstrain, &defgrd, gp);

    // return gp stresses
    if (elestress != NULL)
    { // return 2nd Piola-Kirchhoff stresses
      if (!cauchy)
      {
        for (int i = 0; i < NUMSTR_SOTET4; ++i)
          (*elestress)(gp,i) = stress(i);
      }
      else
      {
        // return Cauchy stresses
        double detF = defgrd(0,0)*defgrd(1,1)*defgrd(2,2) +
                      defgrd(0,1)*defgrd(1,2)*defgrd(2,0) +
                      defgrd(0,2)*defgrd(1,0)*defgrd(2,1) -
                      defgrd(0,2)*defgrd(1,1)*defgrd(2,0) -
                      defgrd(0,0)*defgrd(1,2)*defgrd(2,1) -
                      defgrd(0,1)*defgrd(1,0)*defgrd(2,2);

        LINALG::SerialDenseMatrix pkstress(NUMDIM_SOTET4,NUMDIM_SOTET4);
        pkstress(0,0) = stress(0);
        pkstress(0,1) = stress(3);
        pkstress(0,2) = stress(5);
        pkstress(1,0) = pkstress(0,1);
        pkstress(1,1) = stress(1);
        pkstress(1,2) = stress(4);
        pkstress(2,0) = pkstress(0,2);
        pkstress(2,1) = pkstress(1,2);
        pkstress(2,2) = stress(2);

        LINALG::SerialDenseMatrix temp(NUMDIM_SOTET4,NUMDIM_SOTET4);
        LINALG::SerialDenseMatrix cauchystress(NUMDIM_SOTET4,NUMDIM_SOTET4);
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

    if (force != NULL && stiffmatrix != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).Multiply('T','N',detJ * (tet4_dis.weights)(gp),bop,stress,1.0);

      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      // size is 6x12
      LINALG::SerialDenseMatrix cb(NUMSTR_SOTET4,NUMDOF_SOTET4);
      cb.Multiply('N','N',1.0,cmat,bop,0.0);          // temporary C . B
      // size is 12x12
      stiffmatrix->Multiply('T','N',detJ * (tet4_dis.weights)(gp),bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu
      // auxiliary integrated stress
      Epetra_SerialDenseVector sfac(stress);
      // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      sfac.Scale(detJ * (tet4_dis.weights)(gp));
      // intermediate Sm.B_L
      vector<double> SmB_L(NUMDIM_SOTET4);
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)
      // with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOTET4; ++inod)
      {
        SmB_L[0] = sfac(0) * nxyz(inod,0) + sfac(3) * nxyz(inod,1) + sfac(5) * nxyz(inod,2);
        SmB_L[1] = sfac(3) * nxyz(inod,0) + sfac(1) * nxyz(inod,1) + sfac(4) * nxyz(inod,2);
        SmB_L[2] = sfac(5) * nxyz(inod,0) + sfac(4) * nxyz(inod,1) + sfac(2) * nxyz(inod,2);
        for (int jnod=0; jnod<NUMNOD_SOTET4; ++jnod)
        {
          double bopstrbop = 0.0;            // intermediate value
          for (int idim=0; idim<NUMDIM_SOTET4; ++idim)
            bopstrbop += nxyz(jnod,idim)* SmB_L[idim];
          (*stiffmatrix)(NUMDIM_SOTET4*inod+0,NUMDIM_SOTET4*jnod+0) += bopstrbop;
          (*stiffmatrix)(NUMDIM_SOTET4*inod+1,NUMDIM_SOTET4*jnod+1) += bopstrbop;
          (*stiffmatrix)(NUMDIM_SOTET4*inod+2,NUMDIM_SOTET4*jnod+2) += bopstrbop;
        }
      }
    }
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/


  // static integrator created in any case to safe "if-case"
  const static DRT::ELEMENTS::Integrator_tet4_4point tet4_mass;
  // evaluate mass matrix
  if (massmatrix != NULL)
  {
    //consistent mass matrix evaluated using a 4-point rule
    for (int gp=0; gp<tet4_mass.num_gp; gp++)
    {
      for (int inod=0; inod<NUMNOD_SOTET4; ++inod)
      {
        for (int jnod=0; jnod<NUMNOD_SOTET4; ++jnod)
        {
          double massfactor = (tet4_mass.shapefct_gp[gp])(inod) * density *
                              (tet4_mass.shapefct_gp[gp])(jnod) * detJ *
                              (tet4_mass.weights)(gp);
          (*massmatrix)(NUMDIM_SOTET4*inod+0,NUMDIM_SOTET4*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4*inod+1,NUMDIM_SOTET4*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4*inod+2,NUMDIM_SOTET4*jnod+2) += massfactor;
        }
      }
    }
  }// end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

  return;
} // DRT::ELEMENTS::So_tet4::so_tet4_nlnstiffmass

/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                               bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_lumpmass(Epetra_SerialDenseMatrix* emass)
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

/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Sotet4Register::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_so_tet4) continue;
    DRT::ELEMENTS::So_tet4* actele = dynamic_cast<DRT::ELEMENTS::So_tet4*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}


#if defined(PRESTRESS) || defined(POSTSTRESS)
/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)   gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::DefGradient(const vector<double>& disp,
                                         Epetra_SerialDenseMatrix& gpdefgrd,
                                         DRT::ELEMENTS::PreStress& prestress)
{
  // update element geometry
  LINALG::SerialDenseMatrix xdisp(NUMNOD_SOTET4,NUMDIM_SOTET4);
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOTET4+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET4+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET4+2];
  }

  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    // get derivatives wrt to last spatial configuration
    LINALG::SerialDenseMatrix N_xyz(NUMNOD_SOTET4,NUMDIM_SOTET4);
    prestress_->StoragetoMatrix(gp,N_xyz,prestress_->JHistory());

    // build multiplicative incremental defgrd
    LINALG::SerialDenseMatrix defgrd(NUMDIM_SOTET4,NUMDIM_SOTET4);
    defgrd.Multiply('T','N',1.0,xdisp,N_xyz,0.0);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;

    prestress.MatrixtoStorage(gp,defgrd,gpdefgrd);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  compute Jac.mapping wrt deformed configuration (protected) gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::UpdateJacobianMapping(
                                            const vector<double>& disp,
                                            DRT::ELEMENTS::PreStress& prestress)
{
  // get incremental disp
  LINALG::SerialDenseMatrix xdisp(NUMNOD_SOTET4,NUMDIM_SOTET4);
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOTET4+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET4+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET4+2];
  }

  LINALG::SerialDenseMatrix nxyzhist(NUMNOD_SOTET4,NUMDIM_SOTET4);
  LINALG::SerialDenseMatrix nxyznew(NUMNOD_SOTET4,NUMDIM_SOTET4);
  LINALG::SerialDenseMatrix defgrd(NUMDIM_SOTET4,NUMDIM_SOTET4);

  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    // get the nxyz old state
    prestress.StoragetoMatrix(gp,nxyzhist,prestress.JHistory());
    // build multiplicative incremental defgrd
    defgrd.Multiply('T','N',1.0,xdisp,nxyzhist,0.0);
    defgrd(0,0) += 1.0;
    defgrd(1,1) += 1.0;
    defgrd(2,2) += 1.0;
    // make inverse of this defgrd
    LINALG::NonsymInverse3x3(defgrd);

    // push-forward of nxyz
    nxyznew.Multiply('N','N',1.0,nxyzhist,defgrd,0.0);
    // store new reference configuration
    prestress.MatrixtoStorage(gp,nxyznew,prestress.JHistory());

  } // for (int gp=0; gp<NUMGPT_WEG6; ++gp)

  return;
}
#endif // #if defined(PRESTRESS) || defined(POSTSTRESS)


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
