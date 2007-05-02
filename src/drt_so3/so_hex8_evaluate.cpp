/*!----------------------------------------------------------------------
\file so_hex8_evaluate.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_hex8.H"
#include "../discret/drt_discret.H"
#include "../discret/drt_utils.H"
#include "../discret/drt_exporter.H"
#include "../discret/drt_dserror.H"
#include "../discret/linalg_utils.H"

extern "C" 
{
#include "../headers/standardtypes.h"
// see if we can avoid this #include "../shell8/shell8.h"
}
#include "../discret/dstrc.H"
using namespace std; // cout etc.

/*----------------------------------------------------------------------*
 |                                                         maf 04/07    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::So_hex8::Evaluate(ParameterList& params, 
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DSTraceHelper dst("So_hex8::Evaluate");  
  DRT::Elements::So_hex8::ActionType act = So_hex8::none;
  
  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = So_hex8::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = So_hex8::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = So_hex8::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = So_hex8::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = So_hex8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = So_hex8::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = So_hex8::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = So_hex8::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = So_hex8::calc_struct_update_istep;
  else dserror("Unknown type of action for So_hex8");
  
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
      soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
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
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);
      soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
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
    case calc_struct_nlnstiffmass:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);
      soh8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
    }
    break;
    // evaluate stresses
    case calc_struct_stress:
    {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      soh8_stress(actmat,mydisp);
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
      ;// there is nothing to do here at the moment
    }
    break;
    default:
      dserror("Unknown type of action for Solid3");
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  Do stress calculation (private)                            maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_hex8::soh8_stress(struct _MATERIAL* material, 
                                    vector<double>& mydisp)
{
    dserror("Stress evaluation not yet ready");
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)    maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::So_hex8::EvaluateNeumann(ParameterList& params, 
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_hex8::soh8_nlnstiffmass(vector<int>&              lm, 
                                               vector<double>&           disp, 
                                               vector<double>&           residual,
                                               Epetra_SerialDenseMatrix* stiffmatrix,
                                               Epetra_SerialDenseMatrix* massmatrix,
                                               Epetra_SerialDenseVector* force,
                                               MATERIAL*                 material)
{
  DSTraceHelper dst("So_hex8::soh8_nlnstiffmass");  

/* ======================================================================*
 * SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
 * ======================================================================*/
/* pointer to (static) shape function array 
 * for each node, evaluated at each gp*/
  Epetra_SerialDenseMatrix* shapefct; //[NUMNOD_SOH8][NUMGPT_SOH8]
/* pointer to (static) shape function derivatives array 
 * for each node wrt to each direction, evaluated at each gp*/
  Epetra_SerialDenseMatrix* deriv;    //[NUMGPT_SOH8*NUMDIM][NUMNOD_SOH8]
/* pointer to (static) weight factors at each gp */  
  Epetra_SerialDenseVector* weights;  //[NUMGPT_SOH8]
/* ======================================================================*/
  Epetra_SerialDenseVector internalforce(NUMDOF_SOH8);
  Epetra_SerialDenseMatrix keu(NUMDOF_SOH8,NUMDOF_SOH8);
  Epetra_SerialDenseMatrix kgeo(NUMDOF_SOH8,NUMDOF_SOH8);
 
  soh8_shapederiv(&shapefct,&deriv,&weights);
  
  // update geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  Epetra_SerialDenseMatrix xcurr(NUMNOD_SOH8,NUMDIM_SOH8);  // curr. element displacements
  
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
    
    xcurr(i,0) = xrefe(i,0) + disp[i*NUMDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NUMDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NUMDOF_SOH8+2];
  }
  /*// testing ************
  double delta=0.1905;
  xcurr(1,0) += delta;
  xcurr(2,0) += delta;
  xcurr(5,0) += delta;
  xcurr(6,0) += delta;
  // testing ************/
  
  /* compute the Jacobian matrix which looks like:
   *         [ x_,r  y_,r  z_,r ]
   *     J = [ x_,s  y_,s  z_,s ]
   *         [ x_,t  y_,t  z_,t ]
   * for all GP: J(GP)=dN,i(GP) * X
   * therefore every 3 rows belong to one GP -> 3*8 rows */ 
  Epetra_SerialDenseMatrix jac(NUMDOF_SOH8,NUMDIM_SOH8);
  jac.Multiply('N','N',1.0,*deriv,xrefe,1.0);
  
  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    // compute determinant of J(GP) by Sarrus' rule
    double detJ=0.0;
    detJ= jac(3*gp+0,0) * jac(3*gp+1,1) * jac(3*gp+2,2)
        + jac(3*gp+0,1) * jac(3*gp+1,2) * jac(3*gp+2,0)
        + jac(3*gp+0,2) * jac(3*gp+1,0) * jac(3*gp+2,1)
        - jac(3*gp+0,0) * jac(3*gp+1,2) * jac(3*gp+2,1)
        - jac(3*gp+0,1) * jac(3*gp+1,0) * jac(3*gp+2,2)
        - jac(3*gp+0,2) * jac(3*gp+1,1) * jac(3*gp+2,0);
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");
    
    // compute inverse of J(GP)[3][3]
    Epetra_SerialDenseMatrix invjac(3,3);
    invjac(0,0) = jac(3*gp+1,1)*jac(3*gp+2,2) - jac(3*gp+2,1)*jac(3*gp+1,2);
    invjac(0,1) = jac(3*gp+0,2)*jac(3*gp+2,1) - jac(3*gp+0,1)*jac(3*gp+2,2);
    invjac(0,2) = jac(3*gp+0,1)*jac(3*gp+1,2) - jac(3*gp+0,2)*jac(3*gp+1,1);
    invjac(1,0) = jac(3*gp+1,2)*jac(3*gp+2,0) - jac(3*gp+2,2)*jac(3*gp+1,0);
    invjac(1,1) = jac(3*gp+0,0)*jac(3*gp+2,2) - jac(3*gp+0,2)*jac(3*gp+2,0);
    invjac(1,2) = jac(3*gp+0,2)*jac(3*gp+1,0) - jac(3*gp+0,0)*jac(3*gp+1,2);
    invjac(2,0) = jac(3*gp+1,0)*jac(3*gp+2,1) - jac(3*gp+2,0)*jac(3*gp+1,1);
    invjac(2,1) = jac(3*gp+0,1)*jac(3*gp+2,0) - jac(3*gp+0,0)*jac(3*gp+2,1);
    invjac(2,2) = jac(3*gp+0,0)*jac(3*gp+1,1) - jac(3*gp+0,1)*jac(3*gp+1,0);
    // Scalarmultiply with 1/detJ later (at N_XYZ) !
    
    // build submatrix of deriv at actual gp
    Epetra_SerialDenseMatrix deriv_gp(NUMDIM_SOH8,NUMGPT_SOH8);
    for (int m=0; m<NUMDIM_SOH8; ++m)
    {
        for (int n=0; n<NUMGPT_SOH8; ++n)
        {
            deriv_gp(m,n)=(*deriv)(3*gp+m,n);
        }
    }
    // Derivatives at gp w.r.t. material coordinates
    Epetra_SerialDenseMatrix N_XYZ(NUMDIM_SOH8,NUMNOD_SOH8);
    N_XYZ.Multiply('N','N',1/detJ,invjac,deriv_gp,1.0);
    // (material) deformation gradient F = d xxurr / d xrefe
    Epetra_SerialDenseMatrix defgrd(NUMDIM_SOH8,NUMDIM_SOH8);
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
        defgrd(0,0) += xcurr(i,0) * N_XYZ(0,i);     // dx/dX
        defgrd(1,1) += xcurr(i,1) * N_XYZ(1,i);     // dy/dY
        defgrd(2,2) += xcurr(i,2) * N_XYZ(2,i);     // dz/dZ
        defgrd(0,1) += xcurr(i,0) * N_XYZ(1,i);     // dx/dY
        defgrd(0,2) += xcurr(i,0) * N_XYZ(2,i);     // dx/dZ
        defgrd(1,0) += xcurr(i,1) * N_XYZ(0,i);     // dy/dX
        defgrd(1,2) += xcurr(i,1) * N_XYZ(2,i);     // dy/dZ
        defgrd(2,0) += xcurr(i,2) * N_XYZ(0,i);     // dz/dX
        defgrd(2,1) += xcurr(i,2) * N_XYZ(1,i);     // dz/dY
    }
    // Cauchy-Green Tensor = F^T * F
    Epetra_SerialDenseMatrix cauchygreen(NUMDIM_SOH8,NUMDIM_SOH8);
    cauchygreen.Multiply('T','N',1.0,defgrd,defgrd,1.0);
    // Green-Lagrange strains matrix E = 0.5 * (cauchygreen - Identity)
    // GL strains Vector strain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector strain(NUMSTR_SOH8);
    strain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    strain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    strain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    strain(3) = cauchygreen(0,1);
    strain(4) = cauchygreen(1,2);
    strain(5) = cauchygreen(2,0);
    /* non-linear B-operator (may so be called, meaning
     * of B-operator is not so sharp in the non-linear realm) *
     * B = F . Bl *
     *
     *      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
     *      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
     *      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
     * B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
     *      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
     *      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
     *      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
     *      [                                                         ]
     *      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
     *      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
     *      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
     *      [                                                         ]
     *      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
     *      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
     *      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
     */
    Epetra_SerialDenseMatrix bop(NUMSTR_SOH8,NUMDOF_SOH8);
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
    // call material law
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOH8,NUMSTR_SOH8);
    soh8_mat_sel(&cmat);
    
    // evaluate stresses
    Epetra_SerialDenseVector stress(NUMSTR_SOH8);
    cmat.Multiply('N',strain,stress);               // sigma = C * epsilon
    // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
    internalforce.Multiply('T','N',detJ * (*weights)(gp),bop,stress,1.0); // local intforce
    (*force).Multiply('T','N',detJ * (*weights)(gp),bop,stress,1.0);
    
    // integrate `elastic' and `initial-displacement' stiffness matrix
    // keu = keu + (B^T . C . B) * detJ * w(gp)
    Epetra_SerialDenseMatrix cb(NUMSTR_SOH8,NUMDOF_SOH8);
    cb.Multiply('N','N',1.0,cmat,bop,1.0);          // C . B
    keu.Multiply('T','N',detJ * (*weights)(gp),bop,cb,1.0);  // local keu
    (*stiffmatrix).Multiply('T','N',detJ * (*weights)(gp),bop,cb,1.0);
    
    // intergrate `geometric' stiffness matrix and add to keu
    // kgeo = kgeo + (BLin^T . sigma . BLin) * detJ * w(gp)  with BLin = Ni,Xj see NiliFEM-Skript
    Epetra_SerialDenseVector sfac(stress);     // auxiliary integrated stress
    sfac.Scale(detJ * (*weights)(gp));         // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
    vector<double> SmBlin(NUMDIM_SOH8);        // intermediate Sm.Blin
    for (int inod=0; inod<NUMNOD_SOH8; ++inod)
    {
      SmBlin[0] = sfac(0) * N_XYZ(0,inod) + sfac(3) * N_XYZ(1,inod) + sfac(5) * N_XYZ(2,inod); 
      SmBlin[1] = sfac(3) * N_XYZ(0,inod) + sfac(2) * N_XYZ(1,inod) + sfac(4) * N_XYZ(2,inod); 
      SmBlin[2] = sfac(5) * N_XYZ(0,inod) + sfac(4) * N_XYZ(1,inod) + sfac(3) * N_XYZ(2,inod);
      for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod)
      {
        double bopstrbop = 0.0;
        for (int idim=0; idim<NUMDIM_SOH8; ++idim) bopstrbop += N_XYZ(idim,jnod) * SmBlin[idim];
        kgeo(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += bopstrbop;   // local kgeo
        kgeo(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += bopstrbop;   // local kgeo
        kgeo(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += bopstrbop;   // local kgeo
        (*stiffmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += bopstrbop;
        (*stiffmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += bopstrbop;
        (*stiffmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += bopstrbop;
      }
    } // end of intergrate `geometric' stiffness matrix
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/
  
  return;
} // DRT::Elements::Shell8::s8_nlnstiffmass



/*----------------------------------------------------------------------*
 |  shape functions and derivatives for So_hex8                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::So_hex8::soh8_shapederiv(
                             Epetra_SerialDenseMatrix** shapefct, //pointer to pointer
                             Epetra_SerialDenseMatrix** deriv,    //pointer to pointer
                             Epetra_SerialDenseVector** weights)  //pointer to pointer
{
  DSTraceHelper dst("So_hex8::soh8_shapederiv");
  
  //static matrix objects 
  static Epetra_SerialDenseMatrix  f(NUMNOD_SOH8,NUMGPT_SOH8);
  static Epetra_SerialDenseMatrix df(NUMDOF_SOH8,NUMNOD_SOH8);
  static Epetra_SerialDenseVector weightfactors(NUMGPT_SOH8);
  static int fdf_eval;                        // flag for re-evaluation
  
  const double gploc    = 1.0/sqrt(3.0);      // gp sampling point value for linear fct
  const double gpw      = 1.0;                // weight at every gp for linear fct
  
  if (fdf_eval!=0)              // if true f,df already evaluated
  {
    *shapefct = &f;             // return adress of static object to target of pointer
    *deriv    = &df;            // return adress of static object to target of pointer
    *weights  = &weightfactors; // return adress of static object to target of pointer
    return;
  }
  else 
  {
    // (r,s,t) gp-locations of fully integrated linear 8-node Hex
    const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
    const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};
    const double t[NUMGPT_SOH8] = {-gploc,-gploc,-gploc,-gploc, gploc, gploc, gploc, gploc};
    const double w[NUMGPT_SOH8] = {   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw,   gpw};
    
    for (int i=0; i<NUMGPT_SOH8; ++i)  //fill up nodal f at each gp
    {
        f(0,i) = (1.0-r[i])*(1.0-s[i])*(1.0-t[i])*0.125;  
        f(1,i) = (1.0+r[i])*(1.0-s[i])*(1.0-t[i])*0.125;
        f(2,i) = (1.0+r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
        f(3,i) = (1.0-r[i])*(1.0+s[i])*(1.0-t[i])*0.125;
        f(4,i) = (1.0-r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
        f(5,i) = (1.0+r[i])*(1.0-s[i])*(1.0+t[i])*0.125;
        f(6,i) = (1.0+r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
        f(7,i) = (1.0-r[i])*(1.0+s[i])*(1.0+t[i])*0.125;
        
        weightfactors[i] = w[i]*w[i]*w[i]; // just for clarity how to get weight factors 
    }
    // fill up df wrt to 3 directions (NUMDIM) at each gp 
    for (int i=0; i<NUMGPT_SOH8; ++i)
    {
        // df wrt to r(+0) for each node(0..7) at each gp [i]
        df(NUMDIM_SOH8*i+0,0) = -(1.0-s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,1) =  (1.0-s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,2) =  (1.0+s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,3) = -(1.0+s[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+0,4) = -(1.0-s[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+0,5) =  (1.0-s[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+0,6) =  (1.0+s[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+0,7) = -(1.0+s[i])*(1.0+t[i])*0.125;
        
        // df wrt to s(+1) for each node(0..7) at each gp [i]
        df(NUMDIM_SOH8*i+1,0) = -(1.0-r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,1) = -(1.0+r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,2) =  (1.0+r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,3) =  (1.0-r[i])*(1.0-t[i])*0.125;
        df(NUMDIM_SOH8*i+1,4) = -(1.0-r[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+1,5) = -(1.0+r[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+1,6) =  (1.0+r[i])*(1.0+t[i])*0.125;
        df(NUMDIM_SOH8*i+1,7) =  (1.0-r[i])*(1.0+t[i])*0.125;
        
        // df wrt to t(+2) for each node(0..7) at each gp [i]
        df(NUMDIM_SOH8*i+2,0) = -(1.0-r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,1) = -(1.0+r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,2) = -(1.0+r[i])*(1.0+s[i])*0.125;
        df(NUMDIM_SOH8*i+2,3) = -(1.0-r[i])*(1.0+s[i])*0.125;
        df(NUMDIM_SOH8*i+2,4) =  (1.0-r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,5) =  (1.0+r[i])*(1.0-s[i])*0.125;
        df(NUMDIM_SOH8*i+2,6) =  (1.0+r[i])*(1.0+s[i])*0.125;
        df(NUMDIM_SOH8*i+2,7) =  (1.0-r[i])*(1.0+s[i])*0.125;
    }
    *shapefct = &f;             // return adress of static object to target of pointer
    *deriv = &df;               // return adress of static object to target of pointer
    *weights  = &weightfactors; // return adress of static object to target of pointer
    fdf_eval = 1;               // now all arrays are filled statically
  }
  return;
}  // of soh8_shapederiv

void DRT::Elements::So_hex8::soh8_mat_sel(Epetra_SerialDenseMatrix* cmat)
{
  DSTraceHelper dst("So_hex8::soh8_mat_sel");
  /* Young's modulus (modulus of elasticity */
  double Emod = mat->m.stvenant->youngs;
  /* Poisson's ratio */
  double nu = mat->m.stvenant->possionratio;
  /*--------------------------------------------------------------------*/
  /* isotropic elasticity tensor C in matrix notion */
  /*                       [ 1-nu     nu     nu |          0    0    0 ]
   *                       [        1-nu     nu |          0    0    0 ]
   *           E           [               1-nu |          0    0    0 ]
   *   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
   *       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
   *                       [                    |      (1-2*nu)/2    0 ]
   *                       [ symmetric          |           (1-2*nu)/2 ]
   */
  double mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */
  /* write non-zero components */
  (*cmat)(0,0) = mfac*(1.0-nu);
  (*cmat)(0,1) = mfac*nu;
  (*cmat)(0,2) = mfac*nu;
  (*cmat)(1,0) = mfac*nu;
  (*cmat)(1,1) = mfac*(1.0-nu);
  (*cmat)(1,2) = mfac*nu;
  (*cmat)(2,0) = mfac*nu;
  (*cmat)(2,1) = mfac*nu;
  (*cmat)(2,2) = mfac*(1.0-nu);
  /* ~~~ */
  (*cmat)(3,3) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(4,4) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(5,5) = mfac*0.5*(1.0-2.0*nu);
  
  return;
}  // of soh8_mat_sel
                                            
                             

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
