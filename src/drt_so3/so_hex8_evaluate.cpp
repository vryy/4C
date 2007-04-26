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
  Epetra_SerialDenseMatrix* deriv;    //[NUMNOD_SOH8*NUMDIM][NUMGPT_SOH8]
/* pointer to (static) weight factors at each gp */  
  Epetra_SerialDenseVector* weights;  //[NUMGPT_SOH8]
 
  soh8_shapederiv(&shapefct,&deriv,&weights);

  // update geometry
  Epetra_SerialDenseMatrix ex(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  Epetra_SerialDenseMatrix edis(NUMNOD_SOH8,NUMDIM_SOH8);// curr. element displacements
  
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    ex(i,0) = Nodes()[i]->X()[0];
    ex(i,1) = Nodes()[i]->X()[1];
    ex(i,2) = Nodes()[i]->X()[2];
    
    edis(i,0) = ex(i,0) + disp[i*NUMDOF_SOH8+0];
    edis(i,1) = ex(i,1) + disp[i*NUMDOF_SOH8+1];
    edis(i,2) = ex(i,2) + disp[i*NUMDOF_SOH8+2];
  }
 
 
  (*shapefct)(8,1);
  
  
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
  
  //static matrix objects with dimensions such that access of gp(=column) should be fast 
  static Epetra_SerialDenseMatrix  f(NUMNOD_SOH8,NUMGPT_SOH8);
  static Epetra_SerialDenseMatrix df(NUMDOF_SOH8,NUMGPT_SOH8);
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
        df(NUMDIM_SOH8*i+2,7) = -(1.0-r[i])*(1.0+s[i])*0.125;
    }
    *shapefct = &f;             // return adress of static object to target of pointer
    *deriv = &df;               // return adress of static object to target of pointer
    *weights  = &weightfactors; // return adress of static object to target of pointer
    fdf_eval = 1;               // now all arrays are filled statically
  }
  return;
}  // of soh8_shapederiv
        

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
