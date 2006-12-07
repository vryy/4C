/*!----------------------------------------------------------------------
\file shell8.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "shell8.H"
#include "drt_discret.H"
#include "drt_utils.H"
#include "drt_dserror.H"
extern "C" 
{
#include "../headers/standardtypes.h"
#include "../shell8/shell8.h"
}
#include "dstrc.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Elements::Shell8::Evaluate(ParameterList& params, 
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DSTraceHelper dst("Shell8::Evaluate");  
  DRT::Elements::Shell8::ActionType act = Shell8::none;
  
  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = Shell8::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Shell8::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Shell8::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = Shell8::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = Shell8::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = Shell8::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Shell8::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Shell8::calc_struct_fsiload;
  else dserror("Unknown type of action for Shell8");
  
  // get the material law
  MATERIAL* actmat = &(mat[material_-1]);
  
  // see whether elements where initialized
  if (!Init())
  {
    s8_initialize(discretization);
    SetInit(true);
  }
  
  switch(act)
  {
    case calc_struct_linstiff:
      dserror("Case not yet implemented");
    break;
    case calc_struct_nlnstiff:
      dserror("Case not yet implemented");
    break;
    case calc_struct_internalforce:
      dserror("Case not yet implemented");
    break;
    case calc_struct_linstiffmass:
      dserror("Case not yet implemented");
    break;
    case calc_struct_nlnstiffmass: // do mass, stiffness and internal forces
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual");
      if (disp==null || res==null) dserror("Cannot get state vectors displacement and/or residual");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);
      s8_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
    }
    break;
    case calc_struct_stress:
      dserror("Case not yet implemented");
    break;
    case calc_struct_eleload:
      dserror("Case not yet implemented");
    break;
    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;
    default:
      dserror("Unknown type of action for Shell8");
  }
  
  
  return 0;
}

/*----------------------------------------------------------------------*
 |  init the element (private)                               mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8::s8_initialize(DRT::Discretization& actdis)
{
  DSTraceHelper dst("Shell8::s8_initialize");  
  
  
  
  
  
  
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8::s8_nlnstiffmass(vector<int>&              lm, 
                                            vector<double>&           disp, 
                                            vector<double>&           residual,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* intforce,
                                            struct _MATERIAL*         material)
{
  DSTraceHelper dst("Shell8::s8_nlnstiffmass");  

  // for eas
  Epetra_SerialDenseMatrix P;
  Epetra_SerialDenseMatrix T;
  Epetra_SerialDenseMatrix Lt;
  Epetra_SerialDenseMatrix Dtild;
  Epetra_SerialDenseMatrix Dtildinv;
  Epetra_SerialDenseVector Rtild;

  const int numnode = NumNode();
  const int numdf   = 6;

  // gaussian points
  S8_DATA data;
  s8_integration_points(data);
  
  // eas
  if (nhyb_)
  {
    // init to zero
    P.Shape(12,nhyb_);
    T.Shape(12,12);
    Lt.Shape(nhyb_,numnode*numdf);
    Dtild.Shape(nhyb_,nhyb_);
    Dtildinv.Shape(nhyb_,nhyb_);
    Rtild.Size(nhyb_);
    // access history stuff stored in element
    vector<double>&           alfa        = alfa_;
    Epetra_SerialDenseMatrix& oldDtildinv = Dtildinv_;
    Epetra_SerialDenseMatrix& oldLt       = Lt_;
    vector<double>&           oldRtild    = Rtild_;
    /*---------------- make multiplication eashelp = oldLt * disp[kstep] */
    vector<double> eashelp(nhyb_);
    s8_YpluseqAx(eashelp,oldLt,disp,1.0,true);
    /*---------------------------------------- add old Rtilde to eashelp */
    for (int i=0; i<nhyb_; ++i) eashelp[i] += oldRtild[i];
    /*----------------- make multiplication alfa -= olDtildinv * eashelp */
    s8_YpluseqAx(alfa,oldDtildinv,eashelp,-1.0,false);
  } // if (nhyb_)
  
  // ------------------------------------ check calculation of mass matrix
  int imass=0;
  double density=0.0;
  if (massmatrix)
  {
    imass=1;
    s8_getdensity(mat,&density);
  }
  












  return;
} // DRT::Elements::Shell8::s8_nlnstiffmass













/*----------------------------------------------------------------------*
 |  y(I) = A(I,K)*x(K)*factor -----  y = A*x*factor         m.gee 12/06 |
 |  or                                                                  |
 |  y(I) += A(I,K)*x(K)*factor                                          |
 | (private)                                                            |
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8::s8_YpluseqAx(vector<double>& y,
                                         const Epetra_SerialDenseMatrix& A,
                                         const vector<double>& x, 
                                         const double factor, 
                                         const bool init)
{
  const int rdim = (int)y.size();
  const int ddim = (int)x.size();
  if (A.M()<rdim || A.N()<ddim) dserror("Mismatch in dimensions");
  
  if (init)
    for (int i=0; i<rdim; ++i) y[i] = 0.0;
  for (int i=0; i<rdim; ++i)
  {
    double sum = 0.0;
    for (int k=0; k<ddim; ++k) sum += A(i,k)*x[k];
    y[i] += sum*factor;
  }
  return;
}                    

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8::s8_integration_points(struct _S8_DATA& data)
{
  DSTraceHelper dst("Shell8::s8_integration_points");  
  
  const int numnode = NumNode();

  const double invsqrtthree = 1./sqrt(3.);
  const double sqrtthreeinvfive = sqrt(3./5.);
  const double wgt  = 5.0/9.0;
  const double wgt0 = 8.0/9.0;
  
  switch(ngp_[2])/*---------------- thickness direction t */
  {
    case 2:
      data.xgpt[0] = -invsqrtthree;
      data.xgpt[1] = invsqrtthree;
      data.xgpt[2] = 0.0;
      data.wgtt[0] =  1.0;
      data.wgtt[1] =  1.0;
      data.wgtt[2] =  0.0;
    break;
    default:
      dserror("Unknown no. of gaussian points in thickness direction");
    break;
  }
  
  // quad elements
  if (numnode==4 || numnode==8 || numnode==9)
  {
    switch(ngp_[0]) // r direction
    {
      case 1:
        data.xgpr[0] = 0.0;
        data.xgpr[1] = 0.0;
        data.xgpr[2] = 0.0;
        data.wgtr[0] = 2.0;
        data.wgtr[1] = 0.0;
        data.wgtr[2] = 0.0;
      break;
      case 2:
        data.xgpr[0] = -invsqrtthree;
        data.xgpr[1] =  invsqrtthree;
        data.xgpr[2] =  0.0;
        data.wgtr[0] =  1.0;
        data.wgtr[1] =  1.0;
        data.wgtr[2] =  0.0;
      break;
      case 3:
        data.xgpr[0] = -sqrtthreeinvfive;
        data.xgpr[1] =  0.0;
        data.xgpr[2] =  sqrtthreeinvfive;
        data.wgtr[0] =  wgt;
        data.wgtr[1] =  wgt0;
        data.wgtr[2] =  wgt;
      break;
      default:
        dserror("Unknown no. of gaussian points in r-direction");
      break;
    } // switch(ngp_[0]) // r direction
    
    switch(ngp_[1]) // s direction
    {
      case 1:
        data.xgps[0] = 0.0;
        data.xgps[1] = 0.0;
        data.xgps[2] = 0.0;
        data.wgts[0] = 2.0;
        data.wgts[1] = 0.0;
        data.wgts[2] = 0.0;
      break;
      case 2:
        data.xgps[0] = -invsqrtthree;
        data.xgps[1] =  invsqrtthree;
        data.xgps[2] =  0.0;
        data.wgts[0] =  1.0;
        data.wgts[1] =  1.0;
        data.wgts[2] =  0.0;
      break;
      case 3:
        data.xgps[0] = -sqrtthreeinvfive;
        data.xgps[1] =  0.0;
        data.xgps[2] =  sqrtthreeinvfive;
        data.wgts[0] =  wgt;
        data.wgts[1] =  wgt0;
        data.wgts[2] =  wgt;
      break;
      default:
        dserror("Unknown no. of gaussian points in s-direction");
      break;
    } // switch(ngp_[0]) // s direction
    
  } // if (numnode==4 || numnode==8 || numnode==9)

  else if (numnode==3 || numnode==6) // triangle elements
  {
    switch(ngptri_)
    {
      case 1:
      {
        const double third = 1.0/3.0;
        data.xgpr[0] =  third;
        data.xgpr[1] =  0.0;
        data.xgpr[2] =  0.0;
        data.xgps[0] =  third;
        data.xgps[1] =  0.0;
        data.xgps[2] =  0.0;
        data.wgtr[0] =  0.5;
        data.wgtr[1] =  0.0;
        data.wgtr[2] =  0.0;
        data.wgts[0] =  0.5;
        data.wgts[1] =  0.0;
        data.wgts[2] =  0.0;
      }
      break;
      case 3:
      {
        const double wgt = 1.0/6.0;
        data.xgpr[0] =  0.5;
        data.xgpr[1] =  0.5;
        data.xgpr[2] =  0.0;
        data.xgps[0] =  0.0;
        data.xgps[1] =  0.5;
        data.xgps[2] =  0.5;
        data.wgtr[0] =  wgt;
        data.wgtr[1] =  wgt;
        data.wgtr[2] =  wgt;
        data.wgts[0] =  wgt;
        data.wgts[1] =  wgt;
        data.wgts[2] =  wgt;
      }
      break;
      default:
        dserror("Unknown no. of gaussian points for triangle");
      break;
    } 
  } // else if (numnode==3 || numnode==6)
  
  return;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
