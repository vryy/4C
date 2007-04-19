/*!----------------------------------------------------------------------
\file wall1_evaluate.cpp
\brief

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "wall1.H"
#include "drt_discret.H"
#include "drt_utils.H"
#include "drt_exporter.H"
#include "drt_dserror.H"
#include "linalg_utils.H"

extern "C" 
{
#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
#include "../wall1/wall1_prototypes.h"
}
#include "dstrc.H"

/*----------------------------------------------------------------------*
 |                                                        mgit 03/07    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Elements::Wall1::Evaluate(ParameterList& params, 
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DSTraceHelper dst("Wall1::Evaluate");  
  DRT::Elements::Wall1::ActionType act = Wall1::calc_none;
  
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = Wall1::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Wall1::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Wall1::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = Wall1::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = Wall1::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = Wall1::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Wall1::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Wall1::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Wall1::calc_struct_update_istep;
  else dserror("Unknown type of action for Wall1");

  // get the material law
  MATERIAL* actmat = &(mat[material_-1]);
  switch(act)
  {
    case Wall1::calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      w1_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
    }
    break;
    case Wall1::calc_struct_nlnstiffmass:
    case Wall1::calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::Utils::ExtractMyValues(*res,myres,lm);
      w1_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,actmat);
    }
    break;
    default:
      dserror("Unknown type of action for Wall1 %d", act);
  }
  return 0;
 
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  mgit 04/07|
 *----------------------------------------------------------------------*/

int DRT::Elements::Wall1::EvaluateNeumann(ParameterList& params, 
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::w1_nlnstiffmass(vector<int>&               lm, 
                                            vector<double>&           disp, 
                                            vector<double>&           residual,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force,
                                            struct _MATERIAL*         material)
{
  DSTraceHelper dst("Wall1::w1_nlnstiffmass");  
  const int numnode = NumNode();
  const int numdf   = 2;
  int       ngauss  = 0;
  const int nd      = numnode*numdf;

  // general arrays
  vector<double>           funct(numnode);        
  Epetra_SerialDenseMatrix deriv;
  deriv.Shape(2,numnode);
  Epetra_SerialDenseMatrix xjm;
  xjm.Shape(2,2);
  Epetra_SerialDenseMatrix bop;
  bop.Shape(12,nd);//????
  Epetra_SerialDenseVector intforce;
  intforce.Size(nd);//????
  double det; 
  double xrefe[2][MAXNOD_WALL1];
  double xcure[2][MAXNOD_WALL1];

  // gaussian points
  W1_DATA w1data;
  w1_integration_points(w1data);
  

  // ------------------------------------ check calculation of mass matrix
  int imass=0;
  double density=0.0;
  if (massmatrix)
  {
    imass=1;
    w1_getdensity(material,&density);
  }

  const int nir = ngp_[0];
  const int nis = ngp_[1];
  const int iel = numnode;

  /*----------------------------------------------------- geometry update */
  for (int k=0; k<iel; ++k)
  {
    
    xrefe[0][k] = Nodes()[k]->X()[0];
    xrefe[1][k] = Nodes()[k]->X()[1];
    
    xcure[0][k] = xrefe[0][k] + disp[k*numdf+0];
    xcure[1][k] = xrefe[1][k] + disp[k*numdf+1];

  }

  /*=================================================== integration loops */
  for (int lr=0; lr<nir; ++lr)
  {
    /*================================== gaussian point and weight at it */
    const double e1   = w1data.xgrr[lr];
      for (int ls=0; ls<nis; ++ls)
    {
      const double e2   = w1data.xgss[ls];
      /*-------------------- shape functions at gp e1,e2 on mid surface */
      w1_shapefunctions(funct,deriv,e1,e2,iel,1);
      /*--------------------------------------- compute jacobian Matrix */ 
      w1_jacobianmatrix(xrefe,deriv,xjm,&det,iel);
      
       ngauss++;
    } // for (int ls=0; ls<nis; ++ls)
  } // for (int lr=0; lr<nir; ++lr)

return;
}
/*----------------------------------------------------------------------*
 |  evaluate the element integration points (private)        mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::w1_integration_points(struct _W1_DATA& data)
{
  DSTraceHelper dst("Wall1::w1_integration_points");  
  
  const int numnode = NumNode();

  const double invsqrtthree = 1./sqrt(3.);
  const double sqrtthreeinvfive = sqrt(3./5.);
  const double wgt  = 5.0/9.0;
  const double wgt0 = 8.0/9.0;
  
 
  // quad elements
  if (numnode==4 || numnode==8 || numnode==9)
  {
    switch(ngp_[0]) // r direction
    {
      case 1:
        data.xgrr[0] = 0.0;
        data.xgrr[1] = 0.0;
        data.xgrr[2] = 0.0;
        data.wgtr[0] = 2.0;
        data.wgtr[1] = 0.0;
        data.wgtr[2] = 0.0;
      break;
      case 2:
        data.xgrr[0] = -invsqrtthree;
        data.xgrr[1] =  invsqrtthree;
        data.xgrr[2] =  0.0;
        data.wgtr[0] =  1.0;
        data.wgtr[1] =  1.0;
        data.wgtr[2] =  0.0;
      break;
      case 3:
        data.xgrr[0] = -sqrtthreeinvfive;
        data.xgrr[1] =  0.0;
        data.xgrr[2] =  sqrtthreeinvfive;
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
        data.xgss[0] = 0.0;
        data.xgss[1] = 0.0;
        data.xgss[2] = 0.0;
        data.wgts[0] = 2.0;
        data.wgts[1] = 0.0;
        data.wgts[2] = 0.0;
      break;
      case 2:
        data.xgss[0] = -invsqrtthree;
        data.xgss[1] =  invsqrtthree;
        data.xgss[2] =  0.0;
        data.wgts[0] =  1.0;
        data.wgts[1] =  1.0;
        data.wgts[2] =  0.0;
      break;
      case 3:
        data.xgss[0] = -sqrtthreeinvfive;
        data.xgss[1] =  0.0;
        data.xgss[2] =  sqrtthreeinvfive;
        data.wgts[0] =  wgt;
        data.wgts[1] =  wgt0;
        data.wgts[2] =  wgt;
      break;
      default:
        dserror("Unknown no. of gaussian points in s-direction");
      break;
    } // switch(ngp_[0]) // s direction
    
  } // if (numnode==4 || numnode==8 || numnode==9)

//  else if (numnode==3 || numnode==6) // triangle elements
//  {
//    switch(ngptri_)
//    {
//      case 1:
//      {
//        const double third = 1.0/3.0;
//        data.xgrr[0] =  third;
//        data.xgrr[1] =  0.0;
//        data.xgrr[2] =  0.0;
//        data.xgss[0] =  third;
//        data.xgss[1] =  0.0;
//        data.xgss[2] =  0.0;
//        data.wgtr[0] =  0.5;
//        data.wgtr[1] =  0.0;
//        data.wgtr[2] =  0.0;
//        data.wgts[0] =  0.5;
//        data.wgts[1] =  0.0;
//        data.wgts[2] =  0.0;
//      }
//      break;
//      case 3:
//      {
//        const double wgt = 1.0/6.0;
//        data.xgrr[0] =  0.5;
//        data.xgrr[1] =  0.5;
//        data.xgrr[2] =  0.0;
//        data.xgss[0] =  0.0;
//        data.xgss[1] =  0.5;
//        data.xgss[2] =  0.5;
//        data.wgtr[0] =  wgt;
//        data.wgtr[1] =  wgt;
//        data.wgtr[2] =  wgt;
//        data.wgts[0] =  wgt;
//        data.wgts[1] =  wgt;
//        data.wgts[2] =  wgt;
//      }
//      break;
//      default:
//        dserror("Unknown no. of gaussian points for triangle");
//      break;
//    } 
//  } // else if (numnode==3 || numnode==6)
  
  return;
}

/*----------------------------------------------------------------------*
 |  shape functions and derivatives (private)                mgit 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Wall1::w1_shapefunctions(
                             vector<double>& funct, 
                             Epetra_SerialDenseMatrix& deriv,
                             const double r, const double s, const int numnode, 
                             const int doderiv) const
{
  DSTraceHelper dst("Wall1::w1_shapefunctions");  

  const double q12 = 0.5;
  const double q14 = 0.25;
  const double rr = r*r;
  const double ss = s*s;
  const double rp = 1.0+r;
  const double rm = 1.0-r;
  const double sp = 1.0+s;
  const double sm = 1.0-s;
  const double r2 = 1.0-rr;
  const double s2 = 1.0-ss;
  int i;
  int ii; 

  //const double t;
   
  switch(numnode)
  {
    case 4:
    {
      funct[0] = q14*rp*sp;
      funct[1] = q14*rm*sp;
      funct[2] = q14*rm*sm;
      funct[3] = q14*rp*sm;
      if (doderiv)
      {
        deriv(0,0)= q14*sp;
        deriv(0,1)=-q14*sp;
        deriv(0,2)=-q14*sm;
        deriv(0,3)= q14*sm;
        deriv(1,0)= q14*rp;
        deriv(1,1)= q14*rm;
        deriv(1,2)=-q14*rm;
        deriv(1,3)=-q14*rp;
      }
      return;
    }
    break;
    case 8:
    {
      funct[0] = q14*rp*sp;
      funct[1] = q14*rm*sp;
      funct[2] = q14*rm*sm;
      funct[3] = q14*rp*sm;
      funct[4] = q12*r2*sp;
      funct[5] = q12*rm*s2;
      funct[6] = q12*r2*sm;
      funct[7] = q12*rp*s2;
      funct[0] = funct[0] - q12*(funct[4] + funct[7]);
      if (doderiv)
      {
         deriv(0,0)= q14*sp;
         deriv(0,1)=-q14*sp;
         deriv(0,2)=-q14*sm;
         deriv(0,3)= q14*sm;
         deriv(1,0)= q14*rp;
         deriv(1,1)= q14*rm;
         deriv(1,2)=-q14*rm;
         deriv(1,3)=-q14*rp;
         deriv(0,4)=-1*r*sp;
         deriv(0,5)=-q12*  s2;
         deriv(0,6)=-1*r*sm;
         deriv(0,7)= q12*  s2;
         deriv(1,4)= q12*r2  ;
         deriv(1,5)=-1*rm*s;
         deriv(1,6)=-q12*r2  ;
         deriv(1,7)=-1*rp*s;

         deriv[0][0]=deriv[0][0] - q12*(deriv[0][4] + deriv[0][7]);
         deriv[1][0]=deriv[1][0] - q12*(deriv[1][4] + deriv[1][7]);
      }
      for (i=1; i<=3; i++)
      {
         ii=i + 3;
         funct[i]=funct[i] - q12*(funct[ii] + funct[ii+1]);
         if (doderiv)              /*--- check for derivative evaluation ---*/
         {
             deriv(0,i)=deriv(0,i) - q12*(deriv(0,ii) + deriv(0,ii+1));
             deriv(1,i)=deriv(1,i) - q12*(deriv(1,ii) + deriv(1,ii+1));
         }
      }      
      return;
    }
    break;
    case 9:
    { 
      const double rh  = q12*r;
      const double sh  = q12*s;
      const double rs  = rh*sh;
      const double rhp = r+q12;
      const double rhm = r-q12;
      const double shp = s+q12;
      const double shm = s-q12;
      funct[0] = rs*rp*sp;
      funct[1] =-rs*rm*sp;
      funct[2] = rs*rm*sm;
      funct[3] =-rs*rp*sm;
      funct[4] = sh*sp*r2;
      funct[5] =-rh*rm*s2;
      funct[6] =-sh*sm*r2;
      funct[7] = rh*rp*s2;
      funct[8] = r2*s2;
      if (doderiv==1)
      {
         deriv(0,0)= rhp*sh*sp;
         deriv(0,1)= rhm*sh*sp;
         deriv(0,2)=-rhm*sh*sm;
         deriv(0,3)=-rhp*sh*sm;
         deriv(0,4)=-2.0*r*sh*sp;
         deriv(0,5)= rhm*s2;
         deriv(0,6)= 2.0*r*sh*sm;
         deriv(0,7)= rhp*s2;
         deriv(0,8)=-2.0*r*s2;
         deriv(1,0)= shp*rh*rp;
         deriv(1,1)=-shp*rh*rm;
         deriv(1,2)=-shm*rh*rm;
         deriv(1,3)= shm*rh*rp;
         deriv(1,4)= shp*r2;
         deriv(1,5)= 2.0*s*rh*rm;
         deriv(1,6)= shm*r2;
         deriv(1,7)=-2.0*s*rh*rp;
         deriv(1,8)=-2.0*s*r2;
      }
      return;      
   }
   break;

//*------------------------------------------------- triangular elements */
//case tri3: /* LINEAR shape functions and their natural derivatives -----*/
///*----------------------------------------------------------------------*/
//   funct[0]=ONE-r-s;
//   funct[1]=r;
//   funct[2]=s;
//
//   if(option==1) /* --> first derivative evaluation */
//   {
//      deriv[0][0]=-ONE;
//      deriv[1][0]=-ONE;
//      deriv[0][1]= ONE;
//      deriv[1][1]=ZERO;
//      deriv[0][2]=ZERO;
//      deriv[1][2]= ONE;
//   } /* endif (option==1) */
//break;
///*-------------------------------------------------------------------------*/
//case tri6: /* Quadratic shape functions and their natural derivatives -----*/
//    t = ONE-r-s;
//
//   funct[0] = t*(TWO*t-ONE);
//    funct[1] = r*(TWO*r-ONE);
//    funct[2] = s*(TWO*s-ONE);
//    funct[3] = FOUR*r*t;
//    funct[4] = FOUR*r*s;
//    funct[5] = FOUR*s*t;
//    
//    if (option == 1) /* --> first derivative evaluation */
//    {
//        /* first natural derivative of funct[0] with respect to r */
//        deriv[0][0] = -FOUR*t + ONE;
//        /* first natural derivative of funct[0] with respect to s */
//        deriv[1][0] = -FOUR*t + ONE;
//        deriv[0][1] = FOUR*r - ONE;
//        deriv[1][1] = ZERO;
//        deriv[0][2] = ZERO;
//        deriv[1][2] = FOUR*s - ONE;
//        deriv[0][3] = FOUR*t - FOUR*r;
//        deriv[1][3] = -FOUR*r;
//        deriv[0][4] = FOUR*s;
//        deriv[1][4] = FOUR*r;
//        deriv[0][5] = -FOUR*s;
//        deriv[1][5] = FOUR*t - FOUR*s;
//    } /* end if (option==1) */
//break;
default:
   dserror("Unknown no. of nodes %d to wall1 element",numnode);
break;
} /* end of switch typ */
/*----------------------------------------------------------------------*/

  return;  

} /* DRT::Elements::Wall1::w1_shapefunctions */

/*----------------------------------------------------------------------*
 |  jacobian matrix (private)                                  mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::Elements::Wall1::w1_jacobianmatrix(double xrefe[2][MAXNOD_WALL1],
                          const Epetra_SerialDenseMatrix& deriv,
                          Epetra_SerialDenseMatrix& xjm,
			  double* det,
                          const int iel) 
{

   memset(xjm.A(),0,xjm.N()*xjm.M()*sizeof(double)); 
 
   for (int k=0; k<iel; k++)
   {
        xjm(0,0) += deriv(0,k) * xrefe[0][k];
        xjm(0,1) += deriv(0,k) * xrefe[1][k];
        xjm(1,0) += deriv(1,k) * xrefe[0][k];
        xjm(1,1) += deriv(1,k) * xrefe[1][k];
   }
 cout << xjm; exit(0);

   return;
} // DRT::Elements::Wall1::w1_jacobianmatrix

  







#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
