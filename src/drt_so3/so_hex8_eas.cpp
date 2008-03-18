/*!----------------------------------------------------------------------
\file so_hex8_eas.cpp
\brief Everything concerning EAS technology for so_hex8

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"


using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 |  initialize EAS data (private)                              maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_easinit()
{
  // EAS enhanced strain parameters at currently investigated load/time step
  Epetra_SerialDenseMatrix alpha(neas_,1);
  // EAS enhanced strain parameters of last converged load/time step
  Epetra_SerialDenseMatrix alphao(neas_,1);
  // EAS portion of internal forces, also called enhacement vector s or Rtilde
  Epetra_SerialDenseVector feas(neas_);
  // EAS matrix K_{alpha alpha}, also called Dtilde
  Epetra_SerialDenseMatrix invKaa(neas_,neas_);
  // EAS matrix K_{d alpha}
  Epetra_SerialDenseMatrix Kda(neas_,NUMDOF_SOH8);
  
  // save EAS data into element container easdata_
  data_.Add("alpha",alpha);
  data_.Add("alphao",alphao);
  data_.Add("feas",feas);
  data_.Add("invKaa",invKaa);
  data_.Add("Kda",Kda);
  
  return;
}

/*----------------------------------------------------------------------*
 |  setup of constant EAS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_eassetup(
          Epetra_SerialDenseMatrix** M_GP,    // M-matrix evaluated at GPs
          double& detJ0,                      // det of Jacobian at origin
          Epetra_SerialDenseMatrix& T0invT,   // maps M(origin) local to global
          const Epetra_SerialDenseMatrix& xrefe)    // material element coords
{
  // vector of df(origin)
  double df0_vector[NUMDOF_SOH8*NUMNOD_SOH8] =
               {-0.125,-0.125,-0.125,
                +0.125,-0.125,-0.125,
                +0.125,+0.125,-0.125,
                -0.125,+0.125,-0.125,
                -0.125,-0.125,+0.125,
                +0.125,-0.125,+0.125,
                +0.125,+0.125,+0.125,
                -0.125,+0.125,+0.125};
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  Epetra_DataAccess CV = Copy;
  Epetra_SerialDenseMatrix df0(CV,df0_vector,NUMDIM_SOH8,NUMDIM_SOH8,NUMNOD_SOH8);
  
  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  Epetra_SerialDenseMatrix jac0(NUMDIM_SOH8,NUMDIM_SOH8);
  jac0.Multiply('N','N',1.0,df0,xrefe,1.0);
  
  // compute determinant of Jacobian at origin by Sarrus' rule
  //double detJ0loc;
  detJ0 = jac0(0,0) * jac0(1,1) * jac0(2,2)
           + jac0(0,1) * jac0(1,2) * jac0(2,0)
           + jac0(0,2) * jac0(1,0) * jac0(2,1)
           - jac0(0,0) * jac0(1,2) * jac0(2,1)
           - jac0(0,1) * jac0(1,0) * jac0(2,2)
           - jac0(0,2) * jac0(1,1) * jac0(2,0);
  
  // first, build T0^T transformation matrix which maps the M-matrix
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  T0invT(0,0) = jac0(0,0) * jac0(0,0);
  T0invT(1,0) = jac0(1,0) * jac0(1,0);
  T0invT(2,0) = jac0(2,0) * jac0(2,0);
  T0invT(3,0) = 2 * jac0(0,0) * jac0(1,0);
  T0invT(4,0) = 2 * jac0(1,0) * jac0(2,0);
  T0invT(5,0) = 2 * jac0(0,0) * jac0(2,0);
  
  T0invT(0,1) = jac0(0,1) * jac0(0,1);
  T0invT(1,1) = jac0(1,1) * jac0(1,1);
  T0invT(2,1) = jac0(2,1) * jac0(2,1);
  T0invT(3,1) = 2 * jac0(0,1) * jac0(1,1);
  T0invT(4,1) = 2 * jac0(1,1) * jac0(2,1);
  T0invT(5,1) = 2 * jac0(0,1) * jac0(2,1);

  T0invT(0,2) = jac0(0,2) * jac0(0,2);
  T0invT(1,2) = jac0(1,2) * jac0(1,2);
  T0invT(2,2) = jac0(2,2) * jac0(2,2);
  T0invT(3,2) = 2 * jac0(0,2) * jac0(1,2);
  T0invT(4,2) = 2 * jac0(1,2) * jac0(2,2);
  T0invT(5,2) = 2 * jac0(0,2) * jac0(2,2);
  
  T0invT(0,3) = jac0(0,0) * jac0(0,1);
  T0invT(1,3) = jac0(1,0) * jac0(1,1);
  T0invT(2,3) = jac0(2,0) * jac0(2,1);
  T0invT(3,3) = jac0(0,0) * jac0(1,1) + jac0(1,0) * jac0(0,1);
  T0invT(4,3) = jac0(1,0) * jac0(2,1) + jac0(2,0) * jac0(1,1);
  T0invT(5,3) = jac0(0,0) * jac0(2,1) + jac0(2,0) * jac0(0,1);


  T0invT(0,4) = jac0(0,1) * jac0(0,2);
  T0invT(1,4) = jac0(1,1) * jac0(1,2);
  T0invT(2,4) = jac0(2,1) * jac0(2,2);
  T0invT(3,4) = jac0(0,1) * jac0(1,2) + jac0(1,1) * jac0(0,2);
  T0invT(4,4) = jac0(1,1) * jac0(2,2) + jac0(2,1) * jac0(1,2);
  T0invT(5,4) = jac0(0,1) * jac0(2,2) + jac0(2,1) * jac0(0,2);
  
  T0invT(0,5) = jac0(0,0) * jac0(0,2);
  T0invT(1,5) = jac0(1,0) * jac0(1,2);
  T0invT(2,5) = jac0(2,0) * jac0(2,2);
  T0invT(3,5) = jac0(0,0) * jac0(1,2) + jac0(1,0) * jac0(0,2);
  T0invT(4,5) = jac0(1,0) * jac0(2,2) + jac0(2,0) * jac0(1,2);
  T0invT(5,5) = jac0(0,0) * jac0(2,2) + jac0(2,0) * jac0(0,2);

  // now evaluate T0^{-T} with solver
  Epetra_SerialDenseSolver solve_for_inverseT0;
  solve_for_inverseT0.SetMatrix(T0invT);
  int err2 = solve_for_inverseT0.Factor();        
  int err = solve_for_inverseT0.Invert();
  if ((err != 0) && (err2!=0)) dserror("Inversion of T0inv (Jacobian0) failed");
  
  // build EAS interpolation matrix M, evaluated at the 8 GPs of so_hex8
  static Epetra_SerialDenseMatrix M(NUMSTR_SOH8*NUMGPT_SOH8,neas_);
  static bool M_eval;

  if (M_eval==true){          // if true M already evaluated
      *M_GP = &M;             // return adress of static object to target of pointer
    return;
  } else {
    // (r,s,t) gp-locations of fully integrated linear 8-node Hex
    const double gploc    = 1.0/sqrt(3.0);    // gp sampling point value for linear fct
    const double r[NUMGPT_SOH8] = {-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc,-gploc};
    const double s[NUMGPT_SOH8] = {-gploc,-gploc, gploc, gploc,-gploc,-gploc, gploc, gploc};
    const double t[NUMGPT_SOH8] = {-gploc,-gploc,-gploc,-gploc, gploc, gploc, gploc, gploc};

    // fill up M at each gp
    if (eastype_ == soh8_easmild) {
      /* easmild is the EAS interpolation of 9 modes, based on
      **            r 0 0   0 0 0 0 0 0
      **            0 s 0   0 0 0 0 0 0
      **    M =     0 0 t   0 0 0 0 0 0
      **            0 0 0   r s 0 0 0 0
      **            0 0 0   0 0 s t 0 0 
      **            0 0 0   0 0 0 0 r t 
      */
      for (int i=0; i<NUMGPT_SOH8; ++i) {
        M(i*NUMSTR_SOH8+0,0) = r[i];
        M(i*NUMSTR_SOH8+1,1) = s[i];
        M(i*NUMSTR_SOH8+2,2) = t[i];
        
        M(i*NUMSTR_SOH8+3,3) = r[i]; M(i*NUMSTR_SOH8+3,4) = s[i];
        M(i*NUMSTR_SOH8+4,5) = s[i]; M(i*NUMSTR_SOH8+4,6) = t[i];
        M(i*NUMSTR_SOH8+5,7) = r[i]; M(i*NUMSTR_SOH8+5,8) = t[i];
      }
    
      // return adress of just evaluated matrix
      *M_GP = &M;            // return adress of static object to target of pointer
      M_eval = true;         // now the array is filled statically
    } else if (eastype_ == soh8_easfull) {
      /* easfull is the EAS interpolation of 21 modes, based on
      **            r 0 0   0 0 0 0 0 0   0  0  0  0  0  0   rs rt 0  0  0  0  
      **            0 s 0   0 0 0 0 0 0   0  0  0  0  0  0   0  0  rs st 0  0  
      **    M =     0 0 t   0 0 0 0 0 0   0  0  0  0  0  0   0  0  0  0  rt st 
      **            0 0 0   r s 0 0 0 0   rt st 0  0  0  0   0  0  0  0  0  0
      **            0 0 0   0 0 s t 0 0   0  0  rs rt 0  0   0  0  0  0  0  0
      **            0 0 0   0 0 0 0 r t   0  0  0  0  rs st  0  0  0  0  0  0 
      */
      for (int i=0; i<NUMGPT_SOH8; ++i) {
        M(i*NUMSTR_SOH8+0,0) = r[i];        M(i*NUMSTR_SOH8+0,15) = r[i]*s[i]; M(i*NUMSTR_SOH8+0,16) = r[i]*t[i];
        M(i*NUMSTR_SOH8+1,1) = s[i];        M(i*NUMSTR_SOH8+1,17) = r[i]*s[i]; M(i*NUMSTR_SOH8+1,18) = s[i]*t[i];
        M(i*NUMSTR_SOH8+2,2) = t[i];        M(i*NUMSTR_SOH8+2,19) = r[i]*t[i]; M(i*NUMSTR_SOH8+2,20) = s[i]*t[i];
        
        M(i*NUMSTR_SOH8+3,3) = r[i]; M(i*NUMSTR_SOH8+3,4) = s[i];   M(i*NUMSTR_SOH8+3, 9) = r[i]*t[i]; M(i*NUMSTR_SOH8+3,10) = s[i]*t[i];
        M(i*NUMSTR_SOH8+4,5) = s[i]; M(i*NUMSTR_SOH8+4,6) = t[i];   M(i*NUMSTR_SOH8+4,11) = r[i]*s[i]; M(i*NUMSTR_SOH8+4,12) = r[i]*t[i];
        M(i*NUMSTR_SOH8+5,7) = r[i]; M(i*NUMSTR_SOH8+5,8) = t[i];   M(i*NUMSTR_SOH8+5,13) = r[i]*s[i]; M(i*NUMSTR_SOH8+5,14) = s[i]*t[i];
      }
      // return adress of just evaluated matrix      
      *M_GP = &M;            // return adress of static object to target of pointer
      M_eval = true;         // now the array is filled statically
    } else if (eastype_ == soh8_eassosh8) {
      /* eassosh8 is the EAS interpolation for the Solid-Shell with t=thickness dir.
      ** consisting of 7 modes, based on
      **            r 0 0   0 0 0  0 
      **            0 s 0   0 0 0  0 
      **    M =     0 0 t   0 0 rt st
      **            0 0 0   r s 0  0
      **            0 0 0   0 0 0  0 
      **            0 0 0   0 0 0  0 
      */
      for (int i=0; i<NUMGPT_SOH8; ++i) {
        M(i*NUMSTR_SOH8+0,0) = r[i];
        M(i*NUMSTR_SOH8+1,1) = s[i];
        M(i*NUMSTR_SOH8+2,2) = t[i]; M(i*NUMSTR_SOH8+2,5) = r[i]*t[i]; M(i*NUMSTR_SOH8+2,6) = s[i]*t[i];
        
        M(i*NUMSTR_SOH8+3,3) = r[i]; M(i*NUMSTR_SOH8+3,4) = s[i];
      }
      // return adress of just evaluated matrix      
      *M_GP = &M;            // return adress of static object to target of pointer
      M_eval = true;         // now the array is filled statically
    } else {
    dserror("eastype not implemented");
    }
  }
} // end of soh8_eassetup
  
      
  
  

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
