/*!----------------------------------------------------------------------
\file visconeohooke.cpp
\brief

<pre>
Maintainer: Moritz Frenzel & Thomas Kloeppel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "visconeohooke.H"

extern struct _MATERIAL *mat;


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoNeoHooke::ViscoNeoHooke()
  : matdata_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoNeoHooke::ViscoNeoHooke(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];     // unpack pointer to my specific matdata_

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


/*----------------------------------------------------------------------*
 |  Return density                                (public)         05/08|
 *----------------------------------------------------------------------*/
double MAT::ViscoNeoHooke::Density()
{
  return matdata_->m.visconeohooke->density;  // density, returned to evaluate mass matrix
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         05/08|
 *----------------------------------------------------------------------*

*/

void MAT::ViscoNeoHooke::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                       Epetra_SerialDenseMatrix* cmat,
                                       Epetra_SerialDenseVector* stress,
                                       const double time)

{
  // get material parameters
  double E_s  = matdata_->m.visconeohooke->youngs_slow;
  double nue  = matdata_->m.visconeohooke->poisson;
  double E_f  = matdata_->m.visconeohooke->youngs_fast;  
  double tau  = matdata_->m.visconeohooke->relax;
  double theta= 0.0;//matdata_->m.visconeohooke->theta;
  
  // evaluate "alpha" factors which distribute stress or stiffnes between parallel springs
  // sum_0^i alpha_j = 1
  if (E_f < E_s) dserror("Wrong ratio between fast and slow Young's modulus");
  const double alpha0 = E_s / E_f;
  const double alpha1 = 1.0 - alpha0;
  
  // evaluate Lame constants, bulk modulus
  const double lambda = nue*E_f / ((1.0+nue)*(1.0-2*nue));
  const double mue = E_f / (2*(1.0+nue));
  const double kappa = lambda + 2.0/3.0 * mue;
  
  // remark: in the following we use Voigt-notation
  // and keep the strain-factor-2 at the strains
  
  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  Epetra_SerialDenseVector Id(6);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  Epetra_SerialDenseVector C(*glstrain);
  C.Scale(2.0);
  C += Id;

  // invariants
  const double I1 = C(0) + C(1) + C(2);  // 1st invariant, trace
  const double I3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  const double J = sqrt(I3);
  const double I3invcubroot = pow(I3,-1.0/3.0);
  
  // invert C
  Epetra_SerialDenseVector Cinv(6);

  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);

  Cinv.Scale(1.0/I3);
  
  // Split into volumetric and deviatoric parts. Viscosity affects only deviatoric part
  
  //Volumetric part of PK2 stress
  Epetra_SerialDenseVector SVol(Cinv);
  SVol.Scale(kappa*(J-1.0)*J);
  *stress+=SVol;
  
  //Deviatoric elastic part (2 d W^dev/d C)
  Epetra_SerialDenseVector SDevEla(Cinv);
  SDevEla.Scale(-1.0/3.0*I1);
  SDevEla+=Id;
  SDevEla.Scale(mue*I3invcubroot);  //mue*I3^(-1/3) (Id-1/3*I1*Cinv)
    
  //visco
  return;
}


#endif
