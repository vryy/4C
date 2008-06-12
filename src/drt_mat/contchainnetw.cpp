/*!----------------------------------------------------------------------
\file contchainnetw.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
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
#include "contchainnetw.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_utils.H"


extern struct _MATERIAL *mat;


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         06/08|
 *----------------------------------------------------------------------*/
MAT::ContChainNetw::ContChainNetw()
  : matdata_(NULL)
{
  isinit_=false;
  l1_=rcp(new vector<double>);
  l2_=rcp(new vector<double>);
  l3_=rcp(new vector<double>);
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          06/08|
 *----------------------------------------------------------------------*/
MAT::ContChainNetw::ContChainNetw(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);
  AddtoPack(data,l1_);
  AddtoPack(data,l2_);
  AddtoPack(data,l3_);
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Unpack(const vector<char>& data)
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

  ExtractfromPack(position,data,l1_);
  ExtractfromPack(position,data,l2_);
  ExtractfromPack(position,data,l3_);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Initialize(const int numgp) 
{
  const double isotropy  = 1/sqrt(3.0) * matdata_->m.contchainnetw->r0;

  l1_= rcp(new vector<double> (numgp,isotropy));
  l2_= rcp(new vector<double> (numgp,isotropy));
  l3_= rcp(new vector<double> (numgp,isotropy));
  isinit_ = true;
  
  return ;
  
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         06/08|
 *----------------------------------------------------------------------*

*/

void MAT::ContChainNetw::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                  const int gp,
                                  Teuchos::ParameterList& params,
                                  Epetra_SerialDenseMatrix* cmat,
                                  Epetra_SerialDenseVector* stress)

{
  // bulk (isotropic) NeoHooke material parameters (Lame constants)
  const double lambda = matdata_->m.contchainnetw->lambda;
  const double mue = matdata_->m.contchainnetw->mue;
  // chain network unit cell material parameters
  const double nchain = matdata_->m.contchainnetw->nchain; // chain density ~= cell stiffness
  const double abstemp = matdata_->m.contchainnetw->abstemp; // absolute temperature (K)
  const double L = matdata_->m.contchainnetw->contl_l;  // chain contour length
  const double A = matdata_->m.contchainnetw->persl_a;  // chain persistence length
  const double r0 = matdata_->m.contchainnetw->r0;      // initial chain length
  const double boltzmann = 1.3806503E-23;
  
  // the chain stiffness factor
  const double chn_stiffact = boltzmann * abstemp * nchain / (4*A);
  
  // scalar to arrive at stressfree reference conf
  const double stressfree = - chn_stiffact * ( 1.0/L + 1.0/(4.0*r0*(1.0-r0/L)*(1.0-r0/L)) - 1.0/(4.0*r0) );
  
  const double time = params.get("total time",-1.0);
  const double kappa = matdata_->m.contchainnetw->relax; // relaxation time for remodeling
  const double decay = exp(-kappa*time);

  // initial cell dimensions (isotropy)
  const double isotropy  = 1/sqrt(3.0) * r0;
  const double l10 = isotropy;
  const double l20 = isotropy;
  const double l30 = isotropy;

  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  Epetra_SerialDenseVector Id(6);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  Epetra_SerialDenseVector C(*glstrain);
  C.Scale(2.0);
  C += Id;

  // invariants
  const double IC3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  const double J = sqrt(IC3);
  const double lJ = log(J);
  // 'non-standard' invariants representing stretch^2 in n0_i direction
  double I1 = C(0);
  double I2 = C(1);
  double I3 = C(2);
  
  // invert C
  Epetra_SerialDenseVector Cinv(6);

  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);

  Cinv.Scale(1.0/IC3);

  // isotropic part: NeoHooke  ************************************************
  // W = 1/2 lambda ln^2(J) + 1/2 mue (I1-3) - mue ln(J)
  // S = (lambda ln(J) - mue) Cinv + mue Id
  // Elasticity = lambda (Cinv x Cinv) + 2(mue - lambda ln(J))(Cinv o Cinv) 
  Epetra_SerialDenseVector Siso1(Cinv);
  Siso1.Scale(lambda*lJ-mue);
  *stress += Siso1;
  Siso1 = Id;
  Siso1.Scale(mue);
  *stress += Siso1;
  
  AddtoCmatHolzapfelProduct((*cmat),Cinv,2*(mue-lambda*lJ));
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*cmat)(i,j) += lambda * Cinv(i) * Cinv(j); // add lambda Cinv x Cinv
    }
  }
  // end of isotropic part ****************************************************
  
  // trial S to compute eigenvalues *******************************************
  Epetra_SerialDenseMatrix Strial(3,3);
  double l1sq = l1_->at(gp)*l1_->at(gp);
  double l2sq = l2_->at(gp)*l2_->at(gp);
  double l3sq = l3_->at(gp)*l3_->at(gp);
  
  // structural tensors
  Epetra_SerialDenseVector N01(6); N01(0) = 1.0;
  Epetra_SerialDenseVector N02(6); N02(1) = 1.0;
  Epetra_SerialDenseVector N03(6); N03(2) = 1.0;
  
  // current chain length
  double r = sqrt(I1*l1sq + I2*l2sq + I3*l3sq);
  
  double s_chn = chn_stiffact*(4.0/L + 1.0/(r*(1.0-r/L)*(1.0-r/L)) - 1.0/r);
  
  Strial(0,0) = (*stress)(0) + l1sq*s_chn + l1sq*I1*4*stressfree;
  Strial(1,1) = (*stress)(1) + l2sq*s_chn + l2sq*I2*4*stressfree;
  Strial(2,2) = (*stress)(2) + l3sq*s_chn + l3sq*I3*4*stressfree;
  Strial(0,1) = (*stress)(3); Strial(1,0) = (*stress)(3);
  Strial(1,2) = (*stress)(4); Strial(2,1) = (*stress)(4);
  Strial(0,2) = (*stress)(5); Strial(2,0) = (*stress)(5);
  
  Epetra_SerialDenseVector eig_sp(3);  // lambda^(sigma+)
  LINALG::SymmetricEigenValues(Strial,eig_sp);
  for (int i = 0; i < 3; ++i) {
    if (eig_sp(i) > 0.0) eig_sp(i) = 0.0;
    else eig_sp(i) = 1.0;
  }
  
  // update cell dimensions (remodeling!)
  const double l1 = (eig_sp(0) - l10/r0)*(1 - decay)*r0 + l10;
  const double l2 = (eig_sp(1) - l20/r0)*(1 - decay)*r0 + l20;
  const double l3 = (eig_sp(2) - l30/r0)*(1 - decay)*r0 + l30;
  
  // evaluate chain stress and add to isotropic stress
  l1sq = l1*l1;
  l2sq = l2*l2;
  l3sq = l3*l3;
  r = sqrt(I1*l1sq + I2*l2sq + I3*l3sq);
  s_chn = chn_stiffact*(4.0/L + 1.0/(r*(1.0-r/L)*(1.0-r/L)) - 1.0/r);
  (*stress)(0) += l1sq*s_chn + l1sq*I1*4*stressfree;
  (*stress)(1) += l2sq*s_chn + l2sq*I2*4*stressfree;
  (*stress)(2) += l3sq*s_chn + l3sq*I3*4*stressfree;
  
  //evaluate chain tangent and add to isotropic cmat
  double c_chn = s_chn/(r*r*r) * (1.0 - 1.0/((1.0-r/L)*(1.0-r/L)) + 2.0*r/(L*(1.0-r/L)*(1.0-r/L)*(1.0-r/L)) );
  (*cmat)(0,0) += c_chn*(l1*l1*l1*l1) - 8.0*stressfree*l1*l1/(I1*I1);
  (*cmat)(1,1) += c_chn*(l2*l2*l2*l2) - 8.0*stressfree*l2*l2/(I2*I2);
  (*cmat)(2,2) += c_chn*(l3*l3*l3*l3) - 8.0*stressfree*l3*l3/(I3*I3);
  
  // store history values
  l1_->at(gp) = l1;
  l2_->at(gp) = l2;
  l3_->at(gp) = l3;
  
  return;
}


#endif
