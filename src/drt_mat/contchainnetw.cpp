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
  li_ = rcp(new vector<vector<double> >);
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
  int histsize;
  if (!Initialized())
  {
    histsize=0;
  }
  else 
  {
    histsize = li_->size();
  }
  AddtoPack(data,histsize);  // lenght of history vector(s)
  for (int var = 0; var < histsize; ++var) 
  {
    AddtoPack(data,li_->at(var));
  }
  
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

  // history data
  isinit_ = true;
  int histsize;
  ExtractfromPack(position,data,histsize);
  
  if (histsize == 0) isinit_=false;
  li_ = rcp(new vector<vector<double> >);
  for (int var = 0; var < histsize; ++var) {
    vector<double> li;
    ExtractfromPack(position,data,li);
    li_->push_back(li);
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Initialize(const int numgp) 
{
  const double isotropy  = 1/sqrt(3.0) * matdata_->m.contchainnetw->r0;
  srand ( time(NULL) );

  li0_ = rcp(new vector<vector<double> > (numgp));
  li_ = rcp(new vector<vector<double> > (numgp));
  for(int gp=0; gp<numgp; ++gp){
    li0_->at(gp).resize(3);
    li_->at(gp).resize(3);
//    for (int i = 0; i < 3; ++i) {
//      li0_->at(gp)[i] = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
//    }
    li0_->at(gp)[0] = isotropy;
    li0_->at(gp)[1] = 0.0;
    li0_->at(gp)[2] = 0.0;
    for (int i = 0; i < 3; ++i) li_->at(gp)[i] = li0_->at(gp)[i];
  }

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
                                  Epetra_SerialDenseVector* stress,
                                  int eleId)

{
  // bulk (isotropic) NeoHooke material parameters (Lame constants)
  const double lambda = matdata_->m.contchainnetw->lambda;
  const double mue = matdata_->m.contchainnetw->mue;
  // chain network unit cell material parameters
  const double nchain = matdata_->m.contchainnetw->nchain; // chain density ~= cell stiffness
  const double abstemp = matdata_->m.contchainnetw->abstemp; // absolute temperature (K)
  const double L = matdata_->m.contchainnetw->contl_l;  // chain contour length
  const double A = matdata_->m.contchainnetw->persl_a;  // chain persistence length
  double r0 = matdata_->m.contchainnetw->r0;      // initial chain length
  const double boltzmann = 1.3806503E-23;
  
  // the chain stiffness factor
  const double chn_stiffact = boltzmann * abstemp * nchain / (4.0*A);
  
  const double time = params.get("total time",-1.0);
  double rem_time = time-1.0;
  const double kappa = matdata_->m.contchainnetw->relax; // relaxation time for remodeling
  const double decay = exp(-kappa*rem_time);
  
  const int dim = 3;

  // initial cell dimensions (isotropy)
  vector<double> li0sq(dim);
//  li0[1] = isotropy;
//  li0[0] = 0.0;
//  li0[2] = 0.0;
  for (int i=0; i<dim; ++i) li0sq[i] = li0_->at(gp)[i]*li0_->at(gp)[i];
  r0 = sqrt(li0sq[0] + li0sq[1] + li0sq[2]);

  // scalar to arrive at stressfree reference conf
  const double stressfree = - chn_stiffact * ( 1.0/L + 1.0/(4.0*r0*(1.0-r0/L)*(1.0-r0/L)) - 1.0/(4.0*r0) );

  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  Epetra_SerialDenseVector Id(6);
  for (int i = 0; i < dim; i++) Id(i) = 1.0;
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
//  if (gp==0){
//    cout << "before isotrop" << endl <<(*stress);
//  }
  *stress += Siso1;
  Siso1 = Id;
  Siso1.Scale(mue);
  *stress += Siso1;
//  if (gp==0){
//    cout << "after isotrop" << endl <<(*stress);
//  }
  
  AddtoCmatHolzapfelProduct((*cmat),Cinv,2*(mue-lambda*lJ));
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*cmat)(i,j) += lambda * Cinv(i) * Cinv(j); // add lambda Cinv x Cinv
    }
  }
  // end of isotropic part ****************************************************
  
  // 'non-standard' invariants representing stretch^2 in n0_i direction
  vector<double> I(dim);
  for (int i=0; i<dim; ++i) I[i] = C(i);
  // trial S to compute eigenvalues *******************************************
//  for(int j=0; j<8; ++j){
//    for (int i=0; i<3; ++i){
//      cout << li_->at(j)[i] << ",";
//    }
//    cout << "; gp="<<j<<endl;
//  }
  Epetra_SerialDenseMatrix Strial(dim,dim);
  vector<double> lisq(dim);
  for (int i = 0; i < dim; ++i) lisq[i] = li_->at(gp)[i] * li_->at(gp)[i]; 
  
  // current chain length
  double r = sqrt(I[0]*lisq[0] + I[1]*lisq[1] + I[2]*lisq[2]);
  
  double s_chn = chn_stiffact*(4.0/L + 1.0/(r*(1.0-r/L)*(1.0-r/L)) - 1.0/r);
  
  for (int i = 0; i < dim; ++i) Strial(i,i) = (*stress)(i) + lisq[i]*s_chn + lisq[i]/I[i]*4*stressfree;
  Strial(0,1) = (*stress)(3); Strial(1,0) = (*stress)(3);
  Strial(1,2) = (*stress)(4); Strial(2,1) = (*stress)(4);
  Strial(0,2) = (*stress)(5); Strial(2,0) = (*stress)(5);
  
  Epetra_SerialDenseVector eig_sp(dim);  // lambda^(sigma+)
  if (Strial.NormOne() > 1.0E-13)
//    if (gp==0) cout << Strial;
    LINALG::SymmetricEigenValues(Strial,eig_sp);
//  if ((params.get("total time",-1.0)>=0) && (eleId==0) && (gp==0)){
//    cout << "rem_time: " << rem_time; //params.get("total time",-1.0);
//    cout << ",eigvs: " << eig_sp(0) << "," << eig_sp(1) << "," << eig_sp(2) << endl;
//  }
  vector<double>eig_sort;
  for (int i = 0; i < dim; ++i) eig_sort.push_back(eig_sp(i));
  std::sort(eig_sort.begin(),eig_sort.end());
  std::reverse(eig_sort.begin(),eig_sort.end());
  
  
  for (int i = 0; i < dim; ++i) {
//    if (eig_sort[i] > 1.0E-13) eig_sp(i) = 1.0;
//    else eig_sp(i) = 0.0;
    if (eig_sp(i) > 1.0E-12) eig_sp(i) = 1.0;
    else eig_sp(i) = 0.0;
  }
  // end of trial S to compute eigenvalues ************************************
  
  // update cell dimensions (remodeling!)
  if ((Strial.NormOne() > 1.0E-13) && rem_time > 0.0){
    for (int i = 0; i < dim; ++i){
//      if (gp==0) cout << li_->at(gp)[i] << "->";
      li_->at(gp)[i] = (eig_sp(i) - li0_->at(gp)[i]/r0)*(1-decay)*r0 + li0_->at(gp)[i];
//      if (gp==0) cout << li_->at(gp)[i] << endl;
      lisq[i] = li_->at(gp)[i] * li_->at(gp)[i]; 
    }
  }
  
//  if (gp==0){
//    cout << "isotrop" << endl <<(*stress);
//    for (int i = 0; i < dim; ++i) cout << lisq[i] << "," << I[i] << ";";
//    cout << endl;
//  }
  
  // evaluate chain stress and add to isotropic stress
  r = sqrt(I[0]*lisq[0] + I[1]*lisq[1] + I[2]*lisq[2]);
  s_chn = chn_stiffact*(4.0/L + 1.0/(r*(1.0-r/L)*(1.0-r/L)) - 1.0/r);
  for (int i = 0; i < dim; ++i){
    (*stress)(i) += lisq[i]*s_chn + lisq[i]/I[i]*4*stressfree;
//    if (gp==0){
//      cout << "lisq=" << lisq[i] << " I=" << I[i] << "; " << lisq[i]*s_chn << " zu " << lisq[i]/I[i]*4*stressfree << endl;
//    }
  }
  
//  if (gp==0){
//    cout << "anisotrop" << endl << (*glstrain) <<(*stress) << endl;
//  }
  
  //evaluate chain tangent and add to isotropic cmat
//  if (gp==0){
//    cout << "isotrop" << endl <<(*cmat);
//  }
  double c_chn = chn_stiffact/(r*r*r) * (1.0 - 1.0/((1.0-r/L)*(1.0-r/L)) + 2.0*r/(L*(1.0-r/L)*(1.0-r/L)*(1.0-r/L)) );
  // chain part only affects upper left block of cmat (no shear)
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) (*cmat)(i,j) += c_chn * lisq[i] * lisq[j];
    (*cmat)(i,i) +=  - 8.0*stressfree*lisq[i]/(I[i]*I[i]);
  }
//  if (gp==0){
//    cout << "anisotrop" << endl <<(*cmat) << endl;
//  }
  
  return;
}


#endif
