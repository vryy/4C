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


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "contchainnetw.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::ContChainNetw::ContChainNetw(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  lambda_(matdata->GetDouble("LAMBDA")),
  mue_(matdata->GetDouble("MUE")),
  density_(matdata->GetDouble("DENS")),
  nchain_(matdata->GetDouble("NCHAIN")),
  abstemp_(matdata->GetDouble("ABSTEMP")),
  contl_l_(matdata->GetDouble("CONTL_L")),
  persl_a_(matdata->GetDouble("PERSL_A")),
  r0_(matdata->GetDouble("R0")),
  relax_(matdata->GetDouble("RELAX")),
  initran_(matdata->GetInt("INITRAN")),
  rembegt_(matdata->GetDouble("REMBEGT")),
  updrate_(matdata->GetInt("UPDRATE")),
  difftol_(matdata->GetDouble("DIFFTOL"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::ContChainNetw::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ContChainNetw(this));
}


MAT::ContChainNetwType MAT::ContChainNetwType::instance_;


DRT::ParObject* MAT::ContChainNetwType::Create( const std::vector<char> & data )
{
  MAT::ContChainNetw* chain = new MAT::ContChainNetw();
  chain->Unpack(data);
  return chain;
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         06/08|
 *----------------------------------------------------------------------*/
MAT::ContChainNetw::ContChainNetw()
  : params_(NULL)
{
  isinit_=false;
  //mytime_=0.0;
  li_ = Teuchos::rcp(new vector<vector<double> >);
  li0_ = Teuchos::rcp(new vector<vector<double> >);
  lambda_ = Teuchos::rcp(new vector<vector<double> >);
  ni_ = Teuchos::rcp(new vector<LINALG::Matrix<3,3> >);
  stresses_ = Teuchos::rcp(new vector<LINALG::Matrix<3,3> >);
  mytime_ = Teuchos::rcp(new vector<double>);
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          06/08|
 *----------------------------------------------------------------------*/
MAT::ContChainNetw::ContChainNetw(MAT::PAR::ContChainNetw* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
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
    AddtoPack(data,li0_->at(var));
    AddtoPack(data,ni_->at(var));
    AddtoPack(data,stresses_->at(var));
    AddtoPack(data,mytime_->at(var));
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ContChainNetw*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  isinit_ = true;
  int histsize;
  ExtractfromPack(position,data,histsize);

  if (histsize == 0) isinit_=false;
  li_ = Teuchos::rcp(new vector<vector<double> >);
  li0_ = Teuchos::rcp(new vector<vector<double> >);
  ni_ = Teuchos::rcp(new vector<LINALG::Matrix<3,3> >);
  stresses_ = Teuchos::rcp(new vector<LINALG::Matrix<3,3> >);
  mytime_ = Teuchos::rcp(new vector<double>);
  for (int var = 0; var < histsize; ++var) {
    vector<double> li;
    vector<double> li0;
    LINALG::Matrix<3,3> tmp;
    ExtractfromPack(position,data,li);
    ExtractfromPack(position,data,li0);
    ExtractfromPack(position,data,tmp);
    li_->push_back(li);
    li0_->push_back(li);
    ni_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    stresses_->push_back(tmp);
    double mytime;
    ExtractfromPack(position,data,mytime);
    mytime_->push_back(mytime);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Initialize(const int numgp, const int eleid)
{
  const double r0 = params_->r0_;
  const double isotropy  = 1/sqrt(3.0) * r0;
  srand ( time(NULL) + 5 + eleid*numgp );

  li0_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  li_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  lambda_ = Teuchos::rcp(new vector<vector<double> > (numgp));
  ni_ = Teuchos::rcp(new vector<LINALG::Matrix<3,3> >);
  stresses_ = Teuchos::rcp(new vector<LINALG::Matrix<3,3> >);
  mytime_ = Teuchos::rcp(new vector<double>);
  // initial basis is identity
  LINALG::Matrix<3,3> id(true);
  for (int i=0; i<3; ++i) id(i,i) = 1.0;
  LINALG::Matrix<3,3> initstress(true);

  vector<double> randominit(3);
  double rescale = 0.0;
  for (int i = 0; i < 3; ++i) {
    randominit[i] = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
    rescale += randominit[i]*randominit[i];
  }
  rescale = r0 / sqrt(rescale);

  // initialize cell dimensions
  for(int gp=0; gp<numgp; ++gp){
    li0_->at(gp).resize(3);
    li_->at(gp).resize(3);
    lambda_->at(gp).resize(3);
    if (params_->initran_ == 1){
      // random init, however sum of li^2 has to be r0^2
      li0_->at(gp)[0] = randominit[0]*rescale;
      li0_->at(gp)[1] = randominit[1]*rescale;
      li0_->at(gp)[2] = randominit[2]*rescale;
    } else if (params_->initran_ == 0){
      // pseudo-isotropic init
      li0_->at(gp)[0] = 1.0*isotropy;
      li0_->at(gp)[1] = 1.0*isotropy;  // 2nd dir off! check update!
      li0_->at(gp)[2] = 1.0*isotropy;
//      li0_->at(gp)[0] = 0.95;
//      li0_->at(gp)[1] = 0.95;  // 2nd dir off! check update!
//      li0_->at(gp)[2] = 0.05;
    } else dserror("Unknown remodeling initialization");
    for (int i = 0; i < 3; ++i){
      li_->at(gp)[i] = li0_->at(gp)[i];
      lambda_->at(gp)[i] = 0.0;
    }
    ni_->push_back(id);
    stresses_->push_back(initstress);
    mytime_->push_back(params_->rembegt_);
  }

  //mytime_ = 1.0;  // carefull!
  isinit_ = true;

  return ;

}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         06/08|
 *----------------------------------------------------------------------*

*/

void MAT::ContChainNetw::Evaluate(const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
                                  const int gp,
                                  Teuchos::ParameterList& params,
                                  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
                                  LINALG::Matrix<NUM_STRESS_3D,1> * stress,
                                  int eleId)
{
  // bulk (isotropic) NeoHooke material parameters (Lame constants)
  const double lambda = params_->lambda_;
  const double mue = params_->mue_;
  // chain network unit cell material parameters
  const double nchain = params_->nchain_; // chain density ~= cell stiffness
  const double abstemp = params_->abstemp_; // absolute temperature (K)
  const double L = params_->contl_l_;  // chain contour length
  const double A = params_->persl_a_;  // chain persistence length
  double r0 = params_->r0_;      // initial chain length
  const double boltzmann = 1.3806503E-23;
  const int dim = 3;

  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < dim; i++) Id(i) = 1.0;
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;
  //C.Update(2.0,Id,1.0);

  // we need the 3 by 3 matrix as well later on -> needs improvement
  LINALG::Matrix<dim,dim> CG(false);
  CG(0,0) = C(0); CG(1,1) = C(1); CG(2,2) = C(2);
  CG(0,1) = 0.5*C(3); CG(1,0) = 0.5*C(3);
  CG(1,2) = 0.5*C(4); CG(2,1) = 0.5*C(4);
  CG(0,2) = 0.5*C(5); CG(2,0) = 0.5*C(5);

  // invariants
  const double IC3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  const double J = sqrt(IC3);
  const double lJ = log(J);

  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv(false);

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
//  LINALG::Matrix<NUM_STRESS_3D,1> Siso1(Cinv);
//  Siso1.Scale(lambda*lJ-mue);
//  *stress += Siso1;
//  Siso1 = Id;
//  Siso1.Scale(mue);
//  *stress += Siso1;

  (*stress).Update((lambda*lJ-mue),Cinv,0.0);
  (*stress).Update(mue,Id,1.0);

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*cmat)(i,j) = lambda * Cinv(i) * Cinv(j); // add lambda Cinv x Cinv
    }
  }
  AddtoCmatHolzapfelProduct((*cmat),Cinv,2*(mue-lambda*lJ));

  // anisotropic part *********************************************************

  // the chain stiffness factor
  const double chn_stiffact = boltzmann * abstemp * nchain / (4.0*A);

  // initial cell dimensions
  vector<double> li0sq(dim);
  for (int i=0; i<dim; ++i) li0sq[i] = li0_->at(gp)[i]*li0_->at(gp)[i];
  r0 = sqrt(li0sq[0] + li0sq[1] + li0sq[2]);
  // scalar to arrive at stressfree reference conf
  double stressfree = - chn_stiffact * ( 1.0/L + 1.0/(4.0*r0*(1.0-r0/L)*(1.0-r0/L)) - 1.0/(4.0*r0) );

  // structural tensors Ni0
  vector<LINALG::Matrix<dim,dim> > Ni = EvaluateStructTensors(gp);
  // 'non-standard' invariants representing stretch^2 in n0_i direction
  vector<double> I = EvaluateInvariants(CG,Ni);
  // current cell dimensions
  vector<double> lisq(dim);
  for (int i = 0; i < dim; ++i) lisq[i] = li_->at(gp)[i] * li_->at(gp)[i];
  double r = sqrt(I[0]*lisq[0] + I[1]*lisq[1] + I[2]*lisq[2]);
  double s_chn = chn_stiffact*(4.0/L + 1.0/(r*(1.0-r/L)*(1.0-r/L)) - 1.0/r);

  // evaluate current stress (including isotropic and anisotropic part)
  LINALG::Matrix<dim,dim> S(false);
  StressVoigt2Mat(stress,S);
  UpdateStress(S,Ni,lisq,I,s_chn,stressfree);

  // do remodeling only if we are at a new time step and based on last stress
  const double kappa = params_->relax_; // relaxation time for remodeling
  const double time = params.get("total time",-1.0);
  const double dt = params.get("delta time",-1.0);
  const double lambda_tol = params_->difftol_; //1.0E-12;
  if ( (kappa >= 0.0)  && (time > mytime_->at(gp)) ){
    double rem_time = time - params_->rembegt_;
    mytime_->at(gp) = time;
    const double decay = min(1.0,exp(-kappa*rem_time));

    // evaluate eigenproblem
    Epetra_SerialDenseVector lambda(dim);

    //needs improvement:
    Epetra_SerialDenseMatrix Phi(3,3);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        Phi(i,j) = stresses_->at(gp)(i,j);
      }
    }

    //Epetra_SerialDenseMatrix Phi = S;
    LINALG::SymmetricEigenProblem(Phi,lambda);
    //Epetra_SerialDenseMatrix strain = MAT::StrainVoigt2Mat(glstrain);
    //LINALG::SymmetricEigenProblem(strain,lambda);

    // initialise rem_toggle which means per default diminuish cell lengths
    vector<double> rem_toggle(3,0.0);

    // decide on remodeling strategy and update cell dimensions and history
    if (params_->updrate_ == 0){
      for (int i=0; i<3; ++i){
        if (lambda(i)/lambda(2) > lambda_tol) rem_toggle[i] = 1.0;
      }
      Update(Phi,lambda,rem_toggle,gp,r0,decay);
    } else if (params_->updrate_ == 1){
      for (int i=0; i<3; ++i){
        if (lambda(i)/lambda(2) > lambda_tol) rem_toggle[i] = 1.0;
      }
      UpdateRate(Phi,lambda,rem_toggle,gp,r0,decay,kappa,dt);
    } else if (params_->updrate_ == 2){  // my new strategy
      for (int i=0; i<3; ++i){
        if (lambda(i)/lambda(2) > lambda_tol) rem_toggle[i] = lambda(i)/lambda(2);
      }
      Update(Phi,lambda,rem_toggle,gp,r0,decay);
    } else if (params_->updrate_ == 3){  // my new strategy
      for (int i=0; i<3; ++i){
        if (lambda(i)/lambda(2) > lambda_tol) rem_toggle[i] = lambda(i)/lambda(2);
      }
      UpdateRate(Phi,lambda,rem_toggle,gp,r0,decay,kappa,dt);
    } else dserror("Unknown Update Flag");

#if DEBUG1
    // update initial chain length
    for (int i = 0; i < 3; ++i) li0_->at(gp)[i] = li_->at(gp)[i];
    for (int i=0; i<dim; ++i) li0sq[i] = li_->at(gp)[i]*li_->at(gp)[i];
    double r0new = sqrt(li0sq[0] + li0sq[1] + li0sq[2]);

    if (abs(r0-r0new) > 1.0E-12) dserror("r0 scaling problem");
#endif


    // reevaluate stress with 'remodeled' parameters
    Ni = EvaluateStructTensors(gp);
    I = EvaluateInvariants(CG,Ni);
    for (int i = 0; i < dim; ++i) lisq[i] = li_->at(gp)[i] * li_->at(gp)[i];
    double r = sqrt(I[0]*lisq[0] + I[1]*lisq[1] + I[2]*lisq[2]);

    stressfree = - chn_stiffact * ( 1.0/L + 1.0/(4.0*r0*(1.0-r0/L)*(1.0-r0/L)) - 1.0/(4.0*r0) );

    double s_chn = chn_stiffact*(4.0/L + 1.0/(r*(1.0-r/L)*(1.0-r/L)) - 1.0/r);
    //S = EvaluateStress(MAT::StressVoigt2Mat(stress),Ni,lisq,I,s_chn,stressfree);
    //cout << S << endl;
    LINALG::Matrix<dim,dim> Snew(false);
    StressVoigt2Mat(stress,Snew);
    UpdateStress(Snew,Ni,lisq,I,s_chn,stressfree);
    S.Update(Snew);
 }

  // remember current stress
  stresses_->at(gp) = S;
  // return stress
  StressMat2Voigt(*stress,S);

  // chain stiffness factor
  double c_chn = chn_stiffact/(r*r*r) * (1.0 - 1.0/((1.0-r/L)*(1.0-r/L)) + 2.0*r/(L*(1.0-r/L)*(1.0-r/L)*(1.0-r/L)) );
#if DEBUG1
  //cout << "glstrain" << endl << (*glstrain);
  //cout << "stress" << endl << (*stress);
  cout << "li0: " << PrintVec(li0_->at(gp)) << "; li: " << PrintVec(li_->at(gp)) << endl;
  //cout << PrintAnisoVects(gp) << endl << PrintStructTens(Ni) << endl;
  cout << PrintAnisoCmat((*cmat),Ni,lisq,I,c_chn,stressfree) << endl;
#endif
  EvaluateCmat((*cmat),Ni,lisq,I,c_chn,stressfree);

#if DEBUG1
  //cout << *cmat;
#endif



  return;
}

LINALG::Matrix<3,3> MAT::ContChainNetw::EvaluateStress(
    const LINALG::Matrix<3,3>& isostress,
    const vector<LINALG::Matrix<3,3> >& Ni,
    const vector<double>& cell_li,
    const vector<double>& cell_Inv,
    const double s_chn_scalar,
    const double stressfree)
{
#if DEBUG
  if ((cell_li.size()!=3)||(cell_Inv.size()!=3)||(Ni.size()!=3)||(isostress.M()!=3)||(isostress.N()!=3))
    dserror("Wrong dimensions in stress eval");
#endif
  LINALG::Matrix<3,3> stress(isostress);
  for (int i_fib=0; i_fib<3; ++i_fib){
    double f1 = cell_li[i_fib] * s_chn_scalar;
    double f2 = 4 * stressfree * cell_li[i_fib] / cell_Inv[i_fib];
    for (int i=0; i<3; ++i){
      for (int j=0; j<3; ++j){
        stress(i,j) += f1 * Ni[i_fib](i,j);
        stress(i,j) += f2 * Ni[i_fib](i,j);
      }
    }
  }
  return stress;
}

void MAT::ContChainNetw::UpdateStress(
    LINALG::Matrix<3,3>& stress,
    const vector<LINALG::Matrix<3,3> >& Ni,
    const vector<double>& cell_li,
    const vector<double>& cell_Inv,
    const double s_chn_scalar,
    const double stressfree)
{
#if DEBUG
  if ((cell_li.size()!=3)||(cell_Inv.size()!=3)||(Ni.size()!=3)||(stress.M()!=3)||(stress.N()!=3))
    dserror("Wrong dimensions in stress eval");
#endif
  for (int i_fib=0; i_fib<3; ++i_fib){
    double f1 = cell_li[i_fib] * s_chn_scalar;
    double f2 = 4 * stressfree * cell_li[i_fib] / cell_Inv[i_fib];
    for (int i=0; i<3; ++i){
      for (int j=0; j<3; ++j){
        stress(i,j) += f1 * Ni[i_fib](i,j);
        stress(i,j) += f2 * Ni[i_fib](i,j);
      }
    }
  }
  return;
}

void MAT::ContChainNetw::EvaluateCmat(LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,
    const vector<LINALG::Matrix<3,3> >& Ni,
    const vector<double>& cell_li,
    const vector<double>& cell_Inv,
    const double c_chn_scalar,
    const double stressfree)
{
  LINALG::Matrix<3,3> sumN0i(true);
  for (int i_fib=0; i_fib<3; ++i_fib)
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) sumN0i(i,j) += cell_li[i_fib] * Ni[i_fib](i,j);

  ElastSymTensorMultiply(cmat,c_chn_scalar,sumN0i,sumN0i,1.0);  // fac sumN0i x sumN0i

  for (int i_fib=0; i_fib<3; ++i_fib){
    double f = -8.0*stressfree*cell_li[i_fib] / (cell_Inv[i_fib] * cell_Inv[i_fib]);
    ElastSymTensorMultiply(cmat,f,Ni[i_fib],Ni[i_fib],1.0);
  }
  return;
}

vector<LINALG::Matrix<3,3> > MAT::ContChainNetw::EvaluateStructTensors(const int gp)
{
  vector<LINALG::Matrix<3,3> > Ni;
  for (int i_fib=0; i_fib<3; ++i_fib){
    LINALG::Matrix<3,3> N0(false);
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        N0(i,j) = (ni_->at(gp)(i,i_fib)) * (ni_->at(gp)(j,i_fib));
    Ni.push_back(N0);
  }
  return Ni;
}

vector<double> MAT::ContChainNetw::EvaluateInvariants(
    const LINALG::Matrix<3,3> & CG,
    const vector<LINALG::Matrix<3,3> >& Ni)
{
  vector<double> Inv(3);
  for (int i_fib=0; i_fib<3; ++i_fib){
    LINALG::Matrix<3,3> CNi0(false);
    CNi0.Multiply(CG,Ni.at(i_fib));
    Inv[i_fib] = CNi0(0,0) + CNi0(1,1) + CNi0(2,2); // trace(C:Ni0)
  }
  return Inv;
}

void MAT::ContChainNetw::Update(const Epetra_SerialDenseMatrix& Phi,
    const Epetra_SerialDenseVector& lambda,
    const vector<double> rem_toggle,
    const int gp, const double r0, const double decay)
{
  double rescale = 0.0;
  for (int i = 0; i < 3; ++i){
    lambda_->at(gp)[i] = lambda(i);
    li_->at(gp)[i] = (rem_toggle[i] - li0_->at(gp)[i]/r0)*(1.0-decay)*r0 + li0_->at(gp)[i];
    rescale += li_->at(gp)[i]*li_->at(gp)[i];
  }
  rescale = r0 / sqrt(rescale);
  for (int i = 0; i < 3; ++i){
    li_->at(gp)[i] = rescale * li_->at(gp)[i];
  }
  //li_->at(gp)[1] = 0.0;  // 2nd dir off!
  // needs improvement:
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ni_->at(gp)(i,j) = Phi(i,j);
    }
  }
}

void MAT::ContChainNetw::UpdateRate(const Epetra_SerialDenseMatrix& Phi,
    const Epetra_SerialDenseVector& lambda,
    const vector<double> rem_toggle,
    const int gp, const double r0, const double decay, const double kappa, const double dt)
{
  double rescale = 0.0;
  for (int i = 0; i < 3; ++i){
    lambda_->at(gp)[i] = lambda(i);
    li_->at(gp)[i] = li_->at(gp)[i] + kappa*(rem_toggle[i] - li0_->at(gp)[i]/r0)*decay*r0 *dt;
    rescale += li_->at(gp)[i]*li_->at(gp)[i];
  }
  rescale = r0 / sqrt(rescale);
  for (int i = 0; i < 3; ++i){
    li_->at(gp)[i] = rescale * li_->at(gp)[i];
  }
  // needs improvement:
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ni_->at(gp)(i,j) = Phi(i,j);
    }
  }
}

std::string MAT::ContChainNetw::PrintStructTens(const vector<Epetra_SerialDenseMatrix>& Ni)
{
  std::stringstream out;
  for (int i=0;i<3;++i){
    for (int i_fib=0; i_fib<3; ++i_fib){
      for (int j=0;j<3;++j){
        out << Ni[i_fib](i,j) << " ";
      }
      out << "| ";
    }
    out << endl;
  }
  return out.str();
}

std::string MAT::ContChainNetw::PrintAnisoVects(const int gp)
{
  std::stringstream out;
  for (int i=0;i<3;++i){
    for (int i_fib=0; i_fib<3; ++i_fib){
      out << ni_->at(gp)(i,i_fib) << " | ";
    }
    out << endl;
  }
  return out.str();
}

std::string MAT::ContChainNetw::PrintAnisoCmat(const LINALG::Matrix<6,6>& cmat,
    const vector<LINALG::Matrix<3,3> >& Ni,
    const vector<double>& cell_li,
    const vector<double>& cell_Inv,
    const double c_chn_scalar,
    const double stressfree)
{
  std::stringstream out;
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> newCmat;
  EvaluateCmat(newCmat,Ni,cell_li,cell_Inv,c_chn_scalar,stressfree);
  for (int i=0;i<3;++i){
    for (int j=0;j<3;++j){
      out << cmat(i,j) << " ";
    }
    out << "| ";
    for (int j=0;j<3;++j){
      out << newCmat(i,j) << " ";
    }
    out << endl;
  }
  return out.str();
}

std::string MAT::ContChainNetw::PrintVec(const vector<double> actvec)
{
  std::stringstream out;
  vector<double>::const_iterator i;
  for (i=actvec.begin(); i<actvec.end(); ++i) {
    out << *i << " ";
  }
  return out.str();
}


/// transform Voigt stress vector to matrix
Epetra_SerialDenseMatrix MAT::StressVoigt2Mat(Epetra_SerialDenseVector* stress)
{
#if DEBUG
  if ( ((*stress).M() != 6) || ((*stress).N()!=1) )
    dserror("Wrong dimensions");
#endif
  Epetra_SerialDenseMatrix mat(3,3);
  for (int i = 0; i < 3; ++i) mat(i,i) = (*stress)(i);
  mat(0,1) = (*stress)(3); mat(1,0) = (*stress)(3);
  mat(1,2) = (*stress)(4); mat(2,1) = (*stress)(4);
  mat(0,2) = (*stress)(5); mat(2,0) = (*stress)(5);
  return mat;
}
/// transform Voigt stress vector to matrix
LINALG::Matrix<3,3> MAT::StressVoigt2Mat(const LINALG::Matrix<6,1>* stress)
{
#if DEBUG
  if ( ((*stress).M() != 6) || ((*stress).N()!=1) )
    dserror("Wrong dimensions");
#endif
  LINALG::Matrix<3,3> mat(false);
  for (int i = 0; i < 3; ++i) mat(i,i) = (*stress)(i);
  mat(0,1) = (*stress)(3); mat(1,0) = (*stress)(3);
  mat(1,2) = (*stress)(4); mat(2,1) = (*stress)(4);
  mat(0,2) = (*stress)(5); mat(2,0) = (*stress)(5);
  return mat;
}
/// transform Voigt stress vector to matrix
void MAT::StressVoigt2Mat(const LINALG::Matrix<6,1>* stress, LINALG::Matrix<3,3>& mat)
{
#if DEBUG
  if ( ((*stress).M() != 6) || ((*stress).N()!=1) )
    dserror("Wrong dimensions");
#endif
  for (int i = 0; i < 3; ++i) mat(i,i) = (*stress)(i);
  mat(0,1) = (*stress)(3); mat(1,0) = (*stress)(3);
  mat(1,2) = (*stress)(4); mat(2,1) = (*stress)(4);
  mat(0,2) = (*stress)(5); mat(2,0) = (*stress)(5);
  return;
}

/// transform Voigt strain vector to matrix respecting factor 2 for shear
Epetra_SerialDenseMatrix MAT::StrainVoigt2Mat(const Epetra_SerialDenseVector* strain)
{
#if DEBUG
  if ( ((*strain).M() != 6) || ((*strain).N()!=1) )
    dserror("Wrong dimensions");
#endif
  Epetra_SerialDenseMatrix mat(3,3);
  for (int i = 0; i < 3; ++i) mat(i,i) = (*strain)(i);
  mat(0,1) = 0.5*(*strain)(3); mat(1,0) = 0.5*(*strain)(3);
  mat(1,2) = 0.5*(*strain)(4); mat(2,1) = 0.5*(*strain)(4);
  mat(0,2) = 0.5*(*strain)(5); mat(2,0) = 0.5*(*strain)(5);
  return mat;
}

/// transform stress matrix to Voigt vector
Epetra_SerialDenseVector MAT::StressMat2Voigt(Epetra_SerialDenseMatrix& stressmat)
{
#if DEBUG
  if ( (stressmat.M() != 3) || (stressmat.N()!=3) )
    dserror("Wrong dimensions");
#endif
  Epetra_SerialDenseVector s(6);
  for (int i = 0; i < 3; ++i) s(i) = stressmat(i,i);
  s(3) = stressmat(0,1); s(4) = stressmat(1,2); s(5) = stressmat(0,2);
  return s;
}
LINALG::Matrix<6,1> MAT::StressMat2Voigt(LINALG::Matrix<3,3>& stressmat)
{
#if DEBUG
  if ( (stressmat.M() != 3) || (stressmat.N()!=3) )
    dserror("Wrong dimensions");
#endif
  LINALG::Matrix<6,1> s(false);
  for (int i = 0; i < 3; ++i) s(i) = stressmat(i,i);
  s(3) = stressmat(0,1); s(4) = stressmat(1,2); s(5) = stressmat(0,2);
  return s;
}
void MAT::StressMat2Voigt(LINALG::Matrix<6,1>& s,const LINALG::Matrix<3,3>& stressmat)
{
#if DEBUG
  if ( (stressmat.M() != 3) || (stressmat.N()!=3) )
    dserror("Wrong dimensions");
#endif
  for (int i = 0; i < 3; ++i) s(i) = stressmat(i,i);
  s(3) = stressmat(0,1); s(4) = stressmat(1,2); s(5) = stressmat(0,2);
  return;
}



/// Debug output to txt-file
void MAT::ChainOutputToTxt(const Teuchos::RCP<DRT::Discretization> dis,
    const double time,
    const int iter)
{
    std::stringstream filename;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename << filebase << "_rem" << ".txt";
    ofstream outfile;
    outfile.open(filename.str().c_str(),ios_base::app);
    //int nele = dis->NumMyColElements();
    int endele = 100; //nele;
    for (int iele=0; iele<endele; iele+=12) //++iele) iele+=10)
    {
      const DRT::Element* actele = dis->lColElement(iele);
      RefCountPtr<MAT::Material> mat = actele->Material();
      if (mat->MaterialType() != INPAR::MAT::m_contchainnetw) return;
      MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
      //int ngp = chain->Getni()->size();
      int endgp = 1; //ngp;
      for (int gp = 0; gp < endgp; ++gp){
        vector<double> li = chain->Getli()->at(gp);
        vector<double> lamb = chain->Getlambdas()->at(gp);
        LINALG::Matrix<3,3> ni0 = chain->Getni()->at(gp);

        // time
        outfile << time << ",";
        // iter
        outfile << iter << ",";
        // eleId
        outfile << iele << ",";
        // gp
        outfile << gp << ",";
        // cell dimensions
        for (int i=0;i<3;++i) outfile << li[i] << ",";
        // eigenvalues
        for (int i=0;i<3;++i) outfile << lamb[i] << ",";
        // eigenvectors/cell basis
        for (int i=0;i<3;++i)
          for (int j=0;j<3;++j)
            outfile << ni0(j,i) << ",";
        // end
        outfile << endl;
      }
    }
    outfile.close();
    return;
}


/// Debug output to gmsh-file
void MAT::ChainOutputToGmsh(const Teuchos::RCP<DRT::Discretization> dis,
                                      const double time,
                                      const int iter)
{
  std::stringstream filename;
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  filename << filebase << "_ContChainMat" << std::setw(3) << std::setfill('0') << time << std::setw(2) << std::setfill('0') << iter << ".pos";
  std::ofstream f_system(filename.str().c_str());

  stringstream gmshfilecontent;
  gmshfilecontent << "View \" Time: " << time << " Iter: " << iter << " \" {" << endl;
  for (int iele=0; iele<dis->NumMyColElements(); ++iele)
  {
    const DRT::Element* actele = dis->lColElement(iele);

    // build current configuration
    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
    actele->LocationVector(*dis,lm,lmowner,lmstride);
    RCP<const Epetra_Vector> disp = dis->GetState("displacement");
    vector<double> mydisp(lm.size(),0);
    //DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    const int numnode = actele->NumNode();
    const int numdof = 3;
    LINALG::SerialDenseMatrix xyze(3, numnode);
    for (int inode = 0; inode < numnode; ++inode)
    {
      xyze(0, inode) = actele->Nodes()[inode]->X()[0]+ mydisp[inode*numdof+0];
      xyze(1, inode) = actele->Nodes()[inode]->X()[1]+ mydisp[inode*numdof+1];
      xyze(2, inode) = actele->Nodes()[inode]->X()[2]+ mydisp[inode*numdof+2];
    }
    gmshfilecontent << IO::GMSH::cellWithScalarToString(actele->Shape(),
        1.0, xyze) << endl;

    vector<double> elecenter = MAT::MatPointCoords(actele,mydisp);
    RefCountPtr<MAT::Material> mat = actele->Material();
    MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
    LINALG::Matrix<3,3> ni0 = chain->Getni()->at(0);
    vector<double> lamb0 = chain->Getlambdas()->at(0);
    RCP<vector<vector<double> > > gplis = chain->Getli();
    RCP<vector<vector<double> > > gpli0s = chain->Getli0();
    RCP<vector<LINALG::Matrix<3,3> > > gpnis = chain->Getni();

    vector<double> centerli (3,0.0);
    vector<double> centerli_0 (3,0.0);

//    // material plot at element center
//    const int dim=3;
//    for (int k=0; k<dim; ++k){
//      gmshfilecontent << "VP(" << std::scientific << elecenter[0] << ",";
//      gmshfilecontent << std::scientific << elecenter[1] << ",";
//      gmshfilecontent << std::scientific << elecenter[2] << ")";
//      gmshfilecontent << "{" << std::scientific <<
//      ni0(0,k) * lamb0[k]
//      << "," << ni0(1,k) * lamb0[k] << "," << ni0(2,k) * lamb0[k] << "};" << endl;
//    }

    // material plot at gauss points
    int ngp = chain->Getni()->size();
    for (int gp = 0; gp < ngp; ++gp){
      vector<double> point = MAT::MatPointCoords(actele,mydisp,gp);
      double scalar = 1;
      vector<double> length(3);
      vector<double> gpli =  chain->Getli()->at(gp);
      vector<double> gpli0 = chain->Getli0()->at(gp);
//      length[0] = scalar * gpli[0]/gpli0[0];
//      length[1] = scalar * gpli[1]/gpli0[1];
//      length[2] = scalar * gpli[2]/gpli0[2];
      length[0] = scalar * chain->Getlambdas()->at(gp)[0];
      length[1] = scalar * chain->Getlambdas()->at(gp)[1];
      length[2] = scalar * chain->Getlambdas()->at(gp)[2];
//      for (int i=0; i<3; ++i){
//        cout << gpli[i] << ":" << gpli0[i] << endl;
//      }

      LINALG::Matrix<3,1> loc(&(gplis->at(gp)[0]));
      LINALG::Matrix<3,1> glo(false);
      glo.Multiply(gpnis->at(gp),loc);

      for (int k=0; k<3; ++k){
//        // draw eigenvectors
//        gmshfilecontent << "VP(" << std::scientific << point[0] << ",";
//        gmshfilecontent << std::scientific << point[1] << ",";
//        gmshfilecontent << std::scientific << point[2] << ")";
//        gmshfilecontent << "{" << std::scientific
//        << ((chain->Getni())->at(gp))(0,k)
//        << "," << ((chain->Getni())->at(gp))(1,k)
//        << "," << ((chain->Getni())->at(gp))(2,k) << "};" << endl;

        // draw fiber cell vectors
        LINALG::Matrix<3,1> e(false);
        e(k) = gpli[k];
        glo.Multiply(gpnis->at(gp),e);
        gmshfilecontent << "VP(" << std::scientific << point[0] << ",";
        gmshfilecontent << std::scientific << point[1] << ",";
        gmshfilecontent << std::scientific << point[2] << ")";
        gmshfilecontent << "{" << std::scientific
        <<        glo(0)
        << "," << glo(1)
        << "," << glo(2)
        << "};" << endl;
      }
    }
  }
  gmshfilecontent << "};" << endl;

  f_system << gmshfilecontent.str();
  f_system.close();

  return;
}

/// gmsh-debug: calculate gausspoint coordinates
const vector<double> MAT::MatPointCoords(const DRT::Element* actele,const vector<double>& mydisp, int gp)
{
  // update element geometry
  const int numnode = actele->NumNode();
  const int numdof = 3;
  Epetra_SerialDenseMatrix xrefe(numnode,3);  // material coord. of element
  for (int i=0; i<numnode; ++i){
    xrefe(i,0) = actele->Nodes()[i]->X()[0]+ mydisp[i*numdof+0];
    xrefe(i,1) = actele->Nodes()[i]->X()[1]+ mydisp[i*numdof+1];
    xrefe(i,2) = actele->Nodes()[i]->X()[2]+ mydisp[i*numdof+2];
  }
  const DRT::Element::DiscretizationType distype = actele->Shape();
  Epetra_SerialDenseVector funct(numnode);
  // Element midpoint at r=s=t=0.0
  if (gp==-1) DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  else{
    const DRT::UTILS::GaussRule3D gaussrule_ = DRT::UTILS::intrule_hex_8point;
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule_);
    DRT::UTILS::shape_function_3D(funct, intpoints.qxg[gp][0], intpoints.qxg[gp][1], intpoints.qxg[gp][2], distype);
  }
  Epetra_SerialDenseMatrix point(1,3);
  point.Multiply('T','N',1.0,funct,xrefe,0.0);
  vector<double> coords(3);
  coords[0] = point(0,0);
  coords[1] = point(0,1);
  coords[2] = point(0,2);
  return coords;
}


