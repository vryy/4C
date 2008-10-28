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
#include "visconeohooke.H"
#include "../drt_lib/linalg_serialdensevector.H"

extern struct _MATERIAL *mat;  ///< C-style material struct


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoNeoHooke::ViscoNeoHooke()
  : matdata_(NULL)
{
  isinit_=false;
  histstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  artstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  histstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  artstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
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
  //Pack history data
  int histsize;
  if (!Initialized())
  {
    histsize=0;
  }
  else
  {
    histsize = histstresslast_->size();
  }
  AddtoPack(data,2*histsize);  // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    AddtoPack(data,histstresslast_->at(var));
    AddtoPack(data,artstresslast_->at(var));
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Unpack(const vector<char>& data)
{
  isinit_=true;
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
  int twicehistsize;
  ExtractfromPack(position,data,twicehistsize);

  if (twicehistsize == 0) isinit_=false;

  histstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  artstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  histstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  artstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  for (int var=0; var<twicehistsize; var+=2)
  {
    LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> tmp(true);
    histstresscurr_->push_back(tmp);
    artstresscurr_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstresslast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstresslast_->push_back(tmp);
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Initialise/allocate internal stress variables (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Setup(const int numgp)
{
  histstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  artstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  histstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  artstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  const LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> emptyvec(true);
  histstresscurr_->resize(numgp);
  histstresslast_->resize(numgp);
  artstresscurr_->resize(numgp);
  artstresslast_->resize(numgp);
  for (int j=0; j<numgp; ++j)
  {
    histstresscurr_->at(j) = emptyvec;
    histstresslast_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
    artstresslast_->at(j) = emptyvec;
  }
 
  const double E_s  = matdata_->m.visconeohooke->youngs_slow;
  double E_f  = matdata_->m.visconeohooke->youngs_fast;
  double tau  = matdata_->m.visconeohooke->relax;
 
  if (E_f < E_s) dserror("Wrong ratio between fast and slow Young's modulus");
  if (tau<=0.0) dserror("Relaxation time tau has to be positive!");
  isinit_=true;
  return ;

}

/*----------------------------------------------------------------------*
 |  Update internal stress variables              (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Update()
{
  histstresslast_=histstresscurr_;
  artstresslast_=artstresscurr_;
  const LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> emptyvec(true);
  histstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  artstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  const int numgp=histstresslast_->size();
  histstresscurr_->resize(numgp);
  artstresscurr_->resize(numgp);
  for (int j=0; j<numgp; ++j)
  {
    histstresscurr_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Reset internal stress variables               (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Reset()
{
  // do nothing,
  // because #histstresscurr_ and #artstresscurr_ are recomputed anyway at every iteration
  // based upon #histstresslast_ and #artstresslast_ untouched within time step
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Evaluate(const LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1>* glstrain,
                                  const int gp,
                                  Teuchos::ParameterList& params,
                                  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
                                  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1>* stress)

{
  // get material parameters
  const double E_s  = matdata_->m.visconeohooke->youngs_slow;
  const double nue  = matdata_->m.visconeohooke->poisson;
  double E_f  = matdata_->m.visconeohooke->youngs_fast;
  double tau  = matdata_->m.visconeohooke->relax;
  const double theta= matdata_->m.visconeohooke->theta;

  // get time algorithmic parameters
  double dt = params.get("delta time",-1.0);

  double tau1=tau;
  //check for meaningful values
  if (E_f>E_s)
  {
    tau1=tau*E_s/(E_f-E_s);
  }

  //initialize scalars
  double alpha0;
  double alpha1;
  double lambda;
  double mue;
  double kappa;
  double artscalar1;
  double artscalar2;
  double scalarvisco;

#define GEN_MAXWELL
#ifdef GEN_MAXWELL 
  tau=tau1;
  // evaluate "alpha" factors which distribute stress or stiffness between parallel springs
  // sum_0^i alpha_j = 1
  alpha0 = E_s / E_f;
  alpha1 = 1.0 - alpha0;

  // evaluate Lame constants, bulk modulus
  lambda = nue*E_f / ((1.0+nue)*(1.0-2.0*nue));
  mue = E_f / (2.0*(1.0+nue));
  kappa = lambda + 2.0/3.0 * mue;

  // evaluate scalars to compute
  // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + S^(n+1) - S^n]
  artscalar1=(tau - dt + theta*dt)/tau;
  artscalar2=tau/(tau + theta*dt);

  // factor to calculate visco stiffness matrix from elastic stiffness matrix
  scalarvisco = alpha0+alpha1*tau/(tau+theta*dt);
  
#else
  //in this case we have a parallel layout of a spring and a dashpot,
  //so no stress distribution between parallel springs
  alpha0 = 1.;
  alpha1 = 1.;

  // evaluate Lame constants, bulk modulus
  lambda = nue*E_s / ((1.0+nue)*(1.0-2.0*nue));
  mue = E_s / (2.0*(1.0+nue));
  kappa = lambda + 2.0/3.0 * mue;

  // evaluate sclars to compute
  // Q^(n+1) = tau/(theta*dt) [(-dt+theta*dt)/tau Q + S^(n+1) - S^n]
  artscalar1=(-dt+theta*dt)/tau;
  artscalar2=tau/(theta*dt);

  // factor to calculate visco stiffness matrix from elastic stiffness matrix
  scalarvisco = 1.0+tau/(theta*dt);
  
#endif

  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> Id;
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  for (int i =3; i<6;i++) Id(i)=0.0;
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> C(*glstrain);
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
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> Cinv;
  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);

  // elastic part: NeoHooke  ************************************************
  // NeoHooke with penalty W = W^dev(C) + U(J)
  // W = 1/2 mue (^I1-3) + 1/2 kappa (J-1)^2

  // Split into volumetric and deviatoric parts. Viscosity affects only deviatoric part
  // Volumetric part of PK2 stress
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> SVol(Cinv);
  SVol.Scale(kappa*(J-1.0)*J);
  *stress+=SVol;

  // Deviatoric elastic part (2 d W^dev/d C)
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> SDevEla(Cinv);
  SDevEla.Scale(-1.0/3.0*I1);
  SDevEla+=Id;
  SDevEla.Scale(mue*I3invcubroot);  //mue*I3^(-1/3) (Id-1/3*I1*Cinv)

  // visco part
  // read history
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> S_n (histstresslast_->at(gp));
  S_n.Scale(-1.0);
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> Q_n (artstresslast_->at(gp));

  // artificial visco stresses
  LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> Q(Q_n);
  Q.Scale(artscalar1);
  Q += SDevEla;
  Q += S_n;
  Q.Scale(artscalar2);  // Q^(n+1) = artscalar2* [artscalar1* Q + S^(n+1) - S^n]

  // update history
  histstresscurr_->at(gp) = SDevEla;
  artstresscurr_->at(gp) = Q;

  // add visco PK2 stress, weighted with alphas
  SDevEla.Scale(alpha0);
  *stress += SDevEla;
  Q.Scale(alpha1);
  *stress += Q;
  // elasticity matrix
  double scalar1 = 2.0*kappa*J*J - kappa*J;
  double scalar2 = -2.0*kappa*J*J + 2.0*kappa*J;
  double scalar3 = 2.0/3.0*mue*I3invcubroot*I1;
  double scalar4 = 2.0/3.0*mue*I3invcubroot;

  // add volumetric elastic part 1
  // add scalar2 Cinv o Cinv (see Holzapfel p. 254)
  AddtoCmatHolzapfelProduct((*cmat),Cinv,scalar2);

  // add visco-elastic deviatoric part 1
  AddtoCmatHolzapfelProduct(*cmat,Cinv,scalarvisco*scalar3);

  for (int i=0; i<6; ++i)
  {
     for (int j=0; j<6; ++j)
     {
       // add volumetric elastic part 2
       (*cmat)(i,j) += scalar1 * Cinv(i) * Cinv(j) // add scalar Cinv x Cinv
       // add visco-elastic deviatoric part 2
           + scalarvisco*(-scalar4)*Id(i)*Cinv(j)// add scalar Id x Cinv
           + scalarvisco*(-scalar4)*Id(j)*Cinv(i)// add scalar Cinv x Id
           + scalarvisco*(scalar3)*Cinv(i)*Cinv(j)/3.0;// add scalar Cinv x Cinv
     }
  }
  return;
}

#endif
