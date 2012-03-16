/*!----------------------------------------------------------------------
\file viscopolyC20compressible.cpp
\brief

<pre>
Maintainer: Aline Brunon
            brunon@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include "viscopolyC20compressible.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::ViscoPolyC20Compressible::ViscoPolyC20Compressible(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  //kappa_(matdata->GetDouble("VOLUMETRIC_INF_OGDEN")),
  //tau_volum_(matdata->GetDouble("RELAX_TIME_VOLUME")),
  //beta_volum_(matdata->GetDouble("RELAX_CONST_VOLUME")),
  //density_(matdata->GetDouble("DENS")),
  //C10_(matdata->GetDouble("CONST_NEOHOOKE")),
  //tau_isochoric1_(matdata->GetDouble("RELAX_TIME_ISO")),
  //beta_isochoric1_(matdata->GetDouble("RELAX_CONST_ISO"))
  youngs_slow_(matdata->GetDouble("YOUNGS_SLOW")),
  poisson_(matdata->GetDouble("POISSON")),
  density_(matdata->GetDouble("DENS")),
  youngs_fast_(matdata->GetDouble("YOUNGS_FAST")),
  relax_(matdata->GetDouble("RELAX")),
  theta_(matdata->GetDouble("THETA"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ViscoPolyC20Compressible::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ViscoPolyC20Compressible(this));
}


MAT::ViscoPolyC20CompressibleType MAT::ViscoPolyC20CompressibleType::instance_;


DRT::ParObject* MAT::ViscoPolyC20CompressibleType::Create( const std::vector<char> & data )
{
  MAT::ViscoPolyC20Compressible* visco = new MAT::ViscoPolyC20Compressible();
  visco->Unpack(data);
  return visco;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoPolyC20Compressible::ViscoPolyC20Compressible()
  : params_(NULL)
{
  isinit_=false;
  histstresscurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstresscurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstresslast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstresslast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoPolyC20Compressible::ViscoPolyC20Compressible(MAT::PAR::ViscoPolyC20Compressible* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoPolyC20Compressible::Pack(DRT::PackBuffer& data) const
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

  //  pack history data
  int histsize;
  if (!Initialized())
  {
    histsize=0;
  }
  else
  {
    histsize = histstresslast_->size();
  }
  AddtoPack(data,histsize);  // Length of history vector(s)
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
void MAT::ViscoPolyC20Compressible::Unpack(const vector<char>& data)
{
  isinit_=true;
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
        params_ = static_cast<MAT::PAR::ViscoPolyC20Compressible*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  int histsize;
  ExtractfromPack(position,data,histsize);

  if (histsize == 0) isinit_=false;

  histstresscurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstresscurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstresslast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstresslast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  for (int var=0; var<histsize; var+=1)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> tmp(true);
    histstresscurr_->push_back(tmp);
    artstresscurr_->push_back(tmp);
    histstressvolcurr_->push_back(tmp);
    artstressvolcurr_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstresslast_->push_back(tmp);
    histstressvollast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstresslast_->push_back(tmp);
    artstressvollast_->push_back(tmp);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Initialise/allocate internal stress variables (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoPolyC20Compressible::Setup(const int numgp)
{
  histstresscurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstresscurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstresslast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstresslast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  histstresscurr_->resize(numgp);
  histstresslast_->resize(numgp);
  artstresscurr_->resize(numgp);
  artstresslast_->resize(numgp);
  histstressvolcurr_->resize(numgp);
  histstressvollast_->resize(numgp);
  artstressvolcurr_->resize(numgp);
  artstressvollast_->resize(numgp);
  for (int j=0; j<numgp; ++j)
  {
    histstresscurr_->at(j) = emptyvec;
    histstresslast_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
    artstresslast_->at(j) = emptyvec;
    histstressvolcurr_->at(j) = emptyvec;
    histstressvollast_->at(j) = emptyvec;
    artstressvolcurr_->at(j) = emptyvec;
    artstressvollast_->at(j) = emptyvec;
  }

  const double E_s  = params_->youngs_slow_;
  double E_f  = params_->youngs_fast_;
  double tau  = params_->relax_;

  if (E_f < E_s) dserror("Wrong ratio between fast and slow Young's modulus");
  if (tau<=0.0) dserror("Relaxation time tau has to be positive!");
  isinit_=true;
  return ;

}

/*----------------------------------------------------------------------*
 |  Update internal stress variables              (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoPolyC20Compressible::Update()
{
  histstresslast_=histstresscurr_;
  artstresslast_=artstresscurr_;
  histstressvollast_=histstressvolcurr_;
  artstressvollast_=artstressvolcurr_;
  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  histstresscurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstresscurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  const int numgp=histstresslast_->size();
  histstresscurr_->resize(numgp);
  artstresscurr_->resize(numgp);
  histstressvolcurr_->resize(numgp);
  artstressvolcurr_->resize(numgp);
  for (int j=0; j<numgp; ++j)
  {
    histstresscurr_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
    histstressvolcurr_->at(j) = emptyvec;
    artstressvolcurr_->at(j) = emptyvec;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Reset internal stress variables               (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoPolyC20Compressible::Reset()
{
  // do nothing,
  // because #histstresscurr_ and #artstresscurr_ are recomputed anyway at every iteration
  // based upon #histstresslast_ and #artstresslast_ untouched within time step
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoPolyC20Compressible::Evaluate(const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
                                  const int gp,
                                  Teuchos::ParameterList& params,
                                  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
                                  LINALG::Matrix<NUM_STRESS_3D,1>* stress)

{
  // get material parameters
  const double E_s  = params_->youngs_slow_;
  const double nue  = params_->poisson_;
  double E_f  = params_->youngs_fast_;
  double tau  = params_->relax_;
  const double theta= params_->theta_;

  // get time algorithmic parameters
  // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
  // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
  double dt = params.get<double>("delta time");

  // compute algorithmic relaxation time
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

  // do we have to propagate in time?
  if (dt >0.0)
  {
    // evaluate scalars to compute
    // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + S^(n+1) - S^n]
    artscalar1=(tau-dt+theta*dt)/tau;
    artscalar2=tau/(tau+theta*dt);

    // factor to calculate visco stiffness matrix from elastic stiffness matrix
    scalarvisco = 1.0+tau/(theta*dt);//1+exp(-dt/(2*tau));//1.0+tau/(theta*dt);
  }
  else
  {
    // in case we do not want to propagate in time, Q^{n+1} = Q^{n}
    artscalar1 = 1.0;
    artscalar2 = 1.0;
    // factor to calculate visco stiffness matrix from elastic stiffness matrix
    scalarvisco = 2.0;
  }

#endif

  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id;
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  for (int i =3; i<6;i++) Id(i)=0.0;
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
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
  const double I3invcubroot = pow(I3,-1.0/3.0); // J^{-2/3}

  // Modified first invariant : ^I1 = J^{-2/3}*I1
  const double I1b = I1*I3invcubroot ;

  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv;
  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);

  // elastic part: polynome reduit ordre 2  ************************************************
  // NeoHooke with penalty type Ogden W = W^dev(C) + U_ogd(J)
  // W = 1/2 mue (^I1-3)^2 + (kappa/gamma^2)(gamma.lnJ+J^(-gamma)-1) with gamma = -2

  // Split into volumetric and deviatoric parts. Viscosity affects only deviatoric part
  // Volumetric part of PK2 stress
  LINALG::Matrix<NUM_STRESS_3D,1> SVol(Cinv);
  SVol.Scale(-kappa*(1.0-(J*J))/2.0);//(kappa*(J-1.0)*J);
  *stress+=SVol;

  // Viscous volumetric part
  // read history
  LINALG::Matrix<NUM_STRESS_3D,1> Svol_n (histstressvollast_->at(gp));
  Svol_n.Scale(-1.0);
  LINALG::Matrix<NUM_STRESS_3D,1> Qvol_n (artstressvollast_->at(gp));
  LINALG::Matrix<NUM_STRESS_3D,1> Qvol(Qvol_n);
  Qvol.Scale(artscalar1);
  Qvol += SVol;
  Qvol += Svol_n;
  Qvol.Scale(artscalar2);
  // update history
  histstressvolcurr_->at(gp) = SVol;
  artstressvolcurr_->at(gp) = Qvol;

  // Deviatoric elastic part (2 d W^dev/d C)
  LINALG::Matrix<NUM_STRESS_3D,1> SDevEla(Cinv);
  SDevEla.Scale(-1.0/3.0*I1);
  SDevEla+=Id;
  SDevEla.Scale(2.0*mue*I3invcubroot*(I1b-3.0));  // 2*mue*J^{-2/3}*(^I1-3) (Id-1/3*I1*Cinv)

  // visco part
  // read history
  LINALG::Matrix<NUM_STRESS_3D,1> S_n (histstresslast_->at(gp));
  S_n.Scale(-1.0);
  LINALG::Matrix<NUM_STRESS_3D,1> Q_n (artstresslast_->at(gp));

  // artificial visco stresses
  LINALG::Matrix<NUM_STRESS_3D,1> Q(Q_n);
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
  double scalar1 = kappa*J*J;//2.0*kappa*J*J - kappa*J;
  double scalar2 = kappa*(1-(J*J));//-2.0*kappa*J*J + 2.0*kappa*J;

  double scalara = 4.0*I3invcubroot*mue*(I1b-3.0)*I1/3.0 ;
  double scalarb = 2.0*mue*(I1b-3.0)*I3invcubroot;
  double scalarc = -2.0*I1/3.0;

  // add volumetric elastic part 1
  // add scalar2 Cinv o Cinv (see Holzapfel p. 254)
  AddtoCmatHolzapfelProduct((*cmat),Cinv,scalarvisco*scalar2);

  // add visco-elastic deviatoric part 1
  AddtoCmatHolzapfelProduct(*cmat,Cinv,scalarvisco*scalara);

  for (int i=0; i<6; ++i)
  {
     for (int j=0; j<6; ++j)
     {
       // add volumetric elastic part 2
       (*cmat)(i,j) += scalarvisco*scalar1 * Cinv(i) * Cinv(j) // add scalar Cinv x Cinv
       // add visco-elastic deviatoric part 2
           - 2.0*scalarvisco*scalarb*Id(i)*Cinv(j)/3.0 // add scalar Id x Cinv
           - 2.0*scalarvisco*scalarb*Id(j)*Cinv(i)/3.0 // add scalar Cinv x Id
           - scalarvisco*(scalara+2.0*scalarb*scalarc)*Cinv(i)*Cinv(j)/3.0;// add scalar Cinv x Cinv
     }
  }
  return;
}

#endif
