/*!----------------------------------------------------------------------
\brief
Viscohyperelastic material model containing the following parts:
IsoNeohooke + VolSussmannBathe with Generalized Maxwell (just on isochoric part).
The input line should read
MAT 1 MAT_VISCONEOHOOKE YOUNGS_SLOW 1.0 POISSON 0.499 DENS 0.1 YOUNGS_FAST 100.0 RELAX 10.0 THETA
0.5

\level 2

<pre>
\maintainer Fabian Braeu
</pre>
*----------------------------------------------------------------------*/


#include <vector>
#include "visconeohooke.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::ViscoNeoHooke::ViscoNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_slow_(matdata->GetDouble("YOUNGS_SLOW")),
      poisson_(matdata->GetDouble("POISSON")),
      density_(matdata->GetDouble("DENS")),
      youngs_fast_(matdata->GetDouble("YOUNGS_FAST")),
      relax_(matdata->GetDouble("RELAX")),
      theta_(matdata->GetDouble("THETA"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ViscoNeoHooke::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ViscoNeoHooke(this));
}


MAT::ViscoNeoHookeType MAT::ViscoNeoHookeType::instance_;


DRT::ParObject* MAT::ViscoNeoHookeType::Create(const std::vector<char>& data)
{
  MAT::ViscoNeoHooke* visco = new MAT::ViscoNeoHooke();
  visco->Unpack(data);
  return visco;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoNeoHooke::ViscoNeoHooke() : params_(NULL)
{
  isinit_ = false;
  histstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  artstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  histstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  artstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoNeoHooke::ViscoNeoHooke(MAT::PAR::ViscoNeoHooke* params) : params_(params) {}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  //  pack history data
  int histsize;
  if (!Initialized())
  {
    histsize = 0;
  }
  else
  {
    histsize = histstresslast_->size();
  }
  AddtoPack(data, 2 * histsize);  // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    AddtoPack(data, histstresslast_->at(var));
    AddtoPack(data, artstresslast_->at(var));
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ViscoNeoHooke*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // history data
  int twicehistsize;
  ExtractfromPack(position, data, twicehistsize);

  if (twicehistsize == 0) isinit_ = false;

  histstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  artstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  histstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  artstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  for (int var = 0; var < twicehistsize; var += 2)
  {
    LINALG::Matrix<NUM_STRESS_3D, 1> tmp(true);
    histstresscurr_->push_back(tmp);
    artstresscurr_->push_back(tmp);
    ExtractfromPack(position, data, tmp);
    histstresslast_->push_back(tmp);
    ExtractfromPack(position, data, tmp);
    artstresslast_->push_back(tmp);
  }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  Initialise/allocate internal stress variables (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  histstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  artstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  histstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  artstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  const LINALG::Matrix<NUM_STRESS_3D, 1> emptyvec(true);
  histstresscurr_->resize(numgp);
  histstresslast_->resize(numgp);
  artstresscurr_->resize(numgp);
  artstresslast_->resize(numgp);
  for (int j = 0; j < numgp; ++j)
  {
    histstresscurr_->at(j) = emptyvec;
    histstresslast_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
    artstresslast_->at(j) = emptyvec;
  }

  const double E_s = params_->youngs_slow_;
  double E_f = params_->youngs_fast_;
  double tau = params_->relax_;

  if (E_f < E_s) dserror("Wrong ratio between fast and slow Young's modulus");
  if (tau <= 0.0) dserror("Relaxation time tau has to be positive!");
  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
 |  Update internal stress variables              (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Update()
{
  histstresslast_ = histstresscurr_;
  artstresslast_ = artstresscurr_;
  const LINALG::Matrix<NUM_STRESS_3D, 1> emptyvec(true);
  histstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  artstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  const int numgp = histstresslast_->size();
  histstresscurr_->resize(numgp);
  artstresscurr_->resize(numgp);
  for (int j = 0; j < numgp; ++j)
  {
    histstresscurr_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoNeoHooke::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<NUM_STRESS_3D, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<NUM_STRESS_3D, 1>* stress, LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,
    const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  // get material parameters
  const double E_s = params_->youngs_slow_;
  const double nue = params_->poisson_;
  double E_f = params_->youngs_fast_;
  double tau = params_->relax_;
  const double theta = params_->theta_;

  // get time algorithmic parameters
  // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
  // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
  double dt = params.get<double>("delta time");

  // compute algorithmic relaxation time
  double tau1 = tau;
  // check for meaningful values
  if (E_f > E_s)
  {
    tau1 = tau * E_s / (E_f - E_s);
  }

  // initialize scalars
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
  tau = tau1;
  // evaluate "alpha" factors which distribute stress or stiffness between parallel springs
  // sum_0^i alpha_j = 1
  alpha0 = E_s / E_f;
  alpha1 = 1.0 - alpha0;

  // evaluate Lame constants, bulk modulus
  lambda = nue * E_f / ((1.0 + nue) * (1.0 - 2.0 * nue));
  mue = E_f / (2.0 * (1.0 + nue));
  kappa = lambda + 2.0 / 3.0 * mue;

  // evaluate scalars to compute
  // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + S^(n+1) - S^n]
  artscalar1 = (tau - dt + theta * dt) / tau;
  artscalar2 = tau / (tau + theta * dt);

  // factor to calculate visco stiffness matrix from elastic stiffness matrix
  scalarvisco = alpha0 + alpha1 * tau / (tau + theta * dt);

#else
  // in this case we have a parallel layout of a spring and a dashpot,
  // so no stress distribution between parallel springs
  alpha0 = 1.;
  alpha1 = 1.;

  // do we have to propagate in time?
  if (dt > 0.0)
  {
    // evaluate scalars to compute
    // Q^(n+1) = tau/(theta*dt) [(-dt+theta*dt)/tau Q + S^(n+1) - S^n]
    artscalar1 = (-dt + theta * dt) / tau;
    artscalar2 = tau / (theta * dt);

    // factor to calculate visco stiffness matrix from elastic stiffness matrix
    scalarvisco = 1.0 + tau / (theta * dt);
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
  LINALG::Matrix<NUM_STRESS_3D, 1> Id;
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  for (int i = 3; i < 6; i++) Id(i) = 0.0;
  LINALG::Matrix<NUM_STRESS_3D, 1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  // invariants
  const double I1 = C(0) + C(1) + C(2);  // 1st invariant, trace
  const double I3 = C(0) * C(1) * C(2) + 0.25 * C(3) * C(4) * C(5) - 0.25 * C(1) * C(5) * C(5) -
                    0.25 * C(2) * C(3) * C(3) -
                    0.25 * C(0) * C(4) * C(4);  // 3rd invariant, determinant
  const double J = sqrt(I3);
  const double I3invcubroot = std::pow(I3, -1.0 / 3.0);

  // invert C
  LINALG::Matrix<NUM_STRESS_3D, 1> Cinv;
  Cinv(0) = C(1) * C(2) - 0.25 * C(4) * C(4);
  Cinv(1) = C(0) * C(2) - 0.25 * C(5) * C(5);
  Cinv(2) = C(0) * C(1) - 0.25 * C(3) * C(3);
  Cinv(3) = 0.25 * C(5) * C(4) - 0.5 * C(3) * C(2);
  Cinv(4) = 0.25 * C(3) * C(5) - 0.5 * C(0) * C(4);
  Cinv(5) = 0.25 * C(3) * C(4) - 0.5 * C(5) * C(1);
  Cinv.Scale(1.0 / I3);

  // elastic part: NeoHooke  ************************************************
  // NeoHooke with penalty W = W^dev(C) + U(J)
  // W = 1/2 mue (^I1-3) + 1/2 kappa (J-1)^2

  // Split into volumetric and deviatoric parts. Viscosity affects only deviatoric part
  // Volumetric part of PK2 stress
  LINALG::Matrix<NUM_STRESS_3D, 1> SVol(Cinv);
  SVol.Scale(kappa * (J - 1.0) * J);
  *stress += SVol;

  // Deviatoric elastic part (2 d W^dev/d C)
  LINALG::Matrix<NUM_STRESS_3D, 1> SDevEla(Cinv);
  SDevEla.Scale(-1.0 / 3.0 * I1);
  SDevEla += Id;
  SDevEla.Scale(mue * I3invcubroot);  // mue*I3^(-1/3) (Id-1/3*I1*Cinv)

  // visco part
  // read history
  LINALG::Matrix<NUM_STRESS_3D, 1> S_n(histstresslast_->at(gp));
  S_n.Scale(-1.0);
  LINALG::Matrix<NUM_STRESS_3D, 1> Q_n(artstresslast_->at(gp));

  // artificial visco stresses
  LINALG::Matrix<NUM_STRESS_3D, 1> Q(Q_n);
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
  double scalar1 = 2.0 * kappa * J * J - kappa * J;
  double scalar2 = -2.0 * kappa * J * J + 2.0 * kappa * J;
  double scalar3 = 2.0 / 3.0 * mue * I3invcubroot * I1;
  double scalar4 = 2.0 / 3.0 * mue * I3invcubroot;

  // add volumetric elastic part 1
  // add scalar2 Cinv o Cinv (see Holzapfel p. 254)
  AddtoCmatHolzapfelProduct((*cmat), Cinv, scalar2);

  // add visco-elastic deviatoric part 1
  AddtoCmatHolzapfelProduct(*cmat, Cinv, scalarvisco * scalar3);

  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 6; ++j)
    {
      // add volumetric elastic part 2
      (*cmat)(i, j) += scalar1 * Cinv(i) * Cinv(j)  // add scalar Cinv x Cinv
                                                    // add visco-elastic deviatoric part 2
                       + scalarvisco * (-scalar4) * Id(i) * Cinv(j)        // add scalar Id x Cinv
                       + scalarvisco * (-scalar4) * Id(j) * Cinv(i)        // add scalar Cinv x Id
                       + scalarvisco * (scalar3)*Cinv(i) * Cinv(j) / 3.0;  // add scalar Cinv x Cinv
    }
  }
  return;
}
