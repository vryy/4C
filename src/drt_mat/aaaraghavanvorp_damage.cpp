/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required for aneurysmatic artery wall
damaging following Simo approach explained in Holzapfel (p298).
SEF by Raghavan and Vorp[2000](Isochoric) and Simo&Miehe version of
Ogden (1972) (Volumetric) with beta=-2

The material is a special case of a generalised pover law neo-Hookean material
with a damage model inside relative only to the isochoric part.

the input line should read
MAT 3 MAT_Raghavan_Damage BULK 0.120755 ALPHA 0.068632  BETA 5.799445 EQSTRMIN 0.206141 EQSTRMAX
0.6444974  A 2.4816557 B  0.478626783 DENS 0.001

\level 3

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/

#include <vector>
#include "aaaraghavanvorp_damage.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::AAAraghavanvorp_damage::AAAraghavanvorp_damage(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),

      bulk_(matdata->GetDouble("BULK")),          /// Bulk's modulus (Volumetric)
      alpha_(matdata->GetDouble("ALPHA")),        /// 1st parameter, alpha (Isochoric)
      beta_(matdata->GetDouble("BETA")),          /// 2nd parameter, beta (Isochoric)
      eqstrmin_(matdata->GetDouble("EQSTRMIN")),  /// equivalent strain initial damage
      a_(matdata->GetDouble("A")),                /// 1st parameter, a
      b_(matdata->GetDouble("B")),                /// 2nd parameter, b
      density_(matdata->GetDouble("DENS"))        /// Density
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::AAAraghavanvorp_damage::CreateMaterial()
{
  return Teuchos::rcp(new MAT::AAAraghavanvorp_damage(this));
}


MAT::AAAraghavanvorp_damageType MAT::AAAraghavanvorp_damageType::instance_;


DRT::ParObject* MAT::AAAraghavanvorp_damageType::Create(const std::vector<char>& data)
{
  MAT::AAAraghavanvorp_damage* aaadamage = new MAT::AAAraghavanvorp_damage();
  aaadamage->Unpack(data);
  return aaadamage;  // aaadam;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  ^_^gm 05/09 |
 *----------------------------------------------------------------------*/
MAT::AAAraghavanvorp_damage::AAAraghavanvorp_damage() : params_(NULL)
{
  isinit_ = false;  ///< indicates if material is initialized by calling the #Initialized routine
  // damage history parameters
  histgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);  ///< current damage parameter
  histglast_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);  ///< damage of last converged state
  histeqstrmaxcurr_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);  ///< current damage parameter
  histeqstrmaxlast_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);  ///< damage of last converged state
  elstrength_ = Teuchos::rcp(new double);                   ///< element strength
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)   ^_^gm 05/09 |
 *----------------------------------------------------------------------*/
MAT::AAAraghavanvorp_damage::AAAraghavanvorp_damage(MAT::PAR::AAAraghavanvorp_damage* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  ^_^gm 05/09 |
 *----------------------------------------------------------------------*/
void MAT::AAAraghavanvorp_damage::Pack(DRT::PackBuffer& data) const
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
    histsize = histglast_->size();
  }
  AddtoPack(data, histsize);  // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    AddtoPack(data, histglast_->at(var));
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)   ^_^gm 05/09 |
 *----------------------------------------------------------------------*/
void MAT::AAAraghavanvorp_damage::Unpack(const std::vector<char>& data)
{
  std::cout << "UNPACK \n";
  // isinit_=true; // corretto ore 14:27
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
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
        params_ = static_cast<MAT::PAR::AAAraghavanvorp_damage*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // history data
  int twicehistsize;
  ExtractfromPack(position, data, twicehistsize);

  if (twicehistsize == 0) isinit_ = false;

  histgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);
  histglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);
  //  histgcurr_= Teuchos::rcp(new std::vector<double > );   ///< current damage parameter
  //  histglast_= Teuchos::rcp(new std::vector<double> );   ///< damage of last converged state


  for (int var = 0; var < twicehistsize; ++var)
  {
    LINALG::Matrix<1, 1> tmp(true);
    // double tmp= 0.0;
    histgcurr_->push_back(tmp);            // current vectors have to be initialized
    ExtractfromPack(position, data, tmp);  // last vectors are unpacked
    histglast_->push_back(tmp);
  }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}
/*----------------------------------------------------------------------*
 |  Setup: Initialize/allocate internal g variables (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::AAAraghavanvorp_damage::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  double strength = 0.0;  // section for extracting the element strength
  linedef->ExtractDouble("STRENGTH", strength);

  std::cout << "SETUP \n";
  histgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);
  histglast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);
  //  histgcurr_= Teuchos::rcp(new std::vector<double> );   ///< current damage parameter
  //  histglast_= Teuchos::rcp(new std::vector<double> );   ///< damage of last converged state

  histeqstrmaxcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);
  histeqstrmaxlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);
  //  histeqstrmaxcurr_=Teuchos::rcp(new std::vector<double> );
  //  histeqstrmaxlast_=Teuchos::rcp(new std::vector<double> );

  elstrength_ = Teuchos::rcp(new double);
  //  double emptyvec=0.0;
  //  double emptyvec1=1.0;

  LINALG::Matrix<1, 1> emptyvec(true);
  LINALG::Matrix<1, 1> emptyvec1(true);
  emptyvec1(0, 0) = 1.0;
  histgcurr_->resize(numgp);
  histglast_->resize(numgp);
  histeqstrmaxcurr_->resize(numgp);
  histeqstrmaxlast_->resize(numgp);

  for (int j = 0; j < numgp; ++j)
  {
    histgcurr_->at(j) = emptyvec1;
    histglast_->at(j) = emptyvec1;
    histeqstrmaxcurr_->at(j) = 0.0;
    histeqstrmaxlast_->at(j) = params_->eqstrmin_;
  }

  *elstrength_ = strength;  // uncomment when everything is working
                            //  *elstrength_ = 1.0;
  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
 |  Update: internal stress variables              (public)         05/08|
 *----------------------------------------------------------------------*/

void MAT::AAAraghavanvorp_damage::Update()
{
  // std::cout << " fcg"; // std::string to declare when the function is called,normally after
  // convergence.
  histglast_ = histgcurr_;
  histeqstrmaxlast_ = histeqstrmaxcurr_;

  const LINALG::Matrix<1, 1> emptyvec(true);
  //  double emptyvec=0.0;
  //  histgcurr_= Teuchos::rcp(new std::vector<double > );   ///< current damage parameter
  histgcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);
  histeqstrmaxcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1, 1>>);
  //  histeqstrmaxcurr_=Teuchos::rcp(new std::vector<double>);
  const int numgp = histglast_->size();
  histgcurr_->resize(numgp);
  histeqstrmaxcurr_->resize(numgp);
  for (int j = 0; j < numgp; ++j)
  {
    histgcurr_->at(j) = emptyvec;
    histeqstrmaxcurr_->at(j) = emptyvec;
    //     histeqstrmaxcurr_->at(j)=0.0;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Transform Second Piola-Kirchhoff Stress Tensor in the Cauchy one(public)         05/08|
 *----------------------------------------------------------------------*/

void MAT::AAAraghavanvorp_damage::StressTensTransfSPKtoCauchy(
    LINALG::Matrix<NUM_STRESS_3D, 1>& f,       ///< deformation gradient tensor
    const double detf,                         ///< determinant of deformation gradient tensor
    LINALG::Matrix<NUM_STRESS_3D, 1>& pktwo,   ///< 2nd PK-stress
    LINALG::Matrix<NUM_STRESS_3D, 1>& cstress  ///< Cauchy-stress
)
{
  cstress(0) = (f(0) * pktwo(0) + f(3) * pktwo(3) + f(5) * pktwo(5)) * f(0) +
               (f(0) * pktwo(3) + f(3) * pktwo(1) + f(5) * pktwo(4)) * f(3) +
               (f(0) * pktwo(5) + f(3) * pktwo(4) + f(5) * pktwo(2)) * f(5);  // cstress11

  cstress(3) = (f(0) * pktwo(0) + f(3) * pktwo(3) + f(5) * pktwo(5)) * f(3) +
               (f(0) * pktwo(3) + f(3) * pktwo(1) + f(5) * pktwo(4)) * f(1) +
               (f(0) * pktwo(5) + f(3) * pktwo(4) + f(5) * pktwo(2)) * f(4);  // cstress12

  cstress(5) = (f(0) * pktwo(0) + f(3) * pktwo(3) + f(5) * pktwo(5)) * f(5) +
               (f(0) * pktwo(3) + f(3) * pktwo(1) + f(5) * pktwo(4)) * f(4) +
               (f(0) * pktwo(5) + f(3) * pktwo(4) + f(5) * pktwo(2)) * f(2);  // cstress13

  cstress(1) = (f(3) * pktwo(0) + f(1) * pktwo(3) + f(4) * pktwo(5)) * f(3) +
               (f(3) * pktwo(3) + f(1) * pktwo(1) + f(4) * pktwo(4)) * f(1) +
               (f(3) * pktwo(5) + f(1) * pktwo(4) + f(4) * pktwo(2)) * f(4);  // cstress22

  cstress(4) = (f(3) * pktwo(0) + f(1) * pktwo(3) + f(4) * pktwo(5)) * f(5) +
               (f(3) * pktwo(3) + f(1) * pktwo(1) + f(4) * pktwo(4)) * f(4) +
               (f(3) * pktwo(5) + f(1) * pktwo(4) + f(4) * pktwo(2)) * f(2);  // cstress23

  cstress(2) = (f(5) * pktwo(0) + f(4) * pktwo(3) + f(2) * pktwo(5)) * f(5) +
               (f(5) * pktwo(3) + f(4) * pktwo(1) + f(2) * pktwo(4)) * f(4) +
               (f(5) * pktwo(5) + f(4) * pktwo(4) + f(2) * pktwo(2)) * f(2);  // cstress33

  cstress.Scale(1.0 / detf);
  return;
}

/*----------------------------------------------------------------------*
 |  Transform Cauchy Stress Tensor in the Second Piola-Kirchhoff one(public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::AAAraghavanvorp_damage::StressTensTransfCauchytoSPK(
    LINALG::Matrix<NUM_STRESS_3D, 1>& invf,     ///< deformation gradient tensor
    const double detf,                          ///< determinant of deformation gradient tensor
    LINALG::Matrix<NUM_STRESS_3D, 1>& cstress,  ///< Cauchy-stress
    LINALG::Matrix<NUM_STRESS_3D, 1>& pktwo     ///< 2nd PK-stress
)
{
  pktwo(0) = (invf(0) * cstress(0) + invf(3) * cstress(3) + invf(5) * cstress(5)) * invf(0) +
             (invf(0) * cstress(3) + invf(3) * cstress(1) + invf(5) * cstress(4)) * invf(3) +
             (invf(0) * cstress(5) + invf(3) * cstress(4) + invf(5) * cstress(2)) * invf(5);  // p11

  pktwo(3) = (invf(0) * cstress(0) + invf(3) * cstress(3) + invf(5) * cstress(5)) * invf(3) +
             (invf(0) * cstress(3) + invf(3) * cstress(1) + invf(5) * cstress(4)) * invf(1) +
             (invf(0) * cstress(5) + invf(3) * cstress(4) + invf(5) * cstress(2)) * invf(4);  // p12

  pktwo(5) = (invf(0) * cstress(0) + invf(3) * cstress(3) + invf(5) * cstress(5)) * invf(5) +
             (invf(0) * cstress(3) + invf(3) * cstress(1) + invf(5) * cstress(4)) * invf(4) +
             (invf(0) * cstress(5) + invf(3) * cstress(4) + invf(5) * cstress(2)) * invf(2);  // p13

  pktwo(1) = (invf(3) * cstress(0) + invf(1) * cstress(3) + invf(4) * cstress(5)) * invf(3) +
             (invf(3) * cstress(3) + invf(1) * cstress(1) + invf(4) * cstress(4)) * invf(1) +
             (invf(3) * cstress(5) + invf(1) * cstress(4) + invf(4) * cstress(2)) * invf(4);  // p22

  pktwo(4) = (invf(3) * cstress(0) + invf(1) * cstress(3) + invf(4) * cstress(5)) * invf(5) +
             (invf(3) * cstress(3) + invf(1) * cstress(1) + invf(4) * cstress(4)) * invf(4) +
             (invf(3) * cstress(5) + invf(1) * cstress(4) + invf(4) * cstress(2)) * invf(2);  // p23

  pktwo(2) = (invf(5) * cstress(0) + invf(4) * cstress(3) + invf(2) * cstress(5)) * invf(5) +
             (invf(5) * cstress(3) + invf(4) * cstress(1) + invf(2) * cstress(4)) * invf(4) +
             (invf(5) * cstress(5) + invf(4) * cstress(4) + invf(2) * cstress(2)) * invf(2);  // p33

  pktwo.Scale(detf);
  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                   (public)  gm^_^ 03/08 gee 10/08 |
 *----------------------------------------------------------------------*

 plain strain energy function

 W    = alpha (Ic*IIIc^(-1/3) -3) + beta (Ic*IIIc^(-1/3)-3)^2

 taken from
 M.L. Raghavan, D.A. Vorp: Toward a biomechanical tool to evaluate rupture potential
 of abdominal aortic aneurysm: identification of a finite strain constitutive model
 and evaluation of its applicability, J. of Biomechanics 33 (2000) 475-482.

 and modified to slight compressibility

 here

 Ic   .. first invariant of right Cauchy-Green tensor C
 IIIc .. third invariant of right Cauchy-Green tensor C

 The volumetric part is done by a volumetric strain energy function taken from
 Holzapfel(p244), Simo & Miehe (1991)

 W_vol = K beta2^(-2) ( beta2 ln (J) + J^(-beta2) -1 )

 where

 K    .. bulk modulus
 beta2 = -2 a parameter according to Holzapfel (p244) ,Simo & Miehe (1991)
 J    .. det(F) determinante of the Jacobian matrix, IIIc^0.5

  */
void MAT::AAAraghavanvorp_damage::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<NUM_STRESS_3D, 1>* glstrain,  ///< green lagrange strain
    Teuchos::ParameterList& params,                    ///< parameter list for communication
    LINALG::Matrix<NUM_STRESS_3D, 1>* stress, LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,
    const int eleGID)
{
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  // material parameters for volumetric part
  const double bulk = params_->bulk_;  // Bulk's modulus(Volumetric)
  const double beta2 = -2;             // parameter from Holzapfel (p244), Simo & Miehe (1991)

  // material parameters for isochoric part
  const double alpha = params_->alpha_;  // 1st parameter (Isochoric)
  const double beta = params_->beta_;    // 2nd parameter (Isochoric)

  // material damage parameters
  const double eqstrmin = params_->eqstrmin_;
  const double a = params_->a_;  /// 1st parameter, a
  const double b = params_->b_;  /// 2nd parameter, b
  const double strength = *elstrength_;

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D, 1> identity(true);
  for (int i = 0; i < 3; i++) identity(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D, 1> rcg(*glstrain);
  rcg.Scale(2.0);
  rcg += identity;

  //  if (!gp) printf("strXX %15.10e",sqrt(rcg(0,0)));

  for (int i = 0; i < NUM_STRESS_3D; i++)  // in the first step could be problem without it
    if ((rcg(i) < 0) | (rcg(i) < 1.0e-12)) rcg(i) = 0.0;

  // invariants
  double inv = rcg(0) + rcg(1) + rcg(2);  // 1st invariant, trace
  double iiinv = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                 0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                 0.25 * rcg(0) * rcg(4) * rcg(4);  // 3rd invariant, determinante

  double detf = 0.0;  // J
  if (iiinv < 0.0)
    dserror("fatal failure in aneurysmatic artery wall material");
  else
    detf = sqrt(iiinv);  // determinant of deformation gradient

  //--------------------------------------------------------------------------------------
  // invert C
  LINALG::Matrix<NUM_STRESS_3D, 1> invc(true);

  double invdet = 1. / iiinv;

  invc(0) = rcg(1) * rcg(2) - 0.25 * rcg(4) * rcg(4);
  invc(1) = rcg(0) * rcg(2) - 0.25 * rcg(5) * rcg(5);
  invc(2) = rcg(0) * rcg(1) - 0.25 * rcg(3) * rcg(3);
  invc(3) = 0.25 * rcg(5) * rcg(4) - 0.5 * rcg(3) * rcg(2);
  invc(4) = 0.25 * rcg(3) * rcg(5) - 0.5 * rcg(0) * rcg(4);
  invc(5) = 0.25 * rcg(3) * rcg(4) - 0.5 * rcg(5) * rcg(1);

  invc.Scale(invdet);

  //--- prepare some constants -----------------------------------------------------------
  const double third = 1.0 / 3.0;
  const double twthi = 2.0 / 3.0;


  LINALG::Matrix<NUM_STRESS_3D, 1> f(true);  // deformation gradient tensor F, init to 0

  f(0) = sqrt(rcg(0));        // F11
  f(1) = sqrt(rcg(1));        // F22
  f(2) = sqrt(rcg(2));        // F33
  f(3) = sqrt(rcg(3) / 2.0);  // F12
  f(4) = sqrt(rcg(4) / 2.0);  // F23
  f(5) = sqrt(rcg(5) / 2.0);  // F31

  LINALG::Matrix<NUM_STRESS_3D, 1> invf(
      true);  // inverse of deformation gradient tensor F, init to 0

  invf(0) = f(1) * f(2) - f(4) * f(4);  // invF11
  invf(1) = f(0) * f(2) - f(5) * f(5);  // invF22
  invf(2) = f(0) * f(1) - f(3) * f(3);  // invF33
  invf(3) = f(5) * f(4) - f(3) * f(2);  // invF12
  invf(4) = f(3) * f(5) - f(0) * f(4);  // invF23
  invf(5) = f(3) * f(4) - f(5) * f(1);  // invF31

  // LINALG::Matrix<NUM_STRESS_3D,1> fvol(true); // deformation gradient tensor F, init to 0
  // for(int i=0;i<3;i++)
  //  fvol(i)=pow(detf,third);
  // LINALG::Matrix<NUM_STRESS_3D,1> invfvol(true); // deformation gradient tensor F, init to 0
  // for(int i=0;i<3;i++)
  //  fvol(i)=pow(detf,-third);

  // LINALG::Matrix<NUM_STRESS_3D,1> fiso(f); // deformation gradient tensor F, init to 0
  // fiso.Scale(pow(invdet,-third));
  // LINALG::Matrix<NUM_STRESS_3D,1> invfiso(invf); // deformation gradient tensor F, init to 0
  // fiso.Scale(pow(invdet,third));



  //--- determine 2nd Piola Kirchhoff stresses pktwo -------------------------------------

  // 1st step: volumetric part
  //==========================
  LINALG::Matrix<NUM_STRESS_3D, 1> pktwovol(invc);
  double scalar = bulk / beta2 * (1.0 - pow(detf, -beta2));

  pktwovol.Scale(scalar);  // initialise PKtwo with volumetric part


  // LINALG::Matrix<NUM_STRESS_3D,1> sigmvol(true); // Cauchy stress tensor, init to 0
  // StressTensTransfSPKtoCauchy(fvol,detf,pktwovol,sigmvol);
  // sigmvol.Scale(strength);
  // StressTensTransfCauchytoSPK(invfvol,detf,sigmvol,pktwovol);


  // 2nd step: isochoric part
  //=========================
  double isochor1 =
      2.0 * (alpha * pow(iiinv, third) + 2.0 * beta * inv - 6.0 * beta * pow(iiinv, third)) *
      pow(iiinv, -twthi);

  double isochor2 =
      -twthi * inv *
      (alpha * pow(iiinv, third) + 2.0 * beta * inv - 6.0 * beta * pow(iiinv, third)) *
      pow(iiinv, -twthi);
  // contribution: Cinv
  LINALG::Matrix<NUM_STRESS_3D, 1> pktwoiso(invc);
  pktwoiso.Scale(isochor2);

  // contribution: I
  for (int i = 0; i < 3; i++) pktwoiso(i) += isochor1;

  // LINALG::Matrix<NUM_STRESS_3D,1> sigmiso(true); // Cauchy stress tensor, init to 0
  // StressTensTransfSPKtoCauchy(fvol,1.0,pktwoiso,sigmiso);
  // sigmiso.Scale(strength);
  // StressTensTransfCauchytoSPK(invfiso,1.0,sigmvol,pktwoiso);


  //--- do elasticity matrix -------------------------------------------------------------
  // ensure that cmat is zero when it enters the computation
  // It is an implicit law that cmat is zero upon input
  // cmat.PutScalar(0.0);

  // 3rd step: volumetric part
  //==========================
  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatvol(true);  // volumetric elasticity tensor
  double coeff1 = bulk * pow(detf, -beta2);                    // coefficients
  double coeff2 = -2.0 * scalar;

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++) cmatvol(i, j) = coeff1 * invc(i) * invc(j);

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct(cmatvol, invc, coeff2);

  // 4th step: isochoric part
  //=========================
  LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatiso(true);  // isochoric elasticity tensor

  // deltas (see also Holzapfel p.261)
  // note that these deltas serve for the isochoric part only
  double delta1 = 8.0 * beta * pow(iiinv, -twthi);
  double delta3 =
      -4. / 3 * (alpha + 4. * beta * inv * pow(iiinv, -third) - 6 * beta) * pow(iiinv, -third);
  double delta6 =
      4. / 9 * inv * (alpha + 4. * beta * inv * pow(iiinv, -third) - 6 * beta) * pow(iiinv, -third);
  double delta7 =
      4. / 3 * inv * (alpha + 2. * beta * inv * pow(iiinv, -third) - 6 * beta) * pow(iiinv, -third);

  // contribution: I \obtimes I
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) cmatiso(i, j) = delta1;

  // contribution: Cinv \otimes Cinv
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
    {
      cmatiso(i, j) +=
          delta3 * (identity(i) * invc(j) +
                       invc(i) * identity(j));      // contribution: Cinv \otimes I + I \otimes Cinv
      cmatiso(i, j) += delta6 * invc(i) * invc(j);  // contribution: Cinv \otimes Cinv
    }

  // contribution: boeppel-product
  AddtoCmatHolzapfelProduct(cmatiso, invc, delta7);

  // 5th step: damage evaluation
  //============================
  //  double sefvol=bulk*pow(beta2,-2)*(beta2*log(detf)+pow(detf,-beta2)-1);
  double sefiso =
      alpha * (inv * pow(iiinv, -third) - 3) + beta * pow((inv * pow(iiinv, -third) - 3), 2);

  if (sefiso < 0) sefiso = 0.0;  // for the first step

  double eqstr = sqrt(2.0 * sefiso);
  double hgl = histglast_->at(gp)(0, 0);
  double hgc = 1.0;
  double hesml = histeqstrmaxlast_->at(gp)(0, 0);  // eqstrmin-log(1-(1-hgl)/a)/b;

  if (eqstr < hesml)
  {
    hgc = hgl;

    histgcurr_->at(gp) = hgl;           // storing of actual damage parameter value
    histeqstrmaxcurr_->at(gp) = hesml;  // storing of actual maximum equivalent strain

    cmatiso.Scale(hgc);
  }
  else
  {
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisodam(
        true);  // damage part of isochoric elasticity tensor
    double gdot = 0.0;

    hgc = 1 - a * (1 - exp(-b * (eqstr - eqstrmin)));  // damage parameter evaluation
    gdot = -a * b * exp(-b * (eqstr - eqstrmin)) /
           eqstr;  // derivative of damage function, it�s used to scale cmatisodam

    histgcurr_->at(gp)(0, 0) = hgc;     // storing of actual damage parameter value
    histeqstrmaxcurr_->at(gp) = eqstr;  // storing of actual maximum equivalent strain

    for (int i = 0; i < 6; i++)  // dyadic product
      for (int j = 0; j < 6; j++) cmatisodam(i, j) = pktwoiso(i) * pktwoiso(j);

    cmatiso.Scale(hgc);
    cmatisodam.Scale(gdot);
    cmatiso += cmatisodam;
  }

  // 6th step: final step
  //==========================

  pktwovol.Scale(strength);
  pktwoiso.Scale(hgc * strength);

  *stress = pktwovol;
  *stress += pktwoiso;
  cmatvol.Scale(strength);
  cmatiso.Scale(strength);
  *cmat = cmatvol;
  *cmat += cmatiso;

  /*
    // 7th step: rescaling step - DISASTER ZONE
    //==========================
    LINALG::Matrix<NUM_STRESS_3D,1> pktwo(*stress); // Cauchy stress tensor, init to 0

    LINALG::Matrix<NUM_STRESS_3D,1> f(true); // deformation gradient tensor F, init to 0

    f(0)=sqrt(rcg(0)); 		// F11
    f(1)=sqrt(rcg(1)); 		// F22
    f(2)=sqrt(rcg(2));		// F33
    f(3)=sqrt(rcg(3)/2.0);	// F12
    f(4)=sqrt(rcg(4)/2.0);	// F23
    f(5)=sqrt(rcg(5)/2.0);	// F31

    LINALG::Matrix<NUM_STRESS_3D,1> invf(true); // inverse of deformation gradient tensor F, init to
    0

    invf(0) = f(1)*f(2) - f(4)*f(4);     // invF11
    invf(1) = f(0)*f(2) - f(5)*f(5);     // invF22
    invf(2) = f(0)*f(1) - f(3)*f(3);     // invF33
    invf(3) = f(5)*f(4) - f(3)*f(2);     // invF12
    invf(4) = f(3)*f(5) - f(0)*f(4);     // invF23
    invf(5) = f(3)*f(4) - f(5)*f(1);     // invF31

    LINALG::Matrix<NUM_STRESS_3D,1> fvol(true); // deformation gradient tensor F, init to 0
    for(int =0;i<3;i++)
     fvol(i)=pow(detf,third);
    LINALG::Matrix<NUM_STRESS_3D,1> invfvol(true); // deformation gradient tensor F, init to 0
    for(int =0;i<3;i++)
     fvol(i)=pow(detf,-third);

    LINALG::Matrix<NUM_STRESS_3D,1> fiso(f); // deformation gradient tensor F, init to 0
    fiso.Scale(pow(invdet,-third));
    LINALG::Matrix<NUM_STRESS_3D,1> invfiso(invf); // deformation gradient tensor F, init to 0
    fiso.Scale(pow(invdet,third));


    LINALG::Matrix<NUM_STRESS_3D,1> sigm(true); // Cauchy stress tensor, init to 0

    StressTensTransfSPKtoCauchy(f,detf,pktwo,sigm);

    sigm.Scale(strength);

    LINALG::Matrix<NUM_STRESS_3D,1> realpktwo(true); // Cauchy stress tensor, init to 0

    StressTensTransfCauchytoSPK(invf,detf,sigm,realpktwo);
   */
  return;
}
