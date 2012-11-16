/*----------------------------------------------------------------------*/
/*!
\file logneohooke.cpp
\brief


the input line should read
  MAT 1 MAT_Struct_LogNeoHooke YOUNG 1.044E7 NUE 0.3 DENS 1.0

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "logneohooke.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::LogNeoHooke::LogNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  nue_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  model_(matdata->GetInt("MODEL"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::LogNeoHooke::CreateMaterial()
{
  return Teuchos::rcp(new MAT::LogNeoHooke(this));
}


MAT::LogNeoHookeType MAT::LogNeoHookeType::instance_;


DRT::ParObject* MAT::LogNeoHookeType::Create( const std::vector<char> & data )
{
  MAT::LogNeoHooke* logneo = new MAT::LogNeoHooke();
  logneo->Unpack(data);
  return logneo;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::LogNeoHooke::LogNeoHooke()
  : params_(NULL)
{
  dserror("This material law - LOGNEOHOOKE - is maintained only inside the Elasthyper Toolbox.\n"
    "If you want to use this law, the material input line should read :\n"
    "MAT 1   MAT_ElastHyper   NUMMAT 1 MATIDS 2 DENS 0\n"
    "MAT 2 ELAST_CoupLogNeoHooke MODE YN C1 1.0 C2 0.3\n"
    "or\n"
    "MAT 2 ELAST_CoupLogNeoHooke MODE Lame C1 1.0 C2 1.0");
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::LogNeoHooke::LogNeoHooke(MAT::PAR::LogNeoHooke* params)
  : params_(params)
{
  dserror("This material law - LOGNEOHOOKE - is maintained only inside the Elasthyper Toolbox.\n"
    "If you want to use this law, the material input line should read :\n"
    "MAT 1   MAT_ElastHyper   NUMMAT 1 MATIDS 2 DENS 0\n"
    "MAT 2 ELAST_CoupLogNeoHooke MODE YN C1 1.0 C2 0.3\n"
    "or\n"
    "MAT 2 ELAST_CoupLogNeoHooke MODE Lame C1 1.0 C2 1.0");
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::LogNeoHooke::Pack(DRT::PackBuffer& data) const
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
}

/*----------------------------------------------------------------------*
 |  Unpack                                        (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
void MAT::LogNeoHooke::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::LogNeoHooke*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::LogNeoHooke::Evaluate(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress
  )
{
  const int model = params_->model_;
  if (model == 0)
    EvaluateModelZero(glstrain, cmat, stress);
  else if (model == 1)
    EvaluateModelOne(glstrain, cmat, stress);
  else
    dserror("Cannot handle submodel of type %d", model);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::LogNeoHooke::EvaluateModelZero(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress
  )
{
  // material parameters for isochoric part
  const double youngs = params_->youngs_;  // Young's modulus
  const double nue = params_->nue_;  // Poisson's ratio
  const double lambda = (nue==0.5) ? 0.0 : youngs*nue/((1.0+nue)*(1.0-2.0*nue));  // Lame coeff.
  const double mue = youngs/(2.0*(1.0+nue));  // shear modulus

  // build identity tensor I
  LINALG::Matrix<6,1> identity(true);
  for (int i = 0; i < 3; i++)
    identity(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<6,1> rcg(glstrain);
  rcg.Scale(2.0);
  rcg.Update(1.0, identity, 1.0);

  // invariants
  // 1st invariant, trace
  //const double inv = rcg(0) + rcg(1) + rcg(2);
  // 3rd invariant, determinant
  const double iiinv = rcg(0)*rcg(1)*rcg(2)
                     + 0.25 * rcg(3)*rcg(4)*rcg(5)
                     - 0.25 * rcg(1)*rcg(5)*rcg(5)
                     - 0.25 * rcg(2)*rcg(3)*rcg(3)
                     - 0.25 * rcg(0)*rcg(4)*rcg(4);
  if (iiinv < 0.0)
    dserror("fatal failure in logarithmic neo-Hooke material");
  // determinant of deformation gradient
  const double detf = std::sqrt(iiinv);

  // invert right Cauchy-Green tensor
  LINALG::Matrix<6,1> invc(false);
  {
    invc(0) = ( rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4) ) / iiinv;
    invc(1) = ( rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5) ) / iiinv;
    invc(2) = ( rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3) ) / iiinv;
    invc(3) = ( 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2) ) / iiinv;
    invc(4) = ( 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4) ) / iiinv;
    invc(5) = ( 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1) ) / iiinv;
  }

  // 2nd Piola Kirchhoff stresses
  {
    LINALG::Matrix<6,1> pk2(identity);
    pk2.Scale(mue);
    pk2.Update(-mue+lambda*std::log(detf), invc, 1.0);
    stress.Update(pk2);
  }

  // constitutive tensor
  // It is an implicit law that cmat is zero upon input
  {
    // deltas (see also Holzapfel [2] at p.261)
    const double delta6 = lambda;
    const double delta7 = 2.0*(mue - lambda*std::log(detf));

    // contribution: Cinv \otimes Cinv
    cmat.MultiplyNT(delta6, invc, invc);

    // contribution: Cinv \odot Cinv
    AddtoCmatHolzapfelProduct(cmat, invc, delta7);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::LogNeoHooke::EvaluateModelOne(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress
  )
{
  // material parameters for isochoric part
  const double youngs = params_->youngs_;  // Young's modulus
  const double nue = params_->nue_;  // Poisson's ratio
  const double lambda = (nue==0.5) ? 0.0 : youngs*nue/((1.0+nue)*(1.0-2.0*nue));  // Lame coeff.
  const double mue1 = youngs/(2.0*(1.0+nue));  // shear modulus at I_C-proportional term
  const double mue3 = (nue==0.5) ? 0.0 : mue1;  // shear modulus at III_C-proportional term

  // build identity tensor I
  LINALG::Matrix<6,1> identity(true);
  for (int i = 0; i < 3; i++)
    identity(i) = 1.0;

  // right Cauchy-Green Tensor  C = 2 * E + I
  // REMARK: strain-like 6-Voigt vector
  LINALG::Matrix<6,1> rcg(glstrain);
  rcg.Scale(2.0);
  rcg.Update(1.0, identity, 1.0);

  // invariants
  // 1st invariant, trace
  const double inv = rcg(0) + rcg(1) + rcg(2);
  // 3rd invariant, determinant
  const double iiinv = rcg(0)*rcg(1)*rcg(2)
                     + 0.25 * rcg(3)*rcg(4)*rcg(5)
                     - 0.25 * rcg(1)*rcg(5)*rcg(5)
                     - 0.25 * rcg(2)*rcg(3)*rcg(3)
                     - 0.25 * rcg(0)*rcg(4)*rcg(4);
  if (iiinv < 0.0)
    dserror("fatal failure in logarithmic neo-Hooke material");

  // convenience
  const double iiinvpowthird = std::pow(iiinv,1.0/3.0);

  // invert right Cauchy-Green tensor
  // REMARK: stress-like 6-Voigt vector
  LINALG::Matrix<6,1> icg(false);
  {
    icg(0) = ( rcg(1)*rcg(2) - 0.25*rcg(4)*rcg(4) ) / iiinv;
    icg(1) = ( rcg(0)*rcg(2) - 0.25*rcg(5)*rcg(5) ) / iiinv;
    icg(2) = ( rcg(0)*rcg(1) - 0.25*rcg(3)*rcg(3) ) / iiinv;
    icg(3) = ( 0.25*rcg(5)*rcg(4) - 0.5*rcg(3)*rcg(2) ) / iiinv;
    icg(4) = ( 0.25*rcg(3)*rcg(5) - 0.5*rcg(0)*rcg(4) ) / iiinv;
    icg(5) = ( 0.25*rcg(3)*rcg(4) - 0.5*rcg(5)*rcg(1) ) / iiinv;
  }

  // 2nd Piola Kirchhoff stresses
  {
    LINALG::Matrix<6,1> pk2(false);
    const double gamma1 = mue1/iiinvpowthird;
    pk2.Update(gamma1, identity);
    //pk2.Update(-mue3+lambda*std::log(detf), icg, 1.0);
    const double gamma3 = -mue1*inv/(3.0*iiinvpowthird)
                        - mue3
                        + lambda/2.0*std::log(iiinv);
    pk2.Update(gamma3, icg, 1.0);
    stress.Update(pk2);
  }

  // constitutive tensor
  // It is an implicit law that cmat is zero upon input
  {
    // deltas (see also Holzapfel [2] at p.261)
    const double delta3 = -2.0*mue1/(3.0*iiinvpowthird);
    const double delta6 = 2.0*mue1*inv/(9.0*iiinvpowthird)
                        + lambda;
    const double delta7 = 2.0*mue1*inv/(3.0*iiinvpowthird)
                        + 2.0*mue3
                        - lambda*std::log(iiinv);

    // contribution Id \otimes Cinv + Cinv \otimes Id
    cmat.MultiplyNT(delta3, identity, icg);
    cmat.MultiplyNT(delta3, icg, identity, 1.0);

    // contribution: Cinv \otimes Cinv
    cmat.MultiplyNT(delta6, icg, icg, 1.0);

    // contribution: Cinv \odot Cinv
    AddtoCmatHolzapfelProduct(cmat, icg, delta7);
  }

  return;
}


/*----------------------------------------------------------------------*/
