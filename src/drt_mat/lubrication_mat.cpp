/*--------------------------------------------------------------------------*/
/*!
\file lubrication_mat.cpp

\brief Material model for the lubrication film

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/


#include <vector>
#include "lubrication_mat.H"
#include "lubrication_law.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::LubricationMat::LubricationMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      density_(matdata->GetDouble("DENSITY")),
      lubricationlawID_(matdata->GetInt("LUBRICATIONLAWID")),
      lubricationlaw_(NULL)
{
  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(lubricationlawID_);

  switch (curmat->Type())
  {
    case INPAR::MAT::m_lubrication_law_constant:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::LubricationLawConstant(curmat));
      lubricationlaw_ = static_cast<MAT::PAR::LubricationLaw*>(curmat->Parameter());
      break;
    }
    default:
      dserror("invalid material for porosity law %d", curmat->Type());
      break;
  }
}


Teuchos::RCP<MAT::Material> MAT::PAR::LubricationMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::LubricationMat(this));
}

MAT::LubricationMatType MAT::LubricationMatType::instance_;

DRT::ParObject* MAT::LubricationMatType::Create(const std::vector<char>& data)
{
  MAT::LubricationMat* lubrication_mat = new MAT::LubricationMat();
  lubrication_mat->Unpack(data);
  return lubrication_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::LubricationMat::LubricationMat() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::LubricationMat::LubricationMat(MAT::PAR::LubricationMat* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::LubricationMat::Pack(DRT::PackBuffer& data) const
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::LubricationMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror(
        "wrong instance type data. type = %d, UniqueParObjectId()=%d", type, UniqueParObjectId());

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
        params_ = static_cast<MAT::PAR::LubricationMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
                                                              wirtz 09/16|
*----------------------------------------------------------------------*/
double MAT::LubricationMat::ComputeViscosity(const double press)
{
  //  viscosity = params_->viscosity_;
  //  return;

  double visc = -1.;
  params_->lubricationlaw_->ComputeViscosity(press, visc);

  return visc;
}


/*----------------------------------------------------------------------*
                                                              wirtz 09/16|
*----------------------------------------------------------------------*/
double MAT::LubricationMat::ComputeViscosityDeriv(const double press)
{
  double visc = -1.;
  double visc_dp = -1;
  params_->lubricationlaw_->ConstitutiveDerivatives(press, visc, visc_dp);

  return visc_dp;
}
