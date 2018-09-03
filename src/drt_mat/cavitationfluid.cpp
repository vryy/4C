/*----------------------------------------------------------------------*/
/*!
\file cavitationfluid.cpp

\brief cavitation fluid material

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/


#include <vector>
#include "cavitationfluid.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::CavitationFluid::CavitationFluid(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->GetDouble("DYNVISCOSITY")),
      density_(matdata->GetDouble("DENSITY")),
      gamma_(matdata->GetDouble("GAMMA")),
      p_vapor_(matdata->GetDouble("PVAPOR"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::CavitationFluid::CreateMaterial()
{
  return Teuchos::rcp(new MAT::CavitationFluid(this));
}


MAT::CavitationFluidType MAT::CavitationFluidType::instance_;


DRT::ParObject* MAT::CavitationFluidType::Create(const std::vector<char>& data)
{
  MAT::CavitationFluid* fluid = new MAT::CavitationFluid();
  fluid->Unpack(data);
  return fluid;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::CavitationFluid::CavitationFluid() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::CavitationFluid::CavitationFluid(MAT::PAR::CavitationFluid* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::CavitationFluid::Pack(DRT::PackBuffer& data) const
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
void MAT::CavitationFluid::Unpack(const std::vector<char>& data)
{
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
        params_ = static_cast<MAT::PAR::CavitationFluid*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
