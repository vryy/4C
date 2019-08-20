/*----------------------------------------------------------------------*/
/*! \file
\brief Newtonian fluid material

\level 1

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include "newtonianfluid.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::NewtonianFluid::NewtonianFluid(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->GetDouble("DYNVISCOSITY")),
      density_(matdata->GetDouble("DENSITY")),
      gamma_(matdata->GetDouble("GAMMA"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::NewtonianFluid::CreateMaterial()
{
  return Teuchos::rcp(new MAT::NewtonianFluid(this));
}


MAT::NewtonianFluidType MAT::NewtonianFluidType::instance_;


DRT::ParObject* MAT::NewtonianFluidType::Create(const std::vector<char>& data)
{
  MAT::NewtonianFluid* fluid = new MAT::NewtonianFluid();
  fluid->Unpack(data);
  return fluid;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::NewtonianFluid::NewtonianFluid() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::NewtonianFluid::NewtonianFluid(MAT::PAR::NewtonianFluid* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::NewtonianFluid::Pack(DRT::PackBuffer& data) const
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
void MAT::NewtonianFluid::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::NewtonianFluid*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
