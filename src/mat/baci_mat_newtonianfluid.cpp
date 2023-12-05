/*----------------------------------------------------------------------*/
/*! \file
\brief Newtonian fluid material

\level 1

*/
/*----------------------------------------------------------------------*/


#include "baci_mat_newtonianfluid.H"

#include "baci_lib_globalproblem.H"
#include "baci_mat_par_bundle.H"

#include <vector>



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


CORE::COMM::ParObject* MAT::NewtonianFluidType::Create(const std::vector<char>& data)
{
  MAT::NewtonianFluid* fluid = new MAT::NewtonianFluid();
  fluid->Unpack(data);
  return fluid;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::NewtonianFluid::NewtonianFluid() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::NewtonianFluid::NewtonianFluid(MAT::PAR::NewtonianFluid* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::NewtonianFluid::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::NewtonianFluid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
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
