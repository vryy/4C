/*----------------------------------------------------------------------*/
/*!
\file permeablefluid.cpp

<pre>
Maintainer: ???
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include "permeablefluid.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::PermeableFluid::PermeableFluid(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  viscosity_(matdata->GetDouble("VISCOSITY")),
  density_(matdata->GetDouble("DENSITY")),
  permeability_(matdata->GetDouble("PERMEABILITY"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PermeableFluid::CreateMaterial()
{
  return Teuchos::rcp(new MAT::PermeableFluid(this));
}


MAT::PermeableFluidType MAT::PermeableFluidType::instance_;


DRT::ParObject* MAT::PermeableFluidType::Create( const std::vector<char> & data )
{
  MAT::PermeableFluid* permeable_fluid = new MAT::PermeableFluid();
  permeable_fluid->Unpack(data);
  return permeable_fluid;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PermeableFluid::PermeableFluid()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PermeableFluid::PermeableFluid(MAT::PAR::PermeableFluid* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PermeableFluid::Pack(DRT::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PermeableFluid::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::PermeableFluid*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

#endif
