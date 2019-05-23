/*!------------------------------------------------------------------------------------------------*
\brief
Former file of Martin Winklmaier

\level 3

\maintainer Martin Kronbichler
 *------------------------------------------------------------------------------------------------*/


#include <vector>
#include "optimization_density.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::TopOptDens::TopOptDens(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      poro_bd_down_(matdata->GetDouble("MINPORO")),
      poro_bd_up_(matdata->GetDouble("MAXPORO")),
      smear_fac_(matdata->GetDouble("SMEARFAC"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::TopOptDens::CreateMaterial()
{
  return Teuchos::rcp(new MAT::TopOptDens(this));
}


MAT::TopOptDensType MAT::TopOptDensType::instance_;


DRT::ParObject* MAT::TopOptDensType::Create(const std::vector<char>& data)
{
  MAT::TopOptDens* topoptdens_mat = new MAT::TopOptDens();
  topoptdens_mat->Unpack(data);
  return topoptdens_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::TopOptDens::TopOptDens() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::TopOptDens::TopOptDens(MAT::PAR::TopOptDens* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::TopOptDens::Pack(DRT::PackBuffer& data) const
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
void MAT::TopOptDens::Unpack(const std::vector<char>& data)
{
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
        params_ = static_cast<MAT::PAR::TopOptDens*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
