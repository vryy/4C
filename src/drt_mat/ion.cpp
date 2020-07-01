/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores parameters of ion species in electrolyte solution. Former file of Georg
Bauer

\level 2


*/
/*----------------------------------------------------------------------------*/

#include <vector>
#include "ion.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Ion::Ion(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      valence_(matdata->GetDouble("VALENCE")),
      diffusivity_(matdata->GetDouble("DIFFUSIVITY")),
      densification_(matdata->GetDouble("DENSIFICATION")),
      elimvalence_(matdata->GetDouble("ELIM_VALENCE")),
      elimdiffusivity_(matdata->GetDouble("ELIM_DIFFUSIVITY"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::Ion::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Ion(this));
}

MAT::IonType MAT::IonType::instance_;


DRT::ParObject* MAT::IonType::Create(const std::vector<char>& data)
{
  MAT::Ion* ion = new MAT::Ion();
  ion->Unpack(data);
  return ion;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Ion::Ion() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Ion::Ion(MAT::PAR::Ion* params) : params_(params) {}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Ion::Pack(DRT::PackBuffer& data) const
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

  /*
  for (unsigned i=0;i<data().size();i++)
  std::cout<<"Pack ION: pb["<<i<<"] = "<<(data())[i]<<std::endl;
*/
  /*
  // extract type
  std::vector<char>::size_type posit = 0;
  std::vector<char> pbtest;
  int typio = 0;
  ExtractfromPack(posit,data(),typio);
  std::cout<<"ION Pack: Type will be "<<typio<<std::endl;
*/
  // std::cout<<"Ion Pack: "<<data().size()<<std::endl;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Ion::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::Ion*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}
