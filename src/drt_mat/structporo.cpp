/*!-----------------------------------------------------------------------*
 \file structporo.cpp

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/



#include <vector>
#include "structporo.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::StructPoro::StructPoro(Teuchos::RCP<MAT::PAR::Material> matdata) :
  Parameter(matdata),
  matid_(matdata->GetInt("MATID")),
  bulkmodulus_(matdata->GetDouble("BULKMODULUS")),
  penaltyparameter_(matdata->GetDouble("PENALTYPARAMETER")),
  initporosity_(matdata->GetDouble("INITPOROSITY"))
{
  mat_ = MAT::Material::Factory(matid_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::StructPoro::CreateMaterial()
{
  return Teuchos::rcp(new MAT::StructPoro(this));
}

MAT::StructPoroType MAT::StructPoroType::instance_;

DRT::ParObject* MAT::StructPoroType::Create(const std::vector<char> & data)
{
  MAT::StructPoro* struct_poro = new MAT::StructPoro();
  struct_poro->Unpack(data);
  return struct_poro;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoro::StructPoro() :
  params_(NULL)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoro::StructPoro(MAT::PAR::StructPoro* params) :
  params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL)
    matid = params_->Id(); // in case we are in post-process mode
  AddtoPack(data, matid);

  // porosity_
  int size=0;
  size = (int)porosity_.size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    AddtoPack(data,(porosity_)[i]);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data");

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
        params_ = static_cast<MAT::PAR::StructPoro*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // porosity_
  int size = 0;
  ExtractfromPack(position,data,size);
  porosity_.resize(size);
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,(porosity_)[i]);

  if (position != data.size())
  dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
double MAT::StructPoro::ComputePorosity(double press, double J,
    const double initporosity, int gp) const
{
  // this function is not called yet!!

  const double bulkmodulus = Bulkmodulus();
  const double penalty = Penaltyparameter();

  const double a = (bulkmodulus / (1 - initporosity) + press - penalty
      / initporosity) * J;
  const double b = -a + bulkmodulus + penalty;
  const double c = (b / a) * (b / a) + 4 * penalty / a;
  double d = sqrt(c) * a;
  double sign = 1.0;

  double test = 1 / (2 * a) * (-b + d);
  if (test >= 1.0 or test < 0.0)
  {
    sign = -1.0;
    d = sign * d;
  }

  const double porosity = 1 / (2 * a) * (-b + d);

  if (porosity >= 1.0 or porosity < 0.0)
  {
    dserror("invalid porosity!");
  }

  return porosity;
}

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
double MAT::StructPoro::PorosityAv() const
{
  double porosityav = 0.0;

  std::vector<double>::const_iterator m;
  for (m = porosity_.begin(); m != porosity_.end(); ++m)
  {
    porosityav += *m;
  }
  porosityav = porosityav / (porosity_.size());

  return porosityav;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::SetPorosityAtGP(std::vector<double> porosity_gp)
{
  int numgp = porosity_gp.size();

  porosity_.resize(numgp);

  //set porosity values
  std::vector<double>::iterator m = porosity_gp.begin();
  for (int i = 0; m != porosity_gp.end(); ++m, ++i)
  {
    double porosity = *m;
    porosity_[i]=porosity;
  }

  return;
}

