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
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::StructPoro::CreateMaterial()
{
  return Teuchos::rcp(new MAT::StructPoro(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
  params_(NULL),
  mat_(Teuchos::null),
  porosity_(Teuchos::null),
  dporodt_(Teuchos::null),
  gradporosity_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoro::StructPoro(MAT::PAR::StructPoro* params) :
  params_(params),
  porosity_(Teuchos::null),
  dporodt_(Teuchos::null),
  gradporosity_(Teuchos::null)
{
  mat_ = MAT::Material::Factory(params_->matid_);
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
  size = (int)porosity_->size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
  {
    AddtoPack(data,(*porosity_)[i]);
  }

  // gradporosity_
  size = (int)gradporosity_->size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
  {
    AddtoPack(data,(*gradporosity_)[i]);
  }

  // dporodt_
  size = (int)dporodt_->size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
  {
    AddtoPack(data,(*dporodt_)[i]);
  }

  mat_->Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::Unpack(const vector<char>& data)
{
  std::vector<char>::size_type position = 0;
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
  porosity_=rcp(new std::vector<double >);
  double tmp1 = 0.0;
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp1);
    porosity_->push_back(tmp1);
  }

  // gradporosity_
  gradporosity_=rcp(new std::vector<LINALG::Matrix<3,1> >);
  ExtractfromPack(position,data,size);
  LINALG::Matrix<3,1> tmp2(true);
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp2);
    gradporosity_->push_back(tmp2);
  }

  // dporodt_
  ExtractfromPack(position,data,size);
  dporodt_=rcp(new std::vector<double >);
  double tmp3 = 0.0;
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp3);
    dporodt_->push_back(tmp3);
  }

  // unpack the sub-material
  {
    // we create the sub-material using the factory, because this way we can use the known matid.
    // This way however the material is fully setup according to the dat-file and in order to
    // get a correct material the called unpack method must overwrite all sub-material members.

    if (params_ != NULL)// material are not accessible in postprocessing mode
    {
      mat_ = MAT::Material::Factory(params_->matid_);

      vector<char> basedata(0);
      ExtractfromPack(position,data,basedata);
      mat_->Unpack(basedata);

      if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
    }
  }
}

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputePorosity(double press, double J,
    int gp,double& porosity ,double& dphi_dp,
    double& dphi_dJ,double& dphi_dJdp,double& dphi_dJJ,double& dphi_dpp) const
{

  const double bulkmodulus = params_->bulkmodulus_;
  const double penalty = params_->penaltyparameter_;
  const double initporosity = params_->initporosity_;

  const double a = (bulkmodulus / (1 - initporosity) + press - penalty / initporosity) * J;
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

  double phi = 1 / (2 * a) * (-b + d);

  if (phi >= 1.0 or phi < 0.0)
  {
    dserror("invalid porosity!");
  }

  double d_p = J * (-b+2*penalty)/d;
  double d_p_p = ( d * J + d_p * (b - 2*penalty) ) / (d * d) * J;
  double d_J = a/J * ( -b + 2*penalty ) / d;
  double d_J_p = (d_p / J + ( 1-d_p*d_p/(J*J) ) / d *a);
  double d_J_J = ( a*a/(J*J)-d_J*d_J )/ d;

  //d(porosity) / d(p)
  double tmp1= - J * phi/a + (J+d_p)/(2*a);
  dphi_dp = tmp1;

  //d(porosity) / d(J)
  double tmp2= -phi/J+ 1/(2*J) + d_J / (2*a);
  dphi_dJ = tmp2;

  //d(porosity) / d(J)d(pressure)
  dphi_dJdp= -1/J* tmp1+ d_J_p/(2*a) - d_J*J/(2*a*a);

  //d^2(porosity) / d(J)^2
  dphi_dJJ= phi/(J*J) - tmp2/J - 1/(2*J*J) - d_J/(2*a*J) + d_J_J/(2*a);

  //d^2(porosity) / d(pressure)^2
  dphi_dpp= -J/a* tmp1 + phi*J*J/(a*a) - J/(2*a*a)*(J+d_p) + d_p_p/(2*a);

  porosity= phi;

  //save porosity
  porosity_->at(gp) = phi;

  return;
}

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
double MAT::StructPoro::PorosityAv() const
{
  double porosityav = 0.0;

  std::vector<double>::const_iterator m;
  for (m = porosity_->begin(); m != porosity_->end(); ++m)
  {
    porosityav += *m;
  }
  porosityav = porosityav / (porosity_->size());

  return porosityav;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void MAT::StructPoro::SetPorosityAtGP(std::vector<double> porosity_gp)
{
  //int numgp = porosity_gp.size();

  //porosity_.resize(numgp);

  //set porosity values
  //cout<<"length1: "<<porosity_gp.size()<<endl;
  //cout<<"length2: "<<porosity_->size()<<endl;
  std::vector<double>::iterator m = porosity_gp.begin();
  for (int i = 0; m != porosity_gp.end(); ++m, ++i)
  {
    double porosity = *m;
    porosity_->at(i)=porosity;
  }

  return;
}*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::SetDPoroDtAtGP(std::vector<double> dporodt_gp)
{
  //set dporodt values
  std::vector<double>::iterator m = dporodt_gp.begin();
  for (int i = 0; m != dporodt_gp.end(); ++m, ++i)
  {
    double dporodt = *m;
    dporodt_->at(i) = dporodt;
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::SetGradPorosityAtGP(std::vector<  LINALG::Matrix<3,1> > gradporosity_gp)
{
    std::vector<LINALG::Matrix<3,1> >::iterator m = gradporosity_gp.begin();
    for (int i = 0; m != gradporosity_gp.end(); ++m, ++i)
    {
      LINALG::Matrix<3,1> gradporo  = *m;
      gradporosity_->at(i) = gradporo;
    }
    return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoro::Setup(const int numgp)
{
  porosity_ = rcp(new std::vector< double > (numgp,params_->initporosity_));
  gradporosity_ = rcp(new std::vector< LINALG::Matrix<3,1> > (numgp));
  dporodt_ = rcp(new std::vector<double> (numgp));

  const LINALG::Matrix<3,1> emptyvec(true);
  for (int j=0; j<numgp; ++j)
  {
    gradporosity_->at(j) = emptyvec;
  }
}
