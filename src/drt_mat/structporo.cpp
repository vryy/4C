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
  surfporosity_(Teuchos::null),
  dporodt_(Teuchos::null),
  gradporosity_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::StructPoro::StructPoro(MAT::PAR::StructPoro* params) :
  params_(params),
  porosity_(Teuchos::null),
  surfporosity_(Teuchos::null),
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

  // surfporosity_
  size = (int) surfporosity_->size();
  AddtoPack(data,size);

  // iterator
  std::map<int,std::vector<double> >::const_iterator iter;

  int i=0;
  for(iter=surfporosity_->begin();iter!=surfporosity_->end();++iter)
  {
    AddtoPack(data,iter->first);
    AddtoPack(data,iter->second);
    ++i;
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

  // gradJ_
  size = (int)gradJ_->size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
  {
    AddtoPack(data,(*gradJ_)[i]);
  }

  mat_->Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::Unpack(const std::vector<char>& data)
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
  porosity_=Teuchos::rcp(new std::vector<double >);
  double tmp = 0.0;
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp);
    porosity_->push_back(tmp);
  }

  ExtractfromPack(position,data,size);

  surfporosity_ = Teuchos::rcp(new std::map<int, std::vector< double > >);

  for(int i=0;i<size;i++)
  {
    int dof;
    std::vector<double > value;
    ExtractfromPack(position,data,dof);
    ExtractfromPack(position,data,value);

    //add to map
    surfporosity_->insert(std::pair<int,std::vector<double > >(dof,value));

  }

  // gradporosity_
  gradporosity_=Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> >);
  ExtractfromPack(position,data,size);
  LINALG::Matrix<3,1> tmp2(true);
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp2);
    gradporosity_->push_back(tmp2);
  }

  // dporodt_
  ExtractfromPack(position,data,size);
  dporodt_=Teuchos::rcp(new std::vector<double >);
  tmp = 0.0;
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp);
    dporodt_->push_back(tmp);
  }

  // gradporosity_
  gradJ_=Teuchos::rcp(new std::vector<LINALG::Matrix<1,3> >);
  ExtractfromPack(position,data,size);
  LINALG::Matrix<1,3> tmp3(true);
  for (int i=0; i<size; ++i)
  {
    ExtractfromPack(position,data,tmp3);
    gradJ_->push_back(tmp3);
  }

  // unpack the sub-material
  {
    // we create the sub-material using the factory, because this way we can use the known matid.
    // This way however the material is fully setup according to the dat-file and in order to
    // get a correct material the called unpack method must overwrite all sub-material members.

    if (params_ != NULL)// material are not accessible in postprocessing mode
    {
      mat_ = MAT::Material::Factory(params_->matid_);

      std::vector<char> basedata(0);
      ExtractfromPack(position,data,basedata);
      mat_->Unpack(basedata);

      if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputePorosity( double press,
                                       double J,
                                       int gp,
                                       double& porosity,
                                       double& dphi_dp,
                                       double& dphi_dJ,
                                       double& dphi_dJdp,
                                       double& dphi_dJJ,
                                       double& dphi_dpp) const
{

  const double & bulkmodulus  = params_->bulkmodulus_;
  const double & penalty      = params_->penaltyparameter_;
  const double & initporosity = params_->initporosity_;

  const double a = (bulkmodulus / (1 - initporosity) + press - penalty / initporosity) * J;
  const double b = -a + bulkmodulus + penalty;
  const double c = (b / a) * (b / a) + 4.0 * penalty / a;
  double d = sqrt(c) * a;
  double sign = 1.0;

  double test = 1 / (2.0 * a) * (-b + d);
  if (test >= 1.0 or test < 0.0)
  {
    sign = -1.0;
    d = sign * d;
  }

  const double phi = 1 / (2 * a) * (-b + d);

  if (phi >= 1.0 or phi < 0.0)
    dserror("invalid porosity: %f", porosity);

  const double d_p = J * (-b+2.0*penalty)/d;
  const double d_p_p = ( d * J + d_p * (b - 2.0*penalty) ) / (d * d) * J;
  const double d_J = a/J * ( -b + 2.0*penalty ) / d;
  const double d_J_p = (d_p / J + ( 1-d_p*d_p/(J*J) ) / d *a);
  const double d_J_J = ( a*a/(J*J)-d_J*d_J )/ d;

  //d(porosity) / d(p)
  double tmp1= - J * phi/a + (J+d_p)/(2.0*a);
  dphi_dp = tmp1;

  //d(porosity) / d(J)
  double tmp2= -phi/J+ 1/(2*J) + d_J / (2.0*a);
  dphi_dJ = tmp2;

  //d(porosity) / d(J)d(pressure)
  dphi_dJdp= -1/J* tmp1+ d_J_p/(2*a) - d_J*J/(2.0*a*a);

  //d^2(porosity) / d(J)^2
  dphi_dJJ= phi/(J*J) - tmp2/J - 1/(2.0*J*J) - d_J/(2*a*J) + d_J_J/(2.0*a);

  //d^2(porosity) / d(pressure)^2
  dphi_dpp= -J/a* tmp1 + phi*J*J/(a*a) - J/(2.0*a*a)*(J+d_p) + d_p_p/(2.0*a);

  /*
  dphi_dp = 0.0;
  dphi_dJ = 0.0;
  dphi_dJdp = 0.0;
  dphi_dJJ = 0.0;
  dphi_dpp = 0.0;
  phi = 0.1;
*/
  /*
  double phi = initporosity + (1.0-initporosity) * (1.0-initporosity) / bulkmodulus * press
                + (1.0-initporosity) * (J-1.0);

  dphi_dp = (1.0-initporosity) * (1.0-initporosity) / bulkmodulus;
  dphi_dJ = (1.0-initporosity);
  dphi_dJdp = 0.0;
  dphi_dJJ = 0.0;
  dphi_dpp = 0.0;
  */

  porosity= phi;

  //save porosity
  porosity_->at(gp) = phi;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputeSurfPorosity( double     press,
                                           double     J,
                                           const int  surfnum,
                                           int        gp,
                                           double&    porosity,
                                           double&    dphi_dp,
                                           double&    dphi_dJ,
                                           double&    dphi_dJdp,
                                           double&    dphi_dJJ,
                                           double&    dphi_dpp) const
{
  double phi= 0.0;
  ComputePorosity(press,
                  J,
                  gp,
                  phi,
                  dphi_dp,
                  dphi_dJ,
                  dphi_dJdp,
                  dphi_dJJ,
                  dphi_dpp);

  if(gp==0)
   ( (*surfporosity_)[surfnum] ).clear();

  ( (*surfporosity_)[surfnum] ).push_back(phi);

  porosity = phi;

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
void MAT::StructPoro::PorosityGradientAv(LINALG::Matrix<3,1>& porosityav) const
{
  std::vector<LINALG::Matrix<3,1> >::const_iterator m;
  for (m = gradporosity_->begin(); m != gradporosity_->end(); ++m)
  {
    porosityav.Update(1.0, *m ,1.0) ;
  }
  double size = gradporosity_->size();

  porosityav.Scale(1.0/size);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::PorosityGradientAv(LINALG::Matrix<2,1>& porosityav) const
{
  dserror("not implemented!");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
LINALG::Matrix<1,3> MAT::StructPoro::JGradientAv() const
{
  LINALG::Matrix<1,3> porositygradientav(true);

  std::vector<LINALG::Matrix<1,3> >::const_iterator m;
  for (m = gradJ_->begin(); m != gradJ_->end(); ++m)
  {
    porositygradientav.Update(1.0, *m ,1.0) ;
  }
  double size = gradJ_->size();
  //porositygradientav.Scale(1/(gradporosity_->size()));
  porositygradientav.Scale(1.0/size);

  return porositygradientav;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::StructPoro::DPorosityDtAv() const
{
  double av = 0.0;

  std::vector<double>::const_iterator m;
  for (m = dporodt_->begin(); m != dporodt_->end(); ++m)
  {
    av += *m;
  }
  av = av / (dporodt_->size());

  return av;
}

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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::SetGradPorosityAtGP(std::vector<  LINALG::Matrix<2,1> > gradporosity_gp)
{
 dserror("not implemented");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::SetGradJAtGP(std::vector<  LINALG::Matrix<1,3> > gradporosity_gp)
{
  std::vector<LINALG::Matrix<1,3> >::iterator m = gradporosity_gp.begin();
  for (int i = 0; m != gradporosity_gp.end(); ++m, ++i)
  {
    LINALG::Matrix<1,3> gradporo  = *m;
    gradJ_->at(i) = gradporo;
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::GetGradPorosityAtGP(LINALG::Matrix<3,1>& gradporosity, int gp) const
{
  gradporosity = gradporosity_->at(gp);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::StructPoro::GetGradPorosityAtGP(LINALG::Matrix<2,1>& gradporosity, int gp) const
{
  dserror("not implemented");
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::StructPoro::Setup(const int numgp)
{
  porosity_ = Teuchos::rcp(new std::vector< double > (numgp,params_->initporosity_));
  //surfporosity_ = Teuchos::rcp(new std::vector< double > (numgp,-1.0));
  surfporosity_ = Teuchos::rcp(new std::map<int, std::vector< double > >);
  gradporosity_ = Teuchos::rcp(new std::vector< LINALG::Matrix<3,1> > (numgp));
  dporodt_ = Teuchos::rcp(new std::vector<double> (numgp));
  gradJ_ = Teuchos::rcp(new std::vector< LINALG::Matrix<1,3> > (numgp));

  const LINALG::Matrix<3,1> emptyvec(true);
  const LINALG::Matrix<1,3> emptyvecT(true);
  for (int j=0; j<numgp; ++j)
  {
    gradporosity_->at(j) = emptyvec;
    gradJ_->at(j) = emptyvecT;
  }
}

