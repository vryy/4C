/*!-----------------------------------------------------------------------*
 \file structporo.cpp

 \brief wrapper for structure material of porous media

 <pre>
  \level 2

   \maintainer Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *-----------------------------------------------------------------------*/

#include "structporo.H"
#include "poro_law.H"

#include <vector>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_factory.H"  // for function Factory in Unpack

#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/so3_material.H"


/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
MAT::PAR::StructPoro::StructPoro(Teuchos::RCP<MAT::PAR::Material> matdata) :
  Parameter(matdata),
  matid_(matdata->GetInt("MATID")),
  porolawID_(matdata->GetInt("POROLAWID")),
  initporosity_(matdata->GetDouble("INITPOROSITY"))
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
  Teuchos::RCP<MAT::PAR::Material> curmat = DRT::Problem::Instance(probinst)->Materials()->ById(porolawID_);

  switch (curmat->Type())
  {
  case INPAR::MAT::m_poro_law_linear:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::PoroLawLinear(curmat));
    porolaw_ = static_cast<MAT::PAR::PoroLaw*>(curmat->Parameter());
    break;
  }
  case INPAR::MAT::m_poro_law_constant:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::PoroLawConstant(curmat));
    porolaw_ = static_cast<MAT::PAR::PoroLaw*>(curmat->Parameter());
    break;
  }
  case INPAR::MAT::m_poro_law_logNeoHooke_Penalty:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PoroLawNeoHooke(curmat));
      porolaw_ = static_cast<MAT::PAR::PoroLaw*>(curmat->Parameter());
      break;
    }
  case INPAR::MAT::m_poro_law_incompr_skeleton:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::PoroLawIncompSkeleton(curmat));
    porolaw_ = static_cast<MAT::PAR::PoroLawIncompSkeleton*>(curmat->Parameter());
    break;
  }
  case INPAR::MAT::m_poro_law_linear_biot:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PoroLawLinBiot(curmat));
      porolaw_ = static_cast<MAT::PAR::PoroLawLinBiot*>(curmat->Parameter());
      break;
    }
  case INPAR::MAT::m_poro_law_density_dependent:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PoroLawDensityDependent(curmat));
      porolaw_ = static_cast<MAT::PAR::PoroLawDensityDependent*>(curmat->Parameter());
      break;
    }
  default:
    dserror("invalid material for porosity law %d", curmat->Type());
    break;
  }
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::StructPoro::CreateMaterial()
{
  return Teuchos::rcp(new MAT::StructPoro(this));
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
MAT::StructPoroType MAT::StructPoroType::instance_;

DRT::ParObject* MAT::StructPoroType::Create(const std::vector<char> & data)
{
  MAT::StructPoro* struct_poro = new MAT::StructPoro();
  struct_poro->Unpack(data);
  return struct_poro;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
MAT::StructPoro::StructPoro() :
  params_(NULL),
  mat_(Teuchos::null),
  porosity_(Teuchos::null),
  surfporosity_(Teuchos::null),
  isinitialized_(false)
{
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
MAT::StructPoro::StructPoro(MAT::PAR::StructPoro* params) :
  params_(params),
  porosity_(Teuchos::null),
  surfporosity_(Teuchos::null),
  isinitialized_(false)
{
  mat_ = Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->matid_));
  if (mat_ == Teuchos::null) dserror("MAT::StructPoro: underlying material should be of type MAT::So3Material");
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::PoroSetup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  porosity_ = Teuchos::rcp(new std::vector< double > (numgp,params_->initporosity_));
  surfporosity_ = Teuchos::rcp(new std::map<int, std::vector< double > >);

  isinitialized_=true;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
inline INPAR::MAT::MaterialType MAT::StructPoro::PoroLawType() const
{
  return params_->porolaw_->Type();
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
double MAT::StructPoro::InvBulkmodulus() const
{
  return params_->porolaw_->InvBulkmodulus();
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
double MAT::StructPoro::Density() const
{
  if(params_->initporosity_==1.0)
    return mat_->Density();
  else
    return ( (1.0-params_->initporosity_) * mat_->Density());
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::Pack(DRT::PackBuffer& data) const
{
  if(not isinitialized_)
    dserror("poro material not initialized. Not a poro element?");

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

  // surfporosity_ (i think it is not necessary to pack/unpack this...)
  size = (int) surfporosity_->size();
  AddtoPack(data,size);
  // iterator
  std::map<int,std::vector<double> >::const_iterator iter;
  for(iter=surfporosity_->begin();iter!=surfporosity_->end();++iter)
  {
    AddtoPack(data,iter->first);
    AddtoPack(data,iter->second);
  }

  // Pack data of underlying material
  if (mat_!=Teuchos::null)
    mat_->Pack(data);
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
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

  // surface porosity (i think it is not necessary to pack/unpack this...)
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

  // Unpack data of sub material (these lines are copied from drt_element.cpp)
  std::vector<char> datamat;
  ExtractfromPack(position,data,datamat);
  if (datamat.size()>0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(datamat);  // Unpack is done here
    MAT::So3Material* mat = dynamic_cast<MAT::So3Material*>(o);
    if (mat==NULL)
      dserror("failed to unpack elastic material");
    mat_ = Teuchos::rcp(mat);
  }
  else mat_ = Teuchos::null;

  isinitialized_=true;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputePorosity( const double& refporosity,
                                       const double& press,
                                       const double& J,
                                       const int& gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       double* dphi_dphiref,
                                       bool save)
{
  if(refporosity == 1.0)
  {
    //this is pure fluid. The porosity does not change

    porosity = refporosity;
    if(dphi_dp)      *dphi_dp      = 0.0;
    if(dphi_dJ)      *dphi_dJ      = 0.0;
    if(dphi_dJdp)    *dphi_dJdp    = 0.0;
    if(dphi_dJJ)     *dphi_dJJ     = 0.0;
    if(dphi_dpp)     *dphi_dpp     = 0.0;
    if(dphi_dphiref) *dphi_dphiref = 0.0;
    return;
  }

  params_->porolaw_->ComputePorosity(
      refporosity,
      press,
      J,
      gp,
      porosity,
      dphi_dp,
      dphi_dJ,
      dphi_dJdp,
      dphi_dJJ,
      dphi_dpp,
      dphi_dphiref);

  //save porosity
  if(save)
    porosity_->at(gp) = porosity;

  return;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputePorosity( Teuchos::ParameterList& params,
                                       double press,
                                       double J,
                                       int gp,
                                       double& porosity,
                                       double* dphi_dp,
                                       double* dphi_dJ,
                                       double* dphi_dJdp,
                                       double* dphi_dJJ,
                                       double* dphi_dpp,
                                       bool save)
{

  ComputePorosity( params_->initporosity_, //reference porosity equals initial porosity for non reactive material
                   press,
                   J,
                   gp,
                   porosity,
                   dphi_dp,
                   dphi_dJ,
                   dphi_dJdp,
                   dphi_dJJ,
                   dphi_dpp,
                   NULL, //reference porosity is constant (non reactive) -> derivative not needed
                   save);

  return;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputePorosity( Teuchos::ParameterList& params,
                                       double press,
                                       double J,
                                       int gp,
                                       double& porosity,
                                       bool save)
{

  ComputePorosity( params,
                   press,
                   J,
                   gp,
                   porosity,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   save);

  return;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputeSurfPorosity( Teuchos::ParameterList& params,
                                           double     press,
                                           double     J,
                                           const int  surfnum,
                                           int        gp,
                                           double&    porosity,
                                           double*    dphi_dp,
                                           double*    dphi_dJ,
                                           double*    dphi_dJdp,
                                           double*    dphi_dJJ,
                                           double*    dphi_dpp,
                                           bool save)
{
  ComputePorosity(params,
                  press,
                  J,
                  gp,
                  porosity,
                  dphi_dp,
                  dphi_dJ,
                  dphi_dJdp,
                  dphi_dJJ,
                  dphi_dpp,
                  save);

  if(save)
  {
    if(gp==0)  //it's a new iteration, so old values are not needed any more
     ( (*surfporosity_)[surfnum] ).clear();

    ( (*surfporosity_)[surfnum] ).push_back(porosity);
  }

  return;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::ComputeSurfPorosity( Teuchos::ParameterList& params,
                                           double     press,
                                           double     J,
                                           const int  surfnum,
                                           int        gp,
                                           double&    porosity,
                                           bool save)
{

  ComputeSurfPorosity( params,
                       press,
                       J,
                       surfnum,
                       gp,
                       porosity,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       NULL,
                       save);

  return;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::CouplStress(  const LINALG::Matrix<3,3>& defgrd,
                                    const LINALG::Matrix<3,1>& fluidvel,
                                    const double& press,
                                    LINALG::Matrix<6,1>& couplstress
                                    ) const
{
  const double J = defgrd.Determinant();

  // Right Cauchy-Green tensor = F^T * F
  LINALG::Matrix<3,3> cauchygreen;
  cauchygreen.MultiplyTN(defgrd,defgrd);

  // inverse Right Cauchy-Green tensor
  LINALG::Matrix<3,3> C_inv;
  C_inv.Invert(cauchygreen);

  //inverse Right Cauchy-Green tensor as vector
  LINALG::Matrix<6,1> C_inv_vec;
  for(int i =0, k=0;i<3; i++)
    for(int j =0;j<3-i; j++,k++)
      C_inv_vec(k)=C_inv(i+j,j);

  for(int i=0; i<6 ; i++)
    couplstress(i)= -1.0*J*press*C_inv_vec(i);
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::CouplStress(  const LINALG::Matrix<2,2>& defgrd,
                                    const LINALG::Matrix<2,1>& fluidvel,
                                    const double& press,
                                    LINALG::Matrix<3,1>& couplstress
                                    ) const
{
  const double J = defgrd.Determinant();

  // Right Cauchy-Green tensor = F^T * F
  LINALG::Matrix<2,2> cauchygreen;
  cauchygreen.MultiplyTN(defgrd,defgrd);

  // inverse Right Cauchy-Green tensor
  LINALG::Matrix<2,2> C_inv;
  C_inv.Invert(cauchygreen);

  //inverse Right Cauchy-Green tensor as vector
  LINALG::Matrix<3,1> C_inv_vec;
  for(int i =0, k=0;i<2; i++)
    for(int j =0;j<2-i; j++,k++)
      C_inv_vec(k)=C_inv(i+j,j);

  for(int i=0; i<3 ; i++)
    couplstress(i)= -1.0*J*press*C_inv_vec(i);
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::ConstitutiveDerivatives(Teuchos::ParameterList& params,
                                              double     press,
                                              double     J,
                                              double     porosity,
                                              double*    dW_dp,
                                              double*    dW_dphi,
                                              double*    dW_dJ,
                                              double*    dW_dphiref,
                                              double*    W)
{
  if(porosity == 0.0)
    dserror("porosity equals zero!! Wrong initial porosity?");

  ConstitutiveDerivatives(params,
                         press,
                         J,
                         porosity,
                         params_->initporosity_,
                         dW_dp,
                         dW_dphi,
                         dW_dJ,
                         dW_dphiref,
                         W);

  return;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::ConstitutiveDerivatives(Teuchos::ParameterList& params,
                                              double     press,
                                              double     J,
                                              double     porosity,
                                              double     refporosity,
                                              double*    dW_dp,
                                              double*    dW_dphi,
                                              double*    dW_dJ,
                                              double*    dW_dphiref,
                                              double*    W)
{
  params_->porolaw_->ConstitutiveDerivatives(
    params,
    press,
    J,
    porosity,
    refporosity,
    dW_dp,
    dW_dphi,
    dW_dJ,
    dW_dphiref,
    W);

  return;
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
void MAT::StructPoro::VisNames(std::map<std::string,int>& names)
{
  mat_->VisNames(names);
  std::string porosity = "porosity";
  names[porosity] = 1; // scalar
}

/*----------------------------------------------------------------------*
                                                              vuong 06/11|
*----------------------------------------------------------------------*/
bool MAT::StructPoro::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (mat_->VisData(name,data,numgp))
    return true;
  if (name=="porosity")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    data[0] = PorosityAv();
    return true;
  }
  return false;
}
