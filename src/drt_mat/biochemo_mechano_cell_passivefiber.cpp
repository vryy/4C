/*!----------------------------------------------------------------------
\file biochemo_mechano_cell_passivefiber.cpp

\brief Implementation of Biochemo-Mechano Coupled passive, viscoelastic material model for the cell.

\level 3

\maintainer Andreas Rauch

*----------------------------------------------------------------------*/

#include "biochemo_mechano_cell_passivefiber.H"

#include "../drt_mat/matpar_bundle.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::BioChemoMechanoCellPassiveFiber::BioChemoMechanoCellPassiveFiber(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  idmatelast_(matdata->GetInt("IDMATELAST")),
  mu_(matdata->GetDouble("VISC")),
  analyticalmaterialtangent_(true)
{
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CellMigrationParams().sublist("STRUCTURAL DYNAMIC"),"MATERIALTANGENT"))
    analyticalmaterialtangent_ = false;
  else
    analyticalmaterialtangent_ = true;
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::BioChemoMechanoCellPassiveFiber::CreateMaterial()
{
  return Teuchos::rcp(new MAT::BioChemoMechanoCellPassiveFiber(this));
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::BioChemoMechanoCellPassiveFiberType MAT::BioChemoMechanoCellPassiveFiberType::instance_;


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::BioChemoMechanoCellPassiveFiberType::Create(
    const std::vector<char> & data )
{
  MAT::BioChemoMechanoCellPassiveFiber* material = new MAT::BioChemoMechanoCellPassiveFiber();
  material->Unpack(data);
  return material;
}


/*----------------------------------------------------------------------*
 |  Constructor                                            rauch  01/17 |
 *----------------------------------------------------------------------*/
MAT::BioChemoMechanoCellPassiveFiber::BioChemoMechanoCellPassiveFiber()
  : params_(NULL),
    isinit_(false)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                                       rauch  01/17 |
 *----------------------------------------------------------------------*/
MAT::BioChemoMechanoCellPassiveFiber::BioChemoMechanoCellPassiveFiber(
    MAT::PAR::BioChemoMechanoCellPassiveFiber* params)
  : params_(params),
    isinit_(false)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                                   rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
    AddtoPack(data,matid);

    // pack history data
    int histsize = -1;
    // if material is not initialized, i.e. start simulation, nothing to pack
    if (!Initialized())
    {
      histsize = 0;
    }
    else
    {
      // if material is initialized (restart): size equates number of gausspoints
      histsize = histdefgrdlast_->size();
    }

    AddtoPack(data,histsize);  // length of history vector(s)

    for (int var=0; var<histsize; ++var)
    {
      // insert history vectors to AddtoPack
      AddtoPack(data,histdefgrdlast_->at(var));
    }

    // Pack data of elastic material
    if (matelast_!=Teuchos::null) {
      matelast_->Pack(data);
    }

  return;

}   // Pack()


/*----------------------------------------------------------------------*
 |  Unpack                                                 rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Unpack(const std::vector<char>& data)
{
  // construct current deformation gradient history variable
  histdefgrdcurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );

  isinit_=true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid = -1;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::BioChemoMechanoCellPassiveFiber*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  int histsize = -1;
  ExtractfromPack(position,data,histsize);

  // if system is not yet initialized, the history vectors have to be initialized
  if (histsize == 0){
    isinit_=false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
    return;
  }

  // construct vector of old deformation gradient matrices
  histdefgrdlast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );

  for (int var=0; var<histsize; ++var)
  {
    // initialize
    LINALG::Matrix<3,3> tmp_matrix3x3(true);

    // matrices of last converged state are unpacked
    ExtractfromPack(position,data,tmp_matrix3x3);
    histdefgrdlast_->push_back(tmp_matrix3x3);

    // initialize current deformation gradient
    histdefgrdcurr_->push_back(tmp_matrix3x3);
  }

  // Unpack data of passive elastic material
  std::vector<char> dataelastic;
  ExtractfromPack(position,data,dataelastic);
  if (dataelastic.size()>0) {
    DRT::ParObject* o = DRT::UTILS::Factory(dataelastic);  // Unpack is done here
    MAT::So3Material* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel==NULL)
      dserror("failed to unpack passive material");
    matelast_ = Teuchos::rcp(matel);
  } else matelast_ = Teuchos::null;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

}   // Unpack()


/*----------------------------------------------------------------------*
 | initialize / allocate internal variables (public)       rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Setup(
    int numgp,
    DRT::INPUT::LineDefinition* linedef)
{
  // construct history variables
  histdefgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdlast_->resize(numgp);

  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdcurr_->resize(numgp);

  LINALG::Matrix<3,3> emptymat3x3(true);
  for (int i=0; i<3; i++)
    emptymat3x3(i,i) = 1.0;

  for (int i=0; i<numgp; i++)
  {
    histdefgrdlast_->at(i) = emptymat3x3;
    histdefgrdcurr_->at(i) = emptymat3x3;
  }

  // Setup of passive material
  matelast_ = Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->idmatelast_));
  matelast_->Setup(numgp, linedef);

  isinit_ = true;
  return;
}   // Setup()


/*----------------------------------------------------------------------*
 |  ResetAll                                               rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::ResetAll(const int numgp)
{
  // construct history variables
  histdefgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdlast_->resize(numgp);
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdcurr_->resize(numgp);

  LINALG::Matrix<3,3> emptymat3x3(true);
  for (int i=0; i<3; i++)
    emptymat3x3(i,i) = 1.0;

  for (int i=0; i<numgp; i++)
  {
    histdefgrdlast_->at(i) = emptymat3x3;
    histdefgrdcurr_->at(i) = emptymat3x3;
  }

  matelast_->ResetAll(numgp);
  isinit_ = false;
  return;
}   // ResetAll()


/*----------------------------------------------------------------------*
 |  Update internal variables                              rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Update()
{
  // make current values at time step t^{n+1} to old values at time  t^n
  histdefgrdlast_ = histdefgrdcurr_;
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdcurr_->resize(histdefgrdlast_->size());

  matelast_->Update();
  return;
}   // Update()


/*----------------------------------------------------------------------*
 |  Reset internal variables                               rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::ResetStep()
{
  matelast_->ResetStep();
  return;
}  //ResetStep()


/*----------------------------------------------------------------------*
 |  Evaluate Material                                      rauch  01/17 |
 *--------------------------------------------------------------------- */
void MAT::BioChemoMechanoCellPassiveFiber::Evaluate(
                           const LINALG::Matrix<3,3>* defgrd,
                           const LINALG::Matrix<6,1>* glstrain,
                           Teuchos::ParameterList& params,
                           LINALG::Matrix<6,1>* stress,
                           LINALG::Matrix<6,6>* cmat,
                           const int eleGID)
{

  // get gauss point number
  const int gp = params.get<int>("gp",-1);
  if (gp == -1)   dserror("no Gauss point number provided in material");
  // get time algorithmic parameters
  double dt = params.get<double>("delta time",-1.0);
  if (dt == -1.0) dserror("no time step size provided in material");

  // copy deformation gradient in histroy variable
  histdefgrdcurr_->at(gp) = *defgrd;

  // setup inverse of deformation gradient
  LINALG::Matrix<3,3> invdefgrd(*defgrd);
  invdefgrd.Invert();
  // setup deformation gradient rate, rotation tensor, strain rate and rotation rate
  // \dot{F} = \frac {F^n - F^{n-1}} {\Delta t}
  LINALG::Matrix<3,3> defgrdrate(true);
  // R = F * U^{-1}
  LINALG::Matrix<3,3> R(true);
  // \dot{\epsilon} = d = 0.5 * (\dot{F}F^{-1} + (\dot{F}F^{-1})^{T}
  LINALG::Matrix<6,1> strainrate(true);
  // \dot{R} = \frac {R^n - R^{n-1}} {\Delta t}
  LINALG::Matrix<3,3> rotationrate(true);
  // calc the rates
  SetupRates(*defgrd,invdefgrd,params,defgrdrate,R,strainrate,rotationrate,gp,dt);

  // calc the viscosity tensor
  const double viscosity = params_->mu_;
  // temporary viscous stress tensor
  LINALG::Matrix<6,1> visc_stress_cauchy;
  // calc viscous stress contribution
  visc_stress_cauchy(0) = strainrate(0);
  visc_stress_cauchy(1) = strainrate(1);
  visc_stress_cauchy(2) = strainrate(2);
  visc_stress_cauchy(3) = strainrate(3);
  visc_stress_cauchy(4) = strainrate(4);
  visc_stress_cauchy(5) = strainrate(5);
  visc_stress_cauchy.Scale(2.0*viscosity);

  // Transform Cauchy stress to PK2 stress
  // S = J * F^{-1} \sigma F^{-T}
  LINALG::Matrix<NUM_STRESS_3D,1> visc_stress_PK2(true); //6x1
  LINALG::Matrix<3,3> visc_stress_cauchy_mat(true); //3x3
  CauchytoPK2(visc_stress_PK2,visc_stress_cauchy_mat,*defgrd,invdefgrd,visc_stress_cauchy);

  // temporary viscous material tangent
  LINALG::Matrix<6,6> visc_cmat;

  // Evaluate elastic stress and material tangent
  matelast_->Evaluate(defgrd,glstrain,params,stress,cmat,eleGID);

  // sum up total stress
stress->Update(1.0,visc_stress_PK2,1.0);

  return;
} // MAT::BioChemoMechanoCellPassiveFiber::Evaluate


/*----------------------------------------------------------------------------------*
 | Calculation of deformation gradient rate, rotation tensor, strain rate and       |
 | rotation rate (finite difference scheme)                            rauch  01/17 |
 *----------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::SetupRates(
    LINALG::Matrix<3,3> defgrd,
    LINALG::Matrix<3,3> invdefgrd,
    Teuchos::ParameterList& params,
    LINALG::Matrix<3,3>& defgrdrate,
    LINALG::Matrix<3,3>& R,
    LINALG::Matrix<6,1>& strainrate,
    LINALG::Matrix<3,3>& rotationrate,
    const int& gp,
    const double& dt)
{
  // Read history
  LINALG::Matrix<3,3> defgrdlast;
  defgrdlast = histdefgrdlast_->at(gp);

  // Rate of deformation gradient: \dot{F} = \frac {F^{n+1} - F^{n}} {\Delta t}
  defgrdrate.Update(1.0,defgrd,0.0);
  defgrdrate.Update(-1.0,defgrdlast,1.0);
  defgrdrate.Scale(1.0/dt);

  // Calculate velocity gradient l = \dot{F}F^{-1}
  LINALG::Matrix<3,3> velgradient(true);
  velgradient.MultiplyNN(defgrdrate,invdefgrd);

  // Rate of strain/symmetric part of velocity gradient
  // d = 0.5 * (l + l^{T}) = 0.5 * (\dot{F}F^{-1} + (\dot{F}F^{-1})^{T})
  // Remark: strain-like 6-Voigt vector
  strainrate(0) = velgradient(0,0) + velgradient(0,0);
  strainrate(1) = velgradient(1,1) + velgradient(1,1);
  strainrate(2) = velgradient(2,2) + velgradient(2,2);
  strainrate(3) = velgradient(0,1) + velgradient(1,0);
  strainrate(4) = velgradient(1,2) + velgradient(2,1);
  strainrate(5) = velgradient(0,2) + velgradient(2,0);
  strainrate.Scale(0.5);

}  // SetupRates


/*----------------------------------------------------------------------*
 | pull back of spatial stresses                           rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::CauchytoPK2(
  LINALG::Matrix<6,1>& Sactive,
  LINALG::Matrix<3,3>& cauchystress,
  LINALG::Matrix<3,3> defgrd,
  LINALG::Matrix<3,3> invdefgrd,
  LINALG::Matrix<6,1> sigma)
{
  // calculate the Jacobi-determinant
  double detF = defgrd.Determinant();   // const???

  // Convert stress like 6x1-Voigt vector to 3x3 matrix
  cauchystress(0,0) = sigma(0);
  cauchystress(0,1) = sigma(3);
  cauchystress(0,2) = sigma(5);
  cauchystress(1,0) = cauchystress(0,1);
  cauchystress(1,1) = sigma(1);
  cauchystress(1,2) = sigma(4);
  cauchystress(2,0) = cauchystress(0,2);
  cauchystress(2,1) = cauchystress(1,2);
  cauchystress(2,2) = sigma(2);

  // S = J * F^{-1} * sigma * F^{-T}
  LINALG::Matrix<3,3> temp(true);
  LINALG::Matrix<3,3> S(true);
  temp.MultiplyNN(invdefgrd,cauchystress);
  S.MultiplyNT(temp,invdefgrd);
  S.Scale(detF);

  // Sactive is stress like 6x1-Voigt vector
  Sactive(0) = S(0,0);
  Sactive(1) = S(1,1);
  Sactive(2) = S(2,2);
  Sactive(3) = S(0,1);
  Sactive(4) = S(1,2);
  Sactive(5) = S(0,2);

}  // CauchytoPK2()


/*----------------------------------------------------------------------*
 |  Names of gp data to be visualized                      rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::VisNames(std::map<std::string,int>& names)
{
  matelast_->VisNames(names);
  return;
} // VisNames()


/*----------------------------------------------------------------------*
 |  gp data to be visualized                               rauch  01/17 |
 *----------------------------------------------------------------------*/
bool MAT::BioChemoMechanoCellPassiveFiber::VisData(
    const std::string& name,
    std::vector<double>& data,
    int numgp ,
    int eleID)
{
  return matelast_->VisData(name, data, numgp, eleID);
}
