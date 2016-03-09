/*----------------------------------------------------------------------*/
/*!
\file scalardepinterp.cpp
\brief
This file contains a framework, which allows for the linear interpolation between two different
strain energy functions, where the interpolation ratio come from 'outside'. For now this is supposed to
come from the idea, that there is a scalar quantity (e.g. foam cells) which induces growth. The grown volume
fraction thereby consists of another material (e.g. softer or stiffer material parameters or totally different
strain energy function).

The input line should read
MAT 1 MAT_ScalarDepInterp CONC_0_MAT 2 CONC_INFTY_MAT 3

<pre>
Maintainer: Moritz Thon
            thon@mhpc.mw.tum.de
            089 - 289-10364
</pre>
*/

/*----------------------------------------------------------------------*/

#include "scalardepinterp.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_utils_factory.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScalarDepInterp::ScalarDepInterp(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  id_zero_conc_mat_(matdata->GetInt("IDZEROCONCMAT")),
  id_infty_conc_mat_(matdata->GetInt("IDINFTYCONCMAT"))
{

};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ScalarDepInterp::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScalarDepInterp(this));
}


MAT::ScalarDepInterpType MAT::ScalarDepInterpType::instance_;


DRT::ParObject* MAT::ScalarDepInterpType::Create( const std::vector<char> & data )
{
  MAT::ScalarDepInterp* ScalarDepInterp = new MAT::ScalarDepInterp();
  ScalarDepInterp->Unpack(data);
  return ScalarDepInterp;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScalarDepInterp::ScalarDepInterp()
  : params_(NULL),
    isinit_(false),
    zero_conc_mat_(Teuchos::null),
    infty_conc_mat_(Teuchos::null),
    zero_conc_ratio_(Teuchos::null)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScalarDepInterp::ScalarDepInterp(MAT::PAR::ScalarDepInterp* params)
  : params_(params),
    isinit_(false),
    zero_conc_mat_(Teuchos::null),
    infty_conc_mat_(Teuchos::null),
    zero_conc_ratio_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  if (isinit_)
    dserror("This function should just be called, if the material is jet not initialized.");

  // Setup of elastic material for zero concentration
  zero_conc_mat_ = Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->id_zero_conc_mat_));
  zero_conc_mat_->Setup(numgp, linedef);

  // Setup of elastic material for zero concentration
  infty_conc_mat_ = Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->id_infty_conc_mat_));
  infty_conc_mat_->Setup(numgp, linedef);

  //Some safety check
  const double density1 = zero_conc_mat_->Density();
  const double density2 = infty_conc_mat_->Density();
  if (abs(density1-density2)>1e-14)
    dserror("The densities of the materials specified in IDZEROCONCMAT and IDINFTYCONCMAT must be equal!");

  zero_conc_ratio_ = std::vector<double >(numgp,1.0);

  //initialization done
  isinit_ = true;

  return;
  }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Evaluate(const LINALG::Matrix<3,3>* defgrd,
                               const LINALG::Matrix<6,1>* glstrain,
                               Teuchos::ParameterList& params,
                               LINALG::Matrix<6,1>* stress,
                               LINALG::Matrix<6,6>* cmat,
                               const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1)
    dserror("no Gauss point number provided in material");

  //evaluate elastic material corresponding to zero concentration
  LINALG::Matrix<6,1> stress_zero_conc = *stress;
  LINALG::Matrix<6,6> cmat_zero_conc = *cmat;

  zero_conc_mat_->Evaluate(defgrd, glstrain, params, &stress_zero_conc, &cmat_zero_conc,eleGID);

  //evaluate elastic material corresponding to infinite concentration
  LINALG::Matrix<6,1> stress_infty_conc = *stress;
  LINALG::Matrix<6,6> cmat_infty_conc = *cmat;

  infty_conc_mat_->Evaluate(defgrd, glstrain, params, &stress_infty_conc, &cmat_infty_conc,eleGID);

  double conc_zero_ratio;
  //get the ratio of interpolation
  // NOTE: if no ratio is available, we use the conc_zero material!
  if (params.isParameter("conc_zero_ratio"))
  {
    conc_zero_ratio = params.get< double >("conc_zero_ratio");

    // NOTE: this would be nice, but since negative concentrations can occur,
    // we have to catch 'unnatural' cases different...
    //  if ( conc_zero_ratio < -1.0e-14 or conc_zero_ratio > (1.0+1.0e-14) )
    //      dserror("The conc_zero_ratio must be in [0,1]!");

    // e.g. like that:
    if ( conc_zero_ratio < -1.0e-14 )
      conc_zero_ratio = 0.0;
    if ( conc_zero_ratio > (1.0+1.0e-14) )
      conc_zero_ratio = 1.0;

    zero_conc_ratio_.at(gp)=conc_zero_ratio;
  }
  else
  {
    conc_zero_ratio=zero_conc_ratio_.at(gp);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // \mym S = 2 \frac{\partial}{\partial \mym C} \left( \Psi(\mym C) \right) = ...
  // \mym S = 2 \frac{\partial}{\partial \mym C} \left( \gamma(J) * \Psi_1(\mym C) + (1-\gamma(J)) * \Psi_2(\mym C) \right) = ...
  ///////////////////////////////////////////////////////////////////////////////////////////////

  //do the linear interpolation between stresses:
  // ... = \gamma * 2 * \frac{\partial}{\partial \mym C} \Psi_1 + (1-\gamma) * 2* \frac{\partial}{\partial \mym C} \Psi_2)
  stress->Update(conc_zero_ratio,stress_zero_conc,1-conc_zero_ratio,stress_infty_conc,0.0);

  cmat->Update(conc_zero_ratio,cmat_zero_conc,1-conc_zero_ratio,cmat_infty_conc,0.0);

  if (params.isParameter("dconc_zero_ratio_dC"))
  {
    //get derivative of interpolation ratio w.r.t. glstrain
    Teuchos::RCP<LINALG::Matrix<6,1> > dconc_zero_ratio_dC =
        params.get< Teuchos::RCP<LINALG::Matrix<6,1> > >( "dconc_zero_ratio_dC");

    //evaluate strain energy functions
    double psi_zero_conc = 0.0;
    zero_conc_mat_->StrainEnergy(*glstrain,psi_zero_conc,eleGID);

    double psi_infty_conc = 0.0;
    infty_conc_mat_->StrainEnergy(*glstrain,psi_infty_conc,eleGID);

    //...and add the stresses due to possible dependency of the ratio w.r.t. to C
    // ... + * 2 * \Psi_1 * \frac{\partial}{\partial \mym C} \gamma - 2 * \Psi_2 * \frac{\partial}{\partial \mym C} \gamma )
    stress->Update(2.0*psi_zero_conc,*dconc_zero_ratio_dC,-2.0*psi_infty_conc,*dconc_zero_ratio_dC,1.0);

    //Note: for the linearization we do neglect the derivatives of the ratio w.r.t. glstrain
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Pack(DRT::PackBuffer& data) const
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

  int numgp=0;
  if (isinit_)
  {
    numgp = zero_conc_ratio_.size();;   // size is number of gausspoints
  }
  AddtoPack(data,numgp);

  for (int gp=0; gp<numgp; gp++)
  {
    AddtoPack(data,zero_conc_ratio_.at(gp));
  }

  // Pack data of both elastic materials
  if (zero_conc_mat_!=Teuchos::null and infty_conc_mat_!=Teuchos::null)
  {
    zero_conc_mat_->Pack(data);
    infty_conc_mat_->Pack(data);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Unpack(const std::vector<char>& data)
{
  isinit_=true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::ScalarDepInterp*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  int numgp;
  ExtractfromPack(position,data,numgp);
  if (not (numgp == 0))
  {
    zero_conc_ratio_ = std::vector<double>(numgp,1.0);

    for (int gp=0; gp<numgp; gp++)
    {
      double zero_conc_ratio=1.0;
      ExtractfromPack(position,data,zero_conc_ratio);
      zero_conc_ratio_.at(gp)=zero_conc_ratio;
    }
  }

  // Unpack data of elastic material (these lines are copied from drt_element.cpp)
  std::vector<char> dataelastic;
  ExtractfromPack(position,data,dataelastic);
  if (dataelastic.size()>0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(dataelastic);  // Unpack is done here
    MAT::So3Material* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel==NULL)
      dserror("failed to unpack elastic material");
    zero_conc_mat_ = Teuchos::rcp(matel);
  }
  else
    zero_conc_mat_ = Teuchos::null;

  // Unpack data of elastic material (these lines are copied from drt_element.cpp)
  std::vector<char> dataelastic2;
  ExtractfromPack(position,data,dataelastic2);
  if (dataelastic2.size()>0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(dataelastic2);  // Unpack is done here
    MAT::So3Material* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel==NULL)
      dserror("failed to unpack elastic material");
    infty_conc_mat_ = Teuchos::rcp(matel);
  }
  else
    infty_conc_mat_ = Teuchos::null;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::ScalarDepInterp::ResetAll(const int numgp)
{
  zero_conc_mat_->ResetAll(numgp);
  infty_conc_mat_->ResetAll(numgp);
}

/*----------------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Update()
{
  zero_conc_mat_->Update();
  infty_conc_mat_->Update();
}

/*----------------------------------------------------------------------------*/
void MAT::ScalarDepInterp::ResetStep()
{
  zero_conc_mat_->ResetStep();
  infty_conc_mat_->ResetStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::VisNames(std::map<std::string,int>& names)
{
  std::string fiber = "zero_conc_ratio";
  names[fiber] = 1; // 1-dim vector

  zero_conc_mat_->VisNames(names);
  infty_conc_mat_->VisNames(names);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::ScalarDepInterp::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "zero_conc_ratio")
  {
    if ((int)data.size()!=1)
      dserror("size mismatch");

    double temp = 0.0;
    for (int gp=0; gp<numgp; gp++)
    {
      temp += zero_conc_ratio_.at(gp);
    }

    data[0] = temp/((double)numgp);
  }

  bool tmp1 = zero_conc_mat_->VisData(name, data, numgp, eleID);
  bool tmp2 = infty_conc_mat_->VisData(name, data, numgp, eleID);
  return (tmp1 and tmp2);
}
