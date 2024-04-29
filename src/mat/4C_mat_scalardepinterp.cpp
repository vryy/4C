/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains a framework, which allows for the linear interpolation between two different
strain energy functions, where the interpolation ratio come from 'outside'. For now this is supposed
to come from the idea, that there is a scalar quantity (e.g. foam cells) which induces growth. The
grown volume fraction thereby consists of another material (e.g. softer or stiffer material
parameters or totally different strain energy function).

\level 2

The input line should read
MAT 1 MAT_ScalarDepInterp IDMATZEROSC 2 IDMATUNITSC 3



*/

/*----------------------------------------------------------------------*/

#include "4C_mat_scalardepinterp.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScalarDepInterp::ScalarDepInterp(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      id_lambda_zero_(matdata->Get<int>("IDMATZEROSC")),
      id_lambda_unit_(matdata->Get<int>("IDMATUNITSC")){

      };

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ScalarDepInterp::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScalarDepInterp(this));
}


MAT::ScalarDepInterpType MAT::ScalarDepInterpType::instance_;


CORE::COMM::ParObject* MAT::ScalarDepInterpType::Create(const std::vector<char>& data)
{
  MAT::ScalarDepInterp* ScalarDepInterp = new MAT::ScalarDepInterp();
  ScalarDepInterp->Unpack(data);
  return ScalarDepInterp;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScalarDepInterp::ScalarDepInterp()
    : params_(nullptr),
      isinit_(false),
      lambda_zero_mat_(Teuchos::null),
      lambda_unit_mat_(Teuchos::null),
      lambda_(Teuchos::null)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScalarDepInterp::ScalarDepInterp(MAT::PAR::ScalarDepInterp* params)
    : params_(params),
      isinit_(false),
      lambda_zero_mat_(Teuchos::null),
      lambda_unit_mat_(Teuchos::null),
      lambda_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Setup(int numgp, INPUT::LineDefinition* linedef)
{
  if (isinit_)
    FOUR_C_THROW("This function should just be called, if the material is jet not initialized.");

  // Setup of elastic material for zero concentration
  lambda_zero_mat_ =
      Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Factory(params_->id_lambda_zero_));
  lambda_zero_mat_->Setup(numgp, linedef);

  // Setup of elastic material for zero concentration
  lambda_unit_mat_ =
      Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Factory(params_->id_lambda_unit_));
  lambda_unit_mat_->Setup(numgp, linedef);

  // Some safety check
  const double density1 = lambda_zero_mat_->Density();
  const double density2 = lambda_unit_mat_->Density();
  if (abs(density1 - density2) > 1e-14)
    FOUR_C_THROW(
        "The densities of the materials specified in IDMATZEROSC and IDMATUNITSC must be equal!");


  double lambda = 1.0;
  // Read lambda from input file, if available
  if (linedef->HaveNamed("lambda")) ReadLambda(linedef, "lambda", lambda);

  lambda_ = std::vector<double>(numgp, lambda);

  // initialization done
  isinit_ = true;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
    const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  // evaluate elastic material corresponding to zero concentration
  CORE::LINALG::Matrix<6, 1> stress_lambda_zero = *stress;
  CORE::LINALG::Matrix<6, 6> cmat_zero_conc = *cmat;

  lambda_zero_mat_->Evaluate(
      defgrd, glstrain, params, &stress_lambda_zero, &cmat_zero_conc, gp, eleGID);

  // evaluate elastic material corresponding to infinite concentration
  CORE::LINALG::Matrix<6, 1> stress_lambda_unit = *stress;
  CORE::LINALG::Matrix<6, 6> cmat_infty_conc = *cmat;

  lambda_unit_mat_->Evaluate(
      defgrd, glstrain, params, &stress_lambda_unit, &cmat_infty_conc, gp, eleGID);

  double lambda;
  // get the ratio of interpolation
  // NOTE: if no ratio is available, we use the IDMATZEROSC material!
  if (params.isParameter("lambda"))
  {
    lambda = params.get<double>("lambda");

    // NOTE: this would be nice, but since negative concentrations can occur,
    // we have to catch 'unnatural' cases different...
    //  if ( lambda < -1.0e-14 or lambda > (1.0+1.0e-14) )
    //      FOUR_C_THROW("The lambda must be in [0,1]!");

    // e.g. like that:
    if (lambda < -1.0e-14) lambda = 0.0;
    if (lambda > (1.0 + 1.0e-14)) lambda = 1.0;

    lambda_.at(gp) = lambda;
  }
  else
  {
    lambda = lambda_.at(gp);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // \mym Psi = (1-\lambda(C)) * \Psi_0(\mym C) + \lambda(C) * \Psi_1(\mym C) +  = ...
  // \mym S = 2 \frac{\partial}{\partial \mym C} \left( \Psi(\mym C) \right) = ...
  // \mym S = 2 \frac{\partial}{\partial \mym C} \left( (1-\lambda(C)) * \Psi_0(\mym C) + \lambda(C)
  // * \Psi_1(\mym C) \right) = ...
  ///////////////////////////////////////////////////////////////////////////////////////////////

  // do the linear interpolation between stresses:
  // ... = (1-\lambda(C)) * 2* \frac{\partial}{\partial \mym C} \Psi_0) + \lambda(C) * 2 *
  // \frac{\partial}{\partial \mym C} \Psi_1
  stress->Update(1.0 - lambda, stress_lambda_zero, lambda, stress_lambda_unit, 0.0);

  cmat->Update(1.0 - lambda, cmat_zero_conc, lambda, cmat_infty_conc, 0.0);

  if (params.isParameter("dlambda_dC"))
  {
    // get derivative of interpolation ratio w.r.t. glstrain
    Teuchos::RCP<CORE::LINALG::Matrix<6, 1>> dlambda_dC =
        params.get<Teuchos::RCP<CORE::LINALG::Matrix<6, 1>>>("dlambda_dC");

    // evaluate strain energy functions
    double psi_lambda_zero = 0.0;
    lambda_zero_mat_->StrainEnergy(*glstrain, psi_lambda_zero, gp, eleGID);

    double psi_lambda_unit = 0.0;
    lambda_unit_mat_->StrainEnergy(*glstrain, psi_lambda_unit, gp, eleGID);

    // and add the stresses due to possible dependency of the ratio w.r.t. to C
    // ... - 2 * \Psi_0 * \frac{\partial}{\partial \mym C} \lambda(C) ) + * 2 * \Psi_1 *
    // \frac{\partial}{\partial \mym C} \lambda(C)
    stress->Update(-2.0 * psi_lambda_zero, *dlambda_dC, +2.0 * psi_lambda_unit, *dlambda_dC, 1.0);

    // Note: for the linearization we do neglect the derivatives of the ratio w.r.t. glstrain
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Pack(CORE::COMM::PackBuffer& data) const
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

  int numgp = 0;
  if (isinit_)
  {
    numgp = lambda_.size();
    ;  // size is number of gausspoints
  }
  AddtoPack(data, numgp);

  for (int gp = 0; gp < numgp; gp++)
  {
    AddtoPack(data, lambda_.at(gp));
  }

  // Pack data of both elastic materials
  if (lambda_zero_mat_ != Teuchos::null and lambda_unit_mat_ != Teuchos::null)
  {
    lambda_zero_mat_->Pack(data);
    lambda_unit_mat_->Pack(data);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::ScalarDepInterp*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  int numgp;
  ExtractfromPack(position, data, numgp);
  if (not(numgp == 0))
  {
    lambda_ = std::vector<double>(numgp, 1.0);

    for (int gp = 0; gp < numgp; gp++)
    {
      double lambda = 1.0;
      ExtractfromPack(position, data, lambda);
      lambda_.at(gp) = lambda;
    }
  }

  // Unpack data of elastic material (these lines are copied from element.cpp)
  std::vector<char> dataelastic;
  ExtractfromPack(position, data, dataelastic);
  if (dataelastic.size() > 0)
  {
    CORE::COMM::ParObject* o = CORE::COMM::Factory(dataelastic);  // Unpack is done here
    MAT::So3Material* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel == nullptr) FOUR_C_THROW("failed to unpack elastic material");
    lambda_zero_mat_ = Teuchos::rcp(matel);
  }
  else
    lambda_zero_mat_ = Teuchos::null;

  // Unpack data of elastic material (these lines are copied from element.cpp)
  std::vector<char> dataelastic2;
  ExtractfromPack(position, data, dataelastic2);
  if (dataelastic2.size() > 0)
  {
    CORE::COMM::ParObject* o = CORE::COMM::Factory(dataelastic2);  // Unpack is done here
    MAT::So3Material* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel == nullptr) FOUR_C_THROW("failed to unpack elastic material");
    lambda_unit_mat_ = Teuchos::rcp(matel);
  }
  else
    lambda_unit_mat_ = Teuchos::null;

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::ScalarDepInterp::Update()
{
  lambda_zero_mat_->Update();
  lambda_unit_mat_->Update();
}

/*----------------------------------------------------------------------------*/
void MAT::ScalarDepInterp::ResetStep()
{
  lambda_zero_mat_->ResetStep();
  lambda_unit_mat_->ResetStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::VisNames(std::map<std::string, int>& names)
{
  std::string fiber = "lambda";
  names[fiber] = 1;  // 1-dim vector

  lambda_zero_mat_->VisNames(names);
  lambda_unit_mat_->VisNames(names);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool MAT::ScalarDepInterp::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "lambda")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");

    double temp = 0.0;
    for (int gp = 0; gp < numgp; gp++)
    {
      temp += lambda_.at(gp);
    }

    data[0] = temp / ((double)numgp);
  }

  bool tmp1 = lambda_zero_mat_->VisData(name, data, numgp, eleID);
  bool tmp2 = lambda_unit_mat_->VisData(name, data, numgp, eleID);
  return (tmp1 and tmp2);
}

void MAT::ScalarDepInterp::ReadLambda(
    INPUT::LineDefinition* linedef, std::string specifier, double& lambda)
{
  linedef->ExtractDouble(specifier, lambda);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScalarDepInterp::StrainEnergy(
    const CORE::LINALG::Matrix<6, 1>& glstrain, double& psi, const int gp, const int eleGID)
{
  // evaluate strain energy functions
  double psi_lambda_zero = 0.0;
  lambda_zero_mat_->StrainEnergy(glstrain, psi_lambda_zero, gp, eleGID);

  double psi_lambda_unit = 0.0;
  lambda_unit_mat_->StrainEnergy(glstrain, psi_lambda_unit, gp, eleGID);

  double lambda = 0.0;
  for (unsigned gp = 0; gp < lambda_.size(); gp++)
  {
    lambda += lambda_.at(gp);
  }
  lambda = lambda / (lambda_.size());

  psi += (1 - lambda) * psi_lambda_zero + lambda * psi_lambda_unit;
  return;
}

FOUR_C_NAMESPACE_CLOSE
