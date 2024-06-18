/*----------------------------------------------------------------------*/
/*! \file
 \brief wrapper for structure material of porous media including reactive reference porosity


\level 3
 *----------------------------------------------------------------------*/

#include "4C_mat_structporo_reaction.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::StructPoroReaction::StructPoroReaction(const Core::Mat::PAR::Parameter::Data& matdata)
    : StructPoro(matdata), dofIDReacScalar_(matdata.parameters.get<int>("DOFIDREACSCALAR"))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::StructPoroReaction::create_material()
{
  return Teuchos::rcp(new Mat::StructPoroReaction(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::StructPoroReactionType Mat::StructPoroReactionType::instance_;

Core::Communication::ParObject* Mat::StructPoroReactionType::Create(const std::vector<char>& data)
{
  Mat::StructPoroReaction* struct_poro = new Mat::StructPoroReaction();
  struct_poro->Unpack(data);
  return struct_poro;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::StructPoroReaction::StructPoroReaction()
    : params_(nullptr), refporosity_(-1.0), dphiDphiref_(0.0), refporositydot_(0.0)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::StructPoroReaction::StructPoroReaction(Mat::PAR::StructPoroReaction* params)
    : StructPoro(params),
      params_(params),
      refporosity_(-1.0),
      dphiDphiref_(0.0),
      refporositydot_(0.0)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroReaction::setup(int numgp, Input::LineDefinition* linedef)
{
  StructPoro::setup(numgp, linedef);
  refporosity_ = params_->init_porosity_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::StructPoroReaction::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // refporosity_
  add_to_pack(data, refporosity_);

  // add base class material
  StructPoro::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::StructPoroReaction::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::StructPoroReaction*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // refporosity_
  extract_from_pack(position, data, refporosity_);

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  StructPoro::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::StructPoroReaction::compute_porosity(Teuchos::ParameterList& params, double press,
    double J, int gp, double& porosity, double* dphi_dp, double* dphi_dJ, double* dphi_dJdp,
    double* dphi_dJJ, double* dphi_dpp, bool save)
{
  // evaluate change of reference porosity due to reaction

  // TODO: do not read from parameter list!
  if (params.isParameter("scalar"))
  {
    Teuchos::RCP<std::vector<double>> scalars =
        params.get<Teuchos::RCP<std::vector<double>>>("scalar");
    reaction(porosity, J, scalars, params);
  }

  // call base class to compute porosity
  StructPoro::compute_porosity(refporosity_, press, J, gp, porosity, dphi_dp, dphi_dJ, dphi_dJdp,
      dphi_dJJ, dphi_dpp, &dphiDphiref_, save);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroReaction::constitutive_derivatives(Teuchos::ParameterList& params, double press,
    double J, double porosity, double* dW_dp, double* dW_dphi, double* dW_dJ, double* dW_dphiref,
    double* W)
{
  if (porosity == 0.0)
    FOUR_C_THROW(
        "porosity equals zero!! Wrong initial porosity? (or wrong collagen density for ecm "
        "material)");

  // evaluate change of reference porosity due to reaction

  // TODO: do not read from parameter list!
  if (params.isParameter("scalar"))
  {
    Teuchos::RCP<std::vector<double>> scalars =
        params.get<Teuchos::RCP<std::vector<double>>>("scalar");
    reaction(porosity, J, scalars, params);
  }

  // call base class
  StructPoro::constitutive_derivatives(
      params, press, J, porosity, refporosity_, dW_dp, dW_dphi, dW_dJ, dW_dphiref, W);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroReaction::reaction(const double porosity, const double J,
    Teuchos::RCP<std::vector<double>> scalars, Teuchos::ParameterList& params)
{
  // double dt = params.get<double>("delta time",-1.0);
  double time = params.get<double>("total time", -1.0);

  // use scalar for this type of reaction
  double cnp = scalars->at(params_->dofIDReacScalar_);

  if (time != -1.0)
  {
    // double k = 1.0;
    double tau = 200.0 * cnp;     ///(cnp+k); 20.0/(20*cnp+1.0)
    double limitporosity = 0.45;  // 0.8;

    refporosity_ =
        limitporosity - (limitporosity - params_->init_porosity_) * exp(-1.0 * time / tau);
    refporositydot_ = (limitporosity - params_->init_porosity_) / tau * exp(-1.0 * time / tau);
  }
  else  //(time==-1.0) -> time not set (this happens during setup-> no reaction)
  {
    // do nothing
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroReaction::evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  // call base class
  StructPoro::evaluate(defgrd, glstrain, params, stress, cmat, gp, eleGID);

  // scale stresses and cmat
  stress->Scale((1.0 - refporosity_) / (1.0 - params_->init_porosity_));
  cmat->Scale((1.0 - refporosity_) / (1.0 - params_->init_porosity_));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::StructPoroReaction::RefPorosityAv() const { return refporosity_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroReaction::VisNames(std::map<std::string, int>& names)
{
  // call base class
  StructPoro::VisNames(names);
  std::string name = "reference_porosity";
  names[name] = 1;  // scalar
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Mat::StructPoroReaction::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  // call base class
  if (StructPoro::VisData(name, data, numgp, eleID)) return true;
  if (name == "reference_porosity")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
    data[0] = RefPorosityAv();
    return true;
  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
