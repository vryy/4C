/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation on an active stress material

\level 2
*/
/*---------------------------------------------------------------------*/

#include "4C_matelast_anisoactivestress_evolution.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_anisotropy_extension.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::AnisoActiveStressEvolution::AnisoActiveStressEvolution(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : ParameterAniso(matdata),
      sigma_(matdata.parameters.Get<double>("SIGMA")),
      tauc0_(matdata.parameters.Get<double>("TAUC0")),
      maxactiv_(matdata.parameters.Get<double>("MAX_ACTIVATION")),
      minactiv_(matdata.parameters.Get<double>("MIN_ACTIVATION")),
      activationthreshold_(matdata.parameters.Get<double>("ACTIVATION_THRES")),
      sourceactiv_(matdata.parameters.Get<int>("SOURCE_ACTIVATION")),
      strain_dep_(matdata.parameters.Get<bool>("STRAIN_DEPENDENCY")),
      lambda_lower_(matdata.parameters.Get<double>("LAMBDA_LOWER")),
      lambda_upper_(matdata.parameters.Get<double>("LAMBDA_UPPER")),
      gamma_(matdata.parameters.Get<double>("GAMMA")),
      init_(matdata.parameters.Get<int>("INIT")),
      adapt_angle_(matdata.parameters.Get<bool>("ADAPT_ANGLE"))
{
}

Mat::Elastic::AnisoActiveStressEvolution::AnisoActiveStressEvolution(
    Mat::Elastic::PAR::AnisoActiveStressEvolution* params)
    : params_(params),
      tauc_np_(0.0),
      tauc_n_(0.0),
      anisotropy_extension_(params->init_, params_->gamma_, params_->adapt_angle_ != 0,
          params_->structural_tensor_strategy(), {0})
{
  anisotropy_extension_.register_needed_tensors(
      FiberAnisotropyExtension<1>::FIBER_VECTORS |
      FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS);
}

void Mat::Elastic::AnisoActiveStressEvolution::PackSummand(
    Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, tauc_n_);
  anisotropy_extension_.pack_anisotropy(data);
}

void Mat::Elastic::AnisoActiveStressEvolution::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  extract_from_pack(position, data, tauc_n_);
  anisotropy_extension_.unpack_anisotropy(data, position);
}

void Mat::Elastic::AnisoActiveStressEvolution::register_anisotropy_extensions(
    Mat::Anisotropy& anisotropy)
{
  anisotropy.register_anisotropy_extension(anisotropy_extension_);
}

void Mat::Elastic::AnisoActiveStressEvolution::Setup(int numgp, Input::LineDefinition* linedef)
{
  // Setup of active stress model
  tauc_n_ = params_->tauc0_;
  tauc_np_ = params_->tauc0_;

  // reasonability check...
  if (params_->strain_dep_)
  {
    if (params_->lambda_lower_ >= params_->lambda_upper_)
    {
      FOUR_C_THROW(
          "LAMBDA_LOWER should be lower than LAMBDA_UPPER. Seems reasonable, doesn't it? Dude...");
    }
  }
}

void Mat::Elastic::AnisoActiveStressEvolution::post_setup(Teuchos::ParameterList& params)
{
  Summand::post_setup(params);
}

void Mat::Elastic::AnisoActiveStressEvolution::add_stress_aniso_principal(
    const Core::LinAlg::Matrix<6, 1>& rcg, Core::LinAlg::Matrix<6, 6>& cmat,
    Core::LinAlg::Matrix<6, 1>& stress, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  // Virtual GP (is zero for element fibers, otherwise it is the current GP)
  Core::LinAlg::Matrix<6, 1> A = anisotropy_extension_.get_structural_tensor_stress(gp, 0);

  double dt = params.get("delta time", -1.0);
  if (dt < 0.0)
  {
    FOUR_C_THROW("Parameter 'delta time' could not be read!");
  }

  double activationFunction = 0.0;
  Teuchos::RCP<Core::Mat::Material> scatramat;
  if (params_->sourceactiv_ == 0)
  {
    activationFunction = params.get<double>("scalar", 0.0);
  }
  else if (params_->sourceactiv_ == -1)
  {
    if (params.isParameter("gp_conc"))
    {
      // get pointer to vector containing the scalar values at the Gauss points
      Teuchos::RCP<std::vector<std::vector<double>>> conc =
          params.get<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc");
      // safety check
      if (conc == Teuchos::null)
      {
        FOUR_C_THROW("No concentration from scatra provided for the activation");
      }
      // safety check
      if (gp == -1)
      {
        FOUR_C_THROW("No Gauss point number provided in material");
      }
      // check if activation is above the given threshold
      if (conc->at(gp)[0] >= params_->activationthreshold_)
      {
        activationFunction = 1.0;
      }
    }
  }
  else
  {
    double totaltime = params.get<double>("total time", -1.0);
    if (totaltime < 0.0)
    {
      FOUR_C_THROW("Parameter 'total time' could not be read!");
    }
    const auto& element_center_coordinates_ref =
        params.get<Core::LinAlg::Matrix<3, 1>>("elecenter_coords_ref");
    activationFunction =
        Global::Problem::Instance()
            ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(params_->sourceactiv_ - 1)
            .Evaluate(element_center_coordinates_ref.A(), totaltime, 0);
  }

  double lambda = 0.0;
  double scale = 0.0;
  double n0 = 1.0;
  if (params_->strain_dep_)
  {
    // squared stretch along fiber direction
    double I4;
    I4 = A(0) * rcg(0) + A(1) * rcg(1) + A(2) * rcg(2) + A(3) * rcg(3) + A(4) * rcg(4) +
         A(5) * rcg(5);
    // stretch along fiber direction
    lambda = sqrt(I4);
    scale = -4. / (pow((params_->lambda_lower_ - params_->lambda_upper_), 2.0));
    // Frank-Starling factor n0 for strain-dependent contractility
    if (lambda >= params_->lambda_lower_ and lambda <= params_->lambda_upper_)
    {
      n0 = scale * (lambda - params_->lambda_lower_) * (lambda - params_->lambda_upper_);
    }
    else
    {
      n0 = 0.0;
    }
  }

  activationFunction =
      activationFunction * (params_->maxactiv_ - params_->minactiv_) + params_->minactiv_;
  double abs_u_ = abs(activationFunction);
  double absplus_u_ = abs_u_ * static_cast<double>(activationFunction > 0.0);
  tauc_np_ = (tauc_n_ / dt + n0 * params_->sigma_ * absplus_u_) / (1 / dt + abs_u_);
  stress.Update(tauc_np_, A, 1.0);

  // only contribution to cmat if we have strain dependency!
  if (params_->strain_dep_)
  {
    double dtauc_np_dC;
    // Cmat_active = 2 * dS_active / dC = 2 * dtau(t) / dC \otimes (a_0 \otimes a_0)
    if (lambda >= params_->lambda_lower_ and lambda <= params_->lambda_upper_)
    {
      dtauc_np_dC = 2. * scale *
                    (1. - (params_->lambda_upper_ + params_->lambda_lower_) / (2. * lambda)) *
                    params_->sigma_ * absplus_u_ / (1 / dt + abs_u_);
    }
    else
    {
      dtauc_np_dC = 0.0;
    }
    cmat.MultiplyNT(dtauc_np_dC, A, A, 1.0);
  }
}

void Mat::Elastic::AnisoActiveStressEvolution::GetFiberVecs(
    std::vector<Core::LinAlg::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  if (params_->init_ == DefaultAnisotropyExtension<1>::INIT_MODE_NODAL_FIBERS)
  {
    // This method expects constant fibers within this element but the init mode is such that
    // fibers are defined on the Gauss points
    // We therefore cannot return sth here.

    // ToDo: This may needs improvements later on if needed!
    return;
  }

  fibervecs.push_back(anisotropy_extension_.get_fiber(BaseAnisotropyExtension::GPDEFAULT, 0));
}

// Update internal stress variables
void Mat::Elastic::AnisoActiveStressEvolution::Update() { tauc_n_ = tauc_np_; }

void Mat::Elastic::AnisoActiveStressEvolution::SetFiberVecs(const double newgamma,
    const Core::LinAlg::Matrix<3, 3>& locsys, const Core::LinAlg::Matrix<3, 3>& defgrd)
{
  anisotropy_extension_.set_fiber_vecs(newgamma, locsys, defgrd);
}

FOUR_C_NAMESPACE_CLOSE
