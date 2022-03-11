/*----------------------------------------------------------------------*/
/*! \file
\brief  fluid material for poroelasticity problems


\level 2
 *-----------------------------------------------------------------------*/

#include <vector>
#include "fluidporo.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/matpar_bundle.H"

namespace MAT
{
  namespace FLUIDPORO
  {
    /// compute structure tensor from given anisotropy direction - 2D
    void CreateStructureTensorFromAnisotropyDirection(
        const std::vector<double>& anisotropy_direction, LINALG::Matrix<2, 2>& structure_tensor)
    {
      structure_tensor.Clear();

      // Factor to normalize the structure tensor
      const double square_length = (anisotropy_direction[0] * anisotropy_direction[0] +
                                    anisotropy_direction[1] * anisotropy_direction[1]);

      if (square_length > EPS8)
      {
        for (int i = 0; i < 2; i++)
          for (int j = 0; j < 2; j++)
          {
            structure_tensor(i, j) = anisotropy_direction[i] * anisotropy_direction[j];
          }

        structure_tensor.Scale(1. / square_length);
      }
    }

    /// compute structure tensor from given anisotropy direction - 3D
    void CreateStructureTensorFromAnisotropyDirection(
        const std::vector<double>& anisotropy_direction, LINALG::Matrix<3, 3>& structure_tensor)
    {
      structure_tensor.Clear();

      // Factor to normalize the structure tensor
      const double square_length = (anisotropy_direction[0] * anisotropy_direction[0] +
                                    anisotropy_direction[1] * anisotropy_direction[1] +
                                    anisotropy_direction[2] * anisotropy_direction[2]);

      if (square_length > EPS8)
      {
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
          {
            structure_tensor(i, j) = anisotropy_direction[i] * anisotropy_direction[j];
          }

        structure_tensor.Scale(1. / square_length);
      }
    }

    /*!
     * @brief  Base class for the anisotropy strategy in the FluidPoro material.
     *
     * For materials with anisotropic permeability properties, the reaction tensor will represent
     * the type of anisotropy and therefore will be calculated in different ways. This base class
     * defines the methods that are affected by the anisotropy. These methods are responsible for
     * computations regarding the reaction tensor. Some methods are already implemented by the base
     * class. They can be overridden depending on the needs of the specific type of anisotropy.
     */
    class PoroAnisotropyStrategyBase
    {
     public:
      /// Simple constructor
      PoroAnisotropyStrategyBase(const MAT::PAR::FluidPoro* params) : params_(params){};

      /// Virtual default destructor
      virtual ~PoroAnisotropyStrategyBase() = default;

      /// compute reaction coefficient
      virtual double ComputeReactionCoeff() const = 0;

      /// compute reaction tensor - 2D
      virtual void ComputeReactionTensor(LINALG::Matrix<2, 2>& reaction_tensor, const double J,
          const double porosity,
          const std::vector<std::vector<double>>& anisotropic_permeability_directions) const = 0;

      /// compute reaction tensor - 3D
      virtual void ComputeReactionTensor(LINALG::Matrix<3, 3>& reaction_tensor, const double J,
          const double porosity,
          const std::vector<std::vector<double>>& anisotropic_permeability_directions) const = 0;

      /// compute linearization of the reaction tensor - 2D
      virtual void ComputeLinMatReactionTensor(LINALG::Matrix<2, 2>& linreac_dphi,
          LINALG::Matrix<2, 2>& linreac_dJ, const double J, const double porosity) const
      {
        linreac_dphi.Clear();
        linreac_dJ.Clear();
      };

      /// compute linearization of the reaction tensor - 3D
      virtual void ComputeLinMatReactionTensor(LINALG::Matrix<3, 3>& linreac_dphi,
          LINALG::Matrix<3, 3>& linreac_dJ, const double J, const double porosity) const
      {
        linreac_dphi.Clear();
        linreac_dJ.Clear();
      };

     protected:
      /// Material parameters
      const MAT::PAR::FluidPoro* params_;
    };

    /*!
     * @brief Reaction tensor (anisotropy) strategy for isotropy.
     *
     * In the case of isotropy, both 2D and 3D simulations are allowed. Furthermore, the isotropic
     * case implements two different strategies for the permeability-porosity dependency which are
     *  -# Constant material permeability, and
     *  -# Kozeny-Carman-Equation.
     */
    class PoroIsotropyStrategy : public PoroAnisotropyStrategyBase
    {
     public:
      /// Simple constructor
      PoroIsotropyStrategy(const MAT::PAR::FluidPoro* params)
          : PoroAnisotropyStrategyBase(params){};

      /// compute isotropy reaction coefficient
      double ComputeReactionCoeff() const override
      {
        const double viscosity = params_->viscosity_;
        const double permeability = params_->permeability_;

        if (viscosity < EPS15) dserror("zero or negative viscosity");
        if (permeability < EPS15) dserror("zero or negative permeability");

        // trace of the reaction tensor divided by the dimension
        const double reacoeff = viscosity / permeability;

        return reacoeff;
      };

      /// compute isotropic reaction tensor - 2D
      void ComputeReactionTensor(LINALG::Matrix<2, 2>& reaction_tensor, const double J,
          const double porosity,
          const std::vector<std::vector<double>>& anisotropic_permeability_directions)
          const override
      {
        double reacoeff = ComputeReactionCoeff();
        const double permeability_correction_factor = params_->permeabilitycorrectionfactor_;
        const auto permeability_function = params_->permeabilityfunc_;

        reaction_tensor.Clear();

        if (permeability_function == MAT::PAR::kozeny_karman)
        {
          reacoeff *= (1 - porosity * porosity * J * J) /
                      (porosity * porosity * porosity * J * J * J) * permeability_correction_factor;
        }

        for (int i = 0; i < 2; i++) reaction_tensor(i, i) = reacoeff;
      };

      /// compute isotropic reaction tensor - 3D
      void ComputeReactionTensor(LINALG::Matrix<3, 3>& reaction_tensor, const double J,
          const double porosity,
          const std::vector<std::vector<double>>& anisotropic_permeability_directions)
          const override
      {
        double reacoeff = ComputeReactionCoeff();
        const double permeability_correction_factor = params_->permeabilitycorrectionfactor_;
        const auto permeability_function = params_->permeabilityfunc_;

        reaction_tensor.Clear();

        if (permeability_function == MAT::PAR::kozeny_karman)
        {
          reacoeff *= (1 - porosity * porosity * J * J) /
                      (porosity * porosity * porosity * J * J * J) * permeability_correction_factor;
        }

        for (int i = 0; i < 3; i++) reaction_tensor(i, i) = reacoeff;
      };

      /// compute linearization of isotropic reaction tensor - 2D
      void ComputeLinMatReactionTensor(LINALG::Matrix<2, 2>& linreac_dphi,
          LINALG::Matrix<2, 2>& linreac_dJ, const double J, const double porosity) const override
      {
        const double reacoeff = ComputeReactionCoeff();
        const double permeability_correction_factor = params_->permeabilitycorrectionfactor_;
        const auto permeability_function = params_->permeabilityfunc_;

        linreac_dphi.Clear();
        linreac_dJ.Clear();

        if (permeability_function == MAT::PAR::const_)
          return;  // Permeability is not a function of porosity or J
        else if (permeability_function == MAT::PAR::kozeny_karman)
        {
          // d(isotropic_mat_reactiontensor)/d(phi) = reacoeff * [(J * phi)^2 - 3] / ( J^3 * phi^4 )
          // d(isotropic_mat_reactiontensor)/d(J) = reacoeff * [(J * phi)^2 - 3] / ( J^4 * phi^3 )

          double linreac_tmp = reacoeff * ((J * J * porosity * porosity) - 3.0) /
                               (J * J * J * porosity * porosity * porosity) *
                               permeability_correction_factor;
          linreac_dphi(0, 0) = linreac_tmp / porosity;
          linreac_dJ(0, 0) = linreac_tmp / J;

          for (int i = 1; i < 2; i++)
          {
            linreac_dphi(i, i) = linreac_dphi(0, 0);
            linreac_dJ(i, i) = linreac_dJ(0, 0);
          }
        }
      };

      /// compute linearization of isotropic reaction tensor - 3D
      void ComputeLinMatReactionTensor(LINALG::Matrix<3, 3>& linreac_dphi,
          LINALG::Matrix<3, 3>& linreac_dJ, const double J, const double porosity) const override
      {
        const double reacoeff = ComputeReactionCoeff();
        const double permeability_correction_factor = params_->permeabilitycorrectionfactor_;
        const auto permeability_function = params_->permeabilityfunc_;

        linreac_dphi.Clear();
        linreac_dJ.Clear();

        if (permeability_function == MAT::PAR::const_)
          return;  // Permeability is not a function of porosity or J
        else if (permeability_function == MAT::PAR::kozeny_karman)
        {
          // d(isotropic_mat_reactiontensor)/d(phi) = reacoeff * [(J * phi)^2 - 3] / ( J^3 * phi^4 )
          // d(isotropic_mat_reactiontensor)/d(J) = reacoeff * [(J * phi)^2 - 3] / ( J^4 * phi^3 )

          double linreac_tmp = reacoeff * ((J * J * porosity * porosity) - 3.0) /
                               (J * J * J * porosity * porosity * porosity) *
                               permeability_correction_factor;
          linreac_dphi(0, 0) = linreac_tmp / porosity;
          linreac_dJ(0, 0) = linreac_tmp / J;

          for (int i = 1; i < 3; i++)
          {
            linreac_dphi(i, i) = linreac_dphi(0, 0);
            linreac_dJ(i, i) = linreac_dJ(0, 0);
          }
        }
      };
    };

    /*!
     * @brief Reaction tensor (anisotropy) strategy for constant material transverse isotropy.
     *
     * In the case of transverse isotropy, both 2D and 3D simulations can be done. For now,
     * only constant material permeabilities are allowed, i.e. there is no dependency between the
     * permeability and porosity in the material reaction tensor.
     */
    class PoroConstantMaterialTransverseIsotropyStrategy : public PoroAnisotropyStrategyBase
    {
     public:
      PoroConstantMaterialTransverseIsotropyStrategy(const MAT::PAR::FluidPoro* params)
          : PoroAnisotropyStrategyBase(params){};

      /// compute reaction coefficient for constant material transverse isotropy
      double ComputeReactionCoeff() const override
      {
        const double viscosity = params_->viscosity_;
        const double permeability = params_->permeability_;
        const double axial_permeability = params_->axial_permeability_;

        if (viscosity < EPS15) dserror("zero or negative viscosity");
        if (permeability < EPS15) dserror("zero or negative permeability");
        if (axial_permeability < EPS15) dserror("zero or negative axial permeability");

        // trace of the 3D reaction tensor divided by 3 (even if the problem is in 2D)
        const double reacoeff = viscosity / 3. * (1. / axial_permeability + 2. / permeability);

        return reacoeff;
      };

      /// compute reaction tensor for constant material transverse isotropy - 2D
      void ComputeReactionTensor(LINALG::Matrix<2, 2>& reaction_tensor, const double J,
          const double porosity,
          const std::vector<std::vector<double>>& anisotropic_permeability_directions)
          const override
      {
        const double dynamic_viscosity = params_->viscosity_;
        const double permeability = params_->permeability_;
        const double axial_permeability = params_->axial_permeability_;

        if (dynamic_viscosity < EPS15) dserror("zero or negative viscosity");
        if (permeability < EPS15) dserror("zero or negative permeability");
        if (axial_permeability < EPS15) dserror("zero or negative axial permeability");

        reaction_tensor.Clear();

        LINALG::Matrix<2, 2> permeability_tensor(true);
        LINALG::Matrix<2, 2> structure_tensor(true);
        CreateStructureTensorFromAnisotropyDirection(
            anisotropic_permeability_directions[0], structure_tensor);

        // Transverse component of the transversely isotropic permeability tensor
        for (int i = 0; i < 2; ++i) permeability_tensor(i, i) += permeability;

        // Axial component of the transversely isotropic permeability tensor
        permeability_tensor.Update(axial_permeability - permeability, structure_tensor, 1.0);

        reaction_tensor.Invert(permeability_tensor);
        reaction_tensor.Scale(dynamic_viscosity);
      };

      /// compute reaction tensor for constant material transverse isotropy - 3D
      void ComputeReactionTensor(LINALG::Matrix<3, 3>& reaction_tensor, const double J,
          const double porosity,
          const std::vector<std::vector<double>>& anisotropic_permeability_directions)
          const override
      {
        const double dynamic_viscosity = params_->viscosity_;
        const double permeability = params_->permeability_;
        const double axial_permeability = params_->axial_permeability_;

        if (dynamic_viscosity < EPS15) dserror("zero or negative viscosity");
        if (permeability < EPS15) dserror("zero or negative permeability");
        if (axial_permeability < EPS15) dserror("zero or negative axial permeability");

        reaction_tensor.Clear();

        LINALG::Matrix<3, 3> permeability_tensor(true);
        LINALG::Matrix<3, 3> structure_tensor(true);
        CreateStructureTensorFromAnisotropyDirection(
            anisotropic_permeability_directions[0], structure_tensor);

        // Transverse component of the transversely isotropic permeability tensor
        for (int i = 0; i < 3; ++i) permeability_tensor(i, i) += permeability;

        // Axial component of the transversely isotropic permeability tensor
        permeability_tensor.Update(axial_permeability - permeability, structure_tensor, 1.0);

        reaction_tensor.Invert(permeability_tensor);
        reaction_tensor.Scale(dynamic_viscosity);
      };
    };

    /*!
     * @brief Reaction tensor (anisotropy) strategy for constant material orthotropy.
     *
     * Orthotropy can only be used in 3D simulations. For now, only constant material permeabilities
     * are allowed, i.e. there is no dependency between the permeability and porosity in the
     * material reaction tensor.
     */
    class PoroConstantMaterialOrthotropyStrategy : public PoroAnisotropyStrategyBase
    {
     public:
      /// Simple constructor
      PoroConstantMaterialOrthotropyStrategy(const MAT::PAR::FluidPoro* params)
          : PoroAnisotropyStrategyBase(params){};

      /// compute reaction coefficient for constant material orthotropy
      double ComputeReactionCoeff() const override
      {
        const double viscosity = params_->viscosity_;
        const std::vector<double> orthotropic_permeabilities = params_->orthotropic_permeabilities_;

        if (viscosity < EPS15) dserror("zero or negative viscosity");
        for (int dim = 0; dim < 3; ++dim)
          if (orthotropic_permeabilities[dim] < EPS15)
            dserror("zero or negative permeability in direction " + std::to_string(dim + 1));

        // trace of the reaction tensor divided by 3
        const double reacoeff =
            viscosity / 3. *
            (1. / orthotropic_permeabilities[0] + 1. / orthotropic_permeabilities[1] +
                1. / orthotropic_permeabilities[2]);

        return reacoeff;
      };

      /// constant material orthotropy in 2D is not allowed
      void ComputeReactionTensor(LINALG::Matrix<2, 2>& reaction_tensor, const double J,
          const double porosity,
          const std::vector<std::vector<double>>& anisotropic_permeability_directions)
          const override
      {
        dserror("Use transverse isotropy in 2D!");
      };

      /// compute reaction tensor for constant material orthotropy - 3D
      void ComputeReactionTensor(LINALG::Matrix<3, 3>& reaction_tensor, const double J,
          const double porosity,
          const std::vector<std::vector<double>>& anisotropic_permeability_directions)
          const override
      {
        const double dynamic_viscosity = params_->viscosity_;
        const std::vector<double> orthotropic_permeabilities = params_->orthotropic_permeabilities_;

        if (dynamic_viscosity < EPS15) dserror("zero or negative viscosity");
        for (int dim = 0; dim < 3; ++dim)
          if (orthotropic_permeabilities[dim] < EPS15)
            dserror("zero or negative permeability in direction " + std::to_string(dim + 1));

        reaction_tensor.Clear();

        LINALG::Matrix<3, 3> permeability_tensor(true);
        LINALG::Matrix<3, 3> structure_tensor(true);

        for (int dim = 0; dim < 3; ++dim)
        {
          CreateStructureTensorFromAnisotropyDirection(
              anisotropic_permeability_directions[dim], structure_tensor);
          permeability_tensor.Update(orthotropic_permeabilities[dim], structure_tensor, 1.0);
        }

        reaction_tensor.Invert(permeability_tensor);
        reaction_tensor.Scale(dynamic_viscosity);
      };

      /// constant material orthotropy in 2D is not allowed
      void ComputeLinMatReactionTensor(LINALG::Matrix<2, 2>& linreac_dphi,
          LINALG::Matrix<2, 2>& linreac_dJ, const double J, const double porosity) const override
      {
        dserror("Use transverse isotropy in 2D!");
      };
    };

    Teuchos::RCP<MAT::FLUIDPORO::PoroAnisotropyStrategyBase> CreateAnisotropyStrategy(
        const MAT::PAR::FluidPoro* params)
    {
      switch (params->permeabilityfunc_)
      {
        case MAT::PAR::const_:
        case MAT::PAR::kozeny_karman:
          return Teuchos::rcp(new MAT::FLUIDPORO::PoroIsotropyStrategy(params));
        case MAT::PAR::const_material_transverse:
          return Teuchos::rcp(
              new MAT::FLUIDPORO::PoroConstantMaterialTransverseIsotropyStrategy(params));
        case MAT::PAR::const_material_orthotropic:
          return Teuchos::rcp(new MAT::FLUIDPORO::PoroConstantMaterialOrthotropyStrategy(params));
        default:
          return Teuchos::null;
      }
    }

  }  // namespace FLUIDPORO
}  // namespace MAT

/*----------------------------------------------------------------------*
 *  constructor (public)                               vuong 06/11      |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoro::FluidPoro(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->GetDouble("DYNVISCOSITY")),
      density_(matdata->GetDouble("DENSITY")),
      permeability_(matdata->GetDouble("PERMEABILITY")),
      axial_permeability_(matdata->GetDouble("AXIALPERMEABILITY")),
      type_(undefined),
      varyingpermeability_(false),
      permeabilityfunc_(MAT::PAR::pf_undefined),
      permeabilitycorrectionfactor_(1.0),
      initialporosity_(1.0)
{
  const std::string* typestring = matdata->Get<std::string>("TYPE");

  if (*typestring == "Darcy")
    type_ = darcy;
  else if (*typestring == "Darcy-Brinkman")
    type_ = darcy_brinkman;

  const std::string* pfuncstring = matdata->Get<std::string>("PERMEABILITYFUNCTION");

  if (*pfuncstring == "Const")
    permeabilityfunc_ = MAT::PAR::const_;
  else if (*pfuncstring == "Kozeny_Carman")
    permeabilityfunc_ = MAT::PAR::kozeny_karman;
  else if (*pfuncstring == "Const_Material_Transverse")
    permeabilityfunc_ = MAT::PAR::const_material_transverse;
  else if (*pfuncstring == "Const_Material_Orthotropy")
    permeabilityfunc_ = MAT::PAR::const_material_orthotropic;
  else
    dserror("Unknown permeability function: %s", pfuncstring->c_str());

  orthotropic_permeabilities_.resize(3, 0.0);
  if (permeabilityfunc_ == MAT::PAR::const_material_orthotropic)
  {
    for (int dim = 0; dim < 3; ++dim)
      orthotropic_permeabilities_[dim] =
          matdata->GetDouble("ORTHOPERMEABILITY" + std::to_string(dim + 1));
  }
}

/*----------------------------------------------------------------------*
 *  Create Material (public)                             vuong 06/11      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoro::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoro(this));
}

/*----------------------------------------------------------------------*
  Set Initial Porosity (public)                          vuong 06/11     |
*----------------------------------------------------------------------*/
void MAT::PAR::FluidPoro::SetInitialPorosity(double initialporosity)
{
  initialporosity_ = initialporosity;

  if (permeabilityfunc_ == MAT::PAR::const_)
  {
    permeabilitycorrectionfactor_ = 1.0;
  }
  else if (permeabilityfunc_ == MAT::PAR::kozeny_karman)
  {
    // c = (phi0^3 / (1 - phi0^2))
    permeabilitycorrectionfactor_ = initialporosity_ * initialporosity_ * initialporosity_ /
                                    (1 - initialporosity_ * initialporosity_);
  }
  return;
}

/*----------------------------------------------------------------------*
                                                          vuong 06/11     |
*----------------------------------------------------------------------*/
MAT::FluidPoroType MAT::FluidPoroType::instance_;

/*----------------------------------------------------------------------*
 *                                                           vuong 06/11 |
 *----------------------------------------------------------------------*/

DRT::ParObject* MAT::FluidPoroType::Create(const std::vector<char>& data)
{
  MAT::FluidPoro* fluid_poro = new MAT::FluidPoro();
  fluid_poro->Unpack(data);
  return fluid_poro;
}

MAT::FluidPoro::FluidPoro() : params_(nullptr), anisotropy_strategy_(Teuchos::null) {}

MAT::FluidPoro::FluidPoro(MAT::PAR::FluidPoro* params) : params_(params)
{
  anisotropy_strategy_ = MAT::FLUIDPORO::CreateAnisotropyStrategy(params);
}

/*----------------------------------------------------------------------*
 *                                                          vuong 06/11 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoro::Pack(DRT::PackBuffer& data) const
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
}

void MAT::FluidPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
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
        params_ = static_cast<MAT::PAR::FluidPoro*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // Only execute if not in post-process mode
  if (params_ != nullptr) anisotropy_strategy_ = MAT::FLUIDPORO::CreateAnisotropyStrategy(params_);

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

double MAT::FluidPoro::ComputeReactionCoeff() const
{
  return anisotropy_strategy_->ComputeReactionCoeff();
}

void MAT::FluidPoro::ComputeReactionTensor(LINALG::Matrix<2, 2>& reaction_tensor, const double J,
    const double porosity,
    const std::vector<std::vector<double>>& anisotropic_permeability_directions) const
{
  anisotropy_strategy_->ComputeReactionTensor(
      reaction_tensor, J, porosity, anisotropic_permeability_directions);
}

void MAT::FluidPoro::ComputeReactionTensor(LINALG::Matrix<3, 3>& reaction_tensor, const double J,
    const double porosity,
    const std::vector<std::vector<double>>& anisotropic_permeability_directions) const
{
  anisotropy_strategy_->ComputeReactionTensor(
      reaction_tensor, J, porosity, anisotropic_permeability_directions);
}

void MAT::FluidPoro::ComputeLinMatReactionTensor(LINALG::Matrix<2, 2>& linreac_dphi,
    LINALG::Matrix<2, 2>& linreac_dJ, const double J, const double porosity) const
{
  anisotropy_strategy_->ComputeLinMatReactionTensor(linreac_dphi, linreac_dJ, J, porosity);
}

void MAT::FluidPoro::ComputeLinMatReactionTensor(LINALG::Matrix<3, 3>& linreac_dphi,
    LINALG::Matrix<3, 3>& linreac_dJ, const double J, const double porosity) const
{
  anisotropy_strategy_->ComputeLinMatReactionTensor(linreac_dphi, linreac_dJ, J, porosity);
}

/*----------------------------------------------------------------------*
 *                                                           vuong 06/11 |
 *----------------------------------------------------------------------*/
double MAT::FluidPoro::EffectiveViscosity() const
{
  // set zero viscosity and only modify it for Darcy-Stokes problems
  double viscosity = -1.0;
  if (Type() == PAR::darcy)
    viscosity = 0.0;
  else if (Type() == PAR::darcy_brinkman)
    viscosity = Viscosity();
  else
    dserror("Unknown flow type for porous flow");

  return viscosity;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)     vuong 05/12|
 *----------------------------------------------------------------------*/
void MAT::FluidPoro::EvaluateViscStress(LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
    const LINALG::Matrix<6, 1>* glstrain, const int gp, Teuchos::ParameterList& params)

{
  dserror("macroscopic viscous stress not yet implemented for poroelasticity");
  return;
}
