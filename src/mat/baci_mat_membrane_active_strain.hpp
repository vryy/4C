/*----------------------------------------------------------------------*/
/*! \file
 \brief Active strain membrane material for gastric electromechanics

 The input line should read
 MAT 0 MAT_Membrane_ActiveStrain

 \level 3



 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                     brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_MEMBRANE_ACTIVE_STRAIN_HPP
#define FOUR_C_MAT_MEMBRANE_ACTIVE_STRAIN_HPP

/*----------------------------------------------------------------------*
 | headers                                         brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_membrane_material_interfaces.hpp"
#include "baci_mat_par_parameter.hpp"
#include "baci_mat_so3_material.hpp"

BACI_NAMESPACE_OPEN
/*----------------------------------------------------------------------*
 | active strain membrane material                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
namespace MAT
{
  // forward declaration
  class Membrane_ActiveStrain;

  namespace PAR
  {
    class Membrane_ActiveStrain : public Parameter
    {
      friend class MAT::Membrane_ActiveStrain;

     public:
      /// standard constructor
      Membrane_ActiveStrain(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// Number of the material that describes the elastic behavior
      const int matid_passive_;

      /// Position of the transmembrane voltage of the SMC cell
      const int scalid_voltage_;

      /// density
      const double density_;

      /// beta (parameter dynamics of the VDCC)
      const double beta1_;

      /// beta (parameter dynamics of the Ca-ions)
      const double beta2_;

      /// voltage level for the activation of the Ca-influx
      const double voltage_threshold_;

      /// alpha (parameter for the intensity of the contraction in fiber 1 direction)
      const double alpha1_;

      /// alpha (parameter for the intensity of the contraction in fiber 2 direction)
      const double alpha2_;

      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;
    };
    // class Membrane_ActiveStrain

  }  // namespace PAR

  class Membrane_ActiveStrainType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "Membrane_ActiveStrainType"; }

    static Membrane_ActiveStrainType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Membrane_ActiveStrainType instance_;
  };
  // class Membrane_ActiveStrainType

  /*----------------------------------------------------------------------*/
  /// \author brandstaeter
  /// \date 05/18
  // forward declaration
  class Material;

  class Membrane_ActiveStrain : public So3Material, public MAT::MembraneMaterialLocalCoordinates
  {
   public:
    /// construct empty material object
    Membrane_ActiveStrain();

    /// construct the material object given material parameters
    explicit Membrane_ActiveStrain(MAT::PAR::Membrane_ActiveStrain* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return Membrane_ActiveStrainType::Instance().UniqueParObjectId();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::KinemType::nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_membrane_activestrain;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new Membrane_ActiveStrain(*this));
    }

    /// material mass density
    double Density() const override { return params_->density_; }

    /// setup
    void Setup(int numgp, INPUT::LineDefinition* linedef) override;

    /// Standard SO3 evaluate (not meant to be used)
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const CORE::LINALG::Matrix<6, 1>* glstrain,          ///< Green-Lagrange strain
        Teuchos::ParameterList& params,      ///< Container for additional information
        CORE::LINALG::Matrix<6, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses
        CORE::LINALG::Matrix<6, 6>* cmat,    ///< Constitutive matrix
        int gp,                              ///< Gauss point
        int eleGID) override                 ///< Element ID
    {
      dserror("This a membrane material. Calling So3 evaluate does not make sense.");
    };

    void UpdateMembrane(const CORE::LINALG::Matrix<3, 3>& defgrd, Teuchos::ParameterList& params,
        const CORE::LINALG::Matrix<3, 3>& Q_trafo, int gp, int eleGID) override
    {
      // nothing to do
    }

    void EvaluateMembrane(const CORE::LINALG::Matrix<3, 3>& defgrd,
        const CORE::LINALG::Matrix<3, 3>& cauchygreen, Teuchos::ParameterList& params,
        const CORE::LINALG::Matrix<3, 3>& Q_trafo, CORE::LINALG::Matrix<3, 1>& stress,
        CORE::LINALG::Matrix<3, 3>& cmat, int gp, int eleGID) override;

    /// Update internal variables
    void Update() override;

    /// Reset internal variables
    void ResetStep() override;

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

   private:
    /// My material parameters
    MAT::PAR::Membrane_ActiveStrain* params_;

    /// passive material
    Teuchos::RCP<MAT::So3Material> matpassive_;

    /// (tansmembrane) voltage at every gp
    Teuchos::RCP<std::vector<double>> voltage_;

    /// activation parameter at every gp
    Teuchos::RCP<std::vector<double>> activation_;

    /// indicates if material is initialized
    bool isinit_;

    // setup fiber vectors
    void SetupFiberVectors(int numgp, INPUT::LineDefinition* linedef);

    // read RAD-AXI-CIR
    void ReadDir(
        INPUT::LineDefinition* linedef, std::string specifier, CORE::LINALG::Matrix<3, 1>& dir);

    // calculate normal direction from FIBER1 and FIBER2
    void SetupNormalDirection();

    // pullback of the tangent from intermediate to reference configuration
    void Pullback4thTensorVoigt(const CORE::LINALG::Matrix<2, 2>& defgrd_active_inv_red,
        const CORE::LINALG::Matrix<3, 3>& cmat_passive_intermediate,
        CORE::LINALG::Matrix<3, 3>& cmat_reference);

    // transform voigt to tensor notation
    void Tensor2x2Indices(int p, int* i, int* j);

    // transform tensor to voigt notation
    void Voigt3Index(int i, int j, int* p);

   protected:
    /// vector of fiber vectors
    std::vector<CORE::LINALG::Matrix<3, 1>> fibervecs_;
  };  // class Membrane_ActiveStrain

}  // namespace MAT
// namespace MAT
BACI_NAMESPACE_CLOSE

#endif  // MAT_MEMBRANE_ACTIVE_STRAIN_H
