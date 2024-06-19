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
#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_membrane_material_interfaces.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*
 | active strain membrane material                 brandstaeter 05/2018 |
 *----------------------------------------------------------------------*/
namespace Mat
{
  // forward declaration
  class MembraneActiveStrain;

  namespace PAR
  {
    class MembraneActiveStrain : public Core::Mat::PAR::Parameter
    {
      friend class Mat::MembraneActiveStrain;

     public:
      /// standard constructor
      MembraneActiveStrain(const Core::Mat::PAR::Parameter::Data& matdata);

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
      Teuchos::RCP<Core::Mat::Material> create_material() override;
    };
    // class MembraneActiveStrain

  }  // namespace PAR

  class MembraneActiveStrainType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "Membrane_ActiveStrainType"; }

    static MembraneActiveStrainType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MembraneActiveStrainType instance_;
  };
  // class Membrane_ActiveStrainType

  /*----------------------------------------------------------------------*/
  /// \author brandstaeter
  /// \date 05/18
  // forward declaration
  class Material;

  class MembraneActiveStrain : public So3Material, public Mat::MembraneMaterialLocalCoordinates
  {
   public:
    /// construct empty material object
    MembraneActiveStrain();

    /// construct the material object given material parameters
    explicit MembraneActiveStrain(Mat::PAR::MembraneActiveStrain* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return MembraneActiveStrainType::Instance().UniqueParObjectId();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    void pack(Core::Communication::PackBuffer& data) const override;

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
    void unpack(const std::vector<char>& data) override;

    //@}

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(Inpar::STR::KinemType kinem) override
    {
      if (!(kinem == Inpar::STR::KinemType::nonlinearTotLag))
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    /// material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_membrane_activestrain;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new MembraneActiveStrain(*this));
    }

    /// material mass density
    double Density() const override { return params_->density_; }

    /// setup
    void setup(int numgp, Input::LineDefinition* linedef) override;

    /// Standard SO3 evaluate (not meant to be used)
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const Core::LinAlg::Matrix<6, 1>* glstrain,          ///< Green-Lagrange strain
        Teuchos::ParameterList& params,      ///< Container for additional information
        Core::LinAlg::Matrix<6, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses
        Core::LinAlg::Matrix<6, 6>* cmat,    ///< Constitutive matrix
        int gp,                              ///< Gauss point
        int eleGID) override                 ///< Element ID
    {
      FOUR_C_THROW("This a membrane material. Calling So3 evaluate does not make sense.");
    };

    void UpdateMembrane(const Core::LinAlg::Matrix<3, 3>& defgrd, Teuchos::ParameterList& params,
        const Core::LinAlg::Matrix<3, 3>& Q_trafo, int gp, int eleGID) override
    {
      // nothing to do
    }

    void EvaluateMembrane(const Core::LinAlg::Matrix<3, 3>& defgrd,
        const Core::LinAlg::Matrix<3, 3>& cauchygreen, Teuchos::ParameterList& params,
        const Core::LinAlg::Matrix<3, 3>& Q_trafo, Core::LinAlg::Matrix<3, 1>& stress,
        Core::LinAlg::Matrix<3, 3>& cmat, int gp, int eleGID) override;

    /// Update internal variables
    void update() override;

    /// Reset internal variables
    void reset_step() override;

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

   private:
    /// My material parameters
    Mat::PAR::MembraneActiveStrain* params_;

    /// passive material
    Teuchos::RCP<Mat::So3Material> matpassive_;

    /// (tansmembrane) voltage at every gp
    Teuchos::RCP<std::vector<double>> voltage_;

    /// activation parameter at every gp
    Teuchos::RCP<std::vector<double>> activation_;

    /// indicates if material is initialized
    bool isinit_;

    // setup fiber vectors
    void setup_fiber_vectors(int numgp, Input::LineDefinition* linedef);

    // read RAD-AXI-CIR
    void read_dir(
        Input::LineDefinition* linedef, std::string specifier, Core::LinAlg::Matrix<3, 1>& dir);

    // calculate normal direction from FIBER1 and FIBER2
    void setup_normal_direction();

    // pullback of the tangent from intermediate to reference configuration
    void pullback4th_tensor_voigt(const Core::LinAlg::Matrix<2, 2>& defgrd_active_inv_red,
        const Core::LinAlg::Matrix<3, 3>& cmat_passive_intermediate,
        Core::LinAlg::Matrix<3, 3>& cmat_reference);

    // transform voigt to tensor notation
    void tensor2x2_indices(int p, int* i, int* j);

    // transform tensor to voigt notation
    void voigt3_index(int i, int j, int* p);

   protected:
    /// vector of fiber vectors
    std::vector<Core::LinAlg::Matrix<3, 1>> fibervecs_;
  };  // class MembraneActiveStrain

}  // namespace Mat
// namespace Mat
FOUR_C_NAMESPACE_CLOSE

#endif
