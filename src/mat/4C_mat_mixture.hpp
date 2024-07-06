/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of the mixture material holding a general mixturerule and mixture constituents

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MIXTURE_HPP
#define FOUR_C_MAT_MIXTURE_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_matelast_summand.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_mixture_constituent.hpp"
#include "4C_mixture_rule.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  // forward declaration
  class Mixture;

  namespace PAR
  {
    class Mixture : public Core::Mat::PAR::Parameter
    {
      friend class Mat::Mixture;

     public:
      /// Standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      explicit Mixture(const Core::Mat::PAR::Parameter::Data& matdata);


      /// @name material parameters
      /// @{
      /// list of the references to the constituents
      std::vector<MIXTURE::PAR::MixtureConstituent*> constituents_;

      /// rule of the mixture (contains the physics)
      MIXTURE::PAR::MixtureRule* mixture_rule_;

      /// @}

      Teuchos::RCP<Core::Mat::Material> create_material() override;
    };
  }  // namespace PAR


  class MixtureType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "MixtureType"; }

    static MixtureType& instance() { return instance_; }

    Core::Communication::ParObject* create(const std::vector<char>& data) override;

   private:
    static MixtureType instance_;
  };

  /*!
   * \brief Material class holding a general mixture material
   *
   * This class has to be paired with a mixture rule (MIXTURE::MixtureRule) containing the physics
   * and constituents (MIXTURE::MixtureConstituent). This class manages the interplay between
   * the mixture rule and the constituents defining an clear interface if an extension of the
   * mixture framework is needed.
   *
   * Example input lines:
   * MAT 1 MAT_Mixture NUMCONST 2 MATIDSCONST 11 12 MATIDMIXTURERULE 10
   *  DENS 0.1
   * MAT 10 MIX_Rule_Base NUMCONST 2 MASSFRAC 0.4 0.6
   * MAT 11 MIX_Constituent_ElastHyper NUMMAT 1 MATIDS 101 MAT 12
   * MAT 101 ELAST_CoupLogNeoHooke MODE YN C1 2.5e4 C2 0.27
   * MAT 102 ELAST_CoupAnisoExpo K1 1.666666666666e4 K2 10.0 GAMMA 0.0 K1COMP
   *  0.833333333333e4 K2COMP 10.0 ADAPT_ANGLE No INIT 0 STR_TENS_ID 1000
   * MAT 1000 ELAST_StructuralTensor STRATEGY Standard
   *
   */
  class Mixture : public So3Material
  {
   public:
    /// Constructor for an empty material object
    Mixture();

    /// Constructor for the material given the material parameters
    explicit Mixture(Mat::PAR::Mixture* params);


    /// @name Packing and Unpacking
    /// @{

    /// \brief Return unique ParObject id
    int unique_par_object_id() const override
    {
      return MixtureType::instance().unique_par_object_id();
    }

    /*!
     * \brief Pack this class so it can be communicated
     *
     * Resizes the vector data and stores all information of a class in it. The first information
     * to be stored in data has to be the unique parobject id delivered by unique_par_object_id()
     * which will then identify the exact class on the receiving processor.
     *
     * @param data (in/out): char vector to store class information
     */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
     * \brief Unpack data from a char vector into this class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of this file and delivered by UniqueParObjetId().
     *
     * @param data (in) : vector storing all data to be unpacked into this instance
     */
    void unpack(const std::vector<char>& data) override;

    /// @)

    /// check if element kinematics and material kinematics are compatible
    void valid_kinematics(Inpar::Solid::KinemType kinem) override
    {
      if (!(kinem == Inpar::Solid::KinemType::nonlinearTotLag))
      {
        FOUR_C_THROW(
            "element and material kinematics are not compatible. Use Nonlinear total lagrangian"
            "kinematics (KINEM nonlinear) in your element definition.");
      }
    }

    /// Return material type
    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_mixture;
    }

    /// Create a copy of this material
    /// \return copy of this material
    Teuchos::RCP<Core::Mat::Material> clone() const override
    {
      return Teuchos::rcp(new Mixture(*this));
    }


    /*!
     * \brief Quick access to the material parameters
     *
     * @return Material parameters
     */
    Core::Mat::PAR::Parameter* parameter() const override { return params_; }

    /*!
     * \brief Setup of the material (Read the input line definition of the element)
     *
     * This method will be called during reading of the elements. Here, we create all local tensors.
     *
     * @param numgp Number of Gauss-points
     * @param linedef Line definition of the element
     */
    void setup(int numgp, Input::LineDefinition* linedef) override;

    /*!
     * \brief Post setup routine that will be called before the first Evaluate call
     *
     * @param params Container for additional information
     */
    void post_setup(Teuchos::ParameterList& params, const int eleGID) override;

    /*!
     * \brief Update of the material
     *
     * This method will be called between each timestep.
     *
     * @param defgrd Deformation gradient
     * @param gp Gauss point
     * @param params Container for additional information
     * @param eleGID Global element id
     */
    void update(Core::LinAlg::Matrix<3, 3> const& defgrd, int gp, Teuchos::ParameterList& params,
        int eleGID) override;

    /// \brief This material law uses the extended update method
    bool uses_extended_update() override { return true; }

    /*!
     * \brief Evaluation of the material
     *
     * This method will compute the 2. Piola-Kirchhoff stress and the linearization of the material.
     *
     * @param defgrd (in) Deformation gradient
     * @param glstrain (in) Green-Lagrange strain in strain-like Voigt notation
     * @param params (in/out) Parameter list for additional information
     * @param stress (out) 2nd Piola-Kirchhoff stress tensor in stress-like Voigt notation
     * @param cmat (out) Linearization of the material law in Voigt notation
     * @param eleGID (in) Global element id
     */
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
        const int eleGID) final;

    /// \brief Return material mass density given by mixture rule
    double density() const override { return mixture_rule_->return_mass_density(); };

    void register_output_data_names(
        std::unordered_map<std::string, int>& names_and_size) const override;

    bool evaluate_output_data(
        const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const override;

   private:
    /// Material parameters
    Mat::PAR::Mixture* params_;

    /// list of the references to the constituents
    std::shared_ptr<std::vector<std::unique_ptr<MIXTURE::MixtureConstituent>>> constituents_;

    /// Reference to the mixturerule
    std::shared_ptr<MIXTURE::MixtureRule> mixture_rule_;

    /// Flag whether constituents are already set up.
    bool setup_;

    /// Flag for each gauss point whether the mixture rule and constituents are pre-evaluated
    std::vector<bool> is_pre_evaluated_;

    /// Holder for anisotropic materials
    Anisotropy anisotropy_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
