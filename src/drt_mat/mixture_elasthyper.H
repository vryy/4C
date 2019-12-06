/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of the mixture material holding a general mixture law and mixture constituents

\level 3

\maintainer Amadeus Gebauer
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_ELASTHYPER_H
#define BACI_MIXTURE_ELASTHYPER_H


#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "so3_material.H"
#include "../drt_mixture/mixture_constituent.H"
#include "../drt_mixture/mixture_rule.H"

namespace MAT
{
  // forward declaraction
  class Mixture_ElastHyper;


  namespace PAR
  {
    class Mixture_ElastHyper : public Parameter
    {
      friend class MAT::Mixture_ElastHyper;

     public:
      /// Standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      explicit Mixture_ElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      ~Mixture_ElastHyper() override = default;

      /// @name material parameters
      /// @{

      /// list of mass fractions in the reference configuration at the beginning of the simulation
      const std::vector<double>* mass_fractions_;

      /// material mass density
      const double density_;

      /// list of the references to the constituents
      std::vector<MIXTURE::PAR::MixtureConstituent*> constituents_;

      /// rule of the mixture (contains the physics)
      MIXTURE::PAR::MixtureRule* mixture_rule_;

      /// @}

      Teuchos::RCP<MAT::Material> CreateMaterial() override;
    };
  }  // namespace PAR


  class Mixture_ElastHyperType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "Mixture_ElastHyperType"; }

    static Mixture_ElastHyperType& Instance() { return instance_; }

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static Mixture_ElastHyperType instance_;
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
   * MAT 1 MAT_MixtureElastHyper NUMCONST 2 MATIDSCONST 11 12 MASSFRAC 0.4 0.6 MATIDMIXTURELAW 10
   *  DENS 0.1
   * MAT 10 MIX_Rule_Base
   * MAT 11 MIX_Constituent_ElastHyper NUMMAT 1 MATIDS 101 MAT 12
   * MIX_Constituent_ElastHyper NUMMAT 1 MATIDS 102
   * MAT 101 ELAST_CoupLogNeoHooke MODE YN C1 2.5e4 C2 0.27
   * MAT 102 ELAST_CoupAnisoExpo K1 1.666666666666e4 K2 10.0 GAMMA 0.0 K1COMP
   *  0.833333333333e4 K2COMP 10.0 ADAPT_ANGLE No INIT 0 STR_TENS_ID 1000
   * MAT 1000 ELAST_StructuralTensor STRATEGY Standard
   *
   */
  class Mixture_ElastHyper : public So3Material
  {
   public:
    /// Constructor for an empty material object
    Mixture_ElastHyper();

    /// Constructor for the materiak given the material parameters
    explicit Mixture_ElastHyper(MAT::PAR::Mixture_ElastHyper* params);


    /// @name Packing and Unpacking
    /// @{

    /// \brief Return unique ParObject id
    int UniqueParObjectId() const override
    {
      return Mixture_ElastHyperType::Instance().UniqueParObjectId();
    }

    /*!
     * \brief Pack this class so it can be communicated
     *
     * Resizes the vector data and stores all information of a class in it. The first information
     * to be stored in data has to be the unique parobject id delivered by UniqueParObjectId() which
     * will then identify the exact class on the receiving processor.
     *
     * @param data (in/out): char vector to store class information
     */
    void Pack(DRT::PackBuffer& data) const override;

    /*!
     * \brief Unpack data from a char vector into this class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of this file and delivered by UniqueParObjetId().
     *
     * @param data (in) : vector storing all data to be unpacked into this instance
     */
    void Unpack(const std::vector<char>& data) override;

    /// @)

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
      {
        dserror(
            "element and material kinematics are not compatible. Use Nonlinear total lagrangian"
            "kinematics (KINEM nonlinear) in your element definition.");
      }
    }

    /// Return material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_mixture_elasthyper;
    }

    /// Create a copy of this material
    /// \return copy of this material
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new Mixture_ElastHyper(*this));
    }


    /*!
     * \brief Quick access to the material parameters
     *
     * @return Material parameters
     */
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /*!
     * \brief Return material mass density
     *
     * @return material mass density
     */
    double Density() const override { return params_->density_; }

    /*!
     * \brief Setup of the material (Read the input line definition of the element)
     *
     * This method will be called during reading of the elements. Here, we create all local tensors.
     *
     * @param numgp Number of Gauß-points
     * @param linedef Line definition of the element
     */
    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    /*!
     * \brief Post setup routine that will be called before the first Evaluate call
     *
     * @param params Container for additional information
     */
    void PostSetup(Teuchos::ParameterList& params) override;

    /*!
     * \brief Update of the material
     *
     * This method will be called between each timestep.
     *
     * @param defgrd Deformation gradient
     * @param gp Gauß point
     * @param params Container for additional information
     * @param eleGID Global element id
     */
    void Update(LINALG::Matrix<3, 3> const& defgrd, int gp, Teuchos::ParameterList& params,
        int eleGID) override;

    /// \brief This material law uses the extended update method
    bool UsesExtendedUpdate() override { return true; }

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
    void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        int eleGID) final;

    /*!
     * \brief Get the names and dimension of the quantities to visualize
     *
     * @param names (out) Tuple of names and the corresponding dimension of the quantity
     */
    void VisNames(std::map<std::string, int>& names) override;

    /*!
     * \brief
     * @param name (in) Name of the quantity
     * @param data (out) Values of the quantity
     * @param numgp (in) number of Gauß points
     * @param eleGID (in) Global element id
     */
    bool VisData(
        const std::string& name, std::vector<double>& data, int numgp, int eleGID) override;

   private:
    /// Material parameters
    MAT::PAR::Mixture_ElastHyper* params_;

    /// list of the references to the constituents
    Teuchos::RCP<std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>> constituents_;

    /// Reference to the mixture law
    Teuchos::RCP<MIXTURE::MixtureRule> mixture_rule_;

    /// Flag whether constituents are already set up.
    bool setup_;
  };

}  // namespace MAT

#endif  // BACI_MIXTURE_ELASTHYPER_H
