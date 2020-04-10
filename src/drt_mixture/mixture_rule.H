/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a base mixture rule

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_RULE_H
#define BACI_MIXTURE_RULE_H

#include "../drt_mat/matpar_parameter.H"
#include "../drt_lib/drt_linedefinition.H"
#include "mixture_constituent.H"

namespace MIXTURE
{
  // forward declaration
  class MixtureRule;

  namespace PAR
  {
    class MixtureRule : public MAT::PAR::Parameter
    {
      friend class MIXTURE::MixtureRule;

     public:
      /// constructor
      explicit MixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// default destructor
      ~MixtureRule() override = default;

      /// Override this method and throw error, as only the CreateRule() should be used.
      Teuchos::RCP<MAT::Material> CreateMaterial() final
      {
        dserror("Cannot create mixture rule from this method. Use CreateRule() instead.");
        return Teuchos::null;
      }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MIXTURE::MixtureRule> CreateRule();

      /*!
       * \brief Factory of the mixture rule parameters
       *
       * This static method generates the specific class of the mixture rule defined in the datfile
       * at the corresponding material id
       *
       * @param matid Material id of the mixturerule
       * @return Parameters of the referenced mixture rule
       */
      static MIXTURE::PAR::MixtureRule* Factory(int matid);

      /// @name parameters of the mixture rule
      /// @{
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief Mixture rule containing the physics behind the mixture
   *
   * This class should be used within the Mixture_ElastHyper framework as a mixture rule. This class
   * contains the whole physics. This is the base class defining the simplest possible physics,
   * i.e. constituents all deforming with the same deformation gradient and a homogenized stress
   * response using the mass density of each constituent.
   *
   * Example input lines:
   * MAT 1 MAT_MixtureElastHyper NUMCONST 2 MATIDSCONST 11 12 MASSFRAC 0.4 0.6 MATIDMIXTURERULE 10
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
  class MixtureRule
  {
   public:
    /// Constructor for the material given the material parameters
    explicit MixtureRule(MIXTURE::PAR::MixtureRule* params);

    virtual void PackMixtureRule(DRT::PackBuffer& data) const;

    /*!
     * \brief Unpack data from a char vector into this class to be called from a derived class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of the file and delivered by UniqueParObjectId().
     *
     * @param data (in) : vector storing all data to be unpacked into this instance.
     * @param position (in/out) : current position to unpack data
     */
    virtual void UnpackMixtureRule(
        std::vector<char>::size_type& position, const std::vector<char>& data);

    /*!
     * This method should be called after creation of the constituents
     *
     * @param constituents (in) List of constituents
     */
    void SetConstituents(
        const Teuchos::RCP<std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>>& constituents)
    {
      constituents_ = constituents;
    }

    /*!
     * @brief This method sets the material mass denisty as part of setting up the rule
     * @param density material mass density as defined in the DAT file
     */
    void SetMaterialMassDensity(double density);

    /*!
     * @brief Gets the material mass density as defined in the dat file
     * @return Material mass density
     */
    double GetMaterialMassDensity() const;

    virtual void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
    {
      // do nothing in the default case
    }

    /*!
     * Initialize the mixturerule with the element parameters of the input line
     *
     * @param numgp (in) Number of Gauß-points
     * @param params (in/out) : Parameter list for exchange of parameters
     */
    virtual void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef);

    /*!
     * Returns whether the constituent is already set up
     * @return true if the constituent is already set up, otherwise false
     */
    bool IsSetup() { return is_setup_; }

    /*!
     * \brief Setup the constituent
     *
     * This method is called once for each element. The constituent should setup all his internal
     * variables and materials.
     *
     * @param params (in/out) : Container for additional information
     * @param eleGID (in) : global element id
     */
    virtual void Setup(Teuchos::ParameterList& params, int eleGID);

    /*!
     * \brief Update of the material law
     *
     * This simple mixture rule does not need to update anything, so this method is kept empty
     *
     * @param F (in) : Deformation gradient
     * @param params (in/out) : Container for additional information
     * @param gp (in) : Gauss point
     * @param eleGID (in) : Global element id
     */
    virtual void Update(LINALG::Matrix<3, 3> const& F, Teuchos::ParameterList& params, const int gp,
        const int eleGID)
    {
      // Nothing needs to be updated in this simple mixture rule
    }

    /*!
     *
     * @brief Method that is executed before the first evaluate call, once for each Gauss point
     *
     * @param params (in) : Container for additional information
     * @param gp (in) : Gauss point
     * @param eleGID (in) : Global element id
     */
    virtual void PreEvaluate(Teuchos::ParameterList& params, const int gp, const int eleGID)
    {
      // do nothing in the default case
    }

    /*!
     * Evaluates the constituents. Needs to compute the stress contribution of the constituent out
     * of the displacements. Will be called for each Gauss point
     *
     * @param F (in) : Deformation gradient
     * @param E_strain (in) : Green-Lagrange strain tensor in strain-like Voigt notation
     * @param params (in/out) : Container for additional parameters
     * @param S_stress (out) : 2nd Piola Kirchhoff stress tensor in stress like Voigt-notation
     * @param cmat (out) : Constitutive tensor in Voigt notation
     * @param gp (in) : Gauß point
     * @param eleGID (in) : Global element id
     */
    virtual void Evaluate(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID);

    /*!
     * \brief Get the names and dimension of the quantities to visualize
     *
     * @param names (out) : Tuple of names and the corresponding dimension of the quantity
     */
    virtual void VisNames(std::map<std::string, int>& names) {}

    /*!
     * \brief
     * @param name (in) : Name of the quantity
     * @param data (out) : Values of the quantity
     * @param numgp (in) : number of Gauß points
     * @param eleGID (in) : Global element id
     */
    virtual bool VisData(
        const std::string& name, std::vector<double>& data, const int numgp, const int eleGID)
    {
      return false;
    }

   protected:
    /*!
     * \brief Returns a reference to the constituents
     * @return
     */
    Teuchos::RCP<std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>> Constituents() const
    {
      return constituents_;
    }

    /*!
     * Get number of Gauß points used
     *
     * @return Number of Gauß points
     */
    int NumGP() const { return numgp_; }

   private:
    ///! list of the references to the constituents
    Teuchos::RCP<std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>> constituents_;

    ///! Number of Gauß points
    int numgp_;

    ///! Indicator, whether the constituent has already read the element definition
    bool has_read_element_;

    ///! Indicator, whether the constituent is already set up
    bool is_setup_;

    ///! Material mass density as defined in the input file
    double material_mass_density_;
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_RULE_H
