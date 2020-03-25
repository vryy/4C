/*----------------------------------------------------------------------*/
/*! \file

\brief This file holds the definition of a mixture constituent interface

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_H
#define BACI_MIXTURE_CONSTITUENT_H

#include "../drt_mat/matpar_parameter.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_pack_buffer.H"
#include "../drt_lib/drt_linedefinition.H"

// forward declaration
namespace MAT
{
  class Anisotropy;
}

/*!
 * \brief The mixture namespace holds all mixture specific classes.
 *
 * The idea behind mixtures is that multiple materials share the same
 * deformation and the stress response is a mass fraction weighted
 * average of the stresses of each constituent.
 */
namespace MIXTURE
{
  // forward declaration
  class MixtureConstituent;
  class MixtureRule;

  namespace PAR
  {
    class MixtureConstituent : public MAT::PAR::Parameter
    {
     public:
      explicit MixtureConstituent(
          const Teuchos::RCP<MAT::PAR::Material>& matdata, double ref_mass_fraction);

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MIXTURE::MixtureConstituent> CreateConstituent() = 0;

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() final;

      static MIXTURE::PAR::MixtureConstituent* Factory(int matnum, double ref_mass_fraction);

      double RefMassFraction() const { return ref_mass_fraction_; }

     private:
      const double ref_mass_fraction_;
    };
  }  // namespace PAR

  /*!
   * \brief This is the base class of a constituent in a mixture defining the interface to the
   * holder class
   *
   * This abstract class defines the interface that a constituents needs to implement. It has to be
   * paired with MAT::Mixture_ElastHyper and MIXTURE::MixtureRule.
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
  class MixtureConstituent
  {
   public:
    MixtureConstituent();


    /*!
     * \brief Pack data into a char vector from this class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of the file and delivered by UniqueParObjectId().
     *
     * @param data (in/put) : vector storing all data to be packed into this instance.
     */
    virtual void PackConstituent(DRT::PackBuffer& data) const;

    /*!
     * \brief Unpack data from a char vector into this class to be called from a derived class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of the file and delivered by UniqueParObjectId().
     *
     * @param position (in/out) : current position to unpack data
     * @param data (in) : vector storing all data to be unpacked into this instance.
     */
    virtual void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data);

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

    /*!
     * \brief Register anisotropy extensions of all sub-materials of the constituent
     *
     * \param anisotropy Reference to the global anisotropy manager
     */
    virtual void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
    {
      // do nothing in the default case
    }

    /*!
     * Initialize the constituent with the parameters of the input line
     *
     * @param numgp (in) Number of Gauß-points
     * @param params (in/out) Parameter list for exchange of parameters
     */
    virtual void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef);

    /*!
     * Set the mixture parameters
     * @param mixtureParams parameter of the mixture
     */
    void SetInitialReferenceDensity(const double initialReferenceDensity)
    {
      initialReferenceDensity_ = initialReferenceDensity;
    }

    /*!
     * Returns whether the constituent is already set up
     * @return true if the constituent is already set up, otherwise false
     */
    bool IsSetup() { return is_setup_; }

    /*!
     * \brief Setup the constituent
     *
     * This method is called once for Gauß point at the beginning of the simulation.
     * The constituent should setup all internal variables and materials.
     *
     * @param params Container for additional information
     * @param eleGID Global element id
     */
    virtual void Setup(Teuchos::ParameterList& params, int eleGID);

    /*!
     * \brief Update of the internal variables
     *
     * This method is called once per Gauß point between each time step to update the internal
     * variables. (Not needed for simple constituents)
     *
     * @param defgrd Deformation gradient of the previous timestep
     * @param params Container for additional information
     * @param gp Gauß point
     * @param eleGID Global element identifier
     */
    virtual void Update(LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params,
        const int gp, const int eleGID)
    {
    }

    /*!
     * Evaluates the constituents. Needs to compute the stress contribution of the constituent out
     * of the displacements. Will be called for each Gauß point
     *
     * @param F (in) : Deformation gradient
     * @param E_strain (in) : Green-Lagrange strain in strain-like Voigt notation
     * @param params Container for additional information
     * @param S_stress (out) : 2nd Piola Kirchhoff stress tensor in stress like Voigt-notation
     * @param cmat Constitutive tensor in Voigt notation
     * @param gp (in) : Number of Gauß point
     * @param eleGID (in) : Global element id
     */
    virtual void Evaluate(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) = 0;

    /// Returns the refenrece mass density. Needs to be implemented by the deriving class.
    virtual double CurrentRefDensity() const = 0;

    /*!
     * \brief Get the names and dimension of the quantities to visualize
     *
     * @param names (out) Tuple of names and the corresponding dimension of the quantity
     */
    virtual void VisNames(std::map<std::string, int>& names) {}

    /*!
     * \brief
     * @param name (in) Name of the quantity
     * @param data (out) Values of the quantity
     * @param numgp (in) number of Gauß points
     * @param eleGID (in) Global element id
     */
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleGID)
    {
      return false;
    }

   protected:
    /*!
     * Get number of Gauß points used
     *
     * @return Number of Gauß points
     */
    int NumGP() const { return numgp_; }

    /*!
     * Returns the initial reference density of the consitituent
     * @return
     */
    double InitialRefDensity() const { return initialReferenceDensity_; }

   private:
    ///! Mixture parameters
    double initialReferenceDensity_;

    ///! Number of Gauß points
    int numgp_;

    ///! Indicator, whether the constituent has already read the element
    bool has_read_element_;

    ///! Indicator, whether the constituent is already set up
    bool is_setup_;
  };

}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_H
