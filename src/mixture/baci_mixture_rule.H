/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a base mixture rule

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_RULE_H
#define BACI_MIXTURE_RULE_H

#include <memory>
#include <unordered_map>

#include "baci_mat_par_parameter.H"
#include <Teuchos_ENull.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <string>
#include <utility>
#include <vector>
#include "baci_utils_exceptions.H"
#include "baci_linalg_fixedsizematrix.H"

// forward declarations
namespace LINALG
{
  class SerialDenseMatrix;
}
namespace DRT
{
  class PackBuffer;
  namespace INPUT
  {
    class LineDefinition;
  }
}  // namespace DRT
namespace MAT
{
  class Anisotropy;
  class Material;
  namespace PAR
  {
    class Material;
  }
}  // namespace MAT
namespace Teuchos
{
  class ParameterList;
}

namespace MIXTURE
{
  // forward declaration
  class MixtureRule;
  class MixtureConstituent;

  namespace PAR
  {
    class MixtureRule : public MAT::PAR::Parameter
    {
      friend class MIXTURE::MixtureRule;

     public:
      /// constructor
      explicit MixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// Override this method and throw error, as only the CreateRule() should be used.
      Teuchos::RCP<MAT::Material> CreateMaterial() final
      {
        dserror("Cannot create mixture rule from this method. Use CreateRule() instead.");
        return Teuchos::null;
      }

      /// create material instance of matching type with my parameters
      virtual std::unique_ptr<MIXTURE::MixtureRule> CreateRule() = 0;

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
    };
  }  // namespace PAR

  /*!
   * \brief Mixture rule containing the physics behind the mixture
   *
   * This class should be used within the Mixture framework as a mixture rule. This class
   * contains the whole physics. This is the base class defining the simplest possible physics,
   * i.e. constituents all deforming with the same deformation gradient and a homogenized stress
   * response using the mass density of each constituent.
   *
   * Example input lines:
   * MAT 1 MAT_Mixture NUMCONST 2 MATIDSCONST 11 12 MATIDMIXTURERULE 10
   * MAT 10 MIX_Rule_Simple DENS 0.1 NUMCONST 2 MASSFRAC 0.4 0.6
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

    virtual ~MixtureRule() = default;

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
        std::shared_ptr<std::vector<std::unique_ptr<MIXTURE::MixtureConstituent>>> constituents)
    {
      constituents_ = std::move(constituents);
    }

    virtual void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
    {
      // do nothing in the default case
    }

    /*!
     * Initialize the mixturerule with the element parameters of the input line
     *
     * @param numgp (in) Number of Gauss-points
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
    virtual void Update(CORE::LINALG::Matrix<3, 3> const& F, Teuchos::ParameterList& params,
        const int gp, const int eleGID)
    {
      // Nothing needs to be updated in this simple mixture rule
    }

    /*!
     * \brief Update of prestretch during the internal prestretching mechanism
     *
     * \param defgrd Deformation gradient that ensures global equilibrium in the deformed
     * configuration
     *
     * \param params Container for additional information
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    virtual void UpdatePrestress(CORE::LINALG::Matrix<3, 3> const& defgrd,
        Teuchos::ParameterList& params, const int gp, const int eleGID)
    {
      // Do nothing in the default case
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
     * @param gp (in) : Gauss point
     * @param eleGID (in) : Global element id
     */
    virtual void Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
        const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp,
        int eleGID) = 0;

    /*!
     * @brief Returns the material mass density
     *
     * @return material mass density
     */
    virtual double ReturnMassDensity() const
    {
      dserror("Rule does not provide the evaluation of a material mass density.");
      return 0;
    }

    /*!
     * \brief Register names of the internal data that should be saved during runtime vtk output
     *
     * \param name_and_size [out] : unordered map of names of the data with the respective vector
     * size
     */
    virtual void RegisterVtkOutputDataNames(
        std::unordered_map<std::string, int>& names_and_size) const
    {
      // do nothing for simple mixture rules
    }

    /*!
     * \brief Evaluate internal data for every Gauss point saved for output during runtime vtk
     * output
     *
     * \param name [in] : Name of the data to export
     * \param data [out] : NUMGPxNUMDATA Matrix holding the data
     *
     * \return true if data is set by the material, otherwise false
     */
    virtual bool EvaluateVtkOutputData(
        const std::string& name, CORE::LINALG::SerialDenseMatrix& data) const
    {
      return false;
    }

   protected:
    /*!
     * \brief Returns a reference to the constituents
     * @return
     */
    std::vector<std::unique_ptr<MIXTURE::MixtureConstituent>>& Constituents() const
    {
      return *constituents_;
    }

    /*!
     * Get number of Gauss points used
     *
     * @return Number of Gauss points
     */
    int NumGP() const { return numgp_; }

   private:
    ///! list of the references to the constituents
    std::shared_ptr<std::vector<std::unique_ptr<MIXTURE::MixtureConstituent>>> constituents_;

    ///! Number of Gauss points
    int numgp_;

    ///! Indicator, whether the constituent has already read the element definition
    bool has_read_element_;

    ///! Indicator, whether the constituent is already set up
    bool is_setup_;
  };
}  // namespace MIXTURE

#endif  // MIXTURE_RULE_H
