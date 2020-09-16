/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a hyperelastic constituent basis

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_ELASTHYPER_BASE_H
#define BACI_MIXTURE_CONSTITUENT_ELASTHYPER_BASE_H

#include "../drt_mat/matpar_parameter.H"
#include "mixture_constituent.H"
#include "../drt_mat/elasthyper_service.H"
#include "mixture_prestress_strategy.H"
#include "../drt_mat/anisotropy_extension_cylinder_cosy.H"

namespace MAT
{
  class Anisotropy;
}

namespace MIXTURE
{
  class MixtureConstituent_ElastHyperBase;

  namespace PAR
  {
    class MixtureConstituent_ElastHyperBase : public MIXTURE::PAR::MixtureConstituent
    {
     public:
      explicit MixtureConstituent_ElastHyperBase(
          const Teuchos::RCP<MAT::PAR::Material>& matdata, double ref_mass_fraction);

      ~MixtureConstituent_ElastHyperBase() override = default;

      int GetPrestressingMatId() const { return matid_prestress_strategy_; }

      Teuchos::RCP<MIXTURE::PrestressStrategy> PrestressStrategy() const
      {
        return prestressStrategy_;
      }

      /// @name material parameters
      /// @{

      /// Material id of the prestress strategy
      const int matid_prestress_strategy_;

      /// number of summands
      const int nummat_;

      /// List of material ids of the summands
      const std::vector<int>* matids_;

      /// Strategy for prestressing the constituent
      Teuchos::RCP<MIXTURE::PrestressStrategy> prestressStrategy_;
      /// @}
    };
  }  // namespace PAR

  /*!
   * \brief Constituent for any hyperelastic material
   *
   * This constituent represents any hyperelastic material from the elasthyper toolbox. It has to
   * be paired with the MAT::Mixture_ElastHyper material and a MIXTURE::MixtureRule.
   */
  class MixtureConstituent_ElastHyperBase : public MIXTURE::MixtureConstituent
  {
   public:
    /// Constructor for the materiak given the material parameters
    explicit MixtureConstituent_ElastHyperBase(
        MIXTURE::PAR::MixtureConstituent_ElastHyperBase* params);

    /*!
     * \brief Pack data into a char vector from this class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of the file and delivered by UniqueParObjectId().
     *
     * @param data (in/put) : vector storing all data to be packed into this instance.
     */
    void PackConstituent(DRT::PackBuffer& data) const override;

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
    void UnpackConstituent(
        std::vector<char>::size_type& position, const std::vector<char>& data) override;

    /*!
     * \brief Register all anisotropy extensions also for the sub-summands
     *
     * \param anisotropy Reference to the global anisotropy manager
     */
    void RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy) override;

    /*!
     * Initialize the constituent with the parameters of the input line
     *
     * @param numgp (in) Number of Gauß-points
     * @param params (in/out) Parameter list for exchange of parameters
     */
    void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    /*!
     * \brief Updates the material and all its summands
     *
     * This method is called once between each timestep after convergence.
     *
     * @param defgrd Deformation gradient
     * @param params Container for additional information
     * @param gp Gauß point
     * @param eleGID Global element identifier
     */
    void Update(LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    void UpdatePrestress(LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params, int gp,
        int eleGID) override;

    /*!
     * \brief Returns the current reference mass density of the constituent
     *
     * @param gp (in) : Gauss point
     * @return
     */
    double CurrentRefDensity(int gp) const override;

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

    /*!
     * \brief Returns a reference to all summands
     *
     * \return const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& Reference to the summands
     */
    const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& Summands() const { return potsum_; }

    /*!
     * \brief Returns a reference to all summand properties
     *
     * \return const MAT::SummandProperties& Reference to the summand properties
     */
    const MAT::SummandProperties& SummandProperties() { return summandProperties_; }

    /*!
     * \brief Method that is called to setup the constituent once before the start of the simulation
     *
     * \param params Container for additional information
     * \param eleGID Global element id
     */
    void Setup(Teuchos::ParameterList& params, int eleGID) override;

    /*!
     * \brief Method that is called once for each Gauss point before the first evaluate call
     *
     * \param mixtureRule Reference to the mixture rule
     * \param params Container for additional information
     * \param gp Gauss point
     * \param eleGID Global element id
     */
    void PreEvaluate(
        MixtureRule& mixtureRule, Teuchos::ParameterList& params, int gp, int eleGID) override;

   protected:
    /*!
     * \brief Returns a reference to the prestretch tensor at the Gauss point
     *
     * \param gp Gauss point
     * \return const LINALG::Matrix<3, 3>& Reference to the prestretch tensor
     */
    const LINALG::Matrix<3, 3>& PrestretchTensor(const int gp) const { return prestretch_[gp]; }

    /*!
     * \brief Returns a reference to the cylinder coordinate system
     *
     * \return const MAT::CylinderCoordinateSystemProvider&
     */
    const MAT::CylinderCoordinateSystemAnisotropyExtension&
    CylinderCoordinateSystemAnisotropyExtension() const
    {
      return cosyAnisotropyExtension_;
    }

   private:
    /// @name Flags to specify the elastic formulations (initialize with false)
    //@{
    MAT::SummandProperties summandProperties_;  ///< holder for formulation specification
    //@}

    /// my material parameters
    MIXTURE::PAR::MixtureConstituent_ElastHyperBase* params_;

    /// map to materials/potential summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsum_;

    /// Prestretch of the constituent
    std::vector<LINALG::Matrix<3, 3>> prestretch_;

    /// AnisotropyExtension that handles the management of cylinder coordinate systems
    MAT::CylinderCoordinateSystemAnisotropyExtension cosyAnisotropyExtension_;
  };

}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_ELASTHYPER_BASE_H
