/*----------------------------------------------------------------------*/
/*!

\brief Definition of a hyperelastic constituent

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_ELASTHYPER_H
#define BACI_MIXTURE_CONSTITUENT_ELASTHYPER_H

#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "mixture_constituent.H"
#include "../drt_matelast/elast_summand.H"
#include "../drt_mat/elasthyper_service.H"

namespace MIXTURE
{
  class MixtureConstituent_ElastHyper;

  namespace PAR
  {
    class MixtureConstituent_ElastHyper : public MIXTURE::PAR::MixtureConstituent
    {
      friend class MIXTURE::MixtureConstituent_ElastHyper;

     public:
      explicit MixtureConstituent_ElastHyper(
          const Teuchos::RCP<MAT::PAR::Material>& matdata, double ref_mass_fraction);

      ~MixtureConstituent_ElastHyper() override = default;

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MIXTURE::MixtureConstituent> CreateConstituent() override;

      /// @name material parameters
      /// @{

      /// number of summands
      const int nummat_;

      /// List of material ids of the summands
      const std::vector<int>* matids_;
      /// @}
    };
  }  // namespace PAR

  class MixtureConstituent_ElastHyperType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "MixtureConstituent_ElastHyperType"; }

    static MixtureConstituent_ElastHyperType& Instance() { return instance_; }

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MixtureConstituent_ElastHyperType instance_;
  };

  /*!
   * \brief Constituent for any hyperelastic material
   *
   * This constituent represents any hyperelastic material from the elasthyper toolbox. It has to
   * be paired with the MAT::Mixture_ElastHyper material and a MIXTURE::MixtureRule.
   */
  class MixtureConstituent_ElastHyper : public MIXTURE::MixtureConstituent
  {
   public:
    MixtureConstituent_ElastHyper();

    /// Constructor for the materiak given the material parameters
    explicit MixtureConstituent_ElastHyper(MIXTURE::PAR::MixtureConstituent_ElastHyper* params);

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

    /// Returns the material type enum
    INPAR::MAT::MaterialType MaterialType() const override;


    /*!
     * Initialize the constituent with the parameters of the input line
     *
     * @param numgp (in) Number of Gauß-points
     * @param params (in/out) Parameter list for exchange of parameters
     */
    void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    /*!
     * Evaluates the constituents. Needs to compute the stress contribution of the constituent out
     * of the displacements. Will be called for each Gauß point
     *
     * @param F Deformation gradient
     * @param E_strain Green-Lagrange strain in strain-like Voigt notation
     * @param params Container for additional information
     * @param S_stress 2nd Piola Kirchhoff stress tensor in stress like Voigt-notation
     * @param cmat Constitutive tensor in Voigt notation
     */
    void Evaluate(const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<6, 1>& E_strain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
        int gp, int eleGID) override;

    /*!
     * \brief Updates the material and all its summands
     *
     * This method is called once between each timestep after convergence.
     *
     * @param defgrd Deformation gradient
     * @param params Container for additional information
     * @param gp Gauß point
     */
    void Update(
        LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params, int gp) override;

    /*!
     * \brief Returns the current reference mass density of the constituent
     * @return
     */
    double CurrentRefDensity() const override;

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
    /// @name Flags to specify the elastic formulations (initialize with false)
    //@{
    MAT::SummandProperties summandProperties_;  ///< holder for formulation specification
    //@}

    /// my material parameters
    MIXTURE::PAR::MixtureConstituent_ElastHyper* params_;

    /// map to materials/potential summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsum_;
  };

}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_ELASTHYPER_H
