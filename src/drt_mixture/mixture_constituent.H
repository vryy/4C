/*----------------------------------------------------------------------*/
/*!

\brief This file holds the definition of a mixture constituent interface

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_CONSTITUENT_H
#define BACI_MIXTURE_CONSTITUENT_H
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"

/*!
 * \brief The mixture namespace holds all mixture specific classes.
 *
 * The idea behind mixtures is that multiple materials share the same
 * deformation and the stress response is a mass fraction weighted
 * average of the stresses of each constituent.
 */
namespace MIXTURE
{
  class MixtureConstituent;

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
   */
  class MixtureConstituent : public DRT::ParObject
  {
   public:
    MixtureConstituent();

    /*!
     * \brief Return unique ParObject id
     *
     * Every class implementing ParObject needs a unique id defined at the top of drt_parobject.H
     * and should return it in this method.
     *
     * @return Return unique ParObject id
     */
    int UniqueParObjectId() const override = 0;


    /*!
     * \brief Pack this class so it can be communicated
     *
     * Resizes the vector data and stores all infromation of a class in it. The first information
     * to be stored in the data has to be unique parobect id delivered by UniqueParObjectId() which
     * will then identify the exact class on the receiving processor.
     *
     * @param data (in/out) : char vector to store class information
     */
    void Pack(DRT::PackBuffer& data) const override;

    /*!
     * \brief Unpack data from a char vector into this class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of the file and delivered by UniqueParObjectId().
     *
     * @param data (in) : vector storing all data to be unpacked into this instance.
     */
    void Unpack(const std::vector<char>& data) override;

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
     * Initialize the constituent with the parameters of the input line
     *
     * @param numgp (in) Number of Gauß-points
     * @param params (in/out) Parameter list for exchange of parameters
     */
    virtual void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef);

    /*!
     * This method should be called after ReadElement to initialize all quantities. It is called
     * once at the beginning of the simulation for each element
     *
     * @param params Parameter list containing different parameters
     */
    virtual void Init(Teuchos::ParameterList& params);

    /*!
     * Returns whether the constituent is already initialized
     * @return true if the constituent is already initialized, otherwise false
     */
    bool IsInit() { return is_init_; }

    /*!
     * \brief Setup the constituent
     *
     * This method is called once for Gauß point at the beginning of the simulation.
     * The constituent should setup all internal variables and materials.
     * @param gp Number of the Gauß points
     */
    virtual void Setup(int gp, Teuchos::ParameterList& params);

    /*!
     * \brief Update of the internal variables
     *
     * This method is called once per Gauß point between each time step to update the internal
     * variables. (Not needed for simple constituents)
     *
     * @param defgrd Deformation gradient of the previous timestep
     * @param params Container for additional information
     * @param gp Gauß point
     */
    virtual void Update(
        LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params, const int gp)
    {
    }

    /*!
     * Evaluates the constituents. Needs to compute the stress contribution of the constituent out
     * of the displacements. Will be called for each Gauß point
     *
     * @param defgrd Deformation gradient
     * @param glstrain Green-Lagrange strain in strain-like Voigt notation
     * @param params Container for additional information
     * @param stress 2nd Piola Kirchhoff stress tensor in stress like Voigt-notation
     * @param cmat Constitutive tensor in Voigt notation
     */
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        int gp, int eleGID) = 0;

    /// Returns the refenrece mass fraction. Needs to be implemented by the deriving class.
    virtual double RefMassFraction() const = 0;

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
     * @param eleGID (in) element id
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

   private:
    ///! Number of Gauß points
    int numgp_;

    ///! Indicator, whether the constituent is already initialized
    bool is_init_;

    ///! Indicator, whether the constituent is already set up
    std::vector<char> is_setup_;
  };

}  // namespace MIXTURE

#endif  // BACI_MIXTURE_CONSTITUENT_H
