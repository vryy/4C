/*----------------------------------------------------------------------*/
/*!

\brief Definition of a base mixture rule

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MIXTURE_RULE_H
#define BACI_MIXTURE_RULE_H

#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_linedefinition.H"
#include "mixture_constituent.H"
#include "../drt_mat/matpar_parameter.H"
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
      Teuchos::RCP<MAT::Material> CreateMaterial() override
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
       * @param matid Material id of the mixture law
       * @return Parameters of the referenced mixture rule
       */
      static MIXTURE::PAR::MixtureRule* Factory(int matid);

      /// @name material parameters
      /// @{
      /// @}
    };
  }  // namespace PAR

  class MixtureRuleType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "MixtureRuleType"; }

    static MixtureRuleType& Instance() { return instance_; }

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static MixtureRuleType instance_;
  };

  /*!
   * \brief Mixture rule containing the physics behind the mixture
   *
   * This class should be used within the Mixture_ElastHyper framework as a mixture rule. This class
   * contains the whole physics. This is the base class defining the simplest possible physics,
   * i.e. constituents all deforming with the same deformation gradient and a homogenized stress
   * response using the mass density of each constituent.
   */
  class MixtureRule : public DRT::ParObject
  {
   public:
    /// Constructor for an empty material object
    MixtureRule();

    /// Constructor for the materiak given the material parameters
    explicit MixtureRule(MIXTURE::PAR::MixtureRule* params);

    /*!
     * \brief Return unique ParObject id
     *
     * Every class implementing ParObject needs a unique id defined at the top of drt_parobject.H
     * and should return it in this method.
     *
     * @return Return unique ParObject id
     */
    int UniqueParObjectId() const override
    {
      return MixtureRuleType::Instance().UniqueParObjectId();
    }


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
     * @param data (in/out) : vector storing all data to be packed into this instance.
     */
    void PackMixtureLaw(DRT::PackBuffer& data) const;

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
    void UnpackMixtureLaw(std::vector<char>::size_type& position, const std::vector<char>& data);

    /*!
     * Initialize the mixture law with the element parameters of the input line
     *
     * @param numgp (in) Number of Gauß-points
     * @param params (in/out) : Parameter list for exchange of parameters
     */
    void ReadElement(int numgp, DRT::INPUT::LineDefinition* linedef);

    /*!
     * This method should be called after ReadElement and stores the ParameterList internally
     *
     * @param params (in/out) : Parameter list containing different parameters
     * @param constituents (in) List of constituents
     */
    void Init(Teuchos::ParameterList& params,
        const Teuchos::RCP<std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>>& constituents);

    /*!
     * Returns whether the constituent is already initialized
     * @return true if the constituent is already initialized, otherwise false
     */
    bool IsInit() { return is_init_; }

    /*!
     * \brief Setup the constituent
     *
     * This method is called once for each element. The constituent should setup all his internal
     * variables and materials.
     * @param gp (in) : Gauß point
     * @param params (in/out) : Container for additional information
     */
    virtual void Setup(int gp, Teuchos::ParameterList& params);

    /*!
     * \brief Update of the material law
     *
     * This simple mixture rule does not need to update anything, so this method is kept empty
     *
     * @param defgrd (in) : Deformation gradient
     * @param params (in/out) : Container for additional information
     * @param gp (in) : Gauß point
     */
    virtual void Update(
        LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params, const int gp)
    {
      // Nothing needs to be updated in this simple mixture law
    }

    /*!
     * Evaluates the constituents. Needs to compute the stress contribution of the constituent out
     * of the displacements. Will be called for each Gauß point
     *
     * @param defgrd (in) : Deformation gradient
     * @param stress (out) :  2nd Piola Kirchhoff stress tensor in stress like Voigt-notation
     * @param cmat (out) : Constitutive tensor in Voigt notation
     */
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
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
     * @param eleGID (in) : element id
     */
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleGID)
    {
      return false;
    }

   protected:
    /*!
     * \brief Returns a reference to the constituents
     * @return
     */
    Teuchos::RCP<std::vector<Teuchos::RCP<MIXTURE::MixtureConstituent>>> Constituents()
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

    ///! Indicator, whether the constituent is already initialized
    bool is_init_;

    ///! Indicator, whether the constituent is already set up
    std::vector<char> is_setup_;
  };
}  // namespace MIXTURE

#endif  // BACI_MIXTURE_RULE_H
