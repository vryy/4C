/*----------------------------------------------------------------------*/
/*! \file
 \brief porous material for dissolution reaction in ECM model


\level 3
 *----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_STRUCTPORO_REACTION_ECM_HPP
#define FOUR_C_MAT_STRUCTPORO_REACTION_ECM_HPP

#include "4C_config.hpp"

#include "4C_mat_structporo_reaction.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  // forward declaration
  class StructPoroReactionECM;

  namespace PAR
  {
    class StructPoroReactionECM : public PAR::StructPoroReaction
    {
      friend class Mat::StructPoroReactionECM;

     public:
      /// standard constructor
      StructPoroReactionECM(const Core::Mat::PAR::Parameter::Data& matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      /// @name material parameters
      //@{

      /// density of collagen
      double densCollagen_;
      //@}
    };
    // class StructPoroReactionECM

  }  // namespace PAR

  class StructPoroReactionECMType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "StructPoroReactionECMType"; }

    static StructPoroReactionECMType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static StructPoroReactionECMType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for StructPoroReactionECM material
  ///
  /// This object exists (several times) at every element
  class StructPoroReactionECM : public StructPoroReaction
  {
   public:
    /// construct empty material object
    StructPoroReactionECM();

    /// construct the material object given material parameters
    explicit StructPoroReactionECM(Mat::PAR::StructPoroReactionECM* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of parobject.H (this file) and should return it in this method.
     */
    int UniqueParObjectId() const override
    {
      return StructPoroReactionECMType::Instance().UniqueParObjectId();
    }

    /*!
     \brief Pack this class so it can be communicated

     Resizes the vector data and stores all information of a class in it.
     The first information to be stored in data has to be the
     unique parobject id delivered by UniqueParObjectId() which will then
     identify the exact class on the receiving processor.

     \param data (in/out): char vector to store class information
     */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
     \brief Unpack data from a char vector into this class

     The vector data contains all information to rebuild the
     exact copy of an instance of a class on a different processor.
     The first entry in data has to be an integer which is the unique
     parobject id defined at the top of this file and delivered by
     UniqueParObjectId().

     \param data (in) : vector storing all data to be unpacked into this
     instance.
     */
    void unpack(const std::vector<char>& data) override;

    //@}

    /// material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_structpororeactionECM;
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new StructPoroReactionECM(*this));
    }

    /// Initialize internal variables
    void setup(int numgp,  ///< number of Gauss points
        Input::LineDefinition* linedef) override;

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    /// evaluate chemical potential
    virtual void ChemPotential(
        const Core::LinAlg::Matrix<6, 1>& glstrain,  ///< (i) green lagrange strain
        const double porosity,                       ///< (i) porosity
        const double press,                          ///< (i) pressure at gauss point
        const double J,                              ///< (i) determinant of jacobian at gauss point
        int EleID,                                   ///< (i) element GID
        double& pot,                                 ///< (o) chemical potential
        const int gp);

    /// Update of GP data (e.g., history variables)
    void Update() override;

    //! @name Visualization methods

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

   protected:
    void reaction(const double porosity, const double J, Teuchos::RCP<std::vector<double>> scalars,
        Teuchos::ParameterList& params) override;

    /// reference porosity at time step n
    double refporosity_old_;

    /// time derivative of reference porosity at time step n
    double refporositydot_old_;

    /// chemical potential
    std::vector<double> chempot_;

    /// initial chemical potential
    std::vector<double> chempot_init_;

    /// my material parameters
    Mat::PAR::StructPoroReactionECM* params_;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
