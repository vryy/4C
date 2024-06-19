/*----------------------------------------------------------------------*/
/*! \file
 \brief wrapper for structure material of porous media including reactive reference porosity


\level 3
 *----------------------------------------------------------------------*/


#ifndef FOUR_C_MAT_STRUCTPORO_REACTION_HPP
#define FOUR_C_MAT_STRUCTPORO_REACTION_HPP

#include "4C_config.hpp"

#include "4C_mat_structporo.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  // forward declaration
  class StructPoroReaction;

  namespace PAR
  {
    class StructPoroReaction : public PAR::StructPoro
    {
      friend class Mat::StructPoroReaction;

     public:
      /// standard constructor
      StructPoroReaction(const Core::Mat::PAR::Parameter::Data& matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      /// @name material parameters
      //@{

      int dofIDReacScalar_;
      //@}

    };  // class StructPoroReaction

  }  // namespace PAR

  class StructPoroReactionType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "StructPoroReactionType"; }

    static StructPoroReactionType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static StructPoroReactionType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for StructPoro material
  ///
  /// This object exists (several times) at every element
  class StructPoroReaction : public StructPoro
  {
   public:
    /// construct empty material object
    StructPoroReaction();

    /// construct the material object given material parameters
    explicit StructPoroReaction(Mat::PAR::StructPoroReaction* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of parobject.H (this file) and should return it in this method.
     */
    int UniqueParObjectId() const override
    {
      return StructPoroReactionType::Instance().UniqueParObjectId();
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
      return Core::Materials::m_structpororeaction;
    }

    /// return initial porosity
    double ref_porosity_time_deriv() const override { return dphiDphiref_ * refporositydot_; }

    /// compute current porosity and save it
    void compute_porosity(Teuchos::ParameterList& params,  ///< (i) element parameter list
        double press,                                      ///< (i) pressure at gauss point
        double J,          ///< (i) determinant of jacobian at gauss point
        int gp,            ///< (i) number of current gauss point
        double& porosity,  ///< (o) porosity at gauss point
        double* dphi_dp,   ///< (o) first derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dJ,   ///< (o) first derivative of porosity w.r.t. jacobian at gauss point
        double*
            dphi_dJdp,  ///< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
        double* dphi_dJJ,  ///< (o) second derivative of porosity w.r.t. jacobian at gauss point
        double* dphi_dpp,  ///< (o) second derivative of porosity w.r.t. pressure at gauss point
        bool save = true) override;

    //! evaluate constitutive relation for porosity and compute derivatives
    void constitutive_derivatives(Teuchos::ParameterList& params,  ///< (i) parameter list
        double press,        ///< (i) fluid pressure at gauss point
        double J,            ///< (i) Jacobian determinant at gauss point
        double porosity,     ///< (i) porosity at gauss point
        double* dW_dp,       ///< (o) derivative of potential w.r.t. pressure
        double* dW_dphi,     ///< (o) derivative of potential w.r.t. porosity
        double* dW_dJ,       ///< (o) derivative of potential w.r.t. jacobian
        double* dW_dphiref,  ///< (o) derivative of potential w.r.t. reference porosity
        double* W            ///< (o) inner potential
        ) override;

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new StructPoroReaction(*this));
    }

    /// Initialize internal variables
    void setup(int numgp,  ///< number of Gauss points
        Input::LineDefinition* linedef) override;

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    /// return reference porosity average
    double RefPorosityAv() const;

    //! @name Evaluation methods

    /// evaluate material law
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,  ///< (i) deformation gradient
        const Core::LinAlg::Matrix<6, 1>* glstrain,          ///< (i) green lagrange strain
        Teuchos::ParameterList& params,                      ///< (i) parameter list
        Core::LinAlg::Matrix<6, 1>* stress,                  ///< (o) second piola kirchhoff stress
        Core::LinAlg::Matrix<6, 6>* cmat,                    ///< (o) constitutive matrix
        int gp,                                              ///< (i) Gauss point
        int eleGID) override;

    //@}

    //! @name Visualization methods

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

    //@}

   protected:
    virtual void reaction(const double porosity, const double J,
        Teuchos::RCP<std::vector<double>> scalars, Teuchos::ParameterList& params);

    /// my material parameters
    Mat::PAR::StructPoroReaction* params_;

    /// reference porosity
    double refporosity_;

    /// derivative of porosity w.r.t. reference porosity
    double dphiDphiref_;

    /// time derivative of reference porosity
    double refporositydot_;
  };

}  // namespace Mat


FOUR_C_NAMESPACE_CLOSE

#endif
