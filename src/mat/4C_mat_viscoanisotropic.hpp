/*----------------------------------------------------------------------*/
/*! \file
\brief
Anisotropic viscohyperelastic material
The input line should read
MAT 1 MAT_VISCOANISO KAPPA 1.6667E04 MUE 33.3556 DENS 0.0000000001 K1 100.0 K2 1.0 GAMMA 20.0
BETA_ISO 1.E4 BETA_ANISO 1.E4 RELAX_ISO 0.0010001 RELAX_ANISO 0

\level 2

*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_VISCOANISOTROPIC_HPP
#define FOUR_C_MAT_VISCOANISOTROPIC_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters
    class ViscoAnisotropic : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      ViscoAnisotropic(const Core::Mat::PAR::Parameter::Data& matdata);

      /// @name material parameters
      //@{
      const double kappa_;
      const double mue_;
      const double density_;
      const double k1_;
      const double k2_;
      const double gamma_;
      const int numstresstypes_;
      double beta_[2];
      double relax_[2];
      const double minstretch_;
      const int elethick_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class ViscoAnisotropic

  }  // namespace PAR

  class ViscoAnisotropicType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "ViscoAnisotropicType"; }

    static ViscoAnisotropicType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ViscoAnisotropicType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for Visco-NeoHooke material
  class ViscoAnisotropic : public So3Material
  {
   public:
    /// construct empty material object
    ViscoAnisotropic();

    /// construct the material object given material parameters
    explicit ViscoAnisotropic(Mat::PAR::ViscoAnisotropic* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ViscoAnisotropicType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
      This material contains history variables, which are packed for restart purposes.

      \param data (in/out): char vector to store class information
    */
    void Pack(Core::Communication::PackBuffer& data) const override;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().
      History data is unpacked in restart.

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_viscoanisotropic;
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(Inpar::STR::KinemType kinem) override
    {
      if (!(kinem == Inpar::STR::KinemType::nonlinearTotLag))
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new ViscoAnisotropic(*this));
    }

    /// Setup and Initialize internal stress variables
    void Setup(int numgp,  ///< number of Gauss points
        Input::LineDefinition* linedef) override;

    /// Setup and Initialize internal stress variables and align fibers based on a given vector
    void Setup(const int numgp,             ///< number of Gauss points
        const std::vector<double> thickvec  ///< direction fibers should be oriented in
    );

    /// Update internal stress variables
    void Update() override;

    void UpdateFiberDirs(const int numgp, Core::LinAlg::Matrix<3, 3>* defgrad);

    /// Evaluate material
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,      ///< deformation gradient
        const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>* glstrain,  ///< green lagrange strain
        Teuchos::ParameterList& params,                  ///< parameter list for communication
        Core::LinAlg::Matrix<NUM_STRESS_3D, 1>* stress,  ///< 2nd PK-stress
        Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
        int gp,                                                    ///< Gauss point
        int eleGID                                                 ///< element GID
        ) override;

    /// Return density
    double Density() const override { return params_->density_; };

    /// Return shear modulus
    double shear_mod() const { return params_->mue_; };


    /// Check if history variables are already initialized
    bool Initialized() const { return isinit_ && (histstresscurr_ != Teuchos::null); }

    /// return a1s
    Teuchos::RCP<std::vector<std::vector<double>>> Geta1() const { return ca1_; }

    /// return a2s
    Teuchos::RCP<std::vector<std::vector<double>>> Geta2() const { return ca2_; }


    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

   private:
    /// my material parameters
    Mat::PAR::ViscoAnisotropic* params_;

    // internal variables for fibers
    Teuchos::RCP<std::vector<std::vector<double>>> a1_;  ///< first fiber vector per gp (reference)
    Teuchos::RCP<std::vector<std::vector<double>>> a2_;  ///< second fiber vector per gp (reference)
    Teuchos::RCP<std::vector<std::vector<double>>>
        ca1_;  ///< first fiber vector per gp (spatial config)
    Teuchos::RCP<std::vector<std::vector<double>>>
        ca2_;  ///< second fiber vector per gp (spatial config)

    // visco history stresses for every gausspoint and every stress type
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>>>
        histstresscurr_;  ///< current stress
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>>>
        histstresslast_;  ///< stress of last converged state
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>>>
        artstresscurr_;  ///< current artificial stress
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>>>
        artstresslast_;  ///< artificial stress in last converged state

    bool isinit_;  ///< indicates if material is initialized
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
