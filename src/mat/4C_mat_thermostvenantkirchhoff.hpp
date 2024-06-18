/*----------------------------------------------------------------------*/
/*! \file
\brief St.Venant Kirchhoff with an additional temperature dependent term
       describing heat expansion

       example input line:
       MAT 1   MAT_Struct_ThrStVenantK YOUNGNUM 2 YOUNG 1.48e8 1.48e5 NUE 0.3 DENS
         9.130e-6 THEXPANS 1.72e-5 INITTEMP 293.15

\level 2

*/
/*----------------------------------------------------------------------*
 | definitions                                               dano 02/10 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_THERMOSTVENANTKIRCHHOFF_HPP
#define FOUR_C_MAT_THERMOSTVENANTKIRCHHOFF_HPP


/*----------------------------------------------------------------------*
 | headers                                                   dano 02/10 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_thermomechanical.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    //! material parameters for de St. Venant--Kirchhoff with temperature
    //! dependent term
    //!
    //! <h3>Input line</h3>
    //! MAT 1 MAT_Struct_ThrStVenantK YOUNG 400 NUE 0.3 DENS 1 THEXPANS 1 INITTEMP 20
    class ThermoStVenantKirchhoff : public Core::Mat::PAR::Parameter
    {
     public:
      //! standard constructor
      explicit ThermoStVenantKirchhoff(const Core::Mat::PAR::Parameter::Data& matdata);

      //! @name material parameters
      //@{

      //! Young's modulus (temperature dependent --> polynomial expression)
      const std::vector<double> youngs_;
      //! Possion's ratio \f$ \nu \f$
      const double poissonratio_;
      //! mass density \f$ \rho \f$
      const double density_;
      //! linear coefficient of thermal expansion  \f$ \alpha_T \f$
      const double thermexpans_;
      //! heat capacity \f$ C_V \f$
      const double capa_;
      //! heat conductivity \f$ k \f$
      const double conduct_;
      //! initial temperature (constant) \f$ \theta_0 \f$
      const double thetainit_;
      //! thermal material id, -1 if not used (old interface)
      const int thermomat_;
      //@}

      //! create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class ThermoStVenantKirchhoff

  }  // namespace PAR

  class ThermoStVenantKirchhoffType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "ThermoStVenantKirchhoffType"; }

    static ThermoStVenantKirchhoffType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ThermoStVenantKirchhoffType instance_;
  };

  /*----------------------------------------------------------------------*/
  //! Wrapper for St.-Venant-Kirchhoff material with temperature term
  class ThermoStVenantKirchhoff : public ThermoMechanicalMaterial
  {
   public:
    //! construct empty material object
    ThermoStVenantKirchhoff();

    //! construct the material object given material parameters
    explicit ThermoStVenantKirchhoff(Mat::PAR::ThermoStVenantKirchhoff* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ThermoStVenantKirchhoffType::Instance().UniqueParObjectId();
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(Inpar::STR::KinemType kinem) override
    {
      if (kinem != Inpar::STR::KinemType::linear && kinem != Inpar::STR::KinemType::nonlinearTotLag)
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
    */
    void Pack(
        Core::Communication::PackBuffer& data  //!< (i/o): char vector to store class information
    ) const override;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().
    */
    void Unpack(const std::vector<char>&
            data  //!< (i) : vector storing all data to be unpacked into this instance.
        ) override;

    //@}

    //! material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_thermostvenant;
    }

    //! return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new ThermoStVenantKirchhoff(*this));
    }

    //! evaluates stresses for 3d
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,  //!< deformation gradient
        const Core::LinAlg::Matrix<6, 1>* glstrain,          //!< Green-Lagrange strain
        Teuchos::ParameterList& params,                      //!< parameter list
        Core::LinAlg::Matrix<6, 1>* stress,                  //!< stress
        Core::LinAlg::Matrix<6, 6>* cmat,                    //!< elastic material tangent
        int gp,                                              ///< Gauss point
        int eleGID                                           //!< element GID
        ) override;

    /// add strain energy
    void StrainEnergy(const Core::LinAlg::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
        double& psi,                                               ///< strain energy functions
        int gp,                                                    ///< Gauss point
        int eleGID                                                 ///< element GID
        ) override;

    //! return true if Young's modulus is temperature dependent
    bool youngs_is_temp_dependent() const { return this->params_->youngs_.size() > 1; }

    //! density \f$ \rho \f$
    double Density() const override { return params_->density_; }

    //! conductivity \f$ k \f$
    double Conductivity() const { return params_->conduct_; }

    //! material capacity \f$ C_V \f$
    double Capacity() const override { return params_->capa_; }

    //! Initial temperature \f$ \theta_0 \f$
    double InitTemp() const { return params_->thetainit_; }

    //! Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    void evaluate(const Core::LinAlg::Matrix<3, 1>& gradtemp, Core::LinAlg::Matrix<3, 3>& cmat,
        Core::LinAlg::Matrix<3, 1>& heatflux) const override;

    void evaluate(const Core::LinAlg::Matrix<2, 1>& gradtemp, Core::LinAlg::Matrix<2, 2>& cmat,
        Core::LinAlg::Matrix<2, 1>& heatflux) const override;

    void evaluate(const Core::LinAlg::Matrix<1, 1>& gradtemp, Core::LinAlg::Matrix<1, 1>& cmat,
        Core::LinAlg::Matrix<1, 1>& heatflux) const override;

    void ConductivityDerivT(Core::LinAlg::Matrix<3, 3>& dCondDT) const override;

    void ConductivityDerivT(Core::LinAlg::Matrix<2, 2>& dCondDT) const override;

    void ConductivityDerivT(Core::LinAlg::Matrix<1, 1>& dCondDT) const override;

    double CapacityDerivT() const override;

    void Reinit(double temperature, unsigned gp) override;

    void ResetCurrentState() override;

    void CommitCurrentState() override;

    void Reinit(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, double temperature, unsigned gp) override;

    void GetdSdT(Core::LinAlg::Matrix<6, 1>* dS_dT) override;

    void stress_temperature_modulus_and_deriv(
        Core::LinAlg::Matrix<6, 1>& stm, Core::LinAlg::Matrix<6, 1>& stm_dT) override;

    //! general thermal tangent of material law depending on stress-temperature modulus
    static void FillCthermo(Core::LinAlg::Matrix<6, 1>& ctemp, double m);

   private:
    //! computes isotropic elasticity tensor in matrix notion for 3d
    void setup_cmat(Core::LinAlg::Matrix<6, 6>& cmat);

    //! computes temperature dependent isotropic elasticity tensor in matrix
    //! notion for 3d
    void setup_cthermo(Core::LinAlg::Matrix<6, 1>& ctemp);

    //! calculates stress-temperature modulus
    double st_modulus() const;

    //! calculates stress-temperature modulus
    double get_st_modulus_t() const;

    //! calculates derivative of Cmat with respect to current temperatures
    //! only in case of temperature-dependent material parameters
    void get_cmat_at_tempnp_t(Core::LinAlg::Matrix<6, 6>& derivcmat);

    //! calculates derivative of Cmat with respect to current temperatures
    //! only in case of temperature-dependent material parameters
    void get_cthermo_at_tempnp_t(
        Core::LinAlg::Matrix<6, 1>& derivctemp  //!< linearisation of ctemp w.r.t. T
    );

    //! calculate temperature dependent material parameter and return value
    double get_mat_parameter_at_tempnp(
        const std::vector<double>* paramvector,  //!< (i) given parameter is a vector
        const double& tempnp                     // tmpr (i) current temperature
    ) const;

    //! calculate temperature dependent material parameter and return value
    double get_mat_parameter_at_tempnp_t(
        const std::vector<double>* paramvector,  //!< (i) given parameter is a vector
        const double& tempnp                     // tmpr (i) current temperature
    ) const;

    //! create thermo material object if specified in input (!= -1)
    void create_thermo_material_if_set();

    //! my material parameters
    Mat::PAR::ThermoStVenantKirchhoff* params_;

    //! pointer to the internal thermal material
    Teuchos::RCP<Mat::Trait::Thermo> thermo_;

    //! current temperature (set by Reinit())
    double current_temperature_{};

    //! current Green-Lagrange strain
    const Core::LinAlg::Matrix<6, 1>* current_glstrain_{};

  };  // ThermoStVenantKirchhoff
}  // namespace Mat


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
