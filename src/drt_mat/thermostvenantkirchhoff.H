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
#ifndef THERMOSTVENANTKIRCHHOFF_H_
#define THERMOSTVENANTKIRCHHOFF_H_


/*----------------------------------------------------------------------*
 | headers                                                   dano 02/10 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "material_thermomechanical.H"


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    //! material parameters for de St. Venant--Kirchhoff with temperature
    //! dependent term
    //!
    //! <h3>Input line</h3>
    //! MAT 1 MAT_Struct_ThrStVenantK YOUNG 400 NUE 0.3 DENS 1 THEXPANS 1 INITTEMP 20
    class ThermoStVenantKirchhoff : public Parameter
    {
     public:
      //! standard constructor
      explicit ThermoStVenantKirchhoff(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      ~ThermoStVenantKirchhoff() override { ; }

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
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class ThermoStVenantKirchhoff

  }  // namespace PAR

  class ThermoStVenantKirchhoffType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "ThermoStVenantKirchhoffType"; }

    static ThermoStVenantKirchhoffType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

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
    explicit ThermoStVenantKirchhoff(MAT::PAR::ThermoStVenantKirchhoff* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ThermoStVenantKirchhoffType::Instance().UniqueParObjectId();
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_linear || kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
    */
    void Pack(DRT::PackBuffer& data  //!< (i/o): char vector to store class information
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
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_thermostvenant; }

    //! return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new ThermoStVenantKirchhoff(*this));
    }

    //! evaluates stresses for 3d
    void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  //!< deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,          //!< Green-Lagrange strain
        Teuchos::ParameterList& params,                //!< parameter list
        LINALG::Matrix<6, 1>* stress,                  //!< stress
        LINALG::Matrix<6, 6>* cmat,                    //!< elastic material tangent
        int gp,                                        ///< Gauss point
        int eleGID                                     //!< element GID
        ) override;

    /// add strain energy
    void StrainEnergy(const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
        double& psi,                                         ///< strain energy functions
        int gp,                                              ///< Gauss point
        int eleGID                                           ///< element GID
        ) override;

    //! return true if Young's modulus is temperature dependent
    bool YoungsIsTempDependent() const { return this->params_->youngs_.size() > 1; }

    //! density \f$ \rho \f$
    double Density() const override { return params_->density_; }

    //! conductivity \f$ k \f$
    double Conductivity() const { return params_->conduct_; }

    //! material capacity \f$ C_V \f$
    double Capacity() const override { return params_->capa_; }

    //! Initial temperature \f$ \theta_0 \f$
    double InitTemp() const { return params_->thetainit_; }

    //! Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    void Evaluate(const LINALG::Matrix<3, 1>& gradtemp, LINALG::Matrix<3, 3>& cmat,
        LINALG::Matrix<3, 1>& heatflux) const override;

    void Evaluate(const LINALG::Matrix<2, 1>& gradtemp, LINALG::Matrix<2, 2>& cmat,
        LINALG::Matrix<2, 1>& heatflux) const override;

    void Evaluate(const LINALG::Matrix<1, 1>& gradtemp, LINALG::Matrix<1, 1>& cmat,
        LINALG::Matrix<1, 1>& heatflux) const override;

    void ConductivityDerivT(LINALG::Matrix<3, 3>& dCondDT) const override;

    void ConductivityDerivT(LINALG::Matrix<2, 2>& dCondDT) const override;

    void ConductivityDerivT(LINALG::Matrix<1, 1>& dCondDT) const override;

    double CapacityDerivT() const override;

    void Reinit(double temperature, unsigned gp) override;

    void ResetCurrentState() override;

    void CommitCurrentState() override;

    void Reinit(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        double temperature, unsigned gp) override;

    void GetdSdT(LINALG::Matrix<6, 1>* dS_dT) override;

    void StressTemperatureModulusAndDeriv(
        LINALG::Matrix<6, 1>& stm, LINALG::Matrix<6, 1>& stm_dT) override;

    //! general thermal tangent of material law depending on stress-temperature modulus
    static void FillCthermo(LINALG::Matrix<6, 1>& ctemp, double m);

   private:
    //! computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(LINALG::Matrix<6, 6>& cmat);

    //! computes temperature dependent isotropic elasticity tensor in matrix
    //! notion for 3d
    void SetupCthermo(LINALG::Matrix<6, 1>& ctemp);

    //! calculates stress-temperature modulus
    double STModulus() const;

    //! calculates stress-temperature modulus
    double GetSTModulus_T() const;

    //! calculates derivative of Cmat with respect to current temperatures
    //! only in case of temperature-dependent material parameters
    void GetCmatAtTempnp_T(LINALG::Matrix<6, 6>& derivcmat);

    //! calculates derivative of Cmat with respect to current temperatures
    //! only in case of temperature-dependent material parameters
    void GetCthermoAtTempnp_T(LINALG::Matrix<6, 1>& derivctemp  //!< linearisation of ctemp w.r.t. T
    );

    //! calculate temperature dependent material parameter and return value
    double GetMatParameterAtTempnp(
        const std::vector<double>* paramvector,  //!< (i) given parameter is a vector
        const double& tempnp                     // tmpr (i) current temperature
        ) const;

    //! calculate temperature dependent material parameter and return value
    double GetMatParameterAtTempnp_T(
        const std::vector<double>* paramvector,  //!< (i) given parameter is a vector
        const double& tempnp                     // tmpr (i) current temperature
        ) const;

    //! create thermo material object if specified in input (!= -1)
    void CreateThermoMaterialIfSet();

    //! my material parameters
    MAT::PAR::ThermoStVenantKirchhoff* params_;

    //! pointer to the internal thermal material
    Teuchos::RCP<MAT::TRAIT::Thermo> thermo_;

    //! current temperature (set by Reinit())
    double currentTemperature_{};

    //! current Green-Lagrange strain
    const LINALG::Matrix<6, 1>* currentGlstrain_{};

  };  // ThermoStVenantKirchhoff
}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // THERMOSTVENANTKIRCHHOFF_H
