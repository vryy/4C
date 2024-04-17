/*----------------------------------------------------------------------*/
/*! \file
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following an associative von Mises plasticity and a linear thermoelastic
       material law (Thermo St.Venant Kirchhoff).

       Combination of the two material model
       - PlastiLinElastic
       - ThermoStVenantKirchhoff

       Thermoelastoplastic material with mixed hardening
        - so far:
          - linear isotropic hardening with constant hardening modulus Hiso
            - yield stress sigma_y = sigma_y0 + Hiso . strainplbar
          - linear kinematic hardening with constant hardening modulus Hkin
            - 'linear Armstrong-Frederick kinematic hardening'
            - describing the Bauschinger effect via \f Hkin \,=\, const.\f

       small strains including temperature dependency

       strain-energy potential
       \f \rho \psi \,=\,\rho (\psi^e \,+\, \psi^p)
          \,=\, mu strain^e : strain^e + 1/2 . lambda (strain^e : I)^2 +
                + m (T - T_0) strain^e : I - rho C_V (T ln T/T_o - (T - T_0))
                + 1/2 strain^p : strain^p . (2/3 Hkin)
                + 1/2 strainbar^p : strainbar^p . Hiso \f

       example input line:
       MAT 1 MAT_Struct_ThrPlasticLinElast YOUNG 206.9 NUE 0.29 DENS 0.0
         THEXPANS 1.72e-5 INITTEMP 293.15 YIELD 0.45 ISOHARD 0.0 KINHARD 2.0 TOL 1.0e-6

\level 2

*/
/*----------------------------------------------------------------------*
 | definitions                                               dano 08/11 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_THERMOPLASTICLINELAST_HPP
#define FOUR_C_MAT_THERMOPLASTICLINELAST_HPP


/*----------------------------------------------------------------------*
 | headers                                                   dano 08/11 |
 *----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_mat_par_parameter.hpp"
#include "baci_mat_so3_material.hpp"

FOUR_C_NAMESPACE_OPEN


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    //! material parameters for plastic lin elastic Material
    class ThermoPlasticLinElast : public Parameter
    {
     public:
      //! standard constructor
      ThermoPlasticLinElast(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! @name material parameters
      //@{

      //! Young's modulus
      const double youngs_;
      //! Possion's ratio
      const double poissonratio_;
      //! mass density
      const double density_;
      //! linear coefficient of thermal expansion \f$ \alpha_T \f$
      const double thermexpans_;
      //! initial temperature (constant) \f$ \theta_0 \f$
      const double thetainit_;
      //! yield stress (constant for perfect plasticity)  \f$ \sigma_y \f$
      const double yield_;
      //! isotropic hardening modulus \f$ H_i \f$
      const double isohard_;
      //! kinematic hardening modulus \f$ H_k \f$
      const double kinhard_;
      //! isotropic hardening curve
      const std::vector<double> sigma_y_;
      //! accumulated plastic strain
      const std::vector<double> strainbar_p_ref_;
      //! tolerance for local Newton iteration
      const double abstol_;

      //@}

      //! create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

    };  // class ThermoThermoPlasticLinElast

  }  // namespace PAR


  class ThermoPlasticLinElastType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "ThermoPlasticLinElastType"; }

    static ThermoPlasticLinElastType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ThermoPlasticLinElastType instance_;

  };  // ThermoThermoPlasticLinElastType


  /*----------------------------------------------------------------------*/
  //! wrapper for linear thermo-elasto-plastic material
  class ThermoPlasticLinElast : public So3Material
  {
   public:
    //! construct empty material object
    ThermoPlasticLinElast();

    //! construct the material object given material parameters
    explicit ThermoPlasticLinElast(MAT::PAR::ThermoPlasticLinElast* params);

    //! @name Packing and Unpacking

    //!  \brief return unique ParObject id
    //!
    //!  every class implementing ParObject needs a unique id defined at the
    //!  top of parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return ThermoPlasticLinElastType::Instance().UniqueParObjectId();
    }

    //!  \brief Pack this class so it can be communicated
    //!
    //!  Resizes the vector data and stores all information of a class in it.
    //!  The first information to be stored in data has to be the
    //!  unique parobject id delivered by UniqueParObjectId() which will then
    //!  identify the exact class on the receiving processor.
    void Pack(CORE::COMM::PackBuffer& data  //!<  data (i/o): char vector to store class information
    ) const override;

    //!  \brief Unpack data from a char vector into this class
    //!
    //!  The vector data contains all information to rebuild the
    //!  exact copy of an instance of a class on a different processor.
    //!  The first entry in data has to be an integer which is the unique
    //!  parobject id defined at the top of this file and delivered by
    //!  UniqueParObjectId().
    void Unpack(const std::vector<char>&
            data  //!< (i) : vector storing all data to be unpacked into this instance.
        ) override;

    //@}

    //! material type
    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_thermopllinelast;
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::KinemType::linear))
        dserror("element and material kinematics are not compatible");
    }

    //! return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new ThermoPlasticLinElast(*this));
    }

    //! initialise internal stress variables
    void Setup(int numgp, INPUT::LineDefinition* linedef) override;

    //! update internal stress variables
    void Update() override;

    //! reset internal stress variables
    void Reset();

    //! evaluate material
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>* glstrain,
        Teuchos::ParameterList& params,                  //!< parameter list for communication
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>* stress,  //!< 2nd PK-stress
        CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  //!< material stiffness matrix
        int gp,                                                    ///< Gauss point
        int eleGID) override;

    // computes stress
    void Stress(const double p,                                   //!< volumetric stress tensor
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& devstress,  //!< deviatoric stress tensor
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& stress            //!< 2nd PK-stress
    );

    //! calculate relative/over strees
    void RelDevStress(
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& devstress,  //!< deviatoric stress tensor
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& beta,       //!< back stress tensor
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& eta               //!< relative stress
    );

    //! computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(
        CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat  //!< elastic material tangent
    );

    //! computes isotropic elastoplastic tensor in matrix notion for 3d
    void SetupCmatElastoPlastic(CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
                                    cmat,  //!< elasto-plastic material tangent
        double Dgamma,                     //!< plastic multiplier
        double G,                          //!< shear modulus
        double q,                          //!< relative effective stress
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1> flowvector,  //!< flow vector
        const CORE::LINALG::Matrix<NUM_STRESS_3D, 1> Nbar,        //!< unit vector
        double heaviside,                                         //!< Heaviside function
        double Hiso,                                              //!< isotropic hardening modulus
        double Hkin                                               //!< kinematic hardening modulus
    );

    //! computes isotropic elastoplastic tensor in matrix notation for 3d
    void SetupCmatElastoPlastic2(CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        double Dgamma, double q, CORE::LINALG::Matrix<NUM_STRESS_3D, 1> unitflow);

    //! computes continuum elastoplastic tensor in matrix notation for 3d
    void SetupContinuumCmatElastoPlastic(CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,
        double Dgamma, double q, CORE::LINALG::Matrix<NUM_STRESS_3D, 1> unitflow);

    //! calculates the derivative of GetSigmaYAtStrainbarnp() w.r.t. astrain_{n+1}
    //! and returns the isotropic hardening modulus
    double GetIsoHardAtStrainbarnp(
        const double strainbar_p  //!< current accumulated strain
                                  //!< or if damage!=0: isotropic hardening internal variable
    );

    //! calculates current yield stress from (sigma_y-astrain)-samples which
    //! describe a piecewise constant function
    double GetSigmaYAtStrainbarnp(
        const double strainbar_p  //!< current accumulated strain, in case of dependent hardening
                                  //!< or if damage!=0: isotropic hardening internal variable
    );

    //! return density
    double Density() const override { return params_->density_; }

    //! return scalar-valued accumulated strain at Gauss points
    //! use method to return values at END of time step
    double AccumulatedStrain(int gp  //!< current Gauss point
    ) const
    {
      //! use the old vector (last_) for postprocessing
      //! Output is called after(!!) Update, so that newest values are included in
      //! the old history vectors last_
      //! in contrast: the current history vectors curr_are reset to zero
      return (strainbarpllast_->at(gp));
    }

    // return current scalar-valued accumulated strain
    double AccumulatedStrainCurr(int gp  //!< current Gauss point
    ) const
    {
      return (strainbarplcurr_->at(gp));
    }

    //! check if history variables are already initialised
    bool Initialized() const { return (isinit_ and (strainplcurr_ != Teuchos::null)); }

    //! return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    //! @name temperature specific methods
    //@{

    //! main 3D material call to determine stress and constitutive tensor ctemp
    //  originally method of fourieriso with const!!!
    void Evaluate(const CORE::LINALG::Matrix<1, 1>& Ntemp,  //!< temperature of element
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& ctemp,  //!< temperature dependent material tangent
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& stresstemp  //!< temperature dependent stress term
    );

    //! computes temperature dependent isotropic elasticity tensor in matrix
    //! notion for 3d
    void SetupCthermo(
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& ctemp  //!< temperature dependent material tangent
    );

    //! calculates stress-temperature modulus
    double STModulus();

    //! initial temperature
    double InitTemp() const { return params_->thetainit_; }

    //@}

    //! @name specific methods for TSI and plastic material
    //@{

    //! calculate elastic strain rate, using additive split of strains
    //! (o) save strain^e' in strainelrate_
    void StrainRateSplit(int gp,                            //!< (i): current Gauss point
        const double stepsize,                              //!< (i): stepsize
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& strainrate  //!< (i): total strain rate ( B d')
    );

    //! return current plastic strain vector
    //! \f${\boldsymbol \varepsilon}^p_{n+1}\f$
    CORE::LINALG::Matrix<NUM_STRESS_3D, 1> PlasticStrain(int gp) const
    {
      return strainplcurr_->at(gp);
    }

    //! return current elastic strain rate vector
    //! \f$\dot{\boldsymbol \varepsilon}^e_{n+1}\f$
    CORE::LINALG::Matrix<NUM_STRESS_3D, 1> ElasticStrainRate(int gp) const
    {
      return strainelrate_->at(gp);
    }

    //! compute internal dissipation terms
    void Dissipation(int gp,                           // current Gauss point
        double sigma_yiso,                             // isotropic work hardening von Mises stress
        double Dgamma,                                 // plastic multiplier/increment
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1> N,      // flow vector
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1> stress  // total mechanical stress
    );

    //! compute linearisation of internal dissipation for k_Td
    void DissipationCouplCond(CORE::LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
                                  cmat,                 // elasto-plastic tangent modulus (out)
        int gp,                                         // current Gauss point
        double G,                                       // shear modulus
        double Hiso,                                    // isotropic hardening modulus
        double Hkin,                                    // kinematic hardening modulus
        double heaviside,                               // Heaviside function
        double etanorm,                                 // norm of eta^{trial}_{n+1}
        double Dgamma,                                  // plastic multiplier
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& N,      // flow vector
        CORE::LINALG::Matrix<NUM_STRESS_3D, 1>& stress  // flow vector
    );

    //! return mechanical dissipation term due to kinematic hardening
    //! \f$ D_{\rm mech,kin}\f$
    double MechanicalKinematicDissipation(int gp) const { return Dmech_->at(gp); }

    //! return linearisation of mechanical dissipation w.r.t. displacements
    //! \f$ k_Td += \dfrac{\partial D_{\rm mech}}{\partial d_{n+1}}\f$
    CORE::LINALG::Matrix<NUM_STRESS_3D, 1> DissipationLinearisedForCouplCond(int gp) const
    {
      return Dmech_d_->at(gp);
    }

    //@}

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;


   private:
    //! my material parameters
    MAT::PAR::ThermoPlasticLinElast* params_;

    //! plastic history vector
    //! old plastic strain at t_n
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        strainpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    //! current plastic strain at t_n+1
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        strainplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    //! old back stress at t_n
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        backstresslast_;  //!< \f${\beta}_{n}\f$
    //! current back stress at t_n+1
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        backstresscurr_;  //!< \f${\beta}_{n+1}\f$
    //! old accumulated plastic strain at t_n
    Teuchos::RCP<std::vector<double>> strainbarpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    //! current accumulated plastic strain at t_n+1
    Teuchos::RCP<std::vector<double>> strainbarplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$

    //! @name specific methods for the combination TSI and plastic material
    //@{

    //! mechanical dissipation term due to kinematic and isotropic hardening at t_n+1
    Teuchos::RCP<std::vector<double>> Dmech_;
    //! linearisation of mechanical dissipation term w.r.t. to d_{n+1}
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>> Dmech_d_;
    //! save plastic strain increment: Incstrainpl = strain^p_{n+1} - strain^p_n
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>> Incstrainpl_;
    //! elastic strain rate required for thermoelastic heating term
    //! use additive split: strain^e' = strain' - strain^p'
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>> strainelrate_;

    //@}

    //! indicator if #Initialize routine has been called
    bool isinit_;
    //! indicator if material has started to be plastic
    bool plastic_step_;
    //! element ID, in which first plasticity occurs
    int plastic_eleID_;

  };  // class ThermoThermoPlasticLinElast : public Material
}  // namespace MAT


/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
