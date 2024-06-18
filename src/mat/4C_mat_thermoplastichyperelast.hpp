/*----------------------------------------------------------------------*/
/*! \file
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following finite strain
       von-Mises plasticity with nonlinear isotropic hardening and general
       hyperelasticity (for the time being: NeoHooke).

       implementation is based on
       Simo and Miehe: "Associative coupled thermoplasticity at finite strains:
       Formulation, numerical analysis and implementation", in Computer Methods
       in Applied Mechanics and Engineering, 98:41-104, 1992.

       geometrically nonlinear, for finite strains, rate-independent, thermo-plastic

       example input line:
       [mm,ms,kg,K,GPa]
       MAT 1 MAT_Struct_ThrPlasticHyperElast YOUNG 206.9 NUE 0.29 DENS 7.8e-6
         CTE 1e-5 INITTEMP 293 YIELD 0.45 ISOHARD 0.12924 SATHARDENING 0.715
         HARDEXPO 16.93 YIELDSOFT 0.002 HARDSOFT 0.002 DISSFACT 0.9 TOL 1.0e-06

         Seitz 11/16: There are still some linearizations for the plastic
         heating terms missing, use with caution.

\level 3

*/
/*----------------------------------------------------------------------*
 | definitions                                               dano 03/13 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_THERMOPLASTICHYPERELAST_HPP
#define FOUR_C_MAT_THERMOPLASTICHYPERELAST_HPP

/*----------------------------------------------------------------------*
 | headers                                                   dano 03/13 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    //! material parameters for neo-Hooke
    class ThermoPlasticHyperElast : public Core::Mat::PAR::Parameter
    {
     public:
      //! standard constructor
      ThermoPlasticHyperElast(const Core::Mat::PAR::Parameter::Data& matdata);

      //! @name material parameters
      //@{

      //! Young's modulus
      const double youngs_;
      //! Possion's ratio
      const double poissonratio_;
      //! mass density
      const double density_;
      //! coefficient of thermal expansion
      const double cte_;
      //! initial, reference temperature
      const double inittemp_;
      //! initial yield stress (constant) at reference temperature
      const double yield_;
      //! linear isotropic hardening modulus at reference temperature
      const double isohard_;
      // saturation hardening at reference temperature
      const double sathardening_;
      //! hardening exponent
      const double hardexpo_;
      //! yield stress softening
      const double yieldsoft_;
      //! hardening softening
      const double hardsoft_;
      //! tolerance for local Newton iteration
      const double abstol_;

      //@}

      //! create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

    };  // class ThermoPlasticHyperElast

  }  // namespace PAR


  class ThermoPlasticHyperElastType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "ThermoPlasticHyperElastType"; }

    static ThermoPlasticHyperElastType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ThermoPlasticHyperElastType instance_;

  };  // class ThermoPlasticHyperElastType

  /*----------------------------------------------------------------------*/
  //! wrapper for finite strain elasto-plastic material
  class ThermoPlasticHyperElast : public So3Material
  {
   public:
    //! construct empty material object
    ThermoPlasticHyperElast();

    //! construct the material object given material parameters
    explicit ThermoPlasticHyperElast(Mat::PAR::ThermoPlasticHyperElast* params);

    //! @name Packing and Unpacking

    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return ThermoPlasticHyperElastType::Instance().UniqueParObjectId();
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

    //! @name Access methods

    //! material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_thermoplhyperelast;
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(Inpar::STR::KinemType kinem) override
    {
      if (!(kinem == Inpar::STR::KinemType::nonlinearTotLag))
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    //! return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new ThermoPlasticHyperElast(*this));
    }

    //! Young's modulus
    double Youngs() const { return params_->youngs_; }

    //! Poisson's ratio
    double PoissonRatio() const { return params_->poissonratio_; }

    //! density
    double Density() const override { return params_->density_; }

    //! shear modulus
    double shear_mod() const { return 0.5 * (params_->youngs_) / (1.0 + params_->poissonratio_); }

    //! yield stress
    virtual double Yield() const { return params_->yield_; }

    //! isotropic hardening modulus
    virtual double IsoHard() const { return params_->isohard_; }

    //! flow stress softening
    virtual double FlowStressSoft() const { return params_->yieldsoft_; }

    //! hardening softening
    virtual double HardSoft() const { return params_->hardsoft_; }

    //! saturation hardening
    virtual double SatHardening() const { return params_->sathardening_; }

    //! coefficient of thermal expansion
    virtual double CTE() const { return params_->cte_; }

    //! coefficient of thermal expansion
    virtual double HardExpo() const { return params_->hardexpo_; }

    //! initial, reference temperature
    virtual double InitTemp() const { return params_->inittemp_; }

    //! return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    //! return accumulated strain at Gauss points
    //! use the old vector (last_) for postprocessing
    //! Output is called after(!!) Update, so that newest values are included in
    //! the old history vectors last_, while the current history vectors curr_
    //! are reset
    double AccumulatedStrain(int gp) const { return (accplstrainlast_->at(gp)); }

    //! mechanical dissipation
    double MechDiss(int gp) const { return (mechdiss_->at(gp)); }

    //! linearisation of Dmech w.r.t. temperatures T_{n+1}
    //! contribution to K_TT
    double MechDiss_kTT(int gp) const { return (mechdiss_k_tt_->at(gp)); }

    //! linearisation of the mechanical dissipation w.r.t. displacements d_{n+1}
    //! contribution to K_Td
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> MechDiss_kTd(int gp) const
    {
      return (mechdiss_k_td_->at(gp));
    }

    //! thermoplastic heating
    double ThermoPlastHeating(int gp) const { return (thrplheat_->at(gp)); }

    //! linearisation of thermoplastic heating w.r.t. temperatures T_{n+1}
    //! contribution to K_TT
    double thermo_plast_heating_k_tt(int gp) const { return (thrplheat_k_tt_->at(gp)); }

    //! linearisation of thermoplastic heating w.r.t. temperatures T_{n+1}
    //! contribution to K_Td
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> thermo_plast_heating_k_td(int gp) const
    {
      return (thrplheat_k_td_->at(gp));
    }

    //! linearisation of material tangent w.r.t. temperatures T_{n+1}
    //! contribution to K_dT
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> CMat_kdT(int gp) const { return (cmat_kd_t_->at(gp)); }

    //! check if history variables are already initialised
    bool Initialized() const { return (isinit_ and (defgrdcurr_ != Teuchos::null)); }

    //! return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    //! return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

    //@}

    //! @name Evaluation methods

    //! initialise internal stress variables
    void setup(int numgp, Input::LineDefinition* linedef) override;

    //! update internal stress variables
    void Update() override;

    //! evaluate material law
    void evaluate(const Core::LinAlg::Matrix<3, 3>*
                      defgrd,  //!< input deformation gradient for multiplicative sp
        const Core::LinAlg::Matrix<NUM_STRESS_3D, 1>*
            glstrain,                    //!< input Green-Lagrange strain (redundant with defo
                                         //   but used for neo-hooke evaluation; maybe remove
        Teuchos::ParameterList& params,  //!< input paramter list (e.g. Young's, ...)
        Core::LinAlg::Matrix<NUM_STRESS_3D, 1>*
            stress,  //!< output (mandatory) second Piola-Kirchhoff stress
        Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>*
            cmat,  //!< output (mandatory) material stiffness matrix
        int gp,    ///< Gauss point
        int eleGID) override;

    //! evaluate the elasto-plastic tangent
    void setup_cmat_elasto_plastic(Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
                                       cmat,  //!< elasto-plastic tangent modulus (out)
        double Dgamma,                        //!< plastic multiplier
        double Hiso_temp,           //!< temperature-dependent linear isotropic hardening modulus
        double sigma_y0infty_temp,  //!< temperature-dependent saturation hardening stress
        double sigma_y0_temp,       //!< temperature-dependent flow/yield stress
        double mubar,               //!< deformation-dependent shear modulus
        double q_trial,             //!< trial von Mises equivalent stress
        const Core::LinAlg::Matrix<3, 3>& defgrd,  //!< F_{n+1}
        Core::LinAlg::Matrix<3, 3> invdefgrdcurr,  //!< inverse of F_{n+1}
        Core::LinAlg::Matrix<3, 3> n,              //!< spatial flow vector
        double bulk,                               //!< bulk modulus
        int gp                                     //!< current Gauss-point
    );

    //! calculate updated value of bebar_{n+1}
    void calculate_current_bebar(const Core::LinAlg::Matrix<3, 3>& devtau,  //!< s_{n+1}
        double G,                                                           //!< shear modulus
        const Core::LinAlg::Matrix<3, 3>& id2,  //!< second-order identity
        int gp                                  //!< current Gauss-point
    );


    //! main 3D material call to determine stress and constitutive tensor ctemp
    //  originally method of fourieriso with const!!!
    void evaluate(const Core::LinAlg::Matrix<1, 1>& Ntemp,  //!< temperature of element
        Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& ctemp,  //!< temperature-dependent material tangent
        Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
            cmat_T,                                          //!< temperature-dependent mechanical
                                                             //!< material tangent
        Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& stresstemp,  //!< temperature-dependent stress term
        Teuchos::ParameterList& params                       //!< parameter list
    );

    //! computes temperature-dependent isotropic thermal elasticity tensor in
    //! matrix notion for 3d
    void setup_cthermo(
        Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& ctemp,  //!< temperature dependent material tangent
        Teuchos::ParameterList& params                  //!< parameter list
    );

    //! computes temperature-dependent isotropic mechanical elasticity tensor in
    //! matrix notion for 3d
    void SetupCmatThermo(const Core::LinAlg::Matrix<1, 1>& Ntemp,
        Core::LinAlg::Matrix<6, 6>& cmat_T, Teuchos::ParameterList& params);

    //! calculates stress-temperature modulus
    double st_modulus();

    //! finite difference check of material tangent
    void fd_check(Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& stress,  // updated stress sigma_n+1
        Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
            cmat,  // material tangent calculated with FD of stresses
        Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
            cmatFD,                               // material tangent calculated with FD of stresses
        const Core::LinAlg::Matrix<1, 1>& Ntemp,  // scalar-valued current temperature
        Teuchos::ParameterList& params            // parameter list including F,C^{-1},...
    );

    /// Return whether the material requires the deformation gradient for its evaluation
    bool needs_defgrd() override { return true; };

    //@}

   private:
    //! my material parameters
    Mat::PAR::ThermoPlasticHyperElast* params_;

    //! @name Internal / history variables

    //! plastic history variables
    //! old (i.e. at t_n)  deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<3, 3>>> defgrdlast_;
    //! current (i.e. at t_n+1) deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<3, 3>>> defgrdcurr_;

    //! old (i.e. at t_n) elastic, isochoric right Cauchy-Green tensor
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<3, 3>>> bebarlast_;
    //! current (i.e. at t_n+1) elastic, isochoric right Cauchy-Green tensor
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<3, 3>>> bebarcurr_;

    //! old (i.e. at t_n) accumulated plastic strain
    Teuchos::RCP<std::vector<double>> accplstrainlast_;
    //! current (i.e. at t_n+1) accumulated plastic strain
    Teuchos::RCP<std::vector<double>> accplstraincurr_;

    //@}

    //! @name Linearisation terms for thermal equation

    //! current (i.e. at t_n+1) mechanical dissipation
    Teuchos::RCP<std::vector<double>> mechdiss_;
    //! current (i.e. at t_n+1) linearised mechanical dissipation w.r.t. T_{n+1}
    Teuchos::RCP<std::vector<double>> mechdiss_k_tt_;
    //! current (i.e. at t_n+1) linearised mechanical dissipation w.r.t. d_{n+1}
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<6, 1>>> mechdiss_k_td_;
    //! current (i.e. at t_n+1) thermoplastic heating term
    Teuchos::RCP<std::vector<double>> thrplheat_;
    //! current (i.e. at t_n+1) thermoplastic heating term
    Teuchos::RCP<std::vector<double>> thrplheat_k_tt_;
    //! current (i.e. at t_n+1) thermoplastic heating term w.r.t. d_{n+1}
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<6, 1>>> thrplheat_k_td_;

    //@}

    //! @name Linearisation terms for structural equation

    //! current (i.e. at t_n+1) linearised material tangent w.r.t. T_{n+1}
    Teuchos::RCP<std::vector<Core::LinAlg::Matrix<6, 1>>> cmat_kd_t_;

    //@}

    //! indicator if #Initialize routine has been called
    bool isinit_;
    //! indicator if material has started to be plastic
    bool plastic_step_;
    //! element ID, in which first plasticity occurs
    int plastic_ele_id_;

  };  // class ThermoPlasticHyperElast

}  // namespace Mat

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
