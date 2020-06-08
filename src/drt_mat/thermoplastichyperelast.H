/*----------------------------------------------------------------------*/
/*! \file
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following finite strain
       von-Mises plasticity with nonlinear isotropic hardening and general
       hyperelasticity (for the time being: NeoHooke).

       implementation is based on
       Simo and Miehe: "Associative coupled thermoplasticity at finite strains:
       Formulation, numerical analysis and implementation", in Computer Methods
       in Applied Mechanics and Engineering, 98:41–104, 1992.

       geometrically nonlinear, for finite strains, rate-independent, thermo-plastic

       example input line:
       [mm,ms,kg,K,GPa]
       MAT 1 MAT_Struct_ThrPlasticHyperElast YOUNG 206.9 NUE 0.29 DENS 7.8e-6
         CTE 1e-5 INITTEMP 293 YIELD 0.45 ISOHARD 0.12924 SATHARDENING 0.715
         HARDEXPO 16.93 YIELDSOFT 0.002 HARDSOFT 0.002 DISSFACT 0.9 TOL 1.0e-06

         Seitz 11/16: There are still some linearizations for the plastic
         heating terms missing, use with caution.

\level 3

\maintainer Jan Schnabel
*/
/*----------------------------------------------------------------------*
 | definitions                                               dano 03/13 |
 *----------------------------------------------------------------------*/
#ifndef THERMO_PLAST_HYPERELAST_H_
#define THERMO_PLAST_HYPERELAST_H_

/*----------------------------------------------------------------------*
 | headers                                                   dano 03/13 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    //! material parameters for neo-Hooke
    class ThermoPlasticHyperElast : public Parameter
    {
     public:
      //! standard constructor
      ThermoPlasticHyperElast(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~ThermoPlasticHyperElast() { ; }

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ThermoPlasticHyperElast

  }  // namespace PAR


  class ThermoPlasticHyperElastType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "ThermoPlasticHyperElastType"; }

    static ThermoPlasticHyperElastType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

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
    explicit ThermoPlasticHyperElast(MAT::PAR::ThermoPlasticHyperElast* params);

    //! @name Packing and Unpacking

    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
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
    virtual void Pack(DRT::PackBuffer& data) const;

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
    virtual void Unpack(const std::vector<char>& data);

    //@}

    //! @name Access methods

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_thermoplhyperelast;
    }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new ThermoPlasticHyperElast(*this));
    }

    //! Young's modulus
    double Youngs() const { return params_->youngs_; }

    //! Poisson's ratio
    double PoissonRatio() const { return params_->poissonratio_; }

    //! density
    virtual double Density() const { return params_->density_; }

    //! shear modulus
    double ShearMod() const { return 0.5 * (params_->youngs_) / (1.0 + params_->poissonratio_); }

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
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

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
    double MechDiss_kTT(int gp) const { return (mechdiss_kTT_->at(gp)); }

    //! linearisation of the mechanical dissipation w.r.t. displacements d_{n+1}
    //! contribution to K_Td
    LINALG::Matrix<NUM_STRESS_3D, 1> MechDiss_kTd(int gp) const { return (mechdiss_kTd_->at(gp)); }

    //! thermoplastic heating
    double ThermoPlastHeating(int gp) const { return (thrplheat_->at(gp)); }

    //! linearisation of thermoplastic heating w.r.t. temperatures T_{n+1}
    //! contribution to K_TT
    double ThermoPlastHeating_kTT(int gp) const { return (thrplheat_kTT_->at(gp)); }

    //! linearisation of thermoplastic heating w.r.t. temperatures T_{n+1}
    //! contribution to K_Td
    LINALG::Matrix<NUM_STRESS_3D, 1> ThermoPlastHeating_kTd(int gp) const
    {
      return (thrplheat_kTd_->at(gp));
    }

    //! linearisation of material tangent w.r.t. temperatures T_{n+1}
    //! contribution to K_dT
    LINALG::Matrix<NUM_STRESS_3D, 1> CMat_kdT(int gp) const { return (Cmat_kdT_->at(gp)); }

    //! check if history variables are already initialised
    bool Initialized() const { return (isinit_ and (defgrdcurr_ != Teuchos::null)); }

    //! return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    //! return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    //@}

    //! @name Evaluation methods

    //! initialise internal stress variables
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    //! update internal stress variables
    virtual void Update();

    //! evaluate material law
    virtual void Evaluate(
        const LINALG::Matrix<3, 3>* defgrd,  //!< input deformation gradient for multiplicative sp
        const LINALG::Matrix<NUM_STRESS_3D, 1>*
            glstrain,                    //!< input Green-Lagrange strain (redundant with defo
                                         //   but used for neo-hooke evaluation; maybe remove
        Teuchos::ParameterList& params,  //!< input paramter list (e.g. Young's, ...)
        LINALG::Matrix<NUM_STRESS_3D, 1>*
            stress,  //!< output (mandatory) second Piola-Kirchhoff stress
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>*
            cmat,  //!< output (mandatory) material stiffness matrix
        int gp,    ///< Gauss point
        int eleGID);

    //! evaluate the elasto-plastic tangent
    void SetupCmatElastoPlastic(LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
                                    cmat,  //!< elasto-plastic tangent modulus (out)
        double Dgamma,                     //!< plastic multiplier
        double Hiso_temp,           //!< temperature-dependent linear isotropic hardening modulus
        double sigma_y0infty_temp,  //!< temperature-dependent saturation hardening stress
        double sigma_y0_temp,       //!< temperature-dependent flow/yield stress
        double mubar,               //!< deformation-dependent shear modulus
        double q_trial,             //!< trial von Mises equivalent stress
        const LINALG::Matrix<3, 3>& defgrd,  //!< F_{n+1}
        LINALG::Matrix<3, 3> invdefgrdcurr,  //!< inverse of F_{n+1}
        LINALG::Matrix<3, 3> n,              //!< spatial flow vector
        double bulk,                         //!< bulk modulus
        int gp                               //!< current Gauss-point
    );

    //! calculate updated value of bebar_{n+1}
    void CalculateCurrentBebar(const LINALG::Matrix<3, 3>& devtau,  //!< s_{n+1}
        double G,                                                   //!< shear modulus
        const LINALG::Matrix<3, 3>& id2,                            //!< second-order identity
        int gp                                                      //!< current Gauss-point
    );


    //! main 3D material call to determine stress and constitutive tensor ctemp
    //  originally method of fourieriso with const!!!
    void Evaluate(const LINALG::Matrix<1, 1>& Ntemp,  //!< temperature of element
        LINALG::Matrix<NUM_STRESS_3D, 1>& ctemp,      //!< temperature-dependent material tangent
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat_T,  //!< temperature-dependent mechanical
                                                               //!< material tangent
        LINALG::Matrix<NUM_STRESS_3D, 1>& stresstemp,  //!< temperature-dependent stress term
        Teuchos::ParameterList& params                 //!< parameter list
    );

    //! computes temperature-dependent isotropic thermal elasticity tensor in
    //! matrix notion for 3d
    void SetupCthermo(
        LINALG::Matrix<NUM_STRESS_3D, 1>& ctemp,  //!< temperature dependent material tangent
        Teuchos::ParameterList& params            //!< parameter list
    );

    //! computes temperature-dependent isotropic mechanical elasticity tensor in
    //! matrix notion for 3d
    void SetupCmatThermo(const LINALG::Matrix<1, 1>& Ntemp, LINALG::Matrix<6, 6>& cmat_T,
        Teuchos::ParameterList& params);

    //! calculates stress-temperature modulus
    double STModulus();

    //! finite difference check of material tangent
    void FDCheck(LINALG::Matrix<NUM_STRESS_3D, 1>& stress,  // updated stress sigma_n+1
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
            cmat,  // material tangent calculated with FD of stresses
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
            cmatFD,                         // material tangent calculated with FD of stresses
        const LINALG::Matrix<1, 1>& Ntemp,  // scalar-valued current temperature
        Teuchos::ParameterList& params      // parameter list including F,C^{-1},...
    );

    /// Return whether the material requires the deformation gradient for its evaluation
    virtual bool NeedsDefgrd() { return true; };

    //@}

   private:
    //! my material parameters
    MAT::PAR::ThermoPlasticHyperElast* params_;

    //! @name Internal / history variables

    //! plastic history variables
    //! old (i.e. at t_n)  deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> defgrdlast_;
    //! current (i.e. at t_n+1) deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> defgrdcurr_;

    //! old (i.e. at t_n) elastic, isochoric right Cauchy-Green tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> bebarlast_;
    //! current (i.e. at t_n+1) elastic, isochoric right Cauchy-Green tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> bebarcurr_;

    //! old (i.e. at t_n) accumulated plastic strain
    Teuchos::RCP<std::vector<double>> accplstrainlast_;
    //! current (i.e. at t_n+1) accumulated plastic strain
    Teuchos::RCP<std::vector<double>> accplstraincurr_;

    //@}

    //! @name Linearisation terms for thermal equation

    //! current (i.e. at t_n+1) mechanical dissipation
    Teuchos::RCP<std::vector<double>> mechdiss_;
    //! current (i.e. at t_n+1) linearised mechanical dissipation w.r.t. T_{n+1}
    Teuchos::RCP<std::vector<double>> mechdiss_kTT_;
    //! current (i.e. at t_n+1) linearised mechanical dissipation w.r.t. d_{n+1}
    Teuchos::RCP<std::vector<LINALG::Matrix<6, 1>>> mechdiss_kTd_;
    //! current (i.e. at t_n+1) thermoplastic heating term
    Teuchos::RCP<std::vector<double>> thrplheat_;
    //! current (i.e. at t_n+1) thermoplastic heating term
    Teuchos::RCP<std::vector<double>> thrplheat_kTT_;
    //! current (i.e. at t_n+1) thermoplastic heating term w.r.t. d_{n+1}
    Teuchos::RCP<std::vector<LINALG::Matrix<6, 1>>> thrplheat_kTd_;

    //@}

    //! @name Linearisation terms for structural equation

    //! current (i.e. at t_n+1) linearised material tangent w.r.t. T_{n+1}
    Teuchos::RCP<std::vector<LINALG::Matrix<6, 1>>> Cmat_kdT_;

    //@}

    //! indicator if #Initialize routine has been called
    bool isinit_;
    //! indicator if material has started to be plastic
    bool plastic_step_;
    //! element ID, in which first plasticity occurs
    int plastic_eleID_;

  };  // class ThermoPlasticHyperElast

}  // namespace MAT

/*----------------------------------------------------------------------*/
#endif  // THERMO_PLAST_HYPERELAST_H_
