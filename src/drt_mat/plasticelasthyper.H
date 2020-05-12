/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains the hyperelastic toolbox with application to finite
strain plasticity using a semi-smooth Newton method. It allows summing up
several summands of isotropic non-splitted type to build
a hyperelastic strain energy function.

The input line should read
MAT 1 MAT_PlasticElastHyper NUMMAT 1 MATIDS 2 DENS 1.0 INITYIELD 0.45 ISOHARD 0.12924
EXPISOHARD 16.93 INFYIELD 0.715 KINHARD 0.0 CTE 1.0e-5 INITTEMP 293 YIELDSOFT 0.002 HARDSOFT 0.002
VISC 1e-4 VISC_TEMP 0.003 PL_SPIN_CHI -50 rY_11 1.0 rY_22 0.9 rY_33 0.9 rY_12 0.7 rY_23 0.57385
rY_13 0.7

\maintainer Jan Schnabel
\level 2
*/

/*----------------------------------------------------------------------*/
#ifndef PLASTICELASTHYPER_H
#define PLASTICELASTHYPER_H



#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "elasthyper.H"
#include "../drt_inpar/inpar_tsi.H"

#define AS_CONVERGENCE_TOL 1.e-12

// forward declaration due to avoid header definition
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
  }

  // forward declaration
  class PlasticElastHyper;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of hyperelastic materials
    ///
    /// Storage map of hyperelastic summands.
    class PlasticElastHyper : public MAT::PAR::ElastHyper
    {
      friend class MAT::PlasticElastHyper;

     public:
      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      PlasticElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~PlasticElastHyper() { ; }

      /// @name material parameters
      //@{

      /// provide access to material/summand by its ID
      Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(
          const int id  ///< ID to look for in collection of summands
          ) const;

      // elasthyper needs this for stat inverse analysis
      virtual void OptParams(std::map<std::string, int>* pnames)
      {
        dserror("This material is not implemented in OptParams.");
      };

      /// initial yield stress
      const double inityield_;
      /// linear isotropic hardening parameter
      const double isohard_;
      /// exponent for nonlinear isotropic hardening
      const double expisohard_;
      /// saturation yield stress for nonlinear isotropic hardening
      const double infyield_;
      /// linear kinematic hardening parameter
      const double kinhard_;
      /// Visco-Plasticity parameter 'eta' in Perzyna model
      const double visc_;
      /// rate-dependency in Perzyna model
      const double rate_dependency_;
      /// Visco-Plasticity parameter 'eta' in Perzyna model softening
      const double visc_soft_;
      //! coefficient of thermal expansion
      const double cte_;
      //! initial, reference temperature
      const double inittemp_;
      //! yield stress softening
      const double yieldsoft_;
      //! hardening softening
      const double hardsoft_;
      //! taylor quinney factor
      const double taylor_quinney_;
      /// plastic spin parameter
      const double plspin_chi_;
      /// parameters of Hill yield criterion (relative yield stresses)
      const double rY_11_;
      const double rY_22_;
      const double rY_33_;
      const double rY_12_;
      const double rY_23_;
      const double rY_13_;
      /// complementarity parameter
      double cpl_;
      /// stabilization parameter "s" controlling the shape of the NCP function
      double stab_s_;
      /// method to calculate plastic dissipation in TSI problem
      INPAR::TSI::DissipationMode dis_mode_;


      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();
      //@}

    };  // class PlasticElastHyper

  }  // namespace PAR

  class PlasticElastHyperType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "PlasticElastHyperType"; }

    static PlasticElastHyperType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static PlasticElastHyperType instance_;
  };


  /*----------------------------------------------------------------------*/
  /// Collection of hyperelastic materials
  class Material;

  class PlasticElastHyper : public MAT::ElastHyper
  {
   public:
    /// construct empty material object
    PlasticElastHyper();

    /// construct the material object given material parameters
    explicit PlasticElastHyper(MAT::PAR::PlasticElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const
    {
      return PlasticElastHyperType::Instance().UniqueParObjectId();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    virtual void Pack(DRT::PackBuffer& data) const;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_plelasthyper; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new PlasticElastHyper(*this));
    }

    /// material mass density
    virtual double Density() const { return MatParams()->density_; }

    /// initial yield stress
    virtual double Inityield() const { return MatParams()->inityield_; }

    /// linear isotropic hardening modulus
    virtual double Isohard() const { return MatParams()->isohard_; }

    /// exponent for nonlinear isotropic hardening
    virtual double Expisohard() const { return MatParams()->expisohard_; }

    /// saturation yield stress for nonlinear isotropic hardening
    virtual double Infyield() const { return MatParams()->infyield_; }

    /// linear kinematic hardening modulus
    virtual double Kinhard() const { return MatParams()->kinhard_; }

    /// complementarity parameter
    virtual double cpl() const
    {
      if (MatParams() != NULL)
        return MatParams()->cpl_;
      else
        return 0.;
    }

    /// stabilization parameter
    virtual double s() const
    {
      if (MatParams() != NULL)
        return MatParams()->stab_s_;
      else
        return 0.;
    }

    /// plastic spin unequal zero
    virtual bool HavePlasticSpin() const
    {
      return (MatParams()->plspin_chi_ != 0. &&
              (MatParams()->kinhard_ != 0 || MatParams()->rY_11_ != 0.));
    }

    /// plastic spin parameter chi
    virtual double PlSpinChi() const { return MatParams()->plspin_chi_; }

    /// viscosity (of the plastic flow, no visco-elasticity)
    virtual double Visc() const { return MatParams()->visc_; }

    /// rate dependency of visco-plasticity
    virtual double ViscRate() const { return MatParams()->rate_dependency_; }

    /// viscosity (of the plastic flow, no visco-elasticity)
    virtual double ViscSoft() const { return MatParams()->visc_soft_; }

    /// thermal expansion coefficient
    virtual double Cte() const { return MatParams()->cte_; }

    /// initial temperature
    virtual double InitTemp() const { return MatParams()->inittemp_; }

    /// initial yield stress softening with temperature
    virtual double YieldSoft() const { return MatParams()->yieldsoft_; }

    /// isotropic Hardening softening with temperature
    virtual double HardSoft() const { return MatParams()->hardsoft_; }

    /// Taylor Quinney factor
    virtual double TaylorQuinney() const { return MatParams()->taylor_quinney_; }

    /// set dissipation mode
    virtual void SetDissipationMode(INPAR::TSI::DissipationMode mode)
    {
      if (MatParams() != NULL) MatParams()->dis_mode_ = mode;
    }

    /// get dissipation mode
    virtual INPAR::TSI::DissipationMode DisMode() const
    {
      if (MatParams() != NULL)
        return MatParams()->dis_mode_;
      else
        return (INPAR::TSI::DissipationMode)0;
    }

    /// evaluate quantities for elastic stiffness matrix
    /// in consideration of plastic history/deformation
    virtual void EvaluateElast(const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<3, 3>* deltaLp, LINALG::Matrix<6, 1>* pk2, LINALG::Matrix<6, 6>* cmat,
        const int gp,
        const int eleGID);  ///< global ID of element

    /// evaluate stresses and stiffness contribution
    /// due to thermal expansion
    virtual void EvaluateThermalStress(const LINALG::Matrix<3, 3>* defgrd, const double temp,
        LINALG::Matrix<6, 1>* pk2, LINALG::Matrix<6, 6>* cmat, const int gp,
        const int eleGID);  ///< global ID of element

    /// evaluate stresses and stiffness contribution
    /// due to thermal expansion
    virtual void EvaluateCTvol(const LINALG::Matrix<3, 3>* defgrd, LINALG::Matrix<6, 1>* cTvol,
        LINALG::Matrix<6, 6>* dCTvoldE, const int gp,
        const int eleGID);  ///< global ID of element

    /// evaluate the Gough Joule Effect
    virtual void EvaluateGoughJoule(
        double j, int gp, int eleGID, double& he_fac, double& he_fac_deriv);

    /// evaluate everything needed for the condensation of the plastic deformation
    /// at element level. (with zero plastic spin)
    virtual void EvaluatePlast(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>* deltaDp,  ///< symmetric part of plastic flow increment
        const double* temp,                   ///< current temperature
        Teuchos::ParameterList& params,       ///< Container for additional information
        LINALG::Matrix<6, 6>* dPK2dDp,        ///< derivative of PK2 w.r.t. F_p^{-1}
        LINALG::Matrix<6, 1>* NCP,            ///< NCP function
        LINALG::Matrix<6, 6>* dNCPdC,         ///< derivative of NCP function w.r.t. RCG
        LINALG::Matrix<6, 6>* dNCPdDp,        ///< derivative of NCP function w.r.t. deltaLp
        bool* active,                         ///< gauss point is active
        bool* elast,                          ///< gauss point needs condensation if it is not elast
        bool* as_converged,  ///< convergence of active set (false, if as has changed)
        const int gp,        ///< gauss point
        LINALG::Matrix<6, 1>*
            dNCPdT,  ///< derivative of NCP function w.r.t. temperature (only in TSI case)
        LINALG::Matrix<6, 1>* dHdC,   ///< derivative of Heating w.r.t. RCG (only in TSI case)
        LINALG::Matrix<6, 1>* dHdDp,  ///< derivative of Heating w.r.t. deltaLp (only in TSI case)
        const double dt,              ///< time step size
        const int eleGID,             ///< global ID of element
        LINALG::Matrix<6, 1>* cauchy = NULL, LINALG::Matrix<6, 6>* d_cauchy_ddp = NULL,
        LINALG::Matrix<6, 6>* d_cauchy_dC = NULL, LINALG::Matrix<6, 9>* d_cauchy_dF = NULL,
        LINALG::Matrix<6, 1>* d_cauchy_dT = NULL);

    /// evaluate everything needed for the condensation of the plastic deformation
    /// at element level. (with plastic spin)
    virtual void EvaluatePlast(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>* deltaLp,  ///< plastic deformation gradient (non-symmetric)
        const double* temp,                   ///< current temperature
        Teuchos::ParameterList& params,       ///< Container for additional information
        LINALG::Matrix<6, 9>* dPK2dLp,        ///< derivative of PK2 w.r.t. F_p^{-1}
        LINALG::Matrix<9, 1>* NCP,            ///< NCP function
        LINALG::Matrix<9, 6>* dNCPdC,         ///< derivative of NCP function w.r.t. RCG
        LINALG::Matrix<9, 9>* dNCPdLp,        ///< derivative of NCP function w.r.t. deltaLp
        bool* active,                         ///< gauss point is active
        bool* elast,                          ///< gauss point needs condensation if it is not elast
        bool* as_converged,  ///< convergence of active set (false, if as has changed)
        const int gp,        ///< gauss point
        LINALG::Matrix<9, 1>*
            dNCPdT,  ///< derivative of NCP function w.r.t. temperature (only in TSI case)
        LINALG::Matrix<6, 1>* dHdC,   ///< derivative of Heating w.r.t. RCG (only in TSI case)
        LINALG::Matrix<9, 1>* dHdLp,  ///< derivative of Heating w.r.t. deltaLp (only in TSI case)
        const double dt,              ///< time step size
        const int eleGID,             ///< global ID of element
        LINALG::Matrix<6, 1>* cauchy = NULL, LINALG::Matrix<6, 9>* d_cauchy_ddp = NULL,
        LINALG::Matrix<6, 6>* d_cauchy_dC = NULL, LINALG::Matrix<6, 9>* d_cauchy_dF = NULL,
        LINALG::Matrix<6, 1>* d_cauchy_dT = NULL);

    virtual void EvaluateCauchyPlast(const LINALG::Matrix<3, 1>& dPI,
        const LINALG::Matrix<6, 1>& ddPII, const LINALG::Matrix<3, 3>* defgrd,
        LINALG::Matrix<6, 1>& cauchy, LINALG::Matrix<6, 9>& d_cauchy_dFpi,
        LINALG::Matrix<6, 6>& d_cauchy_dC, LINALG::Matrix<6, 9>& d_cauchy_dF,
        LINALG::Matrix<6, 1>* d_cauchy_dT = NULL);

    /// hyperelastic stress response plus elasticity tensor
    /// (pure virtual in material base class. Not allowed here)
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,                  ///< Green-Lagrange strain
        Teuchos::ParameterList& params,  ///< Container for additional information
        LINALG::Matrix<6, 1>* stress,    ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,      ///< Constitutive matrix
        int gp,                          ///< Gauss point
        int eleGID)                      ///< Element GID
    {
      ElastHyper::Evaluate(defgrd, glstrain, params, stress, cmat, gp, eleGID);
      return;
    }

    virtual double StrainEnergy(const LINALG::Matrix<3, 3>& defgrd, const int gp, const int eleGID)
    {
      return StrainEnergyTSI(defgrd, gp, eleGID, MatParams()->inittemp_);
    }
    virtual double StrainEnergyTSI(
        const LINALG::Matrix<3, 3>& defgrd, const int gp, const int eleGID, const double temp);

    /// setup material data
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /*!
     * \brief Post setup routine, will be called after the complete inout is a´ready befire the
     * first evaluate call
     *
     * \param params Container for additional information
     * \param eleGID Global element id
     */
    void PostSetup(const Teuchos::ParameterList params, int eleGID);

    /// setup material TSI data
    virtual void SetupTSI(const int numgp, const int numdofperelement, const bool eas,
        const INPAR::TSI::DissipationMode mode);

    /// setup plastic orthotropy tensor H
    virtual void SetupHillPlasticity(DRT::INPUT::LineDefinition* linedef);

    /// update sumands
    virtual void Update()
    {
      dserror(
          "Elastic summands in PlasticHyperElast are not allowed to have internal variables"
          "that would need an Update-routine!");
    };

    /// update plastic history variables
    virtual void UpdateGP(const int gp, const LINALG::Matrix<3, 3>* deltaDp);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return MatParams(); }

    /// Access to material params
    virtual MAT::PAR::PlasticElastHyper* MatParams() const { return params_; }

    /// Return whether the material requires the deformation gradient for its evaluation
    virtual bool NeedsDefgrd() { return true; };

    /// get plastic algorithm parameters
    virtual void GetParams(double s, double cpl)
    {
      if (MatParams() == NULL) return;  // ... from post processor
      MatParams()->stab_s_ = s;
      MatParams()->cpl_ = cpl;
    };

    /// return accumulated plastic strain at GP
    virtual double AccumulatedStrain(int gp) { return last_alpha_isotropic_[gp]; }

    /// is this GP active
    virtual bool Active(int gp) { return activity_state_[gp]; }

    /// heating at this gp
    virtual double& HepDiss(int gp) { return (*HepDiss_)[gp]; }

    /// derivative of heating at this gp
    virtual Epetra_SerialDenseVector& dHepDissDd(int gp) { return (*dHepDissdd_)[gp]; }

    // derivative of heating w.r.t. temperature
    virtual double& dHepDT(int gp) { return (*dHepDissdT_)[gp]; }

    // derivative of heating at each gp w.r.t. nodal temperature vector
    // (only EAS contribution)
    virtual Teuchos::RCP<std::vector<Epetra_SerialDenseVector>> dHepDTeas()
    {
      return dHepDissdTeas_;
    }

    //! return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    //! return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    //! purely elastic material (input via negative initial yield stress
    bool AllElastic() { return Inityield() < 0.; }

   protected:
    /// calculates the kinematic quantities and tensors used afterwards
    virtual int EvaluateKinQuantPlast(const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<3, 3>* deltaLp, const int gp, Teuchos::ParameterList& params);

    virtual void EvaluateKinQuantElast(
        const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<3, 3>* deltaLp, const int gp);

    virtual double NormStressLike(const LINALG::Matrix<6, 1>& stress);

    /// calculates the isotropic stress and elasticity tensor for coupled configuration
    virtual void EvaluateIsotropicPrincPlast(LINALG::Matrix<6, 9>& dPK2dFpinvIsoprinc,
        LINALG::Matrix<3, 3>& MandelStressIsoprinc, LINALG::Matrix<6, 6>& dMdCisoprinc,
        LINALG::Matrix<6, 9>& dMdFpinvIsoprinc, const LINALG::Matrix<3, 1>& gamma,
        const LINALG::Matrix<8, 1>& delta);

    /// calculates the isotropic stress and elasticity tensor for coupled configuration
    virtual void EvaluateIsotropicPrincElast(LINALG::Matrix<6, 1>& stressisoprinc,
        LINALG::Matrix<6, 6>& cmatisoprinc, LINALG::Matrix<3, 1> dPI, LINALG::Matrix<6, 1> ddPII);

    virtual void EvaluateNCP(const LINALG::Matrix<3, 3>* mStr, const LINALG::Matrix<6, 6>* dMdC,
        const LINALG::Matrix<6, 9>* dMdFpinv, const LINALG::Matrix<6, 9>* dPK2dFpinv,
        const LINALG::Matrix<3, 3>* deltaDp, const int gp, const double* temp,
        LINALG::Matrix<6, 1>* NCP, LINALG::Matrix<6, 6>* dNCPdC, LINALG::Matrix<6, 6>* dNCPdDp,
        LINALG::Matrix<6, 1>* dNCPdT, LINALG::Matrix<6, 6>* dPK2dDp, bool* active, bool* elast,
        bool* as_converged, LINALG::Matrix<6, 1>* dHdC, LINALG::Matrix<6, 1>* dHdDp,
        Teuchos::ParameterList& params, const double dt, const LINALG::Matrix<6, 9>* d_cauchy_dFpi,
        LINALG::Matrix<6, 6>* d_cauchy_ddp);

    virtual void EvaluateNCPandSpin(const LINALG::Matrix<3, 3>* mStr,
        const LINALG::Matrix<6, 6>* dMdC, const LINALG::Matrix<6, 9>* dMdFpinv,
        const LINALG::Matrix<6, 9>* dPK2dFpinv, const LINALG::Matrix<3, 3>* deltaLp, const int gp,
        LINALG::Matrix<9, 1>* NCP, LINALG::Matrix<9, 6>* dNCPdC, LINALG::Matrix<9, 9>* dNCPdLp,
        LINALG::Matrix<6, 9>* dPK2dLp, bool* active, bool* elast, bool* as_converged,
        const double dt, const LINALG::Matrix<6, 9>* d_cauchy_dFpi,
        LINALG::Matrix<6, 9>* d_cauchy_ddp);

    /**
     * \brief The derivatives of the strain energy function w.r.t. the principle invariants are
     * filled up to the third derivative
     *
     * \param prinv[in]       principle invariants of the tensor
     * \param gp[in]          Gauss point
     * \param eleGID[in]      global ID of the element
     * \param dPI[in,out]     first derivative of strain energy function w.r.t. principle invariants
     * \param ddPII[in,out]   second derivative of strain energy function w.r.t. principle
                              invariants
     * \param dddPIII[in,out] third derivative of strain energy function w.r.t. principle invariants
     * \param temp[in]        temperature
     */
    void EvaluateCauchyDerivs(const LINALG::Matrix<3, 1>& prinv, const int gp, int eleGID,
        LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII, LINALG::Matrix<10, 1>& dddPIII,
        const double* temp = nullptr) override;

    /**
     * \brief Evaluate the thermal dependency of the linearizations of the cauchy stress
     *
     * \param prinv[in]      principle invariants of the tensor
     * \param ndt[in]        dot product of vectors n and t (\f[\bold{t} \cdot \bold{t}\f])
     * \param bdndt[in]      left Cauchy-Green tensor contracted using the vector n and t (\f[
                             \bold{b} \cdot \bold{n} \cdot \bold{t} \f])
     * \param ibdndt[in]     inverse of left Cauchy-Green tensor contracted using the vector n and t
                             (\f[ \bold{b}^{-1} \cdot \bold{n} \cdot \bold{t} \f])
     * \param temp[in]       temperature
     * \param DsntDT[out]    derivative of cauchy stress contracted with vectors n and t w.r.t.
                             temperature (\f[\frac{\mathrm{d} \bold{\sigma} \cdot \bold{n} \cdot
                             \bold{t} }{ \mathrm{d} T }\f])
     * \param iFTV[in]       inverse transposed of the deformation gradient
     * \param DbdndtDFV[in]  derivative of bdndt w.r.t. deformation gradient (\f[\frac{ \mathrm{d}
                             \bold{b} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{F} }\f])
     * \param DibdndtDFV[in] derivative of ibdndt w.r.t. deformation gradient (\f[\frac{ \mathrm{d}
                             \bold{b}^{-1} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{F} }\f])
     * \param DI1DF[in]      derivative of first invariant w.r.t. deformation gradient
     * \param DI2DF[in]      derivative of second invariant w.r.t. deformation gradient
     * \param DI3DF[in]      derivative of third invariant w.r.t. deformation gradient
     * \param D2sntDFDT[out] second derivative of cauchy stress contracted with vectors n and t
                             w.r.t. temperature and deformation gradient (\f[ \frac{\mathrm{d}^2
                             \bold{\sigma} \cdot \bold{n} \cdot \bold{t}} {\mathrm{d} \bold{F}
                             \mathrm{d} T } \f])
     */
    void EvaluateCauchyTempDeriv(const LINALG::Matrix<3, 1>& prinv, const double ndt,
        const double bdndt, const double ibdndt, const double* temp, double* DsntDT,
        const LINALG::Matrix<9, 1>& iFTV, const LINALG::Matrix<9, 1>& DbdndtDFV,
        const LINALG::Matrix<9, 1>& DibdndtDFV, const LINALG::Matrix<9, 1>& DI1DF,
        const LINALG::Matrix<9, 1>& DI2DF, const LINALG::Matrix<9, 1>& DI3DF,
        LINALG::Matrix<9, 1>* D2sntDFDT) override;

    virtual void AddThermalExpansionDerivs(const LINALG::Matrix<3, 1>& prinv,
        LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII, int gp, int eleGID,
        const double& temp);

    //! calculate the exponential of a 3x3 matrix (symmetric or non-symmetric)
    virtual void MatrixExponential3x3(LINALG::Matrix<3, 3>& MatrixInOut);

    //! calculate the derivative of the exponential of a symmetric 3x3 matrix
    virtual void MatrixExponentialDerivativeSym3x3(
        const LINALG::Matrix<3, 3> MatrixIn, LINALG::Matrix<6, 6>& MatrixExpDeriv);

    //! calculate the derivative of the exponential of a non-symmetric 3x3 matrix
    virtual void MatrixExponentialDerivative3x3(
        const LINALG::Matrix<3, 3> MatrixIn, LINALG::Matrix<9, 9>& MatrixExpDeriv);

    /// my material parameters
    MAT::PAR::PlasticElastHyper* params_;

    /// plastic anisotropy tensor for Hill-plasticity
    /// Von Mises plasticity is included for H = H^-1 = P_dev
    LINALG::Matrix<6, 6> PlAniso_full_;
    LINALG::Matrix<6, 6> InvPlAniso_full_;

    /// inverse plastic deformation gradient for each Gauss point at last converged state
    std::vector<LINALG::Matrix<3, 3>> last_plastic_defgrd_inverse_;

    /// accumulated plastic strain for each Gauss point at last converged state
    std::vector<double> last_alpha_isotropic_;

    /// accumulated plastic strain for each Gauss point at last converged state
    std::vector<LINALG::Matrix<3, 3>> last_alpha_kinematic_;

    /// classification, if the Gauss point is currently in the active (true) or inactive (false) set
    std::vector<bool> activity_state_;

    /// isotropic hardening increment over this time step
    std::vector<double> delta_alpha_i_;

    /// TSI infos ***************************************************
    /// use the material to transfer linearization from the structural to the thermo element

    /// Elasto-plastic heating and mechanical dissipation at each gp
    Teuchos::RCP<std::vector<double>> HepDiss_;

    /// derivative of Elasto-plastic heating and mechanical dissipation at each gp w.r.t. nodal
    /// displacements compute the complete derivative w.r.t. nodal displacements (not only RCG) to
    /// make sure, that the same element technology is used.
    Teuchos::RCP<std::vector<Epetra_SerialDenseVector>> dHepDissdd_;

    /// derivative of Elasto-plastic heating and mechanical dissipation at each gp w.r.t. gp
    /// temperature
    Teuchos::RCP<std::vector<double>> dHepDissdT_;

    /// derivative of Elasto-plastic heating and mechanical dissipation at each gp w.r.t. to element
    /// temperature this is an additional term to dHepDissdT_ that only appears in EAS elements
    Teuchos::RCP<std::vector<Epetra_SerialDenseVector>> dHepDissdTeas_;

    //  private:

    /// kinematic quatities
    static LINALG::Matrix<6, 1> Cpi_, CpiCCpi_, ircg_, Ce_, Ce2_, id2V_, bev_, be2v_;
    static LINALG::Matrix<3, 3> invpldefgrd_, CeM_, id2_, CpiC_, FpiCe_, FpiTC_, CeFpiTC_, be_,
        be2_, Fe_, beF_, beFe_, FCpi_, beFCpi_;
    static LINALG::Matrix<9, 1> CFpiCei_, CFpi_, CFpiCe_;
    static LINALG::Matrix<3, 1> prinv_;
  };

}  // namespace MAT

#endif
