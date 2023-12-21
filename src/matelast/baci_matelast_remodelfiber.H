/*----------------------------------------------------------------------*/
/*! \file
\brief Base class for fiber materials which remodel. Has a pointer to the single fiber families.

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MATELAST_REMODELFIBER_H
#define BACI_MATELAST_REMODELFIBER_H

#include "baci_config.H"

#include "baci_mat_par_parameter.H"
#include "baci_matelast_coupanisoexpo.H"
#include "baci_matelast_coupanisoexpoactive.H"
#include "baci_matelast_summand.H"

#include <utility>

BACI_NAMESPACE_OPEN

namespace MAT
{
  /// Growth evolution equation
  struct GrowthEvolution
  {
    /// Constructor
    GrowthEvolution(double const k, double const sig) : sig_h(sig), k_sig(k){};

    /// Set cauchy prestress of new mass which is deposited during G&R
    void SetSigH(double const& sig) { sig_h = sig; };

    /// Evaluate residual of growth evolution equation
    void EvaluateFunc(double& r,     ///< Residual
        double const& sig,           ///< Current Cauchy stress
        double const& cur_rho_col,   ///< Current collagen mass density
        double const& last_rho_col,  ///< "Old" (last time step) collagen mass density
        double const& dt,            ///< Time step size
        int const eleGID) const      ///< Element id
    {
      r = cur_rho_col * (1.0 - dt * k_sig * (sig - sig_h) / sig_h) - last_rho_col;
    };

    /// Evaluate derivative of growth evolution equation of the i-th fiber family w.r.t. the mass
    /// density of the i-th fiber family
    void EvaluatedFuncidrhoi(double& dridrhoi,  ///< Output
        double const& sig,                      ///< Current Cauchy stress
        double const& dsigdrho,     ///< Derivative of Cauchy stress w.r.t. the mass density
        double const& cur_rho_col,  ///< Current collagen mass density
        double const& dt,           ///< Time step size
        int const eleGID) const     ///< Element id
    {
      dridrhoi = 1.0 - dt * k_sig * (sig - sig_h + cur_rho_col * dsigdrho) / sig_h;
    };

    /// Evaluate derivative of growth evolution equation of the i-th fiber family w.r.t. the mass
    /// density of the j-th fiber family
    void EvaluatedFuncidrhoj(double& dridrhoj,  ///< Output
        double const& sig,                      ///< Current Cauchy stress
        double const& dsigdrho,     ///< Derivative of Cauchy stress w.r.t. the mass density
        double const& cur_rho_col,  ///< Current collagen mass density
        double const& dt,           ///< Time step size
        int const eleGID) const     ///< Element id
    {
      dridrhoj = -dt * k_sig * cur_rho_col * dsigdrho / sig_h;
    };

    /// Evaluate derivative of growth evolution equation w.r.t. the remodel stretch
    void EvaluatedFuncidlambr(double& drdlamb,  ///< Output
        double const& dsigdlamb,    ///< Derivative of Cauchy stress w.r.t. the mass density
        double const& cur_rho_col,  ///< Current collagen mass density
        double const& dt,           ///< Time step size
        int const eleGID) const     ///< Element id
    {
      drdlamb = -dt * k_sig * cur_rho_col * dsigdlamb / sig_h;
    };

    /// Evaluate derivative of growth evolution equation w.r.t. the right Cauchy-Green tensor
    void EvaluatedFuncidC(CORE::LINALG::Matrix<1, 6>& drdC,  ///< Output
        CORE::LINALG::Matrix<6, 1> const&
            dsigdCv,  ///< Derivative of Cauchy stress w.r.t. the right Cauchy-Green tensor
        double const& cur_rho_col,  ///< Current collagen mass density
        double const& dt,           ///< Time step size
        int const eleGID) const     ///< Element id
    {
      drdC.UpdateT(-cur_rho_col * dt * k_sig / sig_h, dsigdCv, 0.0);
    };

    /// Evaluate time derivative of the mass density of the i-th fiber family
    void Evaluatedrhodt(double& drhodt,  ///< Output
        double const& sig,               ///< Current Cauchy stress
        double const& cur_rho_col,       ///< Current collagen mass density
        int const eleGID) const          ///< Element id
    {
      drhodt = cur_rho_col * k_sig * (sig - sig_h) / sig_h;
    };

    double sig_h;        ///< Cauchy prestress of new mass which is deposited during G&R
    double const k_sig;  ///< Growth parameter
  };

  /// Remodel evolution equation
  struct RemodelEvolution
  {
    /// Constructor
    RemodelEvolution(double const k, double const sig, double const t)
        : sig_h(sig), k_sig(k), t_decay(t){};

    /// Set cauchy prestress of new mass which is deposited during G&R
    void SetSigH(double const& sig) { sig_h = sig; };

    /// Evaluate residual of remodel evolution equation
    void EvaluateFunc(double& r,                     ///< Residual
        double const& sig,                           ///< Current Cauchy stress
        CORE::LINALG::Matrix<3, 3> const& YM,        ///< Factor in remodel evolution equation
        CORE::LINALG::Matrix<3, 3> const& dsigdCeM,  ///< Derivative of local Cauchy stress w.r.t.
                                                     ///< to the elastic right Cauchy-Green tensor
        int const eleGID) const                      ///< Element id
    {
      r = (k_sig * (sig - sig_h) / sig_h + 1. / t_decay) * (sig - sig_h) - 2.0 * dsigdCeM.Dot(YM);
    };

    /// Evaluate derivative of remodel evolution equation of the i-th fiber family w.r.t. the
    /// inelastic (scalar) remodel stretch of the i-th fiber family
    void EvaluatedFuncidlambri(double& drdlambr,  ///< Output
        double const& sig,                        ///< Current Cauchy stress
        double const& dsigdlambr,  ///< Derivative of Cauchy stress w.r.t. the remodel stretch
        CORE::LINALG::Matrix<3, 3> const& YM,         ///< Factor in remodel evolution equation
        CORE::LINALG::Matrix<3, 3> const& dYdlambrM,  ///< Derivative of factor in remodel evolution
                                                      ///< equation w.r.t. the remodel stretch
        CORE::LINALG::Matrix<3, 3> const& dsigdCeM,   ///< Derivative of local Cauchy stress w.r.t.
                                                      ///< the elastic right Cauchy-Green tensor
        CORE::LINALG::Matrix<3, 3> const&
            dsigdCedlambrM,      ///< Derivative of local Cauchy stress w.r.t. the remodel stretch
        int const eleGID) const  ///< Element id
    {
      drdlambr = 2. * (k_sig / sig_h) * (sig - sig_h) * dsigdlambr;
      drdlambr += dsigdlambr / t_decay;
      drdlambr -= 2.0 * dsigdCedlambrM.Dot(YM);
      drdlambr -= 2.0 * dsigdCeM.Dot(dYdlambrM);
    };

    /// Evaluate derivative of remodel evolution equation of the i-th fiber family w.r.t. the
    /// inelastic (scalar) remodel stretch of the i-th fiber family
    void EvaluatedFuncidrho(double& drdrho,  ///< Output
        double const& sig,                   ///< Current Cauchy stress
        double const& dsigdrho,  ///< Derivative of Cauchy stress w.r.t. the remodel stretch
        CORE::LINALG::Matrix<3, 3> const& YM,        ///< Factor in remodel evolution equation
        CORE::LINALG::Matrix<3, 3> const& dYdrhoM,   ///< Derivative of factor in remodel evolution
                                                     ///< equation w.r.t. the remodel stretch
        CORE::LINALG::Matrix<3, 3> const& dsigdCeM,  ///< Derivative of local Cauchy stress w.r.t.
                                                     ///< the elastic right Cauchy-Green tensor
        CORE::LINALG::Matrix<3, 3> const&
            dsigdCedrhoM,        ///< Derivative of local Cauchy stress w.r.t. the remodel stretch
        int const eleGID) const  ///< Element id
    {
      drdrho = 2. * (k_sig / sig_h) * (sig - sig_h) * dsigdrho;
      drdrho += dsigdrho / t_decay;
      drdrho -= 2.0 * dsigdCedrhoM.Dot(YM);
      drdrho -= 2.0 * dsigdCeM.Dot(dYdrhoM);
    };

    /// Evaluate derivative of remodel evolution equation of the i-th fiber family w.r.t. the
    /// inelastic (scalar) remodel stretch of the i-th fiber family
    void EvaluatedFuncidC(CORE::LINALG::Matrix<1, 6>& drdC,  ///< Output
        double const& sig,                                   ///< Current Cauchy stress
        CORE::LINALG::Matrix<6, 1> const&
            dsigdCv,  ///< Derivative of Cauchy stress w.r.t. the right Cauchy-Green tensor
        CORE::LINALG::Matrix<6, 1> const&
            Y_strain,  ///< Factor in remodel evolution equation (strain-like Voigt notation)
        CORE::LINALG::Matrix<9, 6> const& dYdC,  ///< Derivative of factor in remodel evolution
                                                 ///< equation w.r.t. the right Cauchy-Green tensor
        CORE::LINALG::Matrix<9, 1> const&
            dsigdCe9x1,  ///< Derivative of local Cauchy stress w.r.t. to the elastic right
                         ///< Cauchy-Green tensor (9x1 vector notation)
        CORE::LINALG::Matrix<6, 6> const&
            dsigdCedC,           ///< Derivative of local Cauchy stress w.r.t. to the
                                 ///< elastic and total right Cauchy-Green tensor
        int const eleGID) const  ///< Element id
    {
      drdC.UpdateT(2. * (k_sig / sig_h) * (sig - sig_h), dsigdCv, 0.0);
      drdC.UpdateT(1.0 / t_decay, dsigdCv, 1.0);
      drdC.MultiplyTN(-2.0, dsigdCe9x1, dYdC, 1.0);
      drdC.MultiplyTN(-2.0, Y_strain, dsigdCedC, 1.0);
    };

    /// Evaluate derivative of remodel evolution equation w.r.t. time derivative of the remodel
    /// stretch
    void Evaluatedlambrdt(double& dlambrdt,       ///< Output
        double const& sig,                        ///< Current Cauchy stress
        CORE::LINALG::Matrix<3, 3> const& YredM,  ///< Reduced factor in remodel evolution equation
                                                  ///< (\dot{\lambda}_r is factored out)
        CORE::LINALG::Matrix<3, 3> const& dsigdCeM,  ///< Derivative of local Cauchy stress w.r.t.
                                                     ///< to the elastic right Cauchy-Green tensor
        int const eleGID) const                      ///< Element id
    {
      dlambrdt = ((k_sig * (sig - sig_h) / sig_h + 1. / t_decay) * (sig - sig_h)) /
                 (2.0 * dsigdCeM.Dot(YredM));
    };

    double sig_h;          ///< Cauchy prestress of new mass which is deposited during G&R
    double const k_sig;    ///< Growth parameter
    double const t_decay;  ///< Decay time of collagen
  };

  namespace ELASTIC
  {
    /// Container with all fiber specific data
    struct FiberData
    {
      /// Constructor
      FiberData(Teuchos::RCP<MAT::ELASTIC::Summand> sum) : fiber(std::move(sum)), G(0){};

      /// Update of internal data
      void UpdateNewton(int const gp, double const dt)
      {
        iFrM[gp].Update(G / cur_lambr[gp], AM, 0.0);
        iFrM[gp].Update(1. / std::sqrt(G / cur_lambr[gp]), AM_orth, 1.0);
        FrnM[gp].Update(last_lambr[gp] / G, AM, 0.0);
        FrnM[gp].Update(1. / std::sqrt(last_lambr[gp] / G), AM_orth, 1.0);
        diFrdlambrM[gp].Update(-std::pow(cur_lambr[gp], -2.0) * G, AM, 0.0);
        diFrdlambrM[gp].Update(0.5 * std::pow(cur_lambr[gp] * G, -0.5), AM_orth, 1.0);
        dFrdlambrM[gp].Update(1.0 / G, AM, 0.0);
        dFrdlambrM[gp].Update(
            -0.5 * std::pow(cur_lambr[gp], -1.5) * std::pow(G, 0.5), AM_orth, 1.0);
        FrdotM[gp].Update((cur_lambr[gp] - last_lambr[gp]) / (dt * G), AM, 0.0);
        FrdotM[gp].Update(-0.5 * std::pow(cur_lambr[gp], -1.5) * std::pow(G, 0.5) *
                              (cur_lambr[gp] - last_lambr[gp]) / dt,
            AM_orth, 1.0);
        dFrdotdlambrM[gp].Update(1.0 / (dt * G), AM, 0.0);
        dFrdotdlambrM[gp].Update(0.5 * std::pow(cur_lambr[gp], -1.5) * std::pow(G, 0.5) / dt *
                                     (1.5 * (1.0 - last_lambr[gp] / cur_lambr[gp]) - 1.0),
            AM_orth, 1.0);
      };
      /// Update of internal data
      void UpdateHistory(int const gp)
      {
        last_lambr[gp] = cur_lambr[gp];
        last_rho[gp] = cur_rho[gp];
      };

      Teuchos::RCP<MAT::ELASTIC::Summand> fiber;  ///< Type of remodel fiber
      double G;                                   ///< Prestretch of a fiber family
      std::vector<double> cur_rho;    ///< Mass density (in reference configuration) at time t_{n+1}
      std::vector<double> last_rho;   ///< Mass density (in reference configuration) at time t_{n}
      std::vector<double> cur_lambr;  ///< Inelastic stretch in fiber direction at time t_{n+1}
      std::vector<double> last_lambr;      ///< Inelastic stretch in fiber direction at time t_{n}
      CORE::LINALG::Matrix<3, 3> AM;       ///< Structural tensor: represents fiber direction
      CORE::LINALG::Matrix<3, 3> AM_orth;  ///< Orthogonal structural tensor ( 1_{ij} - A_{ij} )
      std::vector<CORE::LINALG::Matrix<3, 3>>
          iFrM;  ///< Inverse inelastic remodel deformation gradient
      std::vector<CORE::LINALG::Matrix<3, 3>>
          FrnM;  ///< Inelastic remodel deformation gradient at time t_{n}
      std::vector<CORE::LINALG::Matrix<3, 3>>
          dFrdlambrM;  ///< Derivative of inelastic remodel deformation gradient w.r.t. the
                       ///< inelastic stretch
      std::vector<CORE::LINALG::Matrix<3, 3>>
          diFrdlambrM;  ///< Derivative of inverse inelastic remodel deformation gradient w.r.t. the
                        ///< inelastic stretch
      std::vector<CORE::LINALG::Matrix<3, 3>>
          FrdotM;  ///< Time derivative of inelastic remodel deformation gradient
      std::vector<CORE::LINALG::Matrix<3, 3>>
          dFrdotdlambrM;  ///< Time derivative of inelastic remodel deformation gradient
      Teuchos::RCP<RemodelEvolution> remodel;  ///< Remodel evolution equation
      Teuchos::RCP<GrowthEvolution> growth;    ///< Growth evolution equation

      friend std::ostream& operator<<(std::ostream& output, MAT::ELASTIC::FiberData const& F)
      {
        output << "Attention: Only the first element of each vector is written to the screen! -> "
                  "gp=0 \n"
               << "Prestretch : " << F.G << "\n cur_rho: " << F.cur_rho[0]
               << "\n last_rho: " << F.last_rho[0] << "\n cur_lambr: " << F.cur_lambr[0]
               << "\n last_lambr: " << F.last_lambr[0] << "\n Structural tensor: " << F.AM
               << "\n Orthogonal structural tensor: " << F.AM_orth
               << "\n Inverse remodel deformation gradient: " << F.iFrM[0]
               << "\n Remodel deformation gradient time t_n: " << F.FrnM[0]
               << "\n Derivative of remodel deformation gradient w.r.t. lambr: " << F.dFrdlambrM[0]
               << "\n Derivative of inverse remodel deformation gradient w.r.t. lambr: "
               << F.diFrdlambrM[0]
               << "\n Time derivative of remodel deformation gradient: " << F.FrdotM[0]
               << "\n Derivative of time derivative of remodel deformation gradient w.r.t. lambr: "
               << F.dFrdotdlambrM[0] << "\n\n"
               << "RemodelEvolution: \n"
               << "Cauchy prestress: " << F.remodel->sig_h
               << "\nGrowthParameter: " << F.remodel->k_sig
               << "\n Decay time: " << F.remodel->t_decay << "\n\n"
               << "GrowthEvolution: \n"
               << "Cauchy prestress: " << F.growth->sig_h
               << "\nGrowthParameter: " << F.growth->k_sig << "\n\n";
        return output;
      };
    };


    namespace PAR
    {
      /*!
       * @brief <h3>Input line</h3>
       *
       * MAT 1 ELAST_RemodelFiber NUMMAT 1 MATIDS 100 TDECAY 1.0 SIGMAPRE 1.0
       */
      class RemodelFiber : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        RemodelFiber(const Teuchos::RCP<MAT::PAR::Material>& matdata);

        /// length of material list
        const int nummat_;

        /// the list of material IDs
        const std::vector<int>* matids_;

        /// @name material parameters
        //@{

        /// decay time of Poisson (degradation) process
        const double t_decay_;

        /// time constant for collagen growth
        const double k_growth_;

        /// initial mass fraction of each fiber family in constraint mixture
        std::vector<double>* init_w_col_;

        /// deposition stretch of collagen fibers
        const double G_;

        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<MAT::Material> CreateMaterial() override
        {
          dserror(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };

      };  // class RemodelFiber

    }  // namespace PAR

    /*!
     * @brief General interface for fibers which remodel
     */
    class RemodelFiber : public Summand
    {
     public:
      /// constructor with given material parameters
      RemodelFiber(MAT::ELASTIC::PAR::RemodelFiber* params);

      ///@name Packing and Unpacking
      //@{

      void PackSummand(CORE::COMM::PackBuffer& data) const override;

      void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position) override;

      //@}

      /// @name Access material constants
      //@{

      /// material type
      INPAR::MAT::MaterialType MaterialType() const override
      {
        return INPAR::MAT::mes_remodelfiber;
      };

      /// Return number of fiber families
      inline unsigned GetNumFibers() const { return potsumfiber_.size(); };

      /// Return current mass density of k-th specific fiber
      inline double GetCurMassDensity(unsigned const k, int const gp) const
      {
        return potsumfiber_[k]->cur_rho[gp];
      };

      //@}

      /// Update current mass density and inelastic stretch of k-th specific fiber
      inline void UpdateGrowthRemodelParameter(
          double const& drho, double const& dlambr, unsigned const k, int const gp)
      {
        potsumfiber_[k]->cur_rho[gp] += drho;
        potsumfiber_[k]->cur_lambr[gp] += dlambr;
      };

      /// Update of summand
      void Update() override;

      /*!
       * \brief Register anisotropy extensions to be passed to all summands
       *
       * \param anisotropy Reference to the anisotropy
       */
      void RegisterAnisotropyExtensions(Anisotropy& anisotropy) override;

      /// Setup of summand
      virtual void Setup(int numgp, double rho_tot, INPUT::LineDefinition* linedef);

      /// Update fiber directions with new local coordinate system (radaxicirc_)
      void UpdateFiberDirs(CORE::LINALG::Matrix<3, 3> const& locsys,  ///< local coordinate system
          const double& dt);                                          ///< time step size

      /// Update cauchy prestress of new mass which is deposited during G&R
      void UpdateSigH();

      /// Returns the necessary derivations of the growth and remodel evolution equations and
      /// their residuals
      void EvaluateDerivativesInternalNewton(
          CORE::LINALG::Matrix<3, 3> const* const defgrd,  ///< Deformation gradient
          int const nr_grf_proc,  ///< Global number of fibers which were already processed
          int const nr_grf_tot,   ///< Total number of remodel fibers
          int const gp,           ///< Current gp
          double const& dt,       ///< Time step size
          int const eleGID,       ///< Element ID
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse inelastic growth deformation gradient
          CORE::LINALG::Matrix<3, 3> const&
              dFgdrhoM,  ///< Derivative of inelastic growth deformation
                         ///< gradient w.r.t. the mass density
          CORE::LINALG::Matrix<3, 3> const&
              diFgdrhoM,  ///< Derivative of inverse inelastic growth
                          ///< deformation gradient w.r.t. the mass density
          std::vector<std::vector<double>>&
              dWdrho,  ///< Derivative of the growth evolution eq. w.r.t. the current mass density
          std::vector<std::vector<double>>&
              dWdlambr,  ///< Derivative of the growth evolution eq. w.r.t. the inelastic remodel
                         ///< fiber stretch
          std::vector<double>& W,                    ///< Residual of growth evolution eq.
          std::vector<std::vector<double>>& dEdrho,  ///< Derivative of the remodel evolution eq.
                                                     ///< w.r.t. the current mass density
          std::vector<std::vector<double>>&
              dEdlambr,  ///< Derivative of the remodel evolution eq. w.r.t. the inelastic remodel
                         ///< fiber stretch
          std::vector<double>& E);  ///< Residual of remodel evolution eq.

      /// Returns the derivations of the growth and remodel evolution equations w.r.t. the right
      /// Cauchy Green tensor
      void EvaluateDerivativesCauchyGreen(
          CORE::LINALG::Matrix<3, 3> const* const defgrd,  ///< Deformation gradient
          int const nr_grf_proc,  ///< Global number of fibers which were already processed
          int const gp,           ///< Current gp
          double const& dt,       ///< Time step size
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse inelastic growth deformation gradient
          std::vector<CORE::LINALG::Matrix<1, 6>>&
              dWdC,  ///< Derivation of the growth evolution eq. w.r.t. right Cauchy Green
          std::vector<CORE::LINALG::Matrix<1, 6>>&
              dEdC,           ///< Derivation of the remodel evolution eq. w.r.t. right Cauchy Green
          int const eleGID);  ///< Element ID

      /// Retrieve stress and cmat
      void EvaluateAdditionalGrowthRemodelCmat(
          CORE::LINALG::Matrix<3, 3> const* const defgrd,  ///< Deformation gradient
          int const nr_grf_proc,  ///< Global number of fibers which were already processed
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse inelastic growth deformation gradient
          CORE::LINALG::Matrix<3, 3> const&
              diFgdrhoM,  ///< Derivative of inverse inelastic growth deformation gradient w.r.t.
                          ///< the current mass density
          std::vector<CORE::LINALG::Matrix<1, 6>> const&
              drhodC,  ///< Derivation of the current mass density w.r.t. the right Cauchy Green
                       ///< tensor
          std::vector<CORE::LINALG::Matrix<1, 6>> const&
              dlambrdC,  ///< Derivation of the remodel fiber stretch w.r.t. the right Cauchy
                         ///< Green tensor
          CORE::LINALG::Matrix<6, 6>& cmat,  ///< Material stiffness matrix due to G&R
          int const gp,                      ///< Current gp
          int const eleGID) const;           ///< Element ID

      /// prestressing step: evaluate fibers in a "normal" way (no growth and/or remodeling)
      void EvaluateAnisotropicStressCmat(
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,                            ///< Inverse inelastic growth deformation gradient
          CORE::LINALG::Matrix<6, 6>& cmat,    ///< Material stiffness matrix
          CORE::LINALG::Matrix<6, 1>& stress,  ///< 2nd PK-stress
          int const gp,                        ///< Current gp
          double const& dt,                    ///< Time step size
          int const eleGID);                   ///< Element ID

      /// Explicit local time integration: updates current collagen density and inelastic fiber
      /// stretch
      void EvaluateGrowthAndRemodelingExpl(
          CORE::LINALG::Matrix<3, 3> const& defgrd,  ///< Deformation gradient
          double const& dt,                          ///< Time step size
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,           ///< Inverse inelastic growth deformation gradient
          const int gp,       ///< Current gp
          const int eleGID);  ///< Element ID

      /// Return names of visualization data
      virtual void VisNames(std::map<std::string, int>& names, unsigned int p);

      /// Return visualization data
      bool VisData(
          const std::string& name, std::vector<double>& data, int numgp, int eleId) override;

      /// Indicator for the chosen formulations
      void SpecifyFormulation(
          bool& isoprinc,    ///< global indicator for isotropic principal formulation
          bool& isomod,      ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,  ///< global indicator for anisotropic principal formulation
          bool& anisomod,    ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral)
          override  ///< global indicator, if one viscoelastic formulation is used
      {
        anisoprinc = true;
        return;
      };

     private:
      /// Add contribution to elasticity tensor
      void AddStressCmat(CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,                   ///< Inverse inelastic deformation gradient (used in G&R)
          FiberData const& fiberdat,  ///< Type of remodel fiber
          int const gp,               ///< Current Gauss-Point
          int const eleGID,           ///< Element Id
          CORE::LINALG::Matrix<6, 1>& stress,       ///< 2nd Piola Kirchhoff stress
          CORE::LINALG::Matrix<6, 6>& cmat) const;  ///< Elasticity tenor

      /// @name Automatic differentiation
      //@{
      /// Note: The analytical solution is always the default implementation as it is much faster

      /// First derivative w.r.t. to the right Cauchy-Green tensor (automatic differentiation)
      template <class FUNC, typename T, typename ForceAnalytical>
      void DerivdC(CORE::LINALG::Matrix<3, 3, T> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3, T> const&
              iFinM,  ///< Inverse inelastic deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3, T> const&
              AM,                        ///< Structural tensor of associated anisotropic direction
          FUNC const& func,              ///< Function which is differentiated
          int gp,                        ///< Gauss point
          ForceAnalytical const eleGID,  ///< Element Id
          CORE::LINALG::Matrix<3, 3, T>& dfuncdC) const;  ///< First derivative

      /// First derivative w.r.t. to the right Cauchy-Green tensor (analytical differentiation)
      template <class FUNC, typename T>
      void DerivdC(CORE::LINALG::Matrix<3, 3, T> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3, T> const&
              iFinM,  ///< Inverse inelastic deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3, T> const&
              AM,            ///< Structural tensor of associated anisotropic direction
          FUNC const& func,  ///< Function which is differentiated
          int gp,            ///< Gauss point
          int const eleGID,  ///< Element Id
          CORE::LINALG::Matrix<3, 3, T>& dfuncdC) const;  ///< First derivative

      /// Second derivative w.r.t. to the right Cauchy-Green tensor (automatic differentiation)
      template <class FUNC, typename T, typename ForceAnalytical>
      void DerivdCdC(CORE::LINALG::Matrix<3, 3, T> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3, T> const&
              iFinM,  ///< Inverse inelastic deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3, T> const&
              AM,                        ///< Structural tensor of associated anisotropic direction
          FUNC const& func,              ///< Function which is differentiated
          int gp,                        ///< Gauss point
          ForceAnalytical const eleGID,  ///< Element Id
          CORE::LINALG::Matrix<6, 6, T>& dfuncdCdC) const;  ///< Second derivative

      /// Second derivative w.r.t. to the right Cauchy-Green tensor (analytical differentiation)
      template <class FUNC, typename T>
      void DerivdCdC(CORE::LINALG::Matrix<3, 3, T> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3, T> const&
              iFinM,  ///< Inverse inelastic deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3, T> const&
              AM,            ///< Structural tensor of associated anisotropic direction
          FUNC const& func,  ///< Function which is differentiated
          int gp,            ///< Gauss point
          int const eleGID,  ///< Element Id
          CORE::LINALG::Matrix<6, 6, T>& dfuncdCdC) const;  ///< Second derivative

      /// Evaluate local Cauchy stress in current fiber direction (analytical differentiation)
      template <typename T>
      void EvaluateLocalCauchyStress(
          CORE::LINALG::Matrix<3, 3, T> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3, T> const&
              iFinM,  ///< Inverse inelastic deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3, T> const&
              AM,  ///< Structural tensor of associated anisotropic direction
          Teuchos::RCP<MAT::ELASTIC::Summand> const& fiber,  ///< Type of remodel fiber
          int gp,                                            ///< Gauss point
          int eleGID,                                        ///< Element Id
          T& sig) const;                                     ///< Local Cauchy stress

      /// Evaluate derivative of local Cauchy stress in current fiber direction w.r.t. the elastic
      /// right Cauchy Green tensor (analytical differentiation)
      template <typename T>
      void EvaluatedsigdCe(CORE::LINALG::Matrix<3, 3, T> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3, T> const&
              iFgM,  ///< Inverse growth deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3, T> const&
              iFrM,  ///< Inverse remodeling deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3, T> const&
              AM,  ///< Structural tensor of associated anisotropic direction
          Teuchos::RCP<MAT::ELASTIC::Summand> const& fiber,  ///< Type of remodel fiber
          int gp,                                            ///< Gauss point
          int eleGID,                                        ///< Element Id
          CORE::LINALG::Matrix<3, 3, T>& dsigdCeM)
          const;  ///< Derivative of Cauchy stress w.r.t. elastic right Cauchy-Green tensor

      /// Evaluate derivative of the derivative of the local Cauchy stress in current fiber
      /// direction w.r.t. the elastic right Cauchy Green tensor w.r.t. the right Cauchy Green
      /// tensor (automatic differentiation)
      template <typename ForceAnalytical>
      void EvaluatedsigdCedC(CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse growth deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3> const&
              iFrM,  ///< Inverse remodeling deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3> const&
              AM,  ///< Structural tensor of associated anisotropic direction
          Teuchos::RCP<MAT::ELASTIC::Summand> const fiber,  ///< Type of remodel fiber
          int gp,                                           ///< Gauss point
          ForceAnalytical const eleGID,                     ///< Element Id
          CORE::LINALG::Matrix<6, 6>& dsigdCedC) const;     ///< Output

      /// Evaluate derivative of the derivative of the local Cauchy stress in current fiber
      /// direction w.r.t. the elastic right Cauchy Green tensor w.r.t. the right Cauchy Green
      /// tensor (analytical differentiation)
      void EvaluatedsigdCedC(CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse growth deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3> const&
              iFrM,  ///< Inverse remodeling deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3> const&
              AM,  ///< Structural tensor of associated anisotropic direction
          Teuchos::RCP<MAT::ELASTIC::Summand> const&
              fiber,                                     ///< Type of remodel fiberEvaluatedsigdCedC
          int gp,                                        ///< Gauss point
          int eleGID,                                    ///< Element Id
          CORE::LINALG::Matrix<6, 6>& dsigdCedC) const;  ///< Output

      /// Derivative of local Cauchy stress w.r.t. to the right Cauchy-Green tensor (automatic
      /// differentiation)
      template <typename T, typename ForceAnalytical>
      void EvaluatedsigdC(CORE::LINALG::Matrix<3, 3, T> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3, T> const&
              iFinM,  ///< Inverse inelastic deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3, T> const&
              AM,  ///< Structural tensor of associated anisotropic direction
          Teuchos::RCP<MAT::ELASTIC::Summand> const fiber,  ///< Type of remodel fiber
          int gp,                                           ///< Gauss point
          ForceAnalytical const eleGID,                     ///< Element Id
          CORE::LINALG::Matrix<3, 3, T>& dsigdC)
          const;  ///< First derivative of local Cauchy stress
                  ///< w.r.t. the right Cauchy-Green tensor

      /// Derivative of local Cauchy stress w.r.t. to the right Cauchy-Green tensor (analytical
      /// differentiation)
      template <typename T>
      void EvaluatedsigdC(CORE::LINALG::Matrix<3, 3, T> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3, T> const&
              iFinM,  ///< Inverse inelastic deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3, T> const&
              AM,  ///< Structural tensor of associated anisotropic direction
          Teuchos::RCP<MAT::ELASTIC::Summand> const& fiber,  ///< Type of remodel fiber
          int gp,
          int eleGID,  ///< Element Id
          CORE::LINALG::Matrix<3, 3, T>& dsigdC)
          const;  ///< First derivative of local Cauchy stress
                  ///< w.r.t. the right Cauchy-Green tensor

      /// Evaluate derivative of local Cauchy stress/ derivative of the local Cauchy stress w.r.t.
      /// the elastic right Cauchy Green tensor w.r.t. the mass density (automatic
      /// differentiation)
      template <typename ForceAnalytical>
      void EvaluateDerivativesCauchyGrowth(
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse growth deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3> const&
              dFgdrhoM,  ///< Derivative of inelastic growth deformation
                         ///< gradient w.r.t. the mass density
          CORE::LINALG::Matrix<3, 3> const&
              diFgdrhoM,                 ///< Derivative of inverse inelastic growth
                                         ///< deformation gradient w.r.t. the mass density
          FiberData const& fiberdat,     ///< Data Container of fibers
          int const gp,                  ///< Current Gauss-Point
          ForceAnalytical const eleGID,  ///< Element Id
          double& dsigdrho,              ///< Derivative of Cauchy stress w.r.t. the mass density
          CORE::LINALG::Matrix<3, 3>& dsigdCedrhoM)
          const;  ///< Derivative of Cauchy stress w.r.t. elastic right Cauchy Green and w.r.t.
                  ///< the mass density

      /// Evaluate derivative of local Cauchy stress/ derivative of the local Cauchy stress w.r.t.
      /// the elastic right Cauchy Green tensor w.r.t. the mass density (analytical
      /// differentiation)
      void EvaluateDerivativesCauchyGrowth(
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse growth deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3> const&
              dFgdrhoM,  ///< Derivative of inelastic growth deformation
                         ///< gradient w.r.t. the mass density
          CORE::LINALG::Matrix<3, 3> const&
              diFgdrhoM,              ///< Derivative of inverse inelastic growth
                                      ///< deformation gradient w.r.t. the mass density
          FiberData const& fiberdat,  ///< Data Container of fibers
          int const gp,               ///< Current Gauss-Point
          int const eleGID,           ///< Element Id
          double& dsigdrho,           ///< Derivative of Cauchy stress w.r.t. the mass density
          CORE::LINALG::Matrix<3, 3>& dsigdCedrhoM)
          const;  ///< Derivative of Cauchy stress w.r.t. elastic right Cauchy Green and w.r.t.
                  ///< the mass density

      /// Evaluate derivative of local Cauchy stress/ derivative of the local Cauchy stress w.r.t.
      /// the elastic right Cauchy Green tensor w.r.t. inverse remodeling stretch (automatic
      /// differentiation)
      template <typename ForceAnalytical>
      void EvaluateDerivativesCauchyRemodel(
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,                      ///< Inverse growth deformation gradient (used in G&R)
          FiberData const& fiberdat,     ///< Data Container of fibers
          int const gp,                  ///< Current Gauss-Point
          ForceAnalytical const eleGID,  ///< Element Id
          double& dsigdlambr,  ///< Derivative of Cauchy stress w.r.t. inverse remodel deformation
                               ///< gradient
          CORE::LINALG::Matrix<3, 3>& dsigdCedlambrM)
          const;  ///< Derivative of Cauchy stress w.r.t. elastic right Cauchy Green and w.r.t.
                  ///< the inelastic remodeling stretch

      /// Evaluate derivative of local Cauchy stress/ derivative of the local Cauchy stress w.r.t.
      /// the elastic right Cauchy Green tensor w.r.t. inverse remodeling stretch (analytical
      /// differentiation)
      void EvaluateDerivativesCauchyRemodel(
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,                   ///< Inverse growth deformation gradient (used in G&R)
          FiberData const& fiberdat,  ///< Data Container of fibers
          int const gp,               ///< Current Gauss-Point
          int const eleGID,           ///< Element Id
          double& dsigdlambr,  ///< Derivative of Cauchy stress w.r.t. inverse remodel deformation
                               ///< gradient
          CORE::LINALG::Matrix<3, 3>& dsigdCedlambrM)
          const;  ///< Derivative of Cauchy stress w.r.t. elastic right Cauchy Green and w.r.t.
                  ///< the inelastic remodeling stretch

      /// Evaluates growth and remodeling evolution equation
      void EvaluateEvolutionEquation(double& rg,  ///< Residual of growth evolution equation
          double& rr,                             ///< Residual of remodeling evolution equation
          CORE::LINALG::Matrix<3, 3> const& CM,   ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,                   ///< Inverse growth deformation gradient (used in G&R)
          double const& dt,           ///< Time step size
          FiberData const& fiberdat,  ///< Container with all fiber specific data
          int const gp,               ///< Current Gauss-Point
          int const eleGID) const;    ///< Element Id

      /// Evaluates growth and remodeling evolution equation
      void EvaluateDerivativeEvolutionEquation(
          double& dWidrhoi,  ///< Derivative of growth evolution equation of the i-th fiber family
                             ///< w.r.t. the mass density of the i-th fiber family
          double& dWidrhoj,  ///< Derivative of growth evolution equation of the i-th fiber family
                             ///< w.r.t. the mass density of the j-th fiber family
          double& dWdlambr,  ///< Derivative of growth evolution equation w.r.t. the remodel stretch
          double& dEdrho,    ///< Derivative of growth evolution equation w.r.t. the mass density
          double& dEdlambr,  ///< Derivative of growth evolution equation w.r.t. the inelastic
                             ///< remodel stretch
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse growth deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3> const&
              dFgdrhoM,  ///< Derivative of inelastic growth deformation
                         ///< gradient w.r.t. the mass density
          CORE::LINALG::Matrix<3, 3> const&
              diFgdrhoM,              ///< Derivative of inverse inelastic growth
                                      ///< deformation gradient w.r.t. the mass density
          double const& dt,           ///< Time step size
          FiberData const& fiberdat,  ///< Container with all fiber specific data
          int const gp,               ///< Current Gauss-Point
          int const eleGID) const;    ///< Element Id

      /// Derivative of growth and remodel evolution equation w.r.t. the right Cauchy-Green tensor
      void EvaluatedEvolutionEquationdC(
          CORE::LINALG::Matrix<1, 6>& dWdC,      ///< Derivative of growth evolution equation
          CORE::LINALG::Matrix<1, 6>& dEdC,      ///< Derivative of remodel evolution equation
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,                   ///< Inverse growth deformation gradient (used in G&R)
          double const& dt,           ///< Time step size
          FiberData const& fiberdat,  ///< Container with all fiber specific data
          int const k,                ///< k-th fiber family
          int const gp,               ///< Current Gauss-Point
          int const eleGID);          ///< Element Id

      /// Time derivative of mass density and remodel stretch
      void EvaluatedEvolutionEquationdt(double& drhodt,  ///< Output
          double& dlambrdt,                              ///< Output
          CORE::LINALG::Matrix<3, 3> const& CM,          ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,                   ///< Inverse growth deformation gradient (used in G&R)
          FiberData const& fiberdat,  ///< Container with all fiber specific data
          int const k,                ///< k-th fiber family
          int const gp,               ///< Current Gauss-Point
          int const eleGID);          ///< Element Id

      /// Derivative of 2nd Piola Kirchhoff stress w.r.t. the current mass density and the
      /// inelastic remodel stretch (automatic differentiation)
      template <typename ForceAnalytical>
      void EvaluateDerivatives2ndPiolaKirchhoffGrowthRemodel(
          CORE::LINALG::Matrix<6, 1>& dSidrhoi,  ///< Output
          CORE::LINALG::Matrix<6, 1>& dSidrhoj,  ///< Output
          CORE::LINALG::Matrix<6, 1>& dSdlambr,  ///< Output
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse growth deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3> const&
              diFgdrhoM,  ///< Derivative of inverse inelastic growth deformation gradient w.r.t.
                          ///< the current mass density
          FiberData const& fiberdat,            ///< Container with all fiber specific data
          int const gp,                         ///< Current Gauss-Point
          ForceAnalytical const eleGID) const;  ///< Element Id

      /// Derivative of 2nd Piola Kirchhoff stress w.r.t. the current mass density and the
      /// inelastic remodel stretch (analytical differentiation)
      void EvaluateDerivatives2ndPiolaKirchhoffGrowthRemodel(
          CORE::LINALG::Matrix<6, 1>& dSidrhoi,  ///< Output
          CORE::LINALG::Matrix<6, 1>& dSidrhoj,  ///< Output
          CORE::LINALG::Matrix<6, 1>& dSdlambr,  ///< Output
          CORE::LINALG::Matrix<3, 3> const& CM,  ///< Right Cauchy-Green tensor
          CORE::LINALG::Matrix<3, 3> const&
              iFgM,  ///< Inverse growth deformation gradient (used in G&R)
          CORE::LINALG::Matrix<3, 3> const&
              diFgdrhoM,  ///< Derivative of inverse inelastic growth deformation gradient w.r.t.
                          ///< the current mass density
          FiberData const& fiberdat,  ///< Container with all fiber specific data
          int const gp,               ///< Current Gauss-Point
          int const eleGID) const;    ///< Element Id


      /// Initialize all structural tensors necessary in subsequent calculations
      virtual void SetupStructuralTensorsGR();
      //@}

      /// My material parameters
      MAT::ELASTIC::PAR::RemodelFiber* params_;

      /// Map of data containers for each fiber family including material summand
      std::vector<Teuchos::RCP<MAT::ELASTIC::FiberData>> potsumfiber_;

      /// Initial individual mass density of each fiber family at the last time step
      std::vector<double> init_rho_col_;

      /// Current fiber Cauchy stress
      std::vector<std::vector<double>> cauchystress_;
    };  // namespace PAR

  }  // namespace ELASTIC
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif  // MATELAST_REMODELFIBER_H
