/*----------------------------------------------------------------------*/
/*!
\file elast_remodelfiber.H
\brief the input line should read MAT 1 ELAST_RemodelFiber NUMMAT 1 MATIDS 100 TDECAY 1.0 SIGMAPRE 1.0 COLMASSFRAC 0.5

\level 3

\maintainer Fabian Braeu

*----------------------------------------------------------------------*/
#ifndef ELAST_REMODELFIBER_H_
#define ELAST_REMODELFIBER_H_

#include "elast_summand.H"
#include "elast_coupanisoexpoactive.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_RemodelFiber NUMMAT 1 MATIDS 100 TDECAY 1.0 SIGMAPRE 1.0
    class RemodelFiber : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      RemodelFiber(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~RemodelFiber() { ; }

      /// length of material list
      const int nummat_;

      /// the list of material IDs
      const std::vector<int>* matids_;

      /// @name material parameters
      //@{

      /// decay time of Poisson (degradation) process
      const double tdecay_;

      /// time constant for collagen growth
      const double k_growth_;

      /// initial mass fraction of each fiber family in constraint mixture
      std::vector<double>* init_w_col_;

      /// deposition stretch of collagen fibers
      const double G_;

      //@}

      /// create material instance of matching type with my parameters
      ///
      /// Remark:                                                  AN 03/2014
      /// Function has to be implemented, since it is a pure virtual
      /// function in the MAT::PAR::Material class.
      /// Function stays empty, since the call within the material factory is done only
      /// respective the elasthyper class
      ///
      virtual Teuchos::RCP<MAT::Material> CreateMaterial(){ return Teuchos::null; };

    };  // class RemodelFiber

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// General interface for fibers which remodel

  class RemodelFiber : public Summand
  {

  public:

    /// constructor with given material parameters
    RemodelFiber(MAT::ELASTIC::PAR::RemodelFiber* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_remodelfiber; };

    /// Return number of fiber families
    virtual unsigned GetNumFibers() { return (potsumfiberpas_.size()+potsumfiberact_.size()); };

    /// Return current mass density of one specific fiber
    virtual double GetCurMassDensity(const unsigned& k,
                                     const int& gp)
    { return cur_rho_col_[k][gp]; };

    //@}

    /// Set current mass density of one specific fiber
    virtual void UpdateCurMassDensity(const double& drho,
                                      const unsigned& k,
                                      const int& gp)
    { cur_rho_col_[k][gp] += drho; };

    /// Set current inelastic stretch in fiberdirection
    virtual void UpdateCurRemodelStretch(const double& dlamb,
                                         const unsigned& k,
                                         const int& gp)
    { cur_lambda_r_[k][gp] += dlamb; };

    /// Update of summand
    virtual void Update();

    /// Setup of summand
    virtual void Setup(int numgp,double rho_tot,DRT::INPUT::LineDefinition* linedef);

    /// Returns the necessary derivations of the growth and remodel evolution equations and their residuals
    virtual void EvaluateDerivativesInternalNewton(
      const LINALG::Matrix<3,3>* defgrd,            ///< deformation gradient
      const LINALG::Matrix<3,3> id,                 ///< identity tensor
      const int nr_grf_proc,                        ///< global number of fibers which were already processed
      const double density,                         ///< total density of the constraint mixture
      const int gp,                                 ///< current gp
      const double dt,                              ///< time step size
      const double v,                               ///< current volume change
      LINALG::Matrix<3,3> FgM,                      ///< inelastic growth deformation gradient
      LINALG::Matrix<3,3> iFgM,                     ///< inverse inelastic growth deformation gradient
      LINALG::Matrix<3,3> AgM,                      ///< structural tensor of growth direction
      LINALG::Matrix<3,3> AcirM,                    ///< structural tensor of circumferential direction
      LINALG::Matrix<3,3> AradM,                    ///< structural tensor of radial direction
      LINALG::Matrix<3,3> AaxM,                     ///< structural tensor of axial direction
      std::vector<std::vector<double> >& dWdrho,    ///< derivation of the growth evolution eq. w.r.t. the current mass density
      std::vector<std::vector<double> >& dWdlamb,   ///< derivation of the growth evolution eq. w.r.t. the inelastic remodel fiber stretch
      std::vector<double>& W,                       ///< growth evolution eq.
      std::vector<std::vector<double> >& dEdrho,    ///< derivation of the remodel evolution eq. w.r.t. the current mass density
      std::vector<std::vector<double> >& dEdlamb,   ///< derivation of the remodel evolution eq. w.r.t. the inelastic remodel fiber stretch
      std::vector<double>& E,                       ///< remodel evolution eq.
      const int eleGID,                             ///< element ID
      const int growthtype                          ///< flag to decide what type of collagen growth is used: 1: anisotropic growth; 2: isotropic growth
      );

    /// Returns the derivations of the growth and remodel evolution equations w.r.t. the right Cauchy Green tensor
    virtual void EvaluateDerivativesCauchyGreen(
      const LINALG::Matrix<3,3>* defgrd,            ///< deformation gradient
      const LINALG::Matrix<3,3> id,                 ///< identity tensor
      const int nr_grf_proc,                        ///< global number of fibers which were already processed
      const int gp,                                 ///< current gp
      const double dt,                              ///< time step size
      LINALG::Matrix<3,3> FgM,                      ///< inelastic growth deformation gradient
      LINALG::Matrix<3,3> iFgM,                     ///< inverse inelastic growth deformation gradient
      std::vector<LINALG::Matrix<1,6> >& dWdC,      ///< derivation of the growth evolution eq. w.r.t. right Cauchy Green
      std::vector<LINALG::Matrix<1,6> >& dEdC,      ///< derivation of the remodel evolution eq. w.r.t. right Cauchy Green
      const int eleGID                              ///< element ID
      );

    /// Retrieve stress and cmat
    virtual void AddStressCmatGrowthRemodel(
      const LINALG::Matrix<3,3>* defgrd,            ///< deformation gradient
      const LINALG::Matrix<3,3> id,                 ///< identity tensor
      const int nr_grf_tot,                         ///< total number of fibers which grow and remodel
      const int nr_grf_proc,                        ///< global number of fibers which were already processed
      const int gp,                                 ///< current gp
      const double v,                               ///< current volume change
      const double density,                         ///< total density of the constraint mixture
      LINALG::Matrix<3,3> iFgM,                     ///< inverse inelastic growth deformation gradient
      LINALG::Matrix<3,3> AgM,                      ///< structural tensor of growth direction
      LINALG::Matrix<3,3> AcirM,                    ///< structural tensor of circumferential direction
      LINALG::Matrix<3,3> AradM,                    ///< structural tensor of radial direction
      LINALG::Matrix<3,3> AaxM,                     ///< structural tensor of axial direction
      std::vector<LINALG::Matrix<1,6> > drhodC,     ///< derivation of the current mass density w.r.t. the right Cauchy Green tensor
      std::vector<LINALG::Matrix<1,6> > dlambdC,    ///< derivation of the remodel fiber stretch w.r.t. the right Cauchy Green tensor
      LINALG::Matrix<6,1>& stress,                  ///< 2nd PK-stress
      LINALG::Matrix<6,6>& cmat,                    ///< material stiffness matrix
      const int eleGID,                             ///< element ID
      const int growthtype                          ///< flag to decide what type of collagen growth is used: 1: anisotropic growth; 2: isotropic growth
      );

    /// Explicit local time integration: updates current collagen density and inelastic fiber stretch
    virtual void EvaluateGrowthAndRemodelingExpl(
      const LINALG::Matrix<3,3> defgrd,             ///< deformation gradient
      const double dt,                              ///< time step size
      LINALG::Matrix<3,3> iFgM,                     ///< inverse inelastic growth deformation gradient
      const int gp,                                 ///< current gp
      const int eleGID                              ///< element ID
      );

    /// prestressing step: evaluate fibers in a "normal" way (no growth and/or remodeling)
    virtual void EvaluateAnisotropicStressCmat(
      const LINALG::Matrix<3,3>* defgrd,   ///< deformation gradient
      const LINALG::Matrix<3,3> iFg,       ///< inverse inelastic growth deformation gradient
      LINALG::Matrix<6,6>& cmat,           ///< material stiffness matrix
      LINALG::Matrix<6,1>& stress,         ///< 2nd PK-stress
      const int gp,                        ///< current gp
      const int eleGID                     ///< element ID
      );

    /// Converts matrix in "matrix notation" to "stress-like" Voigt notation
    void MatrixtoStressVoigtNotationVector(const LINALG::Matrix<3,3>& in,
                                           LINALG::Matrix<6,1>& out);

    /// This function evaluates the following equation
    /// \f$out_{abcd} = fac * ( in_{ca} in_{db} + in_{da} in_{cb} )\f$
    /// minor symmetries in ab and cd
    void EvaldCedC(LINALG::Matrix<6,6>& out,
        const LINALG::Matrix<3,3>& in,
        const double fac);

    /// This function avaluates the following equation
    /// \f$out_{abcd} = fac * ( in1_{ca} in2_{db} + in1_{da} in2_{cb} + in1_{db} in2_{ca} + in1_{cb} in2_{da} )\f$
    /// minor symmetries in ab and cd
    void EvaldYdC(LINALG::Matrix<6,6>& out,
        const LINALG::Matrix<3,3>& in1,
        const LINALG::Matrix<3,3>& in2,
        const double fac);

    /// add a contribution to the derivative of the PK2 stress with respect to the inelastic growth deformation gradient
    /// in index notation this method does:
    /// \f$dPK2dFpinv_{IJKL} += fac*( A_{IK} * B_{JL} + A_{JK} * B_{IL} )\f$
    /// Note that the result is symmetric within the first two indices
    /// Thus the first two indices are stored in stress-like Voigt notation.
    /// The 2nd index pair has no symmetries and is therefore stored as 9-vector
    /// Be aware that this is NOT symmetric w.r.t. A and B !!!
    void AddtodPK2diFg(LINALG::Matrix<6,9>& out,
        LINALG::Matrix<3,3> A,
        LINALG::Matrix<3,3> B,
        double fac);

    /// This function adds the following contribution to the matrix out_ABCD (9x6)
    /// in index notation this method does:
    /// \f$out_{ABCD} = fac*( A_{AC} * B_{BD} +A_{AD} * B_{BC} )\f$
    /// Note that the result is symmetric within the last two indices
    /// Thus the last two indices are stored in stress-like Voigt notation.
    /// The 1st index pair has no symmetries and is therefore stored as 9-vector
    /// Be aware that this is NOT symmetric w.r.t. A and B !!!
    void AddtoMatrixHalfSymProd(const double& fac,                    ///< scalar factor
                                const LINALG::Matrix<3,3>& A,         ///< Matrix in the above mentioned equation
                                const LINALG::Matrix<3,3>& B,         ///< Matrix in the above mentioned equation
                                LINALG::Matrix<9,6>& X_ABCD);         ///< Output

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names, unsigned int p);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleId);

    /// Indicator for the chosen formulations
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& viscogeneral       ///< global indicator, if one viscoelastic formulation is used
    )
    {
      anisoprinc = true;
      return;
    }

    /// Initialize all structural tensors necessary in subsequent calculations
    virtual void SetupStructuralTensorsGR();


  private:

    /// my material parameters
    MAT::ELASTIC::PAR::RemodelFiber* params_;

    /// last inelastic stretch in fiberdirection for each Gauß point
    std::vector<std::vector<double> > last_lambda_r_;

    /// current inelastic stretch in fiberdirection for each Gauß point
    std::vector<std::vector<double> > cur_lambda_r_;

    /// current individual volume of each fiber family (initially 1)
    std::vector<std::vector<double> > cur_rho_col_;

    /// individual mass density of each fiber family at the last time step
    std::vector<std::vector<double> > last_rho_col_;

    /// initial individual mass density of each fiber family at the last time step
    std::vector<double> init_rho_col_;

    /// current fiber Cauchy stress
    std::vector<std::vector<double> > stress_;

    /// map to materials/potential summands (passive! fibers) which remodel: fibers with no active component
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand> > potsumfiberpas_;

    /// map to materials/potential summands (active! fibers) which remodel: fibers have in addition to their passive component an active one
    std::vector<Teuchos::RCP<MAT::ELASTIC::CoupAnisoExpoActive> > potsumfiberact_;

    /// structural tensor (matrix, stress-like and strain-like Voigt notation)
    std::vector<LINALG::Matrix<3,3> > AM_;
    std::vector<LINALG::Matrix<6,1> > Av_;
    std::vector<LINALG::Matrix<6,1> > A_strain_;

    /// structural tensor in 9x1 vector-like notation
    std::vector<LINALG::Matrix<9,1> > A9x1_;

    /// orthogonal structural tensor ( 1_{ij} - A_{ij} )
    std::vector<LINALG::Matrix<3,3> > AM_orth_;

    /// Cauchy prestress of new mass which is deposited during G&R
    std::vector<double> sigmapre_;
  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ELAST_REMODELFIBER_H_
