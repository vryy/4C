/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_NP.H

\brief evaluation of ScaTra elements for Nernst-Planck ion-transport equations

\level 2

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_ELCH_NP_H
#define SCATRA_ELE_CALC_ELCH_NP_H

#include "scatra_ele_calc_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchNP;

    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcElchNP : public ScaTraEleCalcElch<distype>
    {
     public:
      virtual ~ScaTraEleCalcElchNP() { return; };

      /// Singleton access method
      static ScaTraEleCalcElchNP<distype>* Instance(const int numdofpernode, const int numscal,
          const std::string& disname, const ScaTraEleCalcElchNP* delete_me = NULL);

      /// called upon destruction
      virtual void Done();

     private:
      typedef ScaTraEleCalc<distype> my;
      typedef ScaTraEleCalcElch<distype> myelch;

      /// private constructor, since we are a Singleton.
      ScaTraEleCalcElchNP(const int numdofpernode, const int numscal, const std::string& disname);

      /*========================================================================*/
      //! @name general framework
      /*========================================================================*/

      //! calculate contributions to matrix and rhs (inside loop over all scalars)
      void CalcMatAndRhs(Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,                 //!< element rhs to calculate+
          const int k,                                    //!< index of current scalar
          const double fac,                               //!< domain-integration factor
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double taufac,  //!< tau times domain-integration factor
          const double
              timetaufac,  //!< domain-integration factor times tau times time-integration factor
          const double rhstaufac,  //!< time-integration factor for rhs times tau times
                                   //!< domain-integration factor
          LINALG::Matrix<my::nen_, 1>&
              tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
          double& rhsint  //!< rhs at Gauss point
      );

      //! calculate contributions to matrix and rhs (outside loop over all scalars)
      void CalcMatAndRhsOutsideScalarLoop(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,  //!< element rhs to calculate
          const double fac,                //!< domain-integration factor
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double rhsfac  //!< time-integration factor for rhs times domain-integration factor
      );

      //! Correction for additional flux terms / currents across Dirichlet boundaries
      void CorrectionForFluxAcrossDC(DRT::Discretization& discretization,  //!< discretization
          const std::vector<int>& lm,                                      //!< location vector
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs   //!< element rhs to calculate
      );

      /*========================================================================*/
      //! @name material and related functions
      /*========================================================================*/

      //! get the material parameters
      void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                  //!< density at t_(n)
          std::vector<double>& densnp,                 //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,                 //!< density at t_(n+alpha_M)
          double& visc,                                //!< fluid viscosity
          const int iquad = -1                         //!< id of current gauss point (default = -1)
      );


      //! evaluate material
      void Materials(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< index of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      /*========================================================================*/
      //! @name stabilization and related functions
      /*========================================================================*/

      //! Calculate quantities used for stabilization
      void PrepareStabilization(
          std::vector<double>& tau,  //!< stabilization parameters (one per transported scalar)
          std::vector<LINALG::Matrix<my::nen_, 1>>&
              tauderpot,  //!< derivatives of stabilization parameters w.r.t. electric potential
          const std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_f)
          const double vol                    //!< element volume
      );

      //! Calculate derivative of tau w.r.t. electric potential according to Taylor, Hughes and
      //! Zarins
      void CalcTauDerPotTaylorHughesZarins(
          LINALG::Matrix<my::nen_, 1>&
              tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
          double& tau,    //!< stabilization parameter
          const double densnp,         //!< density at t_(n+1)
          const double frt,            //!< F/(RT)
          const double diffusvalence,  //!< diffusion coefficient times valence
          const LINALG::Matrix<my::nsd_, 1>&
              veleff  //!< effective convective velocity (fluid velocity plus migration velocity if
                      //!< applicable)
      );

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! CalcRes: Residual of Nernst-Planck equation in strong form
      double CalcRes(const int k,  //!< index of current scalar
          const double conint,     //!< concentration at GP
          const double hist,       //!< history value at GP
          const double convphi,    //!< convective term (without convective part of migration term)
          const double frt,        //!< F/(RT)
          const LINALG::Matrix<my::nen_, 1>&
              migconv,         //!< migration operator: -F/(RT) \grad{\Phi} * \grad{N}
          const double rhsint  //!< rhs of Nernst-Planck equation (not of Newton-Raphson scheme) at
                               //!< Gauss point
      );

      //! CalcMat: SUPG Stabilization of convective term due to fluid flow and migration
      void CalcMatConvStab(Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          const int k,                                      //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double taufac,      //!< stabilization parameter tau times domain-integration factor
          const double
              timetaufac,  //!< domain-integration factor times tau times time-integration factor
          LINALG::Matrix<my::nen_, 1>&
              tauderpot,     //!< derivatives of stabilization parameter w.r.t. electric potential
          const double frt,  //!< F/(RT)
          const LINALG::Matrix<my::nen_, 1>&
              conv,  //!< convection operator: u_x*N,x + u_y*N,y + u_z*N,z
          const LINALG::Matrix<my::nen_, 1>&
              migconv,          //!< migration operator: -F/(RT) \grad{\Phi} * \grad{N}
          const double conint,  //!< concentration at GP
          const LINALG::Matrix<my::nsd_, 1>& gradphi,  //!< gradient of concentration at GP
          const double residual  //!< residual of Nernst-Planck equation in strong form
      );

      //! CalcMat: Migration term
      void CalcMatMigr(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                  //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double frt,         //!< F/RT
          const LINALG::Matrix<my::nen_, 1>& migconv,  //!< migration operator
          const double conint                          //!< concentration at GP
      );

      //! CalcMat: Electroneutrality condition in PDE form
      void CalcMatPotEquENCPDE(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                          //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double frt,         //!< F/RT
          const LINALG::Matrix<my::nen_, 1>& migconv,  //!< migration operator
          const double conint                          //!< concentration at GP
      );

      //! CalcMat: Electroneutrality condition in PDE form with Nernst-Planck equation for species m
      //! eliminated
      void CalcMatPotEquENCPDEElim(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                              //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double frt,         //!< F/RT
          const LINALG::Matrix<my::nen_, 1>& migconv,  //!< migration operator
          const double conint                          //!< concentration at GP
      );

      //! CalcMat: Poisson equation for electric potential
      void CalcMatPotEquPoisson(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                           //!< index of current scalar
          const double fac,                                      //!< domain-integration factor
          const double epsilon,                                  //!< dielectric constant
          const double faraday                                   //!< Faraday constant
      );

      //! CalcMat: Laplace equation for electric potential
      void CalcMatPotEquLaplace(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const double fac                                       //!< domain-integration factor
      );

      //! CalcRhs: Additional contributions from conservative formulation of Nernst-Planck equations
      void CalcRhsConvAddCons(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                         //!< index of current scalar
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double conint,  //!< concentration at GP
          const double vdiv     //!< velocity divergence
      );

      //! CalcRhs: SUPG Stabilization of convective term due to fluid flow and migration
      void CalcRhsConvStab(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                      //!< index of current scalar
          const double rhstaufac,  //!< time-integration factor for rhs times tau times
                                   //!< domain-integration factor
          const LINALG::Matrix<my::nen_, 1>&
              conv,  //!< convection operator: u_x*N,x + u_y*N,y + u_z*N,z
          const LINALG::Matrix<my::nen_, 1>&
              migconv,           //!< migration operator: -F/(RT) \grad{\Phi} * \grad{N}
          const double residual  //!< residual of Nernst-Planck equation in strong form
      );

      //! CalcRhs: Migration term
      void CalcRhsMigr(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                  //!< index of current scalar
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const LINALG::Matrix<my::nen_, 1>& migconv,  //!< migration operator
          const double conint                          //!< concentration at GP
      );

      //! CalcRhs: Electroneutrality condition in PDE form
      void CalcRhsPotEquENCPDE(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                          //!< index of current scalar
          const double fac,                                     //!< domain-integration factor
          const LINALG::Matrix<my::nen_, 1>& migconv,           //!< migration operator
          const double conint,                                  //!< concentration at GP
          const LINALG::Matrix<my::nsd_, 1>& gradphi            //!< gradient of concentration at GP
      );

      //! CalcRhs: Electroneutrality condition in PDE form with Nernst-Planck equation for species m
      //! eliminated
      void CalcRhsPotEquENCPDEElim(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                              //!< index of current scalar
          const double fac,                                         //!< domain-integration factor
          const LINALG::Matrix<my::nen_, 1>& migconv,               //!< migration operator
          const double conint,                                      //!< concentration at GP
          const LINALG::Matrix<my::nsd_, 1>& gradphi  //!< gradient of concentration at GP
      );

      //! CalcRhs: Poisson equation for electric potential
      void CalcRhsPotEquPoisson(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                           //!< index of current scalar
          const double fac,                                      //!< domain-integration factor
          const double epsilon,                                  //!< dielectric constant
          const double faraday,                                  //!< Faraday constant
          const double conint,                                   //!< concentration at GP
          const LINALG::Matrix<my::nsd_, 1>& gradpot             //!< gradient of potential at GP
      );

      //! CalcRhs: Laplace equation for electric potential
      void CalcRhsPotEquLaplace(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const double fac,                                      //!< domain-integration factor
          const LINALG::Matrix<my::nsd_, 1>& gradpot             //!< gradient of potential at GP
      );

      /*========================================================================*/
      //! @name additional service routines
      /*========================================================================*/

      //! validity check with respect to input parameters, degrees of freedom, number of scalars
      //! etc.
      void CheckElchElementParameter(DRT::Element* ele  //!< current element
      );

      //! evaluate an electrode boundary kinetics point condition
      void EvaluateElchBoundaryKineticsPoint(const DRT::Element* ele,  ///< current element
          Epetra_SerialDenseMatrix& emat,                              ///< element matrix
          Epetra_SerialDenseVector& erhs,  ///< element right-hand side vector
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
              ephinp,  ///< state variables at element nodes
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
              ehist,                          ///< history variables at element nodes
          double timefac,                     ///< time factor
          Teuchos::RCP<DRT::Condition> cond,  ///< electrode kinetics boundary condition
          const int nume,                     ///< number of transferred electrons
          const std::vector<int> stoich,      ///< stoichiometry of the reaction
          const int kinetics,                 ///< desired electrode kinetics model
          const double pot0,                  ///< electrode potential on metal side
          const double frt,                   ///< factor F/RT
          const double
              scalar  ///< scaling factor for element matrix and right-hand side contributions
      );

      // Get conductivity from material
      void GetConductivity(const enum INPAR::ELCH::EquPot
                               equpot,  //!< type of closing equation for electric potential
          double& sigma_all,            //!< conductivity of electrolyte solution
          std::vector<double>&
              sigma,    //!< conductivity of a single ion + overall electrolyte solution
          bool effCond  //!< flag if effective conductivity should be calculated
      );

      //!  calculate weighted mass flux (no reactive flux so far) -> elch-specific implementation
      void CalculateFlux(LINALG::Matrix<my::nsd_, 1>& q,  //!< flux of species k
          const INPAR::SCATRA::FluxType fluxtype,         //!< type fo flux
          const int k                                     //!< index of current scalar
      );

      //! calculate error of numerical solution with respect to analytical solution
      void CalErrorComparedToAnalytSolution(
          const DRT::Element* ele,          //!< the element we are dealing with
          Teuchos::ParameterList& params,   //!< parameter list
          Epetra_SerialDenseVector& errors  //!< vector containing L2-error norm
      );

      //! set internal variables for Nernst-Planck formulation
      void SetInternalVariablesForMatAndRHS();

      //! get internal variable manager for Nernst-Planck formulation
      Teuchos::RCP<ScaTraEleInternalVariableManagerElchNP<my::nsd_, my::nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerElchNP<my::nsd_, my::nen_>>(
            my::scatravarmanager_);
      };

      /*========================================================================*/
      //! @name flags and enums
      /*========================================================================*/

      //! flag for migration operator in the stabilization terms
      /* The migration term in the ion transport equation can be split into a convective and a
         reactive part. Combining the convective part of the migration term with the convective term
         in the ion transport equation yields an "extended" convective term including not only the
         fluid velocity, but also the migration velocity. In case migrationstab_ == true, SUPG
         stabilization is performed for the extended convective term. In case migrationstab_ ==
         false, SUPG stabilization is performed only for the original convective term.
       */
      bool migrationstab_;

      /*========================================================================*/
      //! @name scalar degrees of freedom and related
      /*========================================================================*/
    };


    /// ScaTraEleInternalVariableManagerElchNP implementation
    /*!
      This class keeps all internal variables for the Nernst-Planck formulation.
    */
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchNP
        : public ScaTraEleInternalVariableManagerElch<NSD, NEN>
    {
     public:
      typedef ScaTraEleInternalVariableManager<NSD, NEN> vm;
      typedef ScaTraEleInternalVariableManagerElch<NSD, NEN> vmelch;

      ScaTraEleInternalVariableManagerElchNP(
          int numscal, const DRT::ELEMENTS::ScaTraEleParameterElch* elchpara)
          : ScaTraEleInternalVariableManagerElch<NSD, NEN>(numscal, elchpara),
            // constant internal variables
            // empty

            //! internal variables evaluated at the Gauss point
            migvelint_(true),
            migconv_(true)
      {
        return;
      }

      virtual ~ScaTraEleInternalVariableManagerElchNP() { return; };

      //! compute and set internal variables for the Nernst-Planck formulation
      void SetInternalVariablesElchNP(
          const LINALG::Matrix<NEN, 1>& funct,  //!< array for shape functions
          const LINALG::Matrix<NSD, NEN>&
              derxy,  //!< global derivatives of shape functions w.r.t x,y,z
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ephinp,  //!< nodal state variables at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>& ephin,  //!< nodal state variables at t_(n)
          const LINALG::Matrix<NSD, NEN>&
              econvelnp,  //!< nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
          const std::vector<LINALG::Matrix<NEN, 1>>&
              ehist  //!< history vector of transported scalars
      )
      {
        // set internal variables in base variable manager
        vmelch::SetInternalVariablesElch(funct, derxy, ephinp, ephin, econvelnp, ehist);

        // migration velocity vector (divided by D_k*z_k) at t_(n+1) or t_(n+alpha_F): -F/(RT)
        // \grad{\Phi}
        migvelint_.Multiply(-vmelch::frt_, derxy, ephinp[vm::numscal_]);

        // convective part of migration term (divided by D_k*z_k): -F/(RT) \grad{\Phi}*\grad{N}
        migconv_.MultiplyTN(-vmelch::frt_, derxy, vmelch::gradpot_);

        return;
      };

      /*========================================================================*/
      //! @name return constant internal variables
      /*========================================================================*/

      // empty

      /*========================================================================*/
      //! @name return methods for internal variables
      /*========================================================================*/

      //! return migration velocity vector (divided by D_k*z_k): -F/(RT) \grad{\Phi}
      const LINALG::Matrix<NSD, 1>& MigVelInt() const { return migvelint_; };

      //! return convective part of migration term (divided by D_k*z_k): -F/(RT) \grad{\Phi} *
      //! \grad{N}
      const LINALG::Matrix<NEN, 1> MigConv() { return migconv_; };

     private:
      /*========================================================================*/
      //! @name constant internal variables
      /*========================================================================*/


      /*========================================================================*/
      //! @name internal variables evaluated at element center or Gauss point
      /*========================================================================*/

      //! migration velocity vector (divided by D_k*z_k): -F/(RT) \grad{\Phi}
      LINALG::Matrix<NSD, 1> migvelint_;

      //! convective part of migration term (divided by D_k*z_k): -F/(RT) \grad{\Phi} * \grad{N}
      LINALG::Matrix<NEN, 1> migconv_;

    };  // class ScaTraEleInternalVariableManagerElchNP

  }  // namespace ELEMENTS

}  // namespace DRT

#endif
