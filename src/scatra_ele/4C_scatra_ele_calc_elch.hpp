/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for ion-transport equation

\level 2

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_ELCH_HPP
#define FOUR_C_SCATRA_ELE_CALC_ELCH_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"
#include "4C_scatra_ele_parameter_elch.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElch;
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElch;
    template <CORE::FE::CellType distype>
    class ScaTraEleUtilsElch;

    // class implementation
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
    class ScaTraEleCalcElch : public ScaTraEleCalc<distype, probdim>
    {
     protected:
      /// protected constructor, since we are a singleton
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalcElch(const int numdofpernode, const int numscal, const std::string& disname);

      using my = ScaTraEleCalc<distype, probdim>;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      /// In this class we do not define a static ScaTraEle...* Instance
      /// since only derived child classes are free to be allocated!!

      //! evaluate the element
      int Evaluate(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

      //! evaluate action
      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

     protected:
      /*========================================================================*/
      //! @name general framework
      /*========================================================================*/

      //! Prepare everything what is needed in CallMatAndRhs() to calculate the sysmat and the rhs
      void Sysmat(DRT::Element* ele,                  //!< the element we are dealing with
          CORE::LINALG::SerialDenseMatrix& emat,      //!< element matrix to calculate
          CORE::LINALG::SerialDenseVector& erhs,      //!< element rhs to calculate
          CORE::LINALG::SerialDenseVector& subgrdiff  //!< subgrid-diff.-scaling vector
          ) override;

      //! calculate contributions to matrix and rhs (inside of loop over all scalars)
      virtual void CalcMatAndRhs(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to calculate
          CORE::LINALG::SerialDenseVector& erhs,  //!< element rhs to calculate+
          const int k,                            //!< index of current scalar
          const double fac,                       //!< domain-integration factor
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double taufac,  //!< tau times domain-integration factor
          const double
              timetaufac,  //!< domain-integration factor times tau times time-integration factor
          const double rhstaufac,  //!< time-integration factor for rhs times tau times
                                   //!< domain-integration factor
          CORE::LINALG::Matrix<nen_, 1>&
              tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
          double& rhsint  //!< rhs at Gauss point
          ) = 0;

      //! calculate contributions to matrix and rhs (outside of loop over all scalars)
      virtual void calc_mat_and_rhs_outside_scalar_loop(
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to calculate
          CORE::LINALG::SerialDenseVector& erhs,  //!< element rhs to calculate
          const double fac,                       //!< domain-integration factor
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double rhsfac  //!< time-integration factor for rhs times domain-integration factor
          ) = 0;

      //! Correction for additional flux terms / currents across Dirichlet boundaries
      virtual void correction_for_flux_across_dc(
          DRT::Discretization& discretization,    //!< discretization
          const std::vector<int>& lm,             //!< location vector
          CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to calculate
          CORE::LINALG::SerialDenseVector& erhs   //!< element rhs to calculate
          ) = 0;

      //! finite difference check for debugging purposes
      void fd_check(DRT::Element* ele,                //!< the element we are dealing with
          CORE::LINALG::SerialDenseMatrix& emat,      //!< element matrix to calculate
          CORE::LINALG::SerialDenseVector& erhs,      //!< element rhs to calculate
          CORE::LINALG::SerialDenseVector& subgrdiff  //!< subgrid-diff.-scaling vector
          ) override;

      //! get material parameters
      void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                  //!< density at t_(n)
          std::vector<double>& densnp,                 //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,                 //!< density at t_(n+alpha_M)
          double& visc,                                //!< fluid viscosity
          const int iquad = -1                         //!< id of current gauss point (default = -1)
          ) override = 0;

      /*========================================================================*/
      //! @name stabilization and related functions
      /*========================================================================*/

      //! Calculate quantities used for stabilization
      virtual void prepare_stabilization(
          std::vector<double>& tau,  //!< stabilization parameters (one per transported scalar)
          std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              tauderpot,  //!< derivatives of stabilization parameters w.r.t. electric potential
          const std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_f)
          const double vol                    //!< element volume
      ){};

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! Potential equation ENC
      void CalcMatPotEquENC(CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                              //!< index of current scalar
          const double fac,                                         //!< domain-integration factor
          const double alphaf                                       //!< time factor for ENC
      );

      //! CalcRhs: Potential equation ENC
      void CalcRhsPotEquENC(CORE::LINALG::SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                              //!< index of current scalar
          const double fac,                                         //!< domain-integration factor
          const double conint                                       //!< concentration at GP
      );

      //! process an electrode boundary kinetics point condition
      void calc_elch_boundary_kinetics_point(DRT::Element* ele,  ///< current element
          Teuchos::ParameterList& params,                        ///< parameter list
          DRT::Discretization& discretization,                   ///< discretization
          std::vector<int>& lm,                                  ///< location vector
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,       ///< element matrix
          CORE::LINALG::SerialDenseVector& elevec1_epetra,       ///< element right-hand side vector
          const double
              scalar  ///< scaling factor for element matrix and right-hand side contributions
      );

      //! evaluate an electrode boundary kinetics point condition
      virtual void evaluate_elch_boundary_kinetics_point(
          const DRT::Element* ele,                ///< current element
          CORE::LINALG::SerialDenseMatrix& emat,  ///< element matrix
          CORE::LINALG::SerialDenseVector& erhs,  ///< element right-hand side vector
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< state variables at element nodes
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ehist,       ///< history variables at element nodes
          double timefac,  ///< time factor
          Teuchos::RCP<CORE::Conditions::Condition>
              cond,                       ///< electrode kinetics boundary condition
          const int nume,                 ///< number of transferred electrons
          const std::vector<int> stoich,  ///< stoichiometry of the reaction
          const int kinetics,             ///< desired electrode kinetics model
          const double pot0,              ///< electrode potential on metal side
          const double frt,               ///< factor F/RT
          const double
              scalar  ///< scaling factor for element matrix and right-hand side contributions
      );

      //! evaluate status information on point electrode
      void evaluate_electrode_status_point(const DRT::Element* ele,  ///< current element
          CORE::LINALG::SerialDenseVector& scalars,                  ///< scalars to be integrated
          Teuchos::ParameterList& params,                            ///< parameter list
          Teuchos::RCP<CORE::Conditions::Condition> cond,            ///< condition
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< state variables at element nodes
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephidtnp,                   ///< nodal time derivative vector
          const int kinetics,             ///< desired electrode kinetics model
          const std::vector<int> stoich,  ///< stoichiometry of the reaction
          const int nume,                 ///< number of transferred electrons
          const double pot0,              ///< electrode potential on metal side
          const double frt,               ///< factor F/RT
          const double timefac,           ///< time factor
          const double scalar             ///< scaling factor for current related quantities
      );

      /*========================================================================*/
      //! @name additional service routines
      /*========================================================================*/

      //! validity check with respect to input parameters, degrees of freedom, number of scalars
      //! etc.
      virtual void check_elch_element_parameter(DRT::Element* ele  //!< current element
          ) = 0;

      //! calculate weighted mass flux (no reactive flux so far) -> elch-specific implementation
      virtual void CalculateFlux(CORE::LINALG::Matrix<nsd_, 1>& q,  //!< flux of species k
          const INPAR::SCATRA::FluxType fluxtype,                   //!< type fo flux
          const int k                                               //!< index of current scalar
          ) = 0;

      //! calculate weighted current flux (no reactive flux so far) -> elch-specific implementation
      virtual void CalculateCurrent(CORE::LINALG::Matrix<nsd_, 1>& q,  //!< flux of species k
          const INPAR::SCATRA::FluxType fluxtype,                      //!< type fo flux
          const double fac                                             //!< integration factor
      ){};

      //! calculate error of numerical solution with respect to analytical solution
      void cal_error_compared_to_analyt_solution(const DRT::Element* ele,  //!< element
          Teuchos::ParameterList& params,                                  //!< parameter list
          CORE::LINALG::SerialDenseVector& errors  //!< vector containing L2 and H1 error norms
          ) override;

      //! calculate conductivity of electrolyte solution
      void calculate_conductivity(const DRT::Element* ele,  //!< the element we are dealing with
          const enum INPAR::ELCH::EquPot
              equpot,  //!< type of closing equation for electric potential
          CORE::LINALG::SerialDenseVector&
              sigma,       //!< conductivity of all single ions + overall electrolyte solution
          bool effCond,    //!< flag if effective conductivity should be calculated
          bool specresist  //!< flag if inverse of conductivity should be calculated -> specific
                           //!< resistance
      );

      // Get conductivity from material
      virtual void GetConductivity(const enum INPAR::ELCH::EquPot
                                       equpot,  //!< type of closing equation for electric potential
          double& sigma_all,                    //!< conductivity of electrolyte solution
          std::vector<double>&
              sigma,  //!< conductivity of all single ions + overall electrolyte solution
          bool effCond) = 0;

      //! set internal variables
      void set_internal_variables_for_mat_and_rhs() override = 0;

      //! get elch diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerElch> diff_manager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);
      };

      //! parameter class for electrochemistry problems
      const ScaTraEleParameterElch* elchparams_;

      //! utility class supporting element evaluation
      DRT::ELEMENTS::ScaTraEleUtilsElch<distype>* utils_;
    };


    /// ScaTraEleDiffManagerElch implementation
    /*!
      This class keeps all elch-specific transport parameter needed for the evaluation of an
      element. The ScaTraEleDiffManagerElch is derived from the standard ScaTraEleDiffManager.
    */
    class ScaTraEleDiffManagerElch : public ScaTraEleDiffManager
    {
     public:
      ScaTraEleDiffManagerElch(int numscal) : ScaTraEleDiffManager(numscal), valence_(numscal, 0.0)
      {
      }
      //! Set valence of the single ionic species
      virtual void SetValence(const double valence, const int k) { valence_[k] = valence; };

      //! Access routine for valence of all ionic species or of single ionic species k
      std::vector<double> GetValence() { return valence_; };
      double GetValence(const int k) { return valence_[k]; };

      //! the length of the vector containing the diffusion coefficients is increased by 1
      //! application: ENC with eliminated species
      void increase_length_vector(const int k, const int numscal)
      {
        if (diff_.size() == (unsigned)numscal)
        {
          valence_.push_back(0.0);
          diff_.push_back(0.0);
        }
      };

     protected:
      //! valence of the single ionic species
      std::vector<double> valence_;
    };

    /// ScaTraEleInternalVariableManagerElch implementation
    /*!
      This class manages all internal variables needed for the evaluation of an element.
      The internal variables stored in this class are used by the Nernst-Planck formulation
      as well as the diffusion-conduction formulation.
      All formulation-specific internal variables are stored and managed by a class derived from
      this class (class ScaTraEleInternalVariableManagerElchNP, class
      ScaTraEleInternalVariableManagerElchElectrode, and class
      ScaTraEleInternalVariableManagerElchDiffCond).
    */
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElch : public ScaTraEleInternalVariableManager<NSD, NEN>
    {
      using my = ScaTraEleInternalVariableManager<NSD, NEN>;

     public:
      ScaTraEleInternalVariableManagerElch(
          int numscal, const DRT::ELEMENTS::ScaTraEleParameterElch* elchpara)
          : ScaTraEleInternalVariableManager<NSD, NEN>(numscal),
            parameters_(elchpara),
            frt_(0.),
            // internal variables evaluated at the Gauss point
            conintinv_(numscal),
            gradpot_(true),
            sgconv_(true)
      {
      }
      // compute and set internal variables used by Nernst-Planck formulation and the
      // Diffusion-Conduction formulation

      /*!
       * \brief ompute and set internal variables used by Nernst-Planck formulation and
       * the Diffusion-Conduction formulation
       *
       * @param funct          array for shape functions
       * @param derxy          global derivatives of shape functions w.r.t x,y,z
       * @param ephinp         nodal state variables at t_(n+1) or t_(n+alpha_F)
       * @param ephin          nodal state variables at t_(n)
       * @param econvelnp      nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
       * @param ehist          history vector of transported scalars
       * @param do_setfrt      should FRT be set?
       */
      void set_internal_variables_elch(const CORE::LINALG::Matrix<NEN, 1>& funct,
          const CORE::LINALG::Matrix<NSD, NEN>& derxy,
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>& ephinp,
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>& ephin,
          const CORE::LINALG::Matrix<NSD, NEN>& econvelnp,
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>& ehist, bool do_setfrt = true)
      {
        // call base class (scatra)
        const CORE::LINALG::Matrix<NSD, NEN> eforcevelocity(true);
        my::set_internal_variables(funct, derxy, ephinp, ephin, econvelnp, ehist, eforcevelocity);

        // loop over all transported scalars
        // get concentration of transported scalar k at integration point
        // evaluation of all concentrations is necessary at this point since
        // -> homogeneous reactions of scalar k may depend on all concentrations
        // -> concentration depending material parameters for the diffusion-convection formulation
        // -> avoiding of possible errors (concentration was always defined as a vector where only
        // one
        //    entry was filled)
        for (int k = 0; k < my::numscal_; ++k)
          // calculate 1/concentration at GP at t_(n+1) or t_(n+alpha_F)
          if (my::phinp_[k] > 1e-16) conintinv_[k] = 1 / my::phinp_[k];

        // calculate gradient of electric potential at GP at t_(n+1) or t_(n+alpha_F)
        gradpot_.Multiply(derxy, ephinp[my::numscal_]);

        // set factor F/RT
        if (do_setfrt) SetFRT();
      }

      /*========================================================================*/
      //! @name return methods for internal variables
      /*========================================================================*/

      //! return factor F/RT
      double FRT() const { return frt_; };

      //! return the homogeneous temperature in the scatra field (can be time dependent)
      double Temperature() const { return parameters_->Temperature(); }

      //! return 1/concentration of species k
      const double& ConIntInv(const int k) const { return conintinv_[k]; };

      //! return 1/concentration of all species in a vector
      const std::vector<double>& ConIntInv() const { return conintinv_; };

      //! return gradient of electric potential
      const CORE::LINALG::Matrix<NSD, 1>& GradPot() const { return gradpot_; };

      //! return subgrid velocity
      const CORE::LINALG::Matrix<NEN, 1>& SGConv() const { return sgconv_; };

      //! set factor F/RT
      virtual void SetFRT() { frt_ = parameters_->FRT(); }

      //! return parameter class
      const DRT::ELEMENTS::ScaTraEleParameterElch* ElchParams() const { return parameters_; };

     protected:
      //! parameter class for electrochemistry problems
      const DRT::ELEMENTS::ScaTraEleParameterElch* parameters_;

      /*========================================================================*/
      //! @name constant internal variables
      /*========================================================================*/

      //! pre-calculation of regularly used constant F/RT
      //! (a division is much more expensive than a multiplication)
      double frt_;

      /*========================================================================*/
      //! @name internal variables evaluated at element center or Gauss point
      /*========================================================================*/

      //! 1/concentration at GP
      std::vector<double> conintinv_;
      //! gradient of electric potential
      CORE::LINALG::Matrix<NSD, 1> gradpot_;
      // subgrid velocity
      CORE::LINALG::Matrix<NEN, 1> sgconv_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
