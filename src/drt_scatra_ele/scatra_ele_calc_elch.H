/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch.H

\brief evaluation of ScaTra elements for ion-transport equation

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_ELCH_H
#define SCATRA_ELE_CALC_ELCH_H

#include "scatra_ele_calc.H"
#include "scatra_ele_parameter_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElch;
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElch;
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleUtilsElch;

    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcElch : public ScaTraEleCalc<distype>
    {
     protected:
      /// protected constructor, since we are a singleton
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalcElch(const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;

     public:
      virtual ~ScaTraEleCalcElch() {}

      /// called upon destruction
      virtual void Done() = 0;

      /// In this class we do not define a static ScaTraEle...* Instance
      /// since only derived child classes are free to be allocated!!

      //! evaluate the element
      int Evaluate(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra);

      //! evaluate action
      virtual int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

     protected:
      /*========================================================================*/
      //! @name general framework
      /*========================================================================*/

      //! Prepare everything what is needed in CallMatAndRhs() to calculate the sysmat and the rhs
      virtual void Sysmat(DRT::Element* ele,   //!< the element we are dealing with
          Epetra_SerialDenseMatrix& emat,      //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,      //!< element rhs to calculate
          Epetra_SerialDenseVector& subgrdiff  //!< subgrid-diff.-scaling vector
      );

      //! calculate contributions to matrix and rhs (inside of loop over all scalars)
      virtual void CalcMatAndRhs(Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,                         //!< element rhs to calculate+
          const int k,                                            //!< index of current scalar
          const double fac,                                       //!< domain-integration factor
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
          ) = 0;

      //! calculate contributions to matrix and rhs (outside of loop over all scalars)
      virtual void CalcMatAndRhsOutsideScalarLoop(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,  //!< element rhs to calculate
          const double fac,                //!< domain-integration factor
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double rhsfac  //!< time-integration factor for rhs times domain-integration factor
          ) = 0;

      //! Correction for additional flux terms / currents across Dirichlet boundaries
      virtual void CorrectionForFluxAcrossDC(
          DRT::Discretization& discretization,  //!< discretization
          const std::vector<int>& lm,           //!< location vector
          Epetra_SerialDenseMatrix& emat,       //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs        //!< element rhs to calculate
          ) = 0;

      //! finite difference check for debugging purposes
      void FDCheck(DRT::Element* ele,          //!< the element we are dealing with
          Epetra_SerialDenseMatrix& emat,      //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,      //!< element rhs to calculate
          Epetra_SerialDenseVector& subgrdiff  //!< subgrid-diff.-scaling vector
      );

      //! get material parameters
      virtual void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                          //!< density at t_(n)
          std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,  //!< density at t_(n+alpha_M)
          double& visc,                 //!< fluid viscosity
          const int iquad = -1          //!< id of current gauss point (default = -1)
          ) = 0;

      /*========================================================================*/
      //! @name stabilization and related functions
      /*========================================================================*/

      //! Calculate quantities used for stabilization
      virtual void PrepareStabilization(
          std::vector<double>& tau,  //!< stabilization parameters (one per transported scalar)
          std::vector<LINALG::Matrix<my::nen_, 1>>&
              tauderpot,  //!< derivatives of stabilization parameters w.r.t. electric potential
          const std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_f)
          const double vol                    //!< element volume
      )
      {
        return;
      };

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! Potential equation ENC
      void CalcMatPotEquENC(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                       //!< index of current scalar
          const double fac,                                  //!< domain-integration factor
          const double alphaf                                //!< time factor for ENC
      );

      //! CalcRhs: Potential equation ENC
      void CalcRhsPotEquENC(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                       //!< index of current scalar
          const double fac,                                  //!< domain-integration factor
          const double conint                                //!< concentration at GP
      );

      //! process an electrode boundary kinetics point condition
      void CalcElchBoundaryKineticsPoint(DRT::Element* ele,  ///< current element
          Teuchos::ParameterList& params,                    ///< parameter list
          DRT::Discretization& discretization,               ///< discretization
          std::vector<int>& lm,                              ///< location vector
          Epetra_SerialDenseMatrix& elemat1_epetra,          ///< element matrix
          Epetra_SerialDenseVector& elevec1_epetra,          ///< element right-hand side vector
          const double
              scalar  ///< scaling factor for element matrix and right-hand side contributions
      );

      //! evaluate an electrode boundary kinetics point condition
      virtual void EvaluateElchBoundaryKineticsPoint(const DRT::Element* ele,  ///< current element
          Epetra_SerialDenseMatrix& emat,                                      ///< element matrix
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

      //! evaluate status information on point electrode
      void EvaluateElectrodeStatusPoint(const DRT::Element* ele,  ///< current element
          Epetra_SerialDenseVector& scalars,                      ///< scalars to be integrated
          Teuchos::ParameterList& params,                         ///< parameter list
          Teuchos::RCP<DRT::Condition> cond,                      ///< condition
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
              ephinp,  ///< state variables at element nodes
          const std::vector<LINALG::Matrix<my::nen_, 1>>&
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
      virtual void CheckElchElementParameter(DRT::Element* ele  //!< current element
          ) = 0;

      //! calculate weighted mass flux (no reactive flux so far) -> elch-specific implementation
      virtual void CalculateFlux(LINALG::Matrix<my::nsd_, 1>& q,  //!< flux of species k
          const INPAR::SCATRA::FluxType fluxtype,                 //!< type fo flux
          const int k                                             //!< index of current scalar
          ) = 0;

      //! calculate weighted current flux (no reactive flux so far) -> elch-specific implementation
      virtual void CalculateCurrent(LINALG::Matrix<my::nsd_, 1>& q,  //!< flux of species k
          const INPAR::SCATRA::FluxType fluxtype,                    //!< type fo flux
          const double fac                                           //!< integration factor
      )
      {
        return;
      };

      //! calculate error of numerical solution with respect to analytical solution
      virtual void CalErrorComparedToAnalytSolution(const DRT::Element* ele,  //!< element
          Teuchos::ParameterList& params,                                     //!< parameter list
          Epetra_SerialDenseVector& errors  //!< vector containing L2 and H1 error norms
      );

      //! calculate conductivity of electrolyte solution
      void CalculateConductivity(const DRT::Element* ele,  //!< the element we are dealing with
          const enum INPAR::ELCH::EquPot
              equpot,  //!< type of closing equation for electric potential
          Epetra_SerialDenseVector&
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
      virtual void SetInternalVariablesForMatAndRHS() = 0;

      //! get elch diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerElch> DiffManager()
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
        return;
      }

      ~ScaTraEleDiffManagerElch(){};

      //! Set valence of the single ionic species
      virtual void SetValence(const double valence, const int k)
      {
        valence_[k] = valence;
        return;
      };

      //! Access routine for valence of all ionic species or of single ionic species k
      const std::vector<double> GetValence() { return valence_; };
      double GetValence(const int k) { return valence_[k]; };

      //! the length of the vector containing the diffusion coefficients is increased by 1
      //! application: ENC with eliminated species
      void IncreaseLengthVector(const int k, const int numscal)
      {
        if (diff_.size() == (unsigned)numscal)
        {
          valence_.push_back(0.0);
          diff_.push_back(0.0);
        }

        return;
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
      typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

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
        return;
      }

      virtual ~ScaTraEleInternalVariableManagerElch(){};

      // compute and set internal variables used by Nernst-Planck formulation and the
      // Diffusion-Conduction formulation
      void SetInternalVariablesElch(
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
        // call base class (scatra)
        my::SetInternalVariables(funct, derxy, ephinp, ephin, econvelnp, ehist);

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
          conintinv_[k] = 1 / my::phinp_[k];

        // calculate gradient of electric potential at GP at t_(n+1) or t_(n+alpha_F)
        gradpot_.Multiply(derxy, ephinp[my::numscal_]);

        // set factor F/RT
        SetFRT();

        return;
      };

      /*========================================================================*/
      //! @name return methods for internal variables
      /*========================================================================*/

      //! return factor F/RT
      double FRT() const { return frt_; };

      //! return 1/concentration of species k
      const double& ConIntInv(const int k) const { return conintinv_[k]; };

      //! return 1/concentration of all species in a vector
      const std::vector<double>& ConIntInv() const { return conintinv_; };

      //! return gradient of electric potential
      const LINALG::Matrix<NSD, 1>& GradPot() const { return gradpot_; };

      //! return subgrid velocity
      const LINALG::Matrix<NEN, 1>& SGConv() const { return sgconv_; };

      //! set factor F/RT
      virtual void SetFRT()
      {
        frt_ = parameters_->FRT();

        return;
      }

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
      LINALG::Matrix<NSD, 1> gradpot_;
      // subgrid velocity
      LINALG::Matrix<NEN, 1> sgconv_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT
#endif
