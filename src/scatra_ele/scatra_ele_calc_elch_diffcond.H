/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations

\level 2

*/
/*--------------------------------------------------------------------------*/
#ifndef SCATRA_ELE_CALC_ELCH_DIFFCOND_H
#define SCATRA_ELE_CALC_ELCH_DIFFCOND_H

#include "scatra_ele_parameter_elch_diffcond.H"
#include "scatra_ele_calc_elch_electrode.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElchDiffCond;
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchDiffCond;
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleUtilsElchDiffCond;

    // class implementation
    template <DRT::Element::DiscretizationType distype,
        int probdim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim>
    class ScaTraEleCalcElchDiffCond : public ScaTraEleCalcElchElectrode<distype, probdim>
    {
     protected:
      /// protected constructor, since we are a Singleton, but a derived class exists
      ScaTraEleCalcElchDiffCond(
          const int numdofpernode, const int numscal, const std::string& disname);

      using my = ScaTraEleCalc<distype, probdim>;
      using myelch = ScaTraEleCalcElch<distype, probdim>;
      using myelectrode = ScaTraEleCalcElchElectrode<distype, probdim>;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      /// Singleton access method
      static ScaTraEleCalcElchDiffCond<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      /*========================================================================*/
      //! @name general framework
      /*========================================================================*/

      //! calculate contributions to matrix and rhs (inside of loop over all scalars)
      virtual void CalcMatAndRhs(Epetra_SerialDenseMatrix& emat,  //!< element matrix to calculate
          Epetra_SerialDenseVector& erhs,                         //!< element rhs to calculate
          const int k,                                            //!< index of current scalar
          const double fac,                                       //!< domain-integration factor
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
      );

      //! calculate contributions to matrix and rhs (outside of loop over all scalars)
      virtual void CalcMatAndRhsOutsideScalarLoop(
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
      //! @name material and related and related functions
      /*========================================================================*/

      //! get material parameters
      virtual void GetMaterialParams(const DRT::Element* ele,  //!< current element
          std::vector<double>& densn,                          //!< density at t_(n)
          std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,  //!< density at t_(n+alpha_M)
          double& visc,                 //!< fluid viscosity
          const int iquad = -1          //!< ID of current integration point
      );

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! CalcMat: Conduction term with inserted current - ohmic overpotential
      virtual void CalcMatCondOhm(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                             //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double invfval,     //!< 1/(F z_k)
          const CORE::LINALG::Matrix<nsd_, 1>& gradpot  //!< gradient of potential at GP
      );

      //! CalcMat: Conduction term with inserted current - conc. overpotential
      void CalcMatCondConc(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                      //!< index of current scalar
          const double timefacfac,      //!< domain-integration factor times time-integration factor
          const double rtffcval,        //!< RT/F^2/Newman_const_c/z_k
          const double newman_const_a,  //!< Newman constant a
          const double newman_const_b,  //!< Newman constant b
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi,  //!< gradient of concentration at GP
          const std::vector<double>& conintinv           //!< inverted concentration at GP
      );

      //! CalcMat: Conduction term without inserted current
      virtual void CalcMatCond(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                          //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double invfval,     //!< 1/(F z_k)
          const CORE::LINALG::Matrix<nsd_, 1>& curint  //!< current at GP
      );

      //! CalcMat: Additional diffusion term without inserted current
      void CalcMatCondDiff(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                      //!< index of current scalar
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double invfval,     //!< 1/(F z_k)
          const std::vector<CORE::LINALG::Matrix<nsd_, 1>>&
              gradphi  //!< gradient of concentration at GP
      );

      //! Potential equation div i inserted current - concentration overpotential
      void CalcMatPotEquDiviConc(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const int k,                                            //!< index of current scalar
          const double timefacfac,      //!< domain-integration factor times time-integration factor
          const double rtffc,           //!< RT/(F^2 Newman_const_c)
          const double rtf,             //!< RT/F
          const double invf,            //!< 1/F
          const double newman_const_a,  //!< Newman constant a
          const double newman_const_b,  //!< Newman constant b
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi,  //!< gradient of concentration at GP
          const double conintinv                         //!< inverted concentration at GP
      );

      //! Potential equation div i without inserted current
      void CalcMatPotEquDivi(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double invf         //!< 1/F
      );

      //! CalcMat: Current equation current
      virtual void CalcMatCurEquCur(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double invf         //!< 1/F
      );

      //! CalcMat: Current equation ohmic overpotential
      virtual void CalcMatCurEquOhm(
          Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double invf,        //!< 1/F
          const CORE::LINALG::Matrix<nsd_, 1>& gradpot  //!< gradient of potenial at GP
      );

      //! CalcMat: Current equation concentration overpotential
      void CalcMatCurEquConc(Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
          const double timefacfac,  //!< domain-integration factor times time-integration factor
          const double rtf,         //!< RT/F
          const double rtffc,       //!< RT/(F^2 Newman_const_c)
          const std::vector<double>& invfval,  //!< 1/(F z_k)
          const double newman_const_a,         //!< Newman constant a
          const double newman_const_b,         //!< Newman constant b
          const std::vector<CORE::LINALG::Matrix<nsd_, 1>>&
              gradphi,                          //!< gradient of concentration at GP
          const std::vector<double>& conintinv  //!< inverted concentration at GP
      );

      //! CalcRhs: Conduction term with inserted current - ohmic overpotential
      virtual void CalcRhsCondOhm(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                             //!< index of current scalar
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double invfval,                         //!< 1/(F z_k)
          const CORE::LINALG::Matrix<nsd_, 1>& gradpot  //!< gradient of potenial at GP
      );

      //! CalcRhs: Conduction term with inserted current - conc. overpotential
      void CalcRhsCondConc(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                      //!< index of current scalar
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double rtffcval,                         //!< RT/(F^2 Newman_const_c z_k)
          const double newman_const_a,                   //!< Newman constant a
          const double newman_const_b,                   //!< Newman constant b
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi,  //!< gradient of concentration at GP
          const std::vector<double>& conintinv           //!< inverted concentration at GP
      );

      //! CalcRhs: Conduction term without inserted current
      virtual void CalcRhsCond(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                          //!< index of current scalar
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double invfval,                        //!< 1/(F z_k)
          const CORE::LINALG::Matrix<nsd_, 1>& curint  //!< current at GP
      );

      //! CalcRhs: Additional diffusion term without inserted current
      void CalcRhsCondDiff(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                      //!< index of current scalar
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const std::vector<CORE::LINALG::Matrix<nsd_, 1>>&
              gradphi  //!< gradient of concentration at GP
      );

      //! CalcRhs: Potential equation div i inserted current - conc. overpotential
      void CalcRhsPotEquDiviConc(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                                            //!< index of current scalar
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double rtf,     //!< RT/F
          const std::vector<double>& invfval,            //!< 1/(F z_k)
          const double rtffc,                            //!< RT/(F^2 Newman_const_c)
          const double newman_const_a,                   //!< Newman constant a
          const double newman_const_b,                   //!< Newman constant b
          const CORE::LINALG::Matrix<nsd_, 1>& gradphi,  //!< gradient of concentration at GP
          const double conintinv                         //!< inverted concentration at GP
      );

      //! CalcRhs: Potential equation divi without inserted current
      void CalcRhsPotEquDivi(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double invf,    //!< 1/F
          const CORE::LINALG::Matrix<nsd_, 1>& curint  //!< current at GP
      );

      //! CalcRhs: Current equation - current
      virtual void CalcRhsCurEquCur(
          Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double invf,    //!< 1/F
          const CORE::LINALG::Matrix<nsd_, 1>& curint  //!< current at GP
      );

      //! CalcRhs: Current equation - ohmic overpotential
      virtual void CalcRhsCurEquOhm(
          Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double invf,    //!< 1/F
          const CORE::LINALG::Matrix<nsd_, 1>& gradpot  //!< gradient of potenial at GP
      );

      //! Current equation - concentration overpotential
      void CalcRhsCurEquConc(Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
          const double rtf,     //!< RT/F
          const std::vector<double>& invfval,  //!< 1/(F z_k)
          const double rtffc,                  //!< RT/(F^2 Newman_const_c)
          const double newman_const_a,         //!< Newman constant a
          const double newman_const_b,         //!< Newman constant b
          const std::vector<CORE::LINALG::Matrix<nsd_, 1>>&
              gradphi,                          //!< vector of gradient of concentration at GP
          const std::vector<double>& conintinv  //!< inverted concentration at GP
      );


      /*========================================================================*/
      //! @name additional service routines
      /*========================================================================*/

      //! evaluate action
      virtual int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

      //! evaluate an electrode boundary kinetics point condition
      void EvaluateElchBoundaryKineticsPoint(const DRT::Element* ele,  ///< current element
          Epetra_SerialDenseMatrix& emat,                              ///< element matrix
          Epetra_SerialDenseVector& erhs,  ///< element right-hand side vector
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< state variables at element nodes
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
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

      //! evaluate electrode kinetics domain condition
      void EvaluateElchDomainKinetics(const DRT::Element* ele,  ///< the actual boundary element
          Epetra_SerialDenseMatrix& emat,                       ///< element-matrix
          Epetra_SerialDenseVector& erhs,                       ///< element-rhs
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ehist,  ///< nodal history vector
          double timefac,                                           ///< time factor
          Teuchos::RCP<DRT::Condition> cond,                        ///< the condition
          const int nume,                 ///< number of transferred electrons
          const std::vector<int> stoich,  ///< stoichiometry of the reaction
          const int kinetics,             ///< desired electrode kinetics model
          const double pot0               ///< actual electrode potential on metal side
      );

      //! validity check for all elements with respect to
      //! formulation-specific parameter, degree's of freedom, number of scalars, ...
      void CheckElchElementParameter(DRT::Element* ele  //!< the element we are dealing with
      );

      //! calculate element mass matrix and element residual for initial time derivative
      void CalcInitialTimeDerivative(DRT::Element* ele,  //!< current element
          Epetra_SerialDenseMatrix& emat,                //!< element matrix
          Epetra_SerialDenseVector& erhs,                //!< element residual
          Teuchos::ParameterList& params,                //!< parameter list
          DRT::Discretization& discretization,           //!< discretization
          DRT::Element::LocationArray& la                //!< location array
      );

      //! Correct RHS calculated from CalcRHSLinMass() for the linearized mass term
      void CorrectRHSFromCalcRHSLinMass(
          Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
          const int k,                     //!< index of current scalar
          const double fac,                //!< domain-integration factor
          const double densnp,             //!< density at time_(n+1)
          const double phinp               //!< scalar at time_(n+1)
      );

      //!  calculate weighted mass flux (no reactive flux so far) -> elch-specific implementation
      void CalculateFlux(CORE::LINALG::Matrix<nsd_, 1>& q,  //!< flux of species k
          const INPAR::SCATRA::FluxType fluxtype,           //!< type fo flux
          const int k                                       //!< index of current scalar
          ) override;

      //!  calculate weighted current flux (no reactive flux so far) -> elch-specific implementation
      void CalculateCurrent(CORE::LINALG::Matrix<nsd_, 1>& q,  //!< flux of species k
          const INPAR::SCATRA::FluxType fluxtype,              //!< type fo flux
          const double fac                                     //!< integration factor
          ) override;

      //! calculate error of numerical solution with respect to analytical solution
      void CalErrorComparedToAnalytSolution(
          const DRT::Element* ele,          //!< the element we are dealing with
          Teuchos::ParameterList& params,   //!< parameter list
          Epetra_SerialDenseVector& errors  //!< vector containing L2-error norm
      );

      void CalcElchDomainKinetics(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseVector& elevec1_epetra);

      void EvaluateElectrodeStatus(const DRT::Element* ele,  ///< the actual boundary element
          Epetra_SerialDenseVector& scalars,                 ///< scalars to be computed
          Teuchos::ParameterList& params,                    ///< the parameter list
          Teuchos::RCP<DRT::Condition> cond,                 ///< the condition
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephinp,  ///< nodal values of concentration and electric potential
          const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
              ephidtnp,                   ///< nodal time derivative vector
          const int kinetics,             ///< desired electrode kinetics model
          const std::vector<int> stoich,  ///< stoichiometry of the reaction
          const int nume,                 ///<  number of transferred electrons
          const double pot0,              ///< actual electrode potential on metal side at t_{n+1}
          const double timefac            ///< factor due to time discretization
      );

      //! set internal variables for diffusion-conduction formulation
      void SetInternalVariablesForMatAndRHS();

      //! get diffusion manager for diffusion-conduction formulation
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> DiffManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerElchDiffCond>(my::diffmanager_);
      }

      //! get internal variable manager for diffusion-conduction formulation
      Teuchos::RCP<ScaTraEleInternalVariableManagerElchDiffCond<nsd_, nen_>> VarManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerElchDiffCond<nsd_, nen_>>(
            my::scatravarmanager_);
      }

      //! get utility class supporting element evaluation for diffusion-conduction formulation
      ScaTraEleUtilsElchDiffCond<distype>* Utils()
      {
        return static_cast<ScaTraEleUtilsElchDiffCond<distype>*>(myelch::utils_);
      }

      void CalculateMeanElectrodeConcentration(const DRT::Element* const& ele,
          const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          Epetra_SerialDenseVector& conc) override;

      //! flag for used element formulation (material based)
      INPAR::ELCH::DiffCondMat diffcondmat_;

      //! parameter class for diffusion-conduction formulation
      const ScaTraEleParameterElchDiffCond* diffcondparams_;
    };


    /// ScaTraEleDiffManagerElchDiffCond implementation
    /*!
      This class keeps all Diffusion-Conduction-specific transport parameter needed for the
      evaluation of an element.
    */
    class ScaTraEleDiffManagerElchDiffCond : public ScaTraEleDiffManagerElchElectrode
    {
     public:
      using dmelch = ScaTraEleDiffManagerElch;
      using dmelectrode = ScaTraEleDiffManagerElchElectrode;

      ScaTraEleDiffManagerElchDiffCond(int numscal)
          : ScaTraEleDiffManagerElchElectrode(numscal),
            transnum_(numscal, 0.0),
            derivtransnum_(numscal, std::vector<double>(numscal, 0.0)),
            thermfac_(1.0),
            derivthermfac_(numscal, 0.0),
            invval_(numscal, 0.),
            invfval_(numscal, 0.),
            // In the moment, we have only one phase but the framework is flexible
            eps_(1, 1.0),
            tort_(1, 1.0)
      {
      }

      /*========================================================================*/
      //! @name access methods for transport parameter
      /*========================================================================*/

      //! set valence and related parameters for single ionic species
      void SetValence(const double valence, const int k) override
      {
        // call base class routine
        dmelch::SetValence(valence, k);

        const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();

        // set additional parameters involving inverse of valence
        invval_[k] = 1. / valence;
        invfval_[k] = 1. / faraday / valence;
      };

      //! Set transference numbers with respect to single ionic species
      virtual void SetTransNum(const double transnum, const int k) { transnum_[k] = transnum; };

      //! Access routine for transference numbers with respect to single ionic species
      double GetTransNum(const int k) { return transnum_[k]; };

      //! Set derivative of transference numbers with respect to concentrations
      virtual void SetDerivTransNum(const double derivtransnum, const int k, const int iscal)
      {
        (derivtransnum_[k])[iscal] = derivtransnum;
      };

      //! Access routine for derivative of transference numbers with respect to concentrations
      double GetDerivTransNum(const int k, const int iscal) { return (derivtransnum_[k])[iscal]; };

      //! Set thermodynamic factor for a specific electrolyte solution
      void SetThermFac(const double thermfac) { thermfac_ = thermfac; };

      //! Access routine for the thermodynamic factor for a specific electrolyte solution
      double GetThermFac() { return thermfac_; };

      //! Set derivative of thermodynamic factor with respect to concentrations
      void SetDerivThermFac(const double derivthermfac, const int k)
      {
        derivthermfac_[k] = derivthermfac;
      };

      //! Access routine for derivative of thermodynamic factor with respect to concentrations
      double GetDerivThermFac(const int k) { return derivthermfac_[k]; };

      //! Calculate conductivity based on diffusion coefficients
      void CalcConductivity(const int numscal, const double ffrt, const std::vector<double>& conint)
      {
        double cond = 0.0;
        for (int ispec = 0; ispec < numscal; ++ispec)
        {
          // conductivity
          cond += ffrt * valence_[ispec] * valence_[ispec] * diff_[ispec] * conint[ispec];
          // derivation of conductivity wrt concentrations
          concderivcond_[ispec] = ffrt * valence_[ispec] * valence_[ispec] * diff_[ispec];
        }

        cond_ = cond;
      };

      //! Calculate transference numbers based on diffusion coefficients
      void CalcTransNum(const int numscal, const std::vector<double>& conint)
      {
        // conductivity without ffrt
        double sum = 0.0;
        for (int k = 0; k < numscal; ++k)
        {
          sum += valence_[k] * valence_[k] * diff_[k] * conint[k];
        }
        double denomin = sum * sum;

        for (int k = 0; k < numscal; ++k)
        {
          transnum_[k] = valence_[k] * valence_[k] * diff_[k] * conint[k] / sum;

          for (int iscal = 0; iscal < numscal; ++iscal)
          {
            if (k == iscal)
            {
              (derivtransnum_[k])[iscal] =
                  (valence_[k] * valence_[k] * diff_[k] *
                      (sum - valence_[k] * valence_[k] * diff_[k] * conint[k])) /
                  denomin;
            }
            else
            {
              (derivtransnum_[k])[iscal] = (-valence_[k] * valence_[k] * diff_[k] * conint[k] *
                                               valence_[iscal] * valence_[iscal] * diff_[iscal]) /
                                           denomin;
            }
          }
        }
      };

      double InvFVal(const int k) const { return invfval_[k]; };
      const std::vector<double>& InvFVal() const { return invfval_; };

      /*========================================================================*/
      //! @name access methods for geometrical parameter
      /*========================================================================*/

      //! Set transference numbers with respect to single ionic species
      virtual void SetPhasePoro(const double eps, const int phase) { eps_[phase] = eps; }

      double GetPhasePoro(const int phase) const { return eps_[phase]; };

      //! Set transference numbers with respect to single ionic species
      virtual void SetPhaseTort(const double tort, const int phase) { tort_[phase] = tort; }

      // get geometrical parameter: porosity*tortuosity
      double GetPhasePoroTort(const int phase) const
      {
        double epstort = eps_[phase] * tort_[phase];
        return epstort;
      };

      /*========================================================================*/
      //! @name output
      /*========================================================================*/

      //! Output of transport parameter (to screen)
      void OutputTransportParams(const int numscal) override
      {
        // call base class routine
        dmelectrode::OutputTransportParams(numscal);

        // additional outputs
        for (int k = 0; k < numscal; ++k)
          std::cout << "valence " << k << ":   " << valence_[k] << std::endl;

        for (int k = 0; k < numscal; ++k)
          std::cout << "transference number " << k << ":   " << transnum_[k] << std::endl;

        for (int k = 0; k < numscal; ++k)
        {
          for (int iscal = 0; iscal < numscal; ++iscal)
          {
            std::cout << "derivation transference number (" << k << "," << iscal
                      << "):  " << (derivtransnum_[k])[iscal] << std::endl;
          }
        }
        std::cout << std::endl;

        std::cout << "thermodynamic factor:   " << thermfac_ << std::endl;

        for (int k = 0; k < numscal; ++k)
          std::cout << "derivation of thermodynamic factor " << k << ":   " << derivthermfac_[k]
                    << std::endl;
        std::cout << std::endl;

        std::cout << "porosity species" << 0 << ":   " << eps_[0] << std::endl;
        std::cout << std::endl;

        std::cout << "tortuosity species" << 0 << ":   " << tort_[0] << std::endl;
        std::cout << std::endl;
      };

     protected:
      /*========================================================================*/
      //! @name transport parameter
      /*========================================================================*/

      //! transference numbers for single ionic species
      std::vector<double> transnum_;

      //! derivative of transference numbers with respect to concentrations
      std::vector<std::vector<double>> derivtransnum_;

      //! thermodynamic factor for a specific electrolyte solution
      //! transport parameter used only in the Newman model
      double thermfac_;

      //! derivative of thermodynamic factor with respect to concentrations
      std::vector<double> derivthermfac_;

      //! constant parameters 1/(z_k)
      std::vector<double> invval_;

      //! constant parameters 1/(F z_k)
      std::vector<double> invfval_;

      /*========================================================================*/
      //! @name geometrical parameters of the porous medium
      /*========================================================================*/

      //! porosity
      std::vector<double> eps_;
      //! tortuosity
      std::vector<double> tort_;
    };


    /// ScaTraEleInternalVariableManagerElchDiffCond implementation
    /*!
      This class keeps all internal variables needed for the diffusion-conduction formulation.
    */
    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerElchDiffCond
        : public ScaTraEleInternalVariableManagerElchElectrode<NSD, NEN>
    {
     public:
      using vm = ScaTraEleInternalVariableManager<NSD, NEN>;
      using vmelch = ScaTraEleInternalVariableManagerElch<NSD, NEN>;
      using vmelectrode = ScaTraEleInternalVariableManagerElchElectrode<NSD, NEN>;

      ScaTraEleInternalVariableManagerElchDiffCond(int numscal,
          const DRT::ELEMENTS::ScaTraEleParameterElch* elchparams,
          const DRT::ELEMENTS::ScaTraEleParameterElchDiffCond* diffcondparams)
          : ScaTraEleInternalVariableManagerElchElectrode<NSD, NEN>(numscal, elchparams),
            diffcondparams_(diffcondparams),
            rtf_(0.),
            rtffc_(0.),
            curint_(true)
      {
      }

      //! compute and set internal variables only used by the Diffusion-Conduction formulation
      void SetInternalVariablesElchDiffCond(
          const CORE::LINALG::Matrix<NEN, 1>& funct,  //!< array for shape functions
          const CORE::LINALG::Matrix<NSD, NEN>&
              derxy,  //!< global derivatives of shape functions w.r.t x,y,z
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ephinp,  //!< nodal state variables at t_(n+1) or t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ephin,  //!< nodal state variables at t_(n)
          const CORE::LINALG::Matrix<NSD, NEN>&
              econvelnp,  //!< nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
          const std::vector<CORE::LINALG::Matrix<NEN, 1>>&
              ehist  //!< history vector of transported scalars
      )
      {
        // set internal variables in base variable manager
        vmelectrode::SetInternalVariablesElchElectrode(
            funct, derxy, ephinp, ephin, econvelnp, ehist);

        // set constant parameter RT/F
        rtf_ = 1.0 / vmelch::frt_;

        // set constant parameter RT/(F^2*Newman_const_C)
        rtffc_ = rtf_ * vmelectrode::invf_ / diffcondparams_->NewmanConstC();

        if (diffcondparams_->CurSolVar())
          for (unsigned idim = 0; idim < NSD; ++idim)
            curint_(idim, 0) = ephinp[vm::numscal_ + 1 + idim].Dot(funct);
      };

      /*========================================================================*/
      //! @name return constant internal variables
      /*========================================================================*/

      double RTF() const { return rtf_; };
      double RTFFC() const { return rtffc_; };

      /*========================================================================*/
      //! @name return methods for GP values
      /*========================================================================*/

      //! return current density at GP
      const CORE::LINALG::Matrix<NSD, 1>& CurInt() { return curint_; };

     protected:
      //! parameter class for diffusion-conduction formulation
      const DRT::ELEMENTS::ScaTraEleParameterElchDiffCond* diffcondparams_;

      /*========================================================================*/
      //! @name constant internal variables
      /*========================================================================*/

      //! constant parameter RT/F
      double rtf_;
      //! constant parameter RT/(F^2*Newman_const_C)
      //! attention: this is a newman specific parameter
      double rtffc_;

      /*========================================================================*/
      //! @name internal variables evaluated at the Gauss point
      /*========================================================================*/

      //! current density at Gauss point
      CORE::LINALG::Matrix<NSD, 1> curint_;
    };  // class ScaTraEleInternalVariableManagerElchDiffCond
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
