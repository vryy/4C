/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch.H

\brief evaluation of scatra boundary terms at integration points

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_H

#include "scatra_ele_boundary_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleParameterElch;

    // class implementation
    template<DRT::Element::DiscretizationType distype> class ScaTraEleBoundaryCalcElch : public ScaTraEleBoundaryCalc<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;

    public:

      //! singleton access method
      // not needed, since class is purely virtual

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElch(){return;};

      //! singleton destruction
      virtual void Done() = 0;

    protected:

      //! protected constructor for singletons
      ScaTraEleBoundaryCalcElch(const int numdofpernode,const int numscal);

      //! evaluate action
      int EvaluateAction(
          DRT::ELEMENTS::TransportBoundary*   ele,
          Teuchos::ParameterList&             params,
          DRT::Discretization&                discretization,
          SCATRA::BoundaryAction              action,
          std::vector<int>&                   lm,
          Epetra_SerialDenseMatrix&           elemat1_epetra,
          Epetra_SerialDenseMatrix&           elemat2_epetra,
          Epetra_SerialDenseVector&           elevec1_epetra,
          Epetra_SerialDenseVector&           elevec2_epetra,
          Epetra_SerialDenseVector&           elevec3_epetra
          );

      //! evaluate electrode kinetics boundary condition
      virtual void EvaluateElchBoundaryKinetics(
         const DRT::Element*                 ele,        ///< the actual boundary element
         Epetra_SerialDenseMatrix&           emat,       ///< element-matrix
         Epetra_SerialDenseVector&           erhs,       ///< element-rhs
         const std::vector<double>&          ephinp,     ///< actual conc. and pot. values
         const std::vector<double>&          ehist,      ///< element history vector
         double                              timefac,    ///< time factor
         Teuchos::RCP<const MAT::Material>   material,   ///< the material
         Teuchos::RCP<DRT::Condition>        cond,       ///< the condition
         const int                           nume,       ///< number of transferred electrons
         const std::vector<int>              stoich,     ///< stoichiometry of the reaction
         const int                           kinetics,   ///< desired electrode kinetics model
         const double                        pot0,       ///< actual electrode potential on metal side
         const double                        frt,        ///< factor F/RT
         const double                        scalar      ///< scaling factor for element matrix and residual contributions
      );

      //! evaluate electrode kinetics boundary condition
      virtual void CalcElchBoundaryKinetics(
          DRT::ELEMENTS::TransportBoundary* ele,
          Teuchos::ParameterList&           params,
          DRT::Discretization&              discretization,
          std::vector<int>&                 lm,
          Epetra_SerialDenseMatrix&         elemat1_epetra,
          Epetra_SerialDenseVector&         elevec1_epetra,
          const double                      scalar
      );

      //! provide electrode kinetics status information
      void ElectrodeStatus(
          const DRT::Element*           ele,        ///< the actual boundary element
          Epetra_SerialDenseVector&     scalars,    ///< scalars to be computed
          Teuchos::ParameterList&       params,     ///< the parameter list
          Teuchos::RCP<DRT::Condition>  cond,       ///< the condition
          const std::vector<double>&    ephinp,     ///< current conc. and potential values
          const std::vector<double>&    ephidtnp,   ///< time derivative vector evaluated at t_{n+1}
          const int                     kinetics,   ///< desired electrode kinetics model
          const std::vector<int>        stoich,     ///< stoichiometry of the reaction
          const int                     nume,       ///<  number of transferred electrons
          const double                  pot0,       ///< actual electrode potential on metal side at t_{n+1}
          const double                  frt,        ///< factor F/RT
          const double                  timefac,    ///< factor due to time discretization
          const double                  scalar      ///< scaling factor for current related quantities
      );

      //! evaluate linearization of nernst equation
      void CalcNernstLinearization(
          DRT::ELEMENTS::TransportBoundary* ele,
          Teuchos::ParameterList&           params,
          DRT::Discretization&              discretization,
          std::vector<int>&                 lm,
          Epetra_SerialDenseMatrix&         elemat1_epetra,
          Epetra_SerialDenseVector&         elevec1_epetra
      );

      //! calculate cell voltage
      void CalcCellVoltage(
          const DRT::Element*               ele,              //!< the element we are dealing with
          Teuchos::ParameterList&           params,           //!< parameter list
          DRT::Discretization&              discretization,   //!< discretization
          const std::vector<int>&           lm,               //!< location vector
          Epetra_SerialDenseVector&         scalars           //!< result vector for scalar integrals to be computed
      );

      //! extract valence of species k from element material
      virtual const double GetValence(
          const Teuchos::RCP<const MAT::Material>&   material,   //! element material
          const int                                  k           //! species number
      ) const = 0;

      //! return parameter class for electrochemistry problems
      const DRT::ELEMENTS::ScaTraEleParameterElch* ElchParams() const {return static_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatraparams_);};
    }; // class ScaTraEleBoundaryCalcElch
  } // namespace ELEMENTS
} // namespace DRT
#endif
