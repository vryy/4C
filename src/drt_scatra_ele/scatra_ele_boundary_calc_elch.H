/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch.H

\brief evaluation of scatra boundary terms at integration points

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
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

template<DRT::Element::DiscretizationType distype>
class ScaTraEleBoundaryCalcElch : public ScaTraBoundaryImpl<distype>
{
  typedef DRT::ELEMENTS::ScaTraBoundaryImpl<distype> my;

public:

  /// Singleton access method
  static ScaTraEleBoundaryCalcElch<distype> * Instance(const int numdofpernode, const int numscal, bool create=true );

  /// Empty destructor
  virtual ~ScaTraEleBoundaryCalcElch() {};

  /// called upon destruction
  void Done();

private:
  /// private constructor since we are singleton
  ScaTraEleBoundaryCalcElch(const int numdofpernode, const int numscal);

  int EvaluateAction(
      DRT::ELEMENTS::TransportBoundary* ele,
      Teuchos::ParameterList&           params,
      DRT::Discretization&              discretization,
      std::vector<int>&                 lm,
      Epetra_SerialDenseMatrix&         elemat1_epetra,
      Epetra_SerialDenseMatrix&         elemat2_epetra,
      Epetra_SerialDenseVector&         elevec1_epetra,
      Epetra_SerialDenseVector&         elevec2_epetra,
      Epetra_SerialDenseVector&         elevec3_epetra
  );


private:

  //! evaluate neumann condition
  int EvaluateNeumann(
      DRT::ELEMENTS::TransportBoundary*   ele,
      Teuchos::ParameterList&             params,
      DRT::Discretization&                discretization,
      DRT::Condition&                     condition,
      std::vector<int>&                   lm,
      Epetra_SerialDenseVector&           elevec1);

  //! evaluate an electrode kinetics boundary condition
  void EvaluateElectrodeKinetics(
     const DRT::Element*          ele,              ///< the actual boundary element
     Epetra_SerialDenseMatrix&    emat,             ///< element-matrix
     Epetra_SerialDenseVector&    erhs,             ///< element-rhs
     const std::vector<double>&   ephinp,           ///< actual conc. and pot. values
     const std::vector<double>&   ehist,            ///< element history vector
     double                       timefac,          ///< time factor
     Teuchos::RCP<const MAT::Material>  material,   ///< the material
     Teuchos::RCP<DRT::Condition>       cond,       ///< the condition
     const int                 nume,                ///< number of transfered electrons
     const std::vector<int>  stoich,                ///< stoichometry of the reaction
     const int             kinetics,                ///< desired electrode kinetics model
     const double              pot0,                ///< actual electrode potential on metal side
     const double               frt,                ///< factor F/RT
     const bool              iselch,                ///< problemtype elch or scatra
     const INPAR::ELCH::ElchType elchtype           ///< scatra problem type
     );

  //! evaluate elch electrode kinetics
  void CalcElchElectrodeKinetics(
      DRT::ELEMENTS::TransportBoundary* ele,
      Teuchos::ParameterList&           params,
      DRT::Discretization&              discretization,
      std::vector<int>&                 lm,
      Epetra_SerialDenseMatrix&         elemat1_epetra,
      Epetra_SerialDenseVector&         elevec1_epetra
      );

  //! provide electrode kinetics status information
  void ElectrodeStatus(
      const DRT::Element*           ele,      ///< the actual boundary element
      Teuchos::ParameterList&       params,   ///< the parameter list
      Teuchos::RCP<DRT::Condition>  cond,     ///< the condition
      const std::vector<double>&    ephinp,   ///< current conc. and potential values
      const std::vector<double>&    ephidtnp, ///< time derivative vector evaluated at t_{n+1}
      const int                     kinetics, ///< desired electrode kinetics model
      const std::vector<int>        stoich,   ///< stoichometry of the reaction
      const int                     nume,     ///<  number of transfered electrons
      const double                  pot0,     ///< actual electrode potential on metal side at t_{n+1}
      const double                  frt,      ///< factor F/RT
      const bool                    iselch,   ///< problemtype elch or scatra
      const double                  timefac   ///< factor due to time discretization
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

  //! flag for closure equation
  enum INPAR::ELCH::EquPot equpot_;

}; // class ScaTraEleBoundaryCalcElch

} //namsepace ELEMENTS
} // namespace DRT

#endif

