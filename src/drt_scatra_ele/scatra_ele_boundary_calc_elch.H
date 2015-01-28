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
// forward declarations
class ScaTraEleParameterElch;
class ScaTraEleDiffManagerElchDiffCond;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleBoundaryCalcElch : public ScaTraBoundaryImpl<distype>
{
  typedef DRT::ELEMENTS::ScaTraBoundaryImpl<distype> my;

public:

  /// Singleton access method
  static ScaTraEleBoundaryCalcElch<distype> * Instance(const int numdofpernode, const int numscal, bool create=true );

  /// Empty destructor
  virtual ~ScaTraEleBoundaryCalcElch(){return;};

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

  //! evaluate neumann condition
  int EvaluateNeumann(
      DRT::ELEMENTS::TransportBoundary*   ele,
      Teuchos::ParameterList&             params,
      DRT::Discretization&                discretization,
      DRT::Condition&                     condition,
      std::vector<int>&                   lm,
      Epetra_SerialDenseVector&           elevec1
  );

  //! evaluate an electrode kinetics boundary condition
  void EvaluateElchBoundaryKinetics(
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
     const double                        frt         ///< factor F/RT
  );

  //! evaluate scatra-scatra interface coupling condition
  void EvaluateS2ICoupling(
      const DRT::Element*         ele,              ///< current boundary element
      Teuchos::ParameterList&     params,           ///< parameter list
      DRT::Discretization&        discretization,   ///< discretization
      std::vector<int>&           lm,               ///< location vector
      Epetra_SerialDenseMatrix&   eslavematrix,     ///< element matrix for slave side
      Epetra_SerialDenseMatrix&   emastermatrix,    ///< element matrix for master side
      Epetra_SerialDenseVector&   eslaveresidual    ///< element residual for slave side
  );

  //! evaluate elch electrode kinetics
  void CalcElchBoundaryKinetics(
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
      const std::vector<int>        stoich,   ///< stoichiometry of the reaction
      const int                     nume,     ///<  number of transferred electrons
      const double                  pot0,     ///< actual electrode potential on metal side at t_{n+1}
      const double                  frt,      ///< factor F/RT
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

  //! extract valence of species k from element material
  const double GetValence(
      const Teuchos::RCP<const MAT::Material>&   material,   //! element material
      const int                                  k           //! species number
  ) const;

  //! equilibrium electric potential difference at electrode-electrolyte interface
  void EquilibriumPotentialDifference(
      const Teuchos::RCP<DRT::Condition>&   condition,       //! boundary condition
      const double&                         emasterphiint,   //! concentration of intercalated Lithium at electrode surface at current Gauss point
      double&                               epd,             //! equilibrium electric potential difference to be computed at current Gauss point
      double&                               epdderiv         //! derivative of equilibrium electric potential difference to be computed at current Gauss point
  ) const;

  /*========================================================================*/
  //! @name flags and enums
  /*========================================================================*/

  //! pointer to parameter class elch
  DRT::ELEMENTS::ScaTraEleParameterElch* elchpara_;

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  //! flag for closure equation
  enum INPAR::ELCH::EquPot equpot_;

  //! manager for diffusion
  Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc_;
}; // class ScaTraEleBoundaryCalcElch
} //namsepace ELEMENTS
} // namespace DRT
#endif
