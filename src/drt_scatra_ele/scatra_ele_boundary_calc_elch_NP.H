/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_NP.H

\brief evaluation of ScaTra boundary elements for Nernst-Planck formulation

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_NP_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_NP_H

#include "scatra_ele_boundary_calc_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template<DRT::Element::DiscretizationType distype> class ScaTraEleBoundaryCalcElchNP : public ScaTraEleBoundaryCalcElch<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype> myelch;

    public:

      //! singleton access method
      static ScaTraEleBoundaryCalcElchNP<distype>* Instance(
          const int numdofpernode,
          const int numscal,
          const std::string& disname,
          bool create=true
          );

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElchNP(){return;};

      //! singleton destruction
      virtual void Done();

    private:

      //! private constructor for singletons
      ScaTraEleBoundaryCalcElchNP(const int numdofpernode,const int numscal, const std::string& disname);

      //! evaluate Neumann boundary condition
      int EvaluateNeumann(
          DRT::FaceElement*                   ele,
          Teuchos::ParameterList&             params,
          DRT::Discretization&                discretization,
          DRT::Condition&                     condition,
          DRT::Element::LocationArray&        la,
          Epetra_SerialDenseVector&           elevec1,
          const double                        scalar
      );

      //! evaluate electrode kinetics boundary condition
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
          const double                        frt,        ///< factor F/RT
          const double                        scalar      ///< scaling factor for element matrix and residual contributions
      );

      //! extract valence of species k from element material
      const double GetValence(
          const Teuchos::RCP<const MAT::Material>&   material,   //! element material
          const int                                  k           //! species number
      ) const;
    }; // class ScaTraEleBoundaryCalcElchNP
  } // namespace ELEMENTS
} // namespace DRT
#endif
