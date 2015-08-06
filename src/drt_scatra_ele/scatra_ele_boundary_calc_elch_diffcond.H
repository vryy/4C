/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_diffcond.H

\brief evaluation of ScaTra boundary elements for diffusion-conduction formulation

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_BOUNDARY_CALC_ELCH_DIFFCOND_H
#define SCATRA_ELE_BOUNDARY_CALC_ELCH_DIFFCOND_H

#include "scatra_ele_boundary_calc_elch_electrode.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    class ScaTraEleDiffManagerElchDiffCond;
    class ScaTraEleParameterElchDiffCond;

    // class implementation
    template<DRT::Element::DiscretizationType distype> class ScaTraEleBoundaryCalcElchDiffCond : public ScaTraEleBoundaryCalcElchElectrode<distype>
    {
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalc<distype> my;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype> myelch;
      typedef DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype> myelectrode;

    public:

      //! singleton access method
      static ScaTraEleBoundaryCalcElchDiffCond<distype>* Instance(
          const int numdofpernode,
          const int numscal,
          const std::string& disname,
          bool create=true
          );

      //! destructor
      virtual ~ScaTraEleBoundaryCalcElchDiffCond(){return;};

      //! singleton destruction
      virtual void Done();

    private:

      /// private constructor for singletons
      ScaTraEleBoundaryCalcElchDiffCond(const int numdofpernode,const int numscal,const std::string& disname);

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

      //! evaluate scatra-scatra interface coupling condition
      void EvaluateS2ICoupling(
          const DRT::FaceElement*     ele,              ///< current boundary element
          Teuchos::ParameterList&     params,           ///< parameter list
          DRT::Discretization&        discretization,   ///< discretization
          std::vector<int>&           lm,               ///< location vector
          Epetra_SerialDenseMatrix&   eslavematrix,     ///< element matrix for slave side
          Epetra_SerialDenseMatrix&   emastermatrix,    ///< element matrix for master side
          Epetra_SerialDenseVector&   eslaveresidual    ///< element residual for slave side
      );

      //! evaluate electrode kinetics boundary condition
      void CalcElchBoundaryKinetics(
          DRT::FaceElement*                 ele,
          Teuchos::ParameterList&           params,
          DRT::Discretization&              discretization,
          std::vector<int>&                 lm,
          Epetra_SerialDenseMatrix&         elemat1_epetra,
          Epetra_SerialDenseVector&         elevec1_epetra,
          const double                      scalar
      );

      //! extract valence of species k from element material
      const double GetValence(
          const Teuchos::RCP<const MAT::Material>&   material,   //! element material
          const int                                  k           //! species number
      ) const;

      //! diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond> dmedc_;
    }; // class ScaTraEleBoundaryCalcElchDiffCond
  } // namespace ELEMENTS
} // namespace DRT
#endif
