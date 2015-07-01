/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_utils_elch.H

\brief utility class supporting element evaluation for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_UTILS_ELCH_H
#define SCATRA_ELE_UTILS_ELCH_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // class implementation
    template<DRT::Element::DiscretizationType distype> class ScaTraEleUtilsElch
    {
      private:

        //! number of element nodes
        static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      public:

        //! singleton access method
        static ScaTraEleUtilsElch<distype>* Instance(
            const int numdofpernode,
            const int numscal,
            bool create=true
            );

        //! singleton destruction
        virtual void Done();

        //! evaluation of electrochemistry kinetics at integration point on domain or boundary element
        void EvaluateElchKineticsAtIntegrationPoint(
            const DRT::Element*                           ele,        ///< current element
            Epetra_SerialDenseMatrix&                     emat,       ///< element matrix
            Epetra_SerialDenseVector&                     erhs,       ///< element right-hand side vector
            const std::vector<LINALG::Matrix<nen_,1> >&   conreact,   ///< concentration values of reactive species at element nodes
            const LINALG::Matrix<nen_,1>&                 pot,        ///< el. potential values at element nodes
            const LINALG::Matrix<nen_,1>&                 phihist,    ///< element history vector for potential at electrode
            const double                                  timefac,    ///< time factor
            const double                                  fac,        ///< Gauss integration factor
            const LINALG::Matrix<nen_,1>&                 funct,      ///< shape functions at int. point
            const Teuchos::RCP<const MAT::Material>&      material,   ///< material
            const Teuchos::RCP<DRT::Condition>&           cond,       ///< condition
            const int                                     nume,       ///< number of transferred electrons
            const std::vector<int>&                       stoich,     ///< stoichiometry of the reaction
            const double                                  valence_k,  ///< valence of the single reactant
            const int                                     kinetics,   ///< desired electrode kinetics model
            const double                                  pot0,       ///< actual electrode potential on metal side
            const double                                  frt,        ///< factor F/RT
            const double                                  fns,        ///< factor fns = s_k / (nume * faraday * (-1))
            const double                                  scalar,     ///< scaling factor for element matrix and right-hand side vector contributions
            const int                                     k           ///< index of evaluated scalar
        ) const;

      private:

        //! private constructor for singletons
        ScaTraEleUtilsElch(const int numdofpernode, const int numscal);

        //! destructor
        virtual ~ScaTraEleUtilsElch(){return;};

        //! number of degrees of freedom per node
        const int numdofpernode_;

        //! number of transported scalars
        const int numscal_;
    }; // class ScaTraEleUtilsElch
  } // namespace ELEMENTS
} // namespace DRT
#endif
