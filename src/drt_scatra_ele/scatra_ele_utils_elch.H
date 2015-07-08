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

#include "../drt_inpar/inpar_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleDiffManagerElch;

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

        //! evaluate electrode kinetics status information at integration point on domain or boundary element
        void EvaluateElectrodeStatusAtIntegrationPoint(
            const DRT::Element*                           ele,        ///< current element
            Epetra_SerialDenseVector&                     scalars,    ///< scalars to be computed
            const Teuchos::ParameterList&                 params,     ///< parameter list
            const Teuchos::RCP<DRT::Condition>&           cond,       ///< condition
            const std::vector<LINALG::Matrix<nen_,1> >&   conreact,   ///< concentration values of reactive species at element nodes
            const LINALG::Matrix<nen_,1>&                 pot,        ///< potential vector
            const LINALG::Matrix<nen_,1>&                 potdtnp,    ///< potential derivative vector
            const LINALG::Matrix<nen_,1>&                 funct,      ///< shape functions at integration point
            const int                                     zerocur,    ///< flag for zero current
            const int                                     kinetics,   ///< desired electrode kinetics model
            const std::vector<int>&                       stoich,     ///< stoichiometry of the reaction
            const int                                     nume,       ///< number of transferred electrons
            const double                                  pot0,       ///< actual electrode potential on metal side at t_{n+1}
            const double                                  frt,        ///< factor F/RT
            const double                                  timefac,    ///< factor due to time discretization
            const double                                  fac,        ///< integration factor
            const double                                  scalar,     ///< scaling factor for current related quantities
            const int                                     k           ///< index of evaluated scalar
        ) const;

        //! evaluate ion material
        void MatIon(
            const Teuchos::RCP<const MAT::Material>         material,     //!< ion material
            const int                                       k,            //!< ID of ion material
            const INPAR::ELCH::EquPot                       equpot,       //!< type of closing equation for electric potential
            const Teuchos::RCP<ScaTraEleDiffManagerElch>&   diffmanager   //!< diffusion manager
        );

      protected:

        //! protected constructor for singletons
        ScaTraEleUtilsElch(const int numdofpernode, const int numscal);

        //! destructor
        virtual ~ScaTraEleUtilsElch(){return;};

      private:

        //! number of degrees of freedom per node
        const int numdofpernode_;

        //! number of transported scalars
        const int numscal_;
    }; // class ScaTraEleUtilsElch
  } // namespace ELEMENTS
} // namespace DRT
#endif
