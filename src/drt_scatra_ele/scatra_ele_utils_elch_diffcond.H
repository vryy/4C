/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_utils_elch_diffcond.H

\brief utility class supporting element evaluation for concentrated electrolytes

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_UTILS_ELCH_DIFFCOND_H
#define SCATRA_ELE_UTILS_ELCH_DIFFCOND_H

#include "scatra_ele_utils_elch_electrode.H"

#include "../drt_inpar/inpar_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declarations
    template<int NSD, int NEN> class ScaTraEleInternalVariableManagerElchDiffCond;
    class ScaTraEleDiffManagerElchDiffCond;

    // class implementation
    template<DRT::Element::DiscretizationType distype> class ScaTraEleUtilsElchDiffCond : public ScaTraEleUtilsElchElectrode<distype>
    {
        //! abbreviations
        typedef ScaTraEleUtilsElch<distype> myelch;
        typedef ScaTraEleUtilsElchElectrode<distype> myelectrode;

      public:

        //! singleton access method
        static ScaTraEleUtilsElchDiffCond<distype>* Instance(
            const int numdofpernode,      ///< number of degrees of freedom per node
            const int numscal,            ///< number of transported scalars per node
            const std::string& disname,   ///< name of discretization
            bool create=true              ///< creation/destruction flag
            );

        //! singleton destruction
        virtual void Done();

        //! evaluate electrolyte material
        void MatElchMat(
            const Teuchos::RCP<const MAT::Material>                                                         material,      //!< electrolyte material
            const Teuchos::RCP<ScaTraEleInternalVariableManagerElchDiffCond<myelch::nsd_,myelch::nen_> >&   varmanager,    //!< internal variable manager
            const Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&                                           diffmanager,   //!< diffusion manager
            INPAR::ELCH::DiffCondMat&                                                                       diffcondmat    //!< ion type
            );

        //! evaluate electrolyte phase
        void MatElchPhase(
            const Teuchos::RCP<const MAT::Material>                                                         material,      //!< electrolyte phase
            const Teuchos::RCP<ScaTraEleInternalVariableManagerElchDiffCond<myelch::nsd_,myelch::nen_> >&   varmanager,    //!< internal variable manager
            const Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&                                           diffmanager,   //!< diffusion manager
            INPAR::ELCH::DiffCondMat&                                                                       diffcondmat    //!< ion type
            );

        //! evaluate Newman material
        void MatNewman(
            const Teuchos::RCP<const MAT::Material>                 material,        //!< Newman material
            const double                                            concentration,   //!< concentration
            const Teuchos::RCP<ScaTraEleDiffManagerElchDiffCond>&   diffmanager      //!< diffusion manager
            );

      private:

        //! private constructor for singletons
        ScaTraEleUtilsElchDiffCond(
            const int numdofpernode,     ///< number of degrees of freedom per node
            const int numscal,           ///< number of transported scalars per node
            const std::string& disname   ///< name of discretization
            );

        //! destructor
        virtual ~ScaTraEleUtilsElchDiffCond(){return;};
    }; // class ScaTraEleUtilsElchDiffCond
  } // namespace ELEMENTS
} // namespace DRT
#endif
