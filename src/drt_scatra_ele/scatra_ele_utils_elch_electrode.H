/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_utils_elch_electrode.H

\brief utility class supporting element evaluation for electrodes

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251

\level 2
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_UTILS_ELCH_ELECTRODE_H
#define SCATRA_ELE_UTILS_ELCH_ELECTRODE_H

#include "scatra_ele_utils_elch.H"

namespace DRT
{
  namespace ELEMENTS
  {
    // forward declaration
    class ScaTraEleDiffManagerElchElectrode;

    // class implementation
    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleUtilsElchElectrode : public ScaTraEleUtilsElch<distype>
    {
      //! abbreviation
      typedef ScaTraEleUtilsElch<distype> myelch;

     public:
      //! singleton access method
      static ScaTraEleUtilsElchElectrode<distype>* Instance(
          const int numdofpernode,     ///< number of degrees of freedom per node
          const int numscal,           ///< number of transported scalars per node
          const std::string& disname,  ///< name of discretization
          const ScaTraEleUtilsElchElectrode* delete_me = NULL  ///< creation/destruction indication
      );

      //! singleton destruction
      virtual void Done();

      //! evaluate electrode material
      void MatElectrode(const Teuchos::RCP<const MAT::Material> material,  //!< electrode material
          const double concentration,                                      //!< concentration
          const Teuchos::RCP<ScaTraEleDiffManagerElchElectrode>& diffmanager  //!< diffusion manager
      );

     protected:
      //! protected constructor for singletons
      ScaTraEleUtilsElchElectrode(
          const int numdofpernode,    ///< number of degrees of freedom per node
          const int numscal,          ///< number of transported scalars per node
          const std::string& disname  ///< name of discretization
      );

      //! destructor
      virtual ~ScaTraEleUtilsElchElectrode() { return; };
    };  // class ScaTraEleUtilsElchElectrode
  }     // namespace ELEMENTS
}  // namespace DRT
#endif
