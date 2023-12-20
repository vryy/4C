/*----------------------------------------------------------------------*/
/*! \file

\brief transparent independent dofset

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_XFEM_DOFSET_TRANSPARENT_INDEPENDENT_H
#define BACI_XFEM_DOFSET_TRANSPARENT_INDEPENDENT_H

#include "baci_config.H"

#include "baci_lib_dofset_transparent_independent.H"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace CORE::GEO
{
  class CutWizard;
}

namespace XFEM
{
  /// Alias dofset that shares dof numbers with another dofset
  /*!
    A special set of degrees of freedom, implemented in order to assign the same degrees of freedom
    to nodes belonging to two discretizations. This way two discretizations can assemble into the
    same position of the system matrix. As internal variable it holds a source discretization
    (Constructor). If such a nodeset is assigned to a sub-discretization, its dofs are assigned
    according to the dofs of the source. The source discretization can be a xfem discretization. In
    this case this  should be called with a Fluidwizard not equal to zero to determine the  number
    of xfem dofs.

   */
  class XFEMTransparentIndependentDofSet : public virtual DRT::TransparentIndependentDofSet
  {
   public:
    /*!
      \brief Standard Constructor
     */
    explicit XFEMTransparentIndependentDofSet(Teuchos::RCP<DRT::Discretization> sourcedis,
        bool parallel, Teuchos::RCP<CORE::GEO::CutWizard> wizard);



   protected:
    int NumDofPerNode(const DRT::Node& node) const override;


   private:
    Teuchos::RCP<CORE::GEO::CutWizard> wizard_;
  };
}  // namespace XFEM

BACI_NAMESPACE_CLOSE

#endif
