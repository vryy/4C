/*----------------------------------------------------------------------*/
/*!
\file xfem_dofset_transparent_independent.H

\brief transparent independent dofset

\level 1

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
*/
/*----------------------------------------------------------------------*/

#ifndef XFEMTRANSPARENTINDEPENDENTDOFSET_H
#define XFEMTRANSPARENTINDEPENDENTDOFSET_H

#include "../drt_lib/drt_dofset_transparent_independent.H"

namespace DRT
{
  class Discretization;
}

namespace GEO
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
        bool parallel, Teuchos::RCP<GEO::CutWizard> wizard);

    /*!
      \brief Destructor
     */
    virtual ~XFEMTransparentIndependentDofSet() { return; }

   protected:
    virtual int NumDofPerNode(const DRT::Node& node) const;


   private:
    Teuchos::RCP<GEO::CutWizard> wizard_;
  };
}  // namespace XFEM

#endif
