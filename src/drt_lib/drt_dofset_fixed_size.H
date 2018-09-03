/*!----------------------------------------------------------------------
\file drt_dofset_fixed_size.H

\brief A modified set of degrees of freedom that always pretends to be
       of a certain size in order to reserve space for fields that vary
       in size, i.e. XFEM field.

\level 1

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/
#ifndef DOFSET_FIXED_SIZE_H
#define DOFSET_FIXED_SIZE_H

#include <vector>
#include <map>

#include <Epetra_IntVector.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_dofset.H"


// forward declarations
namespace DRT
{
  class Discretization;

  /*!
  \brief A set of degrees of freedom

  \author
  */
  class FixedSizeDofSet : virtual public DofSet
  {
   public:
    /*!
    \brief Standard Constructor


    <pre>

    create a dofset that reserves a certain amount of dofGIDs

    </pre>

    \return void

    */
    FixedSizeDofSet(int numreservedofpernode, int nodeindexrange)
        : DofSet(),
          numMyReservedDofs_(numreservedofpernode * nodeindexrange),
          numMyReservedDofsperNode_(numreservedofpernode)
    {
      minGID_ = -1;
      return;
    }

    /*!
    \brief Destructor

    */
    virtual ~FixedSizeDofSet() { return; }

    /// create a copy of this object
    virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new FixedSizeDofSet(*this)); }

    /// Get maximum GID of degree of freedom row map
    virtual int MaxAllGID() const { return MinAllGID() + numMyReservedDofs_; }

    /// Get minimum GID of degree of freedom row map
    virtual int MinAllGID() const
    {
      int mymindof;
      if (minGID_ == -1)
        mymindof = dofrowmap_->MinAllGID();
      else
      {
        int minall = dofrowmap_->MinAllGID();
        mymindof = (minGID_ <= minall) ? minGID_ : minall;
      }
      return mymindof;
    }

    /// set the minimal global id
    virtual void SetMinGID(int mingid) { minGID_ = mingid; }

    /// Get Reserved Max Number Dofs per Node
    virtual void GetReservedMaxNumDofperNode(int& maxnodenumdf)
    {
      if (maxnodenumdf > numMyReservedDofsperNode_)
      {
        dserror("FixedSizeDofSet::GetReservedMaxNumDofperNode: Not enough Dofs reserved!!!");
        return;
      }
      maxnodenumdf = numMyReservedDofsperNode_;

      return;
    }

    /// get the number of reserved DoF's (see also NumMyReservedDofsPerNode())
    int NumMyReservedDofs() const { return numMyReservedDofs_; }

    /// get the number of reserved DoF's per node
    int NumMyReservedDofsPerNode() const { return numMyReservedDofsperNode_; }

    /// minimal global dof id
    int minGID_;

   protected:
    /// Number of reserved dofs
    int numMyReservedDofs_;

    /// Number of reserved dofs per node
    int numMyReservedDofsperNode_;

   private:
  };  // class FixedSizeDofSet

}  // namespace DRT

#endif  // #ifndef DOFSET_FIXED_SIZE_H
