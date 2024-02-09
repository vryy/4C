/*---------------------------------------------------------------------*/
/*! \file

\brief A modified set of degrees of freedom that always pretends to be
       of a certain size in order to reserve space for fields that vary
       in size, i.e. XFEM field.

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_DOFSET_FIXED_SIZE_HPP
#define BACI_LIB_DOFSET_FIXED_SIZE_HPP

#include "baci_config.hpp"

#include "baci_lib_dofset.hpp"

#include <Epetra_Comm.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include <map>
#include <vector>

BACI_NAMESPACE_OPEN

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



    create a dofset that reserves a certain amount of dofGIDs


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



    /// create a copy of this object
    Teuchos::RCP<DofSet> Clone() override { return Teuchos::rcp(new FixedSizeDofSet(*this)); }

    /// Get maximum GID of degree of freedom row map
    int MaxAllGID() const override { return MinAllGID() + numMyReservedDofs_; }

    /// Get minimum GID of degree of freedom row map
    int MinAllGID() const override
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
    void GetReservedMaxNumDofperNode(int& maxnodenumdf) override
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

BACI_NAMESPACE_CLOSE

#endif  // LIB_DOFSET_FIXED_SIZE_H
