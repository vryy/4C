/*---------------------------------------------------------------------*/
/*! \file

\brief A modified set of degrees of freedom that always pretends to be
       of a certain size in order to reserve space for fields that vary
       in size, i.e. XFEM field.

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_DOFSET_FIXED_SIZE_HPP
#define FOUR_C_DISCRETIZATION_DOFSET_FIXED_SIZE_HPP

#include "4C_config.hpp"

#include "4C_discretization_dofset.hpp"

#include <Epetra_Comm.h>
#include <Epetra_IntVector.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Discret
{
  class Discretization;

}

namespace Core::DOFSets
{
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
    void get_reserved_max_num_dofper_node(int& maxnodenumdf) override
    {
      if (maxnodenumdf > numMyReservedDofsperNode_)
      {
        FOUR_C_THROW(
            "FixedSizeDofSet::get_reserved_max_num_dofper_node: Not enough Dofs reserved!!!");
        return;
      }
      maxnodenumdf = numMyReservedDofsperNode_;

      return;
    }

    /// get the number of reserved DoF's (see also num_my_reserved_dofs_per_node())
    int NumMyReservedDofs() const { return numMyReservedDofs_; }

    /// get the number of reserved DoF's per node
    int num_my_reserved_dofs_per_node() const { return numMyReservedDofsperNode_; }

    /// minimal global dof id
    int minGID_;

   protected:
    /// Number of reserved dofs
    int numMyReservedDofs_;

    /// Number of reserved dofs per node
    int numMyReservedDofsperNode_;

   private:
  };  // class FixedSizeDofSet

}  // namespace Core::DOFSets

FOUR_C_NAMESPACE_CLOSE

#endif
