/*!----------------------------------------------------------------------
\file drt_dofset_fixed_size_pbc.H

\brief A modified set of degrees of freedom that always pretends to be
       of a certain size in order to reserve space for fields that vary
       in size, i.e. XFEM fluid. PBC version!

\level 1

\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/
#ifndef DOFSET_FIXED_SIZE_PBC_H
#define DOFSET_FIXED_SIZE_PBC_H

#include <vector>
#include <map>

#include <Epetra_IntVector.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

#include "../drt_lib/drt_dofset_fixed_size.H"
#include "../drt_lib/drt_dofset_pbc.H"


// forward declarations
namespace DRT
{
  class Discretization;

  /*!
  \brief A set of degrees of freedom

  \author
  */
  class FixedSizePBCDofSet : public FixedSizeDofSet, public PBCDofSet
  {
   public:
    /*!
    \brief Standard Constructor


    <pre>

    create a dofset that reserves a certain amount of dofGIDs

    </pre>

    \return void

    */
    FixedSizePBCDofSet(Teuchos::RCP<map<int, std::vector<int>>> couplednodes,
        int numreservedofpernode, int nodeindexrange)
        : FixedSizeDofSet(numreservedofpernode, nodeindexrange), PBCDofSet(couplednodes)
    {
      return;
    }

    /*!
    \brief Destructor

    */
    virtual ~FixedSizePBCDofSet() { return; }

    /// create a copy of this object
    virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new FixedSizePBCDofSet(*this)); }

    /// Get maximum GID of degree of freedom row map
    virtual int MaxAllGID() const { return PBCDofSet::MinAllGID() + numMyReservedDofs_; }

   protected:
   private:
  };  // class StdDofSet

}  // namespace DRT

#endif  // #ifndef DOFSET_INDEPENDENT_H
