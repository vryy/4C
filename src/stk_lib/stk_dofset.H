#ifndef STK_DOFSET_H
#define STK_DOFSET_H

#ifdef STKADAPTIVE

#include <map>
#include <vector>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_IntVector.h>

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>


namespace STK
{
  class Discretization;

  /// Mapping between dof gids and entity,field,position tuples
  /*!
    Finite elements lead to linear systems. We use Epetra objects to describe
    those, thus we need a map of global dof ids and the mapping between the
    stk mesh and the dof ids.
   */
  class DofSet
  {
  public:

    /// Setup of global dof numbers and reverse lookup
    void Setup( STK::Discretization& dis, const std::vector<stk::mesh::FieldBase*>& fields );

    /// Get number of dofs for given node
    int NumDof( stk::mesh::EntityKey key ) const
    {
      if ( key.rank()!=stk::mesh::Node )
        throw std::runtime_error( "unsupported rank" );
      Epetra_BlockMap noderowmap = numdfcolnodes_->Map();
      int gid = static_cast<int>( key.id() );
      int lid = noderowmap.LID( gid );
      if (lid==-1)
        return 0;
      return (*numdfcolnodes_)[lid];
    }

    /// Get the gid of all dofs of a node
    void Dof( stk::mesh::EntityKey key, std::vector<int> & dof ) const
    {
      if ( key.rank()!=stk::mesh::Node )
        throw std::runtime_error( "unsupported rank" );
      Epetra_BlockMap noderowmap = numdfcolnodes_->Map();
      int gid = static_cast<int>( key.id() );
      int lid = noderowmap.LID( gid );
      if (lid==-1)
      {
        dof.clear();
        return;
      }
      int idx = (*idxcolnodes_)[lid];
      int base = dofcolmap_->GID(idx);
      int count = (*numdfcolnodes_)[lid];
      dof.resize( count );
      for (int i=0; i<count; ++i)
        dof[i] = base+i;        // one node always has consecutive dof ids
    }

    /// Reverse lookup
    bool Reverse( int dofgid, stk::mesh::EntityKey & key, stk::mesh::FieldBase *& field, int & pos )
    {
      std::map<int, DofName>::iterator i = reverse_.find( dofgid );
      if ( i!=reverse_.end() )
      {
        DofName & dn = i->second;
        key   = dn.key;
        field = dn.field;
        pos   = dn.pos;
        return true;
      }
      return false;
    }

    /// Get degree of freedom row map
    const Epetra_Map & DofRowMap() const { return *dofrowmap_; }

    /// Get degree of freedom column map
    const Epetra_Map & DofColMap() const { return *dofcolmap_; }

  private:

    struct DofName
    {
      stk::mesh::EntityKey key;
      stk::mesh::FieldBase* field;
      int pos;
    };

    /// reverse lookup
    std::map<int, DofName> reverse_;

    /// unique row map of degrees of freedom (node and element dofs))
    Teuchos::RCP<Epetra_Map>                       dofrowmap_;

    /// unique column map of degrees of freedom (node and element dofs))
    Teuchos::RCP<Epetra_Map>                       dofcolmap_;

    /// number of dofs for each node
    Teuchos::RCP<Epetra_IntVector>                 numdfcolnodes_;

    /// column map lid of first dof for each node
    Teuchos::RCP<Epetra_IntVector>                 idxcolnodes_;

  };
}

#endif
#endif
