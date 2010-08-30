#ifndef STK_FEI_H
#define STK_FEI_H

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

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_dirichletextractor.H"

namespace STK
{
  class Algorithm;
  class Discretization;
  class Mesh;

  /// Finite Element Interface
  namespace FEI
  {

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


  /// Manage assembling inclusive Dirichlet conditions
  class AssembleStrategy : public DRT::AssembleStrategy
  {
  public:

    AssembleStrategy( STK::Discretization & dis,
                      Teuchos::RCP<const Epetra_Map> dbcmap,
                      Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                      Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                      Teuchos::RCP<Epetra_Vector>          systemvector1,
                      Teuchos::RCP<Epetra_Vector>          systemvector2,
                      Teuchos::RCP<Epetra_Vector>          systemvector3 );

    virtual Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( DRT::Discretization & dis, Teuchos::RCP<const Epetra_Map> dbcmap );

    Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( STK::Discretization & dis, Teuchos::RCP<const Epetra_Map> dbcmap );

    virtual void Assemble(LINALG::SparseOperator& sysmat,
                          int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner);

    virtual void Assemble(LINALG::SparseOperator& sysmat,
                          int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lmrow,
                          const std::vector<int>& lmrowowner,
                          const std::vector<int>& lmcol);

    virtual void Assemble(LINALG::SparseOperator& sysmat, double val, int rgid, int cgid);

    virtual void Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
                          const std::vector<int>& lm, const std::vector<int>& lmowner);

    virtual void Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
                          const std::vector<int>& lm, const std::vector<int>& lmowner);

  private:

    /// Dirichlet map
    Teuchos::RCP<const Epetra_Map> dbcmap_;

    struct DofData
    {
      std::set<int> realdofmap;
      double fact;
    };

    /// hn dof to real dof and factor
    std::map<int, DofData> hndofmap_;
  };


  /// State of a discretization at a given time
  /*!
    Links the stk mesh to a set of Epetra_Maps.

    State changes if mesh refinement or parallel distribution changes.

    State is initialized and used. If a state needs changing, it is replaced
    by a new version.

    This state has one dofset. There could be states with more than one.
   */
  class DiscretizationState
  {
  public:

    DiscretizationState( STK::Discretization& dis );

    void Setup( const std::vector<stk::mesh::FieldBase*>& fields );

    const DRT::DirichletExtractor & DirichletExtractor() const { return dirichlet_; }

    void EvaluateDirichlet( double time,
                            const std::vector<stk::mesh::FieldBase*> * v,
                            const std::vector<stk::mesh::FieldBase*> * dv,
                            const std::vector<stk::mesh::FieldBase*> * ddv );

    Teuchos::RCP<Epetra_CrsGraph> NodeGraph();

    /// Get number of dofs for given node
    int NumDof( stk::mesh::EntityKey key ) const
    { return dofset_.NumDof( key ); }

    /// Get the gid of all dofs of a node
    void Dof( stk::mesh::EntityKey key, std::vector<int> & dof ) const
    { dofset_.Dof( key, dof ); }

    /// Reverse lookup
    bool Reverse( int dofgid, stk::mesh::EntityKey & key, stk::mesh::FieldBase *& field, int & pos )
    { return dofset_.Reverse( dofgid, key, field, pos ); }

    /// Get degree of freedom row map
    const Epetra_Map & DofRowMap() const { return dofset_.DofRowMap(); }

    /// Get degree of freedom column map
    const Epetra_Map & DofColMap() const { return dofset_.DofColMap(); }

    /// Get node row map
    const Epetra_Map & NodeRowMap() const { return *rownodemap_; }

    /// Get node column map
    const Epetra_Map & NodeColMap() const { return *colnodemap_; }

    void LocationVector( stk::mesh::Entity & e,
                         std::vector<int> & lm,
                         std::vector<int> & lmowner );

  private:

    void CreateNodeMaps();
    void CreateNodeMap( const stk::mesh::Selector & selector, std::vector<int> & vnids );

    STK::Discretization& dis_;
    DofSet dofset_;
    DRT::DirichletExtractor dirichlet_;

    Teuchos::RCP<Epetra_Map> rownodemap_;
    Teuchos::RCP<Epetra_Map> colnodemap_;
  };

  }
}

#endif

#endif
