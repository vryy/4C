#ifndef STK_DISCRET_STATE_H
#define STK_DISCRET_STATE_H

#ifdef STKADAPTIVE

#include <map>
#include <vector>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_IntVector.h>
#include <Epetra_CrsGraph.h>

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>

#include "stk_dofset.H"


namespace STK
{
  class Algorithm;
  class Discretization;
  class Mesh;

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

    void EvaluateDirichlet( double time,
                            const std::vector<stk::mesh::FieldBase*> * v,
                            const std::vector<stk::mesh::FieldBase*> * dv,
                            const std::vector<stk::mesh::FieldBase*> * ddv );

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

#if 0
    /// Get element row map
    const Epetra_Map & ElementRowMap() const { return *rowelementmap_; }

    /// Get element column map
    const Epetra_Map & ElementColMap() const { return *colelementmap_; }
#endif

    /// Get dof row dirichlet map
    const Epetra_Map & DofRowDirichletMap() const { return *rowdbcmap_; }

    /// Get dof column dirichlet map
    const Epetra_Map & DofColDirichletMap() const { return *coldbcmap_; }

    void LocationVector( stk::mesh::Entity & e,
                         std::vector<int> & lm,
                         std::vector<int> & lmowner );

  private:

    void CreateEpetraMaps();
    void CreateEntityMap( const stk::mesh::Selector & selector, stk::mesh::EntityRank rank, std::vector<int> & vnids );

    STK::Discretization& dis_;
    STK::DofSet dofset_;

    Teuchos::RCP<Epetra_Map> rownodemap_;
    Teuchos::RCP<Epetra_Map> colnodemap_;

#if 0
    Teuchos::RCP<Epetra_Map> rowelementmap_;
    Teuchos::RCP<Epetra_Map> colelementmap_;
#endif

    Teuchos::RCP<Epetra_Map> rowdbcmap_;
    Teuchos::RCP<Epetra_Map> coldbcmap_;
  };
}

#endif
#endif
