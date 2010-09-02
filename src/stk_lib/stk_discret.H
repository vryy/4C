#ifndef STK_DISCRET_H
#define STK_DISCRET_H

#ifdef STKADAPTIVE

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_mesh/fem/FieldDeclarations.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Comm.h>

#include "stk_fei.H"

#include "../drt_lib/drt_condition.H"

namespace MAT
{
  namespace PAR
  {
    class Bundle;
    class Parameter;
  }
}

namespace DRT
{
  class Discretization;
  class Condition;
  class Element;
  class DirichletExtractor;
}

namespace STK
{
  class MetaMesh;
  class Mesh;
  class Algorithm;

  namespace FEI
  {
    class AssembleStrategy;
    class DofSet;
  }

  /// STK based discretization
  /*!
    Contains the mesh and boundary (and other) conditions.

    Provides the connection to the linear system solver.
   */
  class Discretization
  {
  public:

    /// construct empty discretization
    explicit Discretization( const Epetra_Comm & comm );

    /// setup discretization to work with given algorithm
    void Setup( DRT::Discretization & dis, STK::Algorithm & algo );

    const Epetra_Comm & Comm() const { return comm_; }

    Mesh & GetMesh() { return *mesh_; }

    const std::map<int, Teuchos::ParameterList> * Condition( DRT::Condition::ConditionType gtype ) const;

    void AdaptMesh( const std::vector<stk::mesh::EntityKey> & refine,
                    const std::vector<stk::mesh::EntityKey> & unrefine );

    void ScatterFieldData( const Epetra_Vector & v, const std::vector<stk::mesh::FieldBase*> & fields );

    void GatherFieldData( const std::vector<stk::mesh::FieldBase*> & fields, Epetra_Vector & v );

    const DRT::DirichletExtractor & DirichletExtractor() const { return state_->DirichletExtractor(); }

    void EvaluateDirichlet( double time,
                            const std::vector<stk::mesh::FieldBase*> * v,
                            const std::vector<stk::mesh::FieldBase*> * dv=NULL,
                            const std::vector<stk::mesh::FieldBase*> * ddv=NULL )
    { state_->EvaluateDirichlet( time, v, dv, ddv ); }

    Teuchos::RCP<Epetra_CrsGraph> NodeGraph() { return state_->NodeGraph(); }

    /// Get number of dofs for given node
    int NumDof( stk::mesh::EntityKey key ) const
    { return state_->NumDof( key ); }

    /// Get the gid of all dofs of a node
    void Dof( stk::mesh::EntityKey key, std::vector<int> & dof ) const
    { state_->Dof( key, dof ); }

    /// Reverse lookup
    bool Reverse( int dofgid, stk::mesh::EntityKey & key, stk::mesh::FieldBase *& field, int & pos )
    { return state_->Reverse( dofgid, key, field, pos ); }

    /// Get degree of freedom row map
    const Epetra_Map & DofRowMap() const { return state_->DofRowMap(); }

    /// Get degree of freedom column map
    const Epetra_Map & DofColMap() const { return state_->DofColMap(); }

    void LocationVector( stk::mesh::Entity & e,
                         std::vector<int> & lm,
                         std::vector<int> & lmowner )
    { state_->LocationVector( e, lm, lmowner ); }

    MAT::PAR::Parameter* MaterialParameter( stk::mesh::Bucket & bucket );

  private:

    void Refine( const std::vector<stk::mesh::EntityKey> & eids );

    void Unrefine( const std::vector<stk::mesh::EntityKey> & eids );

    void CreateState();

    const Epetra_Comm & comm_;
    STK::Algorithm * algo_;

    Teuchos::RCP<MetaMesh> meta_;
    Teuchos::RCP<Mesh> mesh_;
    Teuchos::RCP<STK::FEI::DiscretizationState> state_;

    std::map<DRT::Condition::ConditionType, std::map<int, Teuchos::ParameterList> > conditions_;

    Teuchos::RCP<MAT::PAR::Bundle> materials_;

    std::map<stk::mesh::Part*, MAT::PAR::Parameter*> matpar_;
  };
}

#endif

#endif
