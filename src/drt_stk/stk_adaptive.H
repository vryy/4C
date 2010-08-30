#ifdef STKADAPTIVE

#ifndef STK_ADAPTIVE_H
#define STK_ADAPTIVE_H

#include <map>
#include <vector>
#include <fstream>

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <Epetra_Vector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace LINALG
{
  class SparseOperator;
  class SparseMatrix;
}

namespace DRT
{
  class Discretization;
  class Condition;
  class Element;
  class DirichletExtractor;
}

/// adaptive meshes
namespace STK
{
  class MetaMesh;
  class Mesh;

class Adaptive
{
public:

  explicit Adaptive( DRT::Discretization & dis );

  virtual ~Adaptive();

  virtual void SetupSTKMesh();

  void RefineAll();

  void RefineHalf();

  void RefineFirst();

  void Refine( const std::vector<stk::mesh::EntityKey> & eids );
  void Unrefine( const std::vector<stk::mesh::EntityKey> & eids );

  Mesh & GetMesh() { return *mesh_; }

  DRT::Discretization & Discretization() { return dis_; }

  Teuchos::RCP<Epetra_Vector> GatherFieldData( const std::vector<stk::mesh::FieldBase*> & fields );

  void ScatterFieldData( Teuchos::RCP<Epetra_Vector> v, const std::vector<stk::mesh::FieldBase*> & fields );

#if 0
  void Evaluate( Teuchos::ParameterList&              params,
                 Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                 Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                 Teuchos::RCP<Epetra_Vector>          systemvector1,
                 Teuchos::RCP<Epetra_Vector>          systemvector2,
                 Teuchos::RCP<Epetra_Vector>          systemvector3 );
#endif

#if 0
  void EvaluateDirichlet( double time,
                          const std::vector<stk::mesh::FieldBase*> & v,
                          const std::vector<stk::mesh::FieldBase*> * dv  = NULL,
                          const std::vector<stk::mesh::FieldBase*> * ddv = NULL );
#endif

#if 0
  void Assemble( stk::mesh::Entity & e,
                 DRT::Element * actele,
                 Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
                 const Epetra_SerialDenseMatrix & Aele,
                 const std::vector<int> & lm,
                 const std::vector<int> & lmowner );
#endif

  //Teuchos::RCP<STK::FixedSparseMatrix> SetupMatrix();

protected:

  virtual void DefineFields( stk::mesh::MetaData & meta_data ) = 0;

private:

  void SyncDRT( bool refine );

  void SetupDirichletNodeMap();

  DRT::Discretization & dis_;
  Teuchos::RCP<MetaMesh> meta_;
  Teuchos::RCP<Mesh> mesh_;

  std::ofstream refine_log_;

  std::map<stk::mesh::Part*, DRT::Condition*> part_condition_;

  std::map<stk::mesh::Entity*, stk::mesh::Part*> dirichlet_condition_;
};

}

#endif
#endif
