#ifdef STKADAPTIVE

#ifndef STK_MESH_H
#define STK_MESH_H

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include <Epetra_MpiComm.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_Map.h>

#include <Teuchos_RCP.hpp>

#include "stk_io.H"

namespace STK
{

class MetaMesh
{
public:

  MetaMesh();

  void Commit();

  typedef phdmesh::exodus::FileSchema::AttributeField AttributeField ;

  stk::mesh::MetaData & MetaData() { return mesh_meta_data_; }

  stk::mesh::VectorField & Coordinates() { return *coordinates_field_; }

  stk::mesh::Part & UniversalPart() { return mesh_meta_data_.universal_part(); }
  stk::mesh::Part & OwnedPart()     { return mesh_meta_data_.locally_owned_part(); }
  stk::mesh::Part & SharedPart()    { return mesh_meta_data_.globally_shared_part(); }
  stk::mesh::Part & Active()   { return * active_; }
  stk::mesh::Part & Hanging()  { return * hanging_; }

  stk::mesh::Part & Quad4()    { return * quad4_; }
  stk::mesh::Part & Quad8()    { return * quad8_; }
  stk::mesh::Part & Quad9()    { return * quad9_; }
  stk::mesh::Part & Tri3()     { return * tri3_; }
  stk::mesh::Part & Tri6()     { return * tri6_; }
  stk::mesh::Part & Hex8()     { return * hex8_; }
  stk::mesh::Part & Hex20()    { return * hex20_; }
  stk::mesh::Part & Hex27()    { return * hex27_; }
  stk::mesh::Part & Tet4()     { return * tet4_; }
  stk::mesh::Part & Tet10()    { return * tet10_; }
  stk::mesh::Part & Wedge6()   { return * wedge6_; }
  stk::mesh::Part & Wedge15()  { return * wedge15_; }
  stk::mesh::Part & Pyramid5() { return * pyramid5_; }
  stk::mesh::Part & Line2()    { return * line2_; }
  stk::mesh::Part & Line3()    { return * line3_; }

  phdmesh::exodus::FileSchema & FileSchema() { return * file_schema_; }

private:

  stk::mesh::MetaData mesh_meta_data_;

  stk::mesh::Part * active_;
  stk::mesh::Part * hanging_;

  stk::mesh::Part * quad4_;     // 4 noded quadrilateral
  stk::mesh::Part * quad8_;     // 8 noded quadrilateral
  stk::mesh::Part * quad9_;     // 9 noded quadrilateral
  stk::mesh::Part * tri3_;      // 3 noded triangle
  stk::mesh::Part * tri6_;      // 6 noded triangle
  stk::mesh::Part * hex8_;      // 8 noded hexahedra
  stk::mesh::Part * hex20_;     // 20 noded hexahedra
  stk::mesh::Part * hex27_;     // 27 noded hexahedra
  stk::mesh::Part * tet4_;      // 4 noded tetrahedra
  stk::mesh::Part * tet10_;     // 10 noded tetrahedra
  stk::mesh::Part * wedge6_;    // 6 noded wedge
  stk::mesh::Part * wedge15_;   // 15 noded wedge
  stk::mesh::Part * pyramid5_;  // 5 noded pyramid
  stk::mesh::Part * line2_;     // 2 noded line
  stk::mesh::Part * line3_;     // 3 noded line

//   stk::mesh::Part * face_quad_;
//   stk::mesh::Part * edge_line_;

  stk::mesh::VectorField * coordinates_field_;
  AttributeField* element_attr_;
  phdmesh::exodus::FileSchema* file_schema_;
};


class Mesh
{
public:

  Mesh( MetaMesh & meta, stk::ParallelMachine pm );

  void Output( std::string filename, std::string title );
  void Output( std::string filename, std::string title, std::vector< const stk::mesh::FieldBase * > & out_fields );
  Teuchos::RCP<phdmesh::exodus::FileOutput> OutputContext( std::string filename, std::string title, std::vector< const stk::mesh::FieldBase * > & out_fields );

  void Modify();
  void Commit();

  stk::mesh::MetaData & MetaData() { return meta_.MetaData(); }
  stk::mesh::BulkData & BulkData() { return mesh_bulk_data_; }

  stk::mesh::VectorField & Coordinates() { return meta_.Coordinates(); }

  stk::mesh::Part & UniversalPart() { return meta_.UniversalPart(); }
  stk::mesh::Part & OwnedPart()     { return meta_.OwnedPart(); }
  stk::mesh::Part & SharedPart()    { return meta_.SharedPart(); }
  stk::mesh::Part & ActivePart()    { return meta_.Active(); }
  stk::mesh::Part & HangingPart()   { return meta_.Hanging(); }

  stk::mesh::Part & Quad4()    { return meta_.Quad4(); }
  stk::mesh::Part & Quad8()    { return meta_.Quad8(); }
  stk::mesh::Part & Quad9()    { return meta_.Quad9(); }
  stk::mesh::Part & Tri3()     { return meta_.Tri3(); }
  stk::mesh::Part & Tri6()     { return meta_.Tri6(); }
  stk::mesh::Part & Hex8()     { return meta_.Hex8(); }
  stk::mesh::Part & Hex20()    { return meta_.Hex20(); }
  stk::mesh::Part & Hex27()    { return meta_.Hex27(); }
  stk::mesh::Part & Tet4()     { return meta_.Tet4(); }
  stk::mesh::Part & Tet10()    { return meta_.Tet10(); }
  stk::mesh::Part & Wedge6()   { return meta_.Wedge6(); }
  stk::mesh::Part & Wedge15()  { return meta_.Wedge15(); }
  stk::mesh::Part & Pyramid5() { return meta_.Pyramid5(); }
  stk::mesh::Part & Line2()    { return meta_.Line2(); }
  stk::mesh::Part & Line3()    { return meta_.Line3(); }

  /** \brief  The parallel machine */
  stk::ParallelMachine parallel() const { return mesh_bulk_data_.parallel() ; }

  /** \brief  Size of the parallel machine */
  unsigned parallel_size()   const { return mesh_bulk_data_.parallel_size() ; }

  /** \brief  Rank of the parallel machine's local processor */
  unsigned parallel_rank()   const { return mesh_bulk_data_.parallel_rank() ; }

  void Print(std::ostream& stream);
  void Dump(std::string filename);

private:

  MetaMesh & meta_;
  stk::mesh::BulkData mesh_bulk_data_;

};

}

#endif
#endif
