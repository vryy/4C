#ifdef STKADAPTIVE

#ifndef FLUID_ERRORESTIMATE_H
#define FLUID_ERRORESTIMATE_H

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

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_Vector.h>

#include <Teuchos_RCP.hpp>

namespace STK
{
  class Discretization;

  namespace FLD
  {

    class FluidEvaluateHelper;
    class FluidSideEvaluateHelper;

    class FluidJumpIntegrator
    {
    public:

      FluidJumpIntegrator( STK::Discretization & dis );

      void Integrate( stk::mesh::VectorField & coords,
                      stk::mesh::VectorField & velnp,
                      stk::mesh::ScalarField & pressure,
                      stk::mesh::ScalarField & error,
                      stk::mesh::ScalarField & volume,
                      stk::mesh::VectorField & resvel,
                      stk::mesh::ScalarField & respres );

    private:

      STK::Discretization & dis_;
    };


    class FluidEvaluateHelper
    {
    public:

      FluidEvaluateHelper( STK::Discretization & dis,
                           stk::mesh::Entity & e,
                           const CellTopologyData & topology,
                           stk::mesh::VectorField & coords,
                           stk::mesh::VectorField & velnp,
                           stk::mesh::ScalarField & pressure,
                           stk::mesh::VectorField * resvel=NULL,
                           stk::mesh::ScalarField * respres=NULL );

      void Integrate( stk::mesh::ScalarField & error, stk::mesh::ScalarField & volume );

      const CellTopologyData & topology() { return topology_; }

      const Epetra_SerialDenseMatrix & xyz()   const { return xyz_; }
      const Epetra_SerialDenseMatrix & vel()   const { return vel_; }
      const Epetra_SerialDenseVector & press() const { return press_; }

      void LocalCoordinates( Epetra_SerialDenseMatrix & rst, unsigned side_number );

      stk::mesh::Entity & element() const { return e_; }

      int NumNode() const { return xyz().N(); }

      double nu() const { return nu_; }

    private:

      void ExtractMyValues( Epetra_SerialDenseMatrix & emat,
                            stk::mesh::FieldBase & f,
                            const stk::mesh::PairIterRelation & nodes,
                            int dim );

      void ExtractMyValues( Epetra_SerialDenseVector & evec,
                            stk::mesh::FieldBase & f,
                            const stk::mesh::PairIterRelation & nodes );

      double Area();

      stk::mesh::Entity & e_;
      const CellTopologyData & topology_;

      double nu_;

      Epetra_SerialDenseMatrix xyz_;
      Epetra_SerialDenseMatrix vel_;
      Epetra_SerialDenseVector press_;

      Epetra_SerialDenseMatrix resvel_;
      Epetra_SerialDenseVector respress_;
    };


    class FluidSideEvaluateHelper
    {
    public:

      FluidSideEvaluateHelper( STK::FLD::FluidEvaluateHelper & e,
                               unsigned side_number );

      void Integrate( stk::mesh::ScalarField & error,
                      STK::FLD::FluidEvaluateHelper & le,
                      unsigned l_side_number,
                      STK::FLD::FluidEvaluateHelper & re,
                      unsigned r_side_number,
                      bool hanging );

      const Epetra_SerialDenseMatrix & xyz()   const { return xyz_; }
      const Epetra_SerialDenseMatrix & vel()   const { return vel_; }
      const Epetra_SerialDenseVector & press() const { return press_; }

      int NumNode() const { return xyz().N(); }

    private:

      void ExtractMyValues( const Epetra_SerialDenseMatrix & emat,
                            Epetra_SerialDenseMatrix & sidemat,
                            const CellTopologyData_Subcell & side,
                            const CellTopologyData * side_topology );

      void ExtractMyValues( const Epetra_SerialDenseVector & evec,
                            Epetra_SerialDenseVector & sidevec,
                            const CellTopologyData_Subcell & side,
                            const CellTopologyData * side_topology );

      Epetra_SerialDenseMatrix xyz_;
      Epetra_SerialDenseMatrix vel_;
      Epetra_SerialDenseVector press_;
    };

  }
}

#endif
#endif
