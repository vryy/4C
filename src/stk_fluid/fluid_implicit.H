#ifndef FLUID_IMPLICIT_H
#define FLUID_IMPLICIT_H

#ifdef STKADAPTIVE

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

#include <Epetra_Vector.h>

#include "../stk_lib/stk_algorithm.H"
#include "../stk_lib/stk_io.H"

#include "../drt_inpar/inpar_fluid.H"

namespace LINALG
{
  class FixedSparseMatrix;
  class Solver;
}

namespace STK
{
  class Discretization;

  namespace FEI
  {
    class AssembleStrategy;
  }

  /// stk based fluid algorithms
  namespace FLD
  {
    class Fluid;

    struct FluidState
    {
      FluidState( int counter, Fluid & fluid );

      int Counter() const { return counter_; }

      void Output( Fluid & fluid );

      const int counter_;

      /// (standard) system matrix
      Teuchos::RCP<LINALG::FixedSparseMatrix> sysmat_;

      Teuchos::RCP<Epetra_Vector> residual_;

      Teuchos::RCP<STK::FEI::AssembleStrategy> assemblestrategy_;

      Teuchos::RCP<phdmesh::exodus::FileOutput> m_exo;
    };

    /// central fluid algorithm
    class Fluid : public STK::Algorithm
    {
      friend struct FluidState;
    public:

      Fluid( STK::Discretization & dis, Teuchos::RCP<LINALG::Solver> solver );

      STK::Discretization & Discretization() const { return dis_; }

      void Integrate();

      virtual void declare_fields( stk::mesh::MetaData & meta );

      virtual void collect_unknowns( std::vector<stk::mesh::FieldBase*> & fields );

      virtual void notify_state_changed();

    private:

      void PrepareTimeStep();
      void LinearSolve();
      void NonlinearSolve();
      void TimeUpdate();
      void StatisticsAndOutput();

      void AdaptiveNonlinearSolve();
      void ErrorEstimate( std::vector<stk::mesh::EntityKey> & refine,
                          std::vector<stk::mesh::EntityKey> & unrefine );

      void ScatterFieldData( const Epetra_Vector & v,
                             stk::mesh::VectorField * vel,
                             stk::mesh::ScalarField * pres );

      void GatherFieldData( stk::mesh::VectorField * vel,
                            stk::mesh::ScalarField * pres,
                            Epetra_Vector & v );

      STK::Discretization & dis_;

      stk::mesh::VectorField * incvel_;
      stk::mesh::VectorField * resvel_;
      stk::mesh::VectorField * velnp_;
      stk::mesh::VectorField * veln_;
      stk::mesh::VectorField * velnm_;
      stk::mesh::VectorField * accnp_;
      stk::mesh::VectorField * accn_;
      stk::mesh::VectorField * hist_;

      stk::mesh::ScalarField * incpres_;
      stk::mesh::ScalarField * respres_;
      stk::mesh::ScalarField * pressure_;
      stk::mesh::ScalarField * constrained_;
      stk::mesh::ScalarField * hanging_;

      stk::mesh::ScalarField * error_;
      stk::mesh::ScalarField * volume_;

      //! @name time stepping variables
      double  time_;        ///< physical time
      int     step_;        ///< timestep
      int     stepmax_;     ///< maximal number of timesteps
      double  maxtime_;     ///< maximal physical computation time
      //@}

      //! @name time step sizes
      double dta_;
      double dtp_;
      //@}

      //! @name time-integration-scheme factors
      double theta_;
      //@}

      int refinestep_;
      int maxrefine_;

      double min_error_;
      double max_error_;
      double min_volume_;

      /// time algorithm flag
      INPAR::FLUID::TimeIntegrationScheme timealgo_;

      /// flag for physical type of fluid flow (standard: incompressible)
      INPAR::FLUID::PhysicalType physicaltype_;

      int dyntype_;
      bool alefluid_;

      //! use (or not) linearisation of reactive terms on the element
      INPAR::FLUID::LinearisationAction newton_;

      /// form of convective term
      std::string convform_;

      Teuchos::RCP<LINALG::Solver> solver_;

      Teuchos::RCP<FluidState> state_;
    };
  }
}

#endif
#endif
