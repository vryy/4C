#ifdef STKADAPTIVE

#ifndef STK_FLUID_H
#define STK_FLUID_H

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

#include "stk_adaptive.H"
#include "../stk_lib/stk_io.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_inpar/inpar_fluid.H"

namespace LINALG
{
  class FixedSparseMatrix;
  class Solver;
}

namespace FLD
{
  namespace UTILS
  {
    class VelPressExtractor;
  }
}

namespace DRT
{
  class ResultTest;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace STK
{
  class AssembleStrategy;
  class Fluid;

/// connection to one DRT discretization
class FluidDRT
{
public:

  FluidDRT( int counter, Fluid & fluid );

  int Counter() const { return counter_; }

  void Output( Fluid & fluid );

  const int counter_;

  /// (standard) system matrix
  Teuchos::RCP<LINALG::FixedSparseMatrix> sysmat_;

  Teuchos::RCP<Epetra_Vector> residual_;

  Teuchos::RCP<STK::AssembleStrategy> assemblestrategy_;

  Teuchos::RCP<FLD::UTILS::VelPressExtractor> velpressplitter_;

//   Teuchos::RCP<IO::DiscretizationWriter> output_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<DRT::DirichletExtractor> dbcmaps_;

  Teuchos::RCP<phdmesh::exodus::FileOutput> m_exo;
};


/// main fluid algorithm
class Fluid : public Adaptive
{
  friend class FluidDRT;

public:

  Fluid( DRT::Discretization & dis, Teuchos::RCP<LINALG::Solver> solver );

  virtual void SetupDRTMesh();

  void Integrate();

  Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

protected:

  virtual void DefineFields( stk::mesh::MetaData & meta_data );

private:

  void PrepareTimeStep();
  void LinearSolve();
  void NonlinearSolve();
  void TimeUpdate();
  void StatisticsAndOutput();

  void AdaptiveNonlinearSolve();
  void ErrorEstimate( std::vector<stk::mesh::EntityKey> & refine,
                      std::vector<stk::mesh::EntityKey> & unrefine );

  stk::mesh::VectorField * velnp_;
  stk::mesh::VectorField * veln_;
  stk::mesh::VectorField * velnm_;
  stk::mesh::VectorField * accnp_;
  stk::mesh::VectorField * accn_;
  stk::mesh::VectorField * hist_;

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

  /// connection to drt discretization objects
  Teuchos::RCP<FluidDRT> drt_;
};

}

#endif
#endif
