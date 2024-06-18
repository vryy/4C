/*----------------------------------------------------------------------*/
/*! \file
 \brief adapter for artery network

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_ART_NET_HPP
#define FOUR_C_ADAPTER_ART_NET_HPP

#include "4C_config.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_result_test.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Adapter
{
  /// basic artery network adapter
  class ArtNet
  {
   public:
    /// constructor
    ArtNet(){};

    /// virtual destructor to support polymorph destruction
    virtual ~ArtNet() = default;

    /// initialization
    virtual void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname) = 0;

    /// initialization
    virtual void InitSaveState() = 0;

    // restart
    virtual void read_restart(int step, bool CoupledTo3D = false) = 0;

    // time integration
    virtual void Integrate(
        bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams) = 0;

    // test results
    virtual void TestResults() = 0;

    // create field test
    virtual Teuchos::RCP<Core::UTILS::ResultTest> CreateFieldTest() = 0;

    //! get discretization
    virtual Teuchos::RCP<Core::FE::Discretization> discretization() = 0;

    // get time step size
    virtual double Dt() const = 0;

    // output
    virtual void Output(bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams) = 0;

    // update of variables from n --> n+1
    virtual void TimeUpdate() = 0;

    // save state
    virtual void SaveState() = 0;

    // load state
    virtual void LoadState() = 0;

    // prepare the loop
    virtual void prepare_time_loop() = 0;

    // prepare step
    virtual void prepare_time_step() = 0;

    // solve
    virtual void Solve(Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams) = 0;

    virtual void assemble_mat_and_rhs() = 0;

    virtual void PrepareLinearSolve() = 0;

    /// direct access to system matrix
    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> SystemMatrix() = 0;

    //! right-hand side alias the dynamic force residual
    virtual Teuchos::RCP<const Epetra_Vector> RHS() const = 0;

    //! return pressure field at time n+1
    virtual Teuchos::RCP<const Epetra_Vector> Pressurenp() const = 0;

    //! iterative update of primary variable
    virtual void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc) = 0;

    // solve scalar transport in arteries
    virtual void SolveScatra() = 0;

    // set solve scalar transport-flag
    virtual void SetSolveScatra(const bool solvescatra) = 0;

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual Teuchos::RCP<const Core::LinAlg::MapExtractor> GetDBCMapExtractor() const = 0;

  };  // class ArtNet

}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
