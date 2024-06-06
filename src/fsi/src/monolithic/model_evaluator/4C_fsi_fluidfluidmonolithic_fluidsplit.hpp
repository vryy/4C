/*------------------------------------------------------*/
/*! \file
\brief Control routine for monolithic fluid-fluid-fsi (fluidsplit)
 using XFEM and NOX.

\level 3

*/
/*------------------------------------------------------*/

#ifndef FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_HPP
#define FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_HPP

#include "4C_config.hpp"

#include "4C_fsi_monolithicfluidsplit.hpp"
#include "4C_inpar_xfem.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class FluidFluidFSI;
  class AleXFFsiWrapper;
}  // namespace Adapter

namespace FSI
{
  /// monolithic hybrid FSI algorithm with overlapping interface equations
  /*!
   * Monolithic fluid-fluid FSI with structure-handled interface motion, employing XFEM and NOX.
   * Fluid interface velocities are condensed.
   * \author kruse
   * \date 05/14
   */
  class FluidFluidMonolithicFluidSplit : public MonolithicFluidSplit
  {
    friend class FSI::FSIResultTest;

   public:
    /// constructor
    explicit FluidFluidMonolithicFluidSplit(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /// update subsequent fields, recover the Lagrange multiplier and relax the ALE-mesh
    void update() override;

    /// start a new time step
    void prepare_time_step() override;

    /// output routine accounting for Lagrange multiplier at the interface
    void output() override;

    /// read restart data (requires distinguation between fluid discretizations)
    void read_restart(int step) override;

   private:
    /// access type-cast pointer to problem-specific fluid-wrapper
    const Teuchos::RCP<Adapter::FluidFluidFSI>& fluid_field() { return fluid_; }

    /// access type-cast pointer to problem-specific ALE-wrapper
    const Teuchos::RCP<Adapter::AleXFFsiWrapper>& ale_field() { return ale_; }

    /// setup of extractor for merged Dirichlet maps
    void setup_dbc_map_extractor() override;


    /// type-cast pointer to problem-specific fluid-wrapper
    Teuchos::RCP<Adapter::FluidFluidFSI> fluid_;

    /// type-cast pointer to problem-specific ALE-wrapper
    Teuchos::RCP<Adapter::AleXFFsiWrapper> ale_;
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
