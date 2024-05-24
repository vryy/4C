/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for monolithic fluid-fluid-fsi
(structuresplit) using XFEM and NOX

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_STRUCTURESPLIT_HPP
#define FOUR_C_FSI_FLUIDFLUIDMONOLITHIC_STRUCTURESPLIT_HPP

#include "4C_config.hpp"

#include "4C_fsi_monolithicstructuresplit.hpp"
#include "4C_inpar_xfem.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  class FluidFluidFSI;
  class AleXFFsiWrapper;
}  // namespace ADAPTER

namespace FSI
{
  /// monolithic hybrid FSI algorithm with overlapping interface equations
  /*!
   * Monolithic fluid-fluid FSI with fluid-handled interface motion, employing XFEM and NOX.
   * Structural interface displacements are condensed.
   * \author kruse
   * \date   05/14
   */
  class FluidFluidMonolithicStructureSplit : public MonolithicStructureSplit
  {
    friend class FSI::FSIResultTest;

   public:
    /// constructor
    explicit FluidFluidMonolithicStructureSplit(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /// update subsequent fields, recover the Lagrange multiplier and relax the ALE-mesh
    void Update() override;

    /// start a new time step
    void prepare_time_step() override;

   private:
    /// setup of extractor for merged Dirichlet maps
    void setup_dbc_map_extractor() override;

    /// access type-cast pointer to problem-specific fluid-wrapper
    const Teuchos::RCP<ADAPTER::FluidFluidFSI>& fluid_field() { return fluid_; }

    /// access type-cast pointer to problem-specific ALE-wrapper
    const Teuchos::RCP<ADAPTER::AleXFFsiWrapper>& ale_field() { return ale_; }

    /// type-cast pointer to problem-specific fluid-wrapper
    Teuchos::RCP<ADAPTER::FluidFluidFSI> fluid_;

    /// type-cast pointer to problem-specific ALE-wrapper
    Teuchos::RCP<ADAPTER::AleXFFsiWrapper> ale_;
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
