/*------------------------------------------------------*/
/*! \file
\brief Control routine for monolithic fluid-fluid-fsi (fluidsplit)
 using XFEM and NOX.

\level 3

*/
/*------------------------------------------------------*/

#ifndef BACI_FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_H
#define BACI_FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_H

#include "baci_config.H"

#include "baci_fsi_monolithicfluidsplit.H"
#include "baci_inpar_xfem.H"

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class FluidFluidFSI;
  class AleXFFsiWrapper;
}  // namespace ADAPTER

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
    void Update() override;

    /// start a new time step
    void PrepareTimeStep() override;

    /// output routine accounting for Lagrange multiplier at the interface
    void Output() override;

    /// read restart data (requires distinguation between fluid discretizations)
    void ReadRestart(int step) override;

   private:
    /// access type-cast pointer to problem-specific fluid-wrapper
    const Teuchos::RCP<ADAPTER::FluidFluidFSI>& FluidField() { return fluid_; }

    /// access type-cast pointer to problem-specific ALE-wrapper
    const Teuchos::RCP<ADAPTER::AleXFFsiWrapper>& AleField() { return ale_; }

    /// setup of extractor for merged Dirichlet maps
    void SetupDBCMapExtractor() override;


    /// type-cast pointer to problem-specific fluid-wrapper
    Teuchos::RCP<ADAPTER::FluidFluidFSI> fluid_;

    /// type-cast pointer to problem-specific ALE-wrapper
    Teuchos::RCP<ADAPTER::AleXFFsiWrapper> ale_;
  };
}  // namespace FSI

BACI_NAMESPACE_CLOSE

#endif  // FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_H
