/*------------------------------------------------------*/
/*! \file
\brief Control routine for monolithic fluid-fluid-fsi (fluidsplit)
 using XFEM and NOX.

\level 3

\maintainer  Matthias Mayr
*/
/*------------------------------------------------------*/

#ifndef FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_H_
#define FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_H_

#include "fsi_monolithicfluidsplit.H"
#include "../drt_inpar/inpar_xfem.H"

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
    virtual void Update();

    /// start a new time step
    virtual void PrepareTimeStep();

    /// output routine accounting for Lagrange multiplier at the interface
    virtual void Output();

    /// read restart data (requires distinguation between fluid discretizations)
    virtual void ReadRestart(int step);

   private:
    /// access type-cast pointer to problem-specific fluid-wrapper
    const Teuchos::RCP<ADAPTER::FluidFluidFSI>& FluidField() { return fluid_; }

    /// access type-cast pointer to problem-specific ALE-wrapper
    const Teuchos::RCP<ADAPTER::AleXFFsiWrapper>& AleField() { return ale_; }

    /// setup of extractor for merged Dirichlet maps
    void SetupDBCMapExtractor();


    /// type-cast pointer to problem-specific fluid-wrapper
    Teuchos::RCP<ADAPTER::FluidFluidFSI> fluid_;

    /// type-cast pointer to problem-specific ALE-wrapper
    Teuchos::RCP<ADAPTER::AleXFFsiWrapper> ale_;
  };
}  // namespace FSI

#endif /* FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_H_ */
