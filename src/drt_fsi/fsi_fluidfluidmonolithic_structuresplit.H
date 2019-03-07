/*----------------------------------------------------------------------*/
/*!
\file fsi_fluidfluidmonolithic_structuresplit.H
\brief Control routine for monolithic fluid-fluid-fsi
(structuresplit) using XFEM and NOX

\level 2

<pre>
\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_FLUIDFLUIDMONOLITHICSTRUCTURESPLIT_H
#define FSI_FLUIDFLUIDMONOLITHICSTRUCTURESPLIT_H

#include "fsi_monolithicstructuresplit.H"
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
    virtual void Update();

    /// start a new time step
    virtual void PrepareTimeStep();

   private:
    /// setup of extractor for merged Dirichlet maps
    void SetupDBCMapExtractor();

    /// access type-cast pointer to problem-specific fluid-wrapper
    const Teuchos::RCP<ADAPTER::FluidFluidFSI>& FluidField() { return fluid_; }

    /// access type-cast pointer to problem-specific ALE-wrapper
    const Teuchos::RCP<ADAPTER::AleXFFsiWrapper>& AleField() { return ale_; }

    /// type-cast pointer to problem-specific fluid-wrapper
    Teuchos::RCP<ADAPTER::FluidFluidFSI> fluid_;

    /// type-cast pointer to problem-specific ALE-wrapper
    Teuchos::RCP<ADAPTER::AleXFFsiWrapper> ale_;
  };
}  // namespace FSI

#endif /* FSI_FLUIDFLUIDMONOLITHIC_STRUCTURESPLIT_H_ */
