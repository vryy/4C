/*----------------------------------------------------------------------*/
/*!
\file fsi_fluidfluidmonolithic_structuresplit.cpp
\brief Control routine for monolithic fluid-fluid-fsi (structuresplit)
 using XFEM and NOX

<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265

            Raffaela Kruse
            kruse@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_FLUIDFLUIDMONOLITHICSTRUCTURESPLIT_H
#define FSI_FLUIDFLUIDMONOLITHICSTRUCTURESPLIT_H

#include "fsi_monolithicstructuresplit.H"
#include "../drt_inpar/inpar_xfem.H"

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
    explicit FluidFluidMonolithicStructureSplit(const Epetra_Comm& comm,
                                            const Teuchos::ParameterList& timeparams);

    /// update subsequent fields, recover the Lagrange multiplier and relax the ALE-mesh
    virtual void Update();

    /// start a new time step
    virtual void PrepareTimeStep();

  private:

    /// setup of extractor for merged Dirichlet maps
    void SetupDBCMapExtractor();

    /// flag, that indicates ALE-relaxation
    bool relaxing_ale_;

    /// number of timesteps in between ALE-relaxation
    int relaxing_ale_every_;
  };
}

#endif /* FSI_FLUIDFLUIDMONOLITHIC_STRUCTURESPLIT_H_ */
