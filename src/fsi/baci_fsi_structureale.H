/*----------------------------------------------------------------------*/
/*! \file

\brief Solve structure only using FSI framework

\level 1


*/
/*----------------------------------------------------------------------*/



#ifndef BACI_FSI_STRUCTUREALE_H
#define BACI_FSI_STRUCTUREALE_H

#include "baci_coupling_adapter_mortar.H"
#include "baci_fsi_algorithm.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// forward declarations

namespace ADAPTER
{
  class CouplingMortar;
}



/// FSI: Fluid Structure Interaction
/*!
  This is the place where we bring fluid and structure together.
 */
namespace FSI
{
  /// Coupling structure and ALE, while keeping fluid velocities at zero
  /*!

    \author tk
    \date 06/07
   */
  class StructureALE : public Algorithm
  {
   public:
    /// setup FSI algorithm
    StructureALE(const Epetra_Comm& comm);

    /// outer level FSI time loop
    void Timeloop();

   protected:
   private:
    void Solve();

    Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) override;

    /// interface structural operator

    /// coupling of structure and fluid at the interface, with mortar. Experimental!
    Teuchos::RCP<CORE::ADAPTER::CouplingMortar> coupsfm_;

    /// nodes at the fluid-structure interface match
    bool matchingnodes_;
  };

}  // namespace FSI

BACI_NAMESPACE_CLOSE

#endif
