/*----------------------------------------------------------------------*/
/*! \file

\brief Solve structure only using FSI framework

\level 1


*/
/*----------------------------------------------------------------------*/



#ifndef FSI_PSEUDOFSISTRUCTURE_H
#define FSI_PSEUDOFSISTRUCTURE_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

#include "fsi_algorithm.H"

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

    Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv);

    /// interface structural operator

    /// coupling of structure and fluid at the interface, with mortar. Experimental!
    Teuchos::RCP<::ADAPTER::CouplingMortar> coupsfm_;

    /// nodes at the fluid-structure interface match
    bool matchingnodes_;
  };

}  // namespace FSI

#endif
