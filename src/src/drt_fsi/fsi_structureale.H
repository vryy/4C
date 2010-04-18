/*----------------------------------------------------------------------*/
/*!
\file fsi_structureale.H

\brief Solve structure only using FSI framework

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef FSI_PSEUDOFSISTRUCTURE_H
#define FSI_PSEUDOFSISTRUCTURE_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <AztecOO.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Map.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>

#include <vector>

#include "fsi_algorithm.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_mortar.H"


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
  StructureALE(Epetra_Comm& comm);

  /// outer level FSI time loop
  void Timeloop();

protected:

private:

  void Solve();

  Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv);

  /// interface structural operator

  /// coupling of structure and fluid at the interface, with mortar. Experimental!
  ADAPTER::CouplingMortar         coupsfm_;

  /// nodes at the fluid-structure interface match
  bool matchingnodes_;
};

}

#endif
#endif
