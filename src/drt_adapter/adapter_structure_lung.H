/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_constrained.H

\brief Structure field adapter for fsi airway simulations with
attached parenchyma balloon

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

#ifndef ADAPTER_STRUCTURE_LUNG_H
#define ADAPTER_STRUCTURE_LUNG_H
/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Vector.h>

#include "adapter_structure.H"
#include "adapter_structure_wrapper.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_structure/structure_utils_mapextractor.H"


namespace ADAPTER {

class StructureLung : public StructureWrapper
{
public:

  /// Constructor
  StructureLung(Teuchos::RCP<Structure> stru);

  /// List of fluid-structure volume constraints
  void ListLungVolCons(std::set<int>& LungVolConIDs,
                       int& MinLungVolConID);

  /// Initialize structural part of lung volume constraint
  void InitializeVolCon(Teuchos::RCP<Epetra_Vector> initvol,     ///< vector of initial volumes
                        Teuchos::RCP<Epetra_Vector> signvol,     ///< vector of signs of initial volumes
                        const int offsetID);                     ///< ID of first volume constraint -> offset

  /// Evaluate structural part of lung volume constraint
  void EvaluateVolCon(Teuchos::RCP<LINALG::BlockSparseMatrixBase> StructMatrix,
                      Teuchos::RCP<Epetra_Vector> StructRHS,
                      Teuchos::RCP<Epetra_Vector> CurrVols,
                      Teuchos::RCP<Epetra_Vector> SignVols,
                      Teuchos::RCP<Epetra_Vector> lagrMultVecRed,
                      const int offsetID);

  /// Write additional volume constraint stuff
  void WriteVolConRestart(Teuchos::RCP<Epetra_Vector> OldFlowRatesRed,
                          Teuchos::RCP<Epetra_Vector> OldVolsRed,
                          Teuchos::RCP<Epetra_Vector> OldLagrMultRed);

  /// Read additional volume constraint stuff
  void ReadVolConRestart(const int step,
                         Teuchos::RCP<Epetra_Vector> OldFlowRatesRed,
                         Teuchos::RCP<Epetra_Vector> OldVolsRed,
                         Teuchos::RCP<Epetra_Vector> OldLagrMultRed);

  /// Get MapExtractor for fsi <-> full map
  LINALG::MapExtractor& FSIInterface() { return fsiinterface_;  }

private:

  /// conditions that define the lung volume constraints
  std::vector<DRT::Condition*> constrcond_;

  /// map extractor for fsi <-> full map
  /// this is needed since otherwise "OtherMap" contains only dofs
  /// which are not part of a condition. however, asi dofs are of
  /// course also "inner" dofs for the structural field.
  LINALG::MapExtractor fsiinterface_;

};

}
#endif
#endif
