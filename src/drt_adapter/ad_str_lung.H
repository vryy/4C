/*----------------------------------------------------------------------*/
/*! \file

\brief Structure field adapter for fsi airway simulations with
attached parenchyma balloon

\maintainer Carolin Geitner

\level 3

*/

/*----------------------------------------------------------------------*/
/* macros */


#ifndef ADAPTER_STRUCTURE_LUNG_H
#define ADAPTER_STRUCTURE_LUNG_H
/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

#include <set>

#include "ad_str_fsiwrapper.H"


// forward declarations
namespace LINALG
{
  class MapExtractor;
}

namespace DRT
{
  class Condition;
}


namespace ADAPTER
{
  class StructureLung : public FSIStructureWrapper
  {
   public:
    /// Constructor
    StructureLung(Teuchos::RCP<Structure> stru);

    /// List of fluid-structure volume constraints
    void ListLungVolCons(std::set<int>& LungVolConIDs, int& MinLungVolConID);

    /// Initialize structural part of lung volume constraint
    void InitializeVolCon(Teuchos::RCP<Epetra_Vector> initvol,  ///< vector of initial volumes
        Teuchos::RCP<Epetra_Vector> signvol,  ///< vector of signs of initial volumes
        const int offsetID);                  ///< ID of first volume constraint -> offset

    /// Evaluate structural part of lung volume constraint
    void EvaluateVolCon(Teuchos::RCP<LINALG::BlockSparseMatrixBase> StructMatrix,
        Teuchos::RCP<Epetra_Vector> StructRHS, Teuchos::RCP<Epetra_Vector> CurrVols,
        Teuchos::RCP<Epetra_Vector> SignVols, Teuchos::RCP<Epetra_Vector> lagrMultVecRed,
        const int offsetID);

    /// Write additional forces due to volume constraint
    void OutputForces(Teuchos::RCP<Epetra_Vector> Forces);

    /// Write additional volume constraint stuff
    void WriteVolConRestart(Teuchos::RCP<Epetra_Vector> OldFlowRatesRed,
        Teuchos::RCP<Epetra_Vector> OldVolsRed, Teuchos::RCP<Epetra_Vector> OldLagrMultRed);

    /// Read additional volume constraint stuff
    void ReadVolConRestart(const int step, Teuchos::RCP<Epetra_Vector> OldFlowRatesRed,
        Teuchos::RCP<Epetra_Vector> OldVolsRed, Teuchos::RCP<Epetra_Vector> OldLagrMultRed);

    /// Get MapExtractor for fsi <-> full map
    Teuchos::RCP<const LINALG::MapExtractor> FSIInterface() { return fsiinterface_; }

    /// Get map of volume coupling dofs
    const Teuchos::RCP<const Epetra_Map> LungConstrMap() { return lungconstraintmap_; }

   private:
    /// conditions that define the lung volume constraints
    std::vector<DRT::Condition*> constrcond_;

    /// conditions that define the structure ale coupling at the outlets
    std::vector<DRT::Condition*> asicond_;

    /// map containing all dofs related to volume coupling (i.e. dofs of the
    /// enclosing boundary)
    Teuchos::RCP<const Epetra_Map> lungconstraintmap_;

    /// map extractor for fsi <-> full map
    /// this is needed since otherwise "OtherMap" contains only dofs
    /// which are not part of a condition. however, asi dofs are of
    /// course also "inner" dofs for the structural field.
    Teuchos::RCP<LINALG::MapExtractor> fsiinterface_;
  };

}  // namespace ADAPTER
#endif
