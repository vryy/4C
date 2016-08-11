/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_utils.H

 \brief helper methods/classes for multiphase porous fluid flow problems

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_UTILS_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_UTILS_H_


#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>

namespace DRT
{
  class Discretization;
}

namespace POROFLUIDMULTIPHASE
{
/// POROFLUIDMULTIPHASE::UTILS: Random stuff that might be helpful when dealing with poromultiphase problems
namespace UTILS
{

/// setup second materials for porosity evaluation within solid phase
void SetupMaterial(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& fluid_disname);


/// convert a dof based vector to a node based multi vector
/*!
  For postprocessing, only vectors based on the primary dof set
  of the discretization can be used. Hence, all other EpetraVectors
  based on secondary dof sets are copied to EpetraMultiVectors with one
  node based vector for each component.

  This method can be deleted, if the post processors would be adapted to
  handle secondary dof sets.

  \param dis              : discretization, the vector is based on
  \param vector           : vector to convert
  \param nds              : number of the dof set the map of the vector corresponds to
  \param numdofpernode    : number of dofs per node of the vector (assumed to be equal for all nodes)
 */
Teuchos::RCP<Epetra_MultiVector> ConvertDofVectorToNodeBasedMultiVector(
    const DRT::Discretization& dis,
    const Epetra_Vector& vector,
    const int nds,
    const int numdofpernode);

}
}


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_UTILS_H_ */
