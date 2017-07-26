/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_utils.H

 \brief helper functions/classes for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_UTILS_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_UTILS_H_

#include "../drt_inpar/inpar_poromultiphase_scatra.H"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

namespace POROMULTIPHASESCATRA
{
class PoroMultiPhaseScaTraBase;

/// POROMULTIPHASE::UTILS: Random stuff that might be helpful when dealing with poromultiphase problems
namespace UTILS
{

/// setup discretizations and dofsets
void SetupDiscretizationsAndFieldCoupling(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    const std::string& scatra_disname,
    int& ndsporo_disp,
    int& ndsporo_vel,
    int& ndsporo_solidpressure,
    int& ndsporofluid_scatra);

//! exchange material pointers of discretizations
void AssignMaterialPointers(
    const std::string& struct_disname,
    const std::string& fluid_disname,
    const std::string& scatra_disname
    );

/// create solution algorithm depending on input file
Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase>
CreatePoroMultiPhaseScatraAlgorithm(
    INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields solscheme, //!< solution scheme to build (i)
    const Teuchos::ParameterList& timeparams, //!< problem parameters (i)
    const Epetra_Comm& comm //!< communicator(i)
    );

//! Determine norm of vector
double CalculateVectorNorm(const enum INPAR::POROMULTIPHASESCATRA::VectorNorm norm, //!< norm to use
    const Teuchos::RCP<const Epetra_Vector> vect //!< the vector of interest
    );

} // UTILS
// Print the logo
void PrintLogo();
}



#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_UTILS_H_ */
