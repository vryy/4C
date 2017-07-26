/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_utils.H

 \brief utils methods for for porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_UTILS_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_UTILS_H_

#include "../drt_inpar/inpar_poromultiphase.H"

#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>

namespace DRT
{
  class Discretization;
}

namespace ADAPTER
{
  class PoroMultiPhase;
}

namespace POROMULTIPHASE
{

/// POROMULTIPHASE::UTILS: Random stuff that might be helpful when dealing with poromultiphase problems
namespace UTILS
{

/// setup discretizations and dofsets
void SetupDiscretizationsAndFieldCoupling(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    int& nds_disp,
    int& nds_vel,
    int& nds_solidpressure);

//! exchange material pointers of both discretizations
void AssignMaterialPointers(
    const std::string& struct_disname,
    const std::string& fluid_disname
    );

/// create solution algorithm depending on input file
Teuchos::RCP<ADAPTER::PoroMultiPhase> CreatePoroMultiPhaseAlgorithm(
    INPAR::POROMULTIPHASE::SolutionSchemeOverFields solscheme, //!< solution scheme to build (i)
    const Teuchos::ParameterList& timeparams, //!< problem parameters (i)
    const Epetra_Comm& comm //!< communicator(i)
    );

//! Determine norm of vector
double CalculateVectorNorm(const enum INPAR::POROMULTIPHASE::VectorNorm norm, //!< norm to use
    const Teuchos::RCP<const Epetra_Vector> vect //!< the vector of interest
    );

} // UTILS
// Print the logo
void PrintLogo();
}



#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_UTILS_H_ */
