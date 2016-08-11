/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_utils.H

 \brief utils methods for for porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_UTILS_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_UTILS_H_

#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>

namespace DRT
{
  class Discretization;
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

}
}



#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_UTILS_H_ */
