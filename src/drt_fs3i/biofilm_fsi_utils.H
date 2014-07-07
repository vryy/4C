
#ifndef BIOFILMFSIUTILS_H_
#define BIOFILMFSIUTILS_H_

#include "../drt_lib/drt_discret.H"

namespace FS3I {
namespace Biofilm {
namespace UTILS {

  //! Update the material configuration of discretization with given displacement
  void updateMaterialConfigWithALE_Disp(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<const Epetra_Vector> disp);

  void ScatraChangeConfig(Teuchos::RCP<DRT::Discretization> scatradis, Teuchos::RCP<DRT::Discretization> dis,
                          Teuchos::RCP<Epetra_Vector> disp);

} /* namespace UTILS */
} /* namespace BiofilmFSI */
} /* namespace FS3I */

#endif
