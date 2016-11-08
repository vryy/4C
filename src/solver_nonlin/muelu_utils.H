/*----------------------------------------------------------------------------*/
/*!
\file muelu_utils.H

\brief Collection of utilities for MueLu

\level 2

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef MUELU_UTILS_H
#define MUELU_UTILS_H

/*----------------------------------------------------------------------------*/
/* headers */

// MueLu
#ifdef HAVE_MueLu
//Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

// MueLu
#include <MueLu_MapTransferFactory_fwd.hpp>
#include <MueLu_Utilities.hpp>

/* MueLu typedefs: header files for default types, must be included after all
 * other MueLu/Xpetra headers */
#include <MueLu_UseDefaultTypes.hpp> // => Scalar = double, LocalOrdinal = GlobalOrdinal = int

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */

// Teuchos
namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of methods */
namespace NLNSOL
{
namespace UTILS
{
  /*! \brief Convert Baci null space (based on ML) to Xpetra::MultiVector
   *
   *  \note This is not specific to the NLNSOL framework. It can also be used at
   *  other places in Baci like in the AMGnxn preconditioner.
   *
   *  ToDo (mayr) Move this to another UTILS place and use it for all MueLu
   *  applications in order to reduce code redundancy.
   *
   *  \author mayr.mt \date 10/2015
   */
  Teuchos::RCP<Xpetra::MultiVector<double,int,int,Node> >
  GetXpetraNullSpaceFromBaci(
      const Teuchos::ParameterList& mlparams, ///< parameter list with null space data created by Baci/ML
      const Teuchos::RCP<const Xpetra::Map<int,int,Node> > rowmap ///< dof row map
      );

} // namespace UTILS
} // namespace NLNSOL

#endif

#endif /* NLN_UTILS_H */
