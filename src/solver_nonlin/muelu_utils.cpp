/*----------------------------------------------------------------------------*/
/*!
\file muelu_utils.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// standard
#include <string>

// Teuchos
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_ParameterList.hpp>

// baci
#include "muelu_utils.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

#ifdef HAVE_MueLu
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Xpetra::MultiVector<double, int, int, Node>> NLNSOL::UTILS::GetXpetraNullSpaceFromBaci(
    const Teuchos::ParameterList& mlparams,
    const Teuchos::RCP<const Xpetra::Map<int, int, Node>> rowmap)
{
  // extract null space data from parameter list
  const int numdof = mlparams.get<int>("PDE equations");          // number of DOFs per node
  const int numnsv = mlparams.get<int>("null space: dimension");  // number of null space vectors
  Teuchos::RCP<std::vector<double>> nsdata = mlparams.get<Teuchos::RCP<std::vector<double>>>(
      "nullspace");  // null space vectors (all in one vector)

  // some safety checks
  if (numdof < 1 or numnsv < 1) dserror("PDE equations or null space dimension wrong.");
  if (nsdata.is_null()) dserror("Null space data is empty");

  // Prepare null space vector for MueLu
  Teuchos::RCP<Xpetra::MultiVector<double, int, int, Node>> nsp =
      Xpetra::MultiVectorFactory<double, int, int, Node>::Build(rowmap, numnsv, true);

  // copy null space vectors from 'nsdata' to 'nsp'
  for (unsigned int i = 0; i < Teuchos::as<size_t>(numnsv); ++i)
  {
    Teuchos::ArrayRCP<double> nspvector = nsp->getDataNonConst(i);
    const int mylength = nsp->getLocalLength();
    for (int j = 0; j < mylength; ++j)
    {
      nspvector[j] = (*nsdata)[i * mylength + j];
    }
  }

  return nsp;
}
#endif
