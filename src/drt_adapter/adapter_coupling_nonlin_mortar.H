/*!----------------------------------------------------------------------
\file adapter_coupling_nonlin_mortar.H

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/

#ifndef ADAPTER_COUPLING_NONLIN_MORTAR_H_
#define ADAPTER_COUPLING_NONLIN_MORTAR_H_

/*---------------------------------------------------------------------*
 | headers                                                 farah 10/14 |
 *---------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 10/14 |
 *---------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
}

namespace CONTACT
{
  class CoInterface;
}

namespace ADAPTER
{
  class CouplingNonLinMortar
  {
  public:

    /*!
    \brief Empty constructor

    */
    CouplingNonLinMortar();

    /*!
    \brief setup routine

    */
    void Setup(Teuchos::RCP<DRT::Discretization>   masterdis,
               Teuchos::RCP<DRT::Discretization>   slavedis,
               std::vector<int>                    coupleddof,
               const std::string&                  couplingcond);


  private:

    // communicator
    Teuchos::RCP<Epetra_Comm>          comm_;            /// communicator
    int                                myrank_;          /// my proc id

    Teuchos::RCP<Epetra_Map>           masterdofrowmap_; ///< map of master row dofs
    Teuchos::RCP<Epetra_Map>           slavedofrowmap_;  ///< map of slave row dofs

    Teuchos::RCP<CONTACT::CoInterface> interface_;       ///< interface

  };
}

#endif /* ADAPTER_COUPLING_NONLIN_MORTAR_H_ */
