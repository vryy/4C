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
#include <Epetra_Vector.h>

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

namespace LINALG
{
  class SparseMatrix;
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
    \brief Virtual destructor

    */
    virtual ~CouplingNonLinMortar(){};

    /*!
    \brief setup routine

    */
    void Setup(Teuchos::RCP<DRT::Discretization>   masterdis,
               Teuchos::RCP<DRT::Discretization>   slavedis,
               std::vector<int>                    coupleddof,
               const std::string&                  couplingcond);

    virtual void Evaluate(const std::string& statename,
                          const Teuchos::RCP<Epetra_Vector> vec,
                          bool onlyD);

    virtual void IntegrateD(const std::string& statename,
                            const Teuchos::RCP<Epetra_Vector> vec,
                            const Teuchos::RCP<Epetra_Vector> veclm);

    virtual void PrintInterface(std::ostream& os);

    virtual Teuchos::RCP<LINALG::SparseMatrix> DMatrix(){return D_;};
    virtual Teuchos::RCP<LINALG::SparseMatrix> MMatrix(){return M_;};

    virtual Teuchos::RCP<LINALG::SparseMatrix> DLinMatrix(){return DLin_;};
    virtual Teuchos::RCP<LINALG::SparseMatrix> MLinMatrix(){return MLin_;};

  private:

    // communicator
    Teuchos::RCP<Epetra_Comm>          comm_;            /// communicator
    int                                myrank_;          /// my proc id

    Teuchos::RCP<Epetra_Map>           masterdofrowmap_; ///< map of master row dofs
    Teuchos::RCP<Epetra_Map>           slavedofrowmap_;  ///< map of slave row dofs

    Teuchos::RCP<LINALG::SparseMatrix> D_;
    Teuchos::RCP<LINALG::SparseMatrix> DLin_;
    Teuchos::RCP<LINALG::SparseMatrix> M_;
    Teuchos::RCP<LINALG::SparseMatrix> MLin_;

    Teuchos::RCP<CONTACT::CoInterface> interface_;       ///< interface

  };
}

#endif /* ADAPTER_COUPLING_NONLIN_MORTAR_H_ */
