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
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 10/14 |
 *---------------------------------------------------------------------*/
namespace DRT
{
  class Condition;
  class Discretization;
  class Node;
  class Element;
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
    \brief initialize routine

    */
    void Setup(
               Teuchos::RCP<DRT::Discretization>   masterdis,
               Teuchos::RCP<DRT::Discretization>   slavedis,
               std::vector<int>                    coupleddof,
               const std::string&                  couplingcond);

    /*!
    \brief Read Mortar Condition

    */
    virtual void ReadMortarCondition(
        Teuchos::RCP<DRT::Discretization>   masterdis,
        Teuchos::RCP<DRT::Discretization>   slavedis,
        std::vector<int>                    coupleddof,
        const std::string&                  couplingcond,
        Teuchos::ParameterList&             input,
        std::map<int, DRT::Node*>& mastergnodes,
        std::map<int, DRT::Node*>& slavegnodes,
        std::map<int, Teuchos::RCP<DRT::Element> >& masterelements,
        std::map<int, Teuchos::RCP<DRT::Element> >& slaveelements);

    /*!
    \brief Add Mortar Nodes

    */
    virtual void AddMortarNodes(
              Teuchos::RCP<DRT::Discretization>   masterdis,
              Teuchos::RCP<DRT::Discretization>   slavedis,
              std::vector<int>                    coupleddof,
              Teuchos::ParameterList&             input,
              std::map<int, DRT::Node*>& mastergnodes,
              std::map<int, DRT::Node*>& slavegnodes,
              std::map<int, Teuchos::RCP<DRT::Element> >& masterelements,
              std::map<int, Teuchos::RCP<DRT::Element> >& slaveelements,
              Teuchos::RCP<CONTACT::CoInterface>& interface,
              int numcoupleddof);

    /*!
    \brief Add Mortar Elements

    */
    virtual void AddMortarElements(
              Teuchos::RCP<DRT::Discretization>   masterdis,
              Teuchos::RCP<DRT::Discretization>   slavedis,
              Teuchos::ParameterList&             input,
              std::map<int, Teuchos::RCP<DRT::Element> >& masterelements,
              std::map<int, Teuchos::RCP<DRT::Element> >& slaveelements,
              Teuchos::RCP<CONTACT::CoInterface>& interface,
              int numcoupleddof);

    virtual void Evaluate(const std::string& statename,
                          const Teuchos::RCP<Epetra_Vector> vec,
                          bool onlyD);

    virtual void SetupSpringDashpot(
        Teuchos::RCP<DRT::Discretization> masterdis,
        Teuchos::RCP<DRT::Discretization> slavedis,
        Teuchos::RCP<DRT::Condition> spring,
        const int coupling_id,
        const Epetra_Comm& comm);

    virtual void IntegrateD(const std::string& statename,
                            const Teuchos::RCP<Epetra_Vector> vec,
                            const Teuchos::RCP<Epetra_Vector> veclm);

    virtual void PrintInterface(std::ostream& os);

    virtual Teuchos::RCP<LINALG::SparseMatrix> DMatrix(){return D_;};
    virtual Teuchos::RCP<LINALG::SparseMatrix> MMatrix(){return M_;};

    virtual Teuchos::RCP<LINALG::SparseMatrix> DLinMatrix(){return DLin_;};
    virtual Teuchos::RCP<LINALG::SparseMatrix> MLinMatrix(){return MLin_;};

    /// the mortar interface itself
    Teuchos::RCP<CONTACT::CoInterface> Interface() const { return interface_;}

  protected:
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
