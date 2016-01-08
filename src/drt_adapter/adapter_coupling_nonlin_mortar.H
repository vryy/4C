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
#include "../drt_lib/drt_dserror.H"

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

    virtual void SetupSpringDashpot(
        Teuchos::RCP<DRT::Discretization> masterdis,
        Teuchos::RCP<DRT::Discretization> slavedis,
        Teuchos::RCP<DRT::Condition> spring,
        const int coupling_id,
        const Epetra_Comm& comm);

    virtual void IntegrateLinD(const std::string& statename,
                            const Teuchos::RCP<Epetra_Vector> vec,
                            const Teuchos::RCP<Epetra_Vector> veclm);

    virtual void IntegrateLinDM(const std::string& statename,
                            const Teuchos::RCP<Epetra_Vector> vec,
                            const Teuchos::RCP<Epetra_Vector> veclm);

    virtual void IntegrateAll(const std::string& statename,
                            const Teuchos::RCP<Epetra_Vector> vec,
                            const Teuchos::RCP<Epetra_Vector> veclm);

    virtual void PrintInterface(std::ostream& os);

    virtual void CreateP();

    virtual Teuchos::RCP<LINALG::SparseMatrix> DMatrix()
    {
      if(D_==Teuchos::null)
        dserror("ERROR: D Matrix is null pointer!");
      return D_;
    };

    virtual Teuchos::RCP<LINALG::SparseMatrix> DInvMatrix()
    {
      if(DInv_==Teuchos::null)
        dserror("ERROR: DInv Matrix is null pointer!");
      return DInv_;
    };

    virtual Teuchos::RCP<LINALG::SparseMatrix> MMatrix()
    {
      if(M_==Teuchos::null)
        dserror("ERROR: M Matrix is null pointer!");
      return M_;
    };

    virtual Teuchos::RCP<LINALG::SparseMatrix> PMatrix()
    {
      if(P_==Teuchos::null)
        dserror("ERROR: P Matrix is null pointer!");
      return P_;
    };

    virtual Teuchos::RCP<LINALG::SparseMatrix> DLinMatrix()
    {
      if(DLin_==Teuchos::null)
        dserror("ERROR: DLin Matrix is null pointer!");
      return DLin_;
    };

    virtual Teuchos::RCP<LINALG::SparseMatrix> MLinMatrix()
    {
      if(MLin_==Teuchos::null)
        dserror("ERROR: MLin Matrix is null pointer!");
      return MLin_;
    };

    virtual Teuchos::RCP<Epetra_Vector> Gap()
    {
      if(gap_==Teuchos::null)
      dserror("ERROR: gap vector is null pointer!");
      return gap_;
    };

    /// the mortar interface itself
    Teuchos::RCP<CONTACT::CoInterface> Interface() const { return interface_;}

  protected:
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

    /// check setup call
    const bool& IsSetup() const { return issetup_; };

    /// check init and setup call
    void CheckSetup() const
    {
      if ( !IsSetup())
        dserror("ERROR: Call Setup() first!");
    }

    bool                               issetup_;         ///< check for setup
    Teuchos::RCP<Epetra_Comm>          comm_;            ///< communicator
    int                                myrank_;          ///< my proc id

    Teuchos::RCP<Epetra_Map>           masterdofrowmap_; ///< map of master row dofs
    Teuchos::RCP<Epetra_Map>           slavedofrowmap_;  ///< map of slave row dofs
    Teuchos::RCP<Epetra_Map>           slavenoderowmap_; ///< map of slave row nodes

    Teuchos::RCP<LINALG::SparseMatrix> D_;               ///< D matrix
    Teuchos::RCP<LINALG::SparseMatrix> DInv_;            ///< inverted D matrix
    Teuchos::RCP<LINALG::SparseMatrix> DLin_;            ///< linearization of D matrix
    Teuchos::RCP<LINALG::SparseMatrix> M_;               ///< M matrix (master side)
    Teuchos::RCP<LINALG::SparseMatrix> MLin_;            ///< linearization of M matrix
    Teuchos::RCP<LINALG::SparseMatrix> P_;               ///< projection operator P = DInv * M

    Teuchos::RCP<Epetra_Vector>        gap_;             ///< gap vector

    Teuchos::RCP<CONTACT::CoInterface> interface_;       ///< interface
  };
}

#endif /* ADAPTER_COUPLING_NONLIN_MORTAR_H_ */
