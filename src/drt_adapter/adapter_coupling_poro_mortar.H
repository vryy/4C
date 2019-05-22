/*!----------------------------------------------------------------------


 \brief coupling adapter for porous meshtying

// Masterthesis of h.Willmann under supervision of Anh-Tu Vuong and Christoph Ager
// Originates from ADAPTER::CouplingNonLinMortar

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
*----------------------------------------------------------------------*/

#ifndef ADAPTER_COUPLING_PORO_MORTAR_H_
#define ADAPTER_COUPLING_PORO_MORTAR_H_

/*---------------------------------------------------------------------*
 | headers                                                  ager 10/15 |
 *---------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "../drt_contact/contact_poro_lagrange_strategy.H"
#include "../drt_adapter/adapter_coupling_nonlin_mortar.H"
/*---------------------------------------------------------------------*
 | forward declarations                                     ager 10/15 |
 *---------------------------------------------------------------------*/
namespace DRT
{
  class Condition;
  class Discretization;
}  // namespace DRT

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
  class CouplingPoroMortar : public CouplingNonLinMortar
  {
   public:
    /*!
    \brief Empty constructor

    */
    CouplingPoroMortar();

    /*!
    \brief Virtual destructor

    */
    virtual ~CouplingPoroMortar(){};

    virtual void EvaluatePoroMt(Teuchos::RCP<Epetra_Vector> fvel, Teuchos::RCP<Epetra_Vector> svel,
        Teuchos::RCP<Epetra_Vector> fpres, Teuchos::RCP<Epetra_Vector> sdisp,
        const Teuchos::RCP<DRT::Discretization> sdis, Teuchos::RCP<LINALG::SparseMatrix>& f,
        Teuchos::RCP<LINALG::SparseMatrix>& k_fs, Teuchos::RCP<Epetra_Vector>& frhs,
        ADAPTER::Coupling& coupfs, Teuchos::RCP<const Epetra_Map> fdofrowmap);

    void UpdatePoroMt();

    void RecoverFluidLMPoroMt(
        Teuchos::RCP<Epetra_Vector> disi, Teuchos::RCP<Epetra_Vector> veli);  // h.Willmann

    // return the used poro lagrange strategy
    Teuchos::RCP<CONTACT::PoroLagrangeStrategy> GetPoroStrategy()
    {
      if (porolagstrategy_ == Teuchos::null) dserror("GetPoroStrategy(): No strategy set!");
      return porolagstrategy_;
    };

   protected:
    /*!
    \brief Read Mortar Condition

    */
    virtual void ReadMortarCondition(Teuchos::RCP<DRT::Discretization> masterdis,
        Teuchos::RCP<DRT::Discretization> slavedis, std::vector<int> coupleddof,
        const std::string& couplingcond, Teuchos::ParameterList& input,
        std::map<int, DRT::Node*>& mastergnodes, std::map<int, DRT::Node*>& slavegnodes,
        std::map<int, Teuchos::RCP<DRT::Element>>& masterelements,
        std::map<int, Teuchos::RCP<DRT::Element>>& slaveelements);

    /*!
    \brief Add Mortar Elments

    */
    virtual void AddMortarElements(Teuchos::RCP<DRT::Discretization> masterdis,
        Teuchos::RCP<DRT::Discretization> slavedis, Teuchos::ParameterList& input,
        std::map<int, Teuchos::RCP<DRT::Element>>& masterelements,
        std::map<int, Teuchos::RCP<DRT::Element>>& slaveelements,
        Teuchos::RCP<CONTACT::CoInterface>& interface, int numcoupleddof);

    /*!
    \brief complete interface, store as internal variable
           store maps as internal variable and do parallel redist.

    */
    virtual void CompleteInterface(
        Teuchos::RCP<DRT::Discretization> masterdis, Teuchos::RCP<CONTACT::CoInterface>& interface);

    /*!
    \brief create strategy object if required

    */
    virtual void CreateStrategy(Teuchos::RCP<DRT::Discretization> masterdis,
        Teuchos::RCP<DRT::Discretization> slavedis, Teuchos::ParameterList& input,
        int numcoupleddof);

   private:
    // poro lagrange strategy
    Teuchos::RCP<CONTACT::PoroLagrangeStrategy> porolagstrategy_;

    // firstinit
    bool firstinit_;

    int slavetype_;   // 1 poro, 0 struct, -1 default
    int mastertype_;  // 1 poro, 0 struct, -1 default
  };
}  // namespace ADAPTER

#endif /* ADAPTER_COUPLING_PORO_MORTAR_H_ */
