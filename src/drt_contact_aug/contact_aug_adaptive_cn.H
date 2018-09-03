/*----------------------------------------------------------------------------*/
/*!
\file contact_aug_adaptive_cn.H

\brief class which contains different strategies to adapt the regularization
       parameter for special scenarios

\maintainer Michael Hiermeier

\date Mar 6, 2018

\level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_ADAPTIVE_CN_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_ADAPTIVE_CN_H_

#include "contact_aug_utils.H"
#include <unordered_map>
#include "../linalg/linalg_serialdensevector.H"

class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
namespace LINALG
{
  class SerialDenseVector;
}
namespace DRT
{
  class DiscretizationInterface;
  class Node;
}  // namespace DRT
namespace STR
{
  namespace MODELEVALUATOR
  {
    class Contact;
  }
}  // namespace STR
namespace CONTACT
{
  class ParamsInterface;
  namespace AUG
  {
    class DataContainer;
    class Strategy;

    enum class AdaptiveCnType
    {
      vague,
      init_nbc
    };

    class Body
    {
      typedef std::vector<const DRT::DiscretizationInterface*> discret_set;

     public:
      /// constructor
      Body() = default;

      /// destructor
      ~Body() = default;

      inline Epetra_Map* NodeMap() { return node_map_.get(); };
      inline Epetra_Map* DofMap() { return dof_map_.get(); };

      void addNodeMapPtr(const Teuchos::RCP<Epetra_Map>& node_map) { node_map_ = node_map; }

      void addDofMapPtr(const Teuchos::RCP<Epetra_Map>& dof_map) { dof_map_ = dof_map; }

      void addContactInterfaceMaps(const plain_interface_set& interfaces, const Epetra_Comm& comm);

      void identifyAndVerifyNumDofsPerNode(
          const DRT::DiscretizationInterface& str_discret, const plain_interface_set& interfaces);

      void computeResultants(const DRT::DiscretizationInterface& str_discret,
          const plain_interface_set& interfaces, const Epetra_Vector& fext,
          const Epetra_Vector& unit_gap_force);

      double leastSquaresForceBalance(const double scale = 1.0) const;

     private:
      void resultantVector(const Epetra_Vector& vec, const Epetra_Map& node_map,
          const discret_set& discrets, LINALG::SerialDenseVector& resultant) const;

      const DRT::Node* findRowNode(int ngid, const discret_set& discrets, int& ilid) const;

      Teuchos::RCP<Epetra_Map> extractContactInterfaceMap(
          const Epetra_Map& body_map, const Epetra_Map& map, const Epetra_Map* pmap) const;

      int numDofsPerNode(const DRT::DiscretizationInterface& discret, int dofset_id,
          const int* num_dofs_per_node = NULL) const;

     private:
      Teuchos::RCP<Epetra_Map> contact_inodes_ = Teuchos::null;

      Teuchos::RCP<Epetra_Map> contact_idofs_ = Teuchos::null;

      std::set<int> interface_lids_;

      /// node map of this body
      Teuchos::RCP<Epetra_Map> node_map_ = Teuchos::null;

      /// dof map of this body
      Teuchos::RCP<Epetra_Map> dof_map_ = Teuchos::null;

      /// resultant of all external forces acting on this body
      LINALG::SerialDenseVector rfext_;

      /// resultant of all contac forces acting on this body
      LINALG::SerialDenseVector rcontact_;

      int num_dofs_per_node_ = -1;
    };

    /*--------------------------------------------------------------------------*/
    class AdaptiveCn
    {
      typedef std::vector<Teuchos::RCP<Body>> nbc_body_set;

     public:
      AdaptiveCn(Strategy& strat, const plain_interface_set& interfaces, DataContainer& data);

      AdaptiveCn(const AdaptiveCn& adcn) = delete;

      ~AdaptiveCn() = default;

      void Execute(const AdaptiveCnType adaptive_cn_type, CONTACT::ParamsInterface& cparams,
          std::ostream* os);

     private:
      class InitNbc;

     private:
      Strategy& strategy_;
      const plain_interface_set& interfaces_;
      DataContainer& data_;

      AdaptiveCnType type_ = AdaptiveCnType::vague;
    };

    /*--------------------------------------------------------------------------*/
    class AdaptiveCn::InitNbc
    {
     public:
      explicit InitNbc(AdaptiveCn& parent) : parent_(parent){};

      void operator()(CONTACT::ParamsInterface& cparams, std::ostream* os = NULL);

     private:
      Teuchos::RCP<Epetra_Vector> UnitGapForce() const;

      const STR::MODELEVALUATOR::Contact& CModel(CONTACT::ParamsInterface& cparams) const;

      void SplitStructureIntoDistinctNeumannBodies(nbc_body_set& neumann_bodies,
          const Epetra_Vector& fext, const DRT::DiscretizationInterface& full_discret) const;

      void PrintUpdate(std::ostream* os_ptr) const;

     private:
      AdaptiveCn& parent_;
      double dt_ = 0.0;
      double cn_init_ = 0.0;
    };
  }  // namespace AUG
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_ADAPTIVE_CN_H_ */
