/*----------------------------------------------------------------------*/
/*! \file
\brief A class for a rough contact node
\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef BACI_CONTACT_ROUGH_NODE_HPP
#define BACI_CONTACT_ROUGH_NODE_HPP

#include "baci_config.hpp"

#include "baci_contact_node.hpp"

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  class RoughNodeType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const final { return "RoughNodeType"; }

    static RoughNodeType& Instance() { return instance_; };

    virtual CORE::COMM::ParObject* Create(const std::vector<char>& data);

   private:
    static RoughNodeType instance_;
  };


  class RoughNode : public Node
  {
   public:
    //! @name Enums and Friends

    /*!
     \brief The Discretization is a friend of RoughNode
     */
    friend class DRT::Discretization;

    //@}

    //! @name Constructors and destructors and related methods

    /*!
     \brief Standard Constructor

     \param id     (in): A globally unique node id
     \param coords (in): vector of nodal coordinates, length 3
     \param owner  (in): Owner of this node.
     \param dofs   (in): list of global degrees of freedom
     \param isslave(in): flag indicating whether node is slave or master
     \param initactive (in): flag indicating whether initially set to active
     \param hurstexponentfunction (in): function id for roughness parameter
     \\ add remaining params here

     */
    RoughNode(int id, const std::vector<double>& coords, const int owner,
        const std::vector<int>& dofs, const bool isslave, const bool initactive,
        const int hurstexponentfunction, int initialtopologystddeviationfunction, int resolution,
        int randomtopologyflag, int randomseedflag, int randomgeneratorseed);

    CORE::LINALG::SerialDenseMatrix* GetTopology() { return &topology_; };
    double GetMaxTopologyHeight() { return maxTopologyHeight_; };

   protected:
    double hurstExponent_ = 0;
    double initialTopologyStdDeviation_ = 0;
    CORE::LINALG::SerialDenseMatrix topology_;
    double maxTopologyHeight_;
  };
}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif  // CONTACT_ROUGH_NODE_HPP
