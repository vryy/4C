/*----------------------------------------------------------------------*/
/*! \file
\brief A class for a rough contact node
\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_ROUGH_NODE_HPP
#define FOUR_C_CONTACT_ROUGH_NODE_HPP

#include "4C_config.hpp"

#include "4C_contact_node.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  class RoughNodeType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const final { return "RoughNodeType"; }

    static RoughNodeType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static RoughNodeType instance_;
  };


  class RoughNode : public Node
  {
   public:
    //! @name Enums and Friends

    /*!
     \brief The discretization is a friend of RoughNode
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
     \param hurstexponentfunction (in): Function for Hurst exponent of the surface
     \param initialtopologystddeviationfunction (in): function for topology standard deviation
     \param resolution (in): resolution of the surface
     \param randomtopologyflag (in): Use random midpoint generator if true
     \param randomseedflag (in): Use random seed for the random midpoint generator
     \param randomgeneratorseed (in): Seed for the random midpoint generator
     \\ add remaining params here

     */
    RoughNode(int id, const std::vector<double>& coords, const int owner,
        const std::vector<int>& dofs, const bool isslave, const bool initactive,
        const int hurstexponentfunction, int initialtopologystddeviationfunction, int resolution,
        bool randomtopologyflag, bool randomseedflag, int randomgeneratorseed);

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of lib/parobject.H.

     */
    int UniqueParObjectId() const override { return RoughNodeType::Instance().UniqueParObjectId(); }

    /*!
     \brief Pack this class so it can be communicated

     \ref Pack and \ref Unpack are used to communicate this node

    */
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /*!
     \brief Unpack data from a char vector into this class

     \ref Pack and \ref Unpack are used to communicate this node

     */
    void Unpack(const std::vector<char>& data) override;

    // //! @name Access methods

    CORE::LINALG::SerialDenseMatrix* GetTopology() { return &topology_; };
    double get_max_topology_height() { return maxTopologyHeight_; };

   protected:
    int hurstexponentfunction_;
    int initialtopologystddeviationfunction_;
    int resolution_;
    bool randomtopologyflag_;
    bool randomseedflag_;
    int randomgeneratorseed_;

    double hurstExponent_ = 0;
    double initialTopologyStdDeviation_ = 0;
    CORE::LINALG::SerialDenseMatrix topology_;
    double maxTopologyHeight_;
  };
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
