/*-----------------------------------------------------------*/
/*! \file

\brief A class for a crosslinker node


\date Oct, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_CROSSLINKER_NODE_HPP
#define FOUR_C_BEAMINTERACTION_CROSSLINKER_NODE_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_fem_general_node.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace Core::Mat
{
  class Material;
}  // namespace Core::Mat

namespace Mat
{
  class CrosslinkerMat;
}

namespace CrossLinking
{
  /*!
  \brief A class for a crosslinker derived from Core::Nodes::Node
  */
  class CrosslinkerNodeType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "CrosslinkerNodeType"; };

    static CrosslinkerNodeType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static CrosslinkerNodeType instance_;
  };

  /*!
   \brief A class containing additional data from crosslinker nodes

   This class contains additional information from crosslinker nodes which are
   needed for correct crosslinking in a biopolymer network simulation. Note they are only
   available on the node's processor (ColMap). The class CrosslinkerNodeDataContainer
   must be declared before the Mortar::Node itself.

   \author eichinger
   */
  class CrosslinkerNodeDataContainer
  {
   public:
    //! @name Constructors and destructors and related methods

    /*!
     \brief Standard Constructor

     */
    CrosslinkerNodeDataContainer();

    /*!
     \brief Destructor

     */
    virtual ~CrosslinkerNodeDataContainer() = default;
    /*!
     \brief Pack this class so that it can be communicated

     This function packs the datacontainer. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    void Pack(Core::Communication::PackBuffer& data) const;

    /*!
     \brief Unpack data from a vector into this class

     This function unpacks the datacontainer. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

    //@}

    //! @name Access methods

    /*!
     \brief Get current binding spot status of linker
     */
    const std::vector<std::pair<int, int>>& GetClBSpotStatus()
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      // safety check
      if ((int)clbspots_.size() != 2) FOUR_C_THROW("crosslinker has wrong bspot size");
#endif
      return clbspots_;
    }

    /*!
     \brief Get
     */
    void SetClBSpotStatus(std::vector<std::pair<int, int>> clbspots)
    {
      clbspots_ = clbspots;
      return;
    }

    /*!
    \brief Get current number of bonds of crosslinker
    */
    const int& GetNumberOfBonds() { return numbond_; }

    /*!
    \brief Set current number of bonds of crosslinker
    */
    void SetNumberOfBonds(int numbond)
    {
      numbond_ = numbond;
      return;
    }

    //@}

   protected:
    // don't want = operator and cctor
    CrosslinkerNodeDataContainer operator=(const CrosslinkerNodeDataContainer& old);
    CrosslinkerNodeDataContainer(const CrosslinkerNodeDataContainer& old);

    /// gid of element to local number of bspot, [0] and [1] first and second bspot of cl
    std::vector<std::pair<int, int>> clbspots_;
    /// number of active bonds
    int numbond_;
  };
  // class CrosslinkerNodeDataContainer

  /*!
   \brief A class for a crosslinker node derived from Core::Nodes::Node

  This class represents a single crosslinker involved in a biopolymer network simulation.
  * note:
  *
  *
  *
  *
  \author eichinger
   */

  class CrosslinkerNode : public Core::Nodes::Node
  {
   public:
    //! @name Enums and Friends

    /*!
     \brief The discretization is a friend of Mortar::Node
     */
    friend class Discret::Discretization;

    //@}

    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param id     (in): A globally unique node id
    \param coords (in): vector of nodal coordinates, length 3
    \param owner  (in): Owner of this node.

    */
    CrosslinkerNode(int id, const std::vector<double>& coords, const int owner);

    /*!
    \brief Copy Constructor

    Makes a deep copy of a CrosslinkerNode

    */
    CrosslinkerNode(const CrossLinking::CrosslinkerNode& old);

    /*!
     \brief Deep copy the derived class and return pointer to it

     */
    CrossLinking::CrosslinkerNode* Clone() const override;



    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of lib/parobject.H.

     */
    int UniqueParObjectId() const override
    {
      return CrosslinkerNodeType::Instance().UniqueParObjectId();
    }

    /*!
     \brief Pack this class so it can be communicated

     \ref Pack and \ref Unpack are used to communicate this node

     */
    void Pack(Core::Communication::PackBuffer& data) const override;

    /*!
     \brief Unpack data from a char vector into this class

     \ref Pack and \ref Unpack are used to communicate this node

     */
    void Unpack(const std::vector<char>& data) override;

    //@}

    //! @name Access methods

    /*!
     \brief Print this Mortar::Node
     */
    void Print(std::ostream& os) const override;


    /*!
    \brief Return material of this node

    This method returns the material associated with this crosslinker node

    */
    inline Teuchos::RCP<Mat::CrosslinkerMat> GetMaterial() const
    {
      if (mat_ == Teuchos::null) FOUR_C_THROW("No crosslinker material attached.");
      return mat_;
    }

    //  /*!
    //   \brief Initializes the data container of the node
    //
    //   With this function, the container with crosslinker binding specific quantities/information
    //   is initialized.
    //
    //   */
    //  virtual void initialize_data_container();

    /*!
     \brief Set material for crosslinker node

     Matnum needs to be assigned to a crosslinker type in the input file

     */
    virtual void SetMaterial(int const matnum);


    /*!
     \brief Set material for crosslinker node
     */
    virtual void SetMaterial(Teuchos::RCP<Core::Mat::Material> material);

    //  /*!
    //   \brief Resets the data container of the node
    //
    //   With this function, the container with crosslinker binding specific quantities/information
    //   is deleted / reset to Teuchos::null pointer
    //
    //   */
    //  virtual void ResetDataContainer();

    //@}


   protected:
    //  /// information of crosslinker binding status, this is different for each crosslinker
    //  //  and may change each time step
    //  Teuchos::RCP<CrossLinking::CrosslinkerNodeDataContainer> cldata_;

    /// this contains information that does not change during the simulation time and is
    //  the same for a subset of crosslinker, we only need one object for each subset
    Teuchos::RCP<Mat::CrosslinkerMat> mat_;


  };  // class CrosslinkerNode
}  // namespace CrossLinking

// << operator
std::ostream& operator<<(std::ostream& os, const CrossLinking::CrosslinkerNode& crosslinker_node);

FOUR_C_NAMESPACE_CLOSE

#endif
