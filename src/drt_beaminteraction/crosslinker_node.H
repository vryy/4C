/*-----------------------------------------------------------*/
/*!
\file crosslinker_node.H

\brief A class for a crosslinker node

\maintainer Jonas Eichinger

\date Oct, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef CROSSLINKER_NODE_H_
#define CROSSLINKER_NODE_H_

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_parobjectfactory.H"

// forward declaration ...
namespace MAT
{
  class Material;
  class CrosslinkerMat;
}  // namespace MAT

namespace CROSSLINKING
{
  /*!
  \brief A class for a crosslinker derived from DRT::Node
  */
  class CrosslinkerNodeType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "CrosslinkerNodeType"; };

    static CrosslinkerNodeType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static CrosslinkerNodeType instance_;
  };

  /*!
   \brief A class containing additional data from crosslinker nodes

   This class contains additional information from crosslinker nodes which are
   needed for correct crosslinking in a biopolymer network simulation. Note they are only
   available on the node's processor (ColMap). The class CrosslinkerNodeDataContainer
   must be declared before the MortarNode itself.

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
    virtual ~CrosslinkerNodeDataContainer() { return; }

    /*!
     \brief Pack this class so that it can be communicated

     This function packs the datacontainer. This is only called
     when the class has been initialized and the pointer to this
     class exists.

     */
    void Pack(DRT::PackBuffer& data) const;

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
#ifdef DEBUG
      // safety check
      if ((int)clbspots_.size() != 2) dserror("crosslinker has wrong bspot size");
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
   \brief A class for a crosslinker node derived from DRT::Node

  This class represents a single crosslinker involved in a biopolymer network simulation.
  * note:
  *
  *
  *
  *
  \author eichinger
   */

  class CrosslinkerNode : public DRT::Node
  {
   public:
    //! @name Enums and Friends

    /*!
     \brief The Discretization is a friend of MortarNode
     */
    friend class DRT::Discretization;

    //@}

    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param id     (in): A globally unique node id
    \param coords (in): vector of nodal coordinates, length 3
    \param owner  (in): Owner of this node.

    */
    CrosslinkerNode(int id, const double* coords, const int owner);

    /*!
    \brief Copy Constructor

    Makes a deep copy of a CrosslinkerNode

    */
    CrosslinkerNode(const CROSSLINKING::CrosslinkerNode& old);

    /*!
     \brief Deep copy the derived class and return pointer to it

     */
    virtual CROSSLINKING::CrosslinkerNode* Clone() const;

    /*!
    \brief Destructor

    */
    virtual ~CrosslinkerNode() { return; }

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_lib/drt_parobject.H.

     */
    virtual int UniqueParObjectId() const
    {
      return CrosslinkerNodeType::Instance().UniqueParObjectId();
    }

    /*!
     \brief Pack this class so it can be communicated

     \ref Pack and \ref Unpack are used to communicate this node

     */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
     \brief Unpack data from a char vector into this class

     \ref Pack and \ref Unpack are used to communicate this node

     */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    //! @name Access methods

    /*!
     \brief Print this MortarNode
     */
    virtual void Print(std::ostream& os) const;


    /*!
    \brief Return material of this node

    This method returns the material associated with this crosslinker node

    */
    inline Teuchos::RCP<MAT::CrosslinkerMat> GetMaterial() const
    {
      if (mat_ == Teuchos::null) dserror("No crosslinker material attached.");
      return mat_;
    }

    //  /*!
    //   \brief Initializes the data container of the node
    //
    //   With this function, the container with crosslinker binding specific quantities/information
    //   is initialized.
    //
    //   */
    //  virtual void InitializeDataContainer();

    /*!
     \brief Set material for crosslinker node

     Matnum needs to be assigned to a crosslinker type in the input file

     */
    virtual void SetMaterial(int const matnum);


    /*!
     \brief Set material for crosslinker node
     */
    virtual void SetMaterial(Teuchos::RCP<MAT::Material> material);

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
    //  Teuchos::RCP<CROSSLINKING::CrosslinkerNodeDataContainer> cldata_;

    /// this contains information that does not change during the simulation time and is
    //  the same for a subset of crosslinker, we only need one object for each subset
    Teuchos::RCP<MAT::CrosslinkerMat> mat_;


  };  // class CrosslinkerNode
}  // namespace CROSSLINKING

// << operator
std::ostream& operator<<(std::ostream& os, const CROSSLINKING::CrosslinkerNode& cnode);

#endif /* CROSSLINKING_NODE_H_ */
