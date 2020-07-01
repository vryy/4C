/*----------------------------------------------------------------------*/
/*! \file

\brief


\level 2

*---------------------------------------------------------------------------*/
#ifndef DRT_MESHFREE_MESHFREEBIN_H
#define DRT_MESHFREE_MESHFREEBIN_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"

namespace DRT
{
  namespace MESHFREE
  {
    /*--------------------------------------------------------------------------*/
    /*!
     * \brief A meshfree bin adds the possibility to add and delete nodes
     *        from elements dynamically
     *
     * Rather than deriving this class directly from DRT::Element, we set the base
     * element as a template argument in order to allow for both standard elements
     * and face elements, i.e., elements that have a parent element. In
     * particular, it allows for a clean interface to the meshfree types with
     * single inheritance for making the meshfree boundary elements aware of the
     * parent relations.
     *
     *
     * \date November, 2012
     */
    /*--------------------------------------------------------------------------*/
    template <typename ELEMENT>
    class MeshfreeBin : public ELEMENT
    {
     public:
      /*========================================================================*/
      //! @name Constructors and destructors and related methods
      /*========================================================================*/

      /*------------------------------------------------------------------------*/
      /*!
       * \brief Standard Constructor
       *///                                                  (public) ghamm 11/12
      /*------------------------------------------------------------------------*/
      MeshfreeBin(int id,  //!< (in): A globally unique bin id
          int owner        //!< (in): owner processor of the meshfree bin
      );

      /*------------------------------------------------------------------------*/
      /*!
       * \brief Copy Constructor
       *
       * Makes a deep copy of a meshfree bin
       *///                                                  (public) ghamm 11/12
      /*------------------------------------------------------------------------*/
      MeshfreeBin(const DRT::MESHFREE::MeshfreeBin<ELEMENT>& old);

      /*------------------------------------------------------------------------*/
      /*!
       * \brief Destructor
       *///                                                  (public) ghamm 11/12
      /*------------------------------------------------------------------------*/
      virtual ~MeshfreeBin();

      /*------------------------------------------------------------------------*/
      /*!
       * \brief Deep copy the derived class and return pointer to it
       *
       * This method is sort of a copy constructor for a class derived from
       * DRT::Element. It allows to copy construct the derived class without
       * knowing what it actually is using the base class Element.
       *
       *///                                                  (public) ghamm 11/12
      /*------------------------------------------------------------------------*/
      virtual DRT::Element* Clone() const = 0;

      /*------------------------------------------------------------------------*/
      /*!
       * \brief Return unique ParObject id
       *
       * Every class implementing ParObject needs a unique id defined at the
       * top of drt_parobject.H
       *///                                                  (public) ghamm 11/12
      /*------------------------------------------------------------------------*/
      virtual int UniqueParObjectId() const = 0;

      /*========================================================================*/
      //! @name Query methods
      /*========================================================================*/

      /*------------------------------------------------------------------------*/
      /*!
       * \brief Add a node to the meshfree bin
       *
       * \warning It is your own responsibility to make sure that there will not
       *          be any double entries! (This would be disastrous!)
       *
       * Adds entry at the end of nodeid_ and node_ pointers
       *///                                                  (public) ghamm 11/12
      /*------------------------------------------------------------------------*/
      virtual void AddNode(DRT::Node* nodeptr  //!< (in): pointer to node to be added
      )
      {
        ELEMENT::nodeid_.push_back(nodeptr->Id());
        ELEMENT::node_.push_back(nodeptr);
        return;
      }

      /*------------------------------------------------------------------------*/
      /*!
       * \brief Add a node to the meshfree bin
       *
       * \warning It is your own responsibility to make sure that there will not
       *          be any double entries! (This would be disastrous!)
       *
       * Adds entry at the end of nodeid_ and node_ pointers
       *///                                                  (public) ghamm 11/12
      /*------------------------------------------------------------------------*/
      virtual void AddNode(const int gid,  //!< (in): global id of node to be added
          DRT::Node* nodeptr               //!< (in): pointer to node to be added
      )
      {
        ELEMENT::nodeid_.push_back(gid);
        ELEMENT::node_.push_back(nodeptr);
        return;
      }

      /*------------------------------------------------------------------------*/
      /*!
       * \brief Delete a node from the meshfree bin
       *
       * Searches for position of node with specified gid and deletes entry in
       * vectors node_ and nodeid_
       *///                                                  (public) ghamm 11/12
      /*------------------------------------------------------------------------*/
      virtual void DeleteNode(int gid  //!< (in): global id of node to be deleted
      );

      /*------------------------------------------------------------------------*/
      /*!
       * \brief Delete all nodes from the meshfree bin
       *///                                                  (public) ghamm 11/12
      /*------------------------------------------------------------------------*/
      virtual inline void DeleteNodes()
      {
        ELEMENT::nodeid_.clear();
        ELEMENT::node_.clear();
      }

     private:
    };  // class MeshfreeBin

  }  // namespace MESHFREE
}  // namespace DRT

#endif  // #ifndef DRT_MESHFREE_MESHFREEBIN_H
