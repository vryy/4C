/*---------------------------------------------------------------------*/
/*!
\file drt_meshfree_multibin.H

\brief element type class of meshfree multibin, creating the same

\maintainer Jonas Eichinger

\level 2

*/
/*---------------------------------------------------------------------*/
#ifndef DRT_MESHFREE_MESHFREEMULTIBIN_H
#define DRT_MESHFREE_MESHFREEMULTIBIN_H

#include "../drt_meshfree_discret/drt_meshfree_bin.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_node.H"
#include "../drt_particle/binning_strategy_utils.H"

namespace DRT
{
namespace MESHFREE
{

/*--------------------------------------------------------------------------*/
/*!
 * \brief The element type class of meshfree multibin, creating the same
 *
 * \author Georg Hammerl (hammerl@lnm.mw.tum.de)
 *
 * \date April, 2013
 */
/*--------------------------------------------------------------------------*/
class MeshfreeMultiBinType : public DRT::ElementType
{
public:

  //!< name of specific element type
  std::string Name() const { return "MeshfreeMultiBinType"; }

  //!< returning instance of specific element type
  static MeshfreeMultiBinType & Instance();

  //!< create parallel object of this element type
  virtual DRT::ParObject* Create( const std::vector<char> & data );

  //!< create element of this element type
  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  //!< create element of this element type
  virtual Teuchos::RCP<DRT::Element> Create( const int id,
                                             const int owner);

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ){}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ){}

private:

  //!< static instance of this element type (for self-instantiation)
  static MeshfreeMultiBinType instance_;
};

/*--------------------------------------------------------------------------*/
/*!
 * \brief A meshfree multibin holds information about associated fluid and wall
 *        elements which can be added and deleted dynamically
 *
 * \author Georg Hammerl (hammerl@lnm.mw.tum.de)
 *
 * \date April, 2013
 */
/*--------------------------------------------------------------------------*/
class MeshfreeMultiBin : public MeshfreeBin<DRT::Element>
{
public:

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  MeshfreeMultiBin(int id,   //!< (in): A globally unique bin id
                   int owner //!< (in): owner processor of the meshfree bin
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Copy Constructor
   *
   * Makes a deep copy of a meshfree multibin
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  MeshfreeMultiBin(const DRT::MESHFREE::MeshfreeMultiBin& old);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Destructor
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  virtual ~MeshfreeMultiBin();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Deep copy the derived class and return pointer to it
   *
   * This method is sort of a copy constructor for a class derived from
   * DRT::Element. It allows to copy construct the derived class without
   * knowing what it actually is using the base class Element.
   *
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  virtual DRT::Element* Clone() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return unique ParObject id
   *
   * Every class implementing ParObject needs a unique id defined at the
   * top of drt_parobject.H
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  virtual int UniqueParObjectId() const { return MeshfreeMultiBinType::Instance().UniqueParObjectId(); };

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Pack this class so it can be communicated
   *
   * \ref Pack and \ref Unpack are used to communicate this meshfree multibin
   *
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  virtual void Pack(DRT::PackBuffer& data) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Unpack data from a char vector into this class
   *
   * \ref Pack and \ref Unpack are used to communicate this meshfree multibin
   *
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  virtual void Unpack(const std::vector<char>& data);

  // return meshfree bin type instance
  virtual DRT::ElementType & ElementType() const { return MeshfreeMultiBinType::Instance(); }

  /*!
  \brief Get number of degrees of freedom of a certain node
  */
  virtual int NumDofPerNode(const DRT::Node& node) const { return 3; }

  /*------------------------------------------------------------------------*/
  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)
   *
   * Meshfree bins do not have dofs
   *///
  /*------------------------------------------------------------------------*/
  virtual int NumDofPerElement() const { return 0; }

  /*========================================================================*/
  //! @name Query methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of associated  elements to this bin
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  inline int NumAssociatedEle(BINSTRATEGY::UTILS::BinContentType bin_content) const
  {
    return associatedeleid_[bin_content].size();
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return id's of associated  elements to this bin
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  inline const int* AssociatedEleIds(BINSTRATEGY::UTILS::BinContentType bin_content) const
  {
    if (associatedeleid_[bin_content].size()) return &associatedeleid_[bin_content][0];
    else return NULL;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of ptrs to associated  elements
   *
   * \return Ptr to pointers to nodes of the element in local nodal ordering.
   *      Returns NULL if pointers to not exist.
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  inline virtual Element** AssociatedEles(BINSTRATEGY::UTILS::BinContentType bin_content)
  {
    if (associatedeleid_[bin_content].size()) return &associatedele_[bin_content][0];
    else return NULL;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Add an associated element to this bin
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  inline virtual void AddAssociatedEle(
    BINSTRATEGY::UTILS::BinContentType bin_content,   //!< (in): type of element to be added
    const int                      gid,           //!< (in): global id of element to be added
    DRT::Element*                  eleptr         //!< (in): pointer to element to be added
  )
  {
    associatedeleid_[bin_content].push_back(gid);
    associatedele_[bin_content].push_back(eleptr);
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Delete an associated element from this bin
   *
   * Searches for position of element with specified gid and deletes entry in
   * vectors associatedeleid_ and associatedele_
   *///                                                  (public) ghamm 04/13
  /*------------------------------------------------------------------------*/
  virtual void DeleteAssociatedEle(BINSTRATEGY::UTILS::BinContentType bin_content, //!< (in): type of element to be deleteted
                                   int                            gid          //!< (in): global id of element to be deleted
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Remove all associated elements of specific kind from this bin

   *///                                                  (public) ghamm 09/13
  /*------------------------------------------------------------------------*/
  void RemoveSpecificAssociatedEles(BINSTRATEGY::UTILS::BinContentType bin_content);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Remove all associated elements of all kind from this bin

   *///                                                  (public) ghamm 09/13
  /*------------------------------------------------------------------------*/
  void RemoveAllAssociatedEles();

  /*!
  \brief Build pointer vector from vector of element pointers

  The method is used to build the variable associatedele_ in this bin.
  It can be used to explicitly pass the element pointers to the bin.

  \param eles (in): Pointer to array of pointers to elements. The array of pointers
                    is implicitly expected to be of length NumAssociatedEle() and contain
                    pointers to elements in the correct element local ordering scheme.
  */
  bool BuildElePointers(BINSTRATEGY::UTILS::BinContentType bin_content, DRT::Element** eles);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this meshfree bin
   *
   * Prints basic information about this bin to std::ostream. This method would
   * usually be called by the print method of a derived class.
   *///                                                  (public) ghamm 11/12
  /*------------------------------------------------------------------------*/
  virtual void Print(std::ostream& os) const;

  /*!
  \brief Evaluate a Neumann boundary condition dummy

  An element derived from this class uses the EvaluateNeumann method to receive commands
  and parameters from some control routine in params and evaluates a Neumann boundary condition
  given in condition

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL)
  { return 0; }

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const
  {
      return DRT::Element::dis_none;
  };

private:

  //! \brief element ids contained in this bin
  std::vector<int> associatedeleid_[BINSTRATEGY::UTILS::enumsize];

  //! \brief pointers to elements contained in this bin
  std::vector<DRT::Element*> associatedele_[BINSTRATEGY::UTILS::enumsize];


}; // class MeshfreeMultiBin

} // namespace MESHFREE
} // namespace DRT

//! << operator
std::ostream& operator << (std::ostream& os, const DRT::MESHFREE::MeshfreeMultiBin& ele);

#endif  // #ifndef DRT_MESHFREE_MESHFREEMULTIBIN_H
