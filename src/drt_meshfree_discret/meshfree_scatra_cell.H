/*!---------------------------------------------------------------------------
\file meshfree_scatra_cell.H

\brief scatra cell for meshfree discretisations

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#ifndef MESHFREE_SCATRA_CELL_H
#define MESHFREE_SCATRA_CELL_H

#include "drt_meshfree_cell.H"

#include "../drt_lib/drt_elementtype.H"

#include "../drt_inpar/inpar_scatra.H"

namespace DRT
{

namespace ELEMENTS
{

/*==========================================================================*/
// forward declarations
/*==========================================================================*/
// none

/*--------------------------------------------------------------------------*/
/*!
 * \brief The element type class of meshfree transport cells, creating the
 *        same
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date March, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeTransportType : public DRT::ElementType
{
public:

  //!< name of specific element type
  std::string Name() const { return "MeshfreeTransportType"; }

  //!< returning instance of specific element type
  static MeshfreeTransportType & Instance();

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

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

  //!< pre-evaluation
  virtual void PreEvaluate(
    DRT::Discretization&                 dis,
    Teuchos::ParameterList&              p,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector>          systemvector1,
    Teuchos::RCP<Epetra_Vector>          systemvector2,
    Teuchos::RCP<Epetra_Vector>          systemvector3);

private:

  //!< static instance of this element type (for self-instantiation)
  static MeshfreeTransportType instance_;
};

/*--------------------------------------------------------------------------*/
/*!
 * \brief A C++ wrapper for the MeshfreeTransport cell
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date February, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeTransport : public DRT::MESHFREE::Cell<DRT::Element>
{
public:

  //! @name Friends
  friend class MeshfreeTransportBoundary;

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  MeshfreeTransport(
    int id,    //!< A unique global id
    int owner  //!< processor id who owns a certain instance of this class
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Copy Constructor
   *
   * Makes a deep copy of a DRT::Element
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  MeshfreeTransport(const MeshfreeTransport& old);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Deep copy this instance of MeshfreeTransport and return pointer to
   *        the copy
   *
   *   The Clone() method is used from the virtual base class Element in cases
   *   where the type of the derived class is unknown and a copy-ctor is
   *   needed
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  DRT::Element* Clone() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Destructor
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual ~MeshfreeTransport();

  /*========================================================================*/
  //! @name access methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Set cell material
   *
   *   Material numbers are read from the input file. The element stores a
   *   corresponding material object. These material objects can be anything
   *   from very simple (just a little calculation) to highly sophisticated
   *   with history data. The material is packed and unpacked along with its
   *   element.
   *
   *  \note reimplementation of this method, due to initialising
   *        numdofpernode_, since the material is known now.
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual void SetMaterial(int matnum //!< material number from input file
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get shape type of cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual DiscretizationType Shape() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief set discretization type of cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual void SetDisType(DiscretizationType shape //!< shape determining distype
    )
  {
    distype_ = shape;
    return;
  };

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of lines of this cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int NumLine() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of surfaces of this cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int NumSurface() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of volumes of this cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int NumVolume() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the lines of this cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the surfaces of this cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the volume of this cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Volumes();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return unique ParObject id
   *
   * Every class implementing ParObject needs a unique id defined at the top
   * of this file.
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int UniqueParObjectId() const
  {
    return MeshfreeTransportType::Instance().UniqueParObjectId();
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Pack this class so it can be communicated
   *
   * \ref Pack and \ref Unpack are used to communicate this cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual void Pack(DRT::PackBuffer& data) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Unpack data from a char vector into this class
   *
   * \ref Pack and \ref Unpack are used to communicate this cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual void Unpack(const std::vector<char>& data);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get number of degrees of freedom of a certain node
   *        (implements pure virtual DRT::Element)
   *
   * The element decides how many degrees of freedom its nodes must have. As
   * this may vary along a simulation, the element can redecide the number of
   * degrees of freedom per node along the way for each of it's nodes
   * separately.
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    dsassert(numdofpernode_!=-1,"NumDofPerNode is -1");
    return numdofpernode_;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get number of degrees of freedom per element
   *        (implements pure virtual DRT::Element)
   *
   * The element decides how many element degrees of freedom it has. It can
   * redecide along the way of a simulation.
   *
   * \note Element degrees of freedom mentioned here are dofs that are visible
   * at the level of the total system of equations. Purely internal element
   * dofs that are condensed internally should NOT be considered.
   *
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int NumDofPerElement() const
  {
    return 0;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this element
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  void Print(std::ostream& os) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief TODO
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual DRT::ElementType & ElementType() const
  {
    return MeshfreeTransportType::Instance();
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Query names of element data to be visualized using BINIO
   *
   * The element fills the provided map with key names of visualization data
   * the element wants to visualize AT THE CENTER of the element geometry. The
   * values is supposed to be dimension of the data to be visualized. It can
   * either be 1 (scalar), 3 (vector), 6 (sym. tensor) or 9 (nonsym. tensor)
   *
   * Example:
   * \code
   *   // Name of data is 'Owner', dimension is 1 (scalar value)
   *   names.insert(std::pair<std::string,int>("Owner",1));
   *   // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
   *   names.insert(std::pair<std::string,int>("StressesXYZ",6));
   * \endcode
   *
   * \param names (out): On return, the derived class has filled names with
   *                     key names of data it wants to visualize and with int
   *                     dimensions of that data.
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual void VisNames(std::map<std::string,int>& names);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Query data to be visualized using BINIO of a given name
   *
   * The method is supposed to call this base method to visualize the owner of
   * the element.
   * If the derived method recognizes a supported data name, it shall fill it
   * with corresponding data.
   * If it does NOT recognizes the name, it shall do nothing.
   *
   * \warning The method must not change size of data
   *
   * \param name (in):   Name of data that is currently processed for visualization
   * \param data (out):  data to be filled by element if element recognizes the name
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual bool VisData(const std::string& name, std::vector<double>& data);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Add values to the data_ container
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual void AddToData(const std::string& name, std::vector<double> data)
  {
    data_.Add(name,data);
  };

  /*========================================================================*/
  //! @name Input and Creation
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Read input for this element
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);

  /*========================================================================*/
  //! @name Evaluation
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate an element, i.e. call the implementation to evaluate MeshfreeTransport
   * element stiffness, mass, internal forces or evaluate errors, statistics or updates
   * etc. directly.
   *
   * Following implementations of the element are allowed:
   * <pre>
   *   o Evaluation of condif system matrix and residual for the One-Step-Theta or BDF2
   *     method
   *
   *   o Evaluation of condif system matrix and residual for the stationary fluid solver
   * </pre>
   *
   * \param params (in/out)    : ParameterList for communication between
   *                             control routine and elements
   * \param discretization (in): A reference to the underlying discretization
   * \param lm (in)            : location vector of this element
   * \param elemat1 (out)      : matrix to be filled by element. If NULL on
   *                             input, the controling method does not epxect
   *                             the element to fill this matrix.
   * \param elemat2 (out)      : matrix to be filled by element. If NULL on
   *                             input, the controling method does not epxect
   *                             the element to fill this matrix.
   * \param elevec1 (out)      : vector to be filled by element. If NULL on input,
   *                             input, the controling method does not epxect
   *                             the element to fill this vector.
   * \param elevec2 (out)      : vector to be filled by element. If NULL on input,
   *                             input, the controling method does not epxect
   *                             the element to fill this vector.
   * \param elevec3 (out)      : vector to be filled by element. If NULL on input,
   *                             input, the controling method does not epxect
   *                             the element to fill this vector.
   * \return 0 if successful, negative otherwise
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  int Evaluate(Teuchos::ParameterList&   params,
               DRT::Discretization&      discretization,
               std::vector<int>&         lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);


  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate a Neumann boundary condition
   *
   * this method evaluates a surfaces Neumann condition on this element
   *
   * \param params (in/out)    : ParameterList for communication between control routine
   *                             and elements
   * \param discretization (in): A reference to the underlying discretization
   * \param condition (in)     : The condition to be evaluated
   * \param lm (in)            : location vector of this element
   * \param elevec1 (out)      : vector to be filled by element. If NULL on input,
   *
   * \return 0 if successful, negative otherwise
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);


private:

  //! don't want = operator
  MeshfreeTransport& operator = (const MeshfreeTransport& old);

  //! container containing all kind of element thisandthat
  DRT::Container data_;

  //! (storage neccessary because we dont know the material in the post filters anymore)
  int numdofpernode_;

  //! the element discretization type
  DRT::Element::DiscretizationType distype_;

  //! implementation type (physics)
  INPAR::SCATRA::ImplType impltype_;
}; // class MeshfreeTransport


/*--------------------------------------------------------------------------*/
/*!
 * \brief The element type class of meshfree transport boundary cells
 *
 * TODO: Are boundary type classes really necessary? Create doesn't create
 *       anything!
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date March, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeTransportBoundaryType : public DRT::ElementType
{
public:

  //!< name of specific parallel object type
  std::string Name() const { return "MeshfreeTransportBoundaryType"; }

  //!< returning instance of specific parallel object type
  static MeshfreeTransportBoundaryType & Instance();

  //!< create parallel object of this parallel object type
  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  //!< static instance of this parallel object type (for self-instantiation)
  static MeshfreeTransportBoundaryType instance_;
};


/*--------------------------------------------------------------------------*/
/*!
 * \brief An element representing a boundary cell of a MeshfreeTransport cell
 *
 * \note This is a pure boundary condition cell. It's only purpose is to
 *       evaluate certain boundary conditions that might be adjacent to a
 *       parent MeshfreeTransport cell.
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date February, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeTransportBoundary : public DRT::MESHFREE::Cell<DRT::FaceElement>
{
public:

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  MeshfreeTransportBoundary(
    int id,                    //!< A unique global id
    int owner,                 //!< Processor owning this surface
    int npoint,                //!< Number of nodes attached to this element
    const int* pointids,       //!< global ids of nodes attached to this element
    DRT::Node** points,        //!< the discretizations map of points to build ptrs to points
    MeshfreeTransport* parent, //!< The parent fluid element of this surface
    const int lbeleid          //!< the local boundary element number of this surface w.r.t. the parent element
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Copy Constructor
   *
   * Makes a deep copy of a MeshfreeTransportBoundary
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  MeshfreeTransportBoundary(const MeshfreeTransportBoundary& old);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Deep copy the derived class and return pointer to it
   *
   * This method is sort of a copy constructor for a class derived from
   * DRT::Element. It allows to copy construct the derived class without
   * knowing what it actually is using the base class DRT::Element.
   *
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  DRT::Element* Clone() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Destructor
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual ~MeshfreeTransportBoundary();

  /*========================================================================*/
  //! @name Query methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get shape type of element
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual DiscretizationType Shape() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of lines of boundary cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int NumLine() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of surfaces of boundary cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int NumSurface() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the lines of this cell
   *        Dummy - not provided.
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the surfaces of this cell
   *        Dummy - not provided.
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Pack dummy class - not supported
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual void Pack(DRT::PackBuffer& data) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Unpack dummy class - not supported
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual void Unpack(const std::vector<char>& data);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get number of degrees of freedom of a certain node
   *        (implements pure virtual DRT::Element)
   *
   * The element decides how many degrees of freedom its nodes must have.
   * As this may vary along a simulation, the element can redecide the
   * number of degrees of freedom per node along the way for each of it's nodes
   * separately.
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    return ParentElement()->NumDofPerNode(node);
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return a pointer to the parent cell of this boundary cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  /*
  virtual DRT::ELEMENTS::MeshfreeTransport* ParentElement()
  {
    return parent_;
  }
  */

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return a boundary element id of this boundary cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  /*
  virtual int BeleNumber()
  {
    return lbeleid_;
  }
  */

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return unique ParObject id
   *
   * every class implementing ParObject needs a unique id defined at the
   * top of the drt_parobject.H file.
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int UniqueParObjectId() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this cell
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  void Print(std::ostream& os) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief TODO
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual DRT::ElementType & ElementType() const
  {
    return MeshfreeTransportBoundaryType::Instance();
  }

  /*========================================================================*/
  //! @name Evaluate methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate an element
   *
   * Evaluate MeshfreeTransport element stiffness, mass, internal forces etc
   *
   * \param params (in/out): ParameterList for communication between control routine
   *                         and elements
   * \param discretization (in): A reference to the underlying discretization
   * \param lm (in):         location vector of this element
   * \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
   *                         the controling method does not epxect the element to fill
   *                         this matrix.
   * \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
   *                         the controling method does not epxect the element to fill
   *                         this matrix.
   * \param elevec1 (out)  : vector to be filled by element. If NULL on input,
   *                         the controlling method does not epxect the element
   *                         to fill this vector
   * \param elevec2 (out)  : vector to be filled by element. If NULL on input,
   *                         the controlling method does not epxect the element
   *                         to fill this vector
   * \param elevec3 (out)  : vector to be filled by element. If NULL on input,
   *                         the controlling method does not epxect the element
   *                         to fill this vector
   * \return 0 if successful, negative otherwise
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  int Evaluate(Teuchos::ParameterList&   params,
               DRT::Discretization&      discretization,
               std::vector<int>&         lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate a Neumann boundary condition
   *
   * this method evaluates a surface Neumann condition on the MeshfreeTransport element
   *
   * \param params (in/out)    : ParameterList for communication between control routine
   *                             and elements
   * \param discretization (in): A reference to the underlying discretization
   * \param condition (in)     : The condition to be evaluated
   * \param lm (in)            : location vector of this element
   * \param elevec1 (out)      : vector to be filled by element. If NULL on input,
   *
   * \return 0 if successful, negative otherwise
   *///                                                   (public) nis Feb12
  /*------------------------------------------------------------------------*/
  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*------------------------------------------------------------------------*/

private:

  /*------------------------------------------------------------------------*/
  // don't want = operator
  /*------------------------------------------------------------------------*/
  MeshfreeTransportBoundary& operator = (const MeshfreeTransportBoundary& old);

}; // class MeshfreeTransportBoundary

} // namespace ELEMENTS
} // namespace DRT


#endif // end MESHFREE_SCATRA_CELL_H

