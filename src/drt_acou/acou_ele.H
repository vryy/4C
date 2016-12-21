/*!----------------------------------------------------------------------
\file acou_ele.H

\brief A C++ wrapper for the acoustic element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ACOU_ELE_H
#define ACOU_ELE_H

#include "Teuchos_RCP.hpp"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_fem_general/drt_utils_polynomial.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{
class AcouIntFace;

class AcouType : public DRT::ElementType
{
public:

  std::string Name() const { return "AcouType"; }

  static AcouType & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

  /// pre-evaluation
  virtual void PreEvaluate(DRT::Discretization&                dis,
                          Teuchos::ParameterList&              p,
                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                          Teuchos::RCP<Epetra_Vector>          systemvector1,
                          Teuchos::RCP<Epetra_Vector>          systemvector2,
                          Teuchos::RCP<Epetra_Vector>          systemvector3);

private:

  static AcouType instance_;

};

/*!
\brief A C++ wrapper for the acoustic element
*/
class Acou : public DRT::Element
{
public:

  //@}
  //! @name constructors and destructors and related methods

  /*!
  \brief standard constructor
  */
  Acou(
    int id,    ///< A unique global id
    int owner  ///< ???
    );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Acou(const Acou& old);

  /*!
  \brief Deep copy this instance and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Destructor
  */
  virtual ~Acou();

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const
  {
    return distype_;
  };

  /*!
  \brief set discretization type of element
  */
  virtual void SetDisType(DRT::Element::DiscretizationType shape)
  {
    distype_ = shape;
    return;
  };

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const
  {
    return DRT::UTILS::getNumberOfElementLines(distype_);
  }

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const
  {
    return DRT::UTILS::getNumberOfElementSurfaces(distype_);
  }

  /*!
  \brief Return number of volumes of this element (always 1)
  */
  virtual int NumVolume() const
  {
    return DRT::UTILS::getNumberOfElementVolumes(distype_);
  }

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of Teuchos::RCPs to the surfaces of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*!
  \brief Get vector of Teuchos::RCPs to the volume of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Volumes();

  /*!
  \brief Get Teuchos::RCP to the internal face adjacent to this element as master element and the parent_slave element
  */
  virtual Teuchos::RCP<DRT::Element> CreateFaceElement(
                                                DRT::Element* parent_slave,          //!< parent slave fluid3 element
                                                int nnode,                           //!< number of surface nodes
                                                const int* nodeids,                  //!< node ids of surface element
                                                DRT::Node** nodes,                   //!< nodes of surface element
                                                const int lsurface_master,           //!< local surface number w.r.t master parent element
                                                const int lsurface_slave,            //!< local surface number w.r.t slave parent element
                                                const std::vector<int>& localtrafomap//! local trafo map
      );

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return AcouType::Instance().UniqueParObjectId(); }
  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element
  */
  virtual void Unpack(const std::vector<char>& data);

  //@}

  //! @name Geometry related methods

  //@}

  //! @name Acess methods

  /*!
  \brief Get number of degrees of freedom per node

  HDG element: No dofs are associated with nodes
  */
  virtual int NumDofPerNode(const DRT::Node&) const { return 0; }

  /*!
  \brief Get number of degrees of freedom per face
  */
  virtual int NumDofPerFace(const unsigned face) const
  {
    return NumDofPerComponent(face); // seems redundant, but this is different for fluid_ele_hdg and acou_sol_ele!
  }

  /*!
  \brief Get number of dofs per component per face
  */
  virtual int NumDofPerComponent(const unsigned face) const
  {
    return DRT::UTILS::getBasisSize(DRT::UTILS::getEleFaceShapeType(distype_), (this->Faces()[face])->Degree(), completepol_);
  }

  virtual int NumDofPerElementAuxiliary() const
  {
    return (DRT::UTILS::getDimension(distype_) + 1) *
            DRT::UTILS::getBasisSize(distype_, degree_, completepol_);
  }

  virtual int NumDofPerElement() const { return 0; }

  /*!
   \brief Returns the degree of the element
   */
  virtual int Degree() const
  {
    return degree_;
  }

  /*!
   \brief Allows to set the element degree
   */
  void SetDegree(int degree)
  {
    degree_ = degree;
    return;
  }

  /*!
   \brief Returns the degree of the element
   */
  int UsesCompletePolynomialSpace() const
  {
    return completepol_;
  }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return AcouType::Instance(); }

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string&          eletype,
                   const std::string&          distype,
                   DRT::INPUT::LineDefinition* linedef);

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element, that is, call the element routines to evaluate acou
  element matrices and vectors or evaluate errors, statistics or updates etc. directly.

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \return 0 if successful, negative otherwise
  */
  int Evaluate(Teuchos::ParameterList&   params,
               DRT::Discretization&      discretization,
               std::vector<int>&         lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);


  /*!
  \brief Evaluate Neumann boundary condition

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): reference to the underlying discretization
  \param condition (in)     : condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  //@}

  // stores the nodal trace values in case of p-adaptivity
  Epetra_SerialDenseVector elenodeTrace_;

  // state vectors which do not have to be stored as global distributed vector, since
  // we are in a discontinuous setting here
  //! @name interior stored states

  /// local values from interior solution vector (gradient, velocity) at n+1
  Epetra_SerialDenseVector eleinteriorVelnp_;
  Epetra_SerialDenseVector eleinteriorPressnp_;

  Epetra_SerialDenseVector eleADERimprovedGrad_;
  Epetra_SerialDenseVector eleADERimprovedDiv_;

  /// local values for acoustical parameter optimization
  std::vector<Epetra_SerialDenseVector> eleintPressForward_;
  std::vector<Epetra_SerialDenseVector> eleintVelForward_;
  //@}

  double GetDensityGradient() const { return grad_rho_; }
  double GetSoSGradient()     const { return grad_c_;   }

  void AddToDensityGradient(double toadd) { grad_rho_ += toadd; return; }
  void AddToSoSGradient(double toadd)     { grad_c_   += toadd; return; }

  void InitDensityAndSoSGradient() { grad_rho_=0.0; grad_c_=0.0; return; }

protected:

  //! discretization type
  DRT::Element::DiscretizationType distype_;
  DRT::Element::DiscretizationType facedis_;

private:

  int fixedPower (const int number, const int pow) const
  {
    int result = 1;
    for (int i=0; i<pow; ++i)
      result *= number;
    return result;
  }

  // don't want = operator
  Acou& operator = (const Acou& old);

  // stores the degree of the element
  unsigned char degree_;

  // stores the polynomial type (tensor product or complete polynomial)
  bool completepol_;

  /// local contributions to the gradient in case of acoustical parameter optimization
  double grad_c_;
  double grad_rho_;

}; // class Acou


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief An element representing a boundary element of an acou element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate surface Neumann boundary conditions that might be
      adjacent to a parent acou element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

*/

class AcouBoundaryType : public DRT::ElementType
{
public:

  std::string Name() const { return "AcouBoundaryType"; }

  static AcouBoundaryType & Instance();

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static AcouBoundaryType instance_;
};

// class AcouBoundary

class AcouBoundary : public DRT::FaceElement
{
public:

  //! @name Constructors and destructors and related methods

  //! number of space dimensions
  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner: Processor owning this surface
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  \param parent: The parent acou element of this surface
  \param lsurface: the local surface number of this surface w.r.t. the parent element
  */
  AcouBoundary(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::Acou* parent, const int lsurface);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  AcouBoundary(const AcouBoundary& old);

  /*!
  \brief Deep copy this instance of an element and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const
  {
    return DRT::UTILS::getNumberOfElementLines(Shape());
  }

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const
  {
    return DRT::UTILS::getNumberOfElementSurfaces(Shape());
  }

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element

  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of Teuchos::RCPs to the surfaces of this element
  */

  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return AcouBoundaryType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const std::vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~AcouBoundary();

  //@}

  //! @name Acess methods


  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    return ParentElement()->NumDofPerNode(node);
  }

  virtual int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return AcouBoundaryType::Instance(); }

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate element

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \return 0 if successful, negative otherwise
  */
  int Evaluate(Teuchos::ParameterList&   params,
               DRT::Discretization&      discretization,
               std::vector<int>&         lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  //@}

  //! @name Evaluate methods

  /*!
  \brief Evaluate Neumann boundary condition

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): reference to the underlying discretization
  \param condition (in)     : condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  //@}

  /*!
  \brief Return the location vector of this element

  The method computes degrees of freedom this element adresses.
  Degree of freedom ordering is as follows:<br>
  First all degrees of freedom of adjacent nodes are numbered in
  local nodal order, then the element internal degrees of freedom are
  given if present.<br>
  If a derived element has to use a different ordering scheme,
  it is welcome to overload this method as the assembly routines actually
  don't care as long as matrices and vectors evaluated by the element
  match the ordering, which is implicitly assumed.<br>
  Length of the output vector matches number of degrees of freedom
  exactly.<br>
  This version is intended to fill the LocationArray with the dofs
  the element will assemble into. In the standard case these dofs are
  the dofs of the element itself. For some special conditions (e.g.
  the weak dirichlet boundary condtion) a surface element will assemble
  into the dofs of a volume element.<br>

  \note The degrees of freedom returned are not neccessarily only nodal dofs.
        Depending on the element implementation, output might also include
        element dofs.

  \param dis (in)      : the discretization this element belongs to
  \param la (out)      : location data for all dofsets of the discretization
  \param doDirichlet (in): whether to get the Dirichlet flags
  \param condstring (in): Name of condition to be evaluated
  \param condstring (in):  List of parameters for use at element level
  */
  virtual void LocationVector(const Discretization&   dis,
                              LocationArray&          la,
                              bool                    doDirichlet,
                              const std::string&      condstring,
                              Teuchos::ParameterList& params
                              ) const;

private:


  // don't want = operator
  AcouBoundary& operator = (const AcouBoundary& old);

}; // class AcouBoundary


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief An element representing an internal face element between two acou elements
*/
class AcouIntFaceType : public DRT::ElementType
{
public:

  std::string Name() const { return "AcouIntFaceType"; }

  static AcouIntFaceType & Instance();

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static AcouIntFaceType instance_;
};


// class AcouIntFace

class AcouIntFace : public DRT::FaceElement
{
public:

  //! @name Constructors and destructors and related methods

  //! number of space dimensions
  /*!
  \brief Standard Constructor

  \param id: A unique global id
  \param owner: Processor owning this surface
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  \param master_parent: The master parent acou element of this surface
  \param slave_parent: The slave parent acou element of this surface
  \param lsurface_master: the local surface number of this surface w.r.t. the master parent element
  \param lsurface_slave: the local surface number of this surface w.r.t. the slave parent element
  \param localtrafomap: transformation map between the local coordinate systems of the face w.r.t the master parent element's face's coordinate system and the slave element's face's coordinate system
  */
  AcouIntFace( int id,
                int owner,
                int nnode,
                const int* nodeids,
                DRT::Node** nodes,
                DRT::ELEMENTS::Acou* parent_master,
                DRT::ELEMENTS::Acou* parent_slave,
                const int lsurface_master,
                const int lsurface_slave,
                const std::vector<int> localtrafomap
                );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  AcouIntFace(const AcouIntFace& old);

  /*!
  \brief Deep copy this instance of an element and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;

  /*!
  \brief Return number of lines of this element
  */
  virtual int NumLine() const
  {
    return DRT::UTILS::getNumberOfElementLines(Shape());
  }

  /*!
  \brief Return number of surfaces of this element
  */
  virtual int NumSurface() const
  {
    return DRT::UTILS::getNumberOfElementSurfaces(Shape());
  }

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return AcouIntFaceType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element
  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element
  */
  virtual void Unpack(const std::vector<char>& data);

  /*!
  \brief Destructor
  */
  virtual ~AcouIntFace();

  //@}

  //! @name Acess methods

  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    return std::max(ParentMasterElement()->NumDofPerNode(node),
                    ParentSlaveElement()->NumDofPerNode(node));
  }

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  virtual int NumDofPerElement() const { return 0; }

  /*!
   \brief Returns the degree of the element
   */
  virtual int Degree() const
  {
    return degree_;
  }

  void SetDegree(int degree)
  {
    degree_ = degree;
    return;
  }

  /*!
  \brief create the location vector for patch of master and slave element

  \note All dofs shared by master and slave element are contained only once. Dofs from interface nodes are also included.
  */
  void PatchLocationVector(
      DRT::Discretization & discretization,       ///< discretization
      std::vector<int>&     nds_master,           ///< nodal dofset w.r.t master parent element
      std::vector<int>&     nds_slave,            ///< nodal dofset w.r.t slave parent element
      std::vector<int>&     patchlm,              ///< local map for gdof ids for patch of elements
      std::vector<int>&     master_lm,            ///< local map for gdof ids for master element
      std::vector<int>&     slave_lm,             ///< local map for gdof ids for slave element
      std::vector<int>&     face_lm,              ///< local map for gdof ids for face element
      std::vector<int>&     lm_masterToPatch,     ///< local map between lm_master and lm_patch
      std::vector<int>&     lm_slaveToPatch,      ///< local map between lm_slave and lm_patch
      std::vector<int>&     lm_faceToPatch,       ///< local map between lm_face and lm_patch
      std::vector<int>&     lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
      std::vector<int>&     lm_slaveNodeToPatch   ///< local map between slave nodes and nodes in patch
      );

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return AcouIntFaceType::Instance(); }

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate element

  \param params (in/out): ParameterList for communication between control routine
                          and elements
  \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                          the controling method does not epxect the element to fill
                          this matrix.
  \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                          the controlling method does not epxect the element
                          to fill this vector
  \return 0 if successful, negative otherwise
  */
  int Evaluate(Teuchos::ParameterList&   params,
               DRT::Discretization&      discretization,
               std::vector<int>&         lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  //@}

  //! @name Evaluate methods

  /*!
  \brief Evaluate Neumann boundary condition

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): reference to the underlying discretization
  \param condition (in)     : condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief return the master parent acou element
  */
  DRT::ELEMENTS::Acou* ParentMasterElement() const
  {
    DRT::Element* parent = this->DRT::FaceElement::ParentMasterElement();
    // make sure the static cast below is really valid
    dsassert(dynamic_cast<DRT::ELEMENTS::Acou*>(parent) != NULL, "Master element is no acou element");
    return static_cast<DRT::ELEMENTS::Acou*>(parent);
  }

  /*!
  \brief return the slave parent acou element
  */
  DRT::ELEMENTS::Acou* ParentSlaveElement() const
  {
    DRT::Element* parent = this->DRT::FaceElement::ParentSlaveElement();
    // make sure the static cast below is really valid
    dsassert(dynamic_cast<DRT::ELEMENTS::Acou*>(parent) != NULL, "Slave element is no acou element");
    return static_cast<DRT::ELEMENTS::Acou*>(parent);
  }

  //@}

private:

  // don't want = operator
  AcouIntFace& operator = (const AcouIntFace& old);

  // degree of this face element
  int degree_;

}; // class AcouIntFace



} // namespace ELEMENTS
} // namespace DRT


#endif /* ACOU_ELE_H */
