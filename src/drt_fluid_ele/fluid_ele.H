/*!
\file fluid_ele.H

\brief A C++ wrapper for the fluid element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

In addition to that, it contains the interface between element call
and Gauss point loop (depending on the fluid implementation)
as well as some additional service routines (for the evaluation
of errors, turbulence statistics etc.).

<pre>
Maintainer: Volker Gravemeier & Andreas Ehrl
            {vgravem,ehrl}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15245/15252
</pre>
*/
#ifndef DRT_FLUID_ELE_H
#define DRT_FLUID_ELE_H

#include "Teuchos_RefCountPtr.hpp"


#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

//=================================================================
//
// This file still contains some deactivated time-dependent subgrid-
// scales stuff. As soon as the time-dependent subgrid-scales are
// fixed, these members and methods should be removed from the element!
// Ursula
//
//=================================================================


// forward declarations
namespace DRT
{

class Discretization;

namespace ELEMENTS
{

class FluidEleParameter;
template<DRT::Element::DiscretizationType distype> class FluidEleCalc;
template<DRT::Element::DiscretizationType distype> class FluidGenalphaResVMM;
class FluidBoundary;
class FluidIntFace;

class FluidType : public DRT::ElementType
{
public:

  std::string Name() const { return "FluidType"; }

  static FluidType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const string eletype,
                                             const string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np );

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

  static FluidType instance_;
};

/*!
\brief A C++ wrapper for the fluid element
*/
class Fluid : public DRT::Element
{
public:

  //! @name friends
  friend class FluidBoundary;
  friend class FluidIntFace;

  friend class FluidSystemEvaluator;

  //@}
  //! @name constructors and destructors and related methods

  /*!
  \brief standard constructor
  */
  Fluid(
    int id,    ///< A unique global id
    int owner  ///< ???
    );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Fluid(const Fluid& old);

  /*!
  \brief Deep copy this instance of fluid and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

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

    if (distype_ == DRT::Element::tet10)
    {
      cout << endl;
      cout << "Warning:The integration rule of a tet10 element was changed from tet-5-point to tet-11-point integration rule, "
      " since the deg. of precision 3 is to small as a standard"<< endl << endl;
    }
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
  \brief Get vector of RefCountPtrs to the lines of this element
  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of RefCountPtrs to the surfaces of this element
  */
  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Get vector of RefCountPtrs to the volume of this element
  */
  virtual vector<RCP<DRT::Element> > Volumes();

  /*!
  \brief Get RefCountPtr to the internal face adjacent to this element as master element and the parent_slave element
  */
  virtual RCP<DRT::Element> CreateInternalFaces(DRT::Element* parent_slave,   //!< parent slave fluid3 element
                                                int nnode,                    //!< number of surface nodes
                                                const int* nodeids,           //!< node ids of surface element
                                                DRT::Node** nodes,            //!< nodes of surface element
                                                const int lsurface_master     //!< local surface number w.r.t master parent element
      );

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return FluidType::Instance().UniqueParObjectId(); }

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
  virtual ~Fluid();

  //@}

  //! @name Geometry related methods

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
    // number of Dof's is fluid-specific.
    // Therefore, it is not pushed into drt_utilis_local_connectivity
    switch(distype_)
    {
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    case DRT::Element::tet4:
    case DRT::Element::tet10:
    case DRT::Element::wedge6:
    case DRT::Element::wedge15:
    case DRT::Element::pyramid5:
    case DRT::Element::nurbs8:
    case DRT::Element::nurbs27:
      return 4;
      break;
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::nurbs4:
    case DRT::Element::nurbs9:
    case DRT::Element::tri3:
    case DRT::Element::tri6:
      return 3;
      break;
    case DRT::Element::line2:
    case DRT::Element::line3:
      dserror("1D Fluid elements are not supported");
    default:
      dserror("discretization type %s not yet implemented", (DRT::DistypeToString(distype_)).c_str());
    }
    return 0;
  }

  /*!
  \brief Get number of degrees of freedom of a certain node in case of multiple
         dofsets, i.e. monolithic solution approaches! (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.

  */
  virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node) const;

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
  \brief Print this element
  */
  void Print(ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return FluidType::Instance(); }

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
  \brief Evaluate an element, that is, call the element routines to evaluate fluid
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
  int Evaluate(ParameterList&            params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
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
  virtual int EvaluateNeumann(ParameterList&            params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

#if 0
  //@}

  //! @name Time-dependent subgrid scales
  /*!
  \brief Memory allocation for subgrid-scale arrays
  */
  void ActivateTDS(int      nquad,
                   int      nsd,
                   double** saccn=NULL,
                   double** sveln=NULL,
                   double** svelnp=NULL);

  /*!
  \brief Nonlinear update for current subgrid-scale velocities according to the current
         residual (reduced version for afgenalpha and one-step-theta)
  */
  void UpdateSvelnpInOneDirection(const double  fac1,
                                  const double  fac2,
                                  const double  fac3,
                                  const double  resM,
                                  const double  alphaF,
                                  const int     dim,
                                  const int     iquad,
                                  double&       svelaf);

  /*!
  \brief Nonlinear update for current subgrid-scale velocities according to the current
         residual (svelnp as additional return value)
  */
  void UpdateSvelnpInOneDirection(const double  fac1,
                                  const double  fac2,
                                  const double  fac3,
                                  const double  resM,
                                  const double  alphaF,
                                  const int     dim,
                                  const int     iquad,
                                  double&       svelnp,
                                  double&       svelaf);

  //@}
#endif


  //! @name Other

  /*!
  \brief Flag for ALE form of equations
  */
  bool IsAle() const { return is_ale_; }

  virtual void SetIsAle(bool is_ale)
  {
	  is_ale_ = is_ale;

    return;
  };


#if 0
  /*!
  \brief Returns the subgrid velocity at time n (sveln_)
  */
  Epetra_SerialDenseMatrix Sveln() const {return sveln_; }

  /*!
  \brief Returns the subgrid velocity at time n+1 (svelnp_)
  */
  Epetra_SerialDenseMatrix Svelnp() const {return svelnp_; }
#endif

private:

  //! discretization type
  DRT::Element::DiscretizationType distype_;

  //! flag for euler/ale net algorithm
  bool is_ale_;

#if 0
  //! matrices of subgrid-scale acceleration values at integration points of this element
  Epetra_SerialDenseMatrix saccn_;
  //! matrices of subgrid-scale velocity values, current iteration value, at integration points of this element
  Epetra_SerialDenseMatrix svelnp_;
  //! matrices of subgrid-scale velocity values, last timestep, at integration points of this element
  Epetra_SerialDenseMatrix sveln_;
#endif

  // internal calculation methods

  // don't want = operator
  Fluid& operator = (const Fluid& old);

}; // class Fluid


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief An element representing a boundary element of a fluid element

\note This is a pure Neumann boundary condition element. It's only
      purpose is to evaluate surface Neumann boundary conditions that might be
      adjacent to a parent fluid element. It therefore does not implement
      the DRT::Element::Evaluate method and does not have its own ElementRegister class.

*/

class FluidBoundaryType : public DRT::ElementType
{
public:

  std::string Name() const { return "FluidBoundaryType"; }

  static FluidBoundaryType & Instance() { return instance_; };

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static FluidBoundaryType instance_;
};


// class FluidBoundary

class FluidBoundary : public DRT::Element
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
  \param parent: The parent fluid element of this surface
  \param lsurface: the local surface number of this surface w.r.t. the parent element
  */
  FluidBoundary(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::Fluid* parent, const int lsurface);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  FluidBoundary(const FluidBoundary& old);

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
  \brief Get vector of RefCountPtrs to the lines of this element

  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */

  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return FluidBoundaryType::Instance().UniqueParObjectId(); }

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
  virtual ~FluidBoundary();

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
    return parent_->NumDofPerNode(node);
  }

  /*!
  \brief Get number of degrees of freedom of a certain node in case of multiple
         dofsets, i.e. monolithic solution approaches! (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.

  */
  virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node) const
  {
    return parent_->NumDofPerNode(nds,node);
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
  \brief Print this element
  */
  void Print(ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return FluidBoundaryType::Instance(); }

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
  int Evaluate(ParameterList&            params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
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
  virtual int EvaluateNeumann(ParameterList&            params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  virtual DRT::ELEMENTS::Fluid* ParentElement()
  {
    return parent_;
  }

  virtual int SurfaceNumber()
  {
    return lsurface_;
  }

  //@}

private:


  // don't want = operator
  FluidBoundary& operator = (const FluidBoundary& old);

  //! The parent element of this surface
  DRT::ELEMENTS::Fluid* parent_;

  //! The local surface number of this surface w.r.t to the parent_ element
  int lsurface_;

}; // class FluidBoundary


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief An element representing an internal face element between two fluid elements

\note It's only purpose is to evaluate edge based stabilizations for XFEM.
*/
class FluidIntFaceType : public DRT::ElementType
{
public:

  std::string Name() const { return "FluidIntFaceType"; }

  static FluidIntFaceType & Instance() { return instance_; };

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static FluidIntFaceType instance_;
};


// class FluidIntFace

class FluidIntFace : public DRT::Element
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
  \param master_parent: The master parent fluid element of this surface
  \param slave_parent: The slave parent fluid element of this surface
  \param lsurface_master: the local surface number of this surface w.r.t. the master parent element
  */
  FluidIntFace(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::Fluid* parent_master, DRT::ELEMENTS::Fluid* parent_slave, const int lsurface_master);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  FluidIntFace(const FluidIntFace& old);

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
  \brief Get vector of RefCountPtrs to the lines of this element
  */
  virtual vector<RCP<DRT::Element> > Lines();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual vector<RCP<DRT::Element> > Surfaces();

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual int UniqueParObjectId() const { return FluidIntFaceType::Instance().UniqueParObjectId(); }

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
  virtual ~FluidIntFace();

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
    return std::max(parent_master_->NumDofPerNode(node),parent_slave_->NumDofPerNode(node));
  }

  /*!
  \brief Get number of degrees of freedom of a certain node in case of multiple
         dofsets, i.e. monolithic solution approaches! (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.

  */
  virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node) const
  {
    return std::max(parent_master_->NumDofPerNode(nds,node),parent_slave_->NumDofPerNode(nds,node));
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
  void Print(ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return FluidIntFaceType::Instance(); }

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
  int Evaluate(ParameterList&            params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
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
  virtual int EvaluateNeumann(ParameterList&            params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief return the master parent fluid element
  */
  virtual DRT::ELEMENTS::Fluid* ParentMasterElement()
  {
    return parent_master_;
  }

  /*!
  \brief return the slave parent fluid element
  */
  virtual DRT::ELEMENTS::Fluid* ParentSlaveElement()
  {
    return parent_slave_;
  }

  /*!
  \brief return the local surface number w.r.t. the master parent fluid element
  */
  virtual int SurfaceMasterNumber()
  {
    return lsurface_master_;
  }

  //@}

private:

  // don't want = operator
  FluidIntFace& operator = (const FluidIntFace& old);

  //! The master parent element of this surface
  DRT::ELEMENTS::Fluid* parent_master_;

  //! The parent slave element of this surface
  DRT::ELEMENTS::Fluid* parent_slave_;

  //! The local surface number of this surface w.r.t to the parent_master_ element
  int lsurface_master_;


}; // class FluidIntFace



} // namespace ELEMENTS
} // namespace DRT




#endif  // #ifndef FLUID3_H
