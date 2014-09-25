/*!---------------------------------------------------------------------------

\file meshfree_fluid_cell.H

\brief A C++ wrapper for the fluid element

This file contains the cell-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

In addition to that, it contains the interface between element call
and Gauss point loop (depending on the fluid implementation)
as well as some additional service routines (for the evaluation
of errors, turbulence statistics etc.).

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*---------------------------------------------------------------------------*/

#ifndef MESHFREE_FLUID_CELL_H
#define MESHFREE_FLUID_CELL_H

#include "drt_meshfree_cell.H"
#include "../drt_lib/drt_elementtype.H"

namespace DRT
{

class Discretization;

namespace ELEMENTS
{

/*==========================================================================*/
// forward declarations
/*==========================================================================*/
class FluidEleParameter;
class FluidEleParameterStd;
class FluidEleParameterTimInt;

/*--------------------------------------------------------------------------*/
/*!
 * \brief The element type class of meshfree fluid cells, creating the same
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date Januar, 2013
 */
/*--------------------------------------------------------------------------*/
class MeshfreeFluidType : public DRT::ElementType
{
public:

  //!< name of specific element type
  std::string Name() const
  {
    return "MeshfreeFluidType";
  }

  //!< returning instance of specific element type
  static MeshfreeFluidType & Instance()
  {
    return instance_;
  }

  //!< create parallel object of this element type
  virtual DRT::ParObject* Create( const std::vector<char> & data );

  //!< create element of this element type
  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  //!< create element of this element type
  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

  /// pre-evaluation ???
  virtual void PreEvaluate(DRT::Discretization&                dis,
                          Teuchos::ParameterList&              p,
                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                          Teuchos::RCP<Epetra_Vector>          systemvector1,
                          Teuchos::RCP<Epetra_Vector>          systemvector2,
                          Teuchos::RCP<Epetra_Vector>          systemvector3);

private:

  //!< static instance of this element type (for self-instantiation)
  static MeshfreeFluidType instance_;
};

/*--------------------------------------------------------------------------*/
/*!
 * \brief A C++ wrapper for the MeshfreeFluid cell
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date January, 2013
 */
/*--------------------------------------------------------------------------*/
class MeshfreeFluid : public DRT::MESHFREE::Cell
{
public:

  //! @name friends
  friend class MeshfreeFluidBoundary;
//  friend class FluidSystemEvaluator; // ???

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  MeshfreeFluid(
    int id,    //!< A unique global id
    int owner  //!< processor id who owns a certain instance of this class
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Copy Constructor
   *
   * Makes a deep copy of a DRT::Element
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  MeshfreeFluid(const MeshfreeFluid& old);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Deep copy this instance of FluidTransport and return pointer to
   *        the copy
   *
   *   The Clone() method is used from the virtual base class Element in cases
   *   where the type of the derived class is unknown and a copy-ctor is
   *   needed
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual DRT::Element* Clone() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Destructor
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual ~MeshfreeFluid();

  /*========================================================================*/
  //! @name access methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get shape type of cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual DiscretizationType Shape() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief set discretization type of cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual void SetDisType(DRT::Element::DiscretizationType shape)
  {
    distype_ = shape;

    if (distype_ == DRT::Element::tet10)
    {
      std::cout << std::endl;
      std::cout << "Warning:The integration rule of a tet10 element was changed from tet-5-point to tet-11-point integration rule, "
        " since the deg. of precision 3 is to small as a standard"<< std::endl << std::endl;
    }
    return;
  };

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of lines of this cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int NumLine() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of surfaces of this cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int NumSurface() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of volumes of this cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int NumVolume() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the lines of this cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the surfaces of this cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the volume of this cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Volumes();


  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return unique ParObject id
   *
   * Every class implementing ParObject needs a unique id defined at the top
   * of this file.
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int UniqueParObjectId() const
  {
    return MeshfreeFluidType::Instance().UniqueParObjectId();
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Pack this class so it can be communicated
   *
   * \ref Pack and \ref Unpack are used to communicate this cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual void Pack(DRT::PackBuffer& data) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Unpack data from a char vector into this class
   *
   * \ref Pack and \ref Unpack are used to communicate this cell
   *///                                                   (public) nis Jan13
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
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int NumDofPerNode(const DRT::Node& node) const;

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
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int NumDofPerElement() const
  {
    return 0;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this element
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  void Print(std::ostream& os) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief TODO
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual DRT::ElementType & ElementType() const
  {
    return MeshfreeFluidType::Instance();
  }

  /*========================================================================*/
  //! @name Input and Creation
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Read input for this element
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  bool ReadElement(const std::string&          eletype,
                   const std::string&          distype,
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
   *///                                                   (public) nis Jan13
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
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  //! @name Other

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return flag for ALE form of equations
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  bool IsAle() const
  {
    return is_ale_;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Set flag for ALE form of equations
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual void SetIsAle(bool is_ale)
  {
    is_ale_ = is_ale;

    return;
  };


private:

  //! don't want = operator
  MeshfreeFluid& operator = (const MeshfreeFluid& old);

  //! discretization type
  DRT::Element::DiscretizationType distype_;

  //! flag for euler/ale net algorithm
  bool is_ale_;

}; // class MeshfreeFluid


/*--------------------------------------------------------------------------*/
/*!
 * \brief The element type class of meshfree transport boundary cells
 *
 * TODO: Are boundary type classes really necessary? Create doesn't create
 *       anything!
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date January, 2013
 */
/*--------------------------------------------------------------------------*/
class MeshfreeFluidBoundaryType : public DRT::ElementType
{
public:

  //!< name of specific parallel object type
  std::string Name() const { return "MeshfreeFluidBoundaryType"; }

  //!< returning instance of specific parallel object type
  static MeshfreeFluidBoundaryType & Instance() { return instance_; };

  //!< create parallel object of this parallel object type
  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  //!< static instance of this parallel object type (for self-instantiation)
  static MeshfreeFluidBoundaryType instance_;
};


/*--------------------------------------------------------------------------*/
/*!
 * \brief An element representing a boundary cell of a MeshfreeTransport cell
 *
 * \note This is a pure Neumann boundary condition element. It's only
 *       purpose is to evaluate surface Neumann boundary conditions that might be
 *       adjacent to a parent meshfree fluid element. It therefore does not implement
 *       the DRT::Element::Evaluate method and does not have its own ElementRegister class.
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date January, 2013
 */
/*--------------------------------------------------------------------------*/
class MeshfreeFluidBoundary : public DRT::MESHFREE::Cell
{
public:

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Standard Constructor
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  MeshfreeFluidBoundary(
    int id,                              //!< A unique global id
    int owner,                           //!< Processor owning this surface
    int npoint,                          //!< Number of nodes attached to this element
    const int* pointids,                 //!< global ids of nodes attached to this element
    DRT::Node** points,                  //!< the discretizations map of points to build ptrs to points
    MeshfreeFluid* parent,               //!< The parent fluid element of this surface
    const int lsurface                   //!< the local surface number of this surface w.r.t. the parent element
  );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Copy Constructor
   *
   * Makes a deep copy of a MeshfreeTransportBoundary
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  MeshfreeFluidBoundary(const MeshfreeFluidBoundary& old);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Deep copy the derived class and return pointer to it
   *
   * This method is sort of a copy constructor for a class derived from
   * DRT::Element. It allows to copy construct the derived class without
   * knowing what it actually is using the base class DRT::Element.
   *
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  DRT::Element* Clone() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Destructor
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual ~MeshfreeFluidBoundary();

  /*========================================================================*/
  //! @name Query methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get shape type of element
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual DiscretizationType Shape() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of lines of boundary cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int NumLine() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return number of surfaces of boundary cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int NumSurface() const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the lines of this cell
   *        Dummy - not provided.
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get vector of Teuchos::RCPs to the surfaces of this cell
   *        Dummy - not provided.
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Pack dummy class for - not supported
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual void Pack(DRT::PackBuffer& data) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Unpack dummy class - not supported
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual void Unpack(const std::vector<char>& data);

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Get number of degrees of freedom of a certain node in case of
   *        multiple dofsets, i.e. monolithic solution approaches!e
   *        (implements pure virtual DRT::Element)
   *
   * The element decides how many degrees of freedom its nodes must have. As
   * this may vary along a simulation, the element can redecide the number of
   * degrees of freedom per node along the way for each of it's nodes
   * separately.
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    return ParentElement()->NumDofPerNode(node);
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
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int NumDofPerElement() const
  {
    return 0;
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return pointer to the parent element
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual DRT::ELEMENTS::MeshfreeFluid* ParentElement() const
  {
    DRT::Element* parent = this->DRT::Element::ParentElement();
    // make sure the static cast below is really valid
    dsassert(dynamic_cast<DRT::ELEMENTS::MeshfreeFluid*>(parent) != NULL, "Parent element is no fluid element");
    return dynamic_cast<DRT::ELEMENTS::MeshfreeFluid*>(parent);
  }


  /*------------------------------------------------------------------------*/
  /*!
   * \brief Return unique ParObject id
   *
   * every class implementing ParObject needs a unique id defined at the
   * top of the drt_parobject.H file.
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int UniqueParObjectId() const
  {
    return MeshfreeFluidBoundaryType::Instance().UniqueParObjectId();
  }

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Print this cell
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  void Print(std::ostream& os) const;

  /*------------------------------------------------------------------------*/
  /*!
   * \brief TODO
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual DRT::ElementType & ElementType() const
  {
    return MeshfreeFluidBoundaryType::Instance();
  }

  /*========================================================================*/
  //! @name Evaluate methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate an element
   *
   * Evaluate MeshfreeFluid element stiffness, mass, internal forces etc
   *
   * \param params (in/out): ParameterList for communication between control
   *                         routine and elements
   * \param discretization (in): A reference to the underlying discretization
   * \param lm (in)        : location vector of this element
   * \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
   *                         the controling method does not epxect the element
   *                         to fill this matrix.
   * \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
   *                         the controling method does not epxect the element
   *                         to fill this matrix.
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
   *///                                                   (public) nis Jan13
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
   * this method evaluates a surface Neumann condition on the MeshfreeFluid cell
   *
   * \param params (in/out)    : ParameterList for communication between control routine
   *                             and elements
   * \param discretization (in): A reference to the underlying discretization
   * \param condition (in)     : The condition to be evaluated
   * \param lm (in)            : location vector of this element
   * \param elevec1 (out)      : vector to be filled by element. If NULL on input,
   *
   * \return 0 if successful, negative otherwise
   *///                                                   (public) nis Jan13
  /*------------------------------------------------------------------------*/
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

private:


  // don't want = operator
  MeshfreeFluidBoundary& operator = (const MeshfreeFluidBoundary& old);

}; // class MeshfreeFluidBoundary

} // namespace ELEMENTS
} // namespace DRT

#endif  // #ifndef MESHFREE_FLUID_CELL_H
