/*!----------------------------------------------------------------------
\file so_surface.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#ifndef SO_SURFACE_H
#define SO_SURFACE_H

#include "Teuchos_RefCountPtr.hpp"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_fem_general/drt_utils_integration.H"

#if defined(PRESTRESS) || defined(POSTSTRESS)
// define whether orthopressure shall respect prestressing or not
#define PRESTRESS_ORTHOPRESSURE 0  // 0: off  1: on
#endif

namespace DRT
{
namespace ELEMENTS
{

/*!
\brief A general structural surface element (tri3/tri6/quad4/quad8/quad9)

This element implements all common shapes. It also implements all
common types of loads currently demanded in baci

\author gee (gee@lnm.mw.tum.de)
*/
class StructuralSurface : public DRT::Element
{
public:

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
  explicit StructuralSurface(int id, int owner, int nnode, const int* nodeids,
                             DRT::Node** nodes, DRT::Element* parent, const int lsurface);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  explicit StructuralSurface(const StructuralSurface& old);

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
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of the drt_parobject.H file.
  */
  virtual inline int UniqueParObjectId() const
  { return ParObject_StructuralSurface;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~StructuralSurface() {}

  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  inline virtual int NumDofPerNode(const DRT::Node& node) const
  {return 3;}

  /*!
  \brief Get number of degrees of freedom per element
         (implements pure virtual DRT::Element)

  The element decides how many element degrees of freedom it has.
  It can redecide along the way of a simulation.

  \note Element degrees of freedom mentioned here are dofs that are visible
        at the level of the total system of equations. Purely internal
        element dofs that are condensed internally should NOT be considered.
  */
  inline virtual int NumDofPerElement() const { return 0; }

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  /*!
  \brief Allocate and return a ElementRegister class

  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  */
  virtual RCP<DRT::ElementRegister> ElementRegister() const
  { return rcp(new DRT::ElementRegister(Type())); }


  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a surface Neumann condition on the fluid3 element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(ParameterList&            params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);

  //! Evaluate method for StructuralSurface-Elements
  virtual int Evaluate(ParameterList& params,
                       DRT::Discretization&      discretization,
                       vector<int>&              lm,
                       Epetra_SerialDenseMatrix& elematrix1,
                       Epetra_SerialDenseMatrix& elematrix2,
                       Epetra_SerialDenseVector& elevector1,
                       Epetra_SerialDenseVector& elevector2,
                       Epetra_SerialDenseVector& elevector3);


protected:

  //! action for surface evaluation
  enum ActionType
  {
    none,
    calc_struct_constrvol,
    calc_struct_volconstrstiff,
    calc_struct_constrarea,
    calc_init_vol,
    calc_surfstress_stiff,
    calc_potential_stiff,
    calc_struct_areaconstrstiff,
    calc_struct_monitarea,
    prestress_update,
    prestress_writerestart,
    prestress_readrestart
  };

  // don't want = operator
  StructuralSurface& operator = (const StructuralSurface& old);

  // variables

  //! The parent element of this surface
  DRT::Element*               parent_;
  //! The local surface number of this surface w.r.t to the parent_ element
  int                         lsurface_;
  //! gaussian integration to be used
  DRT::UTILS::GaussRule2D     gaussrule_;

#if defined(PRESTRESS) || defined(POSTSTRESS)
#if PRESTRESS_ORTHOPRESSURE
  LINALG::SerialDenseMatrix   xrefehist_;
#endif
#endif

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates normal and detA at gaussian point

  \param detA   (out) : are at gaussian point
  \param normal (out) : The normal at gaussian point, length is detA!
  \param x      (in)  : nodal coords in either material or spatial frame
  \param deriv  (in)  : derivatives of shape functions
  */
  void SurfaceIntegration(double& detA,
                          vector<double>& normal,
                          const Epetra_SerialDenseMatrix& x,
                          const Epetra_SerialDenseMatrix& deriv);

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates the normal at gaussian point

  \param normal (out) : The normal at gaussian point, length is detA!
  \param x      (in)  : nodal coords in either material or spatial frame
  \param deriv  (in)  : derivatives of shape functions
  */
  void SurfaceIntegration(vector<double>& normal,
                          const Epetra_SerialDenseMatrix& x,
                          const Epetra_SerialDenseMatrix& deriv);

  /*!
  \brief Create matrix with material configuration

  \param x  (out)  : nodal coords in material frame
  */
  inline void MaterialConfiguration(LINALG::SerialDenseMatrix& x) const
  {
    const int numnode = NumNode();
    for (int i=0; i<numnode; ++i)
    {
      x(i,0) = Nodes()[i]->X()[0];
      x(i,1) = Nodes()[i]->X()[1];
      x(i,2) = Nodes()[i]->X()[2];
    }
    return;
  }

  /*!
  \brief Create matrix with spatial configuration

  \param x     (out)  : nodal coords in spatial frame
  \param disp  (int)  : displacements
  */
  inline void SpatialConfiguration(LINALG::SerialDenseMatrix& x,
                                   const vector<double> disp) const
  {
    const int numnode = NumNode();
    for (int i=0; i<numnode; ++i)
    {
      x(i,0) = Nodes()[i]->X()[0] + disp[i*3+0];
      x(i,1) = Nodes()[i]->X()[1] + disp[i*3+1];
      x(i,2) = Nodes()[i]->X()[2] + disp[i*3+2];
    }
    return;
  }

  /*!
  \brief Create matrix with spatial configuration

  \param x     (out)  : nodal coords in spatial frame
  \param xrefe (out)  : nodal coords in material frame
  \param disp  (int)  : displacements
  */
  inline void SpatialConfiguration(LINALG::SerialDenseMatrix& x,
                                   LINALG::SerialDenseMatrix& xrefe,
                                   const vector<double> disp) const
  {
    const int numnode = NumNode();
    for (int i=0; i<numnode; ++i)
    {
      x(i,0) = xrefe(i,0) + disp[i*3+0];
      x(i,1) = xrefe(i,1) + disp[i*3+1];
      x(i,2) = xrefe(i,2) + disp[i*3+2];
    }
    return;
  }

  //! Submethod to compute the enclosed volume for vlume constraint boundary condition
  double ComputeConstrVols(const LINALG::SerialDenseMatrix& xscurr);

  //! Submethod to compute necessary change to stiffness matrix due to the constraints
  void ComputeVolConstrStiff(const LINALG::SerialDenseMatrix& xscurr,
                             Epetra_SerialDenseMatrix& elematrix);

  //! Submethod to compute first derivatives of constraint volumes w.r.t. the displacements
  void ComputeVolConstrDeriv(const LINALG::SerialDenseMatrix& xscurr,
                             Epetra_SerialDenseVector& elevector1);

  //! Submethod to compute interfacial area and its first and second
  //derivatives w.r.t. the displacements (this is needed for surface energy problems)
  void ComputeAreaDeriv(const LINALG::SerialDenseMatrix& x,
                        const int numnode,
                        const int ndof,
                        double& A,
                        RCP<Epetra_SerialDenseVector> Adiff,
                        RCP<Epetra_SerialDenseMatrix> Adiff2);

}; // class StructuralSurface : public DRT::Element


} // namespace ELEMENTS
} // namespace DRT


#endif // #ifndef SO_SURFACE_H
#endif // #ifdef CCADISCRET
#endif // #ifdef D_SOLID3
