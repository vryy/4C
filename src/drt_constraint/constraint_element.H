/*!----------------------------------------------------------------------
\file constraint_element.H

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONSTRAINT_ELEMENT_H
#define CONSTRAINT_ELEMENT_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_serialdensematrix.H"



namespace DRT
{
// forward declarations
class Discretization;
class ConstraintElementRegister;

namespace ELEMENTS
{
// forward declarations
//class ConstraintElementLine;

/*!
\brief A 3D constraint element with no physics attached

\author tk (kloeppel@lnm.mw.tum.de)
*/
class ConstraintElement : public DRT::Element
{
public:

  //! @name Friends
  friend class ConstraintElementRegister;

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  
  */
  ConstraintElement(
    int id,  ///< A unique global id 
    int owner ///< element owner
  );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  ConstraintElement(const ConstraintElement& old);

  /*!
  \brief Deep copy this instance of ConstraintElement and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const
  {
      dserror("ConstraintElement has no shape!");
      return DRT::Element::dis_none;
  };

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_ConstraintElement;}

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
  virtual ~ConstraintElement();

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
  virtual int NumDofPerNode(const DRT::Node& node) const {return numdof_;}

  /*!
  \brief Set number of degrees of freedom of all connected nodes

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  void SetNumDofPerNode(const int numdof) {numdof_ = numdof;}
  
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

  /*!
  \brief Allocate and return a DesignElementRegister class

  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const;

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  bool ReadElement()
  {
      dserror("No reading for this element! Will be created on the fly, not from a .dat file.");
      return false;
  }

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate fluid2 element stiffness, mass, internal forces etc

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
  int Evaluate(ParameterList& params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);


  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a line Neumann condition on the fluid2 element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1);


private:

  //! action parameters recognized by bele3
  enum ActionType
  {
    none,
    calc_MPC3D_stiff,
    calc_MPC3D_state,
    calc_MPC2D_dist_stiff,
    calc_MPC2D_angle_stiff
  };

  //! container containing all kind of fluid2 thisandthat
  DRT::Container                                  data_;

  //! vector of surfaces of this element (length 1)
  vector<DRT::Element*>                           surface_;

  //! Number of dof per node
  int                                             numdof_;
  
  // don't want = operator
  ConstraintElement& operator = (const ConstraintElement& old);
  
  /*!
  \brief Create matrix with material configuration

  \param x  (out)  : nodal coords in material frame
  */
  inline void MaterialConfiguration(LINALG::SerialDenseMatrix& x, const int& numdim) const
  {
    const int numnode = NumNode();
    for (int i=0; i<numnode; ++i)
    {
      for (int j = 0; j < numdim; ++j) 
      {
        x(i,j) = Nodes()[i]->X()[j];
      }
    }
    return;
  }
  
  /*!
  \brief Create matrix with spatial configuration

  \param x     (out)  : nodal coords in spatial frame
  \param disp  (int)  : displacements
  */
  inline void SpatialConfiguration(LINALG::SerialDenseMatrix& x, 
                                   const vector<double> disp,
                                   const int& numdim) const
  {
    const int numnode = NumNode();
    for (int i=0; i<numnode; ++i)
    {
      for (int j = 0; j < numdim; ++j) 
      {
        x(i,j) = Nodes()[i]->X()[j] + disp[i*numdim+j];
      }
    }
    return;
  }
  
  void ComputeElementNormal3D(const LINALG::SerialDenseMatrix& xc,
      LINALG::SerialDenseVector& elenormal);
  
  double ComputeNormalDist3D(const LINALG::SerialDenseMatrix& xc,
      const LINALG::SerialDenseVector& elenormal);
  
  void ComputeFirstDeriv3D(const LINALG::SerialDenseMatrix& xc,
      Epetra_SerialDenseVector& elevector,
      const LINALG::SerialDenseVector& elenormal);
  
  void ComputeSecondDeriv3D(const LINALG::SerialDenseMatrix&,
      Epetra_SerialDenseMatrix& elematrix,
      const LINALG::SerialDenseVector& elenormal);

  void ComputeElementNormal2D(const LINALG::SerialDenseMatrix& xc,
      LINALG::SerialDenseVector& elenormal);
  
  double ComputeNormalDist2D(const LINALG::SerialDenseMatrix& xc,
      const LINALG::SerialDenseVector& elenormal);
  
  void ComputeFirstDerivDist2D(const LINALG::SerialDenseMatrix& xc,
      Epetra_SerialDenseVector& elevector,
      const LINALG::SerialDenseVector& elenormal);
  
  void ComputeSecondDerivDist2D(const LINALG::SerialDenseMatrix&,
      Epetra_SerialDenseMatrix& elematrix,
      const LINALG::SerialDenseVector& elenormal);

  double ComputeElementAngle2D(const LINALG::SerialDenseMatrix& xc);
  
  void ComputeFirstDerivAngle2D(const LINALG::SerialDenseMatrix& xc,
      Epetra_SerialDenseVector& elevector);
  
  void ComputeSecondDerivAngle2D(const LINALG::SerialDenseMatrix& xs,
      Epetra_SerialDenseMatrix& elematrix);


}; // class ConstraintElement


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*!
\brief A register for 3D constraint element

\author tk (kloeppel@lnm.mw.tum.de)
*/
class ConstraintElementRegister : public DRT::ElementRegister
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  ConstraintElementRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an ConstraintElementRegister

  */
  ConstraintElementRegister(const DRT::ELEMENTS::ConstraintElementRegister& old);

  /*!
  \brief Deep copy this instance of ConstraintElementRegister and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::ConstraintElementRegister* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class imploementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_ConstraintElementRegister;}

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
  virtual ~ConstraintElementRegister();

  //@}

  //! @name Acess methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const;

  //@}

  //! @name Construction

  /*!
  \brief Call Element init routines

  The method performs all steps needed to init the fluid2 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:

}; // class ConstraintElementRegister

} // namespace ELEMENTS
} // namespace DRT


#endif
#endif
