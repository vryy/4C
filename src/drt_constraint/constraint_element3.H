/*!----------------------------------------------------------------------
\file constraint_element3.H
\brief A 3D constraint element with no physics attached
\level 2
<pre>

\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#ifndef CONSTRAINT_ELEMENT3_H
#define CONSTRAINT_ELEMENT3_H


#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_node.H"
#include "../linalg/linalg_fixedsizematrix.H"


namespace DRT
{
class Discretization;

namespace ELEMENTS
{
// forward declarations

class ConstraintElement3Type : public DRT::ElementType
{
public:

  std::string Name() const { return "ConstraintElement3Type"; }

  static ConstraintElement3Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

private:

  static ConstraintElement3Type instance_;
};

/*!
*/
class ConstraintElement3 : public DRT::Element
{
public:

  //! @name Friends

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  ConstraintElement3(
    int id,  ///< A unique global id
    int owner ///< element owner
  );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  ConstraintElement3(const ConstraintElement3& old);

  /*!
  \brief Deep copy this instance of ConstraintElement3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const
  {
      dserror("ConstraintElement3 has no shape!");
      return DRT::Element::dis_none;
  };

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ConstraintElement3Type::Instance().UniqueParObjectId(); }

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
  virtual ~ConstraintElement3();

  //@}

  //! @name Access methods


  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const DRT::Node& node) const {return 3;}

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
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return ConstraintElement3Type::Instance(); }

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  \return 0 if successful, negative otherwise
  */
  int Evaluate
  (
    Teuchos::ParameterList& params, ///< ParameterList for communication
    DRT::Discretization& discretization,  ///< discretization
    std::vector<int>& lm,  ///< location vector
    Epetra_SerialDenseMatrix& elemat1, ///< first matrix to be filled by element
    Epetra_SerialDenseMatrix& elemat2,///< second matrix to be filled by element
    Epetra_SerialDenseVector& elevec1,  ///< third matrix to be filled by element
    Epetra_SerialDenseVector& elevec2,  ///< first vector to be filled by element
    Epetra_SerialDenseVector& elevec3  ///< second vector to be filled by element
  );


  /*!
  \brief Evaluate a Neumann boundary condition

  Since the element has no physics attached this method will give a dserror

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann
  (
    Teuchos::ParameterList& params,  ///< ParameterList for communication
    DRT::Discretization& discretization,  ///< discretization
    DRT::Condition& condition,  ///< Neumann condition to evaluate
    std::vector<int>& lm,  ///< location vector
    Epetra_SerialDenseVector& elevec1,   ///< vector to be filled by element
    Epetra_SerialDenseMatrix* elemat1 = NULL
  );


private:

  //! action parameters recognized by constraint element
  enum ActionType
  {
    none,
    calc_MPC_stiff,
    calc_MPC_state
  };

  //! container containing all kind of constraint element data
  DRT::Container                                  data_;

  //! vector of surfaces of this element (length 1)
  std::vector<DRT::Element*>                      surface_;

  // don't want = operator
  ConstraintElement3& operator = (const ConstraintElement3& old);

  /*!
  \brief Create matrix with material configuration for 3 dimensions and 4 nodes
  */
  inline void MaterialConfiguration
  (
    LINALG::Matrix<4,3>& x ///< nodal coords in material frame
  ) const
  {
    const int numnode = 4;
    const int numdim=3;
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
  \brief Create matrix with spatial configuration for 3 dimensions and 4 nodes
  */
  inline void SpatialConfiguration
  (
    LINALG::Matrix<4,3>& x, ///< nodal coords in spatial frame
    const std::vector<double> disp  ///< displacements
  ) const
  {
      const int numnode = 4;
      const int numdim = 3;
      for (int i=0; i<numnode; ++i)
      {
        for (int j = 0; j < numdim; ++j)
        {
          x(i,j) = Nodes()[i]->X()[j] + disp[i*numdim+j];
        }
      }
      return;
  }

  /*!
  \brief Create matrix with spatial configuration for 3 dimensions and 2 nodes
  */
  inline void SpatialConfiguration
  (
    LINALG::Matrix<2,3>& x, ///< nodal coords in spatial frame
    const std::vector<double> disp  ///< displacements
  ) const
  {
      const int numnode = 2;
      const int numdim = 3;
      for (int i=0; i<numnode; ++i)
      {
        for (int j = 0; j < numdim; ++j)
        {
          x(i,j) = Nodes()[i]->X()[j] + disp[i*numdim+j];
        }
      }
      return;
  }

  /// compute normal for 3D case using the first three nodes to specify a plane
  void ComputeNormal
  (
    const LINALG::Matrix<4,3>& xc,  ///< nodal coords in spatial frame
    LINALG::Matrix<3,1>& elenormal  ///< resulting element normal
  );

  /// Compute normal distance between plane and fourth node
  double ComputeNormalDist
  (
    const LINALG::Matrix<4,3>& xc,  ///< nodal coords in spatial frame
    const LINALG::Matrix<3,1>& elenormal  ///< element normal
  );

  /// Compute first derivative of normal distance with respect to the nodal displacements
  void ComputeFirstDeriv
  (
    const LINALG::Matrix<4,3>& xc,  ///< nodal coords in spatial frame
    Epetra_SerialDenseVector& elevector,  ///< vector to store results into
    const LINALG::Matrix<3,1>& elenormal  ///< element normal
  );

  /// Compute first derivative of normal distance with respect to the nodal displacements
  void ComputeSecondDeriv
  (
    const LINALG::Matrix<4,3>& xc,  ///< nodal coords in spatial frame
    Epetra_SerialDenseMatrix& elematrix,  ///< vector to store results into
    const LINALG::Matrix<3,1>& elenormal  ///< element normal
  );

  /// Compute difference of nodal displacement to masternode in given direction
  double ComputeWeightedDistance
  (
    const std::vector<double> disp,    ///< displacement vector of current node and masternode
    const std::vector<double> direct   ///< direction to weight with
  );

  /// Compute first derivatives nodal displacement to masternode in given direction
  void ComputeFirstDerivWeightedDistance
  (
    Epetra_SerialDenseVector& elevector,  ///< vector to store results into
    const std::vector<double> direct   ///< direction to weight with
  );

  /// Compute difference of spatial configuration to masternode in given direction
  double ComputeWeightedDistance
  (
    const LINALG::Matrix<2,3> x,
    const std::vector<double> direct
  );
}; // class ConstraintElement3


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================



} // namespace ELEMENTS
} // namespace DRT


#endif
