/*!
 * \file bele3_4.H
 *
 * <pre>
 * Maintainer: Benedikt Schott
 *             schott@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15241
 * </pre>
 */

#ifndef BELE3_4_H
#define BELE3_4_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"


namespace DRT
{
// forward declarations
class Discretization;


namespace ELEMENTS
{


class Bele3_4Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Bele3_4Type"; }

  static Bele3_4Type & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const string eletype,
                                             const string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

private:

  static Bele3_4Type instance_;
};



/*!
 * A 3D boundary element with no physics attached
 *
 * This element is meant to have no physics. It can be used to have a boundary discretization
 * of surface/boundary elements. They can be of any 2d shape (quad4,quad9,tri3,...)
 *
 * The number of dof per node is set to 4 !!!, so we can define displacement vectors (using the first three components) and (fluid,pressure) vectors by
 * using FillComplete on the boundary discretization.
 *
 * \author schott (schott@lnm.mw.tum.de)
 */

class Bele3_4 : public DRT::Element
{
public:

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor
  */
  explicit Bele3_4(
      int id,     ///< A unique global id
      int owner   ///< proc num that owns this element
      );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  explicit Bele3_4(const Bele3_4& old);

  DRT::Element* Clone() const;
  virtual DiscretizationType Shape() const;
  virtual int NumLine() const
  {
    if (NumNode()==9 || NumNode()==8 || NumNode()==4) return 4;
    else if (NumNode()==3 || NumNode()==6) return 3;
    else
    {
      dserror("Could not determine number of lines");
      return -1;
    }
  }

  virtual int UniqueParObjectId() const { return Bele3_4Type::Instance().UniqueParObjectId(); }
  virtual void Pack(DRT::PackBuffer& data) const;
  virtual void Unpack(const std::vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~Bele3_4();

  //@}

  //! @name Access methods

  virtual int NumDofPerNode(const DRT::Node&) const {return 4;}
  virtual int NumDofPerElement() const { return 0; }
  void Print(ostream& os) const;
  virtual DRT::ElementType & ElementType() const { return Bele3_4Type::Instance(); }

  //@}

  //! @name Evaluation

  int Evaluate(ParameterList& params,
               DRT::Discretization&      discretization,
               vector<int>&              lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);


  //@}

  //! @name Other
  //! does this element have non-zero displacements or not
//  bool IsMoving() const { return is_moving_; }

  //@}


private:

  //! action parameters recognized by bele3
  enum ActionType
  {
    none
  };

  /*!
   * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
   * \return boolean indicating higher order status
   */
  bool isHigherOrderElement(
    const DRT::Element::DiscretizationType  distype
    ) const
    {
      bool hoel = true;
      switch (distype)
      {
          case quad4: case quad8: case quad9: case tri6:
            hoel = true;
            break;
          case tri3:
            hoel = false;
            break;
          default:
            dserror("distype unknown!");
      }
      return hoel;
    };



  //! don't want = operator
  Bele3_4& operator = (const Bele3_4& old);

  //! set number of gauss points to element shape default
  DRT::UTILS::GaussRule2D getOptimalGaussrule(const DiscretizationType& distype) const;

}; // class Bele3_4



} // namespace ELEMENTS
} // namespace DRT


#endif
