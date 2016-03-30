/*!----------------------------------------------------------------------
\file bele3.H

\maintainer Martin Pfaller

\brief dummy 3D boundary element without any physics
 *----------------------------------------------------------------------*/

#ifndef BELE3_H
#define BELE3_H


#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_element.H"
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

class Bele3Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Bele3Type"; }

  static Bele3Type& Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Bele3Type instance_;
};

/*!
 * A 3D boundary element with no physics attached
 *
 * This element is meant to have no physics. It can be used to have a boundary discretization
 * of surface/boundary elements. They can be of any 2d shape (quad4,quad9,tri3,...)
 *
 * The number of dof per node is set to 3 per default, so we can define displacement vectors by
 * using FillComplete on the boundary discretization. Furthermore numdofpernode can be adapted
 * if necessary.
 *
 * \author kruse (kruse@lnm.mw.tum.de)
 */
class Bele3 : public DRT::Element
{
  // friend class to fill number of dofs per node exclusively during creation
  friend class Bele3Type;

public:

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor
  */
  explicit Bele3(
      int id,     ///< A unique global id
      int owner   ///< proc num that owns this element
      );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  explicit Bele3(const Bele3& old);

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
  virtual int NumSurface() const {return 1;}
  virtual int NumVolume() const {return -1;}
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();
  virtual int UniqueParObjectId() const { return Bele3Type::Instance().UniqueParObjectId(); }
  virtual void Pack(DRT::PackBuffer& data) const;
  virtual void Unpack(const std::vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~Bele3();

  //@}

  //! @name Access methods

  virtual int NumDofPerNode(const DRT::Node&) const {return numdofpernode_;}
  virtual int NumDofPerElement() const { return 0; }
  void Print(std::ostream& os) const;
  virtual DRT::ElementType & ElementType() const { return Bele3Type::Instance(); }

  //@}

  //! @name Evaluation

  int Evaluate(Teuchos::ParameterList& params,
               DRT::Discretization&      discretization,
               std::vector<int>&         lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /// Read input for this element
  virtual bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);
  //@}

  //! @name Other
  //! does this element have non-zero displacements or not
//  bool IsMoving() const { return is_moving_; }

  //@}


private:

  /*!
    \brief Set number of dofs

    \param numdofpernode: number of degress of freedom for one node
   */
  virtual void SetNumDofPerNode(int numdofpernode) { numdofpernode_ = numdofpernode; }

  int numdofpernode_; ///< number of degrees of freedom

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
            break;
      }
      return hoel;
    };

  /*!
    \brief Create matrix with spatial configuration

    \param x     (out)  : nodal coords in spatial frame
    \param disp  (int)  : displacements
  */
  inline void SpatialConfiguration(LINALG::SerialDenseMatrix& x,
                                     const std::vector<double> disp) const
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

  //! Submethod to compute the enclosed volume for volume constraint boundary condition
  double ComputeConstrVols
  (
      const LINALG::SerialDenseMatrix& xscurr, ///< current configuration
      const int numnode ///< num nodes
  );

  //! vector with line elements
//  std::vector<Teuchos::RCP<DRT::Element> >                      lines_;

  //! flag for fixed or moving boundary
//  const bool                                      is_moving_;

  //! don't want = operator
  Bele3& operator = (const Bele3& old);

  //! set number of gauss points to element shape default
  DRT::UTILS::GaussRule2D getOptimalGaussrule(const DiscretizationType& distype) const;

}; // class Bele3



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

class Bele3LineType : public DRT::ElementType
{
public:

  std::string Name() const { return "Bele3LineType"; }

  static Bele3LineType & Instance();

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static Bele3LineType instance_;
};


/*!
\brief An element representing a line of a bele3 element

\author gerstenberger (gerstenberger@lnm.mw.tum.de)
*/
class Bele3Line : public DRT::FaceElement
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner: Processor owning this line
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  \param parent: The parent fluid element of this line
  \param lline: the local line number of this line w.r.t. the parent element
  */
  Bele3Line(int id, int owner, int nnode, const int* nodeids,
             DRT::Node** nodes, DRT::ELEMENTS::Bele3* parent, const int lline);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  Bele3Line(const Bele3Line& old);

  DRT::Element* Clone() const;
  virtual DiscretizationType Shape() const;
  virtual int UniqueParObjectId() const { return Bele3LineType::Instance().UniqueParObjectId(); }
  virtual void Pack(DRT::PackBuffer& data) const;
  virtual void Unpack(const std::vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~Bele3Line();

  //@}

  //! @name Access methods


  /*!
  \brief Get number of degrees of freedom of a certain node
         (implements pure virtual DRT::Element)

  For this 3D boundary element, we have 3 displacements, if needed
  */
  virtual int NumDofPerNode(const DRT::Node&) const { return numdofpernode_; }

  virtual int NumDofPerElement() const { return 0; }

  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return Bele3LineType::Instance(); }

  //@}

    //! @name Evaluation

  virtual int Evaluate(Teuchos::ParameterList& params,
               DRT::Discretization&      discretization,
               std::vector<int>&         lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  //! @name Evaluate methods

  virtual int EvaluateNeumann(Teuchos::ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  //@}

private:

  /*!
    \brief Set number of dofs

    \param numdofpernode: number of degress of freedom for one node
   */
  virtual void SetNumDofPerNode(int numdofpernode) { numdofpernode_ = numdofpernode; }

  int numdofpernode_; ///< number of degrees of freedom

  //! action parameters recognized by Bele3Line
  enum ActionType
  {
     none,
     integrate_Shapefunction
  };

  //! don't want = operator
  Bele3Line& operator = (const Bele3Line& old);


  //! compute infintesimal line element dr for integration along the line
  double  f2_substitution(
          const Epetra_SerialDenseMatrix  xyze,
          const Epetra_SerialDenseMatrix  deriv,
          const int iel);

  //! Get Rule for Gaussintegration according to DRT::UTIL
  DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

  //! integrate shape functions over a line
  void IntegrateShapeFunction(
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      const std::vector<int>&    lm,
      Epetra_SerialDenseVector&  elevec1,
      const std::vector<double>& edispnp
      );


}; // class Bele3Line




} // namespace ELEMENTS
} // namespace DRT


#endif
