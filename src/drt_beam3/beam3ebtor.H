/*----------------------------------------------------------------------------*/
/*!
\file beam3ebtor.H

\brief three dimensional nonlinear rod based on a C1 curve

\level 2

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------------*/

/* 3D nonlinear Euler-Bernoulli-like beam element based on chapter 5 of the diploma thesis "Development of a finite element for
 * nonlinear beams based on the formulas of Frenet-Serret" by Christoph Meier. In contrary to the beam3eb element this element can
 * also display torsion and is therefore a complete, geometrically exact Euler-Bernoulli beam element for initially (i.e. in the stress-free
 * configuration) straight beams. Due to the Euler-Bernoulli theory there is no shear deformation.
 * To be able to use this element correctly structural dynamic parameters need to be set to:
 *
 *   LOADLIN     Yes
 *
 * since due to this beam formulation external point loads are being linearized and have an effect on the stiffness matrix. For this reason a special
 * type of Neumann conditions, namely DESIGN POINT MOMENT EB CONDITIONS, are needed for this element.
 *
 * As the beam curve has to be C1 it is interpolated with hermitien polynomials of order 3. Therefore each node has 7 dofs: the position
 * vector of the node (3 dofs), the tangent vector to the curve at the node (3 dofs) and the nodal twist angle (1 dof). If Dirichlet BC are applied one has to make sure that
 * the flags 3 -6 refer to the tangent at the node (rotational defree of freedom). The flag of the tangent in beam centerline direction must
 * be set to 0, otherwise the axial tension at the boundary would be prescribed. The 7th flag is related to the twist angle. */

//header file only included if not yet included!
#ifndef BEAM3EBTOR_H
#define BEAM3EBTOR_H


#include "beam3_base.H"

#define ScaleFactorLinetor 1;
#define ScaleFactorColumntor 1;

#define ANS_BEA3EBTOR// Flag, to apply Assumed Natural Strain approach for axial tension epsilon
//Default: On

#define mygaussruleebtor intrule_line_4point //define gauss rule; intrule_line_1point - intrule_line_10point is implemented
//Default: 4

//! \brief DRT: namespace of the ccarat discretization module
namespace DRT
{

namespace ELEMENTS
{

class Beam3ebtorType : public DRT::ElementType
{
public:

  std::string Name() const { return "Beam3ebtorType"; }

  static Beam3ebtorType & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis);

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Beam3ebtorType instance_;
};

/*!
\brief 3D nonlinear Euler-Bernoulli-like beam element based on chapter 5 of the diploma thesis "Development of a finite element for
nonlinear beams based on the formulas of Frenet-Serret" by Christoph Meier. The current formulation is only able to display axial
tension and bending curvature based on the curve describing the centerline of an initially (i.e. stress free) straight beam.
There is no shear deformation and no torsion. The expansion of the model to a full Euler Bernoulli beam (inclusive torsion) is possible.

author meier (meier@lnm.mw.tum.de)
*/
class Beam3ebtor : public Beam3Base
{
public:

  //! @name Friends
  friend class Beam3ebtorType;


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3ebtor(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3ebtor(const Beam3ebtor& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3ebtor();

  /*!
  \brief Deep copy this instance of Beam3ebtor and return pointer to the copy

  The Clone() method is used by the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
.
  */
  virtual DRT::Element* Clone() const;

   /*!
  \brief Get shape type of element
  */
  virtual DiscretizationType Shape() const;


  /*!
  \brief Return unique ParObject id

  Every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H
  */
  virtual int UniqueParObjectId() const { return Beam3ebtorType::Instance().UniqueParObjectId(); }

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

  virtual DRT::ElementType & ElementType() const
  { return Beam3ebtorType::Instance(); }

  //@}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 1;}

  /** \brief get number of nodes used for centerline interpolation
   *
   *  \author grill
   *  \date 05/16 */
  virtual inline int NumCenterlineNodes() const { return 2; }

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();


  /*!
  \brief Get number of degrees of freedom of a single node
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    /*note: this is not necessarily the number of DOF assigned to this node by the discretization finally, but only
    *the number of DOF requested for this node by this element; the discretization will finally assign the maximal
    *number of DOF to this node requested by any element connected to this node*/
    return 7;
  }

  /*!
  \brief Get number of degrees of freedom per element not including nodal degrees of freedom
  */
  virtual int NumDofPerElement() const { return 0;}

  /** \brief return flag whether Hermite polynomials are applied for centerline interpolation
   *
   *  \author grill
   *  \date 05/16 */
  virtual inline bool HermiteCenterlineInterpolation() const { return true; }

  /*!
  \brief Print this element
  */
  virtual void Print(std::ostream& os) const;

  /** \brief get Jacobi factor ds/dxi(xi) at xi \in [-1;1]
   *
   *  \author grill
   *  \date 06/16 */
  inline virtual double GetJacobiFacAtXi(const double& xi) const
  {
    dserror("current implementation only valid for initially straight beams,"
            " where the Jacobi factor is constant for all xi");
    return jacobi_;
  }

  //@}

  //! @name Construction


  /*!
  \brief Read input for this element

  This class implements a dummy of this method that prints a warning and
  returns false. A derived class would read one line from the input file and
  store all necessary information.

  */
  //virtual bool ReadElement();

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);

  //@}


  //! @name Evaluation methods


  /*!
  \brief Evaluate an element

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates element matrices and
  vectors accoring to the command in params.

  \note This class implements a dummy of this method that prints a warning and
        returns false.

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elemat2 (out)      : matrix to be filled by element depending on commands
                              given in params
  \param elevec1 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec2 (out)      : vector to be filled by element depending on commands
                              given in params
  \param elevec3 (out)      : vector to be filled by element depending on commands
                              given in params
  \return 0 if successful, negative otherwise
  */
  virtual int Evaluate(Teuchos::ParameterList&   params,
                       DRT::Discretization&      discretization,
                       std::vector<int>&         lm,
                       Epetra_SerialDenseMatrix& elemat1,
                       Epetra_SerialDenseMatrix& elemat2,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseVector& elevec2,
                       Epetra_SerialDenseVector& elevec3);

  /*!
  \brief Evaluate a Neumann boundary condition

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
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief get moment of intertia of area Iyy_
  */
  double Iyy() const
  {
    return Iyy_;
  }

  /*!
  \brief get nodal tangents in reference configuration
  */
  std::vector<LINALG::Matrix<3,1> > Tref() const
  {
    return Tref_;
  }

  /** \brief get access to the reference length
   *
   *  \author grill
   *  \date 05/16 */
  virtual inline double GetRefLength() const { return 2.0*jacobi_; }

  /*!
   \brief Get centerline position at given parameter coordinate xi
  */
  LINALG::Matrix<3,1> GetPos(double& xi, LINALG::Matrix<12,1>& disp_totlag) const;


  //@}

  /** \brief add indices of those DOFs of a given node that are positions
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void PositionDofIndices(std::vector<int>& posdofs, const DRT::Node& node) const
  {
    posdofs.push_back(0);  posdofs.push_back(1);  posdofs.push_back(2);
    return;
  }

  /** \brief add indices of those DOFs of a given node that are tangents (in the case of Hermite interpolation)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void TangentDofIndices(std::vector<int>& tangdofs, const DRT::Node& node) const
  {
    tangdofs.push_back(3);  tangdofs.push_back(4);  tangdofs.push_back(5);
    return;
  }

  /** \brief add indices of those DOFs of a given node that are rotation DOFs (non-additive rotation vectors)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void RotationVecDofIndices(std::vector<int>& rotvecdofs, const DRT::Node& node) const
  {
    return;
  }

  /** \brief add indices of those DOFs of a given node that are 1D rotation DOFs
   *         (planar rotations are additive, e.g. in case of relative twist DOF of beam3k with rotvec=false)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void Rotation1DDofIndices(std::vector<int>& twistdofs, const DRT::Node& node) const
  {
    twistdofs.push_back(6);
    return;
  }

  /** \brief add indices of those DOFs of a given node that represent norm of tangent vector
   *         (additive, e.g. in case of beam3k with rotvec=true)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void TangentLengthDofIndices(std::vector<int>& twistdofs, const DRT::Node& node) const
  {
    return;
  }


private:


  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;

  //! cross section area
  double crosssec_;

  //!moment of inertia of area around y-axis(orthognoal to rod axis in yz-plane)
  double Iyy_;

  //!moment of inertia of area around z-axis(forming cartesian system with rod- and y-axis)
  double Izz_;

  //!polar moment of inertia of area
  double Irr_;

  //! Variable holding the constant Jacobi determinant jacobi of the element
  double jacobi_;

  //! Matrix holding tangents at each node in reference configuration
  std::vector<LINALG::Matrix<3,1> > Tref_;

  //! @name methods for initilization of the element

  //! Computes inner product of two FAD vectors
  double ScalarProduct(LINALG::Matrix<3,1> first_vector, LINALG::Matrix<3,1> second_vector);

  //! Computes outer product of two FAD vectors
  LINALG::Matrix<3,1> VectorProduct(LINALG::Matrix<3,1> first_vector, LINALG::Matrix<3,1> second_vector);

  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration)
  void SetUpReferenceGeometry(const std::vector<double>& xrefe,const bool secondinit = false);

  //@}

  //! @name auxiliary computation methods for non-additive, large rotation variables




  //@}

  //! @name Internal calculation methods

  //!calculation of nonlinear stiffness and mass matrix
  void eb_nlnstiffmass( Teuchos::ParameterList&   params,
                        std::vector<double>&      vel,
                        std::vector<double>&      disp,
                        Epetra_SerialDenseMatrix* stiffmatrix,
                        Epetra_SerialDenseMatrix* massmatrix,
                        Epetra_SerialDenseVector* force);

  //!check of stiffness matrix via FAD
  void FADCheckStiffMatrix(std::vector<double>& disp,
                           Epetra_SerialDenseMatrix* stiffmatrix,
                           Epetra_SerialDenseVector* force);

  //!check of stiffness terms stemming from Neumann Boundary Conditions via FAD
  void FADCheckNeumann(Teuchos::ParameterList& params,
                       DRT::Discretization& discretization,
                       DRT::Condition& condition,
                       std::vector<int>& lm,
                       Epetra_SerialDenseVector& elevec1,
                       Epetra_SerialDenseMatrix* elemat1);

  //! lump mass matrix
  void lumpmass(Epetra_SerialDenseMatrix* emass);


  // don't want = operator
  Beam3ebtor& operator = (const Beam3ebtor& old);


}; // class Beam3ebtor



// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);


} // namespace Elements
} // namespace DRT


#endif  // #ifndef BEAM3EBTOR_H
