/*!----------------------------------------------------------------------
\file beam3k.H

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262


3D nonlinear Kirchhoff-like beam element. It can be switched between a variant with weak enforcement of the Kirchhoff
constraint and variant with strong enforcement of the Kirchhoff constraint. This variant with weak constraint enforcement is based
on a rotation interpolation that is similar to beam3r. As the beam curve has to be C^1-continous, it is interpolated
with Hermite polynomials of order 3. Therefore each of the two boundary nodes has 7 dofs. With the flag rotvec_ is can be
switched between to sets of degrees of freedom on the boundary node. The first set (rotvec_==true) is
[\v{d}_1, \v{theta}_1, t_1, \v{d}_2, \v{theta}_2, t_2, \alpha_3], where \v{d}_i is the vector of nodal positions
on the boundary nodes,\v{theta}_i is a pseudo rotation vector describing the nodal triad orientation on the boundary
nodes (and therewith also the orientation of the boundary tangent vectors), t_i is the length of the boudary
tangent vectors, and alpha_3 is the scalar relative rotation angle between reference and matrial triad at the
interior node. In contrary, the second variant (rotvec_==false) of this element has the dofs:
[\v{d}_1, \v{t}_1, \alpha_1, \v{d}_2, \v{t}_2, \alpha_2, \alpha_3], where \v{t}_i is the nodal tangent vector
(orientation and length) at the boundary nodes, and alpha_i are the scalar relative rotation angles between
reference and matrial triad at the boundary nodes and the interior node. Besides these two boundary nodes
the element has COLLOCATION_POINTS-2 interior nodes which one scalar DoF alpha_i, respectively.


Attention: Since so far linearizations are calculated with FAD, the rotation increments in the case (rotvec_==true)
are of an additive nature, which is in strong contrast to the beam3r implementation, where the iterative
rotation increments are multiplicative. Consequently, the inhomogeneous rotational Dirichlet conditions of beam3k
can be interpreted as additive increments added to the initial values (i.e. if the initial value is zero, the
Dirichlet values in the input file are the total nodal rotation angles). This is not true for beam3r,
where prescribed 3D rotation values have no direct physical interpretation. For 2D rotations both variants are identical.
*/

//header file only included if not yet included!
#ifndef BEAM3K_H
#define BEAM3K_H


#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/largerotations.H"
#include <Sacado.hpp>
#include "../headers/FAD_utils.H"

#include "beam3_base.H"

#include <unistd.h>
#include <sys/times.h>

typedef Sacado::Fad::DFad<double> FAD;

#define REFERENCE_NODE 2      //number of reference-node for the calculation of the rotation of the other triads (local numbering: 0 2 3 4 1)
//DEFAULT: REFERENCE_NODE 2   //REFERENCE_NODE=2 represents the midpoint node in the case of 3 COLLOCATION_POINTS -> standard choice!

#define mygaussrulebeam3k intrule_line_4point //define gauss rule; intrule_line_1point - intrule_line_10point is implemented.
//DEFAULT: intrule_line_4point

#define LENGTHCALCNEWTONTOL 10e-12 //Tolerance for iterative determination of element length.
//DEFAULT: 10e-12

#define COLLOCATION_POINTS 3      //!defines type of element. 2,3,4 are supported. A value of 3 or 4 means that further inner nodes are introduced
//DEFAULT: COLLOCATION_POINTS 3   //to interpolate the torsional degree of freedom alpha_. Furthermore, it specifies the number of collocation
                                  //points defining the number of material triads used to interpolate the triad field.

//#define CONSISTENTSPINSK          //Apply variationally consistent variant of first spin vector component as test function for strong Kirchhoff
                                  //-> this interpolation enables exact energy and angular momentum balances

//! \brief DRT: namespace of the ccarat discretization module
namespace DRT
{
namespace ELEMENTS
{

class Beam3kType : public DRT::ElementType
{
public:

  std::string Name() const { return "Beam3kType"; }

  static Beam3kType & Instance();

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

  static Beam3kType instance_;
};

/*!
\brief 3D nonlinear Kirchhoff-like beam element that can display initially curved beams.

author meier (meier@lnm.mw.tum.de)
*/
class Beam3k : public Beam3Base
{
public:

  //! @name Friends
  friend class Beam3kType;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3k(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3k(const Beam3k& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3k();

  /*!
  \brief Deep copy this instance of Beam3k and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return Beam3kType::Instance().UniqueParObjectId(); }

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
  { return Beam3kType::Instance(); }

  //@}

  /*!
   \brief Get bool inidcating usage of rotation vectors
   */
  const bool& RotVec() const
  {
    return rotvec_;
  }

  /*!
  \brief get reference rotation vector i.e. theta0_
  */
  std::vector<LINALG::Matrix<3,1> > Theta0() const
  {
    return theta0_;
  }

  /*!
  \brief get moment of inertia Iyy
  */
  const double& Iyy() const
  {
    return Iyy_;
  }

  /*!
   \brief Get centerline position at given parameter coordinate xi
  */
  LINALG::Matrix<3,1> GetPos(double& xi, LINALG::Matrix<12,1>& disp_totlag) const;

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 1;}


  /*!
  \brief Get vector of RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();


  /*!
  \brief Get number of degrees of freedom of a single node
  */
  virtual int NumDofPerNode(const DRT::Node & node) const
  {
    /*note: this is not necessarily the number of DOF assigned to this node by the discretization finally, but only
    *the number of DOF requested for this node by this element; the discretization will finally assign the maximal
    *number of DOF to this node requested by any element connected to this node*/

    int dofpn_aux=0;

    if (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id())
    {
      dofpn_aux = 7;
    }
    else
    {
      dofpn_aux = 1;
    }

    const int dofpn = dofpn_aux;
    return dofpn;
  }

  /*!
  \brief Get number of degrees of freedom per element not including nodal degrees of freedom
  */
  virtual int NumDofPerElement() const { return 0;}

  /*!
  \brief Print this element
  */
  virtual void Print(std::ostream& os) const;

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
                              std::vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1);


  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration)
  void SetUpReferenceGeometry(const std::vector<LINALG::Matrix<3,1> >& xrefe,const bool secondinit = false);

  /*!
  \brief get nodal tangents in reference configuration
  */
  std::vector<LINALG::Matrix<3,1> > Tref() const
  {
    return Tref_;
  }
  //@}

private:

  //Variables

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

  //!quaternion describing the nodal reference triads (for the case BEAM3EK_ROT ==false) of the last time step
  std::vector<LINALG::Matrix<4,1> > qrefconv_;
  //!tangent vectors at the two boundary nodes in the reference configuration
  std::vector<LINALG::Matrix<3,1> > Tref_;
  //Variables needed by all local triads
  //! Matrix holding pseudo rotation vectors describing the material triads in the reference configuration at each collocation point (CP)
  std::vector<LINALG::Matrix<3,1> > theta0_;
  //! Matrix with the material curvature in the initial configuration at each gp
  std::vector<LINALG::Matrix<3,1> > K0_;

  //! Length of the element
  double length_;
  //! Jacobi determinant of the element at the Gauss points
  std::vector<double> jacobi_;
  //! Additional jacobi factor appearing in the second derivatives (required in cae of strong Kirchhoff constraint)
  std::vector<double> jacobi2_;

  //! Jacobi determinant at collocation points
  std::vector<double> jacobi_cp_;

  //! bool indicating, if the DoFs at the element boundary are described by rotation vectors or tangent vectors plus relative angle
  bool rotvec_;

  //! bool indicating, if Kirchhoff constraint is enforced in a strong (weakkirchhoff_==false) or weak manner (weakkirchhoff_==true)
  bool weakkirchhoff_;

  //! bool recognizing first element call
  bool firstcall_;

  //! internal energy
  double Eint_;

  //******************************Begin: Class variables required for time integration**************************************************//
  //! triads at Gauss points for exact integration in quaternion at the end of the preceeding time step (required for computation of angular velocity)
  std::vector<LINALG::Matrix<4,1> > Qconvmass_;
  //! current triads at Gauss points for exact integration in quaternion (required for computation of angular velocity)
  std::vector<LINALG::Matrix<4,1> > Qnewmass_;
  //! spatial angular velocity vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > wconvmass_;
  //! current spatial angular velocity vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > wnewmass_;
  //! spatial angular acceleration vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > aconvmass_;
  //! current spatial angular acceleration vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > anewmass_;
  //! modified spatial angular acceleration vector (according to gen-alpha time integration) at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > amodconvmass_;
  //! current modified spatial angular acceleration vector (according to gen-alpha time integration) at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > amodnewmass_;
  //! translational acceleration vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttconvmass_;
  //! current translational acceleration vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttnewmass_;
  //! modified translational acceleration vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttmodconvmass_;
  //! modified current translational acceleration vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttmodnewmass_;
  //! translational velocity vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rtconvmass_;
  //! current translational velocity vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rtnewmass_;
  //! translational displacement vector at Gauss points for exact integration at the end of the preceeding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rconvmass_;
  //! current translational displacement vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rnewmass_;

  //! scale factor for translational inertia terms
  double inertscaletrans_;
  //! scale factor for tangential component of rotational inertia terms
  double inertscalerot1_;
  //! scale factor for remaining components of rotational inertia terms
  double inertscalerot2_;
  //******************************End: Class variables required for time integration**************************************************//

  //Functions

  //! @name methods for initilization of the element

  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration) in case of a Weak Kirchhoff Constraint
  void SetUpReferenceGeometryWK(const std::vector<LINALG::Matrix<3,1> >& xrefe,const bool secondinit);

  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration) in case of a Strong Kirchhoff Constraint
  void SetUpReferenceGeometrySK(const std::vector<LINALG::Matrix<3,1> >& xrefe,const bool secondinit);

  //@}

  //! @name auxiliary computation methods for non-additive, large rotation variables

  //@}

  //! @name Internal calculation methods

  //! Calculate internal forces and stiffness matrix via FAD in case of a Weak Kirchhoff Constraint
  void CalculateInternalForcesWK( Teuchos::ParameterList& params,
                                std::vector<double>&      disp,
                                Epetra_SerialDenseMatrix* stiffmatrix,
                                Epetra_SerialDenseMatrix* massmatrix,
                                Epetra_SerialDenseVector* force,
                                Epetra_SerialDenseVector* inertia_force,
                                bool update);

  //! Calculate internal forces and stiffness matrix via FAD in case of a Strong Kirchhoff Constraint
  void CalculateInternalForcesSK(Teuchos::ParameterList& params,
                                std::vector<double>&      disp,
                                Epetra_SerialDenseMatrix* stiffmatrix,
                                Epetra_SerialDenseMatrix* massmatrix,
                                Epetra_SerialDenseVector* force,
                                Epetra_SerialDenseVector* inertia_force,
                                bool update);

  //! Calculate internal forces and stiffness matrix via FAD in case of a Weak Kirchhoff Constraint
  void CalculateInteriaForces(Teuchos::ParameterList& params,
                              std::vector<LINALG::TMatrix<FAD,3,3> >& triad_mat,
                              std::vector<FAD>& disp_totlag_centerline,
                              std::vector<LINALG::TMatrix<FAD,6*2+COLLOCATION_POINTS,3> >& v_theta,
                              LINALG::TMatrix<FAD,6*2+COLLOCATION_POINTS,1>& f_inert);

   //!compute the relative angle theta between triad_ref and triad so that exp(theta) = triad_ref^T*triad -> inversion of the right translation triad=triad_ref*exp(theta)
  template<typename T>
  void triadtoangleright(LINALG::TMatrix<T,3,1>& theta,LINALG::TMatrix<T,3,3>& triad_ref,LINALG::TMatrix<T,3,3>& triad)
   {
     LINALG::TMatrix<T,3,3> rotation_matrix;
     LINALG::TMatrix<T,4,1> quaternion;

     rotation_matrix.Clear();
     quaternion.Clear();

     rotation_matrix.MultiplyTN(triad_ref,triad);
     LARGEROTATIONS::triadtoquaternion(rotation_matrix,quaternion);
     LARGEROTATIONS::quaterniontoangle(quaternion,theta);
   }

  //!compute the relative angle theta between triad_ref and triad so that exp(theta) = triad*triad_ref^T -> inversion of the right translation triad=exp(theta)*triad_ref
  template<typename T>
  void triadtoangleleft(LINALG::TMatrix<T,3,1>& theta,LINALG::TMatrix<T,3,3>& triad_ref,LINALG::TMatrix<T,3,3>& triad)
  {
    LINALG::TMatrix<T,3,3> rotation_matrix;
    LINALG::TMatrix<T,4,1> quaternion;

    rotation_matrix.Clear();
    quaternion.Clear();

    rotation_matrix.MultiplyNT(triad,triad_ref);
    LARGEROTATIONS::triadtoquaternion(rotation_matrix,quaternion);
    LARGEROTATIONS::quaterniontoangle(quaternion,theta);
  }

   //!compute a triad with an (interpolated) relative angle between triad_ref and triad
   template<typename T>
   void angletotriad(LINALG::TMatrix<T,3,1>& theta,LINALG::TMatrix<T,3,3>& triad_ref,LINALG::TMatrix<T,3,3>& triad)
   {
     LINALG::TMatrix<T,3,3> rotation_matrix;
     LINALG::TMatrix<T,4,1> quaternion;

     rotation_matrix.Clear();
     quaternion.Clear();
     triad.Clear();

     LARGEROTATIONS::angletoquaternion(theta,quaternion);
     LARGEROTATIONS::quaterniontotriad(quaternion,rotation_matrix);
     triad.Multiply(triad_ref,rotation_matrix);// ???
   }

   //!compute (material) strain K
   template<typename T>
   void computestrain(LINALG::TMatrix<T,3,1>& theta, LINALG::TMatrix<T,3,1>& theta_deriv, LINALG::TMatrix<T,3,1>& K)
   {
     LINALG::TMatrix<T,3,3> Tinv;
     K.Clear();
     Tinv.Clear();

     Tinv=LARGEROTATIONS::Tinvmatrix(theta);
     K.MultiplyTN(Tinv,theta_deriv);

     return;
   }

   //! calculate material stress resultants M,N from material strain resultants K, epsilon
   void straintostress(LINALG::TMatrix<FAD,3,1>& Omega, FAD epsilon, LINALG::TMatrix<FAD,3,1>& M, FAD& f_par);

   //! Rotate an arbitrary triad around its first base vector (tangent)
   template<typename T>
   void RotateTriad(LINALG::TMatrix<T,3,3>& triad, T& alpha, LINALG::TMatrix<T,3,3>& triad_rot)
   {
     for (int i=0; i<3; i++)
     {
     triad_rot(i,0)=triad(i,0);
     triad_rot(i,1)=triad(i,1)*cos(alpha)+triad(i,2)*sin(alpha);
     triad_rot(i,2)=triad(i,2)*cos(alpha)-triad(i,1)*sin(alpha);
     }

     return;
   }

  //! Calculate the SR mapping for a given reference system triad_ref and a given tangent vector r_s
  template<typename T>
  void CalculateSRTriads(LINALG::TMatrix<T,3,1>& r_s, LINALG::TMatrix<T,3,3>& triad_ref,LINALG::TMatrix<T,3,3>& triad)
  {

    //In this calculation, r_s does not necessarily have to be a unit vector.
    T temp_scalar1=0.0;
    T temp_scalar2=0.0;
    T temp_scalar3=0.0;
    T fac_n0=0.0;
    T fac_b0=0.0;
    T abs_r_s = FADUTILS::Norm<T>(r_s);

    for (int i=0; i<3; i++)
    {
      temp_scalar1+=triad_ref(i,1)*r_s(i);
      temp_scalar2+=triad_ref(i,2)*r_s(i);
      temp_scalar3+=triad_ref(i,0)*r_s(i);
    }

    fac_n0=temp_scalar1/(abs_r_s+temp_scalar3);
    fac_b0=temp_scalar2/(abs_r_s+temp_scalar3);

    for (int i=0; i<3; i++)
    {
      triad(i,0)=r_s(i)/abs_r_s;
      triad(i,1)=triad_ref(i,1)-fac_n0*(r_s(i)/abs_r_s+triad_ref(i,0));
      triad(i,2)=triad_ref(i,2)-fac_b0*(r_s(i)/abs_r_s+triad_ref(i,0));
     }

    return;
  }

  //! Compute the material triad in case of the strong Kirchhoff (SK)beam formulation
  template<typename T>
  void ComputeTriadSK(T& phi, LINALG::TMatrix<T,3,1>& r_s, LINALG::TMatrix<T,3,3>& triad_ref,LINALG::TMatrix<T,3,3>& triad)
  {
    LINALG::TMatrix<T,3,3> triad_bar(true);
    //Compute triad_bar via SR mapping from triad_ref onto r_s
    CalculateSRTriads(r_s,triad_ref,triad_bar);
    //Compute triad via relative rotation of triad_bar
    RotateTriad(triad_bar, phi, triad);
  }

  void AssembleShapefunctionsNss( LINALG::TMatrix<FAD,1,4>& N_i_xi,
                                  LINALG::TMatrix<FAD,1,4>& N_i_xixi,
                                  FAD jacobi,
                                  FAD jacobi2,
                                  LINALG::TMatrix<FAD,3,2*6+COLLOCATION_POINTS>& N_ss);

  void AssembleShapefunctionsNs(LINALG::TMatrix<FAD,1,4>& N_i_xi,
                                FAD jacobi,
                                LINALG::TMatrix<FAD,3,2*6+COLLOCATION_POINTS>& N_s);

  void AssembleShapefunctionsL(LINALG::TMatrix<FAD,1,COLLOCATION_POINTS>& L_i,
                               LINALG::TMatrix<FAD,1,2*6+COLLOCATION_POINTS>& L);

  void AssembleShapefunctionsN( LINALG::TMatrix<FAD,1,4>& N_i,
                                LINALG::TMatrix<FAD,3,2*6+COLLOCATION_POINTS>& N);

  //! Delivers the tangent vector of the reference configuration at the boundary node "node=1,2" (necessary to calculate the centerline of the reference geometry)
  LINALG::Matrix<3,1> GetNodalRefTangent(const int& node)
  {
    if(node!=0 and node!=1)
      dserror("Method can only deliver tangents for boundary nodes with IDs 0 or 1!!!");

    LINALG::Matrix<3,1> nodal_tangent(true);
    LINALG::Matrix<3,3> nodal_triad(true);
    LARGEROTATIONS::angletotriad(theta0_[node],nodal_triad);
    for(int i=0;i<3;i++)
      nodal_tangent(i)=nodal_triad(i,0);

    return nodal_tangent;
  }

  //! Calculates nodal positions, tangents and relative twist angles out of the corresponding displacements
  void UpdateDispTotlag(std::vector<double>& disp, std::vector<FAD>& disp_totlag);

  //!Set positions vectors and tangents at boundary nodes and triads at all CPs
  void SetNodalVariables(std::vector<FAD>& disp_totlag,
                         std::vector<FAD>& disp_totlag_centerline,
                         std::vector<LINALG::TMatrix<FAD,3,3> >& triad_mat_cp);


  //!Pre-multiply trafo matrix if rotvec_==true: \tilde{\vec{f}_int}=\mat{T}^T*\vec{f}_int
  void ApplyRotVecTrafo( std::vector<FAD>& disp_totlag_centerline,
                         LINALG::TMatrix<FAD,6*2+COLLOCATION_POINTS,1>& f_int);

  //!Update the reference triads qrefconv_
  void UpdateRefTriads(std::vector<FAD>& disp_totlag_centerline);

  //!Update the triads at the end of the time-step
  void UpdateTriads( Teuchos::ParameterList& params,
                     std::vector<double>& vel,
                     std::vector<double>& disp);

  //! lump mass matrix
  void Lumpmass(Epetra_SerialDenseMatrix* emass);

  //! Computes the length of the element and its Jacobi based upon a hermitian interpolation of the geometry
  void Calculate_length(const std::vector<LINALG::Matrix<3,1> >& xrefe, const std::vector<LINALG::Matrix<3,1> >& trefe, double tolerance);

  template<typename T>
  void Calculate_clcurvature(LINALG::TMatrix<T,3,1>& r_s, LINALG::TMatrix<T,3,1>& r_ss, LINALG::TMatrix<T,3,1>& kappacl)
  {
    //spinmatrix Sr' = r'x
    LINALG::TMatrix<T,3,3> Srs(true);
    LARGEROTATIONS::computespin(Srs,r_s);

    //cross-product r'xr''
    LINALG::TMatrix<T,3,1> Srsrss(true);
    Srsrss.Multiply(Srs,r_ss);
    T rsTrs=0.0;

    for (int i=0;i<3;i++)
      rsTrs+=r_s(i)*r_s(i);

    for (int i=0;i<3;i++)
    {
      kappacl(i)=Srsrss(i)/rsTrs;
    }

    return;
  }

  template<typename T>
  void computestrainSK(T& phi_s, LINALG::TMatrix<T,3,1>& kappacl, LINALG::TMatrix<T,3,3>& triadref,LINALG::TMatrix<T,3,3>& triad_mat,LINALG::TMatrix<T,3,1>& K)
  {
    LINALG::TMatrix<T,1,1> scalar_aux(true);
    LINALG::TMatrix<T,3,1> g1(true);
    LINALG::TMatrix<T,3,1> g2(true);
    LINALG::TMatrix<T,3,1> g3(true);
    LINALG::TMatrix<T,3,1> gref1(true);
    T KR1=0.0;

    for (int i=0;i<3;i++)
    {
      g1(i)=triad_mat(i,0);
      g2(i)=triad_mat(i,1);
      g3(i)=triad_mat(i,2);
      gref1(i)=triadref(i,0);
    }

    scalar_aux.MultiplyTN(kappacl,gref1);
    KR1=-scalar_aux(0,0);
    scalar_aux.Clear();
    scalar_aux.MultiplyTN(g1,gref1);
    KR1=KR1/(1.0+scalar_aux(0,0));
    K(0)=KR1+phi_s;

    scalar_aux.Clear();
    scalar_aux.MultiplyTN(kappacl,g2);
    K(1)=scalar_aux(0,0);

    scalar_aux.Clear();
    scalar_aux.MultiplyTN(kappacl,g3);
    K(2)=scalar_aux(0,0);

    return;
  }

  void ResizeClassVariables(const int& n)
  {
    qrefconv_.resize(COLLOCATION_POINTS);
    K0_.resize(n);
    jacobi_.resize(n);
    jacobi2_.resize(n);
    jacobi_cp_.resize(COLLOCATION_POINTS);
    Qconvmass_.resize(n);
    Qnewmass_.resize(n);
    wconvmass_.resize(n);
    wnewmass_.resize(n);
    aconvmass_.resize(n);
    anewmass_.resize(n);
    amodconvmass_.resize(n);
    amodnewmass_.resize(n);
    rttconvmass_.resize(n);
    rttnewmass_.resize(n);
    rttmodconvmass_.resize(n);
    rttmodnewmass_.resize(n);
    rtconvmass_.resize(n);
    rtnewmass_.resize(n);
    rconvmass_.resize(n);
    rnewmass_.resize(n);

    return;
  }

  void SetInitialDynamicClassVariables(const int& num, const LINALG::Matrix<3,3>& triad_mat, const LINALG::Matrix<3,1>& r)
  {
    Qconvmass_[num].Clear();
    Qnewmass_[num].Clear();
    rconvmass_[num].Clear();
    rnewmass_[num].Clear();
    wconvmass_[num].Clear();
    wnewmass_[num].Clear();
    aconvmass_[num].Clear();
    anewmass_[num].Clear();
    amodconvmass_[num].Clear();
    amodnewmass_[num].Clear();
    rtconvmass_[num].Clear();
    rtnewmass_[num].Clear();
    rttconvmass_[num].Clear();
    rttnewmass_[num].Clear();
    rttmodconvmass_[num].Clear();
    rttmodnewmass_[num].Clear();

    LARGEROTATIONS::triadtoquaternion(triad_mat,Qconvmass_[num]);
    Qnewmass_[num]=Qconvmass_[num];
    rconvmass_[num]=r;
    rnewmass_[num]=r;

    return;
  }

  // don't want = operator
  Beam3k& operator = (const Beam3k& old);

  template<int dim>
  void ComputeTripleProduct(LINALG::TMatrix<FAD,3,dim>& mat1, LINALG::TMatrix<FAD,3,1>& vec1, LINALG::TMatrix<FAD,3,1>& vec2, LINALG::TMatrix<FAD,dim,1>& vec_out)
  {
    LINALG::TMatrix<FAD,3,3> auxmatrix1(true);
    LINALG::TMatrix<FAD,3,1> auxvec1(true);
    LARGEROTATIONS::computespin(auxmatrix1,vec1);
    auxvec1.Multiply(auxmatrix1,vec2);
    vec_out.MultiplyTN(mat1,auxvec1);

    return;
  }

}; // class Beam3k

// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);

} // namespace Elements
} // namespace DRT

#endif  // #ifndef BEAM3K_H
