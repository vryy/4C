/*----------------------------------------------------------------------------*/
/*!
\file beam3r.H

\brief 3D nonlinear Reissner beam element

\level 2

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------------*/

/* 3D nonlinear Reissner beam element of type II (according to "The interpolation
 * of rotations and its application to finite element models of geometrically exact
 * rods", Romero 2004)
 *
 * Attention: For this implementation, prescribed 3D rotation values have no
 * direct physical interpretation so far because DBC handling is always additive in BACI.
 * For 2D rotations, multiplicative and additive increments are identical and
 * rotations can be prescribed without problems. */

#ifndef BEAM3R_H
#define BEAM3R_H


#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "Sacado.hpp"
#include "../headers/FAD_utils.H"

#include "beam3_base.H"


//#define BEAM3RAUTOMATICDIFF  //default: off  //Decide whether automatic differentiation via Sacado is used or not for linearization of f_int

#ifdef BEAM3RAUTOMATICDIFF
  typedef Sacado::Fad::DFad<double> FADordouble;

#else
  typedef double FADordouble;
#endif

#define BEAM3RLENGTHCALCNEWTONTOL 10e-12 //Tolerance for iterative determination of element length.
//DEFAULT: 10e-12

//#define BEAM3RCONSTSTOCHFORCE  //Flag in order to hold stochastic forces constant over the element length
                                //and to only provide random numbers for the 3 translational DoFs
                                //(needed in order to compare with beam3eb)

namespace DRT
{

namespace ELEMENTS
{

class Beam3rType : public DRT::ElementType
{
public:

  std::string Name() const { return "Beam3rType"; }

  static Beam3rType & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int         id,
                                             const int         owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis);

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Beam3rType instance_;
};

/*!
\brief 3D nonlinear Reissner beam element implemented according to the following sources:
Jelenic, Crisfield, 1999, "Geometrically exact 3D beam theory: implementations of a strain-invariant finite element for statics and dynamics",
Crisfield, Jelenic, 1999, "Objectivity of strain measures in the geometrically exact three dimensional beam theory and its finite element implementation",
Romero, 2004, "The interpolation of rotations and its application to finite element models of geometrically exact rods",
Crisfield, 2003, "Non-linear Finite Element Analysis of Solids and Structures", Volume 2

author cyron (cyron @lnm.mw.tum.de)
*/
class Beam3r : public Beam3Base
{
public:

  //! @name Friends
  friend class Beam3rType;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3r(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3r(const Beam3r& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3r();

  /*!
  \brief Deep copy this instance of Beam3r and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return Beam3rType::Instance().UniqueParObjectId(); }

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

  virtual DRT::ElementType & ElementType() const { return Beam3rType::Instance(); }

  //@}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 1;}

  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /** \brief Get number of nodes used for centerline interpolation
   *
   *  \author grill
   *  \date 05/16 */
  virtual inline int NumCenterlineNodes() const
  {
    return centerline_hermite_ ? 2 : this->NumNode();
  }

  /*!
  \brief Get number of degrees of freedom of a single node
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    /* note: this is not necessarily the number of DOF assigned to this node by the discretization finally, but only
     * the number of DOF requested for this node by this element; the discretization will finally assign the maximal
     * number of DOF to this node requested by any element connected to this node*/
    if (!centerline_hermite_)
      return 6;
    else
    {
      /* in case of Hermite centerline interpolation (so far always 3rd order = 2nodes), we have
       * 6 translational DOFs for the first two nodes and additionally 3 rotational DOFs for each node */
      int dofpn_aux=0;

      if (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id())
      {
        dofpn_aux = 9;
      }
      else
      {
        dofpn_aux = 3;
      }

      const int dofpn = dofpn_aux;
      return dofpn;
    }
  }

  /*!
  \brief Get number of degrees of freedom per element not including nodal degrees of freedom
  */
  virtual int NumDofPerElement() const { return 0;}

  /*!
  \brief Print this element
  */
  virtual void Print(std::ostream& os) const;

  /*!
  \brief Get moment of inertia of area Iyy_ of this element
  */
  const double& Iyy() const { return Iyy_;}

  /*!
  \brief Get moment of inertia of area Izz_ of this element
  */
  const double& Izz() const { return Izz_;}

  /*!
  \brief Get polar moment of inertia of area Irr_ of this element
  */
  const double& Irr() const { return Irr_;}

  /*!
  \brief return flag whether Hermite polynomials are applied for centerline interpolation
  */
  virtual inline bool HermiteCenterlineInterpolation() const { return centerline_hermite_;}

  /** \brief get centerline position at xi \in [-1,1] (element parameter space)
   *
   *  \author grill
   *  \date 06/16 */
  virtual void GetPosAtXi(LINALG::Matrix<3,1>&       pos,
                          const double&              xi,
                          const std::vector<double>& disp_totlag) const;

  /** \brief get triad at xi \in [-1,1] (element parameter space)
   *
   *  \author grill
   *  \date 07/16 */
  void GetTriadAtXi(LINALG::Matrix<3,3>&       triad,
                    const double&              xi,
                    const std::vector<double>& psi_totlag) const;

  /** \brief get unit tangent vector in reference configuration at i-th node of beam element (element-internal numbering)
   *
   *  \author grill
   *  \date 06/16 */
  virtual inline void GetRefTangentAtNode(LINALG::Matrix<3,1>& Tref_i,
                                          const int&           i) const
  {
    if (not ((unsigned)i<Tref().size()))
      dserror("asked for tangent at node index %d, but only %d centerline nodes existing",i,Tref().size());

    Tref_i = Tref()[i];
  }

  /*!
  \brief compute and get current tangent at specified node i.e. Tcurrnode_
  Attention: strictly speaking, this returns the first vector of the nodal triad which is NOT identical to the tangent of the centerline
             see also corresponding notes in implementation
  */
  LINALG::Matrix<3,1> Tcurr(const int NodeID);

  /*!
  \brief get tangent of centerline at first node in reference configuration (constant over element length for linear Lagrange interpolation)
  */
  LINALG::Matrix<3,1> Treffirst() const;

  /*!
  \brief get tangent of centerline at all nodes in reference configuration
  */
  std::vector<LINALG::Matrix<3,1> > Tref() const
  {
    return Trefnode_;
  }

  /*!
  \brief Get jacobiGP_ factor of first Gauss point for under-integration (constant over element length for linear Lagrange interpolation)
  */
  const double& GetJacobi() const
  {
    return jacobiGPelastf_[0];
  }

  /** \brief get Jacobi factor ds/dxi(xi) at xi \in [-1;1]
   *
   *  \author grill
   *  \date 06/16 */
  virtual double GetJacobiFacAtXi(const double& xi) const;

  /*!
  \brief Get maximal bending curvature occurring in this element
  */
  const double& GetKappaMax() const
  {
    return Kmax_;
  }

  /** \brief get access to the reference length
   *
   *  \author grill
   *  \date 05/16 */
  virtual inline double GetRefLength() const { return reflength_; }

  /*!
  \brief Get current nodal triads as quaternions
  */
  const std::vector<LINALG::Matrix<4,1> >& Qnew() const {return Qnewnode_;}

  /*!
  \brief Get initial nodal rotation vectors
  */
  const std::vector<LINALG::Matrix<3,1> >& InitialNodalRotVecs() const {return theta0node_;}

  /*!
  \brief Set moment of inertia of area Iyy_ of this element
  */
  void SetIyy(const double& Iyy);

  /*!
  \brief Set moment of inertia of area Iyy_ of this element
  */
  void SetIzz(const double& Izz);

  /*!
  \brief Set moment of inertia of area Iyy_ of this element
  */
  void SetIrr(const double& Irr);

  /*!
  \brief Set cross section area of this element
  */
  void SetCrossSec(const double& crosssec);

  /*!
  \brief Set cross section scaled with shear correction factor of this element
  */
  void SetCrossSecShear(const double& crosssecshear);

  //!computes the number of different random numbers required in each time step for generation of stochastic forces
  int HowManyRandomNumbersINeed();

  //! \brief internal force vector
  const Epetra_SerialDenseVector& InternalForceVector() const {return f_;}

  //! \brief return material forces at Gauss point
  const LINALG::Matrix<3,1>& MatForceGp() const {return Ngp_;}

  //! brief epsilonnormal
  const double& EpsilonSgn() const { return eps_;}

  //!sets up all geometric parameters (considering current position as reference configuration)
  /* nnodetriad: number of nodes used for interpolation of triad field
   * nnodecl: number of nodes used for interpolation of centerline
   * vpernode: interpolated values per centerline node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))*/
  template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
  void SetUpReferenceGeometry(const std::vector<double>& xrefe,
                              const std::vector<double>& rotrefe,
                              const bool secondinit = false);

  //@}

  //! @name Construction

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string&          eletype,
                   const std::string&          distype,
                   DRT::INPUT::LineDefinition* linedef);

  //@}


  //! @name Evaluation methods


  /*!
  \brief Evaluate this element

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

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              std::vector<int>&         lm,
                              Epetra_SerialDenseVector& elevec1,
                              Epetra_SerialDenseMatrix* elemat1 = NULL);

  /*!
  \brief Evaluate PTC addition to stiffness for free Brownian motion

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates a  statistical Neumann boundary condition
  as used in the problem type STATISTICAL MECHANICS

  \param params (in/out)       : ParameterList for communication between control routine and elements
  \param vector<double> mydisp : current nodal displacement
  \param elemat1 (out)         : artificial damping matrix to be filled by element

  \return 0 if successful, negative otherwise
  */
  template<unsigned int nnode>
  void EvaluatePTC(Teuchos::ParameterList&   params,
                   Epetra_SerialDenseMatrix& elemat1);

  //@}

  /** \brief add indices of those DOFs of a given node that are positions
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void PositionDofIndices(std::vector<int>& posdofs, const DRT::Node& node) const
  {
    if ( (not centerline_hermite_) or
         (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id()) )
    {
      posdofs.push_back(0);  posdofs.push_back(1);  posdofs.push_back(2);
    }
    return;
  }

  /** \brief add indices of those DOFs of a given node that are tangents (in the case of Hermite interpolation)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void TangentDofIndices(std::vector<int>& tangdofs, const DRT::Node& node) const
  {
    if ( centerline_hermite_ and
         (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id()) )
    {
      tangdofs.push_back(6);  tangdofs.push_back(7);  tangdofs.push_back(8);
    }
    return;
  }

  /** \brief add indices of those DOFs of a given node that are rotation DOFs (non-additive rotation vectors)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void RotationVecDofIndices(std::vector<int>& rotvecdofs, const DRT::Node& node) const
  {
    if ( (not centerline_hermite_) or
         (node.Id() == this->Nodes()[0]->Id() or node.Id() == this->Nodes()[1]->Id()) )
    {
      rotvecdofs.push_back(3);  rotvecdofs.push_back(4);  rotvecdofs.push_back(5);
    }
    else
    {
      rotvecdofs.push_back(0);  rotvecdofs.push_back(1);  rotvecdofs.push_back(2);
    }
    return;
  }

  /** \brief add indices of those DOFs of a given node that are 1D rotation DOFs
   *         (planar rotations are additive, e.g. in case of relative twist DOF of beam3k with rotvec=false)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void Rotation1DDofIndices(std::vector<int>& twistdofs, const DRT::Node& node) const
  {
    return;
  }

  /** \brief add indices of those DOFs of a given node that represent norm of tangent vector
   *         (additive, e.g. in case of beam3k with rotvec=true)
   *
   *  \author grill
   *  \date 07/16 */
  virtual inline void TangentLengthDofIndices(std::vector<int>& tangnormdofs, const DRT::Node& node) const
  {
    return;
  }

private:

  //! purposes of numerical integration
  enum IntegrationPurpose
  {
    res_elastic_force,
    res_elastic_moment,
    res_inertia,
    res_damp_stoch,
    neumann_lineload
  };

  //! bool storing whether Hermite interpolation of centerline is applied (false => Lagrange interpolation)
  bool centerline_hermite_;
  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;
  //! variable storing whether Statmech functionalities are needed or not
  bool needstatmech_;
  //! storing epsilon(0)    // TODO the following three quantities seem to be solely Statmech-related; check usage and validity; if necessary, adapt corresponding access functions
  double eps_;
  //! vector storing the internal force vector
  Epetra_SerialDenseVector f_;
  //! vector storing internal forces at gauss point in material configuration
  LINALG::Matrix<3,1> Ngp_;
  //! norm of maximal bending curvature occurring in this element
  double Kmax_;

//  //! rotation angle displacement at nodes; relates to the end of the preceding time step
  std::vector<LINALG::Matrix<3,1> > dispthetaconvnode_;
  //! rotation angle displacement at nodes; relates to current configuration
  std::vector<LINALG::Matrix<3,1> > dispthetanewnode_;
//  //! nodal triads in quaternion form at the end of the preceding time step
  std::vector<LINALG::Matrix<4,1> > Qconvnode_;
  //! nodal triads in quaternion during the current iteration step
  std::vector<LINALG::Matrix<4,1> > Qnewnode_;

  //************** begin: Class variables required for element-based Lie-group time integration *******************************
  //! triads at Gauss points for exact integration in quaternion at the end of the preceding time step (required for computation of angular velocity)
  std::vector<LINALG::Matrix<4,1> > QconvGPmass_;
  //! current triads at Gauss points for exact integration in quaternion (required for computation of angular velocity)
  std::vector<LINALG::Matrix<4,1> > QnewGPmass_;
  //! spatial angular velocity vector at Gauss points for exact integration at the end of the preceding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > wconvGPmass_;
  //! current spatial angular velocity vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > wnewGPmass_;
  //! spatial angular acceleration vector at Gauss points for exact integration at the end of the preceding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > aconvGPmass_;
  //! current spatial angular acceleration vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > anewGPmass_;
  //! translational acceleration vector at Gauss points for exact integration at the end of the preceding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttconvGPmass_;
  //! current translational acceleration vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttnewGPmass_;
  //! modified translational acceleration vector at Gauss points for exact integration at the end of the preceding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttmodconvGPmass_;
  //! modified current translational acceleration vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rttmodnewGPmass_;
  //! translational velocity vector at Gauss points for exact integration at the end of the preceding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rtconvGPmass_;
  //! current translational velocity vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rtnewGPmass_;
  //! translational displacement vector at Gauss points for exact integration at the end of the preceding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rconvGPmass_;
  //! current translational displacement vector at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > rnewGPmass_;
  //! modified spatial angular acceleration vector (according to gen-alpha time integration) at Gauss points for exact integration at the end of the preceding time step (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > amodconvGPmass_;
  //! current modified spatial angular acceleration vector (according to gen-alpha time integration) at Gauss points for exact integration (required for computation of inertia terms)
  std::vector<LINALG::Matrix<3,1> > amodnewGPmass_;
  //************** end: Class variables required for element-based Lie-group time integration *******************************//

  //! initial length of the element
  double reflength_;
  //! rotational pseudovectors at nodes in reference configuration
  std::vector<LINALG::Matrix<3,1> > theta0node_;
  //! vector holding current tangent at the centerline nodes
  std::vector<LINALG::Matrix<3,1> > Tcurrnode_;
  //! vector holding reference tangent at the centerline nodes
  std::vector<LINALG::Matrix<3,1> > Trefnode_;

  //! initial material curvature at Gauss points for elasticity (corresponding to \Lambda_0^t \Labmda'_0 in eq. (3.5), Crisfield 1999
  std::vector<LINALG::Matrix<3,1> > KrefGP_;

  //! initial axial tension (always zero) and shear deformation at Gauss points for elasticity (corresponding to \Lambda_0^t rprime_0 - (1,0,0) )
  std::vector<LINALG::Matrix<3,1> > GammarefGP_;

  //! node I for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
  unsigned int nodeI_;

  //! node J for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
  unsigned int nodeJ_;

  //! cross-section area
  double crosssec_;

   //! cross-section area scaled with shear correction factor
  double crosssecshear_;

  //! moment of inertia of area around y-axis (orthogonal to beam axis in yz-plane)
  double Iyy_;

  //! moment of inertia of area around z-axis (forming cartesian system with beam- and y-axis)
  double Izz_;

  //! polar moment of inertia of area
  double Irr_;

  //! Vector holding value of Jacobi determinant for each Gauss point of integration purpose res_elastic_force
  std::vector<double> jacobiGPelastf_;

  //! Vector holding value of Jacobi determinant for each Gauss point of integration purpose res_elastic_moment
  std::vector<double> jacobiGPelastm_;

  //! Vector holding value of Jacobi determinant for each Gauss point of integration purpose res_inertia
  std::vector<double> jacobiGPmass_;

  //! Vector holding value of Jacobi determinant for each Gauss point of integration purpose res_damp_stoch
  std::vector<double> jacobiGPdampstoch_;

  //! Vector holding value of Jacobi determinant for each Gauss point of integration purpose neumann_lineload
  std::vector<double> jacobiGPneumannline_;

  //! kinetic energy of element
  double Ekin_;

  //! internal (elastic) energy of element
  double Eint_;

  //! angular momentum of the element
  LINALG::Matrix<3,1> L_;

  //! linear momentum of the element
  LINALG::Matrix<3,1> P_;

  //! kinetic energy from rotational dofs part1
  double Ekintorsion_;

  //! kinetic energy from rotational dofs part2
  double Ekinbending_;

  //! kinetic energy from translational dofs
  double Ekintrans_;

  //! scale factor for translational inertia terms
  double inertscaletrans_;
  //! scale factor for tangential component of rotational inertia terms
  double inertscalerot1_;
  //! scale factor for remaining components of rotational inertia terms
  double inertscalerot2_;


  //! @name Internal calculation methods

  /*!
   \brief Get centerline position at given parameter coordinate xi
  */
  template<unsigned int nnodecl, unsigned int vpernode>
  LINALG::Matrix<3,1> GetPosAtXi(const double& xi,
                                 const LINALG::Matrix<3*nnodecl*vpernode,1>& disp_totlag) const
  {
    LINALG::Matrix<3,1> r;
    LINALG::Matrix<1,vpernode*nnodecl> N_i;

    this->EvaluateShapeFunctionsAtXi<nnodecl,vpernode>(xi,N_i,this->Shape());
    this->Calc_r<nnodecl,vpernode,double>(disp_totlag,N_i,r);

    return r;
  }

  /** \brief get Jacobi factor ds/dxi(xi) at xi \in [-1;1]
   *
   *  \author grill
   *  \date 08/16 */
  template<unsigned int nnodecl, unsigned int vpernode>
  double GetJacobiFacAtXi(const double& xi) const
  {
    /* ||dr_0/ds(xi)||=1 because s is arc-length parameter => ||dr_0/dxi(xi)|| * dxi/ds(xi) = 1
     * => JacobiFac(xi) = ds/dxi(xi) = ||dr_0/dxi(xi)|| */
    LINALG::Matrix<3,1> r0_xi;
    LINALG::Matrix<1,vpernode*nnodecl> N_i_xi;
    LINALG::Matrix<3*nnodecl*vpernode,1> disp_centerline_ref;

    // fill disp_centerline_ref with reference nodal centerline positions and tangents
    for (unsigned int node=0; node<nnodecl; node++)
    {
      for (unsigned int dim=0; dim<3; ++dim)
      {
        disp_centerline_ref(3*vpernode*node+dim) = Nodes()[node]->X()[dim];
        if (vpernode==2)
          disp_centerline_ref(3*vpernode*node+3+dim) = (Trefnode_[node])(dim);
      }
    }

    this->EvaluateShapeFunctionDerivsAtXi<nnodecl,vpernode>(xi,N_i_xi,this->Shape());
    this->Calc_r_xi<nnodecl,vpernode,double>(disp_centerline_ref,N_i_xi,r0_xi);

    return r0_xi.Norm2();
  }

  /*!
   \brief Get triad (three unit base vectors) at given parameter coordinate xi
  */
  template<unsigned int nnodetriad>
  LINALG::Matrix<3,3> GetTriadAtXi(const double& xi,
                                   const std::vector<LINALG::Matrix<3,1> >& psi_totlag) const
  {
    // nodal triads in form of quaternions
    std::vector<LINALG::Matrix<4,1> > Qnode(nnodetriad);
    for (unsigned int node=0; node<nnodetriad; node++)
    {
      // ToDo
      /* make function call independent from class variable Qnewnode_
         by using absolute coordinate vector psi_totlag (nodal rotation pseudo-vectors)
         passed into function from outside */

//      LARGEROTATIONS::angletoquaternion(psi_totlag[node],Qnode[node]);
      Qnode[node]=Qnewnode_[node];
    }

    // compute reference triad Lambda_r (in form of a quaternion) according to (3.9), Jelenic 1999
    LINALG::Matrix<4,1> Q_r(true);
    LINALG::Matrix<3,1> dummy(true);  // dummy 3x1 vector, quantity not required here
    CalcRefQuaternion<double>(Qnode[nodeI_],Qnode[nodeJ_],Q_r,dummy);

    // compute nodal rotation angles between nodal triads and reference triad according to (3.8), Jelenic 1999
    std::vector<LINALG::Matrix<3,1> > Psi_li(nnodetriad);
    for (unsigned int node=0; node<nnodetriad; node++)
    {
      CalcPsi_li<double>(Qnode[node],Q_r,Psi_li[node]);
    }

    // values of Lagrange shape function at xi used for interpolation
    LINALG::Matrix<1,nnodetriad> I_i;
    this->EvaluateShapeFunctionsAtXi<nnodetriad,1>(xi,I_i,this->Shape());

    // compute interpolated local relative rotation \Psi^l at xi according to (3.11), Jelenic 1999
    LINALG::Matrix<3,1> Psi_l;
    Calc_Psi_l<nnodetriad>(Psi_li, I_i, Psi_l);

    // finally compute interpolated triad Lambda at xi
    LINALG::Matrix<3,3> triad;
    Calc_Lambda<double>(Psi_l,Q_r,triad);

    return triad;
  }

  //! evaluate shape functions at position xi in element parameter space [-1,1]
  template<unsigned int nnode, unsigned int vpernode>
  void EvaluateShapeFunctionsAtXi(const double&                           xi,
                                  LINALG::Matrix<1,vpernode*nnode>&       I_i,
                                  const DRT::Element::DiscretizationType& distype) const
  {
    I_i.Clear();

    switch(vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape functions at xi
        DRT::UTILS::shape_function_1D(I_i,xi,distype);
        break;
      }
      case 2:
      {
        // evaluate Hermite shape functions at xi, for now: always 3rd order, i.e. line2
        DRT::UTILS::shape_function_hermite_1D(I_i,xi,reflength_,line2);
        break;
      }
      default:
        dserror("invalid value for vpernode (number of values per node) specified");
    }

    return;
  }

  //! evaluate shape function derivatives at position xi in element parameter space [-1,1]
  template<unsigned int nnode, unsigned int vpernode>
  void EvaluateShapeFunctionDerivsAtXi(const double&                           xi,
                                       LINALG::Matrix<1,vpernode*nnode>&       I_i_xi,
                                       const DRT::Element::DiscretizationType& distype) const
  {
    I_i_xi.Clear();

    switch(vpernode)
    {
      case 1:
      {
        // evaluate Lagrange shape functions at xi
        DRT::UTILS::shape_function_1D_deriv1(I_i_xi,xi,distype);
        break;
      }
      case 2:
      {
        // evaluate Hermite shape functions at xi, for now: always 3rd order, i.e. line2
        DRT::UTILS::shape_function_hermite_1D_deriv1(I_i_xi,xi,reflength_,line2);
        break;
      }
      default:
        dserror("invalid value for vpernode (number of values per node) specified");
    }

    return;
  }

  //! evaluate shape functions at specified Gauss points
  template<unsigned int nnode, unsigned int vpernode>
  void EvaluateShapeFunctionsAllGPs(const DRT::UTILS::IntegrationPoints1D&          gausspoints,
                                    std::vector<LINALG::Matrix<1,vpernode*nnode> >& I_i,
                                    const DRT::Element::DiscretizationType&         distype) const
  {
    if (I_i.size()!= (unsigned int) gausspoints.nquad)
      dserror("vector for %d individual shape functions has wrong size",gausspoints.nquad);

    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      // Get location of GP in element parameter space xi \in [-1;1]
      const double xi = gausspoints.qxg[numgp][0];
      this->EvaluateShapeFunctionsAtXi<nnode,vpernode>(xi,I_i[numgp],distype);
    }

    return;
  }

  //! evaluate shape function derivatives at specified Gauss points
  template<unsigned int nnode, unsigned int vpernode>
  void EvaluateShapeFunctionDerivsAllGPs(const DRT::UTILS::IntegrationPoints1D&          gausspoints,
                                         std::vector<LINALG::Matrix<1,vpernode*nnode> >& I_i_xi,
                                         const DRT::Element::DiscretizationType&         distype) const
  {
    if (I_i_xi.size()!= (unsigned int) gausspoints.nquad)
      dserror("vector for %d individual shape function derivatives has wrong size",gausspoints.nquad);

    for(int numgp=0; numgp < gausspoints.nquad; numgp++)
    {
      // Get location of GP in element parameter space xi \in [-1;1]
      const double xi = gausspoints.qxg[numgp][0];
      this->EvaluateShapeFunctionDerivsAtXi<nnode,vpernode>(xi,I_i_xi[numgp],distype);
    }

    return;
  }

  //! evaluate shape functions and derivatives at specified Gauss points
  template<unsigned int nnode, unsigned int vpernode>
  void EvaluateShapeFunctionsAndDerivsAllGPs(const DRT::UTILS::IntegrationPoints1D&          gausspoints,
                                             std::vector<LINALG::Matrix<1,vpernode*nnode> >& I_i,
                                             std::vector<LINALG::Matrix<1,vpernode*nnode> >& I_i_xi,
                                             const DRT::Element::DiscretizationType&         distype) const
  {
    EvaluateShapeFunctionsAllGPs<nnode,vpernode>(gausspoints,I_i,distype);
    EvaluateShapeFunctionDerivsAllGPs<nnode,vpernode>(gausspoints,I_i_xi,distype);

    return;
  }

  //! compute quaternion corresponding to reference triad Lambda_r according to (3.9), Jelenic 1999
  template<typename T>
  void CalcRefQuaternion(const LINALG::TMatrix<T,4,1>&  Q_nodeI,
                          const LINALG::TMatrix<T,4,1>& Q_nodeJ,
                          LINALG::TMatrix<T,4,1>&       Q_r,
                          LINALG::TMatrix<T,3,1>&       Phi_IJ) const
  {
    Q_r.Clear();
    Phi_IJ.Clear();

    // quaternion of relative rotation between node I and J according to (3.10), Jelenic 1999
    LINALG::TMatrix<T,4,1> QIJ(true);
    // quaternion of half relative rotation between node I and J according to (3.9), Jelenic 1999
    LINALG::TMatrix<T,4,1> QIJhalf(true);

    // compute quaternion Qr corresponding to reference triad \Lambda_r according to eq. (3.10) and (3.9), Jelenic 1999
    LARGEROTATIONS::quaternionproduct<T>(Q_nodeJ,LARGEROTATIONS::inversequaternion<T>(Q_nodeI),QIJ);
    LARGEROTATIONS::quaterniontoangle<T>(QIJ,Phi_IJ);

    LINALG::TMatrix<T,3,1> Phi_IJhalf(Phi_IJ);
    Phi_IJhalf.Scale(0.5);
    LARGEROTATIONS::angletoquaternion<T>(Phi_IJhalf,QIJhalf);

    LARGEROTATIONS::quaternionproduct<T>(QIJhalf,Q_nodeI,Q_r);

    return;
  }

  //! compute nodal local rotations according to (3.8), Jelenic 1999
  template<typename T>
  void CalcPsi_li(const LINALG::TMatrix<T,4,1>& Q_i,
                  const LINALG::TMatrix<T,4,1>& Q_r,
                  LINALG::TMatrix<T,3,1>&       Psi_li) const
  {
    // quaternion equivalent to Psi_li
    LINALG::TMatrix<T,4,1>  Q_li(true);

    LARGEROTATIONS::quaternionproduct<T>(Q_i,LARGEROTATIONS::inversequaternion<T>(Q_r),Q_li);
    LARGEROTATIONS::quaterniontoangle<T>(Q_li,Psi_li);

    return;
  }

  //! compute beam centerline position at element parameter \xi \in [-1;1]
  template<unsigned int nnode, unsigned int vpernode, typename T>
  void Calc_r(const LINALG::TMatrix<T,3*vpernode*nnode,1>& disp_totlag_centerline,
              const LINALG::Matrix<1,vpernode*nnode>&      funct,
              LINALG::TMatrix<T,3,1>&                      r) const
  {
    r.Clear();

    for (unsigned int dim=0; dim<3; ++dim)
      for (unsigned int i=0; i<vpernode*nnode; ++i)
        r(dim) += funct(i)*disp_totlag_centerline(3*i+dim);

     return;
  }

  //! compute derivative of beam centerline with respect to element parameter \xi \in [-1;1]
  template<unsigned int nnode, unsigned int vpernode, typename T>
  void Calc_r_xi(const LINALG::TMatrix<T,3*vpernode*nnode,1>& disp_totlag_centerline,
                 const LINALG::Matrix<1,vpernode*nnode>&      deriv,
                 LINALG::TMatrix<T,3,1>&                      r_xi) const
  {
    r_xi.Clear();

    for (unsigned int dim=0; dim<3; ++dim)
      for (unsigned int i=0; i<vpernode*nnode; ++i)
        r_xi(dim) += deriv(i)*disp_totlag_centerline(3*i+dim);

     return;
  }

  //! compute derivative of beam centerline with respect to arc-length parameter s
  template<unsigned int nnode, unsigned int vpernode, typename T>
  void Calc_r_s(const LINALG::TMatrix<T,3*vpernode*nnode,1>& disp_totlag_centerline,
                const LINALG::Matrix<1, vpernode*nnode>&     deriv,
                const double&                                jacobi,
                LINALG::TMatrix<T,3,1>&                      r_s) const
  {
    Calc_r_xi<nnode,vpernode,T>(disp_totlag_centerline,deriv,r_s);

    /* at this point we have computed derivative with respect to the element parameter \xi \in [-1;1];
     * as r' in (2.12) from (3.11) are derivatives with respect to the reference arc-length parameter s,
     * we have to divide it by the Jacobi determinant at the respective point*/
    r_s.Scale(1.0/jacobi);

     return;
  }

  //! compute interpolated local relative rotation \Psi^l according to (3.11), Jelenic 1999
  template<unsigned int nnode, typename T>
  void Calc_Psi_l(const std::vector<LINALG::TMatrix<T,3,1> >& Psi_li,
                  const LINALG::Matrix<1,nnode>&              func,
                  LINALG::TMatrix<T,3,1>&                     Psi_l) const
  {
    Psi_l.Clear();

    for (int dof=0; dof<3; ++dof)
      for (unsigned int node=0; node<nnode; ++node)
        Psi_l(dof)+= func(node)*(Psi_li[node])(dof);

    return;
  }

  //! compute interpolated local relative rotation \Psi^l according to (3.11), Jelenic 1999
  template<unsigned int nnode>
  void Calc_Psi_l(const std::vector<LINALG::Matrix<3,1> >& Psi_li,
                  const LINALG::Matrix<1,nnode>&           func,
                  LINALG::Matrix<3,1>&                     Psi_l) const
  {
    Psi_l.Clear();

    for (int dof=0; dof<3; ++dof)
      for (unsigned int node=0; node<nnode; ++node)
        Psi_l(dof)+= func(node)*(Psi_li[node])(dof);

    return;
  }

  //!compute derivative of interpolated local relative rotation \Psi^l with respect to reference arc-length parameter s according to (3.11), Jelenic 1999
  template<unsigned int nnode, typename T>
  void Calc_Psi_l_s(const std::vector<LINALG::TMatrix<T,3,1> >& Psi_li,
                  const LINALG::Matrix<1,nnode>&                deriv,
                  const double&                                 jacobi,
                  LINALG::TMatrix<T,3,1>&                       Psi_l_s) const
  {
    Psi_l_s.Clear();

    for (int dof=0; dof<3; ++dof)
      for (unsigned int node=0; node<nnode; ++node)
        Psi_l_s(dof)+= deriv(node)*(Psi_li[node])(dof);

    /* at this point we have computed derivative with respect to the element parameter \xi \in [-1;1];
     * as we want a derivatives with respect to the reference arc-length parameter s,
     * we have to divide it by the Jacobi determinant at the respective point*/
    Psi_l_s.Scale(1.0/jacobi);

    return;
  }

  //!compute local triad \Lambda from Crisfield 1999, eq. (4.7)
  template <typename T>
  void Calc_Lambda(const LINALG::TMatrix<T,3,1>& Psi_l,
                  const LINALG::TMatrix<T,4,1>&  Q_r,
                  LINALG::TMatrix<T,3,3>&        Lambda) const
  {
    Lambda.Clear();

    LINALG::TMatrix<T,4,1> Ql;
    LINALG::TMatrix<T,4,1> Qgauss;

    //c ompute relative rotation between triad at Gauss point and reference triad Qr
    LARGEROTATIONS::angletoquaternion(Psi_l,Ql);

    // compute rotation at Gauss point, i.e. the quaternion equivalent to \Lambda(s) in Crisfield 1999, eq. (4.7)
    LARGEROTATIONS::quaternionproduct(Ql,Q_r,Qgauss);

    // compute rotation matrix at Gauss point, i.e. \Lambda(s) in Crisfield 1999, eq. (4.7)
    LARGEROTATIONS::quaterniontotriad(Qgauss,Lambda);

    return;
  }

  //! compute quaternion equivalent to local triad \Lambda from Crisfield 1999, eq. (4.7)
  template <typename T>
  void Calc_Qgauss(const LINALG::TMatrix<T,3,1>& Psi_l,
                   const LINALG::TMatrix<T,4,1>& Q_r,
                   LINALG::TMatrix<T,4,1>&       Qgauss) const
  {
    Qgauss.Clear();

    LINALG::TMatrix<T,4,1> Ql;

    // compute relative rotation between triad at Gauss point and reference triad Qr
    LARGEROTATIONS::angletoquaternion(Psi_l,Ql);

    // compute rotation at Gauss point, i.e. the quaternion equivalent to \Lambda(s) in Crisfield 1999, eq. (4.7)
    LARGEROTATIONS::quaternionproduct(Ql,Q_r,Qgauss);

    return;
  }

  //! compute material curvature at certain Gauss point according to Crisfield 1999, eq. (4.9)
  template <typename T>
  void computeK(const LINALG::TMatrix<T,3,1>&      Psi_l,
                     const LINALG::TMatrix<T,3,1>& Psi_l_s,
                     const LINALG::Matrix<3,1>&    Kref,
                     LINALG::TMatrix<T,3,1>&       K) const
  {
    K.Clear();

    /* Calculation of material curvature vector according to Crisfield 1999, eq. (4.2) (this equation has been derived
     * for a different beam element formulation but is also valid for the element type considered here),
     * or Jelenic 1999, paragraph on page 153 between NOTE 5 and NOTE 6*/
    LINALG::TMatrix<T,3,3> Tinv(true);
    Tinv=LARGEROTATIONS::Tinvmatrix<T>(Psi_l);
    // It is important to use the transposed matrix Tinv^T instead of Tinv (these two only differ in one of three terms)
    K.MultiplyTN(Tinv,Psi_l_s);

    // mechanically relevant curvature is current curvature minus curvature in reference position
    for (int i=0; i<3; ++i)
      K(i) -= Kref(i);

  }

  //! compute convected strain at certain Gauss point with according to Crisfield 1999, eq. (3.4)
  template <typename T>
  void computeGamma(const LINALG::TMatrix<T,3,1>& r_s,
                    const LINALG::TMatrix<T,3,3>& Lambda,
                    const LINALG::Matrix<3,1>&    Gammaref,
                    LINALG::TMatrix<T,3,1>&       Gamma) const
  {
    Gamma.Clear();

    // convected strain gamma according to Crisfield 1999, eq. (3.4)
    Gamma.MultiplyTN(Lambda,r_s);

    /* In contrary to Crisfield 1999, eq. (3.4), the current implementation allows for initial values of the vector gammaref which
     * has also a second and a third component, i.e. it allows for initial shear deformation. This is the case, when the initial triad
     * at the evaluation point is not parallel to the centerline tangent vector at this point. The geometrically exact beam theory
     * does in general allow for such initial triads if they are considered consistently in the reference strains. While it is standard
     * to assume vanishing initial shear strains in the space-continuous setting, the possibility of initial shear strains might be advantageous
     * for the spatially discretized problem: For curved initial geometries, the nodal triad had to be determined such that the resulting
     * interpolated triad at the Gauss point would be tangential to the centerline tangent at this point resulting from the centerline interpolation.
     * In order to avoid this additional effort and to allow for an independent prescription of the nodal triads (e.g. prescribed by an analytical
     * geometry definition), the approach of considering arbitrary initial shear angles at the Gauss points is applied here.*/
    for (int i=0; i<3; ++i)
      Gamma(i) -= Gammaref(i);
  }

  //! setup constitutive matrices from material law
  template <typename T>
  void GetConstitutiveMatrices(LINALG::TMatrix<T,3,3>& CN,
                               LINALG::TMatrix<T,3,3>& CM) const;

  //! push forward material stress vector and constitutive matrix to their spatial counterparts by rotation matrix Lambda
  template <typename T>
  void pushforward(const LINALG::TMatrix<T,3,3>& Lambda,
                   const LINALG::TMatrix<T,3,1>& stress_mat,
                   const LINALG::TMatrix<T,3,3>& C_mat,
                   LINALG::TMatrix<T,3,1>& stress_spatial,
                   LINALG::TMatrix<T,3,3>& c_spatial) const;

  //! compute \tilde{I}^i in (3.18), page 152, Jelenic 1999, for all nodes i at a certain Gauss point
  template<unsigned int nnode>
  void computeItilde(const LINALG::TMatrix<double,3,1>&               Psil,
                     std::vector<LINALG::TMatrix<double,3,3> >&       Itilde,
                     const LINALG::TMatrix<double,3,1>&               phiIJ,
                     const LINALG::TMatrix<double,3,3>&               Lambdar,
                     const std::vector<LINALG::TMatrix<double,3,1> >& Psili,
                     const LINALG::Matrix<1,nnode>&                   funct)
  {
    /* note: Itilde is always be of type double because it is only needed for computation of analytical stiffness/mass matrix
     *       thus no need to template for FADordouble here*/

    // auxiliary matrices for storing intermediate results
    LINALG::TMatrix<double,3,3> auxmatrix;
    auxmatrix.PutScalar(0);
    LINALG::TMatrix<double,3,3> auxmatrix2;
    auxmatrix2.PutScalar(0);
    LINALG::TMatrix<double,3,3> Tinv_Psil = LARGEROTATIONS::Tinvmatrix<double>(Psil);

    // make sure that Itilde has proper dimensions
    Itilde.resize(nnode);

    // compute squared brackets term in (3.18), Jelenic 1999
    LINALG::TMatrix<double,3,3> squaredbrackets;
    squaredbrackets.PutScalar(0);
    for(unsigned int node=0; node<nnode; ++node)
    {
      auxmatrix.Clear();

      auxmatrix = LARGEROTATIONS::Tmatrix<double>(Psili[node]);
      auxmatrix.Scale(funct(node));
      auxmatrix2.Update(-1.0,auxmatrix,1.0);
    }

    squaredbrackets.Multiply(Tinv_Psil,auxmatrix2);

    for(int i=0; i<3; i++)
      squaredbrackets(i,i) += 1;

    // loop through all nodes i
    for (unsigned int node=0; node<nnode; ++node)
    {
      // compute rightmost term in curley brackets in (3.18), Jelenic 1999
      Itilde[node].Clear();
      Itilde[node].Multiply(Tinv_Psil,LARGEROTATIONS::Tmatrix<double>(Psili[node]));
      Itilde[node].Scale(funct(node));

      // if node i is node I then add squared brackets term times v_I
      if(node == nodeI_)
      {
        auxmatrix.Multiply(squaredbrackets,vI(phiIJ));
        Itilde[node] += auxmatrix;
      }

      // if node i is node J then add squared brackets term times v_J
      if(node == nodeJ_)
      {
        auxmatrix.Multiply(squaredbrackets,vJ(phiIJ));
        Itilde[node] += auxmatrix;
      }

      // now the term in the curly brackets has been computed and has to be rotated by \Lambda_r and \Lambda_r^t
      auxmatrix.MultiplyNT(Itilde[node],Lambdar);
      Itilde[node].MultiplyNN(Lambdar,auxmatrix);
    }

     return;
  }


  //! compute \tilde{I}^i if pre-calculated values are of type FAD
  template<unsigned int nnode>
  void computeItilde(const LINALG::TMatrix<FAD,3,1>&               Psil,
                     std::vector<LINALG::TMatrix<double,3,3> >&    Itilde,
                     const LINALG::TMatrix<FAD,3,1>&               phiIJ,
                     const LINALG::TMatrix<FAD,3,3>&               Lambdar,
                     const std::vector<LINALG::TMatrix<FAD,3,1> >& Psili,
                     const LINALG::Matrix<1,nnode>&                funct)
  {
    /* note: Itilde should always be of type double because it is only needed for computation of analytical stiffness/mass matrix
     * therefore no templated function used here; cast pre-calculated FAD values to double and call 'regular' function afterwards*/

    LINALG::TMatrix<double,3,1> Psil_double = FADUTILS::CastToDouble<FAD,3,1>(Psil);
    LINALG::TMatrix<double,3,1> phiIJ_double = FADUTILS::CastToDouble<FAD,3,1>(phiIJ);
    LINALG::TMatrix<double,3,3> Lambdar_double = FADUTILS::CastToDouble<FAD,3,3>(Lambdar);

    std::vector<LINALG::TMatrix<double,3,1> > Psili_double;
    for (unsigned int i=0; i<Psili.size(); ++i)
      Psili_double.push_back(FADUTILS::CastToDouble<FAD,3,1>(Psili[i]));

    computeItilde(Psil_double,Itilde,phiIJ_double,Lambdar_double,Psili_double,funct);

     return;
  }

  //! compute \tilde{I}^{i'} in (3.19), page 152, Jelenic 1999 for all nodes i at a certain Gauss point
  template<unsigned int nnode, typename T>
  void computeItildeprime(const LINALG::TMatrix<T,3,1>&               Psil,
                          const LINALG::TMatrix<T,3,1>&               Psilprime,
                          std::vector<LINALG::TMatrix<T,3,3> >&       Itildeprime,
                          const LINALG::TMatrix<T,3,1>&               phiIJ,
                          const LINALG::TMatrix<T,3,3>&               Lambdar,
                          const std::vector<LINALG::TMatrix<T,3,1> >& Psili,
                          const LINALG::Matrix<1,nnode>&              funct,
                          const LINALG::Matrix<1,nnode>&              deriv,
                          const double&                               jacobi)
  {
    // auxiliary matrices for storing intermediate results
    LINALG::TMatrix<T,3,3> auxmatrix;

    // make sure that Itildeprime has proper dimensions
    Itildeprime.resize(nnode);

    // matrix d(T^{-1})/dx
    LINALG::TMatrix<T,3,3> dTinvdx;
    LARGEROTATIONS::computedTinvdx<T>(Psil,Psilprime,dTinvdx);

    // compute T^{~} according to remark subsequent to (3.19), Jelenic 1999
    LINALG::TMatrix<T,3,3> Ttilde;
    Ttilde.PutScalar(0);
    for(unsigned int node=0; node<nnode; ++node)
    {
      auxmatrix = LARGEROTATIONS::Tmatrix<T>(Psili[node]);
      auxmatrix.Scale(funct(node));
      Ttilde += auxmatrix;
    }

    // compute T^{~'} according to remark subsequent to (3.19), Jelenic 1999
    LINALG::TMatrix<T,3,3> Ttildeprime;
    Ttildeprime.PutScalar(0);
    for(unsigned int node=0; node<nnode; ++node)
    {
      auxmatrix = LARGEROTATIONS::Tmatrix<T>(Psili[node]);
      auxmatrix.Scale(deriv(node)/jacobi);
      Ttildeprime += auxmatrix;
    }

    // compute first squared brackets term in (3.18), Jelenic 1999
    LINALG::TMatrix<T,3,3> squaredbrackets;
    squaredbrackets.PutScalar(0);
    squaredbrackets.Multiply(dTinvdx,Ttilde);
    auxmatrix.Multiply(LARGEROTATIONS::Tinvmatrix(Psil),Ttildeprime);
    squaredbrackets += auxmatrix;

    // loop through all nodes i
    for (unsigned int node=0; node<nnode; ++node)
    {

      // compute first term in second squared brackets
      Itildeprime[node] = dTinvdx;
      Itildeprime[node].Scale(funct(node));

      // compute second term in second squared brackets
      auxmatrix.PutScalar(0);
      auxmatrix += LARGEROTATIONS::Tinvmatrix(Psil);
      auxmatrix.Scale(deriv(node)/jacobi);

      // compute second squared brackets
      auxmatrix += Itildeprime[node];

      // compute second squared brackets time T(\Psi^l_j)
      Itildeprime[node].Multiply(auxmatrix,LARGEROTATIONS::Tmatrix<T>(Psili[node]));

      // if node i is node I then add first squared brackets term times v_I
      if(node == nodeI_)
      {
        auxmatrix.Multiply(squaredbrackets,vI(phiIJ));
        Itildeprime[node] -= auxmatrix;
      }

      // if node i is node J then add first squared brackets term times v_J
      if(node == nodeJ_)
      {
        auxmatrix.Multiply(squaredbrackets,vJ(phiIJ));
        Itildeprime[node] -= auxmatrix;
      }

      // now the term in the curly brackets has been computed and has to be rotated by \Lambda_r and \Lambda_r^t
      auxmatrix.MultiplyNT(Itildeprime[node],Lambdar);
      Itildeprime[node].MultiplyNN(Lambdar,auxmatrix);
    }

     return;
  }

  //! compute matrix v_I as outlined in the equations above (3.15) on page 152 of Jelenic 1999
  LINALG::TMatrix<double,3,3> vI(const LINALG::TMatrix<double,3,1>& phiIJ)
  {
    // matrix v_I
    LINALG::TMatrix<double,3,3> result;

    LARGEROTATIONS::computespin(result,phiIJ);
    if(phiIJ.Norm2() == 0)
      result.Scale(0.25);
    else
      result.Scale(tan(phiIJ.Norm2()/4.0)/phiIJ.Norm2());

    for(int i=0; i<3; i++)
      result(i,i) +=1;

    result.Scale(0.5);

    return result;
  }

  //! compute matrix v_J as outlined in the equations above (3.15) on page 152 of Jelenic 1999
  LINALG::TMatrix<double,3,3> vJ(const LINALG::TMatrix<double,3,1>& phiIJ)
  {
    // matrix v_J
    LINALG::TMatrix<double,3,3> result;

    LARGEROTATIONS::computespin(result,phiIJ);
    if(phiIJ.Norm2() == 0)
      result.Scale(-0.25);
    else
      result.Scale(-tan(phiIJ.Norm2()/4.0)/phiIJ.Norm2());

    for(int i=0; i<3; i++)
      result(i,i) +=1;

    result.Scale(0.5);

    return result;
  }

  //! calculation of nonlinear stiffness and mass matrix
  template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
  void CalcInternalAndInertiaForcesAndStiff(Teuchos::ParameterList&   params,
                                                 std::vector<double>&      acc,
                                                 std::vector<double>&      vel,
                                                 std::vector<double>&      disp,
                                                 Epetra_SerialDenseMatrix* stiffmatrix,
                                                 Epetra_SerialDenseMatrix* massmatrix,
                                                 Epetra_SerialDenseVector* force,
                                                 Epetra_SerialDenseVector* inertia_force);

  //! calculation of thermal (i.e. stochastic) and damping forces according to Brownian dynamics
  template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
  void CalcBrownianForcesAndStiff(Teuchos::ParameterList&   params,
                                  std::vector<double>&      vel,
                                  std::vector<double>&      disp,
                                  Epetra_SerialDenseMatrix* stiffmatrix,
                                  Epetra_SerialDenseVector* force);

  //! update (total) displacement vector and set nodal triads (as quaternions)
  template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
  void UpdateDispTotLagAndNodalTriads(const std::vector<double>&                         disp,
                                      LINALG::TMatrix<FADordouble,3*vpernode*nnodecl,1>& disp_totlag_centerline,
                                      std::vector<LINALG::TMatrix<FADordouble,4,1> >&    Q_i);

  //! lump mass matrix
  template<unsigned int nnode>
  void lumpmass(Epetra_SerialDenseMatrix* emass);

  //! determine Gauss rule from required type of integration and parameter list
  DRT::UTILS::GaussRule1D MyGaussRule(const Teuchos::ParameterList&  params,
                                      const IntegrationPurpose       intpurpose) const;

  //! determine Gauss rule from required type of integration
  inline DRT::UTILS::GaussRule1D MyGaussRule(const IntegrationPurpose intpurpose) const
  {
    Teuchos::ParameterList p;
    return MyGaussRule(p,intpurpose);
  };

  //! calculates the element length in reference configuration
  template<unsigned int nnode, unsigned int vpernode>
  void Calculate_reflength(const LINALG::Matrix<3*vpernode*nnode,1>& disp_totlag_centerline,
                           const double                              tolerance);

  //@}

  //! @name Methods for Brownian dynamics simulations

  //! compute damping matrix and statistical forces according to fluctuation dissipation theorem
  template<unsigned int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
  void CalcBrownian(Teuchos::ParameterList&                          params,  //!<parameter list
                   const std::vector<double>&                       vel,  //!< element velocity vector
                   const std::vector<double>&                       disp, //!<element disp vector
                   Epetra_SerialDenseMatrix*                        stiffmatrix,  //!< element stiffness matrix
                   Epetra_SerialDenseVector*                        force,
                   std::vector<LINALG::Matrix<1,nnode> >&           Imass,
                   std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildemass);//!< element internal force vector

  //! computes damping coefficients for employed friction model
  void MyDampingConstants(Teuchos::ParameterList& params,
                         LINALG::Matrix<3,1>&    gamma);

  //! computes rotational damping forces and stiffness
  template<unsigned int nnode> //number of nodes
  void MyRotationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                         const std::vector<double>&                             vel,  //!< element velocity vector
                         const std::vector<double>&                             disp, //!<element disp vector
                         Epetra_SerialDenseMatrix*                              stiffmatrix,//!< element stiffness matrix
                         Epetra_SerialDenseVector*                              force,      //!< element internal force vector
                         const DRT::UTILS::IntegrationPoints1D&                 gausspointsdamping,
                         const std::vector<LINALG::Matrix<1,nnode> >&           Idamping,
                         const std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildedamping,
                         const std::vector<LINALG::Matrix<4,1> >&               Qconvdamping,
                         const std::vector<LINALG::Matrix<4,1> >&               Qnewdamping);

  //! computes translational damping forces and stiffness
  template<unsigned int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
  void MyTranslationalDamping(Teuchos::ParameterList&  params,  //!<parameter list
                           const std::vector<double>& vel,  //!< element velocity vector
                           const std::vector<double>& disp, //!<element disp vector
                           Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                           Epetra_SerialDenseVector*  force);//!< element internal force vector

  //! computes stochastic translational forces and resulting stiffness
  template<unsigned int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
  void MyStochasticForces(Teuchos::ParameterList&  params,  //!<parameter list
                       const std::vector<double>& vel,  //!< element velocity vector
                       const std::vector<double>& disp, //!<element disp vector
                       Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                       Epetra_SerialDenseVector*  force);//!< element internal force vector

  //! computes stochastic rotational moments and resulting stiffness
  template<unsigned int nnode, int randompergauss> //number of nodes, number of random numbers required per Gauss point
  void MyStochasticMoments(Teuchos::ParameterList&                             params,  //!<parameter list
                       const std::vector<double>&                             vel,  //!< element velocity vector
                       const std::vector<double>&                             disp, //!<element disp vector
                       Epetra_SerialDenseMatrix*                              stiffmatrix,//!< element stiffness matrix
                       Epetra_SerialDenseVector*                              force,      //!< element internal force vector
                       const DRT::UTILS::IntegrationPoints1D&                 gausspointsdamping,
                       const std::vector<LINALG::Matrix<1,nnode> >&           Idamping,
                       const std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildedamping,
                       const std::vector<LINALG::Matrix<4,1> >&               Qconvdamping,
                       const std::vector<LINALG::Matrix<4,1> >&               Qnewdamping);

  //! computes velocity of background fluid and gradient of that velocity at a certain evaluation point in the physical space and adds respective terms to internal forces and damping matrix
  template<int ndim> //number of dimensions of embedding space
  void MyBackgroundVelocity(Teuchos::ParameterList&       params,  //!<parameter list
                           const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                           LINALG::Matrix<ndim,1>&       velbackground,  //!< velocity of background fluid
                           LINALG::Matrix<ndim,ndim>&    velbackgroundgrad); //!<gradient of velocity of background fluid

  //! shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions
  template<unsigned int nnode, int ndim> //number of nodes, number of dimensions
  void NodeShift(Teuchos::ParameterList& params,  //!<parameter list
                std::vector<double>&    disp); //!<element disp vector


  //@}

  // don't want = operator
  Beam3r& operator = (const Beam3r& old);

}; // class Beam3r


// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);


} // namespace Elements
} // namespace DRT


#endif  // #ifndef BEAM3R_H
