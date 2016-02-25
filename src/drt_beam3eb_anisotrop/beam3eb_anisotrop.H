/*!----------------------------------------------------------------------
\file beam3eb_anisotrop.H

\brief three dimensional nonlinear rod based on a C1 curve

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15301
</pre>

3D nonlinear Kirchhoff-like beam element. In contrary to the beam3ebtor element this element can
also display initially curved beams. To be able to use this element correctly structural dynamic
parameters need to be set to:

  LOADLIN       Yes
  MASSLIN      Yes (only necessary for dynamics)

since due to this beam formulation external point loads are being linearized and have an effect on the stiffness matrix. For this reason a special
type of Neumann conditions, namely DESIGN POINT MOMENT EB CONDITIONS, are needed for this element.

As the beam curve has to be C1 it is interpolated with hermitian polynomials of order 3. Therefore each node has 7 dofs: the position
vector of the node (3 dofs), the tangent vector to the curve at the node (3 dofs) and a relative angle (1 dof) describing the rotation between
the material and a reference triad.
*----------------------------------------------------------------------*/


//header file only included if not yet included!
#ifndef BEAM3EBANISOTROP_H
#define BEAM3EBANISOTROP_H

#include "Teuchos_RCP.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_node.H"

#include "Sacado.hpp"
typedef Sacado::Fad::DFad<double> FAD;

//The default values of the following define flags have been used to create the corresponding test case beam3ebansiotrop_helix.dat

//ONE of the following three reference systems have to be chosen. DEFAULT: NSRISR
#define NSRISR  //Nodal Smallest Rotation Interpolation via Smallest Rotation
//#define VP //Vector product Reference System
//#define SR1 //Smallest Rotation Reference System
//DEFAULT: NSRISR

#define ANS_BEAM3EB_ANISOTROP 1 //Apply Assumed Natural Strain approach. 1=quadratic interpolation, 2=constant average value.
//DEFAULT: ON
//#define CONSISTENTANS //Decide wether the variational correct or the simplified ANS approach should be applied
//DEFAULT: Off

//#define ANSKAPPA 2 //Apply Assumed Natural Strain approach for curvature kappa. DEFAULT: OFF

//#define MATERIALREF  //decide whether the old material or the old reference triad is used to define the new reference triad.
//DEFAULT: Off

#define TWISTDOFS 4 //degrees of freedom per element related to twist angle; only the values 2,3,4 are valid.
// DEFAULT: 4

//#define NODALALPHAT 2 //Needed for dynamics --> 1: alpha_t = 0; 2: old reference triad constant; 3: old reference triad non-constant.
//DEFAULT: OFF

#define mygaussruleebanisotrop intrule_line_4point //define gauss rule; intrule_line_1point - intrule_line_10point is implemented.
//DEFAULT: intrule_line_3point

//#define BEAM3EBANISOTROPORTHOPRESSURE 3.1415926535897938   //Apply orthopressure load. DEFAULT: off

//#define OPTCURVE
//#define OPTCURVEAX
//#define NUMELE 2


#define ROTATIONTOL 1.0e-12 //Tolerance that ensures a rotation per load/time step, that is smaller than 180°: 1+t*t_old !< ROTATIONTOL. DEFAULT: 1.0e-12

#define LENGTHCALCNEWTONTOL 10e-12 //Tolerance for iterative determination of element length. DEFAULT: 10e-12

#if !(defined NSRISR) and !(defined VP) and !(defined SR1)
  #error "Beam3ebanisotrop: No local triad defined!!!"
#endif

#if ((defined NSRISR) and ((defined VP) or (defined SR1))) or ((defined VP) and ((defined SR1) or (defined NSRISR)))
  #error "Only chose one reference triad, SR, VP or NSRISR!!!"
#endif

//! \brief DRT: namespace of the ccarat discretization module
namespace DRT
{

namespace ELEMENTS
{

class Beam3ebanisotropType : public DRT::ElementType
{
public:

  std::string Name() const { return "Beam3ebanisotropType"; }

  static Beam3ebanisotropType & Instance();

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

  static Beam3ebanisotropType instance_;
};

/*!
\brief 3D nonlinear Kirchhoff-like beam element. In contrary to the beam3ebtor element this element can
also display initially curved beams.

author meier (meier@lnm.mw.tum.de)
*/
class Beam3ebanisotrop : public DRT::Element
{
public:

  //! @name Friends
  friend class Beam3ebanisotropType;


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3ebanisotrop(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3ebanisotrop(const Beam3ebanisotrop& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3ebanisotrop();

  /*!
  \brief Deep copy this instance of Beam3ebanisotrop and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return Beam3ebanisotropType::Instance().UniqueParObjectId(); }

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
  { return Beam3ebanisotropType::Instance(); }

  //@}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 1;}


  /*!
  \brief Get vector of Teuchos::RCPs to the lines of this element
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

    int dofpn_aux=0.0;

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
  \brief Get number of detriadgrees of freedom per element not including nodal degrees of freedom
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

  /*!
  \brief get nodal tangents in reference configuration
  */
  std::vector<LINALG::Matrix<3,1> > Tref() const
  {
    return Tref_;
  }


  //@}


private:

  //! action parameters recognized by Beam3ebanisotrop
  enum ActionType
  {
    calc_none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_reset_istep,
    calc_struct_ptcstiff, //calculate artificial stiffness due to PTC solution strategy
    calc_struct_energy, // calculate strain energy
    calc_struct_update_istep,    // two actions necessary to update element data at the end of the
    calc_struct_update_imrlike,  // step. The difference between the two methods is only relevant for EAS technology.
  };

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

  //! Matrices holding the basis vectors of the previous time-step at the gauss points
  std::vector<LINALG::Matrix<3,1> > t0_;
  std::vector<LINALG::Matrix<3,1> > n0_;
  std::vector<LINALG::Matrix<3,1> > b0_;
  //! Matrices holding the derivative of the basis vectors of the previous time-step at the gauss points
  std::vector<LINALG::Matrix<3,1> > dt0ds_;
  std::vector<LINALG::Matrix<3,1> > db0ds_;
  //! Matrices holding the basis vectors of the previous time-step at the element nodes
  std::vector<LINALG::Matrix<3,1> > t0_nodes_;
  std::vector<LINALG::Matrix<3,1> > n0_nodes_;
  std::vector<LINALG::Matrix<3,1> > b0_nodes_;
  //! Matrices holding the derivative of the basis vectors of the previous time-step at the element nodes
  std::vector<LINALG::Matrix<3,1> > dt0ds_nodes_;
  std::vector<LINALG::Matrix<3,1> > db0ds_nodes_;
  //! Matrix holding the angular velocity of the reference system of the previous time-step at the element nodes
  std::vector<LINALG::Matrix<3,1> > w0_nodes_;
  //! d(w0_nodes-(w0_nodes*g1)g1)/dt
  std::vector<LINALG::Matrix<3,1> > dw0perpdt_nodes_;
  //! d((w0_nodes*g1))/dt
  std::vector<double> dw0paralleldt_nodes_;
  //! Matrix holding the angular velocity of the reference system of the current time-step at the element nodes
  std::vector<LINALG::Matrix<3,1> > w_nodes_;
  //! d(w0_nodes-(w0_nodes*g1)g1)/dt
  std::vector<LINALG::Matrix<3,1> > dwperpdt_nodes_;
  //! d((w0_nodes*g1))/dt
  std::vector<double> dwparalleldt_nodes_;
  std::vector<double> alphat_nodes_;
  std::vector<double> alphatt_nodes_;

  //Variables needed by all local triads
  //! Matrix holding tangents at each node in reference configuration
  std::vector<LINALG::Matrix<3,1> > Tref_;
  //! Matrix holding material normal vectors at each node in reference configuration
    std::vector<LINALG::Matrix<3,1> > G2ref_;
  //! Matrix with the curvature in the initial configuration at each gp
  std::vector<double> kappa0_;
  //! Matrizes with the projection of the initial curvature kappa0 onto the initial material base vectors g20 and g30
  std::vector<double> kappa0g20_;
  std::vector<double> kappa0g30_;
  //! vector that holds the relative angle between nodal FS and SR system at the old time step
  std::vector<double> theta_nodes_old_;
  //! relative angle between the two SR systems at the second node at the old time step
  double sr_theta_nodes_old_;
  //! Matrix with the torsion in the initial configuration
  std::vector<double> tau0_;
  //! Length of the element
  double length_;
  //! optimal Hermite constant, calculated via least square optimization
  double copt_;
  //! Jacobi determinant of the element ds/dxi=abs(dr/dxi)=|r'0(xi)|
  std::vector<double> jacobi_;
  //! Variable holding a factor needed for the second derivative r'^T(xi)r''(xi)
  std::vector<double> jacobi2_;
  //! Variable holding a factor needed for the third derivative r''^T(xi)r''(xi)+r'^T(x)r'''(xi)
  std::vector<double> jacobi3_;
  //! variable holding the total internal energy
  double int_energy_;
  //! bool recognizing first element call
  bool firstcall_;

  int timestepcount_;

  //Functions

  //! @name methods for initilization of the element

  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration)
  void SetUpReferenceGeometry(const std::vector<LINALG::Matrix<3,1> >& xrefe,const bool secondinit = false);

  //@}

  //! @name auxiliary computation methods for non-additive, large rotation variables

  //@}

  //! @name Internal calculation methods

  //! Calculate stiffness matrix via FAD
  void CalculateInternalForces(Teuchos::ParameterList& params,
                               std::vector<double>& vel,
                               std::vector<double>& disp,
                               Epetra_SerialDenseMatrix* stiffmatrix,
                               Epetra_SerialDenseVector* force,
                               bool update);

  //! Calculate inertia forces
  void CalculateInertialForces(Teuchos::ParameterList& params,
                               std::vector<double>& acc,
                               std::vector<double>& vel,
                               std::vector<double>& disp,
                               Epetra_SerialDenseMatrix* massmatrix,
                               Epetra_SerialDenseVector* force_inertia,
                               bool update);

  //! Calculate elastic energy of the element
   void CalculateEnergy(Epetra_SerialDenseVector& energy);

   //! Update reference triads of the element

   void UpdateTriads(Teuchos::ParameterList& params,
                     std::vector<double>& vel,
                     std::vector<double>& disp);

   //! Shift definition interval of an given angle theta to prescribed interval [theta_min;theta_min + 360°] with theta_min \in ]-360°;0°]
   FAD ShiftAngleIntervall(FAD theta, double theta_min, bool min);

   //! Integrate the torsion over the element length in order to get the twist angle
   void IntegrateTorsion(std::vector<FAD> disp,
                         std::vector<FAD>& delta_alpha,
                         FAD& alpha_ges,
                         std::vector<FAD>& delta_s,
                         int sumgp,
                         LINALG::TMatrix<FAD,3,6>& nodal_triads,
                         LINALG::TMatrix<FAD,3,2>& nodal_curvature,
                         std::vector<FAD>& gamma_s_nodes,
                         std::vector<FAD>& jacobi_nodes,
                         std::vector<FAD>& theta_nodes,
                         bool refcalc,
                         std::vector<LINALG::TMatrix<FAD,2,3> >& sr_triads_nodes,
                         std::vector<LINALG::TMatrix<FAD,2,3> >& sr_triads_gp,
                         std::vector<FAD>& sr_torsion_gp,
                         FAD& sr_theta_nodes);

   //! Calculate the nodal reference triads as well as the difference angle theta_nodes at the element nodes
   void DetermineNodalTriads(std::vector<FAD> disp_totlag,
                             LINALG::TMatrix<FAD,3,6>& triads_mat_nodes,
                             std::vector<FAD>& theta_nodes,
                             bool refcalc,
                             std::vector<LINALG::TMatrix<FAD,2,3> >& triads_ref_nodes,
                             FAD& taubar0);

   //! Calculate the angular velocity of the reference system at the element nodes
   void DetermineWrefNodes(std::vector<FAD> disp_totlag,
                                std::vector<FAD> vel_totlag,
                                std::vector<FAD> acc_totlag,
                                LINALG::TMatrix<FAD,2,3>& w_ref_nodes,
                                std::vector<LINALG::TMatrix<FAD,3,1> >& w_ref_perp_t_nodes,
                                std::vector<FAD>& w_ref_parallel_t_nodes,
                                FAD& theta_t,
                                FAD& theta_tt);

   //! Calculate material triad and the mechanical curvature at the gauss points
   void CalculateMaterialTriad(double xi,
                               int numgp,
                               std::vector<FAD> theta_nodes,
                               LINALG::TMatrix<FAD,2,3> triad_bar_gp,
                               FAD tau_bar,
                               FAD gamma_gp,
                               LINALG::TMatrix<FAD,2,3>& triad_mat_gp,
                               FAD& tau_gp);

  //! Calculate the SR mapping for a given reference system triad_ref and a given tangent vector r_s
  void CalculateSRTriads( LINALG::TMatrix<FAD,3,1> r_s,
                          LINALG::TMatrix<FAD,3,3> triad_ref,
                          LINALG::TMatrix<FAD,3,3>& triad);

  //! Calculate the VP mapping for a given reference system triad_ref and a given tangent vector r_s
  void CalculateVPTriads( LINALG::TMatrix<FAD,3,1> r_s,
                          LINALG::TMatrix<FAD,3,3> triad_ref,
                          LINALG::TMatrix<FAD,3,3>& triad);

  //! Calculate the tangential component WISR of an angular velocity of a reference system created with a reference triad (tangent_ref with w_ref) onto a given tangent with w_perp
  FAD CalculateWISR(LINALG::TMatrix<FAD,3,1> tangent,
                    LINALG::TMatrix<FAD,3,1> w_perp,
                    LINALG::TMatrix<FAD,3,1> tangent_ref,
                    LINALG::TMatrix<FAD,3,1> w_ref,
                    bool nodal);

  //! Calculate the time derivative of the angular velocity d(WISR)/dt
  FAD CalculateWdotISR( LINALG::TMatrix<FAD,3,1> tangent,
                        LINALG::TMatrix<FAD,3,1> tangent_t,
                        LINALG::TMatrix<FAD,3,1>  w_perp,
                        LINALG::TMatrix<FAD,3,1>  w_perp_t,
                        LINALG::TMatrix<FAD,3,1> tangent_ref,
                        LINALG::TMatrix<FAD,3,1>  tangent_ref_t,
                        LINALG::TMatrix<FAD,3,1> w_ref,
                        LINALG::TMatrix<FAD,3,1> w_ref_perp_t,
                        FAD w_ref_parallel_t,
                        bool nodal);

  //! Assembles the four shape functions of an element into a 3*dofs_per_element block matrix.
  //The same is done for the twist shape functions and the corresponding derivatives.
  void AssembleShapefunctions(LINALG::TMatrix<FAD,1,4> N_i,
                              LINALG::TMatrix<FAD,1,4> N_i_xi,
                              LINALG::TMatrix<FAD,1,4> N_i_xixi,
                              LINALG::TMatrix<FAD,1,4> N_i_xixixi,
                              LINALG::TMatrix<FAD,1,TWISTDOFS> C_i,
                              LINALG::TMatrix<FAD,1,TWISTDOFS> C_i_xi,
                              FAD jacobi_local,
                              FAD jacobi2_local,
                              FAD jacobi3_local,
                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_s,
                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_ss,
                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_sss,
                              LINALG::TMatrix<FAD,1,2*6+TWISTDOFS>& C,
                              LINALG::TMatrix<FAD,1,2*6+TWISTDOFS>& C_s);

  void AssembleShapefunctions( LINALG::Matrix<1,4> N_i_xi,
                               LINALG::Matrix<1,4> N_i_xixi,
                               double jacobi_local,
                               double jacobi2_local,
                               LINALG::Matrix<3,2*6+TWISTDOFS>& N_s,
                               LINALG::Matrix<3,2*6+TWISTDOFS>& N_ss);

  //! Assembles the four first shape function derivatives of an element into a 3*dofs_per_element block matrix.
  void AssembleShapefunctions( LINALG::TMatrix<FAD,1,4> N_i,
                               LINALG::TMatrix<FAD,1,4> N_i_xi,
                               LINALG::TMatrix<FAD,1,TWISTDOFS> C_i,
                               FAD jacobi_local,
                               LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N,
                               LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_s,
                               LINALG::TMatrix<FAD,1,2*6+TWISTDOFS>& C);

  void AssembleShapefunctions(LINALG::TMatrix<FAD,1,4> N_i_xi,
                              FAD jacobi_local,
                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N_s);

  void AssembleShapefunctions(LINALG::TMatrix<FAD,1,4> N_i,
                              LINALG::TMatrix<FAD,3,2*6+TWISTDOFS>& N);

  //! Calculate the intermediate triad triad_bar_gp and the corresponding torsion tau_bar_gp at the gauss points
  void CalculateIntermediateTriad(LINALG::TMatrix<FAD,3,1> r_s,
                                  LINALG::TMatrix<FAD,3,1> r_ss,
                                  const int numgp,
                                  std::vector<LINALG::TMatrix<FAD,2,3> > triads_ref_nodes,
                                  LINALG::TMatrix<FAD,2,3>& triad_bar_gp,
                                  FAD& tau_bar_gp);

  //! Calculates the Norm of a FAD vector, since .Norm2() is not available for FAD vectors
  FAD Norm(LINALG::TMatrix<FAD,3,1> v);

  //! Completes two orthonormal base-vectors to a right-handed orthonormal basis
  LINALG::TMatrix<FAD,3,1> CompleteTriad(LINALG::TMatrix<FAD,3,1> first_basevector, LINALG::TMatrix<FAD,3,1> second_basevector);

  //! Scales a vector with a factor
  LINALG::TMatrix<FAD,3,1> ScaleVector(LINALG::TMatrix<FAD,3,1> vector, FAD factor);

  //! Computes inner product of two FAD vectors
  FAD ScalarProduct(LINALG::TMatrix<FAD,3,1> first_vector, LINALG::TMatrix<FAD,3,1> second_vector);

  //! Computes inner product of two double vectors
  double ScalarProduct(LINALG::Matrix<3,1> first_vector, LINALG::Matrix<3,1> second_vector);

  //! Computes outer product of two FAD vectors
  LINALG::TMatrix<FAD,3,1> VectorProduct(LINALG::TMatrix<FAD,3,1> first_vector, LINALG::TMatrix<FAD,3,1> second_vector);

  //! Calculates nodal positions, tangents and relative twist angles out of the corresponding displacements
  void UpdateDispTotlag(std::vector<double> disp, std::vector<FAD>& disp_totlag);

  //! Calculates nodal positions, tangents and relative twist angles as well as velocities and accelerations out of the corresponding displacements, velocities and accelerations
  void UpdateDispVelAccTotlag(std::vector<double> disp,
                              std::vector<double> vel,
                              std::vector<double> acc,
                              std::vector<FAD>& disp_totlag,
                              std::vector<FAD>& vel_totlag,
                              std::vector<FAD>& acc_totlag);

  //!Calculates tension at collocation points (CPs)
  void CalculateEpsilonCP( std::vector<FAD> disp_totlag,
                           LINALG::TMatrix<FAD,3,1> r_s_midpoint,
                           LINALG::TMatrix<FAD,3,1>& epsilonCP,
                           LINALG::TMatrix<FAD,12+TWISTDOFS,3>& LINepsilonCP_T);


  //!Calculates derivative r' at certain point xi
  void Calculate_r_s( std::vector<FAD> disp_totlag,
                      double xi,
                      LINALG::TMatrix<FAD,3,1>& r_s_midpoint);

  //!Calculates derivatives r' and r'' at certain point xi
  void Calculate_r_s_and_r_ss( std::vector<FAD> disp_totlag,
                             double xi,
                             LINALG::TMatrix<FAD,3,1>& r_s,
                             LINALG::TMatrix<FAD,3,1>& r_ss);

  //!Store the gauss point reference triads of the converged configuration in the corresponding class variables t0_, n0_, b0_, dt0ds_ and db0ds_
  void UpdateGPRefTriads( LINALG::TMatrix<FAD,3,1> tangent,
                             LINALG::TMatrix<FAD,2,3> triad_mat_gp,
                             LINALG::TMatrix<FAD,2,3> triad_bar_gp,
                             LINALG::TMatrix<FAD,3,1> kappa_vec,
                             FAD tau_gp,
                             FAD gamma_s,
                             int numgp);
  //!Store the nodal reference triads of the converged configuration in the corresponding class variables t0_nodes_, n0_nodes_, b0_nodes_
  void UpdateNodalRefTriads(LINALG::TMatrix<FAD,3,1> tangent,
                            LINALG::TMatrix<FAD,3,6> triads_mat_nodes,
                            std::vector<LINALG::TMatrix<FAD,2,3> > triads_ref_nodes,
                            std::vector<FAD> theta_nodes);

  //!Update of the old reference angular velocities and angular accelerations at the nodes
  void UpdateNodalAngularVelocities(Epetra_SerialDenseVector* update_disp,
                                    Epetra_SerialDenseVector* update_vel,
                                    Epetra_SerialDenseVector* update_acc);

  //! lump mass matrix
  void lumpmass(Epetra_SerialDenseMatrix* emass);

  //TODO: Modification
  //! Computes the length of the element and its Jacobi based upon a hermitian interpolation of the geometry
  void calculate_length(const std::vector<LINALG::Matrix<3,1> >& xrefe, double tolerance);
  //! Computes the initial local triads (=Frenet Serret triads) at s specified, returns a vector of the base vectors: t - n - b
  std::vector<LINALG::TMatrix<FAD,3,1> > calculate_fs_triad(LINALG::TMatrix<FAD,3,1>& r_s, LINALG::TMatrix<FAD,3,1>& r_ss);
  //! Computes the initial curvature at a specified point
  LINALG::Matrix<3,1> calculate_curvature(const LINALG::Matrix<3,1>& r_s, const LINALG::Matrix<3,1>& r_ss);

  LINALG::TMatrix<FAD,3,1> calculate_curvature(LINALG::TMatrix<FAD,3,1>& r_s, LINALG::TMatrix<FAD,3,1>& r_ss);

  //! Computes the initial torsion at a specified point
  double calculate_fstorsion(const LINALG::Matrix<3,1>& r_s, const LINALG::Matrix<3,1>& r_ss, const LINALG::Matrix<3,1>& r_sss, bool straight);
  //! Computes the jacobi at a specified point
  std::vector<double> calculate_jacobi(const LINALG::Matrix<3,1>& r_xi, const LINALG::Matrix<3,1>& r_xixi, const LINALG::Matrix<3,1>& r_xixixi);
  //! Calculates Signum function of a given "value"
  double Signum(FAD value);

  // don't want = operator
  Beam3ebanisotrop& operator = (const Beam3ebanisotrop& old);

}; // class Beam3ebanisotrop

// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);

} // namespace Elements
} // namespace DRT

#endif  // #ifndef BEAM3EBANISOTROP_H
