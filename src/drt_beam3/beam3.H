/*!----------------------------------------------------------------------
\file beam3.H

\brief three dimensional nonlinear Reissner beam element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
// line in order to swith on and off by comment in defines.drt
#ifdef D_BEAM3

// administrative line in order to make a difference between BACI and Ccarat files
#ifdef CCADISCRET

//header file only included if not yet included!
#ifndef BEAM3_H
#define BEAM3_H


#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_statmech/statmech_manager.H"
#include "../drt_statmech/statmech_time.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_statmech.H"



//! \brief DRT: namespace of the ccarat discretization module
namespace DRT
{

namespace ELEMENTS
{
/*!
\brief three dimensional nonlinear corotational Timoshenko beam element according to Crisfield, Vol.2,
section 17.2 ("an interpretation of an element due to simo and vu-quoc"), p. 226 - 231; the current
implementation sets the initial director angles and curvature to zero, which can be done without loss
of generality since any initial rotation cancels out in the equations (17.72) and (17.76). As a consequence
the current nodal directors do not appear in the element, but only the values of the current director
angles minus reference director angles, respectively

The element is capable of using higher order functions from linear to quartic. Please make sure you put the nodes in the right order
in the input file.
LIN2  1---2
LIN3	 1---3---2
LIN4  1---4---2---3
LIN5	 1---5---2---3---4

author cyron (cyron @lnm.mw.tum.de)
*/
class Beam3 : public DRT::Element
{
public:

  //! @name Friends
  friend class Beam3Register;
  friend class StatMechManager;
  friend void StatMechManager::SetCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,const std::map<int,LINALG::Matrix<3,1> >& currentpositions,const std::map<int,LINALG::Matrix<3,1> >& currentrotations);
  friend void StatMechTime::Integrate();
  friend StatMechTime::StatMechTime(ParameterList& params,DRT::Discretization& dis,LINALG::Solver& solver,IO::DiscretizationWriter& output);
  
  

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3(const Beam3& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3();

  /*!
  \brief Deep copy this instance of Beam3 and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return ParObject_Beam3; }

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
  \brief Allocate and return an ElementRegister class

  Every implementation of DRT::Element must also implement
  an DRT::ElementRegister. This method allocates and returns
  an instance of a matching DRT::ElementRegister.
  The ElementRegister class is then used to hold track of what type of
  elements are in one discretization.
  */
  virtual RefCountPtr<DRT::ElementRegister> ElementRegister() const;

  //@}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 1;}


  /*!
  \brief Get vector of RefCountPtrs to the lines of this element
  */
  virtual vector<RCP<DRT::Element> > Lines();


  /*!
  \brief Get number of degrees of freedom of a single node
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    /*note: this is not necessarily the number of DOF assigned to this node by the discretization finally, but only
    *the number of DOF requested for this node by this element; the discretization will finally assign the maximal
    *number of DOF to this node requested by any element connected to this node*/
    return 6;
  }

  /*!
  \brief Get number of degrees of freedom per element not including nodal degrees of freedom
  */
  virtual int NumDofPerElement() const { return 0;}

  /*!
  \brief Print this element
  */
  virtual void Print(ostream& os) const;


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
  virtual int Evaluate(ParameterList&            params,
                       DRT::Discretization&      discretization,
                       vector<int>&              lm,
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
  virtual int EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
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
  template<int nnode>
  void EvaluatePTC(ParameterList& params,
                          Epetra_SerialDenseMatrix& elemat1);

  //@}


private:

  //! action parameters recognized by Beam3
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
    calc_struct_update_istep,
    calc_struct_update_imrlike,
    calc_struct_reset_istep, //!< reset elementwise internal variables during iteration to last converged state
    calc_struct_ptcstiff //calculate artificial stiffness due to PTC solution strategy
  };
  
  //! kind of integration to be performed
  enum IntegrationType
  {
    gaussexactintegration,
    gaussunderintegration,
    lobattointegration
  };

  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;

  //! quaternions of triads at GPs at the end of the preceeding time step due to Crisfield, Vol. 2, equation (17.65)
  vector<LINALG::Matrix<4,1> > Qconv_;
  //! quaternions of triads at GPs at the end of the preceeding iteration step due to Crisfield, Vol. 2, equation (17.65)
  vector<LINALG::Matrix<4,1> > Qold_;
  //! quaternions of triads at GPs during the current iteration step due to Crisfield, Vol. 2, equation (17.65)
  vector<LINALG::Matrix<4,1> > Qnew_;
  //!curvaturea at the end of the preceeding time step due to Crisfield, Vol. 2, equation (17.72)
  vector<LINALG::Matrix<3,1> > curvconv_;
  //!curvaturea at the end of the preceeding iteration step due to Crisfield, Vol. 2, equation (17.72)
  vector<LINALG::Matrix<3,1> > curvold_;
  //!curvaturea during the current iteration step due to Crisfield, Vol. 2, equation (17.72)
  vector<LINALG::Matrix<3,1> > curvnew_;

  //!angle displacement taken from discretization and interpolated for a specific Gauss point by Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding time step
  vector<LINALG::Matrix<3,1> > thetaconv_;
  //!angle displacement taken from discretization and interpolated for a specific Gauss point by Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding iteration step
  vector<LINALG::Matrix<3,1> > thetaold_;
  //!angle displacement taken from discretization and interpolated for a specific Gauss point by Crisfield, Vol. 2, equation (17.93); relates to current iteration step
  vector<LINALG::Matrix<3,1> > thetanew_;
  //!derivative of angle displacement taken from discretization with respect to curve parameter at a specific Gauss point; derivative due to Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding time step
  vector<LINALG::Matrix<3,1> > thetaprimeconv_;
  //!derivative of angle displacement taken from discretization with respect to curve parameter at a specific Gauss point; derivative due to Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding iteration step
  vector<LINALG::Matrix<3,1> > thetaprimeold_;
  //!derivative of angle displacement taken from discretization with respect to curve parameter at a specific Gauss point; derivative due to Crisfield, Vol. 2, equation (17.93); relates to the current iteration step
  vector<LINALG::Matrix<3,1> > thetaprimenew_;

  //! cross section area
  double crosssec_;

   //! cross section area including shear corretion factor
  double crosssecshear_;

  //!moment of inertia of area around y-axis(orthognoal to beam axis in yz-plane)
  double Iyy_;

  //!moment of inertia of area around z-axis(forming cartesian system with beam- and y-axis)
  double Izz_;

  //!polar moment of inertia of area
  double Irr_;

  //! Vector holding value of Jacobi determinant alpha for each Gauss point for underintegration
  vector<double> alpha_;

  //! Vector holding value of Jacobi determinant alpha for complete integration of massmatrix
  vector<double> alphamass_;
  
  //! Vector holding value of Jacobi determinant alpha at nodes
  vector<double> alphanode_;

  //! @name methods for initilization of the element

  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration)
  //this function is templated and receives nnode from Beam3Register::Initialize
  template<int nnode>
  void SetUpReferenceGeometry(const vector<double>& xrefe,const vector<double>& rotrefe);

  //@}

  //! @name auxiliary computation methods for non-additive, large rotation variables

  //!computing spin matrix out of a rotation vector
  void computespin(LINALG::Matrix<3,3>& spin, LINALG::Matrix<3,1> rotationangle, const double& spinscale);

  //!computing a rotation matrix R from a quaternion q
  void quaterniontotriad(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,3>& R);

  //!computing a quaternion from an angle vector
  void angletoquaternion(const LINALG::Matrix<3,1>& theta, LINALG::Matrix<4,1>& q);

  //!computing an angle vector from a quaternion q
  void quaterniontoangle(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,1>& theta);

  //!computing from a quaternion q the Rodrigues parameters omega, Crisfield, Vol. 2, equation (16.79)
  void quaterniontorodrigues(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,1>& omega);

  //!computing a quaternion q from a rotation matrix R
  void triadtoquaternion(const LINALG::Matrix<3,3>& R, LINALG::Matrix<4,1>& q);

  //!quaternion product q12 = q2*q1, Crisfield, Vol. 2, equation (16.71)
  void quaternionproduct(const LINALG::Matrix<4,1>& q1,const LINALG::Matrix<4,1>& q2,LINALG::Matrix<4,1>& q12);

  //!computes inverse quaternion q^{-1} for input quaternion q
  LINALG::Matrix<4,1> inversequaternion(const LINALG::Matrix<4,1>& q);

  //!inverse of matrix H which make rotation variables additive in the limit, eq. (16.94)
  LINALG::Matrix<3,3> Hinv(LINALG::Matrix<3,1> theta);


  //@}

  //! @name Internal calculation methods

  /*!this function performs an update of the central triad as in principle given in Crisfield, Vol. 2, equation (16.65),
   * but by means of a quaternion product and then multiplication of the rotation matrix in (16.69) with the reference triad Tref*/
  void updatetriad(const LINALG::Matrix<3,1>& deltatheta, LINALG::Matrix<3,3>& Tnew, const int numgp);

  //!updating local curvature according to Crisfield, Vol. 2, pages 209 - 210, basically equation (16.148)
  void updatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1>& deltatheta,LINALG::Matrix<3,1>& deltathetaprime, const int numgp);

  //!updating local curvature approximately according to Crisfield, Vol. 2, eqs. (17.72) & (17.73) (only possible for linear elements)
  void approxupdatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1> deltatheta,LINALG::Matrix<3,1> deltathetaprime, const int numgp);

  //!computing stiffness matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.106)
  template<int nnode>
  void computeKsig1(LINALG::Matrix<6*nnode,6*nnode>& Ksig1, const LINALG::Matrix<3,1>& stressn, const LINALG::Matrix<3,1>& stressm, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv);

  //!computing stiffness matrix Ksigma2 according to Crisfield, Vol. 2, equation (17.107a) and (17.107b)
  template<int nnode>
  void computeKsig2(LINALG::Matrix<6*nnode,6*nnode>& Ksig2, const LINALG::Matrix<3,1>& stressn, const LINALG::Matrix<3,3>& S, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv);

  //!computing basis of stiffness matrix according to Crisfield, Vol. 2, equation (17.105)
  template<int nnode>
  void computestiffbasis(const LINALG::Matrix<3,3>& Tnew, const LINALG::Matrix<3,1>& Cm, const LINALG::Matrix<3,1>& Cb, const LINALG::Matrix<3,3>& S, LINALG::Matrix<6*nnode,6*nnode>& stiffmatrix, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv);

  //!calculation of nonlinear stiffness and mass matrix
  template<int nnode>
  void b3_nlnstiffmass( ParameterList& params,
                        vector<double>&           vel,
                        vector<double>&           disp,
                        Epetra_SerialDenseMatrix* stiffmatrix,
                        Epetra_SerialDenseMatrix* massmatrix,
                        Epetra_SerialDenseVector* force);

  //! lump mass matrix
  template<int nnode>
  void lumpmass(Epetra_SerialDenseMatrix* emass);
  
  //!determine Gauss rule from required type of integration
  DRT::UTILS::GaussRule1D MyGaussRule(int nnode, IntegrationType integrationtype);

  //@}
  
  //! @name Methods for Brownian dynamics simulations

   //!compute damping matrix and statistical forces according to fluctuation dissipation theorem
   template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
   void CalcBrownian(ParameterList& params,  //!<parameter list
                     const vector<double>&           vel,  //!< element velocity vector
                     const vector<double>&           disp, //!<element disp vector
                     Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                     Epetra_SerialDenseVector* force);//!< element internal force vector  
   
   //!computes damping coefficients for employed friction model
   void MyDampingConstants(ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel);
   
   //!computes the number of different random numbers required in each time step for generation of stochastic forces
   int HowManyRandomNumbersINeed();
   
   //!computes rotational damping forces and stiffness
   template<int nnode> //number of nodes
   void MyRotationalDamping(ParameterList& params,  //!<parameter list
                           const vector<double>&           vel,  //!< element velocity vector
                           const vector<double>&           disp, //!<element disp vector
                           Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                           Epetra_SerialDenseVector* force);//!< element internal force vector
   
   //!computes translational damping forces and stiffness
   template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
   void MyTranslationalDamping(ParameterList& params,  //!<parameter list
                             const vector<double>&           vel,  //!< element velocity vector
                             const vector<double>&           disp, //!<element disp vector
                             Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                             Epetra_SerialDenseVector* force);//!< element internal force vector
   
   //!computes stochastic translational forces and resulting stiffness
   template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
   void MyStochasticForces(ParameterList& params,  //!<parameter list
                         const vector<double>&           vel,  //!< element velocity vector
                         const vector<double>&           disp, //!<element disp vector
                         Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                         Epetra_SerialDenseVector* force);//!< element internal force vector
   
   //!computes stochastic rotational moments and resulting stiffness
   template<int nnode, int randompergauss> //number of nodes, number of random numbers required per Gauss point
   void MyStochasticMoments(ParameterList& params,  //!<parameter list
                         const vector<double>&           vel,  //!< element velocity vector
                         const vector<double>&           disp, //!<element disp vector
                         Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                         Epetra_SerialDenseVector* force);//!< element internal force vector
   
   //!computes velocity of background fluid and gradient of that velocity at a certain evaluation point in the physical space and adds respective terms to internal forces and damping matrix
   template<int ndim> //number of dimensions of embedding space
   void MyBackgroundVelocity(ParameterList& params,  //!<parameter list
                             const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                             LINALG::Matrix<ndim,1>& velbackground,  //!< velocity of background fluid
                             LINALG::Matrix<ndim,ndim>& velbackgroundgrad); //!<gradient of velocity of background fluid

   //@}


  // don't want = operator
  Beam3& operator = (const Beam3& old);


}; // class Beam3



// << operator
ostream& operator << (ostream& os, const DRT::Element& ele);


//------------- class Beam3Register: -------------------------------------------

/*!
\brief A register for Beam3 element

\author cyron (cyron@lnm.mw.tum.de)
*/
class Beam3Register : public DRT::ElementRegister
{
public:

  //! @name Friends
  friend class Beam3;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Beam3Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Beam3Register

  */
  Beam3Register(const DRT::ELEMENTS::Beam3Register& old);

  /*!
  \brief Deep copy this instance of Beam3Register and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Beam3Register* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Beam3Register;}

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
  virtual ~Beam3Register();

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

  The method performs all steps needed to init the Beam3 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Beam3Register
} // namespace Elements
} // namespace DRT




#endif  // #ifndef BEAM_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM3
