/*!----------------------------------------------------------------------
\file beam3ii.H

\brief 3D nonlinear Reissner beam element oy type II (according to "The
 interpolation of rotations and its application to finite element models
  of geometrically exact rods", Romero 2004)

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/

// administrative line in order to make a difference between BACI and Ccarat files
#ifdef CCADISCRET

//header file only included if not yet included!
#ifndef BEAM3II_H
#define BEAM3II_H


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
\brief 3D nonlinear Reissner beam element implemented according to the following sources: 
Jelenic, Crisfield, 1999, "Geometrically exact 3D beam theory: implementations of a strain-invariant finite element for statics and dynamics",
Crisfield, Jelenic, 1998, "Objectivity of strain measures in the geometrically exact three dimensional beam theory and its finite element implementation",
Romero, 2004, "The interpolation of rotations and its application to finite element models of geometrically exact rods",
Crisfield, 2003, "Non-linear Finite Element Analysis of Solids and Structures", Volume 2

The element has been implemented with linear basis functions only so far.

author cyron (cyron @lnm.mw.tum.de)
*/
class Beam3ii : public DRT::Element
{
public:

  //! @name Friends
  friend class Beam3iiRegister;
  friend class StatMechManager;
  friend void StatMechManager::SetCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,const std::map<int,LINALG::Matrix<3,1> >& currentpositions,const std::map<int,LINALG::Matrix<3,1> >& currentrotations);
  friend void StatMechManager::PeriodicBoundaryBeam3iiInit(DRT::Element* element);
  friend void StatMechManager::GmshOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step);
  friend void StatMechManager::GmshOutputPeriodicBoundary(const LINALG::SerialDenseMatrix& coord, const double& color, std::stringstream& gmshfilecontent, int eleid);
  friend void StatMechTime::Integrate();
  friend StatMechTime::StatMechTime(ParameterList& params,DRT::Discretization& dis,LINALG::Solver& solver,IO::DiscretizationWriter& output);
  
  

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3ii(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3ii(const Beam3ii& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3ii();

  /*!
  \brief Deep copy this instance of Beam3ii and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return ParObject_Beam3ii; }

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
  
  /*!
  \brief Get moment of inertia of area Iyy_ of this element
  */
  const double& Iyy() const { return Iyy_;};
  
  /*!
  \brief Get moment of inertia of area Izz_ of this element
  */
  const double& Izz() const { return Izz_;}
  
  /*!
  \brief Get polar moment of inertia of area Irr_ of this element
  */
  const double& Irr() const { return Irr_;}


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

  //! action parameters recognized by Beam3ii
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

  //! nodal triads in quaternion form at the end of the preceeding time step due to Crisfield, Vol. 2, equation (17.65)
  vector<LINALG::Matrix<4,1> > Qconv_;
  //! nodal triads in quaternion form at the end of the preceeding iteration step due to Crisfield, Vol. 2, equation (17.65)
  vector<LINALG::Matrix<4,1> > Qold_;
  //! nodal triads in quaternion form during the current iteration step due to Crisfield, Vol. 2, equation (17.65)
  vector<LINALG::Matrix<4,1> > Qnew_;

  //!rotation angle displacement at nodes; relates to the end of the preceeding time step
  vector<LINALG::Matrix<3,1> > dispthetaconv_;
  //!rotation angle displacement at nodes; relates to the end of the preceeding iteration step
  vector<LINALG::Matrix<3,1> > dispthetaold_;
  //!rotation angle displacement at nodes; relates to current configuration
  vector<LINALG::Matrix<3,1> > dispthetanew_;
  
  //!initial curvature at Gauss points for elasticity (corresponding to \Lambda_0^t \Labmda'_0 in eq. (3.5), Crisfield 1999
  vector<LINALG::Matrix<3,1> > kapparef_;
  
  //!node I for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
  int nodeI_;
  
  //!node J for determination of reference triad, eq. (3.9), (3.10), Jelenic 1999
  int nodeJ_;


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

  //! Vector holding value of Jacobi determinant jacobi for each Gauss point for underintegration
  vector<double> jacobi_;

  //! Vector holding value of Jacobi determinant jacobi for complete integration of massmatrix
  vector<double> jacobimass_;
  
  //! Vector holding value of Jacobi determinant jacobi at nodes
  vector<double> jacobinode_;

  //! @name methods for initilization of the element

  //!sets up from current nodal position all geometric parameters (considering current position as reference configuration)
  //this function is templated and receives nnode from Beam3iiRegister::Initialize
  template<int nnode>
  void SetUpReferenceGeometry(const vector<double>& xrefe,const vector<double>& rotrefe, const bool secondinit = false);

  //@}

  //! @name Internal calculation methods

  //!compute derivative r' (rprime) of line of centroids at Gauss point x with respect to Jelenic 1999 , eq. (2.12)
  template<int nnode, int dof>
  void curvederivative(const vector<double>& disp, const LINALG::Matrix<1,nnode> deriv, LINALG::Matrix<3,1>& rprime, const double& jacobi);
  
  //!compute convected strain at a certain Gauss point with triad rotmat according to Crisfield 1999, eq. (3.4) and eq. (4.9) 
  void computestrain(const LINALG::Matrix<3,1>& rprime, const LINALG::Matrix<3,3>& Lambda,
                     LINALG::Matrix<3,1>& gamma, LINALG::Matrix<3,1>& kappa);
  
  //!computing convected stresses from convected strains vector and returing also convected constitutive matrix between both
  void strainstress(const LINALG::Matrix<3,1>& gamma, const LINALG::Matrix<3,1>& kappa,
                    LINALG::Matrix<3,1>& stressN, LINALG::Matrix<3,3>& CN,
                    LINALG::Matrix<3,1>& stressM, LINALG::Matrix<3,3>& CM);
  
  //!push forward strains, stresses and constitutive matrices to their spatial counterparts by rotation matrix Lambda according to Romero 2004, eq. (3.10)
  void pushforward(const LINALG::Matrix<3,3>& Lambda,
                  const LINALG::Matrix<3,1>& stressN , const LINALG::Matrix<3,3>& CN,
                  const LINALG::Matrix<3,1>& stressM , const LINALG::Matrix<3,3>& CM,
                  LINALG::Matrix<3,1>& stressn, LINALG::Matrix<3,3>& cn,
                  LINALG::Matrix<3,1>& stressm, LINALG::Matrix<3,3>& cm);

  //!compute convected strain at certain Gauss point with triad rotmat according to Crisfield 1999, eq. (3.4) and eq. (4.9)
  void computestrain(const LINALG::Matrix<3,1> Psil, const LINALG::Matrix<3,1> Psilprime,
                     const double& xi, const LINALG::Matrix<3,1>& rprime, const LINALG::Matrix<3,3>& rotmat,
                     LINALG::Matrix<3,1>& gamma, LINALG::Matrix<3,1>& kappa);
  
  //!evaluates basis functions for a given type of Gauss rule  and triads at all Gauss points
  template<int nnode>
  void evaluatebasisfunctionsandtriads(const DRT::UTILS::IntegrationPoints1D gausspoints,
                                       vector<LINALG::Matrix<1,nnode> >& I,
                                       vector<LINALG::Matrix<1,nnode> >& Iprime,
                                       vector<vector<LINALG::Matrix<3,3> > >& Itilde,
                                       vector<vector<LINALG::Matrix<3,3> > >& Itildeprime,
                                       vector<LINALG::Matrix<3,3> >& Lambda);
  
  
  
  //!compute \tilde{I}^i in (3.18), page 152, Jelenic 1999, for all nodes i at a certaom Gauss point
  template<int nnode>
  void computeItilde(const LINALG::Matrix<3,1>& Psil, const LINALG::Matrix<3,1>& Psilprime, vector<LINALG::Matrix<3,3> >& Itilde,
                     const LINALG::Matrix<3,1>& phiIJ, const LINALG::Matrix<3,3>& Lambdar,
                     const vector<LINALG::Matrix<3,1> >& Psili, const LINALG::Matrix<1,nnode>& funct);
  
  //!compute \tilde{I}^{i'} in (3.19), page 152, Jelenic 1999, for all nodes i at a certaom Gauss point
  template<int nnode>
  void computeItildeprime(const LINALG::Matrix<3,1>& Psil, const LINALG::Matrix<3,1>& Psilprime, vector<LINALG::Matrix<3,3> >& Itildeprime,
                          const LINALG::Matrix<3,1>& phiIJ, const LINALG::Matrix<3,3>& Lambdar,
                          const vector<LINALG::Matrix<3,1> >& Psili, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv);
  
  //!compute matrix v_I as outlined in the equations above (3.15) on page 152 of Jelenic 1999
  LINALG::Matrix<3,3> vI(const LINALG::Matrix<3,1>& phiIJ);
  
  //!compute matrix v_J as outlined in the equations above (3.15) on page 152 of Jelenic 1999
  LINALG::Matrix<3,3> vJ(const LINALG::Matrix<3,1>& phiIJ);

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
   
   //!shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions
   template<int nnode, int ndim> //number of nodes, number of dimensions
   void NodeShift(ParameterList& params,  //!<parameter list
                 vector<double>&           disp); //!<element disp vector

   //@}


  // don't want = operator
  Beam3ii& operator = (const Beam3ii& old);


}; // class Beam3ii



// << operator
ostream& operator << (ostream& os, const DRT::Element& ele);


//------------- class Beam3iiRegister: -------------------------------------------

/*!
\brief A register for Beam3ii element

\author cyron (cyron@lnm.mw.tum.de)
*/
class Beam3iiRegister : public DRT::ElementRegister
{
public:

  //! @name Friends
  friend class Beam3ii;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Beam3iiRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Beam3iiRegister

  */
  Beam3iiRegister(const DRT::ELEMENTS::Beam3iiRegister& old);

  /*!
  \brief Deep copy this instance of Beam3iiRegister and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Beam3iiRegister* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Beam3iiRegister;}

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
  virtual ~Beam3iiRegister();

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

  The method performs all steps needed to init the Beam3ii elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Beam3iiRegister
} // namespace Elements
} // namespace DRT


#endif  // #ifndef BEAM_H
#endif  // #ifdef CCADISCRET
