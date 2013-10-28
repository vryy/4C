/*!----------------------------------------------------------------------
\file beam2.H

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
// administrative line in order to make a difference between BACI and Ccarat files

//header file only included if not yet included!
#ifndef BEAM2R_H
#define BEAM2R_H


#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_statmech.H"

namespace DRT
{

namespace ELEMENTS
{

class Beam2rType : public DRT::ElementType
{
public:

  std::string Name() const { return "Beam2rType"; }

  static Beam2rType & Instance() { return instance_; };

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

  static Beam2rType instance_;
};

/*!
\brief two dimensional nonlinear beam element using Reissner`s beam theory according to Crisfield, Vol.1, section 7.4
\The element is capable of using higher order functions from linear to quartic. Please make sure you put the nodes in the right order
\in the input file.
\LIN2    1---2
\LIN3  1---3---2
\LIN4    1---4---2---3
\LIN5  1---5---2---3---4
\author cyron (cyron @lnm.mw.tum.de)
*/
class Beam2r : public DRT::Element
{
public:

  //! @name Friends
  friend class Beam2rType;
  //friend StatMechTime::StatMechTime(ParameterList& params,DRT::Discretization& dis,LINALG::Solver& solver,IO::DiscretizationWriter& output);


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam2r(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam2r(const Beam2r& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam2r();

  /*!
  \brief Deep copy this instance of Beam2r and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return Beam2rType::Instance().UniqueParObjectId(); }

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
  { return Beam2rType::Instance(); }

  //@}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return 1;}


  /*!
  \brief Get vector of RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();


  /*!
  \brief Get number of degrees of freedom of a certain node
  */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    /*note: this is not necessarily the number of DOF assigned to this node by the discretization finally, but only
    *the number of DOF requested for this node by this element; the discretization will finally assign the maximal
    *number of DOF to this node requested by any element connected to this node*/
    return 3;
  }

  /*!
  \brief Get number of degrees of freedom per element not including nodal degrees of freedom
  */
  virtual int NumDofPerElement() const { return 0;}

  /*!
  \brief computes the number of different random numbers required in each time step for generation of stochastic forces
  */
  int HowManyRandomNumbersINeed();

  /*!
  \brief Print this element
  */
  virtual void Print(std::ostream& os) const;


  //@}

  //! @name Construction

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
  virtual int Evaluate(Teuchos::ParameterList&            params,
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




  //@}


private:

  //! action parameters recognized by beam2r
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
    calc_struct_reset_istep, //!< reset elementwise internal variables during iteration to last converged state
    calc_struct_ptcstiff
  };

  //! kind of integration to be performed
  enum IntegrationType
  {
    gaussexactintegration,
    gaussunderintegration,
    lobattointegration
  };

  //! cross section area
  double crosssec_;

   //! cross section area including shear correction factor
  double crosssecshear_;

  //! gaussrule_ will be initialized automatically to an underintegration rule in input.cpp
  DRT::UTILS::GaussRule1D gaussrule_;

  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;

  //! moment of inertia of area
  double mominer_;

  //! Vector holding jacobi for each gausspoint
  std::vector<double> jacobi_;

  //! Vector holding jacobi for each gausspoint for complete integration of massmatrix
  std::vector<double> jacobimass_;

  //! Vector holding value of Jacobi determinant jacobi at nodes
  std::vector<double> jacobinode_;

  /*! Vector holding theta0 for each gausspoint. The current implementation sets the initial
  director angles to zero, which can be done without loss of generality since any initial rotation
  cancels out in the equations (7.131). As a consequence the current nodal directors do not appear
  in the element, but only the values of the current directorat the GP*/
  std::vector<double> theta0_;


  //!internal geometric calculations according to Crisfield Vol.1 (7.118)-(7.121)
  //note: this and a few more functions are templated so we can use the same element as a
  //    higher order element. For templating cf. Practical C++ Programming p.419++
  template<int nnode>
  void local_aux(LINALG::Matrix<3,3*nnode>& Bcurr_gp,//rotation matrix
                 LINALG::Matrix<3*nnode,1>& rcurr_gp,//local vector
                 LINALG::Matrix<3*nnode,1>& zcurr_gp,//local vector at GP
                 LINALG::Matrix<3*nnode,1>& scurr_gp,//local vector
                 const double&              theta_gp, //current angle theta at GP
                 const double&              c1,// geometric variable
                 const double&              c2,// geometric variable
                 LINALG::Matrix<1,nnode>&   funct,//holds values of shape functions at GP
                 LINALG::Matrix<1,nnode>&   deriv);//holds derivates of shape functions at GP

  //! calculation of nonlinear stiffness and mass matrix
  template<int nnode>
  void nlnstiffmass(Teuchos::ParameterList&   params,
                    std::vector<double>&      vel,  //!< element velocity vector
                    std::vector<double>&      disp,  //!< element displacement vector
                    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                    Epetra_SerialDenseMatrix* massmatrix,  //!< element mass matrix
                    Epetra_SerialDenseVector* force);  //!< element internal force vector

  //!this function transforms the consistent massmatrix into a lumped massmatrix. Only called when needed
  template<int nnode>
  void lumpedmass(Epetra_SerialDenseMatrix* massmatrix);

  //! determine Gauss rule from required type of integration
  DRT::UTILS::GaussRule1D MyGaussRule(int nnode, IntegrationType integrationtype);



  //!sets up element reference geometry for reference nodal position vector xrefe (may be used also after simulation start)
  template<int nnode>
  void SetUpReferenceGeometry(const std::vector<double>& xrefe);

  //! @name Methods for Brownian dynamics simulations

   //!compute damping matrix and statistical forces according to fluctuation dissipation theorem
   template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
   void CalcBrownian(Teuchos::ParameterList&    params,  //!<parameter list
                     const std::vector<double>& vel,  //!< element velocity vector
                     const std::vector<double>& disp, //!<element disp vector
                     Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                     Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes damping coefficients for employed friction model
   void MyDampingConstants(Teuchos::ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel, const DRT::UTILS::IntegrationPoints1D gausspoints, const std::vector<double> jacobi);

   //!computes translational damping forces and stiffness
   template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
   void MyTranslationalDamping(Teuchos::ParameterList&  params,  //!<parameter list
                             const std::vector<double>& vel,  //!< element velocity vector
                             const std::vector<double>& disp, //!<element disp vector
                             Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                             Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes stochastic translational forces and resulting stiffness
   template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
   void MyStochasticForces(Teuchos::ParameterList&  params,  //!<parameter list
                         const std::vector<double>& vel,  //!< element velocity vector
                         const std::vector<double>& disp, //!<element disp vector
                         Epetra_SerialDenseMatrix*  stiffmatrix,  //!< element stiffness matrix
                         Epetra_SerialDenseVector*  force);//!< element internal force vector

   //!computes velocity of background fluid and gradient of that velocity at a certain evaluation point in the physical space and adds respective terms to internal forces and damping matrix
   template<int ndim> //number of dimensions of embedding space
   void MyBackgroundVelocity(Teuchos::ParameterList&        params,  //!<parameter list
                             const LINALG::Matrix<ndim,1>&  evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                             LINALG::Matrix<ndim,1>&        velbackground,  //!< velocity of background fluid
                             LINALG::Matrix<ndim,ndim>&     velbackgroundgrad); //!<gradient of velocity of background fluid

   //@}


  // don't want = operator
  Beam2r& operator = (const Beam2r& old);

}; // class Beam2r


// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);




} // namespace Elements
} // namespace DRT




#endif  // #ifndef BEAM_H
