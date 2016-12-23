/*-----------------------------------------------------------*/
/*!
 \file spring3.H

\brief pseudo transient solution method

\maintainer Dhrubajyoti Mukherjee

\level 3

*/
/*-----------------------------------------------------------*/
// administrative line in order to make a difference between BACI and Ccarat files
//header file only included if not yet included!
#ifndef SPRING3_H
#define SPRING3_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "Sacado.hpp"

#include "../drt_inpar/inpar_browniandyn.H"
typedef Sacado::Fad::DFad<double> FAD;


/*!
 \brief DRT: namespace of the ccarat discretization module

 */
namespace DRT
{

namespace ELEMENTS
{

class Spring3Type: public DRT::ElementType
{
public:

  std::string Name() const
  {
    return "Spring3Type";
  }

  static Spring3Type & Instance()
  {
    return instance_;
  }
  ;

  virtual DRT::ParObject* Create(const std::vector<char> & data);

  virtual Teuchos::RCP<DRT::Element> Create(const std::string eletype,
      const std::string eledistype, const int id, const int owner);

  virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

  virtual int Initialize(DRT::Discretization& dis);

  virtual void NodalBlockInformation(DRT::Element * dwele, int & numdf,
      int & dimns, int & nv, int & np);

  virtual void ComputeNullSpace(DRT::Discretization & dis,
      std::vector<double> & ns, const double * x0, int numdf, int dimns);

  virtual void SetupElementDefinition(
      std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition> > & definitions);

private:

  static Spring3Type instance_;
};

/*!
 \brief three dimensional total Lagrange truss element

 \author mukherjee (mukherjee @lnm.mw.tum.de)
 */
class Spring3: public DRT::Element
{
public:

  //! @name Friends
  friend class Spring3Type;
  //! @name Constructors and destructors and related methods

  /*!
   \brief Standard Constructor

   \param id    (in): A globally unique element id
   \param etype (in): Type of element
   \param owner (in): owner processor of the element
   */
  Spring3(int id, int owner);

  /*!
   \brief Copy Constructor

   Makes a deep copy of a Element
   */
  Spring3(const Spring3& old);

  /*!
   \brief Destructor
   */
  virtual ~Spring3();

  /*!
   \brief Deep copy this instance of Spring3 and return pointer to the copy

   The Clone() method is used by the virtual base class Element in cases
   where the type of the derived class is unknown and a copy-ctor is needed
   .
   */
  virtual DRT::Element* Clone() const;

  /*!
   \brief Get shape type of element
   */
  virtual DiscretizationType Shape() const;

  virtual void SetInitialTangents(std::vector<LINALG::Matrix<4,1> >& initquaternions);

  virtual void TcurrBeam3r(LINALG::Matrix<3,1>& Tcurr1, LINALG::Matrix<3,1>& Tcurr2);

  virtual void TrefBeam3r(LINALG::Matrix<3,1>& Tref1, LINALG::Matrix<3,1>& Tref2);

  /*Print the change in angle of this element*/

  virtual LINALG::Matrix<1,3> DeltaTheta() const;

  virtual void GetCurrTangents(std::vector<double>& disp, std::vector<LINALG::Matrix<3,1> >& Tcurr);

  /*!
   \brief Return unique ParObject id

   Every class implementing ParObject needs a unique id defined at the
   top of drt_parobject.H
   */
  virtual int UniqueParObjectId() const
  {
    return Spring3Type::Instance().UniqueParObjectId();
  }

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
  {
    return Spring3Type::Instance();
  }

  //@}

  //! @name methods for initilization of the element

  //!sets reference coordinates X_ and refenrence length lrefe_ for elements added to the discretization after simulation start
  void SetUpReferenceGeometry(const std::vector<double>& xrefe, const std::vector<double>& rotrefe,
      const bool secondinit = false, const bool reissner= false);

  //@}

  /*!
   \brief Return number of lines to this element
   */
  virtual int NumLine() const
  {
    return 1;
  }

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
    return 3;
  }

  /*!
   \brief Get number of degrees of freedom per element not including nodal degrees of freedom
   */
  virtual int NumDofPerElement() const
  {
    return 0;
  }

  /*!
   \brief Set cross section area of this element
   */
  void SetCrossSec(const double& crosssec);

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
  bool ReadElement(const std::string& eletype, const std::string& distype,
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
  virtual int Evaluate(Teuchos::ParameterList& params,
      DRT::Discretization& discretization, std::vector<int>& lm,
      Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
      Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
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
      DRT::Discretization& discretization, DRT::Condition& condition,
      std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
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
  template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
  void EvaluatePTC(Teuchos::ParameterList& params,
      Epetra_SerialDenseMatrix& elemat1);

  double ReturnNormMoment() const
  {
    //for now constant, since we only implemented 4-noded interpolated element with linear shape functions
    return NormMoment;
  }

  double ReturnNormForce() const
  {
    //for now constant, since we only implemented 4-noded interpolated element with linear shape functions
    return NormForce;
  }

  double ReturnRatioNormForceMoment() const
  {
    //for now constant, since we only implemented 4-noded interpolated element with linear shape functions
    return RatioNormForceMoment;
  }

  //! \brief return reference length
  const double& L0() const
  {
    return lrefe_;
  }

  //! \brief return reference length
  double Lcurr() const
  {
    return lcurr_;
  }

  //! \brief return change in angle between filaments
  double GetDeltaTheta() const
  {
    return abs(deltatheta_(2));
  }

  //! \brief internal force vector
  const Teuchos::RCP<Epetra_SerialDenseVector> InternalForceVector() const
      {
    return f_;
      }

  //@}

private:

  //! action parameters recognized by Spring3
  enum ActionType
  {
    calc_none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_nlnstiff2,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_nlnstifflmass, //!< internal force, its stiffness and lumped mass matrix
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_reset_istep, //!< reset elementwise internal variables
    //!< during iteration to last converged state
    postprocess_stress,
    calc_struct_ptcstiff, //calculate artificial stiffness due to PTC solution strategy
    calc_struct_energy
    //calculate current elastic energy of element
  };

  //! container containing all kind of Spring3 thisandthat
  DRT::Container data_;

  // boolian indicating if the spring crosslinker binds to a reissner filament
  bool FilamentIsReissner_;

  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;

  //!reference node position
  LINALG::Matrix<6, 1> X_;

  //!reference tangent position
  std::vector<LINALG::Matrix<3,1> >  trefNode_;

  //! Vector holding values of ref angles
  std::vector<double> ThetaRef_;

  //!reference tangent position
  LINALG::Matrix<1,3> diff_disp_ref_;

  //Change is angle
  LINALG::Matrix<1,3> deltatheta_;


  //! nodal triads in quaternion during the current iteration step
  std::vector<LINALG::Matrix<4,1> > Qnew_;
  //! nodal triads in quaternion form at the end of the preceeding iteration step
  std::vector<LINALG::Matrix<4,1> > Qold_;
  //!rotation angle displacement at nodes; relates to current configuration
  std::vector<LINALG::Matrix<3,1> > dispthetanew_;
  //!rotation angle displacement at nodes; relates to the end of the preceeding iteration step
  std::vector<LINALG::Matrix<3,1> > dispthetaold_;

  //! material type
  int material_;

  //! length in reference configuration
  double lrefe_;

  //! length in reference configuration
  double lcurr_;

  //! Vector holding value of Jacobi determinant jacobi for complete integration of massmatrix
  std::vector<double> jacobimass_;

  //! vector holding value of Jacobi determinant jacobi at nodes
  std::vector<double> jacobinode_;

  //! cross section area
  double crosssec_;

  // L2 norm of total moment
  double NormMoment;
  // L2 norm of total force
  double NormForce;
  // Ratio of norm of Force to Moment
  double RatioNormForceMoment;

  //! angle of last time step (necessary for PTC scheme)
  LINALG::Matrix<3,1> Theta0_;
  //! Angle of current time step (necessary for PTC scheme)
  LINALG::Matrix<3,1> Theta_;

  //! internal force vector
  Teuchos::RCP<Epetra_SerialDenseVector> f_;

  //!calculation of nonlinear stiffness and mass matrix switching between total lagrange and enginerring strains
  void t3_nlnstiffmass(Teuchos::ParameterList& params, std::vector<double>& vel,
      std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
      Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force);

  //!calculation of elastic energy
  void t3_energy(Teuchos::ParameterList& params, std::vector<double>& disp,
      Epetra_SerialDenseVector* intenergy);

  // Calculate change in angle from reference configuration
  void CalcDeltaTheta(std::vector<double>&  disp, LINALG::Matrix<1,3>& thetacurr);

  //   Calculate torsional stiffness matrices
  void torsion_stiffmass(Teuchos::ParameterList&   params,
      std::vector<double>&      disp,
      Epetra_SerialDenseMatrix* stiffmatrix,
      Epetra_SerialDenseVector* force);


  // Calculate torsional stiffness matrices using automatic differentiation
  void FADMyTorsionalStiffTangentCos(Teuchos::ParameterList&   params,
      double theta_0,
      LINALG::Matrix<1,3> & tcurrNode1,
      LINALG::Matrix<1,3> & tcurrNode2,
      Epetra_SerialDenseMatrix & TorStiffmatrix,
      Epetra_SerialDenseVector & TorForce);

  void MyTorsionalStiffTangentDot(Teuchos::ParameterList&   params,
      LINALG::Matrix<1,3> & tangentcurr1,
      LINALG::Matrix<1,3> & tangentcurr2,
      Epetra_SerialDenseMatrix & TorStiffmatrix,
      Epetra_SerialDenseVector & TorForce);

  void MyTorsionalStiffatTangentSin(Teuchos::ParameterList&   params,
      double theta,
      double deltatheta,
      LINALG::Matrix<1,3> & tcurr1,
      LINALG::Matrix<1,3> & tcurr2,
      Epetra_SerialDenseMatrix & TorStiffmatrix,
      Epetra_SerialDenseVector & TorForce);

  void FADThetaLinearisation(LINALG::Matrix<1,3> & diff_disp,
      LINALG::Matrix<1,3> & tcurr,
      LINALG::Matrix<1,3> & A,
      LINALG::Matrix<1,3> & B);

  // Calculate torsional stiffness matrices and forces between tangents of beam
  void MyTorsionalStiffTangentCos(Teuchos::ParameterList&   params,
      double theta,
      double deltatheta,
      LINALG::Matrix<1,3> & tcurr1,
      LINALG::Matrix<1,3> & tcurr2,
      Epetra_SerialDenseMatrix & TorStiffmatrix,
      Epetra_SerialDenseVector & TorForce);

  /*!
   * \brief calculate force, nonlinear stiffness and mass matrix using a engineering strain measure.
   *
   * Unlike the fully nonlinear implementation of #t3_nlnstiffmass_totlag, this function uses \f$\varepsilon=\Delta d / d\f$
   * as strain measure.
   *
   */
  void t3_nlnstiffmass_spring(const LINALG::Matrix<1,6>&      DummyDisp,
      Epetra_SerialDenseMatrix& DummyStiffMatrix,
      Epetra_SerialDenseVector& DummyForce);

  //! lump mass matrix
  void t3_lumpmass(Epetra_SerialDenseMatrix* emass);

  //@}

  //! @name Methods for Brownian dynamics simulations

  //!compute damping matrix and statistical forces according to fluctuation dissipation theorem

  //!computes velocity of background fluid and gradient of that velocity at a certain evaluation point in the physical space and adds respective terms to internal forces and damping matrix
  template<int ndim> //number of dimensions of embedding space
  void MyBackgroundVelocity(Teuchos::ParameterList& params, //!<parameter list
      const LINALG::Matrix<ndim, 1>& evaluationpoint, //!<point at which background velocity and its gradient has to be computed
      LINALG::Matrix<ndim, 1>& velbackground, //!< velocity of background fluid
      LINALG::Matrix<ndim, ndim>& velbackgroundgrad); //!<gradient of velocity of background fluid

  //!shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions
  template<int nnode, int ndim> //number of nodes, number of dimensions
  void NodeShift(Teuchos::ParameterList& params, //!<parameter list
      std::vector<double>& disp); //!<element disp vector

  //@}
  //! Calculates Signum function of a given "value"
  double Signum(double value);
  // don't want = operator
  Spring3& operator =(const Spring3& old);

};
// class Spring3

// << operator
std::ostream& operator <<(std::ostream& os, const DRT::Element& ele);

} // namespace Elements
} // namespace DRT

#endif  // #ifndef SPRING3_H
