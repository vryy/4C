/*----------------------------------------------------------------------------*/
/*!
\file dcohesive.H

\brief DcohesiveType

\level 3

\maintainer N. A.
*/
/*----------------------------------------------------------------------------*/

#ifndef DCOHESIVE_H
#define DCOHESIVE_H


#include "Teuchos_RCP.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include "../drt_inpar/inpar_crack.H"

namespace DRT
{
namespace ELEMENTS
{

class DcohesiveType : public DRT::ElementType
{
public:

  std::string Name() const { return "DcohesiveType"; }

  static DcohesiveType & Instance();

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

  static DcohesiveType instance_;
};

/*!
\brief two dimensional nonlinear corotational Timoshenko beam element according to Crisfield, Vol.1, section 7.3

\author cyron (cyron @lnm.mw.tum.de)
*/
class Dcohesive : public DRT::Element
{
public:

  //! @name Friends
  friend class DcohesiveType;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Dcohesive(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Dcohesive(const Dcohesive& old);

  /*!
  \brief Destructor
  */
  virtual ~Dcohesive();

  /*!
  \brief Deep copy this instance of Beam2 and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return DcohesiveType::Instance().UniqueParObjectId(); }

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
  { return DcohesiveType::Instance(); }

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
  /*bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef);*/


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
  //@}

  /*!
  \brief set reference normal vector which is used to calculate normal and tangential separations
   */
  virtual void setRefNormal( LINALG::Matrix<3,1> norm ){refNormal_ = norm;}

  /*!
  \brief set area over which this spring acts
   */
  virtual void setArea( double area ){ area_ = area; }

  virtual void setCohesiveStrengths( std::vector<double> coh ){ coheStrength_ = coh; }

  virtual void setFractureEnergies( std::vector<double> frac ){ fracEnergy_ = frac; }

  virtual void setCohesiveLaw( enum INPAR::CRACK::tractionSeparation tracLaw ){ tracSepLaw_ = tracLaw; }

  bool isFailed() {return failNorm_ or failTang_;}

private:
  void InitializeElement();

  double Macaulay( double a, double b);

  void ForceStiffnessSpring( const std::vector<double>& deflec,
                             std::vector<double>& forcSpr,
                             std::vector<double>& stifSpr);
  /*!
  \brief Evaluate this element with Discrete Damage Zone Modeling (DDZM) concepts
   */
  int EvaluateDDZM(Teuchos::ParameterList&   params,
                   DRT::Discretization&      discretization,
                   std::vector<int>&         lm,
                   Epetra_SerialDenseMatrix& elemat1,
                   Epetra_SerialDenseMatrix& elemat2,
                   Epetra_SerialDenseVector& elevec1,
                   Epetra_SerialDenseVector& elevec2,
                   Epetra_SerialDenseVector& elevec3);

  /*!
  \brief Evaluate this element with Discrete Cohesive Zone Modeling (DDZM) concepts
   */
  int EvaluateDCZM(Teuchos::ParameterList&   params,
                   DRT::Discretization&      discretization,
                   std::vector<int>&         lm,
                   Epetra_SerialDenseMatrix& elemat1,
                   Epetra_SerialDenseMatrix& elemat2,
                   Epetra_SerialDenseVector& elevec1,
                   Epetra_SerialDenseVector& elevec2,
                   Epetra_SerialDenseVector& elevec3);



  //! action parameters recognized by cohesive elements
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

  /*//! kind of integration to be performed
  enum IntegrationType
  {
    gaussexactintegration,
    gaussunderintegration,
    lobattointegration
  };*/

  enum damageLaw
  {
    linear,
    exponential
  };

  //! whether this spring already failed in normal separation
  bool failNorm_;

  //! whether this spring already failed in shear separation
  bool failTang_;

  //! reference normal vector used to calculate normal and tangential displacements
  LINALG::Matrix<3,1> refNormal_;

  //! spring stiffness is a function of damage which denotes whether it is elastic loading or
  // yielding already started
  LINALG::Matrix<2,1> damage_;

  //! damage in the previous time step -- used to enforce irreversibility
  LINALG::Matrix<2,1> damagePrev_;

  //! area over which this spring is assumed to act
  double area_;

  //! cohesive strength of material in normal and shear
  std::vector<double> coheStrength_;

  //! fracture energy of material in normal and shear
  std::vector<double> fracEnergy_;

  //! critical displacements of the spring -- displacement when traction=cohesive_strength
  std::vector<double> critDisp_;

  enum INPAR::CRACK::tractionSeparation model_;

  enum INPAR::CRACK::tractionSeparation tracSepLaw_;

  //! constants that define the traction-separation curve of PPR model
  double alfa_ppr_, beta_ppr_;

  //! initial slope indicator in normal and tangential direction for PPR model
  std::vector<double> iniSlop_;

  //! non-dimensional exponents of PPR model
  std::vector<double> mn_;

  //! maximum displacement at fracture
  std::vector<double> maxDisp_;

  //! energy constants used in PPR model
  std::vector<double> gamma_;

  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;

  //! length in reference configuration
  double lrefe_;

  //! cross section area
  double crosssec_;

   //! cross section area including shear corretion factor
  double crosssecshear_;

  //! moment of inertia of area
  double mominer_;

  //! number of 2*PI rotations between the aboslute element frame rotation and the angle beta in current iteration step
  int numperiodsnew_;

  //! number of 2*PI rotations between the aboslute element frame rotation and the angle beta in last iteration step
  int numperiodsold_;

  //! number of 2*PI rotations between the aboslute element frame rotation and the angle beta at the end of last time step
  int numperiodsconv_;

  //! absolute angle by which beam has rotated out of x-axis including reference rotation by alpha0_ in current iteration step
  double alphanew_;

  //! absolute angle by which beam has rotated out of x-axis including reference rotation by alpha0_ in last iteration step
  double alphaold_;

  //! absolute angle by which beam has rotated out of x-axis including reference rotation by alpha0_ at the end of last time step
  double alphaconv_;

  //! absolute angle between x-axis and element frame in reference configuration
  double alpha0_;

  //! Vector holding jacobi for each gausspoint
  std::vector<double> jacobi_;

  //! Vector holding jacobi for each gausspoint for complete integration of massmatrix
  std::vector<double> jacobimass_;

  //! Vector holding value of Jacobi determinant jacobi at nodes
  std::vector<double> jacobinode_;

  //! gaussrule_ will be initialized automatically to a 2 point integration rule
  DRT::UTILS::GaussRule1D gaussrule_;

  // internal calculation methods

  std::vector<double> calculateDeflections( DRT::Discretization& discret, std::vector<int>& lm  );

  void computeConstants( double * b, double * initStiff, double * critForce, double * critDisp );

  void computeDamage( double * b,
                      double * initStiff,
                      double * critForce,
                      double * critDisp,
                      const std::vector<double>& deflec,
                      std::vector<double>& forceSpr );

  //! calculation of nonlinear stiffness and mass matrix
  void nlnstiffmass(Teuchos::ParameterList&   params,
                    std::vector<int>&         lm,
                    std::vector<double>&      vel,  //!< element velocity vector
                    std::vector<double>&      disp,  //!< element displacement vector
                    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                    Epetra_SerialDenseMatrix* massmatrix,  //!< element mass matrix
                    Epetra_SerialDenseVector* force,  //!< element internal force vector
                    int lumpedmass); //!< if set to one lumped mass matrix is calculated

  //!sets up element reference geomtry for reference nodal position vector xrefe (may be used also after simulation start)
  void SetUpReferenceGeometry(const LINALG::Matrix<4,1>& xrefe);

  //! determine Gauss rule from required type of integration
  //DRT::UTILS::GaussRule1D MyGaussRule(int nnode, IntegrationType integrationtype);



  //!sets up element reference geometry for reference nodal position vector xrefe (may be used also after simulation start)
  template<int nnode>
  void SetUpReferenceGeometry(const std::vector<double>& xrefe);


  // don't want = operator
   Dcohesive& operator = (const Dcohesive& old);

}; // class Beam2


// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);

} // namespace Elements
} // namespace DRT




#endif  // #ifndef BEAM_H
