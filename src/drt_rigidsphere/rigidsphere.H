/*----------------------------------------------------------------------------*/
/*!
\file rigidsphere.H

\brief spherical particle element for brownian dynamics

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#ifndef RIGIDSPHERE_H
#define RIGIDSPHERE_H


#include "Teuchos_RCP.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include "../drt_beaminteraction/beaminteraction_calc_utils.H"

// forward declaration ...
namespace STR {
namespace ELEMENTS {
  class ParamsInterface;
} // namespace ELEMENTS
} // namespace STR

namespace DRT
{
namespace ELEMENTS
{

class RigidsphereType : public DRT::ElementType
{
public:

  std::string Name() const { return ("RigidsphereType"); }

  static RigidsphereType & Instance();

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

  static RigidsphereType instance_;
};

/*!
\brief Spherical particle element for brownian dynamics

author meier (meier@lnm.mw.tum.de)
*/
class Rigidsphere : public DRT::Element
{
public:

  //! @name Friends
  friend class RigidsphereType;


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Rigidsphere(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Rigidsphere(const Rigidsphere& old);

  /*!
  \brief Destructor
  */
  virtual ~Rigidsphere();

  /*!
  \brief Deep copy this instance of Beam3eb and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return (RigidsphereType::Instance().UniqueParObjectId()); }

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
  { return (RigidsphereType::Instance()); }

  //@}

  /*!
  \brief Return number of lines to this element
  */
  virtual int NumLine() const {return (1);}

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
      return (3);
  }

  /*!
  \brief Get number of degrees of freedom per element not including nodal degrees of freedom
  */
  virtual int NumDofPerElement() const { return (0);}

  /*!
  \brief Print this element
  */
  virtual void Print(std::ostream& os) const;

  /*!
  \brief Return radius of sphere
  */
  const double& Radius() const {return (radius_);};

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
                              Epetra_SerialDenseMatrix* elemat1 = NULL){return (0);};

  /*!
  \brief Evaluate mass matrix
  */
  void nlnstiffmass(Teuchos::ParameterList& params,
                    std::vector<double>& acc,
                    std::vector<double>& vel,
                    std::vector<double>& disp,
                    Epetra_SerialDenseMatrix* stiffmatrix,
                    Epetra_SerialDenseMatrix* massmatrix,
                    Epetra_SerialDenseVector* force,
                    Epetra_SerialDenseVector* inertia_force);

  /*! \brief set the parameter interface ptr for the solid elements
   *
   *  \param p (in): Parameter list coming from the time integrator.
   *
   *  \author hiermeier
   *  \date 04/16 */
  virtual void SetParamsInterfacePtr(const Teuchos::ParameterList& p);

  /*! \brief returns true if the parameter interface is defined and initialized, otherwise false
   *
   *  \author hiermeier
   *  \date 04/16 */
  virtual inline bool IsParamsInterface() const
  { return (not interface_ptr_.is_null()); }

  /*! \brief get access to the parameter interface pointer
   *
   *  \author hiermeier
   *  \date 04/16 */
  virtual Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr();
  //@}

  //! @name methods for biopolymer network simulations

  //!computes the number of different random numbers required in each time step for generation of stochastic forces
  int HowManyRandomNumbersINeed();

  /// \brief get generalized interpolation matrix which yields the variation of the position
  virtual void GetGeneralizedInterpolationMatrixVariationsAtXi(
      LINALG::SerialDenseMatrix& Ivar,
      const double&              dummy1,
      const std::vector<double>& dummy2) const;

  /// \brief get generalized interpolation matrix which yields the increments of the position
  virtual void GetGeneralizedInterpolationMatrixIncrementsAtXi(
      LINALG::SerialDenseMatrix&  Iinc,
      const double&               dummy1,
      const std::vector<double>&  dummy2) const;

  /** \brief get linearization of the product of (generalized interpolation matrix for variations
   *         (see above) and applied force vector) with respect to the primary DoFs of this element
   *
   *  \author grill
   *  \date 01/17 */
  virtual void GetStiffmatResultingFromGeneralizedInterpolationMatrixAtXi(
      LINALG::SerialDenseMatrix&       stiffmat,
      const double&                    xi,
      const std::vector<double>&       disp,
      const LINALG::SerialDenseVector& force) const
  {
    // nothing to do here
    stiffmat.Zero();
  }

  /// return binding spot xi function dummy
  double GetBindingSpotXi(int dummy ) const
  {
    return 0.0;
  }

  /** \brief add bond
   *
   *  \author eichinger
   *  \date 06/17 */
  void AddBond( std::pair< int, int > newbondpartner ) { bondstatus_.insert( newbondpartner ); };

   /** \brief get number of bonds
   *
   *  \author eichinger
   *  \date 06/17 */
  int GetNumberOfBonds( int elepartnergid ) const { return static_cast<int>( bondstatus_.size() ); }

  /** \brief get bond status of sphere element
   *
   *  \author eichinger
   *  \date 06/17 */
  std::set< std::pair< int, int > , BEAMINTERACTION::UTILS::StdPairComparatorOrderCounts > const &
      GetBondStatus() const { return bondstatus_;}

  /** \brief check if bond exists
   *
   *  \author eichinger
   *  \date 06/17 */
  bool DoesBondExist( std::pair< int, int > bondpartner ) const
  {
    return ( bondstatus_.find(bondpartner) != bondstatus_.end() ) ? true : false;
  }


  //@}

protected:
  /** \brief get access to the interface
   *
   *  \author hiermeier
   *  \date 04/16 */
  inline STR::ELEMENTS::ParamsInterface& ParamsInterface()
  { if (not IsParamsInterface())
      dserror("The interface ptr is not set!");
    return *interface_ptr_; }

private:
  /*! \brief interface ptr
   *
   *  data exchange between the element and the time integrator. */
  Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

  //! radius of the sphere
  double radius_;

  //! density of the sphere
  double rho_;

  //! @name variables for biopolymer network simulations

  /// holds gids of (beam) elements to which sphere is bonded to (size equals number of bonds)
  std::set< std::pair< int, int > , BEAMINTERACTION::UTILS::StdPairComparatorOrderCounts > bondstatus_;

  //@}

  //! @name methods for initialization of the element

  //@}

  //! @name Internal calculation methods

  //! calculation of thermal (i.e. stochastic) and damping forces according to Brownian dynamics
  void CalcBrownianForcesAndStiff(Teuchos::ParameterList&   params,
                                  std::vector<double>&      vel,
                                  std::vector<double>&      disp,
                                  Epetra_SerialDenseMatrix* stiffmatrix,
                                  Epetra_SerialDenseVector* force);

  //!calculation of drag force and corresponding stiffness contribution
  void CalcDragForce( Teuchos::ParameterList&   params,
                      const std::vector<double>&      vel,
                      const std::vector<double>&      disp,
                      Epetra_SerialDenseMatrix* stiffmatrix,
                      Epetra_SerialDenseVector* force);

  //!calculation of stochastic force and corresponding stiffness contribution
  void CalcStochasticForce( Teuchos::ParameterList&   params,
                       const std::vector<double>&      vel,
                       const std::vector<double>&      disp,
                       Epetra_SerialDenseMatrix* stiffmatrix,
                       Epetra_SerialDenseVector* force);

//  //!calculation of inertia force and corresponding stiffness contribution
//  void CalcInertiaForce( Teuchos::ParameterList&   params,
//                         std::vector<double>&      vel,
//                         std::vector<double>&      disp,
//                         Epetra_SerialDenseMatrix* stiffmatrix,
//                         Epetra_SerialDenseVector* force);

  //!calculation of background fluid velocity and gradient of velocity
  void GetBackgroundVelocity(Teuchos::ParameterList& params,
                            LINALG::Matrix<3,1>& velbackground,
                            LINALG::Matrix<3,3>& velbackgroundgrad);

  //! computes damping coefficient
  double MyDampingConstant();

  //! Compute artificial damping contribution for PTC
  void EvaluatePTC(Teuchos::ParameterList& params,
                     Epetra_SerialDenseMatrix& elemat1);

  // don't want = operator
  Rigidsphere& operator = (const Rigidsphere& old);

}; // class Rigidsphere

// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);

} // namespace Elements
} // namespace DRT

#endif  // #ifndef RIGIDSPHERE_H
