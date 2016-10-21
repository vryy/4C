/*!----------------------------------------------------------------------
\file truss2.H

\brief three dimensional total Lagrange truss element


\maintainer Christian Cyron

\level 2
*----------------------------------------------------------------------*/
// administrative line in order to make a difference between BACI and Ccarat files

//header file only included if not yet included!
#ifndef TRUSS2_H
#define TRUSS2_H


#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

namespace ELEMENTS
{

class Truss2Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Truss2Type"; }

  static Truss2Type & Instance();

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

  static Truss2Type instance_;
};

/*!
\brief three dimensional total Lagrange truss element

\author cyron (cyron @lnm.mw.tum.de)
*/
class Truss2 : public DRT::Element
{
public:

  //! @name Friends
  friend class Truss2Type;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Truss2(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Truss2(const Truss2& old);

  /*!
  \brief Destructor
  */
  virtual ~Truss2();

  /*!
  \brief Deep copy this instance of Truss2 and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return Truss2Type::Instance().UniqueParObjectId(); }

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

  virtual DRT::ElementType & ElementType() const { return Truss2Type::Instance(); }

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
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    /*note: this is not necessarily the number of DOF assigned to this node by the discretization finally, but only
     *the number of DOF requested for this node by this element; the discretization will finally assign the maximal
     *number of DOF to this node requested by any element connected to this node*/
    return 2;
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
  \brief Evaluate PTC addition to stiffness for free Brownian motion

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates a  statistical Neumann boundary condition
  as used in the problem type STATISTICAL MECHANICS

  \param params (in/out)       : ParameterList for communication between control routine and elements
  \param vector<double> mydisp : current nodal displacement
  \param elemat1 (out)         : artificial damping matrix to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluatePTC(Teuchos::ParameterList&   params,
                          Epetra_SerialDenseMatrix& elemat1,
                          int&                      ActNumDof0);
  //@}


private:

  //! action parameters recognized by Truss2
  enum ActionType
  {
    calc_none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_nlnstiff2,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_reset_istep,  //!< reset elementwise internal variables
                             //!< during iteration to last converged state
    postprocess_stress,
    calc_struct_ptcstiff //calculate artificial stiffness due to PTC solution strategy
  };

  //! possible kinematic types
  enum KinematicType
  {
    tr2_totlag,
    tr2_engstrain
  };

  //! container containing all kind of Truss2 thisandthat
  DRT::Container data_;

  //! variable saving whether element has already been initialized (then isinit_ == true)
  bool isinit_;

  //!reference node position
  LINALG::Matrix<4,1> X_;

  //! material type
  int material_;

  //! length in reference configuration
  double lrefe_;

  //! cross section area
  double crosssec_;

  //! Kinematic type
  KinematicType kintype_;

  //! gaussrule_ will be initialized automatically to a 2 point integration rule
  DRT::UTILS::GaussRule1D gaussrule_;

  //! @name methods for initilization of the element

  //!sets reference coordinates X_ and refenrence length lrefe_ for elements added to the discretization after simulation start
  void SetUpReferenceGeometry(const LINALG::Matrix<4,1>& xrefe);

  //@}

  //! @name Internal calculation methods

  //!calculation of nonlinear stiffness and mass matrix switching between total lagrange and enginerring strains
  void t2_nlnstiffmass(std::vector<double>&      disp,
                       Epetra_SerialDenseMatrix* stiffmatrix,
                       Epetra_SerialDenseMatrix* massmatrix,
                       Epetra_SerialDenseVector* force);


  //!calculation of nonlinear stiffness and mass matrix
  void t2_nlnstiffmass_totlag(std::vector<double>&      disp,
                              Epetra_SerialDenseMatrix* stiffmatrix,
                              Epetra_SerialDenseMatrix* massmatrix,
                              Epetra_SerialDenseVector* force);

  /*!
   * \brief calculate force, nonlinear stiffness and mass matrix using a engineering strain measure.
   *
   * Unlike the fully nonlinear implementation of #t2_nlnstiffmass_totlag, this function uses \f$\varepsilon=\Delta d / d\f$
   * as strain measure.
   *
  */
  void t2_nlnstiffmass_engstr(std::vector<double>&      disp,
                              Epetra_SerialDenseMatrix* stiffmatrix,
                              Epetra_SerialDenseMatrix* massmatrix,
                              Epetra_SerialDenseVector* force);

  //! lump mass matrix
  void t2_lumpmass(Epetra_SerialDenseMatrix* emass);

  //@}

  // don't want = operator
  Truss2& operator = (const Truss2& old);


}; // class Truss2



// << operator
std::ostream& operator << (std::ostream& os, const DRT::Element& ele);


} // namespace Elements
} // namespace DRT




#endif  // #ifndef TRUSS2_H
