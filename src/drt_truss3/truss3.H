/*!----------------------------------------------------------------------
\file truss3.H

\brief three dimensional total Lagrange truss element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
// line in order to swith on and off by comment in defines.drt
#ifdef D_TRUSS3

// administrative line in order to make a difference between BACI and Ccarat files
#ifdef CCADISCRET

//header file only included if not yet included!
#ifndef TRUSS_H
#define TRUSS_H


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
#include "../drt_lib/linalg_fixedsizematrix.H"


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

namespace ELEMENTS
{
/*!
\brief three dimensional total Lagrange truss element

\author cyron (cyron @lnm.mw.tum.de)
*/
class Truss3 : public DRT::Element
{
public:

  //! @name Friends
  friend class Truss3Register;
  friend class StatMechManager;
  friend void StatMechManager::StatMechUpdate(const double dt, const Epetra_Vector& dis);
  friend void StatMechTime::Integrate();


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Truss3(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Truss3(const Truss3& old);

  /*!
  \brief Destructor
  */
  virtual ~Truss3();

  /*!
  \brief Deep copy this instance of Truss3 and return pointer to the copy
    
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
  virtual int UniqueParObjectId() const { return ParObject_Truss3; }

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
    return 3;
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
  virtual bool ReadElement();


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
                              Epetra_SerialDenseVector& elevec1);    
  /*!
  \brief Evaluate a statistical Neumann boundary condition and related damping matrix

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates a  statistical Neumann boundary condition
  as used in the problem type STATISTICAL MECHANICS

  \param params (in/out)       : ParameterList for communication between control routine and elements
  \param vector<double> mydisp : current nodal displacement
  \param elevec1 (out)         : Force vector to be filled by element
  \param elemat1 (out)         : damping matrix to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateStatForceDamp(ParameterList& params,
                                    vector<double> mydisp,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    int& ActNumDof0);  
  
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
  virtual int EvaluatePTC(ParameterList& params,
                          Epetra_SerialDenseMatrix& elemat1,
                          int& ActNumDof0);
  //@}

  
private:

  //! action parameters recognized by Truss3
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
    calc_struct_update_imrlike,
    calc_struct_reset_istep,  //!< reset elementwise internal variables
                             //!< during iteration to last converged state
    postprocess_stress,
    calc_stat_force_damp, //calculate statistical forces and damping for Brownian dynamics
    calc_struct_ptcstiff //calculate artificial stiffness due to PTC solution strategy
  };
  
  //! possible kinematic types
  enum KinematicType
  {
    tr3_totlag,
    tr3_engstrain
  };
  
  //! container containing all kind of Truss3 thisandthat
  DRT::Container data_; 
  
  //!reference node position
  LINALG::Matrix<6,1> X_;
  
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
  void SetUpReferenceGeometry(const LINALG::Matrix<6,1>& xrefe);
  
  //@}

  //! @name Internal calculation methods
              
  //!calculation of nonlinear stiffness and mass matrix switching between total lagrange and enginerring strains
  void t3_nlnstiffmass(vector<double>&           disp, 
                       Epetra_SerialDenseMatrix* stiffmatrix,
                       Epetra_SerialDenseMatrix* massmatrix,
                       Epetra_SerialDenseVector* force,
                       int& ActNumDof0);
  
  
  //!calculation of nonlinear stiffness and mass matrix
  void t3_nlnstiffmass_totlag( vector<double>& disp,
      Epetra_SerialDenseMatrix* stiffmatrix,
      Epetra_SerialDenseMatrix* massmatrix,
      Epetra_SerialDenseVector* force,
      int& ActNumDof0);
  
  /*! 
   * \brief calculate force, nonlinear stiffness and mass matrix using a engineering strain measure.
   * 
   * Unlike the fully nonlinear implementation of #t3_nlnstiffmass_totlag, this function uses \f$\varepsilon=\Delta d / d\f$
   * as strain measure. 
   *  
  */ 
  void t3_nlnstiffmass_engstr( vector<double>& disp,
      Epetra_SerialDenseMatrix* stiffmatrix,
      Epetra_SerialDenseMatrix* massmatrix,
      Epetra_SerialDenseVector* force,
      int& ActNumDof0);
  
  //! lump mass matrix
  void t3_lumpmass(Epetra_SerialDenseMatrix* emass);
  
  //@}
			    
  // don't want = operator
  Truss3& operator = (const Truss3& old);
  

}; // class Truss3



// << operator
ostream& operator << (ostream& os, const DRT::Element& ele);


//------------- class Truss3Register: -------------------------------------------

/*!
\brief A register for Truss3 element

\author cyron (cyron@lnm.mw.tum.de)
*/
class Truss3Register : public DRT::ElementRegister
{
public:
	
  //! @name Friends
  friend class Truss3;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor
  
  */
  Truss3Register(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Truss3Register
  
  */
  Truss3Register(const DRT::ELEMENTS::Truss3Register& old);

  /*!
  \brief Deep copy this instance of Truss3Register and return pointer to the copy
  
  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed
  
  */
  DRT::ELEMENTS::Truss3Register* Clone() const;

  /*!
  \brief Return unique ParObject id
  
  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Truss3Register;}

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
  virtual ~Truss3Register();

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
  
  The method performs all steps needed to init the Truss3 elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Truss3Register
} // namespace Elements
} // namespace DRT




#endif  // #ifndef TRUSS_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_TRUSS3
