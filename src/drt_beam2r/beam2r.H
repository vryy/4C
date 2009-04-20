/*!----------------------------------------------------------------------
\file beam2.H

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
// line in order to swith on and off by comment in defines.drt
#ifdef D_BEAM2R

// administrative line in order to make a difference between BACI and Ccarat files
#ifdef CCADISCRET

//header file only included if not yet included!
#ifndef BEAM2R_H
#define BEAM2R_H


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
#include "../drt_lib/linalg_fixedsizematrix.H"
//including random number library of blitz for statistical forces
#include <random/normal.h>




namespace DRT
{

namespace ELEMENTS
{
/*!
\brief two dimensional nonlinear corotational Timoshenko beam element according to Crisfield, Vol.1, section 7.3 

\author cyron (cyron @lnm.mw.tum.de)
*/
class Beam2r : public DRT::Element
{
public:

  //! @name Friends
  friend class Beam2rRegister;


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
  virtual int UniqueParObjectId() const { return ParObject_Beam2r; }

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
  \brief Evaluate damping matrix for Brownian motion

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates a damping matrix
  as used in the problem type STATISTICAL MECHANICS

  \param params (in/out)       : ParameterList for communication between control routine and elements
  \param vector<double> mydisp : current nodal displacement
  \param elemat1 (out)         : damping matrix to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateBrownianDamp(ParameterList& params,
                                    vector<double> mydisp,
                                    double zeta,
                                    Epetra_SerialDenseMatrix& elemat1);  
  
  /*!
  \brief Evaluate forces for Brownian motion

  An element derived from this class uses the Evaluate method to receive commands
  and parameters from some control routine in params and evaluates a  statistical Neumann boundary condition
  as used in the problem type STATISTICAL MECHANICS

  \param params (in/out)       : ParameterList for communication between control routine and elements
  \param vector<double> mydisp : current nodal displacement
  \param elevec1 (out)         : Force vector to be filled by element

  \return 0 if successful, negative otherwise
  */
  virtual int EvaluateBrownianForces(ParameterList& params,
                                    vector<double> mydisp,
                                    double zeta,
                                    Epetra_SerialDenseVector& elevec1);
  
  
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
    calc_struct_update_imrlike,
    calc_struct_reset_istep, //!< reset elementwise internal variables during iteration to last converged state
    calc_brownian_damp,
    calc_struct_ptcstiff
  };

  //! container containing all kind of beam2r thisandthat
  DRT::Container data_;
  
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

  //! gaussrule_ will be initialized automatically to a 2 point integration rule
  DRT::UTILS::GaussRule1D gaussrule_;

  // internal calculation methods

  //!internal geometric calculations
  void local_aux(LINALG::Matrix<3,6>& B_curr,
                    LINALG::Matrix<6,1>& r_curr,
                    LINALG::Matrix<6,1>& z_curr,
                    const double& lcurr_,
                    const double& lrefe_);
  
  //!computes current rotation angle beta relative to x-axis
  void updatealpha(const LINALG::Matrix<3,2>& xcurr,
                   const double& lcurr);
  
  //! calculation of nonlinear stiffness and mass matrix
  void nlnstiffmass(ParameterList& params,
                       vector<double>&           vel,  //!< element velocity vector
                       vector<double>&           disp,  //!< element displacement vector
                       Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                       Epetra_SerialDenseMatrix* massmatrix,  //!< element mass matrix
                       Epetra_SerialDenseVector* force,  //!< element internal force vector
                       int lumpedmass); //!< if set to one lumped mass matrix is calculated
  
  //!sets up element reference geomtry for reference nodal position vector xrefe (may be used also after simulation start)
  void SetUpReferenceGeometry(const LINALG::Matrix<4,1>& xrefe);
  
  
  // don't want = operator
  Beam2r& operator = (const Beam2r& old);

}; // class Beam2r


// << operator
ostream& operator << (ostream& os, const DRT::Element& ele);




//------------- class Beam2rRegister: -------------------------------------------

/*!
\brief A register for Beam2r element

\author cyron (cyron@lnm.mw.tum.de)
*/
class Beam2rRegister : public DRT::ElementRegister
{
public:

  //! @name Friends
  friend class Beam2r;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  Beam2rRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Beam2rRegister

  */
  Beam2rRegister(const DRT::ELEMENTS::Beam2rRegister& old);

  /*!
  \brief Deep copy this instance of Beam2rRegister and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::Beam2rRegister* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_Beam2rRegister;}

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
  virtual ~Beam2rRegister();

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

  The method performs all steps needed to init the Beam2r elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class Beam2rRegister
} // namespace Elements
} // namespace DRT




#endif  // #ifndef BEAM_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2R

