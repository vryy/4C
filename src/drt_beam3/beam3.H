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
//including random number library of blitz for statistical forces
#include <random/normal.h>
#include <blitz/array.h>
#include "../drt_geometry/blitz_tiny_operation.H"

//! generally usefull typedefs, when working with Blitz arrays
typedef blitz::TinyVector<double,3> BlitzVec3; ///< a short-cut for a one-dimensional BLITZ vector of length 3
typedef blitz::TinyVector<double,6> BlitzVec6; ///< a short-cut for a one-dimensional BLITZ vector of length 6
typedef blitz::TinyMatrix<double,3,2> BlitzMat3x2; ///< a short-cut for a two-dimensional BLITZ matrix of size 3x2
typedef blitz::TinyMatrix<double,3,3> BlitzMat3x3; ///< a short-cut for a two-dimensional BLITZ matrix of size 3x3
typedef blitz::TinyMatrix<double,6,6> BlitzMat6x6; ///< a short-cut for a two-dimensional BLITZ matrix of size 6x6
typedef blitz::TinyMatrix<double,6,2> BlitzMat6x2; ///< a short-cut for a two-dimensional BLITZ matrix of size 6x2


// forward declarations
struct _MATERIAL;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

// forward declarations
class Discretization;
class Beam3Register;

namespace ELEMENTS
{
/*!
\brief three dimensional nonlinear corotational Timoshenko beam element according to Crisfield, Vol.2,
section 17.2 ("an interpretation of an element due to simo and vu-quoc"), p. 226 - 231; the current
implementation sets the initial director angles and curvature to zero, which can be done without loss
of generality since any initial rotation cancels out in the equations (17.72) and (17.76). As a consequence
the current nodal directors do not appear in the element, but only the values of the current director
angles minus reference director angles, respectively 

\author cyron (cyron @lnm.mw.tum.de)
*/
class Beam3 : public DRT::Element
{
public:

  //! @name Friends
  friend class Beam3Register;


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
  virtual int NumDofPerNode(const DRT::Node& node) const { return 6;}

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
    calc_struct_reset_istep  //!< reset elementwise internal variables
                             //!< during iteration to last converged state
  };
  
  //! container containing all kind of Beam3 thisandthat
  DRT::Container data_;  
  
  //! material type
  int material_;
  
  //! length in reference configuration
  double lrefe_;
    
  //!central triad at the end of hte preceeding time step due to Crisfield, Vol. 2, equation (17.65)
  BlitzMat3x3 Told_;
  //! central triad during the current iteration step due to Crisfield, Vol. 2, equation (17.65)
  BlitzMat3x3 Tnew_;
  //!midpoint triad for curvature update due to Crisfield, Vol. 2, equation (17.73)
  BlitzMat3x3 Tmid_;
  //!curvature at the end of hte preceeding time step due to Crisfield, Vol. 2, equation (17.72)
  BlitzVec3 curvold_;
  //!curvature during the current iteration step due to Crisfield, Vol. 2, equation (17.72)
  BlitzVec3 curvnew_;
  
  //!auxiliary angle due to Crisfield, Vol. 2, equation (17.70) at the end of hte preceeding time step
  BlitzVec3 betaplusalphaold_;
  //!auxiliary angle due to Crisfield, Vol. 2, equation (17.70) during current iteration step
  BlitzVec3 betaplusalphanew_;
  //!auxiliary angle due to Crisfield, Vol. 2, equation (17.71) at the end of hte preceeding time step
  BlitzVec3 betaminusalphaold_;
  //!auxiliary angle due to Crisfield, Vol. 2, equation (17.71) during current iteration step
  BlitzVec3 betaminusalphanew_;
  
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

  //! thermal energy responsible for statistical forces
  double thermalenergy_;

  //! gaussrule_ will be initialized automatically to a 2 point integration rule 
  DRT::UTILS::GaussRule1D gaussrule_;
  

  //! @name Internal calculation methods

  //!computing spin matrix out of a rotation vector
  void computespin(BlitzMat3x3& spin, BlitzVec3 rotationangle, const double& spinscale);
  
  //!computing rotation matrix out of a rotation vector theta by Crisfield, Vol. 2, equation (16.34) and (16.35)
  void computerotation(BlitzMat3x3& rotationmatrix, BlitzVec3 theta, const double& scale);
  
  //updating curvature
  void updatecurvature(BlitzVec3 deltabetaplusalpha,BlitzVec3 deltabetaminusalpha);
  
  //!computing stiffness matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.83)
  void computeKsig1(Epetra_SerialDenseMatrix& Ksig1, const BlitzVec3& stressn, const BlitzVec3& stressm);
  
  //!computing stiffness matrix Ksigma2 according to Crisfield, Vol. 2, equation (17.87)
  void computeKsig2(Epetra_SerialDenseMatrix& Ksig1, const BlitzVec3& stressn, const BlitzVec3& x21);
  
  //!computing basis of stiffness matrix according to Crisfield, Vol. 2, equation (17.81)
  void computestiffbasis(const BlitzMat3x3& Tnew_, const BlitzVec3& Cm, const BlitzVec3& Cb, const BlitzMat3x3& spinx21, Epetra_SerialDenseMatrix& stiffmatrix);
                    
  //!calculation of nonlinear stiffness and mass matrix
  void b3_nlnstiffmass(vector<double>&           disp, 
                       Epetra_SerialDenseMatrix* stiffmatrix,
                       Epetra_SerialDenseMatrix* massmatrix,
                       Epetra_SerialDenseVector* force);
  
  //calculation of nonlinear stiffness approximated by finite differences; this function should
  //be activated only for bug fixing 
  /*
  Epetra_SerialDenseMatrix b3_nlnstiff_approx(BlitzMat6x2 xcurr, double h_rel, Epetra_SerialDenseVector force);
  */
  
  
  //! lump mass matrix (bborn 07/08)
  void b3_lumpmass(Epetra_SerialDenseMatrix* emass);
  
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
