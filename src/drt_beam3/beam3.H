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
#include "../drt_lib/drt_utils_integration.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_utils.H"
//including random number library of blitz for statistical forces
#include <random/normal.h>


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
\brief three dimensional nonlinear Reissner beam element

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
  \brief Get vector of ptrs to the lines of this element
  */
  virtual Element** Lines();


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
  
  void Arbeit(double& AN,double& AM,double& AQ, double& xv);
  void Thermik(double& kT);
                              
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
    calc_struct_stress,
    calc_struct_eleload,
    calc_struct_fsiload,
    calc_struct_update_istep,
    calc_struct_update_genalpha_imrlike
  };
  
  //! container containing all kind of Beam3 thisandthat
  DRT::Container data_;  
  
  //! material type
  int material_;
  
  //! cross section area
  double crosssec_;
  
   //! cross section area including shear corretion factor
  double crosssecshear_;
  
  //! moments of inertia of area
  double Iyy_;
  double Izz_;
  double Iyz_;
  double Irr_;
  
  //! flag determining if consistent or lumped mass matrix
  int lumpedflag_;
  
  //! thermal energy responsible for statistical forces
  double thermalenergy_;

  double Arbeit_N;
  double Arbeit_M;
  double Arbeit_Q;
  double x_verschiebung;
  
   //! vector of lines of this element (length 1)
  vector<DRT::Element*> lines_; 
 
  //! gaussrule_ will be initialized automatically to a 2 point integration rule 
  DRT::UTILS::GaussRule1D gaussrule_;


  
  // internal calculation methods
  
  
  //internal geometric calculations
  void b3_local_aux(LINALG::SerialDenseMatrix& B_curr,
                    LINALG::SerialDenseVector& r_curr,
                    LINALG::SerialDenseVector& z_curr,
                    double& beta,
                    const LINALG::SerialDenseMatrix& x_curr,
                    const double& length_curr,
                    const double& length_refe);

  //calculation of nonlinear stiffness and mass matrix
  void b3_nlnstiffmass(vector<double>&           disp, 
                       Epetra_SerialDenseMatrix& stiffmatrix,
                       Epetra_SerialDenseMatrix& massmatrix,
                       Epetra_SerialDenseVector& force);
  
  //calculation of nonlinear stiffness approximated by finite differences; this function should
  //be activated only for bug fixing
  
  Epetra_SerialDenseMatrix b3_nlnstiff_approx(vector<double>& disp, double h_rel);
  
			    
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
