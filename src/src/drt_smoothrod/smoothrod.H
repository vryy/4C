/*!----------------------------------------------------------------------
\file smoothrod.H

\brief three dimensional rod based on smooth curve

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
#ifndef SMOOTHROD_H
#define SMOOTHROD_H


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

class Smoothrod : public DRT::Element
{
public:

  //! @name Friends
  friend class SmoothrodRegister;


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Smoothrod(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Smoothrod(const Smoothrod& old);

  /*!
  \brief Destructor
  */
  virtual ~Smoothrod();

  /*!
  \brief Deep copy this instance of Smoothrod and return pointer to the copy

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
  virtual int UniqueParObjectId() const { return ParObject_Smoothrod; }

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


  //@}


private:

  //! action parameters recognized by Smoothrod
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

  //! quaternions of triads at GPs at the end of the preceeding time step due to Crisfield, Vol. 2, equation (17.65)
  vector<LINALG::Matrix<4,1> > Qconv_;
  //! quaternions of triads at GPs at the end of the preceeding iteration step due to Crisfield, Vol. 2, equation (17.65)
  vector<LINALG::Matrix<4,1> > Qold_;
  //! quaternions of triads at GPs during the current iteration step due to Crisfield, Vol. 2, equation (17.65)
  vector<LINALG::Matrix<4,1> > Qnew_;
  //!curvaturea at the end of the preceeding time step due to Crisfield, Vol. 2, equation (17.72)
  vector<LINALG::Matrix<3,1> > curvconv_;
  //!curvaturea at the end of the preceeding iteration step due to Crisfield, Vol. 2, equation (17.72)
  vector<LINALG::Matrix<3,1> > curvold_;
  //!curvaturea during the current iteration step due to Crisfield, Vol. 2, equation (17.72)
  vector<LINALG::Matrix<3,1> > curvnew_;

  //!angle displacement taken from discretization and interpolated for a specific Gauss point by Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding time step
  vector<LINALG::Matrix<3,1> > thetaconv_;
  //!angle displacement taken from discretization and interpolated for a specific Gauss point by Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding iteration step
  vector<LINALG::Matrix<3,1> > thetaold_;
  //!angle displacement taken from discretization and interpolated for a specific Gauss point by Crisfield, Vol. 2, equation (17.93); relates to current iteration step
  vector<LINALG::Matrix<3,1> > thetanew_;
  //!derivative of angle displacement taken from discretization with respect to curve parameter at a specific Gauss point; derivative due to Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding time step
  vector<LINALG::Matrix<3,1> > thetaprimeconv_;
  //!derivative of angle displacement taken from discretization with respect to curve parameter at a specific Gauss point; derivative due to Crisfield, Vol. 2, equation (17.93); relates to the end of the preceeding iteration step
  vector<LINALG::Matrix<3,1> > thetaprimeold_;
  //!derivative of angle displacement taken from discretization with respect to curve parameter at a specific Gauss point; derivative due to Crisfield, Vol. 2, equation (17.93); relates to the current iteration step
  vector<LINALG::Matrix<3,1> > thetaprimenew_;

  //! cross section area
  double crosssec_;

  //!moment of inertia of area around y-axis(orthognoal to rod axis in yz-plane)
  double Iyy_;

  //!moment of inertia of area around z-axis(forming cartesian system with rod- and y-axis)
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
  //this function is templated and receives nnode from SmoothrodRegister::Initialize
  template<int nnode>
  void SetUpReferenceGeometry(const vector<double>& xrefe,const vector<double>& rotrefe, const bool secondinit = false);

  //@}

  //! @name auxiliary computation methods for non-additive, large rotation variables




  //@}

  //! @name Internal calculation methods

  //!calculation of nonlinear stiffness and mass matrix
  template<int nnode>
  void sr_nlnstiffmass( ParameterList& params,
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
  


  // don't want = operator
  Smoothrod& operator = (const Smoothrod& old);


}; // class Smoothrod



// << operator
ostream& operator << (ostream& os, const DRT::Element& ele);


//------------- class SmoothrodRegister: -------------------------------------------

/*!
\brief A register for Smoothrod element

\author cyron (cyron@lnm.mw.tum.de)
*/
class SmoothrodRegister : public DRT::ElementRegister
{
public:

  //! @name Friends
  friend class Smoothrod;

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  SmoothrodRegister(DRT::Element::ElementType etype);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an SmoothrodRegister

  */
  SmoothrodRegister(const DRT::ELEMENTS::SmoothrodRegister& old);

  /*!
  \brief Deep copy this instance of SmoothrodRegister and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::ELEMENTS::SmoothrodRegister* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_SmoothrodRegister;}

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
  virtual ~SmoothrodRegister();

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

  The method performs all steps needed to init the Smoothrod elements. it e.g. calculates the
  field of averaged directors at nodes.
  */
  int Initialize(DRT::Discretization& dis);

  //@}

private:


}; // class smoothrodRegister
} // namespace Elements
} // namespace DRT


#endif  // #ifndef SMOOTHROD_H
#endif  // #ifdef CCADISCRET
