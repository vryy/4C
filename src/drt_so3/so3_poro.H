/*!----------------------------------------------------------------------
\file so3_poro.H

<pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifndef SO3_PORO_H
#define SO3_PORO_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

/*!
\brief A C++ version of a 8-node hex solid element with modifications for porous media

A structural 8-node hexahedral solid displacement element for large deformations
and (near)-incompressibility.

*/
template<DRT::Element::DiscretizationType distype>
class So3_Poro : virtual public DRT::Element
{

protected:
  //So3_Poro();
public:

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owner
  */
  So3_Poro(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So3_Poro(const So3_Poro& old);

  //@}

  //static const DRT::UTILS::GaussIntegration intpoints_ = DRT::UTILS::GaussIntegration(distype);

  //! number of element nodes (
  static const int numnod_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions
  static const int numdim_ = 3;

  //! number of dofs per node
  static const int noddof_ = 3;

  //! total dofs per element
  static const int numdof_ = noddof_ * numnod_;

  //! number of dofs per node
  static const int numstr_ = 6;

  //!total gauss points per element
  int numgpt_;

  //! @name Acess methods

  /*!
  \brief Deep copy this instance of Solid3 and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const = 0 ;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const = 0 ;

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const vector<char>& data);

  /*!
  \brief Destructor

  */
  virtual ~So3_Poro(){};

  /*!
  \brief Get number of degrees of freedom of a certain node in case of multiple
         dofsets (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.
  */
  virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node) const
  {
    if (nds==1)
    {
      // in case of Poroelasticity: nds=1 (second discretization) fluid: 4Dof/Node
      return 4;
    }
    // structure, 3D: 3 Dofs per node
    //    return 3;
    return 3;
  };

  //@}

  //! @name Access methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const = 0 ;

  virtual DRT::ElementType & ElementType() const = 0 ;
 // { return So3_PoroType::Instance();} //actEle_->ElementType(); }

  //RCP<DRT::ELEMENTS::So3> ActEle(){return actEle_;};

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */

  //virtual bool ReadElement();

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& distype,
                   DRT::INPUT::LineDefinition* linedef) = 0 ;

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate So_hex8fbar element stiffness, mass, internal forces, etc.

  If NULL on input, the controlling method does not expect the element
  to fill these matrices or vectors.

  \return 0 if successful, negative otherwise
  */
  int Evaluate
  (
    ParameterList& params,  ///< ParameterList for communication between control routine and elements
    DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
    DRT::Element::LocationArray& la, ///< location array for de-assembly
    Epetra_SerialDenseMatrix& elemat1,  ///< (stiffness-)matrix to be filled by element.
    Epetra_SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
    Epetra_SerialDenseVector& elevec1,  ///< (internal force-)vector to be filled by element
    Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element
    Epetra_SerialDenseVector& elevec3  ///< vector to be filled by element
  );

  /*!
  \brief Evaluate a Neumann boundary condition

  this method evaluates a surface Neumann condition on the solid3 element

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param condition (in)     : The condition to be evaluated
  \param lm (in)            : location vector of this element
  \param elevec1 (out)      : vector to be filled by element. If NULL on input,

  \return 0 if successful, negative otherwise

  virtual int EvaluateNeumann
  (
    ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1 = NULL
  )=0;*/

  //! init the inverse of the jacobian and its determinant in the material configuration
  virtual void InitJacobianMapping();

  //@}

protected:

  //! don't want = operator
  So3_Poro& operator = (const So3_Poro& old);

  //! action parameters recognized by so_hex8
  enum ActionType
  {
    none,
    calc_struct_linstiff,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_linstiffmass,
    calc_struct_nlnstiffmass,
    calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
    calc_struct_stress,
    calc_struct_update_istep,
    calc_struct_reset_istep,  //!< reset elementwise internal variables
                              //!< during iteration to last converged state
    calc_struct_reset_discretization,  //!< reset elementwise internal variables
                                       //!< to state in the beginning of the computation
    calc_struct_energy,       //!< compute internal energy
    calc_struct_errornorms,   //!< compute error norms (L2,H1,energy)
    postprocess_stress,       //!< postprocess stresses
    prestress_update,
    calc_struct_multidofsetcoupling,//!< structure-fluid coupling: internal force, stiffness for poroelasticity (structural part)
  };

  //! vector of inverses of the jacobian in material frame
  vector<LINALG::Matrix<numdim_,numdim_> >                      invJ_;
  //! determinant of Jacobian in material frame
  vector<double>                                                        detJ_;
  //! container containing all kind of solid3 thisandthat
  DRT::Container                                                        data_;
  //! vector of coordinates of current integration point in reference coordinates
  vector<LINALG::Matrix<numdim_,1> > xsi_;

  //! Calculate nonlinear stiffness and internal force for poroelasticity problems
  virtual void nlnstiff_poroelast(
      vector<int>& lm, ///< location matrix
      vector<double>& disp, ///< current displacements
      vector<double>& vel, ///< current velocities
      //    vector<double>&           residual,       ///< current residual displ
      LINALG::Matrix<numdim_,numnod_> & evelnp, //< fluid velocity of element
      LINALG::Matrix<numnod_,1> & epreaf, //< fluid pressure of element
      LINALG::Matrix<numdof_,numdof_>* stiffmatrix, ///< element stiffness matrix
      LINALG::Matrix<numdof_,numdof_>* reamatrix, // element reactive matrix
      LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
      LINALG::Matrix<numdof_,1>* forcerea, // element reactive force vector
      //LINALG::Matrix<numgpt_,numstr_>* elestress, ///< stresses at GP
      //LINALG::Matrix<numgpt_,numstr_>* elestrain, ///< strains at GP
      ParameterList& params, ///< algorithmic parameters e.g. time
      const INPAR::STR::StressType iostress, ///< stress output option
      const INPAR::STR::StrainType iostrain); ///< strain output option

  //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity problems
  virtual void coupling_poroelast(
      vector<int>& lm, ///< location matrix
      vector<double>& disp, ///< current displacements
      vector<double>& vel, ///< current velocities
      LINALG::Matrix<numdim_,numnod_> & evelnp, //< fluid velocity of element
      LINALG::Matrix<numnod_,1> & epreaf, //< fluid pressure of element
      LINALG::Matrix<numdof_,(numdim_+1)*numnod_>* stiffmatrix, ///< element stiffness matrix
      LINALG::Matrix<numdof_,(numdim_+1)*numnod_>* reamatrix, //< element reactive matrix
      LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
      LINALG::Matrix<numdof_,1>* forcerea, //< element reactive force vector
      ParameterList& params); ///< algorithmic parameters e.g. time

private:

  DRT::UTILS::GaussIntegration intpoints_;
}; // class So3_Poro


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

} // namespace ELEMENTS
} // namespace DRT

#endif /* SO3_PORO_H */
