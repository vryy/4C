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

#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "so3_poro_eletypes.H"


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

/*!
\brief A C++ version of a 3 dimensional solid element with modifications for porous media

A structural 3 dimensional solid displacement element for large deformations
and (near)-incompressibility.

*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
class So3_Poro : public so3_ele
{

  //! @name Friends
  friend class So_tet4PoroType;
  friend class So_hex8PoroType;
  friend class So_hex27PoroType;

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

  /*!
  \brief Destructor

  */
  virtual ~So3_Poro(){};

  //@}

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
  DRT::Element* Clone() const ;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  //virtual int UniqueParObjectId() const = 0 ;
  virtual int UniqueParObjectId() const
  {
    switch(distype)
    {
    case DRT::Element::tet4:
      return So_tet4PoroType::Instance().UniqueParObjectId();
      break;
    case DRT::Element::hex8:
      return So_hex8PoroType::Instance().UniqueParObjectId();
      break;
    default: dserror("unknown element type!");
      break;
    }
    return -1;
  };

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
      // in case of porous media: nds=1 (second discretization) fluid: 4Dof/Node
      return 4;
    }
    // structure, 3D: 3 Dofs per node
    return 3;
  };

  //@}

  //! @name Access methods

  /*!
  \brief Print this element
  */
  void Print(ostream& os) const  ;

  //virtual DRT::ElementType & ElementType() const = 0 ;
  virtual DRT::ElementType & ElementType() const
  {
    switch(distype)
    {
    case DRT::Element::tet4:
      return So_tet4PoroType::Instance();
    case DRT::Element::hex8:
      return So_hex8PoroType::Instance();
    default: dserror("unknown element type!");
      break;
    }
    return So_hex8PoroType::Instance();
  };

  bool ReadElement(const std::string& eletype,
                   const std::string& eledistype,
                   DRT::INPUT::LineDefinition* linedef) ;

  //@}

  //! @name Evaluation

  void PreEvaluate
    (
      ParameterList& params,  ///< ParameterList for communication between control routine and elements
      DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
      DRT::Element::LocationArray& la ///< location array for de-assembly
    );

  //! @name Evaluation

  /*!
  \brief Evaluate an element

  Evaluate So_tet4fbar element stiffness, mass, internal forces, etc.

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


  //! init the inverse of the jacobian and its determinant in the material configuration
  virtual void InitJacobianMapping();

  //@}

protected:

  //! don't want = operator
  So3_Poro& operator = (const So3_Poro& old);

  /*!
  \brief Evaluate an element

  Evaluate So3_poro element stiffness, mass, internal forces, etc.
  Templated evaluate routine of element matrixes

  If NULL on input, the controlling method does not expect the element
  to fill these matrices or vectors.

  \return 0 if successful, negative otherwise
  */
  int MyEvaluate
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
      //LINALG::Matrix<numgpt_,numstr_>* elestress, ///< stresses at GP
      //LINALG::Matrix<numgpt_,numstr_>* elestrain, ///< strains at GP
      ParameterList& params ///< algorithmic parameters e.g. time
  //    const INPAR::STR::StressType iostress ///< stress output option
      );

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
      ParameterList& params); ///< algorithmic parameters e.g. time

  //Gauss integration rule
  DRT::UTILS::GaussIntegration intpoints_;

  //!flag indicating higher order element
  bool ishigherorder_;

private:

  DRT::Node** Nodes();

  RCP<MAT::Material> Material() const;

  int Id() const;
}; // class So3_Poro


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


} // namespace ELEMENTS
} // namespace DRT


#endif /* SO3_PORO_H */
