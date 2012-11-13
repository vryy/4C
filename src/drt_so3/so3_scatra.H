/*!----------------------------------------------------------------------
\file So3_scatra.H

<pre>
   Maintainer: Cristobal Bertoglio
               bertoglio@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifndef SO3_SCATRA_H
#define SO3_SCATRA_H

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_inpar/inpar_structure.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "so3_scatra_eletypes.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

/*!
\brief A C++ version of a 3 dimensional solid element with modifications for volume coupling with scatra

A structural 3 dimensional solid displacement element for large deformations
and (near)-incompressibility.

*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
class So3_Scatra : public so3_ele
{

  //! @name Friends
  friend class So_tet4ScatraType;
  friend class So_hex8ScatraType;
  friend class So_weg6ScatraType;

public:

  //@}
  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner : elements owner
  */
  So3_Scatra(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So3_Scatra(const So3_Scatra& old);

  /*!
  \brief Destructor

  */
  virtual ~So3_Scatra(){};

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
  DRT::Element* Clone() const;

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
      return So_tet4ScatraType::Instance().UniqueParObjectId();
    case DRT::Element::hex8:
      return So_hex8ScatraType::Instance().UniqueParObjectId();
    case DRT::Element::wedge6:
      return So_weg6ScatraType::Instance().UniqueParObjectId();
    default: dserror("unknown element type!");
      break;
    }
    // Intel compiler needs a return so
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
  virtual void Unpack(const std::vector<char>& data);

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
      // in case of Scatra coupling: nds=1 (second discretization): 1Dof/Node
      return 1; // TO CHECK (from Cristobal)
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
  void Print(ostream& os) const ;

  virtual DRT::ElementType & ElementType() const
  {
    switch(distype)
    {
    case DRT::Element::tet4:
      return So_tet4ScatraType::Instance();
    case DRT::Element::hex8:
      return So_hex8ScatraType::Instance();
    case DRT::Element::wedge6:
      return So_weg6ScatraType::Instance();
    default:
      dserror("unknown element type!");
      break;
    }
    // Intel compiler needs a return so
    return So_weg6ScatraType::Instance();
  };


  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& eledistype,
                   DRT::INPUT::LineDefinition* linedef)  ;

  //@}

  //! @name Evaluation

  void PreEvaluate
    (
      ParameterList& params,  ///< ParameterList for communication between control routine and elements
      DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
      DRT::Element::LocationArray& la ///< location array for de-assembly
    );

  /*!
  \brief Evaluate an element

  Evaluate So3_poro element stiffness, mass, internal forces, etc.
  Templated evaluate routine of element matrixes

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
   // virtual void InitJacobianMapping();


  //@}

protected:

  //! don't want = operator
  So3_Scatra& operator = (const So3_Scatra& old);

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
  //vector<LINALG::Matrix<numdim_,numdim_> >                      invJ_;
  //! determinant of Jacobian in material frame
  //vector<double>                                                        detJ_;
  //! container containing all kind of solid3 thisandthat
  //DRT::Container                                                        data_;
  //! vector of coordinates of current integration point in reference coordinates
  //vector<LINALG::Matrix<numdim_,1> > xsi_;


private:

  DRT::UTILS::GaussIntegration intpoints_;

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

#endif /* SO3_SCATRA_H */
