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
  static const int numdim_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! number of dofs per node
  static const int noddof_ = numdim_;

  //! total dofs per element
  static const int numdof_ = noddof_ * numnod_;

  //! number of strains per node
  static const int numstr_ = (numdim_*(numdim_+1))/2;

  //! number of components necessary to store second derivatives
  /*!
   1 component  for nsd=1:  (N,xx)

   3 components for nsd=2:  (N,xx ; N,yy ; N,xy)

   6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
  */
  static const int numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

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
  virtual int UniqueParObjectId() const;

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

  //! @name Access methods

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
    case DRT::Element::hex27:
      return So_hex27PoroType::Instance();
    default: dserror("unknown element type!");
      break;
    }
    return So_hex8PoroType::Instance();
  };

  //@}

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
    Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
    DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
    DRT::Element::LocationArray& la, ///< location array for de-assembly
    Epetra_SerialDenseMatrix& elemat1,  ///< (stiffness-)matrix to be filled by element.
    Epetra_SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
    Epetra_SerialDenseVector& elevec1,  ///< (internal force-)vector to be filled by element
    Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element
    Epetra_SerialDenseVector& elevec3  ///< vector to be filled by element
  );


  //! initialize the inverse of the jacobian and its determinant in the material configuration
  void InitElement();

  //@}

  //! @name Input and Creation

  /*!
  \brief Read input for this element
  */
  /*!
  \brief Query names of element data to be visualized using BINIO

  The element fills the provided map with key names of
  visualization data the element wants to visualize AT THE CENTER
  of the element geometry. The values is supposed to be dimension of the
  data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
  or 9 (nonsym. tensor)

  Example:
  \code
    // Name of data is 'Owner', dimension is 1 (scalar value)
    names.insert(std::pair<string,int>("Owner",1));
    // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
    names.insert(std::pair<string,int>("StressesXYZ",6));
  \endcode

  \param names (out): On return, the derived class has filled names with
                      key names of data it wants to visualize and with int dimensions
                      of that data.
  */
  virtual void VisNames(std::map<string,int>& names);

  /*!
  \brief Query data to be visualized using BINIO of a given name

  The method is supposed to call this base method to visualize the owner of
  the element.
  If the derived method recognizes a supported data name, it shall fill it
  with corresponding data.
  If it does NOT recognizes the name, it shall do nothing.

  \warning The method must not change size of data

  \param name (in):   Name of data that is currently processed for visualization
  \param data (out):  data to be filled by element if element recognizes the name
  */
  virtual bool VisData(const string& name, std::vector<double>& data);

  /*!
  \brief Read input for this element
  */
  bool ReadElement(const std::string& eletype,
                   const std::string& eledistype,
                   DRT::INPUT::LineDefinition* linedef) ;

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
    Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
    DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
    DRT::Element::LocationArray& la, ///< location array for de-assembly
    Epetra_SerialDenseMatrix& elemat1,  ///< (stiffness-)matrix to be filled by element.
    Epetra_SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
    Epetra_SerialDenseVector& elevec1,  ///< (internal force-)vector to be filled by element
    Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element
    Epetra_SerialDenseVector& elevec3  ///< vector to be filled by element
  );

  void PreEvaluate
    (
      Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
      DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
      DRT::Element::LocationArray& la ///< location array for de-assembly
    );

  //! action parameters recognized by so_hex8
  enum ActionType
  {
    none,
    calc_struct_nlnstiff,
    calc_struct_internalforce,
    calc_struct_nlnstiffmass,
    calc_struct_stress,
    calc_struct_multidofsetcoupling,//!< structure-fluid coupling: internal force, stiffness for poroelasticity (structural part)
    postprocess_stress
  };

  //! container containing all kind of solid3 thisandthat
  DRT::Container                                                        data_;
  //! vector of inverses of the jacobian in material frame
  std::vector<LINALG::Matrix<numdim_,numdim_> >                      invJ_;
  //! determinant of Jacobian in material frame
  std::vector<double>                                                        detJ_;
  //! vector of coordinates of current integration point in reference coordinates
  std::vector<LINALG::Matrix<numdim_,1> > xsi_;

  //! Calculate nonlinear stiffness and internal force for poroelasticity problems
  virtual void nlnstiff_poroelast(
      std::vector<int>& lm, ///< location matrix
      std::vector<double>& disp, ///< current displacements
      std::vector<double>& vel, ///< current velocities
      //    std::vector<double>&      residual,       ///< current residual displ
      LINALG::Matrix<numdim_,numnod_> & evelnp, //< fluid velocity of element
      LINALG::Matrix<numnod_,1> & epreaf, //< fluid pressure of element
      LINALG::Matrix<numdof_,numdof_>* stiffmatrix, ///< element stiffness matrix
      LINALG::Matrix<numdof_,numdof_>* reamatrix, // element reactive matrix
      LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
      Teuchos::ParameterList& params ///< algorithmic parameters e.g. time
      );

  //! Calculate coupling terms in nonlinear stiffness and internal force for poroelasticity problems
  virtual void coupling_poroelast(
      std::vector<int>& lm, ///< location matrix
      std::vector<double>& disp, ///< current displacements
      std::vector<double>& vel, ///< current velocities
      LINALG::Matrix<numdim_,numnod_> & evelnp, //< fluid velocity of element
      LINALG::Matrix<numnod_,1> & epreaf, //< fluid pressure of element
      LINALG::Matrix<numdof_,(numdim_+1)*numnod_>* stiffmatrix, ///< element stiffness matrix
      LINALG::Matrix<numdof_,(numdim_+1)*numnod_>* reamatrix, //< element reactive matrix
      LINALG::Matrix<numdof_,1>* force, ///< element internal force vector
      Teuchos::ParameterList& params); ///< algorithmic parameters e.g. time

  //! Calculate coupling stress for poroelasticity problems
  virtual void couplstress_poroelast(
      std::vector<double>& disp, ///< current displacements
      LINALG::Matrix<numdim_,numnod_> & evelnp, //< fluid velocity of element
      LINALG::Matrix<numnod_,1> & epreaf, //< fluid pressure of element
      Epetra_SerialDenseMatrix* elestress, ///< stresses at GP
      Epetra_SerialDenseMatrix* elestrain, ///< strains at GP
      Teuchos::ParameterList& params, ///< algorithmic parameters e.g. time
      const INPAR::STR::StressType iostress ///< stress output option
      );

private:

  //! Extrapolate Gauss-point values (e.g. stresses) to nodes
  void stress_expol(Epetra_SerialDenseMatrix& stresses,
                   Epetra_MultiVector& expolstresses);

  //! push forward of material stresses to the current, spatial configuration
  void PK2toCauchy(
    LINALG::Matrix<numstr_,1>& stress,
    LINALG::Matrix<numdim_,numdim_>& defgrd,
    LINALG::Matrix<numdim_,numdim_>& cauchystress
    );

  //!Gauss integration rule
  DRT::UTILS::GaussIntegration intpoints_;

  //!flag indicating higher order element
  bool ishigherorder_;

  //!flag indicating if element has been initialized
  bool init_;

  //! flag for scatra coupling
  bool scatracoupling_;

  //! get nodes of element
  DRT::Node** Nodes();

  //! get material of element
  Teuchos::RCP<MAT::Material> Material() const;

  //! get global id of element
  int Id() const;
}; // class So3_Poro


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


} // namespace ELEMENTS
} // namespace DRT


#endif /* SO3_PORO_H */
