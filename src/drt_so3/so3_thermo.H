/*!----------------------------------------------------------------------
\file so3_poro.H

<pre>
   Maintainer: Caroline Danowski
               danowski@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15253
</pre>

*----------------------------------------------------------------------*/
#ifndef SO3_THERMO_H
#define SO3_THERMO_H

#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

#include "../drt_fem_general/drt_utils_shapefunctions_service.H"
#include "so3_thermo_eletypes.H"


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

  //! A C++ version of a 3 dimensional solid element with modifications for
  //! thermomechanics
  //!
  //! A structural 3 dimensional solid displacement element for large deformations
  //! and so far with small displacments
  template<class so3_ele, DRT::Element::DiscretizationType distype>
  class So3_Thermo : public so3_ele
  {

    //! @name Friends
    friend class So_tet4ThermoType;
    friend class So_hex8ThermoType;

  public:

    //@}
    //! @name Constructors and destructors and related methods


    //! Standard Constructor
    So3_Thermo(
      int id,  //!< (i) this element's global id
      int owner  //!< elements owner
      );

    //! Copy Constructor
    //! Makes a deep copy of a Element
    So3_Thermo(const So3_Thermo& old);

    //! Destructor
    virtual ~So3_Thermo(){};

    //@}

    //! number of element nodes
    static const int numnod_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    //! number of space dimensions
    static const int numdim_ = 3;

    //! number of dofs per node
    static const int noddof_ = 3;

    //! total dofs per element
    static const int numdof_ = noddof_ * numnod_;

    // TODO 2012-07-26 check if thermo has also 6 dofs???
    //! number of dofs per node
    static const int numstr_ = 6;

    //!total gauss points per element
    int numgpt_;

    //! @name Acess methods

    //! Deep copy this instance of Solid3 and return pointer to the copy
    //!
    //! The Clone() method is used from the virtual base class Element in cases
    //! where the type of the derived class is unknown and a copy-ctor is needed
    DRT::Element* Clone() const ;

    //! Return unique ParObject id
    //!
    //! every class implementing ParObject needs a unique id defined at the top of
    //! this file.
//    virtual int UniqueParObjectId() const = 0 ;
    virtual int UniqueParObjectId() const
    {
      switch(distype)
      {
      case DRT::Element::tet4:
        return So_tet4ThermoType::Instance().UniqueParObjectId();
      case DRT::Element::hex8:
        return So_hex8ThermoType::Instance().UniqueParObjectId();
      default:
        dserror("unknown element type!");
      }
    };

    //! Pack this class so it can be communicated
    //! Pack and \ref Unpack are used to communicate this element
    virtual void Pack(DRT::PackBuffer& data) const;

    //! Unpack data from a char vector into this class
    //!Pack and \ref Unpack are used to communicate this element
    virtual void Unpack(const vector<char>& data);

    //! Get number of degrees of freedom of a certain node in case of multiple
    //! dofsets (implements pure virtual DRT::Element)
    //!
    //! The element decides how many degrees of freedom its nodes must have.
    //! As this may vary along a simulation, the element can redecide the number
    //! of degrees of freedom per node along the way for each of it's nodes separately.
    virtual int NumDofPerNode(
      const unsigned nds,
      const DRT::Node& node
      ) const
    {
      if (nds==1)
      {
        // in case of thermomechanics: nds=1 (second discretisation) temperature:
        // 1Dof/Node
        return 1;
      }
      // structure, 3D: 3 Dofs per node
      return 3;
    };

    //@}

    //! @name Access methods

    //! Print this element
    void Print(ostream& os) const  ;

    //virtual DRT::ElementType & ElementType() const = 0 ;
    virtual DRT::ElementType & ElementType() const
    {
      switch(distype)
      {
      case DRT::Element::tet4:
        return So_tet4ThermoType::Instance();

      case DRT::Element::hex8:
        return So_hex8ThermoType::Instance();

      default:
        dserror("unknown element type!");
      }
    };

    bool ReadElement(
      const std::string& eletype,
      const std::string& eledistype,
      DRT::INPUT::LineDefinition* linedef
      );

    //@}

    // TODO 2012-07-26 not sure if required for thermo, PreEvaluate undertakes the task to calculate coupling term of matrix
    //! evaluation
    void PreEvaluate(
      ParameterList& params,  //!< ParameterList for communication between control routine and elements
      DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
      DRT::Element::LocationArray& la //!< location array for de-assembly
      );

    // TODO 2012-07-26???
    //! @name Evaluation
    //! Evaluate an element
    //!  Evaluate So_tet4fbar element stiffness, mass, internal forces, etc.
    //!
    //! If NULL on input, the controlling method does not expect the element
    //!  to fill these matrices or vectors.
    //!
    //!  \return 0 if successful, negative otherwise
    int Evaluate(
      ParameterList& params,  //!< ParameterList for communication between control routine and elements
      DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
      DRT::Element::LocationArray& la,  //!< location array for de-assembly
      Epetra_SerialDenseMatrix& elemat1,  //!< (stiffness-)matrix to be filled by element.
      Epetra_SerialDenseMatrix& elemat2,  //!< (mass-)matrix to be filled by element.
      Epetra_SerialDenseVector& elevec1,  //!< (internal force-)vector to be filled by element
      Epetra_SerialDenseVector& elevec2,  //!< vector to be filled by element
      Epetra_SerialDenseVector& elevec3  //!< vector to be filled by element
      );

    //! init the inverse of the jacobian and its determinant in the material
    //! configuration
    virtual void InitJacobianMapping();

    //@}

  protected:

    //! don't want = operator
    So3_Thermo& operator = (const So3_Thermo& old);

    //! Evaluate an element
    //!
    //! Evaluate So3_Thermo element stiffness, mass, internal forces, etc.
    //! Templated evaluate routine of element matrixes
    //!
    //! If NULL on input, the controlling method does not expect the element
    //! to fill these matrices or vectors.
    //!
    //! \return 0 if successful, negative otherwise
    int CouplThrEvaluate(
      ParameterList& params,  //!< ParameterList for communication between control routine and elements
      DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
      DRT::Element::LocationArray& la,  //!< location array for de-assembly
      Epetra_SerialDenseMatrix& elemat1,  //!< (stiffness-)matrix to be filled by element.
      Epetra_SerialDenseMatrix& elemat2,  //!< (mass-)matrix to be filled by element.
      Epetra_SerialDenseVector& elevec1,  //!< (internal force-)vector to be filled by element
      Epetra_SerialDenseVector& elevec2,  //!< vector to be filled by element
      Epetra_SerialDenseVector& elevec3  //!< vector to be filled by element
      );

    // TODO 2012-07-26 which actions are required for thermo???
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
      calc_struct_stifftemp,
      calc_struct_update_istep,
      calc_struct_reset_istep,  //!< reset elementwise internal variables
                                //!< during iteration to last converged state
      calc_struct_reset_discretization,  //!< reset elementwise internal variables
                                         //!< to state in the beginning of the computation
      calc_struct_energy,       //!< compute internal energy
      calc_struct_errornorms,   //!< compute error norms (L2,H1,energy)
      postprocess_stress,       //!< postprocess stresses
      prestress_update,
    };

    //! vector of inverses of the jacobian in material frame
    vector<LINALG::Matrix<numdim_,numdim_> > invJ_;
    //! determinant of Jacobian in material frame
    vector<double> detJ_;
    //! container containing all kind of solid3 thisandthat
    DRT::Container data_;
    //! vector of coordinates of current integration point in reference coordinates
    vector<LINALG::Matrix<numdim_,1> > xsi_;
    //! gauss integration rule
    DRT::UTILS::GaussIntegration intpoints_;

    //! flag indicating higher order element
    bool ishigherorder_;

  private:
    DRT::Node** Nodes();

    Teuchos::RCP<MAT::Material> Material() const;

    int Id() const;

  }; // class So3_Thermo

} // namespace ELEMENTS
} // namespace DRT


#endif  // SO3_THERMO_H
