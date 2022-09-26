/*! \file

\brief A C++ wrapper for the solid element

This file contains the element-specific service routines such as
Pack, Unpack, NumDofPerNode etc.

In addition to that, it contains the interface between element call
and Gauss point loop as well as some additional service routines.

\level 1
*/

#ifndef SO_ELE_H
#define SO_ELE_H

#include <memory>
#include "../drt_lib/drt_elementtype.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_element.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"
#include "solid_ele_calc.H"
#include "solid_ele_eas_utils.H"

namespace MAT
{
  class So3Material;
}
namespace DRT::ELEMENTS
{
  // forward declaration
  class SolidEleInterface;

  class SolidType : public DRT::ElementType
  {
   public:
    /// setup the dat file input line definitions for this type of element
    void SetupElementDefinition(
        std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
        override;

    /// create an element from a dat file specifier
    Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
        const int id, const int owner) override;

    /// create an empty element
    Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

    // create and unpack element from data
    DRT::ParObject* Create(const std::vector<char>& data) override;

    std::string Name() const override { return "SolidType"; }

    void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

    Epetra_SerialDenseMatrix ComputeNullSpace(
        DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

    static SolidType& Instance();

   private:
    static SolidType instance_;

  };  // class SolidType

  class Solid : public DRT::Element
  {
    friend class SolidType;

   public:
    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor

    \param id    (in): A globally unique element id
    \param owner (in): owner processor of the element
    */
    Solid(int id, int owner);

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element
    */
    Solid(const DRT::ELEMENTS::Solid& old);

    /*!
    \brief Deep copy
    */
    [[nodiscard]] DRT::Element* Clone() const override;


    /*!
    \brief Return unique ParObject id
    */
    [[nodiscard]] int UniqueParObjectId() const override
    {
      return SolidType::Instance().UniqueParObjectId();
    };

    /*!
    \brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Pack(DRT::PackBuffer& data) const override;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this element

    */
    void Unpack(const std::vector<char>& data) override;

    /// return ElementType instance
    [[nodiscard]] DRT::ElementType& ElementType() const override { return SolidType::Instance(); }

    //@}

    //! @name Query methods

    /*!
    \brief Get shape type of element
    */
    [[nodiscard]] DRT::Element::DiscretizationType Shape() const override { return distype_; };

    /*!
    \brief Get kinematic type of element
    */
    [[nodiscard]] INPAR::STR::KinemType KinematicType() const { return kintype_; }

    /*!
    \brief Set kinematic type of element
    */
    void SetKinematicType(INPAR::STR::KinemType kintype) { kintype_ = kintype; }

    [[nodiscard]] virtual Teuchos::RCP<MAT::So3Material> SolidMaterial(int nummat = 0) const;

    /*!
    \brief Return number of lines to this element
    */
    [[nodiscard]] int NumLine() const override;

    /*!
    \brief Return number of surfaces to this element
    */
    [[nodiscard]] int NumSurface() const override;

    /*!
    \brief Return number of volumes to this element
    */
    [[nodiscard]] int NumVolume() const override;

    /*!
    \brief Get vector of Teuchos::RCPs to the lines of this element
    */
    std::vector<Teuchos::RCP<DRT::Element>> Lines() override;

    /*!
    \brief Get vector of Teuchos::RCPs to the surfaces of this element
    */
    std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;

    /*!
    \brief Get vector of Teuchos::RCPs to the volume of this element
    */
    std::vector<Teuchos::RCP<DRT::Element>> Volumes() override;

    /*!
    \brief Get number of degrees of freedom of a certain node
    */
    [[nodiscard]] int NumDofPerNode(const DRT::Node& node) const override
    {
      return 3;
    }  // todo fix this for 2D elements


    /*!
    \brief Get number of degrees of freedom per element
           (implements pure virtual DRT::Element)

    The element decides how many element degrees of freedom it has.
    It can redecide along the way of a simulation.

    \note Element degrees of freedom mentioned here are dofs that are visible
          at the level of the total system of equations. Purely internal
          element dofs that are condensed internally should NOT be considered.
    */
    [[nodiscard]] int NumDofPerElement() const override { return 0; }



    //! @name Input and Creation

    /*!
    \brief Read input for this element
    */
    bool ReadElement(const std::string& eletype, const std::string& distype,
        DRT::INPUT::LineDefinition* linedef) override;
    /// @}

    //! @name Evaluation

    /*!
    \brief Evaluate an element

    \param params (in/out): ParameterList for communication between control routine
                            and elements
    \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                            the controling method does not epxect the element to fill
                            this matrix.
    \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                            the controling method does not epxect the element to fill
                            this matrix.
    \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                            the controlling method does not epxect the element
                            to fill this vector
    \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                            the controlling method does not epxect the element
                            to fill this vector
    \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                            the controlling method does not epxect the element
                            to fill this vector
    \return 0 if successful, negative otherwise
    */
    int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
        Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
        Epetra_SerialDenseVector& elevec3) override;

    /*!
    \brief Evaluate Neumann boundary condition

    \param params (in/out)    : ParameterList for communication between control routine
                                and elements
    \param discretization (in): reference to the underlying discretization
    \param condition (in)     : condition to be evaluated
    \param lm (in)            : location vector of this element
    \param elevec1 (out)      : vector to be filled by element. If NULL on input,

    \return 0 if successful, negative otherwise
    */
    int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
        Epetra_SerialDenseMatrix* elemat1 = nullptr) override;

    //@}

    /** \brief get access to the parameter interface pointer
     *
     *  \author hiermeier
     *  \date 04/16 */
    Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr() override
    {
      return interface_ptr_;
    }

    /** \brief returns true if the parameter interface is defined and initialized, otherwise false
     *
     *  \author hiermeier
     *  \date 04/16 */
    [[nodiscard]] inline bool IsParamsInterface() const override
    {
      return (not interface_ptr_.is_null());
    }

    /** \brief get access to the interface
     *
     *  \author hiermeier
     *  \date 04/16 */
    inline STR::ELEMENTS::ParamsInterface& ParamsInterface()
    {
      if (not IsParamsInterface()) dserror("The interface ptr is not set!");
      return *interface_ptr_;
    }

    /** \brief set the parameter interface ptr for the solid elements
     *
     *  \param p (in): Parameter list coming from the time integrator.
     *
     *  \author hiermeier
     *  \date 04/16 */
    void SetParamsInterfacePtr(const Teuchos::ParameterList& p) override;

    virtual void SetEAStype(::STR::ELEMENTS::EASType type) { eastype_ = type; }

    virtual void SetElementCondensator(Teuchos::RCP<STR::ELEMENTS::CondensatorBase> data)
    {
      element_condensator_ = data;
    }

    virtual Teuchos::RCP<STR::ELEMENTS::CondensatorBase> GetElementCondensator()
    {
      return element_condensator_;
    }

    virtual ::STR::ELEMENTS::EASType GetEAStype() { return eastype_; }
    virtual const std::set<INPAR::STR::EleTech>& GetEleTech() { return eletech_; }
    virtual bool HaveEAS() { return eletech_.find(INPAR::STR::EleTech::eas) != eletech_.end(); }
    void VisNames(std::map<std::string, int>& names) override;
    bool VisData(const std::string& name, std::vector<double>& data) override;

   protected:
    /*!
    \brief Anything that needs to be done before the standard Evaluate.
    Empty function in the base class that may be overloaded in derived elements.
    */
    virtual int PreEvaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
        Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
        Epetra_SerialDenseVector& elevec3)
    {
      return 0;
    }

    /*!
    \brief Anything that needs to be done after the standard Evaluate.
    Empty function in the base class that may be overloaded in derived elements.
    */
    virtual int PostEvaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
        std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
        Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
        Epetra_SerialDenseVector& elevec3);

    //@}

   private:
    //! discretization type
    DRT::Element::DiscretizationType distype_;

    //! kinematic type
    INPAR::STR::KinemType kintype_;

    //! element techonology
    std::set<INPAR::STR::EleTech> eletech_;

    //! specify EAS type in case  contains EAS
    ::STR::ELEMENTS::EASType eastype_;

    /** \brief interface ptr
     *
     *  data exchange between the element and the time integrator. */
    Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

    //! manager for condensation at element level
    Teuchos::RCP<STR::ELEMENTS::CondensatorBase> element_condensator_;

    //! manager for condensation at Gauss point level
    Teuchos::RCP<std::vector<STR::ELEMENTS::CondensatorBase>> gp_condensator_;

    //! element calculation interface
    std::unique_ptr<SolidEleInterface> solid_interface_ = nullptr;

    //! flag, whether the post setup of materials is already called
    bool material_post_setup_;
  };  // class Solid

}  // namespace DRT::ELEMENTS

#endif  // SO_ELE_H
