/*----------------------------------------------------------------------------*/
/*! \file
\brief Weakly Compressible fluid element based on the HDG method

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_FLUID_ELE_HDG_WEAK_COMP_H
#define BACI_FLUID_ELE_HDG_WEAK_COMP_H

#include "baci_discretization_fem_general_utils_polynomial.H"
#include "baci_fluid_ele.H"
#include "baci_lib_dg_element.h"

#include <Teuchos_RCP.hpp>

namespace DRT
{
  namespace INPUT
  {
    class LineDefinition;
  }


  class Discretization;

  namespace ELEMENTS
  {
    class FluidHDGWeakCompType : public FluidType
    {
     public:
      std::string Name() const { return "FluidHDGWeakCompType"; }

      static FluidHDGWeakCompType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static FluidHDGWeakCompType instance_;
    };


    /*!
    \brief HDG weakly compressible fluid element
    */
    class FluidHDGWeakComp : public Fluid, public DG_Element
    {
     public:
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      FluidHDGWeakComp(int id,  ///< A unique global id
          int owner             ///< ???
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidHDGWeakComp(const FluidHDGWeakComp& old);

      /*!
      \brief Deep copy this instance of fluid and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Destructor

      */
      virtual ~FluidHDGWeakComp();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const
      {
        return FluidHDGWeakCompType::Instance().UniqueParObjectId();
      }

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
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      //! @name Access methods
      /*!
      \brief Get number of degrees of freedom per node, zero for the primary
      and the secondary dof set and equal to the given number for the tertiary dof set
      */
      virtual int NumDofPerNode(const DRT::Node&) const { return 0; }

      /*!
      \brief Returns the number of dofs per node for the ALE displacements
       */
      int NumDofPerNodeAuxiliary() const
      {
        return (CORE::DRT::UTILS::getDimension(distype_) + 0);
        // maybe + 1 to mimic the standard fluid dofset with velocity + pressure
        // this is maybe needed for the fluid-ALE coupling
      }

      /*!
      \brief Get number of degrees of freedom per face

      */
      virtual int NumDofPerFace(const unsigned face) const
      {
        return (1 + CORE::DRT::UTILS::getDimension(distype_)) * NumDofPerComponent(face);
      }

      /*!
      \brief Get number of dofs per component per face
      */
      virtual int NumDofPerComponent(const unsigned face) const
      {
        return CORE::DRT::UTILS::getBasisSize(
            CORE::DRT::UTILS::getEleFaceShapeType(distype_), this->Degree(), completepol_);
      }

      /*!
      \brief Get number of degrees of freedom per element, zero for the primary dof set
      and equal to the given number for the secondary dof set
      */
      virtual int NumDofPerElement() const { return 0; }

      /*!
      \brief Returns the number of dofs per element for the interior variables
       */
      int NumDofPerElementAuxiliary() const
      {
        const int nsd_ = CORE::DRT::UTILS::getDimension(distype_);
        const int msd_ = (nsd_ * (nsd_ + 1.0)) / 2.0;
        return (msd_ + nsd_ + 1.0) *
               CORE::DRT::UTILS::getBasisSize(distype_, degree_, completepol_);
      }

      /*!
       \brief Returns the degree of the element
       */
      int Degree() const { return degree_; }

      /*!
       \brief Returns the degree of the element
       */
      int UsesCompletePolynomialSpace() const { return completepol_; }

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element, that is, call the element routines to evaluate fluid
      element matrices and vectors or evaluate errors, statistics or updates etc. directly.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If nullptr on input,
                              the controlling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If nullptr on input,
                              the controlling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector& elevec3);

      //@}

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return FluidHDGWeakCompType::Instance(); }

     private:
      // don't want = operator
      FluidHDGWeakComp& operator=(const FluidHDGWeakComp& old);

      // stores the degree of the element
      unsigned char degree_;

      // stores the polynomial type (tensor product or complete polynomial)
      bool completepol_;
    };  // class Fluid


  }  // namespace ELEMENTS
}  // namespace DRT



#endif
