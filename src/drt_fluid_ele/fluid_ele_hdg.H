/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_hdg.H

\brief Fluid element based on the HDG method

\level 1

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
*/
/*----------------------------------------------------------------------*/

#ifndef DRT_FLUID_ELE_HDG_H
#define DRT_FLUID_ELE_HDG_H

#include "Teuchos_RCP.hpp"


#include "fluid_ele.H"
#include "../drt_fem_general/drt_utils_polynomial.H"

namespace DRT
{
  namespace INPUT
  {
    class LineDefinition;
  }


  class Discretization;

  namespace ELEMENTS
  {
    class FluidHDGType : public FluidType
    {
     public:
      std::string Name() const { return "FluidHDGType"; }

      static FluidHDGType& Instance();

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
      static FluidHDGType instance_;
    };


    /*!
    \brief HDG fluid element
    */
    class FluidHDG : public Fluid
    {
     public:
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      FluidHDG(int id,  ///< A unique global id
          int owner     ///< ???
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidHDG(const FluidHDG& old);

      /*!
      \brief Deep copy this instance of fluid and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Destructor

      */
      virtual ~FluidHDG();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const { return FluidHDGType::Instance().UniqueParObjectId(); }

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
      \brief Get number of degrees of freedom per node

      HDG element: No dofs are associated with nodes
      */
      virtual int NumDofPerNode(const DRT::Node&) const { return 0; }

      /*!
      \brief Get number of degrees of freedom per face

      */
      virtual int NumDofPerFace(const unsigned face) const
      {
        return DRT::UTILS::getDimension(distype_) * NumDofPerComponent(face);
      }

      /*!
      \brief Get number of dofs per component per face
      */
      virtual int NumDofPerComponent(const unsigned face) const
      {
        return DRT::UTILS::getBasisSize(
            DRT::UTILS::getEleFaceShapeType(distype_), this->Degree(), completepol_);
      }

      /*!
      \brief Get number of degrees of freedom per element, zero for the primary dof set
      and equal to the given number for the secondary dof set
      */
      virtual int NumDofPerElement() const { return 1; }

      /*!
       \brief Returns the degree of the element
       */
      int Degree() const { return degree_; }

      /*!
       \brief Returns the degree of the element
       */
      int UsesCompletePolynomialSpace() const { return completepol_; }

      /*!
       \brief Returns the degree of the element for the interior DG space
       */
      int NumDofPerElementAuxiliary() const
      {
        const int nsd_ = DRT::UTILS::getDimension(distype_);
        return (nsd_ * (nsd_ + 1) + 1) * DRT::UTILS::getBasisSize(distype_, degree_, completepol_) +
               1;
      }

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element, that is, call the element routines to evaluate fluid
      element matrices and vectors or evaluate errors, statistics or updates etc. directly.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                              the controlling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                              the controlling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

      //@}

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return FluidHDGType::Instance(); }

     private:
      // don't want = operator
      FluidHDG& operator=(const FluidHDG& old);

      // stores the degree of the element
      unsigned char degree_;

      // stores the polynomial type (tensor product or complete polynomial)
      bool completepol_;
    };  // class Fluid


  }  // namespace ELEMENTS
}  // namespace DRT



#endif
