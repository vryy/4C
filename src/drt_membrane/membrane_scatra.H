/*----------------------------------------------------------------------*/
/*! \file

\level 3

\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262

\brief Nonlinear Membrane Finite Element with ScaTra coupling

*----------------------------------------------------------------------*/
#ifndef MEMBRANE_SCATRA_H
#define MEMBRANE_SCATRA_H

#include "membrane.H"
#include "membrane_scatra_eletypes.H"

#include "../drt_inpar/inpar_scatra.H"


namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype>
    class MembraneScatra : public Membrane<distype>
    {
     public:
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      MembraneScatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      MembraneScatra(const MembraneScatra<distype>& old);

      /*!
      \brief Deep copy this instance of Membrane and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-constructor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const
      {
        switch (distype)
        {
          case DRT::Element::tri3:
          {
            return MembraneScatra_tri3Type::Instance().UniqueParObjectId();
          }
          case DRT::Element::tri6:
          {
            return MembraneScatra_tri6Type::Instance().UniqueParObjectId();
          }
          case DRT::Element::quad4:
          {
            return MembraneScatra_quad4Type::Instance().UniqueParObjectId();
          }
          case DRT::Element::quad9:
          {
            return MembraneScatra_quad9Type::Instance().UniqueParObjectId();
          }
          default:
            dserror("unknown element type!");
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
      \brief Destructor

      */
      virtual ~MembraneScatra(){};

      //@}

      //! @name Acess methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const
      {
        switch (distype)
        {
          case DRT::Element::tri3:
          {
            return MembraneScatra_tri3Type::Instance();
          }
          break;
          case DRT::Element::tri6:
          {
            return MembraneScatra_tri6Type::Instance();
          }
          break;
          case DRT::Element::quad4:
          {
            return MembraneScatra_quad4Type::Instance();
          }
          break;
          case DRT::Element::quad9:
          {
            return MembraneScatra_quad9Type::Instance();
          }
          break;
          default:
            dserror("unknown element type!");
            break;
        }
        // Intel compiler needs a return so
        return Membrane_quad4Type::Instance();
      };

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& eledistype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      //! @name Evaluation

      /*!
      \brief Pre-evaluate an element

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization : pointer to discretization for de-assembly
      \param la (in)        : location array for de-assembly
      */
      void PreEvaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la);

      /*!
      \brief Evaluate an element

      Evaluate Membrane element stiffness, mass, internal forces etc

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization : pointer to discretization for de-assembly
      \param la (in)        : location array for de-assembly
      \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If NULL on input,
                              the controlling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : (mass-)matrix to be filled by element. If NULL on input,
                              the controlling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : (internal force-)vector to be filled by element. If NULL on input,
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
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

      //@}

      //! @name params
      /*!
      \brief return SCATRA::ImplType
      */
      const INPAR::SCATRA::ImplType& ImplType() const { return impltype_; };

      //@}

     private:
      /*!
      \brief Get vector of ptrs to nodes
      */
      DRT::Node** Nodes();

      /*!
      \brief Get shape type of element
      */
      DRT::Element::DiscretizationType Shape() const;

      //! @{
      //! scalar transport implementation type (physics)
      INPAR::SCATRA::ImplType impltype_;
      //@}

     protected:
      //! don't want = operator
      MembraneScatra& operator=(const MembraneScatra& old);
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
