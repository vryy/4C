/*----------------------------------------------------------------------*/
/*! \file
\brief One-dimensional artery element


\level 3

*----------------------------------------------------------------------*/
#ifndef BACI_ART_NET_ARTERY_H
#define BACI_ART_NET_ARTERY_H

#include "baci_config.H"

#include "baci_comm_parobjectfactory.H"
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_fluid_ele_nullspace.H"
#include "baci_inpar_bio.H"
#include "baci_lib_element.H"
#include "baci_lib_elementtype.H"
#include "baci_lib_node.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// forward declarations
struct _MATERIAL;



namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class ArteryType : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "ArteryType"; }

      static ArteryType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override
      {
        numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
        dimns = numdf;
        nv = numdf;
      }

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override
      {
        return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
      }

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static ArteryType instance_;
    };

    /*!
    \brief A C++ wrapper for the artery element

    */
    class Artery : public DRT::Element
    {
     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      Artery(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Artery(const Artery& old);

      /*!
      \brief Deep copy this instance of Artery and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      CORE::FE::CellType Shape() const override;

      /*!
      \brief Return number of lines of this element
      */
      int NumLine() const override
      {
        if (NumNode() == 2)
          return 1;
        else
        {
          dserror("Could not determine number of lines");
          return -1;
        }
      }

      /*!
      \brief Return number of surfaces of this element (always 1)
      */
      int NumSurface() const override { return -1; }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      int NumVolume() const override { return -1; }

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return ArteryType::Instance().UniqueParObjectId(); }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      /*!
       \brief Get vector of Teuchos::RCPs to the lines of this element
       */
      std::vector<Teuchos::RCP<DRT::Element>> Lines() override;


      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const DRT::Node& node) const override
      {
        switch (impltype_)
        {
          case INPAR::ARTDYN::impltype_lin_exp:
          {
            return 2;
            break;
          }
          case INPAR::ARTDYN::impltype_pressure_based:
          {
            return 1;
            break;
          }
          default:
          {
            dserror("Defined problem type %d does not exist!!", impltype_);
            break;
          }
        }

        return 0;
      }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int NumDofPerElement() const override { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      ArteryType& ElementType() const override { return ArteryType::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          INPUT::LineDefinition* linedef) override;

      /**
       * Set diameter in material
       * @param diam: diameter to be set
       */
      void SetDiamInMaterial(const double diam);

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      An element derived from this class uses the Evaluate method to receive commands
      and parameters from some control routine in params and evaluates element matrices and
      vectors accoring to the command in params.

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param la (in)            : location data for all dofsets of the discretization
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
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2,
          CORE::LINALG::SerialDenseVector& elevec3) override;

      int ScatraEvaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector& elevec3);

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
      int EvaluateNeumann(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a line Neumann condition on the artery element

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateDirichlet(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1);

      /*!
      \brief return implementation type (physics)

      */
      INPAR::ARTDYN::ImplType ImplType() { return impltype_; }


      //@}


      //! @name Other

      CORE::FE::GaussRule1D GaussRule() const { return gaussrule_; }


      //@}


     private:
      //! implementation type (physics)
      INPAR::ARTDYN::ImplType impltype_;

      //! Gaussrule
      CORE::FE::GaussRule1D gaussrule_;

      //! container containing all kind of Artery thisandthat
      DRT::Container data_;


      // internal calculation methods

      // don't want = operator
      Artery& operator=(const Artery& old);


      /// set number of gauss points to element shape default
      CORE::FE::GaussRule1D getOptimalGaussrule(const CORE::FE::CellType& distype);

      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool isHigherOrderElement(const CORE::FE::CellType distype  ///< discretization type
      ) const;


    };  // class Artery


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


  }  // namespace ELEMENTS
}  // namespace DRT



BACI_NAMESPACE_CLOSE

#endif  // ART_NET_ARTERY_H
