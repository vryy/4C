/*----------------------------------------------------------------------*/
/*! \file

\brief Base class for the immersed fluid elements


\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_config.hpp"

#include "4C_fluid_ele.hpp"

#ifndef FOUR_C_FLUID_ELE_IMMERSED_BASE_HPP
#define FOUR_C_FLUID_ELE_IMMERSED_BASE_HPP

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    class FluidTypeImmersedBase : public FluidType
    {
     public:
      /*!
      \brief Decide which element type should be created : FLUIDIMMERSED or FLUIDPOROIMMERSED
      */
      CORE::COMM::ParObject* Create(const std::vector<char>& data) override = 0;

      /*!
      \brief Decide which element type should be created : FLUIDIMMERSED or FLUIDPOROIMMERSED

      */
      Teuchos::RCP<CORE::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      /*!
      \brief Decide which element type should be created : FLUIDIMMERSED or FLUIDPOROIMMERSED
      */
      Teuchos::RCP<CORE::Elements::Element> Create(const int id, const int owner) override = 0;

      /*!
      \brief Setup the definition line for this element
      */
      void setup_element_definition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override = 0;


    };  // class FluidTypeImmersed


    class FluidImmersedBase : public virtual Fluid
    {
     public:
      //@}
      //! @name constructors and destructors and related methods

      /*!
      \brief standard constructor
      */
      FluidImmersedBase(int id,  ///< A unique global id
          int owner              ///< ???
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      FluidImmersedBase(const FluidImmersedBase& old);

      /*!
      \brief Deep copy this instance of fluid and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      CORE::Elements::Element* Clone() const override
      {
        FOUR_C_THROW("not implemented in base class");
        return nullptr;
      };

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        FOUR_C_THROW("not implemented in base class");
        return -1234;
      }

      /*!
      \brief Each element whose nodes are all covered by an immersed dis. are set IsImmersed
      */
      virtual void SetIsImmersed(int isimmersed) { FOUR_C_THROW("not implemented in base class"); };

      /*!
      \brief Each element which has nodes covered by the immersed dis. but at least one node that is
      not covered is set BoundaryIsImmersed

       An element may also be set BoundaryIsImmersed, if an integration point of the structural
      surface lies in this element.

      */
      virtual void set_boundary_is_immersed(int IsBoundaryImmersed)
      {
        FOUR_C_THROW("not implemented in base class");
      };

      /*!
      \brief Each element, which is either Immersed or BoundaryImmersed is also set
      has_projected_dirichlet
      */
      virtual void set_has_projected_dirichlet(int has_projected_dirichletvalues)
      {
        FOUR_C_THROW("not implemented in base class");
      };

      /*!
      \brief set if divergence needs to be projected to an integration point
      */
      virtual void set_int_point_has_projected_divergence(
          int gp, int intpoint_has_projected_divergence)
      {
        FOUR_C_THROW("not implemented in base class");
      };

      /*!
      \brief store the projected divergence
      */
      virtual void store_projected_int_point_divergence(
          int gp, double projected_intpoint_divergence)
      {
        FOUR_C_THROW("not implemented in base class");
      };

      /*!
      \brief returns true if element was set IsImmersed
      */
      virtual int IsImmersed()
      {
        FOUR_C_THROW("not implemented in base class");
        return -1234;
      };

      /*!
      \brief returns true if element was set IsBundaryImmersed
      */
      virtual int IsBoundaryImmersed()
      {
        FOUR_C_THROW("not implemented in base class");
        return -1234;
      };

      /*!
      \brief returns true if element needs to get projected Dirichlet values
      */
      virtual int has_projected_dirichlet()
      {
        FOUR_C_THROW("not implemented in base class");
        return -1234;
      };

      /*!
      \brief returns true if element needs to get projected divergence at integration point

      */
      virtual int int_point_has_projected_divergence(int gp)
      {
        FOUR_C_THROW("not implemented in base class");
        return -1234;
      };

      /*!
      \brief returns projected divergence at integration point
      */
      virtual double projected_int_point_divergence(int gp)
      {
        FOUR_C_THROW("not implemented in base class");
        return -1234.0;
      };

      /*!
      \brief returns rcp to vector containing gps with projected divergence
      */
      virtual Teuchos::RCP<std::vector<int>> get_rcp_int_point_has_projected_divergence()
      {
        FOUR_C_THROW("not implemented in base class");
        return Teuchos::null;
      };

      /*!
      \brief returns rcp to vector containing projected divergence values
      */
      virtual Teuchos::RCP<std::vector<double>> get_rcp_projected_int_point_divergence()
      {
        FOUR_C_THROW("not implemented in base class");
        return Teuchos::null;
      };

      /*!
      \brief construct rcp to vector for divergence projection handling
      */
      virtual void ConstructElementRCP(int size)
      {
        FOUR_C_THROW("not implemented in base class");
        return;
      };

      /*!
      \brief Clean up the element rcp
      */
      virtual void DestroyElementRCP()
      {
        FOUR_C_THROW("not implemented in base class");
        return;
      };

      /*!
      \brief Query data to be visualized using BINIO of a given name

      This method is to be overloaded by a derived method.
      The derived method is supposed to call this base method to visualize the owner of
      the element.
      If the derived method recognizes a supported data name, it shall fill it
      with corresponding data.
      If it does NOT recognizes the name, it shall do nothing.

      \warning The method must not change size of variable data

      \param name (in):   Name of data that is currently processed for visualization
      \param data (out):  data to be filled by element if it recognizes the name
      */
      bool VisData(const std::string& name, std::vector<double>& data) override
      {
        FOUR_C_THROW("not implemented in base class");
        return false;
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(CORE::COMM::PackBuffer& data) const override
      {
        FOUR_C_THROW("not implemented in base class");
      };

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element
      */
      void Unpack(const std::vector<char>& data) override
      {
        FOUR_C_THROW("not implemented in base class");
      };



    };  // class FluidImmersedBase
  }     // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
