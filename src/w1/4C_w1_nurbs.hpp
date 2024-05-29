/*----------------------------------------------------------------------------*/
/*! \file
\brief 2D solid-wall elements using NURBS shape functions.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_W1_NURBS_HPP
#define FOUR_C_W1_NURBS_HPP

#include "4C_config.hpp"

#include "4C_w1.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    namespace NURBS
    {
      class Wall1NurbsType : public CORE::Elements::ElementType
      {
       public:
        std::string Name() const override { return "Wall1NurbsType"; }

        static Wall1NurbsType& Instance();

        CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

        Teuchos::RCP<CORE::Elements::Element> Create(const std::string eletype,
            const std::string eledistype, const int id, const int owner) override;

        Teuchos::RCP<CORE::Elements::Element> Create(const int id, const int owner) override;

        void setup_element_definition(
            std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
            override;

        void nodal_block_information(
            CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

        CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
            CORE::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

       private:
        static Wall1NurbsType instance_;
      };


      class Wall1Nurbs : public DRT::ELEMENTS::Wall1
      {
       public:
        /*!
        \brief Standard Constructor

        \param id    : A unique global id
        \param owner : proc id that will own this element
        */
        Wall1Nurbs(int id, int owner);

        /*!
        \brief Copy Constructor

        Makes a deep copy of a Element

        */
        Wall1Nurbs(const Wall1Nurbs& old);


        /*!
        \brief Deep copy this instance of Wall1 and return pointer to it
        */

        CORE::Elements::Element* Clone() const override;


        /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of this file.

        \return my parobject id
        */
        int UniqueParObjectId() const override
        {
          return Wall1NurbsType::Instance().UniqueParObjectId();
        }


        /// Print this element
        void Print(std::ostream& os) const override;

        CORE::Elements::ElementType& ElementType() const override
        {
          return Wall1NurbsType::Instance();
        }

        /*!
        \brief Get shape type of element

        \return nurbs4 or nurbs9

        */
        CORE::FE::CellType Shape() const override;


        /*!
        \brief Return number of lines of this element.
        */
        int NumLine() const override
        {
          if (num_node() == 9 || num_node() == 4)
          {
            return 4;
          }
          else
          {
            FOUR_C_THROW("Could not determine number of lines");
            return -1;
          }
        }


        /*!
        \brief Get vector of Teuchos::RCPs to the lines of this element
        */
        std::vector<Teuchos::RCP<CORE::Elements::Element>> Lines() override;


        /*!
        \brief Get vector of Teuchos::RCPs to the surfaces of this element
        */
        std::vector<Teuchos::RCP<CORE::Elements::Element>> Surfaces() override;


       private:
      };

    }  // namespace NURBS
  }    // namespace ELEMENTS
}  // namespace DRT


FOUR_C_NAMESPACE_CLOSE

#endif
