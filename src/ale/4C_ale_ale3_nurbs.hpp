/*----------------------------------------------------------------------------*/
/*! \file

\brief A nurbs implementation of the ale3 element

\level 2

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef FOUR_C_ALE_ALE3_NURBS_HPP
#define FOUR_C_ALE_ALE3_NURBS_HPP

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_config.hpp"

#include "4C_ale_ale3.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    namespace Nurbs
    {
      class Ale3NurbsType : public Core::Elements::ElementType
      {
       public:
        std::string name() const override { return "Ale3_NurbsType"; }

        static Ale3NurbsType& instance();

        Core::Communication::ParObject* create(const std::vector<char>& data) override;

        Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
            const std::string eledistype, const int id, const int owner) override;

        Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

        void nodal_block_information(
            Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

        Core::LinAlg::SerialDenseMatrix compute_null_space(
            Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

        void setup_element_definition(
            std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
            override
        {
          // do nothing. Definition inserted by normal wall element.
        }

       private:
        static Ale3NurbsType instance_;
      };

      /*----------------------------------------------------------------------------*/
      /*----------------------------------------------------------------------------*/
      class Ale3Nurbs : public Discret::ELEMENTS::Ale3
      {
       public:
        /*!
        \brief Standard Constructor

        \param id    : A unique global id
        \param owner : proc id that will own this element
        */
        Ale3Nurbs(int id, int owner);

        /*!
        \brief Copy Constructor

        Makes a deep copy of a Element

        */
        Ale3Nurbs(const Ale3Nurbs& old);



        /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of this file.

        \return my parobject id
        */
        int unique_par_object_id() const override
        {
          return Ale3NurbsType::instance().unique_par_object_id();
        }

        virtual Core::Elements::ElementType& element_type() { return Ale3NurbsType::instance(); }

        /// Print this element
        void print(std::ostream& os) const override;


        /*!
        \brief Get shape type of element

        \return nurbs4 or nurbs9

        */
        Core::FE::CellType shape() const override;


        /*!
        \brief Return number of lines of this element.
        */
        int num_line() const override
        {
          return (0);
          /*
          if (num_node()==27 || num_node()==8)
          {
            return 12;
          }
          else {
            FOUR_C_THROW("Could not determine number of lines");
            return -1;
          }
          */
        }


        /*!
        \brief Return number of surfaces of this element
        */
        int num_surface() const override
        {
          if (num_node() == 27 || num_node() == 8)
          {
            return 6;
          }
          else
          {
            FOUR_C_THROW("Could not determine number of surfaces");
            return -1;
          }
        }

       private:
      };

    }  // namespace Nurbs
  }    // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
