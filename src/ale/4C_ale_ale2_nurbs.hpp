/*----------------------------------------------------------------------------*/
/*! \file

\brief Nurbs verison of 2D ALE element

\level 3

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifndef FOUR_C_ALE_ALE2_NURBS_HPP
#define FOUR_C_ALE_ALE2_NURBS_HPP

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_config.hpp"

#include "4C_ale_ale2.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    namespace Nurbs
    {
      class Ale2NurbsType : public Ale2Type
      {
       public:
        std::string Name() const override { return "Ale2_NurbsType"; }

        static Ale2NurbsType& Instance();

        Core::Communication::ParObject* Create(const std::vector<char>& data) override;

        Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
            const std::string eledistype, const int id, const int owner) override;

        Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

        void setup_element_definition(
            std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
            override
        {
          // do nothing. Definition inserted by normal wall element.
        }

       private:
        static Ale2NurbsType instance_;
      };


      class Ale2Nurbs : public Discret::ELEMENTS::Ale2
      {
       public:
        /*!
        \brief Standard Constructor

        \param id    : A unique global id
        \param owner : proc id that will own this element
        */
        Ale2Nurbs(int id, int owner);

        /*!
        \brief Copy Constructor

        Makes a deep copy of a Element

        */
        Ale2Nurbs(const Ale2Nurbs& old);



        /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of this file.

        \return my parobject id
        */
        int UniqueParObjectId() const override
        {
          return Ale2NurbsType::Instance().UniqueParObjectId();
        }


        /// Print this element
        void Print(std::ostream& os) const override;

        Core::Elements::ElementType& ElementType() const override
        {
          return Ale2NurbsType::Instance();
        }

        /*!
        \brief Get shape type of element

        \return nurbs4 or nurbs9

        */
        Core::FE::CellType Shape() const override;


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


       private:
      };

    }  // namespace Nurbs
  }    // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
