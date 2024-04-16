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
#include "baci_config.hpp"

#include "baci_ale_ale2.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    namespace NURBS
    {
      class Ale2_NurbsType : public Ale2Type
      {
       public:
        std::string Name() const override { return "Ale2_NurbsType"; }

        static Ale2_NurbsType& Instance();

        CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

        Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
            const int id, const int owner) override;

        Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

        void SetupElementDefinition(
            std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
            override
        {
          // do nothing. Definition inserted by normal wall element.
        }

       private:
        static Ale2_NurbsType instance_;
      };


      class Ale2Nurbs : public DRT::ELEMENTS::Ale2
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
          return Ale2_NurbsType::Instance().UniqueParObjectId();
        }


        /// Print this element
        void Print(std::ostream& os) const override;

        DRT::ElementType& ElementType() const override { return Ale2_NurbsType::Instance(); }

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
          if (NumNode() == 9 || NumNode() == 4)
          {
            return 4;
          }
          else
          {
            dserror("Could not determine number of lines");
            return -1;
          }
        }


       private:
      };

    }  // namespace NURBS
  }    // namespace ELEMENTS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif
