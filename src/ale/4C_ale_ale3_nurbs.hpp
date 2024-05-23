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

namespace DRT
{
  namespace ELEMENTS
  {
    namespace NURBS
    {
      class Ale3NurbsType : public DRT::ElementType
      {
       public:
        std::string Name() const override { return "Ale3_NurbsType"; }

        static Ale3NurbsType& Instance();

        CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

        Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
            const int id, const int owner) override;

        Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

        void nodal_block_information(
            DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

        CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
            DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

        void setup_element_definition(
            std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
            override
        {
          // do nothing. Definition inserted by normal wall element.
        }

       private:
        static Ale3NurbsType instance_;
      };

      /*----------------------------------------------------------------------------*/
      /*----------------------------------------------------------------------------*/
      class Ale3Nurbs : public DRT::ELEMENTS::Ale3
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
        int UniqueParObjectId() const override
        {
          return Ale3NurbsType::Instance().UniqueParObjectId();
        }

        virtual DRT::ElementType& ElementType() { return Ale3NurbsType::Instance(); }

        /// Print this element
        void Print(std::ostream& os) const override;


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
          return (0);
          /*
          if (NumNode()==27 || NumNode()==8)
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
        int NumSurface() const override
        {
          if (NumNode() == 27 || NumNode() == 8)
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

    }  // namespace NURBS
  }    // namespace ELEMENTS
}  // namespace DRT


FOUR_C_NAMESPACE_CLOSE

#endif
