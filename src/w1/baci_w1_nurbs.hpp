/*----------------------------------------------------------------------------*/
/*! \file
\brief 2D solid-wall elements using NURBS shape functions.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef BACI_W1_NURBS_HPP
#define BACI_W1_NURBS_HPP

#include "baci_config.hpp"

#include "baci_w1.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    namespace NURBS
    {
      class Wall1NurbsType : public DRT::ElementType
      {
       public:
        std::string Name() const override { return "Wall1NurbsType"; }

        static Wall1NurbsType& Instance();

        CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

        Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
            const int id, const int owner) override;

        Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

        void SetupElementDefinition(
            std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
            override;

        void NodalBlockInformation(
            DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

        CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
            DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

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

        DRT::Element* Clone() const override;


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

        DRT::ElementType& ElementType() const override { return Wall1NurbsType::Instance(); }

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


        /*!
        \brief Get vector of Teuchos::RCPs to the lines of this element
        */
        std::vector<Teuchos::RCP<DRT::Element>> Lines() override;


        /*!
        \brief Get vector of Teuchos::RCPs to the surfaces of this element
        */
        std::vector<Teuchos::RCP<DRT::Element>> Surfaces() override;


       private:
      };

    }  // namespace NURBS
  }    // namespace ELEMENTS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // W1_NURBS_H
