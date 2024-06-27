/*----------------------------------------------------------------------------*/
/*! \file
\brief 2D wall element for structure part of porous medium including scatra functionality.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_W1_PORO_SCATRA_HPP
#define FOUR_C_W1_PORO_SCATRA_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_w1_poro.hpp"
#include "4C_w1_poro_scatra_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 2 dimensional solid element with modifications for porous media

    */
    template <Core::FE::CellType distype>
    class Wall1PoroScatra : public Wall1Poro<distype>
    {
      typedef Discret::ELEMENTS::Wall1Poro<distype> my;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      Wall1PoroScatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Wall1PoroScatra(const Wall1PoroScatra& old);

      //@}

      /*!
      \brief Deep copy this instance of Wall1_Poro_Scatra and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed
      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        int parobjectid(-1);
        switch (distype)
        {
          case Core::FE::CellType::tri3:
          {
            parobjectid = Discret::ELEMENTS::WallTri3PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          case Core::FE::CellType::quad4:
          {
            parobjectid =
                Discret::ELEMENTS::WallQuad4PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          case Core::FE::CellType::quad9:
          {
            parobjectid =
                Discret::ELEMENTS::WallQuad9PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          case Core::FE::CellType::nurbs4:
          {
            parobjectid =
                Discret::ELEMENTS::WallNurbs4PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          case Core::FE::CellType::nurbs9:
          {
            parobjectid =
                Discret::ELEMENTS::WallNurbs9PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          default:
          {
            FOUR_C_THROW("unknown element type");
            break;
          }
        }
        return parobjectid;
      };

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element
      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element
      */
      void unpack(const std::vector<char>& data) override;

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void print(std::ostream& os) const override;

      /*!
      \brief Return elementtype instance
      */
      Core::Elements::ElementType& ElementType() const override
      {
        switch (distype)
        {
          case Core::FE::CellType::tri3:
            return Discret::ELEMENTS::WallTri3PoroScatraType::Instance();
            break;
          case Core::FE::CellType::quad4:
            return Discret::ELEMENTS::WallQuad4PoroScatraType::Instance();
            break;
          case Core::FE::CellType::quad9:
            return Discret::ELEMENTS::WallQuad9PoroScatraType::Instance();
            break;
          case Core::FE::CellType::nurbs4:
            return Discret::ELEMENTS::WallNurbs4PoroScatraType::Instance();
            break;
          case Core::FE::CellType::nurbs9:
            return Discret::ELEMENTS::WallNurbs9PoroScatraType::Instance();
            break;
          default:
            FOUR_C_THROW("unknown element type");
            break;
        }
        return Discret::ELEMENTS::WallQuad4PoroScatraType::Instance();
      };

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& eledistype,
          Input::LineDefinition* linedef) override;

      //@}

      /// @name params
      /*!
      \brief Return the SCATRA ImplType
      */
      const Inpar::ScaTra::ImplType& ImplType() const { return impltype_; };

     private:
      //! implementation type (physics)
      Inpar::ScaTra::ImplType impltype_;
      //@}

     protected:
      //! don't want = operator
      Wall1PoroScatra& operator=(const Wall1PoroScatra& old);

    };  // class Wall1_Poro_Scatra

  }  // namespace ELEMENTS
}  // namespace Discret
FOUR_C_NAMESPACE_CLOSE

#endif
