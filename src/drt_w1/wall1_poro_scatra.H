/*----------------------------------------------------------------------------*/
/*!
\brief 2D wall element for structure part of porous medium including scatra functionality.

\level 2

\maintainer Christoph Schmidt

*/
/*---------------------------------------------------------------------------*/

#ifndef WALL1_SCATRA_PORO_H_
#define WALL1_SCATRA_PORO_H_

#include "wall1_poro.H"
#include "wall1_poro_scatra_eletypes.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief A C++ version of a 2 dimensional solid element with modifications for porous media

    */
    template <DRT::Element::DiscretizationType distype>
    class Wall1_Poro_Scatra : public Wall1_Poro<distype>
    {
      typedef DRT::ELEMENTS::Wall1_Poro<distype> my;

     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      Wall1_Poro_Scatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Wall1_Poro_Scatra(const Wall1_Poro_Scatra& old);

      /*!
      \brief Destructor
      */
      virtual ~Wall1_Poro_Scatra(){};

      //@}

      /*!
      \brief Deep copy this instance of Wall1_Poro_Scatra and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed
      */
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const
      {
        int parobjectid(-1);
        switch (distype)
        {
          case DRT::Element::tri3:
          {
            parobjectid = DRT::ELEMENTS::WallTri3PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          case DRT::Element::quad4:
          {
            parobjectid = DRT::ELEMENTS::WallQuad4PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          case DRT::Element::quad9:
          {
            parobjectid = DRT::ELEMENTS::WallQuad9PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          case DRT::Element::nurbs4:
          {
            parobjectid = DRT::ELEMENTS::WallNurbs4PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          case DRT::Element::nurbs9:
          {
            parobjectid = DRT::ELEMENTS::WallNurbs9PoroScatraType::Instance().UniqueParObjectId();
            break;
          }
          default:
          {
            dserror("unknown element type");
            break;
          }
        }
        return parobjectid;
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

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      /*!
      \brief Return elementtype instance
      */
      virtual DRT::ElementType& ElementType() const
      {
        switch (distype)
        {
          case DRT::Element::tri3:
            return DRT::ELEMENTS::WallTri3PoroScatraType::Instance();
            break;
          case DRT::Element::quad4:
            return DRT::ELEMENTS::WallQuad4PoroScatraType::Instance();
            break;
          case DRT::Element::quad9:
            return DRT::ELEMENTS::WallQuad9PoroScatraType::Instance();
            break;
          case DRT::Element::nurbs4:
            return DRT::ELEMENTS::WallNurbs4PoroScatraType::Instance();
            break;
          case DRT::Element::nurbs9:
            return DRT::ELEMENTS::WallNurbs9PoroScatraType::Instance();
            break;
          default:
            dserror("unknown element type");
            break;
        }
        return DRT::ELEMENTS::WallQuad4PoroScatraType::Instance();
      };

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      virtual bool ReadElement(const std::string& eletype, const std::string& eledistype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      /// @name params
      /*!
      \brief Return the SCATRA ImplType
      */
      const INPAR::SCATRA::ImplType& ImplType() const { return impltype_; };

     private:
      //! implementation type (physics)
      INPAR::SCATRA::ImplType impltype_;
      //@}

     protected:
      //! don't want = operator
      Wall1_Poro_Scatra& operator=(const Wall1_Poro_Scatra& old);

    };  // class Wall1_Poro_Scatra

  }  // namespace ELEMENTS
}  // namespace DRT
#endif /* WALL1_PORO_SCATRA_H_ */
