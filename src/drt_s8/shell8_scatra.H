/*----------------------------------------------------------------------------*/
/*!
\brief A 2D shell element with scatra functionality

\level 1

\maintainer Christoph Schmidt

*/
/*---------------------------------------------------------------------------*/
#ifndef SHELL8_SCATRA_H_
#define SHELL8_SCATRA_H_

#include "shell8.H"

#include "../drt_inpar/inpar_scatra.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class Shell8ScatraType : public Shell8Type
    {
     public:
      std::string Name() const { return "Shell8ScatraType"; }

      static Shell8ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static Shell8ScatraType instance_;
    };

    class Shell8_Scatra : public Shell8
    {
      typedef DRT::ELEMENTS::Shell8 my;

     public:
      /// @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owner
      */
      Shell8_Scatra(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Shell8_Scatra(const Shell8_Scatra& old);

      /*!
      \brief Destructor
      */
      virtual ~Shell8_Scatra(){};

      /*!
      \brief Deep copy this instance of Shell8_Scatra and return pointer to the copy

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
        return Shell8ScatraType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element
      */
      void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element
      */
      void Unpack(const std::vector<char>& data);

      //@}

      //! @name Access methods

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      /*!
      \brief Return elementtype instance
      */
      virtual DRT::ElementType& ElementType() const { return Shell8ScatraType::Instance(); }

      //@}

      /// @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      /// @name params
      /// return SCATRA::ImplType
      const INPAR::SCATRA::ImplType& ImplType() const { return impltype_; };

     private:
      //! implementation type (physics)
      INPAR::SCATRA::ImplType impltype_;
      //@}

     protected:
      //! don't want = operator
      Shell8_Scatra& operator=(const Shell8_Scatra& old);
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SHELL8_SCATRA_H_ */
