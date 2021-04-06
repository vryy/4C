/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element used for scalar transport coupling

\level 3

*/
/*---------------------------------------------------------------------------*/
#ifndef TRUSS3_SCATRA_H
#define TRUSS3_SCATRA_H

#include "truss3.H"

#include "../drt_inpar/inpar_scatra.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class Truss3ScatraType : public Truss3Type
    {
     public:
      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      static Truss3ScatraType& Instance();

      std::string Name() const override { return "Truss3ScatraType"; }

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Truss3ScatraType instance_;
    };

    /*!
     \brief three dimensional total Lagrange truss element for scalar transport coupling

     */
    class Truss3Scatra : public Truss3
    {
     public:
      /*!
       \brief Standard Constructor

       \param id    (in): A globally unique element id
       \param owner (in): owner processor of the element
       */
      Truss3Scatra(int id, int owner);

      DRT::ElementType& ElementType() const override { return Truss3ScatraType::Instance(); }

      /// return SCATRA::ImplType
      const INPAR::SCATRA::ImplType& ImplType() const { return impltype_; };

      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;


      int UniqueParObjectId() const override
      {
        return Truss3ScatraType::Instance().UniqueParObjectId();
      }

     private:
      //! scalar transport implementation type (physics)
      INPAR::SCATRA::ImplType impltype_;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef TRUSS3_SCATRA_H
