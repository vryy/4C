/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the 3D solid-poro element including scatra functionality

 \level 2

 *----------------------------------------------------------------------*/

#ifndef BACI_SO3_PORO_SCATRA_ELETYPES_H
#define BACI_SO3_PORO_SCATRA_ELETYPES_H

#include "baci_config.H"

#include "baci_so3_poro_eletypes.H"

BACI_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_hex8PoroScatraType : public So_hex8PoroType
    {
     public:
      std::string Name() const override { return "So_hex8PoroScatraType"; }

      static So_hex8PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex8PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8POROSCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_tet4PoroScatraType : public So_tet4PoroType
    {
     public:
      std::string Name() const override { return "So_tet4PoroScatraType"; }

      static So_tet4PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet4PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDT4POROSCATRA"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_hex27PoroScatraType : public So_hex27PoroType
    {
     public:
      std::string Name() const override { return "So_hex27PoroScatraType"; }

      static So_hex27PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex27PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH27POROSCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_tet10PoroScatraType : public So_tet10PoroType
    {
     public:
      std::string Name() const override { return "So_tet10PoroScatraType"; }

      static So_tet10PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet10PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDT10POROSCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  NURBS 27 Element                                      schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_nurbs27PoroScatraType : public So_nurbs27PoroType
    {
     public:
      std::string Name() const override { return "So_nurbs27PoroScatraType"; }

      static So_nurbs27PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_nurbs27PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SONURBS27POROSCATRA"; }
    };


  }  // namespace ELEMENTS
}  // namespace DRT
BACI_NAMESPACE_CLOSE

#endif  // SO3_PORO_SCATRA_ELETYPES_H
