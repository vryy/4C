/*----------------------------------------------------------------------*/
/*! \file

\brief Solid-scatra elements types

\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_SCATRA_ELETYPES_HPP
#define FOUR_C_SO3_SCATRA_ELETYPES_HPP

#include "baci_config.hpp"

#include "baci_so3_hex27.hpp"
#include "baci_so3_hex8.hpp"
#include "baci_so3_hex8fbar.hpp"
#include "baci_so3_tet10.hpp"
#include "baci_so3_tet4.hpp"
#include "baci_so3_weg6.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                       |
     *----------------------------------------------------------------------*/
    class So_hex8ScatraType : public So_hex8Type
    {
     public:
      std::string Name() const override { return "So_hex8ScatraType"; }

      static So_hex8ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex8ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8SCATRA"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 8 fbar Element                                        Thon 12/14 |
     *----------------------------------------------------------------------*/
    class So_hex8fbarScatraType : public So_hex8fbarType
    {
     public:
      std::string Name() const override { return "So_hex8fbarScatraType"; }

      static So_hex8fbarScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex8fbarScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8FBARSCATRA"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Solid Scatra Element                              thon 12/15 |
     *----------------------------------------------------------------------*/
    class So_hex27ScatraType : public So_hex27Type
    {
     public:
      std::string Name() const override { return "So_hex27ScatraType"; }

      static So_hex27ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex27ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH27SCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                       |
     *----------------------------------------------------------------------*/
    class So_tet4ScatraType : public So_tet4Type
    {
     public:
      std::string Name() const override { return "So_tet4ScatraType"; }

      static So_tet4ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet4ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDT4SCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                       |
     *----------------------------------------------------------------------*/
    class So_tet10ScatraType : public So_tet10Type
    {
     public:
      std::string Name() const override { return "So_tet10ScatraType"; }

      static So_tet10ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet10ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDT10SCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  WEDGE 6 Element                                       |
     *----------------------------------------------------------------------*/
    class So_weg6ScatraType : public So_weg6Type
    {
     public:
      std::string Name() const override { return "So_weg6ScatraType"; }

      static So_weg6ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_weg6ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDW6SCATRA"; }
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
