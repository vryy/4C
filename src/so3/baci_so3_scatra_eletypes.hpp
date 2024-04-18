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
    class SoHex8ScatraType : public SoHex8Type
    {
     public:
      std::string Name() const override { return "So_hex8ScatraType"; }

      static SoHex8ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8SCATRA"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 8 fbar Element                                        Thon 12/14 |
     *----------------------------------------------------------------------*/
    class SoHex8fbarScatraType : public SoHex8fbarType
    {
     public:
      std::string Name() const override { return "So_hex8fbarScatraType"; }

      static SoHex8fbarScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8fbarScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8FBARSCATRA"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Solid Scatra Element                              thon 12/15 |
     *----------------------------------------------------------------------*/
    class SoHex27ScatraType : public SoHex27Type
    {
     public:
      std::string Name() const override { return "So_hex27ScatraType"; }

      static SoHex27ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoHex27ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH27SCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                       |
     *----------------------------------------------------------------------*/
    class SoTet4ScatraType : public SoTet4Type
    {
     public:
      std::string Name() const override { return "So_tet4ScatraType"; }

      static SoTet4ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoTet4ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDT4SCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                       |
     *----------------------------------------------------------------------*/
    class SoTet10ScatraType : public SoTet10Type
    {
     public:
      std::string Name() const override { return "So_tet10ScatraType"; }

      static SoTet10ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoTet10ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDT10SCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  WEDGE 6 Element                                       |
     *----------------------------------------------------------------------*/
    class SoWeg6ScatraType : public SoWeg6Type
    {
     public:
      std::string Name() const override { return "So_weg6ScatraType"; }

      static SoWeg6ScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoWeg6ScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDW6SCATRA"; }
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
