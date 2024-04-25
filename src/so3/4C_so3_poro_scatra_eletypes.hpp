/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the 3D solid-poro element including scatra functionality

 \level 2

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_SO3_PORO_SCATRA_ELETYPES_HPP
#define FOUR_C_SO3_PORO_SCATRA_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_so3_poro_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class SoHex8PoroScatraType : public SoHex8PoroType
    {
     public:
      std::string Name() const override { return "So_hex8PoroScatraType"; }

      static SoHex8PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoHex8PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH8POROSCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class SoTet4PoroScatraType : public SoTet4PoroType
    {
     public:
      std::string Name() const override { return "So_tet4PoroScatraType"; }

      static SoTet4PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoTet4PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDT4POROSCATRA"; }
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class SoHex27PoroScatraType : public SoHex27PoroType
    {
     public:
      std::string Name() const override { return "So_hex27PoroScatraType"; }

      static SoHex27PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoHex27PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDH27POROSCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class SoTet10PoroScatraType : public SoTet10PoroType
    {
     public:
      std::string Name() const override { return "So_tet10PoroScatraType"; }

      static SoTet10PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoTet10PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SOLIDT10POROSCATRA"; }
    };

    /*----------------------------------------------------------------------*
     |  NURBS 27 Element                                      schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class SoNurbs27PoroScatraType : public SoNurbs27PoroType
    {
     public:
      std::string Name() const override { return "So_nurbs27PoroScatraType"; }

      static SoNurbs27PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static SoNurbs27PoroScatraType instance_;

      std::string GetElementTypeString() const { return "SONURBS27POROSCATRA"; }
    };


  }  // namespace ELEMENTS
}  // namespace DRT
FOUR_C_NAMESPACE_CLOSE

#endif
