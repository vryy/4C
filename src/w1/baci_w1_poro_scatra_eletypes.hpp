/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element including scatra functionality.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_W1_PORO_SCATRA_ELETYPES_HPP
#define FOUR_C_W1_PORO_SCATRA_ELETYPES_HPP

#include "baci_config.hpp"

#include "baci_w1_poro.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallQuad4PoroScatraType : public DRT::ELEMENTS::WallQuad4PoroType
    {
     public:
      std::string Name() const override { return "WallQuad4PoroScatraType"; }

      static WallQuad4PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad4PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                        schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallQuad9PoroScatraType : public DRT::ELEMENTS::WallQuad9PoroType
    {
     public:
      std::string Name() const override { return "WallQuad9PoroScatraType"; }

      static WallQuad9PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad9PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 4 Element                                       schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallNurbs4PoroScatraType : public DRT::ELEMENTS::WallNurbs4PoroType
    {
     public:
      std::string Name() const override { return "WallNurbs4PoroScatraType"; }

      static WallNurbs4PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallNurbs4PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 9 Element                                       schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallNurbs9PoroScatraType : public DRT::ELEMENTS::WallNurbs9PoroType
    {
     public:
      std::string Name() const override { return "WallNurbs9PoroScatraType"; }

      static WallNurbs9PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallNurbs9PoroScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class WallTri3PoroScatraType : public DRT::ELEMENTS::WallTri3PoroType
    {
     public:
      std::string Name() const override { return "WallTri3PoroScatraType"; }

      static WallTri3PoroScatraType& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static WallTri3PoroScatraType instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // W1_PORO_SCATRA_ELETYPES_H
