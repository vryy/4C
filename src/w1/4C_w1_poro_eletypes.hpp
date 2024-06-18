/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element.

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_W1_PORO_ELETYPES_HPP
#define FOUR_C_W1_PORO_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_w1.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                       |
     *----------------------------------------------------------------------*/
    class WallQuad4PoroType : public Discret::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallQuad4PoroType"; }

      static WallQuad4PoroType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad4PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                       |
     *----------------------------------------------------------------------*/
    class WallQuad9PoroType : public Discret::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallQuad9PoroType"; }

      static WallQuad9PoroType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad9PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 4 Element                                       |
     *----------------------------------------------------------------------*/
    class WallNurbs4PoroType : public Discret::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallNurbs4PoroType"; }

      static WallNurbs4PoroType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static WallNurbs4PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  NURBS 9 Element                                       |
     *----------------------------------------------------------------------*/
    class WallNurbs9PoroType : public Discret::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallNurbs9PoroType"; }

      static WallNurbs9PoroType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static WallNurbs9PoroType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                       |
     *----------------------------------------------------------------------*/
    class WallTri3PoroType : public Discret::ELEMENTS::Wall1Type
    {
     public:
      std::string Name() const override { return "WallTri3PoroType"; }

      static WallTri3PoroType& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static WallTri3PoroType instance_;
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
