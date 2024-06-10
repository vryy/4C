/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element Type with ScaTra coupling

*----------------------------------------------------------------------*/
#ifndef FOUR_C_MEMBRANE_SCATRA_ELETYPES_HPP
#define FOUR_C_MEMBRANE_SCATRA_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_membrane_eletypes.hpp"

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
     |  TRI 3 Element                                          sfuchs 05/18 |
     *----------------------------------------------------------------------*/
    class MembraneScatraTri3Type : public MembraneTri3Type
    {
     public:
      std::string Name() const override { return "MembraneScatra_tri3Type"; }

      static MembraneScatraTri3Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static MembraneScatraTri3Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 6 Element                                          sfuchs 05/18 |
     *----------------------------------------------------------------------*/
    class MembraneScatraTri6Type : public MembraneTri6Type
    {
     public:
      std::string Name() const override { return "MembraneScatra_tri6Type"; }

      static MembraneScatraTri6Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static MembraneScatraTri6Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                         sfuchs 05/18 |
     *----------------------------------------------------------------------*/
    class MembraneScatraQuad4Type : public MembraneQuad4Type
    {
     public:
      std::string Name() const override { return "MembraneScatra_quad4Type"; }

      static MembraneScatraQuad4Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static MembraneScatraQuad4Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                         sfuchs 05/18 |
     *----------------------------------------------------------------------*/
    class MembraneScatraQuad9Type : public MembraneQuad9Type
    {
     public:
      std::string Name() const override { return "MembraneScatra_quad9Type"; }

      static MembraneScatraQuad9Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static MembraneScatraQuad9Type instance_;
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
