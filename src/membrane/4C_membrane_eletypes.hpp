/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element Type

*----------------------------------------------------------------------*/
#ifndef FOUR_C_MEMBRANE_ELETYPES_HPP
#define FOUR_C_MEMBRANE_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_elementtype.hpp"
#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                          fbraeu 06/16 |
     *----------------------------------------------------------------------*/
    class MembraneTri3Type : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "Membrane_tri3Type"; }

      static MembraneTri3Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Discret::Discretization& dis) override { return 0; };

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static MembraneTri3Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 6 Element                                          fbraeu 06/16 |
     *----------------------------------------------------------------------*/
    class MembraneTri6Type : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "Membrane_tri6Type"; }

      static MembraneTri6Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Discret::Discretization& dis) override { return 0; };

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static MembraneTri6Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                         fbraeu 06/16 |
     *----------------------------------------------------------------------*/
    class MembraneQuad4Type : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "Membrane_quad4Type"; }

      static MembraneQuad4Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Discret::Discretization& dis) override { return 0; };

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static MembraneQuad4Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                         fbraeu 06/16 |
     *----------------------------------------------------------------------*/
    class MembraneQuad9Type : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "Membrane_quad9Type"; }

      static MembraneQuad9Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int Initialize(Discret::Discretization& dis) override { return 0; };

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static MembraneQuad9Type instance_;
    };

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
