/*----------------------------------------------------------------------------*/
/*! \file
\brief Element types of the 2D solid-poro element (p1/mixed approach).

\level 2


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_W1_PORO_P1_ELETYPES_HPP
#define FOUR_C_W1_PORO_P1_ELETYPES_HPP

#include "4C_config.hpp"

#include "4C_w1_poro_eletypes.hpp"

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
     |  QUAD 4 Element                                                      |
     *----------------------------------------------------------------------*/
    class WallQuad4PoroP1Type : public Discret::ELEMENTS::WallQuad4PoroType
    {
     public:
      std::string name() const override { return "WallQuad4PoroP1Type"; }

      static WallQuad4PoroP1Type& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad4PoroP1Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                                      |
     *----------------------------------------------------------------------*/
    class WallQuad9PoroP1Type : public Discret::ELEMENTS::WallQuad9PoroType
    {
     public:
      std::string name() const override { return "WallQuad9PoroP1Type"; }

      static WallQuad9PoroP1Type& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static WallQuad9PoroP1Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                                       |
     *----------------------------------------------------------------------*/
    class WallTri3PoroP1Type : public Discret::ELEMENTS::WallTri3PoroType
    {
     public:
      std::string name() const override { return "WallTri3PoroP1Type"; }

      static WallTri3PoroP1Type& instance();

      Core::Communication::ParObject* create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> create(const int id, const int owner) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix compute_null_space(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      int initialize(Core::FE::Discretization& dis) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static WallTri3PoroP1Type instance_;
    };

  }  // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif
