/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the p1 (mixed) solid-poro element

 \level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_SO3_PORO_P1_ELETYPES_HPP
#define FOUR_C_SO3_PORO_P1_ELETYPES_HPP

#include "baci_config.hpp"

#include "baci_so3_poro_eletypes.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                                       |
     *----------------------------------------------------------------------*/
    class So_hex8PoroP1Type : public So_hex8PoroType
    {
     public:
      std::string Name() const override { return "So_hex8PoroP1Type"; }

      static So_hex8PoroP1Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_hex8PoroP1Type instance_;

      std::string GetElementTypeString() const { return "SOLIDH8POROP1"; }
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                                       |
     *----------------------------------------------------------------------*/
    class So_tet4PoroP1Type : public So_tet4PoroType
    {
     public:
      std::string Name() const override { return "So_tet4PoroP1Type"; }

      static So_tet4PoroP1Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static So_tet4PoroP1Type instance_;

      std::string GetElementTypeString() const { return "SOLIDT4POROP1"; }
    };

  }  // namespace ELEMENTS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif  // SO3_PORO_P1_ELETYPES_H
