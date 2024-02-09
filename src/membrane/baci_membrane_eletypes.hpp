/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element Type

*----------------------------------------------------------------------*/
#ifndef BACI_MEMBRANE_ELETYPES_HPP
#define BACI_MEMBRANE_ELETYPES_HPP

#include "baci_config.hpp"

#include "baci_lib_elementtype.hpp"
#include "baci_linalg_serialdensematrix.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                          fbraeu 06/16 |
     *----------------------------------------------------------------------*/
    class Membrane_tri3Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Membrane_tri3Type"; }

      static Membrane_tri3Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override { return 0; };

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Membrane_tri3Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 6 Element                                          fbraeu 06/16 |
     *----------------------------------------------------------------------*/
    class Membrane_tri6Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Membrane_tri6Type"; }

      static Membrane_tri6Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override { return 0; };

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Membrane_tri6Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                         fbraeu 06/16 |
     *----------------------------------------------------------------------*/
    class Membrane_quad4Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Membrane_quad4Type"; }

      static Membrane_quad4Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override { return 0; };

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Membrane_quad4Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                         fbraeu 06/16 |
     *----------------------------------------------------------------------*/
    class Membrane_quad9Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "Membrane_quad9Type"; }

      static Membrane_quad9Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      int Initialize(DRT::Discretization& dis) override { return 0; };

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Membrane_quad9Type instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // MEMBRANE_ELETYPES_H
