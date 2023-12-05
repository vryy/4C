/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element Type with ScaTra coupling

*----------------------------------------------------------------------*/
#ifndef BACI_MEMBRANE_SCATRA_ELETYPES_H
#define BACI_MEMBRANE_SCATRA_ELETYPES_H

#include "baci_membrane_eletypes.H"


namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  TRI 3 Element                                          sfuchs 05/18 |
     *----------------------------------------------------------------------*/
    class MembraneScatra_tri3Type : public Membrane_tri3Type
    {
     public:
      std::string Name() const override { return "MembraneScatra_tri3Type"; }

      static MembraneScatra_tri3Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static MembraneScatra_tri3Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  TRI 6 Element                                          sfuchs 05/18 |
     *----------------------------------------------------------------------*/
    class MembraneScatra_tri6Type : public Membrane_tri6Type
    {
     public:
      std::string Name() const override { return "MembraneScatra_tri6Type"; }

      static MembraneScatra_tri6Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static MembraneScatra_tri6Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 4 Element                                         sfuchs 05/18 |
     *----------------------------------------------------------------------*/
    class MembraneScatra_quad4Type : public Membrane_quad4Type
    {
     public:
      std::string Name() const override { return "MembraneScatra_quad4Type"; }

      static MembraneScatra_quad4Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static MembraneScatra_quad4Type instance_;
    };

    /*----------------------------------------------------------------------*
     |  QUAD 9 Element                                         sfuchs 05/18 |
     *----------------------------------------------------------------------*/
    class MembraneScatra_quad9Type : public Membrane_quad9Type
    {
     public:
      std::string Name() const override { return "MembraneScatra_quad9Type"; }

      static MembraneScatra_quad9Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static MembraneScatra_quad9Type instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
