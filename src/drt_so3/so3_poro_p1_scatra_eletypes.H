/*----------------------------------------------------------------------*/
/*! \file

 \brief element types of the 3D solid-poro element (p1, mixed approach) including scatra
 functionality

 \level 2

 *----------------------------------------------------------------------*/

#ifndef SO3_PORO_P1_SCATRA_ELETYPES_H_
#define SO3_PORO_P1_SCATRA_ELETYPES_H_

#include "so3_poro_p1_eletypes.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_hex8PoroP1ScatraType : public So_hex8PoroP1Type
    {
     public:
      std::string Name() const { return "So_hex8PoroP1ScatraType"; }

      static So_hex8PoroP1ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex8PoroP1ScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                         schmidt 09/17 |
     *----------------------------------------------------------------------*/
    class So_tet4PoroP1ScatraType : public So_tet4PoroP1Type
    {
     public:
      std::string Name() const { return "So_tet4PoroP1ScatraType"; }

      static So_tet4PoroP1ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet4PoroP1ScatraType instance_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SO3_PORO_P1_SCATRA_ELETYPES_H_ */
