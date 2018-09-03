/*!----------------------------------------------------------------------
\file so3_thermo_eletypes.H

<pre>
   Maintainer: Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>

*----------------------------------------------------------------------*/
#ifndef SO3_THERMO_ELETYPES_H_
#define SO3_THERMO_ELETYPES_H_

#include "so_hex8.H"
#include "so_hex8fbar.H"
#include "so_hex27.H"
#include "so_hex20.H"
#include "so_tet4.H"
#include "so_tet10.H"
#include "so_nurbs27.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     * HEX8 element
     *----------------------------------------------------------------------*/
    class So_hex8ThermoType : public So_hex8Type
    {
     public:
      std::string Name() const { return "So_hex8ThermoType"; }

      static So_hex8ThermoType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex8ThermoType instance_;

    };  // class So_hex8ThermoType


    /*----------------------------------------------------------------------*
     * HEX8FBAR element
     *----------------------------------------------------------------------*/
    class So_hex8fbarThermoType : public So_hex8fbarType
    {
     public:
      std::string Name() const { return "So_hex8fbarThermoType"; }

      static So_hex8fbarThermoType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex8fbarThermoType instance_;

    };  // class So_hex8fbarThermoType


    /*----------------------------------------------------------------------------*
     * TET4 Element
     *----------------------------------------------------------------------------*/
    class So_tet4ThermoType : public So_tet4Type
    {
     public:
      std::string Name() const { return "So_tet4ThermoType"; }

      static So_tet4ThermoType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet4ThermoType instance_;

    };  // class So_tet4ThermoType

    /*----------------------------------------------------------------------------*
     * TET10 Element
     *----------------------------------------------------------------------------*/
    class So_tet10ThermoType : public So_tet10Type
    {
     public:
      std::string Name() const { return "So_tet10ThermoType"; }

      static So_tet10ThermoType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet10ThermoType instance_;

    };  // class So_tet10ThermoType

    /*----------------------------------------------------------------------*
     *  HEX 27 Element
     *----------------------------------------------------------------------*/
    class So_hex27ThermoType : public So_hex27Type
    {
     public:
      std::string Name() const { return "So_hex27ThermoType"; }

      static So_hex27ThermoType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex27ThermoType instance_;

    };  // class So_hex27ThermoType

    /*----------------------------------------------------------------------*
     *  HEX 20 Element
     *----------------------------------------------------------------------*/
    class So_hex20ThermoType : public So_hex20Type
    {
     public:
      std::string Name() const { return "So_hex20ThermoType"; }

      static So_hex20ThermoType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex20ThermoType instance_;

    };  // class So_hex20ThermoType

    /*----------------------------------------------------------------------*
     *  NURBS 27 Element
     *----------------------------------------------------------------------*/
    class So_nurbs27ThermoType : public NURBS::So_nurbs27Type
    {
     public:
      std::string Name() const { return "So_nurbs27ThermoType"; }

      static So_nurbs27ThermoType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_nurbs27ThermoType instance_;

    };  // class So_hex20ThermoType

  }  // namespace ELEMENTS

}  // namespace DRT


/*----------------------------------------------------------------------*/
#endif  // SO3_THERMO_ELETYPES_H_
