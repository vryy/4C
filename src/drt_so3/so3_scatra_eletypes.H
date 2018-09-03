/*!----------------------------------------------------------------------
\file so3_scatra_eletypes.H

\brief Solid-scatra elements types

\level 2

<pre>
   \maintainer Thon Moritz
               thon@mhpc.mw.tum.de
               089 - 289-10264
</pre>

*----------------------------------------------------------------------*/
#ifndef SO3_SCATRA_ELETYPES_H_
#define SO3_SCATRA_ELETYPES_H_

#include "so_hex8.H"
#include "so_hex8fbar.H"
#include "so_hex27.H"
#include "so_tet4.H"
#include "so_tet10.H"
#include "so_weg6.H"
// #include "so_nstet5.H"

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*----------------------------------------------------------------------*
     |  HEX 8 Element                                       |
     *----------------------------------------------------------------------*/
    class So_hex8ScatraType : public So_hex8Type
    {
     public:
      std::string Name() const { return "So_hex8ScatraType"; }

      static So_hex8ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex8ScatraType instance_;
    };


    /*----------------------------------------------------------------------*
     |  HEX 8 fbar Element                                        Thon 12/14 |
     *----------------------------------------------------------------------*/
    class So_hex8fbarScatraType : public So_hex8fbarType
    {
     public:
      std::string Name() const { return "So_hex8fbarScatraType"; }

      static So_hex8fbarScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex8fbarScatraType instance_;
    };


    /*----------------------------------------------------------------------*
     |  HEX 27 Solid Scatra Element                              thon 12/15 |
     *----------------------------------------------------------------------*/
    class So_hex27ScatraType : public So_hex27Type
    {
     public:
      std::string Name() const { return "So_hex27ScatraType"; }

      static So_hex27ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_hex27ScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TET 4 Element                                       |
     *----------------------------------------------------------------------*/
    class So_tet4ScatraType : public So_tet4Type
    {
     public:
      std::string Name() const { return "So_tet4ScatraType"; }

      static So_tet4ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet4ScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  TET 10 Element                                       |
     *----------------------------------------------------------------------*/
    class So_tet10ScatraType : public So_tet10Type
    {
     public:
      std::string Name() const { return "So_tet10ScatraType"; }

      static So_tet10ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_tet10ScatraType instance_;
    };

    /*----------------------------------------------------------------------*
     |  WEDGE 6 Element                                       |
     *----------------------------------------------------------------------*/
    class So_weg6ScatraType : public So_weg6Type
    {
     public:
      std::string Name() const { return "So_weg6ScatraType"; }

      static So_weg6ScatraType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static So_weg6ScatraType instance_;
    };



    /*----------------------------------------------------------------------*
     |  NSTET 5 Element                                       |
     *----------------------------------------------------------------------*/

    /*

    class NStet5ScatraType : public NStet5Type
    {
    public:

      std::string Name() const { return "NStet5ScatraType"; }

      static NStet5ScatraType & Instance() { return instance_; };

      virtual DRT::ParObject* Create( const std::vector<char> & data );

      virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                                 const std::string eledistype,
                                                 const int id,
                                                 const int owner );

      virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

      virtual int Initialize(DRT::Discretization& dis);

      virtual void SetupElementDefinition(
    std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

    private:

      static NStet5ScatraType instance_;
    };
    */
  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SO3_SCATRA_ELETYPES_H_ */
