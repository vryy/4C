/*----------------------------------------------------------------------*/
/*! \file
  \brief  Utility Methods For Fluid Porous Structure Interaction Problems
\level 3

\maintainer  Johannes Kremheller
 *-----------------------------------------------------------------------*/
#ifndef FPSI_UTILS_H_
#define FPSI_UTILS_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include <Epetra_Time.h>
#include <Epetra_Vector.h>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../linalg/linalg_mapextractor.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
namespace FPSI
{
  class FPSI_Base;

  class Utils
  {
   public:
    //! Singleton access method
    static Teuchos::RCP<Utils> Instance();

    //! singleton object
    static Teuchos::RCP<FPSI::Utils> instance_;

    //! Setup Discretizations for FPSI problem (clone ALE and porofluid and setup interfaces)
    Teuchos::RCP<FPSI::FPSI_Base> SetupDiscretizations(const Epetra_Comm& comm,
        const Teuchos::ParameterList& fpsidynparams,
        const Teuchos::ParameterList& poroelastdynparams);

    //! redistribute interface for parallel computations
    void RedistributeInterface(Teuchos::RCP<DRT::Discretization> masterdis,
        Teuchos::RCP<const DRT::Discretization> slavedis, const std::string& condname,
        std::map<int, int>& interfacefacingelementmap);

    //! build map for fpsi interface
    void SetupInterfaceMap(const Epetra_Comm& comm, Teuchos::RCP<DRT::Discretization> structdis,
        Teuchos::RCP<DRT::Discretization> porofluiddis, Teuchos::RCP<DRT::Discretization> fluiddis,
        Teuchos::RCP<DRT::Discretization> aledis);

    //! Fills a map that matches the global id of an interface element on the slave side to the
    //! global id of the opposing bulk element. This is done processor locally. Works only for
    //! matching grids.
    /*!
      \param In
             masterdis - Reference of discretization of master field.
      \param In
             slavedis  - Reference of discretization of slave field.
      \param In
             condname  - String with name of condition on interface to be considered.
      \param Out
             interfacefacingelementmap - processor local map to be filled

       See Detailed Description section for further discussion.
    */
    void SetupLocalInterfaceFacingElementMap(DRT::Discretization& masterdis,
        const DRT::Discretization& slavedis, const std::string& condname,
        std::map<int, int>& interfacefacingelementmap);

    //! access methods
    Teuchos::RCP<std::map<int, int>> Get_Fluid_PoroFluid_InterfaceMap()
    {
      return Fluid_PoroFluid_InterfaceMap;
    };
    Teuchos::RCP<std::map<int, int>> Get_PoroFluid_Fluid_InterfaceMap()
    {
      return PoroFluid_Fluid_InterfaceMap;
    };

   private:
    //! interface maps
    Teuchos::RCP<std::map<int, int>> Fluid_PoroFluid_InterfaceMap;
    Teuchos::RCP<std::map<int, int>> PoroFluid_Fluid_InterfaceMap;

  };  // class Utils


  namespace UTILS
  {
    /// specific MultiMapExtractor to handle the fluid field
    class MapExtractor : public LINALG::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_fsi = 1,
        cond_fpsi = 2
      };

      /// setup the whole thing
      void Setup(
          const DRT::Discretization& dis, bool withpressure = false, bool overlapping = false);

      /*!
       * \brief setup from an existing extractor
       * By calling this setup version we create a map extractor from
       * (1) an existing map extractor and
       * (2) a DOF-map from another discretization, which is appended to othermap.
       * We need this in the context of XFFSI.
       * \param (in) additionalothermap : map of additional unconditioned DOF
       * \param (in) extractor : extractor, from which the conditions are cloned
       * \author kruse
       * \date 05/2014
       */
      void Setup(Teuchos::RCP<const Epetra_Map>& additionalothermap,
          const FPSI::UTILS::MapExtractor& extractor);

      /// get all element gids those nodes are touched by any condition
      Teuchos::RCP<std::set<int>> ConditionedElementMap(const DRT::Discretization& dis) const;

      MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(FSICond, cond_fsi)
      MAP_EXTRACTOR_VECTOR_METHODS(FPSICond, cond_fpsi)
    };

  }  // namespace UTILS
}  // namespace FPSI

#endif /* FPSI_UTILS_H_ */
