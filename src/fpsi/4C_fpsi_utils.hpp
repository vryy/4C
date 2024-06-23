/*----------------------------------------------------------------------*/
/*! \file
  \brief  Utility Methods For Fluid Porous Structure Interaction Problems
\level 3

 *-----------------------------------------------------------------------*/
#ifndef FOUR_C_FPSI_UTILS_HPP
#define FOUR_C_FPSI_UTILS_HPP

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_linalg_mapextractor.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
namespace FPSI
{
  class FpsiBase;

  class Utils
  {
   public:
    //! Singleton access method
    static Teuchos::RCP<Utils> Instance();

    //! singleton object
    static Teuchos::RCP<FPSI::Utils> instance_;

    //! Setup Discretizations for FPSI problem (clone ALE and porofluid and setup interfaces)
    Teuchos::RCP<FPSI::FpsiBase> setup_discretizations(const Epetra_Comm& comm,
        const Teuchos::ParameterList& fpsidynparams,
        const Teuchos::ParameterList& poroelastdynparams);

    //! redistribute interface for parallel computations
    void redistribute_interface(Teuchos::RCP<Core::FE::Discretization> masterdis,
        Teuchos::RCP<const Core::FE::Discretization> slavedis, const std::string& condname,
        std::map<int, int>& interfacefacingelementmap);

    //! build map for fpsi interface
    void SetupInterfaceMap(const Epetra_Comm& comm,
        Teuchos::RCP<Core::FE::Discretization> structdis,
        Teuchos::RCP<Core::FE::Discretization> porofluiddis,
        Teuchos::RCP<Core::FE::Discretization> fluiddis,
        Teuchos::RCP<Core::FE::Discretization> aledis);

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
    void setup_local_interface_facing_element_map(Core::FE::Discretization& masterdis,
        const Core::FE::Discretization& slavedis, const std::string& condname,
        std::map<int, int>& interfacefacingelementmap);

    //! access methods
    Teuchos::RCP<std::map<int, int>> get_fluid_poro_fluid_interface_map()
    {
      return fluid_poro_fluid_interface_map_;
    };
    Teuchos::RCP<std::map<int, int>> get_poro_fluid_fluid_interface_map()
    {
      return poro_fluid_fluid_interface_map_;
    };

   private:
    //! interface maps
    Teuchos::RCP<std::map<int, int>> fluid_poro_fluid_interface_map_;
    Teuchos::RCP<std::map<int, int>> poro_fluid_fluid_interface_map_;

  };  // class Utils


  namespace UTILS
  {
    /// specific MultiMapExtractor to handle the fluid field
    class MapExtractor : public Core::LinAlg::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_fsi = 1,
        cond_fpsi = 2
      };

      /// setup the whole thing
      void setup(
          const Core::FE::Discretization& dis, bool withpressure = false, bool overlapping = false);

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
      void setup(Teuchos::RCP<const Epetra_Map>& additionalothermap,
          const FPSI::UTILS::MapExtractor& extractor);

      /// get all element gids those nodes are touched by any condition
      Teuchos::RCP<std::set<int>> conditioned_element_map(
          const Core::FE::Discretization& dis) const;

      MAP_EXTRACTOR_VECTOR_METHODS(other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(fsi_cond, cond_fsi)
      MAP_EXTRACTOR_VECTOR_METHODS(fpsi_cond, cond_fpsi)
    };

  }  // namespace UTILS
}  // namespace FPSI

FOUR_C_NAMESPACE_CLOSE

#endif
