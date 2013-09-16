/*----------------------------------------------------------------------*/
/*!
 \file fpsi_utils.H

 \brief  Utility Methods For Fluid Porous Structure Interaction Problems

 <pre>
   Maintainer: Rauch Andreas
               rauch.andreas@mytum.de
 </pre>
 *-----------------------------------------------------------------------*/
#ifndef FPSI_UTILS_H_
#define FPSI_UTILS_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include <Epetra_Time.h>
#include <Epetra_Vector.h>
#include "Teuchos_RCP.hpp"
#include "../drt_inpar/inpar_fpsi.H"
#include "fpsi.H"
#include "fpsi_monolithic.H"
#include "fpsi_monolithic_plain.H"
#include "../drt_adapter/adapter_algorithmbase.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
namespace FPSI
{

  class UTILS
  {
  public:

  //! Singleton access method
  static Teuchos::RCP<UTILS> Instance();

  static Teuchos::RCP<FPSI::UTILS> instance_;

  Teuchos::RCP<FPSI::FPSI_Base> SetupDiscretizations(const Epetra_Comm& comm,
                                                     const Teuchos::ParameterList& fpsidynparams,
                                                     const Teuchos::ParameterList& poroelastdynparams);

  void RedistributeInterface(Teuchos::RCP<DRT::Discretization> masterdis, const DRT::Discretization& slavedis, const std::string& condname, std::map<int,int>& interfacefacingelementmap);

  //! Fills a map that matches the global id of an interface element on the slave side to the global id of the opposing bulk element.
  //! This is done processor locally. Works only for matching grids.
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
  void SetupLocalInterfaceFacingElementMap(DRT::Discretization& masterdis, const DRT::Discretization& slavedis, const std::string& condname, std::map<int,int>& interfacefacingelementmap);

  Teuchos::RCP<std::map<int,int> > Fluid_PoroFluid_InterfaceMap;
  Teuchos::RCP<std::map<int,int> > PoroFluid_Fluid_InterfaceMap;


  Teuchos::RCP<std::map<int,int> > Get_Fluid_PoroFluid_InterfaceMap(){return Fluid_PoroFluid_InterfaceMap;};
  Teuchos::RCP<std::map<int,int> > Get_PoroFluid_Fluid_InterfaceMap(){return PoroFluid_Fluid_InterfaceMap;};

  }; // class UTILS
} // namespace FPSI

#endif /* FPSI_UTILS_H_ */
