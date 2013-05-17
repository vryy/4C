/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi_serv.cpp

\brief Helper class for coupled simulations (INCA - BACI)
       --> aero_tfsi.H & aero_tfsi.cpp

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                              ghamm 12/11 |
 *----------------------------------------------------------------------*/
#ifndef AERO_TFSI_SERV_H
#define AERO_TFSI_SERV_H


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include <map>
#include <vector>
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_mortar.H"


/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 12/11 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
  class Element;
}

namespace LINALG
{
  class MapExtractor;
  class SparseMatrix;
}

namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}


namespace FS3I
{
  namespace UTILS
  {
    /// Helper class for coupled simulations INCA-BACI,
    /// Builds interface discretizations (thermal & structural) for
    /// the coupling interface.
    /// Contains methods to project fluid loads onto the structure
    /// and to prepare data before sending to INCA
    /*!
      \author ghamm
      \date 12/11
     */
    class AeroCouplingUtils
    {
    public:
      /// constructor initializing internal variables
        AeroCouplingUtils
      (
        Teuchos::RCP<DRT::Discretization> structdis,  ///< structure discretization
        Teuchos::RCP<DRT::Discretization> thermodis  ///< thermo discretization
      );

      /// empty destructor
      virtual ~AeroCouplingUtils(){}

      /// forces from the aero code are projected on the structural surface and distributed to the nodes
      void ProjectForceOnStruct
      (
        int interface,   ///< interface id
        Teuchos::RCP<Epetra_Vector> idispn,  ///< structural interface displacement of the current time step
        std::map<int, LINALG::Matrix<3,1> >& aerocoords,  ///< coordinates where the corresponding aero forces affect
        std::map<int, LINALG::Matrix<4,1> >& aeroforces,  ///< forces from the aero code
        Teuchos::RCP<Epetra_Vector> iforce,  ///< distributed interface force that loads the structure
        Teuchos::RCP<Epetra_Vector> iredutemp  ///< redundant distributed interface temperatures on the thermal interface
      );

      /// mortar matrices to map data from INCA to Baci are set up
      void BuildMortarCoupling
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> idispn,  ///< structural interface displacement of the current time step
        std::map<int, LINALG::Matrix<3,1> >& aerocoords  ///< coordinates where the corresponding aero forces affect
      );

      /// returns type of shape functions used for mortar
      INPAR::MORTAR::ShapeFcn Shapefcn() {return shapefcn_;}

      /// loads from the aero code are transferred via mortar to the structure (for dual shape functions in mortar)
      void TransferFluidLoadsToStructDual
      (
        int interf,   ///< interface id
        std::map<int, LINALG::Matrix<4,1> >& aeroforces,  ///< forces from the aero code
        Teuchos::RCP<Epetra_Vector> iforce,  ///< distributed interface force that loads the structure
        Teuchos::RCP<Epetra_Vector> iredutemp  ///< redundant distributed interface temperatures on the thermal interface
      );

      /// loads from the aero code are transferred via mortar to the structure (for standard shape functions in mortar)
      void TransferFluidLoadsToStructStd
      (
        int interf,   ///< interface id
        std::map<int, LINALG::Matrix<4,1> >& aeroforces,  ///< forces from the aero code
        Teuchos::RCP<Epetra_Vector> iredutemp  ///< redundant distributed interface temperatures on the thermal interface
      );

      Teuchos::RCP<Epetra_Vector> StructToFluid
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> mv  ///< incoming structural/thermal interface values to be transferred
      ) const;

      Teuchos::RCP<Epetra_Vector> FluidToStruct
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> sv  ///< incoming fluid interface values to be transferred
      ) const;

      /// data is gathered on proc 0 and packed in a specific way before transferring it to INCA
      void PackData
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> idispnp,  ///< interface displacements
        Teuchos::RCP<Epetra_Vector> itempnp,  ///< interface temperatures
        std::vector<double>& packeddata, ///< output: packed data, ready for sending
        bool writedata = true   ///< bool whether interface point data is written to file
      );

      /// data is gathered on proc 0 and packed in a specific way before transferring it to INCA (for mortarcoupling)
      void TransferStructValuesToFluidDual
      (
        int interf,   ///< interface id
        int sizeofphysicaldomain,    ///< size of physical domain
        std::map<int, LINALG::Matrix<3,1> >& aerocoords,  ///< interface displacements of fluid
        Teuchos::RCP<Epetra_Vector> itempnp,  ///< interface temperatures
        std::vector<double>& packeddata, ///< output: packed data, ready for sending
        bool writedata = true   ///< bool whether interface point data is written to file
      );

      /// data is gathered on proc 0 and packed in a specific way before transferring it to INCA (for mortarcoupling)
      void TransferStructValuesToFluidStd
      (
        int interf,   ///< interface id
        int sizeofphysicaldomain,    ///< size of physical domain
        std::map<int, LINALG::Matrix<3,1> >& aerocoords,  ///< interface displacements of fluid
        Teuchos::RCP<Epetra_Vector> itempnp,  ///< interface temperatures
        std::vector<double>& packeddata, ///< output: packed data, ready for sending
        bool writedata = true   ///< bool whether interface point data is written to file
      );

      /// get material configuration
//      Teuchos::RCP<Epetra_Vector> GetMatConfig(){ return matconfig_;}

      /// access function for the ith interface discretization
      Teuchos::RCP<DRT::Discretization> GetInterfaceStructDis(int i){ return istructdis_[i]; }

      /// access function for the ith interface discretization
      Teuchos::RCP<DRT::Discretization> GetInterfaceThermoDis(int i){ return ithermodis_[i]; }

      /// get number of coupling interfaces
      int NumInterfaces(){ return (maxid_+1); }

      /// extract structural interface values for a given full field
      Teuchos::RCP<Epetra_Vector> StrExtractInterfaceVal
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> fullvector  ///< state vector of the full structural field
      );

      /// extract thermal interface values for a given full field
      Teuchos::RCP<Epetra_Vector> ThrExtractInterfaceVal
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> fullvector  ///< state vector of the full thermal field
      );

      /// apply interface tractions from the interface to a full field
      void StrApplyInterfaceVal
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> iforce,  ///< state vector of the structural interface
        Teuchos::RCP<Epetra_Vector> force  ///< state vector of the full structural field
      );

      /// apply interface heat flux from the interface to a full field
      void ThrApplyInterfaceVal
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> iforce,  ///< state vector of the thermal interface
        Teuchos::RCP<Epetra_Vector> force  ///< state vector of the full thermal field
      );


    protected:
      /// calculate current position of structural interface nodes
      std::map<int,LINALG::Matrix<3,1> > CurrentInterfacePos
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> redudisp  ///< redundant version of structure displacements
      );

      /// fluid loads (forces & heat fluxes) are distributed to the adjacent nodes of the closest element 
      void DistributeForceToNodes
      (
        int eleid,  ///< element id of the closest element to the fluid load
        int interf,  /// interface id in which the eleid is contained
        std::map<int,LINALG::Matrix<3,1> >& currentpositions,  ///< current position of structural interface nodes
        LINALG::Matrix<3,1>& projpoint,  ///< projection point on the closest element
        LINALG::Matrix<4,1>& aeroforce,  ///< fluid induced loads (3x forces + 1x heat flux)
        Teuchos::RCP<Epetra_Vector> distributedforces,  ///< distributed forces to the nodes of the element
        Teuchos::RCP<Epetra_Vector> ithermoload  ///< distributed heat fluxes to the nodes of the element
      );

    private:
      std::vector<Teuchos::RCP<DRT::Discretization> > istructdis_;  ///< discretization of structural interface
      std::vector<Teuchos::RCP<DRT::Discretization> > ithermodis_;  ///< discretization of thermal interface

      std::vector<Teuchos::RCP<const Epetra_Map> > istructdofredumap_;  ///< dof map of fully redundant structural interface
      std::vector<Teuchos::RCP<const Epetra_Map> > ithermodofredumap_;  ///< dof map of fully redundant structural interface

      std::vector<std::map<int, Teuchos::RCP<DRT::Element> > > structreduelements_; ///< redundant version of coupling elements for each interface

      std::vector<Teuchos::RCP<const LINALG::MapExtractor> > structrowmapext_;  ///< map extractor between struct discret and interface struct dis
      std::vector<Teuchos::RCP<const LINALG::MapExtractor> > thermorowmapext_;  ///< map extractor between thermo discret and interface thermo dis

      std::vector<Teuchos::RCP<LINALG::SparseMatrix> > DinvM_;  ///< projection operator for mortar coupling

      std::vector<Teuchos::RCP<const Epetra_Map> > serialslrownoderestr_;  ///< restricted row node map after mortar coupling

      // the next two are here for testing standard shape functions
      std::vector<Teuchos::RCP<LINALG::SparseMatrix> > D_;  ///< M matrix of mortar coupling
      std::vector<Teuchos::RCP<LINALG::SparseMatrix> > M_;  ///< M matrix of mortar coupling

      INPAR::MORTAR::ShapeFcn shapefcn_;      // employed type of shape function set

//      Teuchos::RCP<Epetra_Vector> matconfig_; ///< material configuration of structural gids

      double maxmindist_; ///< search radius to find close elements on the structural interface

      int maxid_; ///< max id of the interfaces

      int mastermaxdof_; ///< max dof on TSI side

      int mastermaxnodeid_; ///< max node id on TSI side

    }; // class AeroCouplingUtils

  }
}

#endif
