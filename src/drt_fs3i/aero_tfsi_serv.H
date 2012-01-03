/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi_serv.cpp

\brief Helper class for coupled simulations (INCA - BACI)

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*/


#ifdef CCADISCRET

#ifndef AERO_TFSI_SERV_H
#define AERO_TFSI_SERV_H

#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_mapextractor.H"


namespace FS3I
{

  namespace UTILS
  {

    class AeroCouplingUtils
    {
    public:
      /// constructor initializing internal variables
        AeroCouplingUtils
      (
        Teuchos::RCP<DRT::Discretization> structdis,   ///< structure discretization
        Teuchos::RCP<DRT::Discretization> thermodis   ///< structure discretization
      );

      /// empty destructor
      virtual ~AeroCouplingUtils(){}

      /// forces from the aero code are projected on the structural surface and distributed to the nodes
      void ProjectForceOnStruct
      (
          Teuchos::RCP<Epetra_Vector> idispn,    ///< structural interface displacement of the current time step
          map<int, map<int, LINALG::Matrix<3,1> > >& aerocoords,   ///< coordinates where the corresponding aero forces affect
          map<int, map<int, LINALG::Matrix<4,1> > >& aeroforces,   ///< forces from the aero code
          Teuchos::RCP<Epetra_Vector> iforce,   ///< distributed interface force that loads the structure
          Teuchos::RCP<Epetra_Vector> iredutemp   ///< redundant distributed interface temperatures on the thermal interface
      );

      /// data is gathered on proc 0 and packed in a specific way before transferring it to INCA
      void PackData
      (
        Teuchos::RCP<Epetra_Vector> idispnp,
        Teuchos::RCP<Epetra_Vector> itempnp,
        vector<double>& packeddata
      );

      /// get material configuration
      Teuchos::RCP<Epetra_Vector> GetMatConfig(){ return matconfig_;}

      /// get redundant dof map of the thermal field
      Teuchos::RCP<const Epetra_Map> GetReduInterfThermoMap(){ return ithermodofredumap_;}

      /// get map extractor for the data exchange between the thermal field and the interface
      Teuchos::RCP<const LINALG::MapExtractor> GetInterfaceMapExtractor(){ return thermorowmapext_;}

      /// access function for the interface discretization
      RCP<DRT::Discretization> GetInterfaceStructDis(){ return istructdis_; }

      /// access function for the interface discretization
      RCP<DRT::Discretization> GetInterfaceThermoDis(){ return ithermodis_; }

      int NumInterfaces(){ return (maxid_+1); }

    protected:
      /// calculate current position of structure interface nodes
      std::map<int,LINALG::Matrix<3,1> > CurrentInterfacePos
      (
          Teuchos::RCP<Epetra_Vector> redudisp       ///< redundant version of structure displacements
      );

      /// fluid loads (forces & heat fluxes) are distributed to the adjacent nodes of the closest element 
      void DistributeForceToNodes
      (
        int eleid,
        int interface,
        std::map<int,LINALG::Matrix<3,1> >& currentpositions,
        LINALG::Matrix<3,1>& projpoint,
        LINALG::Matrix<4,1>& aeroforce,
        Teuchos::RCP<Epetra_Vector> distributedforces,
        Teuchos::RCP<Epetra_Vector> ithermoload
      );

      /// calculate the normal of an element
      vector<double> Normal
      (
        std::map<int,LINALG::Matrix<3,1> >& currentpositions,
        const int head1,
        const int origin,
        const int head2
      );


    private:

    Teuchos::RCP<DRT::Discretization> istructdis_;    ///< discretization of structural interface
    Teuchos::RCP<DRT::Discretization> ithermodis_;    ///< discretization of thermal interface

    Teuchos::RCP<const Epetra_Map> istructdofredumap_;   ///< dof map of fully redundant structural interface
    Teuchos::RCP<const Epetra_Map> ithermodofredumap_;   ///< dof map of fully redundant structural interface

    map<int, map<int, RCP<DRT::Element> > > structreduelements_;   ///< redundant version of coupling elements for each interface

    Teuchos::RCP<const LINALG::MapExtractor> thermorowmapext_;   ///< map extractor between thermo discret and interface thermo dis

    Teuchos::RCP<Epetra_Vector> matconfig_;    ///< material configuration of structural gids

    double maxmindist_;    ///< search radius to find close elements on the structural interface

    int maxid_;    ///< currently not needed

    }; // class AeroCouplingUtils

  }
}

#endif
#endif