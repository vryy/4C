/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi_serv.H

\brief Helper class for coupled simulations (INCA - BACI)
       --> aero_tfsi.H & aero_tfsi.cpp

\level 3

\maintainer Georg Hammerl

*/


/*----------------------------------------------------------------------*
 | definitions                                              ghamm 12/11 |
 *----------------------------------------------------------------------*/
#ifndef AERO_TFSI_SERV_H
#define AERO_TFSI_SERV_H


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>


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
  template <unsigned int n, unsigned int m> class Matrix;
}

namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}

namespace MORTAR
{
  class MortarInterface;
}

struct sendingdata {
   double x;
   double y;
   double z;
   double T;
};


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
      /// constructor for TFSI coupling
      AeroCouplingUtils
      (
        Teuchos::RCP<DRT::Discretization> structdis,  ///< structure discretization
        Teuchos::RCP<DRT::Discretization> thermodis  ///< thermo discretization
      );

      /// constructor for mechanical or thermal coupling
      AeroCouplingUtils
      (
        Teuchos::RCP<DRT::Discretization> dis,
        bool structure
      );

      /// useful parameters for mortar coupling
      void MortarParams
      (
        Teuchos::ParameterList& mortarparams
      );

      /// empty destructor
      virtual ~AeroCouplingUtils(){}

      /// mortar matrices to map data from INCA to Baci are set up
      void BuildMortarCoupling
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> idispn,  ///< structural interface displacement of the current time step
        std::vector<LINALG::Matrix<3,1> >& aerocoords,  ///< coordinates where the corresponding aero forces affect
        const int nodeoffset   ///< node offset for mortar from proc to proc
      );

      /// transfer load vector from fluid to structure
      void TransferFluidLoadsToStruct
      (
        int interf,   ///< interface id
        std::vector<LINALG::Matrix<4,1> >& aeroforces,   ///< loads given from INCA
        Teuchos::RCP<Epetra_Vector> iforce,   ///< interface load on structure
        Teuchos::RCP<Epetra_Vector> ithermoload   ///< interface load on thermal problem
      );

      /// transfer vector from slave to master
      Teuchos::RCP<Epetra_MultiVector> SlaveToMaster
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> sv  ///< incoming fluid interface values to be transferred
      ) const;

      /// data is gathered on proc 0 and packed in a specific way before transferring it to INCA
      void PackData
      (
        int interf,   ///< interface id
        Teuchos::RCP<const Epetra_Vector> idispnp,  ///< interface displacements
        Teuchos::RCP<const Epetra_Vector> itempnp,  ///< interface temperatures
        std::vector<double>& packeddata, ///< output: packed data, ready for sending
        bool writedata = true   ///< bool whether interface point data is written to file
      );

      /// data is gathered on proc 0 and packed in a specific way before transferring it to INCA using Immersed Boundary Method
      void TransferStructValuesToFluid
      (
        int interf,   ///< interface id
        Teuchos::RCP<const Epetra_Vector> idisp,   ///< interface displacements
        Teuchos::RCP<const Epetra_Vector> ivel,   ///< interface velocities
        Teuchos::RCP<const Epetra_Vector> itemp,   ///< interface temperatures
        std::vector<double>& packeddata,   ///< output: packed data, ready for sending
        bool writedata = true   ///< bool whether interface point data is written to file
      );

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
        Teuchos::RCP<const Epetra_Vector> fullvector  ///< state vector of the full structural field
      );

      /// extract thermal interface values for a given full field
      Teuchos::RCP<Epetra_Vector> ThrExtractInterfaceVal
      (
        int interf,   ///< interface id
        Teuchos::RCP<const Epetra_Vector> fullvector  ///< state vector of the full thermal field
      );

      /// apply interface tractions from the interface to a full field
      void StrInsertInterfaceVal
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> iforce,  ///< state vector of the structural interface
        Teuchos::RCP<Epetra_Vector> force  ///< state vector of the full structural field
      );

      /// apply interface heat flux from the interface to a full field
      void ThrInsertInterfaceVal
      (
        int interf,   ///< interface id
        Teuchos::RCP<Epetra_Vector> iforce,  ///< state vector of the thermal interface
        Teuchos::RCP<Epetra_Vector> force  ///< state vector of the full thermal field
      );

      /// length scaling factor for coupling is returned
      double LengthScaling() { return lengthscaling_; }

      /// time scaling factor for coupling is returned
      double TimeScaling() { return timescaling_; }

      /// returns whether mechanical coupling is used
      bool MechCoupling() { return mechanicalcoupling_; }

      /// returns whether thermal coupling is used
      bool ThermoCoupling() { return thermalcoupling_; }

    private:
      const int dim_;   ///< dimension of the problem
      const int myrank_;  ///< my processor id
      std::vector<Teuchos::RCP<DRT::Discretization> > istructdis_;  ///< discretization of structural interface
      std::vector<Teuchos::RCP<DRT::Discretization> > ithermodis_;  ///< discretization of thermal interface

      std::vector<Teuchos::RCP<MORTAR::MortarInterface> > mortarinterface_;   ///< mortar interfaces

      std::vector<Teuchos::RCP<const Epetra_Map> > istructdofredumap_;  ///< dof map of fully redundant structural interface
      std::vector<Teuchos::RCP<const Epetra_Map> > ithermodofredumap_;  ///< dof map of fully redundant structural interface

      std::vector<Teuchos::RCP<const LINALG::MapExtractor> > structrowmapext_;  ///< map extractor between struct discret and interface struct dis
      std::vector<Teuchos::RCP<const LINALG::MapExtractor> > thermorowmapext_;  ///< map extractor between thermo discret and interface thermo dis

      std::vector<Teuchos::RCP<const Epetra_Import> > istructproc0importer_;   ///< map extractor between distributed struct interface and data on proc 0
      std::vector<Teuchos::RCP<const Epetra_Import> > ithermoproc0importer_;   ///< map extractor between distributed thermo interface and data on proc 0

      std::vector<Teuchos::RCP<LINALG::SparseMatrix> > DinvM_;  ///< projection operator for mortar coupling

      const int dofpernode_; ///< number of dofs per mortar node
      const bool mechanicalcoupling_; ///< flag whethter mechanical coupling is activated
      const bool thermalcoupling_; ///< flag whether thermal coupling is activated

      int maxid_; ///< max id of the interfaces

      int mastermaxdof_; ///< max dof on TSI side

      int mastermaxnodeid_; ///< max node id on TSI side

      int mastermaxeleid_;   /// max ele id on TSI side

      // fluid uses [m] which becomes here e.g. 100.0 [cm] or 1000.0 [mm]
      double lengthscaling_; ///< scaling factor for extensions for the coupling if other than SI unit [m] is used

      // fluid uses [s] which becomes here e.g. 1000.0 [ms] or 1.0e6 [mikros]
      double timescaling_; ///< scaling factor for time for the coupling if other than SI unit [s] is used

    }; // class AeroCouplingUtils

  }
}

#endif
