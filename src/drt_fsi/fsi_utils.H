
#ifdef CCADISCRET

#ifndef FSI_UTILS_H
#define FSI_UTILS_H

#include <set>
#include <string>
#include <functional>
#include <iostream>

#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_inpar/inpar_fsi.H"

//forward declarations
namespace ADAPTER
{
  class CouplingMortar;
  class FluidMovingBoundary;
  class Structure;
  class Fluid;
}

namespace ALE
{
  class Ale;
}

namespace IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
}

namespace FSI
{
  /// FSI::UTILS: Random stuff that might be helpful when dealing with FSI problems
  namespace UTILS
  {

    /// write interface jacobian to file
    void DumpJacobian(NOX::Epetra::Interface::Required& interface,
                      double alpha,
                      double beta,
                      Teuchos::RefCountPtr<Epetra_Vector> soln,
                      std::string filename);

    /// shift given map to have no overlap with maps in vector
    Teuchos::RCP<const Epetra_Map> ShiftMap(Teuchos::RCP<const Epetra_Map> emap,
                                            const std::vector<Teuchos::RCP<const Epetra_Map> >& vecSpaces);

    /// create ale discretization parallel to the fluid one
    //void CreateAleDiscretization();



    /*!
    \brief implementation of special clone strategy for automatic generation
           of an ALE-mesh from a given fluid discretization
     */
    class AleFluidCloneStrategy
    {
    public:

      /// constructor
      explicit AleFluidCloneStrategy(){}
      /// destructor
      virtual ~AleFluidCloneStrategy() {}

    protected:

      /// determine element type string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element* actele,
          const bool ismyele,
          vector<string>& eletype);

      /// set element-specific data (material etc.)
      void SetElementData(
          RCP<DRT::Element> newele,
          DRT::Element* oldele,
          const int matid,
          const bool nurbsdis);

      /// returns conditions names to be copied (source and target name)
      map<string,string> ConditionsToCopy();

      /// check for correct material
      void CheckMaterialType(const int matid);

    private:

    }; // class AleFluidCloneStrategy

    
    /*!
    \brief implementation of sliding ALE stuff
     */
    class SlideAleUtils
    {
    public:
      
      /// constructor initializing internal variables
      SlideAleUtils
      (
          RCP<DRT::Discretization> structdis,   ///< structure discretization
          RCP<DRT::Discretization> fluiddis,    ///< fluid discretization
          ADAPTER::CouplingMortar& coupsf,      ///< mortar adapter
          bool structcoupmaster                 ///< is structure master of adapter coupling?
      );

      /// empty destructor
      virtual ~SlideAleUtils(){}
      
      /// remesh ALE corresponding
      void Remeshing
      (
          ADAPTER::Structure& structure,          ///< structure adapter
          RCP<DRT::Discretization> fluiddis,      ///< fluid discretization
          Teuchos::RCP<Epetra_Vector> idispale,   ///< standard ALE interface displacement
          Teuchos::RCP<Epetra_Vector> iprojdispale,     ///< projected ALE interface displacement
          ADAPTER::CouplingMortar& coupsf,        ///< mortar adapter
          const Epetra_Comm& comm,                ///< communicator
          INPAR::FSI::SlideALEProj aletype        ///< projection enum
      );

      /// Compute new coupling matrices D and M
      void EvaluateMortar
      (
          Teuchos::RCP<Epetra_Vector> idispstruct, ///< displacement of structure
          Teuchos::RCP<Epetra_Vector> idispale,    ///< (proj.) displacement of ale
          ADAPTER::CouplingMortar& coupsf          ///< mortar adapter
      );
      
      /// write history vectors for restart
      void OutputRestart
      (
        IO::DiscretizationWriter& output
      );

      /// read history values for restart
      void ReadRestart
      (
          IO::DiscretizationReader& reader
      );

    protected:

      /// compute average interface displacement
      vector<double> Centerdisp
      (
          Teuchos::RCP<Epetra_Vector> idispn,      ///< structure displacement t_n
          Teuchos::RCP<Epetra_Vector> idispstep,   ///< structure displacement increment
          RCP<DRT::Discretization> structdis,      ///< discretization
          const Epetra_Comm& comm                  ///< communicator
      );

      /// calculate current position of structure interface nodes
      std::map<int,LINALG::Matrix<3,1> >  CurrentStructPos
      (
          Teuchos::RCP<Epetra_Vector> reddisp,     ///< redundant version of structure displacements
          DRT::Discretization& interfacedis,       ///< interface discretization
          RCP<Epetra_Map> msfullelemap             ///< merged map
      );
      
      /// project ALE nodes onto the structure surface
      void SlideProjection
      (
          DRT::Discretization& interfacedis,
          vector<double> centerdisp,
          std::map<int, RCP<DRT::Element> > structreduelements,
          std::map<int,LINALG::Matrix<3,1> > currentpositions,
          Teuchos::RCP<Epetra_Vector> iprojdispale,             ///< projected ALE interface displacement
          ADAPTER::CouplingMortar& coupsf,
          RCP<DRT::Discretization> fluiddis,
          INPAR::FSI::SlideALEProj aletype
      );
      
    private:  
    
      Teuchos::RCP<Epetra_Vector> idispms_;       ///< merged vector of displacements (struct and fluid interface)
      vector<double> centerdisptotal_;            ///< sum over all center displacement increments
      
      map<int, DRT::Node*> istructslidnodes_;        ///< sliding struct nodes in the interface
      map<int, RCP<DRT::Element> > istructslideles_;  ///< sliding struct elements in the interface
      map<int, DRT::Node*> istructdispnodes_;        ///< struct nodes in the interface used for centerdisp calculation
      map<int, RCP<DRT::Element> > istructdispeles_;  ///< struct elements in the interface used for centerdisp calc
      Teuchos::RCP<Epetra_Map>  slideeleredmap_;      ///< redundant version of sliding struct elements
      map<int, DRT::Node*> ifluidslidnodes_;          ///< sliding fluid nodes in the interface
      map<int, DRT::Node*> ifluidconfnodes_;          ///< sticking fluid nodes in the interface
      
      Teuchos::RCP<Epetra_Vector> iprojhist_;         ///< history of final displacements
      
      bool structcoupmaster_; ///< is structure master of coupling?
      
    }; // class SlideAleUtils
    
  }
}

#endif
#endif
