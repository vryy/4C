
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
          RCP<DRT::Discretization> masterdis,
          RCP<DRT::Discretization> slavedis,
          ADAPTER::CouplingMortar& coupsf
      );

      /// empty destructor
      virtual ~SlideAleUtils(){}
      
      /// remesh ALE corresponding
      void Remeshing
      (
          Teuchos::RCP<Epetra_Vector> idisptotal,
          Teuchos::RCP<Epetra_Vector> idispstep,
          RCP<DRT::Discretization> masterdis,
          RCP<DRT::Discretization> slavedis,
          Teuchos::RCP<Epetra_Vector> islave,
          ADAPTER::CouplingMortar& coupsf,
          const Epetra_Comm& comm,
          INPAR::FSI::SlideALEProj aletype
      );
      

    protected:

      /// compute average interface displacement
      vector<double> Centerdisp
      (
          Teuchos::RCP<Epetra_Vector> idisptotal,
          Teuchos::RCP<Epetra_Vector>idispstep,
          RCP<DRT::Discretization> masterdis,
          const Epetra_Comm& comm
      );

      /// Position of master nodes
      std::map<int,LINALG::Matrix<3,1> >  CurrentMasterPos
      (
          Teuchos::RCP<Epetra_Vector> idisptotal,
          Teuchos::RCP<Epetra_Vector> reddisp,
          DRT::Discretization& interfacedis,
          RCP<Epetra_Map> msfullelemap
      );
      
      /// project ALE  nodes to he structure surface
      void SlideProjection
      (
          Teuchos::RCP<Epetra_Vector> islavestep,
          DRT::Discretization& interfacedis,
          vector<double> centerdisp,
          std::map<int, RCP<DRT::Element> > masterreduelements,
          std::map<int,LINALG::Matrix<3,1> > currentpositions,
          Teuchos::RCP<Epetra_Vector> islave,
          ADAPTER::CouplingMortar& coupsf,
          RCP<DRT::Discretization> slavedis,
          INPAR::FSI::SlideALEProj aletype
      );
      
    private:  
    
      Teuchos::RCP<Epetra_Vector> idispms_;       ///<displacement of master and slave interface for mortar
      vector<double> centerdisptotal_;            ///<sum over all center displacements     
//      Teuchos::RCP<Epetra_Vector> islave_;        ///<displacement of slave side of the interface
      
      map<int, DRT::Node*> imastergnodes_;        ///< masternodes in the interface including ghosts
      map<int, RCP<DRT::Element> > imastereles_;  ///< masterlelements in the interface
      map<int, DRT::Node*> islaveslidnodes_;          ///< sliding slavenodes (= fluid) in the interface
      map<int, DRT::Node*> islaveconfnodes_;          ///< non-sliding slavenodes (= fluid) in the interface
      
    }; // class SlideAleUtils
    
  }
}

#endif
#endif
