/*----------------------------------------------------------------------------*/
/*! \file

\brief Utilities for FSI problems

\level 1

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef FSI_UTILS_H
#define FSI_UTILS_H

#include <set>
#include <string>
#include <functional>
#include <iostream>

#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_fsi.H"

// forward declarations
namespace ADAPTER
{
  class CouplingMortar;
  class FSIStructureWrapper;
}  // namespace ADAPTER

namespace ALE
{
  class Ale;
}

namespace IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
}  // namespace IO

namespace FSI
{
  namespace UTILS
  {
    /// write interface jacobian to file
    void DumpJacobian(NOX::Epetra::Interface::Required& interface, double alpha, double beta,
        Teuchos::RCP<Epetra_Vector> soln, std::string filename);

    /// shift given map to have no overlap with maps in vector
    Teuchos::RCP<const Epetra_Map> ShiftMap(Teuchos::RCP<const Epetra_Map> emap,
        const std::vector<Teuchos::RCP<const Epetra_Map>>& vecSpaces);


    /// Check whether fluid node numbers and ALE node numbers are equal.
    bool FluidAleNodesDisjoint(
        Teuchos::RCP<DRT::Discretization> fluiddis,  ///< pointer to fluid discretization
        Teuchos::RCP<DRT::Discretization> aledis     ///< pointer to ALE discretization
    );

    /*!
    \brief implementation of sliding ALE stuff
     */
    class SlideAleUtils
    {
     public:
      /// constructor initializing internal variables
      SlideAleUtils(Teuchos::RCP<DRT::Discretization> structdis,  ///< structure discretization
          Teuchos::RCP<DRT::Discretization> fluiddis,             ///< fluid discretization
          ::ADAPTER::CouplingMortar& coupsf,                      ///< mortar adapter
          bool structcoupmaster,            ///< is structure master of adapter coupling?
          INPAR::FSI::SlideALEProj aleproj  ///< projection enum
      );

      /// empty destructor
      virtual ~SlideAleUtils() {}

      /// remesh ALE corresponding
      void Remeshing(::ADAPTER::FSIStructureWrapper& structure,  ///< structure adapter
          Teuchos::RCP<DRT::Discretization> fluiddis,            ///< fluid discretization
          Teuchos::RCP<Epetra_Vector> idispale,      ///< standard ALE interface displacement
          Teuchos::RCP<Epetra_Vector> iprojdispale,  ///< projected ALE interface displacement
          ::ADAPTER::CouplingMortar& coupsf,         ///< mortar adapter
          const Epetra_Comm& comm                    ///< communicator
      );

      /// Compute new coupling matrices D and M
      void EvaluateMortar(Teuchos::RCP<Epetra_Vector> idispstruct,  ///< displacement of structure
          Teuchos::RCP<Epetra_Vector> idispale,                     ///< (proj.) displacement of ale
          ::ADAPTER::CouplingMortar& coupsf                         ///< mortar adapter
      );

      /// Compute new coupling matrices D and M
      void EvaluateFluidMortar(Teuchos::RCP<Epetra_Vector> idispma,  ///< displacement of structure
          Teuchos::RCP<Epetra_Vector> idispsl  ///< (proj.) displacement of ale
      );

      /// use fluid-fluid mortar interface to interpolate between fluid quantities before and after
      /// sliding
      Teuchos::RCP<Epetra_Vector> InterpolateFluid(Teuchos::RCP<const Epetra_Vector>
              uold  ///< fluid velocity in configuration before sliding
      );

      /// write history vectors for restart
      void OutputRestart(IO::DiscretizationWriter& output);

      /// read history values for restart
      void ReadRestart(IO::DiscretizationReader& reader);

     protected:
      /// compute average interface displacement
      std::vector<double> Centerdisp(
          ::ADAPTER::FSIStructureWrapper& structure,  ///< structure adapter
          const Epetra_Comm& comm                     ///< communicator
      );

      /// compute approximate interface rotation (structuresplit)
      void Rotation(DRT::Discretization& mtrdis,  ///< mtr interface  discretization
          Teuchos::RCP<Epetra_Vector> idispale,   ///< vector of ALE displacements
          const Epetra_Comm& comm,                ///< communicator
          std::map<int, double>& rotrat,          ///< rotation ratio of tangential displacements
          Teuchos::RCP<Epetra_Vector>
              rotfull  ///< vector of full displacements in tangential directions
      );


      /// calculate current position of structure interface nodes
      std::map<int, LINALG::Matrix<3, 1>> CurrentStructPos(
          Teuchos::RCP<Epetra_Vector> reddisp,  ///< redundant version of structure displacements
          DRT::Discretization& interfacedis,    ///< interface discretization
          std::map<int, double>& maxcoord);


      /// project ALE nodes onto the structure surface
      void SlideProjection(::ADAPTER::FSIStructureWrapper& structure,  ///< structure adapter
          Teuchos::RCP<DRT::Discretization> fluiddis,                  ///< fluid discretization
          Teuchos::RCP<Epetra_Vector> idispale,      ///< standard ALE interface displacement
          Teuchos::RCP<Epetra_Vector> iprojdispale,  ///< projected ALE interface displacement
          ::ADAPTER::CouplingMortar& coupsf,         ///< mortar adapter
          const Epetra_Comm& comm                    ///< communicator
      );

      /// Build full redundant structure and fluid elements.
      /// Necessary for search-trees since MORTAR elements do not know about their facets and edges.
      /// Furthermore, this function builds StructuralSurface elements from the fluid outer surface
      /// for rotation.
      void RedundantElements(::ADAPTER::CouplingMortar& coupsf, const Epetra_Comm& comm);

     private:
      const INPAR::FSI::SlideALEProj aletype_;
      Teuchos::RCP<Epetra_Vector>
          idispms_;  ///< merged vector of displacements (struct and fluid interface)
      std::vector<double> centerdisptotal_;  ///< sum over all center displacement increments
      double maxmindist_;                    ///< maximal distance between fluidpairs

      //      std::map<int, Teuchos::RCP<DRT::Element> > istructslideles_;  ///< sliding struct
      //      elements in the interface
      std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>>
          istructslideles_;  ///< sliding struct elements in the interface
      std::map<int, DRT::Node*>
          istructdispnodes_;  ///< struct nodes in the interface used for centerdisp calculation
      std::map<int, Teuchos::RCP<DRT::Element>>
          istructdispeles_;  ///< struct elements in the interface used for centerdisp calc
      //      std::map<int, Teuchos::RCP<Epetra_Map> >  slideeleredmap_;      ///< redundant version
      //      of sliding struct elements
      std::map<int, std::map<int, DRT::Node*>>
          ifluidslidnodes_;                        ///< sliding fluid nodes in the interface
      std::map<int, DRT::Node*> ifluidconfnodes_;  ///< sticking fluid nodes in the interface

      std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>> ifluidslideles_;
      std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>> ifluidslidstructeles_;

      std::map<int, std::map<int, Teuchos::RCP<DRT::Element>>> structreduelements_;

      Teuchos::RCP<const Epetra_Map> structdofrowmap_;
      Teuchos::RCP<const Epetra_Map> fluiddofrowmap_;
      Teuchos::RCP<Epetra_Map> structfullnodemap_;
      Teuchos::RCP<Epetra_Map> structfullelemap_;
      Teuchos::RCP<Epetra_Map> fluidfullnodemap_;
      Teuchos::RCP<Epetra_Map> fluidfullelemap_;

      int maxid_;

      Teuchos::RCP<Epetra_Vector> iprojhist_;  ///< history of final displacements

      bool structcoupmaster_;  ///< is structure master of coupling?

      /// coupling of fluid before and fluid after the sliding
      Teuchos::RCP<::ADAPTER::CouplingMortar> coupff_;

    };  // class SlideAleUtils

  }  // namespace UTILS
}  // namespace FSI

#endif  // FSI_UTILS_H
