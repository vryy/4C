/*----------------------------------------------------------------------*/
/*!
 \file poro_base.H

 \brief  Basis of all porous media algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/
#ifndef PORO_BASE_H_
#define PORO_BASE_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include <Epetra_Time.h>
#include <Epetra_Vector.h>

#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Condition;
  class Discretization;
}

namespace LINALG
{
  class MapExtractor;
}

namespace ADAPTER
{
  class FluidPoro;
  class Coupling;
  class FSIStructureWrapper;
}

namespace POROELAST
{
  /// poro algorithm base
  class PoroBase: public ADAPTER::AlgorithmBase
  {
    public:

      //! create using a Epetra_Comm
      explicit PoroBase(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

      //! virtual destructor to support polymorph destruction
      virtual ~PoroBase();

      /// read restart data
      virtual void ReadRestart(int restart);

      //! outer level time loop
      virtual void TimeLoop();

      /// initialise system
      virtual void SetupSystem() = 0;

      //! perform result tests
      virtual void TestResults(const Epetra_Comm& comm);

      //! access to structural field
      const Teuchos::RCP<ADAPTER::FSIStructureWrapper>& StructureField(){return structure_;}

      //! access to fluid field
      ::ADAPTER::FluidPoro& FluidField(){return *fluid_;}

      //! solve time step (depending on algorithm)
      virtual void Solve() = 0;

      //! @name Time loop building blocks

      //! start a new time step
      virtual void PrepareTimeStep();

      //! take current results for converged and save for next time step
      virtual void Update();

      //! calculate stresses, strains, energies
      virtual void PrepareOutput();

      //! output
      virtual void Output();

      //@}

      /// HACK
      void CalculateSurfPoro(const string& condstring);

      //! export node-based multivector to column map and add it to a parameter
      //list --> copy from scatra!!! move this to drt_lib eventually
      void AddMultiVectorToParameterList(
         Teuchos::ParameterList&          p,    //!< parameter list
         const std::string                name, //!< naming string for added vector
         Teuchos::RCP<Epetra_MultiVector> vec,   //!< Epetra_MultiVector to be added
         Teuchos::RCP<DRT::Discretization> discret //!< discretization
        );

    protected:

      //! @name Transfer helpers

      //! field transform
      Teuchos::RCP<Epetra_Vector> StructureToFluidField(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      //! field transform
      Teuchos::RCP<Epetra_Vector> FluidToStructureField(Teuchos::RCP<
          const Epetra_Vector> iv) const;

      //@}

      //! @name Transfer methods

      //! set fluid solution on structure
      void SetFluidSolution();

      //! set structure solution on fluid
      void SetStructSolution();

      //@}

      //! build map containing dofs with no penetration condition (fluid)
      void BuidNoPenetrationMap();

      //! return coupling object
      ADAPTER::Coupling& FluidStructureCoupling() { return *coupfs_; }

      //! Extractor used for constraint structure
      Teuchos::RCP< LINALG::MapExtractor > consplitter_;

      //! map of gids of fluid  DOFs and respective structure DOFs
      Teuchos::RCP<std::map<int, int> > nopenetrationmap_;

      //! Extractor used for no penetration condition
      Teuchos::RCP< LINALG::MapExtractor > nopenetration_;

      //! Extractor used for no penetration condition
      Teuchos::RCP< LINALG::MapExtractor > psiextractor_;

      //! vector containing no penetration - conditions
      std::vector<DRT::Condition*> nopencond_;

      //! vector containing global IDs of dofs with no penetration condition
      Teuchos::RCP<std::set<int> >condIDs_;

      //Teuchos::RCP<Epetra_Vector> fluidveln_; //!< global fluid velocities and pressures

     // Teuchos::RCP<Epetra_Vector> veln_; //!< global structure velocities
      //!< \f${V}_{n+1}\f$
      //!< at \f$t_{n+1}\f$
      //Teuchos::RCP<Epetra_Vector> dispn_; //!< global structure displacements
      //!< \f${d}_{n+1}\f$
      //!< at \f$t_{n+1}\f$

    private:
      //! check for submeshes and build of subnode and subelement map if necessary
      bool BuildSubMaps();

      //! coupling of fluid and structure (whole field)
      Teuchos::RCP<ADAPTER::Coupling> coupfs_;

      Teuchos::RCP<Epetra_Map> subnodemap_;

      Teuchos::RCP<Epetra_Map> subelemap_;

      //! @name Underlying fields

      //! underlying structure of the FSI problem
      Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    structure_;

      //! underlying fluid of the FSI problem
      Teuchos::RCP< ::ADAPTER::FluidPoro>        fluid_;

      //@}

      //!flag indicating not fully overlapping fluid and structure discretization
      bool submeshes_;

  }; // MonolithicBase
}

#endif /* PORO_BASE_H_ */
