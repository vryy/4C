/*------------------------------------------------------------------------------------------------*/
/*!
\file adapter_fluid_combust.H

\brief Fluid field adapter

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_COMBUST_H
#define ADAPTER_FLUID_COMBUST_H

#include "../drt_combust/combust_fluidimplicitintegration.H"
#include "adapter_fluid.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

namespace DRT
{
  class ResultTest;
}

namespace COMBUST
{
  class FlameFront;
  class InterfaceHandleCombust;
}

namespace ADAPTER
{
  /// adapter to normal fluid algorithm
  /*!
    <h3> title? </h3>

    Documentation of class FluidCombust is still missing!

    \author henke
    \date 08/08
  */
  class FluidCombust : public Fluid
  {
  public:
    //! constructor
    explicit FluidCombust(
             Teuchos::RCP<DRT::Discretization> dis,
             Teuchos::RCP<LINALG::Solver> solver,
             Teuchos::RCP<ParameterList> params,
             Teuchos::RCP<IO::DiscretizationWriter> output);

    /*--------------------------------------------------------------------------------------------*
     * group of member functions which are not implemented or used
     *
     * The idea is to keep all member functions which are not used in a combustion problem, but
     * inherited from class "Fluid", in the header file. Thus adapter_fluid_combust.cpp only
     * contains functions which can be called
     *
     * These are the unused member functions.                                            henke 08/08
     *--------------------------------------------------------------------------------------------*/

    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() { dserror("not implemented");
                                                               return fluid_.InitialGuess(); };
    virtual Teuchos::RCP<const Epetra_Vector> RHS()          { dserror("not implemented");
                                                               return fluid_.Residual(); };
    virtual Teuchos::RCP<const Epetra_Vector> Velaf()        { dserror("not implemented");
                                                               return null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accam()        { dserror("not implemented");
                                                               return null; };
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp()       { dserror("not implemented");
                                                               return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() { dserror("not implemented");
                                                         const Epetra_Map* dofrowmap = fluiddis_->DofRowMap();
                                                         return Teuchos::rcp(dofrowmap, false); };
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix(){ dserror("not implemented");
    // if anything (e.g. monolithic FSI) we give fluid coupling and interface DOF combined back
                                                               return fluid_.SystemMatrix(); };
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix(){ dserror("not implemented");
                                                                             return Teuchos::null; };
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives(){ dserror("not implemented");
                                                                          return Teuchos::null; };
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor(){ dserror("not implemented");
                                                                           return Teuchos::null;}; //return fluid_.DirichMaps();
    virtual const LINALG::MapExtractor& FreeSurface() const { dserror("not implemented!");
            // this is a hack to compile the code! danger of dangling reference!
            LINALG::MapExtractor object;
            LINALG::MapExtractor& reftoobject=object;
            return reftoobject; }

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin()  { dserror("not implemented!");
                                                                         return Teuchos::null; };
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity(){ dserror("not implemented!");
                                                                         return Teuchos::null; };

    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp){ return; };
    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc){dserror("not implemented!");};
    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)      { return; };
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel,
                                          Teuchos::RCP<Epetra_Vector> iforce){ return; };

    // the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm){ return; };

    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx){ return; };
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx){ return; };

    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)  { dserror("not implemented!"); return; };
    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)  { dserror("not implemented!"); return; };

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(){
            dserror("not implemented!");
            return interface_.ExtractFSICondVector(fluid_.IntegrateInterfaceShape("FSICoupling")); };

    virtual void UseBlockMatrix(bool splitmatrix)  { dserror("not implemented!"); return; };

    /// linear fluid solve with just a interface load
    /*!
      The very special solve done in steepest descent relaxation
      calculation (and matrix free Newton Krylov).

      \note Can only be called after a valid fluid solve.
    */
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel){
            dserror("not implemented!");
            const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
            Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
            interface_.InsertFSICondVector(ivel,relax);
            //fluid_.LinearRelaxationSolve(relax);
            return ExtractInterfaceForces(); };

    ///  lift'n'drag forces, statistics time sample and output of solution and statistics
    virtual void StatisticsAndOutput(){ dserror("not implemented!");
                                        return; };

    ///  set fields for low-Mach-number flow within outer iteration loop
    virtual void SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,
                                   RCP<const Epetra_Vector> scalaram,
                                   RCP<const Epetra_Vector> scalardtam,
                                   const double             thermpressaf,
                                   const double             thermpressam,
                                   const double             thermpressdtaf,
                                   const double             thermpressdtam,
                                   Teuchos::RCP<DRT::Discretization> scatradis)
    { dserror("not implemented!");return; };

    /// set fields for low-Mach-number flow at end of time step
    virtual void SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,
                                   const double             thermpressnp,
                                   RCP<const Epetra_Vector> scatraresidual,
                                   Teuchos::RCP<DRT::Discretization> scatradis,
                                   const int                whichscalar = -1)
    { dserror("not implemented!");return; };

    /// set velocity field (separate computation)
    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp)
    { dserror("not implemented!");return; };

    /// provide access to the turbulence statistic manager
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
    {
      return fluid_.TurbulenceStatisticManager();
    }

    /*--------------------------------------------------------------------------------------------*
     * group member functions which are implemented and used
     *
     * These are the used member functions                                               henke 08/08
     *------------------------------------------------------------------------------------------- */

    virtual Teuchos::RCP<const Epetra_Vector> TrueResidual();
    virtual Teuchos::RCP<const Epetra_Vector> Veln();
    virtual Teuchos::RCP<const Epetra_Vector> Velnp();
    virtual Teuchos::RCP<const Epetra_Vector> Velnm(){ dserror("not implemented"); return null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accnp(){ dserror("not implemented"); return null; };
    virtual Teuchos::RCP<const Epetra_Vector> Accn(){ dserror("not implemented"); return null; };
    virtual Teuchos::RCP<const Epetra_Vector> Scaaf() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Scaam() { dserror("not implemented"); return Teuchos::null; };
    virtual Teuchos::RCP<const Epetra_Vector> Hist();
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel();
    virtual Teuchos::RCP<const Epetra_Vector> FsVel(){ dserror("not implemented"); return Teuchos::null; };

    /// set initial flow field
    virtual void SetInitialFlowField(const INPAR::FLUID::InitialField initfield, const int startfuncno);
    virtual void SetInitialFlowField(const INPAR::COMBUST::InitialField initfield, const int initfuncno);
    virtual Teuchos::RCP<DRT::Discretization> Discretization() { return fluiddis_; }
    virtual Teuchos::RCP<const DRT::DofSet> DofSet() { return fluid_.StandardDofSet(); }
    virtual void TimeLoop();
    virtual void PrepareTimeStep();
    virtual void PrepareTimeandStep() { dserror("not implemented!"); return; }
    virtual void ClearTimeInt();
    virtual void PrepareSolve() { dserror("not implemented!"); return; }
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc);
    virtual bool ConvergenceCheck(int itnum, int itmax, const double ittol) { dserror("not implemented!"); return false; }
    virtual void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment) { dserror("not implemented!"); return; }
    virtual void Update();
    virtual void Output();
    virtual void StatisticsOutput() { dserror("not implemented!"); return; }
    virtual IO::DiscretizationWriter& DiscWriter() { return fluid_.DiscWriter();}

    virtual FLD::UTILS::FluidXFluidMapExtractor XFluidFluidMapExtractor() {dserror("not implemented");
           // this is a hack to compile the code! danger of dangling reference!
            FLD::UTILS::FluidXFluidMapExtractor  object;
            FLD::UTILS::FluidXFluidMapExtractor& reftoobject=object;
            return reftoobject;}

    virtual Teuchos::RCP<const Epetra_Vector> Stepinc() {dserror("not implemented");return Teuchos::null;}

    /// Return dirichlet boundray condition (used for fluidfluid)
    virtual Teuchos::RCP<const Epetra_Map> FluidDirichMaps(){dserror("not implemented"); return Teuchos::null;}

    virtual void NonlinearSolve();
    /// predictor
    virtual void Predictor();
    /// (multiple) corrector
    virtual void MultiCorrector();

    /// get the linear solver object used for this field
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    /// the interface map contains all velocity dofs at the interface
    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap() { return Teuchos::null; } // only used for FSI
    virtual Teuchos::RCP<const Epetra_Map> VelocityRowMap();
    virtual Teuchos::RCP<const Epetra_Map> PressureRowMap();

    /// scaling factor needed to convert the residual to real forces
    virtual double ResidualScaling() const;

    virtual double TimeScaling() const;

    /// communication object at the interface
    virtual const FLD::UTILS::MapExtractor& Interface() const { return interface_; }

    virtual void ReadRestart(int step);
    virtual void SetRestart(const int step,
                    const double time,
                    Teuchos::RCP<const Epetra_Vector> readvelnp,
                    Teuchos::RCP<const Epetra_Vector> readveln,
                    Teuchos::RCP<const Epetra_Vector> readvelnm,
                    Teuchos::RCP<const Epetra_Vector> readaccnp,
                    Teuchos::RCP<const Epetra_Vector> readaccn) { dserror("not implemented!"); return; }
    virtual double Time() const;
    virtual int Step() const;
    virtual double Dt() const;
    virtual int  Itemax() const;
    virtual void SetItemax(int itemax);

    /// return type of time integration scheme
    virtual const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const { return fluid_.TimIntScheme(); }

    virtual void LiftDrag();

    virtual const double EvalTime() const { return fluid_.EvalTime(); }
    virtual void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph) { fluid_.Redistribute(nodegraph); return; }
    virtual void TransferVectorsToNewDistribution(const Teuchos::RCP<COMBUST::InterfaceHandleCombust> interfacehandle) { fluid_.TransferVectorsToNewDistribution(interfacehandle); return; }

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();
    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln()  { dserror("not implemented!");
                                                                         return Teuchos::null; };
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

    /// extract the velocity part of a vector of unknowns
    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres);

    virtual void SetTopOptPorosityField(Teuchos::RCP<Epetra_Vector> porosity)
    { dserror("not implemented!"); };

    /// combustion fluid imports level set discretization (that is, a ScaTra Discretization!)
    virtual void ImportDiscretization(Teuchos::RCP<DRT::Discretization> importdis);

    /// combustion fluid imports geometrical interface information from combustion algorithm
    virtual void ImportInterface(
        const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle,
        const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle_old);

    ///
    virtual void ImportFlameFront(const Teuchos::RCP<COMBUST::FlameFront>& flamefront);

    virtual Teuchos::RCP<DRT::ResultTest>  CreateFieldTest();

  private:

    /// implicit fluid time integration scheme
    FLD::CombustFluidImplicitTimeInt             fluid_;

    /// the interface map setup for interface <-> full translation
    FLD::UTILS::MapExtractor               interface_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      fluiddis_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    //@}

    //! @name state vectors using the solid boundary parallel distribution (row vectors)
//    Teuchos::RCP<Epetra_Vector>            itrueresnp_;///< interface residual      n+1
//    Teuchos::RCP<Epetra_Vector>            ivelnp_;    ///< interface velocity      n+1
//    Teuchos::RCP<Epetra_Vector>            iveln_;     ///< interface velocity      n
//    Teuchos::RCP<Epetra_Vector>            ivelnm_;    ///< interface velocity      n-1
//    Teuchos::RCP<Epetra_Vector>            iaccn_;     ///< interface acceleration  n
    //@}

    /// all velocity dofs not at the interface
//    Teuchos::RCP<Epetra_Map>       innervelmap_;

    /// pointer to the level set discretization
    // this is the link to the interface information contained in the level set discretization
    Teuchos::RCP<DRT::Discretization> gfuncdis_;

  };

} // namespace ADAPTER

#endif
#endif
