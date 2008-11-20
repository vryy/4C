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

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

#include "../drt_combust/combust_fluidimplicitintegration.H"
#include "../drt_combust/combust_fluidresulttest.H"
#include "adapter_fluid.H"

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
    virtual Teuchos::RCP<const Epetra_Vector> Velnp()        { dserror("not implemented");
                                                               return fluid_.Velnp(); };
    virtual Teuchos::RCP<const Epetra_Vector> Velaf()        { dserror("not implemented");
                                                               return null; };
    virtual Teuchos::RCP<const Epetra_Vector> Veln()         { dserror("not implemented");
                                                               return fluid_.Veln(); };
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp()       { dserror("not implemented");
    // CombustFluidImplicitTimeInt does not have a member Dispnp()!
                                                               //return fluid_.Dispnp();
                                                               return null; };
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() { dserror("not implemented");
                                                         const Epetra_Map* dofrowmap = dis_->DofRowMap();
                                                         return Teuchos::rcp(dofrowmap, false); };
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix(){ dserror("not implemented");
    // if anything (e.g. monolithic FSI) we give fluid coupling and interface DOF combined back
                                                               return fluid_.SystemMatrix(); };
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { dserror("not implemented");
                                                                              return Teuchos::null; };
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> MeshMoveMatrix(){ dserror("not implemented");
                                                                          return Teuchos::null; };
    virtual const LINALG::MapExtractor& FreeSurface() const {
            dserror("no free surface in combustion problems!");
            // this is a hack to compile the code! danger of dangling reference!
            LINALG::MapExtractor object;
            LINALG::MapExtractor& reftoobject=object;
            return reftoobject; }

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForcesRobin()  { dserror("no Robin coupling in combustion problems");
                                                                         return Teuchos::null; };
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceFluidVelocity(){ dserror("not implemented for combustion problems");
                                                                         return Teuchos::null; };

    virtual void ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp){ return; };
    virtual void ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel)      { return; };
    virtual void ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel,
                                          Teuchos::RCP<Epetra_Vector> iforce){ return; };

    // the mesh map contains all velocity dofs that are covered by an ALE node
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm){ return; };

    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx){ return; };
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx){ return; };

    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres){
            dserror("not implemented!");
            return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres); };

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(){
            dserror("not implemented!");
            return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling")); };

    virtual void UseBlockMatrix(const LINALG::MultiMapExtractor& domainmaps,
                                const LINALG::MultiMapExtractor& rangemaps,
                                bool splitmatrix){ return; };

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
            interface_.InsertCondVector(ivel,relax);
            //fluid_.LinearRelaxationSolve(relax);
            return ExtractInterfaceForces(); };

    ///  lift'n'drag forces, statistics time sample and output of solution and statistics
    virtual void StatisticsAndOutput(){ dserror("not implemented!");
                                        return; };

    /// set initial flow field
    virtual void SetInitialFlowField(int whichinitialfield,int startfuncno){ dserror("not implemented!");
                                                                             return; };

    ///  set time-step-related fields for low-Mach-number flow
    virtual void SetTimeLomaFields(RCP<const Epetra_Vector> densnp,
                                   RCP<const Epetra_Vector> densn,
                                   RCP<const Epetra_Vector> densnm){ dserror("not implemented!");
                                                                     return; };

    /// set outer-iteration-related fields for low-Mach-number flow
    virtual void SetIterLomaFields(RCP<const Epetra_Vector> densnp){ dserror("not implemented!");
                                                                     return; };

    /// set outer-iteration-related fields for low-Mach-number flow
    virtual void SetGenAlphaIterLomaFields(RCP<const Epetra_Vector> densnp,
                                           RCP<const Epetra_Vector> densdtnp){ dserror("not implemented!");
                                                                     return; };


    /*--------------------------------------------------------------------------------------------* 
     * group member functions which are implemented and used
     * 
     * These are the used member functions                                               henke 08/08
     *------------------------------------------------------------------------------------------- */

    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual void PrepareTimeStep();
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);
    virtual void Update();
    virtual void Output();
    
    virtual void NonlinearSolve();
    /// get the linear solver object used for this field

    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return solver_; }

    /// the interface map contains all velocity dofs at the interface
    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap();
    virtual Teuchos::RCP<const Epetra_Map> VelocityRowMap();
    virtual Teuchos::RCP<const Epetra_Map> PressureRowMap();

    /// scaling factor needed to convert the residual to real forces
    virtual double ResidualScaling() const;

    virtual double TimeScaling() const;

    /// communication object at the interface
    virtual const LINALG::MapExtractor& Interface() const { return interface_; }

    virtual void ReadRestart(int step);
    virtual double Time();
    virtual int Step();
    virtual int  Itemax() const;
    virtual void SetItemax(int itemax);

    ///   return type of time integration scheme
    virtual const FLUID_TIMEINTTYPE TimIntScheme() const {return fluid_.TimIntScheme();}
    virtual void LiftDrag();

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln();
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel);

    virtual Teuchos::RCP<DRT::ResultTest>  CreateFieldTest();

  private:

    void PrintInterfaceVectorField(
        const Teuchos::RCP<Epetra_Vector>   displacementfield,
        const Teuchos::RCP<Epetra_Vector>   vectorfield,
        const std::string filestr,
        const std::string name_in_gmsh
        );

    /// the actual fluid time integration scheme
    FLD::CombustFluidImplicitTimeInt             fluid_;

    /// the interface map setup for interface <-> full translation
    LINALG::MapExtractor                   interface_;

    //! @name local copies of input parameters
    Teuchos::RCP<DRT::Discretization>      dis_;
    Teuchos::RCP<DRT::Discretization>      boundarydis_;
    Teuchos::RCP<LINALG::Solver>           solver_;
    Teuchos::RCP<Teuchos::ParameterList>   params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;
    Teuchos::RCP<IO::DiscretizationWriter> boundaryoutput_;
    //@}

    //! @name state vectors using the solid boundary parallel distribution (row vectors)
//    Teuchos::RCP<Epetra_Vector>            idispnp_;   ///< interface displacement  n+1
//    Teuchos::RCP<Epetra_Vector>            idispn_;    ///< interface displacement  n
    Teuchos::RCP<Epetra_Vector>            itrueresnp_;///< interface residual      n+1
    Teuchos::RCP<Epetra_Vector>            ivelnp_;    ///< interface velocity      n+1
    Teuchos::RCP<Epetra_Vector>            iveln_;     ///< interface velocity      n
    Teuchos::RCP<Epetra_Vector>            ivelnm_;    ///< interface velocity      n-1
    Teuchos::RCP<Epetra_Vector>            iaccn_;     ///< interface acceleration  n
    //@}

    /// all velocity dofs not at the interface
    Teuchos::RefCountPtr<Epetra_Map>       innervelmap_;

    /// ALE dof map
//    LINALG::MapExtractor                   meshmap_;
  };

} // namespace ADAPTER

#endif
#endif
