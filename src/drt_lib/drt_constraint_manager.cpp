/*!----------------------------------------------------------------------
\file drt_constraint_manager.cpp

\brief Class controlling constraint and containing the necessary data

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_constraint_manager.H"
#include "iostream"


/*----------------------------------------------------------------------*
 |  ctor (public)                                               tk 11/07|
 *----------------------------------------------------------------------*/
ConstrManager::ConstrManager(DRT::Discretization& discr,
        RCP<Epetra_Vector> disp):
actdisc_(discr)
{
    //Check, what kind of constraining boundary conditions there are
    numConstrID_=0;
    haveareaconstr_=false;
    havevolconstr_=false;
    //Check for volume constraints
    vector<DRT::Condition*> constrcond(0);
    actdisc_.GetCondition("VolumeConstraint_3D",constrcond);
    // Keep ParameterList p alive during initialization, so global information
    // over IDs as well as element results stored here can be used after all
    // constraints are evaluated
    ParameterList p;
    //Deal with volume constraints
    if (constrcond.size())
    {
        p.set("action","calc_struct_constrvol");
        actdisc_.SetState("displacement",disp);
        actdisc_.EvaluateCondition(p,"VolumeConstraint_3D");
        havevolconstr_=true;
    }
    // Check for Area Constraints
    actdisc_.GetCondition("AreaConstraint_3D",constrcond);
    //Deal with area constraints
    if (constrcond.size())
    {
        p.set("action","calc_struct_monitarea");
        actdisc_.SetState("displacement",disp);
        actdisc_.EvaluateCondition(p,"AreaConstraint_3D");
        haveareaconstr_=true;
    }
    //----------------------------------------------------
    //--------------------------possible other constraints
    //----------------------------------------------------
    haveconstraint_= haveareaconstr_||havevolconstr_;
    if (haveconstraint_)
    {
        ManageIDs(p,minConstrID_,maxConstrID_,numConstrID_);
        // sum up initial values
        initialvalues_=rcp(new Epetra_SerialDenseVector(numConstrID_));
        initialvalues_->Scale(0.0);
        SynchronizeSumConstraint(p,initialvalues_,"computed volume",numConstrID_,minConstrID_);
        SynchronizeSumConstraint(p,initialvalues_,"computed area",numConstrID_,minConstrID_);
        //Initialize Lagrange Multiplicators, reference values and errors
        lagrMultVec_=rcp(new Epetra_SerialDenseVector(numConstrID_));
        lagrMultInc_=rcp(new Epetra_SerialDenseVector(numConstrID_));
        lagrMultVec_->Scale(0.0);
        lagrMultInc_->Scale(0.0);
        actdisc_.ClearState();
        referencevalues_=rcp(new Epetra_SerialDenseVector(numConstrID_));
        actvalues_=rcp(new Epetra_SerialDenseVector(numConstrID_));
        actvalues_->Scale(0.0);
        constrainterr_=rcp(new Epetra_SerialDenseVector(numConstrID_));
        const Epetra_Map* dofrowmap = actdisc_.DofRowMap();
        constrVec_=rcp(new Epetra_MultiVector(*dofrowmap,numConstrID_));
        fact_=rcp(new Epetra_SerialDenseVector(numConstrID_));
    }
    havevolmonitor_=false;
    haveareamonitor_=false;
    //Check for Volume Monitors
    actdisc_.GetCondition("VolumeMonitor_3D",constrcond);
    ParameterList p1;
    if (constrcond.size())
    {
        p1.set("action","calc_struct_constrvol");
        actdisc_.SetState("displacement",disp);
        actdisc_.EvaluateCondition(p1,"VolumeMonitor_3D");
        havevolmonitor_=true;
    }
    // Check for Area Monitor
    actdisc_.GetCondition("AreaMonitor_3D",constrcond);
    //Deal with area Monitors
    if (constrcond.size())
    {
        p1.set("action","calc_struct_monitarea");
        actdisc_.SetState("displacement",disp);
        actdisc_.EvaluateCondition(p1,"AreaMonitor_3D");
        haveareamonitor_=true;
    }
    havemonitor_= haveareamonitor_||havevolmonitor_;
    if (havemonitor_)
    {
        ManageIDs(p1,minMonitorID_,maxMonitorID_,numMonitorID_);
        // sum up initial values
        monitorvalues_=rcp(new Epetra_SerialDenseVector(numMonitorID_));
        initialmonvalues_=rcp(new Epetra_SerialDenseVector(numMonitorID_));
        initialmonvalues_->Scale(0.0);
        SynchronizeSumConstraint(p1,initialmonvalues_,"computed volume",numMonitorID_,minMonitorID_);
        SynchronizeSumConstraint(p1,initialmonvalues_,"computed area",numMonitorID_,minMonitorID_);
    }
    return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 11/07|
|Compute difference between current and prescribed values.              |
|Change Stiffnessmatrix and internal force vector                       |
*-----------------------------------------------------------------------*/

void ConstrManager::StiffnessAndInternalForces(
        const double time,
        RCP<Epetra_Vector> disp,
        RCP<Epetra_Vector> fint,
        RCP<LINALG::SparseMatrix> stiff)
{
    // create the parameters for the discretization
    ParameterList p;
    actvalues_->Scale(0.0);
    const Epetra_Map* dofrowmap = actdisc_.DofRowMap();
    constrVec_=rcp(new Epetra_MultiVector(*dofrowmap,numConstrID_));
    constrVec_->PutScalar(0.0);
    //Deal with volume constraints
    if (havevolconstr_)
    {
        //Evaluate volume at predicted ENDpoint D_{n+1}
        // action for elements
        p.set("action","calc_struct_volconstrstiff");
        // choose what to assemble
        p.set("assemble matrix 1",true);
        p.set("assemble matrix 2",false);
        p.set("assemble vector 1",true);
        p.set("assemble vector 2",true);
        p.set("assemble vector 3",false);
        // other parameters that might be needed by the elements
        p.set("total time",time);
        p.set("MinID",minConstrID_);
        p.set("MaxID",maxConstrID_);
        p.set("NumberofID",numConstrID_);
        p.set("LagrMultVector",lagrMultVec_);
        // set vector values needed by elements
        actdisc_.ClearState();
        actdisc_.SetState("displacement",disp);
        actdisc_.EvaluateCondition(p,stiff,fint,"VolumeConstraint_3D",constrVec_);
        SynchronizeSumConstraint(p,actvalues_,"computed volume",numConstrID_,minConstrID_);
    }
    //Deal with volume constraints
    if (haveareaconstr_)
    {
        //Evaluate volume at predicted ENDpoint D_{n+1}
        // action for elements
        p.set("action","calc_struct_areaconstrstiff");
        // choose what to assemble
        p.set("assemble matrix 1",true);
        p.set("assemble matrix 2",false);
        p.set("assemble vector 1",true);
        p.set("assemble vector 2",true);
        p.set("assemble vector 3",false);
        // other parameters that might be needed by the elements
        p.set("total time",time);
        p.set("MinID",minConstrID_);
        p.set("MaxID",maxConstrID_);
        p.set("NumberofID",numConstrID_);
        p.set("LagrMultVector",lagrMultVec_);
        // set vector values needed by elements
        actdisc_.ClearState();
        actdisc_.SetState("displacement",disp);
        actdisc_.EvaluateCondition(p,stiff,fint,"AreaConstraint_3D",constrVec_);
        SynchronizeSumConstraint(p,actvalues_,"computed area",numConstrID_,minConstrID_);
    }
    SynchronizeMinConstraint(p,fact_,"LoadCurveFactor");
    *referencevalues_=*initialvalues_;
    for (int iter = 0; iter < numConstrID_; ++iter)
    {
        (*referencevalues_)[iter]*=(*fact_)[iter];
        (*constrainterr_)[iter]=(*referencevalues_)[iter]-(*actvalues_)[iter];
    }
    actdisc_.ClearState();
    // do NOT finalize the stiffness matrix, add mass and damping to it later
    return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Compute difference between current and prescribed values.              |
*-----------------------------------------------------------------------*/
void ConstrManager::ComputeError(double time,RCP<Epetra_Vector> disp)
{
    actvalues_->Scale(0.0);
    ParameterList p;
    p.set("total time",time);
    actdisc_.SetState("displacement",disp);
    if(havevolconstr_)
    {
        p.set("action","calc_struct_constrvol");
        actdisc_.EvaluateCondition(p,"VolumeConstraint_3D");
        SynchronizeSumConstraint(p, actvalues_,"computed volume",numConstrID_,minConstrID_);
    }
    if(haveareaconstr_)
    {
        p.set("action","calc_struct_monitarea");
        actdisc_.EvaluateCondition(p,"AreaConstraint_3D");
        SynchronizeSumConstraint(p, actvalues_,"computed area",numConstrID_,minConstrID_);
    }
    for (int iter = 0; iter < numConstrID_; ++iter)
    {
      (*constrainterr_)[iter]=(*referencevalues_)[iter]-(*actvalues_)[iter];
    }
    return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Update Lagrange increment according to linear uzawa algorithm.          |
*-----------------------------------------------------------------------*/
void ConstrManager::UpdateLagrIncr(double factor, Epetra_SerialDenseVector vect)
{
    for (int i=0;i < numConstrID_; i++)
    {
        (*lagrMultInc_)[i]+= factor*(vect[i]+(*constrainterr_)[i]);
    }
    return;
}
/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Add scaled error of constraint to Lagrange multiplier.                 |
*-----------------------------------------------------------------------*/
void ConstrManager::UpdateLagrMult(double factor)
{
    for (int i=0;i < numConstrID_; i++)
    {
        (*lagrMultVec_)[i]+=factor*(*constrainterr_)[i] ;
    }
    return;
}

/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Add Lagrange increment to Lagrange multiplier.                         |
*-----------------------------------------------------------------------*/
void ConstrManager::UpdateLagrMult()
{
    for (int i=0;i < numConstrID_; i++)
    {
        (*lagrMultVec_)[i]+=(*lagrMultInc_)[i] ;
    }
    return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Compute matrix-vector-product of matrix of constraint conditions with   |
|lagrange increment as needed for uzawa algorithm                        |
*------------------------------------------------------------------------*/
void ConstrManager::ComputeConstrTimesLagrIncr(RCP<Epetra_Vector> dotprod)
{
    dotprod->PutScalar(0.0);
    for (int i=0;i < numConstrID_; i++)
    {
        dotprod->Update(-(*lagrMultInc_)(i),*((*constrVec_)(i)),1.0);
    }
    return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Compute matrix-vector-product of matrix of constraint conditions with   |
|displacement increment as needed to compute residual of uzawa algorithm |
*-----------------------------------------------------------------------*/
void ConstrManager::ComputeConstrTimesDisi(Epetra_Vector disi, RCP<Epetra_SerialDenseVector> dotprod)
{

    for (int i=0; i < numConstrID_; i++)
    {
        double temp;
        disi.Dot(*((*constrVec_)(i)),&temp);
        (*dotprod)(i)=temp;
    }
    return;
}

/*-----------------------------------------------------------------------*
|(public)                                                        tk 01/08|
|Compute values defined to keep track of.                                |
*-----------------------------------------------------------------------*/
void ConstrManager::ComputeMonitorValues(RCP<Epetra_Vector> disp)
{
    monitorvalues_->Scale(0.0);
    ParameterList p;
    actdisc_.SetState("displacement",disp);
    if(havevolmonitor_)
    {
        p.set("action","calc_struct_constrvol");
        actdisc_.EvaluateCondition(p,"VolumeMonitor_3D");
        SynchronizeSumConstraint(p, monitorvalues_,"computed volume",numMonitorID_,minMonitorID_);
    }
    if(haveareamonitor_)
    {
        p.set("action","calc_struct_monitarea");
        actdisc_.EvaluateCondition(p,"AreaMonitor_3D");
        SynchronizeSumConstraint(p, monitorvalues_,"computed area",numMonitorID_,minMonitorID_);
    }
    return;
}
/*----------------------------------------------------------------------*
|(public)                                                       tk 01/08|
|Print monitored values                                                 |
*-----------------------------------------------------------------------*/
void ConstrManager::PrintMonitorValues()
{
    for (int i = 0; i < numMonitorID_; ++i)
    {
        printf("Monitor value %2d: %10.5e (%5.2f%% of initial value)\n",i+minMonitorID_,(*monitorvalues_)[i],((*monitorvalues_)[i])*100/((*initialmonvalues_)[i]));

    }

    return;
}

/*----------------------------------------------------------------------*
 |(private)                                                 tk 11/07    |
 |small subroutine for initialization purposes to determine ID's        |
 *----------------------------------------------------------------------*/
void ConstrManager::ManageIDs(ParameterList& params, int& minID, int& maxID, int& numID)
{
    actdisc_.Comm().MaxAll(&(params.get("MaxID",0)),&maxID,1);
    actdisc_.Comm().MinAll(&(params.get("MinID",100000)),&minID,1);
    numID=1+maxID-minID;
    return;
}

/*----------------------------------------------------------------------*
 |(private)                                                 tk 11/07    |
 |small subroutine to synchronize processors after evaluating the       |
 |constraints by summing them up                                        |
 *----------------------------------------------------------------------*/
void ConstrManager::SynchronizeSumConstraint(ParameterList& params,
                                RCP<Epetra_SerialDenseVector>& vect,
                                const char* resultstring, const int numID, const int minID)
{
    for (int i = 0; i < numID; ++i)
    {
        char valname[30];
        sprintf(valname,"%s %d",resultstring,i+minID);
        double currval=0.0;
        actdisc_.Comm().SumAll(&(params.get(valname,0.0)),&(currval),1);
        (*vect)[i]+=currval;
    }
    return;
}

/*----------------------------------------------------------------------*
 |(private)                                                 tk 11/07    |
 |small subroutine to synchronize processors after evaluating the       |
 |constraint by finding minimum value                                   |
 *----------------------------------------------------------------------*/
void ConstrManager::SynchronizeMinConstraint(ParameterList& params,
                                RCP<Epetra_SerialDenseVector>& vect,
                                const char* resultstring)
{
    for (int i = 0; i < numConstrID_; ++i)
    {
        char valname[30];
        sprintf(valname,"%s %d",resultstring,i+minConstrID_);
        double currval=1;
        actdisc_.Comm().MinAll(&(params.get(valname,1.0)),&(currval),1);
        (*vect)[i]=currval;
    }
    return;
}

#endif
