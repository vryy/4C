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
 |  ctor (public)                                            	tk 11/07|
 *----------------------------------------------------------------------*/
ConstrManager::ConstrManager(DRT::Discretization& discr,
		RCP<Epetra_Vector> disp):
actdisc_(discr)
{
	//Check, what kind of constraining boundary conditions there are
	numConstrID_=0;
	haveconstraint_=false;
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
		haveconstraint_=true;
	}
	//----------------------------------------------------
	//--------------------------possible other constraints
	//----------------------------------------------------
	if (haveconstraint_)
	{
		ManageIDs(p);
		// sum up initial values
		initialvalues_=rcp(new Epetra_SerialDenseVector(numConstrID_));
		initialvalues_->Scale(0.0);
		SynchronizeConstraint(p,initialvalues_,"computed volume");
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
		fact_=1;
	}
	return;
}

/*----------------------------------------------------------------------*
|(public)														tk 11/07|
|Compute difference between current volume and prescribed volume.  		|
|Change Stiffnessmatrix and internal force vector						|
*-----------------------------------------------------------------------*/

void ConstrManager::StiffnessAndInternalForces(
		const double time, 
		RCP<Epetra_Vector> disp,
		RefCountPtr<Epetra_Vector> fint,
		RefCountPtr<Epetra_CrsMatrix> stiff)
{
	//Evaluate volume at predicted ENDpoint D_{n+1} 
    // create the parameters for the discretization
    ParameterList p;
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
	const Epetra_Map* dofrowmap = actdisc_.DofRowMap();
	constrVec_=rcp(new Epetra_MultiVector(*dofrowmap,numConstrID_));
	constrVec_->PutScalar(0.0);
	actdisc_.EvaluateCondition(p,stiff,fint,"VolumeConstraint_3D",constrVec_);
    actvalues_->Scale(0.0);
    SynchronizeConstraint(p,actvalues_,"computed volume");
    fact_=p.get("LoadCurveFactor",1.0);
    *referencevalues_=*initialvalues_;
    referencevalues_->Scale(fact_);
    for (int iter = 0; iter < numConstrID_; ++iter) 
    {
	  (*constrainterr_)[iter]=(*referencevalues_)[iter]-(*actvalues_)[iter];
    } 
    actdisc_.ClearState();
    // do NOT finalize the stiffness matrix, add mass and damping to it later	
	return;
}

void ConstrManager::ComputeError(double time,RCP<Epetra_Vector> disp)
{
	ParameterList p;
	p.set("total time",time);
    p.set("action","calc_struct_constrvol");
    actdisc_.SetState("displacement",disp);
    actdisc_.EvaluateCondition(p,"VolumeConstraint_3D");
    actvalues_->Scale(0.0);
    SynchronizeConstraint(p, actvalues_,"computed volume");
    for (int iter = 0; iter < numConstrID_; ++iter) 
    {
	  (*constrainterr_)[iter]=(*referencevalues_)[iter]-(*actvalues_)[iter];
    }
    return;
}
		
void ConstrManager::UpdateLagrIncr(double factor, Epetra_SerialDenseVector vect)
{
	for (int i=0;i < numConstrID_; i++)
	{
		(*lagrMultInc_)[i]+= factor*(vect[i]+(*constrainterr_)[i]);
	}
	return;
}

void ConstrManager::UpdateLagrMult(double factor)
{
	for (int i=0;i < numConstrID_; i++)
	{
		(*lagrMultVec_)[i]+=factor*(*constrainterr_)[i] ;
	}
	return;	
}

void ConstrManager::UpdateLagrMult()
{
	for (int i=0;i < numConstrID_; i++)
	{
		(*lagrMultVec_)[i]+=(*lagrMultInc_)[i] ;
	}
	return;
}

void ConstrManager::ComputeConstrTimesLagrIncr(RCP<Epetra_Vector> dotprod)
{
	dotprod->PutScalar(0.0);
	for (int i=0;i < numConstrID_; i++)
	{
		dotprod->Update(-(*lagrMultInc_)(i),*((*constrVec_)(i)),1.0);
	}
	return;
}

void ConstrManager::ComputeConstrTimesDisi(Epetra_Vector disi, RCP<Epetra_SerialDenseVector> dotprod)
{
	
	for (int i=0; i < numConstrID_; i++)
	{
		double* temp;
		disi.Dot(*((*constrVec_)(i)),temp);
		(*dotprod)(i)=*temp;
	}
	return;
}

/*----------------------------------------------------------------------*
 |(private)                                                 tk 11/07    |
 |small subroutine for initialization purposes to determine ID's        |
 *----------------------------------------------------------------------*/
void ConstrManager::ManageIDs(ParameterList& params)
{
	actdisc_.Comm().MaxAll(&(params.get("MaxID",0)),&maxConstrID_,1);
	actdisc_.Comm().MinAll(&(params.get("MinID",100000)),&minConstrID_,1);
	numConstrID_=1+maxConstrID_-minConstrID_;
	
	return;
}

/*----------------------------------------------------------------------*
 |(private)                                                 tk 11/07    |
 |small subroutine to synchronize processors after evaluating the       |
 |constraint volume                                              		|
 *----------------------------------------------------------------------*/
void ConstrManager::SynchronizeConstraint(ParameterList& params,
								RCP<Epetra_SerialDenseVector>& vect,
								const char* resultstring)
{
	for (int i = 0; i < numConstrID_; ++i) 
	{
		char volname[30];	
		sprintf(volname,"%s %d",resultstring,i+minConstrID_);
		double currvol;
		actdisc_.Comm().SumAll(&(params.get(volname,0.0)),&(currvol),1);
		(*vect)[i]=currvol;
	}
	return;
}

#endif
