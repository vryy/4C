/*!----------------------------------------------------------------------
\file fluid_windkessel_optimization.cpp
\brief evaluation of windkessel optimization bc

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/


#include <stdio.h>

#include "fluid_windkessel_optimization.H"

#include "../linalg/linalg_ana.H"



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::UTILS::FluidWkOptimizationWrapper::FluidWkOptimizationWrapper(
  RCP<DRT::Discretization> actdis,
  IO::DiscretizationWriter& output,
  RCP<FLD::UTILS::FluidImpedanceWrapper> ImpWrapper,
  double dta) :
  // call constructor for "nontrivial" objects
  discret_(actdis),
  output_ (output),
  imp_ (ImpWrapper),
  period_(0.0)
{

  vector<DRT::Condition*> impedancecond;
  discret_->GetCondition("ImpedanceCond",impedancecond);

  // the number of lines of impedance boundary conditions found in the input
  // note that some of these lines could belong to the same physical condition
  // which is then marked by the same 'ConditionID'
  // The corresponding resorting of result values has to be done later
  int numcondlines = impedancecond.size();

  map<const int, RCP<DRT::Condition> > wkmap;

  if (numcondlines > 0) // if there is at least one impedance condition
  {
    // -------------------------------------------------------------------
    // get time period length of first condition, this should always be
    // the same!
    // -------------------------------------------------------------------
    period_ = (impedancecond[0])->GetDouble("timeperiod");
    
    // -------------------------------------------------------------------
    // now care for the fact that there could be more than one input line
    // belonging to the same impedance boundary condition
    // -------------------------------------------------------------------
    for (int i=0; i<numcondlines; i++)
    {
      int condid = (impedancecond[i])->GetInt("ConditionID");
      
      double thisperiod = (impedancecond[i])->GetDouble("timeperiod");
      if (thisperiod != period_)
      {
        dserror("all periods of impedance conditions in one problem have to be the same!!!");
        exit(1);
      }
      
      // -----------------------------------------------------------------
      // stack the impedances into a map
      // -----------------------------------------------------------------
      RCP<DRT::Condition> wkCond = Teuchos::rcp(new DRT::Condition(*(impedancecond[i])));
      bool inserted = wkmap.insert( make_pair( condid, wkCond ) ).second;
      
      if ( !inserted )
      {
	dserror("There are more than one impedance condition lines with the same ID. This can not yet be handled.");
        exit(1);
      }
      
    } // end loop over condition lines from input
  } // end if there were conditions
  
  // ---------------------------------------------------------------------
  // Check if there is any optimzation conditions
  // ---------------------------------------------------------------------
  
  vector<DRT::Condition*> wk_optim_cond;
  discret_->GetCondition("Windkessel_Optimization_Cond",wk_optim_cond);
  
  int numOptlines = wk_optim_cond.size();
  if (numOptlines > 0 )
  {
    // -------------------------------------------------------------------
    // test that we have an integer number of time steps per cycle
    // -------------------------------------------------------------------
    // something more intelligent could possibly be found one day ...
    double doublestepnum = period_/dta;
    int cyclesteps = (int)(doublestepnum+0.5)+1;

    // Problems total dimension
    int Dim = 0;
    // -------------------------------------------------------------------
    // Check if all of the optimization conditions are linked to existing
    // windkessel boundary conditions, and initialize x_
    // -------------------------------------------------------------------

    for (int i = 0; i< numOptlines; i++)
    {
      int optID = wk_optim_cond[i]->GetInt("ConditionID");
      
      // check whether optimization condition has an impedance condition
      if (wkmap.find(optID) == wkmap.end())
      {
        dserror("Windkessel optimization condition %d has no corresponding impedance condition",optID);
        exit(1);
      }

      // -----------------------------------------------------------------
      // sort optimization bc's in map and test if one condition ID
      // appears more than once. Currently this case is forbidden.
      // -----------------------------------------------------------------
      RCP<DRT::Condition> wkoptCond = Teuchos::rcp(new DRT::Condition(*(wk_optim_cond[i])));
      bool inserted = optwkmap_.insert( make_pair( optID, wkoptCond ) ).second;
      if ( !inserted )
      {
	dserror("There are more than one impedance condition lines with the same ID. This can not yet be handled.");
        exit(1);
      }
      Dim += GetObjectiveFunctionSize(wkoptCond,optID);
    }

    Jacobian_  = Teuchos::rcp(new Epetra_SerialDenseMatrix(Dim,Dim));
    Jnm_       = Teuchos::rcp(new Epetra_SerialDenseMatrix(Dim,Dim));
    fn_        = Teuchos::rcp(new Epetra_SerialDenseVector(Dim));
    fnm_       = Teuchos::rcp(new Epetra_SerialDenseVector(Dim));
    xn_        = Teuchos::rcp(new Epetra_SerialDenseVector(Dim));
    xnm_       = Teuchos::rcp(new Epetra_SerialDenseVector(Dim));

    dN_du_     = Teuchos::rcp(new Epetra_SerialDenseMatrix(numOptlines*cyclesteps,numOptlines*cyclesteps));
    dN_dphi_   = Teuchos::rcp(new Epetra_SerialDenseMatrix(numOptlines*cyclesteps,Dim));
    dL_du_     = Teuchos::rcp(new Epetra_SerialDenseMatrix(Dim,numOptlines*cyclesteps));
    dJ_dphi_   = Teuchos::rcp(new Epetra_SerialDenseMatrix(Dim,Dim));
    du_dphi_   = Teuchos::rcp(new Epetra_SerialDenseMatrix(numOptlines*cyclesteps,Dim));

    // -----------------------------------------------------------------
    // initial optimization step is zero
    // -----------------------------------------------------------------
    step_ = 0;
    
  }
} // end FluidWkOptimizationWrapper

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Dimension of the design function (public)               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
int FLD::UTILS::FluidWkOptimizationWrapper::GetObjectiveFunctionSize(RCP<DRT::Condition> cond, int condid)
{
  // Get name of objective function
  string ObjFunType = *(cond->Get<string>("ObjectiveFunction"));

  // Get type of design variables
  string DesignVars = *(cond->Get<string>("DesignVariables"));

  int ObjDim = 0;
  int VarDim = 0;

  // -------------------------------------------------------------------
  // Get the dimension of the objective function
  // -------------------------------------------------------------------
  if (ObjFunType == "Psys_Pdia")
  {
    ObjDim = 2;
  }
  else if (ObjFunType == "Psys_Pdia_Pavg")
  {
    ObjDim = 3;
  }
  else
  {
    dserror("[%s]: No such objective function defined. Check windkessel optimization BC [%d]",ObjFunType.c_str(),condid);
    exit(1);
  }

  // -------------------------------------------------------------------
  // Get the dimension of the design variables
  // -------------------------------------------------------------------
  if (DesignVars == "R_C")
  {
    VarDim = 2;
  }
  else if (DesignVars == "R1_R2_C")
  {
    VarDim = 3;
  }
  else
  {
    dserror("[%s]: No such design variables defined. Check windkessel optimization BC [%d]",DesignVars.c_str(),condid);
    exit(1);
  }

  // -------------------------------------------------------------------
  // So far we don't allow situations were the objective function
  // design variables have different size
  // -------------------------------------------------------------------
  if (VarDim != ObjDim)
  {
    dserror("Dimension of design variables and objective functions must agree. Check windkessel optimization BC [%d]",condid);
    exit(1);
  }

  return ObjDim;

}// end GetObjectiveFunctionSize


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Solve the optimization step (public)                    ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::Solve(ParameterList params)
{

  // -------------------------------------------------------------------
  // do nothing if the windkessel optimization bc is not defined
  // -------------------------------------------------------------------
  if (optwkmap_.size()== 0)
  {
    return;
  }

  // -------------------------------------------------------------------
  // Get time step size
  // -------------------------------------------------------------------
  const double dt = params.get<double> ("time step size");

  if (! this->SteadyStateIsObtained(params,dt))
  {
    return;
  }

  // -------------------------------------------------------------------
  // Read in all of the pressures and evaluate the objective function
  // -------------------------------------------------------------------
  map<int, RCP<DRT::Condition> >::iterator itr;
  
  int index = 0;
  int constrain_num = 0;

  // -------------------------------------------------------------------
  // reset dJ_dphi to zero
  // -------------------------------------------------------------------
  dJ_dphi_->Scale(0.0);

  // -------------------------------------------------------------------
  // reset dL_du_ to zero
  // -------------------------------------------------------------------
  dL_du_->Scale(0.0);

  // -------------------------------------------------------------------
  // Reset dN_dphi to zero
  // -------------------------------------------------------------------
  dN_dphi_->Scale(0.0);

  // -------------------------------------------------------------------
  // reset dN_du_ to zero
  // -------------------------------------------------------------------
  dN_du_-> Scale(0.0);


  for (itr = optwkmap_.begin(); itr != optwkmap_.end(); itr++)
  {
    int dim = this->FluidWkOptimizationWrapper::GetObjectiveFunctionSize(itr->second,
                                                                         itr->first);

    // -----------------------------------------------------------------
    // Get the pointers of the corresponding flow and pressure vectors
    // -----------------------------------------------------------------

    std::stringstream pstream, qstream;
    pstream<<"pressures"<<itr->first;
    qstream<<"flowrates"<<itr->first;

    RCP<vector<double> > pressures = params.get<RCP<vector<double> > > (pstream.str());
    RCP<vector<double> > flowrates = params.get<RCP<vector<double> > > (qstream.str());

    double pres_T, flow_T;
    pres_T = (*pressures)[0];
    flow_T = (*flowrates)[0];
    pressures->push_back(pres_T);
    flowrates->push_back(flow_T);
    cout<<"AP: Number of opt conditions: "<<optwkmap_.size()<<endl;

    // -----------------------------------------------------------------
    // do some tests unless someone messed up the code some where
    // -----------------------------------------------------------------
    {
      if (pressures->size() != flowrates->size())
      {
        dserror("[Can't optimise]:Flowrates and pressures have two different dimensions");
      }
      else if ((int) (pressures->size()*optwkmap_.size()) != dN_du_->RowDim())
      {
        dserror("[Can't optimise]: Pressures (%d) and state-variables (%d) size don't match",pressures->size(),dN_du_->RowDim()/optwkmap_.size());
      }
    }

#if 0
    cout<<"COND("<<itr->first<<"): pushing pressure [0]: "<<(*pressures)[0]<<endl;
    cout<<"COND("<<itr->first<<"): pushing flowrate [0]: "<<(*flowrates)[0]<<endl;
    // -----------------------------------------------------------------
    // print pressure for debugging reasons
    // -----------------------------------------------------------------
    if(discret_->Comm().MyPID() == 0)
    {
      cout<<"Printing Pressures and flowrates"<<endl;
      for (unsigned int i = 0; i<pressures->size();i++)
      {
        cout<<"Cond("<<itr->first<<")\t"<<double(i)*dt<<"\t["<<i<<"]\t"<<(*pressures)[i];
        cout<<"\t"<<(*flowrates)[i]<<endl;
      }
    }
#endif

    // -----------------------------------------------------------------
    // evaluate xnm_
    // -----------------------------------------------------------------
    this->GetDesignVariables(itr->first,
                             index);

    // -----------------------------------------------------------------
    // calculate objective function
    // -----------------------------------------------------------------
    this->FluidWkOptimizationWrapper::CalcObjFunction(index,
                                                      itr->second,
                                                      pressures,
                                                      flowrates,
                                                      dt);
#if 1
    if(discret_->Comm().MyPID() == 0)
    {
      cout<<"fn is:"<<endl;
      cout<<(*fn_)<<endl;
    }
#endif

    // -----------------------------------------------------------------
    // calculate partial derivative of objective function w.r.t
    // state variables
    // -----------------------------------------------------------------    
    this->FluidWkOptimizationWrapper::dL_du(index,
                                            constrain_num,
                                            itr->second,
                                            pressures,
                                            flowrates,
                                            params,
                                            dt);
#if 0
    if(discret_->Comm().MyPID() == 0)
    {
      cout<<"dLdu is:"<<endl;
      cout<<(*dL_du_)<<endl;
    }
#endif

    // -----------------------------------------------------------------
    // calculate partial derivative of constrain function w.r.t
    // state variables
    // -----------------------------------------------------------------
    this->FluidWkOptimizationWrapper::dN_du(index,
                                            constrain_num,
                                            itr->second,
                                            pressures,
                                            flowrates,
                                            params,
                                            dt);
#if 0
    if(discret_->Comm().MyPID() == 0)
    {
      cout<<"dNdu is:"<<endl;
      cout<<(*dN_du_)<<endl;
    }
#endif

    // -----------------------------------------------------------------
    // calculate partial derivative of constrain function w.r.t
    // design variables
    // -----------------------------------------------------------------
    this->FluidWkOptimizationWrapper::dN_dphi(index,
                                              constrain_num,
                                              itr->second,
                                              pressures,
                                              flowrates,
                                              params,
                                              dt);
#if 0
    if(discret_->Comm().MyPID() == 0)
    {
      cout<<"dNdphi is:"<<endl;
      cout<<(*dN_dphi_)<<endl;
    }
#endif

    // -----------------------------------------------------------------
    // calculate partial derivative of design function w.r.t
    // design variables
    // -----------------------------------------------------------------    
    this->FluidWkOptimizationWrapper::dJ_dphi(index,
                                              constrain_num,
                                              itr->second,
                                              pressures,
                                              flowrates,
                                              params,
                                              dt);
#if 0
    if(discret_->Comm().MyPID() == 0)
    {
      cout<<"dJdphi is:"<<endl;
      cout<<(*dJ_dphi_)<<endl;
    }
#endif

    // -----------------------------------------------------------------
    // increment the constraint number by one
    // increment the index number by the problem dimension
    // -----------------------------------------------------------------
    constrain_num ++;
    index += dim;

    pressures->pop_back();
    flowrates->pop_back();
  }

  // -------------------------------------------------------------------
  // calculate the jacobian
  // -------------------------------------------------------------------
  this->FluidWkOptimizationWrapper::CalcAdjointJacobian();
  solver_.SetMatrix(*Jacobian_);
  solver_.SetVectors(*xn_,*fn_);
  solver_.FactorWithEquilibration(true);

  // -------------------------------------------------------------------
  // check if du/dphi
  // -------------------------------------------------------------------
  int err2 = solver_.Factor();
  int err  = solver_.Solve();
  
  if (err!=0 || err2!=0)
  {
    dserror("Unable to solve for du/dphi while claculating the adjoint Jacobian");
  }
  
  // -------------------------------------------------------------------
  // since xn_ has the result of dx = xnm_ -xn_
  // then some manipulation of xn_ must be done to get the new results
  // P.S: this must be modified in the future
  // -------------------------------------------------------------------
  (*xn_).Scale(-1.0);
  (*xn_) += *xnm_ ;

  // -------------------------------------------------------------------
  // update the windkessel parameters
  // -------------------------------------------------------------------
  this->UpdateResidual();

}// end Solve


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  calculate objective function (public)                   ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::CalcObjFunction(
  int                   index,
  RCP<DRT::Condition>   cond,
  RCP<vector<double> >  pressures,
  RCP<vector<double> >  flowrates,
  double                dt
  )
{

  // Get name of objective function
  string ObjFunType = *(cond->Get<string>("ObjectiveFunction"));

  // Get condition Id
  const int condid = cond->GetInt("ConditionID");

  // Find Pmax and Pmin
  double Pmax = (*pressures)[0];
  double Pmin = (*pressures)[0];
  //  int index_max = 0;
  //  int index_min = 0;
  for (unsigned int i = 0; i< pressures->size(); i++)
  {
    if (Pmax<(*pressures)[i])
    {
      Pmax = (*pressures)[i];
      //      index_max = i;
    }
    if (Pmin>(*pressures)[i])
    {
      Pmin = (*pressures)[i];
      //      index_min = i;
    }
  }


  // -------------------------------------------------------------------
  // Get the dimension of the objective function
  // -------------------------------------------------------------------
  if (ObjFunType == "Psys_Pdia")
  {
    const double Psys = cond->GetDouble("Psystolic");
    const double Pdia = cond->GetDouble("Pdiastolic");

    (*fn_)[index   ] = fabs(Pmax - Psys);
    (*fn_)[index +1] = fabs(Pmin - Pdia);

  }
  else if (ObjFunType == "Psys_Pdia_Pavg")
  {

  }
  else
  {
    dserror("[%s]: No such objective function defined. Check windkessel optimization BC [%d]",ObjFunType.c_str(),condid);
    exit(1);
  }

}// end CalcObjFunction



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::dN_du(
  int                   index,
  int                   constrain_num,
  RCP<DRT::Condition>   cond,
  RCP<vector<double> >  pressures,
  RCP<vector<double> >  flowrates,
  ParameterList         params,
  double                dt)
{
  //  int    VarDim = 0;
  // Get type of design variables
  string DesignVars = *(cond->Get<string>("DesignVariables"));


  // -------------------------------------------------------------------
  // Get the dimension of the design variables
  // -------------------------------------------------------------------
  if (DesignVars == "R_C")
  {
    const double alfa = cond->GetDouble("R1R2_ratio");
    //    VarDim = 2;
    
    // -----------------------------------------------------------------
    //   (*xn_)[index] = R1 + R2 = R
    //    alfa = R1/R2
    // => R1 = R*alfa/(alfa+1)
    // => R2 = R - R1
    // -----------------------------------------------------------------
    const double R1 = (*xnm_)[index  ]*alfa/(1.0+alfa);
    const double R2 = (*xnm_)[index  ] - R1;
    const double C  = (*xnm_)[index+1];

    // -----------------------------------------------------------------
    // for constrain function:
    //  N:   Q.(1+R1/R2) + C.R1.dQ/dt - C dP/dt - P/R2 = 0
    //
    //      par N
    //  => ------- =
    //      par u
    //  _                                                             _
    // |               .               .               .               |
    // | 0.5 + R2.C/dt .       0       .      ...      . 0.5 + R2.C/dt |
    // |               .               .               .               |              
    // | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . |
    // |               .               .               .               |
    // | 0.5 - R2.C/dt . 0.5 + R2.C/dt .      ...      .      0        |                          
    // |               .               .               .               |
    // | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . |
    // |               . .             . .             .               |
    // |      0        .    `  .       .    `  .       .      0        |
    // |               .          `  . .          `  . .               |
    // | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . |
    // |               .               .               .               |
    // |      0        .       0       . 0.5 - R2.C/dt . 0.5 + R2.C/dt |
    // |               .               .               .               |
    // |_                                                             _| 
    // -----------------------------------------------------------------
    for (unsigned int i = (constrain_num)*(pressures->size()); i<(constrain_num+1)*(pressures->size()); i++)
    {
      unsigned int k = (constrain_num)*(pressures->size());
      (*dN_du_)(i,i) = 0.5 + R2*C/dt;
      if (i> k)
      {
        (*dN_du_)(i,i-1) = 0.5 - R2*C/dt;
      }
      else
      {
        (*dN_du_)(i,(constrain_num+1)*(pressures->size())-1) = 0.5 - R2*C/dt;
      }
    }
    return;
  }
  else if (DesignVars == "R1_R2_C")
  {
    //    VarDim = 3;
  }
  else
  {
    dserror("[%s]: No such design variables defined. Check windkessel optimization BC [%d]",DesignVars.c_str(),cond->GetInt("ConditionID"));
    exit(1);
  }

}// FluidWkOptimizationWrapper::dN_du


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::dN_dphi(
  int                   index,
  int                   constrain_num,
  RCP<DRT::Condition>   cond,
  RCP<vector<double> >  pressures,
  RCP<vector<double> >  flowrates,
  ParameterList         params,
  double                dt)
{
  //  int    VarDim = 0;
  // Get type of design variables
  string DesignVars = *(cond->Get<string>("DesignVariables"));


  // -------------------------------------------------------------------
  // Get the dimension of the design variables
  // -------------------------------------------------------------------
  if (DesignVars == "R_C")
  {
    const double alfa = cond->GetDouble("R1R2_ratio");
    //    VarDim = 2;
    
    // -----------------------------------------------------------------
    //   (*xn_)[index] = R1 + R2 = R
    //    alfa = R1/R2
    // => R1 = R*alfa/(alfa+1)
    // => R2 = R - R1
    // -----------------------------------------------------------------
    const double R1 = (*xnm_)[index  ]*alfa/(1.0+alfa);
    const double R2 = (*xnm_)[index  ] - R1;
    const double R  = (*xnm_)[index  ];
    const double C  = (*xnm_)[index+1];

    // -----------------------------------------------------------------
    // for constrain function:
    //  N:   Q.(1+R1/R2) + C.R1.dQ/dt - C dP/dt - P/R2 = 0
    //
    //  alfa = R1/R2
    //
    //      par  N
    //  => -------- =
    //      par phi
    //  _                                                             .                               _
    // |      C                  1             2.C.alfa.R             .  R2              R2.R1         |
    // | -----------.(P - P ) - ---(Q + Q ) - -------------.(Q - Q )  . ----.(P - P ) - -----.(Q - Q ) |
    // | dt.(1+alfa)   0  -1     2  -1   0    dt.(1+alfa)^2   0  -1   .  dt    0  -1      dt    0   -1 |
    // | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . |
    // |      C                  1             2.C.alfa.R             .  R2              R2.R1         |
    // | -----------.(P - P ) - ---(Q + Q ) - -------------.(Q - Q )  . ----.(P - P ) - -----.(Q - Q ) |
    // | dt.(1+alfa)   1   0     2   0   1    dt.(1+alfa)^2   1   0   .  dt    1   0      dt    1   0  |
    // | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . |
    // |                                                              .                                |
    // |                               .                              .               .                |
    // |                               .                              .               .                |
    // |                               .                              .               .                |
    // |                                                              .                                |
    // | . . . . . . . . . . . . . . . . . . . . . . . . . .. . . . . . . . . . . . . . . . . . . . . .|
    // |      C                  1             2.C.alfa.R             .  R2              R2.R1         |
    // | -----------.(P - P ) - ---(Q + Q ) - -------------.(Q - Q )  . ----.(P - P ) - -----.(Q - Q ) |
    // | dt.(1+alfa)   n  nm     2   nm  n    dt.(1+alfa)^2   n  nm   .  dt    n  nm      dt    n  nm  |
    // |_. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .  _|

    // -----------------------------------------------------------------

    double Qn  = (*flowrates)[0];
    double Qnm = (*flowrates)[flowrates->size()-2];
    double Pn  = (*pressures)[0];
    double Pnm = (*pressures)[pressures->size()-2];

    int k = (constrain_num)*(pressures->size());
    for (unsigned int i = 0; i<pressures->size(); i++)
    {
      // ---------------------------------------------------------------
      // Update the Qn and Pn
      // ---------------------------------------------------------------
      Qn  = (*flowrates)[i]; 
      Pn  = (*pressures)[i]; 

      // ---------------------------------------------------------------
      // Evaluate dN_dphi_
      // ---------------------------------------------------------------
      (*dN_dphi_)(i+k,index  ) = C/(dt*(1.0+alfa))*(Pn-Pnm)
        - 0.5*(Qn+Qnm) - 2.0*C*alfa*R/((pow(1.0+alfa,2.0))*dt)*(Qn-Qnm);

      (*dN_dphi_)(i+k,index+1) = R2/(dt)*(Pn-Pnm)
        - R1*R2/dt*(Qn-Qnm);

      // ---------------------------------------------------------------
      // Update the Qnm and Pnm
      // ---------------------------------------------------------------
      Qnm = Qn;
      Pnm = Pn;
    }
    return;
  }
  else if (DesignVars == "R1_R2_C")
  {
    //    VarDim = 3;
  }
  else
  {
    dserror("[%s]: No such design variables defined. Check windkessel optimization BC [%d]",DesignVars.c_str(),cond->GetInt("ConditionID"));
    exit(1);
  }
}// FluidWkOptimizationWrapper::dN_dphi

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::dL_du(
  int                   index,
  int                   constrain_num,
  RCP<DRT::Condition>   cond,
  RCP<vector<double> >  pressures,
  RCP<vector<double> >  flowrates,
  ParameterList         params,
  double                dt)
{

  //  int    ObjDim = 0;
  // Get name of objective function
  string ObjFunType = *(cond->Get<string>("ObjectiveFunction"));

  // -------------------------------------------------------------------
  // Get the dimension of the objective functions
  // -------------------------------------------------------------------
  if (ObjFunType == "Psys_Pdia")
  {
    //    ObjDim = 2;

    const double Psys = cond->GetDouble("Psystolic");
    const double Pdia = cond->GetDouble("Pdiastolic");

    // Find Pmax and Pmin
    double Pmax = (*pressures)[0];
    double Pmin = (*pressures)[0];
    unsigned int index_max = 0;
    unsigned int index_min = 0;
    for (unsigned int i = 0; i< pressures->size(); i++)
    {
      if (Pmax<(*pressures)[i])
      {
        Pmax = (*pressures)[i];
        index_max = i;
      }
      if (Pmin>(*pressures)[i])
      {
        Pmin = (*pressures)[i];
        index_min = i;
      }
    }

    // -----------------------------------------------------------------
    // for constrain function:
    //        _                _     _   _
    //       |                  |   |     |
    //       |P    - P          |   |  0  |
    //  f:   | max    systolic  | = |     |
    //       |                  |   |     |
    //       |P    - P          |   |  0  |
    //       |_min    diastolic_|   |_   _|
    //
    //      par L     par f
    //  => ------- = ------ =
    //      par u     par u
    //
    //                    index(Pmax)        index(Pmin)
    //                    __________         ___________
    //                        |                   |
    //  _                     v                   v                     _
    // | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . |
    // |  0  .  0  . - - - .  1  .  0  . - - - .  0  .  0  . - - - .  0  |
    // | . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . |
    // |  0  .  0  . - - - .  0  .  0  . - - - .  1  .  0  . - - - .  0  |
    // |_. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ._|
    //
    // -----------------------------------------------------------------

    unsigned int k = (constrain_num)*(pressures->size());
    for (unsigned int i = 0; i< pressures->size() + k; i++)
    {
      if(i == (index_max) )
      {
        if (Pmax - Psys >= 0.0)
        {
          (*dL_du_)(index  ,i+k) =  1.0;
        }
        else
        {
          (*dL_du_)(index  ,i+k) = -1.0;
        }
        (*dL_du_)(index+1,i+k) = 0.0;
      }
      if(i == (index_min) )
      {
        (*dL_du_)(index  ,i+k) = 0.0;
        if (Pmin - Pdia >= 0.0)
        {
          (*dL_du_)(index+1,i+k) = 1.0;
        }
        else
        {
          (*dL_du_)(index+1,i+k) =-1.0;          
        }
      }
    }
    return;
  }
  else if (ObjFunType == "Psys_Pdia_Pavg")
  {
    //    ObjDim = 3;
  }
  else
  {
    dserror("[%s]: No such objective function defined. Check windkessel optimization BC [%d]",ObjFunType.c_str(),cond->GetInt("ConditionID"));
    exit(1);
  }
}// FluidWkOptimizationWrapper::dL_du


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::dJ_dphi(
  int                   index,
  int                   constrain_num,
  RCP<DRT::Condition>   cond,
  RCP<vector<double> >  pressures,
  RCP<vector<double> >  flowrates,
  ParameterList         params,
  double                dt)
{
  // Get name of objective function
  string ObjFunType = *(cond->Get<string>("ObjectiveFunction"));

  // Get type of design variables
  string DesignVars = *(cond->Get<string>("DesignVariables"));

  //  int ObjDim = 0;
  //  int VarDim = 0;

  // -------------------------------------------------------------------
  // Get the dimension of the objective function
  // -------------------------------------------------------------------
  if (ObjFunType == "Psys_Pdia")
  {
    if (DesignVars == "R_C")
    {
      //      VarDim = 2;
      
      //      int k = (constrain_num)*(pressures->size());
      //      for (unsigned int i = k; i<k+pressures->size(); i++)
      //      {
      //        for (unsigned int j = k; j< k+pressures->size();j++)
      //        {
      //          (*dJ_dphi_)[i][j] = 0.0;
      //        }
      //      }
    }
    else
    {
      dserror("[%s]: No such design variables defined for [] objective function. Check windkessel optimization BC [%d]",DesignVars.c_str(), ObjFunType.c_str(),cond->GetInt("ConditionID"));
      exit(1);
    }
    //    ObjDim = 2;
  }
  else if (ObjFunType == "Psys_Pdia_Pavg")
  {
    //    ObjDim = 3;
    if (DesignVars == "R1_R2_C")
    {
      //      VarDim = 2;
    }
    else
    {
      dserror("[%s]: No such design variables defined for [] objective function. Check windkessel optimization BC [%d]",DesignVars.c_str(), ObjFunType.c_str(),cond->GetInt("ConditionID"));
      exit(1);
    }
  }
  else
  {
    dserror("[%s]: No such objective function defined. Check windkessel optimization BC [%d]",ObjFunType.c_str(),cond->GetInt("ConditionID"));
    exit(1);
  }


}// FluidWkOptimizationWrapper::dJ_dphi

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
bool FLD::UTILS::FluidWkOptimizationWrapper::SteadyStateIsObtained(
  ParameterList params,
  double        dt)
{

  // -------------------------------------------------------------------
  // first check if the results are at time = n*CardiacPeriod
  // -------------------------------------------------------------------
  double time = params.get<double> ("total time");


  // -------------------------------------------------------------------
  // if time is less that period, no need to even check for steadiness
  // that is because such a case mean nothing
  // -------------------------------------------------------------------  
  if (time<period_)
  {
    return false;
  }

  //  if (double(n)*period_ != time)
  //  {
    // -----------------------------------------------------------------
    // if results are not at n*CardiacPeriod, then return false
    // -----------------------------------------------------------------  
    //  return false;
  //  }

  // -------------------------------------------------------------------
  // if results are at n*CardiacPeriod, then
  // check if all of the boundaries reached a steady state
  // -------------------------------------------------------------------
  map<const int, RCP<DRT::Condition> >::iterator itr;

  bool converged = false;
  for(itr = optwkmap_.begin(); itr != optwkmap_.end(); itr++)
  {
    // -----------------------------------------------------------------
    // read in pressure of condition i
    // -----------------------------------------------------------------
    std::stringstream stream, stream2, stream3;
    stream<< "pressures"<<itr->first;
    stream2<<"dP"<<itr->first;
    stream3<<"EndOfCycle"<<itr->first;

#if 0
    // -----------------------------------------------------------------
    // print pressure for debugging reasons
    // -----------------------------------------------------------------
    if(discret_->Comm().MyPID() == 0)
    {
      cout<<"Printing Pressures"<<endl;
      for (unsigned int i = 0; i<pressures->size();i++)
      {
        cout<<"PressureCond("<<itr->first<<")\t"<<double(i)*dt<<"\t"<<(*pressures)[i]<<endl;
      }
    }
#endif

    // -----------------------------------------------------------------
    // calculate the steady state error, which is the difference between
    // the last and first value of pressure
    // -----------------------------------------------------------------
    double error = fabs( params.get<double>(stream2.str()));
    
    // -----------------------------------------------------------------
    // This looks silly, but the user must be careful that the tolerance
    // is a positive value
    // -----------------------------------------------------------------
    const double tol = itr->second->GetDouble("Tolerance");

    if (tol<=0.0)
    {
      dserror("windkessel optimization condition (%d) must have tolerance > 0.0 ",itr->first);
      exit(1);
    }

    // -----------------------------------------------------------------
    // if the error is bigger than the tolerance then return false
    // -----------------------------------------------------------------
    if (error > tol)
    {
      return false;
    }
    else if(params.get<bool>(stream3.str()))
    {
      converged = true;
    }
    else
    {
      return false;
    }
  }

  //if (params.get<bool>(stream3.str()))
  //{
    // -----------------------------------------------------------------
    // if results are not at n*CardiacPeriod, then return false
    // -----------------------------------------------------------------  
  //    return false;
      //}

  // -------------------------------------------------------------------
  // if this stage is obtained then all of the windkessel conditions 
  // reached to a steady state periodicity
  // -------------------------------------------------------------------
  return converged;
}// FluidWkOptimizationWrapper::SteadyStateIsObtained


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::CalcAdjointJacobian()
{


  // -------------------------------------------------------------------
  // Jacobian = dL/du * du/dphi + dJ/dphi
  // du/dphi  = inv(dN/du)*(-dN/dphi)
  // -------------------------------------------------------------------


  // -------------------------------------------------------------------
  // solving for du/dphi
  // -------------------------------------------------------------------

  dN_dphi_->Scale(-1.0);
  Jacobian_->Scale(0.0);
  solver_.SetMatrix(*dN_du_);
  solver_.SetVectors(*du_dphi_,*dN_dphi_);
  solver_.FactorWithEquilibration(true);

  // -------------------------------------------------------------------
  // check if du/dphi
  // -------------------------------------------------------------------
  int err2 = solver_.Factor();
  int err  = solver_.Solve();
  
  if (err!=0 || err2!=0)
  {
    dserror("Unable to solve for du/dphi while claculating the adjoint Jacobian");
  }


  (*Jacobian_).Multiply('N','N', 1.0, *dL_du_, *du_dphi_, 0.0);
  (*Jacobian_) += *dJ_dphi_;

#if 0
  if (discret_->Comm().MyPID() == 0)
  {
    cout<<"du_dphi is:"<<endl;
    cout<<(*du_dphi_)<<endl;

    cout<<"Jacobian is:"<<endl;
    cout<<(*Jacobian_)<<endl;
  }
#endif

}// FluidWkOptimizationWrapper::CalcAdjointJacobian()

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::UpdateResidual()
{
  // -------------------------------------------------------------------
  // do nothing if the windkessel optimization bc is not defined
  // -------------------------------------------------------------------
  if (optwkmap_.size()== 0)
  {
    return;
  }

  *fnm_ = *fn_;
  *xnm_ = *xn_;
  
  int index = 0;
  int VarDim= 0;
  map<const int, RCP<DRT::Condition> >::iterator  itr;
  for (itr = optwkmap_.begin(); itr != optwkmap_.end();itr++)
  {
    ParameterList params;

    // -----------------------------------------------------------------
    // Get the dimension of the design variables
    // -----------------------------------------------------------------
    string DesignVars = *(itr->second->Get<string>("DesignVariables"));

    if (DesignVars == "R_C")
    {
      const double alfa = itr->second->GetDouble("R1R2_ratio");
      VarDim = 2;

      // ---------------------------------------------------------------
      //   (*xn_)[index] = R1 + R2 = R
      //    alfa = R1/R2
      // => R1 = R*alfa/(alfa+1)
      // => R2 = R - R1
      // ---------------------------------------------------------------
      const double R1 = (*xn_)(index  )*alfa/(1.0+alfa);
      const double R2 = (*xn_)(index  ) - R1;
      const double C  = (*xn_)(index+1);

      // ---------------------------------------------------------------
      // Update the impedance parameters
      // ---------------------------------------------------------------
      params.set<double>("R1",R1);
      params.set<double>("R2",R2);
      params.set<double>("C" ,C);
      imp_->SetWindkesselParams(params,itr->first);

      // ---------------------------------------------------------------
      // print out results
      // ---------------------------------------------------------------
      int myrank = discret_->Comm().MyPID();

      if (myrank == 0)
      {
        cout<<"Cond(" << itr->first << ") Adjoint Step(" << step_ << "): R1: " << R1 <<endl;
        cout<<"Cond(" << itr->first << ") Adjoint Step(" << step_ << "): R2: " << R2 <<endl;
        cout<<"Cond(" << itr->first << ") Adjoint Step(" << step_ << "):  C: " << C  <<endl;
      }
    }
    else if (DesignVars == "R1_R2_C")
    {
      VarDim = 3;
      const double R1 = (*xn_)(index  );
      const double R2 = (*xn_)(index+1);
      const double C  = (*xn_)(index+2);

      // ---------------------------------------------------------------
      // Update the impedance parameters
      // ---------------------------------------------------------------
      params.set<double>("R1",R1);
      params.set<double>("R2",R2);
      params.set<double>("C" ,C);
      imp_->SetWindkesselParams(params,itr->first);

      // ---------------------------------------------------------------
      // print out results
      // ---------------------------------------------------------------
      int myrank = discret_->Comm().MyPID();

      if (myrank == 0)
      {
        cout<<"Cond(" << itr->first << ") Adjoint Step(" << step_ << "): R1: " << R1 <<endl;
        cout<<"Cond(" << itr->first << ") Adjoint Step(" << step_ << "): R2: " << R2 <<endl;
        cout<<"Cond(" << itr->first << ") Adjoint Step(" << step_ << "):  C: " << C  <<endl;
      }
    }
    else
    {
      dserror("[%s]: No such design variables defined. Check windkessel optimization BC [%d]",DesignVars.c_str(),itr->first);
      exit(1);
    }
    index += VarDim;
  }

  //--------------------------------------------------------------------
  // recalculate the impedance
  //--------------------------------------------------------------------
  imp_->Impedances();

  //--------------------------------------------------------------------
  // update the step numbers
  //--------------------------------------------------------------------
  step_ ++;

}// FluidWkOptimizationWrapper::UpdateResidual


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::GetDesignVariables(
  int condid,
  int index)
{
  // -------------------------------------------------------------------
  // do nothing if the windkessel optimization bc is not defined
  // -------------------------------------------------------------------
  if (optwkmap_.size()== 0)
  {
    return;
  }

  ParameterList params;

  // -----------------------------------------------------------------
  // Get the dimension of the design variables
  // -----------------------------------------------------------------
  string DesignVars = *(optwkmap_[condid]->Get<string>("DesignVariables"));

  if (DesignVars == "R_C")
  {

    // ---------------------------------------------------------------
    // Get the windkessel parameters
    // ---------------------------------------------------------------
    imp_->GetWindkesselParams(params,condid);

    const double R1 = params.get<double>("R1");
    const double R2 = params.get<double>("R2");
    const double C  = params.get<double>("C" );

    // ---------------------------------------------------------------
    //   (*xn_)[index] = R1 + R2 = R
    //    alfa = R1/R2
    // => R1 = R*alfa/(alfa+1)
    // => R2 = R - R1
    // ---------------------------------------------------------------
    (*xnm_)(index  ) = R1 + R2;
    (*xnm_)(index+1) = C;
  }
  else if (DesignVars == "R1_R2_C")
  {
    // ---------------------------------------------------------------
    // Get the windkessel parameters
    // ---------------------------------------------------------------
    imp_->GetWindkesselParams(params,condid);

    const double R1 = params.get<double>("R1");
    const double R2 = params.get<double>("R2");
    const double C  = params.get<double>("C" );

    // ---------------------------------------------------------------
    //   (*xn_)[index] = R1 + R2 = R
    //    alfa = R1/R2
    // => R1 = R*alfa/(alfa+1)
    // => R2 = R - R1
    // ---------------------------------------------------------------
    (*xnm_)(index  ) = R1;
    (*xnm_)(index+1) = R2;
    (*xnm_)(index+2) = C;
  }
  else
  {
    dserror("[%s]: No such design variables defined. Check windkessel optimization BC [%d]",DesignVars.c_str(),condid);
    exit(1);
  }
#if 0
  if(discret_->Comm().MyPID() == 0)
  {
    cout<<"xnm is:"<<endl;
    cout<<(*xnm_);
  }
#endif

}// FluidWkOptimizationWrapper::GetDesignVariables


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::WriteRestart(
  IO::DiscretizationWriter&  output)
{
  // -------------------------------------------------------------------
  // do nothing if the windkessel optimization bc is not defined
  // -------------------------------------------------------------------
  if (optwkmap_.size()== 0)
  {
    return;
  }

  std::stringstream stream1;

  // output optimization step  
  stream1<<"Optimization_Step";
  output.WriteInt(stream1.str(), step_);

}//FluidWkOptimizationWrapper::WriteRestart


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |   (public)                                               ismail 02/10|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void FLD::UTILS::FluidWkOptimizationWrapper::ReadRestart(
  IO::DiscretizationReader& reader)
{

  // -------------------------------------------------------------------
  // do nothing if the windkessel optimization bc is not defined
  // -------------------------------------------------------------------
  if (optwkmap_.size()== 0)
  {
    return;
  }

  std::stringstream stream1;
  
  // read in step number
  stream1<<"Optimization_Step";
  step_ = reader.ReadInt(stream1.str());

}//FluidWkOptimizationWrapper::ReadRestart

