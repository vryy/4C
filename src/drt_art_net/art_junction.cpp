/*!----------------------------------------------------------------------
\file art_junction.cpp
\brief evaluation of 1d-artery junction bc

\maintainer Lena Yoshihara

\level 3

*----------------------------------------------------------------------*/


#include <stdio.h>

#include "art_junction.H"

#include "../linalg/linalg_ana.H"

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 08/09|
 |                                                                      |
 |                                                                      |
 |  Each junction boundary condition is defined as multiple conditions. |
 |  Taking the following example:                                       |
 |                          ______                                      |
 |                         |      \                                     |
 |                         |   4   \     5         6                    |
 |                         |   o---------o---------o                    |
 |                         |  [3]   |             [4]                   |
 |                         |        |(Artery 2)                         |
 |                         |        |                                   |
 |       1         2       | 3      |                                   |
 |       o---------o---------o      |<---(Junction 1)                   |
 |      [1]                |[2]     |                                   |
 |             (Artery 1)  |        |                                   |
 |                         |        |                                   |
 |                         |   7    |    8         9                    |
 |                         |   o---------o---------o                    |
 |                         |  [5]  /              [6]                   |
 |                         |______/  (Artery 3)                         |
 |                                                                      |
 |                                                                      |
 |    The junction (bifurcation in this case) is connected to the       |
 |    design points [2], [3], and [5], usinf the following definition   |
 |    of the junction boundary contion:                                 |
 |       E [design nodes number] - [junction number (ID)] [options]     |
 |    we will end up with the following expression:                     |
 |       E 2 - 1 [other options]                                        |
 |       E 3 - 1 [other options]                                        |
 |       E 5 - 1 [other options]                                        |
 |                                                                      |
 | .................................................................... |
 |                                                                      |
 |                                                                      |
 |    To create a junction bc we first need to separate the junctions   |
 |    by:                                                               |
 |    (1) Read in all junction boundary conditions                      |
 |                                                                      |
 |    (2) Check where the boundary conditions are connected to the      |
 |        inlet or the outlet of an artery element. This is easily done |
 |        by checking whether the design node of an element is locally  |
 |        has and in/outlet condition flag as "inlet" or "outlet"       |
 |                                                                      |
 |    (3) Group all the conditions having the same ID number            |
 |                                                                      |
 |    (4) Create the junctions.                                         |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
ART::UTILS::ArtJunctionWrapper::ArtJunctionWrapper(Teuchos::RCP<DRT::Discretization> actdis,
                                                   IO::DiscretizationWriter& output,
                                                   Teuchos::ParameterList & params,
                                                   double dta):
  discret_(actdis),
  output_ (output)
{

  //----------------------------------------------------------------------
  // Exit if the function accessed by a non-master processor
  //----------------------------------------------------------------------
  if (discret_->Comm().MyPID()==0)
  {
    //----------------------------------------------------------------------
    // (1) Get the junction boundary conditions
    //----------------------------------------------------------------------

    std::vector<DRT::Condition*> myConditions;
    discret_->GetCondition("ArtJunctionCond",myConditions);
    int numofcond = myConditions.size();

    if(numofcond==1)
    {
      dserror("An arterial junction is supposed to have at least two nodes!");
    }
    else if(numofcond>1) // if there is atleast two arteries connected to each other
    {
      //----------------------------------------------------------------------
      // (2) check whether the condition is connected to an inlet(-1) or
      //     to an outlet(1)
      //----------------------------------------------------------------------

      std::vector<int> IOart(numofcond);
      for(int i =0; i<numofcond; i++)
      {
        // get the node number connected to the condition
        const std::vector<int> * nodes = myConditions[i]->Nodes();

        // The junction condition must be connected to one and only one node
        if(nodes->size()!=1)
          dserror("Artery Connection BC should have only one node connected to it!");

        int local_id =  discret_->NodeRowMap()->LID((*nodes)[0]);
        // Get the actual node connected to the condition
        DRT::Node * nd = discret_->lColNode(local_id);

        // find whether the nodes is at the inlet or at the outlet of the element
        std::string terminalType = *(nd->GetCondition("ArtInOutCond")->Get<std::string>("terminaltype"));
        if(terminalType == "inlet")
          IOart[i] = -1;
        else if (terminalType == "outlet")
          IOart[i] = 1;
        else
          dserror("Something is severely wrong! In/Out terminal condition should be either \"outlet\" or \"inlet\"");
      }
      //----------------------------------------------------------------------
      // (3) Group all of the conditions that belong to the same junction
      //----------------------------------------------------------------------

      DRT::Condition * cond_i;

      //first, sort the condition list according to there IDs
      //In this case the bubble sort algorithm is used
      int IO_i;
      for(int i=myConditions.size(); i>1; i--)
      {

        for(int j = 1; j< i; j++)
        {
          // if Id(j-1) > Id(j) then swap the two values
          if(myConditions[j-1]->GetInt("ConditionID")>myConditions[j]->GetInt("ConditionID"))
          {
            cond_i = myConditions[j];
            IO_i   = IOart[j];
            myConditions[j] = myConditions[j-1];
            IOart[j]        = IOart[j-1];
            myConditions[j-1] = cond_i;
            IOart[j-1]        = IO_i;
          }
        }
      }

      // second, group all the similar conditions in one vector
      std::vector<std::vector<DRT::Condition*> > SortedConds;
      std::vector<DRT::Condition *> grouped_cond;

      std::vector<std::vector<int> > SortedIOarts;
      std::vector<int> grouped_IO;

      for(unsigned int i=0; i<myConditions.size();)
      {
        do
        {
          grouped_IO.push_back(IOart[i]);
          grouped_cond.push_back(myConditions[i++]);

          if(i==myConditions.size())
          break;
        }
        while(myConditions[i]->GetInt("ConditionID") == grouped_cond[0]->GetInt("ConditionID"));

        SortedConds.push_back(grouped_cond);
        grouped_cond.erase(grouped_cond.begin(),grouped_cond.end());

        SortedIOarts.push_back(grouped_IO);
        grouped_IO.erase(grouped_IO.begin(),grouped_IO.end());
      }

      // ---------------------------------------------------------------------
      // (4) Create junction boundary conditions
      // ---------------------------------------------------------------------
      int condid;
      Teuchos::RCP<std::map<const int, Teuchos::RCP<JunctionNodeParams> > >  nodalParams;
      nodalParams = params.get<Teuchos::RCP<std::map<const int, Teuchos::RCP<JunctionNodeParams> > > >("Junctions Parameters");

      for(unsigned int i=0; i<SortedConds.size(); i++)
      {
        // -------------------------------------------------------------------
        // allocate the junction bc class members for every case
        // -------------------------------------------------------------------
        condid = SortedConds[i][0]->GetInt("ConditionID");

        // -------------------------------------------------------------------
        // sort junction BCs in map
        // -------------------------------------------------------------------
        Teuchos::RCP<ArtJunctionBc> junbc = Teuchos::rcp(new ArtJunctionBc(discret_, output_, SortedConds[i], SortedIOarts[i],dta, condid, i) );
        ajunmap_.insert( std::make_pair( condid, junbc ) );

        // -------------------------------------------------------------------
        // Creat the nodes' parameters (material prameters, geometric
        // parameters, n-1 values) in map, and export all the values so that
        // elements could access them.
        // Finally check wheather a node has multiple BC, which is not allowed
        // -------------------------------------------------------------------
        bool inserted;
        // create an empty map associated to the Teuchos::RCP nodalParams_
        //      nodalParams = Teuchos::rcp(new std::map<const int, Teuchos::RCP<JunctionNodeParams> >());

        for (unsigned int j=0 ; j< SortedConds[i].size(); j++)
        {
          const std::vector<int> * nodes = SortedConds[i][j]->Nodes();
          Teuchos::RCP<JunctionNodeParams> nodeparams = Teuchos::rcp(new JunctionNodeParams);

          int local_id =  discret_->NodeRowMap()->LID((*nodes)[0]);
          inserted = nodalParams->insert( std::make_pair( local_id, nodeparams) ).second;
          if(!inserted)
            dserror("Node %d has more than one condition", (*nodes)[0]+1);
        }

      }

    }// end if there is a connection
  }

}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Distructor (public)                                     ismail 08/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
ART::UTILS::ArtJunctionWrapper::~ArtJunctionWrapper()
{
  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Solve      (public)                                     ismail 08/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
int ART::UTILS::ArtJunctionWrapper::Solve(Teuchos::ParameterList & params)
{
  //----------------------------------------------------------------------
  // Exit if the function accessed by a non-master processor
  //----------------------------------------------------------------------

  if (discret_->Comm().MyPID()!=0)
    return 0;

  std::map<const int, Teuchos::RCP<class ArtJunctionBc> >::iterator mapiter;

  for (mapiter = ajunmap_.begin(); mapiter != ajunmap_.end(); mapiter++ )
  {
    mapiter->second->ArtJunctionBc::Solve(params);
  }
  return 0;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 08/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
ART::UTILS::ArtJunctionBc::ArtJunctionBc( Teuchos::RCP<DRT::Discretization>  actdis,
                                          IO::DiscretizationWriter& output,
                                          std::vector<DRT::Condition*> conds,
                                          std::vector<int> IOart_flag,
                                          double dta,
                                          int condid,
                                          int numcond):
  condid_(condid),
  discret_(actdis),
  output_ (output),
  IOart_flag_(IOart_flag)
{
  //----------------------------------------------------------------------
  // Check whether all the nodes have simillar flow direction
  // i.e. whether they all are inlets or all are outlets for the junctions
  //----------------------------------------------------------------------
  int IOartFlag = IOart_flag_[0];
  bool IOartFlags_are_fine = false;
  for(unsigned int i=1; i<IOart_flag_.size(); i++)
  {
    if(IOart_flag_[i]!= IOartFlag)
    {
      IOartFlags_are_fine = true;
      break;
    }
  }
  if(!IOartFlags_are_fine)
  {
    if(IOartFlag==1)
      dserror("Junction (%d) has all of its nodes defined as outlets",conds[0]->GetInt("ConditionID"));
    else
      dserror("Junction (%d) has all of its nodes defined as inlets",conds[0]->GetInt("ConditionID"));
  }

  //----------------------------------------------------------------------
  // Find the size of the nonlinear problem. In this case each nodes is
  // supossed to have two degrees of freedom, i.e. a junction with "N"
  // nodes must have 2*N degrees of freedom to be solved, which in turn
  // is the size of the nonlinear problem
  //----------------------------------------------------------------------
  ProbSize_ = 2*IOart_flag_.size();

  //----------------------------------------------------------------------
  // Extracting the nodes to whome the junction is connected
  //----------------------------------------------------------------------
  for(unsigned int i=0; i<conds.size(); i++)
    nodes_.push_back( (*(conds[i]->Nodes()))[0] );

}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Solve       (public)                                    ismail 08/09|
 |                                                                      |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 | Implimenting the junction boundary condition such that:              |
 |                                                                      |
 |                                                                      |
 | Parent 1    ______________             _______________  Daughter 1   |
 | ---------> ()____________ \           / _____________()..--->        |
 |                          \ \         / /                             |
 | Parent 2    ______________\ \_______/ /_______________  Daughter 2   |
 | ---------> ()_____________             ______________()..--->        |
 |             .     .       |           |             .                |
 |             .     .       . Junction  .             .                |
 |             .     .       .           .             .                |
 | Parent p    ______________|           |_______________  Daughter d   |
 | ---------> ()_____________  ________   ______________()..--->        |
 |                          / /        \ \                              |
 | Parent N    ____________/ /          \ \______________  Daughter M   |
 | ---------> ()____________/            \______________()..--->        |
 |                                                                      |
 |                                                                      |
 | * Mass Conservation Equations:                                       |
 |     _____            _____                                           |
 |     \                \                                               |
 |      \    A * U    =  \    A  * U                                    |
 |      /     p   p      /     d    d                                   |
 |     /____            /____                                           |
 |       N                M                                             |
 |                                                                      |
 | * Moment conservation equations:                                     |
 |                                                                      |
 |                                                                      |
 |  rho    2                       ___     ___       rho    2           |
 |  --- * U  * P   + P    + beta(_/A   - _/Ao  )  =  --- * U  * P       |
 |   2     p    p    ext            p        p        2     d    d      |
 |                                                                      |
 |                                             ___     ___              |
 |                            + P    +  beta(_/A   - _/Ao  )  + H       |
 |                               ext            d        d       Lpd    |
 |                                                                      |
 |  Where, H    is the pressure loss due to the bifurcation type        |
 |          Lpd                                                         |
 |         (e.g. due to sharp angles at the Aortic Arch)                |
 |                                                                      |
 | Solving the nonlinear system:                                        |
 | -----------------------------                                        |
 | Solving this nonlinear system is possible with the help of           |
 | Neuton-Raphson method, which can be implimented as following:        |
 |                                                                      |
 |   1- The previous equations can generate the residual                |
 |                                                                      |
 |      function f( U ) = 0                                             |
 |               ~  ~     ~                                             |
 |        ^                                                             |
 |  f(U)  |                                      o                      |
 |    ~   |                                     o                       |
 |        |                                     o/ par f                |
 |        |                                    o/  ----- (U )           |
 |  f(U ) +  -  -  -  -  -  -  -  -  -  -  -  +/   par U   i            |
 |    ~i  |                                 o /                         |
 |        |                               o  /.                         |
 |        |                             o   / .                         |
 |        |                           o    /  .                         |
 |        |                         o     /   .                         |
 |        |                      o       /    .                         |
 |        |                   o         /     .                         |
 |        |                o           /      .                         |
 |        |             o             /       .                         |
 |      --+---------+----------------+--------+---------------------->  |
 |        |      o                  /U        U                    U    |
 |                                   ~i+1     ~i                   ~    |
 |                       -1                        par f                |
 |   U    =  U   -  H(U )  * f(U )      Where  H = ------               |
 |   ~i+1    ~i       ~i       ~i                  par U                |
 |                                                                      |
 |                                                                      |
 |   2- Generate H and find the liniarized residual f                   |
 |                                                                      |
 |                                                                      |
 |   3- Solve for the linearized residual f                             |
 |                                                                      |
 |                                                                      |
 |   4- Update U and go back to step (1) until norm (f) <= Tol          |
 |                                                                      |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
int ART::UTILS::ArtJunctionBc::Solve(Teuchos::ParameterList & params)
{

  //----------------------------------------------------------------------
  // Define the matricese and the vectors that are needed to solve the
  // nonlinear problem at the junction
  //----------------------------------------------------------------------
  Epetra_SerialDenseMatrix Jacobian(ProbSize_,ProbSize_);
  Epetra_SerialDenseVector f(ProbSize_);
  Epetra_SerialDenseVector x(ProbSize_);
  Epetra_SerialDenseVector dx(ProbSize_);

  //----------------------------------------------------------------------
  // Read the element information at the node of the bifurcation
  //----------------------------------------------------------------------
  std::vector<double> A(ProbSize_/2,0.0);
  std::vector<double> Q(ProbSize_/2,0.0);
  std::vector<double> W(ProbSize_/2,0.0);
  std::vector<double> Ao(ProbSize_/2,0.0);
  std::vector<double> rho(ProbSize_/2,0.0);
  std::vector<double> beta(ProbSize_/2,0.0);
  std::vector<double> Pext(ProbSize_/2,0.0);

  // get the map having the junction nodal information from the elements
  Teuchos::RCP<std::map<const int, Teuchos::RCP<JunctionNodeParams> > > nodalMap =  params.get< Teuchos::RCP<std::map<const int, Teuchos::RCP<JunctionNodeParams> > > >("Junctions Parameters");

  // loop over all the nodes and read in the required parameters
  for(unsigned int i = 0; i< nodes_.size(); i++)
  {
    int local_id =  discret_->NodeRowMap()->LID(nodes_[i]);
    A[i]    = (*nodalMap)[local_id]->A_;
    Q[i]    = (*nodalMap)[local_id]->Q_;
    W[i]    = (*nodalMap)[local_id]->W_;
    Ao[i]   = (*nodalMap)[local_id]->Ao_;
    rho[i]  = (*nodalMap)[local_id]->rho_;
    beta[i] = (*nodalMap)[local_id]->beta_;
    Pext[i] = (*nodalMap)[local_id]->Pext_;
    //Intializing x vector; x = [U1 U2 ... Un A1 A2 ... An]^T
    x[i]               = Q[i]/A[i];
    x[i+nodes_.size()] = A[i];

#if 0
    // for debuging reasons
    std::cout<<"A   ["<<i<<"] is: "<<A[i]<<std::endl;
    std::cout<<"Q   ["<<i<<"] is: "<<Q[i]<<std::endl;
    std::cout<<"W   ["<<i<<"] is: "<<W[i]<<std::endl;
    std::cout<<"Ao  ["<<i<<"] is: "<<Ao[i]<<std::endl;
    std::cout<<"rho ["<<i<<"] is: "<<rho[i]<<std::endl;
    std::cout<<"beta["<<i<<"] is: "<<beta[i]<<std::endl;
    std::cout<<"Pext["<<i<<"] is: "<<Pext[i]<<std::endl;
    std::cout<<std::endl;
#endif
  }


  //----------------------------------------------------------------------
  // Fill the Residual vector
  //----------------------------------------------------------------------
  Residual_Eval(f, A, Q, W, Ao, rho, beta, Pext);

  int  itr = 0;

  // a vector specifying the pivots (reordering)
  int *pivot;
  pivot = new int[2*nodes_.size()];

  while(f.Norm2()>0.000001)
  {

    //--------------------------------------------------------------------
    // Fill the Jacobian matrix
    //--------------------------------------------------------------------
    Jacobian_Eval(Jacobian, A, Q, W, Ao, rho, beta, Pext);

    //--------------------------------------------------------------------
    // Solve for dx
    //--------------------------------------------------------------------
    solver_.SetMatrix(Jacobian);
    solver_.SetVectors(dx,f);
    solver_.FactorWithEquilibration(true);
    int err2 = solver_.Factor();
    int err  = solver_.Solve();

    if (err!=0 || err2!=0)
    {
      dserror("Unable to solve for the jacobian in junction %d, error number %d", condid_, err);
    }

    //--------------------------------------------------------------------
    // Update x = x + dx = x - f
    //--------------------------------------------------------------------
    for(unsigned int i = 0; i< nodes_.size(); i++)
    {
      x[i]               = Q[i]/A[i];
      x[i+nodes_.size()] = A[i];
    }
    Update_Result(x,dx);


    //--------------------------------------------------------------------
    // the junction is not converging! exit
    //--------------------------------------------------------------------
    if(itr++ >= 20)
    {
      delete [] pivot;
      dserror("Junction [%d] is not converging!", condid_);
    }

    //--------------------------------------------------------------------
    // Fill the Residual vector
    //--------------------------------------------------------------------
    for(unsigned int i = 0; i< nodes_.size(); i++)
    {
      A[i] = x[i+nodes_.size()];
      Q[i] = x[i+nodes_.size()]*x[i];
    }

    Residual_Eval(f, A, Q, W, Ao, rho, beta, Pext);

  }
  delete [] pivot;
  std::cout<<"Junction "<<condid_<<" is solved in "<<itr;
  if(itr == 1)
    std::cout<<" iteration"<<std::endl;
  else
    std::cout<<" iterations"<<std::endl;

  //----------------------------------------------------------------------
  // Update the final results
  //----------------------------------------------------------------------
  for(unsigned int i = 0; i< nodes_.size(); i++)
  {
    int local_id =  discret_->NodeRowMap()->LID(nodes_[i]);
    (*nodalMap)[local_id]->A_ = A[i];
    (*nodalMap)[local_id]->Q_ = Q[i];
  }

  return 0;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Evaluate Jacobian (public)                              ismail 09/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::UTILS::ArtJunctionBc::Jacobian_Eval( Epetra_SerialDenseMatrix & Jacobian,
                                               std::vector<double> &A,
                                               std::vector<double> &Q,
                                               std::vector<double> &W,
                                               std::vector<double> &Ao,
                                               std::vector<double> &rho,
                                               std::vector<double> &beta,
                                               std::vector<double> &Pext)
{

  // empty the Jacobian
  Jacobian = Epetra_SerialDenseMatrix(ProbSize_,ProbSize_);

  // fill the entities that have to do with forward charachteristic speeds
  for(unsigned int i = 0; i< nodes_.size(); i++)
  {
    Jacobian(i,i              ) = 1.0;
    Jacobian(i,i+nodes_.size()) = double(IOart_flag_[i])*sqrt(beta[i]/(2.0*Ao[i]*rho[i]))/pow(A[i],0.75);
#if 0
    std::cout<<"beta["<<i<<"] : "<<beta[i]<<std::endl;
    std::cout<<"A   ["<<i<<"] : "<<A[i]<<std::endl;
    std::cout<<"Ao  ["<<i<<"] : "<<Ao[i]<<std::endl;
    std::cout<<"rho ["<<i<<"] : "<<rho[i]<<std::endl;
    std::cout<<std::endl;
#endif
  }

  // fill the entities that have to do with the mass conservation
  for(unsigned int i = 0; i< nodes_.size(); i++)
  {
    Jacobian(nodes_.size(),i                ) = double(IOart_flag_[i]) * A[i];
    Jacobian(nodes_.size(),i + nodes_.size()) = double(IOart_flag_[i]) * Q[i]/A[i];
  }

  // fill the entities that have to do with the pressure conservation
  const double P_u = rho[0]*(Q[0]/A[0]);
  const double P_A = 0.5*beta[0]/(Ao[0]*sqrt(A[0]));
  for(unsigned int i=1;i<nodes_.size();i++)
  {
    Jacobian(i+nodes_.size(),0              ) =  P_u;
    Jacobian(i+nodes_.size(),i              ) = -rho[i]*Q[i]/A[i];
    Jacobian(i+nodes_.size(),nodes_.size()  ) =  P_A;
    Jacobian(i+nodes_.size(),nodes_.size()+i) = -0.5*beta[i]/(Ao[i]*sqrt(A[i]));
  }

}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Evaluate Residual (public)                              ismail 09/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::UTILS::ArtJunctionBc::Residual_Eval( Epetra_SerialDenseVector & f,
                                               std::vector<double> &A,
                                               std::vector<double> &Q,
                                               std::vector<double> &W,
                                               std::vector<double> &Ao,
                                               std::vector<double> &rho,
                                               std::vector<double> &beta,
                                               std::vector<double> &Pext)
{

  // initialize the residual
   f = Epetra_SerialDenseVector(f.Length());

  // fill the entities that have to do with forward charachteristic speeds
  for(unsigned int i = 0; i< nodes_.size(); i++)
  {
    f[i] = Q[i]/A[i] + double(IOart_flag_[i])*4.0* sqrt(beta[i]/(2.0*Ao[i]*rho[i])*sqrt(A[i])) - W[i];
  }

  // fill the entities that have to do with the mass conservation
  f[nodes_.size()] = 0.0;
  for(unsigned int i = 0; i< nodes_.size(); i++)
  {
    f[nodes_.size()] += double(IOart_flag_[i])*Q[i];
  }

  // fill the entities that have to do with the pressure conservation
  const double P0 = 0.5*rho[0]*pow(Q[0]/A[0],2) + beta[0]*(sqrt(A[0]) - sqrt(Ao[0]))/Ao[0];
  for(unsigned int i = 1; i< nodes_.size(); i++)
  {
    f[nodes_.size()+i] = P0 - (0.5*rho[i]*pow(Q[i]/A[i],2) + beta[i]*(sqrt(A[i]) - sqrt(Ao[i]))/Ao[i]);
  }

}



//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Update Residual (public)                                ismail 09/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::UTILS::ArtJunctionBc::Update_Result( Epetra_SerialDenseVector & xn,
                                               Epetra_SerialDenseVector & dx)
{
#if DEBUG
  if(xn.Length() != dx.Length())
  {
    dserror("Both, the result and the result change, must have similar size");
  }
#endif

  for (int i=0; i<xn.Length();i++)
  {
    xn[i]-=dx[i];
  }
}


