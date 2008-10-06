/*!----------------------------------------------------------------------
\file statmech_manager.cpp
\brief management and auxiliary functions for statistical mechanics

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15234
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "statmech_manager.H"
#include "../drt_lib/drt_validparameters.H"
#include "../drt_lib/drt_utils.H"

#include <iostream>
#include <iomanip>


/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 09/08|
 *----------------------------------------------------------------------*/
StatMechManager::StatMechManager(ParameterList& params, DRT::Discretization& discret):
  statmechparams_( DRT::Problem::Instance()->StatisticalMechanicsParams() ),
  maxtime_(params.get<double>("max time",0.0)),
  starttimeoutput_(-1.0),
  endtoendref_(0.0),
  istart_(0),
  rlink_(params.get<double>("R_LINK",0.0)),
  unusedids_(0),
  discret_(discret)
{ 
  //if dynamic crosslinkers are used additional variables are initialized
  if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
  {
    /*crosslinkerpartner_ and crosslinkerelements_ are generated based on the discretization
     * as it is generated from the input data file; this discretization does not yet comprise
     * any additional crosslinker elements which is o.K. since both variables are applied in 
     * order to administrate only the nodes of the original discretization*/
    crosslinkerpartner_ = rcp( new Epetra_Vector(*discret_.NodeRowMap()) );
    crosslinkerelements_ = rcp( new Epetra_Vector(*discret_.NodeRowMap()) );
    
    /*initializing crosslinkerpartner_ and crosslinkerelements with -1 for each element by looping
     * over all local indices assigned by the map dis.NodeRowMap() to a processor*/
    for(int i = 0; i < discret_.NumMyRowNodes(); i++)
    {
      (*crosslinkerpartner_)[i] = -1;
      (*crosslinkerelements_)[i] = -1;
    }
  }

  return;
} // StatMechManager::StatMechManager


/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechOutput(const double& time,const int& num_dof,const int& istep, const double& dt, const Epetra_Vector& dis)
{
  switch(Teuchos::getIntegralValue<int>(statmechparams_,"SPECIAL_OUTPUT"))
  {
    case INPUTPARAMS::statout_endtoendlog:
    {
      FILE* fp = NULL; //file pointer for statistical output file
      double endtoend = 0; //end to end length at a certain time step in single filament dynamics
      double DeltaR2 = 0;
       
      //as soon as system is equilibrated (after time START_FACTOR*maxtime_) a new file for storing output is generated
      if ( (time >= maxtime_ * statmechparams_.get<double>("START_FACTOR",0.0))  && (starttimeoutput_ == -1.0) )
      {
         int testnumber = 0; 
         endtoendref_ = pow ( pow((dis)[num_dof-3]+10 - (dis)[0],2) + pow((dis)[num_dof-2] - (dis)[1],2) , 0.5);  
         starttimeoutput_ = time;
         istart_ = istep;
         //file pointer for operating with numbering file
         FILE* fpnumbering = NULL;
         std::ostringstream numberingfilename;
         
         //look for a numbering file where number of already existing output files is stored:
         numberingfilename.str("NumberOfRealizations");
         fpnumbering = fopen(numberingfilename.str().c_str(), "r");
         
         //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
         if(fpnumbering == NULL)
         {
           do
           {
             testnumber++;
             outputfilename_.str("");
             outputfilename_ << "EndToEnd"<< testnumber << ".dat";
             fp = fopen(outputfilename_.str().c_str(), "r");
           } while(fp != NULL);
           
           //set up new file with name "outputfilename" without writing anything into this file
           fp = fopen(outputfilename_.str().c_str(), "w");
           fclose(fp);
        }         
        //if there already exists a numbering file
        else
        {
          fclose(fpnumbering);
          
          //read the number of the next realization out of the file into the variable testnumber
          std::fstream f(numberingfilename.str().c_str());            
          while (f)
          {
            std::string tok;
            f >> tok;
            if (tok=="Next")
            {
              f >> tok; 
              if (tok=="Number:")     
                f >> testnumber;
            }
          } //while(f)
          
          //defining outputfilename by meands of new testnumber
          outputfilename_.str("");
          outputfilename_ << "EndToEnd"<< testnumber << ".dat";
        }
         
        //increasing the number in the numbering file by one
        fpnumbering = fopen(numberingfilename.str().c_str(), "w");
        std::stringstream filecontent;
        filecontent << "Next Number: "<< testnumber + 1;
        fprintf(fpnumbering,filecontent.str().c_str());
        fclose(fpnumbering);
        
      }
      if (time > starttimeoutput_ && starttimeoutput_ > -1.0)
      { 
        endtoend = pow ( pow((dis)[num_dof-3]+10 - (dis)[0],2) + pow((dis)[num_dof-2] - (dis)[1],2) , 0.5);
        //applying in the following a well conditioned substraction formula according to Crisfield, Vol. 1, equ. (7.53)
        DeltaR2 = pow( (endtoend*endtoend - endtoendref_*endtoendref_) / (endtoend + endtoendref_) ,2 );
    
        //writing output: writing Delta(R^2) according to PhD thesis Hallatschek, eq. (4.60), where t=0 corresponds to starttimeoutput_
        if ( (istep - istart_) % int(ceil(pow( 10, floor(log10((time - starttimeoutput_) / (10*dt))))) ) == 0 )
        {
          // open file and append new data line
          fp = fopen(outputfilename_.str().c_str(), "a");
          //defining temporary stringstream variable
          std::stringstream filecontent;
          filecontent << scientific << setprecision(15) << time - starttimeoutput_ << "  " << DeltaR2 << endl;
          // move temporary stringstream to file and close it
          fprintf(fp,filecontent.str().c_str());
          fclose(fp);
        }
      }
    }
    break;
    case INPUTPARAMS::statout_endtoendergodicity:
    {
      FILE* fp = NULL; //file pointer for statistical output file
      double endtoend = 0; //end to end length at a certain time step in single filament dynamics
       
      //as soon as system is equilibrated (after time START_FACTOR*maxtime_) a new file for storing output is generated
      if ( (time >= maxtime_ * statmechparams_.get<double>("START_FACTOR",0.0))  && (starttimeoutput_ == -1.0) )
      {
       endtoendref_ = pow ( pow((dis)[num_dof-3]+10 - (dis)[0],2) + pow((dis)[num_dof-2] - (dis)[1],2) , 0.5);  
       starttimeoutput_ = time;
       istart_ = istep;   
       //setting up an empty file "EndToEndErgo.dat"
       fp = fopen("EndToEndErgo.dat", "w");
       fclose(fp);
      }
       
        if (time > starttimeoutput_ && starttimeoutput_ > -1.0)
        { 
          endtoend = pow ( pow((dis)[num_dof-3]+10 - (dis)[0],2) + pow((dis)[num_dof-2] - (dis)[1],2) , 0.5);
      
          //writing output: current time and end to end distance are stored without any further postprocessing
          if ( (istep - istart_) % 100 == 0 )
            {
            // open file and append new data line
            fp = fopen("EndToEndErgo.dat", "a");
            //defining temporary stringstream variable
            std::stringstream filecontent;
            filecontent << scientific << setprecision(15) << time << "  " << endtoend << endl;
            // move temporary stringstream to file and close it
            fprintf(fp,filecontent.str().c_str());
            fclose(fp);
            }
        }
    }
    break;
    case INPUTPARAMS::statout_none:
    case INPUTPARAMS::statout_viscoelasticity:
    default:
    {
      std::cout<<"\nno special output written\n";
    }
    break;
  }
 
  return;
} // StatMechManager::StatMechOutput()


/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechUpdate()
{  
  //if dynamic crosslinkers are used update comprises adding and deleting crosslinkers
  if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
  {
    //create random numbers in order to decide whether a crosslinker should be established
    RCP<Epetra_Vector>  setcrosslinker;
    setcrosslinker = LINALG::CreateVector(*discret_.NodeRowMap(),true);
    setcrosslinker->Random();
     
    //create random numbers in order to decide whether a crosslinker should be deleted   
    RCP<Epetra_Vector>  delcrosslinker;
    delcrosslinker = LINALG::CreateVector(*discret_.NodeRowMap(),true);
    delcrosslinker->Random();
    
    //probability with which a crosslinker is established between neighbouring nodes
    double plink = 1;
    double punlink = 1;
    
     
    //searching locally with a tree all nodes forming neighbouring couples and establishing crosslinkers between them
    for(int i = 0; i < discret_.NumMyRowNodes(); i++)
    {
      //if the node is already crosslinked the crosslinker is removed with probability punlink
      if ((*crosslinkerpartner_)[i] != -1.0 && (*delcrosslinker)[i] <  -1.0 + 2*punlink)
      {
        //noting that the Id of the now deleted crosslinker will be unused from now on
        unusedids_.push_back((int)(*crosslinkerelements_)[i] );
        //noting that node responsible for the crosslinker to be deleted has from now on no crosslinker
        (*crosslinkerpartner_)[i] = -1;
        //deleting the crosslinker element number
        (*crosslinkerelements_)[i] = -1;
        
        //delete crosslinker element
        
        //???????????????????????????????????????????????????????????????????????????????????????????????????????
        
        
      }
      
      //searching neighbour for the current node (with global Id "neighbour"):
      int neighbour = -1;
      
      //???????????????????????????????????????????????????????????????????????????????????????????????????????
      
      
  
      
      //defining global Id of new crosslinker element
      int crosslinkerid = 0;
      //if there is a formerly used, but now unused global Id the new crosslinker gets this Id:
      if(unusedids_.size() > 0)
      {
        crosslinkerid = unusedids_.back();
        unusedids_.pop_back();     
      }
      //otherwise the new crosslinker element gets the lowest global Id currently not yet in use
      else
      {
        /*note: if there are no gaps (i.e. unused numbers) in the global Id distribution the highest
         *global element Id currently in use is ( discret_.NumGlobalElements() - 1 ) */    
        crosslinkerid = discret_.NumGlobalElements();
      }
    
      /*crosslinkers are always established from the node with lower global Id to the one with higher
       * global Id only, in order to make sure that the crosslinker is not established twice; furthermore
       * crosslinkers are established only with a certain probability, which is accounted for by menas
       * of the randomized condition (*setcrosslinker)[i] <  -1.0 + 2*plink; finally a crosslinker is established
       * only if the same node has not already one, which is checked by (*crosslinkerpartner_)[i] == -1.0*/
      if (discret_.NodeRowMap()->GID(i) < neighbour && (*setcrosslinker)[i] <  -1.0 + 2*plink && (*crosslinkerpartner_)[i] == -1.0)
      {
        
        /*the owner of the newly established crosslinker element is the current processor on which the search
         * is carried out and whose processor Id can be requested by the command discret_.Comm().MyPID(); the
         * global Id of the new element is chosen appropirately by a special algorithm*/     
        discret_.AddElement( DRT::UTILS::Factory("Beam3","Polynomial",crosslinkerid,discret_.Comm().MyPID()) );    
        
        //noting that local node i has no a crosslinkerelement and noting global Id of this element
        (*crosslinkerpartner_)[neighbour] = -1;
        (*crosslinkerelements_)[crosslinkerid] = -1;
      }
    }
    
    //settling administrative stuff in order to make the discretization ready for the next time step
    discret_.FillComplete(true,true,true);
    
    
    //???????????????????????????????????????????????????????????????????????????????????????????????????????
 
  }
  return;
} // StatMechManager::StatMechUpdate()

/*----------------------------------------------------------------------*
 | updates system damping matrix and external force vector according to |
 | influence of thermal bath (public)                        cyron 10/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechForceDamp(ParameterList& params, RCP<Epetra_Vector> dis, RCP<Epetra_Vector> fext, RCP<LINALG::SparseMatrix> damp)
{  
  // zero out damping matrix
  damp->Zero();
  
  /*declaration of a column and row map Epetra_Vector for evaluation of statistical forces; note: zero initilization 
   * mandatory for correct computations later on*/
  RCP<Epetra_Vector>    fstatcol;
  fstatcol = LINALG::CreateVector(*discret_.DofColMap(),true);
  RCP<Epetra_Vector>    fstatrow;
  fstatrow = LINALG::CreateVector(*discret_.DofRowMap(),true);
  
  //defining parameter list passed down to the elements in order to evalute statistical forces down there
  ParameterList pstat;
  pstat.set("action","calc_stat_force_damp");
  pstat.set("delta time",params.get<double>("delta time",0.0));
  pstat.set("KT",statmechparams_.get<double>("KT",0.0));
  pstat.set("ETA",statmechparams_.get<double>("ETA",0.0));
  pstat.set("STOCH_ORDER",statmechparams_.get<int>("STOCH_ORDER",0));
  
  /*note: the column map statistical force vector is passed down via the parameter list and not as a systemvector 
   * so that assembly is not done by the evaluate method itself, but elementwise; this is in order to account for the
   * special assembly needs of randomly evaluated forces: the evaluate method of the discretization uses the LINALG
   * assembly method in which a processor assembles to the element of the global vector only ir the processor is the 
   * row owner of the related DOF; this is efficient for global row map vectors, but does not assemble correctly if a 
   * global column map vector is used*/
  pstat.set("statistical force vector",fstatcol);
  
  //evaluation of statistical forces on column map vecotor
  discret_.SetState("displacement",dis); //during evaluation of statistical forces access to current displacement possible
  discret_.Evaluate(pstat,damp,null,null,null,null);
  discret_.ClearState();

  
  /*exporting col map statistical force vector to a row map vector additively, i.e. in such a way that a 
   * vector element with a certain GID in the final row vector is the sum of all elements of the column 
   * vector related to the same GID*/
  Epetra_Export exporter(*discret_.DofColMap(),*discret_.DofRowMap());
  fstatrow->Export(*fstatcol,exporter,Add);
  
  //adding statistical forces to external forces passed to this method
  fext->Update(1.0,*fstatrow,1.0);
  
  //complete damping matrix
  damp->Complete();

  return;
} // StatMechManager::StatMechForceDamp()




#endif  // #ifdef CCADISCRET
