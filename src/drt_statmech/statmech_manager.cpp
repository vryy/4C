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

#include <iostream>


/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 09/08|
 *----------------------------------------------------------------------*/
StatMechManager::StatMechManager(ParameterList& params):
  statmechparams_( DRT::Problem::Instance()->StatisticalMechanicsParams() ),
  maxtime_(params.get<double>("max time",0.0)),
  starttimeoutput_(-1.0),
  endtoendref_(0.0),
  istart_(0)
{ 
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
       //look for a not yet existing output file name (numbering upwards)
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
       
        if (time > starttimeoutput_ && starttimeoutput_ > -1.0)
        { 
          endtoend = pow ( pow((dis)[num_dof-3]+10 - (dis)[0],2) + pow((dis)[num_dof-2] - (dis)[1],2) , 0.5);
          //applying in the following a well conditioned substraction formula according to Crisfield, Vol. 1, equ. (7.53)
          DeltaR2 = pow( (endtoend*endtoend - endtoendref_*endtoendref_) / (endtoend + endtoendref_) ,2 );
      
          //writing output
          if ( (istep - istart_) % int(ceil(pow( 10, floor(log10((time - starttimeoutput_) / (10*dt))))) ) == 0 )
            {
            // open file and append new data line
            fp = fopen(outputfilename_.str().c_str(), "a");
            //defining temporary stringstream variable
            std::stringstream filecontent;
            filecontent << scientific << time - starttimeoutput_ << "  " << DeltaR2 << endl;
            // move temporary stringstream to file and close it
            fprintf(fp,filecontent.str().c_str());
            fclose(fp);
            }
        }
    }
    break;
    case INPUTPARAMS::statout_none:
    case INPUTPARAMS::statout_endtoendergodicity:
    case INPUTPARAMS::statout_viscoelasticity:
    default:
    {
      std::cout<<"\nno special output written\n";
    }
    break;
  }
  

  return;
} // StatMechManager::StatMechOutput()


#endif  // #ifdef CCADISCRET
