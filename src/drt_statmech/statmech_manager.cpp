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
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

//including random number library of blitz for statistical forces
#include <random/uniform.h>


#ifdef D_BEAM3
#include "../drt_beam3/beam3.H"
#endif  // #ifdef D_BEAM3
#ifdef D_TRUSS3
#include "../drt_truss3/truss3.H"
#endif  // #ifdef D_TRUSS3

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <math.h>

//MEASURETIME activates measurement of computation time for certain parts of the code
//#define MEASURETIME

/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 09/08|
 *----------------------------------------------------------------------*/
StatMechManager::StatMechManager(ParameterList& params, DRT::Discretization& discret, RCP<LINALG::SparseMatrix>& stiff):
  statmechparams_( DRT::Problem::Instance()->StatisticalMechanicsParams() ),
  maxtime_(params.get<double>("max time",0.0)),
  starttimeoutput_(-1.0),
  endtoendref_(0.0),
  istart_(0),
  rlink_(params.get<double>("R_LINK",0.0)),
  basisnodes_(0),
  outputfilenumber_(-1),
  discret_(discret),
  stiff_(stiff)
{ 
  //initialize random generator on this processor
  
  //random generator for seeding only
  ranlib::UniformClosed<double> seedgenerator;
  //seeding random generator independently on each processor
  int seedvariable = time(0)*(discret_.Comm().MyPID() + 1);
  //seedvariable = 1;
  seedgenerator.seed((unsigned int)seedvariable);


  /*setting and deleting dynamic crosslinkers needs a special code structure for parallel search trees
   * here we provide a fully overlapping column map which is required by the search tree to look for each
   * node for all the neighbouring nodes within a certain distance; note: we pass this fully overlapping
   * map to the discretization and adapt it accordingly; this is not efficient in parallel use, however, it
   * makes sure that it leads to at least correct results; as a consequence this implementation may be con-
   * sidered an implementation capable for parallel use, but not optimized for it; note: the discretization's
   * column map is turned into a fully overlapping one even in case that dynamic crosslinkers are deactivated;
   * the reason is that also the method for Gmsh output currently relies on a fully overlapping column map
   * in case of parallel computing (otherwise the output is not written correctly*/
    
    const Epetra_Map noderowmap = *(discret_.NodeRowMap());

    // fill my own row node ids into vector sdata
    vector<int> sdata(noderowmap.NumMyElements());
    for (int i=0; i<noderowmap.NumMyElements(); ++i)
      sdata[i] = noderowmap.GID(i);
    
    //if current processor has elements it writes its number into stproc
    vector<int> stproc(0); 
    
    
    if (noderowmap.NumMyElements())
      stproc.push_back(discret_.Comm().MyPID());
    
    
    //information how many processor work at all
    vector<int> allproc(discret_.Comm().NumProc());
    
    //in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
    for (int i=0; i<discret_.Comm().NumProc(); ++i) allproc[i] = i;
    
    //declaring new variable into which the information of stproc on all processors is gathered
    vector<int> rtproc(0);
    
    /*gathers information of stproc and writes it into rtproc; in the end rtproc is a vector which 
     * contains the numbers of all processors which have elements*/
    LINALG::Gather<int>(stproc,rtproc,discret_.Comm().NumProc(),&allproc[0],discret_.Comm());
      
    /*in analogy to stproc and rtproc the variable rdata gathers all the element numbers which are 
     * stored on different processors in their own variables sdata; thereby each processor gets
     * the information about all the nodes numbers existing in this problem*/
    vector<int> rdata;

    // gather all gids of nodes redundantly from sdata into rdata
    LINALG::Gather<int>(sdata,rdata,(int)rtproc.size(),&rtproc[0],discret_.Comm());

    // build completely overlapping map (on participating processors)
    RCP<Epetra_Map> newnodecolmap = rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,discret_.Comm()));
    sdata.clear();
    stproc.clear();
    rdata.clear();
    allproc.clear();
    // rtproc still in use
    
    //pass new fully overlapping column map to discretization
    discret_.ExportColumnNodes(*newnodecolmap);

    /*rebuild discretization based on the new column map so that each processor creates new ghost elements
     * if necessary; after the following line we have a discretization, where each processor has a fully
     * overlapping column map regardlesse of how overlapping was managed when starting BACI; having ensured
     * this allows convenient and correct (albeit not necessarily efficient) use of search algorithms and
     * crosslinkers in parallel computing*/     
    discret_.FillComplete(true,false,false);
  
  
  //if dynamic crosslinkers are used additional variables are initialized
  if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
  {  
    /* and filamentnumber_ is generated based on a column map vector as each node has to 
     * know about each other node its filament number in order to decide weather a crosslink may be established 
     * or not; vectors is initalized with -1, which state is changed if filament numbering is used, only*/
    filamentnumber_ = rcp( new Epetra_Vector(*(discret_.NodeColMap())) );
    filamentnumber_->PutScalar(-1);
    
    //initialize crosslinkerneighbours_ and crosslinkerpartner with empty vector for each row node
    std::vector<int> emptyvector;
    crosslinkerneighbours_.resize(discret_.NumMyRowNodes(),emptyvector);
    crosslinkerpartner_.resize(discret_.NumMyRowNodes(),emptyvector);
    

    
    /*force sensors can be applied at any degree of freedom of the discretization the list of force sensors should
     * be based on a column map vector so that each processor has not only the information about each node's
     * displacement, but also about whether this has a force sensor; as a consequence each processor can write the
     * complete information gathered by all force sensors into a file of its own without any additional communication
     * with any other processor; initialization with -1 indicates that so far no forcesensors have been set*/
    forcesensor_ = rcp( new Epetra_Vector(*(discret_.DofColMap())) );
    forcesensor_->PutScalar(-1);
    
    //basisnodes_ is the number of nodes existing in the discretization from the very beginning on (should not change)
    basisnodes_ = discret_.NumGlobalNodes();

  
    /*since crosslinkers should be established only between different filaments the number of the filament
     * each node is belonging to is stored in the condition FilamentNumber; if no such conditions have been defined
     * the default -1 value is upkept in filamentnumber_ and crosslinkers between nodes belonging to the same filament
     * are allowed*/
    
    //gettin a vector consisting of pointers to all filament number conditions set
    vector<DRT::Condition*> filamentnumberconditions(0);
    discret_.GetCondition("FilamentNumber",filamentnumberconditions);
      
    //next all the pointers to all the different conditions are looped
    for (int i=0; i<(int)filamentnumberconditions.size(); ++i)
    {
      //get filament number described by the current condition
      int filamentnumber = filamentnumberconditions[i]->GetInt("Filament Number") ;
      
      //get a pointer to nodal cloud coverd by the current condition
      const vector<int>* nodeids = filamentnumberconditions[i]->Nodes();
      
      //loop through all the nodes of the nodal cloud
      for(int j = 0; j < (int)nodeids->size() ; j++)
      {
        //get the node's global id
        int nodenumber = (*nodeids)[j];
        
        //turning global id into local one
        nodenumber = discret_.NodeColMap()->LID(nodenumber);
        
        //if the node does not belong to current processor Id is set by LID() to -1 and the node is ignored
        if(nodenumber > -1)
          (*filamentnumber_)[nodenumber] = filamentnumber;     
      }
     }  
    
    /*Young's modulus and loss modulus are to be measured by means of the reaction forces at certain sensor points; example: if for a
     * an actin network between two rheometer plates the stiffness is to be deterimined this can be done by measuring the forces exerted
     * to the upper plate which is moving forwards and backwards for example; so for measurements of viscoelastic properties of materials
     * whithin a system certain sensor points have to be specified and in order to handle this whithing BACI these points are marked by 
     * means of the condition sensorcondition*/
    
    //gettin a vector consisting of pointers to all filament number conditions set
    vector<DRT::Condition*> forcesensorconditions(0);
    discret_.GetCondition("ForceSensor",forcesensorconditions);
      
    //next all the pointers to all the different conditions are looped
    for (int i=0; i<(int)forcesensorconditions.size(); ++i)
    {
      //get number of nodal dof with respect to which force is to be measured; note: numbering starts with zero
      int nodedofnumber = forcesensorconditions[i]->GetInt("DOF Number") ;
      
      //get a pointer to nodal cloud coverd by the current condition
      const vector<int>* nodeids = forcesensorconditions[i]->Nodes();
      
      //loop through all the nodes of the nodal cloud
      for(int j = 0; j < (int)nodeids->size() ; j++)
      {
        //get the node's global id
        int nodenumber = (*nodeids)[j];
        
        //testing whether current nodedofnumber makes sense for current node
        if( nodedofnumber < 0 || nodedofnumber >= discret_.NumDof(discret_.gNode(nodenumber)) )
          dserror("ForceSensor condition applied with improper local dof number");
        
        //global id of degree of freedom at which force is to be measured
        int dofnumber = discret_.Dof( discret_.gNode(nodenumber), nodedofnumber-1 );
        
        
        /*if the node does not belong to current processor Id is set by LID() to -1 and the node is ignored; otherwise the degrees of
         * freedom affected by this condition are marked in the vector *forcesensor_ by an one entry*/
        if(nodenumber > -1)
          (*forcesensor_)[dofnumber] = 1;   
      }
     } 
  
      //we construct a search tree with maximal depth of 8 (different numbers might be tried out, too)
      octTree_ = rcp(new GEO::SearchTree(8));
      
      /*after having generated a search tree and a discretization with fully overlapping column map we initialize the search tree
       * for accelerated search for each nodes neighbouring nodes; note: the tree is based on a bounding box
       * with respect to the reference positions (which are the positions at the beginning of the simulation;
       * in case of large overall deformations of the fiber network such an initialization would have to be carried 
       * out in each time step with respect to the current node positions*/

      /*currenpositions is a map which relates each LID of any node on this processor to the nodes coordinates*/
      std::map<int,LINALG::Matrix<3,1> > currentpositions;
      
      currentpositions.clear();

      for (int lid = 0; lid <discret_.NumMyColNodes(); ++lid)
      {
        const DRT::Node* node = discret_.lColNode(lid);
        LINALG::Matrix<3,1> currpos;
        currpos(0) = node->X()[0];
        currpos(1) = node->X()[1];
        currpos(2) = node->X()[2] ;
        currentpositions[node->LID()] = currpos;
      }
      
      //find bounding box for search tree
      const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(discret_, currentpositions);
   
      //initialize search tree
      octTree_->initializePointTree(rootBox,currentpositions,GEO::TreeType(GEO::OCTTREE));
      
      /*we finish the initilization by once determining for each row map node on each processor
       * neighbour nodes; as this search is preformed only once in the beginning, this algorithm
       * is suitable for small network deformation only currently*/
#ifdef MEASURETIME
    const double t_search = ds_cputime();
#endif // #ifdef MEASURETIME   
    
    SearchNeighbours(currentpositions);
      
#ifdef MEASURETIME
    cout << "\n***\nsearch time: " << ds_cputime() - t_search<< " seconds\n***\n";
#endif // #ifdef MEASURETIME
    
  }//if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))

  return;
} // StatMechManager::StatMechManager

/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechOutput(ParameterList& params, const int ndim, const double& time,const int& istep, const double& dt, const Epetra_Vector& dis, const Epetra_Vector& fint)
{
  /*in general simulations in statistical mechanics run over so many time steps that the amount of data stored in the error file
   * may exceed the capacity even of a server hard disk; thus, we rewind the error file in each time step so that the amount of data
   * does not increase after the first time step any longer*/  
  bool printerr    = params.get<bool>("print to err",false);
  FILE* errfile    = params.get<FILE*>("err file",NULL);
  if(printerr)
    rewind(errfile);
  
  
  //the following variable makes sense in case of serial computing only; its use is not allowed for parallel computing!
  int num_dof = dis.GlobalLength();
  

  
  switch(Teuchos::getIntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_,"SPECIAL_OUTPUT"))
  {
    case INPAR::STATMECH::statout_endtoendlog:
    {      
      FILE* fp = NULL; //file pointer for statistical output file
      
      //name of output file
      std::ostringstream outputfilename;
      outputfilename << "EndToEnd"<< outputfilenumber_ << ".dat";
      
      double endtoend = 0; //end to end length at a certain time step in single filament dynamics
      double DeltaR2 = 0;
       
      //as soon as system is equilibrated (after time START_FACTOR*maxtime_) a new file for storing output is generated
      if ( (time >= maxtime_ * statmechparams_.get<double>("START_FACTOR",0.0))  && (starttimeoutput_ == -1.0) )
      {
         endtoendref_ = pow ( pow((dis)[num_dof-3]+10 - (dis)[0],2) + pow((dis)[num_dof-2] - (dis)[1],2) , 0.5);  
         starttimeoutput_ = time;
         istart_ = istep;      
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
          fp = fopen(outputfilename.str().c_str(), "a");
          
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
    case INPAR::STATMECH::statout_endtoendconst:
    {          
      /*in the following we assume that there is only a pulling force point Neumann condition of equal absolute
      *value on either filament length; we get the absolute value of the first one of these two conditions */
      double neumannforce;
      vector<DRT::Condition*> pointneumannconditions(0);
      discret_.GetCondition("PointNeumann",pointneumannconditions);     
      if(pointneumannconditions.size() > 0)
      {
        const vector<double>* val  = pointneumannconditions[0]->Get<vector<double> >("val");
        neumannforce = fabs((*val)[0]);
      }
      else
       neumannforce = 0;

   
      FILE* fp = NULL; //file pointer for statistical output file
         
      //name of output file
      std::ostringstream outputfilename;
      outputfilename << "E2E_"<< discret_.NumMyRowElements() <<'_'<< dt <<'_'<< neumannforce << '_'<< outputfilenumber_ << ".dat";
     
      double endtoend = 0.0; //end to end length at a certain time step in single filament dynamics
           
      //as soon as system is equilibrated (after time START_FACTOR*maxtime_) a new file for storing output is generated
      if ( (time >= maxtime_ * statmechparams_.get<double>("START_FACTOR",0.0))  && (starttimeoutput_ == -1.0) )
      {       
       starttimeoutput_ = time;
       istart_ = istep;   
      }
           
      if (time > starttimeoutput_ && starttimeoutput_ > -1.0)
      { 
        //end to end vector
        LINALG::Matrix<3,1> endtoendvector(true);     
        for(int i = 0; i<ndim; i++)
        {
          endtoendvector(i) -= ( discret_.gNode(0)                            )->X()[i] + dis[i];
          endtoendvector(i) += ( discret_.gNode(discret_.NumMyRowNodes() - 1) )->X()[i] + dis[num_dof - discret_.NumDof(discret_.gNode(discret_.NumMyRowNodes() - 1)) + i];
        }
        
        endtoend = endtoendvector.Norm2();
          
        //writing output: current time and end to end distance are stored at each 100th time step
        if ( (istep - istart_) % 100 == 0 )
          {
          // open file and append new data line
          fp = fopen(outputfilename.str().c_str(), "a");
          //defining temporary stringstream variable
          std::stringstream filecontent;
          filecontent << scientific << setprecision(15) << time << "  " << endtoend << " " << fint[num_dof - discret_.NumDof(discret_.gNode(discret_.NumMyRowNodes() - 1))] << endl;
          // move temporary stringstream to file and close it
          fprintf(fp,filecontent.str().c_str());
          fclose(fp);
          }
      }
    }
    break;
    //the following output allows for anisotropic diffusion simulation of a quasi stiff polymer
    case INPAR::STATMECH::statout_anisotropic:
    {            
      FILE* fp = NULL; //file pointer for statistical output file
      
      //name of output file
      std::ostringstream outputfilename;
      outputfilename << "AnisotropicDiffusion"<< outputfilenumber_ << ".dat";
              
      //positions of first and last node in current time step (note: always 3D variables used; in case of 2D third compoenent just constantly zero)   
      LINALG::Matrix<3,1> beginnew;
      LINALG::Matrix<3,1> endnew;
      
      beginnew.PutScalar(0);
      endnew.PutScalar(0);
      std::cout << "ndim: " << ndim << "\n";
      
      for(int i = 0; i<ndim; i++)
      {
        beginnew(i) = ( discret_.gNode(0)                            )->X()[i] + dis[i];
        endnew(i)   = ( discret_.gNode(discret_.NumMyRowNodes() - 1) )->X()[i] + dis[num_dof - discret_.NumDof(discret_.gNode(discret_.NumMyRowNodes() - 1)) + i];
      }
            
      //unit direction vector for filament axis in last time step
      LINALG::Matrix<3,1> axisold;
      axisold  = endold_;
      axisold -= beginold_;
      axisold.Scale(1/axisold.Norm2());
      
      //unit direction vector for filament axis in current time step
      LINALG::Matrix<3,1> axisnew;
      axisnew = endnew;
      axisnew -= beginnew;      
      axisnew.Scale(1/axisnew.Norm2());
     
      //displacement of first and last node between last time step and current time step
      LINALG::Matrix<3,1> dispbegin;
      LINALG::Matrix<3,1> dispend;
      dispbegin  = beginnew;
      dispbegin -= beginold_;
      dispend  = endnew;
      dispend -= endold_;
        
      //displacement of middle point
      LINALG::Matrix<3,1> dispmiddle;
      dispmiddle  = dispbegin;
      dispmiddle += dispend;
      dispmiddle.Scale(0.5);
      sumdispmiddle_ += dispmiddle;
   
      //update sum of square displacement increments of middle point
      double incdispmiddle = dispmiddle.Norm2()*dispmiddle.Norm2();      
      sumsquareincmid_+=incdispmiddle;
   
      //update sum of square displacement increments of middle point parallel to new filament axis (computed by scalar product)    
      double disppar_square = pow(axisnew(0)*dispmiddle(0) + axisnew(1)*dispmiddle(1) + axisnew(2)*dispmiddle(2), 2);       
      sumsquareincpar_+=disppar_square;
      
      //update sum of square displacement increments of middle point orthogonal to new filament axis (via crossproduct)
      LINALG::Matrix<3,1> aux;
      aux(0) = dispmiddle(1)*axisnew(2) - dispmiddle(2)*axisnew(1);
      aux(1) = dispmiddle(2)*axisnew(0) - dispmiddle(0)*axisnew(2);
      aux(2) = dispmiddle(0)*axisnew(1) - dispmiddle(1)*axisnew(0);
      double disport_square = aux.Norm2()*aux.Norm2();
      sumsquareincort_+=disport_square;
      
      
      //total displacement of rotational angle (in 2D only)
      double incangle = 0;
      if(ndim == 2)
      {
        //angle of old axis relative to x-axis
        double phiold = acos(axisold(0)/axisold.Norm2());
        if(axisold(1) < 0)
          phiold *= -1;
        
        //angle of new axis relative to x-axis
        double phinew = acos(axisnew(0)/axisnew.Norm2());
        if(axisnew(1) < 0)
          phinew *= -1;
        
        //angle increment
        incangle = phinew - phiold;
        if(incangle > PI)
        {
          incangle -= 2*PI;
          incangle *= -1;
        }
        if(incangle < -PI)
        {
          incangle += 2*PI;
          incangle *= -1;
        }
        
        //update absolute rotational displacement compared to reference configuration
        sumsquareincrot_ += incangle*incangle;
        sumrotmiddle_ += incangle;
      }
      

      {
	      // open file and append new data line
	      fp = fopen(outputfilename.str().c_str(), "a");    
	    
	      //defining temporary stringstream variable
	      std::stringstream filecontent;
	      filecontent << scientific << setprecision(15)<<dt<<" "<<sumsquareincmid_ <<" "<<sumsquareincpar_<<" "<<sumsquareincort_<<" "<<sumsquareincrot_<<" "<<sumdispmiddle_.Norm2() * sumdispmiddle_.Norm2()<<" "<<sumrotmiddle_*sumrotmiddle_ << endl;
	    
	      // move temporary stringstream to file and close it
	      fprintf(fp,filecontent.str().c_str());
	      fclose(fp);
      }
       
      //new positions in this time step become old positions in last time step     
      beginold_ = beginnew;
      endold_   = endnew;
    }
    break;
    case INPAR::STATMECH::statout_viscoelasticity:
    {   
       //pointer to file into which each processor writes the output related with the dof of which it is the row map owner
       FILE* fp = NULL; 
            
       //content to be written into the output file
       std::stringstream filecontent;
       
       //name of file into which output is written
       std::ostringstream outputfilename; 
       outputfilename << "ViscoElOutputProc"<< discret_.Comm().MyPID() << ".dat";
      
       fp = fopen(outputfilename.str().c_str(), "a");
       
       /*the output to be written consists of internal forces at exactly those degrees of freedom
        * marked in *forcesensor_ by a one entry*/
       
       filecontent << scientific << setprecision(10) << time;//changed
       
 #ifdef DEBUG
       if(forcesensor_ == null)
         dserror("forcesensor_ is NULL pointer; possible reason: dynamic crosslinkers not activated and forcesensor applicable in this case only");
 #endif  // #ifdef DEBUG
       
        double f = 0;//mean value of force
        double d = 0;//Displacement
        for(int i = 0; i < forcesensor_->MyLength(); i++)//changed
        {
          if( (*forcesensor_)[i] == 1)
             {
               f += fint[i];
               d =  dis[i];          
             }
         }
        
        //Putting time, displacement, meanforce  in Filestream
        filecontent << "   "<< d << "   " << f << "   " << discret_.NumMyRowElements() << endl; //changed
        //writing filecontent into output file and closing it
        fprintf(fp,filecontent.str().c_str());
        fclose(fp);

    }
    break;
    //writing data for generating a Gmsh video of the simulation
    case INPAR::STATMECH::statout_gmsh:
    {
      
      /*construct unique filename for gmsh output with two indices: the first one marking the time step number
       * and the second one marking the newton iteration number, where numbers are written with zeros in the front
       * e.g. number one is written as 000001, number fourteen as 000014 and so on;*/
      
      //note: this kind of output is possilbe for serial computing only (otherwise the following method would have to be adapted to parallel use*/   
      if(discret_.Comm().NumProc() > 1)
        dserror("No Gmsh output for parallel computation possible so far");
      
      // first index = time step index
      std::ostringstream filename;

      //creating complete file name dependent on step number with 6 digits and leading zeros
      if (istep<1000000)
        filename << "./GmshOutput/network"<< std::setw(6) << setfill('0') << istep <<".pos";
      else 
        dserror("Gmsh output implemented for a maximum of 999999 steps");
          
      //calling method for writing Gmsh output
      GmshOutput(dis,filename,istep);     
   
    }
    break;
    case INPAR::STATMECH::statout_none:
    default:
    break;
  }
 
  return;
} // StatMechManager::StatMechOutput()


/*----------------------------------------------------------------------*
 | writing Gmsh data for current step                 public)cyron 01/09|
 *----------------------------------------------------------------------*/
void StatMechManager::GmshOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step)
{
  /*the following method writes output data for Gmsh into file with name "filename"; all line elements are written;
   * the nodal displacements are handed over in the variable "dis"; note: in case of parallel computing only 
   * processor 0 writes; it is assumed to have a fully overlapping column map and hence all the information about
   * all the nodal position*/
  if(discret_.Comm().MyPID() == 0)
  {
  
    //we need displacements also of ghost nodes and hence export displacment vector to column map format
    Epetra_Vector  discol(*(discret_.DofColMap()),true);   
    LINALG::Export(disrow,discol);

    // do output to file in c-style
    FILE* fp = NULL;
    
    //open file to write output data into
    fp = fopen(filename.str().c_str(), "w");
  
    // write output to temporary stringstream; 
    std::stringstream gmshfilecontent;
    /*the beginning of the stream is "View \"" to indicate Gmsh that the following data is in order to create an image and
     * this command is followed by the name of that view displayed during it's shown in the video; in the following example
     * this name is for the 100th time step: Step00100; then the data to be presented within this view is written within { ... };
     * in the following example this data consists of scalar lines defined by the coordinates of their end points*/
    gmshfilecontent << "View \" Step " << step << " \" {" << endl;
  
    //looping through all elements on the processor
    for (int i=0; i<discret_.NumMyColElements(); ++i)
    {
      //getting pointer to current element
      DRT::Element* element = discret_.lColElement(i);
      
      //getting number of nodes of current element
      if( element->NumNode() > 2)
        dserror("Gmsh output for two noded elements only");
      
      //preparing variable storing coordinates of all these nodes
      int nnodes = 2;
      LINALG::SerialDenseMatrix coord(3,nnodes);
      
      for(int id = 0; id<3; id++)
       for(int jd = 0; jd<2; jd++)  
       {
         double referenceposition = ((element->Nodes())[jd])->X()[id];
         vector<int> dofnode = discret_.Dof((element->Nodes())[jd]);        
         double displacement = discol[discret_.DofColMap()->LID( dofnode[id] )];
         coord(id,jd) =  referenceposition + displacement;
       }
      
            
      //declaring variable for color of elements
      double color;
      
      //apply different colors for elements representing filaments and those representing dynamics crosslinkers
      if (element->Id() < basisnodes_)
        color = 1.0;
      else
        color = 0.0;
  
      //writing element by nodal coordinates as a scalar line
      gmshfilecontent << "SL(" << scientific;
      gmshfilecontent<< coord(0,0) << "," << coord(1,0) << "," << coord(2,0) << "," 
                     << coord(0,1) << "," << coord(1,1) << "," << coord(2,1) ;
      /*note: for each node there is one color variable for gmsh and gmsh finally plots the line
       * interpolating these two colors between the nodes*/
      gmshfilecontent << ")" << "{" << scientific << color << "," << color << "};" << endl;
  
    }  
    
    //finish data section of this view by closing curley brackets
    gmshfilecontent << "};" << endl;
    
    //write content into file and close it
    fprintf(fp,gmshfilecontent.str().c_str());
    fclose(fp);
    
  }//if(discret_.Comm().MyPID() == 0)
 
  return;
} // StatMechManager::GmshOutput()


/*----------------------------------------------------------------------*
 | initialize special output for statistical mechanics(public)cyron 12/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechInitOutput(const int ndim,const double& dt)
{
  //initializing special output for statistical mechanics by looking for a suitable name of the outputfile and setting up an empty file with this name
  
  switch(Teuchos::getIntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_,"SPECIAL_OUTPUT"))
  {
    case INPAR::STATMECH::statout_endtoendlog:
    {
       FILE* fp = NULL; //file pointer for statistical output file
       
       //defining name of output file
       std::ostringstream outputfilename;

       outputfilenumber_ = 0; 

       //file pointer for operating with numbering file
       FILE* fpnumbering = NULL;
       std::ostringstream numberingfilename;
       
       //look for a numbering file where number of already existing output files is stored:
       numberingfilename.str("NumberOfRealizationsLog");
       fpnumbering = fopen(numberingfilename.str().c_str(), "r");
       
       //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
       if(fpnumbering == NULL)
       {
         do
         {
           outputfilenumber_++;
           outputfilename.str("");
           outputfilename << "EndToEnd"<< outputfilenumber_ << ".dat";
           fp = fopen(outputfilename.str().c_str(), "r");
         } while(fp != NULL);
         
         //set up new file with name "outputfilename" without writing anything into this file
         fp = fopen(outputfilename.str().c_str(), "w");
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
              f >> outputfilenumber_;
          }
        } //while(f)
        
        //defining outputfilename by means of new testnumber
        outputfilename.str("");
        outputfilename << "EndToEnd"<< outputfilenumber_ << ".dat";
      }
       
      //increasing the number in the numbering file by one
      fpnumbering = fopen(numberingfilename.str().c_str(), "w");
      std::stringstream filecontent;
      filecontent << "Next Number: "<< (outputfilenumber_ + 1);
      fprintf(fpnumbering,filecontent.str().c_str());
      fclose(fpnumbering);

    }
    break;
    case INPAR::STATMECH::statout_endtoendconst:
    {
      FILE* fp = NULL; //file pointer for statistical output file
             
       //defining name of output file
       std::ostringstream outputfilename;
       outputfilename.str("");
       outputfilenumber_ = 0; 

       //file pointer for operating with numbering file
       FILE* fpnumbering = NULL;
       std::ostringstream numberingfilename;
       
       //look for a numbering file where number of already existing output files is stored:
       numberingfilename.str("NumberOfRealizationsConst");
       fpnumbering = fopen(numberingfilename.str().c_str(), "r");
       
       /*in the following we assume that there is only a pulling force point Neumann condition of equal absolute
       *value on either filament length; we get the absolute value of the first one of these two conditions */
       double neumannforce;
       vector<DRT::Condition*> pointneumannconditions(0);
       discret_.GetCondition("PointNeumann",pointneumannconditions);     
       if(pointneumannconditions.size() > 0)
       {
         const vector<double>* val  =pointneumannconditions[0]->Get<vector<double> >("val");
         neumannforce = fabs((*val)[0]);
       }
       else
        neumannforce = 0;
       
       //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
       if(fpnumbering == NULL)
       {
         do
         {              
          //defining name of output file        
           outputfilenumber_++;
           outputfilename.str("");
           outputfilename << "E2E_"<< discret_.NumMyRowElements() <<'_'<< dt <<'_'<< neumannforce << '_'<< outputfilenumber_ << ".dat";
           fp = fopen(outputfilename.str().c_str(), "r");
         } while(fp != NULL);
         
         //set up new file with name "outputfilename" without writing anything into this file
         fp = fopen(outputfilename.str().c_str(), "w");
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
              f >> outputfilenumber_;
          }
        } //while(f)
        
        //defining outputfilename by means of new testnumber
        outputfilename.str("");
        outputfilename << "E2E_"<< discret_.NumMyRowElements() <<'_'<< dt <<'_'<< neumannforce << '_'<<outputfilenumber_ << ".dat";
        
        //set up new file with name "outputfilename" without writing anything into this file
        fp = fopen(outputfilename.str().c_str(), "w");
        fclose(fp);
      }

       
      //increasing the number in the numbering file by one
      fpnumbering = fopen(numberingfilename.str().c_str(), "w");
      std::stringstream filecontent;
      filecontent << "Next Number: "<< (outputfilenumber_ + 1);
      //write fileconent into file!
      fprintf(fpnumbering,filecontent.str().c_str());
      //close file
      fclose(fpnumbering);
      
    }
    break;
    //simulating diffusion coefficient for anisotropic friction
    case INPAR::STATMECH::statout_anisotropic:
    {     
       FILE* fp = NULL; //file pointer for statistical output file
         
       //defining name of output file
       std::ostringstream outputfilename;
  
       outputfilenumber_ = 0; 
  
       //file pointer for operating with numbering file
       FILE* fpnumbering = NULL;
       std::ostringstream numberingfilename;
       
       //look for a numbering file where number of already existing output files is stored:
       numberingfilename.str("NumberOfRealizationsAniso");
       fpnumbering = fopen(numberingfilename.str().c_str(), "r");
       
       //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
       if(fpnumbering == NULL)
       {
         do
         {
           outputfilenumber_++;
           outputfilename.str("");
           outputfilename << "AnisotropicDiffusion"<< outputfilenumber_ << ".dat";
           fp = fopen(outputfilename.str().c_str(), "r");
         } while(fp != NULL);
         
         //set up new file with name "outputfilename" without writing anything into this file
         fp = fopen(outputfilename.str().c_str(), "w");
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
              f >> outputfilenumber_;
          }
        } //while(f)
        
        //defining outputfilename by means of new testnumber
        outputfilename.str("");
        outputfilename <<"AnisotropicDiffusion"<< outputfilenumber_ << ".dat";
      }
       
      //increasing the number in the numbering file by one
      fpnumbering = fopen(numberingfilename.str().c_str(), "w");
      std::stringstream filecontent;
      filecontent << "Next Number: "<< (outputfilenumber_ + 1);
      //write fileconent into file!
      fprintf(fpnumbering,filecontent.str().c_str());
      //close file
      fclose(fpnumbering);
  
      //initializing variables for positions of first and last node at the beginning
      beginold_.PutScalar(0);
      endold_.PutScalar(0);
      for(int i = 0; i<ndim; i++)
       {
         beginold_(i) = ( discret_.gNode(0)                            )->X()[i];
         endold_(i)   = ( discret_.gNode(discret_.NumMyRowNodes() - 1) )->X()[i];
       }  
      
      for (int i=0; i<3; i++)
       sumdispmiddle_(i,0)=0.0;
      
      sumsquareincpar_=0.0;
      sumsquareincort_=0.0;
      sumrotmiddle_=0.0;
      sumsquareincmid_=0.0;
      sumsquareincrot_=0.0;

      
    }
    break;
    case INPAR::STATMECH::statout_viscoelasticity:
    {
      //pointer to file into which each processor writes the output related with the dof of which it is the row map owner
      FILE* fp = NULL; 
          
      //content to be written into the output file
      std::stringstream filecontent;
      
      //defining name of output file related to processor Id
      std::ostringstream outputfilename; 
      outputfilename.str("");
      outputfilename << "ViscoElOutputProc"<< discret_.Comm().MyPID() << ".dat";
    
      fp = fopen(outputfilename.str().c_str(), "w");
        
      //filecontent << "Output for measurement of viscoelastic properties written by processor "<< discret_.Comm().MyPID() << endl;
         
      // move temporary stringstream to file and close it
      fprintf(fp,filecontent.str().c_str());
      fclose(fp);
    }
    break;
    case INPAR::STATMECH::statout_none:
    default:
    break;
  }
 
  return;
} // StatMechManager::StatMechInitOutput()


/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechUpdate(const double dt, const Epetra_Vector& disrow)
{  
  #ifdef D_BEAM3
  
#ifdef MEASURETIME
    const double t_start = ds_cputime();
#endif // #ifdef MEASURETIME
  
  //if dynamic crosslinkers are used update comprises adding and deleting crosslinkers
  if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
  {
    /*the following tow rcp pointers are auxiliary variables which are needed in order provide in the very end of the
     * crosslinker administration a node row and column map; these maps have to be taken here before the first modification
     * by deleting and adding elements have been carried out with the discretization since after such modifications the maps
     * cannot be called from the discretization before calling FillComplete() again which should be done only in the very end
     * note: this way of getting node maps is all right only if no nodes are added/ deleted, but elements only*/
    const Epetra_Map noderowmap = *discret_.NodeRowMap();
    const Epetra_Map nodecolmap = *discret_.NodeColMap(); 
    
    /*in preparation for later decision whether a crosslink should be established between two nodes we first store the 
     * current positions of all column map nodes in the map currentpositions; additionally we store the roational displacments
     * analogously in a map currentrotations for later use in setting up reference geometry of crosslinkers; the maps
     * currentpositions and currentrotations relate positions and rotations to a local column map node Id, respectively*/
    std::map<int,LINALG::Matrix<3,1> > currentpositions;   
    std::map<int,LINALG::Matrix<3,1> > currentrotations; 
    currentpositions.clear();
    currentrotations.clear();
    
    /*note: access by ExtractMyValues requires column map vector, whereas displacements on level of time integration are
     * handled as row map vector*/
    Epetra_Vector  discol(*discret_.DofColMap(),true);
    
    LINALG::Export(disrow,discol);
     
    for (int i = 0; i <discret_.NumMyColNodes(); ++i)
    {
      //get pointer at a node
      const DRT::Node* node = discret_.lColNode(i);
      
      //get GIDs of this node's degrees of freedom
      std::vector<int> dofnode = discret_.Dof(node);

      LINALG::Matrix<3,1> currpos;
      LINALG::Matrix<3,1> currrot;   
      
      currpos(0) = node->X()[0] + discol[discret_.DofColMap()->LID(dofnode[0])];
      currpos(1) = node->X()[1] + discol[discret_.DofColMap()->LID(dofnode[1])];;
      currpos(2) = node->X()[2] + discol[discret_.DofColMap()->LID(dofnode[2])];
      currrot(0) =                discol[discret_.DofColMap()->LID(dofnode[3])];
      currrot(1) =                discol[discret_.DofColMap()->LID(dofnode[4])];
      currrot(2) =                discol[discret_.DofColMap()->LID(dofnode[5])];
      currentpositions[node->LID()] = currpos;
      currentrotations[node->LID()] = currrot;
    }



#ifdef MEASURETIME
    const double t_admin = ds_cputime();
#endif // #ifdef MEASURETIME 
    
    //setting the crosslinkers for neighbours in crosslinkerneighbours_ after probability check
    SetCrosslinkers(dt,noderowmap,nodecolmap,currentpositions,currentrotations);
        
    //deleting the crosslinkers in crosslinkerpartner_ after probability check
    DelCrosslinkers(dt,noderowmap,nodecolmap);
     
    /*settling administrative stuff in order to make the discretization ready for the next time step: the following
     * commmand generates or deletes ghost elements if necessary and calls FillCompete() method of discretization; 
     * this is enough as long as only elements, but no nodes are added in a time step; finally Crs matrices stiff_ has
     * to be deleted completely and made ready for new assembly since their graph was changed*/     
    DRT::UTILS::RedistributeWithNewNodalDistribution(discret_,noderowmap,nodecolmap);      
    discret_.FillComplete(true,false,false);
    stiff_->Reset();


    
#ifdef MEASURETIME
    cout << "\n***\nadministration time: " << ds_cputime() - t_admin<< " seconds\n***\n";
#endif // #ifdef MEASURETIME
      

  }//if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
  
#ifdef MEASURETIME
    const double Delta_t = ds_cputime()-t_start;
    cout << "\n***\ntotal time: " << Delta_t<< " seconds\n***\n";
#endif // #ifdef MEASURETIME
  
  #endif
  
  return;
} // StatMechManager::StatMechUpdate()

/*----------------------------------------------------------------------*
 | Searches and saves in variable  crosslinkerneighbours_ neighbours for|
 | each row nod                                              cyron 07/09|
 *----------------------------------------------------------------------*/
void StatMechManager::SearchNeighbours(const std::map<int,LINALG::Matrix<3,1> > currentpositions)
{   
  //maximal distance bridged by a crosslinker
  double rlink = statmechparams_.get<double>("R_LINK",0.0);
  
  //each processor looks for each of its row nodes for neighbours
  for(int i = 0; i < discret_.NumMyRowNodes(); i++)
  { 
    //get GID and column map LID of row node i
    int iGID = (discret_.lRowNode(i))->Id();
    int icolID = (discret_.gNode(iGID))->LID();
    
    /*getting iterator to current position of local node colID from map currentpositions;
     * note that currentpositions is organized with respect to column map LIDs*/
    const map< int,LINALG::Matrix<3,1> >::const_iterator posi = currentpositions.find(icolID);
    
    /* for each row node all comlumn nodes within a range rlink are searched; this search may
     * be performed either as a brute force search (A) or by means of a search tree (B); the
     * column map LIDs of the nodes within rlink are stored in the vector neighboursLID*/
    std::vector<int> neighboursLID;
   
    //(A): brute force search
    for(std::map<int,LINALG::Matrix<3,1> >::const_iterator posj = currentpositions.begin(); posj != currentpositions.end(); posj++)
    {
      //difference vector between row node i and some column node j
      LINALG::Matrix<3,1> difference;                
      for(int k = 0; k<3; k++)
        difference(k) = (posi->second)(k) - (posj->second)(k);
      
      if(difference.Norm2() < rlink)
        neighboursLID.push_back(posj->first);
    }
    
    //(B): seraching neighbours with search tree
    //neighboursLID = octTree_->searchPointsInRadius(currentpositions,(currentpositions.find(iLID))->second,rlink);
    
    /*after having searched all nodes withing distance rlink around row node i we cancel out those
     * neighbours, which do not comply with certain requirements; here we establish the following
     * requirements: first if filament numbering is activated (i.e. not all filament numbers are set to -1) the
     * neighbour node should belong to a filament different from the one the searching node belongs to;
     * second the GID of the searched node should be greater than the GID of the searching node, as we wish
     * that crosslinkers are established from nodes with smaller GIDs to nodes with greater GIDs only; note that
     * using erase you have to be careful to keep your iterator valid despite conditional deleting of elements
     * during iteration. The following algorithms represents a very efficient and simple way to deal with this
     * problem in a correct manner*/   
    vector<int>::iterator iter = neighboursLID.begin();
    while( iter != neighboursLID.end() )
    {
      if( (  (*filamentnumber_)[*iter] == (*filamentnumber_)[icolID]  && (*filamentnumber_)[icolID] >= 0 )
         ||  (discret_.lColNode(*iter))->Id() < iGID )
        iter = neighboursLID.erase(iter);
      else
        ++iter;
    }
    
    //finally the list of column map LIDs in neighboursLID is assigned to the entry of the i-th row node in crosslinkerneighbours_
    crosslinkerneighbours_[i] = neighboursLID;   
  }

}//void SearchNeighbours(const int rlink, const std::map<int,LINALG::Matrix<3,1> > currentpositions)

/*----------------------------------------------------------------------*
 | set crosslinkers between neighbouring nodes after probability check  |
 | (private)                                                 cyron 02/09|
 *----------------------------------------------------------------------*/
void StatMechManager::SetCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,const std::map<int,LINALG::Matrix<3,1> >& currentpositions,const std::map<int,LINALG::Matrix<3,1> >& currentrotations)
{  
  //probability with which a crosslinker is established between neighbouring nodes
  double plink = 1.0 - exp( -dt*statmechparams_.get<double>("K_ON",0.0)*statmechparams_.get<double>("C_CROSSLINKER",0.0) );

  //creating a random generator object which creates uniformly distributed random numbers in [0;1]
  ranlib::UniformClosed<double> UniformGen;
    
  //we loop through all row nodes on each processor
  for(int i = 0; i < noderowmap.NumMyElements(); i++)
  { 
    //get GID of row node i
    int iGID = noderowmap.GID(i);
    
    //getting current position and rotational status of node with column map LID i from map currentpositions
    map< int,LINALG::Matrix<3,1> >::const_iterator posi         = currentpositions.find(nodecolmap.LID(iGID));
    map< int,LINALG::Matrix<3,1> >::const_iterator roti         = currentrotations.find(nodecolmap.LID(iGID));
    
    //we loop through all column map nodes, which are neighbours of row node i and test whether a crosslinked should be established
    for(int j = 0; j < (int)crosslinkerneighbours_[i].size(); j++)
    {
      //getting current position and rotational status of node with column map LID j from map currentpositions
      map< int,LINALG::Matrix<3,1> >::const_iterator posj = currentpositions.find((crosslinkerneighbours_[i])[j]);
      map< int,LINALG::Matrix<3,1> >::const_iterator rotj = currentrotations.find((crosslinkerneighbours_[i])[j]);
      
      //probability check whether crosslink should be established:
      if(UniformGen.random() < plink)
      {
        //check whether a crosslinker has already been established to the same node
        bool notyet = true;
        for(int k = 0; k < (int)crosslinkerpartner_[i].size(); k++)
        {
          if((crosslinkerpartner_[i])[k] == (crosslinkerneighbours_[i])[j] )
            notyet = false;
        }
        
        //only if no crosslinker has already been established to the same node, yet, the crosslinker is established
        if(notyet)
        {
        
          //add new crosslink to list of already established crosslinks for row node i
          crosslinkerpartner_[i].push_back( (crosslinkerneighbours_[i])[j] );
          
          //save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
          std::vector<double> rotrefe(6);
          std::vector<double> xrefe(6);
       
          for(int k = 0; k<3; k++)
          {
            //set nodal positions
            xrefe[k  ] = (posi->second)(k);
            xrefe[k+3] = (posj->second)(k);
            
            //set nodal rotations (not true ones, only those given in the displacment vector)
            rotrefe[k  ] = (roti->second)(k);
            rotrefe[k+3] = (rotj->second)(k);
          }
          
          /*there is the problem of how to assign to a crosslinker element a GID which is certainly not used by any
           * other element; we know that for a network at the beginning (without crosslinkers) each node has a
           * connectivity of 1 or 2 so that the number of all elements used to discretized the filaments is smaller
           * than basisnodes_. Thus we choose crosslinker GIDs >= basisnodes_. To make sure that two crosslinkers
           * never have the same GID we add to basisnodes_ the value basisnodes_*GID1, where GID1 is the GID of the
           * first nodes of the crosslinker element. Then we add GID2, which is the GID of the second node of the
           * crosslinker element. Hence basisnodes_ + GID1*basisnodes_ + GID2 always give a GID which cannot be 
           * used by any other element*/     
         
   #ifdef D_BEAM3         
          RCP<DRT::ELEMENTS::Beam3> newcrosslinker = rcp(new DRT::ELEMENTS::Beam3((noderowmap.GID(i) + 1)*basisnodes_ +  nodecolmap.GID((crosslinkerneighbours_[i])[j]), discret_.Comm().MyPID()) );
          
          //setting up crosslinker element parameters
          newcrosslinker ->crosssec_ = 2.375829e-05;
          newcrosslinker ->crosssecshear_ = 1.1*2.375829e-05;
          newcrosslinker ->Iyy_ = 4.49180e-11;
          newcrosslinker ->Izz_ = 4.49180e-11;
          newcrosslinker ->Irr_ = 8.9836e-11; 

  
          //nodes are assigned to the new crosslinker element by first assigning global node Ids and then assigning nodal pointers
          int globalnodeids[2] = {noderowmap.GID(i),nodecolmap.GID( (crosslinkerneighbours_[i])[j] )}; 
                 
          newcrosslinker->SetNodeIds(2,globalnodeids);
          DRT::Node *nodes[] = {discret_.gNode( globalnodeids[0] ) , discret_.gNode( globalnodeids[1] )};
          newcrosslinker->BuildNodalPointers(&nodes[0]);
                      
          //correct reference configuration data is computed for the new crosslinker element;
          //function SetUpReferenceGeometry is template and here only linear beam elements can be applied as crosslinkers
          newcrosslinker->SetUpReferenceGeometry<2>(xrefe,rotrefe); 
    
          //set material for new element
          newcrosslinker->SetMaterial(2);
          
          //add new element to discretization
          discret_.AddElement(newcrosslinker);  
          
   #endif 
   
         
    /*
   #ifdef D_TRUSS3         
           RCP<DRT::ELEMENTS::Truss3> newcrosslinker = rcp(new DRT::ELEMENTS::Truss3((noderowmap.GID(i) + 1)*basisnodes_ +  nodecolmap.GID((crosslinkerneighbours_[i])[j]), discret_.Comm().MyPID()) );
           
           //setting up crosslinker element parameters
           newcrosslinker ->crosssec_ = 1.9e-08;
           newcrosslinker ->kintype_ = DRT::ELEMENTS::Truss3::tr3_engstrain;       
          
           //nodes are assigned to the new crosslinker element by first assigning global node Ids and then assigning nodal pointers
           int globalnodeids[2] = {noderowmap.GID(i),nodecolmap.GID( (crosslinkerneighbours_[i])[j] )}; 
                  
           newcrosslinker->SetNodeIds(2,globalnodeids);
           DRT::Node *nodes[] = {discret_.gNode( globalnodeids[0] ) , discret_.gNode( globalnodeids[1] )};
           newcrosslinker->BuildNodalPointers(&nodes[0]);
                       
           //correct reference configuration data is computed for the new crosslinker element;
           //function SetUpReferenceGeometry is template and here only linear beam elements can be applied as crosslinkers
           LINALG::Matrix<6,1> xrefematrix;
           for (int node=0; node<2; node++) 
              for(int dof= 0; dof < 3; dof++)// element node has three coordinates x1, x2 and x3
                xrefematrix(node*3 + dof,0) = xrefe[node*3 + dof];
           
           newcrosslinker->SetUpReferenceGeometry(xrefematrix); 
           
           //set material for new element
           newcrosslinker->SetMaterial(1);
           
           //add new element to discretization
           discret_.AddElement(newcrosslinker);  
     
     #endif 
      */
          
         }//if(notyet)
      }//if(UniformGen.random() < plink)     
    }//for(int j = 0; j < (crosslinkerneighbours_[i]).size(); j++)
  }//for(int i = 0; i < discret_.NumMyRowNodes(); i++)
}//void SetCrosslinkers(const Epetra_Vector& setcrosslinkercol)

/*----------------------------------------------------------------------*
 | delete crosslinkers listed in crosslinkerpartner_ after random check |               
 | (private)                                                 cyron 02/09|
 *----------------------------------------------------------------------*/
void StatMechManager::DelCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap)
{ 
  
  //probability with which a crosslink breaks up in the current time step 
  double punlink = 1.0 - exp( -dt*statmechparams_.get<double>("K_OFF",0.0) );
  
  //creating a random generator object which creates uniformly distributed random numbers in [0;1]
  ranlib::UniformClosed<double> UniformGen;
  
  /*this method removes crosslinkers listed in crosslinkerpartner_ after probability check; note that 
   * removing crosslinkers should be done after setting the crosslinkers since the latter one uses 
   * operations only possible on a discretization on which already FillComplete() has been called;
   * note: here we loop through all row map nodes not by means of the discretization, but by means
   * of a node row map stored when the discretization was still in the fill complete state; thereby
   * we make sure that this method is unaffected by any changes in the discretization and its maps
   * in the turn of adding and deleting crosslinkers*/ 
  for(int i = 0; i < noderowmap.NumMyElements(); i++)
  {      
    vector<int>::iterator iter = crosslinkerpartner_[i].begin();
    while( iter != crosslinkerpartner_[i].end() )
    {
      if(UniformGen.random() < punlink)
      {       
        //we compute the GID of the crosslinker to be deleted 
        int crosslinkerid = (noderowmap.GID(i) + 1)*basisnodes_ +  nodecolmap.GID(*iter);
                
        //delete crosslinker from discretization
        if( !discret_.DeleteElement(crosslinkerid) )
          dserror("Deleting crosslinker element failed");
               
        //delete crosslinker from list of established crosslinkers
        iter = crosslinkerpartner_[i].erase(iter);
      }
      else
        ++iter;
    }
   }//for(int i = 0; i < discret_.NumMyRowNodes(); i++)
  
}//void DelCrosslinkers(const Epetra_Vector& delcrosslinkercol)


/*----------------------------------------------------------------------*
 | (public) writing restart information for manager objects   cyron 12/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechWriteRestart(IO::DiscretizationWriter& output)
{   
  output.WriteInt("istart",istart_);
  output.WriteDouble("starttimeoutput",starttimeoutput_);
  output.WriteDouble("endtoendref",endtoendref_);
  output.WriteInt("basisnodes",basisnodes_);
  output.WriteInt("outputfilenumber",outputfilenumber_);
  
  /*note: the variables crosslinkerpartner_ and crosslinkerneighbours_ are not saved here; this means
   * that for crosslinked networks restarts are not possible in general; as crosslinkerneighbours_ is
   * refreshed in the constructur of the statmech_manager class saving this variable is not that important.
   * However, not saving crosslinkerpartner_ makes a reasonable restart impossible*/
 
  return;
} // StatMechManager::StatMechOutput()

/*----------------------------------------------------------------------*
 |read restart information for statistical mechanics (public)cyron 12/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechReadRestart(IO::DiscretizationReader& reader)
{
  // read restart information for statistical mechanics  
  istart_ = reader.ReadInt("istart");
  starttimeoutput_ = reader.ReadDouble("starttimeoutput");
  endtoendref_ = reader.ReadDouble("endtoendref");
  basisnodes_ = reader.ReadInt("basisnodes");
  outputfilenumber_ = reader.ReadInt("outputfilenumber");

  return;
}// StatMechManager::StatMechReadRestart()


#endif  // #ifdef CCADISCRET
