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


#ifdef D_BEAM3
#include "../drt_beam3/beam3.H"
#endif  // #ifdef D_BEAM3
#ifdef D_TRUSS3
#include "../drt_truss3/truss3.H"
#endif  // #ifdef D_TRUSS3

#include <iostream>
#include <iomanip>


/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 09/08|
 *----------------------------------------------------------------------*/
StatMechManager::StatMechManager(ParameterList& params, DRT::Discretization& discret, RCP<LINALG::SparseMatrix>& stiff, RCP<LINALG::SparseMatrix>& damp):
  statmechparams_( DRT::Problem::Instance()->StatisticalMechanicsParams() ),
  maxtime_(params.get<double>("max time",0.0)),
  starttimeoutput_(-1.0),
  endtoendref_(0.0),
  istart_(0),
  rlink_(params.get<double>("R_LINK",0.0)),
  basiselements_(0),
  outputfilenumber_(-1),
  discret_(discret),
  stiff_(stiff),
  damp_(damp)
{ 
  /*there are two ways for Brownian motion: either by statistical forces or by statistical displacement increments: the first one
   * is activated by setting the following parameter to zero, the latter one by setting it uneuqual zero; note: currently only the 
   * first way is working properly; the second one entails to large predictor residuals and thus crashes the program; however, in 
   * principle it could be a nice option, since it allows simple division of a step into several displacement controlled load steps
   * and hence e.g. an arc-length solver in order to prevent from instabilities*/
  statmechparams_.set("FORCE_OR_DISP",0);
  
  
  /*setting and deleting dynamic crosslinkers needs a special code structure for parallel search trees
   * here we provide a fully overlapping column map which is required by the search tree to look for each
   * node for all the neighbouring nodes withing a certain distance; note: we pass this fully overlapping
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

    /*rebuild discretization basd on the new column map so that each processor creates new ghost elements
     * if necessary; after the following line we have a discretization, where each processor has a fully
     * overlapping column map regardlesse of how overlapping was managed when starting BACI; having ensured
     * this allows convenient and correct (albeit not necessarily efficient) use of search algorithms and
     * crosslinkers in parallel computing*/     
    discret_.FillComplete(true,false,false);
  
  
  //if dynamic crosslinkers are used additional variables are initialized
  if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
  {  
    /* crosslinkerpartner_ and filamentnumber_ is generated based on a column map vector as each node has to 
     * know about each other node its filament number in order to decide weather a crosslink may be established 
     * or not; both vectors are initalized with -1 indicating that no crosslinkers have been set so far nor
     * filament numbers*/
    crosslinkerpartner_ = rcp( new Epetra_Vector(*discret_.NodeColMap()) ); 
    filamentnumber_ = rcp( new Epetra_Vector(*discret_.NodeColMap()) );
    crosslinkerpartner_->PutScalar(-1);
    filamentnumber_->PutScalar(-1);
    

    
    /*force sensors can be applied at any degree of freedom of the discretization the list of force sensors should
     * be based on a column map vector so that each processor has not only the information about each node's
     * displacement, but also about whether this has a force sensor; as a consequence each processor can write the
     * complete information gathered by all force sensors into a file of its own without any additional communication
     * with any other processor; initialization with -1 indicates that so far no forcesensors have been set*/
    forcesensor_ = rcp( new Epetra_Vector(*discret_.DofColMap()) );
    forcesensor_->PutScalar(-1);
    
    //basiselements_ is the number of elements existing in the discretization from the very beginning on
    basiselements_ = discret_.NumGlobalElements();

  
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
      int filamentnumber = filamentnumberconditions[i]->Getint("Filament Number") ;
      
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
      int nodedofnumber = forcesensorconditions[i]->Getint("DOF Number") ;
      
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
        int dofnumber = discret_.Dof( discret_.gNode(nodenumber), nodedofnumber );
        
        
        /*if the node does not belong to current processor Id is set by LID() to -1 and the node is ignored; otherwise the degrees of
         * freedom affected by this condition are marked in the vector *forcesensor_ by an one entry*/
        if(nodenumber > -1)
          (*forcesensor_)[dofnumber] = 1;     
      }
     } 
  
      //we construct a search tree with maximal depth of 8 (different numbers might be tried out, too)
      octTree_ = rcp(new GEO::SearchTree(8));
     
      
      /*after having generated a search tree and a discretization with fully overlapping column map we initialize a search tree
       * for accelerated search for each nodes neighbouring nodes; note: the tree is based on a bounding box
       * with respect to the reference positions (which are the current positions at the beginning of the simulation;
       * in case of large overall deformations of the fiber network such an initialization would have to be carried 
       * out in each time step with respect to the current node positions*/

      std::map<int,LINALG::Matrix<3,1> > currentpositions;
      
      currentpositions.clear();

      for (int lid = 0; lid <discret_.NumMyColNodes(); ++lid)
      {
        const DRT::Node* node = discret_.lColNode(lid);
        LINALG::Matrix<3,1> currpos;
        currpos(0) = node->X()[0];
        currpos(1) = node->X()[1];
        currpos(2) = node->X()[2] ;
        currentpositions[node->Id()] = currpos;
      }
      
      //find bounding box for search tree
      const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(discret_, currentpositions);
   
      //initialize search tree
      octTree_->initializeTree(rootBox, discret_, GEO::TreeType(GEO::OCTTREE));
    
  }//if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))

  return;
} // StatMechManager::StatMechManager


/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechOutput(const double& time,const int& num_dof,const int& istep, const double& dt, const Epetra_Vector& dis, const Epetra_Vector& fint)
{
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
    case INPAR::STATMECH::statout_endtoendergodicity:
    {
      FILE* fp = NULL; //file pointer for statistical output file
      double endtoend = 0; //end to end length at a certain time step in single filament dynamics
       
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
    //measurement of viscoelastic properties should be carried out by means of force sensors
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
      
      filecontent << endl << endl << endl << scientific << setprecision(10) << "measurement data in timestep " << istep << ", time = " << time << endl << endl;
      
      for(int i = 0; i < forcesensor_->MyLength(); i++)
        if( (*forcesensor_)[i] == 1)
          filecontent << "displacement = "<< dis[i] << "  internal force = "<< fint[i] << endl;
         
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
       * e.g. number one is written as 00001, number fourteen as 00014 and so on;*/
      
      //note: this kind of output is possilbe for serial computing only (otherwise the following method would have to be adapted to parallel use*/   
      if(discret_.Comm().NumProc() > 1)
        dserror("No Gmsh output for parallel computation possible so far");
      
      // first index = time step index
      std::ostringstream filename;

      //creating complete file name dependent on step number with 5 digits and leading zeros
      if (istep<100000)
        filename << "./GmshOutput/network"<< std::setw(5) << setfill('0') << istep <<".pos";
      else 
        dserror("Gmsh output implemented for a maximum of 99999 steps");
          
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
      if (element->Id() < basiselements_)
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
void StatMechManager::StatMechInitOutput()
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
       numberingfilename.str("NumberOfRealizations");
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
    case INPAR::STATMECH::statout_endtoendergodicity:
    {
      //defining name of output file
      std::ostringstream outputfilename;
      outputfilename.str("");
      outputfilename << "EndToEndErgo.dat";
      
      FILE* fp = NULL; //file pointer for statistical output file
      
      //making sure that there exists an now empty file named by outputfilename_
      fp = fopen(outputfilename.str().c_str(), "w");
      
      fclose(fp);
      
    }
    break;
    //measurement of viscoelastic properties should be carried out by means of force sensors
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
        
      filecontent << "Output for measurement of viscoelastic properties written by processor "<< discret_.Comm().MyPID() << endl;
         
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
     * current positions of all the nodes in the map currentpositions; additionally we store the roational displacments
     * analogously in a map currentrotations for later use in setting up reference geometry of crosslinkers*/
    std::map<int,LINALG::Matrix<3,1> > currentpositions;   
    std::map<int,LINALG::Matrix<3,1> > currentrotations; 
    currentpositions.clear();
    currentrotations.clear();
    
    //note: access by ExtractMyValues requires column map vector
    Epetra_Vector  discol(*discret_.DofColMap(),true);
    
    LINALG::Export(disrow,discol);
     
    for (int lid = 0; lid <discret_.NumMyColNodes(); ++lid)
    {
      //get pointer at a node
      const DRT::Node* node = discret_.lColNode(lid);
      
      //get GIDs of this node's degrees of freedom
      vector<int> dofnode = discret_.Dof(node);

      LINALG::Matrix<3,1> currpos;
      LINALG::Matrix<3,1> currrot;   
      
      currpos(0) = node->X()[0] + discol[discret_.DofColMap()->LID(dofnode[0])];
      currpos(1) = node->X()[1] + discol[discret_.DofColMap()->LID(dofnode[1])];;
      currpos(2) = node->X()[2] + discol[discret_.DofColMap()->LID(dofnode[2])];
      currrot(0) =                discol[discret_.DofColMap()->LID(dofnode[3])];
      currrot(1) =                discol[discret_.DofColMap()->LID(dofnode[4])];
      currrot(2) =                discol[discret_.DofColMap()->LID(dofnode[5])];
      currentpositions[node->Id()] = currpos;
      currentrotations[node->Id()] = currrot;
    }

    
    //define map between the column map LID of a node and the column map LIDs of all its neighbours
    std::map<int,std::set<int> > neighbours;   
       
    //probability with which a crosslinker is established between neighbouring nodes
    double plink = 1.0 - exp( -dt*statmechparams_.get<double>("K_ON",0.0)*statmechparams_.get<double>("C_CROSSLINKER",0.0) );
      
    //probability with which a crosslink breaks up in the current time step 
    double punlink = 1.0 - exp( -dt*statmechparams_.get<double>("K_OFF",0.0) );
        
    //maximal distance bridged by a crosslinker
    double rlink = statmechparams_.get<double>("R_LINK",0.0);
      
    /*create column map vector which stores for each node of the column map the nearest neighbour node's column map LID;
     * the vector is initialized with -1 in each entry*/
    Epetra_Vector nearestneighbour(nodecolmap);
    nearestneighbour.PutScalar(-1);
 

          
    /*create random numbers in order to decide whether a crosslinker should be established; each
     * processor creates random numbers for its own nodes only; all these numbers are exported to
     * a column map vector so that each processor has access to the same information (additive export
     * allowed due to zero initialization before!); finally for the column map vector each element is
     * considered: a node i can get a crosslinker in the following step if it has not yet one and if
     * the random number calculated for it passes the probability test (i.e. is smaller than -1 + 2*plink).
     * In this case the element i is turned into 1 whereas otherwise into 0. As a consequence finally each
     * element i indicates whether node i can get a crosslinker (entry 1) or not (entry 0); note that this
     * postprocessing of the elements has to be done for the column map vector since also the information
     * of crosslinkerpartner_ is related to a column map*/
    Epetra_Vector  setcrosslinkerrow(noderowmap);
    Epetra_Vector  setcrosslinkercol(nodecolmap);

    setcrosslinkercol.PutScalar(0);
    setcrosslinkerrow.Random();
    
    LINALG::Export(setcrosslinkerrow,setcrosslinkercol);
       
    for(int i = 0; i < setcrosslinkercol.MyLength(); i++)
    {
      if( setcrosslinkercol[i] <  -1.0 + 2*plink && (*crosslinkerpartner_)[i] == -1.0)
        setcrosslinkercol[i] = 1;
      else
        setcrosslinkercol[i] = 0;
    }  

    
    /*create random numbers in order to decide whether a crosslinker should be deleted; each
     * processor creates random numbers for its own nodes only; all these numbers are exported to
     * a column map vector so that each processor has access to the same information; finally for the 
     * column map vector each element is considered: a node i can loose its crosslinker in the following 
     * step if it has already one and if the random number calculated for it passes the probability test 
     * (i.e. is smaller than -1 + 2*punlink). Furthermore the node has to have the lower global Id of the
     * two nodes related to the crosslinker; this makes sure that the off-rate is a quantitiy with respect
     * to crosslinkers and not with respect to their nodes; the node with the lower GID is kind of the owner
     * of the crosslinker responsible for its element GID, for creating it and for deleting it. 
     * If the above three conditions are satisfied the element i is turned into 1 whereas otherwise 
     * into 0. As a consequence finally each element i indicates whether node i can get a crosslinker 
     * (entry 1) or not (entry 0); note that this postprocessing of the elements has to be done for the column
     *  map vector since also the information of crosslinkerpartner_ is related to a column map*/
    Epetra_Vector  delcrosslinkerrow(noderowmap);
    Epetra_Vector  delcrosslinkercol(nodecolmap);
    delcrosslinkercol.PutScalar(0);
    delcrosslinkerrow.Random();
    
    LINALG::Export(delcrosslinkerrow,delcrosslinkercol);
    
    for(int i = 0; i < delcrosslinkercol.MyLength(); i++)
    {
      if( delcrosslinkercol[i] <  -1.0 + 2*punlink && (*crosslinkerpartner_)[i] != -1.0 && nodecolmap.GID(i) < nodecolmap.GID( (*crosslinkerpartner_)[i] )  )
        delcrosslinkercol[i] = 1;
      else
        delcrosslinkercol[i] = 0;
    }  
    
    
    //after having determined for which nodes crosslinkers may be set or deleted we search for each node all its neighbours
    
    //the clever way to do so is using a search tree
    //neighbours = octree_->searchNodesInRadius(discret_,currentpositions,xrefe,rlink)
    
    //the brute force way is to consider all the other nodes as potential neighbours
    for(int i = 0; i < discret_.NumMyColNodes(); i++)
    { 
      //for each node the set of all potential neighbours consists of all column map LIDs
      std::set<int> itset;
      for(int id = 0; id < discret_.NumMyColNodes(); id++)
        itset.insert(id);
      
      neighbours.insert(make_pair(i,itset));
    }
    
    //from the above preselection of neighbours we select the nearest neighbour which is not on the same filament
    for(int i = 0; i < discret_.NumMyColNodes(); i++) //iterating through all the elements of the map neighbours
    {  
      /*declaring constant iterator for access to elements of map "neighbours" related with number i; this iterator
       * hast two components which can be addressed by "->first" and "->second", where the first component is the
       * number i itself and the second one the set at which number i points*/
      map< int,std::set<int> >::iterator mapit = neighbours.find(i);
      
      //getting set at which number i points in map neighbours;
      std::set<int> ineighbourset = mapit->second;
      
      //distance of the so far nearest neighbour (initialized with rlink) 
      double rneighbour = rlink;
      
      //getting current position of node i from map currentpositions
      map< int,LINALG::Matrix<3,1> >::iterator posi = currentpositions.find(i);
                      
      for(int j = 0; j < (int)ineighbourset.size(); j++) //iterating through set ineighbourset
      {
        /*declaring constant iterator for access to elements of set "neighbours" related with number j; this iterator 
         * can be handled as a pointer at the number j points at in the set*/
        std::set<int>::const_iterator setit = ineighbourset.find(j);
        
        //getting current position of node *setit from map currentpositions
        map< int,LINALG::Matrix<3,1> >::iterator posj = currentpositions.find(*setit);
        
        
        /*crosslinkers are established only if either two nodes belong to the same filament or if the filament
         * numbering is deactivated (-> all filamentnumbers set to -1); furthermore the node with the higher
         * global Id decides whether a crosslinker is established between two nodes so that we search only the
         * nodes with lower global Ids whether they are close to a certain node */
         if( (*filamentnumber_)[i] != (*filamentnumber_)[*setit] || (*filamentnumber_)[i] == -1)
         {
           //difference vector between node i and j
           LINALG::Matrix<3,1> difference;                
           for(int k = 0; k<3; k++)
             difference(k) = (posi->second)(k) - (posj->second)(k);

           
           if(difference.Norm2() < rneighbour)
           {
             nearestneighbour[i] = *setit;
             rneighbour = difference.Norm2();
           }           
          }
        }
      } //for(int i = 0; i < discret_.NumMyRowNodes(); i++)
    
    //setting the crosslinkers of nodes marked by a one entry in vector delcrosslinkercol
    SetCrosslinkers(setcrosslinkercol,nodecolmap,nearestneighbour,currentpositions,currentrotations);

     
    //deleting the crosslinkers of all nodes marked by a one entry in vector delcrosslinkercol
    DelCrosslinkers(delcrosslinkercol,nodecolmap);
        
  
    /*settling administrative stuff in order to make the discretization ready for the next time step: the following
     * commmand generates or deletes ghost elements if necessary and calls FillCompete() method of discretization; 
     * this is enough as long as only elements, but no nodes are added in a time step; finally Crs matrices stiff_ and
     * damp_ have to be deleted completely and made ready for new assembly since their graph was changed*/        
    DRT::UTILS::RedistributeWithNewNodalDistribution(discret_,noderowmap,nodecolmap);
    discret_.FillComplete(true,false,false);
    stiff_->Reset();
    damp_->Reset();
      

  }//if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
  
  #endif
  
  return;
} // StatMechManager::StatMechUpdate()

/*----------------------------------------------------------------------*
 | set crosslinkers according to setcrosslinkercol (if entry for a node |
 | equals 1 a crosslinker may be set for this node) (private)           |
 |                                                           cyron 02/09|
 *----------------------------------------------------------------------*/
void StatMechManager::SetCrosslinkers(const Epetra_Vector& setcrosslinkercol,const Epetra_Map& nodecolmap,const Epetra_Vector& nearestneighbour,const std::map<int,LINALG::Matrix<3,1> >& currentpositions,const std::map<int,LINALG::Matrix<3,1> >& currentrotations)
{
  /*this method sets crosslinkers of nodes marked by a one entry in vector delcrosslinkercol; the method
   * sets between two nodes with GIDs n and m for a crosslinker the global element Id basiselements_ + min(n,m); note:
   * in this method not only crosslinkers are added to certain nodes by their owner processors, but also all the processor
   * refresh their variable *crosslinkerpartner_ since they know even of crosslinkers whose owner is another processor*/ 
  
  /*in order not to have to set all the parameters for each single crosslinker element a dummy crosslinker element 
   * is set up before setting crosslinkers*/ 
  #ifdef D_BEAM3
      RCP<DRT::ELEMENTS::Beam3> crosslinkerdummy;   
      crosslinkerdummy = rcp(new DRT::ELEMENTS::Beam3(-1,discret_.Comm().MyPID()) );     
      crosslinkerdummy->crosssecshear_ = 19e-10*1.1;
      crosslinkerdummy->Iyy_ = 28.74e-11;
      crosslinkerdummy->Izz_ = 28.74e-11;
      crosslinkerdummy->Irr_ = 57.48e-11;        
  #endif      
   /*
  #ifdef D_TRUSS3
      RCP<DRT::ELEMENTS::Truss3> crosslinkerdummy;   
      crosslinkerdummy = rcp(new DRT::ELEMENTS::Truss3(-1,discret_.Comm().MyPID()) );       
      crosslinkerdummy->crosssec_ = 19e-6;
  #endif
  */

   
  for(int i = 0; i < (int)currentpositions.size(); i++)
  {
    
    //getting current position of node i and its nearesneighbour from map currentpositions
    map< int,LINALG::Matrix<3,1> >::const_iterator posi         = currentpositions.find(i);
    map< int,LINALG::Matrix<3,1> >::const_iterator roti         = currentrotations.find(i);
    map< int,LINALG::Matrix<3,1> >::const_iterator posneighbour = currentpositions.find((int)nearestneighbour[i]);
    map< int,LINALG::Matrix<3,1> >::const_iterator rotneighbour = currentrotations.find((int)nearestneighbour[i]);
   
    
    /*if a neighbour within the prescribed maximal distance has been found and this neighbour node has not yet
    * a crosslinker itself a crosslinker has to be established in case that also the probability check was passed and
    * that the neighbour has a higher GID; note that the latter condition makes sure that for each crosslinker always
    * the node with the lower GID is kind of the owner: it has to establish it, to delete it, and it is respondible 
    * for the element GID of the crosslinker */
    if(setcrosslinkercol[i] && nearestneighbour[i] > -1 && (*crosslinkerpartner_)[ (int)nearestneighbour[i] ] == -1.0 && nodecolmap.GID(i) < nodecolmap.GID(nearestneighbour[i]))
    {
      //the crosslinker to be established is registered in the variable crosslinkerpartner_
      (*crosslinkerpartner_)[i] = (int)nearestneighbour[i];
      (*crosslinkerpartner_)[ (int)nearestneighbour[i] ] = i;
       
      
      //only the owner processor of node i establishes the crosslinker actually
      if( (discret_.gNode(nodecolmap.GID(i)))->Owner() == discret_.Comm().MyPID() )
      {    
       //fixed size variable for storing positions and rotational displacmenets of the two nodes to be crosslinked
       LINALG::Matrix<6,1> xrefe;
       LINALG::Matrix<6,1> rotrefe;
      
       //current position of nodes with LIDs i and nearestneighbour[i]
       for(int k = 0; k<3; k++)
       {
         //set nodal positions
         xrefe(k  ) = (posi->second)(k);
         xrefe(k+3) = (posneighbour->second)(k);
         
         //set nodal rotations (not true ones, only those given in the displacment vector)
         rotrefe(k  ) = (roti->second)(k);
         rotrefe(k+3) = (rotneighbour->second)(k);
       }

       
       /*a new crosslinker element is generated according to a crosslinker dummy defined during construction 
        * of the statmech_manager; note that the dummy has already the proper owner number*/                   
#ifdef D_BEAM3
       RCP<DRT::ELEMENTS::Beam3> newcrosslinker = rcp(new DRT::ELEMENTS::Beam3(*crosslinkerdummy) );
       
      
       //assigning correct global Id to new crosslinker element: since each node can have one crosslinker element
       //only at the same time a unique global Id can be found by taking the number of elemnts in the discretization
       //before starting dealing with crosslinkers and adding to the smaller one of the two involved global nodal Ids     
       newcrosslinker->SetId( basiselements_ + min(nodecolmap.GID(i),nodecolmap.GID((int)nearestneighbour[i])) );          
                
       //nodes are assigned to the new crosslinker element by first assigning global node Ids and then assigning nodal pointers
       int globalnodeids[] = {nodecolmap.GID(i), nodecolmap.GID((int)nearestneighbour[i])};      
       newcrosslinker->SetNodeIds(2, globalnodeids);
       DRT::Node *nodes[] = {discret_.gNode( globalnodeids[0] ) , discret_.gNode( globalnodeids[1] )};
       newcrosslinker->BuildNodalPointers(&nodes[0]);
                   
       //correct reference configuration data is computed for the new crosslinker element
       newcrosslinker->SetUpReferenceGeometry(xrefe,rotrefe); 
       
       //setting drag coefficient for new crosslinker element
       newcrosslinker->zeta_ = 4*PI*newcrosslinker->lrefe_*( statmechparams_.get<double>("ETA",0.0) );
       
       //set material for new element
       newcrosslinker->SetMaterial(1);
       
       //add new element to discretization
       discret_.AddElement(newcrosslinker); 
       
#endif
       
      } //if( (discret_.lColNode(i))->Owner() == discret_.Comm().MyPID() )
    } //if(nearestneighbour[i] > -1 && (*crosslinkerpartner_)[neighbour] == -1.0)
   } //for(int i = 0; i < discret_.NumMyColNodes(); i++)
  
}//void SetCrosslinkers(const Epetra_Vector& setcrosslinkercol)

/*----------------------------------------------------------------------*
 | delete crosslinkers according to delcrosslinkercol (if entry for a   |
 | node equals 1 its crosslinker is deleted) (private)                  |
 |                                                           cyron 02/09|
 *----------------------------------------------------------------------*/
void StatMechManager::DelCrosslinkers(const Epetra_Vector& delcrosslinkercol,const Epetra_Map& nodecolmap)
{
  /*this method removes crosslinkers of nodes marked by a one entry in vector delcrosslinkercol; the method
   * assumes that a crosslinker between two nodes with GIDs n and m has the global element Id basiselements_ + min(n,m)
   * note that removing crosslinkers should be done after setting the crosslinkers since the latter one uses 
   * operations only possible on a discretization on which already FillComplete() has been called*/ 
  for(int i = 0; i < discret_.NumMyColNodes(); i++)
  {
    if(delcrosslinkercol[i])
    {        
      //we compute the global Id of the crosslinker to be deleted 
      int crosslinkerid = basiselements_ + min(nodecolmap.GID(i),nodecolmap.GID((int)(*crosslinkerpartner_)[i]));
      
      //regardless of whether the node belongs to this processor or not it has to be stored that it looses its crosslinker
      (*crosslinkerpartner_)[(int)(*crosslinkerpartner_)[i]] = -1;
      (*crosslinkerpartner_)[i] = -1;
      
      //the owner processor deletes the crosslinker           
      if( discret_.gElement(crosslinkerid)->Owner() == discret_.Comm().MyPID() )
        if( !discret_.DeleteElement(crosslinkerid) )
          dserror("Deleting crosslinker element failed");  
     }
   }//for(int i = 0; i < discret_.NumMyColNodes(); i++)
  
}//void DelCrosslinkers(const Epetra_Vector& delcrosslinkercol)

/*----------------------------------------------------------------------*
 | updates system damping matrix and either external force vector or    |
 | displacement vector according to influence of thermal bath (public)  |
 |                                                           cyron 10/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechBrownian(ParameterList& params, RCP<Epetra_Vector> dis, RCP<Epetra_Vector> fext, RCP<LINALG::SparseMatrix> damp)
{  
  // zero out damping matrix
  damp->Zero();
  
  /*declaration of a column and row map Epetra_Vector for evaluation of statistical forces or Brownian stets in space;
   *  note: zero initilization mandatory for correct computations later on*/
  RCP<Epetra_Vector>    browniancol;
  browniancol = LINALG::CreateVector(*discret_.DofColMap(),true);
  RCP<Epetra_Vector>    brownianrow;
  brownianrow = LINALG::CreateVector(*discret_.DofRowMap(),true);
  
  //defining parameter list passed down to the elements in order to evalute statistical forces down there
  ParameterList pstat;
  pstat.set("action","calc_brownian_damp");
  pstat.set("delta time",params.get<double>("delta time",0.0));
  pstat.set("KT",statmechparams_.get<double>("KT",0.0));
  pstat.set("ETA",statmechparams_.get<double>("ETA",0.0));
  pstat.set("STOCH_ORDER",statmechparams_.get<int>("STOCH_ORDER",0));
  pstat.set("FORCE_OR_DISP",statmechparams_.get<int>("FORCE_OR_DISP",0));
  
  
  /*note: the column map statistical force or movement vector is passed down via the parameter list and not as a systemvector 
   * so that assembly is not done by the evaluate method itself, but elementwise; this is in order to account for the
   * special assembly needs of randomly evaluated variables: the evaluate method of the discretization uses the LINALG
   * assembly method in which a processor assembles to the element of the global vector only if the processor is the 
   * row owner of the related DOF; this is efficient for global row map vectors, but does not assemble correctly if a 
   * global column map vector is used*/
  pstat.set("statistical vector",browniancol);
  
  
  //evaluation of statistical Brownian forces or displacements on column map vecotor
  discret_.SetState("displacement",dis); //during evaluation of statistical forces or steps in space access to current displacement possible
  discret_.Evaluate(pstat,damp,null,null,null,null);
  discret_.ClearState();
  
  
  /*exporting col map statistical force/ displacement vector to a row map vector additively, i.e. in such a way that a 
   * vector element with a certain GID in the final row vector is the sum of all elements of the column 
   * vector related to the same GID*/
  Epetra_Export dofexporter(*discret_.DofColMap(),*discret_.DofRowMap());
  brownianrow->Export(*browniancol,dofexporter,Add);
  
  //adding Brownian forces to external forces passed to this method or adding Brownian displacements to current displacements
  if(statmechparams_.get<int>("FORCE_OR_DISP",0) == 0)
    fext->Update(1.0,*brownianrow,1.0);
  else
    dis->Update(1.0,*brownianrow,1.0);
  
  //complete damping matrix
  damp->Complete();

  return;
} // StatMechManager::StatMechBrownian()

/*----------------------------------------------------------------------*
 | (public) writing restart information for manager objects   cyron 12/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechWriteRestart(IO::DiscretizationWriter& output)
{   
  output.WriteInt("istart",istart_);
  output.WriteDouble("starttimeoutput",starttimeoutput_);
  output.WriteDouble("endtoendref",endtoendref_);
  output.WriteInt("basiselements",basiselements_);
  output.WriteInt("outputfilenumber",outputfilenumber_);
  /*if no dynamic crosslinkers are calculated the variable crosslinkerpartner is never initialized
   * and hence cannot be saved or read*/
  if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
    output.WriteVector("crosslinkerpartner",crosslinkerpartner_,IO::DiscretizationWriter::nodevector);
  
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
  basiselements_ = reader.ReadInt("basiselements");
  outputfilenumber_ = reader.ReadInt("outputfilenumber");
  /*if no dynamic crosslinkers are calculated the variable crosslinkerpartner is never initialized
   * and hence cannot be saved or read*/
  if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
   reader.ReadVector(crosslinkerpartner_,"crosslinkerpartner");
 
  return;
}// StatMechManager::StatMechReadRestart()


#endif  // #ifdef CCADISCRET
