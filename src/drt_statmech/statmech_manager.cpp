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
  //if dynamic crosslinkers are used additional variables are initialized
  if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
  {
    /*crosslinkerpartner_, filamentnumber_ are generated based on the discretization
     * as it is generated from the input data file; this discretization does not yet comprise
     * any additional crosslinker elements which is o.K. since both variables are applied in 
     * order to administrate only the nodes of the original discretization*/
    crosslinkerpartner_ = rcp( new Epetra_Vector(*discret_.NodeRowMap()) );
    filamentnumber_ = rcp( new Epetra_Vector(*discret_.NodeRowMap()) );
    
    //force sensors can be applied at any degree of freedom of the discretization
    forcesensor_ = rcp( new Epetra_Vector(*discret_.DofRowMap()) );
    
    //basiselements_ is the number of elements existing in the discretization from the very beginning on
    basiselements_ = discret_.NumGlobalElements();
    
    
    for(int i = 0; i < discret_.NumMyRowNodes(); i++)
    {
      /*initializing crosslinkerpartner_ and filamentnumber_ with -1 for each node by looping
       * over all local indices assigned by the map dis.NodeRowMap() to a processor*/
      (*crosslinkerpartner_)[i] = -1;
      (*filamentnumber_)[i] = -1;
    }
    
    for(int i = 0; i < forcesensor_->MyLength(); i++)
    {
      /*initializing forcesensor_ with -1 for each degree of freedom by looping
       * over all local indices assigned by the map dis.DofRowMap() to a processor*/
      (*forcesensor_)[i] = -1;
    }
    
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
        nodenumber = discret_.NodeRowMap()->LID(nodenumber);
        
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
    
  }

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
      const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
      if (istep<100000)
        filename << "./GmshOutput/network"<< std::setw(5) << setfill('0') << istep <<".pos";
      else 
        dserror("Gmsh output implemented for a maximum of 99999 time steps");
               
      // do output to file in c-style
      FILE* fp = NULL;
      
      //start Gmsh output
      fp = fopen(filename.str().c_str(), "w");

      // write output to temporary stringstream; 
      std::stringstream gmshfilecontent;
      /*the beginning of the stream is "View \"" to indicate Gmsh that the following data is in order to create an image and
       * this command is followed by the name of that view displayed during it's shown in the video; in the following example
       * this name is for the 100th time step: Step00100; then the data to be presented within this view is written within { ... };
       * in the following example this data consists of scalar lines defined by the coordinates of their end points*/
      gmshfilecontent << "View \" Step " << istep << " \" {" << endl;
      
      
      //writing the origin point as reference for colors
      gmshfilecontent << "SP(" << scientific;
      gmshfilecontent << 0.0 << "," << 0.0 << "," << 0.0 ;
      /*note: for each node there is one color variable for gmsh and gmsh finally plots the line
       * interpolating these two colors between the nodes*/
      gmshfilecontent << ")" << "{" << scientific << 0.0  << "};" << endl;

      //looping through all elements on the processor
      for (int i=0; i<discret_.NumMyRowElements(); ++i)
      {
        //getting pointer to current element
        DRT::Element* element = discret_.lRowElement(i);
        
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
           double displacement = dis[discret_.DofRowMap()->LID( dofnode[id] )];
           coord(id,jd) =  referenceposition + displacement;
         }
        
              
        //declaring variable for color of elements
        double color;
        
        //apply different colors for elements representing filaments and those representing dynamics crosslinkers
        if (element->Id() < basiselements_)
          color = 1.2;
        else
          color = 2.6;

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
      
    }
    break;
    case INPAR::STATMECH::statout_none:
    default:
    break;
  }
 
  return;
} // StatMechManager::StatMechOutput()


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
void StatMechManager::StatMechUpdate(const double dt, const Epetra_Vector& dis)
{  
  #ifdef D_BEAM3
  
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
    double plink = 1.0 - exp( -dt*statmechparams_.get<double>("K_ON",0.0)*statmechparams_.get<double>("C_CROSSLINKER",0.0) );
    
    //probability with which a crosslink breaks up in the current time step 
    double punlink = 1.0 - exp( -dt*statmechparams_.get<double>("K_OFF",0.0) );
    
    //maximal distance bridged by a crosslinker
    double rlink = statmechparams_.get<double>("R_LINK",0.0);
    
    //searching neighbours of the current node (with global Id "neighbour"):
    std::vector<int> neighbours(0);
       
    /*in order not to have to set all the parameters for each single crosslinker element a dummy crosslinker element 
     * is set up before setting and deleting crosslinkers*/

    
#ifdef D_BEAM3
    RCP<DRT::ELEMENTS::Beam3> crosslinkerdummy;   
    crosslinkerdummy = rcp(new DRT::ELEMENTS::Beam3(-1,discret_.Comm().MyPID()) );
    
    crosslinkerdummy->crosssecshear_ = 19e-6*1.1;
    crosslinkerdummy->Iyy_ = 28.74e-12;
    crosslinkerdummy->Izz_ = 28.74e-12;
    crosslinkerdummy->Irr_ = 28.74e-9;   
    crosslinkerdummy->material_ = 1;   
#endif
    
    /*
#ifdef D_TRUSS3
    RCP<DRT::ELEMENTS::Truss3> crosslinkerdummy;   
    crosslinkerdummy = rcp(new DRT::ELEMENTS::Truss3(-1,discret_.Comm().MyPID()) );
    
    crosslinkerdummy->crosssec_ = 19e-6;
    crosslinkerdummy->material_ = 1;
#endif
*/
    
    /*the following tow rcp pointers are auxiliary variables which are needed in order provide in the very end of the
     * crosslinker administration a node row and column map; these maps have to be taken here before the first modification
     * by deleting and adding elements have been carried out with the discretization since after such modifications the maps
     * cannot be called from the discretization before calling FillComplete() again which should be done only in the very end
     * note: this way of getting node maps is all right only if no nodes are added/ deleted, but elements only*/
    RefCountPtr<Epetra_Map> noderowmap = rcp(new Epetra_Map(*(discret_.NodeRowMap())));
    RefCountPtr<Epetra_Map> nodecolmap = rcp(new Epetra_Map(*(discret_.NodeColMap())));
    

    
    
    /*_________________________________________________________________________________________________________
     * note: the following part of the code is suitable for serial use only !!!!
     * _______________________________________________________________________________________________________*/
    
     
    //searching locally with a tree all nodes forming neighbouring couples and establishing crosslinkers between them
    for(int i = 0; i < discret_.NumMyRowNodes(); i++)
    {
      //if the node is already crosslinked the crosslinker is removed with probability punlink
      if ((*crosslinkerpartner_)[i] != -1.0 && (*delcrosslinker)[i] <  -1.0 + 2*punlink)
      {
        //noting that node responsible for the crosslinker to be deleted has from now on no crosslinker
        (*crosslinkerpartner_)[i] = -1;

        /*since each node can have only one crosslinker at the same time a unique Id for the crosslinker 
         * element can be found by adding the lower one of the global Ids of the two nodes to basiselements_;
         * this is the way new crosslinkers are numbered below and hence also the way old crosslinkers can
         * be deleted*/
        int crosslinkerid = basiselements_ + discret_.NodeRowMap()->GID(i);
        
        /*delete crosslinker element; when trying to delete a crosslinker element not administrated by the calling processor
         * an error is issued*/
        if( !discret_.DeleteElement(crosslinkerid) )
        {
          dserror("Deleting crosslinker element failed"); 
        }
      }
         
      
      /*in the following new crosslinkers are established; however this is possibleonly if the current node has not 
       * yet one, i.e. (*crosslinkerpartner_)[i] == -1.0 and if the node has passed a random test by which it's figured
       * out whether it can get a new crosslinker in this time step, i.e. (*setcrosslinker)[i] <  -1.0 + 2*plink*/
     if( (*crosslinkerpartner_)[i] == -1.0 && (*setcrosslinker)[i] <  -1.0 + 2*plink  )
     {    
        //fixed size variable for storing positions of the two nodes to be crosslinked
        LINALG::Matrix<6,1> xrefe;
       
        //current position of node with LID i  
        for(int k = 0; k<3; k++)
          xrefe(k) = (discret_.lRowNode(i))->X()[k] + dis[discret_.DofRowMap()->LID( discret_.Dof(discret_.lRowNode(i),k) )];
       
       //searching nearest neighbour of the current node (with global Id "neighbour"):
        int neighbour = -1;
        double rneighbour = rlink;
        
        //we apply a naive search algorithm which migth be later replaced by means of a tree
        for(int j = 0; j < i; j++)
        {
          /*crosslinkers are established only if either two nodes belong to the same filament or if the filament
           * numbering is deactivated (-> all filamentnumbers set to -1); furthermore the node with the higher
           * global Id decides whether a crosslinker is established between two nodes so that we search only the
           * nodes with lower global Ids whether they are close to a certain node */
           if( (*filamentnumber_)[i] != (*filamentnumber_)[j] || (*filamentnumber_)[i] == -1)
           {
             //variable for convenient storage of coordinates of node j
             LINALG::Matrix<3,1> xloop;
             
              //current position of node with LID j  
             for(int k = 0; k<3; k++)
               xloop(k) = (discret_.lRowNode(j))->X()[k] + dis[discret_.DofRowMap()->LID( discret_.Dof(discret_.lRowNode(j),k) )];
     
              double rcurrent = pow(pow(xrefe(0) - xloop(0),2)+pow(xrefe(1) - xloop(1),2)+pow(xrefe(2) - xloop(2),2),0.5);
            
              if(rcurrent < rneighbour)
              {
                  neighbour = j;
                  rneighbour = rcurrent;
                  for(int k = 0; k<3; k++)
                    xrefe(k+3) = xloop(k);
              } 
            }
          }
        
        //if a neighbour within the prescribed maximal distance has been found a crosslinker is established
        if(neighbour > -1)
        {
          /*a new crosslinker element is generated according to a crosslinker dummy defined during construction 
           * of the statmech_manager; note that the dummy has already the proper owner number*/          
          
#ifdef D_BEAM3
          RCP<DRT::ELEMENTS::Beam3> newcrosslinker = rcp(new DRT::ELEMENTS::Beam3(*crosslinkerdummy) );
#endif
          
/*
#ifdef D_TRUSS3
          RCP<DRT::ELEMENTS::Truss3> newcrosslinker = rcp(new DRT::ELEMENTS::Truss3(*crosslinkerdummy) );
#endif
*/
          
          //assigning correct global Id to new crosslinker element: since each node can have one crosslinker element
          //only at the same time a unique global Id can be found by taking the number of elemnts in the discretization
          //before starting dealing with crosslinkers and adding to this number the global Id of the node currently involved     
          newcrosslinker->SetId( basiselements_ + discret_.NodeRowMap()->GID(i) );
                   
          //nodes are assigned to the new crosslinker element by first assigning global node Ids and then assigning nodal pointers
          int globalnodeids[] = {discret_.NodeRowMap()->GID(i), neighbour};      
          newcrosslinker->SetNodeIds(2, globalnodeids);
          DRT::Node *nodes[] = {discret_.gNode( globalnodeids[0] ) , discret_.gNode( globalnodeids[1] )};
          newcrosslinker->BuildNodalPointers(&nodes[0]);
          
          //setting reference rotation coordinates
          LINALG::Matrix<6,1> rotrefe;
          for(int k = 0; k<3; k++)
          {
            rotrefe(k)   = dis[discret_.DofRowMap()->LID( discret_.Dof(discret_.lRowNode(i),k+3) )];
            rotrefe(k+3) = dis[discret_.DofRowMap()->LID( discret_.Dof(discret_.lRowNode(neighbour),k+3) )];
          }
          
          //correct reference configuration data is computed for the new crosslinker element
          newcrosslinker->SetUpReferenceGeometry(xrefe,rotrefe); 
          //set drag coefficient for new crosslinker element dependent on its length
          newcrosslinker->zeta_ = 4*PI*newcrosslinker->lrefe_*( statmechparams_.get<double>("ETA",0.0) );
          
          //add new element to discretization
          discret_.AddElement(newcrosslinker);  

          
          //noting that local node i has now a crosslinkerelement and noting global Id of this element
          (*crosslinkerpartner_)[i] = neighbour;

        }
      }
      
      /*settling administrative stuff in order to make the discretization ready for the next time step: the following
       * commmand generates ghost elements if necessary and calls FillCompete() method of discretization; note: this is
       * enough as long as only elements, but no nodes are added in a time step*/        
      DRT::UTILS::RedistributeWithNewNodalDistribution(discret_,*noderowmap,*nodecolmap);
      discret_.FillComplete(true,false,false);
      
      /*since graph of Crs matrices has changed stiff_ and damp_ are deleted completely and made ready for completely new assembly*/
      stiff_->Reset();
      damp_->Reset();
    }
    
  }
  
  #endif
  
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
