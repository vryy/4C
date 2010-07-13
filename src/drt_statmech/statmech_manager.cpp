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

#include <Teuchos_Time.hpp>

#include "statmech_manager.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"

//including random number library of blitz for statistical forces
#include <random/uniform.h>
#include <random/normal.h>


#ifdef D_BEAM3
#include "../drt_beam3/beam3.H"
#endif  // #ifdef D_BEAM3
#ifdef D_BEAM3II
#include "../drt_beam3ii/beam3ii.H"
#endif  // #ifdef D_BEAM3II
#ifdef D_TRUSS3
#include "../drt_truss3/truss3.H"
#endif  // #ifdef D_TRUSS3

#include "../drt_torsion3/torsion3.H"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <math.h>

//MEASURETIME activates measurement of computation time for certain parts of the code
//#define MEASURETIME

/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 09/08|
 *----------------------------------------------------------------------*/
StatMechManager::StatMechManager(ParameterList& params, DRT::Discretization& discret):
  statmechparams_( DRT::Problem::Instance()->StatisticalMechanicsParams() ),
  konswitch_(false),
  nsearch_(0),
  starttimeoutput_(-1.0),
  endtoendref_(0.0),
  istart_(0),
  rlink_(params.get<double>("R_LINK",0.0)),
  basisnodes_(discret.NumGlobalNodes()),
  basiselements_(discret.NumGlobalElements()),
  currentelements_(discret.NumGlobalElements()),
  outputfilenumber_(-1),
  discret_(discret)
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


    //information how many processors work at all
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

    //initialize crosslinkerpartner_ as a col map multivector consisting of as many vectors as a node can have crosslinkers at the maximum
    crosslinkerpartner_ = rcp(new Epetra_MultiVector(*(discret_.NodeColMap()),statmechparams_.get<int>("N_CROSSMAX",0)));
    crosslinkerpartner_->PutScalar(-1);
    numcrosslinkerpartner_ = rcp(new Epetra_Vector(*(discret_.NodeColMap())));
    numcrosslinkerpartner_->PutScalar(0);


    /*force sensors can be applied at any degree of freedom of the discretization the list of force sensors should
     * be based on a column map vector so that each processor has not only the information about each node's
     * displacement, but also about whether this has a force sensor; as a consequence each processor can write the
     * complete information gathered by all force sensors into a file of its own without any additional communication
     * with any other processor; initialization with -1 indicates that so far no forcesensors have been set*/
    forcesensor_ = rcp( new Epetra_Vector(*(discret_.DofColMap())) );
    forcesensor_->PutScalar(-1);

    /*since crosslinkers should be established only between different filaments the number of the filament
     * each node is belonging to is stored in the condition FilamentNumber; if no such conditions have been defined
     * the default -1 value is upkept in filamentnumber_ and crosslinkers between nodes belonging to the same filament
     * are allowed*/

    //getting a vector consisting of pointers to all filament number conditions set
    vector<DRT::Condition*> filamentnumberconditions(0);
    discret_.GetCondition("FilamentNumber",filamentnumberconditions);

    //next all the pointers to all the different conditions are looped
    for (int i=0; i<(int)filamentnumberconditions.size(); ++i)
    {
      //get filament number described by the current condition
      int filamentnumber = filamentnumberconditions[i]->GetInt("Filament Number") ;

      //get a pointer to nodal cloud covered by the current condition
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
     * an actin network between two rheometer plates the stiffness is to be determined this can be done by measuring the forces exerted
     * to the upper plate which is moving forwards and backwards for example; so for measurements of viscoelastic properties of materials
     * whithin a system certain sensor points have to be specified and in order to handle this within BACI these points are marked by
     * means of the condition sensorcondition*/

    //gettin a vector consisting of pointers to all filament number conditions set
    vector<DRT::Condition*> forcesensorconditions(0);
    discret_.GetCondition("ForceSensor",forcesensorconditions);

    //next all the pointers to all the different conditions are looped
    for (int i=0; i<(int)forcesensorconditions.size(); ++i)
    {
      //get number of nodal dof with respect to which force is to be measured; note: numbering starts with zero
      int nodedofnumber = forcesensorconditions[i]->GetInt("DOF Number") ;

      //get a pointer to nodal cloud covered by the current condition
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
         * freedom affected by this condition are marked in the vector *forcesensor_ by a one entry*/
        if(nodenumber > -1)
          (*forcesensor_)[dofnumber] = 1.0;
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

      //as soon as system is equilibrated (after time STARTTIME) a new file for storing output is generated
      if ( (time >= statmechparams_.get<double>("STARTTIME",0.0))  && (starttimeoutput_ == -1.0) )
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

      //as soon as system is equilibrated (after time STARTTIME) a new file for storing output is generated
      if ( (time >= statmechparams_.get<double>("STARTTIME",0.0))  && (starttimeoutput_ == -1.0) )
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
        if ( (istep - istart_) % statmechparams_.get<int>("OUTPUTINTERVALS",1) == 0 )
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
    /*computing and writing into file data about correlation of orientation of different elements
     *as it is considered in the context of the persistence length*/
    case INPAR::STATMECH::statout_orientationcorrelation:
    {

      FILE* fp = NULL; //file pointer for statistical output file

      //name of output file
      std::ostringstream outputfilename;
      outputfilename.str("");
      outputfilename << "OrientationCorrelation"<< outputfilenumber_ << ".dat";

      vector<double> arclength(discret_.NumMyRowElements(),0);
      vector<double> cosdiffer(discret_.NumMyRowElements(),0);

      //after initilization time write output cosdiffer in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps,
      //when discret_.NumMyRowNodes()-1 = 0,cosdiffer is always equil to 1!!
      if((time >= statmechparams_.get<double>("STARTTIME",0.0)) && (istep % statmechparams_.get<int>("OUTPUTINTERVALS",1)  == 0) )
      {

        Epetra_SerialDenseMatrix coord;
        coord.Shape(discret_.NumMyRowNodes(),ndim);

        for(int id = 0; id < discret_.NumMyRowNodes(); id++)
          for(int j = 0; j< ndim; j++)
            coord(id,j) = (discret_.gNode(id))->X()[j] + (dis)[(id)*(ndim-1)*3+j];

        for(int id = 0; id < discret_.NumMyRowElements(); id++)
        {

          //calculate the deformed length of every element
          for(int j = 0; j< ndim; j++)
          {
            arclength[id] += pow( (coord(id+1,j) - coord(id,j)), 2);
            cosdiffer[id] += (coord(id+1,j) - coord(id,j))*(coord(0+1,j) - coord(0,j));
          }

          //calculate the cosine difference referring to the first element
          //Dot product of the (id+1)th element with the 1st element and devided by the length of the (id+1)th element and the 1st element
          arclength[id] = pow(arclength[id],0.5);
          cosdiffer[id]= cosdiffer[id]/(arclength[id]*arclength[0]);
        }

        fp = fopen(outputfilename.str().c_str(), "a");
        std::stringstream filecontent;
        filecontent << istep;
        filecontent << scientific << setprecision(10);

        for(int id = 0; id < discret_.NumMyRowElements(); id++)
          filecontent<<" "<<cosdiffer[id];

        filecontent<<endl;
        fprintf(fp,filecontent.str().c_str());
        fclose(fp);

     }

    }
    break;
    //the following output allows for anisotropic diffusion simulation of a quasi stiff polymer
    case INPAR::STATMECH::statout_anisotropic:
    {
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
      if((time >= statmechparams_.get<double>("STARTTIME",0.0)) && (istep % statmechparams_.get<int>("OUTPUTINTERVALS",1)  == 0) )
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
          if(incangle > M_PI)
          {
            incangle -= 2*M_PI;
            incangle *= -1;
          }
          if(incangle < -M_PI)
          {
            incangle += 2*M_PI;
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
    }
    break;
    case INPAR::STATMECH::statout_viscoelasticity:
    {
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
      if( istep % statmechparams_.get<int>("OUTPUTINTERVALS",1)  == 0 )
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
          int count = 0;

          for(int i = 0; i < forcesensor_->MyLength(); i++)//changed
          {
            if( (*forcesensor_)[i] == 1.0)
            {
            	count++;
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
    }
    break;
    //writing data for generating a Gmsh video of the simulation
    case INPAR::STATMECH::statout_gmsh:
    {
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
      if( istep % statmechparams_.get<int>("OUTPUTINTERVALS",1)  == 0 )
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
    }
    break;
    //writing
    case INPAR::STATMECH::statout_structpolymorph:
    {
      //output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
      if( istep % statmechparams_.get<int>("OUTPUTINTERVALS",1)  == 0 )
      {
        // first index = time step index
        std::ostringstream filename;

        filename << "./GmshOutput/StructPolymorph_"<<discret_.Comm().MyPID()<<"_"<<std::setw(6) << setfill('0') << istep <<".dat";

        StructPolymorphOutput(dis,filename);
      }
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
   * all the nodal position; note that no parallel operations have to be carried out in the following block as it
   * is processed by only one processor; note that one could get also a parallel gmsh output if all processors write
   * their output after each other into the gmsh output file as it is realized for example for solid contact problems*/
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
      //if( element->NumNode() > 2)
        //dserror("Gmsh output for two noded elements only");

			//preparing variable storing coordinates of all these nodes
			LINALG::SerialDenseMatrix coord(3,element->NumNode());
			for(int id = 0; id<3; id++)
				for(int jd = 0; jd<element->NumNode(); jd++)
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
				color = 0.5;

      //if no periodic boundary conditions are to be applied, we just plot the current element
      if(statmechparams_.get<double>("PeriodLength",0.0) == 0)
      {
      	// check whether the kinked visualization is to be applied
      	bool kinked = CheckForKinkedVisual(color, element->Id());

        const DRT::ElementType & eot = element->ElementType();
#ifdef D_BEAM3
      	if(eot==DRT::ELEMENTS::Beam3Type::Instance())
      	{
      		if(!kinked | color==1.0)
      		{
						for(int j=0; j<element->NumNode()-1; j++)
						{
							//writing element by nodal coordinates as a scalar line
							gmshfilecontent << "SL(" << scientific;
							gmshfilecontent<< coord(0,j) << "," << coord(1,j) << "," << coord(2,j) << ","
														 << coord(0,j+1) << "," << coord(1,j+1) << "," << coord(2,j+1) ;
							/*note: colors are chosen by values between 0.0 and 1.0. These values refer to a color vector which
							 * is given in a .geo-setup-file. If, for example, 5 colors are given(either in X11 color expressions or RGB),
							 * possible values are 0.0, 0.25, 0.5, 0.75, 1.0.
							 */
							gmshfilecontent << ")" << "{" << scientific << color << "," << color << "};" << endl;
						}
      		}
      		else
      		{
      			for(int j=0; j<element->NumNode()-1; j++)
      			{
							// calculate third point of the
							std::vector<double> p3(3,0.0);
							GmshKinkedVisual(coord, p3, element->Id());

							//visualization
							// (1/2)
							gmshfilecontent << "SL(" << scientific
															<< coord(0,j) << "," << coord(1,j) << "," << coord(2,j) << ","
															<< p3.at(0) << "," << p3.at(1) << "," << p3.at(2)
															<< ")" << "{" << scientific << 0.5 << "," << 0.5 << "};" << endl;
							//(2/2)
							gmshfilecontent << "SL(" << scientific
															<< p3.at(0) << "," << p3.at(1) << "," << p3.at(2)<<","
															<< coord(0,j+1) << "," << coord(1,j+1) << "," << coord(2,j+1)
															<< ")" << "{" << scientific << 0.5 << "," << 0.5 << "};" << endl;
      			}
      		}
      	}
        else
#endif
#ifdef D_BEAM3II
      	if(eot==DRT::ELEMENTS::Beam3iiType::Instance())
      	{
      		if(!kinked | color==1.0)
      		{
						for(int j=0; j<element->NumNode()-1; j++)
						{
							gmshfilecontent << "SL(" << scientific;
							gmshfilecontent<< coord(0,j) << "," << coord(1,j) << "," << coord(2,j) << ","
														 << coord(0,j+1) << "," << coord(1,j+1) << "," << coord(2,j+1) ;
							gmshfilecontent << ")" << "{" << scientific << color << "," << color << "};" << endl;
						}
      		}
      		else
      		{
						for(int j=0; j<element->NumNode()-1; j++)
						{
							std::vector<double> p3(3,0.0);
							GmshKinkedVisual(coord, p3, element->Id());

							gmshfilecontent << "SL(" << scientific
															<< coord(0,j) << "," << coord(1,j) << "," << coord(2,j) << ","
															<< p3.at(0) << "," << p3.at(1) << "," << p3.at(2)
															<< ")" << "{" << scientific << 0.5 << "," << 0.5 << "};" << endl;
							gmshfilecontent << "SL(" << scientific
															<< p3.at(0) << "," << p3.at(1) << "," << p3.at(2)<<","
															<< coord(0,j+1) << "," << coord(1,j+1) << "," << coord(2,j+1)
															<< ")" << "{" << scientific << 0.5 << "," << 0.5 << "};" << endl;
						}
      		}
      	}
        else
#endif
#ifdef D_TRUSS3
        if(eot==DRT::ELEMENTS::Truss3Type::Instance())
      	{
					for(int j=0; j<element->NumNode()-1; j++)
					{
						gmshfilecontent << "SL(" << scientific;
						gmshfilecontent<< coord(0,j) << "," << coord(1,j) << "," << coord(2,j) << ","
													 << coord(0,j+1) << "," << coord(1,j+1) << "," << coord(2,j+1) ;
						gmshfilecontent << ")" << "{" << scientific << color << "," << color << "};" << endl;
					}
      	}
        else
#endif
#ifdef D_TORSION3
        if(eot==DRT::ELEMENTS::Torsion3Type::Instance())
        {
          double beadcolor = 0.75;
          for(int j=0; j<element->NumNode(); j++)
          {
            gmshfilecontent << "SP(" << scientific;
            gmshfilecontent<< coord(0,j) << "," << coord(1,j) << "," << coord(2,j);
            gmshfilecontent << ")" << "{" << scientific << beadcolor << "," << beadcolor << "};" << endl;
          }
        }
        else
#endif
        {
        }
      }
      //in case of periodic boundary conditions we have to take care to plot correctly an element broken at some boundary plane
      else
        GmshOutputPeriodicBoundary(coord, color, gmshfilecontent, element->Id());

    }

    // plot (cubic) periodic box in case of periodic boundary conditions
    if(statmechparams_.get<double>("PeriodLength",0.0) > 0)
    {
    	// get current period length
    	double pl = statmechparams_.get<double>("PeriodLength",0.0);

			double boundarycolor=0.0;

    	// define boundary lines
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << 0.0 << "," << 0.0 << "," << 0.0 << ","
											<< pl << "," << 0.0 << "," << 0.0 ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 2
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << pl << "," << 0.0 << "," << 0.0 << ","
    									<< pl << "," << pl << "," << 0.0 ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 3
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << pl << "," << pl << "," << 0.0 << ","
    									<< pl << "," << pl << "," << pl ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 4
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << pl << "," << pl << "," << pl << ","
    									<< 0.0 << "," << pl << "," << pl ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 5
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << 0.0 << "," << pl << "," << pl << ","
    									<< 0.0 << "," << 0.0 << "," << pl ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 6
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << 0.0 << "," << 0.0 << "," << pl << ","
											<< 0.0 << "," << 0.0 << "," << 0.0 ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 7
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << 0.0 << "," << 0.0 << "," << 0.0 << ","
    									<< 0.0 << "," << pl << "," << 0.0 ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 8
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << 0.0 << "," << pl << "," << 0.0 << ","
    									<< pl << "," << pl << "," << 0.0 ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 9
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << 0.0 << "," << pl << "," << 0.0 << ","
    									<< 0.0 << "," << pl << "," << pl ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 10
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << pl << "," << 0.0 << "," << 0.0 << ","
    									<< pl << "," << 0.0 << "," << pl ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 11
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << pl << "," << 0.0 << "," << pl << ","
    									<< pl << "," << pl << "," << pl ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
    	// line 12
    	gmshfilecontent << "SL(" << scientific;
    	gmshfilecontent << pl << "," << 0.0 << "," << pl << ","
    									<< 0.0 << "," << 0.0 << "," << pl ;
    	gmshfilecontent << ")" << "{" << scientific << boundarycolor << "," << boundarycolor << "};" << endl;
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
 | gmsh output data in case of periodic boundary conditions             |
 |                                                    public)cyron 02/10|
 *----------------------------------------------------------------------*/
void StatMechManager::GmshOutputPeriodicBoundary(const LINALG::SerialDenseMatrix& coord, const double& color, std::stringstream& gmshfilecontent, int eleid)
{
  //number of spatial dimensions
  const int ndim = 3;

  // get Element Type of the first Element to determine the graphics output
  DRT::Element* element = discret_.gElement(eleid);

  // draw colored lines between two nodes of a beam3 or truss3 element (meant for filaments/crosslinks/springs)
  const DRT::ElementType & eot = element->ElementType();
  bool dotline = false;

#ifdef D_BEAM3
  if(element->ElementType().Name()=="Beam3Type")
  	dotline = eot==DRT::ELEMENTS::Beam3Type::Instance();
#endif
#ifdef D_BEAM3II
  if(element->ElementType().Name()=="Beam3iiType")
  	dotline = eot==DRT::ELEMENTS::Beam3iiType::Instance();
#endif
#ifdef D_TRUSS3
  if(element->ElementType().Name()=="Truss3Type")
  	dotline = dotline or eot==DRT::ELEMENTS::Truss3Type::Instance();
#endif

  // determine whether crosslink connects two filaments or occupies two binding spots on the same filament;
  // variable triggers different visualizations
  //bool kinked = CheckForKinkedVisual(color, element->Id());

  if (dotline)
  {
    /*detect and save in vector "cut", at which boundaries the element is broken due to periodic boundary conditions;
     * the entries of cut have the following meaning: 0: element not broken in respective coordinate direction, 1:
     * element broken in respective coordinate direction (node 0 close to zero boundary and node 1 close to boundary
     * at PeriodLength);  2: element broken in respective coordinate direction (node 1 close to zero boundary and node
     * 0 close to boundary at PeriodLength);*/
    LINALG::SerialDenseMatrix cut(3,(int)element->NumNode()-1,true);

    /* "coord" currently holds the shifted set of coordinates.
     * In order to determine the correct vector "dir" of the visualization at the boundaries,
     * a copy of "coord" with adjustments in the proper places is introduced*/
    LINALG::SerialDenseMatrix unshift = coord;

    for(int i=0; i<cut.N(); i++)
    {
      for(int dof=0; dof<ndim; dof++)
      {
        if( fabs(coord(dof,i+1) - statmechparams_.get<double>("PeriodLength",0.0) - coord(dof,i))  < fabs(coord(dof,i+1) - coord(dof,i)) )
        {
          cut(dof,i) = 1;
          unshift(dof,i+1) -= statmechparams_.get<double>("PeriodLength",0.0);
        }
        if( fabs(coord(dof,i+1) + statmechparams_.get<double>("PeriodLength",0.0) - coord(dof,i))  < fabs(coord(dof,i+1) - coord(dof,i)) )
        {
          cut(dof,i) = 2;
          unshift(dof,i+1) += statmechparams_.get<double>("PeriodLength",0.0);
        }
      }
    }

    // write special output for broken elements
    for(int i=0; i<cut.N() ; i++)
    {
      if(cut(0,i) + cut(1,i) + cut(2,i) > 0)
      {
        //compute direction vector between first(i-th) and second(i+1-th) node of element (normed):
        LINALG::Matrix<3,1> dir;
        double mod=0.0;
        for(int dof=0; dof<ndim; dof++)
        {
          dir(dof) = unshift(dof,i+1) - unshift(dof,i);
          mod += dir(dof)*dir(dof);
        }
        for(int dof=0; dof<ndim; dof++)
          dir(dof) /= mod;

        //from node 0 to nearest boundary where element is broken you get by vector X + lambda0*dir
        double lambda0 = dir.Norm2();
        for(int dof=0; dof<ndim; dof++)
        {
          if(cut(dof,i) == 1)
          {
            if(fabs( - coord(dof,i) / dir(dof)) < fabs(lambda0))
              lambda0 = - coord(dof,i) / dir(dof);
          }
          else if(cut(dof,i) == 2)
          {
            if( fabs( (statmechparams_.get<double>("PeriodLength",0.0) - coord(dof,i)) / dir(dof) ) < fabs(lambda0) )
              lambda0 = ( statmechparams_.get<double>("PeriodLength",0.0) - coord(dof,i)) / dir(dof);
          }
        }

        //from node 1 to nearest boundary where element is broken you get by vector X + lambda1*dir
        double lambda1 = dir.Norm2();
        for(int dof=0; dof<ndim; dof++)
        {
          if(cut(dof,i) == 2)
          {
            if(fabs( - coord(dof,i+1) / dir(dof) ) < fabs(lambda1))
              lambda1 = - coord(dof,i+1) / dir(dof);
          }
          else if(cut(dof,i) == 1)
          {
            if(fabs((statmechparams_.get<double>("PeriodLength",0.0) - coord(dof,i+1)) / dir(dof)) < fabs(lambda1))
              lambda1 = (statmechparams_.get<double>("PeriodLength",0.0) - coord(dof,i+1)) / dir(dof);
          }
        }

        //writing element by nodal coordinates as a scalar line
        gmshfilecontent << "SL(" << scientific;
        gmshfilecontent<< coord(0,i)  << "," << coord(1,i) << "," << coord(2,i) << ","
                       << coord(0,i) + lambda0*dir(0)  << "," << coord(1,i) + lambda0*dir(1)  << "," << coord(2,i) + lambda0*dir(2)  ;
        /*note: for each node there is one color variable for gmsh and gmsh finally plots the line
         * interpolating these two colors between the nodes*/
        gmshfilecontent << ")" << "{" << scientific << color << "," << color << "};" << endl;
        //writing element by nodal coordinates as a scalar line
        gmshfilecontent << "SL(" << scientific;
        gmshfilecontent<< coord(0,i+1)  << "," << coord(1,i+1) << "," << coord(2,i+1) << ","
                       << coord(0,i+1) + lambda1*dir(0)  << "," << coord(1,i+1) + lambda1*dir(1)  << "," << coord(2,i+1) + lambda1*dir(2)  ;
        /*note: for each node there is one color variable for gmsh and gmsh finally plots the line
         * interpolating these two colors between the nodes*/
        gmshfilecontent << ")" << "{" << scientific << color << "," << color << "};" << endl;
      }
      else	// output for continuous elements
      {
      	bool kinked = CheckForKinkedVisual(color, element->Id());
      	if(!kinked | color==1.0)
      	{
					// filament or crosslink between two filaments
					//writing element by nodal coordinates as a scalar line
					gmshfilecontent << "SL(" << scientific;
					gmshfilecontent<< coord(0,i) << "," << coord(1,i) << "," << coord(2,i) << ","
												 << coord(0,i+1) << "," << coord(1,i+1) << "," << coord(2,i+1) ;
					/*note: for each node there is one color variable for gmsh and gmsh finally plots the line
					 * interpolating these two colors between the nodes*/
					gmshfilecontent << ")" << "{" << scientific << color << "," << color << "};" << endl;
      	}
    		else
    		{

					std::vector<double> p3(3,0.0);
					GmshKinkedVisual(coord, p3, element->Id());

					gmshfilecontent << "SL(" << scientific
													<< coord(0,i) << "," << coord(1,i) << "," << coord(2,i) << ","
													<< p3.at(0) << "," << p3.at(1) << "," << p3.at(2)
													<< ")" << "{" << scientific << 0.5 << "," << 0.5 << "};" << endl;
					gmshfilecontent << "SL(" << scientific
													<< p3.at(0) << "," << p3.at(1) << "," << p3.at(2)<<","
													<< coord(0,i+1) << "," << coord(1,i+1) << "," << coord(2,i+1)
													<< ")" << "{" << scientific << 0.5 << "," << 0.5 << "};" << endl;
    		}
      }
    }
  }
#ifdef D_TORSION3
  // draw spheres at node positions ("beads" of the bead spring model)
  else if(eot==DRT::ELEMENTS::Torsion3Type::Instance())
  {
    double beadcolor = 0.75;
    for(int i=0; i<element->NumNode(); i++)
    {
      //writing element by nodal coordinates as a scalar line
      gmshfilecontent << "SP(" << scientific;
      gmshfilecontent<< coord(0,i) << "," << coord(1,i) << "," << coord(2,i);
      gmshfilecontent << ")" << "{" << scientific << beadcolor << "," << beadcolor << "};" << endl;
    }
  }
#endif
  return;
} // StatMechManager::GmshOutputPeriodicBoundary()


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
    /*computing and writing into file data about correlation of orientation of different elements
     *as it is considered in the context of the persistence length*/
    case INPAR::STATMECH::statout_orientationcorrelation:
    {
      FILE* fp = NULL; //file pointer for statistical output file

      //defining name of output file
      std::ostringstream outputfilename;

      outputfilenumber_ = 0;

      //file pointer for operating with numbering file
      FILE* fpnumbering = NULL;
      std::ostringstream numberingfilename;

      //look for a numbering file where number of already existing output files is stored:
      numberingfilename.str("NumberOfRealizationsOrientCorr");
      fpnumbering = fopen(numberingfilename.str().c_str(), "r");



      //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
      if(fpnumbering == NULL)
      {
        do
        {
          outputfilenumber_++;
          outputfilename.str("");
          outputfilename << "OrientationCorrelation"<< outputfilenumber_ << ".dat";
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
       outputfilename <<"OrientationCorrelation"<< outputfilenumber_ << ".dat";
     }

     //set up new file with name "outputfilename" without writing anything into this file
     fp = fopen(outputfilename.str().c_str(), "w");
     fclose(fp);

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

      //set up new file with name "outputfilename" without writing anything into this file
      fp = fopen(outputfilename.str().c_str(), "w");
      fclose(fp);

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
void StatMechManager::StatMechUpdate(const double dt, Epetra_Vector& disrow, RCP<LINALG::SparseOperator>& stiff, int ndim)
{

#ifdef MEASURETIME
    const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME

  /*first we modify the displacement vector so that current nodal position at the end of current time step compies with
   * periodic boundary conditions, i.e. no node lies outside a cube of edge length Hperiodic*/
  PeriodicBoundaryShift(disrow,ndim);

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

    /*search neighbours for each row node on this processor within column nodes*/
    const double t_search = Teuchos::Time::wallTime();

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
      currpos(1) = node->X()[1] + discol[discret_.DofColMap()->LID(dofnode[1])];
      currpos(2) = node->X()[2] + discol[discret_.DofColMap()->LID(dofnode[2])];
      //if node has also rotational degrees of freedom
      if(discret_.NumDof(node) == 6)
      {
        currrot(0) =                discol[discret_.DofColMap()->LID(dofnode[3])];
        currrot(1) =                discol[discret_.DofColMap()->LID(dofnode[4])];
        currrot(2) =                discol[discret_.DofColMap()->LID(dofnode[5])];
      }
      currentpositions[node->LID()] = currpos;
      currentrotations[node->LID()] = currrot;
    }


    //if(time_ - nsearch_*statmechparams_.get<double>("Delta_t_search",0.0) > statmechparams_.get<double>("Delta_t_search",0.0) )(for now, call for each time step)
    //new search for neighbour nodes after average time specified in input file
    RCP<Epetra_MultiVector> crosslinkerneighbours;
    SearchNeighbours(currentpositions,crosslinkerneighbours);

    cout << "\n***\nsearch time: " << Teuchos::Time::wallTime() - t_search<< " seconds\n***\n";


#ifdef MEASURETIME
    const double t_admin = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME

    //number of elements in this time step before adding or deleting any elements
    currentelements_ = discret_.NumGlobalElements();
    //set crosslinkers
    SetCrosslinkers(dt,noderowmap,nodecolmap,currentpositions,currentrotations,*crosslinkerneighbours);

    //delete crosslinkers
    DelCrosslinkers(dt,noderowmap,nodecolmap);


    /*settling administrative stuff in order to make the discretization ready for the next time step: synchronize
     *the Filled() state on all processors after having added or deleted elements by ChekcFilledGlobally(); then build
     *new element maps and call FillComplete(); finally Crs matrices stiff_ has to be deleted completely and made ready
     *for new assembly since their graph was changed*/
    discret_.CheckFilledGlobally();
    discret_.FillComplete(true,false,false);
    stiff->Reset();



#ifdef MEASURETIME
    cout << "\n***\nadministration time: " << Teuchos::Time::wallTime() - t_admin<< " seconds\n***\n";
#endif // #ifdef MEASURETIME


  }//if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))

#ifdef MEASURETIME
    const double Delta_t = Teuchos::Time::wallTime()-t_start;
    cout << "\n***\ntotal time: " << Delta_t<< " seconds\n***\n";
#endif // #ifdef MEASURETIME


  return;
} // StatMechManager::StatMechUpdate()

/*----------------------------------------------------------------------*
 | Shifts current position of nodes so that they comply with periodic   |
 | boundary conditions                                       cyron 04/10|
 *----------------------------------------------------------------------*/
void StatMechManager::PeriodicBoundaryShift(Epetra_Vector& disrow, int ndim)
{
	//only if period length >0 has been defined periodic boundary conditions are swithced on
  if(statmechparams_.get<double>("PeriodLength",0.0) > 0.0)
    for(int i = 0; i < discret_.NumMyRowNodes(); i++)
    {
      //get a pointer at i-th row node
      const DRT::Node* node = discret_.lRowNode(i);

      //get GIDs of this node's degrees of freedom
      std::vector<int> dofnode = discret_.Dof(node);

      for(int j = ndim - 1; j > -1; j--)
      {
        /*if node currently has coordinate value greater than statmechparams_.get<double>("PeriodLength",0.0),
         *it is shifted by -statmechparams_.get<double>("PeriodLength",0.0) to lie again in the domain*/
        if(node->X()[j] + disrow[discret_.DofRowMap()->LID(dofnode[j])] > statmechparams_.get<double>("PeriodLength",0.0))
        {
          disrow[discret_.DofRowMap()->LID(dofnode[j])] -= statmechparams_.get<double>("PeriodLength",0.0);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may fixed by DBC. To avoid problmes when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if(j == 2 && statmechparams_.get<int>("CURVENUMBER",-1) >=  1 && time_ > statmechparams_.get<double>("STARTTIME",0.0))
            disrow[discret_.DofRowMap()->LID(dofnode[statmechparams_.get<int>("OSCILLDIR",-1)])] -= statmechparams_.get<double>("SHEARAMPLITUDE",0.0)*DRT::Problem::Instance()->Curve(statmechparams_.get<int>("CURVENUMBER",-1)-1).f(time_);
        }
        /*if node currently has coordinate value smaller than zero, it is shifted by statmechparams_.get<double>("PeriodLength",0.0)
         *to lie again in the domain*/
        if(node->X()[j] + disrow[discret_.DofRowMap()->LID(dofnode[j])]< 0)
        {
          disrow[discret_.DofRowMap()->LID(dofnode[j])] += statmechparams_.get<double>("PeriodLength",0.0);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if(j == 2 && statmechparams_.get<int>("CURVENUMBER",-1) >=  1 && time_ > statmechparams_.get<double>("STARTTIME",0.0))
            disrow[discret_.DofRowMap()->LID(dofnode[statmechparams_.get<int>("OSCILLDIR",-1)])] += statmechparams_.get<double>("SHEARAMPLITUDE",0.0)*DRT::Problem::Instance()->Curve(statmechparams_.get<int>("CURVENUMBER",-1)-1).f(time_);
        }
      }
    }
}


/*------------------------------------------------------------------------*
 | This function loops through all the elements of the discretization and |
 | tests whether beam3 are broken by periodic boundary conditions in the  |
 | reference configuration; if yes initial values of curvature and jacobi |
 | determinants are adapted in a proper way                    cyron 02/10|
 *-----------------------------------------------------------------------*/
void StatMechManager::PeriodicBoundaryBeam3Init(DRT::Element* element)
{

#ifdef D_BEAM3

  DRT::ELEMENTS::Beam3* beam = dynamic_cast<DRT::ELEMENTS::Beam3*>(element);

  //3D beam elements are embeddet into R^3:
  const int ndim = 3;

  /*get reference configuration of beam3 element in proper format for later call of SetUpReferenceGeometry;
   * note that rotrefe for beam3 elements is related to the entry in the global total Lagrange displacement
   * vector related to a certain rotational degree of freedom; as the displacement is initially zero also
   * rotrefe is set to zero here*/
  vector<double> xrefe(beam->NumNode()*ndim,0);
  vector<double> rotrefe(beam->NumNode()*ndim,0);


  for(int i=0;i<beam->NumNode();i++)
    for(int dof=0; dof<ndim; dof++)
    {
      xrefe[3*i+dof] = beam->Nodes()[i]->X()[dof];
      rotrefe[3*i+dof] = 0.0;
    }

  /*loop through all nodes except for the first node which remains fixed as reference node; all other nodes are
   * shifted due to periodic boundary conditions if required*/
  for(int i=1;i<beam->NumNode();i++)
  {
    for(int dof=0; dof<ndim; dof++)
    {
      /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
       * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
       * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
       * is smaller than half the periodic length*/
      if( fabs( (beam->Nodes()[i]->X()[dof]) + statmechparams_.get<double>("PeriodLength",0.0) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] += statmechparams_.get<double>("PeriodLength",0.0);

      if( fabs( (beam->Nodes()[i]->X()[dof]) - statmechparams_.get<double>("PeriodLength",0.0) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] -= statmechparams_.get<double>("PeriodLength",0.0);
    }
  }

  /*SetUpReferenceGeometry is a templated function; note that the third argument "true" is necessary as all beam elements
   * have already been initialized once upon reading input file*/
  switch(beam->NumNode())
  {
    case 2:
    {
      beam->SetUpReferenceGeometry<2>(xrefe,rotrefe,true);
      break;
    }
    case 3:
    {
      beam->SetUpReferenceGeometry<3>(xrefe,rotrefe,true);
      break;
    }
    case 4:
    {
      beam->SetUpReferenceGeometry<4>(xrefe,rotrefe,true);
      break;
    }
    case 5:
    {
      beam->SetUpReferenceGeometry<5>(xrefe,rotrefe,true);
      break;
    }
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
  }

#endif
}

/*------------------------------------------------------------------------*
 | This function loops through all the elements of the discretization and |
 | tests whether beam3 are broken by periodic boundary conditions in the  |
 | reference configuration; if yes initial values of curvature and jacobi |
 | determinants are adapted in a proper way                    cyron 02/10|
 *-----------------------------------------------------------------------*/
void StatMechManager::PeriodicBoundaryBeam3iiInit(DRT::Element* element)
{

#ifdef D_BEAM3II

  DRT::ELEMENTS::Beam3ii* beam = dynamic_cast<DRT::ELEMENTS::Beam3ii*>(element);

  //3D beam elements are embeddet into R^3:
  const int ndim = 3;

  /*get reference configuration of beam3ii element in proper format for later call of SetUpReferenceGeometry;
   * note that rotrefe for beam3ii elements is related to the entry in the global total Lagrange displacement
   * vector related to a certain rotational degree of freedom; as the displacement is initially zero also
   * rotrefe is set to zero here*/
  vector<double> xrefe(beam->NumNode()*ndim,0);
  vector<double> rotrefe(beam->NumNode()*ndim,0);


  for(int i=0;i<beam->NumNode();i++)
    for(int dof=0; dof<ndim; dof++)
    {
      xrefe[3*i+dof] = beam->Nodes()[i]->X()[dof];
      rotrefe[3*i+dof] = 0.0;
    }

  /*loop through all nodes except for the first node which remains fixed as reference node; all other nodes are
   * shifted due to periodic boundary conditions if required*/
  for(int i=1;i<beam->NumNode();i++)
  {
    for(int dof=0; dof<ndim; dof++)
    {
      /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
       * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
       * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
       * is smaller than half the periodic length*/
      if( fabs( (beam->Nodes()[i]->X()[dof]) + statmechparams_.get<double>("PeriodLength",0.0) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] += statmechparams_.get<double>("PeriodLength",0.0);

      if( fabs( (beam->Nodes()[i]->X()[dof]) - statmechparams_.get<double>("PeriodLength",0.0) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] -= statmechparams_.get<double>("PeriodLength",0.0);
    }
  }

  /*SetUpReferenceGeometry is a templated function; note that the third argument "true" is necessary as all beam elements
   * have already been initialized once upon reading input file*/
  switch(beam->NumNode())
  {
    case 2:
    {
      beam->SetUpReferenceGeometry<2>(xrefe,rotrefe,true);
      break;
    }
    case 3:
    {
      beam->SetUpReferenceGeometry<3>(xrefe,rotrefe,true);
      break;
    }
    case 4:
    {
      beam->SetUpReferenceGeometry<4>(xrefe,rotrefe,true);
      break;
    }
    case 5:
    {
      beam->SetUpReferenceGeometry<5>(xrefe,rotrefe,true);
      break;
    }
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
  }

#endif
}

/*------------------------------------------------------------------------*
 | This function loops through all the elements of the discretization and |
 | tests whether truss3 are broken by periodic boundary conditions in the |
 | reference configuration; if yes initial values of jacobi determinants  |
 | are adapted in a proper way                                 cyron 03/10|
 *-----------------------------------------------------------------------*/
void StatMechManager::PeriodicBoundaryTruss3Init(DRT::Element* element)
{

#ifdef D_TRUSS3

  DRT::ELEMENTS::Truss3* truss = dynamic_cast<DRT::ELEMENTS::Truss3*>(element);

  //3D beam elements are embeddet into R^3:
  const int ndim = 3;

  /*get reference configuration of truss3 element in proper format for later call of SetUpReferenceGeometry*/
  vector<double> xrefe(truss->NumNode()*ndim,0);

  for(int i=0;i<truss->NumNode();i++)
    for(int dof=0; dof<ndim; dof++)
      xrefe[3*i+dof] = truss->Nodes()[i]->X()[dof];


  /*loop through all nodes except for the first node which remains fixed as reference node; all other nodes are
   * shifted due to periodic boundary conditions if required*/
  for(int i=1;i<truss->NumNode();i++)
  {
    for(int dof=0; dof<ndim; dof++)
    {
      /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
       * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
       * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
       * is smaller than half the periodic length*/
      if( fabs( (truss->Nodes()[i]->X()[dof]) + statmechparams_.get<double>("PeriodLength",0.0) - (truss->Nodes()[0]->X()[dof]) ) < fabs( (truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] += statmechparams_.get<double>("PeriodLength",0.0);

      if( fabs( (truss->Nodes()[i]->X()[dof]) - statmechparams_.get<double>("PeriodLength",0.0) - (truss->Nodes()[0]->X()[dof]) ) < fabs( (truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] -= statmechparams_.get<double>("PeriodLength",0.0);
    }
  }

  /*note that the third argument "true" is necessary as all truss elements have already been initialized once upon reading input file*/
  truss->SetUpReferenceGeometry(xrefe,true);


#endif
}

/*----------------------------------------------------------------------*
 | Searches and saves in variable  crosslinkerneighbours neighbours for|
 | each col node                                              cyron 07/09|
 *----------------------------------------------------------------------*/
void StatMechManager::SearchNeighbours(const std::map<int,LINALG::Matrix<3,1> > currentpositions, RCP<Epetra_MultiVector>& crosslinkerneighbours)
{
  //update the number of times the function SearchNeighbours has already been called
  nsearch_++;

  //variable for local storage of neighbours of each nodes (initialized with with empty vector for each row node)
  std::vector< std::vector<int> > crosslinkerneighboursloc;
  std::vector<int> emptyvector;
  crosslinkerneighboursloc.resize(discret_.NumMyRowNodes(),emptyvector);

  //mean distance bridged by a crosslinker
  double rlink = statmechparams_.get<double>("R_LINK",0.0);
  //absolute value of difference between maximal/minimal and mean distance bridged by crosslinkers
  double deltarlink = statmechparams_.get<double>("DeltaR_LINK",0.0);

  //maximal number of numbers of any nodes locally on this processor and globally
  int maxneighbourslocal  = 0;
  int maxneighboursglobal = 0;


  //each processor looks for each of its row nodes for neighbours; loop index i is the local row node Id
  for(int i = 0; i < discret_.NumMyRowNodes(); i++)
  {
    //get GID and column map LID of row node i
    int iGID = (discret_.lRowNode(i))->Id();
    int icolID = (discret_.gNode(iGID))->LID();


    /*getting iterator to current position of local node icolID from map currentpositions;
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

      if(difference.Norm2() < rlink + deltarlink && difference.Norm2() > rlink - deltarlink )
        neighboursLID.push_back(posj->first);
    }

    //(B): seraching neighbours with search tree
    //neighboursLID = octTree_->searchPointsInRadius(currentpositions,(currentpositions.find(iLID))->second,rlink);

    /* after having searched all nodes within distance rlink around row node i we delete those
     * neighbours, which do not comply with certain requirements; here we establish the following
     * requirement: if filament numbering is activated (i.e. not all filament numbers are set to -1) the
     * neighbour node should belong to a filament different from the one the searching node belongs to;
     * note:
     * using erase you have to be careful to keep your iterator valid despite conditional deleting of elements
     * during iteration. The following algorithms represents a very efficient and simple way to deal with this
     * problem in a correct manner*/
    /*
    vector<int>::iterator iter = neighboursLID.begin();
    while( iter != neighboursLID.end() )
    {
      if( (  (*filamentnumber_)[*iter] == (*filamentnumber_)[icolID]  && (*filamentnumber_)[icolID] >= 0 )
         ||  (discret_.lColNode(*iter))->Id() < iGID )
        iter = neighboursLID.erase(iter);
      else
        ++iter;

    }
    */

    //finally the list of column map LIDs in neighboursLID is assigned to the entry of the i-th row node in crosslinkerneighbours_
    crosslinkerneighboursloc[i] = neighboursLID;


    //update maximal number of neighbours if this node has now more than any other before on this processor
    maxneighbourslocal = max(maxneighbourslocal,(int)crosslinkerneighboursloc[i].size());
  }


  //get global maximal number of neighbours of any node
  discret_.Comm().MaxAll(&maxneighbourslocal,&maxneighboursglobal,1);


  //write information in crosslinkerneighboursloc into an Epetra_MultiVector
  Epetra_MultiVector crosslinkerneighboursrow(*(discret_.NodeRowMap()),maxneighboursglobal);
  crosslinkerneighboursrow.PutScalar(-1);
  for(int i = 0; i < discret_.NumMyRowNodes(); i++)
    for(int j = 0; j < (int)(crosslinkerneighboursloc[i]).size(); j++)
      crosslinkerneighboursrow[j][i] = (crosslinkerneighboursloc[i])[j];


  //export information in crosslinkerneighboursrow to col map variable crosslinkerneighbours
  crosslinkerneighbours = rcp(new Epetra_MultiVector(*(discret_.NodeColMap()),maxneighboursglobal));


  Epetra_Import importer(*(discret_.NodeColMap()),*(discret_.NodeRowMap()));
  crosslinkerneighbours->Import(crosslinkerneighboursrow,importer,Add);



}//void SearchNeighbours(const int rlink, const std::map<int,LINALG::Matrix<3,1> > currentpositions)

/*----------------------------------------------------------------------*
 | set crosslinkers between neighbouring nodes after probability check  |
 | (private)                                                 cyron 02/09|
 *----------------------------------------------------------------------*/
void StatMechManager::SetCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
                                      const std::map<int,LINALG::Matrix<3,1> >& currentpositions,
                                      const std::map<int,LINALG::Matrix<3,1> >& currentrotations,
                                      Epetra_MultiVector& crosslinkerneighbours)
{
  //get current on-rate for crosslinkers
  double kon = 0;
  if( (currentelements_ - basiselements_) < statmechparams_.get<double>("N_crosslink",0.0) && !konswitch_)
    kon = statmechparams_.get<double>("K_ON_start",0.0);
  else
  {
    kon = statmechparams_.get<double>("K_ON_end",0.0);
    konswitch_=true;
  }

  //probability with which a crosslinker is established between neighbouring nodes
  double plink = 1.0 - exp( -dt*kon*statmechparams_.get<double>("C_CROSSLINKER",0.0) );

  //creating a random generator object which creates uniformly distributed random numbers in [0;1]
  ranlib::UniformClosed<double> UniformGen;


  /*if FIXEDSEED is set to Yes, reproducible simulations are required. Note: Blitz uses a different
   *seed for each type of random generator (e.g. normal, uniform, ...); thus former seed e.g. for
   *normally distributed random generators are not sufficient to guarantee a seed for the generator
   *UniformGen*/
  if(Teuchos::getIntegralValue<int>(statmechparams_,"FIXEDSEED"))
    UniformGen.seed((unsigned int)(time_/dt));


  //create Epetra_MultiVector which saves on proc 0 for each col map node to which nodes (GID) crosslinker is to be set
  Epetra_MultiVector crosslinkstobeadded(nodecolmap,statmechparams_.get<int>("N_CROSSMAX",0),true);
  crosslinkstobeadded.PutScalar(-1);

  //number of crosslinks to be set at this node to other nodes
  Epetra_Vector numcrosslinkstobeadded(nodecolmap,true);


  /*the following part of the code is executed on processor 0 only and leads to the decision, at which nodes crosslinks shall be set
   *processor 0 goes through all its column nodes in random order and checks whether a crosslink is to be set; this works precisely as follows:
   *(1) a random order is generated in which the nodes are went through
   *(2) if a node as well as its neighbor has not yet reached its maximal number of crosslinks, a crosslink may be set
   *(3) a crosslink is set if and only if the node passes a probability check
   *note: a fully overlapping node map on processor 0 is imperative for the following part of the code to work correctly!*/
  if(discret_.Comm().MyPID() == 0)
  {
    //random order is generated in which the nodes are went through
    std::vector<int> order(Permutation(nodecolmap.NumMyElements()));

    for(int i = 0; i< nodecolmap.NumMyElements(); i++)
    {
      int nodenumber = order[i];

      //go through all neighbors of this node in random order and look for suitable ones for crosslinking
      std::vector<int> neighboursorder(Permutation(crosslinkerneighbours.NumVectors()));
      int j = 0;

      //go through all neighbors (i.e. entries unequal -1) as long as current node has not yet reached maximal number of crosslinkers
      while(j < min(crosslinkerneighbours.NumVectors(),statmechparams_.get<int>("N_CROSSMAX",0)) && crosslinkerneighbours[neighboursorder[j]][nodenumber] != -1 && (*numcrosslinkerpartner_)[nodenumber] < statmechparams_.get<int>("N_CROSSMAX",0) )
      {

        //probability check whether at this node a crosslink can be set at all
        if(UniformGen.random() < plink)
        {

          //get LID of current neighbor node
          int neighbourLID = nodecolmap.LID((int)crosslinkerneighbours[neighboursorder[j]][nodenumber]);

          //array with nodal pointers to currently considered two nodes
          DRT::Node* nodes[2] = {discret_.gNode( nodecolmap.GID(nodenumber) ) , discret_.gNode( nodecolmap.GID(neighbourLID) )};

          //unit direction vector between currently considered two nodes
          LINALG::Matrix<3,1> direction((currentpositions.find(neighbourLID))->second);
          direction -= (currentpositions.find(nodenumber))->second;
          direction.Scale(1.0 / direction.Norm2());

          //a crosslinker can be set only if constraints wsith respect to its orientation are satisfied and also second node does not exceed it maximal number of crosslinkers
          if(CheckOrientation(direction,&nodes[0]) && (*numcrosslinkerpartner_)[neighbourLID] < statmechparams_.get<int>("N_CROSSMAX",0))
          {

            //only the node nodenumber stores the node to which a crosslink has to be set up, not vice versa, but for both the crosslinker counter is increased
            (*numcrosslinkerpartner_)[nodenumber]++;
            (*numcrosslinkerpartner_)[neighbourLID]++;
            /*this crosslinker is registered to be added only for one node and not for both nodes so that later only one crosslinker is set between the two nodes;
             *it does not matter which one of the two nodes registers the crosslinker to be setted; this decision was made here arbitrarily*/
            numcrosslinkstobeadded[nodenumber]++;
            crosslinkstobeadded[(int)(numcrosslinkstobeadded[nodenumber] - 1)][nodenumber] = nodecolmap.GID(neighbourLID);

            //we check whether a crosslinker has already been established between the same nodes; if yes nothing is done
            bool notyet = true;
            for(int k = 0; k < statmechparams_.get<int>("N_CROSSMAX",0); k++)
              if(nodecolmap.GID(neighbourLID) == (*crosslinkerpartner_)[k][nodenumber] )
                notyet = false;

            /*these two nodes are not yet crosslinked the link will actually be set and this is registered in crosslinkerpartner_ for
             *the two linked nodes in the first currently unused multivector column, respectively*/
            bool registered1 = false;
            bool registered2 = false;
            if(notyet)
              for(int k = 0; k < statmechparams_.get<int>("N_CROSSMAX",0); k++)
              {
                if((*crosslinkerpartner_)[k][nodenumber] == -1 && registered1 == false)
                {
                  (*crosslinkerpartner_)[k][nodenumber] = nodecolmap.GID(neighbourLID);
                  registered1 = true;
                }

                if((*crosslinkerpartner_)[k][neighbourLID] == -1 && registered2 == false)
                {
                  (*crosslinkerpartner_)[k][neighbourLID] = nodecolmap.GID(nodenumber);
                  registered2 = true;
                }
              }

          }
        }
        j++;
      }

    }
  }//if(discret_.NumProc() == 0)
  //set to zero on all other processors to allow for simple additive export of information generated by processor 0 later
  else
  {
    crosslinkstobeadded.PutScalar(0);
    numcrosslinkstobeadded.PutScalar(0);
    numcrosslinkerpartner_->PutScalar(0);
    crosslinkerpartner_->PutScalar(0);
  }


  //synchronize information about crosslinks to be added by exporting it to row map formant and then reimporting it to column map format
  Epetra_Export exporter(nodecolmap,noderowmap);
  Epetra_Import importer(nodecolmap,noderowmap);
  Epetra_MultiVector crosslinkstobeaddedrow(noderowmap,crosslinkstobeadded.NumVectors(),true);
  Epetra_Vector numcrosslinkstobeaddedrow(noderowmap,true);
  Epetra_MultiVector crosslinkerpartnerrow(noderowmap,crosslinkerpartner_->NumVectors(),true);
  Epetra_Vector numcrosslinkerpartnerrow(noderowmap,true);

  crosslinkstobeaddedrow.Export(crosslinkstobeadded,exporter,Add);
  crosslinkstobeadded.Import(crosslinkstobeaddedrow,importer,Insert);
  numcrosslinkstobeaddedrow.Export(numcrosslinkstobeadded,exporter,Add);
  numcrosslinkstobeadded.Import(numcrosslinkstobeaddedrow,importer,Insert);
  numcrosslinkerpartnerrow.Export(*numcrosslinkerpartner_,exporter,Add);
  numcrosslinkerpartner_->Import(numcrosslinkerpartnerrow,importer,Insert);
  crosslinkerpartnerrow.Export(*crosslinkerpartner_,exporter,Add);
  crosslinkerpartner_->Import(crosslinkerpartnerrow,importer,Insert);


  /*now it has been determined which crosslinkers are to be added; a crosslinker is added on each processor which is row node owner of at least one of its nodes;
   *the processor which is row map owner of the node with the smalles GID will be the owner of the new crosslinker element; on the other processors the new
   *crosslinker element will be present as a ghost element only*/
  for(int i = 0; i < crosslinkstobeadded.MyLength(); i++)
  {
    //getting current position and rotational status of node with column map LID i from map currentpositions
    map< int,LINALG::Matrix<3,1> >::const_iterator posi         = currentpositions.find(i);
    map< int,LINALG::Matrix<3,1> >::const_iterator roti         = currentrotations.find(i);

    //loop through all node GIDs to which a crosslink should be established from column map node i
    for(int j = 0; j < crosslinkstobeadded.NumVectors(); j++)
    {
      /*a crosslinker is to be established if the following two conditions are satisfied:
       *(1) there is a node stored to which it is to be established (i.e. the i-th element of the j-th vector in in the Multi_Vector crosslinkerstobeaddedglobalcol is unequal the default value -1
       *(2) second: one of the two nodes between which the crosslinker should be established is a row node on this processor; otherwise one does not need to add the crosslinker element on this processor*/
      if(crosslinkstobeadded[j][i] > -0.9 &&
         ( noderowmap.LID((nodecolmap.GID(i))) > -0.9 || noderowmap.LID(((int)crosslinkstobeadded[j][i])) > -0.9 ) )
      {

        //getting current position and rotational status of node with GID crosslinkerstobeaddedglobalcol[j][i]
        map< int,LINALG::Matrix<3,1> >::const_iterator posj = currentpositions.find( nodecolmap.LID((int)crosslinkstobeadded[j][i]) );
        map< int,LINALG::Matrix<3,1> >::const_iterator rotj = currentrotations.find( nodecolmap.LID((int)crosslinkstobeadded[j][i]) );



        /*there is the problem of how to assign to a crosslinker element a GID which is certainly not used by any
         * other element; we know that for a network at the beginning (without crosslinkers) each node has a
         * connectivity of 1 or 2 so that the number of all elements used to discretize the filaments is smaller
         * than basisnodes_. Thus we choose crosslinker GIDs >= basisnodes_ for the additionally added crosslinker
         * elements. To make sure that two crosslinkers never have the same GID we add to basisnodes_ the value
         * basisnodes_*GID1, where GID1 is the GID of the first node of the crosslinker element. Then we add GID2
         * (note: GID2 < basisnodes_), which is the GID of the second node of the crosslinker element.
         * Hence basisnodes_ + GID1*basisnodes_ + GID2 always gives a GID which cannot be used by any other element;
         * the first node of the crosslinker element is always assumed to be the one with the greater GID; the owner
         * of the first node is assumed to be the owner of the crosslinker element*/
        int GID2 = min(nodecolmap.GID(i),(int)crosslinkstobeadded[j][i]);
        int GID1 = max(nodecolmap.GID(i),(int)crosslinkstobeadded[j][i]);
        int newcrosslinkerGID = (GID1 + 1)*basisnodes_ + GID2;
        //create nodal pointers for crosslinker element
        int globalnodeids[2] = {GID1,GID2};
        DRT::Node* nodes[2] = {discret_.gNode( globalnodeids[0] ) , discret_.gNode( globalnodeids[1] )};

        //save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
        std::vector<double> rotrefe(6);
        std::vector<double> xrefe(6);

        for(int k = 0; k<3; k++)
        {
          //the first three positions in xrefe and rotrefe are for the data of node GID1
          if(nodecolmap.GID(i) > (int)crosslinkstobeadded[j][i])
          {
            //set nodal positions
            xrefe[k  ] = (posi->second)(k);
            xrefe[k+3] = (posj->second)(k);

            //set nodal rotations (not true ones, only those given in the displacment vector)
            rotrefe[k  ] = (roti->second)(k);
            rotrefe[k+3] = (rotj->second)(k);
          }
          else
          {
            //set nodal positions
            xrefe[k  ] = (posj->second)(k);
            xrefe[k+3] = (posi->second)(k);

            //set nodal rotations (not true ones, only those given in the displacment vector)
            rotrefe[k  ] = (rotj->second)(k);
            rotrefe[k+3] = (roti->second)(k);
          }
        }

 #ifdef D_TRUSS3
 #ifdef D_BEAM3
 #ifdef D_BEAM3II

        if(statmechparams_.get<double>("ILINK",0.0) > 0.0)
        {
          RCP<DRT::ELEMENTS::Beam3> newcrosslinker = rcp(new DRT::ELEMENTS::Beam3(newcrosslinkerGID,(discret_.gNode(GID1))->Owner() ) );

          newcrosslinker->SetNodeIds(2,globalnodeids);
          newcrosslinker->BuildNodalPointers(&nodes[0]);

          //setting up crosslinker element parameters and reference geometry
          newcrosslinker ->crosssec_ = statmechparams_.get<double>("ALINK",0.0);
          newcrosslinker ->crosssecshear_ = 1.1*statmechparams_.get<double>("ALINK",0.0);
          newcrosslinker ->Iyy_ = statmechparams_.get<double>("ILINK",0.0);
          newcrosslinker ->Izz_ = statmechparams_.get<double>("ILINK",0.0);
          newcrosslinker ->Irr_ = statmechparams_.get<double>("IPLINK",0.0);
          newcrosslinker->SetMaterial(2);

          //set up reference configuration of crosslinker
          newcrosslinker->SetUpReferenceGeometry<2>(xrefe,rotrefe);

          //add element to discretization
          discret_.AddElement(newcrosslinker);
        }
        else
        {
          RCP<DRT::ELEMENTS::Truss3> newcrosslinker = rcp(new DRT::ELEMENTS::Truss3(newcrosslinkerGID, (discret_.gNode(GID1))->Owner() ) );

          newcrosslinker->SetNodeIds(2,globalnodeids);
          newcrosslinker->BuildNodalPointers(&nodes[0]);

          //setting up crosslinker element parameters
          newcrosslinker ->crosssec_ = statmechparams_.get<double>("ALINK",0.0);
          newcrosslinker->SetMaterial(2);

          //correct reference configuration data is computed for the new crosslinker element;
          newcrosslinker->SetUpReferenceGeometry(xrefe);

          //add new element to discretization
          discret_.AddElement(newcrosslinker);
        }
 #endif
 #endif
 #endif

      }
    }
  }
}//void SetCrosslinkers(const Epetra_Vector& setcrosslinkercol)

/*----------------------------------------------------------------------*
 | delete crosslinkers listed in crosslinkerpartner_ after random check;|
 | the random check is conducted by the owner processor of the          |
 | crosslinker, respectively, and the result of that check is           |
 | communicated subsequently to all the other processors                |
 | (private)                                                 cyron 11/09|
 *----------------------------------------------------------------------*/
void StatMechManager::DelCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap)
{
  //get current off-rate for crosslinkers
  double koff = 0;
  if( (currentelements_ - basiselements_) < statmechparams_.get<double>("N_crosslink",0.0) )
    koff = statmechparams_.get<double>("K_OFF_start",0.0);
  else
    koff = statmechparams_.get<double>("K_OFF_end",0.0);


  //probability with which a crosslink breaks up in the current time step
  double punlink = 1.0 - exp( -dt*koff );

  //creating a random generator object which creates uniformly distributed random numbers in [0;1]
  ranlib::UniformClosed<double> UniformGen;


  //create Epetra_MultiVector which saves on proc 0 for each col map node to which nodes (GID) crosslinker is to be deleted
  Epetra_MultiVector crosslinkstobedeleted(nodecolmap,statmechparams_.get<int>("N_CROSSMAX",0),true);
  crosslinkstobedeleted.PutScalar(-1);

  //number of crosslinks to be deleted at this node to other nodes
  Epetra_Vector numcrosslinkstobedeleted(nodecolmap,true);


  /*as already the decision of setting crosslinkers the decision of deleting them is made on processor 0 only*/
  if(discret_.Comm().MyPID() == 0)
  {
    /*go through all column map nodes (as deleting one crosslinker does not affect the decision whether another one may be deleted
     *this has not to be done in random order (unlike when setting crosslinkers)*/
    for(int i = 0; i< nodecolmap.NumMyElements(); i++)
      //go through all crosslink-partners
      for(int j = 0; j < statmechparams_.get<int>("N_CROSSMAX",0); j++)
        //is a crosslink registered in (*crosslinkerpartner_)[j][i]? If yes: is the probability check to delete it passed?
        if((*crosslinkerpartner_)[j][i] > -0.9 && UniformGen.random() < punlink)
        {

          //get LID of current neighbor node to which crosslink is to be deleted
          int neighbourLID = nodecolmap.LID((int)(*crosslinkerpartner_)[j][i]);


          //GID of crosslinker to be deleted
          int GID2 = min(nodecolmap.GID(i),(int)(*crosslinkerpartner_)[j][i]);
          int GID1 = max(nodecolmap.GID(i),(int)(*crosslinkerpartner_)[j][i]);
          int crosslinkerid = (GID1 + 1)*basisnodes_ + GID2;


          //decrease number of crosslinkers at both nodes where a crosslinker is deleted by one
          (*numcrosslinkerpartner_)[i]--;
          (*numcrosslinkerpartner_)[neighbourLID]--;

          //delete crosslinker from crosslinkerpartner_
          for(int k = 0; k < statmechparams_.get<int>("N_CROSSMAX",0); k++)
          {
            if((*crosslinkerpartner_)[k][i] == neighbourLID)
              (*crosslinkerpartner_)[k][i] = -1.0;

            if((*crosslinkerpartner_)[k][neighbourLID] == nodecolmap.GID(i))
              (*crosslinkerpartner_)[k][neighbourLID] = -1.0;
          }

          //register at either node that a crosslinker has to be deleted
          crosslinkstobedeleted[ (int)(numcrosslinkstobedeleted[i]) ][i] = crosslinkerid;
          numcrosslinkstobedeleted[i]++;
          crosslinkstobedeleted[ (int)(numcrosslinkstobedeleted[neighbourLID]) ][neighbourLID] = crosslinkerid;
          numcrosslinkstobedeleted[neighbourLID]++;

        }

  }//if(discret_.NumProc() == 0)
  //set col map variables to zero on all other processors to allow for simple additive export of information generated by processor 0 later
  else
  {
    crosslinkstobedeleted.PutScalar(0);
    numcrosslinkstobedeleted.PutScalar(0);
    numcrosslinkerpartner_->PutScalar(0);
    crosslinkerpartner_->PutScalar(0);
  }


  //synchronize information about crosslinks to be added by exporting it to row map formant and then reimporting it to column map format
  Epetra_Export exporter(nodecolmap,noderowmap);
  Epetra_Import importer(nodecolmap,noderowmap);
  Epetra_MultiVector crosslinkstobedeletedrow(noderowmap,crosslinkstobedeleted.NumVectors(),true);
  Epetra_Vector numcrosslinkstobedeletedrow(noderowmap,true);
  Epetra_MultiVector crosslinkerpartnerrow(noderowmap,crosslinkerpartner_->NumVectors(),true);
  Epetra_Vector numcrosslinkerpartnerrow(noderowmap,true);

  crosslinkstobedeletedrow.Export(crosslinkstobedeleted,exporter,Add);
  crosslinkstobedeleted.Import(crosslinkstobedeletedrow,importer,Insert);
  numcrosslinkstobedeletedrow.Export(numcrosslinkstobedeleted,exporter,Add);
  numcrosslinkstobedeleted.Import(numcrosslinkstobedeletedrow,importer,Insert);
  numcrosslinkerpartnerrow.Export(*numcrosslinkerpartner_,exporter,Add);
  numcrosslinkerpartner_->Import(numcrosslinkerpartnerrow,importer,Insert);
  crosslinkerpartnerrow.Export(*crosslinkerpartner_,exporter,Add);
  crosslinkerpartner_->Import(crosslinkerpartnerrow,importer,Insert);


  /*at this point the information which crosslinkers are to be deleted to a certain node is present on each processor which knows a certain node at least in
   *its column map; now each processor loops through all column map nodes and checks for the crosslinker to be deleted whether it exists on this processor (at
   *least as a ghost element); if yes, the crosslinker element is actually deleted*/
  for(int i = 0; i < crosslinkstobedeleted.MyLength(); i++)
  {
    for(int j = 0; j < numcrosslinkstobedeleted[i]; j++)
    {
      //GID of crosslinker to be deleted
      int crosslinkerid = (int)crosslinkstobedeleted[j][i];

      //delete crosslinker if present on this processosr
      if(discret_.HaveGlobalElement(crosslinkerid))
        discret_.DeleteElement(crosslinkerid);

    }
  }



}//void DelCrosslinkers(const Epetra_Vector& delcrosslinkercol)

/*----------------------------------------------------------------------*
 | (public) generate gaussian randomnumbers with mean "meanvalue" and   |
 | standarddeviation "standarddeviation" for parallel use     cyron10/09|
 *----------------------------------------------------------------------*/
void StatMechManager::GenerateGaussianRandomNumbers(RCP<Epetra_MultiVector> randomnumbers,const double meanvalue, const double standarddeviation)
{
  //multivector for stochastic forces evaluated by each element based on row map
  Epetra_MultiVector randomnumbersrow( *(discret_.ElementRowMap()),randomnumbers->NumVectors());

  //creating a random generator object which creates random numbers with zero mean and unit standard deviation using the Blitz routine
  ranlib::Normal<double> normalGen(meanvalue,standarddeviation);

  for(int i=0; i<randomnumbersrow.MyLength(); i++)
    for(int j=0; j<randomnumbersrow.NumVectors(); j++)
      randomnumbersrow[j][i] = normalGen.random();

  //export stochastic forces from row map to column map
  Epetra_Export exporter(*discret_.ElementRowMap(),*discret_.ElementColMap());
  randomnumbers->Export(randomnumbersrow,exporter,Add);

  //now fstoch contains stochastic forces identical on all processors

  return;
} // StatMechManager::SynchronizeRandomForces()


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

  /*note: the variable crosslinkerpartner_  is not saved here; this means
   * that for crosslinked networks restarts are not possible in general;
   * However, not saving crosslinkerpartner_ makes a reasonable restart impossible*/

  return;
} // StatMechManager::StatMechWriteRestart()

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

/*----------------------------------------------------------------------*
 | check for broken element                        (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void StatMechManager::CheckForBrokenElement(LINALG::SerialDenseMatrix& coord, LINALG::SerialDenseMatrix& cut, bool *broken)
{
	// empty cut just in case it was handed over non-empty
	cut.Zero();
	*broken = false;
	int ndim = coord.M();
	// flag "broken" signals a broken element, flag "cut" hints at location of nodes 0 and 1 (of the two-noded beam)
  for(int dof=0; dof<ndim; dof++)
  	//loop through columns of "cut"
  	for(int n=0; n<cut.N(); n++)
  	{
  		// broken element with node_n close to "0.0"-boundary
			if( fabs(coord(dof,n+1) - statmechparams_.get<double>("PeriodLength",0.0) - coord(dof,n))  < fabs(coord(dof,n+1)  - coord(dof,n)) )
			{
				*broken = true;
				// set value for the spatial component in question at n-th cut
				cut(dof,n) = 1.0;
			}
			else if( fabs(coord(dof,n+1) + statmechparams_.get<double>("PeriodLength",0.0) - coord(dof,n))  < fabs(coord(dof,n+1) - coord(dof,n)) )
			{
				*broken = true;
				cut(dof,n) = 2.0;
			}
  	}
  return;
}// StatMechManager::CheckForBrokenElement

/*----------------------------------------------------------------------*
 | get a matrix with node coordinates and their DOF LIDs                |
 |																							   (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void StatMechManager::GetElementNodeCoords(DRT::Element* element, RCP<Epetra_Vector> dis, LINALG::SerialDenseMatrix& coord, vector<int>* lids)
{
	// clear LID vector just in case it was handed over non-empty
	lids->clear();
	for(int j=0; j<element->NumNode();j++)
	{
		int numdof = discret_.Dof(0,element->Nodes()[j]).size();
		for(int k = 0; k<numdof; k++)
	 	{
			// obtain k-th spatial component of the reference position of the j-th node
	 		double referenceposition = ((element->Nodes())[j])->X()[k];
	 		// get the GIDs of the node's DOFs
	 	  vector<int> dofnode = discret_.Dof((element->Nodes())[j]);
	 	  // store the displacement of the k-th spatial component
	 	  double displacement = (*dis)[discret_.DofRowMap()->LID( dofnode[k] )];
	 	  // write updated components into coord (only translational
	 	  coord(k,j) =  referenceposition + displacement;
	 	  // store current lid(s) (3 translational DOFs per node)
	 	  if(lids!=NULL)
	 	  	lids->push_back(discret_.DofRowMap()->LID( dofnode[k] ));
	 	}
	}
	return;
} // StatMechManager::GetElementNodeCoords

/*----------------------------------------------------------------------*
 | update force sensor locations                   (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void StatMechManager::UpdateForceSensors(vector<int>& sensornodes, int oscdir)
{
	// reinitialize forcesensor_, so that tedious looping
	// in order to cancel previous force sensor positions is avoided
	forcesensor_->PutScalar(-1.0);

	// loop over DOFs subjected to oscillation (by DBC)
	for(int i=0; i<(int)sensornodes.size(); i++)
	{
		// check if node is on proc (if not, continue)
		bool havenode = discret_.HaveGlobalNode(sensornodes.at(i));
		if(!havenode)
			continue;

		// get the node
		DRT::Node* actnode = discret_.gNode(sensornodes.at(i));
		// get the GID of the DOF of the oscillatory motion
		int gid = discret_.Dof(0,actnode)[oscdir];
		// now, get the LID
		int lid = discret_.DofRowMap()->LID(gid);
		// activate force sensor at lid-th position
		(*forcesensor_)[lid] = 1.0;
	}
} // StatMechManager::UpdateForceSensors

/*----------------------------------------------------------------------*
 | generates a vector with a random permutation of the numbers between 0|
 | and N - 1                                        (public) cyron 06/10|
 *----------------------------------------------------------------------*/
std::vector<int> StatMechManager::Permutation(const int& N)
{
  //auxiliary variable
  int j = 0;

  //result vector initialized with ordered numbers from 0 to N-1
  std::vector<int> result(N,0);
  for(int i=0; i<(int)result.size(); i++)
    result[i] = i;

  for (int i = 0; i < N; ++i)
  {
    //generate random number between 0 and i
     j = rand() % (i + 1);

     /*exchange values at positions i and j (note: value at position i is i due to above initialization
      *and because so far only positions <=i have been changed*/
     result[i] = result[j];
     result[j] = i;
  }

  return result;

} // StatMechManager::Permutation

/*----------------------------------------------------------------------*
 | checks orientation of crosslinker relative to linked filaments       |
 |                                                  (public) cyron 06/10|
 *----------------------------------------------------------------------*/

bool StatMechManager::CheckOrientation(const LINALG::Matrix<3,1> direction, DRT::Node** nodes)
{

#ifdef D_BEAM3
#ifdef D_BEAM3II
  //creating a random generator object which creates uniformly distributed random numbers in [0;1]
  ranlib::UniformClosed<double> UniformGen;

  //if orientation is not to be checked explicitly, this function always returns true
  if(!Teuchos::getIntegralValue<int>(statmechparams_,"CHECKORIENT"))
    return 1;


  //triad of node on filament which is affected by the new crosslinker
  LINALG::Matrix<3,3> Tfil;

  //direction vectors of filaments projected to plane perpendicular to crosslinker axis and normalized to absolute value one
  vector<LINALG::Matrix<3,1> > filprojected;
  filprojected.resize(2);

  //angles between crosslinker axis and filament axes minus equilibrium angles phi0, respectively
  vector<double> Deltaphi;
  Deltaphi.resize(2);

  //angle between filament axes projected into the plane perpendicular to the crosslinker axis minus equilibrium angles tau0
  double Deltatau;


  //bending and torsion stiffness of crosslinkers
  double EI = statmechparams_.get<double>("ELINK",0.0)*statmechparams_.get<double>("ILINK",0.0);
  double GJ = statmechparams_.get<double>("ELINK",0.0)*statmechparams_.get<double>("IPLINK",0.0);


  for(int i=0; i<2; i++)
  {
    //lowest Id of any connected element (the related element cannot be a crosslinker, but has to belong to the actual filament discretization)
    int lowestid( ( (nodes[i]->Elements())[0] )->Id() );
    int lowestidele(0);
    for(int j=0; j<nodes[i]->NumElement(); j++)
      if( ( (nodes[i]->Elements())[j] )->Id() < lowestid )
      {
        lowestid = ( (nodes[i]->Elements())[j] )->Id();
        lowestidele = j;
      }

    //check whether filaments are discretized with beam3ii elements (otherwise no orientation triads at nodes available)
    DRT::ELEMENTS::Beam3ii* filele;
    DRT::ElementType & eot = ((nodes[i]->Elements())[lowestidele])->ElementType();
    if(eot == DRT::ELEMENTS::Beam3iiType::Instance())
      filele = dynamic_cast<DRT::ELEMENTS::Beam3ii*>(nodes[i]->Elements()[lowestidele]);
    else
      dserror("Filaments have to be discretized with beam3ii elements for orientation check!!!");

    //check whether crosslinker is connected to first or second node of that element
    int nodenumber = 0;
    if(nodes[i]->Id() == ((filele->Nodes())[1])->Id() )
      nodenumber = 1;

    //compute triads from quaterions and unit direction vector of crosslinker
    quaterniontotriad(filele->Qnew_[nodenumber],Tfil);

    //auxiliary variable
    double scalarproduct = (direction(0)*Tfil(0,0) + direction(1)*Tfil(1,0) + direction(2)*Tfil(2,0));

    //computing Deltaphi and Delta at binding site i
    Deltaphi[i] = acos(scalarproduct) - statmechparams_.get<double>("PHI0",0.0);

    //substract from filament direction vector the component parallel to the crosslinker axis and normalize the result
    (filprojected[i])(0) = Tfil(0,0) - scalarproduct*direction(0);
    (filprojected[i])(1) = Tfil(1,0) - scalarproduct*direction(1);
    (filprojected[i])(2) = Tfil(2,0) - scalarproduct*direction(2);
    filprojected[i].Scale(1.0 / filprojected[i].Norm2() );
  }

  //tau is the angle between the projections of the two filament direction vectors into the plane perpendicular to the crosslinker axis
  Deltatau = acos((filprojected[0])(0)*(filprojected[1])(0) + (filprojected[0])(1)*(filprojected[1])(1) + (filprojected[0])(2)*(filprojected[1])(2)) - statmechparams_.get<double>("TAU0",0.0);

  //assuming bending and torsion potentials 0.5*EI*Deltaphi^2, 0.5*GJ*Deltatau^2 and a Boltzmann distribution for the different states of the crosslinker we get
  double p0 = exp(-0.5*EI*Deltaphi[0]*Deltaphi[0] / statmechparams_.get<double>("KT",0.0) ) ;
  double p1 = exp(-0.5*EI*Deltaphi[1]*Deltaphi[1] / statmechparams_.get<double>("KT",0.0) ) ;
  double p2 = exp(-0.5*GJ*Deltatau*Deltatau       / statmechparams_.get<double>("KT",0.0) ) ;

  //crosslinker has to pass three probability checks with respect to orientation
  return(UniformGen.random() < p0 && UniformGen.random() < p1 && UniformGen.random() < p2);

#endif  // #ifdef D_BEAM3
#endif  // #ifdef D_BEAM3II

} // StatMechManager::Permutation


/*----------------------------------------------------------------------*
 | output for structural polymorphism             (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void StatMechManager::StructPolymorphOutput(const Epetra_Vector& disrow, const std::ostringstream& filename)
{
	/* The following code is basically copied from GmshOutput() and therefore remains largely uncommented.
	 * The output consists of a vector v E R³ for each element and the respective node coordinates C.
	 * Since the output is redundant in the case of ghosted elements (parallel use), the GID of the
	 * element is written, too. This ensures an easy adaption of the output for post-processing.
	 *
	 * columns of the output file: element-ID, filament-ID, v_x, v_y, v_z, C_0_x,  C_0_y, C_0_z, C_1_x, C_1_y, C_1_z
	 * */

	FILE* fp = NULL;

	fp = fopen(filename.str().c_str(), "w");

	std::stringstream fiberorientfilecontent;

	for (int i=0; i<discret_.NumMyRowElements(); ++i)
	{
		DRT::Element* element = discret_.lRowElement(i);

		// we only need filament elements, crosslinker elements are neglected with the help of their high element IDs
		if(element->Id()>basisnodes_)
			continue;

		// get filament number conditions
		vector<DRT::Condition*> filamentnumberconditions(0);
		discret_.GetCondition("FilamentNumber",filamentnumberconditions);

		LINALG::SerialDenseMatrix coord(3,element->NumNode());
		// filament number
		int filnumber = 0;
		// a switch to avoid redundant search for filament number
		bool done = false;

		for(int j = 0; j<3; j++)
		{
			for(int k = 0; k<element->NumNode(); k++)
			{
				double referenceposition = ((element->Nodes())[k])->X()[j];
				vector<int> dofnode = discret_.Dof((element->Nodes())[k]);
				double displacement = disrow[discret_.DofRowMap()->LID( dofnode[j] )];
				coord(j,k) =  referenceposition + displacement;

				// get corresponding filament number
				// loop over all filaments
				if(!done)
					for(int l=0; l<(int)filamentnumberconditions.size(); l++)
						for(int m=0; m<(int)filamentnumberconditions[l]->Nodes()->size(); m++)
						 if(element->Nodes()[k]->Id() == filamentnumberconditions[l]->Nodes()->at(m))
						 {
							 filnumber = filamentnumberconditions[l]->Id();
							 done = true;
						 }
			}
		}

		const DRT::ElementType & eot = element->ElementType();

#ifdef D_BEAM3
		if(eot==DRT::ELEMENTS::Beam3Type::Instance())
			for(int j=0; j<element->NumNode()-1; j++)
				fiberorientfilecontent<<element->Id()<<" "<<filnumber<<"   "<<coord(0,j+1)-coord(0,j) << " " <<coord(1,j+1)-coord(1,j) << " " <<coord(2,j+1)-coord(2,j)<<"   "
																										<< coord(0,j) << " " << coord(1,j) << " " << coord(2,j) << " "
																										<< coord(0,j+1) << " " << coord(1,j+1) << " " << coord(2,j+1)<<endl ;
		else
#endif
#ifdef D_BEAM3II
		if(eot==DRT::ELEMENTS::Beam3iiType::Instance())
			for(int j=0; j<element->NumNode()-1; j++)
				fiberorientfilecontent<<element->Id()<<" "<<filnumber<<"   "<<coord(0,j+1)-coord(0,j) << " " <<coord(1,j+1)-coord(1,j) << " " <<coord(2,j+1)-coord(2,j)<<"   "
																										<< coord(0,j) << " " << coord(1,j) << " " << coord(2,j) << " "
																										<< coord(0,j+1) << " " << coord(1,j+1) << " " << coord(2,j+1)<<endl ;
		else
#endif
#ifdef D_TRUSS3
		if(eot==DRT::ELEMENTS::Truss3Type::Instance())
			for(int j=0; j<element->NumNode()-1; j++)
				fiberorientfilecontent<<element->Id()<<" "<<filnumber<<"   "<<coord(0,j+1)-coord(0,j) << " " <<coord(1,j+1)-coord(1,j) << " " <<coord(2,j+1)-coord(2,j)<<"   "
																										<< coord(0,j) << " " << coord(1,j) << " " << coord(2,j) << " "
																										<< coord(0,j+1) << " " << coord(1,j+1) << " " << coord(2,j+1)<<endl ;
		else
#endif
		{
		}
	}
  //write content into file and close it
  fprintf(fp,fiberorientfilecontent.str().c_str());
  fclose(fp);
}// StatMechManager:StructPolymorphOutput()

/*----------------------------------------------------------------------*
 | check the binding mode of a crosslinker        (public) mueller 07/10|
 *----------------------------------------------------------------------*/
bool StatMechManager::CheckForKinkedVisual(const double& color, int eleid)
{
	// if element is of a filament
	if(color != 0.5)
		return false;
	// else, if element is a crosslinker
	else
	{
		DRT::Element *element = discret_.gElement(eleid);

		std::vector<int> filamentids;
		std::vector<DRT::Condition*> filaments(0);
		discret_.GetCondition("FilamentNumber", filaments);
		// determine filament IDs of the nodes belonging to the crosslinker
		for(int i=0; i<(int)filaments.size(); i++)
		{
			// exit loop if a match was found for all crosslinker element nodes
			if((int)filamentids.size()==element->NumNode())
				break;
			else
				for(int j=0; j<(int)filaments[i]->Nodes()->size(); j++)
					for(int k=0; k<element->NumNode(); k++)
						if(element->Nodes()[k]->Id() == filaments[i]->Nodes()->at(j))
							filamentids.push_back(filaments[i]->Id());
		}
		/*cout<<"filament IDs: ";
		for(int i=0; i<(int)filamentids.size(); i++)
			cout<<filamentids.at(i)<<" ";
		cout<<" "<<endl;*/
		// check affiliation of (first and last) crosslinker nodes to filament(s)
		if(filamentids.at(0)==filamentids.at(element->NumNode()-1))
			return true;
		else
			return false;
	}
}

/*----------------------------------------------------------------------*
 | Special Gmsh output for crosslinkers occupying two binding spots on  |
 | the same filament											        (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void StatMechManager::GmshKinkedVisual(const LINALG::SerialDenseMatrix& coord, std::vector<double>& thirdpoint, int eleid)
{
	/* We need a third point in order to visualize the crosslinker.
	 * It marks the location of the kink.
	 */
	// get the element
	DRT::Element* element = discret_.gElement(eleid);

	// calculate tangent
	double ltan = 0.0;
	for(int j=0; j<coord.M(); j++)
		ltan += (coord(j,coord.N()-1)-coord(j,0))*(coord(j,coord.N()-1)-coord(j,0));
	ltan = sqrt(ltan);

	std::vector<double> t(3,0.0);
	for(int j=0; j<coord.M(); j++)
		t.at(j) = ((coord(j,coord.N()-1)-coord(j,0))/ltan);

	// calculate normal via cross product: [0 0 1]x[tx ty tz]
	std::vector<double> n(3,0.0);
	n.at(0) = -t.at(1);
	n.at(1) = t.at(0);
	// norm it since the cross product does not keep the length
	double lnorm = 0.0;
	for(int j=0; j<(int)n.size(); j++)
		lnorm += n.at(j)*n.at(j);
	lnorm = sqrt(lnorm);
	for(int j=0; j<(int)n.size(); j++)
		n.at(j) /= lnorm;


	// obtain angle
	// random angle
	//ranlib::UniformClosed<double> UniformGen;
	//double alpha = 2*M_PI*UniformGen.random();
	// by modulo operation involving the node IDs
	double alpha = fmod( (double)(element->Nodes()[element->NumNode()-1]->Id() + element->Nodes()[0]->Id()), 2*M_PI);

	// rotate the normal by alpha
	LINALG::SerialDenseMatrix RotMat(3,3);
	// build the matrix of rotation
	for(int j=0; j<3; j++)
		RotMat(j,j) = cos(alpha)+t.at(j)*t.at(j)*(1-cos(alpha));
	RotMat(0,1) = t.at(0)*t.at(1)*(1-cos(alpha))-t.at(2)*sin(alpha);
	RotMat(0,2) = t.at(0)*t.at(2)*(1-cos(alpha))+t.at(1)*sin(alpha);
	RotMat(1,0) = t.at(1)*t.at(0)*(1-cos(alpha))+t.at(2)*sin(alpha);
	RotMat(1,2) = t.at(1)*t.at(2)*(1-cos(alpha))-t.at(0)*sin(alpha);
	RotMat(2,0) = t.at(2)*t.at(0)*(1-cos(alpha))-t.at(1)*sin(alpha);
	RotMat(2,1) = t.at(2)*t.at(1)*(1-cos(alpha))+t.at(0)*sin(alpha);

	// rotation
	std::vector<double> nrot(3,0.0);
	for(int j=0; j<3; j++)
		for(int k=0; k<3; k++)
			nrot.at(j) += RotMat(j,k)*n.at(k);

	// calculation of the third point lying in the direction of the rotated normal
	// height of the third point above the filament
	double h = 0.33*(statmechparams_.get<double>("R_LINK",0.0)+statmechparams_.get<double>("DeltaR_LINK",0.0));
	for(int j=0; j<3; j++)
		thirdpoint.at(j) = (coord(j,0) + coord(j,element->NumNode()-1))/2.0 + h*nrot.at(j);

	/*/ cout block
	cout<<"coord  = \n"<<coord<<endl;
	cout<<"ltan   = "<<ltan<<endl;
	cout<<"t      = [ "<<t.at(0)<<" "<<t.at(1)<<" "<<t.at(2)<<" ]"<<endl;
	cout<<"lnorm  = "<<lnorm<<endl;
	cout<<"n      = [ "<<n.at(0)<<" "<<n.at(1)<<" "<<n.at(2)<<" ]"<<endl;
	cout<<"alpha  = "<<alpha<<endl;
	cout<<"RotMat =\n"<<RotMat<<endl;
	cout<<"nrot   = [ "<<nrot.at(0)<<" "<<nrot.at(1)<<" "<<nrot.at(2)<<" ]"<<endl;
	cout<<"thirdp = [ "<<thirdpoint.at(0)<<" "<<thirdpoint.at(1)<<" "<<thirdpoint.at(2)<<" ]\n\n\n"<<endl;*/
}
#endif  // #ifdef CCADISCRET
