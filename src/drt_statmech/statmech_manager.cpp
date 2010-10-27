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
isinit_(false),
konswitch_(false),
nsearch_(0),
starttimeoutput_(-1.0),
endtoendref_(0.0),
istart_(0),
basisnodes_(discret.NumGlobalNodes()),
basiselements_(discret.NumGlobalElements()),
currentelements_(discret.NumGlobalElements()),
outputfilenumber_(-1),
normalgen_(0,1),
discret_(discret)
{
  //if input flag FIXEDSEED == YES: use same random numbers in each program start based on input variable INITIALSEED
	int seedvariable = 0;
  if(Teuchos::getIntegralValue<int>(statmechparams_,"FIXEDSEED"))
  	seedvariable = statmechparams_.get<int>("INITIALSEED", 0);
  //else set seed according to system time and different for each processor
  else
  	seedvariable = time(0)*(discret_.Comm().MyPID() + 1);

	normalgen_.seed( (unsigned int)seedvariable );
	uniformclosedgen_.seed( (unsigned int)seedvariable );
	uniformclosedopengen_.seed( (unsigned int)seedvariable );


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
			for(int j=0; j<(int)nodeids->size() ; j++)
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
			for (int j=0; j<(int)nodeids->size(); j++)
			{
				//get the node's global id
				int nodenumber = (*nodeids)[j];

				//testing whether current nodedofnumber makes sense for current node
				if (nodedofnumber < 0 || nodedofnumber >= discret_.NumDof(discret_.gNode(nodenumber)))
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
			currpos(2) = node->X()[2];
			currentpositions[node->LID()] = currpos;
		}

		//find bounding box for search tree
		const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(discret_, currentpositions);

		//initialize search tree
		octTree_->initializePointTree(rootBox,currentpositions,GEO::TreeType(GEO::OCTTREE));

	/* Initialization of N_CROSSLINK crosslinker molecule REPRESENTATIONS. As long as the molecules do not act as a link
	 * between two filaments, only their positions are calculated. Here, the molecules' initial positions are determined.
	 * Calculations are made on Proc 0, only.*/
	if(Teuchos::getIntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(),"CRSLNKDIFFUSION"))
		CrosslinkerMoleculeInit();

	}//if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))
	return;
}// StatMechManager::StatMechManager

/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechOutput(ParameterList& params, const int ndim,
																		 const double& time, const int& istep, const double& dt,
																		 const Epetra_Vector& dis, const Epetra_Vector& fint)
{
	/*in general simulations in statistical mechanics run over so many time steps that the amount of data stored in the error file
	 * may exceed the capacity even of a server hard disk; thus, we rewind the error file in each time step so that the amount of data
	 * does not increase after the first time step any longer*/
	bool printerr = params.get<bool> ("print to err", false);
	FILE* errfile = params.get<FILE*> ("err file", NULL);
	if (printerr)
		rewind(errfile);

	//the following variable makes sense in case of serial computing only; its use is not allowed for parallel computing!
	int num_dof = dis.GlobalLength();

  double starttime = statmechparams_.get<double>("STARTTIME",0.0);

	switch (Teuchos::getIntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT"))
	{
		case INPAR::STATMECH::statout_endtoendlog:
		{

      double endtoend = 0; //end to end length at a certain time step in single filament dynamics
      double DeltaR2 = 0;

      //as soon as system is equilibrated (after time STARTTIME) a new file for storing output is generated
      if ( (time > starttime && fabs(time-starttime)>dt/1e4) && (starttimeoutput_ == -1.0))
      {
        endtoendref_ = pow(pow((dis)[num_dof - 3] + 10 - (dis)[0], 2) + pow(
            (dis)[num_dof - 2] - (dis)[1], 2), 0.5);
        starttimeoutput_ = time;
        istart_ = istep;
      }
      if (time > starttimeoutput_ && starttimeoutput_ > -1.0)
      {
        endtoend = pow(pow((dis)[num_dof - 3] + 10 - (dis)[0], 2) + pow(
            (dis)[num_dof - 2] - (dis)[1], 2), 0.5);

        //applying in the following a well conditioned substraction formula according to Crisfield, Vol. 1, equ. (7.53)
        DeltaR2 = pow((endtoend * endtoend - endtoendref_ * endtoendref_) / (endtoend + endtoendref_), 2);

        //writing output: writing Delta(R^2) according to PhD thesis Hallatschek, eq. (4.60), where t=0 corresponds to starttimeoutput_
        if ((istep - istart_) % int(ceil(pow(10, floor(log10((time - starttimeoutput_) / (10* dt )))))) == 0)
        {

          //proc 0 write complete output into file, all other proc inactive
          if(!discret_.Comm().MyPID())
          {
            FILE* fp = NULL; //file pointer for statistical output file

            //name of output file
            std::ostringstream outputfilename;
            outputfilename << "EndToEnd" << outputfilenumber_ << ".dat";

            // open file and append new data line
            fp = fopen(outputfilename.str().c_str(), "a");

            //defining temporary stringstream variable
            std::stringstream filecontent;
            filecontent << scientific << setprecision(15) << time - starttimeoutput_ << "  " << DeltaR2 << endl;

            // move temporary stringstream to file and close it
            fprintf(fp, filecontent.str().c_str());
            fclose(fp);
          }
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
			discret_.GetCondition("PointNeumann", pointneumannconditions);
			if (pointneumannconditions.size() > 0)
			{
				const vector<double>* val = pointneumannconditions[0]->Get<vector<double> > ("val");
				neumannforce = fabs((*val)[0]);
			}
			else
				neumannforce = 0;

			double endtoend = 0.0; //end to end length at a certain time step in single filament dynamics

			//as soon as system is equilibrated (after time STARTTIME) a new file for storing output is generated
			if ((time > starttime && fabs(time-starttime)>dt/1e4) && (starttimeoutput_ == -1.0))
			{
				starttimeoutput_ = time;
				istart_ = istep;
			}

			if (time > starttimeoutput_ && starttimeoutput_ > -1.0)
			{
				//end to end vector
				LINALG::Matrix<3, 1> endtoendvector(true);
				for (int i = 0; i < ndim; i++)
				{
					endtoendvector(i) -= (discret_.gNode(0))->X()[i] + dis[i];
					endtoendvector(i) += (discret_.gNode(discret_.NumMyRowNodes() - 1))->X()[i] + dis[num_dof - discret_.NumDof(discret_.gNode(discret_.NumMyRowNodes() - 1)) + i];
				}

				endtoend = endtoendvector.Norm2();

				//writing output: current time and end to end distance are stored at each 100th time step
				if ((istep - istart_) % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0)
				{

		      //proc 0 write complete output into file, all other proc inactive
		      if(!discret_.Comm().MyPID())
		      {
            FILE* fp = NULL; //file pointer for statistical output file

            //name of output file
            std::ostringstream outputfilename;
            outputfilename << "E2E_" << discret_.NumMyRowElements() << '_' << dt<< '_' << neumannforce << '_' << outputfilenumber_ << ".dat";

            // open file and append new data line
            fp = fopen(outputfilename.str().c_str(), "a");
            //defining temporary stringstream variable
            std::stringstream filecontent;
            filecontent << scientific << setprecision(15) << time << "  "
                        << endtoend << " " << fint[num_dof - discret_.NumDof(
                           discret_.gNode(discret_.NumMyRowNodes() - 1))] << endl;
            // move temporary stringstream to file and close it
            fprintf(fp, filecontent.str().c_str());
            fclose(fp);
		      }
				}
			}
		}
		break;
		/*computing and writing into file data about correlation of orientation of different elements
		 *as it is considered in the context of the persistence length*/
		case INPAR::STATMECH::statout_orientationcorrelation:
		{

		  //we need displacements also of ghost nodes and hence export displacment vector to column map format
		  Epetra_Vector discol(*(discret_.DofColMap()), true);
		  LINALG::Export(dis, discol);

			vector<double> arclength(discret_.NumMyColNodes() - 1, 0);
			vector<double> cosdiffer(discret_.NumMyColNodes() - 1, 0);

			//after initilization time write output cosdiffer in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps,
			//when discret_.NumMyRowNodes()-1 = 0,cosdiffer is always equil to 1!!
			if ((time > starttime && fabs(time-starttime)>dt/1e4) && (istep% statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0))
			{
				Epetra_SerialDenseMatrix coord;
				coord.Shape(discret_.NumMyColNodes(), ndim);

				for (int id=0; id<discret_.NumMyColNodes(); id++)
					for (int j=0; j<ndim; j++)
						coord(id, j) = (discret_.lColNode(id))->X()[j] + (discol)[(id) * (ndim - 1) * 3 + j];

				for (int id=0; id < discret_.NumMyColNodes() - 1; id++)
				{

					//calculate the deformed length of every element
					for (int j=0; j<ndim; j++)
					{
						arclength[id] += pow((coord(id + 1, j) - coord(id, j)), 2);
						cosdiffer[id] += (coord(id + 1, j) - coord(id, j)) * (coord(0 + 1, j) - coord(0, j));
					}

					//calculate the cosine difference referring to the first element
					//Dot product of the (id+1)th element with the 1st element and devided by the length of the (id+1)th element and the 1st element
					arclength[id] = pow(arclength[id], 0.5);
					cosdiffer[id] = cosdiffer[id] / (arclength[id] * arclength[0]);
				}

				//proc 0 write complete output into file, all other proc inactive
				if(!discret_.Comm().MyPID())
				{
		      FILE* fp = NULL; //file pointer for statistical output file

		      //name of output file
		      std::ostringstream outputfilename;
		      outputfilename.str("");
		      outputfilename << "OrientationCorrelation" << outputfilenumber_ << ".dat";

          fp = fopen(outputfilename.str().c_str(), "a");
          std::stringstream filecontent;
          filecontent << istep;
          filecontent << scientific << setprecision(10);

          for (int id = 0; id < discret_.NumMyColNodes() - 1; id++)
            filecontent << " " << cosdiffer[id];

          filecontent << endl;
          fprintf(fp, filecontent.str().c_str());
          fclose(fp);
				}

			}

		}
		break;
		//the following output allows for anisotropic diffusion simulation of a quasi stiff polymer
		case INPAR::STATMECH::statout_anisotropic:
		{
			//output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
			if ((time > starttime && fabs(time-starttime)>dt/1e4) && (istep % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0))
			{


				//positions of first and last node in current time step (note: always 3D variables used; in case of 2D third compoenent just constantly zero)
				LINALG::Matrix<3, 1> beginnew;
				LINALG::Matrix<3, 1> endnew;

				beginnew.PutScalar(0);
				endnew.PutScalar(0);
				std::cout << "ndim: " << ndim << "\n";

				for (int i = 0; i < ndim; i++)
				{
					beginnew(i) = (discret_.gNode(0))->X()[i] + dis[i];
					endnew(i) = (discret_.gNode(discret_.NumMyRowNodes() - 1))->X()[i] + dis[num_dof - discret_.NumDof(discret_.gNode(discret_.NumMyRowNodes() - 1)) + i];
				}

				//unit direction vector for filament axis in last time step
				LINALG::Matrix<3, 1> axisold;
				axisold = endold_;
				axisold -= beginold_;
				axisold.Scale(1 / axisold.Norm2());

				//unit direction vector for filament axis in current time step
				LINALG::Matrix<3, 1> axisnew;
				axisnew = endnew;
				axisnew -= beginnew;
				axisnew.Scale(1 / axisnew.Norm2());

				//displacement of first and last node between last time step and current time step
				LINALG::Matrix<3, 1> dispbegin;
				LINALG::Matrix<3, 1> dispend;
				dispbegin = beginnew;
				dispbegin -= beginold_;
				dispend = endnew;
				dispend -= endold_;

				//displacement of middle point
				LINALG::Matrix<3, 1> dispmiddle;
				dispmiddle = dispbegin;
				dispmiddle += dispend;
				dispmiddle.Scale(0.5);
				sumdispmiddle_ += dispmiddle;

				//update sum of square displacement increments of middle point
				double incdispmiddle = dispmiddle.Norm2() * dispmiddle.Norm2();
				sumsquareincmid_ += incdispmiddle;

				//update sum of square displacement increments of middle point parallel to new filament axis (computed by scalar product)
				double disppar_square = pow(axisnew(0) * dispmiddle(0) + axisnew(1) * dispmiddle(1) + axisnew(2) * dispmiddle(2), 2);
				sumsquareincpar_ += disppar_square;

				//update sum of square displacement increments of middle point orthogonal to new filament axis (via crossproduct)
				LINALG::Matrix<3, 1> aux;
				aux(0) = dispmiddle(1) * axisnew(2) - dispmiddle(2) * axisnew(1);
				aux(1) = dispmiddle(2) * axisnew(0) - dispmiddle(0) * axisnew(2);
				aux(2) = dispmiddle(0) * axisnew(1) - dispmiddle(1) * axisnew(0);
				double disport_square = aux.Norm2() * aux.Norm2();
				sumsquareincort_ += disport_square;

				//total displacement of rotational angle (in 2D only)
				double incangle = 0;
				if (ndim == 2)
				{
					//angle of old axis relative to x-axis
					double phiold = acos(axisold(0) / axisold.Norm2());
					if (axisold(1) < 0)
						phiold *= -1;

					//angle of new axis relative to x-axis
					double phinew = acos(axisnew(0) / axisnew.Norm2());
					if (axisnew(1) < 0)
						phinew *= -1;

					//angle increment
					incangle = phinew - phiold;
					if (incangle > M_PI)
					{
						incangle -= 2*M_PI;
						incangle *= -1;
					}
					if (incangle < -M_PI)
					{
						incangle += 2*M_PI;
						incangle *= -1;
					}

					//update absolute rotational displacement compared to reference configuration
					sumsquareincrot_ += incangle * incangle;
					sumrotmiddle_ += incangle;
				}

        //proc 0 write complete output into file, all other proc inactive
        if(!discret_.Comm().MyPID())
				{
	        FILE* fp = NULL; //file pointer for statistical output file

	        //name of output file
	        std::ostringstream outputfilename;
	        outputfilename << "AnisotropicDiffusion" << outputfilenumber_ << ".dat";

					// open file and append new data line
					fp = fopen(outputfilename.str().c_str(), "a");

					//defining temporary stringstream variable
					std::stringstream filecontent;
					filecontent << scientific << setprecision(15) << dt << " "
							<< sumsquareincmid_ << " " << sumsquareincpar_ << " "
							<< sumsquareincort_ << " " << sumsquareincrot_ << " "
							<< sumdispmiddle_.Norm2() * sumdispmiddle_.Norm2() << " "
							<< sumrotmiddle_ * sumrotmiddle_ << endl;

					// move temporary stringstream to file and close it
					fprintf(fp, filecontent.str().c_str());
					fclose(fp);
				}

				//new positions in this time step become old positions in last time step
				beginold_ = beginnew;
				endold_ = endnew;
			}
		}
		break;
		case INPAR::STATMECH::statout_viscoelasticity:
		{
			//output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
			if (istep % statmechparams_.get<int> ("OUTPUTINTERVALS", 1) == 0)
			{
				//pointer to file into which each processor writes the output related with the dof of which it is the row map owner
				FILE* fp = NULL;

				//content to be written into the output file
				std::stringstream filecontent;

				//name of file into which output is written
				std::ostringstream outputfilename;
				outputfilename << "ViscoElOutputProc" << discret_.Comm().MyPID()<< ".dat";

				fp = fopen(outputfilename.str().c_str(), "a");

				/*the output to be written consists of internal forces at exactly those degrees of freedom
				 * marked in *forcesensor_ by a one entry*/

				filecontent << scientific << setprecision(10) << time;//changed

#ifdef DEBUG
				if (forcesensor_ == null)
					dserror("forcesensor_ is NULL pointer; possible reason: dynamic crosslinkers not activated and forcesensor applicable in this case only");
#endif  // #ifdef DEBUG
				double f = 0;//mean value of force
				double d = 0;//Displacement
				int count = 0;

				for(int i=0; i<forcesensor_->MyLength(); i++)//changed

				{
					if((*forcesensor_)[i]>0.9)
					{
						count++;
						f += fint[i];
						d = dis[i];
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
			if( istep % statmechparams_.get<int>("OUTPUTINTERVALS",1) == 0 )
			{
				/*construct unique filename for gmsh output with two indices: the first one marking the time step number
				 * and the second one marking the newton iteration number, where numbers are written with zeros in the front
				 * e.g. number one is written as 000001, number fourteen as 000014 and so on;*/

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
			if( istep % statmechparams_.get<int>("OUTPUTINTERVALS",1) == 0 )
			{
				std::ostringstream filename1;
				std::ostringstream filename2;

				filename1 << "./GmshOutput/StruPolyCoords_Proc"<<discret_.Comm().MyPID()<<"_"<<std::setw(6) << setfill('0') << istep <<".dat";
				filename2 <<"./GmshOutput/StruPolyCrslnks_Proc"<<discret_.Comm().MyPID()<<"_"<<std::setw(6) << setfill('0') << istep <<".dat";
				StructPolymorphOutput(dis,filename1,filename2);
			}
		}
		break;
		case INPAR::STATMECH::statout_densitydensitycorr:
		{
			//output in every statmechparams_.get<int>("OUTPUTINTERVALS",1) timesteps
			if( istep % statmechparams_.get<int>("OUTPUTINTERVALS",1) == 0 )
			{

				std::ostringstream filename;
				filename << "./DensityDensityCorrFunction_"<<std::setw(6) << setfill('0') << istep <<".dat";
				DensityDensityCorrOutput(filename);
			}
		}
		break;
		case INPAR::STATMECH::statout_none:
		default:
		break;
	}

	return;
}	// StatMechManager::StatMechOutput()


/*----------------------------------------------------------------------*
 | writing Gmsh data for current step                 public)cyron 01/09|
 *----------------------------------------------------------------------*/
void StatMechManager::GmshOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step)
{
	/*the following method writes output data for Gmsh into file with name "filename"; all line elements are written;
	 * the nodal displacements are handed over in the variable "dis"; note: in case of parallel computing only
	 * processor 0 writes; it is assumed to have a fully overlapping column map and hence all the information about
	 * all the nodal position; parallel output is now possible with the restriction that the nodes(processors) in question
	 * are of the same machine*/

	//we need displacements also of ghost nodes and hence export displacment vector to column map format
	Epetra_Vector discol(*(discret_.DofColMap()), true);
	LINALG::Export(disrow, discol);

	// do output to file in c-style
	FILE* fp = NULL;

	// first processor starts by opening the file and writing the header, other processors have to wait
	if (discret_.Comm().MyPID() == 0)
	{
		//open file to write output data into
		fp = fopen(filename.str().c_str(), "w");
		// write output to temporary stringstream;
		std::stringstream gmshfileheader;
		/*the beginning of the stream is "View \"" to indicate Gmsh that the following data is in order to create an image and
		 * this command is followed by the name of that view displayed during it's shown in the video; in the following example
		 * this name is for the 100th time step: Step00100; then the data to be presented within this view is written within { ... };
		 * in the following example this data consists of scalar lines defined by the coordinates of their end points*/
		gmshfileheader << "View \" Step " << step << " \" {" << endl;

		//write content into file and close it
		fprintf(fp, gmshfileheader.str().c_str());
		fclose(fp);
	}

	// wait for all processors to arrive at this point
	discret_.Comm().Barrier();

	// loop over the participating processors each of which appends its part of the output to one output file
	for (int proc = 0; proc < discret_.Comm().NumProc(); proc++)
	{
		if (discret_.Comm().MyPID() == proc)
		{
			//open file again to append ("a") output data into
			fp = fopen(filename.str().c_str(), "a");
			// write output to temporary stringstream;
			std::stringstream gmshfilecontent;

			//looping through all elements on the processor
			for (int i=0; i<discret_.NumMyColElements(); ++i)
			{
				//getting pointer to current element
				DRT::Element* element = discret_.lColElement(i);

				//getting number of nodes of current element
				//if( element->NumNode() > 2)
				//dserror("Gmsh output for two noded elements only");

				//preparing variable storing coordinates of all these nodes
				LINALG::SerialDenseMatrix coord(3, element->NumNode());
				for (int id=0; id<3; id++)
					for (int jd=0; jd<element->NumNode(); jd++)
					{
						double referenceposition = ((element->Nodes())[jd])->X()[id];
						vector<int> dofnode = discret_.Dof((element->Nodes())[jd]);
						double displacement =	discol[discret_.DofColMap()->LID(dofnode[id])];
						coord(id, jd) = referenceposition + displacement;
					}

				//declaring variable for color of elements
				double color;

				//apply different colors for elements representing filaments and those representing dynamics crosslinkers
				if (element->Id() < basisnodes_)
					color = 1.0;
				else
					color = 0.5;

				//if no periodic boundary conditions are to be applied, we just plot the current element
				if (statmechparams_.get<double> ("PeriodLength", 0.0) == 0)
				{
					// check whether the kinked visualization is to be applied
					bool kinked = CheckForKinkedVisual(element->Id());

					const DRT::ElementType & eot = element->ElementType();
#ifdef D_BEAM3
					if (eot == DRT::ELEMENTS::Beam3Type::Instance())
					{
						if (!kinked)
						{
							for (int j=0; j<element->NumNode() - 1; j++)
							{
								//writing element by nodal coordinates as a scalar line
								gmshfilecontent << "SL(" << scientific;
								gmshfilecontent << coord(0, j) << "," << coord(1, j) << ","<< coord(2, j) << "," << coord(0, j + 1) << "," << coord(1,j + 1) << "," << coord(2, j + 1);
								/*note: colors are chosen by values between 0.0 and 1.0. These values refer to a color vector which
								 * is given in a .geo-setup-file. If, for example, 5 colors are given(either in X11 color expressions or RGB),
								 * possible values are 0.0, 0.25, 0.5, 0.75, 1.0.
								 */
								gmshfilecontent << ")" << "{" << scientific << color << ","<< color << "};" << endl;
							}
						}
						else
							GmshKinkedVisual(coord, 0.375, element->Id(), gmshfilecontent);
					}
					else
#endif
#ifdef D_BEAM3II
					if(eot==DRT::ELEMENTS::Beam3iiType::Instance())
					{
						if(!kinked)
						{
							for(int j=0; j<element->NumNode()-1; j++)
							{
								gmshfilecontent << "SL(" << scientific;
								gmshfilecontent<< coord(0,j) << "," << coord(1,j) << "," << coord(2,j) << ","
															 << coord(0,j+1) << "," << coord(1,j+1) << "," << coord(2,j+1);
								gmshfilecontent<< ")" << "{" << scientific << color << "," << color << "};" << endl;
							}
						}
						else
							GmshKinkedVisual(coord, 0.375, element->Id(), gmshfilecontent);
					}
					else
#endif
#ifdef D_TRUSS3
					if (eot == DRT::ELEMENTS::Truss3Type::Instance())
					{
						if (!kinked)
						{
							for (int j=0; j<element->NumNode()-1; j++)
							{
								gmshfilecontent << "SL(" << scientific;
								gmshfilecontent << coord(0, j) << "," << coord(1, j) << ","<< coord(2, j) << ","
																<< coord(0, j + 1) << "," << coord(1,j + 1) << "," << coord(2, j + 1);
								gmshfilecontent << ")" << "{" << scientific << color << ","<< color << "};" << endl;
							}
						}
						else
							GmshKinkedVisual(coord, 0.375, element->Id(), gmshfilecontent);
					}
					else
#endif
#ifdef D_TORSION3
					if (eot == DRT::ELEMENTS::Torsion3Type::Instance())
					{
						double beadcolor = 0.75;
						for (int j=0; j<element->NumNode(); j++)
						{
							gmshfilecontent << "SP(" << scientific;
							gmshfilecontent << coord(0, j) << "," << coord(1, j) << ","<< coord(2, j);
							gmshfilecontent << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
						}
					}
					else
#endif
					{
					}
				}
				//in case of periodic boundary conditions we have to take care to plot correctly an element broken at some boundary plane
				else
					GmshOutputPeriodicBoundary(coord, color, gmshfilecontent,element->Id());
			}
			//write content into file and close it (this way we make sure that the output is written serially)
			fprintf(fp, gmshfilecontent.str().c_str());
			fclose(fp);
		}
		discret_.Comm().Barrier();
	}
	// plot the periodic boundary box
	GmshOutputBoundaryBox(0.0, &filename);
	// plot crosslink molecule diffusion and (partial) bonding
	GmshOutputCrosslinkDiffusion(0.125, &filename, disrow);
	// finish data section of this view by closing curly brackets
	if (discret_.Comm().MyPID() == 0)
	{
		fp = fopen(filename.str().c_str(), "a");
		std::stringstream gmshfileend;
		gmshfileend << "};" << endl;
		fprintf(fp, gmshfileend.str().c_str());
		fclose(fp);
	}

	// return simultaneously (not sure if really needed)
	discret_.Comm().Barrier();

	return;
} // StatMechManager::GmshOutput()


/*----------------------------------------------------------------------*
 | gmsh output data in case of periodic boundary conditions             |
 |                                                    public)cyron 02/10|
 *----------------------------------------------------------------------*/
void StatMechManager::GmshOutputPeriodicBoundary(const LINALG::SerialDenseMatrix& coord, const double& color, std::stringstream& gmshfilecontent, int eleid, bool ignoreeleid)
{
	//number of spatial dimensions
	const int ndim = 3;
	// get Element Type of the first Element to determine the graphics output
	DRT::Element* element = discret_.gElement(eleid);

	bool dotline = false;
	bool kinked = false;

	if (ignoreeleid)
		dotline = true;
	else
	{
		// draw colored lines between two nodes of a beam3 or truss3 element (meant for filaments/crosslinks/springs)
		const DRT::ElementType & eot = element->ElementType();
#ifdef D_BEAM3
		if (element->ElementType().Name() == "Beam3Type")
			dotline = eot == DRT::ELEMENTS::Beam3Type::Instance();
#endif
#ifdef D_BEAM3II
		if(element->ElementType().Name()=="Beam3iiType")
		dotline = eot==DRT::ELEMENTS::Beam3iiType::Instance();
#endif
#ifdef D_TRUSS3
		if (element->ElementType().Name() == "Truss3Type")
			dotline = dotline or eot == DRT::ELEMENTS::Truss3Type::Instance();
#endif
#ifdef D_TORSION3
		// draw spheres at node positions ("beads" of the bead spring model)
		if (eot == DRT::ELEMENTS::Torsion3Type::Instance())
		{
			double beadcolor = 0.75;
			for (int i=0; i<element->NumNode(); i++)
			{
				//writing element by nodal coordinates as a sphere
				gmshfilecontent << "SP(" << scientific;
				gmshfilecontent << coord(0,i) << "," << coord(1, i) << "," << coord(2,i);
				gmshfilecontent << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
			}
		}
#endif
		/* determine whether crosslink connects two filaments or occupies two binding spots on the same filament;
		 * variable triggers different visualizations.*/
		kinked = CheckForKinkedVisual(element->Id());
	}

	if (dotline)
	{
		/*detect and save in vector "cut", at which boundaries the element is broken due to periodic boundary conditions;
		 * the entries of cut have the following meaning: 0: element not broken in respective coordinate direction, 1:
		 * element broken in respective coordinate direction (node 0 close to zero boundary and node 1 close to boundary
		 * at PeriodLength);  2: element broken in respective coordinate direction (node 1 close to zero boundary and node
		 * 0 close to boundary at PeriodLength);*/
		LINALG::SerialDenseMatrix cut;
		if (ignoreeleid)
			cut = LINALG::SerialDenseMatrix(3, 1, true);
		else
			cut = LINALG::SerialDenseMatrix(3, (int)element->NumNode()-1, true);

		/* "coord" currently holds the shifted set of coordinates.
		 * In order to determine the correct vector "dir" of the visualization at the boundaries,
		 * a copy of "coord" with adjustments in the proper places is introduced*/
		LINALG::SerialDenseMatrix unshift = coord;

		for (int i=0; i<cut.N(); i++)
		{
			for (int dof=0; dof<ndim; dof++)
			{
				if (fabs(coord(dof,i+1)-statmechparams_.get<double>("PeriodLength",0.0)-coord(dof,i)) < fabs(coord(dof,i+1) - coord(dof,i)))
				{
					cut(dof, i) = 1;
					unshift(dof, i + 1) -= statmechparams_.get<double>("PeriodLength",0.0);
				}
				if (fabs(coord(dof, i+1)+statmechparams_.get<double>("PeriodLength",0.0) - coord(dof,i)) < fabs(coord(dof,i+1)-coord(dof,i)))
				{
					cut(dof,i) = 2;
					unshift(dof,i+1) += statmechparams_.get<double>("PeriodLength",0.0);
				}
			}
		}
		/*for(int i=0;i<coord.M(); i++)
		 {
		 for(int j=0; j<coord.N(); j++)
		 if(coord(i,j)-unshift(i,j)>0.0)
		 cout<<"shift++ : "<<coord(i,j) - unshift(i,j)<<endl;
		 else if(coord(i,j)-unshift(i,j)<0.0)
		 cout<<"shift-- : "<<coord(i,j) - unshift(i,j)<<endl;
		 }*/

		// write special output for broken elements
		for (int i=0; i<cut.N(); i++)
		{
			if (cut(0, i) + cut(1, i) + cut(2, i) > 0)
			{
				//compute direction vector between first(i-th) and second(i+1-th) node of element (normed):
				LINALG::Matrix<3, 1> dir;
				double ldir = 0.0;
				for (int dof = 0; dof < ndim; dof++)
				{
					dir(dof) = unshift(dof, i + 1) - unshift(dof, i);
					ldir += dir(dof) * dir(dof);
				}
				for (int dof = 0; dof < ndim; dof++)
					dir(dof) /= ldir;

				//from node 0 to nearest boundary where element is broken you get by vector X + lambda0*dir
				double lambda0 = dir.Norm2();
				for (int dof = 0; dof < ndim; dof++)
				{
					if (cut(dof, i) == 1)
					{
						if (fabs(-coord(dof, i) / dir(dof)) < fabs(lambda0))
							lambda0 = -coord(dof, i) / dir(dof);
					}
					else if (cut(dof, i) == 2)
					{
						if (fabs((statmechparams_.get<double> ("PeriodLength", 0.0) - coord(dof, i)) / dir(dof)) < fabs(lambda0))
							lambda0 = (statmechparams_.get<double> ("PeriodLength", 0.0) - coord(dof, i)) / dir(dof);
					}
				}

				//from node 1 to nearest boundary where element is broken you get by vector X + lambda1*dir
				double lambda1 = dir.Norm2();
				for (int dof = 0; dof < ndim; dof++)
				{
					if (cut(dof, i) == 2)
					{
						if (fabs(-coord(dof, i + 1) / dir(dof)) < fabs(lambda1))
							lambda1 = -coord(dof, i + 1) / dir(dof);
					}
					else if (cut(dof, i) == 1)
					{
						if (fabs((statmechparams_.get<double> ("PeriodLength", 0.0) - coord(dof, i + 1)) / dir(dof)) < fabs(lambda1))
							lambda1 = (statmechparams_.get<double> ("PeriodLength", 0.0) - coord(dof, i + 1)) / dir(dof);
					}
				}

				//writing element by nodal coordinates as a scalar line
				gmshfilecontent << "SL(" << scientific;
				gmshfilecontent << coord(0, i) << "," << coord(1, i) << "," << coord(2,i) << ","
												<< coord(0, i) + lambda0 * dir(0) << ","<< coord(1, i)+ lambda0 * dir(1) << "," << coord(2, i) + lambda0 * dir(2);
				/*note: for each node there is one color variable for gmsh and gmsh finally plots the line
				 * interpolating these two colors between the nodes*/
				gmshfilecontent << ")" << "{" << scientific << color << "," << color<< "};" << endl;
				//writing element by nodal coordinates as a scalar line
				gmshfilecontent << "SL(" << scientific;
				gmshfilecontent << coord(0, i + 1) << "," << coord(1, i + 1) << "," << coord(2, i + 1) << ","
												<< coord(0, i + 1) + lambda1 * dir(0)<< "," << coord(1, i + 1) + lambda1 * dir(1) << "," << coord(2, i + 1) + lambda1 * dir(2);
				/*note: for each node there is one color variable for gmsh and gmsh finally plots the line
				 * interpolating these two colors between the nodes*/
				gmshfilecontent << ")" << "{" << scientific << color << "," << color << "};" << endl;

				// crosslink molecules with one bond
				if (ignoreeleid)
				{
					double beadcolor = color;
					gmshfilecontent << "SP(" << scientific;
					gmshfilecontent << coord(0, i + 1) << "," << coord(1, i + 1) << ","<< coord(2, i + 1);
					gmshfilecontent << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
				}
			}
			else // output for continuous elements
			{
				if (!kinked)
				{
					// filament or crosslink between two filaments
					//writing element by nodal coordinates as a scalar line
					gmshfilecontent << "SL(" << scientific;
					gmshfilecontent << coord(0, i) << "," << coord(1, i) << "," << coord(2, i) << ","
												  << coord(0, i + 1) << "," << coord(1, i + 1) << ","<< coord(2, i + 1);
					/*note: for each node there is one color variable for gmsh and gmsh finally plots the line
					 * interpolating these two colors between the nodes*/
					gmshfilecontent << ")" << "{" << scientific << color << "," << color<< "};" << endl;

					if (ignoreeleid)
					{
						double beadcolor = color; //blue
						gmshfilecontent << "SP(" << scientific;
						gmshfilecontent << coord(0, i + 1) << "," << coord(1, i + 1) << ","<< coord(2, i + 1);
						gmshfilecontent << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
					}
				}
				else
					GmshKinkedVisual(coord, 0.375, element->Id(), gmshfilecontent);
			}
		}
	}
	return;
} // StatMechManager::GmshOutputPeriodicBoundary()

/*----------------------------------------------------------------------*
 | plot the periodic boundary box                  (public) mueller 7/10|
 *----------------------------------------------------------------------*/
void StatMechManager::GmshOutputBoundaryBox(double boundarycolor,const std::ostringstream *filename)
{
	// plot the periodic box in case of periodic boundary conditions (first processor)
	if (statmechparams_.get<double> ("PeriodLength", 0.0) > 0 && discret_.Comm().MyPID() == 0)
	{
		FILE *fp = fopen(filename->str().c_str(), "a");
		std::stringstream gmshfilefooter;
		// get current period length
		double pl = statmechparams_.get<double> ("PeriodLength", 0.0);

		// define boundary lines
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << 0.0 << "," << 0.0 << "," << 0.0 << "," << pl << ","<< 0.0 << "," << 0.0;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 2
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << pl << "," << 0.0 << "," << 0.0 << "," << pl << "," << pl<< "," << 0.0;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 3
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << pl << "," << pl << "," << 0.0 << "," << pl << "," << pl<< "," << pl;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 4
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << pl << "," << pl << "," << pl << "," << 0.0 << "," << pl<< "," << pl;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 5
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << 0.0 << "," << pl << "," << pl << "," << 0.0 << "," << 0.0<< "," << pl;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 6
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << 0.0 << "," << 0.0 << "," << pl << "," << 0.0 << ","<< 0.0 << "," << 0.0;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 7
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << ","<< pl << "," << 0.0;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 8
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << 0.0 << "," << pl << "," << 0.0 << "," << pl << "," << pl<< "," << 0.0;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 9
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << 0.0 << "," << pl << "," << 0.0 << "," << 0.0 << "," << pl<< "," << pl;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 10
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << pl << "," << 0.0 << "," << 0.0 << "," << pl << "," << 0.0<< "," << pl;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 11
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << pl << "," << 0.0 << "," << pl << "," << pl << "," << pl<< "," << pl;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;
		// line 12
		gmshfilefooter << "SL(" << scientific;
		gmshfilefooter << pl << "," << 0.0 << "," << pl << "," << 0.0 << "," << 0.0<< "," << pl;
		gmshfilefooter << ")" << "{" << scientific << boundarycolor << ","<< boundarycolor << "};" << endl;

		fprintf(fp, gmshfilefooter.str().c_str());
		fclose(fp);
	}
	// wait for Proc 0 to catch up to the others
	discret_.Comm().Barrier();
}// StatMechManager::GmshOutputBoundaryBox

/*----------------------------------------------------------------------*
 | Gmsh output for crosslink molecule diffusion    (public) mueller 7/10|
 *----------------------------------------------------------------------*/
void StatMechManager::GmshOutputCrosslinkDiffusion(double color, const std::ostringstream *filename, const Epetra_Vector& disrow)
{
	// Visualization of crosslink molecule diffusion
	if (Teuchos::getIntegralValue<int>(statmechparams_, "CRSLNKDIFFUSION"))
	{
		// export row displacement to column map format
		Epetra_Vector discol(*(discret_.DofColMap()), true);
		LINALG::Export(disrow, discol);

		if (discret_.Comm().MyPID() == 0)
		{
			FILE *fp = fopen(filename->str().c_str(), "a");
			/*/ visualization of crosslink molecule positions by spheres on Proc 0
			std::stringstream gmshfilecross;
			for(int i=0; i<visualizepositions_->MyLength(); i++)
			{
				if((*numbond_)[i]<0.1)
				{
					double beadcolor = 5*color;
					//writing element by nodal coordinates as a sphere
					gmshfilecross << "SP(" << scientific;
					gmshfilecross<< (*visualizepositions_)[0][i]<< "," << (*visualizepositions_)[1][i] << "," << (*visualizepositions_)[2][i];
					gmshfilecross << ")" << "{" << scientific << beadcolor << "," << beadcolor << "};" << endl;
				}
			}
			fprintf(fp,gmshfilecross.str().c_str());
			fclose(fp);*/

			//special visualization for crosslink molecules with one/two bond(s); going through the Procs

			fp = fopen(filename->str().c_str(), "a");
			std::stringstream gmshfilebonds;

			// first, just update positions: redundant information on all procs
			for (int i=0; i<numbond_->MyLength(); i++)
			{
				switch ((int) (*numbond_)[i])
				{
					// crosslink molecule - filament node pair
					case 1:
					{
						// determine position of nodeGID entry
						int occupied = 0;
						for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
							if ((int) (*crosslinkerbond_)[j][i] != -1)
							{
								occupied = j;
								break;
							}
						int nodeGID = (int) (*crosslinkerbond_)[occupied][i];

						DRT::Node *node = discret_.lColNode(discret_.NodeColMap()->LID(nodeGID));
						LINALG::SerialDenseMatrix coord(3, 2, true);
						for (int j=0; j<coord.M(); j++)
						{
							int dofgid = discret_.Dof(node)[j];
							coord(j, 0) = node->X()[j] + discol[dofgid];
							coord(j, 1) = (*visualizepositions_)[j][i];
							//length += (coord(j,1)-coord(j,0))*(coord(j,1)-coord(j,0));
						}

						double beadcolor = 2*color;
						// in case of periodic boundary conditions
						if (statmechparams_.get<double> ("PeriodLength", 0.0) > 0.0)
						{
							GmshOutputPeriodicBoundary(coord, 2*color, gmshfilebonds, 0, true);
							// visualization of "real" crosslink molecule positions
							//beadcolor = 0.0; //black
							//gmshfilebonds << "SP(" << scientific;
							//gmshfilebonds << (*crosslinkerpositions_)[0][i] << ","<< (*crosslinkerpositions_)[1][i] << ","<< (*crosslinkerpositions_)[2][i];
							//gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
						}
						else
						{
							gmshfilebonds << "SL(" << scientific;
							gmshfilebonds << coord(0, 0) << "," << coord(1, 0) << ","<< coord(2, 0) << "," << coord(0, 1) << "," << coord(1, 1)<< "," << coord(2, 1);
							gmshfilebonds << ")" << "{" << scientific << 2*color << ","<< 2*color << "};" << endl;
							gmshfilebonds << "SP(" << scientific;
							gmshfilebonds << coord(0, 1) << "," << coord(1, 1) << ","<< coord(2, 1);
							gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
						}
					}
						break;
						// crosslinker element: crosslink molecule (representation) position (Proc 0 only)
					case 2:
					{
						// actual crosslinker element connecting two filaments (self-binding kinked crosslinkers are visualized in GmshKinkedVisual())
						if((*searchforneighbours_)[i] > 0.9)
						{
							if((*crosslinkonsamefilament_)[i] < 0.1)
							{
								double beadcolor = 4* color ; //red
								gmshfilebonds << "SP(" << scientific;
								gmshfilebonds << (*visualizepositions_)[0][i] << ","<< (*visualizepositions_)[1][i] << ","<< (*visualizepositions_)[2][i];
								gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
							}
						}
						else	// passive crosslink molecule
						{
							// determine position of nodeGID entry
							int occupied = 0;
							for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
								if ((int) (*crosslinkerbond_)[j][i] != -1)
								{
									occupied = j;
									break;
								}
							int nodeGID = (int) (*crosslinkerbond_)[occupied][i];

							DRT::Node *node = discret_.lColNode(discret_.NodeColMap()->LID(nodeGID));
							LINALG::SerialDenseMatrix coord(3, 2, true);
							for (int j=0; j<coord.M(); j++)
							{
								int dofgid = discret_.Dof(node)[j];
								coord(j, 0) = node->X()[j] + discol[dofgid];
								coord(j, 1) = (*visualizepositions_)[j][i];
							}

							double beadcolor = color;
							// in case of periodic boundary conditions
							if (statmechparams_.get<double> ("PeriodLength", 0.0) > 0.0)
							{
								GmshOutputPeriodicBoundary(coord, color, gmshfilebonds, 0, true);
								// visualization of "real" crosslink molecule positions
								//beadcolor = 0.0; //black
								//gmshfilebonds << "SP(" << scientific;
								//gmshfilebonds << (*crosslinkerpositions_)[0][i] << ","<< (*crosslinkerpositions_)[1][i] << ","<< (*crosslinkerpositions_)[2][i];
								//gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
							}
							else
							{
								gmshfilebonds << "SL(" << scientific;
								gmshfilebonds << coord(0, 0) << "," << coord(1, 0) << ","<< coord(2, 0) << "," << coord(0, 1) << "," << coord(1, 1)<< "," << coord(2, 1);
								gmshfilebonds << ")" << "{" << scientific << color << ","<< color << "};" << endl;
								gmshfilebonds << "SP(" << scientific;
								gmshfilebonds << coord(0, 1) << "," << coord(1, 1) << ","<< coord(2, 1);
								gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << endl;
							}
						}
					}
						break;
					default:
						continue;
				}
			}
			fprintf(fp, gmshfilebonds.str().c_str());
			fclose(fp);
		}
	}
	discret_.Comm().Barrier();
}// GmshOutputCrosslinkDiffusion

/*----------------------------------------------------------------------*
 | Special Gmsh output for crosslinkers occupying two binding spots on  |
 | the same filament											        (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void StatMechManager::GmshKinkedVisual(const LINALG::SerialDenseMatrix& coord, const double& color, int eleid, std::stringstream& gmshfilecontent)
{
	/* We need a third point in order to visualize the crosslinker.
	 * It marks the location of the kink.
	 */
	std::vector<double> thirdpoint(3,0.0);
	// get the element
	DRT::Element* element = discret_.gElement(eleid);

	// calculate tangent
	double ltan = 0.0;
	for (int j=0; j<coord.M(); j++)
		ltan += (coord(j, coord.N() - 1) - coord(j, 0)) * (coord(j, coord.N() - 1) - coord(j, 0));
	ltan = sqrt(ltan);

	std::vector<double> t(3, 0.0);
	for (int j=0; j<coord.M(); j++)
		t.at(j) = ((coord(j, coord.N() - 1) - coord(j, 0)) / ltan);

	// calculate normal via cross product: [0 0 1]x[tx ty tz]
	std::vector<double> n(3, 0.0);
	n.at(0) = -t.at(1);
	n.at(1) = t.at(0);
	// norm it since the cross product does not keep the length
	double lnorm = 0.0;
	for (int j=0; j<(int)n.size(); j++)
		lnorm += n.at(j) * n.at(j);
	lnorm = sqrt(lnorm);
	for (int j=0; j<(int) n.size(); j++)
		n.at(j) /= lnorm;

	// by modulo operation involving the node IDs
	double alpha = fmod((double) (element->Nodes()[element->NumNode() - 1]->Id() + element->Nodes()[0]->Id()), 2*M_PI);

	// rotate the normal by alpha
	LINALG::SerialDenseMatrix RotMat(3, 3);
	// build the matrix of rotation
	for (int j=0; j<3; j++)
		RotMat(j, j) = cos(alpha) + t.at(j) * t.at(j) * (1 - cos(alpha));
	RotMat(0, 1) = t.at(0) * t.at(1) * (1 - cos(alpha)) - t.at(2) * sin(alpha);
	RotMat(0, 2) = t.at(0) * t.at(2) * (1 - cos(alpha)) + t.at(1) * sin(alpha);
	RotMat(1, 0) = t.at(1) * t.at(0) * (1 - cos(alpha)) + t.at(2) * sin(alpha);
	RotMat(1, 2) = t.at(1) * t.at(2) * (1 - cos(alpha)) - t.at(0) * sin(alpha);
	RotMat(2, 0) = t.at(2) * t.at(0) * (1 - cos(alpha)) - t.at(1) * sin(alpha);
	RotMat(2, 1) = t.at(2) * t.at(1) * (1 - cos(alpha)) + t.at(0) * sin(alpha);

	// rotation
	std::vector<double> nrot(3, 0.0);
	for (int j=0; j<3; j++)
		for (int k=0; k<3; k++)
			nrot.at(j) += RotMat(j, k) * n.at(k);

	// calculation of the third point lying in the direction of the rotated normal
	// height of the third point above the filament
	double h = 0.33 * (statmechparams_.get<double> ("R_LINK", 0.0) + statmechparams_.get<double> ("DeltaR_LINK", 0.0));
	for (int j=0; j<3; j++)
		thirdpoint.at(j) = (coord(j, 0) + coord(j, element->NumNode() - 1)) / 2.0	+ h * nrot.at(j);

	gmshfilecontent << "SL(" << scientific << coord(0,0) << ","<< coord(1,0) << "," << coord(2,0) << ","
									<< thirdpoint.at(0) << ","<< thirdpoint.at(1) << "," << thirdpoint.at(2) << ")"
									<< "{" << scientific<< color << "," << color << "};" << endl;
	gmshfilecontent << "SL(" << scientific << thirdpoint.at(0) << "," << thirdpoint.at(1)<< "," << thirdpoint.at(2) << ","
									<< coord(0, 1) << "," << coord(1,1) << "," << coord(2,1) << ")"
									<< "{" << scientific<< color << "," << color << "};" << endl;
	gmshfilecontent << "SP(" << scientific
									<< thirdpoint.at(0) << "," << thirdpoint.at(1) << ","<< thirdpoint.at(2)
									<< ")" << "{" << scientific << color << ","<< color << "};" << endl;

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
	return;
}// StatMechManager::GmshKinkedVisual

/*----------------------------------------------------------------------*
 | prepare visualization vector for Gmsh Output  (private) mueller 08/10|
 *----------------------------------------------------------------------*/
void StatMechManager::GmshPrepareVisualization(const Epetra_Vector& dis)
{
	double ronebond = statmechparams_.get<double> ("R_LINK", 0.0) / 2.0;

	// column map displacement and displacement increment

	Epetra_Vector discol(*(discret_.DofColMap()), true);
	LINALG::Export(dis, discol);

	if (discret_.Comm().MyPID() == 0)
	{
		for (int i=0; i<numbond_->MyLength(); i++)
		{
			switch((int)(*numbond_)[i])
			{
				// diffusion
				case 0:
				{
					for (int j=0; j<visualizepositions_->NumVectors(); j++)
						(*visualizepositions_)[j][i] = (*crosslinkerpositions_)[j][i];
				}
				break;
				// one bonded crosslink molecule
				case 1:
				{
					// determine position of nodeGID entry
					int occupied = -1;
					for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
						if ((*crosslinkerbond_)[j][i] > -0.9)
						{
							occupied = j;
							break;
						}

					int nodeLID = discret_.NodeColMap()->LID((int)(*crosslinkerbond_)[occupied][i]);
					int currfilament = (int)(*filamentnumber_)[nodeLID];
					const DRT::Node *node0 = discret_.lColNode(nodeLID);
					// choose a second (neighbour) node
					const DRT::Node *node1 = NULL;
					if(nodeLID < basisnodes_-1)
						if((*filamentnumber_)[nodeLID+1]==currfilament)
							node1 = discret_.lColNode(nodeLID+1);
						else
							node1 = discret_.lColNode(nodeLID-1);
					if(nodeLID == basisnodes_-1)
						if((*filamentnumber_)[nodeLID-1]==currfilament)
							node1 = discret_.lColNode(nodeLID-1);

					//calculate unit tangent
					LINALG::Matrix<3,1> nodepos0;
					LINALG::Matrix<3,1> tangent;

					for (int j=0; j<3; j++)
					{
						int dofgid0 = discret_.Dof(node0)[j];
						int dofgid1 = discret_.Dof(node1)[j];
						nodepos0(j,0) = node0->X()[j] + discol[discret_.DofColMap()->LID(dofgid0)];
						double nodeposj1 = node1->X()[j] + discol[discret_.DofColMap()->LID(dofgid1)];
						tangent(j) = nodeposj1 - nodepos0(j,0);
					}
					tangent.Scale(1 / tangent.Norm2());

					// calculate normal via cross product: [0 0 1]x[tx ty tz]
					LINALG::Matrix<3, 1> normal;
					normal(0) = -tangent(1);
					normal(1) = tangent(0);
					// norm it since the cross product does not keep the length
					normal.Scale(1 / normal.Norm2());

					// obtain angle
					// random angle
					// by modulo operation involving the crosslink molecule number
					double alpha = fmod((double) i, 2*M_PI);

					// rotate the normal by alpha
					LINALG::Matrix<3, 3> RotMat;
					// build the matrix of rotation
					for (int j=0; j<3; j++)
						RotMat(j, j) = cos(alpha) + tangent(j) * tangent(j) * (1 - cos(alpha));
					RotMat(0, 1) = tangent(0) * tangent(1) * (1 - cos(alpha)) - tangent(2) * sin(alpha);
					RotMat(0, 2) = tangent(0) * tangent(2) * (1 - cos(alpha)) + tangent(1) * sin(alpha);
					RotMat(1, 0) = tangent(1) * tangent(0) * (1 - cos(alpha)) + tangent(2) * sin(alpha);
					RotMat(1, 2) = tangent(1) * tangent(2) * (1 - cos(alpha)) - tangent(0) * sin(alpha);
					RotMat(2, 0) = tangent(2) * tangent(0) * (1 - cos(alpha)) - tangent(1) * sin(alpha);
					RotMat(2, 1) = tangent(2) * tangent(1) * (1 - cos(alpha)) + tangent(0) * sin(alpha);

					// rotation
					LINALG::Matrix<3, 1> rotnormal;
					rotnormal.Multiply(RotMat, normal);

					// calculation of the visualized point lying in the direction of the rotated normal
					for (int j=0; j<visualizepositions_->NumVectors(); j++)
						(*visualizepositions_)[j][i] = nodepos0(j,0) + ronebond*rotnormal(j,0);
				}
				break;
				case 2:
				{
					// actual crosslinker element (not kinked)
					if((*searchforneighbours_)[i]>0.9)
					{
						// loop over filament node components
						for (int j=0; j<visualizepositions_->NumVectors(); j++)
						{
							std::vector<double> dofnodepositions(crosslinkerbond_->NumVectors(), 0.0);
							// loop over filament node GIDs
							for (int k=0; k<crosslinkerbond_->NumVectors(); k++)
							{
								int nodeGID = (int) (*crosslinkerbond_)[k][i];
								DRT::Node *node = discret_.lColNode(discret_.NodeColMap()->LID(nodeGID));
								int dofgid = discret_.Dof(node)[j];
								dofnodepositions.at(k) = node->X()[j] + discol[dofgid];
							}
							/* Check if the crosslinker element is broken/discontinuous; if so, reposition the second nodal value.
							 * It does not matter which value is shifted as long the shift occurs in a consistent way.*/
							if (statmechparams_.get<double> ("PeriodLength", 0.0) > 0.0)
							{
								(*visualizepositions_)[j][i] = dofnodepositions.at(0);
								for (int k=0; k<1; k++)
								{
									// shift position if it is found to be outside the boundary box
									if (fabs(dofnodepositions.at(k+1) - statmechparams_.get<double> ("PeriodLength", 0.0) - dofnodepositions.at(k))< fabs(dofnodepositions.at(k+1) - dofnodepositions.at(k)))
										dofnodepositions.at(k+1) -= statmechparams_.get<double> ("PeriodLength", 0.0);
									if (fabs(dofnodepositions.at(k+1) + statmechparams_.get<double> ("PeriodLength", 0.0) - dofnodepositions.at(k))< fabs(dofnodepositions.at(k+1) - dofnodepositions.at(k)))
										dofnodepositions.at(k+1) += statmechparams_.get<double> ("PeriodLength", 0.0);

									(*visualizepositions_)[j][i] += dofnodepositions.at(k+1);
								}
								// new crosslink molecule position
								(*visualizepositions_)[j][i] /= 2.0;
							}
							else
								(*visualizepositions_)[j][i] /= 2.0;
						}
					}
					else	// passive crosslink molecule
					{
						// determine position of nodeGID entry
						int occupied = -1;
						for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
							if ((*crosslinkerbond_)[j][i] > -0.9)
							{
								occupied = j;
								break;
							}

						int nodeLID = discret_.NodeColMap()->LID((int)(*crosslinkerbond_)[occupied][i]);
						int currfilament = (int)(*filamentnumber_)[nodeLID];
						const DRT::Node *node0 = discret_.lColNode(nodeLID);
						// choose a second (neighbour) node
						const DRT::Node *node1 = NULL;
						if(nodeLID < basisnodes_-1)
							if((*filamentnumber_)[nodeLID+1]==currfilament)
								node1 = discret_.lColNode(nodeLID+1);
							else
								node1 = discret_.lColNode(nodeLID-1);
						if(nodeLID == basisnodes_-1)
							if((*filamentnumber_)[nodeLID-1]==currfilament)
								node1 = discret_.lColNode(nodeLID-1);

						//calculate unit tangent
						LINALG::Matrix<3,1> nodepos0;
						LINALG::Matrix<3,1> tangent;

						for (int j=0; j<3; j++)
						{
							int dofgid0 = discret_.Dof(node0)[j];
							int dofgid1 = discret_.Dof(node1)[j];
							nodepos0(j,0) = node0->X()[j] + discol[discret_.DofColMap()->LID(dofgid0)];
							double nodeposj1 = node1->X()[j] + discol[discret_.DofColMap()->LID(dofgid1)];
							tangent(j) = nodeposj1 - nodepos0(j,0);
						}
						tangent.Scale(1 / tangent.Norm2());

						// calculate normal via cross product: [0 0 1]x[tx ty tz]
						LINALG::Matrix<3, 1> normal;
						normal(0) = -tangent(1);
						normal(1) = tangent(0);
						// norm it since the cross product does not keep the length
						normal.Scale(1 / normal.Norm2());

						// obtain angle
						// random angle
						// by modulo operation involving the crosslink molecule number
						double alpha = fmod((double) i, 2*M_PI);

						// rotate the normal by alpha
						LINALG::Matrix<3, 3> RotMat;
						// build the matrix of rotation
						for (int j=0; j<3; j++)
							RotMat(j, j) = cos(alpha) + tangent(j) * tangent(j) * (1 - cos(alpha));
						RotMat(0, 1) = tangent(0) * tangent(1) * (1 - cos(alpha)) - tangent(2) * sin(alpha);
						RotMat(0, 2) = tangent(0) * tangent(2) * (1 - cos(alpha)) + tangent(1) * sin(alpha);
						RotMat(1, 0) = tangent(1) * tangent(0) * (1 - cos(alpha)) + tangent(2) * sin(alpha);
						RotMat(1, 2) = tangent(1) * tangent(2) * (1 - cos(alpha)) - tangent(0) * sin(alpha);
						RotMat(2, 0) = tangent(2) * tangent(0) * (1 - cos(alpha)) - tangent(1) * sin(alpha);
						RotMat(2, 1) = tangent(2) * tangent(1) * (1 - cos(alpha)) + tangent(0) * sin(alpha);

						// rotation
						LINALG::Matrix<3, 1> rotnormal;
						rotnormal.Multiply(RotMat, normal);

						// calculation of the visualized point lying in the direction of the rotated normal
						for (int j=0; j<visualizepositions_->NumVectors(); j++)
							(*visualizepositions_)[j][i] = nodepos0(j,0) + ronebond*rotnormal(j,0);
					}
				}
				break;
			}
		}
		if (statmechparams_.get<double>("PeriodLength", 0.0) > 0.0)
			CrosslinkerPeriodicBoundaryShift(*visualizepositions_);
	}
	else
	{
		visualizepositions_->PutScalar(0.0);
	}
	// synchronize results
	Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
	Epetra_Import crosslinkimporter(*crosslinkermap_, *transfermap_);
	Epetra_MultiVector visualizepositionstrans(*transfermap_, 3, true);
	visualizepositionstrans.Export(*visualizepositions_, crosslinkexporter, Add);
	visualizepositions_->Import(visualizepositionstrans, crosslinkimporter, Insert);
	// debug couts
	//cout<<*visualizepositions_<<endl;
}//GmshPrepareVisualization

/*----------------------------------------------------------------------*
 | initialize special output for statistical mechanics(public)cyron 12/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechInitOutput(const int ndim, const double& dt)
{
	//initializing special output for statistical mechanics by looking for a suitable name of the outputfile and setting up an empty file with this name

	switch (Teuchos::getIntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT"))
	{
		case INPAR::STATMECH::statout_endtoendlog:
		{

		  //output is written on proc 0 only
      if(!discret_.Comm().MyPID())
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
        if (fpnumbering == NULL)
        {
          do
          {
            outputfilenumber_++;
            outputfilename.str("");
            outputfilename << "EndToEnd" << outputfilenumber_ << ".dat";
            fp = fopen(outputfilename.str().c_str(), "r");
          }
          while (fp != NULL);

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
            if (tok == "Next")
            {
              f >> tok;
              if (tok == "Number:")
                f >> outputfilenumber_;
            }
          } //while(f)

          //defining outputfilename by means of new testnumber
          outputfilename.str("");
          outputfilename << "EndToEnd" << outputfilenumber_ << ".dat";
        }

        //increasing the number in the numbering file by one
        fpnumbering = fopen(numberingfilename.str().c_str(), "w");
        std::stringstream filecontent;
        filecontent << "Next Number: " << (outputfilenumber_ + 1);
        fprintf(fpnumbering, filecontent.str().c_str());
        fclose(fpnumbering);
      }

		}
			break;
		case INPAR::STATMECH::statout_endtoendconst:
		{
      //output is written on proc 0 only
      if(!discret_.Comm().MyPID())
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
        discret_.GetCondition("PointNeumann", pointneumannconditions);
        if (pointneumannconditions.size() > 0)
        {
          const vector<double>* val = pointneumannconditions[0]->Get<vector<double> > ("val");
          neumannforce = fabs((*val)[0]);
        }
        else
          neumannforce = 0;

        //if there is no such numbering file: look for a not yet existing output file name (numbering upwards)
        if (fpnumbering == NULL)
        {
          do
          {
            //defining name of output file
            outputfilenumber_++;
            outputfilename.str("");
            outputfilename << "E2E_" << discret_.NumMyRowElements() << '_' << dt<< '_' << neumannforce << '_' << outputfilenumber_ << ".dat";
            fp = fopen(outputfilename.str().c_str(), "r");
          }
          while (fp != NULL);

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
            if (tok == "Next")
            {
              f >> tok;
              if (tok == "Number:")
                f >> outputfilenumber_;
            }
          } //while(f)

          //defining outputfilename by means of new testnumber
          outputfilename.str("");
          outputfilename << "E2E_" << discret_.NumMyRowElements() << '_' << dt << '_' << neumannforce << '_' << outputfilenumber_ << ".dat";

          //set up new file with name "outputfilename" without writing anything into this file
          fp = fopen(outputfilename.str().c_str(), "w");
          fclose(fp);
        }

        //increasing the number in the numbering file by one
        fpnumbering = fopen(numberingfilename.str().c_str(), "w");
        std::stringstream filecontent;
        filecontent << "Next Number: " << (outputfilenumber_ + 1);
        //write fileconent into file!
        fprintf(fpnumbering, filecontent.str().c_str());
        //close file
        fclose(fpnumbering);
      }

		}
    break;
    /*computing and writing into file data about correlation of orientation of different elements
     *as it is considered in the context of the persistence length*/
		case INPAR::STATMECH::statout_orientationcorrelation:
		{
      //output is written on proc 0 only
      if(!discret_.Comm().MyPID())
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
        if (fpnumbering == NULL)
        {
          do
          {
            outputfilenumber_++;
            outputfilename.str("");
            outputfilename << "OrientationCorrelation" << outputfilenumber_<< ".dat";
            fp = fopen(outputfilename.str().c_str(), "r");
          }
          while (fp != NULL);

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
            if (tok == "Next")
            {
              f >> tok;
              if (tok == "Number:")
                f >> outputfilenumber_;
            }
          } //while(f)

          //defining outputfilename by means of new testnumber
          outputfilename.str("");
          outputfilename << "OrientationCorrelation" << outputfilenumber_<< ".dat";
        }

        //set up new file with name "outputfilename" without writing anything into this file
        fp = fopen(outputfilename.str().c_str(), "w");
        fclose(fp);

        //increasing the number in the numbering file by one
        fpnumbering = fopen(numberingfilename.str().c_str(), "w");
        std::stringstream filecontent;
        filecontent << "Next Number: " << (outputfilenumber_ + 1);
        //write fileconent into file!
        fprintf(fpnumbering, filecontent.str().c_str());
        //close file
        fclose(fpnumbering);
      }
		}
    break;
    //simulating diffusion coefficient for anisotropic friction
		case INPAR::STATMECH::statout_anisotropic:
		{
      //output is written on proc 0 only
      if(!discret_.Comm().MyPID())
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
        if (fpnumbering == NULL)
        {
          do
          {
            outputfilenumber_++;
            outputfilename.str("");
            outputfilename << "AnisotropicDiffusion" << outputfilenumber_
                << ".dat";
            fp = fopen(outputfilename.str().c_str(), "r");
          }
          while (fp != NULL);

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
            if (tok == "Next")
            {
              f >> tok;
              if (tok == "Number:")
                f >> outputfilenumber_;
            }
          } //while(f)

          //defining outputfilename by means of new testnumber
          outputfilename.str("");
          outputfilename << "AnisotropicDiffusion" << outputfilenumber_ << ".dat";
        }

        //set up new file with name "outputfilename" without writing anything into this file
        fp = fopen(outputfilename.str().c_str(), "w");
        fclose(fp);

        //increasing the number in the numbering file by one
        fpnumbering = fopen(numberingfilename.str().c_str(), "w");
        std::stringstream filecontent;
        filecontent << "Next Number: " << (outputfilenumber_ + 1);
        //write fileconent into file!
        fprintf(fpnumbering, filecontent.str().c_str());
        //close file
        fclose(fpnumbering);

        //initializing variables for positions of first and last node at the beginning
        beginold_.PutScalar(0);
        endold_.PutScalar(0);
        for (int i=0; i<ndim; i++)
        {
          beginold_(i) = (discret_.gNode(0))->X()[i];
          endold_(i) = (discret_.gNode(discret_.NumMyRowNodes() - 1))->X()[i];
        }

        for (int i=0; i<3; i++)
          sumdispmiddle_(i, 0) = 0.0;

        sumsquareincpar_ = 0.0;
        sumsquareincort_ = 0.0;
        sumrotmiddle_ = 0.0;
        sumsquareincmid_ = 0.0;
        sumsquareincrot_ = 0.0;
      }
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
			outputfilename << "ViscoElOutputProc" << discret_.Comm().MyPID()<< ".dat";

			fp = fopen(outputfilename.str().c_str(), "w");

			//filecontent << "Output for measurement of viscoelastic properties written by processor "<< discret_.Comm().MyPID() << endl;

			// move temporary stringstream to file and close it
			fprintf(fp, filecontent.str().c_str());
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
 | output for structural polymorphism             (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void StatMechManager::StructPolymorphOutput(const Epetra_Vector& disrow, const std::ostringstream& filename1, const std::ostringstream& filename2)
{
	/* The following code is basically copied from GmshOutput() and therefore remains largely uncommented.
	 * The output consists of a vector v E R³ for each element and the respective node coordinates C.
	 * Since the output is redundant in the case of ghosted elements (parallel use), the GID of the
	 * element is written, too. This ensures an easy adaption of the output for post-processing.
	 *
	 * columns of the output file: element-ID, filament-ID, v_x, v_y, v_z, C_0_x,  C_0_y, C_0_z, C_1_x, C_1_y, C_1_z
	 * */

	FILE* fp = NULL;

	fp = fopen(filename1.str().c_str(), "w");

	std::stringstream filamentcoords;

	// get filament number conditions
	vector<DRT::Condition*> filamentnumbers(0);
	discret_.GetCondition("FilamentNumber", filamentnumbers);

	for (int i=0; i<discret_.NumMyRowElements(); ++i)
	{
		DRT::Element* element = discret_.lRowElement(i);

		// we only need filament elements, crosslinker elements are neglected with the help of their high element IDs
		if (element->Id() > basisnodes_)
			continue;

		LINALG::SerialDenseMatrix coord(3, element->NumNode(), true);
		// filament number
		int filnumber = 0;
		// a switch to avoid redundant search for filament number
		bool done = false;

		for (int j=0; j<3; j++)
		{
			for (int k=0; k<element->NumNode(); k++)
			{
				double referenceposition = ((element->Nodes())[k])->X()[j];
				vector<int> dofnode = discret_.Dof((element->Nodes())[k]);
				double displacement = disrow[discret_.DofRowMap()->LID(dofnode[j])];
				coord(j,k) = referenceposition + displacement;

				// get corresponding filament number
				// loop over all filaments
				if (!done)
					for (int l = 0; l < (int) filamentnumbers.size(); l++)
						for (int m = 0; m < (int) filamentnumbers[l]->Nodes()->size(); m++)
							if (element->Nodes()[k]->Id() == filamentnumbers[l]->Nodes()->at(m))
							{
								filnumber = filamentnumbers[l]->Id();
								done = true;
							}
			}
		}

		const DRT::ElementType & eot = element->ElementType();

#ifdef D_BEAM3
		if (eot == DRT::ELEMENTS::Beam3Type::Instance())
			for (int j=0; j<element->NumNode()-1; j++)
				filamentcoords << element->Id() << " " << filnumber << "   "
											 << std::scientific << std::setprecision(15)//<<coord(0,j+1)-coord(0,j) << " " <<coord(1,j+1)-coord(1,j) << " " <<coord(2,j+1)-coord(2,j)<<"   "
											 << coord(0, j) << " " << coord(1, j) << " " << coord(2, j) << " "
											 << coord(0, j + 1) << " " << coord(1, j + 1) << " " << coord(2, j+ 1) << endl;
		else
#endif
#ifdef D_BEAM3II
		if(eot==DRT::ELEMENTS::Beam3iiType::Instance())
		for(int j=0; j<element->NumNode()-1; j++)
		filamentcoords<<element->Id()<<" "<<filnumber<<"   "<<std::scientific<<std::setprecision(15)//<<coord(0,j+1)-coord(0,j) << " " <<coord(1,j+1)-coord(1,j) << " " <<coord(2,j+1)-coord(2,j)<<"   "
									<< coord(0,j) << " " << coord(1,j) << " " << coord(2,j) << " "
									<< coord(0,j+1) << " " << coord(1,j+1) << " " << coord(2,j+1)<<endl;
		else
#endif
#ifdef D_TRUSS3
		if (eot == DRT::ELEMENTS::Truss3Type::Instance())
			for (int j=0; j<element->NumNode()-1; j++)
				filamentcoords << element->Id() << " " << filnumber << "   "
										   << std::scientific << std::setprecision(15)//<<coord(0,j+1)-coord(0,j) << " " <<coord(1,j+1)-coord(1,j) << " " <<coord(2,j+1)-coord(2,j)<<"   "
											 << coord(0, j) << " " << coord(1, j) << " " << coord(2, j) << " "
											 << coord(0, j + 1) << " " << coord(1, j + 1) << " " << coord(2, j+ 1) << endl;
		else
#endif
		{
		}
	}
	//write content into file and close it
	fprintf(fp, filamentcoords.str().c_str());
	fclose(fp);

	/* A measure to quantify bundling, is hereafter established and written to a separate output file.
	 * It works as follows:
	 * 1) Specify all possible filament pairs P(i,j)
	 * 2) n(i,j) counts the number of crosslinkers of P(i,j)
	 * 3) n_av determines the average number of non-zero n(i,j) [e.g. n_av==1 indicates a perfectly
	 *    homogenous isotropic network]
	 */

	//if(discret_.Comm().MyPID()==0)
	//{
	LINALG::SerialDenseMatrix ncross((int) filamentnumbers.size(), (int) filamentnumbers.size(), true);
	// loop over all possible filament pairs and their nodes (wow!)
	// filament i
	for (int i=0; i<(int)filamentnumbers.size() - 1; i++)
	{
		// filament j
		for (int j=i+1; j<(int)filamentnumbers.size(); j++)
		{
			// number of crosslinkers of filament pair P(i,j)
			int n_ij = 0;
			// node k of filament i
			for (int k=0; k<(int)filamentnumbers[i]->Nodes()->size(); k++)
			{
				// get node k of filament i if available on proc, else: next node in condition
				if (discret_.HaveGlobalNode(filamentnumbers[i]->Nodes()->at(k)))
				{
					//cout<<"Node_k="<<k<<",i="<<i<<" is on Proc "<<discret_.Comm().MyPID()<<endl;
					DRT::Node* node_ki = discret_.gNode(filamentnumbers[i]->Nodes()->at(k));
					// node l of filament j
					for (int l = 0; l < (int) filamentnumbers[j]->Nodes()->size(); l++)
					{
						// get node l of filament j if available on proc, else: next node in condition
						if (discret_.HaveGlobalNode(filamentnumbers[j]->Nodes()->at(l)))
						{
							//cout<<"Node_l="<<l<<",j="<<j<<" is on Proc "<<discret_.Comm().MyPID()<<endl;
							DRT::Node* node_lj = discret_.gNode(filamentnumbers[j]->Nodes()->at(l));
							// identify any crosslinkers with nodes k and l
							DRT::Element** adjelement_k = node_ki->Elements();
							DRT::Element** adjelement_l = node_lj->Elements();

							// loop over Elements adjacent to nodes k and l and increment n_ij if the element
							// GIDs of element m and n match, i.e. the two filament share a crosslinker
							for (int m = 0; m < node_ki->NumElement(); m++)
							{
								if (adjelement_k[m]->Id() > basisnodes_)
									for (int n = 0; n < node_lj->NumElement(); n++)
										if (adjelement_k[m]->Id() == adjelement_l[n]->Id())
											n_ij++;
							}
						}
					}
				}
			}
			// store the number of crosslinker of P(i,j) at the respective position in the matrix
			ncross(i, j) = n_ij;
		}
	}

	// write the matrix ncross into a file via stream
	fp = fopen(filename2.str().c_str(), "w");
	std::stringstream crosslinkercount;

	for (int i=0; i<ncross.M(); i++)
	{
		for (int j=0; j<ncross.N(); j++)
			crosslinkercount << ncross(i, j) << "  ";
		crosslinkercount << endl;
	}

	//write content into file and close it
	fprintf(fp, crosslinkercount.str().c_str());
	fclose(fp);
	//}
}// StatMechManager:StructPolymorphOutput()

/*----------------------------------------------------------------------*
 | output for density-density-correlation-function(public) mueller 07/10|
 *----------------------------------------------------------------------*/
void StatMechManager::DensityDensityCorrOutput(const std::ostringstream& filename)
{
	/* Calculate the distances for all tuples of crosslink molecules.
	 * Each processor calculates the distances between its row map molecules and column map molecules
	 * Since we compare a set to itself, we just calculate one half of the matrix ( (n²-n)/2 calculations for n molecules)*/
	int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 0);
	// number of overall independent combinations
	int numcombinations = (statmechparams_.get<int>("N_crosslink", 0)*statmechparams_.get<int>("N_crosslink", 0)-statmechparams_.get<int>("N_crosslink", 0))/2;
	// combinations on each processor
	int combinationsperproc = (int)floor((double)numcombinations/(double)discret_.Comm().NumProc());
	int remainder = numcombinations%combinationsperproc;

	// get starting index tuples for later use
	if(!isinit_)
	{
		// first entry
		startindex_.push_back(std::vector<int>(2,0));
		for(int mypid=0; mypid<discret_.Comm().NumProc(); mypid++)
		{
			std::vector<int> start(2,0);
			bool continueloop = false;
			bool quitloop = false;
			int counter = 0;
			int appendix = 0;
			if(mypid==discret_.Comm().NumProc()-1)
				appendix = remainder;

			// loop over crosslinker pairs
			for(int i=0; i<crosslinkermap_->NumMyElements(); i++)
			{
				for(int j=0; j<crosslinkermap_->NumMyElements(); j++)
				{
					if(i==startindex_[mypid][0] && j==startindex_[mypid][1])
						continueloop = true;
					if(j>i && continueloop)
					{
						if(counter<combinationsperproc+appendix)
							counter++;
						else
						{
							// new start index j
							if(j==crosslinkermap_->NumMyElements()-1)
								start[1] = 0;
							else
								start[1] = j;
							quitloop = true;
							break;
						}
					}
				}
				if(quitloop)
				{
					// new start index i
					if(start[1]==0)
						start[0] = i+1;
					else
						start[0] = i;
					// new start tuple
					startindex_.push_back(start);
					break;
				}
			}
		}
		isinit_ = true;
	}

	int combicount = 0;
	Epetra_Vector crosslinksperbinrow(*ddcorrrowmap_, true);
	// loop over crosslinkermap_ (column map, same for all procs: maps all crosslink molecules)
	for(int mypid=0; mypid<discret_.Comm().NumProc(); mypid++)
	{
		bool quitloop = false;
		if(mypid==discret_.Comm().MyPID())
		{
			bool continueloop = false;
			int appendix = 0;
			if(mypid==discret_.Comm().NumProc()-1)
				appendix = remainder;

			for(int i=0; i<crosslinkermap_->NumMyElements(); i++)
			{
				for(int j=0; j<crosslinkermap_->NumMyElements(); j++)
				{
					// start adding crosslink from here
					if(i==startindex_[mypid][0] && j==startindex_[mypid][1])
						continueloop = true;
					// only entries above main diagonal and within limits of designated number of crosslink molecules per processor
					if(j>i && continueloop)
						if(combicount<combinationsperproc+appendix)
						{
							combicount++;
							double periodlength = statmechparams_.get<double>("PeriodLength", 0.0);
							double deltaxij = 0.0;

							// calculate the distance between molecule i and molecule j
							for(int k=0; k<crosslinkerpositions_->NumVectors(); k++)
								deltaxij += ((*crosslinkerpositions_)[k][i]-(*crosslinkerpositions_)[k][j])*((*crosslinkerpositions_)[k][i]-(*crosslinkerpositions_)[k][j]);
							deltaxij = sqrt(deltaxij);
							// calculate the actual bin to which the current distance belongs and increment the count for that bin (currbin=LID)
							int currbin = (int)floor(deltaxij/(periodlength*sqrt(3.0))*numbins);
							// in case the distance is exactly periodlength*sqrt(3)
							if(currbin==numbins)
								currbin--;
							crosslinksperbinrow[currbin] += 1.0;
						}
						else
						{
							quitloop = true;
							break;
						}
				}
				if(quitloop)
					break;
			}
			if(quitloop)
				break;
		}
		else
			continue;
	}

	// Export
	Epetra_Vector crosslinksperbincol(*ddcorrcolmap_, true);
	Epetra_Import importer(*ddcorrcolmap_, *ddcorrrowmap_);
	crosslinksperbincol.Import(crosslinksperbinrow, importer, Insert);

	// Add the processor-specific data up
	std::vector<int> crosslinksperbin(numbins, 0);
	for(int i=0; i<numbins; i++)
		for(int pid=0; pid<discret_.Comm().NumProc(); pid++)
			crosslinksperbin[i] += (int)crosslinksperbincol[pid*numbins+i];

	// write data to file
	if(!discret_.Comm().MyPID())
	{
		FILE* fp = NULL;
		fp = fopen(filename.str().c_str(), "w");
		std::stringstream histogram;

		for(int i=0; i<numbins; i++)
			histogram<<i+1<<"    "<<crosslinksperbin[i]<<endl;
		//write content into file and close it
		fprintf(fp, histogram.str().c_str());
		fclose(fp);
	}
}

/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechUpdate(const int& istep, const double dt, Epetra_Vector& disrow,RCP<LINALG::SparseOperator>& stiff, int ndim)
{
#ifdef MEASURETIME
	const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME
	/*first we modify the displacement vector so that current nodal position at the end of current time step complies with
	 * periodic boundary conditions, i.e. no node lies outside a cube of edge length PeriodLength*/

	Epetra_Vector disrowpre = disrow;
	//cout<<"disrow.MyLength()"<<disrow.MyLength()<<endl;
	PeriodicBoundaryShift(disrow, ndim, dt);

	// debuggin cout
	/*for(int i=0; i<discret_.NumMyRowNodes(); i++)
	{
		DRT::Node* node = discret_.lRowNode(i);
		std::vector<int> dofnode = discret_.Dof(node);
		for(int j=0; j<3; j++)
			if(node->X()[j]+disrow[discret_.DofRowMap()->LID(dofnode[j])]<0.0 || node->X()[j]+disrow[discret_.DofRowMap()->LID(dofnode[j])]>statmechparams_.get<double>("PeriodLength", 0.0))
			{
				cout<<"Proc "<<discret_.Comm().MyPID()<<": disrowpre entry, node "<<i<<", dof "<<dofnode[j]<<", Numele "<<node->NumElement()<<": "<<node->X()[j]+disrowpre[discret_.DofRowMap()->LID(dofnode[j])]<<endl;
				cout<<"Proc "<<discret_.Comm().MyPID()<<": disrow entry   , node "<<i<<", dof "<<dofnode[j]<<", Numele "<<node->NumElement()<<": "<<node->X()[j]+disrow[discret_.DofRowMap()->LID(dofnode[j])]<<endl;
			}
	}*/

	//if dynamic crosslinkers are used update comprises adding and deleting crosslinkers
	if (Teuchos::getIntegralValue<int>(statmechparams_, "DYN_CROSSLINKERS"))
	{
		// crosslink molecule diffusion
		if (Teuchos::getIntegralValue<int>(statmechparams_, "CRSLNKDIFFUSION"))
		{
			double standarddev = sqrt(statmechparams_.get<double> ("KT", 0.0) / (2*M_PI * statmechparams_.get<double> ("ETA", 0.0) * statmechparams_.get<double> ("R_LINK", 0.0)) * dt);
			CrosslinkerDiffusion(disrow, 0.0, standarddev, dt);
		}

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
		std::map<int, LINALG::Matrix<3, 1> > currentpositions;
		std::map<int, LINALG::Matrix<3, 1> > currentrotations;
		currentpositions.clear();
		currentrotations.clear();

		/*note: access by ExtractMyValues requires column map vector, whereas displacements on level of time integration are
		 * handled as row map vector*/
		Epetra_Vector discol(*discret_.DofColMap(), true);

		LINALG::Export(disrow, discol);

		for (int i=0; i<discret_.NumMyColNodes(); ++i)
		{

			//get pointer at a node
			const DRT::Node* node = discret_.lColNode(i);

			//get GIDs of this node's degrees of freedom
			std::vector<int> dofnode = discret_.Dof(node);

			LINALG::Matrix<3, 1> currpos;
			LINALG::Matrix<3, 1> currrot;

			currpos(0) = node->X()[0] + discol[discret_.DofColMap()->LID(dofnode[0])];
			currpos(1) = node->X()[1] + discol[discret_.DofColMap()->LID(dofnode[1])];
			currpos(2) = node->X()[2] + discol[discret_.DofColMap()->LID(dofnode[2])];
			//if node has also rotational degrees of freedom
			if (discret_.NumDof(node) == 6)
			{
				currrot(0) = discol[discret_.DofColMap()->LID(dofnode[3])];
				currrot(1) = discol[discret_.DofColMap()->LID(dofnode[4])];
				currrot(2) = discol[discret_.DofColMap()->LID(dofnode[5])];
			}

			/*/ debugging cout
			for(int j=0; j<(int)currpos.M(); j++)
				if(currpos(j)<0.0 || currpos(j)>statmechparams_.get<double>("PeriodLength", 0.0))
					cout<<"Proc "<<discret_.Comm().MyPID()<<": currpos(discol["<<i<<"] = "<<currpos<<endl;*/

			currentpositions[node->LID()] = currpos;
			currentrotations[node->LID()] = currrot;
		}

		//new search for neighbour nodes after average time specified in input file
		RCP<Epetra_MultiVector> crosslinkerneighbours;
		// in case of deactivated crosslink molecule diffusion (i.e. crosslinking is a purely geometric matter: node-to-node-distance of filaments)
		if (!(Teuchos::getIntegralValue<int>(statmechparams_, "CRSLNKDIFFUSION")))
		{
			SearchNeighbours(currentpositions, crosslinkerneighbours);
			cout << "\n***\nsearch time: " << Teuchos::Time::wallTime() - t_search
					<< " seconds\n***\n";
		}

#ifdef MEASURETIME
		const double t_admin = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME
		//number of elements in this time step before adding or deleting any elements
		currentelements_ = discret_.NumGlobalElements();
		//set crosslinkers the "old way", i.e. geometric proximity check of filament nodes
		if (!(Teuchos::getIntegralValue<int>(statmechparams_, "CRSLNKDIFFUSION")))
		{
			SetCrosslinkers(dt, noderowmap, nodecolmap, currentpositions,currentrotations, *crosslinkerneighbours);

			//delete crosslinkers
			DelCrosslinkers(dt, noderowmap, nodecolmap);
		}
		// set crosslinkers the "new way", i.e. considering crosslink molecule diffusion
		else
		{
			SearchAndSetCrosslinkers(istep, dt, noderowmap, nodecolmap, currentpositions,currentrotations);
			// filled in because #DEBUG wanted it
			discret_.CheckFilledGlobally();
			discret_.FillComplete(true, false, false);
			SearchAndDeleteCrosslinkers(dt, noderowmap, nodecolmap, currentpositions);
			discret_.CheckFilledGlobally();
			discret_.FillComplete(true, false, false);
			if (Teuchos::getIntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT") == INPAR::STATMECH::statout_gmsh &&
					istep % statmechparams_.get<int>("OUTPUTINTERVALS",1)==0)
					GmshPrepareVisualization(disrow);
		}

		/*settling administrative stuff in order to make the discretization ready for the next time step: synchronize
		 *the Filled() state on all processors after having added or deleted elements by ChekcFilledGlobally(); then build
		 *new element maps and call FillComplete(); finally Crs matrices stiff_ has to be deleted completely and made ready
		 *for new assembly since their graph was changed*/
		discret_.CheckFilledGlobally();
		discret_.FillComplete(true, false, false);
		stiff->Reset();

#ifdef MEASURETIME
		cout << "\n***\nadministration time: " << Teuchos::Time::wallTime() - t_admin<< " seconds\n***\n";
#endif // #ifdef MEASURETIME

	}//if(Teuchos::getIntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))

#ifdef MEASURETIME
	const double Delta_t = Teuchos::Time::wallTime()-t_start;
	cout << "\n***\ntotal time: " << Delta_t<< " seconds\n***\n";
#endif // #ifdef MEASURETIME

	// test cout
	/*for(int i=0; i<discret_.NumMyColNodes(); i++)
	{
		for(int proc=0; proc<discret_.Comm().NumProc(); proc++)
		{
			if(proc==discret_.Comm().MyPID())
			{
				DRT::Node *node = discret_.lColNode(i);
				if(node->NumElement()>2+statmechparams_.get<int>("N_CROSSMAX",0))
				{
					cout<<"Proc "<<discret_.Comm().MyPID()<<": nodeLID = "<<i<<", "<<node->NumElement()<<" elements: ";
					for(int j=0; j<node->NumElement(); j++)
						cout<<node->Elements()[j]->Id()<<" ";
					cout<<endl;
				}
			}
			discret_.Comm().Barrier();
		}
	}*/
	return;
} // StatMechManager::StatMechUpdate()

/*----------------------------------------------------------------------*
 | Shifts current position of nodes so that they comply with periodic   |
 | boundary conditions                                       cyron 04/10|
 *----------------------------------------------------------------------*/
void StatMechManager::PeriodicBoundaryShift(Epetra_Vector& disrow, int ndim, const double &dt)
{
  double starttime = statmechparams_.get<double>("STARTTIME",0.0);
  //period length of simulated box
  double H = statmechparams_.get<double> ("PeriodLength", 0.0);

	//only if period length >0 has been defined periodic boundary conditions are swithced on
	if (H > 0.0)
		for (int i=0; i<discret_.NumMyRowNodes(); i++)
		{
			//get a pointer at i-th row node
			const DRT::Node* node = discret_.lRowNode(i);

			//get GIDs of this node's degrees of freedom
			std::vector<int> dofnode = discret_.Dof(node);

			for (int j=ndim-1; j>-1; j--)
			{
				//current coordinate value
				double xcurr = node->X()[j] + disrow[discret_.DofRowMap()->LID(dofnode[j])];

				/*if node currently has coordinate value greater than statmechparams_.get<double>("PeriodLength",0.0),
				 *it is shifted by -statmechparams_.get<double>("PeriodLength",0.0) sufficiently often to lie again in the domain*/
				if (xcurr > H)
				{
					disrow[discret_.DofRowMap()->LID(dofnode[j])] -= H*floor(xcurr/H);

					/*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
					 *may fixed by DBC. To avoid problems when nodes exit the domain through the upper z-surface and reenter through the lower
					 *z-surface, the shear has to be substracted from nodal coordinates in that case */
					if (j == 2 && statmechparams_.get<int> ("CURVENUMBER", -1) >= 1 && time_ > starttime && fabs(time_-starttime)>dt/1e4)
						disrow[discret_.DofRowMap()->LID(dofnode[statmechparams_.get<int> ("OSCILLDIR", -1)])] -= statmechparams_.get<double> ("SHEARAMPLITUDE", 0.0) * DRT::Problem::Instance()->Curve(statmechparams_.get<int> ("CURVENUMBER", -1) - 1).f(time_);
				}
				/*if node currently has coordinate value smaller than zero, it is shifted by statmechparams_.get<double>("PeriodLength",0.0) sufficiently often
				 *to lie again in the domain*/
				if (xcurr < 0.0)
				{
					disrow[discret_.DofRowMap()->LID(dofnode[j])] -= H*floor(xcurr/H);

					/*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
					 *may be fixed by DBC. To avoid problems when nodes exit the domain through the lower z-surface and reenter through the upper
					 *z-surface, the shear has to be added to nodal coordinates in that case */
					if (j == 2 && statmechparams_.get<int> ("CURVENUMBER", -1) >= 1 && time_ > starttime && fabs(time_-starttime)>dt/1e4)
						disrow[discret_.DofRowMap()->LID(dofnode[statmechparams_.get<int> ("OSCILLDIR", -1)])] += statmechparams_.get<double> ("SHEARAMPLITUDE", 0.0) * DRT::Problem::Instance()->Curve(statmechparams_.get<int> ("CURVENUMBER", -1) - 1).f(time_);
				}
			}
		}

	return;
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

	DRT::ELEMENTS::Beam3* beam = dynamic_cast<DRT::ELEMENTS::Beam3*> (element);

	//3D beam elements are embeddet into R^3:
	const int ndim = 3;

	/*get reference configuration of beam3 element in proper format for later call of SetUpReferenceGeometry;
	 * note that rotrefe for beam3 elements is related to the entry in the global total Lagrange displacement
	 * vector related to a certain rotational degree of freedom; as the displacement is initially zero also
	 * rotrefe is set to zero here*/
	vector<double> xrefe(beam->NumNode() * ndim, 0);
	vector<double> rotrefe(beam->NumNode() * ndim, 0);

	for (int i=0; i<beam->NumNode(); i++)
		for (int dof = 0; dof < ndim; dof++)
		{
			xrefe[3* i + dof] = beam->Nodes()[i]->X()[dof];
			rotrefe[3* i + dof] = 0.0;
		}

	/*loop through all nodes except for the first node which remains fixed as reference node; all other nodes are
	 * shifted due to periodic boundary conditions if required*/
	for (int i=1; i<beam->NumNode(); i++)
	{
		for (int dof = 0; dof < ndim; dof++)
		{
			/*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
			 * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
			 * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
			 * is smaller than half the periodic length*/
			if (fabs((beam->Nodes()[i]->X()[dof]) + statmechparams_.get<double> ("PeriodLength", 0.0) - (beam->Nodes()[0]->X()[dof])) < fabs((beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof])))
				xrefe[3* i + dof] += statmechparams_.get<double> ("PeriodLength", 0.0);

			if (fabs((beam->Nodes()[i]->X()[dof]) - statmechparams_.get<double> ("PeriodLength", 0.0) - (beam->Nodes()[0]->X()[dof])) < fabs((beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof])))
				xrefe[3* i + dof] -= statmechparams_.get<double> ("PeriodLength", 0.0);
		}
	}

	/*SetUpReferenceGeometry is a templated function; note that the third argument "true" is necessary as all beam elements
	 * have already been initialized once upon reading input file*/
	switch (beam->NumNode())
	{
		case 2:
		{
			beam->SetUpReferenceGeometry<2> (xrefe, rotrefe, true);
			break;
		}
		case 3:
		{
			beam->SetUpReferenceGeometry<3> (xrefe, rotrefe, true);
			break;
		}
		case 4:
		{
			beam->SetUpReferenceGeometry<4> (xrefe, rotrefe, true);
			break;
		}
		case 5:
		{
			beam->SetUpReferenceGeometry<5> (xrefe, rotrefe, true);
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

	DRT::ELEMENTS::Truss3* truss = dynamic_cast<DRT::ELEMENTS::Truss3*> (element);

	//3D beam elements are embeddet into R^3:
	const int ndim = 3;

	/*get reference configuration of truss3 element in proper format for later call of SetUpReferenceGeometry*/
	vector<double> xrefe(truss->NumNode() * ndim, 0);

	for (int i=0; i<truss->NumNode(); i++)
		for (int dof = 0; dof < ndim; dof++)
			xrefe[3* i + dof] = truss->Nodes()[i]->X()[dof];

	/*loop through all nodes except for the first node which remains fixed as reference node; all other nodes are
	 * shifted due to periodic boundary conditions if required*/
	for (int i=1; i<truss->NumNode(); i++)
	{
		for (int dof = 0; dof < ndim; dof++)
		{
			/*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
			 * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
			 * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
			 * is smaller than half the periodic length*/
			if (fabs((truss->Nodes()[i]->X()[dof]) + statmechparams_.get<double> ("PeriodLength", 0.0) - (truss->Nodes()[0]->X()[dof])) < fabs((truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof])))
				xrefe[3* i + dof] += statmechparams_.get<double> ("PeriodLength", 0.0);

			if (fabs((truss->Nodes()[i]->X()[dof]) - statmechparams_.get<double> ("PeriodLength", 0.0) - (truss->Nodes()[0]->X()[dof])) < fabs((truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof])))
				xrefe[3* i + dof] -= statmechparams_.get<double> ("PeriodLength", 0.0);
		}
	}

	/*note that the third argument "true" is necessary as all truss elements have already been initialized once upon reading input file*/
	truss->SetUpReferenceGeometry(xrefe, true);

#endif
}
/*----------------------------------------------------------------------*
 | Searches and saves in variable  crosslinkerneighbours neighbours for|
 | each col node (OLD)                                      cyron 07/09|
 *----------------------------------------------------------------------*/
void StatMechManager::SearchNeighbours(const std::map<int, LINALG::Matrix<3, 1> > currentpositions, RCP<Epetra_MultiVector>& crosslinkerneighbours)
{
	//update the number of times the function SearchNeighbours has already been called
	nsearch_++;

	//variable for local storage of neighbours of each nodes (initialized with with empty vector for each row node)
	std::vector<std::vector<int> > crosslinkerneighboursloc;
	std::vector<int> emptyvector;
	crosslinkerneighboursloc.resize(discret_.NumMyRowNodes(), emptyvector);

	//mean distance bridged by a crosslinker
	double rlink = statmechparams_.get<double> ("R_LINK", 0.0);
	//absolute value of difference between maximal/minimal and mean distance bridged by crosslinkers
	double deltarlink = statmechparams_.get<double> ("DeltaR_LINK", 0.0);

	//maximal number of numbers of any nodes locally on this processor and globally
	int maxneighbourslocal = 0;
	int maxneighboursglobal = 0;

	//each processor looks for each of its row nodes for neighbours; loop index i is the local row node Id
	for (int i=0; i<discret_.NumMyRowNodes(); i++)
	{
		//get GID and column map LID of row node i
		int iGID = (discret_.lRowNode(i))->Id();
		int icolID = (discret_.gNode(iGID))->LID();

		/*getting iterator to current position of local node icolID from map currentpositions;
		 * note that currentpositions is organized with respect to column map LIDs*/
		const map<int, LINALG::Matrix<3, 1> >::const_iterator posi = currentpositions.find(icolID);

		/* for each row node all comlumn nodes within a range rlink are searched; this search may
		 * be performed either as a brute force search (A) or by means of a search tree (B); the
		 * column map LIDs of the nodes within rlink are stored in the vector neighboursLID*/
		std::vector<int> neighboursLID;

		//(A): brute force search
		for (std::map<int, LINALG::Matrix<3, 1> >::const_iterator posj = currentpositions.begin(); posj != currentpositions.end(); posj++)
		{
			//difference vector between row node i and some column node j
			LINALG::Matrix<3, 1> difference;
			for (int k=0; k<3; k++)
				difference(k) = (posi->second)(k) - (posj->second)(k);

			if (difference.Norm2() < rlink + deltarlink && difference.Norm2() > rlink - deltarlink)
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
		maxneighbourslocal = max(maxneighbourslocal,(int) crosslinkerneighboursloc[i].size());
	}

	//get global maximal number of neighbours of any node
	discret_.Comm().MaxAll(&maxneighbourslocal, &maxneighboursglobal, 1);

	//write information in crosslinkerneighboursloc into an Epetra_MultiVector
	Epetra_MultiVector crosslinkerneighboursrow(*(discret_.NodeRowMap()),maxneighboursglobal);
	crosslinkerneighboursrow.PutScalar(-1);
	for (int i=0; i<discret_.NumMyRowNodes(); i++)
		for (int j=0; j<(int) (crosslinkerneighboursloc[i]).size(); j++)
			crosslinkerneighboursrow[j][i] = (crosslinkerneighboursloc[i])[j];

	//export information in crosslinkerneighboursrow to col map variable crosslinkerneighbours
	crosslinkerneighbours = rcp(new Epetra_MultiVector(*(discret_.NodeColMap()),maxneighboursglobal));

	Epetra_Import importer(*(discret_.NodeColMap()), *(discret_.NodeRowMap()));
	crosslinkerneighbours->Import(crosslinkerneighboursrow, importer, Add);

}//void SearchNeighbours(const int rlink, const std::map<int,LINALG::Matrix<3,1> > currentpositions)

/*----------------------------------------------------------------------*
 | Assign crosslink molecules and nodes to volume partitions            |
 |																								(public) mueller 08/10|
 *----------------------------------------------------------------------*/
void StatMechManager::PartitioningAndSearch(const std::map<int,LINALG::Matrix<3,1> >& currentpositions, RCP<Epetra_MultiVector>& neighbourslid)
{
	//filament nodes and crosslink molecules are indexed according to their positions within the boundary box volume
	std::vector<std::vector<std::vector<int> > > nodeinpartition(3, std::vector<std::vector<int> >(statmechparams_.get<int>("SEARCHRES",1), std::vector<int>()));

	// initialize vectors related to volume indexing
	if(statmechparams_.get<int>("SEARCHRES",1)<1)
		dserror("Please give a plausible value for SEARCHRES!");

	double pl = statmechparams_.get<double>("PeriodLength",0.0);
	int N = statmechparams_.get<int>("SEARCHRES", 1);

	/*nodes*/
	//cout<<"Proc "<<discret_.Comm().MyPID()<<" : partitioning of nodes"<<endl;
	// loop over node positions to map their column map LIDs to partitions
	for (std::map<int, LINALG::Matrix<3, 1> >::const_iterator posi = currentpositions.begin(); posi != currentpositions.end(); posi++)
		for(int j=0; j<(int)nodeinpartition.size(); j++) // nodeinpartition.size==3
		{
			int partition = (int)std::floor((posi->second)(j)/pl*(double)N);
			if(partition==N)
				partition--;
			if(partition<0 || partition>N-1)
				cout<<"Proc "<<discret_.Comm().MyPID()<<": pos("<<j<<") = "<<(posi->second)(j)<<", partition = "<<partition<<endl;
			nodeinpartition[j][partition].push_back((int)(posi->first)); //column lid
		}

	/*crosslink molecules*/
	//cout<<"Proc "<<discret_.Comm().MyPID()<<" : Epetra_Export"<<endl;
	// Export crosslinkerpositions_ to transfermap_ format (kind of a row map format for crosslink molecules)
	Epetra_MultiVector crosslinkerpositionstrans(*transfermap_, 3, true);
	Epetra_Vector numbondtrans(*transfermap_, true);
	Epetra_MultiVector crosslinkpartitiontrans(*transfermap_, 3, false);
	Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);

	// preparations
	if(discret_.Comm().MyPID()!=0)
	{
	  numbond_->PutScalar(0.0);
		crosslinkerpositions_->PutScalar(0.0);
	}
	// Export to transfer map format
	numbondtrans.Export(*numbond_, crosslinkexporter, Add);
	crosslinkerpositionstrans.Export(*crosslinkerpositions_, crosslinkexporter, Add);
	//cout<<"Proc "<<discret_.Comm().MyPID()<<" : partitioning of crosslink molecules"<<endl;
	for(int i=0; i<crosslinkpartitiontrans.MyLength(); i++)
	{
	  // mark entries with double-bonded crosslink molecules
	  if(numbondtrans[i]>1.9)
	  {
	    for(int j=0; j<crosslinkpartitiontrans.NumVectors(); j++)
	      crosslinkpartitiontrans[j][i] = -1.0;
	    continue;
	  }
	  else
	  {
			for(int j=0; j<crosslinkpartitiontrans.NumVectors(); j++)
			{
				int partition = (int)std::floor(crosslinkerpositionstrans[j][i]/pl*(double)N);
				if(partition==N)
					partition--;
				crosslinkpartitiontrans[j][i] = partition;
			}
	  }
	}
	//cout<<"Proc "<<discret_.Comm().MyPID()<<" : DetectNeighbourNodes()...";
	// detection of nodes within search proximity of the crosslink molecules
	DetectNeighbourNodes(currentpositions, &nodeinpartition, numbondtrans, crosslinkerpositionstrans, crosslinkpartitiontrans, neighbourslid);
	//cout<<"done!"<<endl;
	return;
}//void StatMechManager::PartitioningAndSearch

/*----------------------------------------------------------------------*
 | detect neighbour nodes to crosslink molecules (public) mueller (7/10)|
 *----------------------------------------------------------------------*/
void StatMechManager::DetectNeighbourNodes(const std::map<int,LINALG::Matrix<3,1> >& currentpositions,
																					 std::vector<std::vector<std::vector<int> > >* nodeinpartition,
																					 Epetra_Vector& numbond,
																					 Epetra_MultiVector& crosslinkerpositions,
																					 Epetra_MultiVector& crosslinkpartitions,
																					 RCP<Epetra_MultiVector>& neighbourslid)
{
	/* Description:
	 * A vector containing the volume partitions of the crosslink molecules is handed over to this method. We loop over the partition
	 * number of the first (x) component and its two neighbouring partition numbers, hence gathering information on all nodes
	 * within these three layers.
	 *
	 * Now, we check, whether or not an LID of the first component matches one of the LIDs in the second component's partition layers
	 * (once again, we check the given partition number and its immediate neighbours).
	 * If so, we head to the third component and repeat this procedure.
	 * If a match is found for all components, we can be sure that the crosslink molecule in question lies within
	 * the 27 partitions encompassing the crosslink molecule partition.
	 *
	 * Eventually, the distance between the crosslink molecule and the filament node is calculated.
	 * If the distance lies beneath the boundaries given by the crosslinker lengths rmin and rmax,
	 * the node's LID is added to a storage vector for further use.
	 *
	 * After having found a match in the next component, we exit the loop to avoid unnecessary computational cost
	 */

	std::vector<std::vector<int> > neighbournodes(crosslinkpartitions.MyLength(), std::vector<int>());

	int maxneighbourslocal = 0;
	int maxneighboursglobal = 0;

	for(int part=0; part<crosslinkpartitions.MyLength(); part++)
	{
		if(crosslinkpartitions[0][part]<-0.9)
			continue;
		// determine search radius in accordance to bonding status
		double rmin, rmax;
		if ((int)numbond[part]<0.1)
		{
			rmin = (statmechparams_.get<double>("R_LINK", 0.0)-statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
			rmax = (statmechparams_.get<double>("R_LINK", 0.0)+statmechparams_.get<double>("DeltaR_LINK", 0.0)) / 2.0;
		}
		else
		{
			rmin = statmechparams_.get<double>("R_LINK", 0.0)-statmechparams_.get<double>("DeltaR_LINK", 0.0);
			rmax = statmechparams_.get<double>("R_LINK", 0.0)+statmechparams_.get<double>("DeltaR_LINK", 0.0);
		}
		// first component
		for(int ilayer=(int)crosslinkpartitions[0][part]-1; ilayer<(int)crosslinkpartitions[0][part]+2; ilayer++)
			if(ilayer>-1 && ilayer<statmechparams_.get<int>("SEARCHRES",1))
				for(int i=0; i<(int)(*nodeinpartition)[0][ilayer].size(); i++)
				{
					int tmplid = (int)(*nodeinpartition)[0][ilayer][i];
					// second component
					for(int jlayer=(int)crosslinkpartitions[1][part]-1; jlayer<(int)crosslinkpartitions[1][part]+2; jlayer++)
						if(jlayer>-1 && jlayer<statmechparams_.get<int>("SEARCHRES",1))
							for(int j=0; j<(int)(*nodeinpartition)[1][jlayer].size(); j++)
								if((*nodeinpartition)[1][jlayer][j]==tmplid)
								{
									//third component
									for(int klayer=(int)crosslinkpartitions[2][part]-1; klayer<(int)crosslinkpartitions[2][part]+2; klayer++)
										if(klayer>-1 && klayer<statmechparams_.get<int>("SEARCHRES",1))
											for(int k=0; k<(int)(*nodeinpartition)[2][klayer].size(); k++)
												if((*nodeinpartition)[2][klayer][k]==tmplid)
												{
													// get the current node position for the node with LID==tmplid
													const map<int, LINALG::Matrix<3, 1> >::const_iterator nodepos = currentpositions.find(tmplid);
													// calculate distance crosslinker-node
													LINALG::Matrix<3, 1> difference;
													for (int l=0; l<(int)difference.M(); l++)
														difference(l) = crosslinkerpositions[l][part]-(nodepos->second)(l);
													// only nodes within the search volume are stored
													if(difference.Norm2()<rmax && difference.Norm2()>rmin)
														neighbournodes[part].push_back(tmplid);
													// exit loop immediately
													break;
												}
									break;
								}
				}
		// "-1" indicates the possibility of a crosslink molecule becoming passive, i.e. hypothetically bonding to the same filament
		if((int)numbond[part]==1)
			neighbournodes[part].push_back(-1);
		// store local maximal number of LIDs per molecule in order to determine neighbourslid->NumVectors()
		maxneighbourslocal = max(maxneighbourslocal, (int)neighbournodes[part].size());
	}

	// get global maximal number of LIDs per molecule
	discret_.Comm().MaxAll(&maxneighbourslocal, &maxneighboursglobal, 1);
	if(maxneighboursglobal==0)
		maxneighboursglobal = 1;
	// copy information to Epetra_MultiVector for communication
	RCP<Epetra_MultiVector> neighbourslidtrans = rcp(new Epetra_MultiVector(*transfermap_, maxneighboursglobal, false));

	Epetra_Import crosslinkimporter(*crosslinkermap_, *transfermap_);
	/* assign "-2" (also to distinguish between 'empty' and passive crosslink molecule "-1"
	 * to be able to determine entries which remain "empty" due to number of LIDs < maxneighboursglobal*/
	neighbourslidtrans->PutScalar(-2.0);
	for(int i=0; i<neighbourslidtrans->MyLength(); i++)
		for(int j=0; j<(int)neighbournodes[i].size(); j++)
			(*neighbourslidtrans)[j][i] = (double)neighbournodes[i][j];

	// make information redundant on al Procs
	neighbourslid = rcp(new Epetra_MultiVector(*crosslinkermap_,maxneighboursglobal,false));
	neighbourslid->Import(*neighbourslidtrans, crosslinkimporter, Insert);

	return;
}// StatMechManager::DetectNeighbourNodes

/*----------------------------------------------------------------------*
 | Searches for crosslink molecule-filament node pairs and adds actual  |
 | crosslinker elements once linking conditions are met.       					|
 | (public)                                               mueller (7/10)|
 *----------------------------------------------------------------------*/
void StatMechManager::SearchAndSetCrosslinkers(const int& istep,const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
																							 const std::map<int, LINALG::Matrix<3, 1> >& currentpositions,
																							 const std::map<int,LINALG::Matrix<3, 1> >& currentrotations)
{
#ifdef D_TRUSS3
#ifdef D_BEAM3
#ifdef D_BEAM3II

	// update the number of times the function SearchNeighbours has already been called
	nsearch_++;

	double t_search = Teuchos::Time::wallTime();
	/*preliminaries*/

	// NODAL TRIAD UPDATE
	//first get triads at all row nodes
	Epetra_MultiVector nodaltriadsrow(noderowmap, 4);
	nodaltriadsrow.PutScalar(0);

	//update nodaltriads_
	for (int i=0; i<noderowmap.NumMyElements(); i++)
	{
		//lowest GID of any connected element (the related element cannot be a crosslinker, but has to belong to the actual filament discretization)
		int lowestid(((discret_.lRowNode(i)->Elements())[0])->Id());
		int lowestidele(0);
		for (int j=0; j<discret_.lRowNode(i)->NumElement(); j++)
			if (((discret_.lRowNode(i)->Elements())[j])->Id() < lowestid)
			{
				lowestid = ((discret_.lRowNode(i)->Elements())[j])->Id();
				lowestidele = j;
			}

		//check whether filaments are discretized with beam3ii elements (otherwise no orientation triads at nodes available)
		DRT::ELEMENTS::Beam3ii* filele = NULL;
		DRT::ElementType & eot = ((discret_.lRowNode(i)->Elements())[lowestidele])->ElementType();
		if (eot == DRT::ELEMENTS::Beam3iiType::Instance())
			filele = dynamic_cast<DRT::ELEMENTS::Beam3ii*> (discret_.lRowNode(i)->Elements()[lowestidele]);
		else
			dserror("Filaments have to be discretized with beam3ii elements for orientation check!!!");

		//check whether crosslinker is connected to first or second node of that element
		int nodenumber = 0;
		if(discret_.lRowNode(i)->Id() == ((filele->Nodes())[1])->Id() )
			nodenumber = 1;

		//save nodal triad of this node in nodaltriadrow
		for(int j=0; j<4; j++)
			nodaltriadsrow[j][i] = (filele->Qnew_[nodenumber])(j);
	}

	//export nodaltriadsrow to col map variable
	Epetra_MultiVector nodaltriadscol(nodecolmap,4,true);
	Epetra_Import importer(nodecolmap,noderowmap);
	nodaltriadscol.Import(nodaltriadsrow,importer,Insert); // NODAL TRIAD UPDATE

	//get current on-rate for crosslinkers
	double kon = 0;
	double starttime = statmechparams_.get<double>("STARTTIME", 0.0);

	if(time_ <= starttime || (time_>starttime && fabs(time_-starttime) < dt/1e4))
		kon = statmechparams_.get<double>("K_ON_start",0.0);
	else
		kon = statmechparams_.get<double>("K_ON_end",0.0);

	//probability with which a crosslinker is established between crosslink molecule and neighbour node
	double plink = 1.0 - exp( -dt*kon*statmechparams_.get<double>("C_CROSSLINKER",0.0) );
	//probability with which a crosslinker blocks two binding spots on the same filament (self-binding)
	double pself = 1.0 - exp( -dt*statmechparams_.get<double>("K_ON_SELF", 0.0) );

	//Volume partitioning, assignment of nodes and molecules to partitions, search for neighbours
	// map filament (column map, i.e. entire node set on each proc) node positions to volume partitions every SEARCHINTERVAL timesteps
	RCP<Epetra_MultiVector> neighbourslid;
	if(statmechparams_.get<int>("SEARCHRES",1)>0)
		PartitioningAndSearch(currentpositions, neighbourslid);
	//cout<<*neighbourslid<<endl;
	/*the following part of the code is executed on processor 0 only and leads to the decision, at which nodes crosslinks shall be set
	 *processor 0 goes through all the crosslink molecules and checks whether a crosslink is to be set; this works precisely as follows:
	 *(1) the crosslink molecules are looped through in a random order
	 *(2) if a node has not yet reached its maximal number of crosslinks, a crosslink may be set
	 *(3) a crosslink is set if and only if the node passes a probability check
	 *note: a fully overlapping node map on processor 0 is imperative for the following part of the code to work correctly!*/

	// a vector indicating the crosslink molecule which is going to constitute a crosslinker element
	Epetra_Vector doublebond(*crosslinkermap_, true);

	if(discret_.Comm().MyPID()==0)
	{
		//cout<<*neighbourslid<<endl;
		// obtain a random order in which the crosslinkers are addressed
		std::vector<int> order = Permutation(statmechparams_.get<int>("N_crosslink", 0));

		for(int i=0; i<numbond_->MyLength(); i++)
		{
			int irandom = order[i];

			// skip this crosslink molecule if it already bonded with two nodes or if it is passive
			if((*numbond_)[irandom]>1.9)
				continue;

			// PROBABILITY CHECK & PREPARATION OF VECTORS
			// obtain a random order of neighboursLID indices
			std::vector<int> neighbourorder = Permutation(neighbourslid->NumVectors());

			// loop over neighbour nodes
			for(int j=0; j<neighbourslid->NumVectors(); j++)
			{
				// random index
				int index = neighbourorder[j];
				// continue, if neighbourslid entry is '-2', meaning empty
				if((*neighbourslid)[index][irandom] < -1.9)
					continue;

				// current neighbour LID
				int nodeLID = (int)(*neighbourslid)[index][irandom];
				//continue in case of N_CROSSMAX crosslinkers at the current node
				if((*numcrossnodes_)[nodeLID]>=statmechparams_.get<int>("N_CROSSMAX",0))
					continue;

				// flag indicating loop break after first new bond has been established between i-th crosslink molecule and j-th neighbour node
				bool bondestablished = false;
				// necessary condition to be fulfilled in order to set a crosslinker
				double probability = 0.0;
				// switch between probability for regular inter-filament binding and self-binding crosslinker
				if(nodeLID>=0)
					probability = plink;
				else
					probability = pself;

				if(uniformclosedgen_.random() < probability)
				{
					int free = 0;
					int occupied = 0;
					bool set = false;

					// check both entries of crosslinkerbond_
					for(int k=0; k<crosslinkerbond_->NumVectors(); k++)
						if((*crosslinkerbond_)[k][irandom]<-0.9 && !set)
						{
							// store free bond position
							free = k;
							set = true;
						}
						else
							occupied = k;

					switch((int)(*numbond_)[irandom])
					{
						// free crosslink molecule creating a bond
						case 0:
						{
							// update of crosslink molecule positions
							LINALG::SerialDenseMatrix LID(1,1,true);
							LID(0,0) = nodeLID;
							(*crosslinkerbond_)[free][irandom] = nodecolmap.GID(nodeLID);
							// increment the number of crosslinkers at this node
							((*numcrossnodes_)[nodeLID]) += 1.0;
							// increment the number of bonds of this crosslinker
							((*numbond_)[irandom]) = 1.0;
							CrosslinkerIntermediateUpdate(currentpositions, LID, irandom);
							bondestablished = true;
						}
						break;
						// one bond already exists -> establish a second bond/passive crosslink molecule
						// Potentially, add an element (later)
						case 1:
						{
							//Col map LIDs of nodes to be crosslinked
							LINALG::Matrix<2,1> LID;
							LID(0) = nodecolmap.LID((int)(*crosslinkerbond_)[occupied][irandom]);
							// distinguish between a real nodeLID and the entry '-1', which indicates a potentially passive crosslink molecule
							if(nodeLID>=0)
								LID(1) = nodeLID;
							else // choose the neighbours node on the same filament as nodeLID as second entry and take basisnodes_ into account
							{
								int currfilament = (int)(*filamentnumber_)[(int)LID(0)];
								if((int)LID(0)<basisnodes_-1)
									if((int)(*filamentnumber_)[(int)LID(0)+1]==currfilament)
										LID(1) = LID(0) + 1.0;
									else
										LID(1) = LID(0) - 1.0;
								if((int)LID(0)==basisnodes_-1)
									if((int)(*filamentnumber_)[(int)LID(0)-1]==currfilament)
										LID(1) = LID(0) - 1.0;
							}

							//unit direction vector between currently considered two nodes
							LINALG::Matrix<3,1> direction((currentpositions.find((int)LID(0)))->second);
							direction -= (currentpositions.find((int)LID(1)))->second;
							direction.Scale(1.0/direction.Norm2());

							/*Orientation Check: only if the two nodes in question pass a check for their
							 * orientation, a marker, indicating an element to be added, will be set*/
							if(CheckOrientation(direction,nodaltriadscol,LID))
							{
								LINALG::SerialDenseMatrix lid(2,1,true);
								lid(0,0) = LID(0);
								lid(1,0) = LID(1);

								((*numbond_)[irandom]) = 2.0;
								// actually existing two bonds
								if(nodeLID>=0)
								{
									(*crosslinkerbond_)[free][irandom] = nodecolmap.GID(nodeLID);
									((*numcrossnodes_)[nodeLID]) += 1.0;
									// insert node GID in order to check the correct setup
									doublebond[irandom] = 1.0;
									// update molecule positions
									CrosslinkerIntermediateUpdate(currentpositions, lid, irandom);

									// consider crosslinkers covering two binding spots of one filament
									if((*filamentnumber_)[(int)LID(0)]==(*filamentnumber_)[nodeLID])
										(*crosslinkonsamefilament_)[irandom] = 1.0;
								}
								else // passive crosslink molecule
								{
									(*searchforneighbours_)[irandom] = 0.0;
									LINALG::SerialDenseMatrix oneLID(1,1,true);
									oneLID(0,0) = LID(0);
									CrosslinkerIntermediateUpdate(currentpositions, oneLID, irandom);
								}
								bondestablished = true;
							}
						}
						break;
					}// switch((int)(*numbond_)[irandom])

					// for now, break after a new bond was established, i.e crosslinker elements cannot be established starting from zero bonds
					if(bondestablished)
						break;
				}// if(uniformclosedgen_.random() < plink)
			}// for(int j=0; j<(int)neighboursLID.size(); j++)

		}// for(int i=0; i<numbond_->MyLength(); i++)
		cout << "\n***\nsearch time: " << Teuchos::Time::wallTime() - t_search<< " seconds\n***\n";
	}// if(discret_.Comm().MypPID==0)
	else
	{
		/* zerofy numcrossnodes at the beginning of each search except for Proc 0
		 * for subsequent export and reimport. This way, we guarantee redundant information
		 * on all processors.*/
		numcrossnodes_->PutScalar(0.0);
		crosslinkerpositions_->PutScalar(0.0);
		crosslinkerbond_->PutScalar(0.0);
		numbond_->PutScalar(0.0);
		doublebond.PutScalar(0.0);
	}

	//synchronize information about number of bonded filament nodes by exporting it to row map format and then reimporting it to column map format
	Epetra_Export exporter(nodecolmap,noderowmap);
	Epetra_Vector numcrossnodesrow(noderowmap,true);
	numcrossnodesrow.Export(*numcrossnodes_,exporter,Add);
	numcrossnodes_->Import(numcrossnodesrow,importer,Insert);

	// exporter and importer for crosslink molecules
	Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
	Epetra_Import crosslinkimporter(*crosslinkermap_, *transfermap_);

	// transfer vectors
	Epetra_MultiVector crosslinkerpositionstrans(*transfermap_,3,true);
	Epetra_MultiVector crosslinkerbondtrans(*transfermap_,2,true);
	Epetra_Vector numbondtrans(*transfermap_, true);
	Epetra_Vector doublebondtrans(*transfermap_, true);

	// exports and reimports
	crosslinkerpositionstrans.Export(*crosslinkerpositions_, crosslinkexporter, Add);
	crosslinkerpositions_->Import(crosslinkerpositionstrans, crosslinkimporter, Insert);
	crosslinkerbondtrans.Export(*crosslinkerbond_, crosslinkexporter, Add);
	crosslinkerbond_->Import(crosslinkerbondtrans, crosslinkimporter, Insert);
	numbondtrans.Export(*numbond_, crosslinkexporter, Add);
	numbond_->Import(numbondtrans, crosslinkimporter, Insert);
	doublebondtrans.Export(doublebond, crosslinkexporter, Add);
	doublebond.Import(doublebondtrans, crosslinkimporter, Insert);

	/*/ check numcrossnodes for correct values
	Epetra_Vector check(*crosslinkermap_, true);
	Epetra_Vector checktrans(*transfermap_, true);

	checktrans.Export(*numbond_,crosslinkexporter, Add);
	check.Import(checktrans, crosslinkimporter, Insert);

	if(discret_.Comm().MyPID()==1)
		for(int i=0; i<check.MyLength(); i++)
			if((int)check[i] % discret_.Comm().NumProc()!=0)
				dserror("wrong communication of numbond");*/

	// debug couts
	//cout<<*crosslinkerpositions_<<endl;
	//cout<<*crosslinkerbond_<<endl;
	//cout<<*numbond_<<endl;
	//cout<<doublebond<<endl;

	std::vector<int> crosslinkerids;
	// ADDING ELEMENTS
	for(int i=0; i<doublebond.MyLength(); i++)
	{
		if(doublebond[i]>0.9)
		{
			// obtain node GID
			int nodeGID[2] = {	(int)(*crosslinkerbond_)[0][i],(int)(*crosslinkerbond_)[1][i]};

			/*a crosslinker is added on each processor which is row node owner of at least one of its nodes;
			 *the processor which is row map owner of the node with the larger GID will be the owner of the new crosslinker element; on the other processors the new
			 *crosslinker element will be present as a ghost element only*/
			if(noderowmap.LID(nodeGID[0]) > -0.9 || noderowmap.LID(nodeGID[1]) > -0.9)
			{
				//getting current position and rotational status of nodes with GID nodeGID[]
				// node 1
				map< int,LINALG::Matrix<3,1> >::const_iterator pos0 = currentpositions.find( nodecolmap.LID(nodeGID[0]) );
				map< int,LINALG::Matrix<3,1> >::const_iterator rot0 = currentrotations.find( nodecolmap.LID(nodeGID[0]) );
				// node 2
				map< int,LINALG::Matrix<3,1> >::const_iterator pos1 = currentpositions.find( nodecolmap.LID(nodeGID[1]) );
				map< int,LINALG::Matrix<3,1> >::const_iterator rot1 = currentrotations.find( nodecolmap.LID(nodeGID[1]) );

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
				int GID2 = min(nodeGID[0],nodeGID[1]);
				int GID1 = max(nodeGID[0],nodeGID[1]);
				int globalnodeids[2] = {GID1, GID2};
				int newcrosslinkerGID = (GID1 + 1)*basisnodes_ + GID2;

				DRT::Node* nodes[2] = {	discret_.gNode( GID1 ) , discret_.gNode( GID2 )};

				//save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
				std::vector<double> rotrefe(6);
				std::vector<double> xrefe(6);

				for(int k=0; k<3; k++)
				{
					//the first three positions in xrefe and rotrefe are for the data of node GID1
					if(nodeGID[0] > nodeGID[1])
					{
						//set nodal positions
						xrefe[k ] = (pos0->second)(k);
						xrefe[k+3] = (pos1->second)(k);

						//set nodal rotations (not true ones, only those given in the displacement vector)
						rotrefe[k ] = (rot0->second)(k);
						rotrefe[k+3] = (rot1->second)(k);
					}
					else
					{
						//set nodal positions
						xrefe[k ] = (pos1->second)(k);
						xrefe[k+3] = (pos0->second)(k);

						//set nodal rotations (not true ones, only those given in the displacement vector)
						rotrefe[k ] = (rot1->second)(k);
						rotrefe[k+3] = (rot0->second)(k);
					}
				}
				if(statmechparams_.get<double>("ILINK",0.0) > 0.0)
				{
					crosslinkerids.push_back(newcrosslinkerGID);
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
					//cout<<"Proc "<<discret_.Comm().MyPID()<<": ADDED GID = "<<newcrosslinkerGID<<", numbond = "<<(*numbond_)[i]<<endl;
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
			}
		}
	} // ADDING ELEMENTS
	// couts
#endif //#ifdef D_TRUSS3
#endif //#ifdef D_BEAM3
#endif //#ifdef D_BEAM3II
}//void StatMechManager::SearchAndSetCrosslinkers

/*----------------------------------------------------------------------*
 | set crosslinkers between neighbouring nodes after probability check  |
 | (OLD)   (public)                                          cyron 02/09|
 *----------------------------------------------------------------------*/
void StatMechManager::SetCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
																			const std::map<int, LINALG::Matrix<3, 1> >& currentpositions,
																			const std::map<int, LINALG::Matrix<3, 1> >& currentrotations,
																			Epetra_MultiVector& crosslinkerneighbours)
{
#ifdef D_TRUSS3
#ifdef D_BEAM3
#ifdef D_BEAM3II

	//first get triads at all row nodes
	Epetra_MultiVector nodaltriadsrow(*(discret_.NodeRowMap()),4);
	nodaltriadsrow.PutScalar(0);

	//update nodaltriads_
	for(int i=0; i<(discret_.NodeRowMap())->NumMyElements(); i++)
	{

		//lowest GID of any connected element (the related element cannot be a crosslinker, but has to belong to the actual filament discretization)
		int lowestid( ( (discret_.lRowNode(i)->Elements())[0] )->Id() );
		// position of the element with the lowest GID within the list of all elements adjacent to row node i
		int lowestconnectedele(0);
		for(int j=0; j<discret_.lRowNode(i)->NumElement(); j++)
			if( ( (discret_.lRowNode(i)->Elements())[j] )->Id() < lowestid )
			{
				lowestid = ( (discret_.lRowNode(i)->Elements())[j] )->Id();
				lowestconnectedele = j;
			}

		//check whether filaments are discretized with beam3ii elements (otherwise no orientation triads at nodes available)
		DRT::ELEMENTS::Beam3ii* filele = NULL;
		DRT::ElementType & eot = ((discret_.lRowNode(i)->Elements())[lowestconnectedele])->ElementType();
		if(eot == DRT::ELEMENTS::Beam3iiType::Instance())
			filele = dynamic_cast<DRT::ELEMENTS::Beam3ii*>(discret_.lRowNode(i)->Elements()[lowestconnectedele]);
		else
			dserror("Filaments have to be discretized with beam3ii elements for orientation check!!!");

		//check whether crosslinker is connected to first or second node of that element
		int nodenumber = 0;
		if(discret_.lRowNode(i)->Id() == ((filele->Nodes())[1])->Id() )
			nodenumber = 1;

		//save nodal triad of this node in nodaltriadrow
		for(int j=0; j<4; j++)
			nodaltriadsrow[j][i] = (filele->Qnew_[nodenumber])(j);

	}

	//export nodaltriadsrow to col map variable
	Epetra_MultiVector nodaltriadscol(*(discret_.NodeColMap()),4,true);
	Epetra_Import importer(nodecolmap,noderowmap);
	nodaltriadscol.Import(nodaltriadsrow,importer,Insert);

	//get current on-rate for crosslinkers
	double kon = 0;
	if( (currentelements_ - basiselements_) < statmechparams_.get<int>("N_crosslink", 0) && !konswitch_)
		kon = statmechparams_.get<double>("K_ON_start",0.0);
	else
	{
		kon = statmechparams_.get<double>("K_ON_end",0.0);
		konswitch_=true;
	}

	//probability with which a crosslinker is established between neighbouring nodes
	double plink = 1.0 - exp( -dt*kon*statmechparams_.get<double>("C_CROSSLINKER",0.0) );

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

		for(int i=0; i< nodecolmap.NumMyElements(); i++)
		{
			int nodenumber = order[i];

			//go through all neighbors of this node in random order and look for suitable ones for crosslinking
			std::vector<int> neighboursorder(Permutation(crosslinkerneighbours.NumVectors()));
			int j = 0;

			//go through all neighbors (i.e. entries unequal -1) as long as current node has not yet reached maximal number of crosslinkers
			while(j < min(crosslinkerneighbours.NumVectors(),statmechparams_.get<int>("N_CROSSMAX",0)) && crosslinkerneighbours[neighboursorder[j]][nodenumber] != -1 && (*numcrosslinkerpartner_)[nodenumber] < statmechparams_.get<int>("N_CROSSMAX",0) )
			{
				//probability check whether at this node a crosslink can be set at all
				if(uniformclosedgen_.random() < plink)
				{

					//get LID of current neighbor node
					int neighbourLID = nodecolmap.LID((int)crosslinkerneighbours[neighboursorder[j]][nodenumber]);

					//unit direction vector between currently considered two nodes
					LINALG::Matrix<3,1> direction((currentpositions.find(neighbourLID))->second);
					direction -= (currentpositions.find(nodenumber))->second;
					direction.Scale(1.0 / direction.Norm2());

					//Col map LIDs of nodes to be crosslinked
					LINALG::Matrix<2,1> LID;
					LID(0) = nodenumber;
					LID(1) = neighbourLID;

					//a crosslinker can be set only if constraints wsith respect to its orientation are satisfied and also second node does not exceed it maximal number of crosslinkers
					if(CheckOrientation(direction,nodaltriadscol,LID) && (*numcrosslinkerpartner_)[neighbourLID] < statmechparams_.get<int>("N_CROSSMAX",0))
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
						for(int k=0; k<statmechparams_.get<int>("N_CROSSMAX",0); k++)
						if(nodecolmap.GID(neighbourLID) == (*crosslinkerpartner_)[k][nodenumber] )
							notyet = false;

						/*these two nodes are not yet crosslinked the link will actually be set and this is registered in crosslinkerpartner_ for
						 *the two linked nodes in the first currently unused multivector column, respectively*/
						bool registered1 = false;
						bool registered2 = false;
						if(notyet)
						for(int k=0; k<statmechparams_.get<int>("N_CROSSMAX",0); k++)
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
	for(int i=0; i<crosslinkstobeadded.MyLength(); i++)
	{
		//getting current position and rotational status of node with column map LID i from map currentpositions
		map< int,LINALG::Matrix<3,1> >::const_iterator posi = currentpositions.find(i);
		map< int,LINALG::Matrix<3,1> >::const_iterator roti = currentrotations.find(i);

		//loop through all node GIDs to which a crosslink should be established from column map node i
		for(int j=0; j<crosslinkstobeadded.NumVectors(); j++)
		{
			/*a crosslinker is to be established if the following two conditions are satisfied:
			 *(1) there is a node stored to which it is to be established (i.e. the i-th element of the j-th vector in in the Multi_Vector crosslinkerstobeaddedglobalcol is unequal the default value -1
			 *(2) second: one of the two nodes between which the crosslinker should be established is a row node on this processor; otherwise one does not need to add the crosslinker element on this processor*/
			if(crosslinkstobeadded[j][i] > -0.9 && ( noderowmap.LID((nodecolmap.GID(i))) > -0.9 || noderowmap.LID(((int)crosslinkstobeadded[j][i])) > -0.9 ) )
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
				int globalnodeids[2] = {	GID1,GID2};
				DRT::Node* nodes[2] = {	discret_.gNode( globalnodeids[0] ) , discret_.gNode( globalnodeids[1] )};

				//save positions of nodes between which a crosslinker has to be established in variables xrefe and rotrefe:
				std::vector<double> rotrefe(6);
				std::vector<double> xrefe(6);

				for(int k=0; k<3; k++)
				{
					//the first three positions in xrefe and rotrefe are for the data of node GID1
					if(nodecolmap.GID(i) > (int)crosslinkstobeadded[j][i])
					{
						//set nodal positions
						xrefe[k ] = (posi->second)(k);
						xrefe[k+3] = (posj->second)(k);

						//set nodal rotations (not true ones, only those given in the displacment vector)
						rotrefe[k ] = (roti->second)(k);
						rotrefe[k+3] = (rotj->second)(k);
					}
					else
					{
						//set nodal positions
						xrefe[k ] = (posj->second)(k);
						xrefe[k+3] = (posi->second)(k);

						//set nodal rotations (not true ones, only those given in the displacment vector)
						rotrefe[k ] = (rotj->second)(k);
						rotrefe[k+3] = (roti->second)(k);
					}
				}

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

			}
		}
	}

#endif  // #ifdef D_TRUSS3
#endif  // #ifdef D_BEAM3
#endif  // #ifdef D_BEAM3II
}//void SetCrosslinkers(const Epetra_Vector& setcrosslinkercol)

/*----------------------------------------------------------------------*
 | searches crosslinkers and deletes them if probability check is passed|
 | (private)																							mueller (7/10)|
 *----------------------------------------------------------------------*/
void StatMechManager::SearchAndDeleteCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,
																									const std::map<int, LINALG::Matrix<3, 1> >& currentpositions)
{
	/* synchronize information about adjacent elements and their respective node GIDs so that it can be used by Proc 0.
	 * This is currently somehow tedious and not too elegant (node map for dealing with elements, redundant element entries etc).*/
	Epetra_Import importer(nodecolmap, noderowmap);
	// vectors holding crosslinker element GIDs
	Epetra_MultiVector crosslinkergids(nodecolmap, statmechparams_.get<int>("N_CROSSMAX",0), true);
	Epetra_MultiVector crosslinkergidsrow(noderowmap, statmechparams_.get<int>("N_CROSSMAX",0), false);
	// node IDs of crosslinker elements; NumVectors()==number of unique node GIDs per filament node with attached crosslinker elements
	Epetra_MultiVector crosslinkernodeids(nodecolmap, statmechparams_.get<int>("N_CROSSMAX",0)+1, true);
	Epetra_MultiVector crosslinkernodeidsrow(noderowmap, statmechparams_.get<int>("N_CROSSMAX",0)+1, false);

	crosslinkergidsrow.PutScalar(-1.0);
	crosslinkernodeidsrow.PutScalar(-1.0);

	for(int i=0; i<discret_.NumMyRowNodes(); i++)
	{
		// insert crosslinker element GIDs
		int gidposition = 0;
		DRT::Node *node = discret_.lRowNode(i);
		for(int j=0; j<node->NumElement(); j++)
		{
			DRT::Element *element = node->Elements()[j];
			// element type either beam3 or truss3 (both with NumNode()==2)
			if(element->Id() > basisnodes_)
			{
				crosslinkergidsrow[gidposition][i] = node->Elements()[j]->Id();
				for(int k=0; k<element->NumNode(); k++)
					for(int l=0; l<crosslinkernodeidsrow.NumVectors(); l++)
						if(element->NodeIds()[k]!=crosslinkernodeidsrow[l][i]) // no identical entry yet
							crosslinkernodeidsrow[gidposition*element->NumNode()+k][i] = element->NodeIds()[k];
						else
							break;
				gidposition++;
			}
		}
	}
	// Import from row map to column map
	crosslinkergids.Import(crosslinkergidsrow, importer, Insert);
	crosslinkernodeids.Import(crosslinkernodeidsrow, importer, Insert);

	// a vector indicating the upcoming deletion of crosslinker elements
	Epetra_Vector delcrosslinkers(*crosslinkermap_, true);
	delcrosslinkers.PutScalar(-1.0);

	//get current off-rate for crosslinkers
	double koff = 0;
	double starttime = statmechparams_.get<double>("STARTTIME", 0.0);

	if (time_ <= starttime || (time_>starttime && fabs(time_-starttime)<dt/1e4))
		koff = statmechparams_.get<double> ("K_OFF_start", 0.0);
	else
		koff = statmechparams_.get<double> ("K_OFF_end", 0.0);

	//probability with which a crosslink breaks up in the current time step
	double punlink = 1.0 - exp(-dt * koff);

	// SEARCH
	// search and setup for the deletion of elements is done by Proc 0
	if (discret_.Comm().MyPID()==0)
	{
		for (int i=0; i<numbond_->MyLength(); i++)
		{
			// take action according to the number of bonds of a crosslink molecule
			switch ((int)(*numbond_)[i])
			{
				case 0:
					continue;
				break;
				// crosslink molecule with one bond
				case 1:
				{
					for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
						if ((*crosslinkerbond_)[j][i]>-0.9)
							if (uniformclosedgen_.random() < punlink)
							{
								// obtain LID and reset crosslinkerbond_ at this position
								int nodeLID = nodecolmap.LID((int) (*crosslinkerbond_)[j][i]);
								((*numcrossnodes_)[nodeLID]) -= 1.0;
								(*numbond_)[i] = 0.0;
								(*crosslinkerbond_)[j][i] = -1.0;

								LINALG::SerialDenseMatrix LID(1, 1, true);
								LID(0, 0) = nodeLID;
								CrosslinkerIntermediateUpdate(currentpositions, LID, i, false);
							}
				}
				break;
				// crosslinker element
				case 2:
				{
					// currently, if an element is deleted, a one-bonded crosslink molecule remains (no simultaneous cut of both bonds)
					for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
						if (uniformclosedgen_.random() < punlink)
						{
							(*numbond_)[i] = 1.0;

							// an actual crosslinker element exists
							if((*searchforneighbours_)[i]>0.9)
							{
								// get the nodal LID and then reset the entry
								int nodeLID = nodecolmap.LID((int)(*crosslinkerbond_)[j][i]);
								// get the crosslinker element GID and store it for later deletion:
								// k<2 maps to first crosslinker element GID, from then on k to the k-1-th crosslinker element (reason above: unique node GIDs)
								for(int k=0; k<crosslinkernodeids.NumVectors(); k++)
									if(crosslinkernodeids[k][nodeLID]==(int)(*crosslinkerbond_)[j][i])
										if(k<2)
											delcrosslinkers[i] = crosslinkergids[0][nodeLID];
										else
											delcrosslinkers[i] = crosslinkergids[k-1][nodeLID];

								((*numcrossnodes_)[nodeLID]) -= 1.0;
								(*crosslinkerbond_)[j][i] = -1.0;
								if((*crosslinkonsamefilament_)[i] > 0.9)
									(*crosslinkonsamefilament_)[i] = 0.0;
							}
							else	// passive crosslink molecule
								(*searchforneighbours_)[i] = 1.0;

							for (int k=0; k<crosslinkerbond_->NumVectors(); k++)
								if ((*crosslinkerbond_)[k][i]>-0.9)
								{
									LINALG::SerialDenseMatrix LID(1, 1, true);
									LID(0,0) = nodecolmap.LID((int)(*crosslinkerbond_)[k][i]);
									CrosslinkerIntermediateUpdate(currentpositions, LID, i);
									break;
								}
							break;
						}
				}
				break;
			}// switch ((int)(*numbond_)[i])
		}// for (int i=0; i<numbond_->MyLength(); i++)
	}// if(discret_.Comm().MyPID()==0)
	else
	{
		numcrossnodes_->PutScalar(0.0);
		// set crosslinkerbond_ to zero, but not to "-1" because of addition
		crosslinkerbond_->PutScalar(0.0);
		numbond_->PutScalar(0.0);
		delcrosslinkers.PutScalar(0.0);
	}

	//synchronize information about number of bonded filament nodes by exporting it to row map format and then reimporting it to column map format
	Epetra_Export exporter(nodecolmap, noderowmap);
	Epetra_Vector numcrossnodesrow(noderowmap, true);
	numcrossnodesrow.Export(*numcrossnodes_, exporter, Add);
	numcrossnodes_->Import(numcrossnodesrow, importer, Insert);

	// exporter and importer
	Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
	Epetra_Import crosslinkimporter(*crosslinkermap_, *transfermap_);
	// transfer vectors
	Epetra_MultiVector crosslinkerpositionstrans(*transfermap_, true);
	Epetra_MultiVector crosslinkerbondtrans(*transfermap_, 2, true);
	Epetra_Vector numbondtrans(*transfermap_, true);
	Epetra_Vector delcrosslinkerstrans(*transfermap_, true);

	//exports and reimports
	crosslinkerpositionstrans.Export(*crosslinkerpositions_, crosslinkexporter, Add);
	crosslinkerpositions_->Import(crosslinkerpositionstrans, crosslinkimporter, Insert);
	crosslinkerbondtrans.Export(*crosslinkerbond_, crosslinkexporter, Add);
	crosslinkerbond_->Import(crosslinkerbondtrans, crosslinkimporter, Insert);
	numbondtrans.Export(*numbond_, crosslinkexporter, Add);
	numbond_->Import(numbondtrans, crosslinkimporter, Insert);
	delcrosslinkerstrans.Export(delcrosslinkers, crosslinkexporter, Add);
	delcrosslinkers.Import(delcrosslinkerstrans, crosslinkimporter, Insert);

	// debug couts
	//cout<<*crosslinkerbond_<<endl;
	//cout<<*numbond_<<endl;
	//cout<<delcrosslinkers<<endl;

	// DELETION OF ELEMENTS
	int delelement = 0;
	for (int i=0; i<delcrosslinkers.MyLength(); i++)
		if (discret_.HaveGlobalElement((int)delcrosslinkers[i]))
		{
			delelement++;
			discret_.DeleteElement( (int)delcrosslinkers[i]);
		}

} //StatMechManager::SearchAndDeleteCrosslinkers()

/*----------------------------------------------------------------------*
 | delete crosslinkers listed in crosslinkerpartner_ after random check;|
 | the random check is conducted by the owner processor of the          |
 | crosslinker, respectively, and the result of that check is           |
 | communicated subsequently to all the other processors   (OLD)        |
 | (private)                                                 cyron 11/09|
 *----------------------------------------------------------------------*/
void StatMechManager::DelCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap)
{
	//get current off-rate for crosslinkers
	double koff = 0;
	if ((currentelements_ - basiselements_) < statmechparams_.get<int> ("N_crosslink", 0))
		koff = statmechparams_.get<double> ("K_OFF_start", 0.0);
	else
		koff = statmechparams_.get<double> ("K_OFF_end", 0.0);

	//probability with which a crosslink breaks up in the current time step
	double punlink = 1.0 - exp(-dt * koff);

	//create Epetra_MultiVector which saves on proc 0 for each col map node to which nodes (GID) crosslinker is to be deleted
	Epetra_MultiVector crosslinkstobedeleted(nodecolmap, statmechparams_.get<int> ("N_CROSSMAX", 0), true);
	crosslinkstobedeleted.PutScalar(-1);

	//number of crosslinks to be deleted at this node to other nodes
	Epetra_Vector numcrosslinkstobedeleted(nodecolmap, true);

	/*as already the decision of setting crosslinkers the decision of deleting them is made on processor 0 only*/
	if (discret_.Comm().MyPID() == 0)
	{
		/*go through all column map nodes (as deleting one crosslinker does not affect the decision whether another one may be deleted
		 *this has not to be done in random order (unlike when setting crosslinkers)*/
		for (int i=0; i<nodecolmap.NumMyElements(); i++)
			//go through all crosslink-partners
			for (int j=0; j<statmechparams_.get<int> ("N_CROSSMAX", 0); j++)
				//is a crosslink registered in (*crosslinkerpartner_)[j][i]? If yes: is the probability check to delete it passed?
				if ((*crosslinkerpartner_)[j][i] > -0.9 && uniformclosedgen_.random() < punlink)
				{

					//get LID of current neighbor node to which crosslink is to be deleted
					int neighbourLID = nodecolmap.LID((int) (*crosslinkerpartner_)[j][i]);

					//GID of crosslinker to be deleted
					int GID2 = min(nodecolmap.GID(i), (int) (*crosslinkerpartner_)[j][i]);
					int GID1 = max(nodecolmap.GID(i), (int) (*crosslinkerpartner_)[j][i]);
					int crosslinkerid = (GID1 + 1) * basisnodes_ + GID2;

					//decrease number of crosslinkers at both nodes where a crosslinker is deleted by one
					(*numcrosslinkerpartner_)[i]--;
					(*numcrosslinkerpartner_)[neighbourLID]--;

					//delete crosslinker from crosslinkerpartner_
					for (int k=0; k<statmechparams_.get<int> ("N_CROSSMAX", 0); k++)
					{
						if ((*crosslinkerpartner_)[k][i] == neighbourLID)
							(*crosslinkerpartner_)[k][i] = -1.0;

						if ((*crosslinkerpartner_)[k][neighbourLID] == nodecolmap.GID(i))
							(*crosslinkerpartner_)[k][neighbourLID] = -1.0;
					}

					//register at either node that a crosslinker has to be deleted
					crosslinkstobedeleted[(int) (numcrosslinkstobedeleted[i])][i] = crosslinkerid;
					numcrosslinkstobedeleted[i]++;
					crosslinkstobedeleted[(int) (numcrosslinkstobedeleted[neighbourLID])][neighbourLID] = crosslinkerid;
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
	Epetra_Export exporter(nodecolmap, noderowmap);
	Epetra_Import importer(nodecolmap, noderowmap);
	Epetra_MultiVector crosslinkstobedeletedrow(noderowmap,crosslinkstobedeleted.NumVectors(), true);
	Epetra_Vector numcrosslinkstobedeletedrow(noderowmap, true);
	Epetra_MultiVector crosslinkerpartnerrow(noderowmap,crosslinkerpartner_->NumVectors(), true);
	Epetra_Vector numcrosslinkerpartnerrow(noderowmap, true);

	crosslinkstobedeletedrow.Export(crosslinkstobedeleted, exporter, Add);
	crosslinkstobedeleted.Import(crosslinkstobedeletedrow, importer, Insert);
	numcrosslinkstobedeletedrow.Export(numcrosslinkstobedeleted, exporter, Add);
	numcrosslinkstobedeleted.Import(numcrosslinkstobedeletedrow, importer, Insert);
	numcrosslinkerpartnerrow.Export(*numcrosslinkerpartner_, exporter, Add);
	numcrosslinkerpartner_->Import(numcrosslinkerpartnerrow, importer, Insert);
	crosslinkerpartnerrow.Export(*crosslinkerpartner_, exporter, Add);
	crosslinkerpartner_->Import(crosslinkerpartnerrow, importer, Insert);

	/*at this point the information which crosslinkers are to be deleted to a certain node is present on each processor which knows a certain node at least in
	 *its column map; now each processor loops through all column map nodes and checks for the crosslinker to be deleted whether it exists on this processor (at
	 *least as a ghost element); if yes, the crosslinker element is actually deleted*/
	for (int i=0; i<crosslinkstobedeleted.MyLength(); i++)
	{
		for (int j=0; j<numcrosslinkstobedeleted[i]; j++)
		{
			//GID of crosslinker to be deleted
			int crosslinkerid = (int) crosslinkstobedeleted[j][i];

			//delete crosslinker if present on this processosr
			if (discret_.HaveGlobalElement(crosslinkerid))
				discret_.DeleteElement(crosslinkerid);

		}
	}

}//void DelCrosslinkers(const Epetra_Vector& delcrosslinkercol)

/*----------------------------------------------------------------------*
 | (public) generate gaussian randomnumbers with mean "meanvalue" and   |
 | standarddeviation "standarddeviation" for parallel use     cyron10/09|
 *----------------------------------------------------------------------*/
void StatMechManager::GenerateGaussianRandomNumbers(RCP<Epetra_MultiVector> randomnumbers, const double meanvalue, const double standarddeviation)
{
	//multivector for stochastic forces evaluated by each element based on row map
	Epetra_MultiVector randomnumbersrow(*(discret_.ElementRowMap()), randomnumbers->NumVectors());

	for (int i=0; i<randomnumbersrow.MyLength(); i++)
		for (int j=0; j<randomnumbersrow.NumVectors(); j++)
			randomnumbersrow[j][i] = standarddeviation*normalgen_.random() + meanvalue;

	//export stochastic forces from row map to column map
	Epetra_Export exporter(*discret_.ElementRowMap(), *discret_.ElementColMap());
	randomnumbers->Export(randomnumbersrow, exporter, Add);

	//now fstoch contains stochastic forces identical on all processors

	return;
} // StatMechManager::SynchronizeRandomForces()

/*----------------------------------------------------------------------*
 | (public) writing restart information for manager objects   cyron 12/08|
 *----------------------------------------------------------------------*/
void StatMechManager::StatMechWriteRestart(IO::DiscretizationWriter& output)
{
	output.WriteInt("istart", istart_);
	output.WriteDouble("starttimeoutput", starttimeoutput_);
	output.WriteDouble("endtoendref", endtoendref_);
	output.WriteInt("basisnodes", basisnodes_);
	output.WriteInt("outputfilenumber", outputfilenumber_);

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
	for (int dof = 0; dof < ndim; dof++)
		//loop through columns of "cut"
		for (int n = 0; n < cut.N(); n++)
		{
			// broken element with node_n close to "0.0"-boundary
			if (fabs(coord(dof, n + 1) - statmechparams_.get<double> ("PeriodLength", 0.0) - coord(dof, n)) < fabs(coord(dof, n + 1) - coord(dof, n)))
			{
				*broken = true;
				// set value for the spatial component in question at n-th cut
				cut(dof, n) = 1.0;
			}
			else if (fabs(coord(dof, n + 1) + statmechparams_.get<double> ("PeriodLength", 0.0) - coord(dof, n)) < fabs(coord(dof, n + 1) - coord(dof, n)))
			{
				*broken = true;
				cut(dof, n) = 2.0;
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
	for (int j=0; j<element->NumNode(); j++)
	{
		for (int k=0; k<3; k++)
		{
			// obtain k-th spatial component of the reference position of the j-th node
			double referenceposition = ((element->Nodes())[j])->X()[k];
			// get the GIDs of the node's DOFs
			vector<int> dofnode = discret_.Dof((element->Nodes())[j]);
			// store the displacement of the k-th spatial component
			double displacement = (*dis)[discret_.DofRowMap()->LID(dofnode[k])];
			// write updated components into coord (only translational
			coord(k, j) = referenceposition + displacement;
			// store current lid(s) (3 translational DOFs per node)
			if (lids != NULL)
				lids->push_back(discret_.DofRowMap()->LID(dofnode[k]));
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
	for (int i=0; i<(int)sensornodes.size(); i++)
	{
		// check if node is on proc (if not, continue)
		bool havenode = discret_.HaveGlobalNode(sensornodes.at(i));
		if (!havenode)
			continue;

		// get the node
		DRT::Node* actnode = discret_.gNode(sensornodes.at(i));
		// get the GID of the DOF of the oscillatory motion
		int gid = discret_.Dof(0, actnode)[oscdir];
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
	std::vector<int> result(N, 0);
	for (int i=0; i<(int)result.size(); i++)
		result[i] = i;

	for (int i=0; i<N; ++i)
	{
		//generate random number between 0 and i
		j = (int)floor((i + 1.0)*uniformclosedopengen_.random());

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

bool StatMechManager::CheckOrientation(const LINALG::Matrix<3, 1> direction, const Epetra_MultiVector& nodaltriadscol, const LINALG::Matrix<2, 1>& LID)
{

  //if orientation is not to be checked explicitly, this function always returns true
  if (!Teuchos::getIntegralValue<int>(statmechparams_, "CHECKORIENT"))
    return true;

  //triads on filaments at the two nodes connected by crosslinkers
  LINALG::Matrix<3, 3> T1;
  LINALG::Matrix<3, 3> T2;

  //auxiliary variable for storing a triad in quaternion form
  LINALG::Matrix<4, 1> qnode;

  //angle between filament axes at crosslinked points, respectively
  double Phi;

  //Deltaphi = Phi - Phi0, where Phi0 is the angle between crosslinked filaments with zero potential energy (i.e. the most likely one)
  double DeltaPhi;

  //triad of node on first filament which is affected by the new crosslinker
  for (int j=0; j<4; j++)
    qnode(j) = nodaltriadscol[j][(int) LID(0)];
  LARGEROTATIONS::quaterniontotriad(qnode, T1);

  //triad of node on second filament which is affected by the new crosslinker
  for (int j=0; j<4; j++)
    qnode(j) = nodaltriadscol[j][(int) LID(1)];
  LARGEROTATIONS::quaterniontotriad(qnode, T2);

  //auxiliary variable
  double scalarproduct = T1(0, 0)*T2(0,0) + T1(1,0)*T2(1,0) + T1(2,0)*T2(2,0);

  Phi = acos(scalarproduct);
  //Phi should be the acute angle between 0° and 90° between the filament axes
  if(Phi > M_PI/2.0)
    Phi = M_PI - Phi;

  DeltaPhi = Phi - statmechparams_.get<double> ("PHI0",0.0);

  //assuming bending and torsion potentials 0.5*EI*Deltaphi^2 and a Boltzmann distribution for the different states of the crosslinker we get
  double pPhi = exp(-0.5 * statmechparams_.get<double> ("CORIENT", 0.0) * DeltaPhi * DeltaPhi / statmechparams_.get<double> ("KT", 0.0));


  //pPhi = 0.0 if DeltaPhi is outside allowed range
  if(Phi < statmechparams_.get<double>("PHI0",0.0) - statmechparams_.get<double>("PHI0DEV",6.28) ||
     Phi > statmechparams_.get<double>("PHI0",0.0) + statmechparams_.get<double>("PHI0DEV",6.28))
     pPhi = 0.0;

  //crosslinker has to pass three probability checks with respect to orientation
  return(uniformclosedgen_.random() < pPhi);
} // StatMechManager::CheckOrientation

/*----------------------------------------------------------------------*
 | check the binding mode of a crosslinker        (public) mueller 07/10|
 *----------------------------------------------------------------------*/
bool StatMechManager::CheckForKinkedVisual(int eleid)
{
	bool kinked = true;
	// if element is a crosslinker
	if(eleid>basisnodes_)
	{
		DRT::Element *element = discret_.gElement(eleid);
		int lid = discret_.NodeColMap()->LID(element->NodeIds()[0]);
		for(int i=0; i<element->NumNode(); i++)
			if((*filamentnumber_)[lid]!=(*filamentnumber_)[element->NodeIds()[i]])
				kinked = false;
		return kinked;
	}
	else
		kinked = false;
	return kinked;
}//StatMechManager::CheckForKinkedVisual

/*----------------------------------------------------------------------*
 | simulation of crosslinker diffusion		        (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void StatMechManager::CrosslinkerDiffusion(const Epetra_Vector& dis, double mean, double standarddev, const double &dt)
{
	/* Here, the diffusion of crosslink molecules is handled.
	 * Depending on the number of occupied binding spots of the molecule, its movement
	 * is calculated differently.
	 */

	// export row displacement to column map format
	Epetra_Vector discol(*(discret_.DofColMap()), true);
	LINALG::Export(dis, discol);
	/*In this section, crosslinker positions are updated*/
	// diffusion processed by Proc 0
	if (discret_.Comm().MyPID()==0)
	{

		// bonding cases
		for (int i=0; i<crosslinkerpositions_->MyLength(); i++)
		{
			// number of crosslink molecule bonds to filament nodes, larger GID in case of two bonds
			bool set = false;
			int marker = 0;
			int largergid = -1;

			for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
			{
				if ((*crosslinkerbond_)[j][i]>-0.9 && !set)
				{
					marker = j;
					set = true;
				}
				if ((int)(*crosslinkerbond_)[j][i] > largergid)
					largergid = (int)(*crosslinkerbond_)[j][i];
			}

			for (int j=0; j<crosslinkerpositions_->NumVectors(); j++)
			{
				switch ((int)(*numbond_)[i])
				{
					// bonding case 1:  no bonds, diffusion
					case 0:
						(*crosslinkerpositions_)[j][i] += standarddev*normalgen_.random() + mean;
					break;
					// bonding case 2: crosslink molecule attached to one filament
					case 1:
					{
						int nodegid = (int)(*crosslinkerbond_)[marker][i];

						DRT::Node *node = discret_.lColNode(discret_.NodeColMap()->LID(nodegid));
						int dofgid = discret_.Dof(node).at(j);
						// set current crosslink position to coordinates of node which it is attached to
						(*crosslinkerpositions_)[j][i] = node->X()[j] + discol[discret_.DofColMap()->LID(dofgid)];
					}
					break;
					// bonding case 3: an actual crosslinker has been established
					case 2:
					{
						// crosslink molecule position is set to position of the node with the larger GID
						DRT::Node *node = discret_.lColNode(discret_.NodeColMap()->LID(largergid));
						int dofgid = discret_.Dof(node).at(j);
						(*crosslinkerpositions_)[j][i] = node->X()[j]	+ discol[discret_.DofColMap()->LID(dofgid)];
					}
					break;
				}
			}// end of j-loop
		}// end of i-loop
		// check for compliance with periodic boundary conditions if existent
		if (statmechparams_.get<double> ("PeriodLength", 0.0) > 0.0)
			CrosslinkerPeriodicBoundaryShift(*crosslinkerpositions_);
	} // if(discret_.Comm().MyPID()==0)
	else
		crosslinkerpositions_->PutScalar(0.0);

	// Update by Broadcast: copy this information to all processors
	Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
	Epetra_Import crosslinkimporter(*crosslinkermap_, *transfermap_);
	Epetra_MultiVector crosslinkerpositionstrans(*transfermap_, 3, true);
	crosslinkerpositionstrans.Export(*crosslinkerpositions_, crosslinkexporter, Add);
	crosslinkerpositions_->Import(crosslinkerpositionstrans, crosslinkimporter, Insert);

	return;
}// StatMechManager::CrosslinkerDiffusion

/*----------------------------------------------------------------------*
 | update crosslink molecule positions	                 								|
 |																				        (public) mueller 08/10|
 *----------------------------------------------------------------------*/
void StatMechManager::CrosslinkerIntermediateUpdate(const std::map<int,
																										LINALG::Matrix<3, 1> >& currentpositions,
																										const LINALG::SerialDenseMatrix& LID, int crosslinkernumber,
																										bool coupledmovement)
{
	// case: one-bonded crosslink molecule (i.e. two cases: +1 bond (starting from 0 bonds) or -1 bond (molecule is free to diffuse again)
	if (LID.M()==1 && LID.N()==1)
	{
		// set molecule position to node position
		map<int, LINALG::Matrix<3, 1> >::const_iterator pos0 = currentpositions.find((int)LID(0, 0));
		if (coupledmovement)
		{
			for (int i=0; i < crosslinkerpositions_->NumVectors(); i++)
				(*crosslinkerpositions_)[i][crosslinkernumber] = (pos0->second)(i);
		}
		else
		{
			// generate vector in random direction of length R_LINK to "reset" crosslink molecule position:
			// it may now reenter or leave the bonding proximity
			LINALG::Matrix<3, 1> deltapos;
			for (int i=0; i<(int)deltapos.M(); i++)
				deltapos(i) = uniformclosedgen_.random();
			deltapos.Scale(statmechparams_.get<double> ("R_LINK", 0.0) / deltapos.Norm2());
			for (int i=0; i<crosslinkerpositions_->NumVectors(); i++)
				(*crosslinkerpositions_)[i][crosslinkernumber] += deltapos(i);
		}
	}
	// case: crosslinker element
	if (LID.M()==2 && LID.N()==1)
	{
		int largerlid = max((int)LID(0,0), (int)LID(1,0));
		map<int,LINALG::Matrix<3,1> >::const_iterator updatedpos = currentpositions.find(largerlid);
		for (int i=0; i<crosslinkerpositions_->NumVectors(); i++)
			(*crosslinkerpositions_)[i][crosslinkernumber] = (updatedpos->second)(i);
	}
	return;
}// StatMechManager::CrosslinkerIntermediateUpdate

/*----------------------------------------------------------------------*
 | Initialize crosslinker positions  			        (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void StatMechManager::CrosslinkerMoleculeInit()
{
	// create crosslinker maps
	std::vector<int> gids;
	for (int i=0; i<statmechparams_.get<int> ("N_crosslink", 0); i++)
		gids.push_back(i);
	// crosslinker column and row map
	crosslinkermap_ = rcp(new Epetra_Map(-1, statmechparams_.get<int> ("N_crosslink", 0), &gids[0], 0, discret_.Comm()));
	transfermap_    = rcp(new Epetra_Map(statmechparams_.get<int> ("N_crosslink", 0), 0, discret_.Comm()));

	// create density-density-correlation-function map with
	std::vector<int> bins;
	for(int i=0; i<discret_.Comm().NumProc()*statmechparams_.get<int>("HISTOGRAMBINS", 0); i++)
		bins.push_back(i);
	ddcorrcolmap_     = rcp(new Epetra_Map(-1, discret_.Comm().NumProc()*statmechparams_.get<int>("HISTOGRAMBINS", 0), &bins[0], 0, discret_.Comm()));
	// create processor-specific density-density-correlation-function map
	ddcorrrowmap_ = rcp(new Epetra_Map(discret_.Comm().NumProc()*statmechparams_.get<int>("HISTOGRAMBINS", 1), 0, discret_.Comm()));

	double upperbound = 0.0;
	// handling both cases: with and without periodic boundary conditions
	if (statmechparams_.get<double> ("PeriodLength", 0.0) > 0.0)
		upperbound = statmechparams_.get<double> ("PeriodLength", 0.0);
	else
		upperbound = statmechparams_.get<double> ("MaxRandValue", 0.0);

	crosslinkerpositions_ = rcp(new Epetra_MultiVector(*crosslinkermap_, 3, true));

	for (int i=0; i<crosslinkerpositions_->MyLength(); i++)
		for (int j=0; j<crosslinkerpositions_->NumVectors(); j++)
			(*crosslinkerpositions_)[j][i] = upperbound * uniformclosedgen_.random();

	// initial bonding status is set (no bonds)
	crosslinkerbond_ = rcp(new Epetra_MultiVector(*crosslinkermap_, 2));
	crosslinkerbond_->PutScalar(-1.0);

	numbond_ = rcp(new Epetra_Vector(*crosslinkermap_, true));

	crosslinkonsamefilament_ = rcp(new Epetra_Vector(*crosslinkermap_, true));

	// initialize the beautiful visuals vector (aka beevee-vector)
	visualizepositions_ = rcp(new Epetra_MultiVector(*crosslinkermap_, 3, false));
	for (int i=0; i<visualizepositions_->MyLength(); i++)
		for (int j=0; j<visualizepositions_->NumVectors(); j++)
			(*visualizepositions_)[j][i] = (*crosslinkerpositions_)[j][i];

	searchforneighbours_ = rcp(new Epetra_Vector(*crosslinkermap_, false));
	searchforneighbours_->PutScalar(1.0);

	numcrossnodes_ = rcp(new Epetra_Vector(*(discret_.NodeColMap()), true));

	return;
}//StatMechManager::CrosslinkerPosInit

/*----------------------------------------------------------------------*
 | Periodic Boundary Shift for crosslinker diffusion simulation					|
 |																						  	(public) mueller 07/10|
 *----------------------------------------------------------------------*/
void StatMechManager::CrosslinkerPeriodicBoundaryShift(Epetra_MultiVector& crosslinkerpositions)
{
	for (int i=0; i<crosslinkerpositions.MyLength(); i++)
		for (int j=0; j<crosslinkerpositions.NumVectors(); j++)
		{
			if (crosslinkerpositions[j][i] > statmechparams_.get<double> ("PeriodLength", 0.0))
				crosslinkerpositions[j][i] -= statmechparams_.get<double> ("PeriodLength", 0.0);
			if (crosslinkerpositions[j][i] < 0.0)
				crosslinkerpositions[j][i] += statmechparams_.get<double> ("PeriodLength", 0.0);
		}
	return;
}// StatMechManager::CrosslinkerPeriodicBoundaryShift

#endif  // #ifdef CCADISCRET
