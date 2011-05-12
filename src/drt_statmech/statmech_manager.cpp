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
#endif
#ifdef D_BEAM3II
#include "../drt_beam3ii/beam3ii.H"
#endif
#ifdef D_TRUSS3
#include "../drt_truss3/truss3.H"
//#include "../drt_trusslm/trusslm.H"
#endif

#include "../drt_torsion3/torsion3.H"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <math.h>

//namespace with utility functions for operations with large rotations used
using namespace LARGEROTATIONS;

//MEASURETIME activates measurement of computation time for certain parts of the code
//#define MEASURETIME

/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 09/08|
 *----------------------------------------------------------------------*/
StatMechManager::StatMechManager(ParameterList& params, DRT::Discretization& discret):
statmechparams_( DRT::Problem::Instance()->StatisticalMechanicsParams() ),
unconvergedsteps_(0),
starttimeoutput_(-1.0),
endtoendref_(0.0),
istart_(0),
basisnodes_(discret.NumGlobalNodes()),
basiselements_(discret.NumGlobalElements()),
outputfilenumber_(-1),
normalgen_(0,1),
discret_(discret)
{
  //initialize random generators
  SeedRandomGenerators(0);


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

  /* and filamentnumber_ is generated based on a column map vector as each node has to
   * know about each other node its filament number in order to decide weather a crosslink may be established
   * or not; vectors is initalized with -1, which state is changed if filament numbering is used, only*/
  filamentnumber_ = rcp( new Epetra_Vector(*(discret_.NodeColMap())) );
  filamentnumber_->PutScalar(-1);


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
	CrosslinkerMoleculeInit();

	// initialize istart_ with step number of the beginning of statistical mechanics output
	istart_ = (int)(statmechparams_.get<double>("STARTTIMEOUT", 0.0) / params.get<double>("delta time", 0.01));

	return;
}// StatMechManager::StatMechManager


/*----------------------------------------------------------------------*
 | seed all random generators of this object with fixed seed if given and|
 | with system time otherwise; seedparameter is used only in the first   |
 | case to calculate the actual seed variable based on some given fixed  |
 | seed value; note that seedparameter may be any integer, but has to be |
 | been set in a deterministic way so that it for a certain call of this |
 | method at a certain point in the program always the same number       |
 | whenever the program is used                               cyron 11/10|
 *----------------------------------------------------------------------*/
void StatMechManager::SeedRandomGenerators(const int seedparameter)
{
  //integer for seeding all random generators
  int seedvariable = 0;

  /*if input flag FIXEDSEED == YES: use same random numbers in each program start;
   *to this end compute seedvariable from given parameter FIXEDSEED and some other
   *deterministic parameter seedparameter given to this method at runtime*/
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"FIXEDSEED"))
  {
    seedvariable = (statmechparams_.get<int>("INITIALSEED", 0) + seedparameter)*(discret_.Comm().MyPID() + 1);

    normalgen_.seed( (unsigned int)seedvariable );
    uniformclosedgen_.seed( (unsigned int)seedvariable );
    uniformclosedopengen_.seed( (unsigned int)seedvariable );
  }
   /*else set seed according to system time and different for each processor
   *(pseudo-random seed) if seedparameter == 0 (this allows for conveniently
   *using a random seed only at certain points in the program, e.g. only once
   *in the beginning; one has just to make sure that seedparameter == 0 does
   *not happen at any other point in the program*/
  else if(seedparameter == 0)
  {
    seedvariable = time(0)*(discret_.Comm().MyPID() + 1);

    normalgen_.seed( (unsigned int)seedvariable );
    uniformclosedgen_.seed( (unsigned int)seedvariable );
    uniformclosedopengen_.seed( (unsigned int)seedvariable );
  }


}


/*----------------------------------------------------------------------*
 | write special output for statistical mechanics (public)    cyron 09/08|
 *----------------------------------------------------------------------*/
void StatMechManager::Update(const int& istep, const double dt, Epetra_Vector& disrow,RCP<LINALG::SparseOperator>& stiff, int ndim)
{
#ifdef MEASURETIME
	const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME
	/*first we modify the displacement vector so that current nodal position at the end of current time step complies with
	 * periodic boundary conditions, i.e. no node lies outside a cube of edge length PeriodLength*/

	//if dynamic crosslinkers are used update comprises adding and deleting crosslinkers
	if (DRT::INPUT::IntegralValue<int>(statmechparams_, "DYN_CROSSLINKERS"))
	{
		// crosslink molecule diffusion
    double standarddev = sqrt(statmechparams_.get<double> ("KT", 0.0) / (2*M_PI * statmechparams_.get<double> ("ETA", 0.0) * statmechparams_.get<double> ("R_LINK", 0.0)) * dt);
    CrosslinkerDiffusion(disrow, 0.0, standarddev, dt);


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

			currentpositions[node->LID()] = currpos;
			currentrotations[node->LID()] = currrot;
		}

		// set crosslinkers, i.e. considering crosslink molecule diffusion after filaments had time to equilibrate
		if(time_>=statmechparams_.get<double>("EQUILIBTIME",0.0) || fabs(time_-statmechparams_.get<double>("EQUILIBTIME",0.0))<(dt/1e3))
		{
			SearchAndSetCrosslinkers(istep, dt, noderowmap, nodecolmap, currentpositions,currentrotations);

			SearchAndDeleteCrosslinkers(dt, noderowmap, nodecolmap, currentpositions);
		}


		// reset thermal energy to new value (simple value change for now, maybe Load Curve later on)
		if(fabs(time_-statmechparams_.get<double>("KTSWITCHTIME",0.0))<(dt/1e3))
			statmechparams_.set("KT",statmechparams_.get<double>("KTACT",statmechparams_.get<double>("KT",0.0)));

		// actions taken when reaching STARTTIMEACT
		if(fabs(time_-statmechparams_.get<double>("STARTTIMEACT",0.0))<(dt/1e3))
		{
			if(!discret_.Comm().MyPID())
			{
				cout<<"\n\n==========================================================="<<endl;
				cout<<"-- "<<crosslinkermap_->NumMyElements()<<" crosslink molecules in volume"<<endl;
				cout<<"-- removing "<<statmechparams_.get<int>("REDUCECROSSLINKSBY",0)<<" crosslinkers...\n"<<endl;
			}
			// Set a certain number of double-bonded crosslinkers free
			ReduceNumOfCrosslinkersBy(statmechparams_.get<int>("REDUCECROSSLINKSBY",0));
			if(!discret_.Comm().MyPID())
			{
				cout<<"\n-- "<<crosslinkermap_->NumMyElements()<<" crosslink molecules left in volume"<<endl;
				cout<<"===========================================================\n"<<endl;
			}
		}

		/*settling administrative stuff in order to make the discretization ready for the next time step:
		 * done in SearchAndeDeleteCrosslinkers():
		 * synchronize the Filled() state on all processors after having added or deleted elements by CheckFilledGlobally();
		 * then build new element maps and call FillComplete();
		 * done here: finally Crs matrices stiff_ has to be deleted completely and made ready
		 * for new assembly since their graph was changed*/
		stiff->Reset();

#ifdef MEASURETIME
		cout << "\n***\nadministration time: " << Teuchos::Time::wallTime() - t_admin<< " seconds\n***\n";
#endif // #ifdef MEASURETIME

	}//if(DRT::INPUT::IntegralValue<int>(statmechparams_,"DYN_CROSSLINKERS"))

#ifdef MEASURETIME
	const double Delta_t = Teuchos::Time::wallTime()-t_start;
	cout << "\n***\ntotal time: " << Delta_t<< " seconds\n***\n";
#endif // #ifdef MEASURETIME

	return;
} // StatMechManager::Update()

/*----------------------------------------------------------------------*
 | Shifts current position of nodes so that they comply with periodic   |
 | boundary conditions                                       cyron 04/10|
 *----------------------------------------------------------------------*/
void StatMechManager::PeriodicBoundaryShift(Epetra_Vector& disrow, int ndim, const double &dt)
{
  double starttime = statmechparams_.get<double>("STARTTIMEACT",0.0);
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

/*------------------------------------------------------------------------*
 | This function loops through all the elements of the discretization and |
 | tests whether truss3 are broken by periodic boundary conditions in the |
 | reference configuration; if yes initial values of jacobi determinants  |
 | are adapted in a proper way                                 cyron 03/10|
 *-----------------------------------------------------------------------*/
/*void StatMechManager::PeriodicBoundaryTrussLmInit(DRT::Element* element)
{

#ifdef D_TRUSS3

	DRT::ELEMENTS::TrussLm* truss = dynamic_cast<DRT::ELEMENTS::TrussLm*> (element);

	//3D beam elements are embeddet into R^3:
	const int ndim = 3;

	//get reference configuration of truss3 element in proper format for later call of SetUpReferenceGeometry
	vector<double> xrefe(truss->NumNode() * ndim, 0);

	for (int i=0; i<truss->NumNode(); i++)
		for (int dof = 0; dof < ndim; dof++)
			xrefe[3* i + dof] = truss->Nodes()[i]->X()[dof];

	// loop through all nodes except for the first node which remains fixed as reference node; all other nodes are
	// shifted due to periodic boundary conditions if required
	for (int i=1; i<truss->NumNode(); i++)
	{
		for (int dof = 0; dof < ndim; dof++)
		{
			//if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
			// the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
			// back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
			// is smaller than half the periodic length
			if (fabs((truss->Nodes()[i]->X()[dof]) + statmechparams_.get<double> ("PeriodLength", 0.0) - (truss->Nodes()[0]->X()[dof])) < fabs((truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof])))
				xrefe[3* i + dof] += statmechparams_.get<double> ("PeriodLength", 0.0);

			if (fabs((truss->Nodes()[i]->X()[dof]) - statmechparams_.get<double> ("PeriodLength", 0.0) - (truss->Nodes()[0]->X()[dof])) < fabs((truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof])))
				xrefe[3* i + dof] -= statmechparams_.get<double> ("PeriodLength", 0.0);
		}
	}

	//note that the third argument "true" is necessary as all truss elements have already been initialized once upon reading input file
	truss->SetUpReferenceGeometry(xrefe, true);

#endif
}*/

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
	// detection of nodes within search proximity of the crosslink molecules
	DetectNeighbourNodes(currentpositions, &nodeinpartition, numbondtrans, crosslinkerpositionstrans, crosslinkpartitiontrans, neighbourslid);
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
	Epetra_MultiVector crosslinkerbond(*transfermap_, 2, true);
	Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
	Epetra_Export crosslinkimporter(*crosslinkermap_, *transfermap_);
	if(discret_.Comm().MyPID()!=0)
		crosslinkerbond_->PutScalar(0.0);
	// forth...
	crosslinkerbond.Export(*crosslinkerbond_, crosslinkexporter, Add);
	//...and back
	crosslinkerbond_->Import(crosslinkerbond, crosslinkimporter, Insert);

	std::vector<std::vector<int> > neighbournodes(crosslinkpartitions.MyLength(), std::vector<int>());

	int maxneighbourslocal = 0;
	int maxneighboursglobal = 0;

	for(int part=0; part<crosslinkpartitions.MyLength(); part++)
	{
		// i.e. numbond==2.0
		if(crosslinkpartitions[0][part]<-0.9)
			continue;
		// determine search radius in accordance to bonding status
		// reason for this: we set the crosslinker position of molecules with numbond>0 on a node, i.e. on a binding spot of a filement/linker
		// So, the maximal radius is now R_LINK-DeltaR_LINK assuming the linker can bond in any direction.
		// When numbond==0 -> crosslinker position is thought to be in the center of a sphere with the adjusted search radii rmin and rmax.
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

	double t_search = Teuchos::Time::wallTime();
	/*preliminaries*/
	// NODAL TRIAD UPDATE
	//first get triads at all row nodes
	Epetra_MultiVector nodaltriadsrow(noderowmap, 4);
	Epetra_MultiVector nodaltriadscol(nodecolmap,4,true);
	Epetra_Import importer(nodecolmap,noderowmap);
	nodaltriadsrow.PutScalar(0);

	if (DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT"))
	{
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

      //check type of element (orientation triads are not for all elements available in the same way
      DRT::ElementType & eot = ((discret_.lRowNode(i)->Elements())[lowestidele])->ElementType();
      //if element is of type beam3ii get nodal triad
      if (eot == DRT::ELEMENTS::Beam3iiType::Instance())
      {
        DRT::ELEMENTS::Beam3ii* filele = NULL;
        filele = dynamic_cast<DRT::ELEMENTS::Beam3ii*> (discret_.lRowNode(i)->Elements()[lowestidele]);

        //check whether crosslinker is connected to first or second node of that element
        int nodenumber = 0;
        if(discret_.lRowNode(i)->Id() == ((filele->Nodes())[1])->Id() )
          nodenumber = 1;

        //save nodal triad of this node in nodaltriadrow
        for(int j=0; j<4; j++)
          nodaltriadsrow[j][i] = (filele->Qnew_[nodenumber])(j);
      }
      else if (eot == DRT::ELEMENTS::Beam3Type::Instance())
      {
        DRT::ELEMENTS::Beam3* filele = NULL;
        filele = dynamic_cast<DRT::ELEMENTS::Beam3*> (discret_.lRowNode(i)->Elements()[lowestidele]);

        //approximate nodal triad by triad at the central element Gauss point (assuming 2-noded beam elements)
        for(int j=0; j<4; j++)
          nodaltriadsrow[j][i] = (filele->Qnew_[0])(j);
      }
      else
        dserror("Filaments have to be discretized with beam3ii elements for orientation check!!!");


    }
    //export nodaltriadsrow to col map variable
    nodaltriadscol.Import(nodaltriadsrow,importer,Insert); // NODAL TRIAD UPDATE
	}

	//get current on-rate for crosslinkers
	double kon = 0;
	double starttime = statmechparams_.get<double>("STARTTIMEACT", 0.0);
	double ktswitchtime = statmechparams_.get<double>("KTSWITCHTIME", starttime);

	if(time_ <= ktswitchtime || (time_>ktswitchtime && fabs(time_-ktswitchtime) < dt/1e4))
		kon = statmechparams_.get<double>("K_ON_start",0.0);
	else
		kon = statmechparams_.get<double>("K_ON_end",0.0);

	//probability with which a crosslinker is established between crosslink molecule and neighbour node
	double plink = 1.0 - exp( -dt*kon );
	//probability with which a crosslinker blocks two binding spots on the same filament (self-binding)
	double pself = 1.0 - exp( -dt*statmechparams_.get<double>("K_ON_SELF", 0.0) );

	//Volume partitioning, assignment of nodes and molecules to partitions, search for neighbours
	// map filament (column map, i.e. entire node set on each proc) node positions to volume partitions every SEARCHINTERVAL timesteps
	RCP<Epetra_MultiVector> neighbourslid;
	if(statmechparams_.get<int>("SEARCHRES",1)>0)
		PartitioningAndSearch(currentpositions, neighbourslid);

	//cout<<"\n\nneighbourslid\n"<<*neighbourslid<<endl;
	/*the following part of the code is executed on processor 0 only and leads to the decision, at which nodes crosslinks shall be set
	 *processor 0 goes through all the crosslink molecules and checks whether a crosslink is to be set; this works precisely as follows:
	 *(1) the crosslink molecules are looped through in a random order
	 *(2) if a node has not yet reached its maximal number of crosslinks, a crosslink may be set
	 *(3) a crosslink is set if and only if the node passes a probability check
	 *note: a fully overlapping node map on processor 0 is imperative for the following part of the code to work correctly!*/

	// a vector indicating the crosslink molecule which is going to constitute a crosslinker element
	Epetra_Vector doublebond(*crosslinkermap_, true);

	int numsetelements = 0;
	if(discret_.Comm().MyPID()==0)
	{
		// obtain a random order in which the crosslinkers are addressed
		std::vector<int> order = Permutation(crosslinkermap_->NumMyElements());

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
				//continue in case of N_CROSSMAX crosslinkers at the current node (only if nodeLID>-1, i.e. no passive crosslinker)
				if(nodeLID>-1)
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

				if( uniformclosedgen_.random() < probability )
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
							// distinguish between a real nodeLID and the entry '-1', which indicates a passive crosslink molecule
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
								numsetelements++;
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
		cout << "\nsearch time: " << Teuchos::Time::wallTime() - t_search<< " seconds";
	}// if(discret_.Comm().MypPID==0)
	else
	{
		/* zerofy numcrossnodes at the beginning of each search except for Proc 0
		 * for subsequent export and reimport. This way, we guarantee redundant information
		 * on all processors.
		 * note: searchforneighbours_ and crosslinkonsamefilament_ are not being communicated
		 * to the other Procs because their information is of concern to Proc 0 only. */
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
			if(noderowmap.LID(nodeGID[0]) > -1 || noderowmap.LID(nodeGID[1]) > -1)
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
					addedelements_.push_back(newcrosslinkerGID);
					discret_.AddElement(newcrosslinker);
				}
				else
				{
					RCP<DRT::ELEMENTS::Truss3> newcrosslinker = rcp(new DRT::ELEMENTS::Truss3(newcrosslinkerGID, (discret_.gNode(GID1))->Owner() ) );
					//RCP<DRT::ELEMENTS::TrussLm> newcrosslinker = rcp(new DRT::ELEMENTS::TrussLm(newcrosslinkerGID, (discret_.gNode(GID1))->Owner() ) );

					newcrosslinker->SetNodeIds(2,globalnodeids);
					newcrosslinker->BuildNodalPointers(&nodes[0]);

					//setting up crosslinker element parameters
					newcrosslinker ->crosssec_ = statmechparams_.get<double>("ALINK",0.0);
					newcrosslinker->SetMaterial(2);

					//correct reference configuration data is computed for the new crosslinker element;
					newcrosslinker->SetUpReferenceGeometry(xrefe);

					//add new element to discretization and list this event in addedelements_
					addedelements_.push_back(newcrosslinkerGID);
					discret_.AddElement(newcrosslinker);
				}
			}
		}
	} // ADDING ELEMENTS
  // synchronization
  discret_.CheckFilledGlobally();
  discret_.FillComplete(true, false, false);
  //couts
  if(!discret_.Comm().MyPID())
  	cout<<"\n\n"<<numsetelements<<" crosslinker element(s) added!"<<endl;
#endif //#ifdef D_TRUSS3
#endif //#ifdef D_BEAM3
#endif //#ifdef D_BEAM3II
}//void StatMechManager::SearchAndSetCrosslinkers


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
	double starttime = statmechparams_.get<double>("STARTTIMEACT", 0.0);

	if (time_ <= starttime || (time_>starttime && fabs(time_-starttime)<dt/1e4))
		koff = statmechparams_.get<double> ("K_OFF_start", 0.0);
	else
		koff = statmechparams_.get<double> ("K_OFF_end", 0.0);

	//probability with which a crosslink breaks up in the current time step
	double punlink = 1.0 - exp(-dt * koff);

	// SEARCH
	// search and setup for the deletion of elements is done by Proc 0
	int numdelelements = 0;
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
								numdelelements++;
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
		// searchforneighbours_ and crosslinkonsamefilament_ are not communicated, hence, no resetting to zero here
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

        unsigned startsize = deletedelements_.size();
        std::vector<DRT::PackBuffer> vdata( startsize );

	// DELETION OF ELEMENTS
	for (int i=0; i<delcrosslinkers.MyLength(); i++)
		if (discret_.HaveGlobalElement((int)delcrosslinkers[i]))
		{
			//save the element by packing before elimination to make it restorable in case that needed
			vdata.push_back( DRT::PackBuffer() );
			discret_.gElement((int)delcrosslinkers[i])->Pack( vdata.back() );
		}
	for ( unsigned i=startsize; i<vdata.size(); ++i )
		vdata[i].StartPacking();
	for (int i=0; i<delcrosslinkers.MyLength(); i++)
		if (discret_.HaveGlobalElement((int)delcrosslinkers[i]))
		{
			//save the element by packing before elimination to make it restorable in case that needed
			deletedelements_.push_back( std::vector<char>() );
			discret_.gElement((int)delcrosslinkers[i])->Pack( vdata[deletedelements_.size()-1] );
			discret_.DeleteElement( (int)delcrosslinkers[i] );
		}
	for ( unsigned i=startsize; i<vdata.size(); ++i )
		swap( deletedelements_[i], vdata[i]() );
	/*synchronize
	 *the Filled() state on all processors after having added or deleted elements by ChekcFilledGlobally(); then build
	 *new element maps and call FillComplete();*/
	discret_.CheckFilledGlobally();
	discret_.FillComplete(true, false, false);
	if(!discret_.Comm().MyPID())
		cout<<numdelelements<<" crosslinker element(s) deleted!"<<endl;
} //StatMechManager::SearchAndDeleteCrosslinkers()

/*----------------------------------------------------------------------*
 | (private) Reduce currently existing crosslinkers by a certain        |
 | percentage.                                              mueller 1/11|
 *----------------------------------------------------------------------*/
void StatMechManager::ReduceNumOfCrosslinkersBy(const int numtoreduce)
{
	int ncrosslink = statmechparams_.get<int>("N_crosslink", 0);

	// check for the correctness of the given input value
	if(numtoreduce>ncrosslink)
		dserror("REDUCECROSSLINKSBY is greater than N_crosslink. Please check your input file!");

	/* synchronize information about adjacent elements and their respective node GIDs so that it can be used by Proc 0.
	* This is currently somehow tedious and not too elegant (node map for dealing with elements, redundant element entries etc).*/
	Epetra_Import importer(*(discret_.NodeColMap()), *(discret_.NodeRowMap()));
	// vectors holding crosslinker element GIDs
	Epetra_MultiVector crosslinkergids(*(discret_.NodeColMap()), statmechparams_.get<int>("N_CROSSMAX",0), true);
	Epetra_MultiVector crosslinkergidsrow(*(discret_.NodeRowMap()), statmechparams_.get<int>("N_CROSSMAX",0), false);
	// node IDs of crosslinker elements; NumVectors()==number of unique node GIDs per filament node with attached crosslinker elements
	Epetra_MultiVector crosslinkernodeids(*(discret_.NodeColMap()), statmechparams_.get<int>("N_CROSSMAX",0)+1, true);
	Epetra_MultiVector crosslinkernodeidsrow(*(discret_.NodeRowMap()), statmechparams_.get<int>("N_CROSSMAX",0)+1, false);

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
	Epetra_Vector deletecrosslinkerelements(*crosslinkermap_);
        Epetra_Vector deletecrosslinkermolecules(*crosslinkermap_, true);
	deletecrosslinkerelements.PutScalar(-1.0);

	// SEARCH
	// search and setup for the deletion of elements is done by Proc 0
	int numdelelements = 0;
	int numdelmolecules = 0;
	if (discret_.Comm().MyPID()==0)
	{
		//create random order in which crosslinkers are addressed
		std::vector<int> randomorder = Permutation(ncrosslink);
		// number of deleted crosslinkers
		int numdelcrosslinks = 0;
		for (int i=0; i<numbond_->MyLength(); i++)
		{
			int irandom = randomorder[i];
			if(numdelcrosslinks<numtoreduce)
			{
				// take action according to the number of bonds of a crosslink molecule
				switch ((int)(*numbond_)[irandom])
				{
					case 0:
					{
						numdelmolecules++;
						deletecrosslinkermolecules[irandom] = 1.0;
					}
					break;
					// crosslink molecule with one bond
					case 1:
					{
						numdelmolecules++;
						for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
							if ((*crosslinkerbond_)[j][irandom]>-0.9)
							{
								// obtain LID and reset crosslinkerbond_ at this position
								int nodeLID = discret_.NodeColMap()->LID((int) (*crosslinkerbond_)[j][irandom]);
								((*numcrossnodes_)[nodeLID]) -= 1.0;
								deletecrosslinkermolecules[irandom] = 1.0;
							}
					}
					break;
					// crosslinker element
					case 2:
					{
						// an actual crosslinker element exists
						if((*searchforneighbours_)[irandom]>0.9)
						{
							numdelelements++;
							// get first nodal LID ( second one not needed here)
							int nodeLID = discret_.NodeColMap()->LID((int)(*crosslinkerbond_)[0][irandom]);
							// get the crosslinker element GID and store it for later deletion:
							// k<2 maps to first crosslinker element GID, from then on k to the k-1-th crosslinker element (reason above: unique node GIDs)
							for(int k=0; k<crosslinkernodeids.NumVectors(); k++)
								if(crosslinkernodeids[k][nodeLID]==(int)(*crosslinkerbond_)[0][irandom])
								{
									if(k<2)
										deletecrosslinkerelements[irandom] = crosslinkergids[0][nodeLID];
									else
										deletecrosslinkerelements[irandom] = crosslinkergids[k-1][nodeLID];
								}

							((*numcrossnodes_)[(int)(*crosslinkerbond_)[0][irandom]]) -= 1.0;
							((*numcrossnodes_)[(int)(*crosslinkerbond_)[1][irandom]]) -= 1.0;
							deletecrosslinkermolecules[irandom] = 1.0;
						}
						else	// passive crosslink molecule
						{
							numdelmolecules++;
							for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
								if ((*crosslinkerbond_)[j][irandom]>-0.9)
								{
									// obtain LID and reset crosslinkerbond_ at this position
									int nodeLID = discret_.NodeColMap()->LID((int)(*crosslinkerbond_)[j][irandom]);
									((*numcrossnodes_)[nodeLID]) -= 1.0;
									deletecrosslinkermolecules[irandom] = 1.0;
								}
						}
					}
				}// switch ((int)(*numbond_)[i])
				numdelcrosslinks++;
			}// if(numdelcrosslinks<numtoreduce)
			else // upon reaching the given number of deleted crosslinks, exit the loop
				break;
		}// for (int i=0; i<numbond_->MyLength(); i++)
	}// if(discret_.Comm().MyPID()==0)
	else
	{
		numcrossnodes_->PutScalar(0.0);
		deletecrosslinkerelements.PutScalar(0.0);
                deletecrosslinkermolecules.PutScalar(0.0);
		// searchforneighbours_ and crosslinkonsamefilament_ are not communicated, hence, no resetting to zero here
	}

	//synchronize information about number of bonded filament nodes by exporting it to row map format and then reimporting it to column map format
	Epetra_Export exporter(*(discret_.NodeColMap()), *(discret_.NodeRowMap()));
	Epetra_Vector numcrossnodesrow(*(discret_.NodeRowMap()), true);
	numcrossnodesrow.Export(*numcrossnodes_, exporter, Add);
	numcrossnodes_->Import(numcrossnodesrow, importer, Insert);
	// exporter and importer
	Epetra_Export crosslinkexporter(*crosslinkermap_, *transfermap_);
	Epetra_Import crosslinkimporter(*crosslinkermap_, *transfermap_);
	// transfer vector
	Epetra_Vector deletecrosslinkerelementstrans(*transfermap_, true);
        Epetra_Vector deletecrosslinkermoleculestrans(*transfermap_, true);
	//export and reimport
	deletecrosslinkerelementstrans.Export(deletecrosslinkerelements, crosslinkexporter, Add);
	deletecrosslinkerelements.Import(deletecrosslinkerelementstrans, crosslinkimporter, Insert);
	deletecrosslinkermoleculestrans.Export(deletecrosslinkermolecules, crosslinkexporter, Add);
	deletecrosslinkermolecules.Import(deletecrosslinkermoleculestrans, crosslinkimporter, Insert);

	//PREPARE REBUILDING OF CROSSLINKER MAPS AND CORRESPONDING VECTORS
	// std::vectors for transfer of data to new maps and vectors
	std::vector<std::vector<double> > newcrosslinkerpositions;
	std::vector<std::vector<double> > newcrosslinkerbond;
	std::vector<std::vector<double> > newvisualizepositions;
	std::vector<int> newcrosslinkonsamefilament;
	std::vector<int> newsearchforneighbours;
	std::vector<int> newnumbond;
	// Build temporary vectors
	for(int i=0; i<deletecrosslinkermolecules.MyLength(); i++)
	{
		if(deletecrosslinkermolecules[i]<0.1)
		{
			// add the crosslinkers which are kept, to temporary vectors
			std::vector<double> crosspos;
			std::vector<double> visualpos;
			std::vector<double> crossbond;
			for(int j=0; j<crosslinkerpositions_->NumVectors(); j++)
				crosspos.push_back((*crosslinkerpositions_)[j][i]);
			for(int j=0; j<visualizepositions_->NumVectors(); j++)
				visualpos.push_back((*visualizepositions_)[j][i]);
			for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
				crossbond.push_back((*crosslinkerbond_)[j][i]);

			newcrosslinkerpositions.push_back(crosspos);
			newvisualizepositions.push_back(visualpos);
			newcrosslinkerbond.push_back(crossbond);
			newcrosslinkonsamefilament.push_back((int)(*crosslinkonsamefilament_)[i]);
			newsearchforneighbours.push_back((int)(*searchforneighbours_)[i]);
			newnumbond.push_back((int)(*numbond_)[i]);
		}
	}

	unsigned startsize = deletedelements_.size();
	std::vector<DRT::PackBuffer> vdata( startsize );

	// DELETION OF ELEMENTS
	for (int i=0; i<deletecrosslinkerelements.MyLength(); i++)
		if (discret_.HaveGlobalElement((int)deletecrosslinkerelements[i]))
		{
			//save the element by packing before elimination to make it restorable in case that needed
			//cout<<"Proc "<<discret_.Comm().MyPID()<<": deleting element
			//"<<(int)deletecrosslinkerelements[i]<<endl;
			vdata.push_back( DRT::PackBuffer() );
			discret_.gElement((int)deletecrosslinkerelements[i])->Pack( vdata.back() );
		}
	for ( unsigned i=startsize; i<vdata.size(); ++i )
		vdata[i].StartPacking();
	for (int i=0; i<deletecrosslinkerelements.MyLength(); i++)
		if (discret_.HaveGlobalElement((int)deletecrosslinkerelements[i]))
		{
			//save the element by packing before elimination to make it restorable in case that needed
			//cout<<"Proc "<<discret_.Comm().MyPID()<<": deleting element "<<(int)deletecrosslinkerelements[i]<<endl;
			deletedelements_.push_back( std::vector<char>() );
			discret_.gElement((int)deletecrosslinkerelements[i])->Pack( vdata[deletedelements_.size()-1] );
			discret_.DeleteElement( (int)deletecrosslinkerelements[i]);
		}
	for ( unsigned i=startsize; i<vdata.size(); ++i )
		swap( deletedelements_[i], vdata[i]() );
	/*synchronize
	*the Filled() state on all processors after having added or deleted elements by ChekcFilledGlobally(); then build
	*new element maps and call FillComplete();*/
	discret_.CheckFilledGlobally();
	discret_.FillComplete(true, false, false);

	// SET UP OF NEW CROSSLINKER MAPS AND VECTORS
	// create crosslinker maps
	std::vector<int> newgids;
	for (int i=0; i<(int)newcrosslinkerpositions.size(); i++)
		newgids.push_back(i);

	// crosslinker column and row maps
	crosslinkermap_ = rcp(new Epetra_Map(-1, (int)newgids.size(), &newgids[0], 0, discret_.Comm()));
	transfermap_    = rcp(new Epetra_Map((int)newgids.size(), 0, discret_.Comm()));
	//vectors
	crosslinkerpositions_ = rcp(new Epetra_MultiVector(*crosslinkermap_, 3, true));
	visualizepositions_ = rcp(new Epetra_MultiVector(*crosslinkermap_, 3, true));
	crosslinkerbond_ = rcp(new Epetra_MultiVector(*crosslinkermap_, 2));
	crosslinkonsamefilament_ = rcp(new Epetra_Vector(*crosslinkermap_, true));
	searchforneighbours_ = rcp(new Epetra_Vector(*crosslinkermap_, false));
	numbond_ = rcp(new Epetra_Vector(*crosslinkermap_, true));

	//copy information from the temporary vectors to the adjusted crosslinker vectors
	for(int i=0; i<crosslinkerpositions_->MyLength(); i++)
	{
		for(int j=0; j<crosslinkerpositions_->NumVectors(); j++)
			(*crosslinkerpositions_)[j][i] = newcrosslinkerpositions[i][j];
		for(int j=0; j<visualizepositions_->NumVectors(); j++)
			(*visualizepositions_)[j][i] = (double)newvisualizepositions[i][j];
		for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
			(*crosslinkerbond_)[j][i] = (double)newcrosslinkerbond[i][j];

		(*crosslinkonsamefilament_)[i] = (double)newcrosslinkonsamefilament[i];
		(*searchforneighbours_)[i] = (double)newsearchforneighbours[i];
		(*numbond_)[i] = (double)newnumbond[i];
	}
	if(!discret_.Comm().MyPID())
	{
		cout<<"-- "<<numdelelements<<" crosslinker elements removed"<<endl;
		cout<<"-- "<<numdelmolecules<<" free/one-bonded crosslinker molecules removed"<<endl;
		cout<<"------------------------------------------------------"<<endl;
		cout<<"-- "<<numdelelements+numdelmolecules<<" crosslinkers removed"<<endl;
	}
	return;
}//ReduceNumberOfCrosslinkersBy()


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
void StatMechManager::WriteRestart(IO::DiscretizationWriter& output)
{
	output.WriteInt("istart", istart_);
	output.WriteInt("unconvergedsteps", unconvergedsteps_);
	output.WriteDouble("starttimeoutput", starttimeoutput_);
	output.WriteDouble("endtoendref", endtoendref_);
	output.WriteInt("basisnodes", basisnodes_);
	output.WriteInt("outputfilenumber", outputfilenumber_);
  //note: beginold_,endold_,sumdispmiddle_ not considered; related methods not restartable
	output.WriteInt("basiselements", basiselements_);
	output.WriteDouble("sumsquareincpar", sumsquareincpar_);
	output.WriteDouble("sumsquareincort", sumsquareincort_);
	output.WriteDouble("sumrotmiddle", sumrotmiddle_);
	output.WriteDouble("sumsquareincmid", sumsquareincmid_);
	output.WriteDouble("sumsquareincrot", sumsquareincrot_);
  /*note: crosslinkermap_, transfermap_, ddcorrrowmap_, ddcorrcolmap_,
   * filamentnumber_, forcesensor_, octTree_ not considered, because generated in constructor*/

	//note: Using WriteVector and ReadMultiVector requires unique map of MultiVector thus export/import for restart to/from row map
	const Epetra_Map noderowmap = *discret_.NodeRowMap();
	const Epetra_Map nodecolmap = *discret_.NodeColMap();

	Epetra_Export exporter(nodecolmap,noderowmap);
	RCP<Epetra_Vector> numcrossnodesrow = rcp(new Epetra_Vector(noderowmap,true));
	numcrossnodesrow->Export(*numcrossnodes_,exporter,Insert);
	output.WriteVector("numcrossnodes",numcrossnodesrow,IO::DiscretizationWriter::nodevector);

  output.WriteRedundantDoubleVector("startindex",startindex_);

  WriteRestartRedundantMultivector(output,"crosslinkerbond",crosslinkerbond_);
  WriteRestartRedundantMultivector(output,"crosslinkerpositions",crosslinkerpositions_);
  WriteRestartRedundantMultivector(output,"numbond",numbond_);
  WriteRestartRedundantMultivector(output,"crosslinkonsamefilament",crosslinkonsamefilament_);
  WriteRestartRedundantMultivector(output,"visualizepositions",visualizepositions_);
  WriteRestartRedundantMultivector(output,"searchforneighbours",searchforneighbours_);


	return;
} // StatMechManager::WriteRestart()

/*----------------------------------------------------------------------------*
 | (public) write restart information for fully redundant   Epetra_Multivector|
 | with name "name"                                                cyron 11/10|
 *----------------------------------------------------------------------------*/
void StatMechManager::WriteRestartRedundantMultivector(IO::DiscretizationWriter& output, const string name, RCP<Epetra_MultiVector> multivector)
{
  //create stl vector to store information in multivector
  RCP<vector<double> > stlvector = rcp(new vector<double>);
  stlvector->resize(multivector->MyLength()*multivector->NumVectors());

  for (int i=0; i<multivector->MyLength(); i++)
    for (int j=0; j<multivector->NumVectors(); j++)
      (*stlvector)[i + j*multivector->MyLength()] = (*multivector)[j][i];

  //write information to output file; note that WriteRedundantDoubleVector is active on proc 0 only
  output.WriteRedundantDoubleVector(name,stlvector);

  return;
} // StatMechManager::WriteRestartRedundantMultivector()



/*----------------------------------------------------------------------*
 |read restart information for statistical mechanics (public)cyron 12/08|
 *----------------------------------------------------------------------*/
void StatMechManager::ReadRestart(IO::DiscretizationReader& reader)
{

  // read restart information for statistical mechanics
	istart_ = reader.ReadInt("istart");
	unconvergedsteps_ = reader.ReadInt("unconvergedsteps");
	starttimeoutput_ = reader.ReadDouble("starttimeoutput");
	endtoendref_ = reader.ReadDouble("endtoendref");
	basisnodes_ = reader.ReadInt("basisnodes");
	outputfilenumber_ = reader.ReadInt("outputfilenumber");
	//note: beginold_,endold_,sumdispmiddle_ not considered; related methods not restartable
	basiselements_ = reader.ReadInt("basiselements");
	sumsquareincpar_ = reader.ReadDouble("sumsquareincpar");
	sumsquareincort_ = reader.ReadDouble("sumsquareincort");
	sumrotmiddle_ = reader.ReadDouble("sumrotmiddle");
	sumsquareincmid_ = reader.ReadDouble("sumsquareincmid");
	sumsquareincrot_ = reader.ReadDouble("sumsquareincrot");
	/*note: crosslinkermap_, transfermap_, ddcorrrowmap_, ddcorrcolmap_,
	 * filamentnumber_, forcesensor_, octTree_ not considered, because generated in constructor*/

	//note: Using WriteVector and ReadMultiVector requires uniquen map of MultiVector thus export/import for restart to/from row map
  const Epetra_Map noderowmap = *(discret_.NodeRowMap());
  const Epetra_Map nodecolmap = *(discret_.NodeColMap());
  Epetra_Import importer(nodecolmap,noderowmap);
  RCP<Epetra_Vector> numcrossnodesrow = rcp(new Epetra_Vector(noderowmap),true);
  reader.ReadVector(numcrossnodesrow,"numcrossnodes");
  numcrossnodes_->Import(*numcrossnodesrow,importer,Insert);


	//Read redundant Epetra_Multivectors and STL vectors
	reader.ReadRedundantDoubleVector(startindex_,"startindex");
	ReadRestartRedundantMultivector(reader,"crosslinkerbond",crosslinkerbond_);
	ReadRestartRedundantMultivector(reader,"crosslinkerpositions",crosslinkerpositions_);
	ReadRestartRedundantMultivector(reader,"numbond",numbond_);
	ReadRestartRedundantMultivector(reader,"crosslinkonsamefilament",crosslinkonsamefilament_);
	ReadRestartRedundantMultivector(reader,"visualizepositions",visualizepositions_);
	ReadRestartRedundantMultivector(reader,"searchforneighbours",searchforneighbours_);


	return;
}// StatMechManager::ReadRestart()

/*----------------------------------------------------------------------------*
 | (public) read restart information for fully redundant Epetra_Multivector   |
 | with name "name"                                                cyron 11/10|
 *----------------------------------------------------------------------------*/
void StatMechManager::ReadRestartRedundantMultivector(IO::DiscretizationReader& reader, const string name, RCP<Epetra_MultiVector> multivector)
{
    //we assume that information was stored like for a redundant stl vector
    RCP<vector<double> > stlvector = rcp(new vector<double>);
    stlvector->resize(multivector->MyLength()*multivector->NumVectors());

    /*ReadRedundantDoubleVector reads information of stlvector on proc 0 and distributes
     *this information then to all other processors so that after the following line
     *the variable stlvector has the same information on all procs*/
    reader.ReadRedundantDoubleVector(stlvector,name);

    //transfer data from stlvector to Epetra_Multivector
    for (int i=0; i<multivector->MyLength(); i++)
      for (int j=0; j<multivector->NumVectors(); j++)
         (*multivector)[j][i] = (*stlvector)[i + j*multivector->MyLength()];

  return;
} // StatMechManager::WriteRestartRedundantMultivector()

/*-----------------------------------------------------------------------*
 | (public) saves all relevant variables *_ as *conv_ to allow  for      |
 | returning to the beginning of a time step                 cyron 11/10 |
 *-----------------------------------------------------------------------*/
void StatMechManager::WriteConv()
{
  //save relevant class variables at the very end of the time step
  crosslinkerbondconv_ = rcp(new Epetra_MultiVector(*crosslinkerbond_));
  crosslinkerpositionsconv_ = rcp(new Epetra_MultiVector(*crosslinkerpositions_));
  numcrossnodesconv_ = rcp(new Epetra_Vector(*numcrossnodes_));
  numbondconv_ = rcp(new Epetra_Vector(*numbond_));
  crosslinkonsamefilamentconv_ = rcp(new Epetra_Vector(*crosslinkonsamefilament_));
  searchforneighboursconv_ = rcp(new Epetra_Vector(*searchforneighbours_));

  //set addedelements_, deletedelements_ empty vectors
  addedelements_.resize(0);
  deletedelements_.resize(0);

  return;
} // StatMechManager::WriteConv()

/*-----------------------------------------------------------------------*
 | (public) restore state at the beginning of this time step cyron 11/10 |
 *-----------------------------------------------------------------------*/
void StatMechManager::RestoreConv(RCP<LINALG::SparseOperator>& stiff)
{
  //restore state at the beginning of time step for relevant class variables
  crosslinkerbond_ = rcp(new Epetra_MultiVector(*crosslinkerbondconv_));
  crosslinkerpositions_ = rcp(new Epetra_MultiVector(*crosslinkerpositionsconv_));
  numcrossnodes_ = rcp(new Epetra_Vector(*numcrossnodesconv_));
  numbond_ = rcp(new Epetra_Vector(*numbondconv_));
  crosslinkonsamefilament_ = rcp(new Epetra_Vector(*crosslinkonsamefilamentconv_));
  searchforneighbours_ = rcp(new Epetra_Vector(*searchforneighboursconv_));

  /*restore state of the discretization at the beginning of this time step; note that to this and
   *adding and deleting crosslinker element has to be undone exactly vice-versa compared to the way
   *it was done first in order to handle also those crosslinkers correctly added and deleted in one
   *and the same time step*/

  //loop through all elements deleted in this time step and restore them in the discretization
  for(int i=0; i<(int)deletedelements_.size(); i++)
  {
    vector<char> tmp;
    vector<char>::size_type position = 0;
    DRT::ParObject::ExtractfromPack(position,deletedelements_[i],tmp);
    DRT::ParObject* o = DRT::UTILS::Factory(tmp);
    DRT::Element* ele = dynamic_cast<DRT::Element*>(o);
    if (ele == NULL)
      dserror("Failed to build an element from the element data");
    discret_.AddElement(rcp(ele));
  }
  deletedelements_.resize(0);

  //loop through addedelements_, delete all these elements and then set addedelements_ an empty vector
  for(int i=0; i<(int)addedelements_.size(); i++)
    discret_.DeleteElement(addedelements_[i]);

  addedelements_.resize(0);

  /*settling administrative stuff in order to make the discretization ready for the next time step: synchronize
   *the Filled() state on all processors after having added or deleted elements by ChekcFilledGlobally(); then build
   *new element maps and call FillComplete(); finally Crs matrices stiff_ has to be deleted completely and made ready
   *for new assembly since their graph was changed*/
  discret_.CheckFilledGlobally();
  discret_.FillComplete(true, false, false);
  stiff->Reset();

  return;
} // StatMechManager::WriteConv()

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
	// reinitialize forcesensor_
	forcesensor_->PutScalar(-1.0);

	// loop over DOFs subjected to oscillation (by DBC)
	for (int i=0; i<(int)sensornodes.size(); i++)
	{
		// check if node is row node (this ensures processor-wise unique forcesensor_ vectors)
		if(discret_.NodeRowMap()->LID(sensornodes[i])>-1)
		{
			// get the node
			DRT::Node* actnode = discret_.gNode(sensornodes.at(i));
			// get the GID of the DOF of the oscillatory motion
			int dofgid = discret_.Dof(0, actnode)[oscdir];
			// now, get the LID
			int collid = discret_.DofColMap()->LID(dofgid);
			// activate force sensor at lid-th position
			(*forcesensor_)[collid] = 1.0;
		}
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
 | Computes current internal energy of discret_ (public)     cyron 12/10|
 *----------------------------------------------------------------------*/
void StatMechManager::ComputeInternalEnergy(const RCP<Epetra_Vector> dis, double& energy,const double& dt, const std::ostringstream& filename)
{
  ParameterList p;
  p.set("action", "calc_struct_energy");

  //add statistical vector to parameter list for statistical forces and damping matrix computation
  p.set("delta time",dt);
  p.set("ETA",statmechparams_.get<double>("ETA",0.0));
  p.set("THERMALBATH",DRT::INPUT::IntegralValue<INPAR::STATMECH::ThermalBathType>(statmechparams_,"THERMALBATH"));
  p.set<int>("FRICTION_MODEL",DRT::INPUT::IntegralValue<INPAR::STATMECH::FrictionModel>(statmechparams_,"FRICTION_MODEL"));
  p.set("SHEARAMPLITUDE",statmechparams_.get<double>("SHEARAMPLITUDE",0.0));
  p.set("CURVENUMBER",statmechparams_.get<int>("CURVENUMBER",-1));
  p.set("OSCILLDIR",statmechparams_.get<int>("OSCILLDIR",-1));
  p.set("PeriodLength",statmechparams_.get<double>("PeriodLength",0.0));


  discret_.ClearState();
  discret_.SetState("displacement", dis);
  RCP<Epetra_SerialDenseVector> energies = Teuchos::rcp(new Epetra_SerialDenseVector(1));
  energies->Scale(0.0);
  discret_.EvaluateScalars(p, energies);
  discret_.ClearState();
  energy = (*energies)(0);

  if(!discret_.Comm().MyPID())
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream writetofile;
    writetofile<<energy;
    for(int i=0; i<16; i++)
    	writetofile<<"    "<<0;
    writetofile<<endl;
		fprintf(fp, writetofile.str().c_str());
		fclose(fp);
  }

  return;

} // StatMechManager::ComputeInternalEnergy

/*----------------------------------------------------------------------*
 | checks orientation of crosslinker relative to linked filaments       |
 |                                                  (public) cyron 06/10|
 *----------------------------------------------------------------------*/

bool StatMechManager::CheckOrientation(const LINALG::Matrix<3, 1> direction, const Epetra_MultiVector& nodaltriadscol, const LINALG::Matrix<2, 1>& LID, RCP<double> phifil)
{

  //if orientation is not to be checked explicitly, this function always returns true
  if (!DRT::INPUT::IntegralValue<int>(statmechparams_, "CHECKORIENT"))
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

  if(phifil!=Teuchos::null)
  	*phifil = Phi;

  DeltaPhi = Phi - statmechparams_.get<double> ("PHI0",0.0);

  //assuming bending and torsion potentials 0.5*EI*Deltaphi^2 and a Boltzmann distribution for the different states of the crosslinker we get
  double pPhi = exp(-0.5 * statmechparams_.get<double> ("CORIENT", 0.0) * DeltaPhi * DeltaPhi / statmechparams_.get<double> ("KT", 0.0));


  //pPhi = 0.0 if DeltaPhi is outside allowed range
  if(Phi < statmechparams_.get<double>("PHI0",0.0) - 0.5*statmechparams_.get<double>("PHI0DEV",6.28) ||
     Phi > statmechparams_.get<double>("PHI0",0.0) + 0.5*statmechparams_.get<double>("PHI0DEV",6.28))
     pPhi = 0.0;

  //crosslinker has to pass three probability checks with respect to orientation
  return(uniformclosedgen_.random() < pPhi);
} // StatMechManager::CheckOrientation

/*----------------------------------------------------------------------*
 | simulation of crosslinker diffusion		        (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void StatMechManager::CrosslinkerDiffusion(const Epetra_Vector& dis, double mean, double standarddev, const double &dt)
{
	/* Here, the diffusion of crosslink molecules is handled.
	 * Depending on the number of occupied binding spots of the molecule, its motion
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
			int nodegid0 = -1;
			int nodegid1 = -1;

			// getting node gid for crosslink molecules with one bond (both active and passive, i.e. numbond={1,2})
			for (int j=0; j<crosslinkerbond_->NumVectors(); j++)
				if ((*crosslinkerbond_)[j][i]>-0.9 && !set)
				{
					nodegid0 = (int)(*crosslinkerbond_)[j][i];
					set = true;
				}
			// case: crosslinker beam element (numbond=2 AND both gid entries>-1)
			if ((*crosslinkerbond_)[0][i]>-0.9 && (*crosslinkerbond_)[1][i]>-0.9)
				nodegid1 = (int)(*crosslinkerbond_)[1][i];

			// determine crosslinker position according to bonding status (only crosslinker representations, actual crosslinker elements handled thereafter)
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
						DRT::Node *node = discret_.lColNode(discret_.NodeColMap()->LID(nodegid0));
						int dofgid = discret_.Dof(node).at(j);
						// set current crosslink position to coordinates of node which it is attached to
						(*crosslinkerpositions_)[j][i] = node->X()[j] + discol[discret_.DofColMap()->LID(dofgid)];
					}
					break;
					// bonding case 3: an actual crosslinker has been established
					case 2:
					{
						// crosslink molecule position is set to position of the node with the larger GID
						if(nodegid1 == -1)
						{
							DRT::Node *node = discret_.lColNode(discret_.NodeColMap()->LID(nodegid0));
							int dofgid = discret_.Dof(node).at(j);
							(*crosslinkerpositions_)[j][i] = node->X()[j]	+ discol[discret_.DofColMap()->LID(dofgid)];
						}
					}
					break;
				}
			}
			// crosslinker position in case of an actual crosslinker element (mid position)
			if(nodegid1 > -1)
			{
				DRT::Node *node0 = discret_.lColNode(discret_.NodeColMap()->LID(nodegid0));
				DRT::Node *node1 = discret_.lColNode(discret_.NodeColMap()->LID(nodegid1));
				for(int j=0; j<crosslinkerpositions_->NumVectors(); j++)
				{
					int dofgid0 = discret_.Dof(node0).at(j);
					int dofgid1 = discret_.Dof(node1).at(j);
					(*crosslinkerpositions_)[j][i] = (node0->X()[j]+discol[discret_.DofColMap()->LID(dofgid0)] + node1->X()[j]+discol[discret_.DofColMap()->LID(dofgid1)])/2.0;
				}
			}
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
	int ncrosslink = statmechparams_.get<int> ("N_crosslink", 0);
	int numbins = statmechparams_.get<int>("HISTOGRAMBINS", 1);
	double periodlength = statmechparams_.get<double> ("PeriodLength", 0.0);

	// create crosslinker maps
	std::vector<int> gids;
	for (int i=0; i<ncrosslink; i++)
		gids.push_back(i);
	// crosslinker column and row map
	crosslinkermap_ = rcp(new Epetra_Map(-1, ncrosslink, &gids[0], 0, discret_.Comm()));
	transfermap_    = rcp(new Epetra_Map(ncrosslink, 0, discret_.Comm()));
	startindex_ = rcp(new std::vector<double>);

	// create density-density-correlation-function map with
	if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_densitydensitycorr ||
		 DRT::INPUT::IntegralValue<int>(statmechparams_, "GMSHOUTPUT"))
	{
		std::vector<int> bins;
		for(int i=0; i<discret_.Comm().NumProc()*numbins; i++)
			bins.push_back(i);
		ddcorrcolmap_ = rcp(new Epetra_Map(-1, discret_.Comm().NumProc()*numbins, &bins[0], 0, discret_.Comm()));
		// create processor-specific density-density-correlation-function map
		ddcorrrowmap_ = rcp(new Epetra_Map(discret_.Comm().NumProc()*numbins, 0, discret_.Comm()));
		// create new trafo matrix (for later use in DDCorr Function where we evaluate in layer directions), initialize with identity matrix
		trafo_ = rcp(new LINALG::SerialDenseMatrix(3,3,true));
		for(int i=0; i<trafo_->M(); i++)
			(*trafo_)(i,i) = 1.0;
		cog_.PutScalar(periodlength/2.0);

		// start indices for parallel handling of orientation correlation etc in NumLinkerSpotsAndOrientation()
		// calculation of start indices for each processor
		// number of overall independent combinations
		int numnodes = discret_.NodeColMap()->NumMyElements();
		int numcombinations = (numnodes*numnodes-numnodes)/2;
		// combinations on each processor
		int combinationsperproc = (int)floor((double)numcombinations/(double)discret_.Comm().NumProc());
		int remainder = numcombinations%combinationsperproc;

		// get starting index tuples for later use
		startindex_->assign(2*discret_.Comm().NumProc(), 0.0);

		for(int mypid=0; mypid<discret_.Comm().NumProc()-1; mypid++)
		{
			std::vector<int> start(2,0);
			bool continueloop = false;
			bool quitloop = false;
			int counter = 0;
			int appendix = 0;
			if(mypid==discret_.Comm().NumProc()-1)
				appendix = remainder;

			// loop over crosslinker pairs
			for(int i=0; i<discret_.NodeColMap()->NumMyElements(); i++)
			{
				for(int j=0; j<discret_.NodeColMap()->NumMyElements(); j++)
				{
					if(i==(*startindex_)[2*mypid] && j==(*startindex_)[2*mypid+1])
						continueloop = true;
					if(j>i && continueloop)
					{
						if(counter<combinationsperproc+appendix)
							counter++;
						else
						{
							// new start index j
							if(j==discret_.NodeColMap()->NumMyElements()-1)
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
					(*startindex_)[2*(mypid+1)] = (double)(start[0]);
					(*startindex_)[2*(mypid+1)+1] = (double)(start[1]);
					break;
				}
			}
		}
	}
	/*cout<<"start indices: ";
	for(int i=0; i<(int)startindex_->size(); i++)
		cout<<(*startindex_)[i]<<" ";
	cout<<endl;*/

	double upperbound = 0.0;
	// handling both cases: with and without periodic boundary conditions
	if (periodlength > 0.0)
		upperbound = periodlength;
	else
		upperbound = statmechparams_.get<double> ("MaxRandValue", 0.0);

	crosslinkerpositions_ = rcp(new Epetra_MultiVector(*crosslinkermap_, 3, true));

	for (int i=0; i<crosslinkerpositions_->MyLength(); i++)
		for (int j=0; j<crosslinkerpositions_->NumVectors(); j++)
			(*crosslinkerpositions_)[j][i] = upperbound * uniformclosedgen_.random();

	// initial bonding status is set (no bonds)
	crosslinkerbond_ = rcp(new Epetra_MultiVector(*crosslinkermap_, 2));
	crosslinkerbond_->PutScalar(-1.0);
	// initial bond counter is set (no bonds)
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
}//StatMechManager::CrosslinkerMoleculeInit

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

/*----------------------------------------------------------------------*
 |  Evaluate DBCs in case of periodic BCs (public)         mueller  2/10|
 *----------------------------------------------------------------------*/
void StatMechManager::EvaluateDirichletPeriodic(ParameterList& params,
																							  RCP<Epetra_Vector> disn,
																							  RCP<Epetra_Vector> dirichtoggle,
																							  RCP<Epetra_Vector> invtoggle)
/*The idea behind EvaluateDirichletPeriodic() is simple:
 * Give Dirichlet values to nodes of an element that is broken in
 * z-direction due to the application of periodic boundary conditions.
 * The motion of the node close to z=0.0 in the cubic volume of edge
 * length l (==PeriodLength in this case) is inhibited in direction of the
 * oscillatory motion. The oscillation is imposed on the node close to z=l.
 * This method is triggered in case of PeriodLength>0.0 (i.e. Periodic BCs
 * exist). Since the DBC setup happens dynamically by checking element
 * positions with each new time step, the static definition of DBCs in the input
 * file is only used to get the direction of the oscillatory motion as well as
 * the time curve. Therefore, only one DBC needs to be specified.
 *
 * How it works:
 * Each time this method is called, the system vector and the toggle vector are
 * modified to fit the current geometric situation.
 * DOFs holdings Dirichlet values are marked by setting the corresponding toggle
 * vector component to 1.0. In case of an element which was broken the step before
 * and is now whole again,just the toggle vector components in question are
 * reset to 0.0.
 * A position vector deltadbc is needed in order to calculate the correct Dirichlet
 * values to be imposed on nodes of an element which has drifted over the boundaries
 * and thus has been broken.
 * These positions are used to calculate the zero position of the oscillation which then can
 * be added to the time curve value in DoDirichletConditionPeriodic().
 */
{
#ifdef MEASURETIME
  const double t_start = Teuchos::Time::wallTime();
#endif // #ifdef MEASURETIME

	if (!(discret_.Filled())) dserror("FillComplete() was not called");
	if (!(discret_.HaveDofs())) dserror("AssignDegreesOfFreedom() was not called");

//----------------------------------- some variables
	// indicates broken element
	bool broken;
	// time trigger
	bool usetime = true;
	// to avoid redundant or wrong actions when filling vectors or deleting the last element of the free nodes vector
	bool alreadydone = false;
	// store node Id of previously handled node
  int tmpid = -1;
  // store LIDs of element nodes
  vector<int>	lids;
  // vectors to manipulate DBC properties
  vector<int> oscillnodes;
  vector<int> fixednodes;
  vector<int> freenodes;

	// get the apmlitude of the oscillation
	double amp = statmechparams_.get<double>("SHEARAMPLITUDE",0.0);
	// retrieve direction of oscillatory motion
	int oscdir = statmechparams_.get<int>("OSCILLDIR",-1);
	// retrieve number of time curve that is to be applied
	int curvenumber = statmechparams_.get<int>("CURVENUMBER",0)-1;

	if(oscdir!=0 && oscdir!=1 && oscdir!=2)
		dserror("Please define the StatMech Parameter OSCILLDIR correctly");

  // get the current time
	const double time = params.get<double>("total time", 0.0);
	double dt = params.get<double>("delta time", 0.01);
	double starttime = statmechparams_.get<double>("STARTTIMEACT",-1.0);
	// check if time has superceeded start. If not, do nothing (i.e. no application of Dirichlet values) and just return!
	if(time <= starttime || (time > starttime && fabs(time-starttime)<dt/1e4) || curvenumber == -1)
		return;
	if (time<0.0)
		usetime = false;

//---------------------------------------------------------- loop over elements
	// increment vector for dbc values
	Epetra_Vector deltadbc(*(discret_.DofRowMap()), true);

  // loop over row elements
	for(int lid=0; lid<discret_.NumMyRowElements(); lid++)
	{
		// An element used to browse through local Row Elements
	  DRT::Element* element = discret_.lRowElement(lid);

	  // skip element if it is a crosslinker element or in addition, in case of the Bead Spring model, Torsion3 elements
	  if(element->Id() > basisnodes_ || element->Id() >= statmechparams_.get<int>("NUM_EVAL_ELEMENTS", basiselements_))
	  	continue;

	  // number of translational DOFs (not elegant but...ah well...!)
	  int numdof = 3;
	  // positions of nodes of an element with n nodes
	  LINALG::SerialDenseMatrix coord(3,(int)discret_.lRowElement(lid)->NumNode(), true);
	  // indicates location, direction and component of a broken element with n nodes->n-1 possible cuts
	  LINALG::SerialDenseMatrix cut(3,(int)discret_.lRowElement(lid)->NumNode()-1,true);
//-------------------------------- obtain nodal coordinates of the current element
	  // get nodal coordinates and LIDs of the nodal DOFs
	  GetElementNodeCoords(element, disn, coord, &lids);
//-----------------------detect broken/fixed/free elements and fill position vector
	  // determine existence and location of broken element
	  CheckForBrokenElement(coord, cut, &broken);

	  // loop over number of cuts (columns)
	  for(int n=0; n<cut.N(); n++)
	  {
	  	// case 1: broken element (in z-dir); node_n+1 oscillates, node_n is fixed in dir. of oscillation
			if(broken && cut(2,n)==1.0)
			{
				// indicates beginning of a new filament (in the very special case that this is needed)
				bool newfilament = false;
				// check for case: last element of filament I as well as first element of filament I+1 broken
				if(tmpid!=element->Nodes()[n]->Id() && alreadydone)
				{
					// in this case, reset alreadydone...
					alreadydone = false;
					// ...and set newfilament to true. Otherwise the last free nodes vector element will be deleted
					newfilament = true;
				}
				// add GID of fixed node to fixed-nodes-vector (to be added to condition later)
				if(!alreadydone)
					fixednodes.push_back(element->Nodes()[n]->Id());
				// add GID of oscillating node to osc.-nodes-vector
				oscillnodes.push_back(element->Nodes()[n+1]->Id());
				/* When an element is cut, there are always two nodes involved: one that is subjected to a fixed
				 * displacement in one particular direction (oscdir), another which oscillates in the same direction.
				 * The following code section calculates increments for both node types and stores this increment in
				 * a vector deltadbc. This increment is later added to the nodes' displacement of the preceding time step.
				 */
				// the new displacement increments
				// incremental displacement for a fixed node...(all DOFs = 0.0)
				for(int i=0; i<numdof; i++)
					if(i==oscdir)
						deltadbc[lids.at(numdof*n+i)] = 0.0;
				// incremental Dirichlet displacement for an oscillating node (all DOFs except oscdir = 0.0)
				// time step size
				double dt = params.get<double>("delta time" ,-1.0);
				// time curve increment
				double tcincrement = 0.0;
				if(curvenumber>-1)
					tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
												DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
				deltadbc[lids.at(numdof*(n+1)+oscdir)] = amp*tcincrement;
				// delete last Id of freenodes if it was previously and falsely added
				if(element->Nodes()[n]->Id()==tmpid && !alreadydone && !newfilament)
					freenodes.pop_back();
				// store gid of the "n+1"-node to avoid overwriting during the following iteration,
				// e.g. oscillating node becomes free if the following CheckForBrokenElement() call yields "!broken".
				tmpid = element->Nodes()[n+1]->Id();
				// Set to true to initiate certain actions if the following element is also broken.
				// If the following element isn't broken, alreadydone will be reset to false (see case: !broken)
				alreadydone=true;
			}// end of case 1

			// case 2: broken element (in z-dir); node_n oscillates, node_n+1 is fixed in dir. of oscillation
			if(broken && cut(2,n)==2.0)
			{
				bool newfilament = false;
				if(tmpid!=element->Nodes()[n]->Id() && alreadydone)
				{
					alreadydone = false;
					newfilament = true;
				}
				if(!alreadydone)
					oscillnodes.push_back(element->Nodes()[n]->Id());
				fixednodes.push_back(element->Nodes()[n+1]->Id());
				// oscillating node
				double dt = params.get<double>("delta time" ,-1.0);
				double tcincrement = 0.0;
				if(curvenumber>-1)
					tcincrement = DRT::Problem::Instance()->Curve(curvenumber).f(time) -
												DRT::Problem::Instance()->Curve(curvenumber).f(time-dt);
				deltadbc[lids.at(numdof*n+oscdir)] = amp*tcincrement;
				// fixed node
				for(int i=0; i<numdof; i++)
					if(i==oscdir)
						deltadbc[lids.at(numdof*(n+1)+i)] = 0.0;
				if(element->Nodes()[n]->Id()==tmpid && !alreadydone && !newfilament)
					freenodes.pop_back();
				tmpid = element->Nodes()[n+1]->Id();
				alreadydone = true;
			} // end of case 2

			// case 3: unbroken element or broken in another than z-direction
			if(cut(2,n)==0.0)
			{
				if(element->Nodes()[n]->Id()!=tmpid)
				{
					freenodes.push_back(element->Nodes()[n]->Id());
					freenodes.push_back(element->Nodes()[n+1]->Id());
				}
				else
					freenodes.push_back(element->Nodes()[n+1]->Id());
				tmpid=element->Nodes()[n+1]->Id();
				// set to false to handle annoying special cases
				alreadydone = false;
			}
	  }
	}

//---------check/set force sensors anew for each time step
	// add DOF LID where a force sensor is to be set
  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->StatisticalMechanicsParams(),"DYN_CROSSLINKERS"))
  	UpdateForceSensors(oscillnodes, oscdir);
  //if(!discret_.Comm().MyPID())
  //	cout<<"\n=========================================="<<endl;
  //for(int pid=0; pid<discret_.Comm().NumProc();pid++)
  //{
	//	if(pid==discret_.Comm().MyPID())
  //		cout<<"UpdateForceSensors_"<<pid<<": "<<oscillnodes.size()<< " nodes @ t="<<time<<endl;
  //	discret_.Comm().Barrier();
  //}
  //cout<<"==========================================\n"<<endl;

//------------------------------------set Dirichlet values
	// preliminary
  DRT::Node* node = discret_.gNode(discret_.NodeRowMap()->GID(0));
	int numdof = (int)discret_.Dof(node).size();
	vector<int>  	 addonoff(numdof, 0);

  // set condition for oscillating nodes
	// inhibit all DOFs (for now, testing)
	for(int i=0; i<numdof; i++)
		if(i==oscdir)
			addonoff.at(i) = 1;

  // do not do anything if vector is empty
  if(!oscillnodes.empty())
  	DoDirichletConditionPeriodic(&oscillnodes, &addonoff, disn, dirichtoggle, invtoggle, deltadbc);

  // set condition for fixed nodes
	if(!fixednodes.empty())
  	DoDirichletConditionPeriodic(&fixednodes, &addonoff, disn, dirichtoggle, invtoggle, deltadbc);

  // set condition for free or recently set free nodes
  for(int i=0; i<numdof; i++)
  	if(i==oscdir)
  		addonoff.at(i) = 0;

	if(!freenodes.empty())
  	DoDirichletConditionPeriodic(&freenodes, &addonoff, disn, dirichtoggle, invtoggle, deltadbc);

#ifdef MEASURETIME
  const double t_end = Teuchos::Time::wallTime();
  if(!discret_.Comm().MyPID())
  	cout<<"DBC Evaluation time: "<<t_end-t_start<<endl;
#endif // #ifdef MEASURETIME

	return;
}

/*----------------------------------------------------------------------*
 |  fill system vector and toggle vector (public)          mueller  3/10|
 *----------------------------------------------------------------------*/
void StatMechManager::DoDirichletConditionPeriodic(vector<int>* nodeids,
																									 vector<int>* onoff,
																									 RCP<Epetra_Vector> disn,
																									 RCP<Epetra_Vector> dirichtoggle,
																									 RCP<Epetra_Vector> invtoggle,
																									 Epetra_Vector& deltadbc)
/*
 * This basically does the same thing as DoDirichletCondition() (to be found in drt_discret_evaluate.cpp),
 * but with the slight difference of taking current displacements into account.
 * Time curve values aren't added to the reference position(s) of the discretization as usual,
 * but to the latest known 0-position(s). These positions are calculated using the deltadbc
 * vector holding the latest incremental Dirichlet displacement. It is added to the displacement
 * at the end of the preceding time step.
 */
{
	/*/ "condition output"
	cout<<"Node Ids: ";
	for(int i=0; i<(int)nodeids->size(); i++)
		cout<<nodeids->at(i)<<" ";
	cout<<"onoff: ";
	for(int i=0; i<(int)discret_.Dof(0,discret_.gNode(nodeids->at(0))).size(); i++)
		cout<<onoff->at(i)<<" ";
	cout<<endl;*/

	// some checks for errors
	if (!nodeids) dserror("No Node IDs were handed over!");

	// get the condition properties
	const int nnode = nodeids->size();

	// loop over all nodes in condition
	for (int i=0; i<nnode; ++i)
	{
		// do only nodes in my row map
		if (!discret_.NodeRowMap()->MyGID(nodeids->at(i))) continue;
		DRT::Node* actnode = discret_.gNode(nodeids->at(i));
		if (!actnode) dserror("Cannot find global node %d",nodeids->at(i));
		// call explicitly the main dofset, i.e. the first column
		vector<int> dofs = discret_.Dof(0,actnode);
		const unsigned numdf = dofs.size();

		// loop over DOFs
		for (unsigned j=0; j<numdf; ++j)
		{
			// get the LID of the currently handled DOF
			const int lid = (*disn).Map().LID(dofs[j]);

			// if DOF in question is not subject to DBCs (anymore)
      if (onoff->at(j)==0)
      {
        if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
        if (dirichtoggle!=Teuchos::null)
        {
        	// turn off application of Dirichlet value
          (*dirichtoggle)[lid] = 0.0;
          // in addition, modify the inverse vector (needed for manipulation of the residual vector)
          (*invtoggle)[lid] = 1.0;
        }
        continue;
      }
//---------------------------------------------Dirichlet Value Assignment
	    // assign value
	    if (lid<0) dserror("Global id %d not on this proc in system vector", dofs[j]);
	    if (disn != Teuchos::null)
	    	(*disn)[lid] += deltadbc[lid];
	    // set toggle vector and the inverse vector
	    if (dirichtoggle != Teuchos::null)
	    {
	      (*dirichtoggle)[lid] = 1.0;
	      (*invtoggle)[lid] = 0.0;
	    }
	  }  // loop over nodal DOFs
	}  // loop over nodes
	return;
}
#endif  // #ifdef CCADISCRET
