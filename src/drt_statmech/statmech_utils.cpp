/*!----------------------------------------------------------------------
\file statmech_utils.cpp
\brief management and auxiliary functions for statistical mechanics

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/

#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_io/io.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3r/beam3r.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_beam3cl/beam3cl.H"
#include "../drt_truss3/truss3.H"
#include "../drt_truss3cl/truss3cl.H"
#include "../drt_torsion3/torsion3.H"

// defines flags for debugging and optimization purposes
//#define MEASURETIME
//#define DEBUGCOUT

/*----------------------------------------------------------------------*
 | seed all random generators of this object with fixed seed if given and|
 | with system time otherwise; seedparameter is used only in the first   |
 | case to calculate the actual seed variable based on some given fixed  |
 | seed value; note that seedparameter may be any integer, but has to be |
 | been set in a deterministic way so that it for a certain call of this |
 | method at a certain point in the program always the same number       |
 | whenever the program is used                               cyron 11/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::SeedRandomGenerators(const int seedparameter, const int seedparameter2)
{
  //integer for seeding all random generators
  int seedvariable = 0;

  Teuchos::ParameterList statmechparams = GetStatMechParams();
  double randnumtimeinc = statmechparams.get<double>("RANDNUMTIMEINT",-1.0);

  /*if input flag FIXEDSEED == YES: use same random numbers in each program start;
   *to this end compute seedvariable from given parameter FIXEDSEED and some other
   *deterministic parameter seedparameter given to this method at runtime*/
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"FIXEDSEED"))
  {
    //Decide if random numbers should change in every time step...
    if(randnumtimeinc==-1.0)
      seedvariable = (statmechparams_.get<int>("INITIALSEED", 0) + seedparameter)*(discret_->Comm().MyPID() + 1);
    //...or not before a prescribed interval RANDNUMTIMEINT
    else
      seedvariable = (statmechparams_.get<int>("INITIALSEED", 0) + seedparameter2)*(discret_->Comm().MyPID() + 1);

    randomnumbergen_.seed((unsigned int)seedvariable);
  }
   /*else set seed according to system time and different for each processor
   *(pseudo-random seed) if seedparameter == 0 (this allows for conveniently
   *using a random seed only at certain points in the program, e.g. only once
   *in the beginning; one has just to make sure that seedparameter == 0 does
   *not happen at any other point in the program*/
  else if(seedparameter == 0)
  {
    seedvariable = time(0)*(discret_->Comm().MyPID() + 1);

    randomnumbergen_.seed((unsigned int)seedvariable);
  }

#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
  uniformgen_ = Teuchos::rcp(new boost::uniform_01<randnumgen&>(randomnumbergen_));
#else
  boost::uniform_01<>           uniformdist;
  uniformgen_ = Teuchos::rcp(new boost::variate_generator<randnumgen&,boost::uniform_01<> >(randomnumbergen_,uniformdist));
#endif
  boost::normal_distribution<>  normaldist(0.0,1.0);
  normalgen_ = Teuchos::rcp(new boost::variate_generator<randnumgen&,boost::normal_distribution<> >(randomnumbergen_,normaldist));

  return;
}

/*----------------------------------------------------------------------*
 | Create fully overlapping node map             (private) mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CreateFullyOverlappingNodeMap()
{
  const Epetra_Map noderowmap = *(discret_->NodeRowMap());

  // fill my own row node ids into vector sdata
  std::vector<int> sdata(noderowmap.NumMyElements());
  for (int i=0; i<noderowmap.NumMyElements(); ++i)
    sdata[i] = noderowmap.GID(i);

  /*
   * If current processor has elements it writes the processor ID into stproc.
   * In this case, elements doesn't mean finite elements, but rather components
   * of a particular array.
   */
  std::vector<int> stproc(0);


  if (noderowmap.NumMyElements())
    stproc.push_back(discret_->Comm().MyPID());


  /*
   * Information how many processors work at all.
   * Size of variable allproc= number of processors
   */
  std::vector<int> allproc(discret_->Comm().NumProc());

  //in case of n processors allproc becomes a std::vector with entries (0,1,...,n-1)
  for (int i=0; i<discret_->Comm().NumProc(); ++i) allproc[i] = i;

  //declaring new variable into which the information of stproc on all processors is gathered
  std::vector<int> rtproc(0);

  /*
   * Gathers information of stproc and writes it into rtproc. In the end rtproc is a vector which
   * contains the numbers of all processors which have elements. Lookup method "LINALG::Gather" to
   * understand more about the variables.
   */
  LINALG::Gather<int>(stproc,rtproc,discret_->Comm().NumProc(),&allproc[0],discret_->Comm());

  /*in analogy to stproc and rtproc the variable rdata gathers all the element numbers which are
   * stored on different processors in their own variables sdata. Therefore, each processor gets
   * the information about all the nodes numbers existing in this problem*/
  std::vector<int> rdata;

  // gather all gids of nodes redundantly from sdata into rdata
  LINALG::Gather<int>(sdata,rdata,(int)rtproc.size(),&rtproc[0],discret_->Comm());

  // build completely overlapping map (on participating processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)rdata.size(),&rdata[0],0,discret_->Comm()));
  sdata.clear();
  stproc.clear();
  rdata.clear();
  allproc.clear();
  // rtproc still in use

  //pass new fully overlapping column map to discretization
  discret_->ExportColumnNodes(*newnodecolmap);
  /*rebuild discretization based on the new column map so that each processor creates new ghost elements
   * if necessary; after the following line we have a discretization, where each processor has a fully
   * overlapping column map regardless of how overlapping was managed when starting BACI; having ensured
   * this allows convenient and correct (albeit not necessarily efficient) use of search algorithms and
   * crosslinkers in parallel computing*/
  discret_->FillComplete(true,false,false);
  return;
}

/*----------------------------------------------------------------------*
 | Checks for row node for a given nodal vector          mueller (01/14)|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::CheckIfRowNodes(const Epetra_Map& noderowmap, Teuchos::RCP<std::vector<int> > nodeids)
{
  bool hasrownode = false;
  for(int k=0; k<(int)nodeids->size(); k++)
    if(noderowmap.LID(nodeids->at(k)) > -1)
    {
      hasrownode = true;
      break;
    }
  return hasrownode;
}

/*----------------------------------------------------------------------*
 | Checks for row node for a given node                mukherjee (03/15)|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::CheckRowNode(const Epetra_Map& noderowmap, int&  nodeid)
{
  bool hasrownode = false;
  if(noderowmap.LID(nodeid) > -1)
    hasrownode = true;

  return hasrownode;
}

/*------------------------------------------------------------------------*
 | This function loops through all the elements of the discretization and |
 | tests whether beam3 are broken by periodic boundary conditions in the  |
 | reference configuration; if yes initial values of curvature and jacobi |
 | determinants are adapted in a proper way                    cyron 02/10|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryBeam3Init(DRT::Element* element)
{
  DRT::ELEMENTS::Beam3* beam = dynamic_cast<DRT::ELEMENTS::Beam3*> (element);

  //3D beam elements are embedded into R^3:
  const int ndim = 3;

  /*get reference configuration of beam3 element in proper format for later call of SetUpReferenceGeometry;
   * note that rotrefe for beam3 elements is related to the entry in the global total Lagrange displacement
   * vector related to a certain rotational degree of freedom; as the displacement is initially zero also
   * rotrefe is set to zero here*/
  std::vector<double> xrefe(beam->NumNode() * ndim, 0);
  std::vector<double> rotrefe(beam->NumNode() * ndim, 0);

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
      if (fabs((beam->Nodes()[i]->X()[dof]) + periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof])) < fabs((beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof])))
        xrefe[3* i + dof] += periodlength_->at(dof);

      if (fabs((beam->Nodes()[i]->X()[dof]) - periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof])) < fabs((beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof])))
        xrefe[3* i + dof] -= periodlength_->at(dof);
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
    break;
  }
  return;
}

/*------------------------------------------------------------------------*
 | Beam3r initialization when periodic BCs are applied        cyron 02/10|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryBeam3rInit(DRT::Element* element)
{
  // note: in analogy to PeriodicBoundaryBeam3Init()

  DRT::ELEMENTS::Beam3r* beam = dynamic_cast<DRT::ELEMENTS::Beam3r*>(element);

  const int ndim = 3;

  std::vector<double> xrefe(beam->NumNode()*ndim,0);
  std::vector<double> rotrefe(beam->NumNode()*ndim,0);

  for(int i=0;i<beam->NumNode();i++)
    for(int dof=0; dof<ndim; dof++)
    {
      xrefe[3*i+dof] = beam->Nodes()[i]->X()[dof];
      rotrefe[3*i+dof] = 0.0;
    }

  for(int i=1;i<beam->NumNode();i++)
  {
    for(int dof=0; dof<ndim; dof++)
    {
      if( fabs( (beam->Nodes()[i]->X()[dof]) + periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] += periodlength_->at(dof);
      if( fabs( (beam->Nodes()[i]->X()[dof]) - periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] -= periodlength_->at(dof);
    }
  }

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
    break;
  }
  return;
}

/*------------------------------------------------------------------------*
 | BeamCL initialization when periodic BCs are applied       mueller 10/12|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryBeamCLInit(DRT::Element* element, bool setuprefgeo)
{
  // note: in analogy to PeriodicBoundaryBeam3Init()

  DRT::ELEMENTS::BeamCL* beam = dynamic_cast<DRT::ELEMENTS::BeamCL*>(element);
  if(beam->NumNode()!= 4)
    dserror("PeriodicBoundaryBeamCLInit() only implemented for 4-noded BeamCL");
  const int ndim = 3;

  std::vector<double> xrefe(6,0);
  std::vector<double> rotrefe(6,0);

  xrefe=beam->XRef();
  for(int i=0;i<6;i++)
   rotrefe[i]=0.0;
  for(int dof=0; dof<ndim; dof++)
  {
    if( fabs( beam->XRef()[3+dof] + periodlength_->at(dof) - beam->XRef()[dof] ) < fabs( beam->XRef()[3+dof] - beam->XRef()[dof] ) )
      xrefe[3+dof] += periodlength_->at(dof);
    if( fabs( beam->XRef()[3+dof] - periodlength_->at(dof) - beam->XRef()[dof] ) < fabs( beam->XRef()[3+dof] - beam->XRef()[dof] ) )
      xrefe[3+dof] -= periodlength_->at(dof);
  }

  beam->SetUpReferenceGeometry<2>(xrefe,rotrefe,true);
}

/*------------------------------------------------------------------------*
 | Beam3eb initialization when periodic BCs are applied      mueller 03/13|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryBeam3ebInit(DRT::Element* element)
{
  // note: in analogy to PeriodicBoundaryBeam3Init()

  DRT::ELEMENTS::Beam3eb* beam = dynamic_cast<DRT::ELEMENTS::Beam3eb*>(element);
  const int ndim = 3;
  std::vector<double> xrefe(beam->NumNode()*ndim,0);

  for(int i=0;i<beam->NumNode();i++)
    for(int dof=0; dof<ndim; dof++)
      xrefe[3*i+dof] = beam->Nodes()[i]->X()[dof];

  for(int i=1;i<beam->NumNode();i++)
  {
    for(int dof=0; dof<ndim; dof++)
    {
      if( fabs( (beam->Nodes()[i]->X()[dof]) + periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] += periodlength_->at(dof);
      if( fabs( (beam->Nodes()[i]->X()[dof]) - periodlength_->at(dof) - (beam->Nodes()[0]->X()[dof]) ) < fabs( (beam->Nodes()[i]->X()[dof]) - (beam->Nodes()[0]->X()[dof]) ) )
        xrefe[3*i+dof] -= periodlength_->at(dof);
    }
  }
  beam->SetUpReferenceGeometry(xrefe,true);
}


/*------------------------------------------------------------------------*
 | This function loops through all the elements of the discretization and |
 | tests whether truss3 are broken by periodic boundary conditions in the |
 | reference configuration; if yes initial values of jacobi determinants  |
 | are adapted in a proper way                                 cyron 03/10|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryTruss3Init(DRT::Element* element)
{
  DRT::ELEMENTS::Truss3* truss = dynamic_cast<DRT::ELEMENTS::Truss3*> (element);

  //3D truss elements are embeddet into R^3:
  const int ndim = 3;

  /*get reference configuration of truss3 element in proper format for later call of SetUpReferenceGeometry*/
  std::vector<double> xrefe(truss->NumNode() * ndim, 0);
  std::vector<double> rotrefe(truss->NumNode() * ndim, 0);

  for (int i=0; i<truss->NumNode(); i++)
    for (int dof = 0; dof < ndim; dof++)
      xrefe[3* i + dof] = truss->Nodes()[i]->X()[dof];

  for(int i=0;i<6;i++)
   rotrefe[i]=0.0;

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
      if (fabs((truss->Nodes()[i]->X()[dof]) + periodlength_->at(dof) - (truss->Nodes()[0]->X()[dof])) < fabs((truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof])))
        xrefe[3* i + dof] += periodlength_->at(dof);

      if (fabs((truss->Nodes()[i]->X()[dof]) - periodlength_->at(dof) - (truss->Nodes()[0]->X()[dof])) < fabs((truss->Nodes()[i]->X()[dof]) - (truss->Nodes()[0]->X()[dof])))
        xrefe[3* i + dof] -= periodlength_->at(dof);
    }
  }
    /*note that the third argument "true" is necessary as all truss elements have already been initialized once upon reading input file*/
  truss->SetUpReferenceGeometry(xrefe,rotrefe,true);
}

/*------------------------------------------------------------------------*
 | Truss3cl initialization when periodic BCs are applied   mukherjee 01/14|
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryTruss3CLInit(DRT::Element* element)
{
  // note: in analogy to PeriodicBoundaryBeam3clInit()

  DRT::ELEMENTS::Truss3CL* truss3cl = dynamic_cast<DRT::ELEMENTS::Truss3CL*>(element);


  if(truss3cl->NumNode()!= 4)
    dserror("PeriodicBoundarytTrussCLInit() only implemented for 4-noded Truss eleem");
  const int ndim = 3;

  /*get reference configuration of truss3 element in proper format for later call of SetUpReferenceGeometry*/
  std::vector<double> xrefe(truss3cl->NumNode() * ndim, 0);

  xrefe=truss3cl->XRef();
  for (int i=0; i<truss3cl->NumNode(); i++)
    for (int dof = 0; dof < ndim; dof++)
      xrefe[3* i + dof] = truss3cl->Nodes()[i]->X()[dof];

  for(int dof=0; dof<ndim; dof++)
  {
    if( fabs( truss3cl->XRef()[3+dof] + periodlength_->at(dof) - truss3cl->XRef()[dof] ) < fabs( truss3cl->XRef()[3+dof] - truss3cl->XRef()[dof] ) )
      xrefe[3+dof] += periodlength_->at(dof);
    if( fabs( truss3cl->XRef()[3+dof] - periodlength_->at(dof) - truss3cl->XRef()[dof] ) < fabs( truss3cl->XRef()[3+dof] - truss3cl->XRef()[dof] ) )
      xrefe[3+dof] -= periodlength_->at(dof);
  }

  truss3cl->SetUpReferenceGeometry(xrefe,true);
}

/*-----------------------------------------------------------------------*
 | (public) saves all relevant variables *_ as *conv_ to allow  for      |
 | returning to the beginning of a time step                 cyron 11/10 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::WriteConv(Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  //save relevant class variables at the very end of the time step
  crosslinkerbondconv_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerbond_));
  crosslinkerpositionsconv_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerpositions_));
  bspotstatusconv_ = Teuchos::rcp(new Epetra_Vector(*bspotstatus_));
  numbondconv_ = Teuchos::rcp(new Epetra_Vector(*numbond_));
  crosslink2elementconv_ = Teuchos::rcp(new Epetra_Vector(*crosslink2element_));
  if(linkermodel_==statmech_linker_myosinthick)
    additionalcross2eleconv_ = Teuchos::rcp(new Epetra_MultiVector(*additionalcross2ele_));
  if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_myosinthick ||
     statmechparams_.get<double>("ACTIVELINKERFRACTION", 0.0)>0.0)
  {
    crosslinkeractlengthconv_ = Teuchos::rcp(new Epetra_Vector(*crosslinkeractlength_));
    crosslinkeractcycletimeconv_ = Teuchos::rcp(new Epetra_Vector(*crosslinkeractcycletime_));
  }

  //set addedelements_, deletedelements_ empty vectors
  addedelements_.clear();
  deletedelements_.clear();
  if(beamcmanager!=Teuchos::null)
  {
    addedcelements_.clear();
    deletedcelements_.clear();
  }

  return;
} // StatMechManager::WriteConv()

/*-----------------------------------------------------------------------*
 | (public) restore state at the beginning of this time step cyron 11/10 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::RestoreConv(Teuchos::RCP<Epetra_Vector>           dis,
                                            Teuchos::RCP<LINALG::SparseOperator>& stiff,
                                            Teuchos::RCP<CONTACT::Beam3cmanager>  beamcmanager,
                                            const bool                            printscreen)
{
  //restore state at the beginning of time step for relevant class variables
  crosslinkerbond_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerbondconv_));
  crosslinkerpositions_ = Teuchos::rcp(new Epetra_MultiVector(*crosslinkerpositionsconv_));
  bspotstatus_ = Teuchos::rcp(new Epetra_Vector(*bspotstatusconv_));
  numbond_ = Teuchos::rcp(new Epetra_Vector(*numbondconv_));
  crosslink2element_ = Teuchos::rcp(new Epetra_Vector(*crosslink2elementconv_));
  if(linkermodel_==statmech_linker_myosinthick)
    additionalcross2ele_ = Teuchos::rcp(new Epetra_MultiVector(*additionalcross2eleconv_));
  if(linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_myosinthick)
  {
    Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::null;
    Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::null;

    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"CROSSBRIDGEMODEL"))
    {
      Teuchos::RCP<Epetra_Vector> discol = Teuchos::rcp(new Epetra_Vector(*discret_->DofColMap(), true));
      LINALG::Export(*dis, *discol);
      //  reverts reference lengths of elements back to the way they were before... (has to be done prior to restoring actlinklength_. Otherwise, fatal loss of information)
      bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
      bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
      GetBindingSpotPositions(*discol, bspotpositions, bspotrotations);
    }
    ChangeActiveLinkerLength(1e9, 1e9, crosslinkeractlength_, crosslinkeractlengthconv_, printscreen, true, bspotpositions, bspotrotations);
    crosslinkeractlength_ = Teuchos::rcp(new Epetra_Vector(*crosslinkeractlengthconv_));
    crosslinkeractcycletime_ = Teuchos::rcp(new Epetra_Vector(*crosslinkeractcycletimeconv_));
  }

  /*restore state of the discretization at the beginning of this time step; note that to this and
   *adding and deleting crosslinker element has to be undone exactly vice-versa compared to the way
   *it was done first in order to handle also those crosslinkers correctly added and deleted in one
   *and the same time step*/

  //loop through all elements deleted in this time step and restore them in the discretization
  for(int i=0; i<(int)deletedelements_.size(); i++)
  {
    std::vector<char> tmp;
    std::vector<char>::size_type position = 0;
    DRT::ParObject::ExtractfromPack(position,deletedelements_[i],tmp);
    DRT::ParObject* o = DRT::UTILS::Factory(tmp);
    DRT::Element* ele = dynamic_cast<DRT::Element*>(o);
    if (ele == NULL)
      dserror("Failed to build an element from the element data");
    discret_->AddElement(Teuchos::rcp(ele));
  }
  deletedelements_.clear();

  //loop through addedelements_, delete all these elements and then set addedelements_ an empty vector
  for(int i=0; i<(int)addedelements_.size(); i++)
    discret_->DeleteElement(addedelements_[i]);
  addedelements_.clear();

  /*settling administrative stuff in order to make the discretization ready for the next time step: synchronize
   *the Filled() state on all processors after having added or deleted elements by ChekcFilledGlobally(); then build
   *new element maps and call FillComplete(); finally Crs matrices stiff_ has to be deleted completely and made ready
   *for new assembly since their graph was changed*/
  discret_->CheckFilledGlobally();
  discret_->FillComplete(true, false, false);
  stiff->Reset();

  // same procedure for contact discretization
  if(beamcmanager!=Teuchos::null)
  {
    for(int i=0; i<(int)deletedcelements_.size(); i++)
    {
      std::vector<char> tmp;
      std::vector<char>::size_type position = 0;
      DRT::ParObject::ExtractfromPack(position,deletedcelements_[i],tmp);
      DRT::ParObject* o = DRT::UTILS::Factory(tmp);
      DRT::Element* ele = dynamic_cast<DRT::Element*>(o);
      if (ele == NULL)
        dserror("Failed to build an element from the element data");
      beamcmanager->BTSolDiscret().AddElement(Teuchos::rcp(ele));
    }
    deletedcelements_.clear();

    for(int i=0; i<(int)addedcelements_.size(); i++)
      beamcmanager->BTSolDiscret().DeleteElement(addedcelements_[i]);
    addedcelements_.clear();

    // contact discretization
    beamcmanager->BTSolDiscret().CheckFilledGlobally();
    beamcmanager->BTSolDiscret().FillComplete();
  }

  return;
} // StatMechManager::RestoreConv()

/*-----------------------------------------------------------------------*
 | communicate Vector to all Processors                    mueller 11/11 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CommunicateVector(Teuchos::RCP<Epetra_Vector> InVec,
                                                  Teuchos::RCP<Epetra_Vector> OutVec,
                                                  bool                        doexport,
                                                  bool                        doimport,
                                                  bool                        zerofy,
                                                  bool                        exportinsert)
{
  /* zerofy InVec at the beginning of each search except for Proc 0
   * for subsequent export and reimport. This way, we guarantee redundant information
   * on all processors. */

  const Epetra_BlockMap& InMap = InVec->Map();
  const Epetra_BlockMap& OutMap = OutVec->Map();

  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutMap, InMap);
  Epetra_Import importer(OutMap, InMap);
  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_->Comm().MyPID()!=0 && zerofy)
      OutVec->PutScalar(0.0);
    if(exportinsert)
      InVec->Export(*OutVec, exporter, Insert);
    else
      InVec->Export(*OutVec, exporter, Add);
  }
  if(doimport)
    OutVec->Import(*InVec,importer,Insert);
  return;
}

/*-----------------------------------------------------------------------*
 | communicate MultiVector to all Processors               mueller 11/11 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CommunicateMultiVector(Teuchos::RCP<Epetra_MultiVector> InVec,
                                                       Teuchos::RCP<Epetra_MultiVector> OutVec,
                                                       bool                             doexport,
                                                       bool                             doimport,
                                                       bool                             zerofy,
                                                       bool                             exportinsert)
{
  // first, export the values of OutVec on Proc 0 to InVecs of all participating processors
  Epetra_Export exporter(OutVec->Map(), InVec->Map());
  Epetra_Import importer(OutVec->Map(), InVec->Map());

  if(doexport)
  {
    // zero out all vectors which are not Proc 0. Then, export Proc 0 data to InVec map.
    if(discret_->Comm().MyPID()!=0 && zerofy)
      OutVec->PutScalar(0.0);
    if(exportinsert)
      InVec->Export(*OutVec, exporter, Insert);
    else
      InVec->Export(*OutVec, exporter, Add);
  }
  if(doimport)
    OutVec->Import(*InVec,importer,Insert);
  return;
}


/*----------------------------------------------------------------------*
 | (public) writing restart information for manager objects   cyron 12/08|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::WriteRestart(Teuchos::RCP<IO::DiscretizationWriter> output, double& dt)
{
  output->WriteInt("istart", istart_);
  output->WriteInt("unconvergedsteps", unconvergedsteps_);
  output->WriteDouble("starttimeoutput", starttimeoutput_);
  output->WriteDouble("endtoendref", endtoendref_);
  output->WriteInt("basisnodes", basisnodes_);
  output->WriteInt("outputfilenumber", outputfilenumber_);
  //note: beginold_,endold_,sumdispmiddle_ not considered; related methods not restartable
  output->WriteInt("basiselements", basiselements_);
  output->WriteDouble("sumsquareincpar", sumsquareincpar_);
  output->WriteDouble("sumsquareincort", sumsquareincort_);
  output->WriteDouble("sumrotmiddle", sumrotmiddle_);
  output->WriteDouble("sumsquareincmid", sumsquareincmid_);
  output->WriteDouble("sumsquareincrot", sumsquareincrot_);
  output->WriteDouble("timestepsize", dt);
  /*note: crosslinkermap_, transfermap_, ddcorrrowmap_, ddcorrcolmap_,
   * filamentnumber_, forcesensor_, not considered, because generated in constructor*/

  //note: Using WriteVector and ReadMultiVector requires unique map of MultiVector thus export/import for restart to/from row map
  Epetra_Export exporter(*bspotcolmap_,*bspotrowmap_);
  Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_,true));
  bspotstatusrow->Export(*bspotstatus_,exporter,Insert);
  output->WriteVector("bspotstatus",bspotstatusrow,IO::DiscretizationWriter::nodevector);

  output->WriteRedundantDoubleVector("startindex",startindex_);

  if(linkermodel_==statmech_linker_active || linkermodel_==statmech_linker_activeintpol || linkermodel_==statmech_linker_myosinthick || statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)>0.0)
  {
    WriteRestartRedundantMultivector(output,"crosslinkertype",crosslinkertype_);
    WriteRestartRedundantMultivector(output,"crosslinkeractlength",crosslinkeractlength_);
    WriteRestartRedundantMultivector(output,"crosslinkeractcycletime", crosslinkeractcycletime_);
  }
  WriteRestartRedundantMultivector(output,"crosslinkerbond",crosslinkerbond_);
  WriteRestartRedundantMultivector(output,"crosslinkerpositions",crosslinkerpositions_);
  WriteRestartRedundantMultivector(output,"numbond",numbond_);
  WriteRestartRedundantMultivector(output,"crosslink2element",crosslink2element_);
  if(linkermodel_==statmech_linker_myosinthick)
    WriteRestartRedundantMultivector(output,"additionalcross2ele",additionalcross2ele_);
  WriteRestartRedundantMultivector(output,"visualizepositions",visualizepositions_);

  if(crosslinkunbindingtimes_!=Teuchos::null)
    WriteRestartRedundantMultivector(output,"crosslinkunbindingtimes",crosslinkunbindingtimes_);

  return;
} // StatMechManager::WriteRestart()

/*----------------------------------------------------------------------------*
 | (public) write restart information for fully redundant   Epetra_Multivector|
 | with name "name"                                                cyron 11/10|
 *----------------------------------------------------------------------------*/
void STATMECH::StatMechManager::WriteRestartRedundantMultivector(Teuchos::RCP<IO::DiscretizationWriter> output, const std::string name, Teuchos::RCP<Epetra_MultiVector> multivector)
{
  //create stl vector to store information in multivector
  Teuchos::RCP<std::vector<double> > stlvector = Teuchos::rcp(new std::vector<double>);
  stlvector->resize(multivector->MyLength()*multivector->NumVectors());

  for (int i=0; i<multivector->MyLength(); i++)
    for (int j=0; j<multivector->NumVectors(); j++)
      (*stlvector)[i + j*multivector->MyLength()] = (*multivector)[j][i];

  //write information to output file; note that WriteRedundantDoubleVector is active on proc 0 only
  output->WriteRedundantDoubleVector(name,stlvector);

  return;
} // StatMechManager::WriteRestartRedundantMultivector()



/*----------------------------------------------------------------------*
 |read restart information for statistical mechanics (public)cyron 12/08|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ReadRestart(IO::DiscretizationReader& reader, double& dt)
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
  dt = reader.ReadDouble("timestepsize");
  /*note: crosslinkermap_, transfermap_, ddcorrrowmap_, ddcorrcolmap_,
   * filamentnumber_, forcesensor_, not considered, because generated in constructor*/

  //note: Using WriteVector and ReadMultiVector requires uniquen map of MultiVector thus export/import for restart to/from row map
  Epetra_Import importer(*bspotcolmap_,*bspotrowmap_);
  Teuchos::RCP<Epetra_Vector> bspotstatusrow = Teuchos::rcp(new Epetra_Vector(*bspotrowmap_),true);
  reader.ReadVector(bspotstatusrow,"bspotstatus");
  bspotstatus_->Import(*bspotstatusrow,importer,Insert);


  //Read redundant Epetra_Multivectors and STL vectors
  reader.ReadRedundantDoubleVector(startindex_,"startindex");

  ReadRestartRedundantMultivector(reader,"crosslinkerbond",crosslinkerbond_);
  ReadRestartRedundantMultivector(reader,"crosslinkerpositions",crosslinkerpositions_);
  ReadRestartRedundantMultivector(reader,"numbond",numbond_);
  ReadRestartRedundantMultivector(reader,"crosslink2element",crosslink2element_);
  if(linkermodel_==statmech_linker_myosinthick)
    ReadRestartRedundantMultivector(reader,"additionalcross2ele",additionalcross2ele_);
  ReadRestartRedundantMultivector(reader,"visualizepositions",visualizepositions_);
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
    ReadRestartRedundantMultivector(reader,"crosslinkunbindingtimes",crosslinkunbindingtimes_);

  if(linkermodel_==statmech_linker_active || linkermodel_==statmech_linker_activeintpol || linkermodel_==statmech_linker_myosinthick || statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)>0.0)
  {
    ReadRestartRedundantMultivector(reader,"crosslinkertype", crosslinkertype_);
    ReadRestartRedundantMultivector(reader,"crosslinkeractlength",crosslinkeractlength_);
    ReadRestartRedundantMultivector(reader,"crosslinkeractcycletime", crosslinkeractcycletime_);
  }

  return;
}// StatMechManager::ReadRestart()

/*----------------------------------------------------------------------------*
 | (public) read restart information for fully redundant Epetra_Multivector   |
 | with name "name"                                                cyron 11/10|
 *----------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ReadRestartRedundantMultivector(IO::DiscretizationReader& reader, const std::string name, Teuchos::RCP<Epetra_MultiVector> multivector)
{
    //we assume that information was stored like for a redundant stl vector
    Teuchos::RCP<std::vector<double> > stlvector = Teuchos::rcp(new std::vector<double>);
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

/*----------------------------------------------------------------------*
 | add parameters to given parameter list          (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::AddStatMechParamsTo(Teuchos::ParameterList& params, Teuchos::RCP<Epetra_MultiVector> randomnumbers)
{
  params.set("ETA",statmechparams_.get<double>("ETA",0.0));
  params.set("THERMALBATH",DRT::INPUT::IntegralValue<INPAR::STATMECH::ThermalBathType>(statmechparams_,"THERMALBATH"));
  params.set<int>("FRICTION_MODEL",DRT::INPUT::IntegralValue<INPAR::STATMECH::FrictionModel>(statmechparams_,"FRICTION_MODEL"));
  params.set<INPAR::STATMECH::DBCType>("DBCTYPE", DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE"));
  if(randomnumbers!=Teuchos::null)
    params.set("RandomNumbers",randomnumbers);
  params.set("SHEARAMPLITUDE",statmechparams_.get<double>("SHEARAMPLITUDE",0.0));
  // note CURVENUMBER and DBCDISPDIR follow the counting convention 1,2,3,... (not C++). They are internally decremented
  params.set("CURVENUMBER",statmechparams_.get<int>("CURVENUMBER",-1));
  params.set("DBCDISPDIR",statmechparams_.get<int>("DBCDISPDIR",-1));
  params.set("STARTTIMEACT",actiontime_->at(bctimeindex_));
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE")==INPAR::STATMECH::dbctype_affineshear)
  params.set("STARTTIMEACT",actiontime_->at(bctimeindex_));
  params.set("DELTA_T_NEW",timestepsizes_->at(bctimeindex_));
  params.set("PERIODLENGTH",GetPeriodLength());

  if(linkermodel_ == statmech_linker_bellseq || linkermodel_ == statmech_linker_bellseqintpol ||
     linkermodel_ == statmech_linker_active || linkermodel_ == statmech_linker_activeintpol ||
     linkermodel_ == statmech_linker_myosinthick || networktype_ == statmech_network_casimir ||
     DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_viscoelasticity)
    params.set<std::string>("internalforces","yes");

  // adds information on the used time integration
  params.set<std::string>("DYNAMICTYP", DRT::Problem::Instance()->StructuralDynamicParams().get<std::string>("DYNAMICTYP"));

  return;
}

/*----------------------------------------------------------------------*
 | Decomposition of vector vec into parallel and orthogonal part with   |
 | respect to a given vector refvec (private)             mueller (3/14)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::DoVectorDecomposition(LINALG::Matrix<3,1>& vec,
                                                      LINALG::Matrix<3,1>& parvec,
                                                      LINALG::Matrix<3,1>& orthovec,
                                                      LINALG::Matrix<3,1>& refvec)
{
  parvec = refvec;
  parvec.Scale(vec.Dot(refvec));
  orthovec = vec;
  orthovec -= parvec;
  return;
}

/*----------------------------------------------------------------------*
 | rotation around a fixed axis by angle phirot          mueller (11/11)|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::RotationAroundFixedAxis(LINALG::Matrix<3,1>& axis, LINALG::Matrix<3,1>& vector, const double& phirot)
{
  // unit length
  axis.Scale(1.0/axis.Norm2());
  vector.Scale(1.0/vector.Norm2());

  LINALG::Matrix<3, 3> RotMat;

  for (int j=0; j<3; j++)
    RotMat(j, j) = cos(phirot) + axis(j) * axis(j) * (1 - cos(phirot));
  RotMat(0, 1) = axis(0) * axis(1) * (1 - cos(phirot)) - axis(2) * sin(phirot);
  RotMat(0, 2) = axis(0) * axis(2) * (1 - cos(phirot)) + axis(1) * sin(phirot);
  RotMat(1, 0) = axis(1) * axis(0) * (1 - cos(phirot)) + axis(2) * sin(phirot);
  RotMat(1, 2) = axis(1) * axis(2) * (1 - cos(phirot)) - axis(0) * sin(phirot);
  RotMat(2, 0) = axis(2) * axis(0) * (1 - cos(phirot)) - axis(1) * sin(phirot);
  RotMat(2, 1) = axis(2) * axis(1) * (1 - cos(phirot)) + axis(0) * sin(phirot);

  LINALG::Matrix<3,1> tmpvec = vector;
  vector.Multiply(RotMat, tmpvec);

  return;
}

/*----------------------------------------------------------------------*
 | Shifts current position of nodes so that they comply with periodic   |
 | boundary conditions                                       cyron 04/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PeriodicBoundaryShift(Epetra_Vector& disrow,
                                                      int            ndim,
                                                      const double   timen,
                                                      const double   dt)
{
  double starttime = actiontime_->at((int)(actiontime_->size()-1));
  double shearamplitude = statmechparams_.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = statmechparams_.get<int> ("CURVENUMBER", -1)-1;
  int oscilldir = statmechparams_.get<int> ("DBCDISPDIR", -1)-1;


  //only if period length >0 has been defined periodic boundary conditions are swithced on
  if (periodlength_->at(0) > 0.0)
  {
    for (int i=0; i<discret_->NumMyRowNodes(); i++)
    {
      //get a pointer at i-th row node
      const DRT::Node* node = discret_->lRowNode(i);

      //get GIDs of this node's degrees of freedom
      std::vector<int> dofnode = discret_->Dof(node);

      for (int j=ndim-1; j>-1; j--)
      {
        //current coordinate value
        double xcurr = node->X()[j] + disrow[discret_->DofRowMap()->LID(dofnode[j])];

        /*if node currently has coordinate value greater than periodlength,
         *it is shifted by -periodlength sufficiently often to lie again in the domain*/
        if (xcurr > (*periodlength_)[j])
        {
          disrow[discret_->DofRowMap()->LID(dofnode[j])] -= (*periodlength_)[j]*floor(xcurr/(*periodlength_)[j]);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may fixed by DBC. To avoid problems when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if (j == 2 && curvenumber >= 0 && timen > starttime && fabs(timen-starttime)>dt/1e4)
            disrow[discret_->DofRowMap()->LID(dofnode[oscilldir])] -= shearamplitude * DRT::Problem::Instance()->Curve(curvenumber).f(timen);
        }
        /*if node currently has coordinate value smaller than zero, it is shifted by periodlength sufficiently often
         *to lie again in the domain*/
        if (node->X()[j] + disrow[discret_->DofRowMap()->LID(dofnode[j])] < 0.0)
        {
          disrow[discret_->DofRowMap()->LID(dofnode[j])] -= (*periodlength_)[j]*floor(xcurr/(*periodlength_)[j]);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problems when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if (j == 2 && curvenumber >= 0 && timen > starttime && fabs(timen-starttime)>dt/1e4)
            disrow[discret_->DofRowMap()->LID(dofnode[oscilldir])] += shearamplitude * DRT::Problem::Instance()->Curve(curvenumber).f(timen);
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Periodic Boundary Shift for crosslinker diffusion simulation        |
 |                                                (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::CrosslinkerPeriodicBoundaryShift(Teuchos::RCP<Epetra_MultiVector> crosslinkerpositions)
{
  for (int i=0; i<crosslinkerpositions->NumVectors(); i++)
    for (int j=0; j<crosslinkerpositions->MyLength(); j++)
    {
      if ((*crosslinkerpositions)[i][j] > (*periodlength_)[i])
        (*crosslinkerpositions)[i][j] -= (floor((*crosslinkerpositions)[i][j]/(*periodlength_)[i]))*(*periodlength_)[i];
      if ((*crosslinkerpositions)[i][j] < 0.0)
        (*crosslinkerpositions)[i][j] -= (floor((*crosslinkerpositions)[i][j]/(*periodlength_)[i]))*(*periodlength_)[i];
    }
  return;
}// StatMechManager::CrosslinkerPeriodicBoundaryShift

/*----------------------------------------------------------------------*
 | check for broken element                        (public) mueller 3/10|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::CheckForBrokenElement(LINALG::SerialDenseMatrix& coord, LINALG::SerialDenseMatrix& cut)
{
  // empty cut just in case it was handed over non-empty
  cut.Zero();
  bool broken = false;
  int ndim = coord.M();
  // flag "broken" signals a broken element, flag "cut" hints at location of nodes 0 and 1 (of the two-noded beam)
  for (int dof = 0; dof < ndim; dof++)
    //loop through columns of "cut"
    for (int n = 0; n < cut.N(); n++)
    {
      // broken element with node_n close to "0.0"-boundary
      if (fabs(coord(dof, n + 1) - periodlength_->at(dof) - coord(dof, n)) < fabs(coord(dof, n + 1) - coord(dof, n)))
      {
        broken = true;
        // set value for the spatial component in question at n-th cut
        cut(dof, n) = 1.0;
      }
      else if (fabs(coord(dof, n + 1) + periodlength_->at(dof) - coord(dof, n)) < fabs(coord(dof, n + 1) - coord(dof, n)))
      {
        broken = true;
        cut(dof, n) = 2.0;
      }
    }
  return broken;
}// StatMechManager::CheckForBrokenElement

/*----------------------------------------------------------------------*
 | check for crosslink that reaches across periodic BCs     mueller 1/14|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::CheckCrossPeriodicBCLink(const Epetra_SerialDenseMatrix& LID,
                                                         const Epetra_Vector&            discol)
{
  bool crossboundarylink = false;

  switch(linkermodel_)
  {
    case statmech_linker_stdintpol:
    case statmech_linker_activeintpol:
    case statmech_linker_bellseqintpol:
    case statmech_linker_myosinthick:
    {

      switch(DRT::INPUT::IntegralValue<INPAR::STATMECH::DBCType>(statmechparams_,"DBCTYPE"))
      {
        case INPAR::STATMECH::dbctype_shearfixed:
        case INPAR::STATMECH::dbctype_affineshear:
        {
          // binding spot loop
          for(int k=0; k<LID.M(); k++)
          {
            LINALG::SerialDenseMatrix coord(3,2);
            LINALG::SerialDenseMatrix cut(3,1);
            std::vector<int> nodegids(coord.N());
            // node loop
            for(int l=0; l<coord.N(); l++)
            {
              nodegids[l] = (int)(*bspot2nodes_)[l][(int)LID(k,0)];
              DRT::Node* node = discret_->lColNode(discret_->NodeColMap()->LID(nodegids[l]));
              for(int m=0; m<coord.M(); m++)
              {
                int dofgid = discret_->Dof(0, node)[m];
                coord(m,l) = node->X()[m]+discol[dofgid];
              }
            }
            // if filament element considered for linking reaches across the z-boundary
            if(CheckForBrokenElement(coord,cut) && cut(2,0)>0)
              crossboundarylink = true;
          }
        }
        break;
        default: { }
      }
    }
    break;
    default: { }
  }
  return crossboundarylink;
}

/*----------------------------------------------------------------------------------------------*
 | Shift vector in accordance to periodic BCs                                    mueller 10/12  |
 *----------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ShiftPositions(std::vector<double>& pos,
                                               const int&           nnode)
{
  for (int i=0; i<nnode; i++)
  {

    for (int k = 0; k < 3; k++)
    {
      // in case of positive overlap of component
      if (pos[3*i+k] > periodlength_->at(k))
          pos[3*i+k] -=  periodlength_->at(k);
      // in case of negative overlap of component
      if (pos[3*i+k] < 0)
          pos[3*i+k] += periodlength_->at(k);
    }
  }
  return;
}


/*----------------------------------------------------------------------------------------------*
 | Reverse the effect of periodic BCs on a vector                                mueller 10/12  |
 *----------------------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::UnshiftPositions(std::vector<double>& pos,
                                                 const int&           nnode,
                                                 bool                 separatefil)
{
  if (periodlength_->at(0) > 0.0)
  {
    if((int)pos.size()%nnode!=0)
      dserror("Check your position vector! It does not have an equal number of entries per node!");
    int ndim = (int)pos.size()/nnode;

    if(separatefil)
    {
      if(nnode==4)
      {
        int nfil = 2;
        for(int h=0; h<nfil; h++)
        {
          for(int i=1; i<nnode/2; i++)
          {
            for(int dof= ndim - 1; dof > -1; dof--)
            {
              if( fabs( pos[nfil*ndim*h+ndim*i+dof] + periodlength_->at(dof) - pos[nfil*ndim*h+ndim*(i-1)+dof] ) < fabs( pos[nfil*ndim*h+ndim*i+dof] - pos[nfil*ndim*h+ndim*(i-1)+dof] ) )
                pos[nfil*ndim*h+ndim*i+dof] += periodlength_->at(dof);
              if( fabs( pos[nfil*ndim*h+ndim*i+dof] - periodlength_->at(dof) - pos[nfil*ndim*h+ndim*(i-1)+dof]) < fabs( pos[nfil*ndim*h+ndim*i+dof] - pos[nfil*ndim*h+ndim*(i-1)+dof] ) )
                pos[nfil*ndim*h+ndim*i+dof] -= periodlength_->at(dof);
            }
          }
        }
      }
      else
        dserror("Only 4-noded element implemented!");
    }
    else
    {
      for(int i=1; i<nnode; i++)
      {
        for(int dof= ndim - 1; dof > -1; dof--)
        {
          if( fabs( pos[3*i+dof] + periodlength_->at(dof) - pos[3*(i-1)+dof] ) < fabs( pos[3*i+dof] - pos[3*(i-1)+dof] ) )
            pos[3*i+dof] += periodlength_->at(dof);
          if( fabs( pos[3*i+dof] - periodlength_->at(dof) - pos[3*(i-1)+dof]) < fabs( pos[3*i+dof] - pos[3*(i-1)+dof] ) )
            pos[3*i+dof] -= periodlength_->at(dof);
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | Computes current internal energy of discret_ (public)     cyron 12/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ComputeInternalEnergy(const Teuchos::RCP<Epetra_Vector> dis,
                                                      std::vector<double>&              energy,
                                                      const double&                     dt,
                                                      const std::ostringstream&         filename,
                                                      bool                              fillzeros,
                                                      bool                              writefile)
{
  Teuchos::ParameterList p;
  p.set("action", "calc_struct_energy");

  //add statistical vector to parameter list for statistical forces and damping matrix computation
  p.set("delta time",dt);
  p.set("ETA",statmechparams_.get<double>("ETA",0.0));
  p.set("THERMALBATH",DRT::INPUT::IntegralValue<INPAR::STATMECH::ThermalBathType>(statmechparams_,"THERMALBATH"));
  p.set<int>("FRICTION_MODEL",DRT::INPUT::IntegralValue<INPAR::STATMECH::FrictionModel>(statmechparams_,"FRICTION_MODEL"));
  p.set("SHEARAMPLITUDE",statmechparams_.get<double>("SHEARAMPLITUDE",0.0));
  p.set("CURVENUMBER",statmechparams_.get<int>("CURVENUMBER",-1));
  p.set("DBCDISPDIR",statmechparams_.get<int>("DBCDISPDIR",-1));
  p.set("PERIODLENGTH", periodlength_);

  discret_->ClearState();
  discret_->SetState("displacement", dis);
  Teuchos::RCP<Epetra_SerialDenseVector> energies = Teuchos::rcp(new Epetra_SerialDenseVector(6));
  energies->Scale(0.0);
  //filaments
  p.set("energyoftype", "beam3r");
  discret_->EvaluateScalars(p, energies);
  for(int i=0; i<energies->M(); i++)
    energy.push_back((*energies)(i));
  //crosslinkers
  energies->Scale(0.0);
  p.remove("energyoftype", true);
  if(linkermodel_ == statmech_linker_stdintpol  ||
     linkermodel_ == statmech_linker_activeintpol ||
     linkermodel_ == statmech_linker_bellseqintpol ||
     linkermodel_ == statmech_linker_myosinthick)
    p.set("energyoftype", "beam3cl");
  else
    p.set("energyoftype", "beam3");
  discret_->EvaluateScalars(p, energies);
  for(int i=0; i<energies->M(); i++)
    energy.push_back((*energies)(i));
  discret_->ClearState();

  if(!discret_->Comm().MyPID() && writefile)
  {
    FILE* fp = NULL;
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream writetofile;
    for(int i=0; i<(int)energy.size(); i++)
      writetofile<<energy[i]<<"    ";
    if(fillzeros)
      for(int i=0; i<17-(int)energy.size(); i++)
        writetofile<<"    "<<0;
    writetofile<<std::endl;
    fprintf(fp, writetofile.str().c_str());
    fclose(fp);
  }
  return;
} // StatMechManager::ComputeInternalEnergy

/*----------------------------------------------------------------------*
 | generates a vector with a random permutation of the numbers between 0|
 | and N - 1                                        (public) cyron 06/10|
 *----------------------------------------------------------------------*/
std::vector<int> STATMECH::StatMechManager::Permutation(const int& N)
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
    j = (int)floor((i + 1.0)*(*uniformgen_)());

    /*exchange values at positions i and j (note: value at position i is i due to above initialization
     *and because so far only positions <=i have been changed*/
    result[i] = result[j];
    result[j] = i;
  }

  return result;
} // StatMechManager::Permutation

/*----------------------------------------------------------------------*
 | get a matrix with node coordinates and their DOF LIDs                |
 |                                                 (public) mueller 3/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetElementNodeCoords(DRT::Element* element, Teuchos::RCP<Epetra_Vector> discol, LINALG::SerialDenseMatrix& coord, std::vector<int>* lids)
{
  // clear LID vector just in case it was handed over non-empty
  lids->clear();
  for (int j=0; j<element->NumNode(); j++)
  {
    int nodegid = element->NodeIds()[j];
    DRT::Node* node = discret_->lColNode(discret_->NodeColMap()->LID(nodegid));
    for (int k=0; k<3; k++)
    {
      // obtain k-th spatial component of the reference position of the j-th node
      double referenceposition = node->X()[k];
      // get the GIDs of the node's DOFs
      std::vector<int> dofnode = discret_->Dof(node);
      // store the displacement of the k-th spatial component
      double displacement = (*discol)[discret_->DofColMap()->LID(dofnode[k])];
      // write updated components into coord (only translational
      coord(k, j) = referenceposition + displacement;
      // store current lid(s) (3 translational DOFs per node)
      if (lids != NULL)
        lids->push_back(discret_->DofRowMap()->LID(dofnode[k]));
    }
  }
  return;
} // StatMechManager::GetElementNodeCoords


/*--------------------------------------------------------------------*
 | Return the inclusive angle between two given vectors               |
 |                                             (public) mukherjee 3/14|
 *--------------------------------------------------------------------*/
double STATMECH::StatMechManager::GetTheta(LINALG::Matrix<3,1> &T1, LINALG::Matrix<3,1> &T2)
{
  // Calculate current angle
  double dotprod=0.0;
  LINALG::Matrix  <1,3> crossprod(true);
  double CosTheta=0.0;
  double SinTheta=0.0;

  double norm_T1 = T1.Norm2();
  double norm_T2=T2.Norm2();
  for (int j=0; j<3; ++j)
    dotprod +=  T1(j) * T2(j);

  CosTheta = dotprod/(norm_T1*norm_T2);

  // cross product
  crossprod(0) = T1(1)*T2(2) - T1(2)*T2(1);
  crossprod(1) = T1(2)*T2(0) - T1(0)*T2(2);
  crossprod(2) = T1(0)*T2(1) - T1(1)*T2(0);

  double norm= crossprod.Norm2();
  SinTheta= norm/(norm_T1*norm_T2);

  double ThetaBoundary1=M_PI/4;
  double ThetaBoundary2=3*M_PI/4;
  double theta=0;
  if (SinTheta >=0)
  {
    if (CosTheta >= cos(ThetaBoundary1))
      theta=asin(SinTheta);
    else if (CosTheta <= cos(ThetaBoundary2))
      theta=M_PI-asin(SinTheta);
    else
      theta=acos(CosTheta);
  }
  else
    dserror("Angle more than 180 degrees!");

  return theta;
} // StatMechManager::GetTheta
