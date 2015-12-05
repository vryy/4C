/*!----------------------------------------------------------------------
\file statmech_bilayer_utils.cpp
\brief management and auxiliary functions for statistical mechanics
       of lipid bilayer

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/
#include "statmech_manager_bilayer.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_discsh3/discsh3.H"


/*-----------------------------------------------------------------------*
 | communicate MultiVector to all Processors             mukherjee 10/15 |
 *-----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::CommunicateMultiVector(Teuchos::RCP<Epetra_MultiVector> InVec,
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
 | Get Spatial postion of the node              (public) mukherjee 09/15|
 *----------------------------------------------------------------------*/
LINALG::Matrix <1,3> STATMECH::StatMechManagerBilayer::GetSpatialPosition(DRT::Node * node,
                                                          Epetra_Vector&  discol)
{
    LINALG::Matrix <1,3> current_pos(true);
    for (int dim=0;dim<3;++dim)
    {
        double refpos = node->X()[dim];
        std::vector<int> dofnode = discret_->Dof(node);
        double displacement = (double)(discol)[discret_->DofColMap()->LID(dofnode[dim])];
        current_pos(dim) =  refpos + displacement;
    }

  return current_pos;
}

/*------------------------------------------------------------------------*
 | (public) generate gaussian randomnumbers with mean "meanvalue" and     |
 | standarddeviation "standarddeviation" for parallel use  mukherjee 10/15|
 *------------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::GenerateGaussianRandomNumbers(Teuchos::RCP<Epetra_MultiVector> randomnumbers, const double meanvalue, const double standarddeviation)
{
  randomnumbers->PutScalar(0.0);

  //multivector for stochastic forces evaluated by each element based on row map
  Teuchos::RCP<Epetra_MultiVector> randomnumbersrow = Teuchos::rcp(new Epetra_MultiVector(*(discret_->ElementRowMap()), randomnumbers->NumVectors()));

  //MAXRANDFORCE is a multiple of the standard deviation
  double maxrandforcefac = statmechBilayerparams_.get<double>("MAXRANDFORCE",-1.0);

  if(maxrandforcefac==-1.0)
  {
    for (int i=0; i<randomnumbersrow->MyLength(); i++)
      for (int j=0; j<randomnumbersrow->NumVectors(); j++)
        (*randomnumbersrow)[j][i] = standarddeviation*(*normalgen_)() + meanvalue;
  }
  else
  {
    for (int i=0; i<randomnumbersrow->MyLength(); i++)
      for (int j=0; j<randomnumbersrow->NumVectors(); j++)
      {
        (*randomnumbersrow)[j][i] = standarddeviation*(*normalgen_)() + meanvalue;
        if((*randomnumbersrow)[j][i]>maxrandforcefac*standarddeviation)
        {
          (*randomnumbersrow)[j][i]=maxrandforcefac*standarddeviation;
        }
        else if((*randomnumbersrow)[j][i]<-maxrandforcefac*standarddeviation)
        {
          (*randomnumbersrow)[j][i]=-maxrandforcefac*standarddeviation;
        }
      }
  }

  //export stochastic forces from row map to column map (unusual CommunicateMultiVector() call but does the job!)
  CommunicateMultiVector(randomnumbers,randomnumbersrow,true,false,false);

#ifdef DEBUGCOUT
  for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
  {
    if(pid==discret_->Comm().MyPID())
      std::cout<<"\n\nProc "<<pid<<": Row\n\n"<<*randomnumbersrow<<std::endl;
    discret_->Comm().Barrier();
  }

  discret_->Comm().Barrier();

  for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
  {
    if(pid==discret_->Comm().MyPID())
      std::cout<<"\n\nProc "<<pid<<": Column\n\n"<<*randomnumbers<<std::endl;
    discret_->Comm().Barrier();
  }
#endif

  return;
} // StatMechManagerBilayer::SynchronizeRandomForces()

/*----------------------------------------------------------------------*
 | seed all random generators of this object with fixed seed if given and|
 | with system time otherwise; seedparameter is used only in the first   |
 | case to calculate the actual seed variable based on some given fixed  |
 | seed value; note that seedparameter may be any integer, but has to be |
 | been set in a deterministic way so that it for a certain call of this |
 | method at a certain point in the program always the same number       |
 | whenever the program is used                               cyron 11/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::SeedRandomGenerators(const int seedparameter, const int seedparameter2)
{
  //integer for seeding all random generators
  int seedvariable = 0;

  Teuchos::ParameterList statmechBilayerparams = GetStatMechBilayerParams();
  double randnumtimeinc = statmechBilayerparams.get<double>("RANDNUMTIMEINT",-1.0);

  /*if input flag FIXEDSEED == YES: use same random numbers in each program start;
   *to this end compute seedvariable from given parameter FIXEDSEED and some other
   *deterministic parameter seedparameter given to this method at runtime*/
  if(DRT::INPUT::IntegralValue<int>(statmechBilayerparams_,"FIXEDSEED"))
  {
    //Decide if random numbers should change in every time step...
    if(randnumtimeinc==-1.0)
      seedvariable = (statmechBilayerparams_.get<int>("INITIALSEED", 0) + seedparameter)*(discret_->Comm().MyPID() + 1);
    //...or not before a prescribed interval RANDNUMTIMEINT
    else
      seedvariable = (statmechBilayerparams_.get<int>("INITIALSEED", 0) + seedparameter2)*(discret_->Comm().MyPID() + 1);

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
} // StatMechManagerBilayer::SeedRandomGenerators()

