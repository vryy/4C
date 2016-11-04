/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_uniform.cpp

\brief Patch-wise uniform distribution of material parameters

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager_uniform.H"

#include "Epetra_Import.h"

#include "invana_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_element.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"

INVANA::MatParManagerUniform::MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManager(discret)
{}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::Setup()
{
  paramlayoutmap_ = Teuchos::rcp(new Epetra_Map(NumParams(),NumParams(),0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));
  int numeleperproc=0;
  if (Discret()->Comm().MyPID()==0) numeleperproc = NumParams();
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(-1,numeleperproc,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));

  // build the mapextractor
  // the partial maps (only one map in case of uniform material parameters)
  std::vector< Teuchos::RCP<const Epetra_Map> > partials;
  partials.push_back(paramlayoutmapunique_);
  paramapextractor_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*paramlayoutmapunique_,partials));

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_,1,true));
  optparams_initial_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_,1,true));

  //initialize parameter vector from material parameters given in the input file
  InitParams();
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  // make optparams redundant on every proc
  Teuchos::RCP<Epetra_MultiVector> optparams = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,1,false));
  Epetra_Import importer(optparams->Map(), optparams_->Map());
  int err = optparams->Import(*optparams_, importer, Insert);
  if (err)
    dserror("Export using exporter returned err=%d", err);

  // now every proc can fill in the material parameters
  for (int i=0; i<NumParams(); i++)
    (*params)(i)->PutScalar((*(*optparams)(0))[i]);
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::InitParameters(int parapos, double val)
{
  // only proc 0 is filled here from the input file material
  optparams_->ReplaceGlobalValue(parapos,0,val);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
    double val, int elepos, int paraposglobal,int paraposlocal)
{
  // every proc can do the 'product rule' on his own
  // summation can be done after the assembly is complete
  dfint->SumIntoGlobalValue(paraposglobal,0,val);

  return;

}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::Finalize(Teuchos::RCP<Epetra_MultiVector> source,
    Teuchos::RCP<Epetra_MultiVector> target)
{
  // sum across processor
  std::vector<double> val(source->MyLength(),0.0);
  Discret()->Comm().SumAll((*source)(0)->Values(),&val[0],source->MyLength());

  // put into the global gradient; procs who dont own
  // the gid dont get contributions
  for (int i=0; i<NumParams(); i++)
    target->SumIntoGlobalValue(i,0,val[i]);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerUniform::FillAdjacencyMatrix(const Epetra_Map& paramrowmap,
    Teuchos::RCP<Epetra_CrsMatrix> graph)
{
  /* there is no particular connectivity for uniform distributions.
   if you want a connectivity between the different patches, do it.
   But for the moment the graph has just a diagonal with a value which is
   zeroed out in the calling routine anyways. So the graph will be left
   just empty here. */

  return;
}
