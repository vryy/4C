/*----------------------------------------------------------------------*/
/*!
\file regularization_totalvariation.cpp

\brief Total variation type regularization

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_totalvariation.H"

#include "matpar_manager.H" //for ConnectivityData
#include "invana_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"

#include "Teuchos_ParameterList.hpp"
#include "Epetra_CombineMode.h"


/*----------------------------------------------------------------------*/
INVANA::RegularizationTotalVariation::RegularizationTotalVariation() :
  RegularizationBase(),
  eps_(0.0)
{}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTotalVariation::Setup(const Teuchos::ParameterList& invp)
{
  adjacency_ = connectivity_->AdjacencyMatrix();
  eps_ = invp.get<double>("TVD_EPS");
  weight_ = invp.get<double>("REG_WEIGHT");
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTotalVariation::Evaluate(const Epetra_MultiVector& theta, double* value)
{
  // check if the map of theta is the row map of the adjacency matrix
  if (not theta.Map().SameAs(adjacency_->RowMap()))
    dserror("Map of the parameter vector is not the row map of the adjacency matrix");

  // check if only one vector in the multivector
  if (theta.NumVectors()>1)
    dserror("tvd regularization not made for use with more than one vector!");

  // communicate theta data from other procs such that every proc can compute sums
  // over adjacent parameters
  Epetra_MultiVector thetacol(adjacency_->ColMap(),theta.NumVectors(),false);
  LINALG::Export(theta,thetacol);

  double functvalue=0.0;
  for (int i=0; i<theta(0)->MyLength(); i++)
  {
    // get weights of neighbouring parameters
    int lenindices = adjacency_->NumMyEntries(i);
    int numindex;
    std::vector<int> indices(lenindices,0);
    std::vector<double> weights(lenindices,0);
    adjacency_->ExtractMyRowCopy(i,lenindices,numindex,&weights[0],&indices[0]);

    double rowsum=0.0;
    double rowval=(*theta(0))[i];
    for (int j=0; j<lenindices; j++)
    {
      if (indices[j]!=i) // skip substracting from itself
      {
        double colval=(*thetacol(0))[indices[j]];
        rowsum += weights[j]*(colval-rowval)*(colval-rowval);
      }
    }

    // sum over all the rows
    functvalue+=sqrt(rowsum+eps_);
  }

  // collect contributions from all procs
  double gfunctvalue=0.0;
  adjacency_->Comm().SumAll(&functvalue,&gfunctvalue,1);

  *value +=weight_*gfunctvalue;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTotalVariation::EvaluateGradient(const Epetra_MultiVector& theta,
    Teuchos::RCP<Epetra_MultiVector> gradient)
{
  // check if the map of theta is the row map of the adjacency matrix
  if (not theta.Map().SameAs(adjacency_->RowMap()))
    dserror("Map of the parameter vector is not the row map of the adjacency matrix");

  // check if only one vector in the multivector
  if (theta.NumVectors()>1 or gradient->NumVectors()>1)
    dserror("tvd regularization not made for use with more than one vector!");

  // communicate theta data from other procs such that every proc can compute sums
  // over adjacent parameters
  Epetra_MultiVector thetacol(adjacency_->ColMap(),theta.NumVectors(),false);
  LINALG::Export(theta,thetacol);

  // a gradient with the column layout of the adjacency matrix to be filled
  // on each proc and then summed up via the Epetra_Export back to the
  // unique layout of the gradient vector coming in
  Epetra_MultiVector gradientcol(adjacency_->ColMap(), gradient->NumVectors(),true);

  for (int i=0; i<theta(0)->MyLength(); i++)
  {
    // get weights of neighbouring parameters
    int lenindices = adjacency_->NumMyEntries(i);
    int numindex;
    std::vector<int> indices(lenindices,0);
    std::vector<double> weights(lenindices,0);
    adjacency_->ExtractMyRowCopy(i,lenindices,numindex,&weights[0],&indices[0]);

    // value of the parameter in this row
    double rowval=(*theta(0))[i];

    // denominator
    double denom=0.0;
    for (int j=0; j<lenindices; j++)
    {
      if (indices[j]!=i) // skip substracting from itself
      {
        double colval=(*thetacol(0))[indices[j]];
        denom += weights[j]*(colval-rowval)*(colval-rowval);
      }
    }
    denom = sqrt(denom+eps_);

    // nominator i
    double nomi=0.0;
    for (int j=0; j<lenindices; j++)
    {
      if (indices[j]!=i) // skip substracting from itself
      {
        double colval=(*thetacol(0))[indices[j]];
        nomi += weights[j]*(colval-rowval);
      }
    }
    nomi=nomi*(-1);

    // insert contributions into dof i of the gradient
    gradientcol.SumIntoMyValue(i,0,nomi/denom);

    // nominator j
    double nomj=0.0;
    for (int j=0; j<lenindices; j++)
    {
      if (indices[j]!=i) // skip substracting from itself
      {
        double colval=(*thetacol(0))[indices[j]];
        nomj = weights[j]*(colval-rowval);

        // insert contributions into dof j of the gradient
        gradientcol.SumIntoMyValue(indices[j],0,nomj/denom);
      }
    }

  }// loop theta

  // bring back gradient to the unique layout he came here with.
  // we have to add up off proc components so we cannot use
  // the LINALG::Export since it only provides CombineMode::Insert
  Epetra_MultiVector tmp(gradient->Map(), gradient->NumVectors(),true);
  Epetra_Export exporter(gradientcol.Map(), tmp.Map());
  int err = tmp.Export(gradientcol, exporter, Add);
  if (err)
    dserror("Export using exporter returned err=%d", err);

  gradient->Update(weight_,tmp,1.0);

  return;
}

