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

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::RegularizationTotalVariation::Evaluate(const Epetra_MultiVector& theta, double* value)
{
  // check if only one vector in the multivector
  if (theta.NumVectors()>1)
    dserror("tv regularization not made for use with more than one vector!");

  // project to elementwise layout
  Epetra_Vector thetaele(adjacency_->RowMap(),true);
  connectivity_->Projector()->Multiply(true,theta,thetaele);

  // check if the map of theta is the row map of the adjacency matrix
  if (not thetaele.Map().SameAs(adjacency_->RowMap()))
    dserror("Map of the parameter vector is not the row map of the adjacency matrix");

  // communicate theta data from other procs such that every proc can compute sums
  // over adjacent parameters
  Epetra_MultiVector thetacol(adjacency_->ColMap(),thetaele.NumVectors(),false);
  LINALG::Export(thetaele,thetacol);

  double functvalue=0.0;
  for (int i=0; i<thetaele(0)->MyLength(); i++)
  {
    // get weights of neighbouring parameters
    int lenindices = adjacency_->NumMyEntries(i);
    int numindex;
    std::vector<int> indices(lenindices,0);
    std::vector<double> weights(lenindices,0);
    adjacency_->ExtractMyRowCopy(i,lenindices,numindex,&weights[0],&indices[0]);

    // row in local index space of the collayout
    int rowi = thetacol.Map().LID(thetaele.Map().GID(i));

    double rowsum=0.0;
    double rowval=(*thetacol(0))[i];
    for (int j=0; j<lenindices; j++)
    {
      if (indices[j]!=rowi) // skip substracting from itself
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
  // check if only one vector in the multivector
  if (theta.NumVectors()>1 or gradient->NumVectors()>1)
    dserror("tv regularization not made for use with more than one vector!");

  // project to elementwise layout
  Epetra_Vector thetaele(adjacency_->RowMap(),true);
  connectivity_->Projector()->Multiply(true,theta,thetaele);

  // check if the map of theta is the row map of the adjacency matrix
  if (not thetaele.Map().SameAs(adjacency_->RowMap()))
    dserror("Map of the parameter vector is not the row map of the adjacency matrix");


  // communicate theta data from other procs such that every proc can compute sums
  // over adjacent parameters
  Epetra_MultiVector thetacol(adjacency_->ColMap(),thetaele.NumVectors(),false);
  LINALG::Export(thetaele,thetacol);

  // a gradient with the column layout of the adjacency matrix to be filled
  // on each proc and then summed up via the Epetra_Export back to the
  // unique layout of the gradient vector coming in
  Epetra_MultiVector gradientcol(adjacency_->ColMap(), gradient->NumVectors(),true);

  for (int i=0; i<thetaele(0)->MyLength(); i++)
  {
    // get weights of neighbouring parameters
    int lenindices = adjacency_->NumMyEntries(i);
    int numindex;
    std::vector<int> indices(lenindices,0);
    std::vector<double> weights(lenindices,0);
    adjacency_->ExtractMyRowCopy(i,lenindices,numindex,&weights[0],&indices[0]);

    if (numindex!=lenindices)
      dserror("bla");

    // value of the parameter in this row
    double rowval=(*thetaele(0))[i];

    // row in local index space of the collayout
    int rowi = thetacol.Map().LID(thetaele.Map().GID(i));

    // denominator
    double denom=0.0;
    for (int j=0; j<lenindices; j++)
    {
      if (indices[j]!=rowi) // skip substracting from itself
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
      if (indices[j]!=rowi) // skip substracting from itself
      {
        double colval=(*thetacol(0))[indices[j]];
        nomi += weights[j]*(colval-rowval);
      }
    }
    nomi=nomi*(-1);

    // insert contributions into dof rowi of the gradient
    gradientcol.SumIntoMyValue(rowi,0,nomi/denom);

    // nominator j
    double nomj=0.0;
    for (int j=0; j<lenindices; j++)
    {
      if (indices[j]!=rowi) // skip substracting from itself
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

  // apply possible linear transformation
  // project to elementwise layout
  Epetra_Vector tmptrans(connectivity_->Projector()->RangeMap(),true);
  connectivity_->Projector()->Multiply(true,tmp,tmptrans);

  gradient->Update(weight_,tmptrans,1.0);

  return;
}

