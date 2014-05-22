/*!----------------------------------------------------------------------
\file drt_utils.cpp
\brief A collection of helper methods for namespace DRT

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/


#include <algorithm>
#include <numeric>
#include <vector>
#include <set>
#include <map>
#include <string>

#include "drt_utils.H"
#include "drt_discret.H"


/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)            mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(const Epetra_Vector& global,
                                 std::vector<double>& local,
                                 const std::vector<int>& lm)
{
  const size_t ldim = lm.size();
  local.resize(ldim);
  for (size_t i=0; i<ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 |  locally extract a subset of values  (public)             henke 12/09|
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyValues(const Epetra_Vector&      global,
                                 Epetra_SerialDenseVector& local,
                                 const std::vector<int>&        lm)
{
  const size_t ldim = lm.size();
  local.Size(ldim);
  for (size_t i=0; i<ldim; ++i)
  {
    const int lid = global.Map().LID(lm[i]);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),lm[i]);
    local[i] = global[lid];
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based (multi) vector           |
 |                                                          henke 06/09 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(
    const DRT::Element* ele,
    std::vector<double>& local,
    const Epetra_MultiVector& global)
{
  const int numnode = ele->NumNode();
  const int numcol = global.NumVectors();
  local.resize(numnode*numcol);

  // loop over element nodes
  for (int i=0; i<numnode; ++i)
  {
    const int nodegid = (ele->Nodes()[i])->Id();
    const int lid = global.Map().LID(nodegid);
    if (lid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",global.Comm().MyPID(),nodegid);

    // loop over multi vector columns (numcol=1 for Epetra_Vector)
    for (int col=0; col<numcol; col++)
    {
      double* globalcolumn = (global)[col];
      local[col+(numcol*i)] = globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract local values from global node-based multi vector   gjb 08/08 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::ExtractMyNodeBasedValues(
    const DRT::Element* ele,
    Epetra_SerialDenseVector& local,
    const Teuchos::RCP<Epetra_MultiVector>& global,
    const int nsd
    )
{
  if (global==Teuchos::null) dserror("received a TEUCHOS::null pointer");
  if (nsd > global->NumVectors())
    dserror("Requested %d of %d available columns", nsd,global->NumVectors());
  const int iel = ele->NumNode(); // number of nodes
  if (local.Length()!=(iel*nsd)) dserror("vector size mismatch.");

  for (int i=0; i<nsd; i++)
  {
    // access actual component column of multi-vector
    double* globalcolumn = (*global)[i];
    // loop over the element nodes
    for (int j=0;j<iel;j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      const int lid = global->Map().LID(nodegid);
      local(i+(nsd*j))=globalcolumn[lid];
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | extract location vector based on numdof of dis      winklmaier 12/12 |
 *----------------------------------------------------------------------*/
void DRT::UTILS::DisBasedLocationVector(
    const DRT::Discretization & dis,
    const DRT::Element& ele,
    std::vector<int>& lm,
    const int num)
{
  lm.clear();
  std::vector<int> giddofs;
  const int numnodes = ele.NumNode();
  for (int i=0;i<numnodes;i++)
  {
    giddofs.clear();
    giddofs = dis.Dof(ele.Nodes()[i]);
    for (int j=0;j<num;j++)
      lm.push_back(giddofs[j]);
  }
}


/*---------------------------------------------------------------------------------------*
 * Equate the values at DOFs of global nodeids, given by original and copy
 * After equating dof_values(copy) = dof_values(original)                         sudhakar 12/13
 * At the moment, this is used when duplicating nodes in crack simulation
 *---------------------------------------------------------------------------------------*/
void DRT::UTILS::EquateValuesAtTheseNodes( Epetra_Vector& vec,
                                           Teuchos::RCP<DRT::Discretization> dis,
                                           int original,
                                           int copy )
{
  if( dis->HaveGlobalNode( copy ) )
  {
    DRT::Node * newnode = dis->gNode( copy );
    if( newnode->Owner() == dis->Comm().MyPID() )
    {
      DRT::Node * oldnode = dis->gNode( original );

      if( not (oldnode->Owner() == dis->Comm().MyPID() ) )
        dserror( "Both nodes should be owned by the same processor" );

      const std::vector<int> dof_new = dis->Dof( newnode );
      const std::vector<int> dof_old = dis->Dof( oldnode );

      if( not (dof_new.size() == dof_old.size() ) )
        dserror( "Both nodes should have same number of DOFs" );

      const int lid_new = vec.Map().LID(dof_new[0]);
      const int lid_old = vec.Map().LID(dof_old[0]);
      if ( lid_new<0 ) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",vec.Comm().MyPID(),dof_new[0]);
      if ( lid_old<0 ) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",vec.Comm().MyPID(),dof_old[0]);

      for( unsigned i=0; i<dof_new.size(); i++ )
      {
        vec[ lid_new + i ] = vec[ lid_old + i ];
      }
    }
  }
}

DRT::UTILS::Random::Random():
      rand_engine_(0),                        //< set random seed
      uni_dist_(-1.0, 1.0),                         //< set range of uniform distributed rnd no
      uni_rand_no_(rand_engine_, uni_dist_),  //< create the actual rnd no generator
      norm_dist_(0.0, 1.0),                   //< set mean and variance for normal distribution
      norm_rand_no_(rand_engine_, norm_dist_) //< create the actual rnd no generator
{}

DRT::UTILS::Random::~Random()
{}

/// get a random number
double DRT::UTILS::Random::Uni()
{
  return uni_rand_no_();
}

/// get a random number
double DRT::UTILS::Random::Normal()
{
  return norm_rand_no_();
}

/// set the random seed
void DRT::UTILS::Random::SetRandSeed(const unsigned int seed)
{
  rand_engine_.seed(seed);
}

/// set the range for the uniform rng
void DRT::UTILS::Random::SetRandRange(const double lower, const double upper)
{
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
  boost::random::uniform_real_distribution<double>::param_type parm(lower, upper);
  uni_rand_no_.distribution().param(parm);
  return;
#else
  dserror("Your outdated boost version does not support changing the range afterwards!");
#endif
}

/// set the mean and variance for the normal rng
void DRT::UTILS::Random::SetMeanVariance(const double mean, const double var)
{
#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
  boost::random::normal_distribution<double>::param_type parm(mean, var);
  norm_rand_no_.distribution().param(parm);
  return;
#else
  dserror("Your outdated boost version does not support changing mean or sigma afterwards!");
#endif
}

