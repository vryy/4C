/*!----------------------------------------------------------------------
\file statmech_neumann.cpp
\brief special Neumann boundary conditions for StatMech problems

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/

#include "statmech_manager.H"
#include "statmech_search.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../drt_inpar/inpar_structure.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | Determine if application of NBCs starts at a given time  mueller 5/12|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::NBCStart(Teuchos::ParameterList& params)
{
  double eps = 2.0e-11;
  // get the current time
  double time = params.get<double>("total time", 0.0);
  double starttime = actiontime_->at(bctimeindex_);
  //double dt = params.get<double>("delta time", 0.01);
  if (time<0.0) dserror("t = %f ! Something is utterly wrong here. The absolute time should be positive!", time);

  if(time < starttime && fabs(starttime-time)>eps)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------*
 | Create Dirichlet DOF maps                    (private)  mueller 09/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::EvaluateNeumannStatMech(Teuchos::ParameterList&              params,
                                                        Teuchos::RCP<Epetra_Vector>          disn,
                                                        Teuchos::RCP<Epetra_Vector>          systemvector,
                                                        Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  // some checks
  if (!discret_->Filled()) dserror("FillComplete() was not called");
  if (!discret_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get load vector
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  bool loadlin = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "LOADLIN");
  INPAR::STATMECH::NBCType nbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::NBCType>(statmechparams_,"NBCTYPE");
  switch(nbctype)
  {
    case INPAR::STATMECH::nbctype_std:
    {
      if (!loadlin)
        discret_->EvaluateNeumann(params, systemvector);
      else
      {
        discret_->SetState(0,"displacement new", disn);
        discret_->EvaluateNeumann(params, systemvector, systemmatrix);
      }
    }
    break;
    case INPAR::STATMECH::nbctype_constcreep:
    {
      if(NBCStart(params))
      {
        GetNBCNodes(params);

        int freedir = statmechparams_.get<int>("DBCDISPDIR",0)-1;
        int curvenum = statmechparams_.get<int>("NBCCURVENUMBER", 0)-1;
        double amp = statmechparams_.get<double>("NBCFORCEAMP",0.0);

        // some checks
        int numnbclocal = (int)nbcnodesets_[0].size();
        int numnbcglobal = -1;
        discret_->Comm().SumAll(&numnbclocal,&numnbcglobal,1);
        if(curvenum<0)
          dserror("Invalid NBC curve number %i! Check NBCCURVENUMBER", curvenum);
        if (!numnbcglobal)
          dserror("PointNeumann condition does not have nodal cloud");
        if(amp<=0.0)
          dserror("NBCFORCEAMP = %d! Choose a meaninful value!",amp);

        // construct Point Neumann BC definition vector
        std::vector<int> curve(6,0);
        std::vector<int> onoff(6,0);
        std::vector<double> val(6,0.0);

        curve.at(freedir) = 1;
        for(int i=0; i<3; i++)
          if(i!=freedir)
            onoff.at(i) = 0;
          else
            onoff.at(i) = 1;
        val.at(freedir) = amp;

        ApplyNeumannValueStatMech(params,onoff,curve,val,systemvector);
      }
    }
    break;
    case INPAR::STATMECH::nbctype_randompointforce:
    {
      if(NBCStart(params))
      {
        GetNBCNodes(params);

        int curvenum = statmechparams_.get<int>("NBCCURVENUMBER", 0)-1;
        double nbcval = statmechparams_.get<double>("NBCFORCEAMP",0.0);

        // some checks
        if(curvenum<0) dserror("Invalid NBC curve number %i! Check NBCCURVENUMBER", curvenum);
        int numdbclocal = (int)dbcnodesets_[0].size();
        int numdbcglobal = -1;
        discret_->Comm().SumAll(&numdbclocal,&numdbcglobal,1);
        if (!numdbcglobal) dserror("PointNeumann condition does not have nodal cloud");

        // construct Point Neumann BC definition vector
        std::vector<int> curve(6,0);
        std::vector<int> onoff(6,0);
        std::vector<double> val(6,0.0);

        // all three spatial direction subject to time-dependent force pattern
        for(int i=0; i<3; i++)
        {
          onoff.at(i) = 1;
          curve.at(i) = 1;
          val.at(i) = nbcval;
        }

        ApplyNeumannValueStatMech(params, onoff, curve, val, systemvector,bctimeindex_-timeintervalstep_);
      }
    }
    break;
    default: break;
  }
  return;
}

/*----------------------------------------------------------------------*
 | Retrieve Neumann nodes                       (private)  mueller 09/14|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GetNBCNodes(Teuchos::ParameterList& params)
{
  INPAR::STATMECH::NBCType nbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::NBCType>(statmechparams_,"NBCTYPE");
  switch(nbctype)
  {
    case INPAR::STATMECH::nbctype_std:
      break;
    case INPAR::STATMECH::nbctype_constcreep:
    {
      // since Point Neumann forces are applied to nodes also affected by Dirichlet BCs, nbcnodeset_=dbcnodeset_
      if(nbcnodesets_.empty() && (int)dbcnodesets_.size()>0)
        nbcnodesets_ = dbcnodesets_;
    }
    break;
    case INPAR::STATMECH::nbctype_randompointforce:
    {
      // only once at the beginning. Each node set will be interpreted according to the current time interval defined by ACTIONTIME
      if(nbcnodesets_.empty())
      {
        // for now, hard-coded number of nodes that will be subjected to a Neumann loads SEQUENTIALLY
        int N = 10;
        for(int i=0; i<N; i++)
        {
          std::vector<int> id(1,0);
          id.at(0) = (int)(floor((*uniformgen_)()*(double)(discret_->NumGlobalNodes())));
          nbcnodesets_.push_back(id);
        }
      }
    }
    break;
    default: dserror("Unknown NBCTYPE %i", nbctype);
  }
  return;
}

/*----------------------------------------------------------------------*
 | Apply Neumann                     (private)  mueller 09/14|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ApplyNeumannValueStatMech(Teuchos::ParameterList&      params,
                                                          std::vector<int>&            onoff,
                                                          std::vector<int>&            curve,
                                                          std::vector<double>&         val,
                                                          Teuchos::RCP<Epetra_Vector>  systemvector,
                                                          int                          nodesetindex)
{
  if(nodesetindex<0)
    return;

  int curvenum = statmechparams_.get<int>("NBCCURVENUMBER", 0)-1;
  // get the current time
  const double time = params.get("total time",-1.0);

  double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
  for (int i=0; i<(int)nbcnodesets_[nodesetindex].size(); ++i)
  {
    if (!discret_->NodeRowMap()->MyGID(nbcnodesets_[nodesetindex][i])) continue;
    DRT::Node* actnode = discret_->gNode(nbcnodesets_[nodesetindex][i]);
    if (!actnode) dserror("Cannot find global node %d",nbcnodesets_[nodesetindex][i]);
    std::vector<int> dofs = discret_->Dof(0,actnode);
    const unsigned numdf = dofs.size();
    for (unsigned j=0; j<numdf; ++j)
    {
      if (onoff[j]==0)
        continue;
      const int gid = dofs[j];
      const int lid = systemvector->Map().LID(gid);
      if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
      (*systemvector)[lid] += val[j]*curvefac;
    }
  }
  return;
}
