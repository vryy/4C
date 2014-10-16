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
        Teuchos::RCP<Epetra_MultiVector> nodepositions = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()), 3, true));
        Epetra_Vector discol(*discret_->DofColMap(),true);
        LINALG::Export(*disn, discol);
        GetNodalBindingSpotPositionsFromDisVec(discol, nodepositions, Teuchos::null);

        GetNBCNodes(params, nodepositions);

        int oscdir = statmechparams_.get<int>("DBCDISPDIR",0)-1;
        int curvenum = statmechparams_.get<int>("NBCCURVENUMBER", 0)-1;
        double nbcval = statmechparams_.get<double>("NBCFORCEAMP",0.0);

        // some checks
        if(curvenum<0) dserror("Invalid NBC curve number %i! Check NBCCURVENUMBER", curvenum);
        int numnbclocal = (int)nbcnodesets_[0].size();
        int numnbcglobal = -1;
        discret_->Comm().SumAll(&numnbclocal,&numnbcglobal,1);
        if (!numnbcglobal) dserror("PointNeumann condition does not have nodal cloud");

        // construct Point Neumann BC definition vector
        std::vector<int> curve(6,0);
        std::vector<int> onoff(6,0);
        std::vector<double> val(6,0.0);

        // all three spatial direction subject to time-dependent force pattern
//        for(int i=0; i<3; i++)
//        {
//          onoff.at(i) = 1;
//          curve.at(i) = 1;
//          val.at(i) = nbcval;
//        }
        // only in one direction
        onoff.at(oscdir) = 1;
        curve.at(oscdir) = 1;
        val.at(oscdir) = nbcval;

        int nodesetindex = timeintervalstep_-bctimeindex_;
        // necessary adjustment beyond the first timeintervalstep_ since timeintervalstep_ does not mark the current position in actiontime_
        // but the next position to reach (hence --).
        if((int)nbcnodesets_.size()==1)
          nodesetindex = 0;
        else
        {
          if(timeintervalstep_>1)// && timeintervalstep_<(int)actiontime_->size())
            nodesetindex--;
        }

        ApplyNeumannValueStatMech(params, onoff, curve, val, systemvector, nodesetindex);
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
void STATMECH::StatMechManager::GetNBCNodes(Teuchos::ParameterList&          params,
                                            Teuchos::RCP<Epetra_MultiVector> nodeposcol)
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
        int N = statmechparams_.get<int>("NUMNBCNODES",0);
        if(N<=0)  dserror("Provide the number of Neumann nodes by means of the input parameter NUMNBCNODES (>0)!");

        // find node closest to COG
        if(N==1)
        {
          // calculate center of gravity
          LINALG::Matrix<3,1> cog(true);
          for(int i=0; i<nodeposcol->MyLength(); i++)
            for(int j=0; j<nodeposcol->NumVectors(); j++)
              cog(j) += (*nodeposcol)[j][i];
          cog.Scale(1.0/(double)(nodeposcol->MyLength()));

          // all processors should arrive at the same node (also a check for consistency!)
          double dmin = 1e9;
          std::vector<int> nodecolid(1,-1);
          if(nodeposcol==Teuchos::null) dserror("No nodal positions supplied!");

          for(int i=0; i<nodeposcol->MyLength(); i++)
          {
            double disti = sqrt(((*nodeposcol)[0][i]-cog(0))*((*nodeposcol)[0][i]-cog(0))+
                                ((*nodeposcol)[1][i]-cog(1))*((*nodeposcol)[1][i]-cog(1))+
                                ((*nodeposcol)[2][i]-cog(2))*((*nodeposcol)[2][i]-cog(2)));
            if(disti<dmin)
            {
              nodecolid[0] = nodeposcol->Map().GID(i);
              dmin = disti;
            }
          }
          // check
          int procid = nodecolid[0];
          int summed = -1;
          discret_->Comm().SumAll(&procid,&summed,1);
          if(summed!=discret_->Comm().NumProc()*procid) dserror("Wrong check sum! Either communication went wrong or your column map node vector is incorrect!");

          nbcnodesets_.push_back(nodecolid);

          if(!discret_->Comm().MyPID())
          {
            std::cout<<"=== Point force location ==="<<std::endl;
            std::cout<<"X: "<<(*nodeposcol)[0][nodecolid[0]]<<std::endl;
            std::cout<<"Y: "<<(*nodeposcol)[1][nodecolid[0]]<<std::endl;
            std::cout<<"Z: "<<(*nodeposcol)[2][nodecolid[0]]<<std::endl;
            std::cout<<"============================"<<std::endl;
          }
          discret_->Comm().Barrier();
        }
        else // choose N nodes randomly
        {
          std::vector<int> nbcnodes;
          nbcnodes.clear();
          if(!discret_->Comm().MyPID())
            for(int i=0; i<N; i++)
              nbcnodes.push_back((int)(floor((*uniformgen_)()*(double)(discret_->NumGlobalNodes()))));
          discret_->Comm().Barrier();
          // communicate
          // gid vector
          std::vector<int> gids(N,0);
          for(int i=0; i<N; i++)
            gids.at(i) = i;
          Epetra_Map fullmap(-1, N, &gids[0], 0, discret_->Comm());
          Epetra_Map partmap(N, 0, discret_->Comm());
          Teuchos::RCP<Epetra_Vector> target(new Epetra_Vector(fullmap,true));
          Teuchos::RCP<Epetra_Vector> part(new Epetra_Vector(partmap, true));
          if(!discret_->Comm().MyPID())
            for(int i=0; i<target->MyLength(); i++)
              (*target)[i] = (double)nbcnodes.at(i);
          CommunicateVector(part, target, true, true);

          // give all processors all nodes (the ordering is important information)
          for(int i=0; i<N; i++)
            nbcnodesets_.push_back(std::vector<int>(1,(*target)[i]));
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
  if(curvenum<0)  dserror("Provide NBCCURVENUMBER>=1 in your input file!");

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
