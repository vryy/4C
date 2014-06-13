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
  double starttime = actiontime_->at(dbctimeindex_);
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
        int freedir = statmechparams_.get<int>("DBCDISPDIR",0)-1;
        int curvenum = statmechparams_.get<int>("NBCCURVENUMBER", 0)-1;
        double nbcval = statmechparams_.get<double>("NBCCREEPFORCE",0.0);

        // some checks
        if (!discret_->Filled()) dserror("FillComplete() was not called");
        if (!discret_->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
        if(curvenum<0) dserror("Invalid NBC curve number %i! Check NBCCURVENUMBER", curvenum);
        int numdbclocal = (int)dbcnodesets_[0].size();
        int numdbcglobal = -1;
        discret_->Comm().SumAll(&numdbclocal,&numdbcglobal,1);
        if (!numdbcglobal) dserror("PointNeumann condition does not have nodal cloud");

        // get the current time
        const double time = params.get("total time",-1.0);

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
        val.at(freedir) = nbcval;

        double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
        for (int i=0; i<(int)dbcnodesets_[0].size(); ++i)
        {
          if (!discret_->NodeRowMap()->MyGID(dbcnodesets_[0][i])) continue;
          DRT::Node* actnode = discret_->gNode(dbcnodesets_[0][i]);
          if (!actnode) dserror("Cannot find global node %d",dbcnodesets_[0][i]);
          std::vector<int> dofs = discret_->Dof(0,actnode);
          const unsigned numdf = dofs.size();
          for (unsigned j=0; j<numdf; ++j)
          {
            if (onoff[j]==0) continue;
            const int gid = dofs[j];
            double value  = val[j];
            value *= curvefac;
            const int lid = systemvector->Map().LID(gid);
            if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
            (*systemvector)[lid] += value;
          }
        }
      }
    }
    break;
    default: break;
  }
  return;
}
