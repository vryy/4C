/*!----------------------------------------------------------------------
\file statmech_mechanaly.cpp
\brief mechanical analysis and output for StatMech network structures

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/

#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../linalg/linalg_utils.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3r/beam3r.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_truss3/truss3.H"
#include "../drt_truss3cl/truss3cl.H"
#include "../drt_beam3cl/beam3cl.H"
#include "../drt_torsion3/torsion3.H"
#include "../drt_lib/drt_globalproblem.H"

/*------------------------------------------------------------------------------*
| Viscoelasticity output                                  (public) mueller 11/12|
*-------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::ViscoelasticityOutput(const double&        time,
                                                      const Epetra_Vector& dis,
                                                      const Epetra_Vector& fint,
                                                      std::ostringstream&  filename)
{
#ifdef DEBUG
  if (forcesensor_ == Teuchos::null)
    dserror("forcesensor_ is NULL pointer; possible reason: dynamic crosslinkers not activated and forcesensor applicable in this case only");
#endif  // #ifdef DEBUG
  double f = 0;// summed force
  double d = 0;//Displacement

  // forcesensor_ is unique on each Proc (row map!) (see UpdateForceSensors() !)
  for(int i=0; i<forcesensor_->MyLength(); i++)
    if((*forcesensor_)[i]>0.9)
    {
      // translate i to DofRowMap LID
      int rowid = discret_->DofRowMap()->LID(discret_->DofColMap()->GID(i));
      if(rowid>-1)
      {
        f += fint[rowid];
        d = dis[rowid];
      }
    }

  //f is the sum of all forces at the top on this processor; compute the sum fglob on all processors all together
  double fglob = 0;
  discret_->Comm().SumAll(&f,&fglob,1);

  if(!discret_->Comm().MyPID())
  {
    //pointer to file into which each processor writes the output related with the dof of which it is the row map owner
    FILE* fp = NULL;

    //content to be written into the output file
    std::stringstream filecontent;

    fp = fopen(filename.str().c_str(), "a");

    /*the output to be written consists of internal forces at exactly those degrees of freedom
     * marked in *forcesensor_ by a one entry*/

    filecontent << std::scientific << std::setprecision(15) << time;//changed

    //Putting time, displacement, meanforce  in Filestream
    filecontent << "   "<< d << "   " << fglob << "   " << discret_->NumGlobalElements() << std::endl; //changed
    //writing filecontent into output file and closing it
    fputs(filecontent.str().c_str(), fp);
    fclose(fp);
  }
}

/*------------------------------------------------------------------------------*
 | output element internal forces                                  mueller 5/12 |
 *------------------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputElementSpatialInternalForces(const std::ostringstream& filename)
{
  std::stringstream elementfint;

  for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
  {
    if(pid==discret_->Comm().MyPID())
    {
      FILE* fp = NULL;
      if(pid==0)
        fp = fopen(filename.str().c_str(), "w");
      else
        fp = fopen(filename.str().c_str(), "a");

      for(int i=0; i<discret_->ElementRowMap()->NumMyElements(); i++)
      {
        DRT::Element* element = discret_->lRowElement(i);
        // element internal force vector
        Epetra_SerialDenseVector force(6*element->NumNode());
        // normal strain
        double eps = 0.0;

        const DRT::ElementType &eot = element->ElementType();
        if(eot == DRT::ELEMENTS::Beam3Type::Instance())
        {
          force = (dynamic_cast<DRT::ELEMENTS::Beam3*>(element))->InternalForceVector();
          eps = (dynamic_cast<DRT::ELEMENTS::Beam3*>(element))->EpsilonSgn();
        }
        else if(eot == DRT::ELEMENTS::Beam3rType::Instance())
        {
          force = (dynamic_cast<DRT::ELEMENTS::Beam3r*>(element))->InternalForceVector();
          eps = (dynamic_cast<DRT::ELEMENTS::Beam3r*>(element))->EpsilonSgn();
        }
        else if(eot == DRT::ELEMENTS::BeamCLType::Instance())
        {
          force = (dynamic_cast<DRT::ELEMENTS::BeamCL*>(element))->InternalForceVector();
          eps = (dynamic_cast<DRT::ELEMENTS::BeamCL*>(element))->EpsilonSgn();
        }
        else
          dserror("No implementation for other Beam elements yet!");

        LINALG::Matrix<3,1> nodalforce;
        for(int j=0; j<(int)nodalforce.M(); j++)
          nodalforce(j) = force(j);
        double fabsolute = nodalforce.Norm2();
        if(eps<0.0)
          fabsolute *= -1.0;
        elementfint << fabsolute <<std::endl;
      }

      fputs(elementfint.str().c_str(), fp);
      fclose(fp);
    }
    discret_->Comm().Barrier();
  }
  return;
}


/*----------------------------------------------------------------------*
 | output of internal forces in material coords  (private) mueller 10/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputElementMaterialInternalForces(const Epetra_Vector&      disrow,
                                                                    const std::ostringstream& filename)
{
  Epetra_Vector discol(*(discret_->DofColMap()),true);
  LINALG::Export(disrow,discol);
  std::map<int, LINALG::Matrix<3, 1> > currentpositions;
  std::map<int, LINALG::Matrix<3, 1> > currentrotations;
  GetNodePositionsFromDisVec(discol,currentpositions,currentrotations,true);

  Teuchos::RCP<Epetra_Vector> element2crosslink(element2crosslink_);
  if(element2crosslink_==Teuchos::null)
  {
    element2crosslink = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementColMap()), true));
    element2crosslink->PutScalar(-1.0);
    ElementToCrosslinkMapping(element2crosslink);
  }

  std::stringstream filelefint;
  std::stringstream crosselefint;

  for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
  {
    if(pid==discret_->Comm().MyPID())
    {
      FILE* fp = NULL;
      if(pid==0)
      {
        fp = fopen(filename.str().c_str(), "w");
        filelefint <<std::scientific<<std::setprecision(8)<<basiselements_<<"\t"<<basisnodes_<<std::endl;
      }
      else
        fp = fopen(filename.str().c_str(), "a");

      for(int i=0; i<discret_->NumMyRowElements(); i++)
      {
        // material forces at element center
        DRT::Element* rowele = discret_->lRowElement(i);
        LINALG::Matrix<3,1> force;
        const DRT::ElementType &eot = rowele->ElementType();
        if(eot == DRT::ELEMENTS::Beam3Type::Instance())
          force = (dynamic_cast<DRT::ELEMENTS::Beam3*>(rowele))->MatForceGp();
        else if(eot == DRT::ELEMENTS::Beam3rType::Instance())
          force = (dynamic_cast<DRT::ELEMENTS::Beam3r*>(rowele))->MatForceGp();
        else if(eot == DRT::ELEMENTS::BeamCLType::Instance())
          force = (dynamic_cast<DRT::ELEMENTS::BeamCL*>(rowele))->MatForceGp();
        else
          dserror("No implementation for other Beam elements yet!");

        if(rowele->Id()<basisnodes_) // filament element
        {
          filelefint <<std::scientific<<std::setprecision(8)<<rowele->Id()<<"\t"<<(*filamentnumber_)[discret_->NodeColMap()->LID(rowele->NodeIds()[0])]
                                                                          <<"\t"<<force(0)<<"\t"<<force(1)<<"\t"<<force(2)<<"\t";
          // get nodal positions
          for(int j=0; j<rowele->NumNode(); j++)
          {
            int nodecollid = discret_->NodeColMap()->LID(rowele->NodeIds()[j]);
            std::map< int,LINALG::Matrix<3,1> >::const_iterator it = currentpositions.find(nodecollid);
            filelefint <<(it->second)(0)<<"\t"<<(it->second)(1)<<"\t"<<(it->second)(2)<<"\t";
          }
          if(crosslinkertype_!=Teuchos::null)
            filelefint<<"-1\t-1";
          filelefint<<std::endl;
        }
        else
        {
          int crosslid = (*element2crosslink)[discret_->ElementColMap()->LID(discret_->ElementRowMap()->GID(i))];
          crosselefint <<std::scientific<<std::setprecision(8)<<rowele->Id()<<"\t"<<(*filamentnumber_)[discret_->NodeColMap()->LID(rowele->NodeIds()[0])]
                                                                            <<"\t"<<force(0)<<"\t"<<force(1)<<"\t"<<force(2)<<"\t";
          // get nodal positions
          for(int j=0; j<rowele->NumNode(); j++) // crosslinker element
          {
            int nodecollid = discret_->NodeColMap()->LID(rowele->NodeIds()[j]);
            std::map< int,LINALG::Matrix<3,1> >::const_iterator it = currentpositions.find(nodecollid);
            crosselefint <<(it->second)(0)<<"\t"<<(it->second)(1)<<"\t"<<(it->second)(2)<<"\t";
          }
          if(crosslinkertype_!=Teuchos::null)
            crosselefint<<(int)(*crosslinkertype_)[crosslid]<<"\t"<<(int)(*crosslinkeractlength_)[crosslid];
          crosselefint<<std::endl;
        }
      }
      // write filament information
      fputs(filelefint.str().c_str(), fp);
      fclose(fp);
    }
    discret_->Comm().Barrier();
  }

  // write crosslinker information at the end
  for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
  {
    if(pid==discret_->Comm().MyPID())
    {
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "a");
      fputs(crosselefint.str().c_str(), fp);
      fclose(fp);
    }
    discret_->Comm().Barrier();
  }

  return;
}

/*----------------------------------------------------------------------*
 | output of internal forces in material coords  (private) mueller 11/14|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::OutputNeumannPointForce(const double&             time,
                                                        const std::ostringstream& filename)
{
  std::stringstream fileneumanptforce;

  int curvenum = statmechparams_.get<int>("NBCCURVENUMBER", 0)-1;
  double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
  double amp = statmechparams_.get<double>("NBCFORCEAMP",0.0);

  for(int pid=0; pid<discret_->Comm().NumProc(); pid++)
  {
    if(pid==discret_->Comm().MyPID())
    {
      FILE* fp = NULL;
      fp = fopen(filename.str().c_str(), "a");

      int nodesetindex = timeintervalstep_-bctimeindex_;
      // necessary adjustment beyond the first timeintervalstep_ since timeintervalstep_ does not mark the current position in actiontime_
      // but the next position to reach (hence --).
      if((int)nbcnodesets_.size()==1)
        nodesetindex = 0;
      else if(timeintervalstep_>1)// && timeintervalstep_<(int)actiontime_->size())
        nodesetindex--;

      for(int i=0; i<(int)nbcnodesets_[nodesetindex].size(); i++)
      {
        if (discret_->NodeRowMap()->MyGID(nbcnodesets_[nodesetindex][0]))
        {
          fileneumanptforce <<std::scientific<<std::setprecision(8)<<time<<"\t"<<nbcnodesets_[nodesetindex][0]<<"\t"<<std::scientific<<std::setprecision(8)<<amp*curvefac<<std::endl;
        }
      }

      fputs(fileneumanptforce.str().c_str(), fp);
      fclose(fp);
    }
    discret_->Comm().Barrier();
  }

  return;
}
