/*!----------------------------------------------------------------------
\file statmech.cpp
\brief time integration for structural problems with statistical mechanics

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15234
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "statmech_time.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 08/08|
 *----------------------------------------------------------------------*/
StatMechTime::StatMechTime(ParameterList& params,
                          DRT::Discretization& dis,
                          LINALG::Solver& solver,
                          IO::DiscretizationWriter& output) :
statmechmanager_(params),
StruGenAlpha(params,dis,solver,output)
{
  return;
} // StatMechTime::StatMechTime


/*----------------------------------------------------------------------*
 |  integrate in time          (static/public)               cyron 08/08|
 *----------------------------------------------------------------------*/
void StatMechTime::Integrate()
{
  int    step    = params_.get<int>   ("step" ,0);
  int    nstep   = params_.get<int>   ("nstep",5);
  double maxtime = params_.get<double>("max time",0.0);
  // can have values "full newton" , "modified newton" , "nonlinear cg"
  string equil = params_.get<string>("equilibrium iteration","full newton");

  // can have values takes values "constant" consistent"
  string pred  = params_.get<string>("predictor","constant");
  int predictor=-1;
  if      (pred=="constant")   predictor = 1;
  else if (pred=="consistent") predictor = 2;
  else dserror("Unknown type of predictor");

  if (equil=="full newton")
  {

  double dt = params_.get<double>("delta time" ,0.01);
  
  int num_dof = (*fext_).GlobalLength();


  Epetra_SerialDenseVector  v0;
  v0.Size(num_dof);
  Epetra_SerialDenseVector d0;
  d0.Size(num_dof);
  Epetra_SerialDenseVector d1_ap;
  d1_ap.Size(num_dof);
  Epetra_SerialDenseVector v1_ap;
  v1_ap.Size(num_dof);
  Epetra_SerialDenseVector fint0;
  fint0.Size(num_dof);
  Epetra_SerialDenseVector relerr_d;
  relerr_d.Size(num_dof);
  Epetra_SerialDenseVector relerr_v;
  relerr_v.Size(num_dof);
  Epetra_SerialDenseVector Delta_d;
  Delta_d.Size(num_dof);
  Epetra_SerialDenseVector Delta_v;
  Delta_v.Size(num_dof);
  int kd = 0;
  double gamma = 0.125663706 / ((num_dof/6) - 1);

    for (int i=step; i<nstep; ++i)
    {
      double time = params_.get<double>("total time",0.0);

      predictor = 2;
      /*
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      */
      //if(i<10000)
        ConsistentPredictor();
      //else
       // BrownianPredictor3D();

      FullNewton();
      UpdateandOutput();

      //Freiheitsgrade längs zur Filamentachse: Da nur geringe axiale Dehnung zu erwarten ist, kann angenommen werden,
      //dass alle Freiheitsgrade in Längsrichtung dieselbe Bewegung Delta_x ausführen, die approximiert werden kann durch:
      // Gamma * Delta_x / dt = fext_axial, wobei fext_axial die Summe der externen Kräfte in Axialrichtung längs des
      //gesamten Filaments ist und Gamma die Gesamtreibung eines Filaments der Länge 10 gegenüber axialer Verschiebung ist
      double fext_axial = 0;
      for(int id = 0; id < num_dof; id = id+6)
      {
        fext_axial += (*fext_)[id];
      }
      //Freiheitsgrade entlang der Filamentachse:
      v1_ap(0) = fext_axial / 0.125663706;
      d1_ap(0) = 0.5*dt*(v0(0) + v1_ap(0)) + d0(0);
      double lrefe = 10.0 / (num_dof/6 - 1);

      for(int jd = 0; jd < (num_dof/6); jd++)
      {
        //Freiheitsgrade entlang der Filamentachse aus Undehnbarkeitsbedingung:
        //v1_ap(jd*6) = fext_axial / 0.125663706;
        //d1_ap(jd*6) = 0.5*dt*(v0(jd*6) + v1_ap(jd*6)) + d0(jd*6);


        //Freiheitsgrade quer zur Filamentachse
        for(int id = 1; id < 3; id++)
        {
          kd = 6*jd + id;
          v1_ap(kd) = ( (*fext_)[kd] - fint0(kd) ) / gamma;
          //an den Randknoten nur jeweils halbes gamma:
          if (jd == 0 || jd == (num_dof/6 -1) )
            v1_ap(kd) = 2*v1_ap(kd);
          d1_ap(kd) = 0.5*dt*(v0(kd) + v1_ap(kd)) + d0(kd);
        }

        if(jd>0)
        {
          double dy = d1_ap(jd*6+1) - d1_ap((jd-1)*6+1);
          double dz = d1_ap(jd*6+2) - d1_ap((jd-1)*6+2);
          d1_ap(jd*6) = pow(lrefe*lrefe - dy*dy - dz*dz  ,0.5) + d1_ap((jd-1)*6) - lrefe;
          v1_ap(jd*6) = 2*(d1_ap(jd*6) - d0(jd*6))/dt - v0(jd*6);
        }

        for(int id = 0; id < 6; id++)
        {
          kd = 6*jd + id;
          //Berechnung des relativen Fehlers im Prädiktorschritt:
          Delta_d(kd) = (*dis_)[kd] - d0(kd);
          Delta_v(kd) = (*velm_)[kd] - v0(kd);
          relerr_d(kd) = ( (d1_ap(kd) - d0(kd) ) - Delta_d(kd) ) / Delta_d(kd);
          relerr_v(kd) = ( (v1_ap(kd) - v0(kd) ) - Delta_v(kd) ) / Delta_v(kd);

          //Zwischenspeichern der Endgrößen im abgeschlossenen Zeitschritt
          d0(kd) = (*dis_)[kd];
          v0(kd) = (*velm_)[kd];
          fint0(kd) = (*fint_)[kd];
        }
      }


      //std::cout<<"\nfext nach update and output"<<*fext_;
      //std::cout<<"\nvelm nach update and output"<<*velm_;
      //std::cout<<"\nfint nach update and output"<<*fint_;
      //std::cout<<"\ndis nach update and output"<< d0;

      //std::cout<<"\n*dis_ nach update and output"<<d0;

      //std::cout<<"\n\nrelerr_d \n" << relerr_d;
      //std::cout<<"\n\nrelerr_v \n" << relerr_v;
      //std::cout<<"\n\nDelta_d \n" << Delta_d;
      //std::cout<<"\n\nDelta_v \n" << Delta_v;


      statmechmanager_.StatMechOutput(time,num_dof,i,dt,*dis_);

      if (time>=maxtime) break;
    }
#if 0
    for (int i=0; i<discret_.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = discret_.lRowNode(i);
      printf("NODE %d COORD ",actnode->Id()+1);
      for (int j=0; j<discret_.NumDof(actnode); ++j)
      {
        const int gdof = discret_.Dof(actnode,j);
        const int lid  = dis_->Map().LID(gdof);
        printf("%20.15f ",actnode->X()[j]+(*dis_)[lid]);
      }
      printf("\n");
    }
#endif
  }
  else if (equil=="line search newton")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      LineSearchNewton();
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
#if 0
    for (int i=0; i<discret_.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = discret_.lRowNode(i);
      printf("NODE %d COORD ",actnode->Id()+1);
      for (int j=0; j<discret_.NumDof(actnode); ++j)
      {
        const int gdof = discret_.Dof(actnode,j);
        const int lid  = dis_->Map().LID(gdof);
        printf("%20.15f ",actnode->X()[j]+(*dis_)[lid]);
      }
      printf("\n");
    }
#endif
  }
  else if (equil=="modified newton")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      ModifiedNewton();
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
  }
  else if (equil=="nonlinear cg")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      NonlinearCG();
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
  }
  else if (equil=="ptc")
  {
    for (int i=step; i<nstep; ++i)
    {
      if      (predictor==1) ConstantPredictor();
      else if (predictor==2) ConsistentPredictor();
      PTC();
      UpdateandOutput();
      double time = params_.get<double>("total time",0.0);
      if (time>=maxtime) break;
    }
  }
  else dserror("Unknown type of equilibrium iteration");

  return;
} // void StatMechTime::Integrate()


#endif  // #ifdef CCADISCRET
