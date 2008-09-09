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

#include <iostream>


/*----------------------------------------------------------------------*
 |  ctor (public)                                             cyron 08/08|
 *----------------------------------------------------------------------*/
StatMechTime::StatMechTime(ParameterList& params,
                          DRT::Discretization& dis,
                          LINALG::Solver& solver,
                          IO::DiscretizationWriter& output) :
StruGenAlpha(params,dis,solver,output)
{
  return;
} // StatMechTime::StatMechTime




/*----------------------------------------------------------------------*
 |  do Newton iteration (public)                             cyron 08/08|
 *----------------------------------------------------------------------*/
void StatMechTime::FullNewton()
{
  // -------------------------------------------------------------------
  // get some parameters from parameter list
  // -------------------------------------------------------------------
  double time      = params_.get<double>("total time"             ,0.0);
  double dt        = params_.get<double>("delta time"             ,0.01);
  double timen     = time + dt;
  int    maxiter   = params_.get<int>   ("max iterations"         ,10);
  bool   damping   = params_.get<bool>  ("damping"                ,false);
  double beta      = params_.get<double>("beta"                   ,0.292);
#ifdef STRUGENALPHA_BE
  double delta     = params_.get<double>("delta"                  ,beta);
#endif
  double gamma     = params_.get<double>("gamma"                  ,0.581);
  double alpham    = params_.get<double>("alpha m"                ,0.378);
  double alphaf    = params_.get<double>("alpha f"                ,0.459);
  string convcheck = params_.get<string>("convcheck"              ,"AbsRes_Or_AbsDis");
  double toldisp   = params_.get<double>("tolerance displacements",1.0e-07);
  double tolres    = params_.get<double>("tolerance residual"     ,1.0e-07);
  bool printscreen = params_.get<bool>  ("print to screen",true);
  bool printerr    = params_.get<bool>  ("print to err",false);
  FILE* errfile    = params_.get<FILE*> ("err file",NULL);
  bool structrobin = params_.get<bool>  ("structrobin"            ,false);
  if (!errfile) printerr = false;
  bool dynkindstat = (params_.get<string>("DYNAMICTYP") == "Static");
  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_.get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_.get<double>("ADAPTCONV_BETTER",0.01);

  // check whether we have a stiffness matrix, that is not filled yet
  // and mass and damping are present
  if (stiff_->Filled()) dserror("stiffness matrix may not be filled here");
  if (!mass_->Filled()) dserror("mass matrix must be filled here");
  if (damping)
    if (!damp_->Filled()) dserror("damping matrix must be filled here");

  //=================================================== equilibrium loop
  int numiter=0;
  double fresmnorm = 1.0e6;
  double disinorm = 1.0e6;
  fresm_->Norm2(&fresmnorm);
  Epetra_Time timer(discret_.Comm());
  timer.ResetStartTime();
  bool print_unconv = true;

  while (!Converged(convcheck, disinorm, fresmnorm, toldisp, tolres) and numiter<=maxiter)
  {
    //------------------------------------------- effective rhs is fresm
    //---------------------------------------------- build effective lhs
    // (using matrix stiff_ as effective matrix)
    if (dynkindstat); // do nothing, we have the ordinary stiffness matrix ready
    else
    {
#ifdef STRUGENALPHA_BE
      stiff_->Add(*mass_,false,(1.-alpham)/(delta*dt*dt),1.-alphaf);
#else
      stiff_->Add(*mass_,false,(1.-alpham)/(beta*dt*dt),1.-alphaf);
#endif
      if (damping)
      {
#ifdef STRUGENALPHA_BE
        stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(delta*dt),1.0);
#else
        stiff_->Add(*damp_,false,(1.-alphaf)*gamma/(beta*dt),1.0);
#endif
      }
    }
    stiff_->Complete();

    //----------------------- apply dirichlet BCs to system of equations
    disi_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(stiff_,disi_,fresm_,zeros_,dirichtoggle_);

    //--------------------------------------------------- solve for disi
    // Solve K_Teffdyn . IncD = -R  ===>  IncD_{n+1}
    if (isadapttol && numiter)
    {
      double worst = fresmnorm;
      double wanted = tolres;
      solver_.AdaptTolerance(wanted,worst,adaptolbetter);
    }
    solver_.Solve(stiff_->EpetraMatrix(),disi_,fresm_,true,numiter==0);
    solver_.ResetTolerance();

    //---------------------------------- update mid configuration values
    // displacements
    // D_{n+1-alpha_f} := D_{n+1-alpha_f} + (1-alpha_f)*IncD_{n+1}
#ifdef STRUGENALPHA_FINTLIKETR
    disn_->Update(1.0,*disi_,1.0);
    dism_->Update(1.-alphaf,*disn_,alphaf,*dis_,0.0);
#else
    dism_->Update(1.-alphaf,*disi_,1.0);
#endif
    // velocities
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // V_{n+1-alpha_f} := V_{n+1-alpha_f}
    //                  + (1-alpha_f)*gamma/beta/dt*IncD_{n+1}
    velm_->Update((1.-alphaf)*gamma/(beta*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    velm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
#ifdef STRUGENALPHA_BE
    velm_->Update((delta-(1.0-alphaf)*gamma)/delta,*vel_,
                  (1.0-alphaf)*(-gamma-2.*delta*gamma+2.*beta*gamma+2.*delta)*dt/(2.*delta),*acc_,
                  gamma/(delta*dt));
#else
    velm_->Update((beta-(1.0-alphaf)*gamma)/beta,*vel_,
                  (1.0-alphaf)*(2.*beta-gamma)*dt/(2.*beta),*acc_,
                  gamma/(beta*dt));
#endif
#endif
    // accelerations
#ifndef STRUGENALPHA_INCRUPDT
    // iterative
    // A_{n+1-alpha_m} := A_{n+1-alpha_m}
    //                  + (1-alpha_m)/beta/dt^2*IncD_{n+1}
    accm_->Update((1.-alpham)/(beta*dt*dt),*disi_,1.0);
#else
    // incremental (required for constant predictor)
    accm_->Update(1.0,*dism_,-1.0,*dis_,0.0);
#ifdef STRUGENALPHA_BE
    accm_->Update(-(1.-alpham)/(delta*dt),*vel_,
                  (2.*beta-1.+alpham-2.*alpham*beta+2.*alpham*delta)/(2.*delta),*acc_,
                  (1.-alpham)/((1.-alphaf)*delta*dt*dt));
#else
    accm_->Update(-(1.-alpham)/(beta*dt),*vel_,
                  (2.*beta-1.+alpham)/(2.*beta),*acc_,
                  (1.-alpham)/((1.-alphaf)*beta*dt*dt));
#endif
#endif

    // zerofy velocity and acceleration in case of statics
    if (dynkindstat)
    {
      velm_->PutScalar(0.0);
      accm_->PutScalar(0.0);
    }

    //---------------------------- compute internal forces and stiffness
    {
      // zero out stiffness
      stiff_->Zero();
      // create the parameters for the discretization
      ParameterList p;
      // action for elements
      p.set("action","calc_struct_nlnstiff");
      // other parameters that might be needed by the elements
      p.set("total time",timen);
      p.set("delta time",dt);
      p.set("alpha f",alphaf);
      // set vector values needed by elements
      discret_.ClearState();
#ifdef STRUGENALPHA_FINTLIKETR
#else
      // scale IncD_{n+1} by (1-alphaf) to obtain mid residual displacements IncD_{n+1-alphaf}
      disi_->Scale(1.-alphaf);
#endif
      discret_.SetState("residual displacement",disi_);
#ifdef STRUGENALPHA_FINTLIKETR
      discret_.SetState("displacement",disn_);
      discret_.SetState("velocity",veln_);
#else
      discret_.SetState("displacement",dism_);
      discret_.SetState("velocity",velm_);
#endif
      //discret_.SetState("velocity",velm_); // not used at the moment
#ifdef STRUGENALPHA_FINTLIKETR
      fintn_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fintn_,null,null);
#else
      fint_->PutScalar(0.0);  // initialise internal force vector
      discret_.Evaluate(p,stiff_,null,fint_,null,null);
#endif
      discret_.ClearState();

      if (pot_man_!=null)
      {
        p.set("pot_man", pot_man_);
        pot_man_->EvaluatePotential(p,dism_,fint_,stiff_);
      }

      // do NOT finalize the stiffness matrix to add masses to it later
      // If we have a robin condition we need to modify both the rhs and the
      // matrix diagonal corresponding to the dofs at the robin interface.
      if (structrobin)
      {
        double alphas = params_.get<double>("alpha s",-1.);

        // Add structural part of Robin force
        fsisurface_->AddCondVector(alphas/dt,
                                   fsisurface_->ExtractCondVector(dism_),
                                   fint_);

        double scale = alphas*(1.-alphaf)/dt;
        const Epetra_Map& robinmap = *fsisurface_->CondMap();
        int numrdofs = robinmap.NumMyElements();
        int* rdofs = robinmap.MyGlobalElements();
        for (int lid=0; lid<numrdofs; ++lid)
        {
          int gid = rdofs[lid];
          // Note: This assemble might fail if we have a block matrix here.
          stiff_->Assemble(scale,gid,gid);
        }
      }
    }

    //------------------------------------------ compute residual forces
    if (dynkindstat)
    {
      // static residual
      // Res = F_int - F_ext
      fresm_->PutScalar(0.0);
    }
    else
    {
      // dynamic residual
      // Res = M . A_{n+1-alpha_m}
      //     + C . V_{n+1-alpha_f}
      //     + F_int(D_{n+1-alpha_f})
      //     - F_{ext;n+1-alpha_f}
      // add mid-inertial force
      mass_->Multiply(false,*accm_,*finert_);
      fresm_->Update(1.0,*finert_,0.0);
      // add mid-viscous damping force
      if (damping)
      {
        //RefCountPtr<Epetra_Vector> fviscm = LINALG::CreateVector(*dofrowmap,true);
        damp_->Multiply(false,*velm_,*fvisc_);
        fresm_->Update(1.0,*fvisc_,1.0);
      }
    }
    // add static mid-balance
#ifdef STRUGENALPHA_FINTLIKETR
    fresm_->Update(1.0,*fextm_,-1.0);
    fresm_->Update(-(1.0-alphaf),*fintn_,-alphaf,*fint_,1.0);
#else
    fresm_->Update(-1.0,*fint_,1.0,*fextm_,-1.0);
#endif
    // blank residual DOFs that are on Dirichlet BC
    {
      Epetra_Vector fresmcopy(*fresm_);
      fresm_->Multiply(1.0,*invtoggle_,fresmcopy,0.0);
    }

    //---------------------------------------------- build residual norm
    disi_->Norm2(&disinorm);

    fresm_->Norm2(&fresmnorm);

    // a short message
    if (!myrank_ and (printscreen or printerr))
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }

    //--------------------------------- increment equilibrium loop index
    ++numiter;

  }
  //=============================================== end equilibrium loop
  print_unconv = false;

  //-------------------------------- test whether max iterations was hit
  if (numiter>=maxiter)
  {
     dserror("Newton unconverged in %d iterations",numiter);
  }
  else
  {
    if (!myrank_ and printscreen)
    {
      PrintNewton(printscreen,printerr,print_unconv,errfile,timer,numiter,maxiter,
                  fresmnorm,disinorm,convcheck);
    }
  }

  params_.set<int>("num iterations",numiter);

  return;
} // StatMechTime::FullNewton()


/*----------------------------------------------------------------------*
 |  integrate in time          (static/public)               cyron 08/08|
 *----------------------------------------------------------------------*/
void StatMechTime::Integrate()
{
  int    step    = params_.get<int>   ("step" ,0);
  int    nstep   = params_.get<int>   ("nstep",5);
  double maxtime = params_.get<double>("max time",0.0);
  //special parameters for problem type STATISTICHAL MECHANICS
  
  

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
  //defining and initializing a file pointer variable
  FILE* fp = NULL;

  //statistical measurements start at time = maxtime * start_factor
  double start_factor = 0.25; // 0.25 für hinged und free ends
  double endzuend=0;
  double Referenzabstand = 0;
  int schalter = 0;
  int zaehler = 0;
  int i_start;
  double startzeit = 0;
  int num_dof = (*fext_).GlobalLength();

  double DeltaR2 = 0;
  int Dateinummer = 0;
  int testflag = 0;
  double dt = params_.get<double>("delta time" ,0.01);
  //variable "filename" is assigned the name of the file where output should be written
  std::ostringstream filename;

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

      std::cout<<"\n\nSimulation with FEM like forces\n\n";

      zaehler ++;
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




      //die Abstandsmessung soll im equilibrierten Zustand beginnen:
     if ( (time >= maxtime *start_factor)  && (schalter == 0) )
      {
        Referenzabstand = pow ( pow((*dis_)[num_dof-3]+10 - (*dis_)[0],2) + pow((*dis_)[num_dof-2] - (*dis_)[1],2) , 0.5);
        schalter = 1;
        i_start = i;
        startzeit = time;
        //suche eine noch nicht existierende Dateinummer
         while(!testflag)
           {
             Dateinummer++;
             FILE* testdatei = NULL;
             std::ostringstream testname;
             testname << "EndToEnd"<< Dateinummer << ".dat";
             testdatei = fopen(testname.str().c_str(), "r");
             if (testdatei == NULL)
               testflag = 1;
           }

       filename << "EndToEnd"<< Dateinummer << ".dat";

      }
      if ( schalter == 1  && i > i_start)
      {

        endzuend = pow ( pow((*dis_)[num_dof-3]+10 - (*dis_)[0],2) + pow((*dis_)[num_dof-2] - (*dis_)[1],2) , 0.5);
        DeltaR2 = pow( endzuend - Referenzabstand ,2 );

        //writing data

        if ( (i - i_start) % int(ceil(pow( 10, floor(log10((time - startzeit) / (10*dt))))) ) == 0 )
          {
          // open file and append new data line
          fp = fopen(filename.str().c_str(), "a");
          //defining temporary stringstream variable
          std::stringstream filecontent;
          filecontent << scientific << time - startzeit << "  " << DeltaR2 << endl;
          // move temporary stringstream to file and close it
          fprintf(fp,filecontent.str().c_str());
          fclose(fp);
          }
      }

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
