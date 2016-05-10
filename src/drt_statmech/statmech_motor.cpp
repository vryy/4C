/*!----------------------------------------------------------------------
\file statmech_motor.cpp
\brief management and auxiliary functions for molecular motors

\maintainer Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276

*----------------------------------------------------------------------*/
#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../drt_inpar/inpar_statmech.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3cl.H"


/*----------------------------------------------------------------------*
 | (private) cahnge reference length of crosslinker elements            |
 |                                                         mueller 08/13|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ChangeActiveLinkerLength(const double&                    timen,
                                                         const double&                    dt,
                                                         Teuchos::RCP<Epetra_Vector>      actlinklengthin,
                                                         Teuchos::RCP<Epetra_Vector>      actlinklengthout,
                                                         const bool                       printscreen,
                                                         const bool                       revertchanges,
                                                         Teuchos::RCP<Epetra_MultiVector> bspotpositions,
                                                         Teuchos::RCP<Epetra_MultiVector> bspotrotations)
{
//  std::cout<<"Pre-ChangeLinkerLegth:\n"<<*crosslinkeractlength_<<std::endl;

  int numprobshort = 0;
  int numproblong = 0;
  if(!revertchanges)
  {
    // get current hydrolysis-rate for crosslinkers to change it's length to long or short
    double kactlinklong = 0;
    double kactlinkshort = 0;
    double starttime = actiontime_->at(1);

    if(timen <= starttime || (timen>starttime && fabs(timen-starttime) < dt/1e4))
    {
      kactlinklong = statmechparams_.get<double>("K_ACT_LONG_start",0.0);
      kactlinkshort = statmechparams_.get<double>("K_ACT_SHORT_start",0.0);
    }
    else
    {
      kactlinklong = statmechparams_.get<double>("K_ACT_LONG_end",0.0);
      kactlinkshort = statmechparams_.get<double>("K_ACT_SHORT_end",0.0);
    }

    // probability with which a crosslinker change it's length to long
    double plong = 1.0 - exp( -dt*kactlinklong );
    // probability with which a crosslinker change it's length to short
    double pshort = 1.0 - exp( -dt*kactlinkshort );

    /*the following part of the code leads to the decision, which active crosslinker should change it's length.
     * This works precisely as follows:
     *(1) the crosslink molecules are looped through
     *(2) long crosslinkers are checked, if probability to change to short is fulfilled --> if:
     *(3)   beam3-pointer on this crosslinker to change linker length of this crosslinker from long to short
     *(4) do also with short crosslinkers
     *note: a fully overlapping node map on processor 0 is imperative for the following part of the code to work correctly!*/

    Teuchos::RCP<Epetra_Vector> actlinklengthtrans = Teuchos::rcp(new Epetra_Vector(*transfermap_,true));
    CommunicateVector(actlinklengthtrans,actlinklengthout,true,false);

    // probability check, if length of active linker should be changed
    for(int i=0; i<actlinklengthtrans->MyLength(); i++)
    {
      // when multiple species of linkers are present, we skip the standard ones.
      if(crosslinkertype_!=Teuchos::null)
        if((*crosslinkertype_)[crosslinkermap_->LID(transfermap_->GID(i))]==0.0) // i.e. standard, non-active linker
          continue;

      // only for doubly-bound crosslinkers
      if((*numbond_)[crosslinkermap_->LID(transfermap_->GID(i))]>1.9)
      {
        switch((int)(*actlinklengthtrans)[i])
        {
          case 0:
          // probability check for long to short
            if((*uniformgen_)() < pshort)
            {
              numprobshort++;
              (*actlinklengthtrans)[i] = 1.0;
            }
          break;
          case 1:
            // probability check for short to long
            if((*uniformgen_)() < plong)
            {
              numproblong++;
              (*actlinklengthtrans)[i] = 0.0;
            }
          break;
          default: dserror("Wrong status %d in actlinklengthtrans", (int)(*actlinklengthtrans)[i]);
        }
      }
    }
    CommunicateVector(actlinklengthtrans, actlinklengthout, false, true);
  }

  int numprobshortglob = 0;
  int numproblongglob = 0;

  discret_->Comm().SumAll(&numprobshort, &numprobshortglob, 1);
  discret_->Comm().SumAll(&numproblong, &numproblongglob, 1);

  // for rigid body rotation of reference frame (only interpolated elements)
  Teuchos::RCP<Epetra_MultiVector> nodalquaternions = Teuchos::null;
  Teuchos::RCP<Epetra_MultiVector> bspotquaternions = Teuchos::null;
  if((linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_myosinthick) &&
      DRT::INPUT::IntegralValue<int>(statmechparams_,"CROSSBRIDGEMODEL"))
  {
    nodalquaternions = Teuchos::rcp(new Epetra_MultiVector(*(discret_->NodeColMap()),4));
    bspotquaternions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4));
    GetElementBindingSpotTriads(nodalquaternions);
    GetInterpolatedBindingSpotTriads(bspotrotations,bspotquaternions);
  }

  //TODO
//  if((int)(round(timen/dt)-1.0)%2!=0)
//  {
//    cout<<"Change back to long"<<endl;
//    (*actlinklengthout)[0] = 0.0;
//  }

  // store linker length change
  int toshort=0;
  int tolong=0;
  // change reference length in actual elements
  for(int i=0; i<actlinklengthout->MyLength(); i++)
  {
    int collid = discret_->ElementColMap()->LID((int)(*crosslink2element_)[i]);

    // if length of the active linker changed in the probability check
    if(fabs((*actlinklengthout)[i]-(*actlinklengthin)[i])>1e-6 && (int)(*crosslink2element_)[i]>-0.9 && collid != -1)
    {
      // just on the processor, where element is located
      int rowlid = discret_->ElementRowMap()->LID((int)(*crosslink2element_)[i]);
      double sca = statmechparams_.get<double>("LINKERSCALEFACTOR", 0.8);

      switch((int)(*actlinklengthout)[i])
      {
        // change linker length to short
        case 1:
        {
          switch(linkermodel_)
          {
            case statmech_linker_active:
              (dynamic_cast<DRT::ELEMENTS::Beam3*>(discret_->lColElement(collid)))->SetReferenceLength(sca);
            break;
            case statmech_linker_activeintpol:
            case statmech_linker_myosinthick:
              ActiveLinkerPowerStroke(collid,i,(int)(*actlinklengthout)[i],revertchanges, bspotpositions, bspotquaternions,nodalquaternions);
            break;
            default: dserror("Unknown active linker beam element!");
          }
          if(rowlid!=-1)
            toshort++;
        }
        break;
        // change linker length to long
        case 0:
        {
          switch(linkermodel_)
          {
            case statmech_linker_active:
              (dynamic_cast<DRT::ELEMENTS::Beam3*>(discret_->lColElement(collid)))->SetReferenceLength(1.0/sca);
            break;
            case statmech_linker_activeintpol:
            case statmech_linker_myosinthick:
              ActiveLinkerPowerStroke(collid, i, (int)(*actlinklengthout)[i], revertchanges,bspotpositions, bspotquaternions,nodalquaternions);
            break;
            default: dserror("Unknown active linker beam element!");
          }
          if(rowlid!=-1)
            tolong++;
        }
        break;
        default: dserror("Wrong status %d in actlinklength_", (int)(*crosslinkeractlength_)[i]);
      }
    }
  }

  int toshortglob = 0;
  int tolongglob = 0;

  discret_->Comm().SumAll(&toshort, &toshortglob, 1);
  discret_->Comm().SumAll(&tolong, &tolongglob,1);

  // number of changed crosslinker lengths
  if(revertchanges)
  {
    if(!discret_->Comm().MyPID() && printscreen)
      std::cout<<"\nRevert reference lengths..."<<std::endl;
  }
  else
  {
    // update cycle time
    for(int i=0; i<crosslinkeractcycletime_->MyLength(); i++)
      (*crosslinkeractcycletime_)[i] += dt;
    Teuchos::RCP<Epetra_Vector> crosslinkeractcycletimetrans = Teuchos::rcp(new Epetra_Vector(*transfermap_,true));
    CommunicateVector(crosslinkeractcycletimetrans,crosslinkeractcycletime_);

    if(!discret_->Comm().MyPID() && (toshortglob+tolongglob)>0 && printscreen)
    {
      std::cout<<"\n"<<"--Crosslinker length changes--"<<std::endl;
      std::cout<<" - long to short: "<<toshortglob<<std::endl;
      std::cout<<"   - due to prob: "<<numprobshortglob<<std::endl;
      std::cout<<" - short to long: "<<tolongglob<<std::endl;
      std::cout<<"   - due to prob: "<<numproblongglob<<std::endl;
      std::cout<<"------------------------------"<<std::endl;
    }
  }
//    std::cout<<"Post-ChangeLinkerLegth:\n"<<*crosslinkeractlength_<<std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | rotate active linker reference position in addition to contraction   |
 |                                               (private) mueller 03/14|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::ActiveLinkerPowerStroke(const int&                       elecollid,
                                                        const int&                       crosslid,
                                                        const int&                       conformation,
                                                        const bool                       revertchanges,
                                                        Teuchos::RCP<Epetra_MultiVector> bspotpositions,
                                                        Teuchos::RCP<Epetra_MultiVector> bspotquaternions,
                                                        Teuchos::RCP<Epetra_MultiVector> nodalquaternions)
{
 /*This method changes the reference configuration of the linker beam element.
 * It scales the reference length by a factor "sca". Additionally, it rotates the
 * reference orientation, which is the result of the latest converged time step,
 * by an angle "thetarot".
 * Both quantities "sca" and "thetarot" are computed such that the ensuing restoring forces
 * and moments lead to a relaxed configuration, which assumes a covered distance delta along the inverse binding spot tangent (due to positive polarity)
 * */
  // default length scale factor (going to be adjusted below)
  double sca = statmechparams_.get<double>("LINKERSCALEFACTOR", 0.8);
  DRT::Element* activeele = discret_->lColElement(elecollid);

  // more complex cross bridge model
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"CROSSBRIDGEMODEL"))
  {
    // difference vector between the two binding sites involved as well as its parallel and orthogonal projection with respect to the tangent of the binding site
    LINALG::Matrix<3,1> linkdirstart;
    LINALG::Matrix<3,1> linkdirend;

    // rotation angle as pseudovector
    LINALG::Matrix<3,1> thetarot;

    // short to long (back swing; usually not required for swinging crossbridge model)
    if(conformation==0)
    {
      // revertchanges is only TRUE, when the iterative scheme has diverged. Then, all values are reset to the latest converged state.
      // We don't have to do anything then, the reset is done by discret_->Evaluate()!!!
      if(!revertchanges)
      {
        // note: this only works because we an active linker always starts out long and then contracts/swings. So, when arriving here, it's
        // supposed to be always short!!!
        // get inverse rotation
        LINALG::Matrix<4,1> invQtheta = LARGEROTATIONS::inversequaternion((dynamic_cast<DRT::ELEMENTS::BeamCL*>(activeele))->Qrot());
        LARGEROTATIONS::quaterniontoangle(invQtheta,thetarot);
        // reversion of the contraction
        sca = 1.0/(dynamic_cast<DRT::ELEMENTS::BeamCL*>(activeele))->LScaleFac();
      }
    }
    else
    {
      // retrieve hypotenuse vector (vector from one binding spot to the other)
      // get binding spot LIDs
      int lid0 = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][crosslid]);
      int lid1 = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][crosslid]);

      for(int j=0;j<(int)linkdirstart.M(); j++)
        linkdirstart(j) = (*bspotpositions)[j][lid1]-(*bspotpositions)[j][lid0];

      // calculate vector decomposition of hypotenuse for subsequent calculation of the rotation pseudo vector theta
      LINALG::Matrix<3,1> bspottangent(true);
      LINALG::Matrix<3,3> bspottriad = ExtractTriadFromGlobalQuaternionVector(lid1, bspotquaternions);

      for(int j=0; j<(int)bspottangent.M(); j++)
        bspottangent(j) = bspottriad(j,0);
      bspottangent.Scale(1.0/bspottangent.Norm2());

      // calculate end position assuming polarity criterion is met (i.e. the stroke direction is -bspottangent
      // calculate remaining distance to be covered
      linkdirend = bspottangent;
      linkdirend.Scale(-statmechparams_.get<double>("STROKEDISTANCE", 0.005));
      linkdirend += linkdirstart;

      // ensure parallel trajectory of filament by an additional corresponding contraction
      sca = linkdirend.Norm2()/linkdirstart.Norm2();
      // security measure in case of no contraction
      if(fabs(sca-1.0)<1e-9)
        sca = 1.0;

      thetarot(0) = linkdirstart(1)*linkdirend(2)-linkdirstart(2)*linkdirend(1);
      thetarot(1) = linkdirstart(2)*linkdirend(0)-linkdirstart(0)*linkdirend(2);
      thetarot(2) = linkdirstart(0)*linkdirend(1)-linkdirstart(1)*linkdirend(0);
      if(thetarot.Norm2()>1e-16)
        thetarot.Scale(1.0/thetarot.Norm2());
      // calculate rotation angle
      //rotation angle in radian
      linkdirstart.Scale(1.0/linkdirstart.Norm2());
      linkdirend.Scale(1.0/linkdirend.Norm2());
      double thetaabs = acos(linkdirstart.Dot(linkdirend));
      // make theta a pseudo vector
      thetarot.Scale(thetaabs);
    }

    // Set converged quaternions to rotated values. Currently only direction long->short, no inverse operation.
    // set new reference length
    if(!revertchanges)
    {
      (dynamic_cast<DRT::ELEMENTS::BeamCL*>(activeele))->SetReferenceLength(sca);
      (dynamic_cast<DRT::ELEMENTS::BeamCL*>(activeele))->SetRotation(thetarot);
    }
  }
  else // standard case without rotation
  {
    //TODO
    if(conformation==1)
      (dynamic_cast<DRT::ELEMENTS::BeamCL*>(activeele))->SetReferenceLength(sca);
    else
      (dynamic_cast<DRT::ELEMENTS::BeamCL*>(activeele))->SetReferenceLength(1.0/sca);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  filament polarity check for active linkers upon bond breaking       |
 |                                               (private) mueller 08/13|
 *----------------------------------------------------------------------*/
int STATMECH::StatMechManager::LinkerPolarityCheckDetach(Teuchos::RCP<Epetra_MultiVector>       punlink,
                                                          const Teuchos::RCP<Epetra_MultiVector> bspotpositions,
                                                          const Teuchos::RCP<Epetra_MultiVector> bspottriadscol)
{
  if(DRT::INPUT::IntegralValue<int>(statmechparams_,"FILAMENTPOLARITY"))
  {

    int numdetach = 0;

    if(!discret_->Comm().MyPID())
    {
      double recoverytime = (statmechparams_.get<double>("ACTIVERECOVERYFRACTION",0.95))*(statmechparams_.get<double>("ACTIVELINKERCYCLE",0.04));
      double detachtime = (1.0 - (statmechparams_.get<double>("ACTIVERECOVERYFRACTION",0.95)))*(statmechparams_.get<double>("ACTIVELINKERCYCLE",0.04));

      for(int i=0; i<crosslinkermap_->NumMyElements(); i++)
      {
        // there exists a crosslinker element
        if((*crosslink2element_)[i]>-0.9)
        {
          // skip this linker if its not an active one
          if(statmechparams_.get<double>("ACTIVELINKERFRACTION",0.0)>0.0)
            if((*crosslinkertype_)[i]==0.0)
              continue;

          std::vector<int> order = Permutation(2);
          // prescribed order because only filament-linker bonds may be broken for motility assays
          if(networktype_==statmech_network_motassay)
          {
            order[0] = 0;
            order[1] = 1;
          }
          for(int j=0; j<crosslinkerbond_->NumVectors(); j++)
          {
            int bspotlid0 = bspotcolmap_->LID((int)(*crosslinkerbond_)[order[j]][i]);
            int bspotlid1 = bspotcolmap_->LID((int)(*crosslinkerbond_)[!order[j]][i]);

            // check for periodic boundary conditions
            std::vector<double> pos(6,0.0);
            for(int k=0; k<3; k++)
            {
              pos[k] = (*bspotpositions)[k][bspotlid0];
              pos[k+3] = (*bspotpositions)[k][bspotlid1];
            }
            UnshiftPositions(pos,2);

            LINALG::Matrix<3,1> linkdir(true);
            for(int k=0; k<(int)linkdir.M(); k++)
              linkdir(k) = pos[k+3]-pos[k];
            linkdir.Scale(1.0/linkdir.Norm2());

            LINALG::Matrix<3,3> bspottriad(true);
            LINALG::Matrix<4, 1> qnode(true);
            LINALG::Matrix<3,1> tangent(true);

            for (int k=0; k<4; k++)
              qnode(k) = (*bspottriadscol)[k][bspotlid1];
            LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);
            for (int k=0; k<(int)bspottriad.M(); k++)
              tangent(k) = bspottriad(k,0);
            tangent.Scale(1/tangent.Norm2());


        //    std::cout<<"direct : "<<direction(0)<<" "<<direction(1)<<" "<<direction(2)<<std::endl;
        //    std::cout<<"linkdir: "<<linkdir(0)<<" "<<linkdir(1)<<" "<<linkdir(2)<<std::endl;
        //    std::cout<<"tangent: "<<tangent(0)<<" "<<tangent(1)<<" "<<tangent(2)<<std::endl;
        //    std::cout<<"  pardir  : "<<linkpardir(0)<<" "<<linkpardir(1)<<" "<<linkpardir(2)<<std::endl;
        //    std::cout<<"  orthodir: "<<linkorthodir(0)<<" "<<linkorthodir(1)<<" "<<linkorthodir(2)<<std::endl;
            double polarityscale = linkdir.Dot(tangent);
            // For now, manipulate j==1 since this is the free filament binding site
            if (polarityscale < 0.0 || ((*crosslinkeractcycletime_)[i] >= recoverytime+detachtime && (*crosslinkeractlength_)[i]>0.9))
            {
              numdetach++;
              (*punlink)[!order[j]][i]=1.0;
            }
            // break after evaluating the first bond
            break;
          }
        }
      }
    }

    Teuchos::RCP<Epetra_MultiVector> punlinktrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_,2,true));
    CommunicateMultiVector(punlinktrans,punlink);

    return numdetach;
  }
  else
    return 0;
}// StatMechManager::LinkerPolarityCheckDetach()

/*----------------------------------------------------------------------*
 |  filament polarity check for active linkers upon bond establishment  |
 |                                             (private)   mueller 08/13|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::LinkerPolarityCheckAttach(Teuchos::RCP<Epetra_MultiVector> bspottriadscol,
                                                          Epetra_SerialDenseMatrix&      LID,
                                                          LINALG::Matrix<3,1>&           direction)
{
  // cycle time check:
  // time during which an active linker is in recovery conformation (default values from myosin cycle)
  double recoverytime = (statmechparams_.get<double>("ACTIVERECOVERYFRACTION",0.95))*(statmechparams_.get<double>("ACTIVELINKERCYCLE",0.04));
  int crosslid = (*bspotstatus_)[(int)LID(0,0)];
  if((*crosslinkeractcycletime_)[crosslid] >= recoverytime)
  {
    // 1. polarity criterion (2D): filaments are not allowed to link from too far away + angle criterion
    // retrieve tangential vector from binding spot quaternions
    LINALG::Matrix<3,3> bspottriad = ExtractTriadFromGlobalQuaternionVector((int)LID(1,0),bspottriadscol);
    // unit tangent corresponding to first column of bspottriad
    LINALG::Matrix<3,1> tangent(true);
    for (int l=0; l<(int)bspottriad.M(); l++)
      tangent(l) = bspottriad(l,0);
    tangent.Scale(1.0/tangent.Norm2());

    // direction vector from linker->binding spot, parallel/orthogonal projection of that direction w. respect to binding spot tangent
    LINALG::Matrix<3,1> linkdir(true);
    LINALG::Matrix<3,1> linkpardir(true);
    LINALG::Matrix<3,1> linkorthodir(true);
    linkdir -= direction; // according to convention X(LID_1)-X(LID_0), note: direction is unscaled, i.e., carries distance info

    DoVectorDecomposition(linkdir,linkpardir,linkorthodir, tangent);

//    std::cout<<"direct : "<<direction(0)<<" "<<direction(1)<<" "<<direction(2)<<std::endl;
//    std::cout<<"linkdir: "<<linkdir(0)<<" "<<linkdir(1)<<" "<<linkdir(2)<<std::endl;
//    std::cout<<"tangent: "<<tangent(0)<<" "<<tangent(1)<<" "<<tangent(2)<<std::endl;
//    std::cout<<"  pardir  : "<<linkpardir(0)<<" "<<linkpardir(1)<<" "<<linkpardir(2)<<std::endl;
//    std::cout<<"  orthodir: "<<linkorthodir(0)<<" "<<linkorthodir(1)<<" "<<linkorthodir(2)<<std::endl;

    // scale to unit length for later use
    direction.Scale(1.0/direction.Norm2());

    // 1. criterion: angle
    double alpha = acos(linkorthodir.Norm2()/linkdir.Norm2());

    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"CROSSBRIDGEMODEL") && alpha > statmechparams_.get<double>("PHIBSPOT", 1.0472))
      return false;
    else if(!(DRT::INPUT::IntegralValue<int>(statmechparams_,"CROSSBRIDGEMODEL")) && (alpha > statmechparams_.get<double>("PHIBSPOT", 1.0472) || alpha<(statmechparams_.get<double>("PHIBSPOT", 1.0472)/10.0)))
      return false;

    // 2. criterion: polarity
    linkdir.Scale(1.0/linkdir.Norm2());
    double polarityscale = linkdir.Dot(tangent);
    // zero-positive scale factor signals parallelity, negative scale factor signals antiparallelity
    if (polarityscale >= 0)
      return true;
    else
      return false;
  }
  else
    return false;
}
