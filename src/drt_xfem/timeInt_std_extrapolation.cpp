/*!------------------------------------------------------------------------------------------------*
\file timeInt_std_extrapolation.cpp

\brief provides the Extrapolation class

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
 *------------------------------------------------------------------------------------------------*/


#include "timeInt_std_extrapolation.H"
#include "../drt_combust/combust_flamefront.H"
#include "../drt_combust/combust_defines.H"


/*------------------------------------------------------------------------------------------------*
 * Extrapolation constructor                                                     winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
XFEM::ExtrapolationOld::ExtrapolationOld(
    XFEM::TIMEINT& timeInt,
    INPAR::COMBUST::XFEMTimeIntegration timeIntType,
    const Teuchos::RCP<Epetra_Vector> veln,
    const double& dt,
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    const double& veljump,
    bool initialize
) : STD(timeInt,timeIntType,veln,dt,flamefront,initialize),
veljump_(veljump)
{
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * call the computation based on an extrapolation                                winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationOld::compute(
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
)
{
  if (FGIType_==FRSNot1_)
    return;

  handleVectors(newRowVectorsn,newRowVectorsnp);

  resetState(TimeIntData::basicStd_,TimeIntData::extrapolateStd_);

  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
#ifndef COMBUST_SETJUMP
    extrapolationMain(&*data);
#else
    setJump(&*data);
#endif
  }

#ifdef PARALLEL
  exportFinalData();
#endif
  setFinalData();

  // fill vectors with the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_!=TimeIntData::doneStd_)
      dserror("All data should be set here, having status 'done'. Thus something is wrong!");
  }
} // end compute



/*------------------------------------------------------------------------------------------------*
 * call the computation based on an extrapolation                                winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationOld::compute(std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectors)
{
  if (oldVectors_.size() != newRowVectors.size())
  {
    std::cout << "sizes are " << oldVectors_.size() << " and " << newRowVectors.size() << std::endl;
    dserror("Number of state-vectors at new and old discretization are different!");
  }

  newVectors_ = newRowVectors;

  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
#ifndef COMBUST_SETJUMP
    extrapolationMain(&*data);
#else
    setJump(&*data);
#endif
  }

#ifdef PARALLEL
  exportFinalData();
#endif
  setFinalData();

  // fill vectors with the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_!=TimeIntData::doneStd_)
      dserror("All data should be set here, having status 'done'. Thus something is wrong!");
  }
}// end compute



/*------------------------------------------------------------------------------------------------*
 * extrapolation of values to data-requiring nodes:                              winklmaier 10/11 *
 * a straight line through the data-requiring node and the nearest node on the correct interface  *
 * side is set up, two appropriate points on this line are used for extrapolation                 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationOld::extrapolationMain(
    TimeIntData* data
)
{
  if (data->state_!=TimeIntData::extrapolateStd_)
    dserror("extrapolation must be called with according state");

  const int nsd = 3;

  DRT::Element* startele = NULL; // element of startpoint
  DRT::Element* midele = NULL; // element of midpoint

  LINALG::Matrix<nsd,1> startpoint(true); // global coordinates of startpoint
  LINALG::Matrix<nsd,1> xistartpoint(true); // local coordinates of startpoint

  LINALG::Matrix<nsd,1> midpoint(true); // global coordinates of midpoint
  LINALG::Matrix<nsd,1> ximidpoint(true); // local coordinates of midpoint

  LINALG::Matrix<nsd,1> endpoint(data->node_.X());

//    std::cout << "searching node = endpoint = " << endpoint << std::endl;
  // identify final element and local coordinates of startpoint and midpoint for extrapolation
  bisection(data,startele,startpoint,xistartpoint,midele,midpoint,ximidpoint);

//    std::cout << std::endl << "startpoint is " << startpoint;
//    std::cout << "midpoint is " << midpoint;

  // compute the constants for the extrapolation:
  // value = c1*valuestartpoint + c2*valuemidpoint
  LINALG::Matrix<nsd,1> dist1;
  LINALG::Matrix<nsd,1> dist2;

  dist1.Update(1.0,midpoint,-1.0,startpoint);
  dist2.Update(1.0,endpoint,-1.0,midpoint);
//  std::cout << "distance from startpoint to midpoint is " << dist1;
//  std::cout << "distance from midpoint to endpoint is " << dist2;

  if (dist1.Norm2()<1e-14 or dist2.Norm2()<1e-14)
    dserror("something wrong in bisection");

  double c1 = - dist2.Norm2()/dist1.Norm2(); // 1 + dist2/dist1
  double c2 = 1.0 + dist2.Norm2()/dist1.Norm2(); // dist2/dist1

  // get the velocities and the pressures at the start- and midpoint
  std::vector<LINALG::Matrix<nsd,1> > velstartpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));
  std::vector<LINALG::Matrix<nsd,1> > velmidpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));

  std::vector<double> presstartpoint(oldVectors_.size(),0.0);
  std::vector<double> presmidpoint(oldVectors_.size(),0.0);

  callInterpolation(startele,xistartpoint,velstartpoint,presstartpoint);
  callInterpolation(midele,ximidpoint,velmidpoint,presmidpoint);

  //  std::cout << "pres at startpoint is " << presstartpoint[0];
  //  std::cout << "pres at midpoint is " << presmidpoint[0];

  // compute the final velocities and pressure due to the extrapolation
  std::vector<LINALG::Matrix<nsd,1> > velendpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));
  std::vector<double> presendpoint(oldVectors_.size(),0.0);

  for (size_t index=0;index<oldVectors_.size();index++)
  {
    velendpoint[index].Update(c1,velstartpoint[index],c2,velmidpoint[index]);
    presendpoint[index] = c1*presstartpoint[index] + c2*presmidpoint[index];

    //if (index == 0)
    //std::cout << "final pressure is " << presendpoint[index] <<
    //" and final velocity is " << velendpoint[index];

  } // loop over vectors to be set

  data->startOwner_ = std::vector<int>(1,myrank_);
  data->velValues_ = velendpoint;
  data->presValues_ = presendpoint;
  data->state_ = TimeIntData::doneStd_;
}



void XFEM::ExtrapolationOld::setJump(
    TimeIntData* data
)
{
  DRT::Node node = data->node_;

  // get the global id for current node
  const int gid = node.Id();

  // get local processor id according to global node id
  const int nodelid = (*gradphi_).Map().LID(gid);
  if (nodelid<0) dserror("Proc %d: Cannot find gid=%d in Epetra_Vector",(*gradphi_).Comm().MyPID(),gid);

  const int numcol = (*gradphi_).NumVectors();
  if( numcol != 3) dserror("number of columns in Epetra_MultiVector is not 3");

  //-------------------------------------------------------------
  // get (smoothed) gradient of the G-function field at this node
  //-------------------------------------------------------------
  LINALG::Matrix<3,1> mygradphi(true);

  // loop over dimensions (= number of columns in multivector)
  for(int col=0; col< numcol; col++)
  {
    // get columns vector of multivector
    double* globalcolumn = (*gradphi_)[col];
    // set smoothed gradient entry of phi into column of global multivector
    mygradphi(col) = globalcolumn[nodelid];
  }

  //------------------------------------
  // compute smoothed unit normal vector
  // n = - grad phi / |grad phi|
  //------------------------------------
  // smoothed normal vector at this node
  LINALG::Matrix<3,1> nvec = mygradphi;
  // compute norm of smoothed normal vector
  const double ngradphi = mygradphi.Norm2();

  // divide vector by its norm to get unit normal vector
  if (fabs(ngradphi < 1.0E-12))// 'ngradnorm' == 0.0
  {
    // length of smoothed normal is zero at this node -> node must be on a singularity of the
    // level set function (e.g. "regular level set cone"); all normals add up to zero normal vector
    // -> The fluid convective velocity 'fluidvel' alone constitutes the flame velocity, since the
    //    relative flame velocity 'flvelrel' turns out to be zero due to the zero average normal vector.
    std::cout << "\n/!\\ phi gradient too small at node " << gid << " -> flame velocity is only the convective velocity" << std::endl;
    nvec.PutScalar(0.0);
  }
  else
  {
    nvec.Scale(-1.0/ngradphi);
  }

  LINALG::Matrix<3,1> veljump(true);
  veljump.Update(veljump_,nvec);

  double wallfac = 1.0;
#ifdef ORACLES
  //--------------------------------------------------------
  // check physical coordinates of node
  // remark: we want to blend the flame speed close to walls
  //--------------------------------------------------------

  //    wall
  // 1.0 |     ______
  //     |    /
  // 0.0 |___/ ,
  //     |     H/6
  const double wallzone = 0.0299/6.0;
  if (node.X()[0] > 0.0) // inside combustion chamber
  {
    if ( (0.0653-abs(node.X()[1])) < wallzone or // close to top or bottom wall
                       node.X()[0] < wallzone )  // close to step
    {
      // wall factor is 0 at the wall and 1 at H/6 or further away from the wall
      wallfac = 6.0/0.0299 * std::min(0.0653-abs(node.X()[1]),node.X()[0]);
      if (wallfac < 0.1) // cut off the last 10% to guarantee a zero jump at the wall
        wallfac = 0.0;
    }
  }
#endif
  veljump.Scale(wallfac);

  // node velocities of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<3,1> > nodeveldata(oldVectors_.size(),LINALG::Matrix<3,1>(true));
  // node pressures of the element nodes for the data that should be changed
  std::vector<double> nodepresdata(oldVectors_.size(),0.0);

  // get nodal velocities and pressures with help of the field set of node
  const std::set<XFEM::FieldEnr>& fieldEnrSet(olddofman_->getNodeDofSet(node.Id()));
  for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
      fieldenr != fieldEnrSet.end();++fieldenr)
  {
    const DofKey olddofkey(gid, *fieldenr);
    const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;
    if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
    {
      if (fieldenr->getField() == XFEM::PHYSICS::Velx)
      {
        for (size_t index=0;index<oldVectors_.size();index=index+2)
          nodeveldata[index](0) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)] - veljump(0)*interfaceSide(data->phiValue_);
      }
      else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
      {
        for (size_t index=0;index<oldVectors_.size();index=index+2)
          nodeveldata[index](1) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)] - veljump(1)*interfaceSide(data->phiValue_);
      }
      else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
      {
        for (size_t index=0;index<oldVectors_.size();index=index+2)
          nodeveldata[index](2) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)] - veljump(2)*interfaceSide(data->phiValue_);
      }
      else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
      {
        for (size_t index=0;index<oldVectors_.size();index=index+2)
          nodepresdata[index] = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
      }
      else
      {
        //std::cout << XFEM::PHYSICS::physVarToString(fieldenr->getField()) << std::endl;
        //dserror("not implemented physical field!");
      }
    }
  } // end loop over fieldenr

  data->velValues_ = nodeveldata;
  data->presValues_ = nodepresdata;
  data->state_ = TimeIntData::doneStd_;
}



/*------------------------------------------------------------------------------------------------*
 * perform a bisection on the line startpoint-endpoint                           winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationOld::bisection(
    TimeIntData* data,
    DRT::Element*& startele,
    LINALG::Matrix<3,1>& startpoint,
    LINALG::Matrix<3,1>& xistartpoint,
    DRT::Element*& midele,
    LINALG::Matrix<3,1>& midpoint,
    LINALG::Matrix<3,1>& ximidpoint
)
{
  // general initialization
  const int nsd = 3;

  DRT::Node* node = discret_->gNode(data->startGid_[0]); // startpoint node on current proc
  LINALG::Matrix<nsd,1> nodecoords(node->X());

  std::vector<const DRT::Element*> eles;
  addPBCelements(node,eles);
  const int numele=eles.size();

  startpoint = nodecoords; // first point used for extrapolation
  midpoint = nodecoords; // second point used for extrapolation
  LINALG::Matrix<nsd,1> endpoint(data->node_.X()); // coordinates of data-requiring node
//    std::cout << std::endl << "initial startpoint and midpoint is " << nodecoords;

  LINALG::Matrix<nsd,1> pointTmp; // temporarily point used for computations and for bisection
  pointTmp.Update(0.5,startpoint,0.5,endpoint); // midpoint of startpoint and endpoint

  LINALG::Matrix<nsd,1> dist;
  dist.Update(1.0,endpoint,-1.0,startpoint);

  LINALG::Matrix<nsd,1> xipointTmp(true);

  int iter = 0;
  const int crit_max_iter = 100; // maximal iteration number for critical cases
  int curr_max_iter = 10; // currently used iteration number
  int std_max_iter = 10; // iteration number usually sufficient
  bool elefound = false;

  DRT::Element* eletmp = NULL;
  bool elefoundtmp = false;

  // prework: search for the element around "startpoint" into "endpoint"-direction
  while(true)
  {
//    std::cout << "potential midpoint is " << pointTmp;
    iter++;
    for (int i=0; i<numele; i++)
    {
      eletmp = (DRT::Element*)eles[i];
      callXToXiCoords(eletmp,pointTmp,xipointTmp,elefoundtmp);
      if (elefoundtmp) // usually loop can stop here, just in special cases of critical cuts problems might happen
      {
        midele = eletmp;
        elefound = true;
        if (intersectionStatus(midele)==XFEM::TIMEINT::cut_) // really cut -> ok
          break;
        else if (intersectionStatus(midele)==XFEM::TIMEINT::uncut_) // really uncut -> ok
          break;
        else // special cases -> try to take another element
          ; // do not break, but potential element is saved
      }
    }

    if (elefound)
      break;
    else // corresponds to CFL > 1 -> not very good, but possible
      pointTmp.Update(-pow(0.5,iter+1),dist,1.0);

    if (iter==crit_max_iter) break;
  }

  curr_max_iter = crit_max_iter;

  // search for midpoint with bisection of the straight line:
  //distances p1-p2 and p2-pend as even as possible
  //-> p2 near midpoint of p1 and pend as far as possible
  //gives back either a point between pend and p1 or p2 = p1
  for (int i=iter;i<=curr_max_iter;i++)
  {
//    std::cout << "potential midpoint is " << pointTmp;
    callXToXiCoords(midele,pointTmp,xipointTmp,elefound);

    if (!elefound)
      dserror("two points of a line segment in one element "
          "-> also any point between them should be in the same element.");

    if (interfaceSideCompare(midele,pointTmp,0,data->phiValue_) == true) // Lagrangian origin and original node on different interface sides
    {
//      std::cout << "current midpoint is " << pointTmp;
      curr_max_iter = std_max_iter; // usable point found
      midpoint = pointTmp; // possible point
      ximidpoint = xipointTmp;
      pointTmp.Update(pow(0.5,i+1),dist,1.0); // nearer to endpoint
    }
    else
      pointTmp.Update(-pow(0.5,i+1),dist,1.0); // nearer to startpoint
  }

  // if element is (nearly) touched, the above bisection might fail
  if (midpoint == nodecoords) // nothing was changed above
  {
    std::cout << std::endl << std::endl << "WARNING: this case should no more happen!" << std::endl << std::endl;
    midele = (DRT::Element*)eles[0];
    callXToXiCoords(midele,midpoint,ximidpoint,elefound);
  }
//  std::cout << "final midpoint is " << midpoint << " in " << *midele;

  // get current distances of the three points
  LINALG::Matrix<nsd,1> dist1; // distance from startpoint to midpoint
  LINALG::Matrix<nsd,1> dist2; // distance from midpoint to endpoint

  dist1.Update(1.0,midpoint,-1.0,startpoint);
  dist2.Update(1.0,endpoint,-1.0,midpoint);

  curr_max_iter = crit_max_iter;

  // compute potentially final startpoint
  if ((dist1.Norm2()/dist.Norm2() > 1.0/3.0) and (dist2.Norm2()/dist.Norm2() > 1.0/3.0)) // distances are ok
  {
    startele = midele;
    callXToXiCoords(startele,startpoint,xistartpoint,elefound);

    if (!elefound)
      dserror("extrapolation should be done within one convex element. Thus this should work");

    if (interfaceSideCompare(startele,startpoint,0,data->phiValue_) == false) // Lagrangian origin and original node on different interface sides
      dserror("node which was evaluated to be the nearest node on correct side is on wrong side???");
  }
  else if (dist2.Norm2()/dist.Norm2() <= 1.0/3.0) // midpoint much nearer at endpoint than at startpoint
  {
    startpoint.Update(1.0,endpoint,-2.0,dist2); // same, moderate distances
    startele = midele;
    callXToXiCoords(startele,startpoint,xistartpoint,elefound);

    if (!elefound)
      dserror("extrapolation should be done within one convex element. Thus this should work");

    if (interfaceSideCompare(startele,startpoint,0,data->phiValue_) == false) // Lagrangian origin and original node on different interface sides
    {
      LINALG::Matrix<nsd,1> startpointLeft = startpoint;
      LINALG::Matrix<nsd,1> xistartpointLeft(true);
      LINALG::Matrix<nsd,1> startpointRight = startpoint;
      LINALG::Matrix<nsd,1> xistartpointRight(true);

      // get the potential startpoint on the "node"-side of the current, not usable startpoint
      dist.Update(1.0,startpoint,-1.0,nodecoords);
      pointTmp.Update(0.5,startpoint,0.5,nodecoords);
      for (int i=1;i<=curr_max_iter;i++)
      {
        callXToXiCoords(startele,pointTmp,xipointTmp,elefound);

        if (interfaceSideCompare(startele,pointTmp,0,data->phiValue_) == true) // Lagrangian origin and original node on different interface sides
        {
          curr_max_iter = std_max_iter;
          startpointLeft = pointTmp; // possible point
          pointTmp.Update(pow(0.5,i+1),dist,1.0); // nearer to optimal startpoint
        }
        else
          pointTmp.Update(-pow(0.5,i+1),dist,1.0); // nearer to not-optimal, but possible node
      }

      curr_max_iter = crit_max_iter;
      // get the potential startpoint on the "midpoint"-side of the current, not usable startpoint
      dist.Update(1.0,midpoint,-1.0,startpoint);
      pointTmp.Update(0.5,midpoint,0.5,startpoint);
      for (int i=1;i<=curr_max_iter;i++)
      {
        callXToXiCoords(startele,pointTmp,xipointTmp,elefound);

        if (interfaceSideCompare(startele,pointTmp,0,data->phiValue_) == true) // Lagrangian origin and original node on different interface sides
        {
          curr_max_iter = std_max_iter;
          startpointRight = pointTmp; // possible point
          pointTmp.Update(-pow(0.5,i+1),dist,1.0); // nearer to original startpoint
        }
        else
          pointTmp.Update(+pow(0.5,i+1),dist,1.0); // nearer to midpoint, worse ratio, but possible
      }

      // check if startpoints changed
      if (startpointLeft == startpoint) // left possible startpoint didn't change
        startpointLeft = nodecoords; // possible, but not very good point due to high ratio
      if (startpointRight == startpoint) // right possible startpoint didn't change
        startpointRight = nodecoords; // possible, but not very good point due to high ratio

      // check which startpoint is the better one
      dist.Update(1.0,midpoint,-1.0,startpointLeft);
      double ratioLeft = dist1.Norm2()/dist.Norm2(); // dist is the greater entry in this case

      dist.Update(1.0,midpoint,-1.0,startpointRight);
      double ratioRight =  dist.Norm2()/dist1.Norm2(); // dist1 is the greater entry in this case

      if (ratioLeft > ratioRight) // nearer-to-one ratio = better
        startpoint = startpointLeft;
      else
        startpoint = startpointRight;

      // compute final local coordinates of the startpoint
      callXToXiCoords(startele,startpoint,xistartpoint,elefound);
    }
  }
  else // distance startpoint - midpoint much small than distance midpoint - endpoint
  {
    pointTmp.Update(1.0,endpoint,-1.5,dist2); // dist2 = 2*dist1 -> moderate factor and moderate total lenght
    dist.Update(1.0,midpoint,-1.0,pointTmp);

    iter = 0;
    // possibly the startpoint is too far away from the "node"
    // -> reduce distance until startpoint is in element adjacent to the "node"
    while(true)
    {
      iter++;
      for (int i=0; i<numele; i++)
      {
        startele = (DRT::Element*)eles[i];
        callXToXiCoords(startele,pointTmp,xipointTmp,elefound);
        if (elefound) break;
      }

      if (elefound) break;
      if (iter==curr_max_iter) break;

      pointTmp.Update(+pow(0.5,iter),dist,1.0);
    }

    if (!elefound)
      dserror("extrapolation should be done within one convex element. Thus this should work");

    // get the potential startpoint on the "node"-side of the current, not usable startpoint
    for (int i=iter;i<=curr_max_iter+1;i++)
    {
      if (!elefound)
        pointTmp.Update(+pow(0.5,i),dist,1.0); // nearer to midpoint
      else
      {
        if (interfaceSideCompare(startele,pointTmp,0,data->phiValue_) == true) // Lagrangian origin and original node on different interface sides
        {
          curr_max_iter = std_max_iter;
          startpoint = pointTmp; // possible point
          xistartpoint = xipointTmp;
          pointTmp.Update(-pow(0.5,i),dist,1.0); // nearer to original startpoint
        }
        else
          pointTmp.Update(+pow(0.5,i),dist,1.0); // nearer to midpoint
      }

      if (i<=curr_max_iter)
      {
        callXToXiCoords(startele,pointTmp,xipointTmp,elefound);
        if (!elefound)
        {
          for (int i=0; i<numele; i++)
          {
            startele = (DRT::Element*)eles[i];
            callXToXiCoords(startele,pointTmp,xipointTmp,elefound);
            if (elefound) break;
          }
        }
      }
    }
  }
//  std::cout << "final startpoint is " << startpoint << " in " << *startele;

} // end bisection



/*------------------------------------------------------------------------------------------------*
 * call the interpolation                                                        winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationOld::callInterpolation(
    DRT::Element* ele,
    LINALG::Matrix<3,1>& xi,
    std::vector<LINALG::Matrix<3,1> >& velValues,
    std::vector<double>& presValues
)
{
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    interpolation<numnode,DRT::Element::hex8>(ele,xi,velValues,presValues);
  }
  break;
  case DRT::Element::hex20:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
    interpolation<numnode,DRT::Element::hex20>(ele,xi,velValues,presValues);
  }
  break;
  default:
    dserror("xfem assembly type not yet implemented in time integration");
  }
}



/*------------------------------------------------------------------------------------------------*
 * perform the interpolation                                                     winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
void XFEM::ExtrapolationOld::interpolation(
    DRT::Element* ele,
    LINALG::Matrix<3,1>& xi,
    std::vector<LINALG::Matrix<3,1> >& velValues,
    std::vector<double>& presValues
)
{
  const int nsd = 3;

  // node velocities of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<nsd,2*numnode> > nodeveldata(oldVectors_.size(),LINALG::Matrix<nsd,2*numnode>(true));
  // node pressures of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<1,2*numnode> > nodepresdata(oldVectors_.size(),LINALG::Matrix<1,2*numnode>(true));

  // required enriched shape functions
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true);

  // dummies for function call
  LINALG::Matrix<nsd,nsd> dummy_jacobi(true);
  LINALG::Matrix<numnode,1> dummy_shpFcn(true);
  LINALG::Matrix<nsd,2*numnode> dummy_enrShapeXYVelDeriv1(true);
  LINALG::Matrix<nsd,2*numnode> dummy_enrShapeXYPresDeriv1(true);

  // evaluate data for the given point
  pointdataXFEM<numnode,DISTYPE>(
      ele,
      xi,
      dummy_jacobi,
      dummy_shpFcn,
      enrShapeFcnVel,
      enrShapeFcnPres,
      dummy_enrShapeXYVelDeriv1,
      dummy_enrShapeXYPresDeriv1,
      false
  );

  const int* elenodeids = ele->NodeIds();

  int dofcounterVelx = 0;
  int dofcounterVely = 0;
  int dofcounterVelz = 0;
  int dofcounterPres = 0;

  for (int nodeid=0;nodeid<ele->NumNode();nodeid++) // loop over element nodes
  {
    // get nodal velocities and pressures with help of the field set of node
    const std::set<XFEM::FieldEnr>& fieldEnrSet(olddofman_->getNodeDofSet(elenodeids[nodeid]));
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
        fieldenr != fieldEnrSet.end();++fieldenr)
    {
      const DofKey olddofkey(elenodeids[nodeid], *fieldenr);
      const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;
      switch (fieldenr->getEnrichment().Type())
      {
      case XFEM::Enrichment::typeStandard :
      case XFEM::Enrichment::typeJump :
      case XFEM::Enrichment::typeVoid :
      case XFEM::Enrichment::typeKink :
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodeveldata[index](0,dofcounterVelx) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterVelx++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodeveldata[index](1,dofcounterVely) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterVely++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodeveldata[index](2,dofcounterVelz) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterVelz++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodepresdata[index](0,dofcounterPres) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterPres++;
        }
        else
        {
          std::cout << XFEM::PHYSICS::physVarToString(fieldenr->getField()) << std::endl;
          dserror("not implemented physical field!");
        }
        break;
      }
      case XFEM::Enrichment::typeUndefined :
      default :
      {
        std::cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << std::endl;
        dserror("unknown enrichment type");
        break;
      }
      } // end switch enrichment
    } // end loop over fieldenr
  } // end loop over element nodes
  // compute shape functions

  LINALG::Matrix<nsd,1> vel(true);
  LINALG::Matrix<1,1> pres(true);

  for (size_t index=0;index<oldVectors_.size();index++)
  {
    vel.Multiply(nodeveldata[index],enrShapeFcnVel); // v = theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n
    velValues[index]=vel;

    pres.Multiply(nodepresdata[index],enrShapeFcnPres); // p = p_n + dt*(theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n)
    presValues[index] = pres(0);
  } // loop over vectors to be set
}



/*------------------------------------------------------------------------------------------------*
 * Extrapol constructor                                                     winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
XFEM::ExtrapolationNew::ExtrapolationNew(
    XFEM::TIMEINT& timeInt,
    INPAR::COMBUST::XFEMTimeIntegration timeIntType,
    const Teuchos::RCP<Epetra_Vector> veln,
    const double& dt,
    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
    bool initialize
) : STD(timeInt,timeIntType,veln,dt,flamefront,initialize)
{
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * call the computation based on an Extrapol                                winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationNew::compute(
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsn,
    std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectorsnp
)
{
  if (FGIType_==FRSNot1_)
    return;

  handleVectors(newRowVectorsn,newRowVectorsnp);

  resetState(TimeIntData::basicStd_,TimeIntData::extrapolateStd_);

  exportDataToNodeProc(); // export data of failed nodes

  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    ExtrapolationMain(&*data);
  }

#ifdef PARALLEL
  exportFinalData();
#endif
  setFinalData();

  // fill vectors with the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_!=TimeIntData::doneStd_)
      dserror("All data should be set here, having status 'done'. Thus something is wrong!");
  }
} // end compute



/*------------------------------------------------------------------------------------------------*
 * call the computation based on an Extrapol                                winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationNew::compute(std::vector<Teuchos::RCP<Epetra_Vector> > newRowVectors)
{
  exportDataToNodeProc(); // export data of failed nodes

  if (oldVectors_.size() != newRowVectors.size())
  {
    std::cout << "sizes are " << oldVectors_.size() << " and " << newRowVectors.size() << std::endl;
    dserror("Number of state-vectors at new and old discretization are different!");
  }

  newVectors_ = newRowVectors;

  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    ExtrapolationMain(&*data);
  }

#ifdef PARALLEL
  exportFinalData();
#endif
  setFinalData();
}// end compute



/*------------------------------------------------------------------------------------------------*
 * Extrapol of values to data-requiring nodes:                              winklmaier 10/11 *
 * a straight line through the data-requiring node and the nearest node on the correct interface  *
 * side is set up, two appropriate points on this line are used for Extrapol                 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationNew::ExtrapolationMain(
    TimeIntData* data
)
{
  if (data->state_!=TimeIntData::extrapolateStd_)
    dserror("Extrapol must be called with according state");

  const int nsd = 3;

  const DRT::Element* ele = NULL; // element of midpoint

  LINALG::Matrix<nsd,1> startpoint(true); // global coordinates of startpoint
  LINALG::Matrix<nsd,1> xistartpoint(true); // local coordinates of startpoint

  LINALG::Matrix<nsd,1> midpoint(true); // global coordinates of midpoint
  LINALG::Matrix<nsd,1> ximidpoint(true); // local coordinates of midpoint

  LINALG::Matrix<nsd,1> endpoint(data->node_.X());

  LINALG::Matrix<3,1> dummynormal(true);

//    std::cout << "searching node = endpoint = " << endpoint << std::endl;
  // identify final element and local coordinates of startpoint and midpoint for Extrapol
  Cases extrapolcase = EvalPoints(data,ele,midpoint,ximidpoint,startpoint,xistartpoint);

  if (extrapolcase==extrapol)
  {
    //    std::cout << std::endl << "startpoint is " << startpoint;
    //    std::cout << "midpoint is " << midpoint;

    // compute the constants for the Extrapol:
    // value = c1*valuestartpoint + c2*valuemidpoint
    LINALG::Matrix<nsd,1> dist1;
    LINALG::Matrix<nsd,1> dist2;

    dist1.Update(1.0,midpoint,-1.0,startpoint);
    dist2.Update(1.0,endpoint,-1.0,midpoint);
    //  std::cout << "distance from startpoint to midpoint is " << dist1;
    //  std::cout << "distance from midpoint to endpoint is " << dist2;

    double c1 = - dist2.Norm2()/dist1.Norm2(); // 1 + dist2/dist1
    double c2 = 1.0 + dist2.Norm2()/dist1.Norm2(); // dist2/dist1

    // get the velocities and the pressures at the start- and midpoint
    std::vector<LINALG::Matrix<nsd,1> > velstartpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));
    std::vector<LINALG::Matrix<nsd,1> > velmidpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));

    std::vector<double> presstartpoint(oldVectors_.size(),0.0);
    std::vector<double> presmidpoint(oldVectors_.size(),0.0);

    int side = interfaceSide(data->phiValue_);
    callInterpolation(ele,xistartpoint,velstartpoint,presstartpoint,side);
    callInterpolation(ele,ximidpoint,velmidpoint,presmidpoint,side);

    //  std::cout << "pres at startpoint is " << presstartpoint[0];
    //  std::cout << "pres at midpoint is " << presmidpoint[0];

    // compute the final velocities and pressure due to the Extrapol
    std::vector<LINALG::Matrix<nsd,1> > velendpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));
    std::vector<double> presendpoint(oldVectors_.size(),0.0);

    for (size_t index=0;index<oldVectors_.size();index++)
    {
      velendpoint[index].Update(c1,velstartpoint[index],c2,velmidpoint[index]);
      presendpoint[index] = c1*presstartpoint[index] + c2*presmidpoint[index];

      //  if (index == 0)
      //    std::cout << "final pressure is " << presendpoint[index] <<
      //        " and final velocity is " << velendpoint[index];

    } // loop over vectors to be set

    data->startOwner_ = std::vector<int>(1,myrank_);
    data->velValues_ = velendpoint;
    data->presValues_ = presendpoint;
    data->state_ = TimeIntData::doneStd_;
  }
  else if (extrapolcase==project)
  {
    // get the velocities and the pressures at the start- and midpoint
    std::vector<LINALG::Matrix<nsd,1> > velmidpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));
    std::vector<double> presmidpoint(oldVectors_.size(),0.0);

    int side = interfaceSide(data->phiValue_);
    callInterpolation(ele,ximidpoint,velmidpoint,presmidpoint,side);

    data->startOwner_ = std::vector<int>(1,myrank_);
    data->velValues_ = velmidpoint;
    data->presValues_ = presmidpoint;
    data->state_ = TimeIntData::doneStd_;
  }
  else if (extrapolcase==failed)
    data->state_ = TimeIntData::failedSL_;
  else
    dserror("extrapolation failed with unknown state");
}



/*------------------------------------------------------------------------------------------------*
 * perform a bisection on the line startpoint-endpoint                           winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
XFEM::ExtrapolationNew::Cases XFEM::ExtrapolationNew::EvalPoints(
    TimeIntData* data,
    const DRT::Element*& ele,
    LINALG::Matrix<3,1>& midpoint,
    LINALG::Matrix<3,1>& ximidpoint,
    LINALG::Matrix<3,1>& startpoint,
    LINALG::Matrix<3,1>& xistartpoint
)
{
  Cases extrapolcase = undefined;

  LINALG::Matrix<3,1> endpoint(data->node_.X());

  std::vector<const DRT::Element*> eles;
  addPBCelements(&data->node_, eles);
  const int numele = eles.size();

  double mindist = 0.0;
  double dist = 0.0;
  LINALG::Matrix<3,1> dummynormal(true);
  LINALG::Matrix<3,1> point(true);

  bool found = false;

  for (int i=0;i<numele;i++)
  {
    if (intersectionStatus(eles[i])!=XFEM::TIMEINT::uncut_)
    {
      bool worked = SignedDistance(&data->node_,eles[i]->Id(),dist,dummynormal,point,true);

      if ((found==false) or (dist<mindist))
      {
        if (worked)
        {
          found = true;
          mindist = dist;
          midpoint = point;
          ele = eles[i];
        }
      }
    }
  }

  if (!found)
  {
    extrapolcase = failed; // not working
  }
  else
  {
    bool pointInDomain = false;
    callXToXiCoords(ele,midpoint,ximidpoint,pointInDomain);
    if (!pointInDomain)
    {
      std::cout.precision(16);
      std::cout << "midpoint is " << midpoint;
      std::cout << "ele is " << *ele << std::endl;
      const DRT::Node*const* nodes = ele->Nodes();
      for (int i=0;i<ele->NumNode();i++)
        std::cout << " ele node is " << *nodes[i] << std::endl;

      std::cout << "local coords are " << ximidpoint;
      dserror("point shall be in element");
    }


    startpoint.Update(2.0,midpoint,-1.0,endpoint);
    callXToXiCoords(ele,startpoint,xistartpoint,pointInDomain);
    if (!pointInDomain)
    {
      startpoint.Update(1.5,midpoint,-1.0,endpoint);
      callXToXiCoords(ele,startpoint,xistartpoint,pointInDomain);
    }

    if (!pointInDomain)
      extrapolcase = failed;
    else
    {
      const int* elenodeids = ele->NodeIds();

      Epetra_SerialDenseVector nodephi(ele->NumNode()); // nodal phivalues
      for (int nodeid=0;nodeid<ele->NumNode();nodeid++) // loop over element nodes
        nodephi(nodeid) = (*phin_)[discret_->gNode(elenodeids[nodeid])->LID()];

      Epetra_SerialDenseVector funct(ele->NumNode()); // nodal phivalues
      DRT::UTILS::shape_function_3D(funct, xistartpoint(0),xistartpoint(1),xistartpoint(2),ele->Shape()); // evaluate shape functions at xi

      double phivalue = funct.Dot(nodephi);

      if (interfaceSideCompare(phivalue,data->phiValue_))
        extrapolcase = extrapol;
      else
        extrapolcase = project;
    }
  }

  return extrapolcase;
} // end EvalPoints



/*------------------------------------------------------------------------------------------------*
 * call the interpolation                                                        winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationNew::callInterpolation(
    const DRT::Element* ele,
    LINALG::Matrix<3,1>& xi,
    std::vector<LINALG::Matrix<3,1> >& velValues,
    std::vector<double>& presValues,
    int side
)
{
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    interpolation<numnode,DRT::Element::hex8>(ele,xi,velValues,presValues,side);
  }
  break;
  case DRT::Element::hex20:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
    interpolation<numnode,DRT::Element::hex20>(ele,xi,velValues,presValues,side);
  }
  break;
  default:
    dserror("xfem assembly type not yet implemented in time integration");
  }
}



/*------------------------------------------------------------------------------------------------*
 * perform the interpolation                                                     winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
void XFEM::ExtrapolationNew::interpolation(
    const DRT::Element* ele,
    LINALG::Matrix<3,1>& xi,
    std::vector<LINALG::Matrix<3,1> >& velValues,
    std::vector<double>& presValues,
    int side
)
{
  const int nsd = 3;

  // node velocities of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<nsd,2*numnode> > nodeveldata(oldVectors_.size(),LINALG::Matrix<nsd,2*numnode>(true));
  // node pressures of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<1,2*numnode> > nodepresdata(oldVectors_.size(),LINALG::Matrix<1,2*numnode>(true));

  // required enriched shape functions
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true);

  // dummies for function call
  LINALG::Matrix<nsd,nsd> dummy_jacobi(true);
  LINALG::Matrix<numnode,1> dummy_shpFcn(true);
  LINALG::Matrix<nsd,2*numnode> dummy_enrShapeXYVelDeriv1(true);
  LINALG::Matrix<nsd,2*numnode> dummy_enrShapeXYPresDeriv1(true);

  // evaluate data for the given point
  pointdataXFEM<numnode,DISTYPE>(
      ele,
      xi,
      dummy_jacobi,
      dummy_shpFcn,
      enrShapeFcnVel,
      enrShapeFcnPres,
      dummy_enrShapeXYVelDeriv1,
      dummy_enrShapeXYPresDeriv1,
      false,
      side
  );

  const int* elenodeids = ele->NodeIds();

  int dofcounterVelx = 0;
  int dofcounterVely = 0;
  int dofcounterVelz = 0;
  int dofcounterPres = 0;

  for (int nodeid=0;nodeid<ele->NumNode();nodeid++) // loop over element nodes
  {
    // get nodal velocities and pressures with help of the field set of node
    const std::set<XFEM::FieldEnr>& fieldEnrSet(olddofman_->getNodeDofSet(elenodeids[nodeid]));
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
        fieldenr != fieldEnrSet.end();++fieldenr)
    {
      const DofKey olddofkey(elenodeids[nodeid], *fieldenr);
      const int olddofpos = oldNodalDofColDistrib_.find(olddofkey)->second;
      switch (fieldenr->getEnrichment().Type())
      {
      case XFEM::Enrichment::typeStandard :
      case XFEM::Enrichment::typeJump :
      case XFEM::Enrichment::typeVoid :
      case XFEM::Enrichment::typeKink :
      {
        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodeveldata[index](0,dofcounterVelx) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterVelx++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodeveldata[index](1,dofcounterVely) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterVely++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodeveldata[index](2,dofcounterVelz) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterVelz++;
        }
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodepresdata[index](0,dofcounterPres) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterPres++;
        }
        else
        {
          std::cout << XFEM::PHYSICS::physVarToString(fieldenr->getField()) << std::endl;
          dserror("not implemented physical field!");
        }
        break;
      }
      case XFEM::Enrichment::typeUndefined :
      default :
      {
        std::cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << std::endl;
        dserror("unknown enrichment type");
        break;
      }
      } // end switch enrichment
    } // end loop over fieldenr
  } // end loop over element nodes
  // compute shape functions

  LINALG::Matrix<nsd,1> vel(true);
  LINALG::Matrix<1,1> pres(true);

  for (size_t index=0;index<oldVectors_.size();index++)
  {
    vel.Multiply(nodeveldata[index],enrShapeFcnVel); // v = theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n
    velValues[index]=vel;

    pres.Multiply(nodepresdata[index],enrShapeFcnPres); // p = p_n + dt*(theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n)
    presValues[index] = pres(0);
  } // loop over vectors to be set
}



#ifdef PARALLEL
/*------------------------------------------------------------------------------------------------*
 * export alternative algo data to neighbour proc                                winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::ExtrapolationNew::exportDataToNodeProc()
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  std::vector<std::vector<TimeIntData> > dataVec(numproc_);

  // fill vectors with the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::extrapolateStd_)
      dataVec[data->node_.Owner()].push_back(*data);
  }

  clearState(TimeIntData::extrapolateStd_);
  timeIntData_->insert(timeIntData_->end(),
      dataVec[myrank_].begin(),
      dataVec[myrank_].end());

  dataVec[myrank_].clear(); // clear the set data from the vector

  // send data to the processor where the point lies (1. nearest higher neighbour 2. 2nd nearest higher neighbour...)
  for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
  {
    // Initialization of sending
    DRT::PackBuffer dataSend; // vector including all data that has to be send to dest proc

    // Initialization
    int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
    if (source<0)
      source+=numproc_;
    else if (source>=numproc_)
      source -=numproc_;

    // pack data to be sent
    for (std::vector<TimeIntData>::iterator data=dataVec[dest].begin();
        data!=dataVec[dest].end(); data++)
    {
      if (data->state_==TimeIntData::extrapolateStd_)
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
        DRT::ParObject::AddtoPack(dataSend,data->startGid_);
        DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
      }
    }

    dataSend.StartPacking();

    for (std::vector<TimeIntData>::iterator data=dataVec[dest].begin();
        data!=dataVec[dest].end(); data++)
    {
      if (data->state_==TimeIntData::extrapolateStd_)
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->phiValue_);
        DRT::ParObject::AddtoPack(dataSend,data->startGid_);
        DRT::ParObject::AddtoPack(dataSend,data->startOwner_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
      }
    }

    // clear the no more needed data
    dataVec[dest].clear();

    std::vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    std::vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      double coords[nsd] = {0.0};
      DRT::Node node(0,(double*)coords,0);
      LINALG::Matrix<nsd,1> vel;
      std::vector<LINALG::Matrix<nsd,nsd> > velDeriv;
      std::vector<LINALG::Matrix<1,nsd> > presDeriv;
      LINALG::Matrix<nsd,1> startpoint;
      double phiValue;
      std::vector<int> startGid;
      std::vector<int> startOwner;
      int newtype;

      unpackNode(posinData,dataRecv,node);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,phiValue);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startGid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startOwner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,newtype);

      timeIntData_->push_back(TimeIntData(
          node,
          vel,
          velDeriv,
          presDeriv,
          startpoint,
          phiValue,
          startGid,
          startOwner,
          (TimeIntData::type)newtype)); // startOwner is current proc
    } // end loop over number of nodes to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end loop over processors

  resetState(TimeIntData::failedSL_,TimeIntData::extrapolateStd_);
} // end exportDataToStartpointProc
#endif


