/*!------------------------------------------------------------------------------------------------*
\file startvalues.cpp

\brief provides the Extrapol class

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "timeInt_std_extrapol.H"


/*------------------------------------------------------------------------------------------------*
 * Extrapol constructor                                                     winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
XFEM::Extrapol::Extrapol(
    XFEM::TIMEINT& timeInt,
    INPAR::COMBUST::XFEMTimeIntegration timeIntType,
    const RCP<Epetra_Vector> veln,
    const double& dt,
    const RCP<COMBUST::FlameFront> flamefront,
    bool initialize
) : STD(timeInt,timeIntType,veln,dt,flamefront,initialize)
{
  return;
} // end constructor



/*------------------------------------------------------------------------------------------------*
 * call the computation based on an Extrapol                                winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Extrapol::compute(
    vector<RCP<Epetra_Vector> > newRowVectorsn,
    vector<RCP<Epetra_Vector> > newRowVectorsnp
)
{
  if (FGIType_==FRSNot1_)
    return;

  handleVectors(newRowVectorsn,newRowVectorsnp);

  resetState(TimeIntData::basicStd_,TimeIntData::extrapolateStd_);

  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    ExtrapolMain(&*data);
  }

#ifdef PARALLEL
  exportFinalData();
#endif
  setFinalData();

  // fill vectors with the data
  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_!=TimeIntData::doneStd_)
      dserror("All data should be set here, having status 'done'. Thus something is wrong!");
  }

} // end compute



/*------------------------------------------------------------------------------------------------*
 * call the computation based on an Extrapol                                winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Extrapol::compute(vector<RCP<Epetra_Vector> > newRowVectors)
{
  if (oldVectors_.size() != newRowVectors.size())
  {
    cout << "sizes are " << oldVectors_.size() << " and " << newRowVectors.size() << endl;
    dserror("Number of state-vectors at new and old discretization are different!");
  }

  newVectors_ = newRowVectors;

  for (vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    ExtrapolMain(&*data);
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
void XFEM::Extrapol::ExtrapolMain(
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

//    cout << "searching node = endpoint = " << endpoint << endl;
  // identify final element and local coordinates of startpoint and midpoint for Extrapol
  bool failed = false;
  EvalPoints(data,ele,midpoint,ximidpoint,startpoint,xistartpoint,failed);

  if (!failed)
  {
    //    cout << endl << "startpoint is " << startpoint;
    //    cout << "midpoint is " << midpoint;

    // compute the constants for the Extrapol:
    // value = c1*valuestartpoint + c2*valuemidpoint
    LINALG::Matrix<nsd,1> dist1;
    LINALG::Matrix<nsd,1> dist2;

    dist1.Update(1.0,midpoint,-1.0,startpoint);
    dist2.Update(1.0,endpoint,-1.0,midpoint);
    //  cout << "distance from startpoint to midpoint is " << dist1;
    //  cout << "distance from midpoint to endpoint is " << dist2;

    double c1 = - dist2.Norm2()/dist1.Norm2(); // 1 + dist2/dist1
    double c2 = 1.0 + dist2.Norm2()/dist1.Norm2(); // dist2/dist1

    // get the velocities and the pressures at the start- and midpoint
    vector<LINALG::Matrix<nsd,1> > velstartpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));
    vector<LINALG::Matrix<nsd,1> > velmidpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));

    vector<double> presstartpoint(oldVectors_.size(),0.0);
    vector<double> presmidpoint(oldVectors_.size(),0.0);

    int side = interfaceSide(data->phiValue_);
    callInterpolation(ele,xistartpoint,velstartpoint,presstartpoint,side);
    callInterpolation(ele,ximidpoint,velmidpoint,presmidpoint,side);

    //  cout << "pres at startpoint is " << presstartpoint[0];
    //  cout << "pres at midpoint is " << presmidpoint[0];

    // compute the final velocities and pressure due to the Extrapol
    vector<LINALG::Matrix<nsd,1> > velendpoint(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));
    vector<double> presendpoint(oldVectors_.size(),0.0);

    for (size_t index=0;index<oldVectors_.size();index++)
    {
      velendpoint[index].Update(c1,velstartpoint[index],c2,velmidpoint[index]);
      presendpoint[index] = c1*presstartpoint[index] + c2*presmidpoint[index];

      //  if (index == 0)
      //    cout << "final pressure is " << presendpoint[index] <<
      //        " and final velocity is " << velendpoint[index];

    } // loop over vectors to be set

    data->startOwner_ = vector<int>(1,myrank_);
    data->velValues_ = velendpoint;
    data->presValues_ = presendpoint;
    data->state_ = TimeIntData::doneStd_;
  }
}



/*------------------------------------------------------------------------------------------------*
 * perform a bisection on the line startpoint-endpoint                           winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Extrapol::EvalPoints(
    TimeIntData* data,
    const DRT::Element*& ele,
    LINALG::Matrix<3,1>& midpoint,
    LINALG::Matrix<3,1>& ximidpoint,
    LINALG::Matrix<3,1>& startpoint,
    LINALG::Matrix<3,1>& xistartpoint,
    bool& failed
)
{
  failed = true;

  LINALG::Matrix<3,1> endpoint(data->node_.X());

  vector<const DRT::Element*> eles;
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

  if (found)
  {
    bool pointInDomain = false;
    callXToXiCoords(ele,midpoint,ximidpoint,pointInDomain);
    if (!pointInDomain)
    {
      cout.precision(16);
      cout << "midpoint is " << midpoint;
      cout << "ele is " << *ele << endl;
      const DRT::Node*const* nodes = ele->Nodes();
      for (int i=0;i<ele->NumNode();i++)
        cout << " ele node is " << *nodes[i] << endl;

      cout << "local coords are " << ximidpoint;
      dserror("point shall be in element");
    }


    startpoint.Update(2.0,midpoint,-1.0,endpoint);
    callXToXiCoords(ele,startpoint,xistartpoint,pointInDomain);
    if (!pointInDomain)
    {
      startpoint.Update(1.5,midpoint,-1.0,endpoint);
      callXToXiCoords(ele,startpoint,xistartpoint,pointInDomain);
    }

    if (pointInDomain)
      failed = false;
    //    cout << "searching node = endpoint = " << endpoint << endl;
    // identify final element and local coordinates of startpoint and midpoint for Extrapol
  }
} // end EvalPoints



/*------------------------------------------------------------------------------------------------*
 * call the interpolation                                                        winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::Extrapol::callInterpolation(
    const DRT::Element* ele,
    LINALG::Matrix<3,1>& xi,
    vector<LINALG::Matrix<3,1> >& velValues,
    vector<double>& presValues,
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
void XFEM::Extrapol::interpolation(
    const DRT::Element* ele,
    LINALG::Matrix<3,1>& xi,
    vector<LINALG::Matrix<3,1> >& velValues,
    vector<double>& presValues,
    int side
)
{
  const int nsd = 3;

  // node velocities of the element nodes for the data that should be changed
  vector<LINALG::Matrix<nsd,2*numnode> > nodeveldata(oldVectors_.size(),LINALG::Matrix<nsd,2*numnode>(true));
  // node pressures of the element nodes for the data that should be changed
  vector<LINALG::Matrix<1,2*numnode> > nodepresdata(oldVectors_.size(),LINALG::Matrix<1,2*numnode>(true));
#ifdef COMBUST_NORMAL_ENRICHMENT
  vector<LINALG::Matrix<1,numnode> > nodevelenrdata(oldVectors_.size(),LINALG::Matrix<1,numnode>(true));
  LINALG::Matrix<1,numnode> nodevelenr(true);
#endif

  // required enriched shape functions
  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
  LINALG::Matrix<2*numnode,1> enrShapeFcnPres(true);

  // dummies for function call
  LINALG::Matrix<nsd,nsd> dummy_jacobi(true);
  LINALG::Matrix<numnode,1> dummy_shpFcn(true);
  LINALG::Matrix<nsd,2*numnode> dummy_enrShapeXYVelDeriv1(true);
  LINALG::Matrix<nsd,2*numnode> dummy_enrShapeXYPresDeriv1(true);

  // evaluate data for the given point
#ifdef COMBUST_NORMAL_ENRICHMENT
#ifdef COLLAPSE_FLAME_NORMAL
  LINALG::Matrix<nsd,1> normal(node->X());
  normal(2) = 0.0;
  normal.Scale(1.0/normal.Norm2());
#endif
  ApproxFuncNormalVector<2,2*numnode> shp;
  pointdataXFEMNormal<numnode,DISTYPE>(
      ele,
#ifdef COLLAPSE_FLAME_NORMAL
      normal,
#endif
      xi,
      xji,
      shapeFcn,
      enrShapeFcnPres,
      enrShapeXYPresDeriv1,
      shp,
      false
  );
#else
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
#endif

  const int* elenodeids = ele->NodeIds();

  int dofcounterVelx = 0;
  int dofcounterVely = 0;
  int dofcounterVelz = 0;
  int dofcounterPres = 0;
#ifdef COMBUST_NORMAL_ENRICHMENT
  int dofcounterVeln = 0;
#endif

  for (int nodeid=0;nodeid<ele->NumNode();nodeid++) // loop over element nodes
  {
    // get nodal velocities and pressures with help of the field set of node
    const std::set<XFEM::FieldEnr>& fieldEnrSet(olddofman_->getNodeDofSet(elenodeids[nodeid]));
    for (set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
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
#ifdef COMBUST_NORMAL_ENRICHMENT
        else if (fieldenr->getField() == XFEM::PHYSICS::Veln)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodevelenrdata[index](0,dofcounterVeln) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterVeln++;
        }
#endif
        else if (fieldenr->getField() == XFEM::PHYSICS::Pres)
        {
          for (size_t index=0;index<oldVectors_.size();index++)
            nodepresdata[index](0,dofcounterPres) = (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
          dofcounterPres++;
        }
        else
        {
          cout << XFEM::PHYSICS::physVarToString(fieldenr->getField()) << endl;
          dserror("not implemented physical field!");
        }
        break;
      }
      case XFEM::Enrichment::typeUndefined :
      default :
      {
        cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << endl;
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
