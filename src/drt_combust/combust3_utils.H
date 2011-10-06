/*!
\file combust3_utils.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
#ifdef CCADISCRET
#ifndef COMBUST3_UTILS_H
#define COMBUST3_UTILS_H

#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_condition_utils.H"

namespace COMBUST
{

  template <size_t numparam>
  struct EnrViscs2
  {
    LINALG::Matrix<numparam,1> xx;
    LINALG::Matrix<numparam,1> xy;
    LINALG::Matrix<numparam,1> xz;
    LINALG::Matrix<numparam,1> yx;
    LINALG::Matrix<numparam,1> yy;
    LINALG::Matrix<numparam,1> yz;
    LINALG::Matrix<numparam,1> zx;
    LINALG::Matrix<numparam,1> zy;
    LINALG::Matrix<numparam,1> zz;
  };

  /*!
   * \brief calculate the body force from the nodes condition
   *
   *  get the body force in the nodes of the element
   *  the Neumann condition associated with the nodes is stored in the
   *  array edeadng only if all nodes have a VolumeNeumann condition
   *
   *  \author rasthofer 08/09
   */
  template <DRT::Element::DiscretizationType DISTYPE,
            class M>
  void BodyForce(
      const DRT::Element* ele,
      const double time,
      const string neumtype,
      M& edead)
  {
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    vector<DRT::Condition*> myneumcond;
    int myindex = -1;

    // check whether all nodes have a unique VolumeNeumann condition
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);

    for (size_t i = 0; i < myneumcond.size(); ++i)
    {
      const string* condtype = myneumcond[i]->Get<string>("type");
      if (*condtype == neumtype or (neumtype == "DeadOrLive" and (*condtype == "neum_dead" or *condtype == "neum_live")))
      {
        myindex = (int)i;
        break;
      }
    }

    if (myindex >= 0)
    {
      // find out whether we will use a time curve
      const vector<int>* curve  = myneumcond[myindex]->Get<vector<int> >("curve");
      int curvenum = -1;

      if (curve) curvenum = (*curve)[0];

      // initialisation
      double curvefac    = 0.0;

      if (curvenum >= 0) // yes, we have a timecurve
      {
        // time factor for the intermediate step
        if(time >= 0.0)
        {
          curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
        }
        else
        {
          // do not compute an "alternative" curvefac here since a negative time value
          // indicates an error.
          dserror("Negative time value in body force calculation: time = %f",time);
          //curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
        }
      }
      else // we do not have a timecurve --- timefactors are constant equal 1
      {
        curvefac = 1.0;
      }

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[myindex]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[myindex]->Get<vector<double> >("val"  );
      const vector<int>*    functions = myneumcond[myindex]->Get<vector<int> >("funct");

      // factor given by spatial function
      double functionfac = 1.0;
      int functnum = -1;

      // set this condition to the edeadng array
      for(int isd=0;isd<3;isd++)
      {
        // get factor given by spatial function
        if (functions) functnum = (*functions)[isd];
        else functnum = -1;

        double num = (*onoff)[isd]*(*val)[isd]*curvefac;

        for (std::size_t jnode=0; jnode<numnode; jnode++)
        {
          if (functnum>0)
          {
            // evaluate function at the position of the current node
            functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,(ele->Nodes()[jnode])->X(),time,NULL);
          }
          else functionfac = 1.0;

          edead(isd,jnode) = num*functionfac;
        }
      }
    }
    else if(myneumcond.size()==2)
    {

    }
  }


  //! compute largest element diameter for Nitsche stabilization parameter
  template<DRT::Element::DiscretizationType DISTYPE, class M1>
  double getEleDiameter(const M1& xyze)
  {
    double elediam = 0.0;

    // number of nodes of this element
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // check all possible connections between nodes of an element
    // node 1 to 2
    //    :
    // node 1 to 8 = numnode
    // node 2 to 3
    //    :
    // node 2 to 8
    //    :
    //    :
    //    :
    // node 7 to 8
    for(size_t i_start=0; i_start< numnode-2; ++i_start)
    {
      for(size_t i_end= i_start+1; i_end < numnode-1; ++i_end)
      {
        LINALG::Matrix<3,1> direction;
        direction.Clear();
        direction(0) = xyze(0, i_start) - xyze(0, i_end);
        direction(1) = xyze(1, i_start) - xyze(1, i_end);
        direction(2) = xyze(2, i_start) - xyze(2, i_end);

        // update elediam
        if (direction.Norm2() > elediam) elediam=direction.Norm2();
      }
    }

    return elediam;
  }

} // namespace COMBUST

#endif
#endif  // #ifdef CCADISCRET
