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
  template <DRT::Element::DiscretizationType DISTYPE>
  LINALG::SerialDenseMatrix BodyForce(
      const DRT::Element* ele,
      const double time)
  {
    const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    LINALG::SerialDenseMatrix edeadng(3,numnode);
    edeadng.Zero();
    vector<DRT::Condition*> myneumcond;

    // check whether all nodes have a unique VolumeNeumann condition
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);

    if (myneumcond.size()>1)
      dserror("more than one VolumeNeumann cond on one node");

    if (myneumcond.size()==1)
    {
      // find out whether we will use a time curve
      const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
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
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );
      const vector<int>*    functions = myneumcond[0]->Get<vector<int> >("funct");

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

          edeadng(isd,jnode) = num*functionfac;
        }
      }
    }

    return edeadng;
  }

} // namespace COMBUST

#endif
#endif  // #ifdef CCADISCRET
