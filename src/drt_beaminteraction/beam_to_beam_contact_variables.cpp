/*----------------------------------------------------------------------------*/
/*!
\file beam_to_beam_contact_variables.cpp

\brief one beam contact segment living on an element pair

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "../drt_beaminteraction/beam_to_beam_contact_variables.H"

#include "beam3contact_utils.H"

#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_structure/strtimint_impl.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"

#include "Teuchos_TimeMonitor.hpp"

#include "../drt_beaminteraction/beam_to_beam_contact_pair.H"
#include "beam3contact_defines.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                     meier 01/14|
 *----------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
BEAMINTERACTION::BeamToBeamContactVariables<numnodes, numnodalvalues>::BeamToBeamContactVariables(
    std::pair<TYPE,TYPE>& closestpoint,
    std::pair<int,int>& segids,
    std::pair<int,int>& intids,
    const double& pp,
    TYPE jacobi):
closestpoint_(closestpoint),
segids_(segids),
intids_(intids),
jacobi_(jacobi),
gap_(0.0),
normal_(LINALG::TMatrix<TYPE,3,1>(true)),
pp_(pp),
ppfac_(0.0),
dppfac_(0.0),
fp_(0.0),
dfp_(0.0),
energy_(0.0),
integratedenergy_(0.0),
angle_(0.0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  end: constructor
 *----------------------------------------------------------------------*/

//Possible template cases: this is necessary for the compiler
template class BEAMINTERACTION::BeamToBeamContactVariables<2,1>;
template class BEAMINTERACTION::BeamToBeamContactVariables<3,1>;
template class BEAMINTERACTION::BeamToBeamContactVariables<4,1>;
template class BEAMINTERACTION::BeamToBeamContactVariables<5,1>;
template class BEAMINTERACTION::BeamToBeamContactVariables<2,2>;
