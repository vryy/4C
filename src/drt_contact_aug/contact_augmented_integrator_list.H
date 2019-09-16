/*---------------------------------------------------------------------*/
/*! \file
\brief List of template combinations for the augmented contact integrator.

\level 2

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------*/

#include "contact_aug_integrator_policy.H"

template class CONTACT::AUG::Integrator<2, DRT::Element::line2, DRT::Element::line2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::line2, DRT::Element::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::line2, DRT::Element::line2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::line2, DRT::Element::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::line2, DRT::Element::nurbs2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::line2, DRT::Element::nurbs2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::line2, DRT::Element::nurbs3,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs3>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::line2, DRT::Element::nurbs3,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs3>>;

template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs2, DRT::Element::nurbs2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs2, DRT::Element::nurbs2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs2, DRT::Element::line2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs2, DRT::Element::line2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs2, DRT::Element::nurbs3,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs3>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs2, DRT::Element::nurbs3,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs3>>;

template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs3, DRT::Element::nurbs3,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs3>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs3, DRT::Element::nurbs3,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs3>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs3, DRT::Element::line2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs3, DRT::Element::line2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs3, DRT::Element::nurbs2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::nurbs3, DRT::Element::nurbs2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs2>>;

template class CONTACT::AUG::Integrator<3, DRT::Element::quad4, DRT::Element::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::quad4, DRT::Element::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::quad4, DRT::Element::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::quad4, DRT::Element::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::quad4, DRT::Element::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::quad4, DRT::Element::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::quad4, DRT::Element::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::quad4, DRT::Element::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs9>>;

template class CONTACT::AUG::Integrator<3, DRT::Element::tri3, DRT::Element::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::tri3, DRT::Element::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::tri3, DRT::Element::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::tri3, DRT::Element::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::tri3, DRT::Element::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::tri3, DRT::Element::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::tri3, DRT::Element::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::tri3, DRT::Element::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs9>>;

template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs4, DRT::Element::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs4, DRT::Element::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs4, DRT::Element::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs4, DRT::Element::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs4, DRT::Element::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs4, DRT::Element::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs4, DRT::Element::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs4, DRT::Element::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs9>>;

template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs9, DRT::Element::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs9, DRT::Element::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs9, DRT::Element::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs9, DRT::Element::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs9, DRT::Element::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs9, DRT::Element::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs9, DRT::Element::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::nurbs9, DRT::Element::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs4>>;
