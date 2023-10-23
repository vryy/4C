/*---------------------------------------------------------------------*/
/*! \file
\brief List of template combinations for the augmented contact integrator.

\level 2

*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_AUG_INTEGRATOR_INST_H
#define BACI_CONTACT_AUG_INTEGRATOR_INST_H

#include "baci_contact_aug_integrator_policy.H"

template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::line2,
    DRT::Element::DiscretizationType::line2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::DiscretizationType::line2,
        DRT::Element::DiscretizationType::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::line2,
    DRT::Element::DiscretizationType::line2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::DiscretizationType::line2,
        DRT::Element::DiscretizationType::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::line2,
    DRT::Element::DiscretizationType::nurbs2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::DiscretizationType::line2,
        DRT::Element::DiscretizationType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::line2,
    DRT::Element::DiscretizationType::nurbs2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::DiscretizationType::line2,
        DRT::Element::DiscretizationType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::line2,
    DRT::Element::DiscretizationType::nurbs3,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::DiscretizationType::line2,
        DRT::Element::DiscretizationType::nurbs3>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::line2,
    DRT::Element::DiscretizationType::nurbs3,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::DiscretizationType::line2,
        DRT::Element::DiscretizationType::nurbs3>>;

template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs2,
    DRT::Element::DiscretizationType::nurbs2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs2,
        DRT::Element::DiscretizationType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs2,
    DRT::Element::DiscretizationType::nurbs2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs2,
        DRT::Element::DiscretizationType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs2,
    DRT::Element::DiscretizationType::line2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs2,
        DRT::Element::DiscretizationType::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs2,
    DRT::Element::DiscretizationType::line2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs2,
        DRT::Element::DiscretizationType::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs2,
    DRT::Element::DiscretizationType::nurbs3,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs2,
        DRT::Element::DiscretizationType::nurbs3>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs2,
    DRT::Element::DiscretizationType::nurbs3,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs2,
        DRT::Element::DiscretizationType::nurbs3>>;

template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs3,
    DRT::Element::DiscretizationType::nurbs3,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs3,
        DRT::Element::DiscretizationType::nurbs3>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs3,
    DRT::Element::DiscretizationType::nurbs3,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs3,
        DRT::Element::DiscretizationType::nurbs3>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs3,
    DRT::Element::DiscretizationType::line2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs3,
        DRT::Element::DiscretizationType::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs3,
    DRT::Element::DiscretizationType::line2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs3,
        DRT::Element::DiscretizationType::line2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs3,
    DRT::Element::DiscretizationType::nurbs2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs3,
        DRT::Element::DiscretizationType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, DRT::Element::DiscretizationType::nurbs3,
    DRT::Element::DiscretizationType::nurbs2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, DRT::Element::DiscretizationType::nurbs3,
        DRT::Element::DiscretizationType::nurbs2>>;

template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::quad4,
    DRT::Element::DiscretizationType::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::quad4,
        DRT::Element::DiscretizationType::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::quad4,
    DRT::Element::DiscretizationType::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::quad4,
        DRT::Element::DiscretizationType::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::quad4,
    DRT::Element::DiscretizationType::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::quad4,
        DRT::Element::DiscretizationType::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::quad4,
    DRT::Element::DiscretizationType::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::quad4,
        DRT::Element::DiscretizationType::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::quad4,
    DRT::Element::DiscretizationType::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::quad4,
        DRT::Element::DiscretizationType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::quad4,
    DRT::Element::DiscretizationType::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::quad4,
        DRT::Element::DiscretizationType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::quad4,
    DRT::Element::DiscretizationType::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::quad4,
        DRT::Element::DiscretizationType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::quad4,
    DRT::Element::DiscretizationType::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::quad4,
        DRT::Element::DiscretizationType::nurbs9>>;

template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::tri3,
    DRT::Element::DiscretizationType::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::tri3,
        DRT::Element::DiscretizationType::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::tri3,
    DRT::Element::DiscretizationType::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::tri3,
        DRT::Element::DiscretizationType::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::tri3,
    DRT::Element::DiscretizationType::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::tri3,
        DRT::Element::DiscretizationType::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::tri3,
    DRT::Element::DiscretizationType::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::tri3,
        DRT::Element::DiscretizationType::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::tri3,
    DRT::Element::DiscretizationType::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::tri3,
        DRT::Element::DiscretizationType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::tri3,
    DRT::Element::DiscretizationType::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::tri3,
        DRT::Element::DiscretizationType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::tri3,
    DRT::Element::DiscretizationType::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::tri3,
        DRT::Element::DiscretizationType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::tri3,
    DRT::Element::DiscretizationType::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::tri3,
        DRT::Element::DiscretizationType::nurbs9>>;

template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs4,
    DRT::Element::DiscretizationType::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs4,
        DRT::Element::DiscretizationType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs4,
    DRT::Element::DiscretizationType::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs4,
        DRT::Element::DiscretizationType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs4,
    DRT::Element::DiscretizationType::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs4,
        DRT::Element::DiscretizationType::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs4,
    DRT::Element::DiscretizationType::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs4,
        DRT::Element::DiscretizationType::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs4,
    DRT::Element::DiscretizationType::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs4,
        DRT::Element::DiscretizationType::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs4,
    DRT::Element::DiscretizationType::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs4,
        DRT::Element::DiscretizationType::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs4,
    DRT::Element::DiscretizationType::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs4,
        DRT::Element::DiscretizationType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs4,
    DRT::Element::DiscretizationType::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs4,
        DRT::Element::DiscretizationType::nurbs9>>;

template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs9,
    DRT::Element::DiscretizationType::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs9,
        DRT::Element::DiscretizationType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs9,
    DRT::Element::DiscretizationType::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs9,
        DRT::Element::DiscretizationType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs9,
    DRT::Element::DiscretizationType::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs9,
        DRT::Element::DiscretizationType::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs9,
    DRT::Element::DiscretizationType::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs9,
        DRT::Element::DiscretizationType::tri3>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs9,
    DRT::Element::DiscretizationType::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs9,
        DRT::Element::DiscretizationType::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs9,
    DRT::Element::DiscretizationType::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs9,
        DRT::Element::DiscretizationType::quad4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs9,
    DRT::Element::DiscretizationType::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs9,
        DRT::Element::DiscretizationType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, DRT::Element::DiscretizationType::nurbs9,
    DRT::Element::DiscretizationType::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, DRT::Element::DiscretizationType::nurbs9,
        DRT::Element::DiscretizationType::nurbs4>>;

#endif