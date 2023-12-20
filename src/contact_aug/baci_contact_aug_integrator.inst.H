/*---------------------------------------------------------------------*/
/*! \file
\brief List of template combinations for the augmented contact integrator.

\level 2

*/
/*---------------------------------------------------------------------*/
#ifndef BACI_CONTACT_AUG_INTEGRATOR_INST_H
#define BACI_CONTACT_AUG_INTEGRATOR_INST_H

#include "baci_config.H"

#include "baci_contact_aug_integrator_policy.H"

BACI_NAMESPACE_OPEN

template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::line2, CORE::FE::CellType::line2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::line2, CORE::FE::CellType::line2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::line2, CORE::FE::CellType::line2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::line2,
        CORE::FE::CellType::line2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::line2, CORE::FE::CellType::nurbs2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::line2, CORE::FE::CellType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::line2, CORE::FE::CellType::nurbs2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::line2,
        CORE::FE::CellType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::line2, CORE::FE::CellType::nurbs3,
    CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::line2, CORE::FE::CellType::nurbs3>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::line2, CORE::FE::CellType::nurbs3,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::line2,
        CORE::FE::CellType::nurbs3>>;

template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs2, CORE::FE::CellType::nurbs2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
        CORE::FE::CellType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs2, CORE::FE::CellType::nurbs2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
        CORE::FE::CellType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs2, CORE::FE::CellType::line2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs2, CORE::FE::CellType::line2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs2, CORE::FE::CellType::line2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
        CORE::FE::CellType::line2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs2, CORE::FE::CellType::nurbs3,
    CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
        CORE::FE::CellType::nurbs3>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs2, CORE::FE::CellType::nurbs3,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs2,
        CORE::FE::CellType::nurbs3>>;

template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs3, CORE::FE::CellType::nurbs3,
    CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
        CORE::FE::CellType::nurbs3>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs3, CORE::FE::CellType::nurbs3,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
        CORE::FE::CellType::nurbs3>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs3, CORE::FE::CellType::line2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs3, CORE::FE::CellType::line2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs3, CORE::FE::CellType::line2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
        CORE::FE::CellType::line2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs3, CORE::FE::CellType::nurbs2,
    CONTACT::AUG::DebugCompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
        CORE::FE::CellType::nurbs2>>;
template class CONTACT::AUG::Integrator<2, CORE::FE::CellType::nurbs3, CORE::FE::CellType::nurbs2,
    CONTACT::AUG::DebugIncompleteIntPolicy<2, CORE::FE::CellType::nurbs3,
        CORE::FE::CellType::nurbs2>>;

template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::quad4, CORE::FE::CellType::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::quad4, CORE::FE::CellType::quad4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::quad4, CORE::FE::CellType::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::quad4,
        CORE::FE::CellType::quad4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::quad4, CORE::FE::CellType::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::quad4, CORE::FE::CellType::tri3>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::quad4, CORE::FE::CellType::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::quad4, CORE::FE::CellType::tri3>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::quad4, CORE::FE::CellType::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::quad4, CORE::FE::CellType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::quad4, CORE::FE::CellType::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::quad4,
        CORE::FE::CellType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::quad4, CORE::FE::CellType::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::quad4, CORE::FE::CellType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::quad4, CORE::FE::CellType::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::quad4,
        CORE::FE::CellType::nurbs9>>;

template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::tri3, CORE::FE::CellType::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::tri3, CORE::FE::CellType::tri3>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::tri3, CORE::FE::CellType::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::tri3, CORE::FE::CellType::tri3>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::tri3, CORE::FE::CellType::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::tri3, CORE::FE::CellType::quad4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::tri3, CORE::FE::CellType::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::tri3, CORE::FE::CellType::quad4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::tri3, CORE::FE::CellType::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::tri3, CORE::FE::CellType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::tri3, CORE::FE::CellType::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::tri3,
        CORE::FE::CellType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::tri3, CORE::FE::CellType::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::tri3, CORE::FE::CellType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::tri3, CORE::FE::CellType::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::tri3,
        CORE::FE::CellType::nurbs9>>;

template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
        CORE::FE::CellType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
        CORE::FE::CellType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::tri3>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
        CORE::FE::CellType::tri3>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::quad4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
        CORE::FE::CellType::quad4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
        CORE::FE::CellType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs4, CORE::FE::CellType::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs4,
        CORE::FE::CellType::nurbs9>>;

template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::nurbs9,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
        CORE::FE::CellType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::nurbs9,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
        CORE::FE::CellType::nurbs9>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::tri3,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::tri3>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::tri3,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
        CORE::FE::CellType::tri3>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::quad4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::quad4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::quad4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
        CORE::FE::CellType::quad4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::nurbs4,
    CONTACT::AUG::DebugCompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
        CORE::FE::CellType::nurbs4>>;
template class CONTACT::AUG::Integrator<3, CORE::FE::CellType::nurbs9, CORE::FE::CellType::nurbs4,
    CONTACT::AUG::DebugIncompleteIntPolicy<3, CORE::FE::CellType::nurbs9,
        CORE::FE::CellType::nurbs4>>;

BACI_NAMESPACE_CLOSE

#endif