/*---------------------------------------------------------------------*/
/*! \file
\brief List of template combinations for the augmented contact integrator.

\level 2

*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_AUG_INTEGRATOR_INST_HPP
#define FOUR_C_CONTACT_AUG_INTEGRATOR_INST_HPP

#include "4C_config.hpp"

#include "4C_contact_aug_integrator_policy.hpp"

FOUR_C_NAMESPACE_OPEN

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

FOUR_C_NAMESPACE_CLOSE

#endif