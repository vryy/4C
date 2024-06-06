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

template class CONTACT::Aug::Integrator<2, Core::FE::CellType::line2, Core::FE::CellType::line2,
    CONTACT::Aug::DebugCompleteIntPolicy<2, Core::FE::CellType::line2, Core::FE::CellType::line2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::line2, Core::FE::CellType::line2,
    CONTACT::Aug::DebugIncompleteIntPolicy<2, Core::FE::CellType::line2,
        Core::FE::CellType::line2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::line2, Core::FE::CellType::nurbs2,
    CONTACT::Aug::DebugCompleteIntPolicy<2, Core::FE::CellType::line2, Core::FE::CellType::nurbs2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::line2, Core::FE::CellType::nurbs2,
    CONTACT::Aug::DebugIncompleteIntPolicy<2, Core::FE::CellType::line2,
        Core::FE::CellType::nurbs2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::line2, Core::FE::CellType::nurbs3,
    CONTACT::Aug::DebugCompleteIntPolicy<2, Core::FE::CellType::line2, Core::FE::CellType::nurbs3>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::line2, Core::FE::CellType::nurbs3,
    CONTACT::Aug::DebugIncompleteIntPolicy<2, Core::FE::CellType::line2,
        Core::FE::CellType::nurbs3>>;

template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs2, Core::FE::CellType::nurbs2,
    CONTACT::Aug::DebugCompleteIntPolicy<2, Core::FE::CellType::nurbs2,
        Core::FE::CellType::nurbs2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs2, Core::FE::CellType::nurbs2,
    CONTACT::Aug::DebugIncompleteIntPolicy<2, Core::FE::CellType::nurbs2,
        Core::FE::CellType::nurbs2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs2, Core::FE::CellType::line2,
    CONTACT::Aug::DebugCompleteIntPolicy<2, Core::FE::CellType::nurbs2, Core::FE::CellType::line2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs2, Core::FE::CellType::line2,
    CONTACT::Aug::DebugIncompleteIntPolicy<2, Core::FE::CellType::nurbs2,
        Core::FE::CellType::line2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs2, Core::FE::CellType::nurbs3,
    CONTACT::Aug::DebugCompleteIntPolicy<2, Core::FE::CellType::nurbs2,
        Core::FE::CellType::nurbs3>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs2, Core::FE::CellType::nurbs3,
    CONTACT::Aug::DebugIncompleteIntPolicy<2, Core::FE::CellType::nurbs2,
        Core::FE::CellType::nurbs3>>;

template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs3, Core::FE::CellType::nurbs3,
    CONTACT::Aug::DebugCompleteIntPolicy<2, Core::FE::CellType::nurbs3,
        Core::FE::CellType::nurbs3>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs3, Core::FE::CellType::nurbs3,
    CONTACT::Aug::DebugIncompleteIntPolicy<2, Core::FE::CellType::nurbs3,
        Core::FE::CellType::nurbs3>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs3, Core::FE::CellType::line2,
    CONTACT::Aug::DebugCompleteIntPolicy<2, Core::FE::CellType::nurbs3, Core::FE::CellType::line2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs3, Core::FE::CellType::line2,
    CONTACT::Aug::DebugIncompleteIntPolicy<2, Core::FE::CellType::nurbs3,
        Core::FE::CellType::line2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs3, Core::FE::CellType::nurbs2,
    CONTACT::Aug::DebugCompleteIntPolicy<2, Core::FE::CellType::nurbs3,
        Core::FE::CellType::nurbs2>>;
template class CONTACT::Aug::Integrator<2, Core::FE::CellType::nurbs3, Core::FE::CellType::nurbs2,
    CONTACT::Aug::DebugIncompleteIntPolicy<2, Core::FE::CellType::nurbs3,
        Core::FE::CellType::nurbs2>>;

template class CONTACT::Aug::Integrator<3, Core::FE::CellType::quad4, Core::FE::CellType::quad4,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::quad4, Core::FE::CellType::quad4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::quad4, Core::FE::CellType::quad4,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::quad4,
        Core::FE::CellType::quad4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::quad4, Core::FE::CellType::tri3,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::quad4, Core::FE::CellType::tri3>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::quad4, Core::FE::CellType::tri3,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::quad4, Core::FE::CellType::tri3>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::quad4, Core::FE::CellType::nurbs4,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::quad4, Core::FE::CellType::nurbs4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::quad4, Core::FE::CellType::nurbs4,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::quad4,
        Core::FE::CellType::nurbs4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::quad4, Core::FE::CellType::nurbs9,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::quad4, Core::FE::CellType::nurbs9>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::quad4, Core::FE::CellType::nurbs9,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::quad4,
        Core::FE::CellType::nurbs9>>;

template class CONTACT::Aug::Integrator<3, Core::FE::CellType::tri3, Core::FE::CellType::tri3,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::tri3, Core::FE::CellType::tri3>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::tri3, Core::FE::CellType::tri3,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::tri3, Core::FE::CellType::tri3>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::tri3, Core::FE::CellType::quad4,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::tri3, Core::FE::CellType::quad4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::tri3, Core::FE::CellType::quad4,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::tri3, Core::FE::CellType::quad4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::tri3, Core::FE::CellType::nurbs4,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::tri3, Core::FE::CellType::nurbs4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::tri3, Core::FE::CellType::nurbs4,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::tri3,
        Core::FE::CellType::nurbs4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::tri3, Core::FE::CellType::nurbs9,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::tri3, Core::FE::CellType::nurbs9>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::tri3, Core::FE::CellType::nurbs9,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::tri3,
        Core::FE::CellType::nurbs9>>;

template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs4, Core::FE::CellType::nurbs4,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::nurbs4,
        Core::FE::CellType::nurbs4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs4, Core::FE::CellType::nurbs4,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::nurbs4,
        Core::FE::CellType::nurbs4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs4, Core::FE::CellType::tri3,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::nurbs4, Core::FE::CellType::tri3>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs4, Core::FE::CellType::tri3,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::nurbs4,
        Core::FE::CellType::tri3>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs4, Core::FE::CellType::quad4,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::nurbs4, Core::FE::CellType::quad4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs4, Core::FE::CellType::quad4,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::nurbs4,
        Core::FE::CellType::quad4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs4, Core::FE::CellType::nurbs9,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::nurbs4,
        Core::FE::CellType::nurbs9>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs4, Core::FE::CellType::nurbs9,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::nurbs4,
        Core::FE::CellType::nurbs9>>;

template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs9, Core::FE::CellType::nurbs9,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::nurbs9,
        Core::FE::CellType::nurbs9>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs9, Core::FE::CellType::nurbs9,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::nurbs9,
        Core::FE::CellType::nurbs9>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs9, Core::FE::CellType::tri3,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::nurbs9, Core::FE::CellType::tri3>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs9, Core::FE::CellType::tri3,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::nurbs9,
        Core::FE::CellType::tri3>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs9, Core::FE::CellType::quad4,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::nurbs9, Core::FE::CellType::quad4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs9, Core::FE::CellType::quad4,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::nurbs9,
        Core::FE::CellType::quad4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs9, Core::FE::CellType::nurbs4,
    CONTACT::Aug::DebugCompleteIntPolicy<3, Core::FE::CellType::nurbs9,
        Core::FE::CellType::nurbs4>>;
template class CONTACT::Aug::Integrator<3, Core::FE::CellType::nurbs9, Core::FE::CellType::nurbs4,
    CONTACT::Aug::DebugIncompleteIntPolicy<3, Core::FE::CellType::nurbs9,
        Core::FE::CellType::nurbs4>>;

FOUR_C_NAMESPACE_CLOSE

#endif