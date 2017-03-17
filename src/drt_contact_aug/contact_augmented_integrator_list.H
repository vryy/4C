/*---------------------------------------------------------------------*/
/*!
\file contact_augmented_integrator_list.H

\brief List of template combinations for the augmented contact integrator.

\level 2

\maintainer Michael Hiermeier

\date Mar 8, 2017

*/
/*---------------------------------------------------------------------*/

template class CONTACT::AUG::Integrator<2,DRT::Element::line2,DRT::Element::line2>;
template class CONTACT::AUG::Integrator<3,DRT::Element::line2,DRT::Element::line2>;
template class CONTACT::AUG::Integrator<3,DRT::Element::line2,DRT::Element::quad4>;
template class CONTACT::AUG::Integrator<3,DRT::Element::line2,DRT::Element::tri3>;

template class CONTACT::AUG::Integrator<3,DRT::Element::quad4,DRT::Element::line2>;
template class CONTACT::AUG::Integrator<3,DRT::Element::quad4,DRT::Element::quad4>;
template class CONTACT::AUG::Integrator<3,DRT::Element::quad4,DRT::Element::tri3>;

template class CONTACT::AUG::Integrator<3,DRT::Element::tri3,DRT::Element::line2>;
template class CONTACT::AUG::Integrator<3,DRT::Element::tri3,DRT::Element::quad4>;
template class CONTACT::AUG::Integrator<3,DRT::Element::tri3,DRT::Element::tri3>;

