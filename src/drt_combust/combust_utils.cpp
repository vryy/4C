/*!---------------------------------------------------------------------*
\file combust_utils.cpp

\brief collection of functions in namespace COMBUST

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <iostream>

#include "combust_utils.H"

/*----------------------------------------------------------------------*
 | print COMBUST module logo on screen                      henke 06/08 |
 *----------------------------------------------------------------------*/
void COMBUST::printlogo()
{
    std::cout<<"     ___            ___    "<<std::endl;
    std::cout<<"    /   \\          /   \\ "<<std::endl;
    std::cout<<"    \\_   \\        /  __/ "<<std::endl;
    std::cout<<"     _\\   \\      /  /__  "<<" Das ist               "<<std::endl;
    std::cout<<"     \\___  \\____/   __/  "<<" das Verbrennungsmodul "<<std::endl;
    std::cout<<"         \\_       _/      "<<" in BACI               "<<std::endl;
    std::cout<<"           | @ @  \\_      "<<"                       "<<std::endl;
    std::cout<<"           |               "<<" Der Elch wird bald    "<<std::endl;
    std::cout<<"         _/     /\\        "<<" ein feuerspeiender    "<<std::endl;
    std::cout<<"        /o)  (o/\\ \\_     "<<" Drache sein!          "<<std::endl;
    std::cout<<"        \\_____/ /         "<<std::endl;
    std::cout<<"          \\____/          "<<std::endl;
    std::cout<<"                           "<<std::endl;
} // void printlogo()

#endif  // #ifdef CCADISCRET
