/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file HomotopicThinning.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2014/11/13
 *
 * Header file for module HomotopicThinning.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(HomotopicThinning_RECURSES)
#error Recursive header files inclusion detected in HomotopicThinning.h
#else // defined(HomotopicThinning_RECURSES)
/** Prevents recursive inclusion of headers. */
#define HomotopicThinning_RECURSES

#if !defined HomotopicThinning_h
/** Prevents repeated inclusion of headers. */
#define HomotopicThinning_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class HomotopicThinning
/**
 * Description of class 'HomotopicThinning' <p>
 * \brief Aim:
 */
template <class ObjectA_B, class FixedSet, class PriorityFunctor>
class HomotopicThinning
{
    // ----------------------- Standard services ------------------------------
public:
  typedef typename ObjectA_B::DigitalSet DigitalSet;
    HomotopicThinning ( PriorityFunctor & );

    /**
     * Destructor.
     */
    ~HomotopicThinning(){}

    // ----------------------- Interface --------------------------------------
public:
  
    void operator()( ObjectA_B &, const FixedSet &, const unsigned int & );
    void operator()( ObjectA_B &, const unsigned int & );

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;
    

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
private:
    // ------------------------- Private Datas --------------------------------
private:

    // ------------------------- Hidden services ------------------------------
protected:

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    HomotopicThinning ( const HomotopicThinning & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    HomotopicThinning & operator= ( const HomotopicThinning & other );

    // ------------------------- Internals ------------------------------------
private:
  PriorityFunctor & functor;

}; // end of class HomotopicThinning


/**
 * Overloads 'operator<<' for displaying objects of class 'HomotopicThinning'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'HomotopicThinning' to write.
 * @return the output stream after the writing.
 */
template <class ObjectA_B, class FixedSet, class PriorityFunctor>
std::ostream&
operator<< ( std::ostream & out, const HomotopicThinning<ObjectA_B, FixedSet, PriorityFunctor> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/topology/HomotopicThinning.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined HomotopicThinning_h

#undef HomotopicThinning_RECURSES
#endif // else defined(HomotopicThinning_RECURSES)
