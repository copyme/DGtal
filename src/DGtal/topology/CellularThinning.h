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
 * @file CellularThinning.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2015/10/30
 *
 * Header file for module CellularThinning.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CellularThinning_RECURSES)
#error Recursive header files inclusion detected in CellularThinning.h
#else // defined(CellularThinning_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CellularThinning_RECURSES

#if !defined CellularThinning_h
/** Prevents repeated inclusion of headers. */
#define CellularThinning_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class CellularThinning
/**
 * Description of class 'CellularThinning' <p>
 * \brief Aim: Implements ParaDirCollapse from Chaussard's paper.
 * 
 * \todo write 2D version for students and ask them for 3D version.
 */
class CellularThinning
{
    // ----------------------- Standard services ------------------------------
public:

    /**
     * Destructor.
     */
    ~CellularThinning();

    // ----------------------- Interface --------------------------------------
public:

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

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    CellularThinning();

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    CellularThinning ( const CellularThinning & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    CellularThinning & operator= ( const CellularThinning & other );

    // ------------------------- Internals ------------------------------------
private:

}; // end of class CellularThinning


/**
 * Overloads 'operator<<' for displaying objects of class 'CellularThinning'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'CellularThinning' to write.
 * @return the output stream after the writing.
 */
std::ostream&
operator<< ( std::ostream & out, const CellularThinning & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/topology/CellularThinning.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CellularThinning_h

#undef CellularThinning_RECURSES
#endif // else defined(CellularThinning_RECURSES)
