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
 * @file DigitalSetToCellular2D.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2015/11/03
 *
 * Header file for module DigitalSetToCellular2D.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DigitalSetToCellular2D_RECURSES)
#error Recursive header files inclusion detected in DigitalSetToCellular2D.h
#else // defined(DigitalSetToCellular2D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DigitalSetToCellular2D_RECURSES

#if !defined DigitalSetToCellular2D_h
/** Prevents repeated inclusion of headers. */
#define DigitalSetToCellular2D_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/topology/CCellularGridSpaceND.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class DigitalSetToCellular2D
/**
 * Description of class 'DigitalSetToCellular2D' <p>
 * \brief Aim: Convert a DigitalSet to cubical complexes
 * 
 * \todo Write a 3D implementation.
 */
template < typename TKSpace, typename TDigitalSet >
class DigitalSetToCellular2D
{
  //Checking concepts
  BOOST_STATIC_ASSERT(( TKSpace::dimension == 2 ));
  BOOST_CONCEPT_ASSERT(( concepts::CCellularGridSpaceND<TKSpace> ));
  // ----------------------- Types ------------------------------
public:
  typedef typename TDigitalSet::ConstIterator ConstIterator;
  
  typedef typename TKSpace::CellSet Cells;
  typedef typename TKSpace::Cell Cell;
  typedef typename Cells::const_iterator ConstCellsIterator;
    // ----------------------- Standard services ------------------------------
public:

    /**
     * Destructor.
     */
    ~DigitalSetToCellular2D();
    DigitalSetToCellular2D ( const TKSpace & kspace );
    
    /**
     * Initialization.
     * @param itb begin iterator
     * @param ite end iterator
     */
    void init ( const ConstIterator & itb, const ConstIterator & ite );
    
    /**
     * \todo there must be the faster way without std::find, try to find it.
     * Converts digital set to three sets which contains different cells.
     * @param cells0d - collection of 0D cells
     * @param cells1d - collections of 1D cells
     * @param cells2d - collections of 2D cells
     */
    void toCellularGrid ( Cells & cells0d, Cells & cells1d, Cells & cells2d );

    // ----------------------- Interface --------------------------------------
public:

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
    DigitalSetToCellular2D();

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    DigitalSetToCellular2D ( const DigitalSetToCellular2D & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    DigitalSetToCellular2D & operator= ( const DigitalSetToCellular2D & other );

    // ------------------------- Internals ------------------------------------
private:
  const TKSpace & K;
  /**
   * Iterator which corresponds to the beginning of a valid range - [myBegin, myEnd)
   */
  ConstIterator myBegin;
  /**
   * Iterator which corresponds to the end of a valid range - [myBegin, myEnd)
   */
  ConstIterator myEnd;

}; // end of class DigitalSetToCellular2D

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/topology/DigitalSetToCellular2D.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DigitalSetToCellular2D_h

#undef DigitalSetToCellular2D_RECURSES
#endif // else defined(DigitalSetToCellular2D_RECURSES)
