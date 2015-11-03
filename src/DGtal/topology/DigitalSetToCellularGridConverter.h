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
 * @file DigitalSetToCellularGridConverter.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2015/11/03
 *
 * Header file for module DigitalSetToCellularGridConverter.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DigitalSetToCellularGridConverter_RECURSES)
#error Recursive header files inclusion detected in DigitalSetToCellularGridConverter.h
#else // defined(DigitalSetToCellularGridConverter_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DigitalSetToCellularGridConverter_RECURSES

#if !defined DigitalSetToCellularGridConverter_h
/** Prevents repeated inclusion of headers. */
#define DigitalSetToCellularGridConverter_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/topology/CCellularGridSpaceND.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class DigitalSetToCellularGridConverter
/**
 * Description of class 'DigitalSetToCellularGridConverter' <p>
 * \brief Aim: Convert a DigitalSet to cubical complexes
 * @tparam TKSpace - Khalimsky space
 * @tparam TDigitalSet - a type of digital set.
 */
template < typename TKSpace, typename TDigitalSet >
class DigitalSetToCellularGridConverter
{
  //Checking concepts
  BOOST_STATIC_ASSERT(( TKSpace::dimension == TDigitalSet::Domain::dimension ));
  BOOST_CONCEPT_ASSERT(( concepts::CCellularGridSpaceND<TKSpace> ));
  // ----------------------- Types ------------------------------
public:
  /// Type of iterator, at least readable and forward
  typedef typename TDigitalSet::ConstIterator ConstIterator;
  /// Type of set which stores cells and guarantee that each cell can be added only once.
  typedef typename TKSpace::CellSet Cells;
  /// Type of cell.
  typedef typename TKSpace::Cell Cell;
  /// Type of iterator over cells.
  typedef typename Cells::const_iterator ConstCellsIterator;
    // ----------------------- Standard services ------------------------------
public:

    /**
     * Constructor
     * @param kspace - constant reference to Khalimsky space.
     */
    DigitalSetToCellularGridConverter ( const TKSpace & kspace );
    
    /**
     * Destructor.
     */
    ~DigitalSetToCellularGridConverter();

    // ----------------------- Interface --------------------------------------
public:
  
    /**
     * Initialization.
     * @param itb begin iterator
     * @param ite end iterator
     */
    void init ( const ConstIterator & itb, const ConstIterator & ite );
    
    /**
     * Extract all faces of spels which have a given dimension.
     * @param cells - collection of nD cells.
     * @param dim   - dimension of cells to extract.
     */
    void extractAllCells ( Cells & cells, typename TKSpace::Integer dim );
    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Hidden services ------------------------------
protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    DigitalSetToCellularGridConverter();

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    DigitalSetToCellularGridConverter ( const DigitalSetToCellularGridConverter & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    DigitalSetToCellularGridConverter & operator= ( const DigitalSetToCellularGridConverter & other );

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

}; // end of class DigitalSetToCellularGridConverter

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/topology/DigitalSetToCellularGridConverter.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DigitalSetToCellularGridConverter_h

#undef DigitalSetToCellularGridConverter_RECURSES
#endif // else defined(DigitalSetToCellularGridConverter_RECURSES)
