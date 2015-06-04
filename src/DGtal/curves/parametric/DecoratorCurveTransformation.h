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
 * @file DecoratorCurveTransformation.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
 *
 * @date 2014/10/10
 *
 * Header file for module DecoratorCurveTransformation.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DecoratorCurveTransformation_RECURSES)
#error Recursive header files inclusion detected in DecoratorCurveTransformation.h
#else // defined(DecoratorCurveTransformation_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DecoratorCurveTransformation_RECURSES

#if !defined DecoratorCurveTransformation_h
/** Prevents repeated inclusion of headers. */
#define DecoratorCurveTransformation_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/curves/parametric/C3DParametricCurve.h"
#include "DGtal/curves/ParametricCurveDigitizer3D.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class DecoratorCurveTransformation
/**
 * Description of class 'DecoratorCurveTransformation' <p>
 * \brief Aim:
 */
template <typename TCurve, typename TTransfromation, typename TInverseTransformation>
class DecoratorCurveTransformation
{
    BOOST_CONCEPT_ASSERT(( concepts::C3DParametricCurve < TCurve > ));
    // ----------------------- Standard services ------------------------------
public:

    typedef TCurve TypeCurve;
    typedef typename TCurve::Space Space;
    typedef typename Space::RealPoint RealPoint;
    typedef typename Space::Point Point;

    DecoratorCurveTransformation ( const TCurve &, const TTransfromation &, const TInverseTransformation & );
    /**
     * Destructor.
     */
    ~DecoratorCurveTransformation() {}

    // ----------------------- Interface --------------------------------------
public:
  
    const TCurve & curve;

    /**
     * @param t any angle between 0 and k*Pi.
     *
     * @return the vector (x(t),y(t), z(t))
     */
    RealPoint x ( double t ) const;

    /**
     * @param t any angle between 0 and k*Pi.
     *
     * @return the vector (x(t)',y(t)', z(t)')
     */
    RealPoint xp ( double t ) const;

    /**
     * @brief inverse function of x
     * @param p = x(t)
     * @return t
     */
    double f ( const RealPoint & p ) const;

    /**
     * @brief inverse function of y
     * @param p = x(t)[
     * @return t
     */
    double g ( const RealPoint & p ) const;

    /**
     * @brief inverse function of z
     * @param p = x(t)
     * @return t
     */
    double h ( const RealPoint & ) const;

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
    DecoratorCurveTransformation();

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    DecoratorCurveTransformation ( const DecoratorCurveTransformation & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    DecoratorCurveTransformation & operator= ( const DecoratorCurveTransformation & other );

    // ------------------------- Internals ------------------------------------
private:
    const TTransfromation & trans;
    const TInverseTransformation & inverse;

}; // end of class DecoratorCurveTransformation


/**
 * Overloads 'operator<<' for displaying objects of class 'DecoratorCurveTransformation'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'DecoratorCurveTransformation' to write.
 * @return the output stream after the writing.
 */
template <typename TCurve, typename TTransfromation, typename TInverseTransformation>
inline
std::ostream&
operator<< ( std::ostream & out, const DecoratorCurveTransformation< TCurve, TTransfromation, TInverseTransformation > & object );


template <typename TParametricCurve>
class ParametricCurveDigitizer3DDecorator : public ParametricCurveDigitizer3D<TParametricCurve>
{
  BOOST_CONCEPT_ASSERT(( concepts::C3DParametricCurveDecorator< TParametricCurve > ));
  
private:
  virtual unsigned char findMainAxis ( const char & blockAxis, const double & i );
};

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/curves/parametric/DecoratorCurveTransformation.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DecoratorCurveTransformation_h

#undef DecoratorCurveTransformation_RECURSES
#endif // else defined(DecoratorCurveTransformation_RECURSES)
