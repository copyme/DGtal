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

/**
 * @file testHomotopicThinning.cpp
 * @ingroup Tests
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2014/11/13
 *
 * Functions for testing class HomotopicThinning.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <functional>
#include <vector>
#include <limits>
#include "DGtal/base/Common.h"
#include "ConfigTest.h"
#include "DGtal/helpers/StdDefs.h"
#include <DGtal/images/ImageSelector.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include "DGtal/topology/HomotopicThinning.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"

#ifdef __GNUC__
  #ifndef NDEBUG
    #include <fenv.h>
  #endif
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace functors;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class HomotopicThinning.
///////////////////////////////////////////////////////////////////////////////
/**
 * Example of a test. To be completed.
 *
 */

class testHomotopicThinning
{
public:
  typedef ImageSelector<Z2i::Domain, unsigned char >::Type Image2D;
  typedef ImageSelector<Z3i::Domain, unsigned char >::Type Image3D;
  typedef IntervalForegroundPredicate<Image2D> Predicate2D;
  typedef DistanceTransformation<Z2i::Space, Predicate2D , Z2i::L2Metric> DistanceTrans2D;
  typedef std::function<int(Z3i::Point)> Zeros;
  typedef HomotopicThinning < Z3i::Object26_6, std::vector < Z3i::Point >, Zeros > HomoThinnSimple3D;
  typedef HomotopicThinning < Z2i::Object8_4, std::vector < Z2i::Point >, DistanceTrans2D > HomoThinnSimple2D;
  testHomotopicThinning () : image3d ( PGMReader<Image3D>::importPGM3D ( testPath + "samples/cat10.pgm3d" ) ),
  image2d ( PGMReader<Image2D>::importPGM ( testPath + "samples/church-small.pgm" ) ),
  zeros ([](Z3i::Point) { return 0; } ),
  aPredicate ( image2d, 125, 255 ),
  disTrans2D ( image2d.domain(),aPredicate, L2metric ),
  homoThinnSimple3D ( zeros ),
  homoThinnSimple2D ( disTrans2D )
  {
  }
  bool Binary3D()
  {
    using namespace Z3i;
    DigitalSet shape_set( image3d.domain() );
    SetFromImage<DigitalSet>::append<Image3D>(shape_set, image3d, 1, 255 );
    Object26_6 shape ( dt26_6, shape_set );
    homoThinnSimple3D ( shape, 3 );
    Point low, up;
    shape.pointSet().computeBoundingBox ( low, up );
    Image3D output (Domain ( low, up ) );
    for ( DigitalSet::ConstIterator it = shape.begin(); it != shape.end(); ++it )
      output.setValue ( *it, 255 );
    PGMWriter<Image3D>::exportPGM3D ( "output_homothinn.pgm3d", output );
    return true;
  }
  bool Gray2D()
  {
    using namespace Z2i;
    DigitalSet shape_set( image2d.domain() );
    SetFromImage<DigitalSet>::append<Image2D>(shape_set, image2d, 125, 255 );
    Object8_4 shape ( dt8_4, shape_set );
    homoThinnSimple2D ( shape, 3 );
    Point low, up;
    shape.pointSet().computeBoundingBox ( low, up );
    Image2D output (Domain ( low, up ) );
    for ( DigitalSet::ConstIterator it = shape.begin(); it != shape.end(); ++it )
      output.setValue ( *it, 255 );
    PGMWriter<Image2D>::exportPGM ( "output_homothinn.pgm", output );
    return true;
  }
  
  bool Gray2DFixSet()
  {
    using namespace Z2i;
    DigitalSet shape_set( image2d.domain() );
    DigitalSet fixedShape( image2d.domain() );
    SetFromImage<DigitalSet>::append<Image2D>(shape_set, image2d, 0, 255 );
    SetFromImage<DigitalSet>::append<Image2D>(fixedShape, image2d, 125, 190 );
    std::vector < Point > fixPoints;
    for ( DigitalSet::ConstIterator it = fixedShape.begin(); it != fixedShape.end(); ++it )
      fixPoints.push_back ( *it );
    Object8_4 shape ( dt8_4, shape_set );
    homoThinnSimple2D ( shape, fixPoints, image2d.domain().size() );
    Point low, up;
    shape.pointSet().computeBoundingBox ( low, up );
    Image2D output (Domain ( low, up ) );
    for ( DigitalSet::ConstIterator it = shape.begin(); it != shape.end(); ++it )
      output.setValue ( *it, 255 );
    PGMWriter<Image2D>::exportPGM ( "output_homothinn_fix.pgm", output );
    return true;
  }
  
private:
  Image3D image3d;
  Image2D image2d;
  Zeros zeros;
  Predicate2D aPredicate;
  Z2i::L2Metric L2metric;
  DistanceTrans2D disTrans2D;
  HomoThinnSimple3D homoThinnSimple3D;
  HomoThinnSimple2D homoThinnSimple2D;
};

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int , char** )
{
  #ifdef __GNUC__
    #ifndef NDEBUG
      feenableexcept ( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW );
    #endif
  #endif
  
  bool res = true;
  testHomotopicThinning testHomoThinn;
  trace.beginBlock ( "Testing homotopic thinning." );
  res &= testHomoThinn.Binary3D();
  res &= testHomoThinn.Gray2D();
  res &= testHomoThinn.Gray2DFixSet();
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
