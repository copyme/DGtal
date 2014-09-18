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
 * @file testtestQuasiAffineTransformation2D.cpp
 * @ingroup Tests
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2014/09/16
 *
 * Functions for testing class testQuasiAffineTransformation2D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "ConfigTest.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/QuasiAffineTransformation2D.h"
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/ConstImageAdapter.h"

#ifndef NDEBUG
#include <fenv.h>
#endif


///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace functors;
///////////////////////////////////////////////////////////////////////////////
// Functions for testing class testQuasiAffineTransformation2D.
///////////////////////////////////////////////////////////////////////////////
/**
 * Example of a test. To be completed.
 *
 */
bool testtestQuasiAffineTransformation2D()
{
  unsigned int nbok = 0;
  unsigned int nb = 0;
  
  trace.beginBlock ( "Testing block ..." );
  nbok += true ? 1 : 0; 
  nb++;
  trace.info() << "(" << nbok << "/" << nb << ") "
	       << "true == true" << std::endl;
  trace.endBlock();
  
  return nbok == nb;
}

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{

#ifndef NDEBUG
    fetestexcept ( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW );
#endif

  typedef ForwardsQuasiAffineTransformation2D<Z2i::Space> QuasiAffineTrans;
  typedef BackwardsQuasiAffineTransformation2D <Z2i::Space> BackwardsNN2D;
  typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image;
  typedef QuasiAffineTransformer2D < Image > QuasiAffineTransformer2D;
  typedef ConstImageAdapter<Image, Z2i::Domain, BackwardsNN2D, Image::Value, Identity > MyImageBackwardAdapter;
  typedef DomainQuasiAffineTransformation2D<Z2i::Domain, QuasiAffineTrans> DomainQuasiAffineTrans;
  typedef DomainQuasiAffineTransformation2D<Z2i::Domain, QuasiAffineTrans>::Bounds Bounds;
  QuasiAffineTrans::Matrix mat;

  mat(0,0) = 3;
  mat(0,1) = 1;
  mat(1,0) = -1;
  mat(1,1) = 3;
  QuasiAffineTrans trans(mat, 6, Z2i::Vector(0,0));
  DomainQuasiAffineTrans domainTrans (  trans );
  Image image = PGMReader<Image>::importPGM(testPath+"samples/church.pgm");
  Bounds bounds = domainTrans ( image.domain() );

  Identity idD;
  BackwardsNN2D backwardsnn ( mat,6, Z2i::Vector(0,0) );

  MyImageBackwardAdapter adapter ( image, Z2i::Domain ( bounds.first, bounds.second ) , backwardsnn, idD );

  adapter >> "test-quasi_backNN.pgm";

  QuasiAffineTransformer2D imageTransformer(mat, 6, Z2i::Vector(0,0));

  Image final (Z2i::Domain ( bounds.first, bounds.second ) );
  Image final2 (Z2i::Domain ( bounds.first, bounds.second ) );

  imageTransformer.BackwardsNearestNeighborInterpolation ( image, final );
  imageTransformer.BackwardsLinearInterpolation ( image, final2 );

  final >> "test-image-trans-backwardsNN.pgm";
  final2 >> "test-image-trans-linearI.pgm";

  trace.beginBlock ( "Testing class testQuasiAffineTransformation2D" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  bool res = testtestQuasiAffineTransformation2D(); // && ... other tests
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
//  return res ? 0 : 1;
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
