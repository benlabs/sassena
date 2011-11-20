/** \file
This file contains helper class to generate a icosahedral geodesic sphere.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/

#ifndef MATH__GEOMETRY_HPP_
#define MATH__GEOMETRY_HPP_

// common header
#include "common.hpp"

// standard header
#include <vector>

// special library headers

// other headers
#include "math/coor3d.hpp"

/** 
Helper class which is used to construct evenly distributed grid points on a sphere. This class is adapted from code on:http://local.wasp.uwa.edu.au/~pbourke/geometry/platonic/
*/
class DrawSphereHelper {
public:
	std::vector<CartesianCoor3D> vectors;
	
    CartesianCoor3D center;
    float radius;
    float maxEdgeLength;

    // "kitchen sink" constructor (allows specifying everything)
    DrawSphereHelper (CartesianCoor3D _center,
                      const double _radius,
                      const double _maxEdgeLength)
        : center (_center),
          radius (_radius),
          maxEdgeLength (_maxEdgeLength)
    {}

    // draw as an icosahedral geodesic sphere
    void draw (void) 
    {
        // Geometry based on Paul Bourke's excellent article:
        //   Platonic Solids (Regular polytopes in 3D)
        //   http://astronomy.swin.edu.au/~pbourke/polyhedra/platonic/
        float sqrt5 = sqrt (5.0f);
        float phi = (1.0f + sqrt5) * 0.5f; // "golden ratio"
        // ratio of edge length to radius
        float ratio = sqrt (10.0f + (2.0f * sqrt5)) / (4.0f * phi);
        float a = (radius / ratio) * 0.5;
        float b = (radius / ratio) / (2.0f * phi);

        // define the icosahedron's 12 vertices:
        CartesianCoor3D v1  = center + CartesianCoor3D( 0,  b, -a);
        CartesianCoor3D v2  = center + CartesianCoor3D( b,  a,  0);
        CartesianCoor3D v3  = center + CartesianCoor3D(-b,  a,  0);
        CartesianCoor3D v4  = center + CartesianCoor3D( 0,  b,  a);
        CartesianCoor3D v5  = center + CartesianCoor3D( 0, -b,  a);
        CartesianCoor3D v6  = center + CartesianCoor3D(-a,  0,  b);
        CartesianCoor3D v7  = center + CartesianCoor3D( 0, -b, -a);
        CartesianCoor3D v8  = center + CartesianCoor3D( a,  0, -b);
        CartesianCoor3D v9  = center + CartesianCoor3D( a,  0,  b);
        CartesianCoor3D v10 = center + CartesianCoor3D(-a,  0, -b);
        CartesianCoor3D v11 = center + CartesianCoor3D( b, -a,  0);
        CartesianCoor3D v12 = center + CartesianCoor3D(-b, -a,  0);

        // draw the icosahedron's 20 triangular faces:
        drawMeshedTriangleOnSphere (v1, v2, v3);
        drawMeshedTriangleOnSphere (v4, v3, v2);
        drawMeshedTriangleOnSphere (v4, v5, v6);
        drawMeshedTriangleOnSphere (v4, v9, v5);
        drawMeshedTriangleOnSphere (v1, v7, v8);
        drawMeshedTriangleOnSphere (v1, v10, v7);
        drawMeshedTriangleOnSphere (v5, v11, v12);
        drawMeshedTriangleOnSphere (v7, v12, v11);
        drawMeshedTriangleOnSphere (v3, v6, v10);
        drawMeshedTriangleOnSphere (v12, v10, v6);
        drawMeshedTriangleOnSphere (v2, v8, v9);
        drawMeshedTriangleOnSphere (v11, v9, v8);
        drawMeshedTriangleOnSphere (v4, v6, v3);
        drawMeshedTriangleOnSphere (v4, v2, v9);
        drawMeshedTriangleOnSphere (v1, v3, v10);
        drawMeshedTriangleOnSphere (v1, v8, v2);
        drawMeshedTriangleOnSphere (v7, v10, v12);
        drawMeshedTriangleOnSphere (v7, v11, v8);
        drawMeshedTriangleOnSphere (v5, v12, v6);
        drawMeshedTriangleOnSphere (v5, v9, v11);
    }

    // given two points, take midpoint and project onto this sphere
    inline CartesianCoor3D midpointOnSphere (CartesianCoor3D& a, CartesianCoor3D& b) 
    {
        CartesianCoor3D midpoint = (a + b) * 0.5f;
        CartesianCoor3D unitRadial = (midpoint - center);

        return center + ((unitRadial * radius)/unitRadial.length());
    }

    // given three points on the surface of this sphere, if the triangle
    // is "small enough" draw it, otherwise subdivide it into four smaller
    // triangles and recursively draw each of them.
    void drawMeshedTriangleOnSphere (CartesianCoor3D& a,
                                     CartesianCoor3D& b,
                                     CartesianCoor3D& c) 
    {
        // if all edges are short enough
        if ((((a - b).length ()) < maxEdgeLength) &&
            (((b - c).length ()) < maxEdgeLength) &&
            (((c - a).length ()) < maxEdgeLength))
        {
            // draw triangle
            drawTriangleOnSphere (a, b, c);
        }
        else // otherwise subdivide and recurse
        {
            // find edge midpoints
             CartesianCoor3D ab = midpointOnSphere (a, b);
             CartesianCoor3D bc = midpointOnSphere (b, c);
             CartesianCoor3D ca = midpointOnSphere (c, a);

            // recurse on four sub-triangles
            drawMeshedTriangleOnSphere ( a, ab, ca);
            drawMeshedTriangleOnSphere (ab,  b, bc);
            drawMeshedTriangleOnSphere (ca, bc,  c);
            drawMeshedTriangleOnSphere (ab, bc, ca);
        }
    }

    // draw one mesh element for drawMeshedTriangleOnSphere
    void drawTriangleOnSphere ( CartesianCoor3D& a,
                                CartesianCoor3D& b,
                                CartesianCoor3D& c) 
    {
        // draw triangle, subject to the camera orientation criteria
        // (according to drawBackFacing and drawFrontFacing)
         CartesianCoor3D triCenter = (a + b + c) / 3.0f;
         CartesianCoor3D triNormal = triCenter - center; // not unit length
		 vectors.push_back( triNormal/triNormal.length() );
    }

};

#endif

// end of file
