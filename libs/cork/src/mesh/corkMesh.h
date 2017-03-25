#include "mesh.h"
#include "../rawmesh/rawMesh.h"

#pragma once

struct CorkTriangle;

struct CorkVertex :
    public MinimalVertexData,
    public RemeshVertexData,
    public IsctVertexData,
    public BoolVertexData
{
    void merge(const CorkVertex &v0, const CorkVertex &v1) {
        double                              a0 = 0.5;
        if(v0.manifold && !v1.manifold)     a0 = 0.0;
        if(!v0.manifold && v1.manifold)     a0 = 1.0;
        double a1 = 1.0 - a0;
        
        pos         = a0 * v0.pos       + a1 * v1.pos;
    }
    void interpolate(const CorkVertex &v0, const CorkVertex &v1) {
        double a0 = 0.5;
        double a1 = 0.5;
        pos         = a0 * v0.pos       + a1 * v1.pos;
    }
    
    
    void isct(IsctVertEdgeTriInput<CorkVertex,CorkTriangle> input)
    {
        Vec2d       a_e     = Vec2d(1,1)/2.0;
        Vec3d       a_t     = Vec3d(1,1,1)/3.0;
        a_e /= 2.0;
        a_t /= 2.0;
    }
    void isct(IsctVertTriTriTriInput<CorkVertex,CorkTriangle> input)
    {
        Vec3d       a[3];
        for(uint k=0; k<3; k++) {
            a[k]    = Vec3d(1,1,1)/3.0;
            a[k] /= 3.0;
        }
        for(uint i=0; i<3; i++) {
          for(uint j=0; j<3; j++) {
        }}
    }
    void isctInterpolate(const CorkVertex &v0, const CorkVertex &v1) {
        double a0 = len(v1.pos - pos);
        double a1 = len(v0.pos - pos);
        if(a0 + a1 == 0.0) a0 = a1 = 0.5; // safety
        double sum = a0+a1;
        a0 /= sum;
        a1 /= sum;
    }
};

struct CorkTriangle :
    public MinimalTriangleData,
    public RemeshTriangleData,
    public IsctTriangleData,
    public BoolTriangleData
{
    void merge(const CorkTriangle &, const CorkTriangle &) {}
    static void split(CorkTriangle &, CorkTriangle &,
                      const CorkTriangle &) {}
    void move(const CorkTriangle &) {}
    void subdivide(SubdivideTriInput<CorkVertex,CorkTriangle> input)
    {
        bool_alg_data = input.pt->bool_alg_data;
    }
};

//using RawCorkMesh = RawMesh<CorkVertex, CorkTriangle>;
//using CorkMesh = Mesh<CorkVertex, CorkTriangle>;
typedef Mesh<CorkVertex, CorkTriangle> CorkMesh;
typedef RawMesh<CorkVertex, CorkTriangle> RawCorkMesh;

