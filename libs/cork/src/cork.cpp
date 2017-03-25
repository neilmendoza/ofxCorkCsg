// +-------------------------------------------------------------------------
// | cork.cpp
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#include "cork.h"

#include "mesh/mesh.h"
#include "mesh/corkMesh.h"

void freeCorkTriMesh(CorkTriMesh *mesh)
{
	if (mesh->triangles)
		delete[] mesh->triangles;
	if (mesh->vertices)
		delete[] mesh->vertices;
    mesh->n_triangles = 0;
    mesh->n_vertices = 0;
}

void corkTriMesh2CorkMesh(
    CorkTriMesh in,
    CorkMesh *mesh_out
) {
    RawCorkMesh raw;
    raw.vertices.resize(in.n_vertices);
    raw.triangles.resize(in.n_triangles);
    if(in.n_vertices == 0 || in.n_triangles == 0) {
        CORK_ERROR("empty mesh input to Cork routine.");
        *mesh_out = CorkMesh(raw);
        return;
    }
    
    uint max_ref_idx = 0;
    for(uint i=0; i<in.n_triangles; i++) {
        raw.triangles[i].a = in.triangles[3*i+0];
        raw.triangles[i].b = in.triangles[3*i+1];
        raw.triangles[i].c = in.triangles[3*i+2];
        max_ref_idx = std::max(
                        std::max(max_ref_idx,
                                 in.triangles[3*i+0]),
                        std::max(in.triangles[3*i+1],
                                 in.triangles[3*i+2])
                      );
    }
    if(max_ref_idx > in.n_vertices) {
        CORK_ERROR("mesh input to Cork routine has an out of range reference "
              "to a vertex.");
        raw.vertices.clear();
        raw.triangles.clear();
        *mesh_out = CorkMesh(raw);
        return;
    }
    
    for(uint i=0; i<in.n_vertices; i++) {
        raw.vertices[i].pos.x = in.vertices[3*i+0];
        raw.vertices[i].pos.y = in.vertices[3*i+1];
        raw.vertices[i].pos.z = in.vertices[3*i+2];
    }
    
    *mesh_out = CorkMesh(raw);
}
void corkMesh2CorkTriMesh(
    CorkMesh *mesh_in,
    CorkTriMesh *out
) {
    RawCorkMesh raw = mesh_in->raw();
    
    out->n_triangles = static_cast<uint>(raw.triangles.size());
    out->n_vertices  = static_cast<uint>(raw.vertices.size());
    
    out->triangles = new uint[(out->n_triangles) * 3];
    out->vertices  = new float[(out->n_vertices) * 3];
    
    for(uint i=0; i<out->n_triangles; i++) {
        (out->triangles)[3*i+0] = raw.triangles[i].a;
        (out->triangles)[3*i+1] = raw.triangles[i].b;
        (out->triangles)[3*i+2] = raw.triangles[i].c;
    }
    
	for(uint i=0; i<out->n_vertices; i++)
	{
		(out->vertices)[3*i+0] = static_cast<float>(raw.vertices[i].pos.x);
		(out->vertices)[3*i+1] = static_cast<float>(raw.vertices[i].pos.y);
		(out->vertices)[3*i+2] = static_cast<float>(raw.vertices[i].pos.z);
	}
}


bool isSolid(CorkTriMesh cmesh)
{
    CorkMesh mesh;
    corkTriMesh2CorkMesh(cmesh, &mesh);
    
    bool solid = true;
    
    if(mesh.isSelfIntersecting()) {
        CORK_ERROR("isSolid() was given a self-intersecting mesh");
        solid = false;
    }
    
    if(!mesh.isClosed()) {
        CORK_ERROR("isSolid() was given a non-closed mesh");
        solid = false;
    }
    
    return solid;
}

void computeUnion(
    CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
) {
    CorkMesh cmIn0, cmIn1;
    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.boolUnion(cmIn1);
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}

void computeDifference(
    CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
) {
    CorkMesh cmIn0, cmIn1;
    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.boolDiff(cmIn1);
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}

void computeIntersection(
    CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
) {
    CorkMesh cmIn0, cmIn1;
    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.boolIsct(cmIn1);
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}

void computeSymmetricDifference(
    CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
) {
    CorkMesh cmIn0, cmIn1;
    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.boolXor(cmIn1);
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}

void resolveIntersections(
    CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
) {
    CorkMesh cmIn0, cmIn1;
    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);
    
    cmIn0.disjointUnion(cmIn1);
    cmIn0.resolveIntersections();
    
    corkMesh2CorkTriMesh(&cmIn0, out);
}