/*
 *  CorkMeshWrapper.cpp
 *
 *  Copyright (c) 2017, Neil Mendoza, http://www.neilmendoza.com
 *  All rights reserved. 
 *  
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met: 
 *  
 *  * Redistributions of source code must retain the above copyright notice, 
 *    this list of conditions and the following disclaimer. 
 *  * Redistributions in binary form must reproduce the above copyright 
 *    notice, this list of conditions and the following disclaimer in the 
 *    documentation and/or other materials provided with the distribution. 
 *  * Neither the name of Neil Mendoza nor the names of its contributors may be used 
 *    to endorse or promote products derived from this software without 
 *    specific prior written permission. 
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE. 
 *
 */
#include "CorkMeshWrapper.h"

namespace nm
{
    CorkMeshWrapper::CorkMeshWrapper(const ofMesh& mesh)
    {
        corkTriMesh.n_triangles = mesh.getNumIndices() / 3;
        corkTriMesh.n_vertices = mesh.getNumVertices();
        corkTriMesh.triangles = new unsigned[corkTriMesh.n_triangles * 3];
        corkTriMesh.vertices = new float[corkTriMesh.n_vertices * 3];
        for (unsigned i = 0; i < mesh.getNumVertices(); ++i)
        {
            corkTriMesh.vertices[i * 3] = mesh.getVertices()[i].x;
            corkTriMesh.vertices[i * 3 + 1] = mesh.getVertices()[i].y;
            corkTriMesh.vertices[i * 3 + 2] = mesh.getVertices()[i].z;
        }
        for (unsigned i = 0; i < mesh.getNumIndices(); ++i)
        {
            corkTriMesh.triangles[i] = mesh.getIndices()[i];
        }
    }
    
    CorkMeshWrapper::~CorkMeshWrapper()
    {
        delete[] corkTriMesh.triangles;
        delete[] corkTriMesh.vertices;
    }
}
