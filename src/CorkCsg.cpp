/*
 *  CorkCsg.cpp
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
#include "CorkCsg.h"

namespace nm
{
    ofMesh CorkCsg::computeUnion(const ofMesh& in0, const ofMesh& in1)
    {
        CorkMeshWrapper corkIn0(in0);
        CorkMeshWrapper corkIn1(in1);
        CorkTriMesh* out;
        computeUnion(corkIn0.corkTriMesh, corkIn1.corkTriMesh, out);
        freeCorkTriMesh(out);
    }
    
    bool CorkCsg::isSolid(const ofMesh& mesh)
    {
        CorkMeshWrapper corkMesh(mesh);
        return isSolid(corkMesh.corkTriMesh);
    }
    
    CorkMeshWrapper CorkCsg::toCork(const ofMesh& mesh)
    {
        return CorkMeshWrapper(mesh);
    }
    
    unsigned CorkCsg::getIndex(const ofMesh& mesh, const unsigned idx)
    {
        if (mesh.getIndices().empty()) return idx;
        else return mesh.getIndices()[idx];
    }
    
    unsigned CorkCsg::getNumVertices(const ofMesh& mesh)
    {
        if (mesh.getIndices().empty()) return mesh.getNumVertices();
        else return mesh.getNumIndices();
    }
    
    ofMesh CorkCsg::unifyVertices(const ofMesh& mesh, float epsilonSq)
    {
        // this is the unified mesh
        ofMesh unified;
        
        if (mesh.getMode() == OF_PRIMITIVE_TRIANGLES)
        {
            // keep track of how many vertices have been unified so we can divide out at the end
            vector<unsigned> divisors;
            
            for (unsigned i = 0; i < getNumVertices(mesh); ++i)
            {
                unsigned oldIdx = getIndex(mesh, i);
                auto& v = mesh.getVertices()[oldIdx];
                unsigned idx = numeric_limits<unsigned>::max();
                for (unsigned i = 0; i < unified.getNumVertices(); ++i)
                {
                    if (glm::distance2(v, unified.getVertices()[i]) < epsilonSq)
                    {
                        // same vertex
                        idx = i;
                        unified.addIndex(idx);
                        divisors[idx]++;
                        
                        if (!mesh.getNormals().empty()) unified.getNormals()[idx] += mesh.getNormal(oldIdx);
                        if (!mesh.getColors().empty()) unified.getColors()[idx] += mesh.getColor(oldIdx);
                        if (!mesh.getTexCoords().empty()) unified.getTexCoords()[idx] += mesh.getTexCoord(oldIdx);
                        
                        break;
                    }
                }
                
                if (idx == numeric_limits<unsigned>::max())
                {
                    // didn't find vertex so add a new one
                    idx = unified.getNumVertices();
                    unified.addVertex(v);
                    unified.addIndex(idx);
                    divisors.push_back(1);
                    
                    if (!mesh.getNormals().empty()) unified.addNormal(mesh.getNormal(oldIdx));
                    if (!mesh.getColors().empty()) unified.addColor(mesh.getColor(oldIdx));
                    if (!mesh.getTexCoords().empty()) unified.addTexCoord(mesh.getTexCoord(oldIdx));
                }
            }
            
            // go through and take averages of non-vertex data
            for (unsigned i = 0; i < unified.getNumVertices(); ++i)
            {
                float divisor = 1.f / divisors[i];
                if (!mesh.getNormals().empty()) unified.getNormals()[i] *= divisor;
                if (!mesh.getColors().empty()) unified.getColors()[i] *= divisor;
                if (!mesh.getTexCoords().empty()) unified.getTexCoords()[i] *= divisor;
            }
        }
        else ofLogError() << "unifyVertices only implemented for OF_PRIMITIVE_TRIANGLES";
        
        return unified;
    }
    
    ofMesh CorkCsg::fastUnifyVertices(const ofMesh& mesh)
    {
        // this is the unified mesh
        ofMesh unified;
        
        if (mesh.getMode() == OF_PRIMITIVE_TRIANGLES)
        {
            auto compare = [](glm::vec3 a, glm::vec3 b)
            {
                if (a.x == b.x)
                {
                    if (a.y == b.y) return a.z < b.z;
                    else return a.y < b.y;
                }
                else return a.x < b.x;
            };
            map<glm::vec3, unsigned, decltype(compare)> vertexLookup(compare);
            
            // keep track of how many vertices have been unified so we can divide out at the end
            vector<unsigned> divisors;
            
            for (unsigned i = 0; i < getNumVertices(mesh); ++i)
            {
                unsigned oldIdx = getIndex(mesh, i);
                auto& v = mesh.getVertices()[oldIdx];
                auto it = vertexLookup.find(v);
                unsigned idx = numeric_limits<unsigned>::max();
                if (it == vertexLookup.end())
                {
                    idx = unified.getNumVertices();
                    unified.addVertex(v);
                    unified.addIndex(idx);
                    vertexLookup[v] = idx;
                    divisors.push_back(1);
                    
                    if (!mesh.getNormals().empty()) unified.addNormal(mesh.getNormal(oldIdx));
                    if (!mesh.getColors().empty()) unified.addColor(mesh.getColor(oldIdx));
                    if (!mesh.getTexCoords().empty()) unified.addTexCoord(mesh.getTexCoord(oldIdx));
                }
                else
                {
                    idx = vertexLookup[v];
                    unified.addIndex(idx);
                    divisors[idx]++;
                    
                    if (!mesh.getNormals().empty()) unified.getNormals()[idx] += mesh.getNormal(oldIdx);
                    if (!mesh.getColors().empty()) unified.getColors()[idx] += mesh.getColor(oldIdx);
                    if (!mesh.getTexCoords().empty()) unified.getTexCoords()[idx] += mesh.getTexCoord(oldIdx);
                }
            }
            
            // go through and take averages of non-vertex data
            for (unsigned i = 0; i < unified.getNumVertices(); ++i)
            {
                float divisor = 1.f / divisors[i];
                if (!mesh.getNormals().empty()) unified.getNormals()[i] *= divisor;
                if (!mesh.getColors().empty()) unified.getColors()[i] *= divisor;
                if (!mesh.getTexCoords().empty()) unified.getTexCoords()[i] *= divisor;
            }
        }
        else ofLogError() << "unifyVertices only implemented for OF_PRIMITIVE_TRIANGLES";
        
        return unified;
    }
}
