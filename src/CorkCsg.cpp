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
    void CorkCsg::computeUnion(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh)
    {
        CorkMeshWrapper corkIn0(in0);
        CorkMeshWrapper corkIn1(in1);
        computeUnion(corkIn0, corkIn1, outMesh);
    }
    
    void CorkCsg::computeUnion(const CorkMeshWrapper& in0, const CorkMeshWrapper& in1, ofMesh& outMesh)
    {
        CorkTriMesh corkOutMesh;
        computeUnion(in0.corkTriMesh, in1.corkTriMesh, &corkOutMesh);
        toOf(corkOutMesh, outMesh);
        freeCorkTriMesh(&corkOutMesh);
    }
    
    void CorkCsg::toOf(const CorkTriMesh& inMesh, ofMesh& outMesh)
    {
        outMesh.clear();
        outMesh.setMode(OF_PRIMITIVE_TRIANGLES);
        for (unsigned i = 0; i < inMesh.n_vertices; ++i)
        {
            glm::vec3 v;
            v.x = inMesh.vertices[i * 3];
            v.y = inMesh.vertices[i * 3 + 1];
            v.z = inMesh.vertices[i * 3 + 2];
            outMesh.addVertex(v);
        }
        for (unsigned i = 0; i < inMesh.n_triangles; ++i)
        {
            outMesh.addIndex(inMesh.triangles[i * 3]);
            outMesh.addIndex(inMesh.triangles[i * 3 + 1]);
            outMesh.addIndex(inMesh.triangles[i * 3 + 2]);
        }
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
    
    void CorkCsg::unifyVertices(const ofMesh& inMesh, ofMesh& outMesh, float epsilonSq)
    {
        outMesh.clear();
        if (inMesh.getMode() == OF_PRIMITIVE_TRIANGLES)
        {
            // keep track of how many vertices have been unified so we can divide out at the end
            vector<unsigned> divisors;
            
            for (unsigned i = 0; i < getNumVertices(inMesh); ++i)
            {
                unsigned oldIdx = getIndex(inMesh, i);
                auto& v = inMesh.getVertices()[oldIdx];
                unsigned idx = numeric_limits<unsigned>::max();
                for (unsigned i = 0; i < outMesh.getNumVertices(); ++i)
                {
                    if (glm::distance2(v, outMesh.getVertices()[i]) < epsilonSq)
                    {
                        // same vertex
                        idx = i;
                        outMesh.addIndex(idx);
                        divisors[idx]++;
                        
                        if (!inMesh.getNormals().empty()) outMesh.getNormals()[idx] += inMesh.getNormal(oldIdx);
                        if (!inMesh.getColors().empty()) outMesh.getColors()[idx] += inMesh.getColor(oldIdx);
                        if (!inMesh.getTexCoords().empty()) outMesh.getTexCoords()[idx] += inMesh.getTexCoord(oldIdx);
                        
                        break;
                    }
                }
                
                if (idx == numeric_limits<unsigned>::max())
                {
                    // didn't find vertex so add a new one
                    idx = outMesh.getNumVertices();
                    outMesh.addVertex(v);
                    outMesh.addIndex(idx);
                    divisors.push_back(1);
                    
                    if (!inMesh.getNormals().empty()) outMesh.addNormal(inMesh.getNormal(oldIdx));
                    if (!inMesh.getColors().empty()) outMesh.addColor(inMesh.getColor(oldIdx));
                    if (!inMesh.getTexCoords().empty()) outMesh.addTexCoord(inMesh.getTexCoord(oldIdx));
                }
            }
            
            // go through and take averages of non-vertex data
            for (unsigned i = 0; i < outMesh.getNumVertices(); ++i)
            {
                float divisor = 1.f / divisors[i];
                if (!inMesh.getNormals().empty()) outMesh.getNormals()[i] *= divisor;
                if (!inMesh.getColors().empty()) outMesh.getColors()[i] *= divisor;
                if (!inMesh.getTexCoords().empty()) outMesh.getTexCoords()[i] *= divisor;
            }
        }
        else ofLogError() << "unifyVertices only implemented for OF_PRIMITIVE_TRIANGLES";
    }
    
    void CorkCsg::fastUnifyVertices(const ofMesh& inMesh, ofMesh& outMesh)
    {
        if (inMesh.getMode() == OF_PRIMITIVE_TRIANGLES)
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
            
            for (unsigned i = 0; i < getNumVertices(inMesh); ++i)
            {
                unsigned oldIdx = getIndex(inMesh, i);
                auto& v = inMesh.getVertices()[oldIdx];
                auto it = vertexLookup.find(v);
                unsigned idx = numeric_limits<unsigned>::max();
                if (it == vertexLookup.end())
                {
                    idx = outMesh.getNumVertices();
                    outMesh.addVertex(v);
                    outMesh.addIndex(idx);
                    vertexLookup[v] = idx;
                    divisors.push_back(1);
                    
                    if (!inMesh.getNormals().empty()) outMesh.addNormal(inMesh.getNormal(oldIdx));
                    if (!inMesh.getColors().empty()) outMesh.addColor(inMesh.getColor(oldIdx));
                    if (!inMesh.getTexCoords().empty()) outMesh.addTexCoord(inMesh.getTexCoord(oldIdx));
                }
                else
                {
                    idx = vertexLookup[v];
                    outMesh.addIndex(idx);
                    divisors[idx]++;
                    
                    if (!inMesh.getNormals().empty()) outMesh.getNormals()[idx] += inMesh.getNormal(oldIdx);
                    if (!inMesh.getColors().empty()) outMesh.getColors()[idx] += inMesh.getColor(oldIdx);
                    if (!inMesh.getTexCoords().empty()) outMesh.getTexCoords()[idx] += inMesh.getTexCoord(oldIdx);
                }
            }
            
            // go through and take averages of non-vertex data
            for (unsigned i = 0; i < outMesh.getNumVertices(); ++i)
            {
                float divisor = 1.f / divisors[i];
                if (!inMesh.getNormals().empty()) outMesh.getNormals()[i] *= divisor;
                if (!inMesh.getColors().empty()) outMesh.getColors()[i] *= divisor;
                if (!inMesh.getTexCoords().empty()) outMesh.getTexCoords()[i] *= divisor;
            }
        }
        else ofLogError() << "unifyVertices only implemented for OF_PRIMITIVE_TRIANGLES";
    }
}
