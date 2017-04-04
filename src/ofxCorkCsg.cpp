/*
 *  ofxCorkCsg.cpp
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
#include "ofxCorkCsg.h"

namespace ofxCorkCsg
{
    void computeUnion(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh)
    {
        computeUnion(MeshWrapper(in0), MeshWrapper(in1), outMesh);
    }
    
    void computeDifference(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh)
    {
        computeDifference(MeshWrapper(in0), MeshWrapper(in1), outMesh);
    }
    
    void computeIntersection(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh)
    {
        computeIntersection(MeshWrapper(in0), MeshWrapper(in1), outMesh);
    }
    
    void computeSymmetricDifference(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh)
    {
        computeSymmetricDifference(MeshWrapper(in0), MeshWrapper(in1), outMesh);
    }
    
    void resolveIntersections(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh)
    {
        resolveIntersections(MeshWrapper(in0), MeshWrapper(in1), outMesh);
    }
    
    void computeUnion(const MeshWrapper& in0, const MeshWrapper& in1, ofMesh& outMesh)
    {
        CorkTriMesh corkOutMesh;
        computeUnion(in0.corkTriMesh, in1.corkTriMesh, &corkOutMesh);
        toOf(corkOutMesh, outMesh);
        freeCorkTriMesh(&corkOutMesh);
    }
    
    void computeDifference(const MeshWrapper& in0, const MeshWrapper& in1, ofMesh& outMesh)
    {
        CorkTriMesh corkOutMesh;
        computeDifference(in0.corkTriMesh, in1.corkTriMesh, &corkOutMesh);
        toOf(corkOutMesh, outMesh);
        freeCorkTriMesh(&corkOutMesh);
    }
    
    void computeIntersection(const MeshWrapper& in0, const MeshWrapper& in1, ofMesh& outMesh)
    {
        CorkTriMesh corkOutMesh;
        computeIntersection(in0.corkTriMesh, in1.corkTriMesh, &corkOutMesh);
        toOf(corkOutMesh, outMesh);
        freeCorkTriMesh(&corkOutMesh);
    }
    
    void computeSymmetricDifference(const MeshWrapper& in0, const MeshWrapper& in1, ofMesh& outMesh)
    {
        CorkTriMesh corkOutMesh;
        computeSymmetricDifference(in0.corkTriMesh, in1.corkTriMesh, &corkOutMesh);
        toOf(corkOutMesh, outMesh);
        freeCorkTriMesh(&corkOutMesh);
    }
    
    void toOf(const CorkTriMesh& inMesh, ofMesh& outMesh)
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
    
    void toCork(const ofMesh& inMesh, MeshWrapper& outMesh)
    {
        outMesh.init(inMesh);
    }
    
    unsigned getIndex(const ofMesh& mesh, const unsigned idx)
    {
        if (mesh.getIndices().empty()) return idx;
        else return mesh.getIndices()[idx];
    }
    
    unsigned getNumVertices(const ofMesh& mesh)
    {
        if (mesh.getIndices().empty()) return mesh.getNumVertices();
        else return mesh.getNumIndices();
    }
    
    void unifyVertices(const ofMesh& inMesh, ofMesh& outMesh, float epsilonSq)
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
    
    void octohedron(ofMesh& mesh, float width, float height)
    {
        const glm::vec3 vertices[6] = {
            glm::vec3(0.f, height, 0.f),
            glm::vec3(width, 0.f, 0.f),
            glm::vec3(0.f, 0.f, width),
            glm::vec3(-width, 0.f, 0.f),
            glm::vec3(0.f, 0.f, -width),
            glm::vec3(0.f, -height, 0.f)
        };
        
        const unsigned indices[24] = {
            0, 2, 1,
            0, 3, 2,
            0, 4, 3,
            0, 1, 4,
            5, 1, 2,
            5, 2, 3,
            5, 3, 4,
            5, 4, 1
        };
        
        for (auto& v : vertices) mesh.addVertex(v);
        for (auto& i : indices) mesh.addIndex(i);
    }
    
    void box(ofMesh& mesh, float width, float height, float depth, int resX, int resY, int resZ)
    {
        unifyVertices(ofMesh::box(width, height, depth, resX, resY, resZ), mesh);
    }
    
    void sphere(ofMesh& mesh, float radius, unsigned resolution)
    {
        unifyVertices(ofMesh::sphere(radius, resolution, OF_PRIMITIVE_TRIANGLES), mesh);
        // wind the vertices correctly
        for (unsigned i = 0; i < mesh.getNumIndices() / 3; ++i)
        {
            unsigned i0 = mesh.getIndices()[i * 3];
            mesh.getIndices()[i * 3] = mesh.getIndices()[i * 3 + 1];
            mesh.getIndices()[i * 3 + 1] = i0;
        }
    }
    
    void cylinder(ofMesh& cylinder,
                  float height,
                  float radius,
                  unsigned segments,
                  unsigned verticalSlices,
                  unsigned radialSlices)
    {
        // slice in the at y = 0;
        ofMesh cap;
        for (unsigned sliceIdx = 0; sliceIdx < radialSlices; ++sliceIdx)
        {
            float r = radius * (sliceIdx + 1) / (float)radialSlices;
            for (unsigned segmentIdx = 0; segmentIdx < segments; ++segmentIdx)
            {
                float theta = TWO_PI * segmentIdx / (float)segments;
                cap.addVertex(glm::vec3(r * sin(theta), 0.f, r * cos(theta)));
            }
        }
        cap.addVertex(glm::vec3());
        
        // inner ring
        for (unsigned segmentIdx = 1; segmentIdx <= segments; ++segmentIdx)
        {
            cap.addIndex(cap.getNumVertices() - 1);
            cap.addIndex(segmentIdx - 1);
            cap.addIndex(segmentIdx % segments);
        }
        
        // outer rings
        for (unsigned sliceIdx = 1; sliceIdx < radialSlices; ++sliceIdx)
        {
            for (unsigned segmentIdx = 1; segmentIdx <= segments; ++ segmentIdx)
            {
                unsigned i0 = (sliceIdx - 1) * segments + segmentIdx - 1;
                unsigned i1 = sliceIdx * segments + segmentIdx - 1;
                unsigned i2 = sliceIdx * segments + (segmentIdx % segments);
                unsigned i3 = (sliceIdx - 1) * segments + (segmentIdx % segments);
                
                cap.addIndex(i0);
                cap.addIndex(i1);
                cap.addIndex(i3);
                
                cap.addIndex(i3);
                cap.addIndex(i1);
                cap.addIndex(i2);
            }
        }
        
        // body
        ofMesh body;
        float sectionHeight = height * (verticalSlices - 2) / (float)verticalSlices;
        for (unsigned sliceIdx = 0; sliceIdx < verticalSlices - 1; ++sliceIdx)
        {
            float y = ofMap(sliceIdx, 0, verticalSlices - 2, -.5f * sectionHeight, .5f * sectionHeight);
            for (unsigned segmentIdx = 0; segmentIdx < segments; ++segmentIdx)
            {
                float theta = TWO_PI * segmentIdx / (float)segments;
                body.addVertex(glm::vec3(radius * sin(theta), y, radius * cos(theta)));
            }
        }
        for (unsigned sliceIdx = 1; sliceIdx < verticalSlices - 1; ++sliceIdx)
        {
            for (unsigned segmentIdx = 1; segmentIdx <= segments; ++segmentIdx)
            {
                unsigned i0 = (sliceIdx - 1) * segments + segmentIdx - 1;
                unsigned i1 = sliceIdx * segments + segmentIdx - 1;
                unsigned i2 = sliceIdx * segments + (segmentIdx % segments);
                unsigned i3 = (sliceIdx - 1) * segments + (segmentIdx % segments);
                
                body.addIndex(i0);
                body.addIndex(i3);
                body.addIndex(i1);
                
                body.addIndex(i3);
                body.addIndex(i2);
                body.addIndex(i1);
            }
        }
        
        // add bottom cap with triangle indices reversed
        const unsigned bottomIdx = 0;
        glm::vec3 bottomOffset(0.f, -.5f * height, 0.f);
        for (auto& v : cap.getVertices()) { cylinder.addVertex(v + bottomOffset); }
        for (unsigned i = 0; i < cap.getNumIndices() / 3; ++i)
        {
            cylinder.addIndex(cap.getIndex(3 * i));
            cylinder.addIndex(cap.getIndex(3 * i + 2));
            cylinder.addIndex(cap.getIndex(3 * i + 1));
        }
        
        // add top cap
        const unsigned topIdx = cylinder.getNumVertices();
        glm::vec3 topOffset(0.f, .5f * height, 0.f);
        for (auto& v : cap.getVertices()) { cylinder.addVertex(v + topOffset); }
        for (auto& i : cap.getIndices()) cylinder.addIndex(i + topIdx);
        
        // add body
        const unsigned bodyIdx = cylinder.getNumVertices();
        cylinder.addVertices(body.getVertices());
        for (auto& i : body.getIndices()) { cylinder.addIndex(i + bodyIdx); }
        
        // stitch together
        const unsigned capBottomOuterIdx = segments * (radialSlices - 1);
        const unsigned bodyBottomOuterIdx = bodyIdx;
        const unsigned capTopOuterIdx = capBottomOuterIdx + topIdx;
        const unsigned bodyTopOuterIdx = bodyIdx + segments * (verticalSlices - 2);
        for (unsigned segmentIdx = 1; segmentIdx <= segments; ++segmentIdx)
        {
            unsigned i0 = capBottomOuterIdx + segmentIdx - 1;
            unsigned i1 = bodyBottomOuterIdx + segmentIdx - 1;
            unsigned i2 = bodyBottomOuterIdx + (segmentIdx % segments);
            unsigned i3 = capBottomOuterIdx + (segmentIdx % segments);
            
            cylinder.addIndex(i0);
            cylinder.addIndex(i3);
            cylinder.addIndex(i1);
            
            cylinder.addIndex(i3);
            cylinder.addIndex(i2);
            cylinder.addIndex(i1);
            
            i0 = bodyTopOuterIdx + segmentIdx - 1;
            i1 = capTopOuterIdx + segmentIdx - 1;
            i2 = capTopOuterIdx + (segmentIdx % segments);
            i3 = bodyTopOuterIdx + (segmentIdx % segments);
            
            cylinder.addIndex(i0);
            cylinder.addIndex(i3);
            cylinder.addIndex(i1);
            
            cylinder.addIndex(i3);
            cylinder.addIndex(i2);
            cylinder.addIndex(i1);
        }
    }
    
    void fastUnifyVertices(const ofMesh& inMesh, ofMesh& outMesh)
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
