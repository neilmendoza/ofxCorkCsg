/*
 *  ofxCorkCsg.h
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
#pragma once

#include "MeshWrapper.h"

namespace ofxCorkCsg
{
    // Boolean operations follow
    // result = A U B
    void computeUnion(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh);
    void computeUnion(const MeshWrapper& in0, const MeshWrapper& in1, ofMesh& outMesh);
    
    // result = A - B
    void computeDifference(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh);
    void computeDifference(const MeshWrapper& in0, const MeshWrapper& in1, ofMesh& outMesh);
    
    // result = A ^ B
    void computeIntersection(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh);
    void computeIntersection(const MeshWrapper& in0, const MeshWrapper& in1, ofMesh& outMesh);
    
    // result = A XOR B
    void computeSymmetricDifference(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh);
    void computeSymmetricDifference(const MeshWrapper& in0, const MeshWrapper& in1, ofMesh& outMesh);
    
    // Not a Boolean operation, but related:
    //  No portion of either surface is deleted.  However, the
    //  curve of intersection between the two surfaces is made explicit,
    //  such that the two surfaces are now connected.
    void resolveIntersections(const ofMesh& in0, const ofMesh& in1, ofMesh& outMesh);
    void resolveIntersections(const MeshWrapper& in0, const MeshWrapper& in1, ofMesh& outMesh);
    
    // share vertices that are close enought to each other
    void unifyVertices(const ofMesh& inMesh, ofMesh& outMesh, float epsilonSq = 1e-8);
    
    // create ofMesh from CorkTriMesh
    void toOf(const CorkTriMesh& inMesh, ofMesh& outMesh);
    
    // helper functions to interate through indexed or non-indexed mesh
    unsigned getIndex(const ofMesh& mesh, const unsigned idx);
    unsigned getNumVertices(const ofMesh& mesh);
    
    // functions to create primitives with shared vertices
    void cylinder(ofMesh& cylinder,
                  float height,
                  float radius,
                  unsigned segments = 10,
                  unsigned verticalSlices = 4,
                  unsigned radialSlices = 2);
    
    // Probably don't use this, only works if there's no floating point errors in the vertices
    void fastUnifyVertices(const ofMesh& inMesh, ofMesh& outMesh);
}
