#include "ofApp.h"
#include "CorkCsg.h"

//--------------------------------------------------------------
void ofApp::setup()
{
    ofSetFrameRate(60);
    ofBackground(0);
    
    //nm::CorkCsg::unifyVertices(ofMesh::cylinder(50, 50, 5, 5, 5, true, OF_PRIMITIVE_TRIANGLES), in0);
    
    nm::CorkCsg::unifyVertices(ofMesh::box(50, 50, 50), in0);
    
    ofMesh sphere = ofMesh::sphere(20, 20, OF_PRIMITIVE_TRIANGLES);
    for (auto& v : sphere.getVertices()) { v += glm::vec3(20, 20, 20); }
    nm::CorkCsg::unifyVertices(sphere, in1);
    
    nm::CorkCsg::computeUnion(in0, in1, outMesh);
}

//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw()
{
    ofSetWindowTitle(ofToString(ofGetFrameRate(), 2));
    
    cam.begin();
    
    outMesh.drawWireframe();
    
    cam.end();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
