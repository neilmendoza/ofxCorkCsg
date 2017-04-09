#include "ofApp.h"
#include "ofxCorkCsg.h"

//--------------------------------------------------------------
void ofApp::setup()
{
    ofSetFrameRate(60);
    ofBackground(255);
    
    // create box
    ofxCorkCsg::box(boxMesh, 150.f, 150.f, 150.f);
    
    // create sphere
    ofxCorkCsg::sphere(sphereMesh, 100.f);
    
    // if we wanted to translate the sphere we could do this
    // for (auto& v : sphereMesh.getVertices()) { v += glm::vec3(60.f, 60.f, 60.f); }
    
    // do union
    operation = UNION;
    operationChanged = true;
}

//--------------------------------------------------------------
void ofApp::update()
{
    if (operationChanged)
    {
        operate();
        operationChanged = false;
    }
}

//--------------------------------------------------------------
void ofApp::draw()
{
    ofSetWindowTitle(ofToString(ofGetFrameRate(), 2));
    
    cam.begin();
    ofPushStyle();
    ofSetColor(0);
    outMesh.drawWireframe();
    ofPopStyle();
    cam.end();
    
    ofPushStyle();
    for (unsigned i = 0; i < NUM_OPERATIONS; ++i)
    {
        if (operation == (Operation)i) ofSetColor(0);
        else ofSetColor(200);
        ostringstream oss;
        oss << i << ": " << toString((Operation)i);
        ofDrawBitmapString(oss.str(), 10.f, 20.f + i * 20.f);
    }
    ofPopStyle();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key)
{
    if (key >= '0' && key <= '3')
    {
        operation = (Operation)(key - '0');
        operationChanged = true;
    }
}

void ofApp::operate()
{
    switch (operation)
    {
        case UNION:
            ofxCorkCsg::computeUnion(boxMesh, sphereMesh, outMesh);
            break;
            
        case DIFFERENCE:
            ofxCorkCsg::computeDifference(boxMesh, sphereMesh, outMesh);
            break;
            
        case INTERSECTION:
            ofxCorkCsg::computeIntersection(boxMesh, sphereMesh, outMesh);
            break;
            
        case SYMETRIC_DIFFERENCE:
            ofxCorkCsg::computeSymmetricDifference(boxMesh, sphereMesh, outMesh);
            break;
            
        default:
            break;
    }
}

//--------------------------------------------------------------
string ofApp::toString(Operation operation) const
{
    switch (operation)
    {
        case UNION:
            return "UNION";
            break;
            
        case DIFFERENCE:
            return "DIFFERENCE";
            break;
            
        case INTERSECTION:
            return "INTERSECTION";
            break;
            
        case SYMETRIC_DIFFERENCE:
            return "SYMETRIC_DIFFERENCE";
            break;
            
        default:
            return "";
            break;
    }
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
