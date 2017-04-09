#pragma once

#include "ofMain.h"

class ofApp : public ofBaseApp
{
public:
    enum Operation
    {
        UNION,
        DIFFERENCE,
        INTERSECTION,
        SYMETRIC_DIFFERENCE,
        
        NUM_OPERATIONS
    };
    
    void setup();
    void update();
    void draw();

    void keyPressed(int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y );
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void mouseEntered(int x, int y);
    void mouseExited(int x, int y);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);

private:
    void operate();
    string toString(Operation operation) const;
    
    ofVboMesh boxMesh, sphereMesh, outMesh;
    Operation operation;
    ofEasyCam cam;
    bool operationChanged;
};
