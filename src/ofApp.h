#pragma once

#include "ofMain.h"
#include "ofxLeapMotion2.h"
#include "threadedObject.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>


#define OF_MAX_FRAME_RATE   60

// include our ThreadedObject class.
class ofApp : public ofBaseApp{
    

public:
    // oF Drawing Functions
    void setup();
    void update();
    void draw();
    void keyPressed  (int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y );
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);
    void exit();
    // Leap
	ofxLeapMotion leap;
	vector <ofxLeapMotionSimpleHand> simpleHands;
	vector <int> fingersFound;
    //
    // Open FrameWorks Camera
    ofEasyCam cam;
    // Threaded Objects
    ThreadedObject threadedObject;
    
};
