#include <mex.h>
#include "GLFW/glfw3.h"

GLFWwindow *windows[100];
mwSize numWindows;
int keyPressed = -1;
int keyReleased = -1;
int modifiers = -1;
int buttonPressed = -1;
int buttonReleased = -1;
int activeWindow = -1;

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS) {
        keyPressed = key;
        modifiers = mods;
    }
    else if (action == GLFW_RELEASE) {
        keyReleased = key;
        modifiers = mods;
    }
}

static void mouse_callback(GLFWwindow* window, int button, int action, int mods)
{
    int i;
    
    if (action == GLFW_PRESS) {
        buttonPressed = button;
        modifiers = mods;
    }
    if (action == GLFW_RELEASE) {
        buttonReleased = button;
        modifiers = mods;
    }
    for (i = 0; i < numWindows; i++) {
        if (windows[i] == window) {
            activeWindow = i;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int command;
    int dummy;
    int width, height;
    int xpos, ypos;
    int antialiasing;
    double aspectRatio;
    double *windowInfo;
    const GLFWvidmode *mode;
    GLdouble *surfaceVertices;
    GLuint *surfaceIndices;
    GLdouble *surfaceColors;
    GLdouble *lineVertices;
    GLuint *lineIndices;
    GLdouble *lineColors;
    mwSize colorSize;
    double *currentKey, *currentKeyReleased, *currentButton, *currentButtonReleased, *currentModifiers, *cursorPosition, *currentWindow;
    double *background;
    double *colorSize3;
    int i;
    int wind, transformation, numVertices,numTriangles;
    double isMac;
    double *data;
    int dims[2];
    
    command = mxGetScalar(prhs[0]);
    
    // Initialize window
    if (command == 0) {
        // Register OpenGL termination to occur on Matlab exit
        mexAtExit(glfwTerminate);
        
        // Create new OpenGL window
        dummy = glfwInit();
        
        // Read in windows information
        windowInfo = mxGetPr(prhs[1]);
        numWindows = mxGetN(prhs[1]);
        
        // Do some things differently if this is a Mac.
        isMac = mxGetScalar(prhs[2]);
        
        for (i = 0; i < numWindows; i++) {
            // Create new windows
            // Set antialiasing
            antialiasing = windowInfo[5*i+4];
            glfwWindowHint(GLFW_SAMPLES, antialiasing);
            
            glfwWindowHint(GLFW_DECORATED, GL_FALSE);
            width = windowInfo[5*i+2];
            height = windowInfo[5*i+3];
            windows[i] = glfwCreateWindow(width, height, "ViRMEn", NULL, NULL);
            
            glfwMakeContextCurrent(windows[i]);
            glfwGetFramebufferSize(windows[i], &width, &height);
            glfwSwapInterval(1);
            glViewport(0, 0, width, height);
            
            xpos = windowInfo[5*i];
            ypos = windowInfo[5*i+1];
            glfwSetWindowPos(windows[i], xpos, ypos);
            
            // Callbacks for keyboard press and mouse clicks
            glfwSetKeyCallback(windows[i], key_callback);
            glfwSetMouseButtonCallback(windows[i], mouse_callback);
            
            // Initialize OpenGL properties
            aspectRatio = (double)width / (double)height;
            glOrtho(-aspectRatio, aspectRatio, -1, 1, -1000, 0);  // orthographic projection
            glEnable(GL_DEPTH_TEST);  // enable depth (for object occlusion)
            glClearDepth(-1.0);
            glDepthFunc(GL_GEQUAL);
            glShadeModel(GL_FLAT);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            
            // Enable the use of arrays to specify coordinates and colors
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);
            
            if (!isMac) {
                glfwIconifyWindow(windows[i]);
                glfwRestoreWindow(windows[i]);
            }
        }
        
        keyPressed = -1;
        keyReleased = -1;
        modifiers = -1;
        buttonPressed = -1;
        buttonReleased = -1;
        activeWindow = -1;
    }
    
    // Render
    else if (command == 1) {
        // Get surface arrays from Matlab
        surfaceVertices = (GLdouble *)mxGetData(prhs[1]);
        surfaceIndices = (GLuint *)mxGetData(prhs[2]);
        surfaceColors = (GLdouble *)mxGetData(prhs[3]);
        
        // Get line arrays from Matlab
        lineVertices = (GLdouble *)mxGetData(prhs[4]);
        lineIndices = (GLuint *)mxGetData(prhs[5]);
        lineColors = (GLdouble *)mxGetData(prhs[6]);
        
        wind = mxGetScalar(prhs[7]);
        transformation = mxGetScalar(prhs[8]);
        numVertices = mxGetScalar(prhs[9]);
        numTriangles = mxGetScalar(prhs[10]);
        
        glfwMakeContextCurrent(windows[wind-1]);
        
        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Determine size of the surface color matrix
        colorSize = mxGetM(prhs[3]);
        
        // Point OpenGL to the surface arrays
        glColorPointer(colorSize, GL_DOUBLE, 0, surfaceColors);
        glVertexPointer(3, GL_DOUBLE, 0, &(surfaceVertices[(transformation-1)*numVertices]));
        
        // Render the surface
        glDrawElements(GL_TRIANGLES,numTriangles, GL_UNSIGNED_INT, &(surfaceIndices[(transformation-1)*numTriangles]));
        
        // Determine size of the line color matrix
        colorSize = mxGetM(prhs[6]);
        
        if (colorSize > 0) { // only if any lines exist
            // Point OpenGL to the line arrays
            glColorPointer(3, GL_DOUBLE, 0, lineColors);
            glVertexPointer(2, GL_DOUBLE, 0, lineVertices);
            
            // Render the lines
            glDrawElements(GL_LINES, mxGetNumberOfElements(prhs[5]), GL_UNSIGNED_INT, lineIndices);
        }
        
        // Display and flush all events
        glFlush();
        glfwSwapBuffers(windows[wind-1]);
        glfwPollEvents();
        
        // Return current pressed key
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        currentKey = mxGetPr(plhs[0]);
        
        currentKey[0] = keyPressed;
        for (i = 0; i < numWindows; i++) {
            if (glfwWindowShouldClose(windows[i])) {
                currentKey[0] = 256;
            }
        }
        keyPressed = -1;
        
        // Return current clicked mouse button
        plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        currentKeyReleased = mxGetPr(plhs[1]);
        currentKeyReleased[0] = keyReleased;
        keyReleased = -1;
        
        // Return current clicked mouse button
        plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
        currentButton = mxGetPr(plhs[2]);
        currentButton[0] = buttonPressed;
        buttonPressed = -1;
        
        // Return current clicked mouse button
        plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
        currentButtonReleased = mxGetPr(plhs[3]);
        currentButtonReleased[0] = buttonReleased;
        buttonReleased = -1;
        
        // Return modifier keys that are currently pressed
        plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
        currentModifiers = mxGetPr(plhs[4]);
        currentModifiers[0] = modifiers;
        modifiers = -1;
        
        // Return current window
        plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
        currentWindow = mxGetPr(plhs[5]);
        currentWindow[0] = activeWindow;
        activeWindow = -1;
        
        // Return cursor position
        plhs[6] = mxCreateDoubleMatrix(1, 2, mxREAL);
        cursorPosition = mxGetPr(plhs[6]);
        glfwGetCursorPos(windows[wind-1], &(cursorPosition[0]), &(cursorPosition[1]));
    }
    
    // Terminate window
    else if (command == 2) {
        // Destroy window
        for (i = 0; i < numWindows; i++) {
            glfwDestroyWindow(windows[i]);
        }
        
        // Terminate GLFW
        glfwTerminate();
    }
    
    // Change transparency
    else if (command == 3) {
        // Get color matrix size (3 or 4) from Matlab
        colorSize3 = mxGetPr(prhs[1]);
        for (i = 0; i < numWindows; i++) {
            glfwMakeContextCurrent(windows[i]);
            if (colorSize3[0] == 3) {
                glDisable(GL_BLEND);
            }
            if (colorSize3[0] == 4) {
                glEnable(GL_BLEND);
            }
        }
        glFlush();
        glfwPollEvents();
    }
    
    // Change background color
    else if (command == 4) {
        // Get background color
        background = mxGetPr(prhs[1]);
        for (i = 0; i < numWindows; i++) {
            glfwMakeContextCurrent(windows[i]);
            glClearColor(background[0], background[1], background[2], 0.0);
        }
        glFlush();
        glfwPollEvents();
    }
    
    // Get pixel data
    else if (command == 5) {
        wind = mxGetScalar(prhs[1]);
        glfwMakeContextCurrent(windows[wind-1]);
        glfwGetFramebufferSize(windows[wind-1], &width, &height);
        
        // allocate space to store the pixels
        dims[0] = 3;
        dims[1] = width;
        dims[2] = height;
        plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
        
        data = mxGetPr(plhs[0]);
        
        glReadPixels(0, 0, width, height, GL_RGB, GL_FLOAT, data);
      
    }
   
}