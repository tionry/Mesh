#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include <GL\glut.h>
#include <stdlib.h>
#include <iostream>
#include "shape.h"

using namespace std;
shape demo;

GLfloat xRotAngle=0.0f;
GLfloat yRotAngle=0.0f;
int mode = 0;
int norm_idx = -1;

void SetupRC(){
    glEnable(GL_COLOR_MATERIAL);  
    glEnable(GL_DEPTH_TEST);

    glDepthMask(GL_TRUE);

    glEnable(GL_LIGHTING); 
    GLfloat ambientLight[] = { 0.3f, 0.3f, 0.3f, 1.0f };
    GLfloat diffuseLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat specularLight[] = { 0.3f, 0.3f, 0.3f, 1.0f };
    GLfloat position1[] = {-5.0 ,5., 1., 0.};

    // Assign created components to GL_LIGHT0
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, position1);

    glShadeModel(GL_SMOOTH);

    glEnable(GL_LIGHT0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(90.0f, 1.0f, 1.0f, 500.0f);
}

void Reshape(int w, int h)
{
    if(h == 0)	h = 1;

    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    GLfloat fAspect;
    fAspect = (float)w/(float)h;
    gluPerspective(45.0, fAspect, 1.0, 500.0);


    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void renderObject(){

    glRotatef(xRotAngle,1.0f,0.0f,0.0f);

    glRotatef(yRotAngle,0.0f,1.0f,0.0f);
    
    demo.Show(mode);
    
	glFlush(); 
}

void Display(void)
{
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    

    glTranslatef(0.0f,0.0f,-3.0f);

    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH,GL_NICEST);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH,GL_NICEST);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_POLYGON_SMOOTH,GL_NICEST);
	
    renderObject();
    
    glPopMatrix();

    glutSwapBuffers();
}

void specialKey(int key,int x,int y){

    if(key == GLUT_KEY_UP){
        xRotAngle -= 8.0f;
    }
    else if(key == GLUT_KEY_DOWN){
        xRotAngle += 8.0f;
    }
    else if(key==GLUT_KEY_LEFT){
        yRotAngle -= 8.0f;
    }
    else if(key==GLUT_KEY_RIGHT){
        yRotAngle += 8.0f;
    }

    glutPostRedisplay();

}

void processMenu(int value){
    
    int* result_v;
    int n_v;
    int* result_f;
    int n_f;

    int n = 30;
    int face_region_ids[10] = {100,5425,7470,26256,6409,7481,30404,37880,1604,26664};
    int *vertex_ids;
	mode = value;
    
    glutPostRedisplay();
}

int main(int argc, char* argv[])
{
    // load obj
    char objfile[] = "./281.obj";
	char offfile[] = "./281.off";
    demo.LoadObjShape(objfile);
	//demo.LoadOffShape(offfile);

	demo.noise_flag = true;
	demo.denoise_flag = true;
	demo.stat_flag = false;
	demo.mse_flag = false;
	demo.average_point_length = 0;
    //original.LoadObjShape(objfile);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1200,800);
    glutCreateWindow("Mesh");
    glutReshapeFunc(Reshape);
    glutDisplayFunc(Display);

    SetupRC();

    // special keys
    glutSpecialFunc(specialKey);

    // menu
    int nMainMenu=glutCreateMenu(processMenu);
    glutAddMenuEntry("Shade Model Continuous Color", 1);
	glutAddMenuEntry("Shade Model Discrete Color", 2);
    glutAddMenuEntry("Show Neighborhood Vertexes", 3);
    glutAddMenuEntry("Show Neighborhood Faces", 4);
    glutAddMenuEntry("Show Neighborhood Faces for Face", 5);
    glutAddMenuEntry("Shade Region", 6);
    glutAddMenuEntry("Show Normal Vector", 7);
    glutAddMenuEntry("Add Noise", 8);
    glutAddMenuEntry("Denoise", 9);
	glutAddMenuEntry("CalculateMSE",10);

    glutAttachMenu(GLUT_RIGHT_BUTTON);

    glutMainLoop();
    
    return 0;
}