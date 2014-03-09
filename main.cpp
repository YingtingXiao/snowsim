#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "fluidsim.h"
#include "array3_utils.h"
#include "GL/glew.h"
#include "gluvi.h"
#include "stb_image_write.h"
#include "obj.h"
#include "objloader.h"

using namespace std;

//Simulation stuff
int grid_resolution = 100;
float timestep = 0.01f;
int frame = 0;
float grid_width = 1;
FluidSim sim;
string outpath;

//Display stuff
bool filming = false;
bool running = true;
std::vector<Vec3f> particle_positions;
float particle_radius;
string frame_number="frame 0";
obj* solidBoundary = new obj();

float sphere_phi(const Vec3f& position, const Vec3f& centre, float radius) {
    return (dist(position,centre) - radius);
}

Vec3f c0(0.5f,0.5f,0.5f);
float rad0 = 0.35f;

float boundary_phi(const Vec3f& position) {
    return -sphere_phi(position, c0, rad0);
}

float liquid_phi(const Vec3f& position) {
    return sphere_phi(position, Vec3f(0.55f, 0.55f, 0.4f), 0.23f);
}

void simulate_frame(int frame);
void export_particles(string path, int frame, const std::vector<Vec3f>& particles, float radius);

//Display fluid
void set_view(Gluvi::Target3D &cam)
{
    cam.target[0] = 0.5;
    cam.target[1] = 0.5;
    cam.target[2] = 0.5;
    cam.dist = 2;
}

void set_lights_and_material(int object)
{
    //Enable blending
    glEnable (GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDisable(GL_LIGHTING);
    glShadeModel(GL_FLAT);

    glColor4f(.4, .4, .7, .3);
}

void timer(int value)
{
    frame++;
    if(filming) {
        simulate_frame(frame);
        if(frame == 0) {
            filming = false;
        }
        glutPostRedisplay();
        glutTimerFunc(200, timer, 0);
    }

    if(running) {
        simulate_frame(frame);
        glutTimerFunc(1, timer, 0);
        glutPostRedisplay();
    }

}

//Draw solid boundary from mesh
void drawSolidBoundary() {
    glColor3f(0,0,0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glBegin(GL_TRIANGLES);
    for (int i=0; i<solidBoundary->faces.size(); ++i) {
        vector<int> face = solidBoundary->faces[i];
        glm::vec3 p0 = glm::vec3(solidBoundary->points[face[0]]);
        for (int j=1; j<face.size()-1; j++) {
            glm::vec3 p1 = glm::vec3(solidBoundary->points[face[j]]);
            glm::vec3 p2 = glm::vec3(solidBoundary->points[face[j+1]]);
            glVertex3f(p0.x, p0.y, p0.z);
            glVertex3f(p1.x, p1.y, p1.z);
            glVertex3f(p2.x, p2.y, p2.z);
        }
    }
    glEnd();
}

void display(void)
{
    glClearColor(0.5f, 0.5f, 0.5f, 1);

    //Coordinate system
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    glColor3f(1,0,0); glVertex3f(0,0,0); glVertex3f(0.1,0,0);
    glColor3f(0,1,0); glVertex3f(0,0,0); glVertex3f(0,0.1,0);
    glColor3f(0,0,1); glVertex3f(0,0,0); glVertex3f(0,0,0.1);
    glEnd();

    glEnable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    set_lights_and_material(1); 

    //Draw the liquid particles as simple spheres for now.
    glPolygonMode(GL_FRONT, GL_LINES);
    GLUquadric* particle_sphere;
    particle_sphere = gluNewQuadric();
    gluQuadricDrawStyle(particle_sphere, GLU_FILL );
    for(unsigned int p = 0; p < sim.particles.size(); ++p) {
        glPushMatrix();
        Vec3f pos = sim.particles[p].v;
        glTranslatef(pos[0], pos[1], pos[2]);
        gluSphere(particle_sphere, sim.particle_radius, 20, 20);
        glPopMatrix();   
    }

    //Draw the bound box for good measure
    glDisable(GL_LIGHTING);
    glColor3f(0,0,0);
    glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(0,0,1);

    glVertex3f(0,0,0);
    glVertex3f(0,1,0);

    glVertex3f(0,0,0);
    glVertex3f(1,0,0);

    glVertex3f(0,1,0);
    glVertex3f(1,1,0);

    glVertex3f(1,1,0);
    glVertex3f(1,0,0);

    glVertex3f(1,0,0);
    glVertex3f(1,0,1);

    glVertex3f(0,1,0);
    glVertex3f(0,1,1);

    glVertex3f(1,1,0);
    glVertex3f(1,1,1);

    glVertex3f(0,1,1);
    glVertex3f(1,1,1);

    glVertex3f(1,0,1);
    glVertex3f(1,1,1);

    glVertex3f(0,0,1);
    glVertex3f(1,0,1);

    glVertex3f(0,0,1);
    glVertex3f(0,1,1);

    glEnd();

    //Draw wireframe sphere geometry (specific to this scene).
    drawSolidBoundary();
}

void keyPress(unsigned char key, int x, int y) {
    if(key == 'r') {
        if(!filming) {
            if(!running) {
                running = true;
                glutTimerFunc(200, timer, 0);
            }
            else {
                running = false;
            }
        }
    }
    else if(key == 'f') {
        if(!running) {
            if(!filming) {
                filming = true;
            }
            else {
                filming = false;
            }
        }
    }
    glutPostRedisplay();
}

//Set up view fluid window
void setupDisplay(int argc, char **argv) {
    Gluvi::init("Liquid Data Viewer", &argc, argv);
    glutKeyboardFunc(keyPress);

    Gluvi::Target3D cam;
    set_view(cam);
    Gluvi::camera=&cam;

    Gluvi::userDisplayFunc=display;

    //Gluvi::StaticText frametext(frame_number.c_str());
    //Gluvi::root.list.push_back(&frametext);

    Gluvi::run();
}


//Main testing code
//-------------
int main(int argc, char **argv)
{
    if(argc != 3){
    	cerr << "Not enough arguments" << endl;
    	return 1;
    }

    string solid_file(argv[1]);
    string liquid_file(argv[2]);

    printf("Initializing data\n");
    sim.initialize(grid_width, grid_resolution, grid_resolution, grid_resolution);

    printf("Initializing boundary\n");
    sim.set_boundary(solid_file);

    printf("Initializing liquid\n");
    sim.set_liquid(liquid_file);

    printf("Exporting initial data\n");
    export_particles(outpath, 0, sim.particles, sim.particle_radius);

    //Load obj for drawing boundary
    objLoader* loader = new objLoader(solid_file, solidBoundary);
    solidBoundary->buildVBOs();
    delete loader;

    setupDisplay(argc, argv);

    return 0;
}

void simulate_frame(int frame) {
    printf("--------------------\nFrame %d\n", frame);

    //Simulate
    printf("Simulating liquid\n");
    sim.advance(timestep);

    //printf("Exporting particle data\n");
    //export_particles(outpath, frame, sim.particles, sim.particle_radius);

    //Save an image
    int recordWidth = 720;
    int recordHeight = 480;
    unsigned char* bitmapData = new unsigned char[3 * recordWidth * recordHeight];

    for (int i=0; i<recordHeight; i++) 
    {
        glReadPixels(0,i,recordWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
            bitmapData + (recordWidth * 3 * ((recordHeight-1)-i)));
    }

    char anim_filename[2048];
    sprintf_s(anim_filename, 2048, "output/fluid_%03d.png", frame); 

    stbi_write_png(anim_filename, recordWidth, recordHeight, 3, bitmapData, recordWidth * 3);

    delete [] bitmapData;
}

void export_particles(string path, int frame, const std::vector<Vec3f>& particles, float radius) {
    //Write the output
    std::stringstream strout;
    strout << path << "particles_" << frame << ".txt";
    string filepath = strout.str();

    ofstream outfile(filepath.c_str());
    //write vertex count and particle radius
    outfile << particles.size() << " " << radius << std::endl;
    //write vertices
    for(unsigned int i = 0; i < particles.size(); ++i)
        outfile << particles[i][0] << " " << particles[i][1] << " " << particles[i][2] << std::endl;
    outfile.close();
}


