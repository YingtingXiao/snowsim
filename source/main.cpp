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

#include "Eigen/Dense"

using namespace std;

//Simulation stuff
int grid_resolution = 0;
float timestep = 0.001f;
int frame = 0;
float grid_width = 1;
FluidSim sim;
string outpath;

//Display stuff
bool filming = false;
bool running = true;
std::vector<Vec3f> particle_positions;
float particle_radius;
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
void save_image(int frame);
void export_grid(int frame);
void export_mesh(int frame);
void export_particles(string path, int frame, const std::vector<Particle>& particles);

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
   //---Draw translucent spheres---//
   //Enable blending
   glEnable (GL_BLEND);
   glEnable(GL_ALPHA_TEST);
   glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   glEnable(GL_DEPTH_TEST);
   glDepthFunc(GL_LEQUAL);
   glDisable(GL_LIGHTING);
   glShadeModel(GL_FLAT);

   glColor4f(.8, .8, .85, .3);

   ////---Draw blinn spheres---//
   //glEnable(GL_LIGHTING);
   //GLfloat global_ambient[4] = {0.1f, 0.1f, 0.1f, 1.0f};
   //glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
   //glShadeModel(GL_SMOOTH);

   ////Light #1
   //GLfloat color[4] = {1.0f, 1.0f, 1.0f, 1.0f};
   //GLfloat position[3] = {1.0f, 1.0f, 1.0f};
   //glLightfv(GL_LIGHT0, GL_SPECULAR, color);
   //glLightfv(GL_LIGHT0, GL_DIFFUSE, color);
   //glLightfv(GL_LIGHT0, GL_POSITION, position);

   ////Light #2
   //GLfloat color2[4] = {1.0f, 1.0f, 1.0f, 1.0f};
   //GLfloat position2[3] = {-1.0f, -1.0f, 1.0f};
   //glLightfv(GL_LIGHT1, GL_SPECULAR, color2);
   //glLightfv(GL_LIGHT1, GL_DIFFUSE, color2);
   //glLightfv(GL_LIGHT1, GL_POSITION, position2);

   //GLfloat obj_color[4] = {.2, .3, .7};
   //glMaterialfv (GL_FRONT, GL_AMBIENT, obj_color);
   //glMaterialfv (GL_FRONT, GL_DIFFUSE, obj_color);

   //GLfloat specular[4] = {.4, .2, .8};
   //glMaterialf (GL_FRONT, GL_SHININESS, 32);
   //glMaterialfv (GL_FRONT, GL_SPECULAR, specular);
   //glEnable(GL_LIGHT0);
   //glEnable(GL_LIGHT1);
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
      Vec3f pos = sim.particles[p].pos;
      glTranslatef(pos[0], pos[1], pos[2]);
      gluSphere(particle_sphere, sim.particles[p].radius/*0.01*/, 20, 20);
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

   Gluvi::run();
}


//Main testing code
//-------------
int main(int argc, char **argv)
{
   if(argc != 4){
      cerr << "Not enough arguments" << endl;
      return 1;
   }

   string solid_file(argv[1]);
   string liquid_file(argv[2]);
   grid_resolution = atoi(argv[3]);

   printf("Initializing data\n");
   sim.initialize(grid_width, grid_resolution, grid_resolution, grid_resolution);

   printf("Initializing boundary\n");
   sim.set_boundary(solid_file);

   printf("Initializing liquid\n");
   sim.set_liquid(liquid_file);

   //printf("Exporting initial data\n");
   //export_particles(outpath, 0, sim.particles);

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
   sim.advance(timestep, frame==1);

   //printf("Exporting particle data\n");
   //export_particles(outpath, frame, sim.particles, sim.particle_radius);
   
   //if (frame % 10 == 1) {
      int f = frame;
      save_image(f);
      //export_grid(f);
      export_mesh(f);
   //}

   //save_image(frame);
}

void save_image(int frame) {
   int recordWidth = 1080;
   int recordHeight = 720;
   unsigned char* bitmapData = new unsigned char[3 * recordWidth * recordHeight];

   for (int i=0; i<recordHeight; i++) 
   {
      glReadPixels(0,i,recordWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
         bitmapData + (recordWidth * 3 * ((recordHeight-1)-i)));
   }

   char anim_filename[2048];
   sprintf_s(anim_filename, 2048, "../output/image/snow_%04d.png", frame); 

   stbi_write_png(anim_filename, recordWidth, recordHeight, 3, bitmapData, recordWidth * 3);

   delete [] bitmapData;

   cout << "Exported image to " << anim_filename << endl;
}

void export_grid(int frame) {
   //Write the output
   char filepath[2048];
   sprintf_s(filepath, 2048, "../output/data/data_%04d.txt", frame);

   ofstream outfile(filepath);

   //Write grid resolution
   outfile << grid_resolution << endl;
   outfile << endl;

   for(int i=0; i<sim.mass.a.size(); ++i) {
      outfile << sim.mass.a[i] << endl;
   }

   cout << "Exported grid data to " << filepath << endl;
}

void export_mesh(int frame) {
   std::vector<Vec3f> position, normal;
	std::vector<unsigned int> indices;
   sim.marching_cube(position, normal, indices);
   std::cout<<"mesh point size: "<<position.size()<<" indices size: "<<indices.size()<<endl;
   
   char filepath[2048];
   sprintf_s(filepath, 2048, "../output/mesh/mesh_%04d.obj", frame);
   std::ofstream outfile(filepath);
   outfile<< "g vertex"<<endl;
   for (int i = 0; i<position.size(); i++) {
	    outfile<< "v " <<position[i][0]<<" "<<position[i][1]<<" "<<position[i][2]<<endl;
   }
   
   outfile<<"g normal"<<endl;
   for (int j = 0;j<normal.size();j++)
	   outfile<< "vn "<<normal[j][0]<<" "<<normal[j][1]<<" "<<normal[j][2]<<endl;
   outfile<<"g faces"<<endl;
   int size = indices.size();
   for(int k = 0; k<size;k+=3)
		outfile<<"f " <<indices[k]+1<<" "<<indices[k+1]+1<<" "<<indices[k+2]+1<<endl;
   outfile.close();

   cout << "Exported mesh to " << filepath << endl;
}

void export_particles(string path, int frame, const std::vector<Particle>& particles) {
   //Write the output
   std::stringstream strout;
   strout << path << "particles_" << frame << ".txt";
   string filepath = strout.str();

   ofstream outfile(filepath.c_str());
   //write vertices
   for(unsigned int i = 0; i < particles.size(); ++i)
      outfile << particles[i].pos[0] << " " << particles[i].pos[1] << " " << particles[i].pos[2] << std::endl;
   outfile.close();
}


