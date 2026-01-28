#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <GL/glut.h>
#include <math.h>
#include <malloc.h>
#include <windows.h>
#include <vector>
#include <GLM/glm/vec3.hpp>
#include <GLM/glm/gtc/type_ptr.hpp>
#include <GLM/glm/gtc/random.hpp>
GLuint window;
GLuint width = 800, height = 500;

typedef struct _Ociste {
	GLdouble	x;
	GLdouble	y;
	GLdouble	z;
} Ociste;
typedef struct {
	double x;
	double y;
} Point2D;
typedef struct {
	float x;
	float y;
	float z;
}Point3D;
typedef struct {
	unsigned int Pcnt; // broj kontrolnih to�aka
	//double* knots; //vektor uzlova, pretpostavka da �e bit onak dobri pa mi ovo ne�e ni trebat
	std::vector<glm::vec3> points; //kontrolne to�ke
}B2spline;

Ociste	ociste = { 0, 2 ,4 };
double daljina_kamere = ociste.z;
double LOD = 20;
void myDisplay();
void myReshape(int width, int height);
//void myObject();
void Gen_curve(std::vector<glm::vec3>& points, B2spline spline);
void myKeyboard(unsigned char theKey, int mouseX, int mouseY);
void Draw_axis();
void Draw_ground();
void Load_object(char* path, std::vector < glm::vec3 >& vertices);
void Draw_object(std::vector<glm::vec3> verticies);
void Draw_wind_dir(void);
std::vector<glm::vec3> offset_object(std::vector<glm::vec3>& verticies, glm::vec3 center, glm::vec3 offset, glm::vec3 deriv1, glm::vec3 deriv2);
void animate(int val);

void redisplay_all(void)
{
	glutSetWindow(window);
	myReshape(width, height);
	glutPostRedisplay();
}



void Gen_curve(std::vector<glm::vec3>& points, B2spline spline)
{
	;
	for (size_t i = 0; i < spline.Pcnt - 3; i++)
	{
		glm::vec3 r1 = spline.points[i];
		glm::vec3 r2 = spline.points[i + 1];
		glm::vec3 r3 = spline.points[i + 2];
		glm::vec3 r4 = spline.points[i + 3];
		//printf("kontrolne tocke -> (%f,%f),(%f,%f),(%f,%f),(%f,%f)\n", r1.x, r1.y, r2.x, r2.y, r3.x, r3.y, r4.x, r4.y);
		for (size_t j = 0; j < LOD; j++)
		{
			double t = j / LOD;
			double B1 = 1.0 / 6 * (-1.0 * pow(t, 3) + 1.0) + 0.5 * (pow(t, 2) - t);
			double B2 = pow(t, 3) * 0.5 - pow(t, 2) + 2.0 / 3;
			double B3 = 0.5 * (-1 * pow(t, 3) + pow(t, 2) + t) + 1.0 / 6;
			double B4 = pow(t, 3) / 6;
			//printf("%f-%f,%f,%f,%f\n",t, B1, B2, B4, B4);
			GLdouble x = B1 * r1.x + B2 * r2.x + B3 * r3.x + B4 * r4.x;
			GLdouble y = B1 * r1.y + B2 * r2.y + B3 * r3.y + B4 * r4.y;
			GLdouble z = B1 * r1.z + B2 * r2.z + B3 * r3.z + B4 * r4.z;
			//printf("%f, %f, %f \n", x, y, z);
			//glVertex3d(x, y,z);
			glm::vec3 point;
			point.x = x; point.y = y; point.z = z;
			points.push_back(point);
		}
	}

}
typedef struct {
	glm::vec3 position;
	glm::vec3 velocity; 
	float mass;
	float size;
} particle;

particle  create_particle(glm::vec3 position, float size, bool zero_vel) {
	particle ret;
	ret.position = position;
	glm::vec3 vel;
	if (!zero_vel) {
		vel.x = 2 * glm::linearRand(-1.0, 1.0);
		vel.y = 1 * glm::linearRand(0.0, 1.0);
		vel.z = 2 * glm::linearRand(-1.0, 1.0);
	}
	else vel = glm::vec3(0.f, 0.f, 0.f);
	ret.velocity = vel;
	ret.mass = 1;
	ret.size = size;
	return ret;
}
particle create_initial_rain(glm::vec3 center, float a) {
	particle ret;
	glm::vec3 position;
	position.x= glm::linearRand(center.x-a/2, center.x + a/2);
	position.y = glm::linearRand(center.z, center.y + a / 2);
	position.z = glm::linearRand(center.z - a / 2, center.z + a / 2);
	ret.position = position;
	glm::vec3 vel= glm::vec3(0.f, 0.f, 0.f);
	ret.velocity = vel;
	ret.mass = 1;
	ret.size = 4;
	return ret;
}

std::vector < glm::vec3 > vertices;
std::vector < glm::vec3 > vertices2;
std::vector <particle> particles;
std::vector<particle> raindrops;
size_t rain_size = 500;
glm::vec3 gravity(0.f, -9.8f, 0.f);
float wind_strength = 3.5f;
glm::vec3 wind(wind_strength, 1.f, 0.f);
int dir = 1;
glm::vec3 spawn_point(0.f, 4.f, 0.f);
int anim_timer = 20;
float min_norm = 0.02;
float dt = 0.02;
size_t max_size = 1256;
GLuint texture;
unsigned int index = 1;

void particles_update(std::vector <particle> &particles, float dt, glm::vec3 gravity, glm::vec3 wind, glm::vec3 spawn_point) {
	
	//for (size_t i = 0; i < particles.size(); i++)
	for (std::vector <particle>::iterator i = particles.begin(); i!=particles.end();)
	{
		//if(particles[i].position.y<0.f)
		//{
		//	particles[i].velocity.y *= -1;
		//	particles[i].velocity *= 0.5;
		//}
		
		//if (glm::length(particles[i].velocity) < min_norm) particles[i] = create_particle(spawn_point,5.0,false);
		glm::vec3 force = (*i).mass * (gravity+wind);
		(*i).velocity += dt * force;
		(*i).position += (*i).velocity * dt;
		(*i).size -= 3 * dt;
		/*glm::vec3 force = particles[i].mass * gravity;
		particles[i].velocity += dt * force;
		particles[i].position += particles[i].velocity * dt;
		particles[i].size -= 4 * dt;*/
		//if (particles[i].size < 0.2 || particles[i].position.y < 0.f) {
		//	particles.erase(particles.begin() + i);
		//	//rem.push_back(i);
		//}
		if ((*i).size < 0.2 || (*i).position.y < 0.f) {
			i=particles.erase(i);
			
		}
		else {
			i++;
		}
		//i++;
		
	}

		
	
	return;
}
void rain_update(std::vector <particle>& particles, std::vector <particle>& particles2, float dt, glm::vec3 gravity, glm::vec3 wind, glm::vec3 center, float a) {
	for (size_t i = 0; i < particles.size(); i++)
	{
		if (particles[i].position.y < 0.f) {
			//printf("thing se dogodio ");
			for (size_t j = 0; j < 10+ i%4 ; j++){
				//printf("%d ", j);
				glm::vec3 pos(particles[i].position.x, 0, particles[i].position.z);
				particle droplet = create_particle(pos, 1, false);
				if (particles2.size() < max_size) {
					particles2.push_back(droplet);
				}
			}
			//printf("\n");

			glm::vec3 pos2;
			pos2.x= glm::linearRand(center.x - a / 2, center.x + a / 2);
			pos2.z = glm::linearRand(center.z - a / 2, center.z + a / 2);
			pos2.y = center.y;
			particles[i] = create_particle(pos2, particles[i].size, true);
		
			}
		

		else {
			glm::vec3 force = particles[i].mass * (gravity + wind);
			particles[i].velocity += dt * force;
			particles[i].position += particles[i].velocity * dt;
		}
	}
}

void Draw_points(std::vector <particle> particles)

{
	
	//glColor4f(0.3, 0.7, 1);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	//glVertex3d(0, 0, 0);
	printf("%d, crtanje \n",particles.size());
	for (size_t i = 0; i < particles.size(); i++)
	{
		//printf("%d ", i);
		//printf("....%f, %f\n", spline.points[i].x, spline.points[i].y);
		//printf("%f,%f,%f,%f,%f,%f\n", particles[i].position.x, particles[i].position.y, particles[i].position.z, particles[i].velocity.x, particles[i].velocity.y, particles[i].velocity.z);
		float p_size = 5 * particles[i].size;
		glColor4f(0.3, 0.7,1, particles[i].size);
		glPointSize(p_size);
		glVertex3fv(glm::value_ptr(particles[i].position));
	}
	glEnd();
}

void Draw_rain(std::vector <particle> raindrops) {
	glColor3f(0.3, 0.7, 1);
	glLineWidth(3);
	glBegin(GL_LINES);
	for (size_t i = 0; i < raindrops.size(); i++)
	{
		//printf("%f,%f,%f,%f,%f,%f\n", raindrops[i].position.x, raindrops[i].position.y, raindrops[i].position.z, raindrops[i].velocity.x, raindrops[i].velocity.y, raindrops[i].velocity.z);
		glVertex3fv(glm::value_ptr(raindrops[i].position));
		glm::vec3 b = raindrops[i].position - 1.5f*raindrops[i].velocity* dt;
		glVertex3fv(glm::value_ptr(b));
	}
	glEnd();
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);

	// postavljanje dvostrukog spremnika za prikaz (zbog titranja)
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(100, 100);

	//glEnable(GL_TEXTURE_2D);
	window = glutCreateWindow("Prozor");
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	//for (size_t i = 0; i < max_size; i++)
	//{
	//	particle par = create_particle(spawn_point);
	//	particles.push_back(par);

	//}
	//for (size_t i = 0; i < 4; i++)
	//{
	//	printf("%f,%f,%f,%f,%f,%f\n", particles[i].position.x, particles[i].position.y, particles[i].position.z, particles[i].velocity.x, particles[i].velocity.y, particles[i].velocity.z);
	//}
	//char filename[] = "lab2_materijali/cestica.bmp";
	//texture = load_texture(filename);

	
	//window = glutCreateWindow("Prozor");

	//ki�a, poku�aj 1

	//init po�etnu ki�u
	for (size_t i = 0; i < rain_size ; i++)
	{
		particle rain_drop = create_initial_rain(spawn_point, 4);
		raindrops.push_back(rain_drop);
	}
	//for (size_t i = 0; i <4; i++)
	//{
	//	printf("%f,%f,%f,%f,%f,%f\n", raindrops[i].position.x, raindrops[i].position.y, raindrops[i].position.z, raindrops[i].velocity.x, raindrops[i].velocity.y, raindrops[i].velocity.z);
	//}

	glutKeyboardFunc(myKeyboard);
	glutReshapeFunc(myReshape);
	glutDisplayFunc(myDisplay);
	redisplay_all();
	glutTimerFunc(anim_timer, animate, 0);
	glutMainLoop();

	return 0;
}



void myReshape(int w, int h)
{
	width = w; height = h;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);        // aktivirana matrica projekcije
	glLoadIdentity();
	glClearColor(0.8f, 0.8f, 0.8f, 1.0f);		         // boja pozadine - bijela
	glClear(GL_COLOR_BUFFER_BIT);				     // brisanje zaslona
	gluPerspective(45.0, (float)width / height, 0.5, 30.0); // kut pogleda, x/y, prednja i straznja ravnina odsjecanja
	gluLookAt(ociste.x, ociste.y, ociste.z, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);	// ociste x,y,z; glediste x,y,z; up vektor x,y,z
	glColor3ub(0, 0, 255);
	glMatrixMode(GL_MODELVIEW);         // aktivirana matrica modela
}



void myDisplay(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	rain_update(raindrops,particles, dt, gravity,wind, spawn_point,4);
	particles_update(particles, dt, gravity,wind, spawn_point);
	//Draw_object(vertices2);
	Draw_axis();
	
	Draw_ground();
	Draw_wind_dir();
	Draw_rain(raindrops);
	Draw_points(particles);
	glutSwapBuffers();      // iscrtavanje iz dvostrukog spemnika (umjesto glFlush)
}



void Draw_object(std::vector<glm::vec3> verticies) {
	glColor3f(0.8, 0.8, 0.8);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_TRIANGLES);

	for (size_t i = 0; i < verticies.size(); i++)
	{
		glVertex3fv(glm::value_ptr(verticies[i]));

	}
	glEnd();
	glColor3f(0.0, 0.0, 0.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_TRIANGLES);

	for (size_t i = 0; i < verticies.size(); i++)
	{
		glVertex3fv(glm::value_ptr(verticies[i]));

	}
	glEnd();
}
std::vector<glm::vec3> offset_object(std::vector<glm::vec3>& verticies, glm::vec3 center, glm::vec3 offset, glm::vec3 deriv1, glm::vec3 deriv2) {
	std::vector<glm::vec3> out;

	glm::vec3 norm = glm::cross(deriv1, deriv2);
	glm::vec3 binorm = glm::cross(deriv1, norm);
	glm::mat3x3 rotation(glm::normalize(deriv1), glm::normalize(norm), glm::normalize(binorm));
	rotation = glm::transpose(rotation);
	glm::mat3x3 inv_rot = glm::inverse(rotation);
	printf("test -> %f %f \n", deriv1[1], rotation[1][0]);
	for (size_t i = 0; i < vertices.size(); i++)
	{
		glm::vec3 in = inv_rot * vertices[i] + offset;
		out.push_back(in);
	}
	return out;
}
void animate(int val) {
	wind.x -= (float)dir * dt ;
	wind.z += (float)dir * dt;
	if (wind.z > wind_strength) dir *= -1;
	if (wind.x > wind_strength) dir *= -1;
	glutPostRedisplay();
	glutTimerFunc(anim_timer, animate, 0);
}
void Draw_axis()
{
	glLineWidth(2.0);
	glBegin(GL_LINES);

	glColor4f(1.0, 1.0, 0.0, 0.2);
	glVertex3f(-6.0, 0.0, 0.0);
	glVertex3f(6.0, 0.0, 0.0);

	glColor4f(1.0, 1.0, 0.0, 0.5);
	glVertex3f(0.0, -6.0, 0.0);
	glVertex3f(0.0, 6.0, 0.0);

	glColor4f(1.0, 1.0, 0.0, 0.5);
	glVertex3f(0.0, 0.0, -6.0);
	glVertex3f(0.0, 0.0, 6.0);
	glEnd();
}

void Draw_wind_dir() {
	glLineWidth(6.0);
	glBegin(GL_LINES);
	glColor4f(1.0, 0.0, 0.0, 0.3);
	glVertex3f(0.0, 0.f, 0.0);
	glVertex3fv(glm::value_ptr(0.5f* glm::vec3(wind.x,0.f,wind.z)));


	glEnd();


}
void Draw_ground() {

	glBegin(GL_QUADS);
	glColor4f(0.0, 1.0, 0.0, 1.0);
	glVertex3f(-6.0, 0.0, 6.0);
	glVertex3f(6.0, 0.0, 6.0);
	glVertex3f(6.0, 0.0, -6.0);
	glVertex3f(-6.0, 0.0, -6.0);
	glEnd();
	return;
}

void myKeyboard(unsigned char theKey, int mouseX, int mouseY)
{
	double phi = atan(ociste.x / ociste.z);
	double rho = acos(ociste.y / ociste.z);
	//printf("kut gledanja -> %f\n", phi);
	//printf("ociste->(%f,%f,%f)\n", ociste.x, ociste.y, ociste.z);
	switch (theKey)
	{
	case 'i':
		ociste.x *= 1.04;
		ociste.y *= 1.04;
		ociste.z *= 1.04;
		daljina_kamere *= 1.04;
		break;
	case 'p':
		ociste.x *= 1 / 1.04;
		ociste.y *= 1 / 1.04;
		ociste.z *= 1 / 1.04;
		daljina_kamere *= 1 / 1.04;
		break;

	case 'l':

		//ociste.x = ociste.x + 0.1;
		ociste.x = daljina_kamere * sin((phi + 0.03));
		ociste.z = daljina_kamere * cos((phi + 0.03));
		break;

	case 'k':
		//ociste.x = ociste.x - 0.1;
		ociste.x = daljina_kamere * sin((phi - 0.03));
		ociste.z = daljina_kamere * cos((phi - 0.03));
		break;

	case 'r': ociste.x = 0.0;
		break;
	case 'o': index = 0;
		break;

	case 27:  exit(0);
		break;
	}

	redisplay_all();
}