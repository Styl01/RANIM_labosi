#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <math.h>
#include <malloc.h>
#include <windows.h>
#include <vector>
#include <GLM/glm/vec3.hpp>
#include <GLM/glm/gtc/type_ptr.hpp>

#define PI 3.1415926
GLuint window;
GLuint width = 800, height = 500;

typedef struct _Ociste {
	GLdouble	x;
	GLdouble	y;
	GLdouble	z;
} Ociste;
Ociste	ociste = { 0, 0 ,6 };

typedef struct _sine {
	float f;
	float A;
	float phase;

} Sine;

void myDisplay();
void myReshape(int width, int height);
void Draw_curve(std::vector<glm::vec2> points,  float r, float g, float b, float a);
void Draw_chain(std::vector<Sine> dtf_out, float t);
void Draw_circles(std::vector<glm::vec2> origin, std::vector<float> radius);

std::vector<glm::vec2> points;
std::vector<glm::vec2> points2;
std::vector<Sine> dft;
int camera_index = 0;
float camera_dist = 6;
int Band_width = 50;
float t = 0.0;
float dt = 0.001;

std::vector<glm::vec2> load_points(char* file_name) {
	//FILE* file;
	//errno_t err = fopen_s(&file, file_name, "r");
	FILE* file=fopen(file_name,"r");
	std::vector<glm::vec2> ret;
	while (1) {
		glm::vec2 buffer;
		int r = fscanf(file, "%f,%f", &buffer.x, &buffer.y);
		if (r == EOF) break;

		else {
			//buffer.y *= -1;
			ret.push_back(buffer);
		}
	}


	return ret;

}

std::vector<Sine> DFT(std::vector<glm::vec2> points, std::vector<int> freqs) {
	std::vector<Sine> ret;
	size_t N = points.size();
	//printf("%d,%d\n", N, freqs.size());
	for (size_t k = 0; k < freqs.size(); k++)
	{
		
		float Re = 0.0;
		float Im = 0.0;
		Sine temp;
		for (size_t i = 0; i < N; i++)
		{
			float a = cos(2 * PI * freqs[k]*i / N);
			float b = -1 * sin(2 * PI * freqs[k] * i / N);

			Re += (points[i].x * a - points[i].y * b)/N;
			Im += (points[i].y * a + points[i].x * b)/N;
			
		}
		temp.f = freqs[k];
		temp.A = sqrtf(powf(Re, 2) + powf(Im, 2));
		//temp.phase = atan(Re/Im);
		//printf("%f,%f -- %f,%f \n", atan(Re / Im), atan(Im / Re), Re, Im);

		temp.phase = atan2(Im , Re);
		ret.push_back(temp);
	}
	return ret;
}
void Draw_curve(std::vector<glm::vec2> points, float r , float g , float b, float a) {
	glColor4f(r,g,b,a);
	glLineWidth(5.0);
	glBegin(GL_LINE_STRIP);
	for (size_t i = 0; i < points.size(); i++)
	{
		glVertex2fv(glm::value_ptr(points[i]));
	}
	
	glEnd();

}
void Draw_chain(std::vector<Sine> dft_out, float t) {
	glColor4f(0,0,0,1.0);
	glLineWidth(2.0);
	glBegin(GL_LINE_STRIP);
	glm::vec2 end = glm::vec2(0, 0);
	std::vector<glm::vec2> origin;
	std::vector<float> radius;
	glVertex2fv(glm::value_ptr(end));
	for (size_t i = 0; i < dft_out.size(); i++)
	{
		origin.push_back(end);
		radius.push_back(dft_out[i].A);
		float phi = 2*PI*t * dft_out[i].f + dft_out[i].phase;
		glm::vec2 temp = glm::vec2(dft_out[i].A * cos(phi), dft_out[i].A *sin( phi));
		//if (dft_out[i].f == 0.0) {
		//	printf("%f,%f,:: %f, %f \n", dft_out[i].A, dft_out[i].phase, dft_out[i].A * cos(phi), dft_out[i].A * sin(phi));
		//}
		//Draw_circle(end, dft_out[i].A);
		if (i == camera_index) {
			ociste.x = end.x;
			ociste.y = end.y;
			ociste.z = camera_dist;
		}
		temp += end;
		end = temp;
		glVertex2fv(glm::value_ptr(temp));
	}
	glEnd();
	Draw_circles(origin, radius);
	points2.push_back(end);
}
void Draw_circles(std::vector<glm::vec2> origin, std::vector<float> radius) {
	int LOD = 100;
	glLineWidth(1.0);
	
	for(int k=0; k<origin.size(); k++){
		glBegin(GL_LINE_LOOP);

		for (size_t i = 0; i < LOD; i++)
		{
			float t = 1.0 * i / LOD;
			glm::vec2 point = glm::vec2(radius[k] * cos(2 * PI * t) + origin[k].x, radius[k] * sin(2 * PI * t) + origin[k].y);
			glVertex2fv(glm::value_ptr(point));
		}
		glEnd();
	}


}

void redisplay_all(void)
{
	glutSetWindow(window);
	myReshape(width, height);
	glutPostRedisplay();
}

void myReshape(int w, int h)
{
	width = w; height = h;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);        // aktivirana matrica projekcije
	glLoadIdentity();
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);		         // boja pozadine - bijela
	glClear(GL_COLOR_BUFFER_BIT);				     // brisanje zaslona
	gluPerspective(45.0, (float)width / height, 0.5, 100.0); // kut pogleda, x/y, prednja i straznja ravnina odsjecanja
	gluLookAt(ociste.x, ociste.y, ociste.z, ociste.x,ociste.y, 0.0, 0.0, 1.0, 0.0);	// ociste x,y,z; glediste x,y,z; up vektor x,y,z
	glColor3ub(0, 0, 255);
	glMatrixMode(GL_MODELVIEW);         // aktivirana matrica modela
}


void myDisplay(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	Draw_curve(points, 0.0,0.0,1.0,0.35);
	Draw_chain(dft, t);
	Draw_curve(points2, 1.0, 0.0, 0.0, 1.0);
	//glMatrixMode(GL_PROJECTION);        // aktivirana matrica projekcije
	//glLoadIdentity();
	//gluLookAt(ociste.x, ociste.y, ociste.z, ociste.x, ociste.y, 0.0, 0.0, 1.0, 0.0);
	//glMatrixMode(GL_MODELVIEW);
	glutSwapBuffers();      // iscrtavanje iz dvostrukog spemnika (umjesto glFlush)
}

void animate(int val) {
	t += dt;
	
	if (t > 1.0){ 
		t = 0.0;
		points2.clear();
	}
	glutPostRedisplay();
	myReshape(width, height);
	glutTimerFunc(5, animate, 0);
}
void myKeyboard(unsigned char theKey, int mouseX, int mouseY)
{

	//printf("kut gledanja -> %f\n", phi);
	//printf("ociste->(%f,%f,%f)\n", ociste.x, ociste.y, ociste.z);
	switch (theKey)
	{
	case 'i':
		camera_index--;
		camera_dist += 0.4;
		break;

	case 'p':
		camera_index++;
		camera_dist -=0.4;
		break;

	
	}
	printf("%f", camera_dist);
	redisplay_all();
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);

	// postavljanje dvostrukog spremnika za prikaz (zbog titranja)
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(100, 100);

	window = glutCreateWindow("Prozor");

	
	std::vector<int> freqs;

	for (int i = 0; i < Band_width+1; i++)
	{
		freqs.push_back(i);
		if (i != 0) freqs.push_back(-1*i);
	}
		
	char file[] = "curves/letterN.txt";
	printf(":%s\n", file);
	points = load_points(file);
	for (size_t i = 0; i < points.size(); i++)
	{
		//printf("%f,%f\n", points[i].x, points[i].y);
	}


	dft = DFT(points, freqs);
	char file2[] = "curves/out2.txt";
	//dft = test(file2);
	//printf("%d--\n", dft.size());
	for (size_t i = 0; i <dft.size(); i++)
	{
		//printf("%f,%f,%f\n", dft[i].f, dft[i].A, dft[i].phase);
	}
	//printf("done\n");

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glutKeyboardFunc(myKeyboard);
	glutReshapeFunc(myReshape);
	glutDisplayFunc(myDisplay);
	redisplay_all();
	glutTimerFunc(20, animate, 0);
	glutMainLoop();


	return 0;
}
