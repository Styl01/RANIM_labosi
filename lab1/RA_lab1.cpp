#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <GL/freeglut.h>
#include <math.h>
#include <malloc.h>
#include <windows.h>
#include <vector>
#include <GLM/glm/vec3.hpp>
#include <GLM/glm/gtc/type_ptr.hpp>
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

Ociste	ociste = { 0, 4 ,12 };
double daljina_kamere = ociste.z;
double LOD = 20;
void myDisplay();
void myReshape(int width, int height);
//void myObject();
void Gen_curve(std::vector<glm::vec3>& points, B2spline spline);
void gen_deriv1(std::vector<glm::vec3>& derivs, B2spline spline);
void gen_deriv2(std::vector<glm::vec3>& derivs, B2spline spline);

void Draw_curve(/*B2spline spline,*/ std::vector<glm::vec3> points);
void Draw_points(B2spline spline);
void myKeyboard(unsigned char theKey, int mouseX, int mouseY);
void Draw_axis();
void Load_object(char* path, std::vector < glm::vec3 >& vertices);
void Draw_object(std::vector<glm::vec3> verticies);
std::vector<glm::vec3> offset_object(std::vector<glm::vec3>& verticies, glm::vec3 center, glm::vec3 offset, glm::vec3 deriv1, glm::vec3 deriv2);
void animate(int val);

void redisplay_all(void)
{
	glutSetWindow(window);
	myReshape(width, height);
	glutPostRedisplay();
}
void Load_object(char* path, std::vector < glm::vec3 >& vertices)
{
	//FILE* file;
	//errno_t err = fopen_s(&file,path, "r");
	FILE* file=fopen(path, "r");
	std::vector< glm::vec3 > temp_vertices;
	std::vector< UINT32> faces;
	glm::vec3 center = glm::vec3(0, 0, 0);
	while (1) {
		char buffer[64];
		int r = fscanf(file, "%s", buffer);
		if (r == EOF) break;
		else if (strcmp(buffer, "v") == 0) {
			glm::vec3 vertex;
			(void)fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
			temp_vertices.push_back(vertex);
			center += vertex;
		}
		else if (strcmp(buffer, "f") == 0) {
			int face[3];
			(void)fscanf(file, "%d %d %d\n", &face[0], &face[1], &face[2]);
			faces.push_back(face[0]);
			faces.push_back(face[1]);
			faces.push_back(face[2]);
		}

	
	}
	center.x/=1.0*temp_vertices.size();
	center.y /= 1.0 * temp_vertices.size();
	center.z /= 1.0 * temp_vertices.size();
	
	for (size_t i = 0; i < faces.size(); i++)
	{
		glm::vec3 vertex = temp_vertices[faces[i] - 1]-center;
		vertices.push_back(vertex);
	}
	printf("->->%d\n", vertices.size());
	//for (size_t i = 0; i < vertices.size(); i++)
	//{

	//}
	fclose(file);
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
void gen_deriv1(std::vector<glm::vec3>& derivs, B2spline spline) 
{
	glm::mat4 Bi = (1.f / 6)* glm::mat4({ {-1,3,-3,1},{3,-6,3,0},{-3,0,3,0},{1,4,1,0} });
	glm::vec4 T;
	for (size_t i = 0; i < spline.Pcnt - 3; i++)
	{
		glm::vec3 r1 = spline.points[i];
		glm::vec3 r2 = spline.points[i + 1];
		glm::vec3 r3 = spline.points[i + 2];
		glm::vec3 r4 = spline.points[i + 3];

		glm::mat4x3 R(r1,r2,r3,r4);
		for (size_t j = 0; j <LOD ; j++)
		{
			double t = j / LOD;
			T = glm::vec4(3 * t * t, 2 * t, 1, 0);
			glm::vec3 point = R * Bi * T;
			derivs.push_back(point);

		}
	}

}
void gen_deriv2(std::vector<glm::vec3>& derivs, B2spline spline)
{
	glm::mat4 Bi = (1.f / 6) * glm::mat4({ {-1,3,-3,1},{3,-6,3,0},{-3,0,3,0},{1,4,1,0} });
	glm::vec4 T;
	for (size_t i = 0; i < spline.Pcnt - 3; i++)
	{
		glm::vec3 r1 = spline.points[i];
		glm::vec3 r2 = spline.points[i + 1];
		glm::vec3 r3 = spline.points[i + 2];
		glm::vec3 r4 = spline.points[i + 3];

		glm::mat4x3 R(r1, r2, r3, r4);
		for (size_t j = 0; j < LOD; j++)
		{
			double t = j / LOD;
			T = glm::vec4(6 * t, 2 , 0, 0);
			glm::vec3 point = R*Bi*T;
			derivs.push_back(point);

		}
	}

}


std::vector < glm::vec3 > vertices;
std::vector < glm::vec3 > vertices2;
glm::vec3 center = glm::vec3(0.5, 0.5, 0.5);
B2spline spline;
std::vector < glm::vec3 > points;
std::vector < glm::vec3 > derivs1;
std::vector < glm::vec3 > derivs2;
unsigned int index = 1;

int main(int argc, char** argv)
{
	glutInit(&argc, argv);

	// postavljanje dvostrukog spremnika za prikaz (zbog titranja)
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(100, 100);


	char file[] = "lab1_materijali/kocka.obj";
	Load_object(file, vertices);
	spline.Pcnt = 10;
	double a[10][3] = { {-4,0.0,0.0},{0,0.5,4},{4,1,0},{0,1.5,-4}, {-4,2,0}, {0,2.5,4},{4,3,0},{0,3.5,-4},{-4,4,0}, {0,4.5,4} };
	for (size_t i = 0; i < spline.Pcnt; i++)
	{
		spline.points.push_back(glm::vec3(a[i][0], a[i][1], a[i][2]));
	
	}
	Gen_curve(points, spline);
	gen_deriv1(derivs1, spline);
	gen_deriv2(derivs2, spline);

	printf("veli�ine ->%d , %d, %d\n",points.size(), derivs1.size(), derivs2.size());
	//glm::vec3 init_offset;
	//init_offset = points[0]- glm::vec3(0.5,0.5,0.5);
	
	offset_object(vertices,center, points[0],derivs1[0],derivs2[0]);

	window = glutCreateWindow("Prozor");

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



void myReshape(int w, int h)
{
	width = w; height = h;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);        // aktivirana matrica projekcije
	glLoadIdentity();
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);		         // boja pozadine - bijela
	glClear(GL_COLOR_BUFFER_BIT);				     // brisanje zaslona
	gluPerspective(45.0, (float)width / height, 0.5, 30.0); // kut pogleda, x/y, prednja i straznja ravnina odsjecanja
	gluLookAt(ociste.x, ociste.y, ociste.z, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);	// ociste x,y,z; glediste x,y,z; up vektor x,y,z
	glColor3ub(0, 0, 255);
	glMatrixMode(GL_MODELVIEW);         // aktivirana matrica modela
}


void myDisplay(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	

	//printf("ide u crtanje\n");
	if (index > 1 && index<points.size()) {
		//glm::vec3 offset = points[index] - points[index - 1];
		vertices2=offset_object(vertices, center, points[index], derivs1[index],derivs2[index]);
	}
	Draw_object(vertices2);
	Draw_curve(/*spline,*/ points);
	Draw_points(spline);
	Draw_axis();
	glutSwapBuffers();      // iscrtavanje iz dvostrukog spemnika (umjesto glFlush)
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

void Draw_curve(/*B2spline spline,*/ std::vector<glm::vec3> points)
{
	glColor4f(1, 0, 0,0.5);
	//double LOD = 20;
	glLineWidth(5.0);
	glBegin(GL_LINE_STRIP);
	for (size_t i = 0; i < points.size(); i++)
	{
		glVertex3fv(glm::value_ptr(points[i]));
	}
	glEnd();

	glColor4f(0, 1, 0,1);
	glLineWidth(4.0);
	glBegin(GL_LINES);
	for (size_t i = 0; i < points.size(); i+=5)
	{
		glVertex3fv(glm::value_ptr(points[i]));
		glm::vec3 b = points[i] + 1.2f * glm::normalize(derivs1[i]);
		glVertex3fv(glm::value_ptr(b));
	}
	glEnd();


}
void Draw_points(B2spline spline)
{
	glColor3f(0, 1, 0);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	glVertex3d(0,0,0);
	for (size_t i = 0; i < spline.Pcnt; i++)
	{	
		//printf("....%f, %f\n", spline.points[i].x, spline.points[i].y);
		glVertex3d((GLdouble)spline.points[i].x, (GLdouble)spline.points[i].y, (GLdouble)spline.points[i].z);
	}
	glEnd();
}


void myKeyboard(unsigned char theKey, int mouseX, int mouseY)
{
	double phi = atan(ociste.x / ociste.z);
	double rho= acos(ociste.y/ociste.z);
	printf("kut gledanja -> %f\n", phi);
	printf("ociste->(%f,%f,%f)\n", ociste.x, ociste.y, ociste.z);
	switch (theKey)
	{
	case 'i':
		ociste.x *= 1.04;
		ociste.y *= 1.04;
		ociste.z *= 1.04;
		daljina_kamere *= 1.04;
		break;
	case 'p':
		ociste.x *= 1/1.04;
		ociste.y *= 1/1.04;
		ociste.z *= 1/1.04;
		daljina_kamere *=1/1.04;
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
	/*glBegin(GL_LINES);

	for (size_t i = 0; i < verticies.size(); i++)
	{
		glVertex3fv(glm::value_ptr(verticies[i]));

	}
	glEnd();*/
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
		glm::vec3 in =inv_rot*vertices[i] + offset;
		out.push_back(in);
	}
	return out;
}
void animate(int val) {
	
	if (index < points.size()) {
		index++;
	}
	glutPostRedisplay();
	glutTimerFunc(20, animate, 0);
}