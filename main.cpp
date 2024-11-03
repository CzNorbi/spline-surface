#define _USE_MATH_DEFINES

#include <iostream>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <unordered_map>
#include <cmath>

std::unordered_map<int, long long> memo;

struct Point {
	float x, y, z;
};

Point** controlPoints = nullptr;
Point** bezierPoints = nullptr;

static int rows, cols;
static float size = 8.0;
static float axisSize = 10.0;
static float pointSize = 5.0;
static const float resolution = 0.1;
static int numPoints = int(1 / resolution) + 1;


static bool showCoordinateSystem = true;
static bool showControlPoints = true;
static bool showControlHull = true;
static bool showSplineSurface = true;
static bool fill = false;

static float angle = 0;

Point** createPointMatrx(int rows, int cols) {
	Point** matrix = new Point * [rows];
	for (int i = 0; i < rows; ++i) {
		matrix[i] = new Point[cols];
	}
	return matrix;
}

void initialieControlPoints(Point** controlPoints, int row, int cols) {
	for (int i = 0; i < rows; ++i) {
		controlPoints[i][0] = Point{ 0.0, 0.0, size * i / (rows - 1) };
		controlPoints[i][cols - 1] = Point{ size, 0.0, size * i / (rows - 1) };
	}

	for (int i = 0; i < cols; ++i) {
		controlPoints[0][i] = { size * i / (cols - 1), 0.0, 0.0 };
		controlPoints[rows - 1][i] = { size * i / (cols - 1), 0.0, size };
	}

	for (int i = 1; i < rows - 1; ++i) {
		for (int j = 1; j < cols - 1; ++j) {
			controlPoints[i][j] = { size * j / (cols - 1), 1.0, size * i / (rows - 1) };
		}
	}
}


void printMatrix() {
	std::cout << "Matrix of Points (" << rows << "x" << cols << "):" << std::endl;
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			std::cout << "(" << controlPoints[i][j].x << ", "
				<< controlPoints[i][j].y << ", "
				<< controlPoints[i][j].z << ") ";
		}
		std::cout << std::endl;
	}
}

void drawControlPoints(Point** controlPoints, int rows, int cols, float pointSize = 5.0) {
	glPointSize(pointSize);
	glBegin(GL_POINTS);
	glColor3f(0.0, 0.0, 0.0);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			glVertex3f(controlPoints[i][j].x, controlPoints[i][j].y, controlPoints[i][j].z);
		}
	}
	glEnd();
}

void drawControlMesh(Point** controlPoints, int rows, int cols) {
	for (int i = 0; i < rows; ++i) {
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < cols; ++j) {
			glVertex3f(controlPoints[i][j].x, controlPoints[i][j].y, controlPoints[i][j].z);
		}
		glEnd();
	}

	for (int j = 0; j < cols; ++j) {
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < rows; ++i) {
			glVertex3f(controlPoints[i][j].x, controlPoints[i][j].y, controlPoints[i][j].z);
		}
		glEnd();
	}
}

long long factorial(int n) {
	if (memo.find(n) != memo.end()) {
		return memo[n];
	}

	if (n == 0 || n == 1) {
		return 1;
	}

	memo[n] = n * factorial(n - 1);
	return memo[n];
}

long double NCR(int n, int i) {
	return (long double) factorial(n) / (factorial(i) * (factorial(n - i)));
}

float bernstein(int n, int i, float t) {
	return NCR(n, i) * pow(t, i) * pow(1.0f - t, n - i);
}


void generateBezierSurface(Point** bezierPoints, int brows, int bcols, Point** controlPoints, int crows, int ccols, float resolution) {
	for (float u = 0; u < 1 + resolution; u += resolution)
	{
		for (float v = 0; v < 1 + resolution; v += resolution)
		{
			Point nextPoint = Point{ 0.0, 0.0, 0.0 };
			for (int i = 0; i < crows; i++) {
				GLfloat B1 = bernstein(crows - 1, i, u);
				for (int j = 0; j < ccols; j++) {
					GLfloat B2 = bernstein(ccols - 1, j, v);
					nextPoint.x += B1 * B2 * controlPoints[i][j].x;
					nextPoint.y += B1 * B2 * controlPoints[i][j].y;
					nextPoint.z += B1 * B2 * controlPoints[i][j].z;
				}
			}
			bezierPoints[int(u * (brows - 1))][int(v * (bcols - 1))] = Point{ nextPoint.x, nextPoint.y, nextPoint.z };
		}
	}
}

void drawBezierSurfave(Point** bezierPoints, int rows, int cols, bool fill)
{
	if (fill) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}

	glColor3f(0.0, 0.0, 0.0);
	for (int i = 0; i < rows - 1; i++) {
		glBegin(GL_TRIANGLE_STRIP);
		for (int j = 0; j < cols; j++) {
			glVertex3f(bezierPoints[i][j].x, bezierPoints[i][j].y, bezierPoints[i][j].z);
			glVertex3f(bezierPoints[i + 1][j].x, bezierPoints[i + 1][j].y, bezierPoints[i + 1][j].z);
		}
		glEnd();
	}
}

void drawCoordinateSystem() {
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(axisSize, 0.0, 0.0);

	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, axisSize, 0.0);

	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, axisSize);
	glEnd();
}


// Drawing routine.
void drawScene(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	glTranslatef(-2.0, -3.0, -20.0);

	// coordinate system
	if (showCoordinateSystem) {
		glPushMatrix();
		glTranslatef(-2.0, 0.0, 0.0);
		glTranslatef(axisSize / 2, 0.0, axisSize / 2);
		glRotatef(angle, 0.0, 1.0, 0.0);
		glTranslatef(-axisSize / 2, 0.0, -axisSize / 2);
		drawCoordinateSystem();
		glPopMatrix();
	}

	// control points and surface
	glPushMatrix();
	glTranslatef(-1.0, 1.0, 0.0);
	glTranslatef(size / 2, 0.0, size / 2);
	glRotatef(angle, 0.0, 1.0, 0.0);
	glTranslatef(-size / 2, 0.0, -size / 2);
	if (showControlPoints) {
		drawControlPoints(controlPoints, rows, cols);
	}
	if (showControlHull) {
		drawControlMesh(controlPoints, rows, cols);
	}
	if (showSplineSurface) {
		generateBezierSurface(bezierPoints, numPoints, numPoints, controlPoints, rows, cols, resolution);
		drawBezierSurfave(bezierPoints, numPoints, numPoints, fill);
	}
	glPopMatrix();

	glFlush();
}

// Initialization routine.
void setup(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
}

// OpenGL window reshape routine.
void resize(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (float)w / (float)h, 1.0, 50.0);
	glMatrixMode(GL_MODELVIEW);
}

// Keyboard input processing routine.
void keyInput(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 27:
		exit(0);
		break;
	case '0':
		showCoordinateSystem = !showCoordinateSystem;
		glutPostRedisplay();
		break;
	case '1':
		showControlPoints = !showControlPoints;
		glutPostRedisplay();
		break;
	case '2':
		showControlHull = !showControlHull;
		glutPostRedisplay();
		break;
	case '3':
		showSplineSurface = !showSplineSurface;
		glutPostRedisplay();
		break;
	case 'f':
		fill = !fill;
		glutPostRedisplay();
		break;
	case 'y': 
		angle -= 5.0;
		glutPostRedisplay();
		break;
	case 'Y':
		angle += 5.0;
		glutPostRedisplay();
		break;
	default:
		break;
	}
}

void printUserManual() {
	std::cout << "Press '0' to hide/show coordinate system." << std::endl;
	std::cout << "Press '1' to hide/show control points." << std::endl;
	std::cout << "Press '2' to hide/show control hull." << std::endl;
	std::cout << "Press '3' to hide/show spline surface." << std::endl;
	std::cout << "Press 'f' to switch between wired and filled surface." << std::endl;
	std::cout << "Press 'y'/'Y' to rotate around 'Y' axis." << std::endl;
}

// Main routine.
int main(int argc, char** argv)
{
	if (argc != 3) {
		std::cerr << "Usage: " << argv[0] << " <n> <m>" << std::endl;
		return 1;
	}

	rows = std::atoi(argv[1]);
	cols = std::atoi(argv[2]);

	if (rows < 3 || rows > 10 || cols < 3 || cols > 10) {
		std::cerr << "Invalid input arguments. The number of rows and columns must be in the [3, 10] range."<< std::endl;
	}

	controlPoints = createPointMatrx(rows, cols);
	initialieControlPoints(controlPoints, rows, cols);
	// printMatrix();
	bezierPoints = createPointMatrx(numPoints, numPoints);

	printUserManual();

	glutInit(&argc, argv);

	glutInitContextVersion(4, 3);
	glutInitContextProfile(GLUT_COMPATIBILITY_PROFILE);

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);

	glutInitWindowSize(500, 500);
	glutInitWindowPosition(100, 100);

	glutCreateWindow("spline-surface.cpp");

	glutDisplayFunc(drawScene);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keyInput);

	glewExperimental = GL_TRUE;
	glewInit();

	setup();

	glutMainLoop();
}
