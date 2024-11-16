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
static bool showBSplineSurface = false;
static bool fill = false;

static float angleY, angleX, angleZ = 0.0;

static float knot_u[7] = { 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0 };
static float knot_v[8] = { 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0 };

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

// -------------------------- Bezier ------------------------ //

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

// --------------------------- End of Bezier ---------------------------- //


// -------------------------- B-Spline ------------------------------------ //

// B-spline basis function definition
float basisFunction(int i, int degree, float t, const float* knots) {
	if (degree == 0) {
		return (t >= knots[i] && t < knots[i + 1]) ? 1.0f : 0.0f;
	}
	else {
		float denom1 = knots[i + degree] - knots[i];
		float denom2 = knots[i + degree + 1] - knots[i + 1];
		float term1 = (denom1 == 0) ? 0 : ((t - knots[i]) / denom1) * basisFunction(i, degree - 1, t, knots);
		float term2 = (denom2 == 0) ? 0 : ((knots[i + degree + 1] - t) / denom2) * basisFunction(i + 1, degree - 1, t, knots);
		return term1 + term2;
	}
}

// Compute a point on the B-spline surface
void computeBSplineSurface(float u, float v, int rows, int cols, Point& nextPoint) {
	if (u >= 1.0f) u = 0.999f;
	if (v >= 1.0f) v = 0.999f;

	//int num_u = 6, num_v = 4;
	int num_u = rows, num_v = cols;
	int degree_u = 2, degree_v = 2;

	nextPoint.x = nextPoint.y = nextPoint.z = 0.0f;
	for (int i = 0; i < num_u; ++i) {
		float bu = basisFunction(i, degree_u, u, knot_u);
		for (int j = 0; j < num_v; ++j) {
			float bv = basisFunction(j, degree_v, v, knot_v);
			nextPoint.x += bu * bv * controlPoints[i][j].x;
			nextPoint.y += bu * bv * controlPoints[i][j].y;
			nextPoint.z += bu * bv * controlPoints[i][j].z;

		}
	}
}

// Draw the surface by rendering line strips for a wireframe
void drawBSplineSurface() {
	float u, v, x, y, z;
	Point nextPoint;
	int num_samples = 10; // Number of samples along u and v directions

	// Wireframe along the u direction (rows)
	for (int j = 0; j <= num_samples; ++j) {
		v = (float)j / num_samples;

		// Start a new line strip for each row of constant v
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i <= num_samples; ++i) {
			u = (float)i / num_samples;
			computeBSplineSurface(u, v, rows, cols, nextPoint);
			glVertex3f(nextPoint.x, nextPoint.y, nextPoint.z);
		}
		glEnd();
	}

	// Wireframe along the v direction (columns)
	for (int i = 0; i <= num_samples; ++i) {
		u = (float)i / num_samples;

		// Start a new line strip for each column of constant u
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j <= num_samples; ++j) {
			v = (float)j / num_samples;
			computeBSplineSurface(u, v, rows, cols, nextPoint);
			glVertex3f(nextPoint.x, nextPoint.y, nextPoint.z);
		}
		glEnd();
	}
}

// --------------------------- End of B-Spline ---------------------------- //

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
		glTranslatef(-3.0, 0.0, 0.0);
		glTranslatef(axisSize / 2, 0.0, axisSize / 3);
		glRotatef(angleX, 1.0, 0.0, 0.0);
		glRotatef(angleY, 0.0, 1.0, 0.0);
		glRotatef(angleZ, 0.0, 0.0, 1.0);
		glTranslatef(-axisSize / 2, 0.0, -axisSize / 3);
		drawCoordinateSystem();
		glPopMatrix();
	}

	// control points and surface
	glPushMatrix();
	glTranslatef(-1.0, 1.0, 0.0);
	glTranslatef(size / 2, 0.0, size / 2);
	glRotatef(angleX, 1.0, 0.0, 0.0);
	glRotatef(angleY, 0.0, 1.0, 0.0);
	glRotatef(angleZ, 0.0, 0.0, 1.0);
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
	if (showBSplineSurface)
	{
		drawBSplineSurface();
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
		showBSplineSurface = false;
		showSplineSurface = !showSplineSurface;
		glutPostRedisplay();
		break;
	case '4':
		showSplineSurface = false;
		showBSplineSurface = !showBSplineSurface;
		glutPostRedisplay();
		break;
	case 'f':
		fill = !fill;
		glutPostRedisplay();
		break;
	case 'y': 
		angleY -= 5.0;
		glutPostRedisplay();
		break;
	case 'Y':
		angleY += 5.0;
		glutPostRedisplay();
		break;
	case 'x': 
		angleX -= 5.0;
		glutPostRedisplay();
		break;
	case 'X':
		angleX += 5.0;
		glutPostRedisplay();
		break;
	case 'c': 
		angleZ -= 5.0;
		glutPostRedisplay();
		break;
	case 'C':
		angleZ += 5.0;
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
	std::cout << "Press '3' to hide/show Bezier surface." << std::endl;
	std::cout << "Press '4' to hide/show B-spline surface." << std::endl;
	std::cout << "Press 'f' to switch between wired and filled surface." << std::endl;
	std::cout << "Press 'x'/'X' to rotate around 'X' axis." << std::endl;
	std::cout << "Press 'y'/'Y' to rotate around 'Y' axis." << std::endl;
	std::cout << "Press 'c'/'C' to rotate around 'Z' axis." << std::endl;
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
	glutInitWindowPosition(700, 100);

	glutCreateWindow("spline-surface.cpp");

	glutDisplayFunc(drawScene);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keyInput);

	glewExperimental = GL_TRUE;
	glewInit();

	setup();

	glutMainLoop();
}
