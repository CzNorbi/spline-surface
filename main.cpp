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
Point** surfacePoints = nullptr;

static int rows, cols;
static float size = 8.0;
static float axisSize = 10.0;
static float pointSize = 5.0;
static const float resolution = 0.1;
static int numSurfaceRows = int(1 / resolution) + 1;
static int numSurfaceCols = int(1 / resolution) + 1;
static int selectedRow = 0, selectedColumn = 0;


static bool showCoordinateSystem = true;
static bool showControlPoints = true;
static bool showControlHull = true;
static bool showBezierSurface = true;
static bool showBSplineSurface = false;
static bool showNurbsSurface = false;
static bool fill = false;

static GLfloat material_black[] = { 0.0f, 0.0f, 0.0f, 1.0f };
static GLfloat material_pink[] = { 1.0f, 0.0f, 1.0f, 1.0f };
static GLfloat material_red[] = { 1.0f, 0.0f, 0.0f, 1.0f };
static GLfloat material_green[] = { 0.0f, 1.0f, 0.0f, 1.0f };
static GLfloat material_blue[] = { 0.0f, 0.0f, 1.0f, 1.0f };

static float angleY, angleX, angleZ = 0.0;

// ---- B-Spline ---- //
static int uLength, vLength;

float* knot_u = nullptr;
float* knot_v = nullptr;

// ---------------- //

void initiateKnotVector(float* knotVectorU, int ulength, float numberOfControlPointsU) {
	if (ulength < 3) {
		std::cerr << "Error: Knot vector length too small!" << std::endl;
		return;
	}
	
	for (int i = 0; i < ulength; i++) {
		knotVectorU[i] = 0.0;
	}

	for (int i = ulength - 3; i < ulength; i++) {
		knotVectorU[i] = 1.0;
	}

	bool still = true;
	int j = 3, index = 1;
	while (still) {
		if (knotVectorU[j] == 1.0) {
			still = false;
		}
		else {
			knotVectorU[j] = floorf((index / (float)(numberOfControlPointsU - 2)) * 100) / 100;
			j++;
			index++;
		}
	}
}

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
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material_black);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_black);
	glPointSize(pointSize);

	glBegin(GL_POINTS);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			glVertex3f(controlPoints[i][j].x, controlPoints[i][j].y, controlPoints[i][j].z);
		}
	}
	glEnd();

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material_pink);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_pink);

	glBegin(GL_POINTS);
	glVertex3f(controlPoints[selectedRow][selectedColumn].x, controlPoints[selectedRow][selectedColumn].y, controlPoints[selectedRow][selectedColumn].z);
	glEnd();
}

void drawControlMesh(Point** controlPoints, int rows, int cols) {
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material_black);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_black);
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

Point calculateNormal(Point** surfacePoints, int i, int j, int rows, int cols) {
	Point normal = { 0.0f, 0.0f, 0.0f };

	// Get neighboring points for tangent vectors
	Point right = (j < cols - 1) ?
		Point{
			surfacePoints[i][j + 1].x - surfacePoints[i][j].x,
			surfacePoints[i][j + 1].y - surfacePoints[i][j].y,
			surfacePoints[i][j + 1].z - surfacePoints[i][j].z
	} :
		Point{
			surfacePoints[i][j].x - surfacePoints[i][j - 1].x,
			surfacePoints[i][j].y - surfacePoints[i][j - 1].y,
			surfacePoints[i][j].z - surfacePoints[i][j - 1].z
	};

	Point down = (i < rows - 1) ?
		Point{
			surfacePoints[i + 1][j].x - surfacePoints[i][j].x,
			surfacePoints[i + 1][j].y - surfacePoints[i][j].y,
			surfacePoints[i + 1][j].z - surfacePoints[i][j].z
	} :
		Point{
			surfacePoints[i][j].x - surfacePoints[i - 1][j].x,
			surfacePoints[i][j].y - surfacePoints[i - 1][j].y,
			surfacePoints[i][j].z - surfacePoints[i - 1][j].z
	};

	// Calculate cross product for normal
	normal.x = right.y * down.z - right.z * down.y;
	normal.y = right.z * down.x - right.x * down.z;
	normal.z = right.x * down.y - right.y * down.x;

	// Normalize the vector
	float length = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
	if (length > 0) {
		normal.x /= length;
		normal.y /= length;
		normal.z /= length;
	}

	normal.x *= -1;
	normal.y *= -1;
	normal.z *= -1;
	return normal;
}

void drawSplineSurface(Point** surfacePoints, int rows, int cols, bool fill) {
	if (fill) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		GLfloat material_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
		GLfloat material_diffuse[] = { 0.7f, 0.7f, 0.7f, 1.0f };
		GLfloat material_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
		GLfloat material_shininess[] = { 50.0f };

		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, material_shininess);
	}
	else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material_black);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_black);
	}

	// Enable smooth shading
	glShadeModel(GL_SMOOTH);

	for (int i = 0; i < rows - 1; i++) {
		glBegin(GL_TRIANGLE_STRIP);
		for (int j = 0; j < cols; j++) {
			// Calculate and set normal for the current vertex
			Point normal1 = calculateNormal(surfacePoints, i, j, rows, cols);
			glNormal3f(normal1.x, normal1.y, normal1.z);
			glVertex3f(surfacePoints[i][j].x, surfacePoints[i][j].y, surfacePoints[i][j].z);

			// Calculate and set normal for the next vertex
			Point normal2 = calculateNormal(surfacePoints, i + 1, j, rows, cols);
			glNormal3f(normal2.x, normal2.y, normal2.z);
			glVertex3f(surfacePoints[i + 1][j].x, surfacePoints[i + 1][j].y, surfacePoints[i + 1][j].z);
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


void generateBezierSurface(Point** surfacePoints, int numSurfaceRows, int numSurfaceCols, Point** controlPoints, int crows, int ccols, float resolution) {
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
			surfacePoints[int(u * (numSurfaceRows - 1))][int(v * (numSurfaceCols - 1))] = Point{ nextPoint.x, nextPoint.y, nextPoint.z };
		}
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

void generateBSplineSurfacePoints(Point**& surfacePoints, int& numSurfaceRows, int& numSurfaceCols) {

	for (float u = 0; u < 1 + resolution; u += resolution) {
		for (float v = 0; v < 1 + resolution; v += resolution) {
			Point nextPoint = { 0.0f, 0.0f, 0.0f };
			computeBSplineSurface(u, v, rows, cols, nextPoint);
			surfacePoints[int(u * (numSurfaceRows - 1))][int(v * (numSurfaceCols - 1))] = nextPoint;
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


// -------------------------- NURBS ------------------------------------ //

int numControlPointsU = 0;
int numControlPointsV = 0;

std::vector<float> knotVectorU;
std::vector<float> knotVectorV;

const int degreeU = 2;
const int degreeV = 2;

float** weights = nullptr;

// Function to dynamically allocate a 2D array of weights
float** createWeights(int numControlPointsU, int numControlPointsV, float initialValue = 1.0f) {
	// Allocate an array of pointers for the rows
	float** weights = new float* [numControlPointsU];

	// Allocate each row and initialize values
	for (int i = 0; i < numControlPointsU; ++i) {
		weights[i] = new float[numControlPointsV];
		for (int j = 0; j < numControlPointsV; ++j) {
			weights[i][j] = initialValue; // Initialize with the given initial value
		}
	}

	return weights;
}

void deleteWeights() {
	if (weights != nullptr) {
		for (int i = 0; i < numControlPointsU; ++i) {
			delete[] weights[i];
		}
		delete[] weights;
		weights = nullptr;
		std::cout << "Weights array deallocated.\n";
	}
}

std::vector<float> generateUniformKnotVector(int numControlPoints, int degree) {
	int numKnots = numControlPoints + degree + 1;
	std::vector<float> knotVector(numKnots);

	for (int i = 0; i < numKnots; ++i) {
		if (i < degree) {
			knotVector[i] = 0.0f;
		}
		else if (i > numKnots - degree - 1) {
			knotVector[i] = 1.0f;
		}
		else {
			knotVector[i] = float(i - degree) / (numKnots - 2 * degree - 1);
		}
	}

	return knotVector;
}

// Helper function to calculate the basis function
float N(int i, int p, float u, const std::vector<float>& knotVector) {
	if (p == 0) {
		return (knotVector[i] <= u && u < knotVector[i + 1]) ? 1.0f : 0.0f;
	}
	float left = 0.0f;
	float right = 0.0f;

	if (knotVector[i + p] - knotVector[i] != 0) {
		left = ((u - knotVector[i]) / (knotVector[i + p] - knotVector[i])) * N(i, p - 1, u, knotVector);
	}
	if (knotVector[i + p + 1] - knotVector[i + 1] != 0) {
		right = ((knotVector[i + p + 1] - u) / (knotVector[i + p + 1] - knotVector[i + 1])) * N(i + 1, p - 1, u, knotVector);
	}
	return left + right;
}

// Function to compute a point on the NURBS surface
void computeSurfacePoint(float u, float v, float& x, float& y, float& z) {
	x = y = z = 0.0f;
	float wSum = 0.0f;

	for (int i = 0; i < numControlPointsU; i++) {
		for (int j = 0; j < numControlPointsV; j++) {
			float Nu = N(i, degreeU, u, knotVectorU);
			float Nv = N(j, degreeV, v, knotVectorV);
			float weight = weights[i][j] * Nu * Nv;

			x += weight * controlPoints[i][j].x;
			y += weight * controlPoints[i][j].y;
			z += weight * controlPoints[i][j].z;
			wSum += weight;
		}
	}

	x /= wSum;
	y /= wSum;
	z /= wSum;
}

void drawNurbsSurface() {
	glColor3f(0.0, 1.0, 0.0);
	glPointSize(3.0);

	float resolution = 0.0999f;
	float uMin = knotVectorU[degreeU];
	float uMax = knotVectorU[knotVectorU.size() - degreeU - 1];
	float vMin = knotVectorV[degreeV];
	float vMax = knotVectorV[knotVectorV.size() - degreeV - 1];

	// Draw u-direction lines
	for (float v = vMin; v <= vMax; v += resolution) {
		glBegin(GL_LINE_STRIP);
		for (float u = uMin; u <= uMax; u += resolution) {
			float x, y, z;
			computeSurfacePoint(u, v, x, y, z);
			glVertex3f(x, y, z);
		}
		glEnd();
	}

	// Draw v-direction lines
	for (float u = uMin; u <= uMax; u += resolution) {
		glBegin(GL_LINE_STRIP);
		for (float v = vMin; v <= vMax; v += resolution) {
			float x, y, z;
			computeSurfacePoint(u, v, x, y, z);
			glVertex3f(x, y, z);
		}
		glEnd();
	}
}

void generateNurbsSurfacePoints(Point**& surfacePoints, int& numSurfaceRows, int& numSurfaceCols) {
	// Get parameter ranges
	float uMin = knotVectorU[degreeU];
	float uMax = knotVectorU[knotVectorU.size() - degreeU - 1];
	float vMin = knotVectorV[degreeV];
	float vMax = knotVectorV[knotVectorV.size() - degreeV - 1];

	// Generate surface points
	for (float u = uMin; u <= uMax; u += resolution) {
		for (float v = vMin; v <= vMax; v += resolution) {
			float x, y, z;
			computeSurfacePoint(u, v, x, y, z);
			Point nextPoint = { x, y, z };
			surfacePoints[int((u - uMin) / (uMax - uMin) * (numSurfaceRows - 1))]
				[int((v - vMin) / (vMax - vMin) * (numSurfaceCols - 1))] = nextPoint;
		}
	}
}

// --------------------------- End of NURBS ---------------------------- //

void cleanupSurfacePoints() {
	if (surfacePoints != nullptr) {
		for (int i = 0; i < numSurfaceRows; ++i) {
			delete[] surfacePoints[i];
		}
		delete[] surfacePoints;
		surfacePoints = nullptr;
	}
}

void cleanupControlPoints() {
	if (controlPoints != nullptr) {
		for (int i = 0; i < rows; ++i) {
			delete[] controlPoints[i];
		}
		delete[] controlPoints;
		controlPoints = nullptr;
	}
}

// Register cleanup function to be called on program exit
void registerCleanup() {
	atexit(deleteWeights);
	atexit(cleanupControlPoints);
	atexit(cleanupSurfacePoints);
}

void drawCoordinateSystem() {
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material_red);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_red);
	glBegin(GL_LINES);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(axisSize, 0.0, 0.0);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material_blue);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_blue);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, axisSize, 0.0);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, material_green);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_green);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, axisSize);
	glEnd();
}


// Drawing routine.
void drawScene(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
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

	if (showControlHull) {
		drawControlMesh(controlPoints, rows, cols);
	}
	if (showBezierSurface) {
		generateBezierSurface(surfacePoints, numSurfaceRows, numSurfaceCols, controlPoints, rows, cols, resolution);
		drawSplineSurface(surfacePoints, numSurfaceRows, numSurfaceCols, fill);
	}
	if (showBSplineSurface)
	{
		generateBSplineSurfacePoints(surfacePoints, numSurfaceRows, numSurfaceCols);
		drawSplineSurface(surfacePoints, numSurfaceRows, numSurfaceCols, fill);
	}
	if (showNurbsSurface)
	{
		generateNurbsSurfacePoints(surfacePoints, numSurfaceRows, numSurfaceCols);
		drawSplineSurface(surfacePoints, numSurfaceRows, numSurfaceCols, fill);
	}
	if (showControlPoints) {
		drawControlPoints(controlPoints, rows, cols);
	}
	glPopMatrix();

	glFlush();
}

void setupLighting() {
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);

	GLfloat material_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	GLfloat material_diffuse[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	GLfloat material_specular[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat material_shininess[] = { 25.0f };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, material_shininess);


	GLfloat light_position[] = { 0.0f, 10.0f, 5.0f, 0.0f };
	GLfloat light_ambient[] = { 0.15f, 0.15f, 0.15f, 1.0f };
	GLfloat light_diffuse[] = { 0.6f, 0.6f, 0.6f, 1.0f };
	GLfloat light_specular[] = { 0.4f, 0.4f, 0.4f, 1.0f };

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
}

// Initialization routine.
void setup(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	setupLighting();
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
		showNurbsSurface = false;
		showBezierSurface = !showBezierSurface;
		glutPostRedisplay();
		break;
	case '4':
		showBezierSurface = false;
		showNurbsSurface = false;
		showBSplineSurface = !showBSplineSurface;
		glutPostRedisplay();
		break;
	case '5':
		showBSplineSurface = false;
		showBezierSurface = false;
		showNurbsSurface = !showNurbsSurface;
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
	case 'r':
		initialieControlPoints(controlPoints, rows, cols);
		deleteWeights();
		weights = createWeights(rows, cols);
		std::cout << "Control points and weights reseted." << std::endl;
		glutPostRedisplay();
		break;
	case 'i':
		weights[selectedRow][selectedColumn] += 0.01;
		glutPostRedisplay();
		break;
	case 'k':
		weights[selectedRow][selectedColumn] -= 0.01;
		glutPostRedisplay();
		break;
	case 9:
		if (selectedRow < rows-1) selectedRow++;
		else selectedRow = 0;
		glutPostRedisplay();
		break;
	case ' ':
		if (selectedColumn < cols-1) selectedColumn++;
		else selectedColumn = 0;
		glutPostRedisplay();
	default:
		break;
	}
}

void specialKeyInput(int key, int x, int y)
{
	if (key == GLUT_KEY_LEFT) controlPoints[selectedRow][selectedColumn].x -= 0.1;
	if (key == GLUT_KEY_RIGHT) controlPoints[selectedRow][selectedColumn].x += 0.1;
	if (key == GLUT_KEY_DOWN) controlPoints[selectedRow][selectedColumn].y -= 0.1;
	if (key == GLUT_KEY_UP) controlPoints[selectedRow][selectedColumn].y += 0.1;
	if (key == GLUT_KEY_PAGE_DOWN) controlPoints[selectedRow][selectedColumn].z -= 0.1;
	if (key == GLUT_KEY_PAGE_UP) controlPoints[selectedRow][selectedColumn].z += 0.1;
	glutPostRedisplay();
}

void printUserManual() {
	std::cout << "Press '0' to hide/show coordinate system." << std::endl;
	std::cout << "Press '1' to hide/show control points." << std::endl;
	std::cout << "Press '2' to hide/show control hull." << std::endl;
	std::cout << "Press '3' to hide/show Bezier surface." << std::endl;
	std::cout << "Press '4' to hide/show B-spline surface." << std::endl;
	std::cout << "Press '5' to hide/show NURBS surface." << std::endl;
	std::cout << "Press 'f' to switch between wired and filled surface." << std::endl;
	std::cout << "Press 'x'/'X' to rotate around 'X' axis." << std::endl;
	std::cout << "Press 'y'/'Y' to rotate around 'Y' axis." << std::endl;
	std::cout << "Press 'c'/'C' to rotate around 'Z' axis." << std::endl;
	std::cout << "Press 'r' to reset control points and weights." << std::endl;
	std::cout << "Press space and tab to select a control point." << std::endl;
	std::cout << "Press the right/left arrow keys to move the control point up/down the x-axis." << std::endl;
	std::cout << "Press the up/down arrow keys to move the control point up/down the y-axis." << std::endl;
	std::cout << "Press the page up/down keys to move the control point up/down the z-axis." << std::endl;
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
		return 1;
	}

	controlPoints = createPointMatrx(rows, cols);
	initialieControlPoints(controlPoints, rows, cols);
	// printMatrix();
	surfacePoints = createPointMatrx(numSurfaceRows, numSurfaceCols);
	
	uLength = rows + 2 + 1; // controlPoint + degree + 1
	vLength = cols + 2 + 1; // same

	knot_u = new float[uLength];
	knot_v = new float[vLength];

	initiateKnotVector(knot_u, uLength, rows);
	initiateKnotVector(knot_v, vLength, cols);

	numControlPointsU = rows;
	numControlPointsV = cols;

	knotVectorU = generateUniformKnotVector(numControlPointsU, degreeU);
	knotVectorV = generateUniformKnotVector(numControlPointsV, degreeV);

	
	registerCleanup();
	weights = createWeights(rows, cols);

	printUserManual();

	glutInit(&argc, argv);

	glutInitContextVersion(4, 3);
	glutInitContextProfile(GLUT_COMPATIBILITY_PROFILE);

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);

	glutInitWindowSize(500, 500);
	glutInitWindowPosition(700, 100);

	glutCreateWindow("spline-surface.cpp");

	glutDisplayFunc(drawScene);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keyInput);
	glutSpecialFunc(specialKeyInput);

	glewExperimental = GL_TRUE;
	atexit(cleanupSurfacePoints);
	glewInit();

	setup();

	glutMainLoop();
}
