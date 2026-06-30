#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>

// ================================================
// Struct definitions
// ================================================

typedef struct {
    double x,y,z;
} vec3;

typedef struct {
    vec3 position, velocity;
    double mass, charge;
    // Stores a unique identifier for each body, to avoid self comparison
    int id;
} body;

// ================================================
// Random generation functions
// ================================================

double* randomContinuousPositive(double maxMag, int N) {
    double* adress = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        double val = maxMag * rand()/RAND_MAX;
        *(adress + i) = val;
    }
    return adress;
}

double* randomContinuous(double maxMag, int N) {
    double *adress = calloc(N, sizeof(double));

    for (int i = 0; i < N; i++) {
        if ((double)rand()/RAND_MAX > 0.5) {
            double val = maxMag * rand()/RAND_MAX;
            *(adress + i) = val;
        }
        else {
            double val = - maxMag * rand()/RAND_MAX;
            *(adress + i) = val;
        }
    }
    return adress;
}

// ================================================
// Progress printing function
// ================================================

void printProgress(int stepVal, int tSteps, int width) {
    double perVal = (double) stepVal / (tSteps) * 100;

    char progBar[] = "========================================================================================================================";
    char emptyBar[] = "------------------------------------------------------------------------------------------------------------------------";
    char empty[] = "";

    // Width of the filled bar
    int filledWidth = (int) (perVal / 100 * width);
    int emptyPad = width - filledWidth;
    printf("%.s", emptyBar);
    // Bunch of mumbo jumbo I took ages to figure out
    // First bit (%5.2f%) just means 'print this number to 2 decimal places, and make it take up at least 5 places'
    // %.*s Basically just means 'Input how many characters to print of a string, input the string to be truncated'
    printf("\r %5.2f%% [%.*s%.*s]", perVal, filledWidth, progBar, emptyPad, emptyBar);
    fflush(stdout);
}

// Rendering helper function
void writeFrame(FILE *dataFile, body *bodies, unsigned int N, float *frameBuffer) {
    for (int i = 0; i < N; i++) {
        frameBuffer[3*i] = (float)bodies[i].position.x;
        frameBuffer[(3*i)+1] = (float)bodies[i].position.y;
        frameBuffer[(3*i)+2] = (float)bodies[i].position.z;
    }
    fwrite(frameBuffer, sizeof(float), 3*N, dataFile);
}

// ================================================
// SPH helper functions
// ================================================


// Kernel function
// This is the cubic spline kernel (mon92)
// if the distance r is greater than or equal to 2, the kernel should have a zero value


double kernel(double r, double normConst, double h) {
    double W = 0.0;
    double q = r/h;
    if (q >= 0 && q <= 2) {
        if (q < 1) {
            W = 2.0/3.0 - (q*q) + 1.0/2.0 * (q*q*q);
        }
        else {
            W = 1.0/6.0 * (2-q) * (2-q) * (2-q);
        }
    }
    else {
        return 0.0;
    }
    return W * normConst;
}

// Distance helper function

vec3 getSeparation(vec3 pos1, vec3 pos2) {
    double dx = pos1.x - pos2.x;
    double dy = pos1.y - pos2.y;
    double dz = pos1.z - pos2.z;
    return (vec3) {dx,dy,dz};
}

double getDistMag(vec3 pos1, vec3 pos2) {
    double dx = pos1.x - pos2.x;
    double dy = pos1.y - pos2.y;
    double dz = pos1.z - pos2.z;  
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    return r;
}

// Calculates the gradient of the smoothing kernel 

double gradKernelScalar(vec3 pos1, vec3 pos2, double h, double normConst) {
    double gradW = 0.0;
    vec3 dis = getSeparation(pos1, pos2);
    double r2 = dis.x * dis.x + dis.y * dis.y + dis.z * dis.z;
    if (r2/(h*h) >= 0.0 && r2/(h*h) <= 4.0) {
        double q = sqrt(r2)/h;
        if (q < 1) {
            gradW = -2*q + 3.0/2.0 * (q*q);
        }
        else {
            gradW = -1.0/2.0 * (2-q)*(2-q);
        }
    }

    return gradW * normConst;
}

double calculateDensity(vec3 pos, body* bodies, int N, double h, double normConst) {
    double density = 0.0;
    // Only particles inside the smoothing radius should be visited, so I will have to optimise this later
    for (int i = 0; i < N; i++) {
        double dist = getDistMag(pos, bodies[i].position);
        double weighting = kernel(dist,normConst,h);
        density += bodies[i].mass * weighting;
    }
    return density;
}


// Calculate prssure
// if n = 1, the power just simplifies to 2 using the polytropic equation of state
inline double getPressure(double density, double kConst, double nConst, body* bodies) {
    //double Pressure = k * pow(density, 1.0 + 1.0/n);
    double Pressure = kConst * density * density;
    return Pressure;
}

void SPHAccelCalc(body* bodies, vec3* accels, int N, double h, double kConst, double nConst, double normConst) {

    // initialise accels array to zero
    for (int i = 0; i < N; i++) {
        accels[i].x = accels[i].y = accels[i].z = 0.0;
    }

    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j ) {
                continue;
            }
            vec3 sep = getSeparation(bodies[i].position,bodies[j].position);
            double r2 = ((sep.x*sep.x) + (sep.y*sep.y) + (sep.z*sep.z));
            if (r2 > 5.0) {
                continue;
            }

            double densityI = calculateDensity(bodies[i].position,bodies,N,h,normConst);
            double densityJ = calculateDensity(bodies[i].position,bodies,N,h,normConst);
            double pressureI = getPressure(densityI, kConst, nConst, bodies);
            double pressureJ = getPressure(densityJ, kConst, nConst, bodies);

            // Gradient values
            double WScalar = gradKernelScalar(bodies[i].position, bodies[j].position, h, normConst);
            
            double constTerm = bodies[i].mass * (pressureI/(densityI*densityI) + pressureJ/(densityJ*densityJ));
            //printf("WScalar term: %f", WScalar);
            accels[i].x -= constTerm * (WScalar * sep.x);
            accels[i].y -= constTerm * (WScalar * sep.y);
            accels[i].z -= constTerm * (WScalar * sep.z);
        }
        //printf("Acceleration of body 0: (%f, %f, %f)", accels[0].x, accels[0].y, accels[0].z);
    }
}




void timeStep(body* bodies, vec3* accels, double dt, int N, double g, double k, double e, double xmax, double ymax, double zmax, double n, double normConst, double h, double kConst) {
    // Half kick and drift
    for (int i = 0; i < N; i++) {
        // The half kick
        bodies[i].velocity.x += accels[i].x * dt / 2;
        bodies[i].velocity.y += accels[i].y * dt / 2;
        bodies[i].velocity.z += accels[i].z * dt / 2;

        // The full drift
        bodies[i].position.x += bodies[i].velocity.x * dt;
        bodies[i].position.y += bodies[i].velocity.y * dt;
        bodies[i].position.z += bodies[i].velocity.z * dt;
    }

    //printf("Body number 0 has position: (%f, %f, %f) \n", bodies[0].position.x,bodies[0].position.y,bodies[0].position.z);
    //printf("Body number 0 has velocity: (%f, %f, %f) \n", bodies[0].velocity.x,bodies[0].velocity.y,bodies[0].velocity.z);
    //printf("Body number 0 has acceleration: (%f, %f, %f) \n", accels[0].x,accels[0].y,accels[0].z);
    // Calculate new accelerations
    SPHAccelCalc(bodies, accels, N, h, kConst, n, normConst);

    // Last half kick

    for (int i = 0; i < N; i++) {
        // The half kick
        bodies[i].velocity.x += accels[i].x * dt / 2;
        bodies[i].velocity.y += accels[i].y * dt / 2;
        bodies[i].velocity.z += accels[i].z * dt / 2;
    }

    for (int i = 0; i < N; i++) {
        // The effects of gravity
        bodies[i].velocity.z -= g * dt;

        // Velocity damping
        bodies[i].velocity.x -= bodies[i].velocity.x * k;
        bodies[i].velocity.y -= bodies[i].velocity.y * k;
        bodies[i].velocity.z -= bodies[i].velocity.z * k;
    }

    // Collisions with the wall
    #pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < N; i++) {
        if (bodies[i].position.x > xmax) {
            //printf("Position before contact: %f, velocity before contact: %f \n", bodies[i].position.x, bodies[i].velocity.x);
            double xDiff = bodies[i].position.x - xmax;
            bodies[i].position.x = xmax - xDiff;
            bodies[i].velocity.x = - e * bodies[i].velocity.x;
            //printf("Position after contact: %f, velocity after contact: %f \n", bodies[i].position.x, bodies[i].velocity.x);
        }
        if (bodies[i].position.x < -xmax) {
            double xDiff = - bodies[i].position.x - xmax;
            bodies[i].position.x = -xmax - xDiff;
            bodies[i].velocity.x = - e * bodies[i].velocity.x;
        }
        if (bodies[i].position.y > ymax) {
            double yDiff = bodies[i].position.y - ymax;
            bodies[i].position.y = ymax - yDiff;
            bodies[i].velocity.y = - e * bodies[i].velocity.y;
        }
        if (bodies[i].position.y < -ymax) {
            double yDiff = - bodies[i].position.y - ymax;
            bodies[i].position.y = -ymax - yDiff;
            bodies[i].velocity.y = - e * bodies[i].velocity.y;
        }
        if (bodies[i].position.z > zmax) {
            double zDiff = bodies[i].position.z - zmax;
            bodies[i].position.z = zmax-zDiff;
            bodies[i].velocity.z = - e * bodies[i].velocity.z;
        }
        if (bodies[i].position.z < -zmax) {
            double zDiff = - bodies[i].position.z - zmax;
            bodies[i].position.z = -zmax-zDiff;
            bodies[i].velocity.z = - e * bodies[i].velocity.z;
        }
    }

}

int main() {

    double xMax = 20;
    double yMax = 20;
    double zMax = 40;

    double dt = 0.005;

    // For 2d the kernel normalisation constant is:
    double h = 1.0; // smoothing radius
    double normConst = 15.0 / (7 * 3.14159 * h * h);

    double n = 1.0; // polytropic index
    double kConst = 1.0; // Equation of state constant
    // Choosing sensible parameters
    /*
    The average atomic spacing of iron atoms is ~2.86 x 10^-10 m. 
    Since r_avg = 2^(1/6) * \sigma, -> /sigma = r_avg / (2^(1/6))
    Iron also has a binding energy of 706.7 eV, so we set epsiolon to 706.7 x 1.61 x 10^(-19)
    */
    //double r_avg = 2.86e-10;
    //double sigma = r_avg/(pow(2,1.0/6.0));
    //printf("Sigma %e \n", sigma);
    //double epsilon = 706.7*1.61e-19; 
    //double m = 1e-26;

    double sigma = 1.0;
    double epsilon = 1.0;
    double m = 100.0;

    int tSteps = 4000;

    double softener = 0.1;

    // force constants
    double g = 9.81; // gravitational acceleration
    double k = 0.001; // velocity damping
    double e = 0.9; // coefficient of restitution

    // Initialise bodies. 
    // Bodies will be spawned along a grid
    int gridX = 20; int gridY = 20; int gridZ = 10;
    double sep = 1.0;
    int N = gridX*gridY*gridZ;
    body* bodies = calloc(N, sizeof(body));
    vec3* accels = calloc(N, sizeof(vec3));
    int id = 0;
    double xOff = gridX * 1.0/2.0 * sep;
    double yOff = gridY * 1.0/2.0 * sep;
    double zOff = gridZ * 1.0/2.0 * sep;
    printf("Initialising bodies \n");
    for (int x = 0; x < gridX; x++) {
        for (int y = 0; y < gridY; y++) {
            for (int z = 0; z < gridZ; z++) {
                bodies[id].mass = m;
                bodies[id].position.x = x * sep - xOff;
                bodies[id].position.y = y * sep - yOff;
                bodies[id].position.z = z * sep - zOff;
                bodies[id].id = id;
                id++;
            }
        }
    }
    printf("Initialised bodies \n");
    // Initialise acceleration to zero
    for (int i = 0; i < N; i++) {
        accels[i].x = accels[i].y = accels[i].z = 0.0;
    }

    // File handling
    FILE* dataFile;
    // Writing parameter data to file
    // Number of bodies, number of time steps, size of the box in (X,Y,Z)
    int params[] = {N,tSteps,h,xMax,yMax,zMax};
    int numParam = sizeof(params) / sizeof(params[0]);
    dataFile = fopen("LiquidBoxData.bin", "wb");
    fwrite(params, sizeof(int), numParam, dataFile);
    // Create the frame buffer
    float * frameBuffer = calloc(3*N, sizeof(float));

    for (int i = 0; i < tSteps; i++) {
        timeStep(bodies, accels,dt,N,g,k,e,xMax,yMax,zMax,n,normConst,h, kConst);
        printProgress(i, tSteps, 50);
        writeFrame(dataFile,bodies,N,frameBuffer);
    }
    printProgress(tSteps, tSteps, 50);

    free(bodies);
    free(accels);
    free(frameBuffer);


    return 0;
}