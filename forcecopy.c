#include <stdio.h>
#include"mpi.h"
#include<stdlib.h>
#include<math.h>

typedef struct Particle {
   int force;
   float x, y, z;
} Particle;

typedef struct Force {
   float x_vec;
   float y_vec;
   float z_vec;
} Force;

void calculateSumForce(Particle *p, Particle *p_set2, int particlesPerCore, Force* forces);
void initializeParticles(Particle *p, int numParticles);
int randomInt(int a, int b);
float randomFloat(float a, float b);
Particle* alloc_array_particle(int number);
Force* alloc_array_force(int number);
void printParticles(Particle *p, int numParticles);
Force addForce(Force f1, Force f2);

int main( int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size < 2) {
        fprintf(stderr, "World size must be greater than 1");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    
    int numParticles = 1000;
    if(argc > 1) {
        numParticles = atoi(argv[1]);
    }
    int modulo = numParticles % world_size;
    if(modulo != 0) {
        numParticles = numParticles + (numParticles - ((numParticles / world_size) * world_size));
    }

    int particlesPerCore = numParticles / world_size; 

    MPI_Request request=MPI_REQUEST_NULL;
    Particle* p;
    Force* f;
    f = alloc_array_force(particlesPerCore);
    Force* total_f;
    f = alloc_array_force(particlesPerCore);
    int i;
    if (world_rank == 0) {
    	p  = alloc_array_particle(numParticles);
        initializeParticles(p, numParticles);	
        //printParticles(p, numParticles);
        for(i = 0; i < world_size; i ++) {
            MPI_Isend(&p[particlesPerCore * i], particlesPerCore * sizeof(Particle), MPI_BYTE, i, 0, MPI_COMM_WORLD, &request);
        }
    }

    double t1, t2; 

    t1 = MPI_Wtime(); 
    
    p  = alloc_array_particle(particlesPerCore);
    
    MPI_Recv(&p[0], numParticles * sizeof(Particle), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    calculateSumForce(p, p, particlesPerCore, f);
    
    Particle* q;	
    q  = alloc_array_particle(particlesPerCore);
    for (i = 1; i < world_size; i++) {
        MPI_Isend(&p[0], particlesPerCore * sizeof(Particle), MPI_BYTE, getMod((world_rank + i), world_size) , 0, MPI_COMM_WORLD, &request);

        MPI_Recv(&q[0], particlesPerCore * sizeof(Particle), MPI_BYTE, getMod((world_rank - i), world_size), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        calculateSumForce(p, q, particlesPerCore, f);
    }

    MPI_Isend(&f[0], particlesPerCore * sizeof(Force), MPI_BYTE, 0, world_rank + 20, MPI_COMM_WORLD, &request);

    
    t2 = MPI_Wtime(); 
    printf( "Elapsed time is %f\n", t2 - t1 ); 
    
    if(world_rank == 0) {
    	Force* total_f;
    	total_f = alloc_array_force(numParticles);
        for(i = 0; i < world_size; i++) {
           MPI_Recv(&total_f[particlesPerCore * i], particlesPerCore * sizeof(Force), MPI_BYTE, i, i + 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
   
        //printf("\nForce per particle: \n");
        //for(i = 0; i< numParticles; i++) {
        //    printf("%g, %g, %g\n", total_f[i].x_vec, total_f[i].y_vec, total_f[i].z_vec);
	//}
        
    }
    MPI_Finalize();
    return 0;
}

int getMod(int a, int b){ return (a%b+b)%b; }

void calculateSumForce(Particle *p, Particle *q, int particlesPerCore, Force* forces) {
   int i, j;
   for(i = 0;i < particlesPerCore; i++) {
       for(j = 0; j < particlesPerCore; j++) {
           float xval = p[i].x - q[j].x;
           float yval = p[i].y - q[j].y;
           float zval = p[i].z - q[j].z;
           float co_ord_mag = ((float)sqrt(pow(xval, 2) + pow(yval, 2) + pow(zval, 2)));
           if(co_ord_mag == 0) continue;
           float mag = (p[i].force * q[j].force)/pow(co_ord_mag, 3);
           forces[i].x_vec += xval * mag;
           forces[i].y_vec += yval * mag;
           forces[i].z_vec += zval * mag;
       } 
   }
}

void printParticles(Particle *p, int numParticles) {
   int i;
   printf("Particles: \n");
   for(i = 0;i < numParticles; i++) {
       printf("%d, %g, %g, %g\n", p[i].force, p[i].x, p[i].y, p[i].z);
       fflush(stdout);
   }
};

Particle* alloc_array_particle(int number) {
    Particle *data = (Particle *)malloc(number*sizeof(Particle));
    return data;
}


Force* alloc_array_force(int number) {
    Force *data = (Force *)calloc(number, sizeof(Force));
    return data;
}

float randomFloat(float a, float b) {
    float random = ((float) rand() + 1) / ((float) RAND_MAX + 1);
    float diff = b - a;
    float r = random * diff;
    return a + r;
}

Force addForce(Force f1, Force f2) {
    f1.x_vec += f2.x_vec;
    f1.y_vec += f2.y_vec;
    f1.z_vec += f2.z_vec;
    return f1;
}

int randomInt(int a, int b) {
    int random = rand() % (b-a);
    return random + a;
}

void initializeParticles(Particle *p, int numParticles) {
   int i;
   for(i = 0;i < numParticles; i++) {
       p[i].force = randomInt(-3, 3);
       p[i].x = randomFloat(0, 10);
       p[i].y = randomFloat(0, 10);
       p[i].z = randomFloat(0, 10);
   }
}
