
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graphics/graphics.h"
#include <pthread.h>
#include <omp.h>

typedef struct QuadTree {
        struct QuadTree *child1;
        struct QuadTree *child2;
        struct QuadTree *child3;
        struct QuadTree *child4;
        double mass;
        double C_O_x;
        double C_O_y;
        double posx;
        double posy;
        double side;
        int nop;
        struct particle ** particles;
} QT;

// Global variables
const double epsilon = 0.001;


// Struct desrcibing a particle
typedef struct particle {
    double x_force;
    double y_force;
    double x;
    double y;
    double mass;
    double x_vel;
    double y_vel;
    double brightness;
} part;


void vel_update(part** particles, int N, double t) {

    for(int i = 0; i < N; i++) {
        particles[i]->x_vel = particles[i]->x_vel + ((particles[i]->x_force)/(particles[i]->mass)) * t;
        particles[i]->y_vel = particles[i]->y_vel + ((particles[i]->y_force)/(particles[i]->mass)) * t;    
    }
}

void pos_update(part** particles, int N, double t) {

        for(int i = 0; i < N; i++) {
            particles[i]->x = particles[i]->x + (particles[i]->x_vel) * t;
            particles[i]->y = particles[i]->y + (particles[i]->y_vel) * t;
            particles[i]->y_force = 0;
            particles[i]->x_force = 0;

    }

}
void centerOfMass(QT * node, int N){
    double mass = 0;
    double x = 0;
    double y = 0;
    if(N == 0) {
        return;
    }
    else if(N==1) {
        node->C_O_x = node->particles[0]->x;
        node->C_O_y = node->particles[0]->y;
        node->mass = node->particles[0]->mass;                          
    }
    else {
        for(int i=0; i<N; i++) {
            mass += node->particles[i]->mass;
            x += node->particles[i]->mass*node->particles[i]->x;
            y += node->particles[i]->mass*node->particles[i]->y;
        }
        node->C_O_x = x/(mass);
        node->C_O_y = y/(mass);
        node->mass = mass;
    } 
}

void create_tree(QT * node, int N) {
        QT * child1 = (QT *)malloc(sizeof(QT));
        QT * child2 = (QT *)malloc(sizeof(QT));
        QT * child3 = (QT *)malloc(sizeof(QT));
        QT * child4 = (QT *)malloc(sizeof(QT)); 
        node->child1 = child1;
        node->child2 = child2;
        node->child3 = child3;
        node->child4 = child4;
        child1->mass = 0;
        child2->mass = 0;
        child3->mass = 0;
        child4->mass = 0;
        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        int count4 = 0;
        double new_side = node->side/2;
        double posx = node->posx;
        double posy = node->posy;
        node->child1->side = new_side;
        node->child1->posx = posx-new_side/2;
        node->child1->posy = posy-new_side/2;
        node->child2->side = new_side;
        node->child2->posx = posx + new_side/2;
        node->child2->posy = posy - new_side/2;
        node->child3->side = new_side;
        node->child3->posx = posx - new_side/2;
        node->child3->posy = posy + new_side/2;
        node->child4->side = new_side;
        node->child4->posx = posx + new_side/2;
        node->child4->posy = posy + new_side/2;
                
        part *part1[N];
        part *part2[N];
        part *part3[N];
        part *part4[N];


        for(int i=0; i<N; i++) {
            if(node->particles[i]->y <= posy) {
                if(node->particles[i]->x <= posx) {
                    part1[count1] = node->particles[i];
                    count1++;
                } else {
                    part2[count2] = node->particles[i];
                    count2++;
                }
            } else {
                if(node->particles[i]->x <= posx) {
                    part3[count3] = node->particles[i];
                    count3++;    
                } else {
                    part4[count4] = node->particles[i];
                    count4++;
                }
            }
        }
        node->child1->particles = part1;
        node->child2->particles = part2;
        node->child3->particles = part3;
        node->child4->particles = part4;
        node->child1->nop = count1;
        node->child2->nop = count2;
        node->child3->nop = count3;
        node->child4->nop = count4;
        
        centerOfMass(node->child1, count1);
        centerOfMass(node->child2, count2);
        centerOfMass(node->child3, count3);
        centerOfMass(node->child4, count4);
        
        if(count1 > 1) {
            create_tree(node->child1, count1);
        }
        if(count2 > 1) {
            create_tree(node->child2, count2);
        }
        if(count3 > 1) {
            create_tree(node->child3, count3);
        }
        if(count4 > 1) {
            create_tree(node->child4, count4);
        }
    
}

void free_tree(QT * node){
    if(node->nop <=1) {
        return;
    }
    else {

        free_tree(node->child1);
        free_tree(node->child2);
        free_tree(node->child3);
        free_tree(node->child4);
        free(node->child1);
        free(node->child2);
        free(node->child3);
        free(node->child4);
        node->child1 = NULL;
        node->child2 = NULL;
        node->child3 = NULL;
        node->child4 = NULL;
        return;
    }
}

void force(part * particle, QT * restrict node, double theta, double G){
    if(node->nop == 0) {
        return;
    }
    else if(node->nop == 1) {

        double denomerator = (sqrt((particle->x-node->C_O_x)*(particle->x-node->C_O_x) + (particle->y-node->C_O_y)*(particle->y-node->C_O_y)) + epsilon)*
        (sqrt((particle->x-node->C_O_x)*(particle->x-node->C_O_x) + (particle->y-node->C_O_y)*(particle->y-node->C_O_y)) + epsilon)*
        (sqrt((particle->x-node->C_O_x)*(particle->x-node->C_O_x) + (particle->y-node->C_O_y)*(particle->y-node->C_O_y)) + epsilon);
        double constant = - G * (particle->mass) * (node->mass)  /denomerator;
        
        double forcex = (particle->x - node->C_O_x) * constant;
        double forcey = (particle->y - node->C_O_y)* constant;

        particle->x_force += forcex;
        particle->y_force += forcey;
        return;
    }
    double theta_comp = (node->side)/sqrt((node->posx -particle->x)*(node->posx -particle->x) + (node->posy -particle->y)*(node->posy -particle->y));
    if(theta_comp <= theta) {
        
        double denomerator = (sqrt((particle->x-node->C_O_x)*(particle->x-node->C_O_x) + (particle->y-node->C_O_y)*(particle->y-node->C_O_y)) + epsilon)*
        (sqrt((particle->x-node->C_O_x)*(particle->x-node->C_O_x) + (particle->y-node->C_O_y)*(particle->y-node->C_O_y)) + epsilon)*
        (sqrt((particle->x-node->C_O_x)*(particle->x-node->C_O_x) + (particle->y-node->C_O_y)*(particle->y-node->C_O_y)) + epsilon);
        
        double forcex = - G * (particle->mass) * (node->mass) * ((particle->x) - (node->C_O_x)) / denomerator;
        double forcey = - G * (particle->mass) * (node->mass) * ((particle->y) - (node->C_O_y)) / denomerator;

        particle->x_force += forcex;
        particle->y_force += forcey;
        return;
    }
    else {
        force(particle, node->child1, theta, G);
        force(particle, node->child2, theta, G);
        force(particle, node->child3, theta, G);
        force(particle, node->child4, theta, G);
        return;
    }
}


int main(int argc, char* argv[]) {

    if (argc != 8) {
        printf("Wrong number of inputs");
        return -1;
    }
    // Read user input
    const int N = atoi(argv[1]);
    const char* filename = argv[2];
    const int nsteps = atoi(argv[3]);
    const double delta_t = atof(argv[4]);
    const double theta_max = atof(argv[5]);
    const int graphics = atoi(argv[6]);
    const int not = atoi(argv[7]);
    omp_set_num_threads(not);


     // Creates array to store particles
    struct particle **array = (part**) malloc(N*sizeof(part*));

    FILE *input = fopen(filename, "r");
    
	for (int i =0; i<N; i++){
		double attr[6] = {0}; //compiler will set all values to zero
		fread(attr,sizeof(double),6,input);
		array[i] = (part*)malloc(sizeof(part));
		array[i]->x = attr[0];
		array[i]->y = attr[1];
		array[i]->mass = attr[2];
		array[i]->x_vel = attr[3];
		array[i]->y_vel = attr[4];
		array[i]->brightness = attr[5];
		array[i]->x_force = 0;
		array[i]->y_force = 0;
	}
		fclose(input);



   


    // Read file



    // Prepare for the loop
    int t = 0;
    const double G = (double) 100/N;
    const float W = 1; 
    const float L = 1;

    // Loop with graphics
    
    if(graphics==1) {
        InitializeGraphics(argv[0], 800, 800);
        SetCAxes(0,1);

        QT * root = (QT*)malloc(sizeof(QT));
        root->posx = 0.5;
        root->posy = 0.5;
        root->side = 1;
        root->particles = array;
        root->nop = N;
        while(t < nsteps){
            ClearScreen();
    
            create_tree(root, N);
           
            for(int i = 0; i < N; i++)  
            {
                DrawCircle(array[i]->x, array[i]->y, W, L, 0.01, 0);
                force(array[i], root, theta_max, G);
            }
            vel_update(array, N, delta_t);
            pos_update(array, N, delta_t); 
            
            Refresh();
            usleep(3000);

            t+=1;
            free_tree(root);
        }
        free(root);
        FlushDisplay();
        CloseDisplay();
    }

    // Loop without graphics
    else {

            QT * root = (QT*)malloc(sizeof(QT));
            root->posx = 0.5;
            root->posy = 0.5;
            root->side = 1;
            root->particles = array;
            root->nop = N;
        while(t < nsteps){
            
            create_tree(root, N);
            #pragma omp parallel
            {
                #pragma omp for
            for(int i = 0; i < N; i++)  
            {    
                force(array[i], root, theta_max, G);
            }
            }
            //printf("First paricle force %lf \n", array[0]->x_force);
            vel_update(array, N, delta_t);
            pos_update(array, N, delta_t);
            
            t +=1;
            free_tree(root);
            
            
        } 
        free(root);                               
    }

    // Writing to new file
    FILE * pw = fopen("result.gal", "w");

    for(int i = 0; i < N; i++)
    {
        fwrite(&(array[i]->x), sizeof(double), 1, pw);
        fwrite(&(array[i]->y), sizeof(double), 1, pw);
        fwrite(&(array[i]->mass), sizeof(double), 1, pw);
        fwrite(&(array[i]->x_vel), sizeof(double), 1, pw);
        fwrite(&(array[i]->y_vel), sizeof(double), 1, pw);
        fwrite(&(array[i]->brightness), sizeof(double), 1, pw);

    }
    
    fclose(pw);

    for(int i = 0; i < N; i++)
    {
        free(array[i]);
    }
    free(array);
    
return 0;
}