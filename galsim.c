
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graphics/graphics.h"
#include <omp.h>

// Struct desrcibing a particle
typedef struct part {
    double x_force;
    double y_force;
    double x;
    double y;
    double mass;
    double x_vel;
    double y_vel;
    double brightness;
} part;

typedef struct QT {
        struct QT *child1;
        struct QT *child2;
        struct QT *child3;
        struct QT *child4;
        double mass;
        double C_O_x;
        double C_O_y;
        double posx;
        double posy;
        double side;
        int nop;
        struct part ** particles;
} QT;






void vel_update(part** particles, int N, double t) {

    #pragma omp parallel
            {
        #pragma omp for 
    for(int i = 0; i < N; i++) {
        particles[i]->x_vel = particles[i]->x_vel + ((particles[i]->x_force)/(particles[i]->mass)) * t;
        particles[i]->y_vel = particles[i]->y_vel + ((particles[i]->y_force)/(particles[i]->mass)) * t;    
    }
    }
}

void pos_update(part** particles, int N, double t) {

        #pragma omp parallel
            {
            #pragma omp for 
        for(int i = 0; i < N; i++) {
            particles[i]->x = particles[i]->x + (particles[i]->x_vel) * t;
            particles[i]->y = particles[i]->y + (particles[i]->y_vel) * t;
            particles[i]->y_force = 0;
            particles[i]->x_force = 0;

    }
            }

}
void centerOfMass(QT * node){
    double mass = 0;
    double x = 0;
    double y = 0;
    for(int i=0; i<node->nop; i++) {
        x += node->particles[i]->mass*node->particles[i]->x;
        y += node->particles[i]->mass*node->particles[i]->y;
        mass += node->particles[i]->mass;
        }
    node->C_O_x = x/(mass);
    node->C_O_y = y/(mass);
    node->mass = mass;
    } 

void create_tree(QT * node) {

        QT * child1 = (QT *)malloc(sizeof(QT));
        QT * child2 = (QT *)malloc(sizeof(QT));
        QT * child3 = (QT *)malloc(sizeof(QT));
        QT * child4 = (QT *)malloc(sizeof(QT));

        double new_side = (node->side)/2;
        double posx = (node->posx);
        double posy = (node->posy);

        node->child1 = child1;
        node->child2 = child2;
        node->child3 = child3;
        node->child4 = child4;

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
                
        
        part **part1=(part**)malloc((node->nop)*sizeof(part*));
        part **part2=(part**)malloc((node->nop)*sizeof(part*));
        part **part3=(part**)malloc((node->nop)*sizeof(part*));
        part **part4=(part**)malloc((node->nop)*sizeof(part*));

        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        int count4 = 0;

        for(int i=0; i<node->nop; i++) {
            if(node->particles[i]->x < posx) {
                if(node->particles[i]->y < posy) {
                    part1[count1] = node->particles[i];
                    count1++;
                } else {
                    part3[count3] = node->particles[i];
                    count3++;
                }
            } else {
                if(node->particles[i]->y < posy) {
                    part2[count2] = node->particles[i]; 
                    count2++;   
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
        
        if(count1 > 1) {
            centerOfMass(node->child1);
            create_tree(node->child1);
        }
        if(count2 > 1) {
            centerOfMass(node->child2);
            create_tree(node->child2);
        }
        if(count3 > 1) {
            centerOfMass(node->child3);
            create_tree(node->child3);
        }
        if(count4 > 1) {
            centerOfMass(node->child4);
            create_tree(node->child4);
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
        return;
    }
}

void force(part * p, QT * restrict n, double theta_max, double G){
    double x = 0;
    double y = 0;
    const double epsilon = 0.001;
    if(n->nop == 0) {
        return;
    }    
    if(n->nop>1){
		double dnodecenter=sqrt((((p->x)-(n->posx))*((p->x)-(n->posx)))+(((p->y)-(n->posy))*((p->y)-(n->posy))));
		double theta=(n->side)/dnodecenter;


		if(theta<=theta_max){
			double dmasscenter= sqrt((((p->x)-(n->C_O_x))*((p->x)-(n->C_O_x)))+(((p->y)-(n->C_O_y))*((p->y)-(n->C_O_y))));
			x= (-G*(p->mass) * (n->mass) * ((p->x)-(n->C_O_x)))/((dmasscenter+epsilon)*(dmasscenter+epsilon)*(dmasscenter+epsilon));
			y= (-G*(p->mass) * (n->mass) * ((p->y)-(n->C_O_y)))/((dmasscenter+epsilon)*(dmasscenter+epsilon)*(dmasscenter+epsilon));
			p->x_force=p->x_force+x;
			p->y_force=p->y_force+y;
			return;
		}else{
			force(p,n->child1,G,theta_max);
			force(p,n->child2,G,theta_max);
			force(p,n->child3,G,theta_max);
			force(p,n->child4,G,theta_max);
		}
	}
	else if (n->nop==1 && p->x !=n->C_O_x){
			double dmasscenter= sqrt(((p->x)-(n->particles[0]->x))*((p->x)-(n->particles[0]->x))+((p->y)-(n->particles[0]->y))*((p->y)-(n->particles[0]->y)));
			x= -G*(p->mass) * (n->particles[0]->mass) * ((p->x)-(n->particles[0]->x))/((dmasscenter+epsilon)*(dmasscenter+epsilon)*(dmasscenter+epsilon));
			y= -G*(p->mass) * (n->particles[0]->mass) * ((p->y)-(n->particles[0]->y))/((dmasscenter+epsilon)*(dmasscenter+epsilon)*(dmasscenter+epsilon));
			p->x_force +=x;
			p->y_force += y;
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
    omp_set_num_threads(atoi(argv[7]));



    // Creates array to store particles
    part **array = (part**) malloc(N*sizeof(part*));
    for(int i = 0; i < N; i++)
    {
        array[i] = (part*) malloc(sizeof(part));
    }



    // Read file
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
        centerOfMass(root);
        while(t < nsteps){
            ClearScreen();
    
            create_tree(root);
           
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
        QT * root = (QT *)malloc(sizeof(QT));
        root->posx = 0.5;
        root->posy = 0.5;
        root->side = 1;
        root->particles = array;
        root->nop = N;
        centerOfMass(root);

        while(t < nsteps){
            
            create_tree(root);
            
            //#pragma omp parallel
            //{
            //    #pragma omp for 
            for(int i = 0; i < N; i++)  
            {    
                force(array[i], root, theta_max, G);
            }
            //partnew
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