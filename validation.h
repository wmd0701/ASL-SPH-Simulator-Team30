#include <stdlib.h>
#include <stdio.h>
#include "data_set.h"


const int N = 60;
const double dx = 0.02916;
const double dy = dx;
const int M = 38;

vector index_from_coord(vector particle_pos){
    double x = particle_pos.first;
    double y = particle_pos.second;
    if(!(x < 1.73 - 0.5*dx && y < 1.1 - 0.5*dy && x > -0.5*dx && y > -0.5*dy)){
        printf("Something wrong with particle position! \n");
    }
    
    vector result;
    result.first = x/dx + 0.5;
    result.second = y/dy + 0.5;
	return result;
}

double displace_to_origin(Particle* all_particle){
    double x_min = all_particle[0].position.first;
    for (int i = 1; i < NUMBER_OF_PARTICLE; ++i){
        if(all_particle[i].tag == repulsive){
            if(all_particle[i].position.first < x_min){
                x_min = all_particle[i].position.first;
            }
        }
    }
    
    for (int i = 1; i < NUMBER_OF_PARTICLE; ++i){
        all_particle[i].position.first -= x_min;
    }
    
    return x_min;
}

void validation(Particle* all_particle){
    //Displace all the particles, such that the left most repulsive particle has position.first = 0
    double x_min = displace_to_origin(all_particle);
    
    //Initialize the array which contains the coordinates in which the level set function is computed
    vector* points = (vector*)malloc(sizeof(vector) * (N+2) * (M+2));
    vector temp;
    for(int j = -1; j < M+1; ++j){
        for(int i = -1; i < N+1; ++i){
            int index = (i+1) + (j+1)*(N+2);
            temp.first = i*dx;
            temp.second = j*dy;
            points[index] = temp;
        }
    }
    
    
    double h = dx;
	//Initialization of plevel_set_, x_avrg_num, r_avrg_num and den
    double* level_set = (double*) malloc(sizeof(double)*(N+2)*(M+2));
    int* is_fluid = (int*) malloc(sizeof(int)*(N+2)*(M+2));
    vector* x_avrg_num = (vector*) malloc(sizeof(vector)*N*M);
    double* r_avrg_num = (double*) malloc(sizeof(double)*N*M);
	double* den = (double*) malloc(sizeof(double)*N*M);
	
	//Compute the values of x_avrg_num, r_avrg_num and den
	for(int p = 0; p < NUMBER_OF_PARTICLE; ++p){
        if(all_particle[p].tag == interior){
            vector particle_pos = all_particle[p].position;
            vector init_cell = index_from_coord(particle_pos);
            is_fluid[(int)init_cell.first + (int)init_cell.second * (N+2)];
        
            //~ for(int j = (int)init_cell.second - 2; j <= (int)init_cell.second + 2; ++j){
                //~ for(int i = (int)init_cell.first - 2; i <= (int)init_cell.first + 2; ++i){
            for(int j = (int)init_cell.second - 1; j <= (int)init_cell.second + 1; ++j){
                for(int i = (int)init_cell.first - 1; i <= (int)init_cell.first + 1; ++i){
                    
                    if (j < 0 || j >= M || i < 0 || i >= N ){
                        continue;
                    }
                    
                    vector cell;
                    cell.first = i*dx;
                    cell.second = j*dy;
                    int index = i + j*N;
                    double temp = vec_distance_vec(vec_sub_vec(cell, particle_pos), vec_sub_vec(cell, particle_pos))/h;
                    
                    if (temp < 1){
                        double W_surf = (1-temp*temp)*(1-temp*temp);
                        vector temp1;
                        temp1.first = particle_pos.first*W_surf;
                        temp1.second = particle_pos.second*W_surf;
                        x_avrg_num[index] = vec_add_vec(x_avrg_num[index], temp1);
                        *(r_avrg_num + index) += W_surf*0.87*dx;
                        *(den + index) += W_surf;
                    }
                }
            }
        }
    }
        
    //Compute the values of level set function
    for(int j = -1; j < M+1; ++j){
        for(int i = -1; i < N+1; ++i){
            int index = (i+1) + (j+1)*(N+2);
            
            if(i == -1 || i == N || j == -1 || j == M){
                level_set[index] = 0.5*dx;
            }
            
            else{
                int index2 = i + j*N;
                double temp = *(den+index2);
                vector x_avrg;
                x_avrg.first = x_avrg_num[index2].first/temp;
                x_avrg.second = x_avrg_num[index2].second/temp;
                double r_avrg = 0.87*dx;
                vector cell_pos;
                cell_pos.first = i*dx;
                cell_pos.second = j*dx;
                
                if(*(den+index2) != 0)
                    level_set[index] = vec_distance_vec(vec_sub_vec(cell_pos, x_avrg),vec_sub_vec(cell_pos, x_avrg)) - r_avrg;
                else{
                    if(is_fluid[index2] == 1)
                        level_set[index] = -1.;
                    else
                        level_set[index] = 1.;
                }
            }
        }
    }
    
    //Output of level set function
	char file_name[22] = "level_set_function.csv";
	FILE *fp = NULL;
	fp = fopen(file_name,"w");
	fprintf(fp, "x coord, y coord, level set\n"); 
	for (int i = 0; i < (N+2)*(M+2); i++) {
		vector this_point = points[i];
		fprintf(fp, "%lf, %lf, %lf, \n",  
				this_point.first + x_min,
				this_point.second,
				level_set[i]);
	}
	fclose(fp);
    
    //Finding the height of the fluid at x = 0.05
    //First linear interpolation to find the values of the level set function at x = 0.05
    double* level_set_x005 = (double*) malloc(sizeof(double)*(M+2));
    for(int j = -1; j < M+1; ++j){
        level_set_x005[j+1] = level_set[(2+1)+(j+1)*(N+2)] + (0.05 - dx)*((level_set[(3+1)+(j+1)*(N+2)] - level_set[(2+1)+(j+1)*(N+2)])/(2*dx - dx));
    }
    //Then finding the zero of the function
    for(int j = -1; j < M+1; ++j){
        if(level_set_x005[j+2] - level_set_x005[j+1] != 0){
            double y0 = dy*j;
            double y1 = dy*(j+1);
            double y = y0 - level_set_x005[j+1]*(y1 - y0)/(level_set_x005[j+2] - level_set_x005[j+1]);
            if ((y >= y0 && y <= y1) && y > 0 && y < 1.1){
                printf("One 0 of the level set function is situated at y = %f (x = 0.05) \n", y);
            }
        }
    }
    
}
