#include "shape.h"
#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <Windows.h>
#include <GL\GL.h>

#include <GL\GLU.h>
#include <GL\glut.h>

using namespace std;
#define MAX_COLOR 13


bool shape::LoadObjShape(char* objshape)
{
	char buffer[512];

	ifstream file_in(objshape);

	if(!file_in.is_open())
	{
		cout<<"Failed to open file "<<objshape<<endl;
		return false;
	}

	int count_v = 0;
	int count_f = 0;
	while(!file_in.eof())
	{
		file_in.getline(buffer,512);
		if(buffer[0] == 'f')	count_f ++;
		if(buffer[0] == 'v')	count_v ++;
	}

	this->shape_face_number = count_f;
	this->shape_vertex_number = count_v;
	this->shape_edge_number = count_f * 3 / 2;

	this->vertexes = new double*[this->shape_vertex_number];
	this->vertex_color = new double[this->shape_vertex_number];
	for (int i = 0; i < this->shape_vertex_number; i++)
	{
		this->vertexes[i] = new double[3];
	}

	this->faces = new int*[this->shape_face_number];
	for (int i = 0; i < this->shape_face_number; i++)
	{
		this->faces[i] = new int[3];
		this->vertex_color[i] = 0;
	}

	file_in.close();
	file_in.open(objshape);

	count_v = 0;
	count_f = 0;
	while(!file_in.eof())
	{
		file_in.getline(buffer,512);

		if(buffer[0] == 'f')
		{
			int* p = this->faces[count_f];
			sscanf(buffer,"f %d %d %d",p , p + 1, p + 2);
			count_f ++;
		}

		if(buffer[0] == 'v')
		{
			double* p = this->vertexes[count_v];
			sscanf(buffer,"v %lf %lf %lf",p, p + 1, p + 2);
			count_v ++;
		}
	}

	//neighbor faces of one vertex
	this->vertex_neighbor_faces = new vector<int> [this->shape_vertex_number];
	list<int> *temp_edges;
	temp_edges = new list<int> [this->shape_vertex_number + 1];
	for (int i = 0; i < this->shape_face_number; i++)
	{
		int a = this->faces[i][0];
		int b =	this->faces[i][1];
		int c = this->faces[i][2];

		if (a < b){temp_edges[a].push_back(b);}else{temp_edges[b].push_back(a);};
		if (a < c){temp_edges[a].push_back(c);}else{temp_edges[c].push_back(a);};
		if (b < c){temp_edges[b].push_back(c);}else{temp_edges[c].push_back(b);};

		this->vertex_neighbor_faces[a-1].push_back(i);
		this->vertex_neighbor_faces[b-1].push_back(i);
		this->vertex_neighbor_faces[c-1].push_back(i);
	}

	this->edges = new int*[this->shape_edge_number];
	for (int i = 0; i < this->shape_edge_number; i++)
	{
		this->edges[i] = new int[2];
	}

	//neighbor vertexes of one vertex
	this->vertex_neighbor_vertexes = new vector<int>[this->shape_vertex_number];

	int count_e = 0;
	int *count_e_array = new int[this->shape_vertex_number + 1];
	for(int i = 1; i <= this->shape_vertex_number; i++){
		temp_edges[i].sort();
		count_e_array[i] = count_e;
		set<int> vset;
		for(list<int>::iterator iter=temp_edges[i].begin(); iter != temp_edges[i].end(); iter++){
			int edge_v = *iter;
			if (vset.find(edge_v) == vset.end()){
				this->edges[count_e][0] = i;
				this->edges[count_e][1] = edge_v;

				this->vertex_neighbor_vertexes[i-1].push_back(edge_v-1);
				this->vertex_neighbor_vertexes[edge_v-1].push_back(i-1);

				count_e++;
				vset.insert(edge_v);
			}
		}
	}

	this -> face_edges = new int*[this->shape_face_number];

	for(int i = 0; i < this->shape_face_number; i++)
	{
		this -> face_edges[i] = new int[3];
	}

	for(int i = 0;i < this->shape_face_number; i++){
		for(int k = 0;k < 3; k++){
			int a = this->faces[i][k];
			int b = this->faces[i][(k+1) % 3];
			if (a > b) {int tmp = a;a = b; b = tmp;};
			int count_b = 0;
			for(list<int>::iterator int_iter = temp_edges[a].begin();int_iter!= temp_edges[a].end();int_iter++){
				if (*int_iter == b)
					break;
				else
					count_b++;
			}
			this->face_edges[i][k] = count_e_array[a] + count_b / 2;
		}
	}

	file_in.close();
	CalculateAreaLenth();
	NormalizeShape();
	//get color
	this->face_colors = new double[this->shape_face_number];
	file_in.open("./281_color.txt");

	if(!file_in.is_open())
	{
		cout<< "Failed to open file 281_color.txt"<<endl;
		return false;
	}

	double original_point[] = {0,0,0};
	for(int i = 0; i < this->shape_edge_number;i++)
	{
		this->average_point_length += this->DistanceEu(this->vertexes[edges[i][0] - 1],original_point);
	}
	this->average_point_length /= this->shape_edge_number;
	cout<<"average edge length is "<<average_point_length<<endl;
	count_f = 0;
	while(!file_in.eof())
	{
		if(count_f >= shape_face_number) break;
		file_in.getline(buffer,512);
		sscanf_s(buffer,"%lf",&(this->face_colors[count_f]));

		/*int a = this->faces[count_f][0] - 1;
		int b = this->faces[count_f][1] - 1;
		int c = this->faces[count_f][2] - 1;
		this->vertex_color[a] = this->face_colors[count_f];
		this->vertex_color[b] = this->face_colors[count_f];
		this->vertex_color[c] = this->face_colors[count_f];*/
		count_f ++;
	}
	file_in.close();

	return true;
}

bool shape::LoadOffShape(char* offshape)
{
	char buffer[512];

	ifstream file_in(offshape);

	if(!file_in.is_open())
	{
		cout<<"Failed to open file "<<offshape<<endl;
		return false;
	}

	int count_v = 0;
	int count_f = 0;
	int line = 0;
	while(!file_in.eof())
	{
		file_in.getline(buffer,512);
		if(line == 1)
		{
			sscanf_s(buffer,"%d %d %d",&this->shape_vertex_number,&this->shape_face_number,&this->shape_edge_number);
			break;
		}
		line++;
	}

	this->vertexes = new double*[this->shape_vertex_number];
	this->vertex_color = new double[this->shape_vertex_number];
	for (int i = 0; i < this->shape_vertex_number; i++)
	{
		this->vertexes[i] = new double[3];
	}

	this->faces = new int*[this->shape_face_number];
	for (int i = 0; i < this->shape_face_number; i++)
	{
		this->faces[i] = new int[3];
		this->vertex_color[i] = 0;
	}

	file_in.close();
	file_in.open(offshape);

	count_v = 0;
	count_f = 0;
	line = 0;
	while(!file_in.eof())
	{
		file_in.getline(buffer,512);

		if(line > 1 && line <= (1+ shape_vertex_number))
		{
			double* p = this->vertexes[count_v];
			sscanf(buffer,"%lf %lf %lf",p, p + 1, p + 2);
			count_v ++;
		}

		if(line > 1 + shape_vertex_number)
		{
			int* p = this->faces[count_f];
			if(buffer[0] != '3') break;
			sscanf(buffer,"3 %d %d %d",p , p + 1, p + 2);
			this->faces[count_f][0] ++;
			this->faces[count_f][1] ++;
			this->faces[count_f][2] ++;
			count_f ++;
		}
		line++;
	}

	//neighbor faces of one vertex
	this->vertex_neighbor_faces = new vector<int> [this->shape_vertex_number];
	list<int> *temp_edges;
	temp_edges = new list<int> [this->shape_vertex_number + 1];
	for (int i = 0; i < this->shape_face_number; i++)
	{
		int a = this->faces[i][0];
		int b =	this->faces[i][1];
		int c = this->faces[i][2];

		if (a < b){temp_edges[a].push_back(b);}else{temp_edges[b].push_back(a);};
		if (a < c){temp_edges[a].push_back(c);}else{temp_edges[c].push_back(a);};
		if (b < c){temp_edges[b].push_back(c);}else{temp_edges[c].push_back(b);};

		this->vertex_neighbor_faces[a-1].push_back(i);
		this->vertex_neighbor_faces[b-1].push_back(i);
		this->vertex_neighbor_faces[c-1].push_back(i);
	}

	this->edges = new int*[this->shape_edge_number];
	for (int i = 0; i < this->shape_edge_number; i++)
	{
		this->edges[i] = new int[2];
	}

	//neighbor vertexes of one vertex
	this->vertex_neighbor_vertexes = new vector<int>[this->shape_vertex_number];

	int count_e = 0;
	int *count_e_array = new int[this->shape_vertex_number + 1];
	for(int i = 1; i <= this->shape_vertex_number; i++){
		temp_edges[i].sort();
		count_e_array[i] = count_e;
		set<int> vset;
		for(list<int>::iterator iter=temp_edges[i].begin(); iter != temp_edges[i].end(); iter++){
			int edge_v = *iter;
			if (vset.find(edge_v) == vset.end()){
				this->edges[count_e][0] = i;
				this->edges[count_e][1] = edge_v;

				this->vertex_neighbor_vertexes[i-1].push_back(edge_v-1);
				this->vertex_neighbor_vertexes[edge_v-1].push_back(i-1);

				count_e++;
				vset.insert(edge_v);
			}
		}
	}

	this -> face_edges = new int*[this->shape_face_number];

	for(int i = 0; i < this->shape_face_number; i++)
	{
		this -> face_edges[i] = new int[3];
	}

	for(int i = 0;i < this->shape_face_number; i++){
		for(int k = 0;k < 3; k++){
			int a = this->faces[i][k];
			int b = this->faces[i][(k+1) % 3];
			if (a > b) {int tmp = a;a = b; b = tmp;};
			int count_b = 0;
			for(list<int>::iterator int_iter = temp_edges[a].begin();int_iter!= temp_edges[a].end();int_iter++){
				if (*int_iter == b)
					break;
				else
					count_b++;
			}
			this->face_edges[i][k] = count_e_array[a] + count_b / 2;
		}
	}

	file_in.close();
	CalculateAreaLenth();
	NormalizeShape();
	//get color
	this->face_colors = new double[this->shape_face_number];
	file_in.open("./281_color.txt");

	if(!file_in.is_open())
	{
		cout<< "Failed to open file 281_color.txt"<<endl;
		return false;
	}

	double original_point[] = {0,0,0};
	for(int i = 0; i < this->shape_edge_number;i++)
	{
		this->average_point_length += this->DistanceEu(this->vertexes[edges[i][0] - 1],original_point);
	}
	this->average_point_length /= this->shape_edge_number;
	cout<<"average point length is "<<average_point_length<<endl;
	count_f = 0;
	while(!file_in.eof())
	{
		if(count_f >= shape_face_number) break;
		file_in.getline(buffer,512);
		sscanf_s(buffer,"%lf",&(this->face_colors[count_f]));

		/*int a = this->faces[count_f][0] - 1;
		int b = this->faces[count_f][1] - 1;
		int c = this->faces[count_f][2] - 1;
		this->vertex_color[a] = this->face_colors[count_f];
		this->vertex_color[b] = this->face_colors[count_f];
		this->vertex_color[c] = this->face_colors[count_f];*/
		count_f ++;
	}
	file_in.close();

	return true;
}
void shape::NormalizeShape()
{
	if(!this->normalized_faces)
	{
		this->normalized_faces = new double*[this->shape_face_number];
		for (int i = 0; i < this->shape_face_number; i++)
		{
			normalized_faces[i] = new double[3];
		}

		for (int i = 0; i < this->shape_face_number; i++)
		{
			int vertex_a = this->faces[i][0] - 1;
			int vertex_b = this->faces[i][1] - 1;
			int vertex_c = this->faces[i][2] - 1;

			double x1 = this->vertexes[vertex_a][0] - this->vertexes[vertex_c][0];
			double x2 = this->vertexes[vertex_b][0] - this->vertexes[vertex_c][0];
			double y1 = this->vertexes[vertex_a][1] - this->vertexes[vertex_c][1];
			double y2 = this->vertexes[vertex_b][1] - this->vertexes[vertex_c][1];
			double z1 = this->vertexes[vertex_a][2] - this->vertexes[vertex_c][2];
			double z2 = this->vertexes[vertex_b][2] - this->vertexes[vertex_c][2];

			double t1 = y1 * z2 - z1 * y2;
			double t2 = z1 * x2 - x1 * z2;
			double t3 = x1 * y2 - y1 * x2; 

			double norm_delta = sqrt(t1*t1+t2*t2+t3*t3);
			this->normalized_faces[i][0] = t1 / norm_delta;
			this->normalized_faces[i][1] = t2 / norm_delta;
			this->normalized_faces[i][2] = t3 / norm_delta;
		}
	}

	if(!this->normalized_vertexes)
	{
		this->normalized_vertexes = new double*[this->shape_vertex_number];
		for (int i = 0; i < this->shape_vertex_number; i++)
		{
			this->normalized_vertexes[i] = new double[3];
		}
	}

	for (int i = 0; i < this->shape_vertex_number; i++)
	{
		double t1 = 0;
        double t2 = 0;
        double t3 = 0;
        double sum_area = 0;

        for(int j=0;j<(this->vertex_neighbor_faces[i]).size();j++){
            int fid = this -> vertex_neighbor_faces[i][j];
            double area = face_area[fid];
            sum_area += area;

            t1 += area * this-> normalized_faces[fid][0];
            t2 += area * this-> normalized_faces[fid][1];
            t3 += area * this-> normalized_faces[fid][2];
        }
        
        double norm_delta = sqrt(t1*t1+t2*t2+t3*t3);

        this->normalized_vertexes[i][0] = t1 / norm_delta;
        this->normalized_vertexes[i][1] = t2 / norm_delta;
        this->normalized_vertexes[i][2] = t3 / norm_delta;
	}
}

void shape::CalculateAreaLenth(){
    if (! this->face_area) this->face_area = new double[this->shape_face_number];
    
    for(int i=0;i < this -> shape_face_number;i++){
        int vertex_a = this->faces[i][0]-1;
        int vertex_b = this->faces[i][1]-1;
        int vertex_c = this->faces[i][2]-1;
        
        double length_c = DistanceEu(this->vertexes[vertex_a], this->vertexes[vertex_b]);
        double length_b = DistanceEu(this->vertexes[vertex_a], this->vertexes[vertex_c]);
        double length_a = DistanceEu(this->vertexes[vertex_c], this->vertexes[vertex_b]);
        
        double r = (length_a + length_b + length_c)/2;
        double area = sqrt(r * (r - length_a) * (r - length_b) * (r - length_c) );
        this->face_area[i] = area;
    }

    if (! this->edge_length ) this->edge_length = new double[this->shape_edge_number];
    
    for(int i=0;i < this-> shape_edge_number;i++){
        int a = this->edges[i][0] - 1;
        int b = this->edges[i][1] - 1;
        this->edge_length[i] = DistanceEu(this->vertexes[a], this->vertexes[b]);
    }
}

double shape::DistanceEu(double* a,double* b)
{
		return sqrt(pow(a[0] - b[0],2) + pow(a[1] - b[1],2) + pow(a[2] - b[2],2));
	
}

void shape::Show(int mode){
    glColor3f(1,1,1);

	if(mode == 8 && this->noise_flag)
	{
		this->AddGaussianNoise(0.005);
		this->NormalizeShape();
		this->noise_flag = !this->noise_flag;
	}

	if(mode == 9 && this->denoise_flag)
	{
		BilateralDenoise(2, 0.02, 0.03, 0.02);
		this->denoise_flag = !this->denoise_flag;
	}

    for(int i=0;i < this->shape_face_number; i++){
		glColor3f(1,1,1);
        int a = this->faces[i][0]-1;
        int b = this->faces[i][1]-1;
        int c = this->faces[i][2]-1;
		if(mode == 0 && !this->stat_flag)
		{
			this->stat_flag = true;
			cout << "The model's face number:"<<this->shape_face_number<<endl;
			cout << "The model's vertex number:"<<this->shape_vertex_number<<endl;
			cout << "The model's edge number:"<<this->shape_edge_number<<endl; 
		}
		if(mode == 1)
		{
			double r,g,b;
			r = (double)i/this->shape_face_number;
			g = 0.6;
			b = 0.8;
			glColor3f(r,g,b);
		}

        if (mode == 2)
		{
			double r,g,b;
			r = int(this->face_colors[i]) % 3;
			g = int(this->face_colors[i]) % 5;
			b = int(this->face_colors[i]) % 7;
			glColor3f(r/3,g/5,b/7);
        }            
		if(mode == 4)
		{
			if(a == 1500 || b == 1500 ||c == 1500)	glColor3f(0,0,0);
		}
		if(mode == 5)
		{
			int face_id = 3500;
			if(FaceRegionOfFace(this->faces[i],this->faces[face_id],3))
			{
				glColor3f(0,0,0);
			}
		}
		if(mode == 6)
		{
			int face_region_ids[10] = {50,354,658,1257,5789,14281,20007,7864,5304,464};
			if(FaceRegion(this->faces[i],face_region_ids,10))
			{
				double color_r,color_g,color_b;
				color_r = int(this->face_colors[i]) % 3;
				color_g = int(this->face_colors[i]) % 5;
				color_b = int(this->face_colors[i]) % 7;
				glColor3f(color_r/3,color_g/5,color_b/7);
			}
		}
        glBegin(GL_TRIANGLES);
        glNormal3f(this->normalized_faces[i][0], this->normalized_faces[i][1],this->normalized_faces[i][2]); 
		if (this->vertex_color[a] > 0) glColor3f(0.4, 1 - vertex_color[a], 0.8);
        glVertex3f(this->vertexes[a][0],this->vertexes[a][1],this->vertexes[a][2]);   
		if (this->vertex_color[b] > 0) glColor3f(0.4, 1 - vertex_color[b], 0.8);
        glVertex3f(this->vertexes[b][0],this->vertexes[b][1],this->vertexes[b][2]);  
		if (this->vertex_color[c] > 0) glColor3f(0.4, 1 - vertex_color[c], 0.8);
        glVertex3f(this->vertexes[c][0],this->vertexes[c][1],this->vertexes[c][2]);
        glEnd();
    }

	if(mode == 3)
	{
		//pick a random vertex
		srand((UINT)GetCurrentTime());
		//int v_id = rand() % this->shape_vertex_number;
		int v_id = 1500;
		glPointSize(3.0f);
		glEnable(GL_POINT_SMOOTH); 
		glBegin(GL_POINTS);  
		glColor3f(1.0f,0,0);
		glVertex3f(this->vertexes[v_id][0],this->vertexes[v_id][1],this->vertexes[v_id][2]);    
		glEnd();		 
		for(int i = 0;i < this->vertex_neighbor_vertexes[v_id].size();i++)
		{
			int neighbor_id = this->vertex_neighbor_vertexes[v_id][i];
			glBegin(GL_POINTS);  
			glColor3f(0,0,1.0);       
			glVertex3f(this->vertexes[neighbor_id][0],this->vertexes[neighbor_id][1],this->vertexes[neighbor_id][2]);    
			glEnd();
		}
	}
	
	if(mode == 7)
	{
		int f_id =1500;
		glColor3f(255, 0, 0);
        glLineWidth(1);
        int a = this->faces[f_id][0]-1;
        int b = this->faces[f_id][1]-1;
        int c = this->faces[f_id][2]-1;
        
        double ox = (this->vertexes[a][0] + this->vertexes[b][0] + this->vertexes[c][0])/3;
        double oy = (this->vertexes[a][1] + this->vertexes[b][1] + this->vertexes[c][1])/3;
        double oz = (this->vertexes[a][2] + this->vertexes[b][2] + this->vertexes[c][2])/3;
        
        double length = 1;

        glBegin(GL_LINES);
            glVertex3f(ox, oy, oz);
            glVertex3f(ox + length * this->normalized_faces[f_id][0],
                        oy + length * this->normalized_faces[f_id][1],
                        oz + length * this->normalized_faces[f_id][2]);
        glEnd();
	}

	if(mode == 10 && !this->mse_flag)
	{
		double mse = this->CalculateMSE();
		cout<<"MSE value is "<<mse<<endl;
		this->mse_flag = true;
	}
    /*if (showLines){
        glColor3f(0, 0, 0);
        glLineWidth(0.1);
        glBegin(GL_LINES);
        for(int i=0;i<this-> shape_edge_number;i++){
            int a = this->edges[i][0]-1;
            int b = this->edges[i][1]-1;
            glVertex3f(this->vertexes[a][0],this->vertexes[a][1],this->vertexes[a][2]);
            glVertex3f(this->vertexes[b][0],this->vertexes[b][1],this->vertexes[b][2]);
        }
        glEnd();
    }*/

}

bool shape::FaceRegion(int* face,int* faces,int size)
{
	int count = 0;
	for(int i = 0;i < 3;i++)
	{
		for (int j = 0; j < size; j++)
		{
			if(face[i] == faces[j])	
			{
				count ++;
				break;
			}
		}
	}
	if(count >= 1) return true;
	return false;
}

bool shape::FaceRegionOfFace(int* face,int* faces,int size)
{
	int count = 0;
	for(int i = 0;i < 3;i++)
	{
		for (int j = 0; j < size; j++)
		{
			if(face[i] == faces[j])	
			{
				count ++;
				break;
			}
		}
	}
	if(count == 2) return true;
	return false;
}

double gaussrand(){
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    do {
        double U1 = (double)rand() / RAND_MAX;
        double U2 = (double)rand() / RAND_MAX;
             
        V1 = 2 * U1 - 1;
        V2 = 2 * U2 - 1;
        S = V1 * V1 + V2 * V2;
    } while(S >= 1 || S == 0);         
    X = V1 * sqrt(-2 * log(S) / S);    
 
    return X;
}

void shape::AddGaussianNoise(double delta){
    for(int i = 0;i< this->shape_vertex_number;i++){
        double d = gaussrand() * delta;
        this->vertexes[i][0] +=  d * this->normalized_vertexes[i][0];
        this->vertexes[i][1] +=  d * this->normalized_vertexes[i][1];
        this->vertexes[i][2] +=  d * this->normalized_vertexes[i][2];
    }
}

void shape::BilateralDenoise(int iter_time, double sigma_c, double sigma_s, double radius){
    
    double **tmp_vertex = new double* [this->shape_vertex_number];
	for(int i=0;i<this->shape_vertex_number;i++)
	{ 
		tmp_vertex[i] = new double[3];
	}
	this->original_vertex = new double* [this->shape_vertex_number];
	for(int i=0;i<this->shape_vertex_number;i++)
	{ 
		original_vertex[i] = new double[3];
		original_vertex[i][0] = this->vertexes[i][0];
		original_vertex[i][1] = this->vertexes[i][1];
		original_vertex[i][2] = this->vertexes[i][2];
	}
    /*ÖÐÖµÂË²¨
	for(int i=0;i<this->shape_vertex_number;i++)
	{ 
		tmp_vertex[i] = new double[3];
		tmp_vertex[i][0] = this->vertexes[i][0];
		tmp_vertex[i][1] = this->vertexes[i][1];
		tmp_vertex[i][2] = this->vertexes[i][2];
	}

	for(int i = 0; i < this->shape_vertex_number;i++)
	{
		double t_x,t_y,t_z;
		vector<int> neighbor_vertex = vertex_neighbor_vertexes[i];
		t_x = t_y = t_z = 0;
		int j;
		for (j = 0; j < this->vertex_neighbor_vertexes[i].size(); j++)
		{
			int neighbor = neighbor_vertex[j];
			t_x += this->vertexes[neighbor][0];
			t_y += this->vertexes[neighbor][1];
			t_z += this->vertexes[neighbor][2];
		}
		vertexes[i][0] = (vertexes[i][0] + t_x)/j;
		vertexes[i][1] = (vertexes[i][1] + t_y)/j;
		vertexes[i][2] = (vertexes[i][2] + t_z)/j;
	}
    this->NormalizeShape();*/  
    for(int iter=0;iter < iter_time;iter++){        
        
        vector<int> neighbor_vertex;
		set<int> neighbor_set;
        for(int i=0;i<this->shape_vertex_number;i++)
		{            
            neighbor_vertex = this->vertex_neighbor_vertexes[i];
			neighbor_set.clear();
            int K = neighbor_vertex.size();
			for (int s = 0; s < K; s++)
			{
				int j = neighbor_vertex[s];
				if(this->DistanceEu(this->vertexes[i],this->vertexes[j]) < 2 * sigma_c) neighbor_set.insert(j);
				vector<int> tmp_neighbor = this->vertex_neighbor_vertexes[j];
				for (int l = 0; l < tmp_neighbor.size(); l++)
				{
					if(!this->DistanceEu(this->vertexes[i],this->vertexes[tmp_neighbor[l]]) < 2 * sigma_c)	continue;
					neighbor_set.insert(tmp_neighbor[l]);
				}
			}
			K = neighbor_set.size();
            double sum = 0; double normalizer = 0;
			set<int>::iterator it;
            for( it = neighbor_set.begin();it != neighbor_set.end();it ++){               
			//for( int k = 0;k < K; k++){      
                int j = *it;
				//int j = neighbor_vertex[k];

                double t = DistanceEu(this->vertexes[i],this->vertexes[j]);
				double temp_h = 0;
				for (int l = 0; l < 3; l++)
				{
					temp_h += (this->vertexes[i][l] - this->vertexes[j][l]) * normalized_vertexes[i][l];
				}
                double h = temp_h;

                double wc = exp(-0.5*t*t/(sigma_c *sigma_c));
                double ws = exp(-0.5*h*h/(sigma_s *sigma_s));
                sum -= wc * ws * h;
                normalizer += wc * ws;
            }
            double d = sum / normalizer;
            
            tmp_vertex[i][0] = this->vertexes[i][0] + this->normalized_vertexes[i][0] * d;
            tmp_vertex[i][1] = this->vertexes[i][1] + this->normalized_vertexes[i][1] * d;
            tmp_vertex[i][2] = this->vertexes[i][2] + this->normalized_vertexes[i][2] * d;
            
        }
        for(int i = 0;i < this-> shape_vertex_number; i++){
            this->vertexes[i][0] = tmp_vertex[i][0];
            this->vertexes[i][1] = tmp_vertex[i][1];
            this->vertexes[i][2] = tmp_vertex[i][2];
        }

        // update face normals and vector normals
        this->NormalizeShape();     
		cout<<"iteration: "<<iter + 1<<endl;
    }
	cout<<"denoised!"<<endl;
}

vector<int> shape::GetNeighbor(int v_id,double rings)
{
	vector<int> neighbors;
	for (int i = 0; i < this->shape_vertex_number; i++)
	{
		if(DistanceEu(this->vertexes[i],this->vertexes[v_id]) < rings)
		{
			neighbors.push_back(i);
		}
	}
	return neighbors;
}

double shape::CalculateMSE()
{
	double mse = 0.0;
    for(int i=0;i < this->shape_vertex_number;i++){
        mse += pow(this->DistanceEu(this->vertexes[i], this->original_vertex[i]),2);
    }
    
    mse = mse / (double)this->shape_vertex_number;	

	double max_d = 0;
    for(int i=0;i< this->shape_vertex_number;i++){
        double d =  this->DistanceEu(this->vertexes[i], this->original_vertex[i]);
        if(d > max_d)	max_d = d;
    }
	for(int i=0;i< this->shape_vertex_number;i++){
        double d =  this->DistanceEu(this->vertexes[i], this->original_vertex[i]);
        this->vertex_color[i] = d / max_d;
    }
	return mse;
}