#ifndef SHAPE_H
#define SHAPE_H
#include <vector>

using namespace std;

class shape
{
public:
	shape(){};
	shape(char* objshape);
	~shape(){};

	int shape_vertex_number;
	int shape_face_number;
	int shape_edge_number;

	double** vertexes;//Matrix vertex_number * 3
	double** original_vertex;
	double*  vertex_color;//Matrix vertex_number * 1
	int**	 faces;//Matrix face_number * 3
	double*  face_colors;//Matrix face_number * 1
	int**	 edges;//Matrix edge_number * 2
	int**	 face_edges;//Matrix face_number * 3
	double** normalized_faces;
	double** normalized_vertexes;
	bool	noise_flag;
	bool	denoise_flag;
	bool	stat_flag;
	bool	mse_flag;
	double	average_point_length;

	double* face_area;
	double* edge_length;

	vector<int> *vertex_neighbor_faces;
	vector<int> *vertex_neighbor_vertexes;

	void shape::NormalizeShape();
	bool LoadObjShape(char*);
	bool LoadOffShape(char*);
	void CalculateAreaLenth();
	double DistanceEu(double*,double*);
	void DistanceD();
	bool FaceRegion(int*,int*,int);
	bool FaceRegionOfFace(int*,int*,int);
	void DrawVertexNeighborVertexes(int vertexId);
	void DrawVertexNeighborFaces(int vertexId);
	void Show(int);//1:read file;2:shade mesh shape;3:find vertex neiborhood vertexes
	void AddGaussianNoise(double);
	void BilateralDenoise(int,double,double,double);
	double CalculateMSE();
	vector<int> GetNeighbor(int,double);
};
#endif