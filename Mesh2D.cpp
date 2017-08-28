#ifndef FILE_MESH_2D_CPP

#include "Mesh2D.h"
#include<iostream>
#include<fstream>
#include <assert.h>

using namespace std;
using namespace Eigen;


/********
 * Edge *
 ********/


// constructeur par defaut
Edge::Edge()
{
  _edg[0] = -1;
  _edg[1] = -1;
  _ref = 0;
  _t1 = -1;
  _t2 = -1;
}

// constructeur prenant en arguments les deux sommets de l'arete et sa reference
Edge::Edge(int a, int b, int ref) : _ref(ref)
{
  // sort
  if (a > b)
    {
      _edg[0] = b;
      _edg[1] = a;
    }
  else
    {
      _edg[0] = a;
      _edg[1] = b;
    }

  _t1 = -1;
  _t2 = -1;
}


// on imprime l'arete
void Edge::Print() const
{
  cout << "[pt1, pt2] = [" << _edg[0] << " " << _edg[1] << "];" << endl;
  cout << "[T1, T2] = [" << _t1 << " " << _t2 << "];" << endl;
  cout << "ref = " << _ref << endl;
}


/************
 * Triangle *
 ************/


// constructeur par defaut
Triangle::Triangle()
{
  _tri[0] = -1; _tri[1] = -1; _tri[2] = -1;
  _ref = 0;
}


// constructeur prenant en arguments les trois sommets du triangle et sa reference
Triangle::Triangle(int a, int b, int c, int ref) : _ref(ref)
{
  _tri[0] = a; _tri[1] = b; _tri[2] = c;
  _ref = 0;
}


// on imprime le triangle
void Triangle::Print() const
{
  cout << "[pt1, pt2, pt3] = [" << _tri[0] << " " << _tri[1] << " " << _tri[2] << "];" << endl;
  cout << "ref = " << _ref << endl;
}


/**********
 * Mesh2D *
 **********/


// constructeur par defaut
Mesh2D::Mesh2D()
{
}


// methode construisant le centre des triangles et leur aire
void Mesh2D::BuildTrianglesCenter()
{
  _tri_center.resize(_triangles.size());
  _tri_area.resize(_triangles.size());

  int p1,p2,p3; //Définition des sommets
  double s1,s2;
  for (int i =0; i<_triangles.size(); i++) //Boucle sur les triangles
  {
    //Récupération des sommets
    p1 = _triangles[i].GetTriangle()(0);
    p2 = _triangles[i].GetTriangle()(1);
    p3 = _triangles[i].GetTriangle()(2);

    //Coordonnées du centre
    _tri_center[i](0) = (_vertices[p1](0)+_vertices[p2](0)+_vertices[p3](0))/3.;
    _tri_center[i](1) = (_vertices[p1](1)+_vertices[p2](1)+_vertices[p3](1))/3.;

    s1=(_vertices[p2](0)-_vertices[p1](0))*(_vertices[p3](1)-_vertices[p1](1));
    s2=(_vertices[p2](1)-_vertices[p1](1))*(_vertices[p3](0)-_vertices[p1](0));
    //Aire du triangle
    _tri_area[i] = 0.5*abs(s1-s2);

  }
}


// methode construisant les normales des aretes
void Mesh2D::BuildEdgesNormal()
{
  _edg_center.resize(_edges.size());
  _edg_normal.resize(_edges.size());
  int A1,A2,T1;
  double ps;
  vector<double> w(2);
  for (int i =0; i<_edges.size(); i++)
  {
    //Récupération des arêtes
    A1 = _edges[i].GetEdge()(0);
    A2 = _edges[i].GetEdge()(1);
    //Centre des arrêtes
    _edg_center[i](0)=0.5*(_vertices[A1](0)+_vertices[A2](0));
    _edg_center[i](1)=0.5*(_vertices[A1](1)+_vertices[A2](1));

    //Coordonnées des normales
    _edg_normal[i](0)=_vertices[A2](1)-_vertices[A1](1);
    _edg_normal[i](1)=_vertices[A1](0)-_vertices[A2](0);

    //Vérification de l'orientation des normales
    T1 = _edges[i].GetT1();
    w[0]=_edg_center[i](0)-_tri_center[T1](0);
    w[1]=_edg_center[i](1)-_tri_center[T1](1);

    if (_edg_normal[i](0)*w[0]+_edg_normal[i](1)*w[1] < 0)
    {
      _edg_normal[i](0)=-1.*_edg_normal[i](0);
      _edg_normal[i](1)=-1.*_edg_normal[i](1);
    }
  }
}


// methode interne qui rajoute une arete
void Mesh2D::AddSingleEdge(const Edge& edge, int ne, vector<int>& head_minv,
			   vector<int>& next_edge, int& nb_edges)
{
  int n1 = edge.GetEdge()(0);
  int n2 = edge.GetEdge()(1);
  int ref = edge.GetReference();

  bool existe = false;
  // we look at the list of edges leaving from n1
  // if we find the same edge than n1->n2 we add the edge
  for (int e = head_minv[n1]; e != -1; e = next_edge[e])
    {
      if (_edges[e].GetEdge()(1) == n2)
	{
	  if (ne >= 0)
	    _edges[e].AddElement(ne);

	  existe = true;
	}
    }

  // if the edge has not been found, we create it
  if (!existe)
    {
      // we initialize the edge
      _edges[nb_edges] = Edge(n1, n2, ref);
      if (ne >= 0)
	_edges[nb_edges].AddElement(ne);

      // we update the arrays next_edge and head_minv
      next_edge[nb_edges] = head_minv[n1];
      head_minv[n1] = nb_edges;
      nb_edges++;
    }
}


// lecture du maillage et creation des aretes, et calcul des aires, centres, normales
void Mesh2D::Read(string name_mesh)
{
  ifstream mesh_file(name_mesh.data());
  if (!mesh_file.is_open())
    {
      cout << "Unable to open file " << name_mesh << endl;
      abort();
    }
  else
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Reading mesh: " << name_mesh << endl;
    }

  string file_line;

  vector<Edge> edges_boundary;
  while (!mesh_file.eof())
    {
      getline(mesh_file, file_line);
      if (file_line.find("Vertices") != std::string::npos)
	{
	  int nb_vertices(0);
	  mesh_file >> nb_vertices;
	  cout << "Number of vertices  (" << nb_vertices << ")" << endl;

	  _vertices.resize(nb_vertices);
	  double x, y, z; int ref;
	  for (int i = 0 ; i < nb_vertices ; ++i)
	    {

	      mesh_file >> x >> y >> z >> ref;
	      _vertices[i](0) = x;
	      _vertices[i](1) = y;
	    }
	}
      else if (file_line.find("Edges") != std::string::npos)
	{
	  int nb_edges(0);
	  mesh_file >> nb_edges;
	  cout << "Number of edges (" << nb_edges << ")" << endl;

	  edges_boundary.resize(nb_edges);
	  int n1, n2, ref;
	  for (int i = 0 ; i < nb_edges ; ++i)
	    {
	      mesh_file >> n1 >> n2 >> ref;
	      n1--; n2--;
	      edges_boundary[i] = Edge(n1, n2, ref);
	    }
	}
      else if (file_line.find("Triangles") != std::string::npos)
	{
	  int nb_triangles(0);
	  mesh_file >> nb_triangles;
	  cout << "Number of triangles (" << nb_triangles << ")" << endl;

	  _triangles.resize(nb_triangles);
	  int n1, n2, n3, ref;
	  for (int i = 0 ; i < nb_triangles ; ++i)
	    {
	      mesh_file >> n1 >> n2 >> n3 >> ref;
	      n1--; n2--; n3--;
	      _triangles[i] = Triangle(n1, n2, n3, ref);
	    }
	}
    }

  cout << "---------Edges and Associated Triangles----------" << endl;

  // Toutes les aretes exterieures du maillage sont presentes
  int nb_edges = (3*_triangles.size() + edges_boundary.size())/2;
  _edges.resize(nb_edges);

  int nb_vertices = _vertices.size();
  vector<int> head_minv(nb_vertices, -1);
  vector<int> next_edge(nb_edges, -1);

  // on rajoute d'abord les aretes du bord
  nb_edges = 0;
  for (int i = 0; i < edges_boundary.size(); i++)
    this->AddSingleEdge(edges_boundary[i], -1, head_minv, next_edge, nb_edges);

  // ensuite les aretes interieures
  for (int i = 0; i < _triangles.size(); i++)
    {
      const Matrix<int, 3, 1>& nv = _triangles[i].GetTriangle();
      for (int j = 0; j < 3; j++)
	{
	  Edge edge(nv(j), nv((j+1)%3), 0);
	  AddSingleEdge(edge, i, head_minv, next_edge, nb_edges);
	}
    }

  cout << "-----------Triangles center and area-------------" << endl;
  BuildTrianglesCenter();

  cout << "------------Edges Normal --------------" << endl;
  BuildEdgesNormal();

  cout << "-------------------------------------------------" << endl;
}

#define FILE_MESH_2D_CPP
#endif
