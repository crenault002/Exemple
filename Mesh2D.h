#ifndef FILE_MESH_2D_H

#include<vector>
#include<string>
#include "Dense"
#include "Sparse"

//! Arete du maillage
class Edge
{
 private:
  // numeros des deux sommets de l'arete
  Eigen::Matrix<int, 2, 1> _edg;

  // reference de l'arete
  int _ref;

  // numeros des triangle T1 et T2
  int _t1, _t2;

 public:
  Edge();
  Edge(int a, int b, int ref);

  void Print() const;

  inline const Eigen::Matrix<int, 2, 1>& GetEdge() const { return _edg;}
  inline int GetT1() const { return _t1; }
  inline int GetT2() const { return _t2; }
  inline int GetReference() const { return _ref; }

  inline void AddElement(int t)
  {
    if (_t1 == -1)
      _t1 = t;
    else
      _t2 = t;
  }

};

inline std::ostream& operator<<(std::ostream& out, const Edge& edge)
{ edge.Print(); }

//! triangle du maillage
class Triangle
{
 private:
  // numeros des trois sommets du triangle
  Eigen::Matrix<int, 3, 1> _tri;

  // arete du triangle
  int _ref;

 public:
  Triangle();
  Triangle(int a, int b, int c, int ref);

  void Print() const;

  inline const Eigen::Matrix<int, 3, 1>& GetTriangle() const { return _tri; }

};

inline std::ostream& operator<<(std::ostream& out, const Triangle& tri)
{ tri.Print(); }

//! classe contenant le maillage
class Mesh2D
{
 private:
  // liste de tous les sommets
  std::vector<Eigen::Matrix<double, 2, 1> > _vertices;

  // liste de tous les triangles
  std::vector<Triangle> _triangles;

  // centre de tous les triangles
  std::vector<Eigen::Matrix<double, 2, 1> > _tri_center;

  // aire de tous les triangles
  std::vector<double> _tri_area;

  // liste de toutes les aretes
  std::vector<Edge> _edges;

  // liste de toutes les normales
  std::vector<Eigen::Matrix<double, 2, 1> > _edg_normal;

  // centre des aretes
  std::vector<Eigen::Matrix<double, 2, 1> > _edg_center;

 public:
  Mesh2D();

  void BuildTrianglesCenter();
  void BuildEdgesNormal();

 protected:
  void AddSingleEdge(const Edge& edge, int ne, std::vector<int>& head_minv,
		     std::vector<int>& next_edge, int& nb_edges);

 public:
  void Read(std::string name_mesh);

  inline const std::vector<Eigen::Matrix<double, 2, 1> >& GetVertices() const {return _vertices;}

  inline const std::vector<Triangle>& GetTriangles() const { return _triangles; }
  inline const std::vector<Eigen::Matrix<double, 2, 1> >& GetTrianglesCenter() const {return _tri_center; }
  inline const std::vector<double>& GetTrianglesArea() const {return _tri_area; }

  inline const std::vector<Edge>& GetEdges() const { return _edges; }
  inline const std::vector<Eigen::Matrix<double, 2, 1> >& GetEdgesNormal() const { return _edg_normal; }
  inline const std::vector<Eigen::Matrix<double, 2, 1> >& GetEdgesCenter() const { return _edg_center; }

};

#define FILE_MESH_2D_H
#endif
