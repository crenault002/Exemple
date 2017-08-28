#ifndef FILE_ADVECTION_CPP

#include "Advection.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;

double square(double x)
{
  return x*x;
}

// Constructeur
Advection::Advection(const Mesh2D& mesh, string type_flux, string type_ecoulement)
  : _mesh(mesh), _triangles(mesh.GetTriangles()), _edges(mesh.GetEdges()),
    _tri_center(mesh.GetTrianglesCenter()), _tri_area(mesh.GetTrianglesArea()),
    _edg_normal(mesh.GetEdgesNormal()), _edg_center(mesh.GetEdgesCenter()),
    _type_flux(type_flux), _type_ecoulement(type_ecoulement)
{
}


// Construit le tableau contenant V pour tous les centres des triangles
void Advection::BuildVelocity(double t)
{
  _V.resize(_triangles.size());

  for (int k=0; k < _triangles.size(); k++){   // Boucle sur les triangles
    if ( _type_ecoulement == "uniform"){
        _V[k](0) = 0.8;
        _V[k](1) = 0.4;
    }
    else if ( _type_ecoulement == "rotational"){
        _V[k](0) = _tri_center[k](1)/5.;
        _V[k](1) = - _tri_center[k](0)/5.;
    }
    else if ( _type_ecoulement == "sinusoidal_esp"){
        _V[k](0) = 0.8*cos(2.*3.1415*_tri_center[k](1));
        _V[k](1) = 0.;
    }
    else if ( _type_ecoulement == "sinusoidal_tps"){
        _V[k](0) = 0;
        _V[k](1) = 0.5*cos(2.*3.1415*t);
    }
  }

}


// Construit _f = g(t, U)
void Advection::BuildF(const double& t, const Eigen::VectorXd& sol)
{
  // Calcule v pour tous les triangles
  BuildVelocity(t);

  _f.resize(_triangles.size());
  _f.setZero();

  int t1,t2; //numéros des triangles 1 et 2 de l'arrête considérée
  double u1,u2,u_donnee; //Valeur de u sur les triangles 1 et 2

  for (int k=0; k<_edges.size();k++)
  {
    t1 = _edges[k].GetT1();
    t2 = _edges[k].GetT2();
    u1 = sol(t1);
    double y=_edg_center[k](1);
    u_donnee = exp(-6*(t-3)*(t-3))*exp(-6*(y+2)*(y+2));

    if (t2 == -1)
    {//On se situe sur une arrête du bord
      if (_V[t1].dot(_edg_normal[k]) > 0.){
        _f(t1) = _f(t1) - _V[t1].dot(_edg_normal[k])*u1/_tri_area[t1];
      }
      else if (_cond_init == "nulle")
      {
        if (abs(_edg_center[k](0)+4.0) < 1e-10)
          _f(t1) -= _V[t1].dot(_edg_normal[k])*u_donnee/_tri_area[t1];
      }
    }

    else
    {//L'arête n'est pas sur le bord
      u2 = sol(t2);
      if (_type_flux == "centre")
        {
          _f(t1) -= (0.5/_tri_area[t1])*(_V[t1].dot(_edg_normal[k])*u1 + _V[t2].dot(_edg_normal[k])*u2);
          _f(t2) += (0.5/_tri_area[t2])*(_V[t1].dot(_edg_normal[k])*u1 + _V[t2].dot(_edg_normal[k])*u2);
        }
      else if (_type_flux == "decentre")
      {
        if (_V[t1].dot(_edg_normal[k]) > 0)
        {
          _f(t1) -= (_V[t1].dot(_edg_normal[k])*u1)/_tri_area[t1];
          _f(t2) += (_V[t1].dot(_edg_normal[k])*u1)/_tri_area[t2];
        }
        else
        {
          _f(t1) -= (_V[t2].dot(_edg_normal[k])*u2)/_tri_area[t1];
          _f(t2) += (_V[t1].dot(_edg_normal[k])*u2)/_tri_area[t2];
        }
    }
  }
}
}


// Construit la condition initiale au centre des triangles
VectorXd Advection::InitialCondition(string initial_condition)
{
  _cond_init = initial_condition;
  VectorXd sol0(_triangles.size());

  //Point initial
  double x0, y0;
  x0=-2.;
  y0=-1.;

  if(_cond_init == "gaussian")
  {
    for (int i = 0; i<_triangles.size(); i++)
    {
      sol0(i)=exp(-7*((_tri_center[i](0)-x0)*(_tri_center[i](0)-x0)+(_tri_center[i](1)-y0)*(_tri_center[i](1)-y0)));
    }
  }
  else if (_cond_init == "creneau")
  {
    sol0.setZero();
    for (int i =0; i<_triangles.size(); i++)
    {
      if ((_tri_center[i](0)-x0)*(_tri_center[i](0)-x0)+(_tri_center[i](1)-y0)*(_tri_center[i](1)-y0)<=1)
      {
        sol0(i)=1.;
      }
    }
  }
  else if (_cond_init == "nulle")
  {
    sol0.setZero();
  }
  return sol0;
}


// Solution exacte au centre des triangles
VectorXd Advection::ExactSolution(double t)
{
  VectorXd exact_sol(_triangles.size());
  double x0(-2.);
  double y0(-1.);

  for(int k = 0; k < _triangles.size(); k++)
  {
    if((_cond_init == "gaussian")||(_cond_init == "creneau"))
      exact_sol[k] = exp(-7*(square(_tri_center[k](0)-_V[k](0)*t-x0)+square(_tri_center[k](1)-_V[k](1)*t-y0)));
  }

  return exact_sol;
}

// Erreur en norme L2
double Advection::Norm(VectorXd& sol, VectorXd& exact_sol)
{
  double norm(0.);
  int nb_tri = _triangles.size();

  for (int i = 0 ; i < nb_tri ; ++i)
    norm+= _tri_area[i]*square(sol[i]-exact_sol[i]);

  return sqrt(norm);
}


// Sauvegarde la solution
void Advection::SaveSol(const VectorXd& sol, const double& t, int n)
{
  string name_file = "solution_" + std::to_string(n) + ".vtk";

  const vector<Matrix<double, 2, 1> >& vertices = _mesh.GetVertices();
  int nb_vert = vertices.size();

  assert((sol.size() == _triangles.size()) && "The size of the solution vector is not the same than the number of _triangles !");


  ofstream solution;
  solution.open(name_file, ios::out);
  solution.precision(7);

  if (!solution.is_open())
    {
      cout << "Unable to open file " << name_file << endl;
      abort();
    }

  solution << "# vtk DataFile Version 3.0 " << '\n';
  solution << "2D Unstructured Grid" << '\n';
  solution << "ASCII" << '\n';
  solution << "DATASET UNSTRUCTURED_GRID" << '\n';

  solution << "POINTS " << nb_vert << " float " << '\n';
  for (int i = 0 ; i < nb_vert ; ++i)
    solution << vertices[i][0] << " " << vertices[i][1] << " 0." << '\n';

  solution << '\n';

  int nb_tri = _triangles.size();
  solution << "CELLS " << nb_tri << " " << nb_tri*4 << '\n';
  for (int i = 0 ; i < nb_tri ; ++i)
    solution << 3 << " " << _triangles[i].GetTriangle()[0] << " " << _triangles[i].GetTriangle()[1]
	     << " " << _triangles[i].GetTriangle()[2] << '\n';

  solution << '\n';

  solution << "CELL_TYPES " << nb_tri << '\n';
  for (int i = 0 ; i < nb_tri ; ++i)
    solution << 5 << '\n';

  solution << '\n';

  solution << "CELL_DATA " << nb_tri << '\n';
  solution << "SCALARS sol float 1" << '\n';
  solution << "LOOKUP_TABLE default" << '\n';
  for (int i = 0 ; i < nb_tri ; ++i)
  {
    solution << float(sol[i]) << '\n';
  }

  solution << '\n';

  solution.close();
}

#define FILE_ADVECTION_CPP
#endif
