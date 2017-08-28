#ifndef FILE_ADVECTION_H

#include <string>
#include "Dense"
#include "Sparse"

#include "Mesh2D.h"
#include "OdeSystem.h"

class Advection : public OdeSystem
{
 private:
  // Reference de la classe Mesh2D (accès à tout ce qui concerne le maillage
  // et les données géométriques)
  const Mesh2D& _mesh;

  // Variables géométriques
  const std::vector<Triangle>& _triangles;
  const std::vector<Edge>& _edges;
  const std::vector<Eigen::Matrix<double, 2, 1> >& _tri_center;
  const std::vector<double>& _tri_area;
  const std::vector<Eigen::Matrix<double, 2, 1> >& _edg_normal;
  const std::vector<Eigen::Matrix<double, 2, 1> >& _edg_center;

  // Le choix du flux
  std::string _type_flux;

  // le choix de l'ecoulement
  std::string _type_ecoulement;

  // condition initiale
  std::string _cond_init;

  // ecoulement V aux centres des triangles
  std::vector<Eigen::Matrix<double, 2, 1> > _V;

 public:
  // Constructeur
  Advection(const Mesh2D& mesh, std::string type_flux, std::string type_ecoulement);

  // calcul de V au centre des triangles
  void BuildVelocity(double t);

  // Construit le vecteur f = g(u) (EDO : du/dt = g(u))
  void BuildF(const double& t, const Eigen::VectorXd& sol);

  // Condition Initiale au centre des triangles
  Eigen::VectorXd InitialCondition(std::string initial_condition);

  // Solution exacte au centre des triangles
  Eigen::VectorXd ExactSolution(double t);

  // Erreur en norme L2
  double Norm(Eigen::VectorXd& approxSol,Eigen::VectorXd& exactSol);

  // Sauvegarde la solution
  void SaveSol(const Eigen::VectorXd& sol, const double& t, int n);

};

#define FILE_ADVECTION_H
#endif
