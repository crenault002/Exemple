#include <iostream>
#include <fstream>
#include <chrono>

#include "Advection.h"
#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
  if (argc < 8)
    {
      cout << "Please, enter the following arguments :." << endl;
      cout << "mesh_file final_time time_scheme dt initial_condition flux velocity" << endl;
      cout << "Example :" << endl;
      cout << "./run Meshes/square_0_02.mesh 10.0 rk4 0.02 gaussian centered uniform" << endl;
      cout << "Current " << argc << " arguments are " << endl;
      for (int i = 0; i < argc; i++)
	cout << "Argument " << i << " = " << argv[i] << endl;

      abort();
    }

  // on recupere les parametres d'entree
  string mesh_file(argv[1]);
  double tf = atof(argv[2]);
  string name_time_scheme(argv[3]);
  double dt = atof(argv[4]);
  string initial_condition(argv[5]);
  string type_flux(argv[6]);
  string type_ecoulement(argv[7]);

  // on les affiche
  cout << endl << endl;
  cout << "Simulation sur maillage " << mesh_file << " avec des flux " << type_flux << endl;
  cout << "Jusqu'au temps " << tf << " avec le schema " << name_time_scheme << " et un pas de temps " << dt << endl;
  cout << "Condition initiale " << initial_condition << " avec une vitesse " << type_ecoulement << endl;
  cout << endl;

  // nombre d'iterations en temps et dt modifie pour tomber pile sur tf
  int nb_iter = ceil(tf / dt);
  dt = tf / nb_iter;

  // on construit le maillage
  Mesh2D mesh;
  mesh.Read(mesh_file);

  // on initialise le probleme d'advection
  Advection adv(mesh, type_flux, type_ecoulement);

  // on initialise le schema en temps
  TimeScheme* time_scheme;

  if (name_time_scheme == "rk4")
    time_scheme = new RungeKuttaScheme(adv);
  else if (name_time_scheme == "euler")
    time_scheme = new EulerScheme(adv);

  // on calcule la condition initiale
  VectorXd sol0 = adv.InitialCondition(initial_condition);
  adv.SaveSol(sol0, 0.0, 0);

  double error2(1);
  int k(0);
  VectorXd approxSol;
  VectorXd exactSol;
  auto start = chrono::high_resolution_clock::now();
  auto finish = chrono::high_resolution_clock::now();
  double t;

  time_scheme->Initialize(0.0, dt, sol0);

  if (((initial_condition == "gaussian")||(initial_condition == "window"))&&(type_ecoulement == "uniform")){
    start = chrono::high_resolution_clock::now();
    while (error2 > 0.5){
      time_scheme->Advance();
      adv.SaveSol(time_scheme->GetIterateSolution(), (k+1)*dt, k+1);
      approxSol = time_scheme->GetIterateSolution();
      exactSol = adv.ExactSolution((k+1)*dt);
      error2 = (approxSol-exactSol).norm();
      k++;
    }
    finish = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
    cout << "Compilation time to reach 10 percent error = " << t << " microseconds" << endl;
    cout << "------------------------------------" << endl;
  }
  else{
    start = chrono::high_resolution_clock::now();
    finish = chrono::high_resolution_clock::now();
    t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
    while (t < 6*pow(10,7)){
      time_scheme->Advance();
      adv.SaveSol(time_scheme->GetIterateSolution(), (k+1)*dt, k+1);
      approxSol = time_scheme->GetIterateSolution();
      exactSol = adv.ExactSolution((k+1)*dt);
      error2 = (approxSol-exactSol).norm();
      k++;
      finish = chrono::high_resolution_clock::now();
      t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
    }
    cout << "Error = " << error2 << endl;
    cout << "------------------------------------" << endl;
  }

  return 0;
}
