#ifndef FILE_TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur avec le systeme d'edo a resoudre
TimeScheme::TimeScheme(OdeSystem& sys) : _sys(sys)
{}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{}

// Initialisation de vos différentes variables
void TimeScheme::Initialize(double t0, double dt, Eigen::VectorXd& rho0)
{
  _dt = dt;
  _t = t0;
  _sol0 = rho0;
  _sol = _sol0;
}

// Renvoie _sol (pratique pour vérifier la résolution)
const VectorXd & TimeScheme::GetIterateSolution() const
{
  return _sol;
}

// Euler Explicite
void EulerScheme::Advance()
{
  _t += _dt;
  _sys.BuildF(_t, _sol);
  _sol += _dt*_sys.GetF();
}

// RungeKutta
void RungeKuttaScheme::Advance()
{
  _t += _dt;
  VectorXd k1, k2, k3, k4;
  _sys.BuildF(_t, _sol);
  k1 = _sys.GetF();
  _sys.BuildF(_t+_dt/2., _sol+_dt/2.*k1);
  k2 = _sys.GetF();
  _sys.BuildF(_t+_dt/2., _sol+_dt/2.*k2);
  k3 = _sys.GetF();
  _sys.BuildF(_t+_dt, _sol+_dt*k3);
  k4 = _sys.GetF();
  _sol += _dt/6.*(k1 + 2.*k2 + 2.*k3 + k4);
}

#define FILE_TIME_SCHEME_CPP
#endif
