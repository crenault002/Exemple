#ifndef FILE_ODE_SYSTEM_CPP

#include "OdeSystem.h"

using namespace Eigen;

// Constructeur par défaut
OdeSystem::OdeSystem()
{
}
    
// Destructeur par défaut
OdeSystem::~OdeSystem()
{
}
    

// Pour récupérer _f
const VectorXd&  OdeSystem::GetF()  const
{
  return _f;
}

#define FILE_ODE_SYSTEM_CPP
#endif

