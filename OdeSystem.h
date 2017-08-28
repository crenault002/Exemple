#ifndef FILE_ODE_SYSTEM_H

#include "Dense"
#include <fstream>

class OdeSystem
{
  protected:
    // vecteur f
    Eigen::VectorXd _f;
    
  public:
    // Constructeur par défaut
    OdeSystem();
    
    // Destructeur par défaut
    virtual ~OdeSystem();
    
    // Pour récupérer _f
    const Eigen::VectorXd&  GetF()  const;

    // Pour construire _f en fonction de votre système
    virtual void BuildF(const double& t, const Eigen::VectorXd& rho) = 0;
    
};

#define FILE_ODE_SYSTEM_H
#endif

