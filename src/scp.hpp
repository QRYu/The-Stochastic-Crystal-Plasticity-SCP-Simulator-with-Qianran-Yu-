//
//  deformationGradient.hpp
//  CP_Algorithm2
//
//  Created by fx on 7/15/20.
//  Copyright © 2020 fx. All rights reserved.
//

#ifndef scp_hpp
#define scp_hpp

#include "slipSystem.hpp"

typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

class SCP{
private:
    /* deformation parameters */
    SlipSystem slip_systems; /* declare necessary constants in this class -- for tungsten*/
    double h; // [cm] kink height
    double xi[12][12]; //dislocation junction strength matrix
    double totalStrainRate; //[s^{-1}]total strain rate, assume constant strain rate -- [total rate]
    double plasticStrainRate[12];
    double elasticStrainRate;
    double TOL;
    double t; // total time
    double Th; // Theta
    double Ph; // Phi
    double dt;
    double dGamma[12], V[12], B_drag[12];
    double Rho[12], dRho[12], Rho_f[12], tau[12];
    double Rho_int[12];
    double one_over_L;
    double rates[13];
    double acctime;
    int round;
    int maxiter;
    int yield;
    /**
     * Rho : [cm^-2] dislocation density
     * dRho : dislocation density rate (small increment of dislocation density)
     * Lambda : dislocation mean free path
     * Rho_f : forest density
     * tau: [GPa] resolved shear stress 1 GPa = 1e+9 Pa
     * Rho_int: interaction dislocation density
     **/
    Matrix3d strain_tot; // total strain matrix
    Matrix3d strain_elastic; // elastic strain matrix
    Matrix3d strain_plastic; // plastic strain matrix
    Matrix3d Cauchy_stress; // cauchy stress tensor
    Vector6d stress_vec; //stress vector
    Vector6d strain_vec; //strain vector
    Matrix6d stiffness; // stiffness matrix for Hookie's law
    fstream fss,fss1;
    fstream fsr;
    fstream fratio;
    /* private functions */
    void computeStiffness(); // compute stiffness matrix
    int sign(double&); //sign function
    int selectReaction();
    Matrix3d computeCauchyStress(Matrix3d&);
    bool computeOneRSS(int&); // compute one resolved applied stress
    void computeOneRho_f(int&); // compute one forest dislocation density
    void computeOneV(int&); // compute one dislocation velocity
    void computeOneDGamma(int&); // compute one slip rate
    void computeOneDRho(int&); // compute one dislocation density rate
    void updateDislocationDensity(int&);
    void computeOneB_drag(int&); // compute one drag coefficient
public:
    SCP(double&, double&);
    ~SCP(){};
    int deformation(); // process whole deformation
    void initializeParameters(double&, double&);
    void plotStressStrainTensors();
    double getMainStrain();
    
};
#endif /* deformationGradient_hpp */
