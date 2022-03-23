//
//  deformationGradient.cpp
//  CP_Algorithm2
//
//  Created by fx on 7/15/20.
//  Copyright © 2020 fx. All rights reserved.
//

#include "scp.hpp"
#include<ctime>
#include <vector>
/* public functions */
SCP::SCP(double& theta, double& phi):slip_systems(theta, phi){
    round = -1;
    h = ALATT*sqrt(3.0/2.0);
    /* initialize total strain, or total rate */
    totalStrainRate = 0.001;
    maxiter = 50;
    computeStiffness();
    srand(time(0));
    xi[0][0] = xi[1][1] = xi[2][2] = xi[3][3] = xi[4][4] =xi[5][5] = xi[6][6] = xi[7][7] = xi[8][8] = xi[9][9] = xi[11][11] = Self;
    xi[1][0] = xi[3][2] = xi[5][4] = xi[7][6] = xi[9][8] = xi[11][10] = xi[0][1] = xi[2][3] = xi[4][5] = xi[6][7] = xi[8][9] = xi[10][11] = Coplanar;
    xi[10][10] = xi[2][0] = xi[2][1] = xi[3][0] = xi[3][1] = xi[6][4] = xi[6][5] = xi[7][4] = xi[7][5] = xi[10][8] = xi[10][9] = xi[11][8] = xi[11][9] = xi[0][2] = xi[1][2] = xi[0][3] = xi[1][3] = xi[4][6] = xi[5][6] = xi[4][7] = xi[5][7] = xi[8][10] = xi[9][10] = xi[8][11] = xi[9][11] = Sessile;
    xi[4][3] = xi[5][1] = xi[6][2] = xi[7][0] = xi[8][3] = xi[8][4] = xi[9][0] = xi[9][7] = xi[10][2] = xi[10][6] = xi[11][1] = xi[11][5] = xi[3][4] = xi[1][5] = xi[2][6] = xi[0][7] = xi[3][8] = xi[4][8] = xi[0][9] = xi[7][9] = xi[2][10] = xi[6][10] = xi[1][11] = xi[5][11] = Collinear;
    xi[4][0] = xi[5][2] = xi[6][1] = xi[7][3] = xi[8][1] = xi[8][6] = xi[9][2] = xi[9][5] = xi[10][0] = xi[10][4] = xi[11][3] = xi[11][7] = xi[0][4] = xi[2][5] = xi[1][6] = xi[3][7] = xi[1][8] = xi[6][8] = xi[2][9] = xi[5][9] = xi[0][10] = xi[4][10] = xi[3][11] = xi[7][11] = Glissile;
    xi[4][1] = xi[1][4] = xi[4][2] = xi[2][4] = xi[5][0] = xi[0][5] = xi[5][3] = xi[3][5] = xi[6][3] = xi[3][6] = xi[6][0] = xi[0][6] = xi[7][1] = xi[1][7] = xi[7][2] = xi[2][7] = xi[8][0] = xi[0][8] = xi[8][2] = xi[2][8] = xi[8][5] = xi[5][8] = xi[8][7] = xi[7][8] = xi[9][1] = xi[1][9] = xi[9][3] = xi[3][9] = xi[9][4] = xi[4][9] = xi[9][6] = xi[6][9] = xi[10][1] = xi[1][10] = xi[10][3] = xi[3][10] = xi[10][5] = xi[5][10] = xi[10][7] = xi[7][10] = xi[11][0] = xi[0][11] = xi[11][2] = xi[2][11] = xi[11][4] = xi[4][11] = xi[11][6] = xi[6][11] = Orthorgonal;
    /* initialize stress/strain tensors */
    initializeParameters(theta, phi);
}

int SCP::deformation(){
    fstream fy;
    fstream ff;
    fstream fd;
    Matrix3d d_strain_tot; // total strain matrix
    Matrix3d d_strain_elastic; // elastic strain matrix
    Matrix3d d_strain_plastic; // plastic strain matrix
    d_strain_tot = Matrix3d::Zero();
    d_strain_elastic = Matrix3d::Zero();
    d_strain_plastic = Matrix3d::Zero();
    double deltaTotStrain = totalStrainRate*dt;
    double deltaPlasticStrain = 0.0;
    double deltaElasticStrain = 0.0;
    TOL = 1e-6*dt*totalStrainRate;
    /* select event */
    
    int event = selectReaction();
    if(event == 12){
        d_strain_elastic(1,1) = rates[event]*dt;
        d_strain_tot(1,1) = rates[event]*dt;
    }else{
        
        for(int i = 0; i<12; i++){
            d_strain_elastic(1,1) -= rates[i]*dt;
            d_strain_plastic(1,1) += rates[i]*dt;
            computeOneDRho(i);
            updateDislocationDensity(i);
        }
        
        //d_strain_elastic(1,1) -= rates[event]*dt;
        //d_strain_plastic(1,1) += rates[event]*dt;
        
        computeOneDRho(event);
        //updateDislocationDensity(event);
        
    }
    strain_elastic += d_strain_elastic;
    strain_plastic += d_strain_plastic;
    strain_tot += d_strain_tot;
    Cauchy_stress = computeCauchyStress(strain_elastic);
    
    t += dt;
    acctime +=dt;
    if(Cauchy_stress(1,1) < -1.0){
        ff.open("fracture_toughness.txt", ios::app);
        ff << round << "    " << Th << " "<< Ph << " " << strain_tot(1,1) <<" " << Cauchy_stress(1,1) << endl;
        ff.close();
        return -1;
    }
    if((Cauchy_stress(1,1)-268.0*(strain_tot(1,1)-0.002))<0.001 && yield == 0){
        //if(((Cauchy_stress(1,1)-268.0*(strain_tot(1,1)-0.00001))<0.001 || n == 50000) && yield == 0){
        yield = 1;
        fy.open("yield_stress.txt", ios::app);
        double x0 = (-cos(Th)+sqrt(cos(Th)*cos(Th)+1))*Ph*4.0/PI;
        fy << round << "    " << x0*cos(Th) << "  " << x0*sin(Th) << " " << strain_tot(1,1) <<" " << Cauchy_stress(1,1) << endl;
        fy.close();
        return -1;
    }
    if(acctime >= pointPlot){
        acctime = 0.0;
        plotStressStrainTensors();
        char strainStressFile[30];
        sprintf(strainStressFile, "strain_stress%d.txt", round);
        fss.open(strainStressFile, ios::app);
        //sprintf(strainStressFile, "stress%d.txt", round);
        //fss1.open(strainStressFile, ios::app);
        fss << strain_tot(1,1) << " " << Cauchy_stress(1,1) << endl;
        fss.close();
        //fss1.close();
        plotStressStrainTensors();
        fd.open("dd.txt", ios::app);
        fd << strain_tot(1,1)<<"    "<<strain_plastic(1,1) << " ";
        fratio.open("ratio.txt", ios::app);
        fratio<<strain_plastic(1,1)<<"  ";
        double sumRho = 0.0;
        for(int i = 0; i < 1; i++){
            fd << Rho[i] << "   ";
            sumRho += Rho[i];
            fratio << abs(dGamma[i]/0.001) << "    ";
            //cout<<abs(dGamma[i]/0.001)<<endl;
        }
        fratio<<endl;
        fratio.close();
        fd << sumRho << endl;
        fd.close();
        //cout << "unstable rate = " << (double)unstable_time/(double)total_step*100 << endl;
    }
    /*
    for(int i = 0; i<12; i++){
        computeOneDRho(i);
        updateDislocationDensity(i);
    }
    */
    totalStrainRate = rates[12];
    for(int i = 0; i<12; i++){
        bool whether = computeOneRSS(i);
        computeOneRho_f(i);
        computeOneV(i);
        computeOneDGamma(i);
        plasticStrainRate[i] = abs((abs(dGamma[i])*slip_systems.getSCCN(i))(1,1));
        rates[i] = plasticStrainRate[i];
        totalStrainRate += rates[i];
    }
    double random = (double)rand()/RAND_MAX;
    fstream fts;
    fd.open("dt.txt", ios::app);
    dt = (-1)/totalStrainRate*log(random)*1e-3;
    fd<<dt<<endl;
    fd.close();
    //cout<<"dt = " << dt<<endl;
}

double SCP::getMainStrain(){
    return strain_tot(1,1);
}

/* private functions */
void SCP::computeStiffness(){
    double prefactor = E/(1+Poisson)/(1-2*Poisson);
    double element = (1.0-2.0*Poisson)/2.0;
    stiffness << prefactor*(1-Poisson), prefactor*Poisson, prefactor*Poisson, 0.0, 0.0, 0.0,
    prefactor*Poisson, prefactor*(1-Poisson), prefactor*Poisson, 0.0, 0.0, 0.0,
    prefactor*Poisson, prefactor*Poisson, prefactor*(1-Poisson), 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, prefactor*element, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, prefactor*element, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, prefactor*element;
}

void SCP::initializeParameters(double& theta, double& phi)
{
    round++;
    yield = 0;
    t = 0.0;
    acctime = 0.0;
    double random = (double)rand()/RAND_MAX;
    dt = (-1)/totalStrainRate*log(random)*1e-4;
    strain_tot = Matrix3d::Zero();
    strain_elastic = Matrix3d::Zero();
    strain_plastic = Matrix3d::Zero();
    Cauchy_stress = Matrix3d::Zero();
    Th = theta;
    Ph = phi;
    slip_systems.updateSlipSystems(theta, phi);
    elasticStrainRate = totalStrainRate; // Initialize elastic strain rate
    rates[12] = totalStrainRate;
    /* initialize dislocation densities, plastic strain rate, etc. */
    for(int i = 0; i<12; i++){
        Rho[i] = 1.0e10/12.0; // [cm^-2] initial dislocation density
        double random = (double)rand()/RAND_MAX;
        double m = random-0.5;
        //Rho[i] = 1.0e10/12.0*(1.0+sign(m)*random*T/Tm);
        //Rho[i] = 1.0e10/12.0*(1.0+sign(m)*random*0.001);
        //cout<<"Rho["<<i<<"] = "<<Rho[i]<<endl;
        Rho_f[i] = 0.0;
        V[i] = 0.0;
        dGamma[i] = 0.0;
        dRho[i] = 0.0;
        plasticStrainRate[i] = 0.0;
        rates[i] = plasticStrainRate[i];
    }
    /* initialize one_over_L */
    one_over_L = 0.0; //0.0
    cout<<"Rho_0 = " << Rho[0] << endl << endl;
    //cin.get();
}

int SCP::sign(double& a){
    if(a >= 0.0){
        return 1;
    }else{
        return -1;
    }
}

int SCP::selectReaction(){
    double randomRate = (double)rand()/RAND_MAX * totalStrainRate;
    double compareRate = 0.0;
    for(int i = 0; i < 13; i++){
        compareRate += rates[i];
        if(randomRate <= compareRate){
            return i;
        }
    }
}

Matrix3d SCP::computeCauchyStress(Matrix3d& strain_elastic){
    Matrix3d CauchyStress;
    strain_vec << strain_elastic(0,0), strain_elastic(1, 1), strain_elastic(2,2), 2.0*strain_elastic(1, 2), 2.0*strain_elastic(0, 2),2.0*strain_elastic(0,1);
    stress_vec = stiffness*strain_vec;
    CauchyStress << stress_vec(0), stress_vec(5), stress_vec(4),
    stress_vec(5), stress_vec(1), stress_vec(3),
    stress_vec(4), stress_vec(3), stress_vec(2);
    return CauchyStress;
}

bool SCP::computeOneRSS(int& count){
    
    Vector3d s = slip_systems.getSlipDirection(count);
    Vector3d n = slip_systems.getSlipPlaneNormal(count);
    double RSS = s.transpose()*Cauchy_stress*n;
    double dtau_f = mu*BURGER*sqrt(Rho_int[count]);
    double ddis_d = mu*BURGER*one_over_L;
    tau[count] = RSS - sign(RSS)*(dtau_f+ddis_d);
    if(abs(RSS)>abs(dtau_f)+abs(ddis_d)){
        return true;
    }else{
        return false;
    }
    
}

void SCP::computeOneRho_f(int& count){
    double sum = 0.0;
    double sum_1 = 0.0;
    for(int i = 0; i < 12; i++){
        sum += Rho[i]*abs(slip_systems.getSlipPlaneNormal(count).dot(slip_systems.getSlipDirection(i)));
        sum_1 += xi[count][i]*Rho[i];
    }
    Rho_f[count] = sum;
    Rho_int[count] = sum_1;
}

void SCP::computeOneV(int& count){
    
    Vector3d s = slip_systems.getSlipDirection(count);
    Vector3d n = slip_systems.getSlipPlaneNormal(count);
    double RSS = s.transpose()*Cauchy_stress*n;
    double w = 11.0*BURGER;
    double lambda_alpha = 1.0/(sqrt(Rho_f[count])+one_over_L);
    double dtau_f = mu*BURGER*sqrt(Rho_int[count]+Self*Rho[count]);
    double v_0 = sign(tau[count])*nu0*h/BURGER*(lambda_alpha-w);
    //V[count] = v_0*exp(-Delta_H0/KB/T*(pow((1-pow(abs(RSS/(tau_p+mu*BURGER*sqrt(Rho_int[count]))),p)),q)));
    //V[count] = v_0*exp(-Delta_H0/KB/T*(pow((1-pow(abs(tau[count]/tau_p),p)),q)));
    /*
    if(abs(tau[count])<=tau_p){
        V[count] = v_0*exp(-Delta_H0/KB/T*(pow((1-pow(abs(tau[count]/tau_p),p)),q)));
    }else{
        computeOneB_drag(count);
        V[count] = abs((tau[count]-sign(RSS)*tau_p)*1e9*BURGER/B_drag[count]);
    }
    */
    if(abs(RSS)>dtau_f){
        V[count] = v_0*exp(-Delta_H0/KB/T*(pow((1-pow((abs(RSS)-dtau_f)/tau_p,p)),q)));
    }else{
        V[count] = 0.0;
    }
    //cout<<"V = "<< V[count]<<endl;
    //cout<<"Rho" << Rho[count]<<endl;
}

void SCP::computeOneDGamma(int& count){
    dGamma[count] = Rho[count]*BURGER*V[count];
}

void SCP::computeOneDRho(int& count){
    double lambda = sqrt(Rho_f[count]) + one_over_L;
    dRho[count] = abs(dGamma[count])*dt/BURGER*(lambda-2*BURGER*Rho[count]);
}

void SCP::updateDislocationDensity(int& count){
    Rho[count] = Rho[count] + dRho[count];
}

void SCP::computeOneB_drag(int& count){
    B_drag[count] = 1e-4; // Pa*s
}

void SCP::plotStressStrainTensors(){
    cout << "sigma = " << Cauchy_stress(1,1)<<endl;
    cout << "e_E = " << strain_elastic(1,1) << endl;
    cout << "e_P = " << strain_plastic(1,1) << endl;
    cout << "e_T = " << strain_tot(1,1) << endl << endl;
}
