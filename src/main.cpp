//
//  main.cpp
//  CP
//
//  Created by fx on 7/7/20.
//  Copyright © 2020 fx. All rights reserved.
//  This code is developed based on Algorithm 1 created by Jaime.
//  Equations of some parameters is not accurate and used the simplified numbers. (i.e. drag coefficient)
//

#include <iostream>
#include "scp.hpp"

int main(int argc, const char * argv[]){
    
    double theta = 0.0;
    double phi = 0.0;
    SCP F(theta, phi);
    double strainT = F.getMainStrain();
    double oneDegree = 1.0/180.0*PI;
    double totTime = 0.0;
    double dtime = 0.0;
    if(Scenario == 64){
        while(theta <= PI/4.0){
            phi = 0.0;
            while(phi <= PI/4.0){
                F.initializeParameters(theta, phi);
                strainT = F.getMainStrain();
                while(strainT <= TotalStrain){
                    F.deformation();
                    strainT = F.getMainStrain();
                }
                phi += oneDegree;
            }
            theta += oneDegree;
        }
    }else{
        dtime = clock();
        while(strainT <= TotalStrain){
            F.deformation();
            strainT = F.getMainStrain();
        }
        cout<<strainT<<endl;
        dtime = clock()-dtime;
        totTime += dtime;
    }
    cout << "totTime = " << totTime/CLOCKS_PER_SEC << endl;
    
    return 0;
}

