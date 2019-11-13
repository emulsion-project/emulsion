#include <iostream>
#include <fstream>
#include <time.h>
#include <limits>
#include <stdlib.h>
#include "Routines.h"
#include "Multilayer.h"
#include "Multilayer1.h"

using namespace std;

int main()
{
  //calcEnergyLossProbability("Configuration/config.txt",1.5,3.,0.005);
  //calcNearFieldPlane("Configuration/config.txt",200,200,-250e-9,250e-9,-250e-9,250e-9,0,true);
  //calcSurfaceCharges("Configuration/config.txt",25,25);
  //calcNearFieldParticleSurface("Configuration/config.txt",10,10,true);
  calcCrossSections("Configuration/config.txt",350e-9,400e-9,2e-9,false,30e-9);
  //calcCextVSn("Configuration/config.txt",1,15,1,false,1e-9);
  //calcEELSMapping("Configuration/config.txt",100,100,-200e-9,200e-9,-200e-9,200e-9);
  //calcFarFieldSurface("Configuration/config.txt",200,200);

    return 0;
}
