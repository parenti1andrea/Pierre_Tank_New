/*----------------------------------------------------------------------------
Syntax : ./TankSim <request_file>
request_file (ASCII file) gives the quantities (numerical values or strings) defining the configuration and the requests
each line contains the name of this quantity (not to be modified) followed by the value chosen by the user
the lines may be given in any order; they are read until the name is END (so the following lines are ignored)

- InputFile: list of incoming particles
    First line (primary): code,energy,theta,phi
    then one line for each ground particles: code,E,x,y,cx,cy,t,weight
- ArrayFile: name of the file describing the array (x,y coordinates of the tanks)
- ArrayGeom: may be GROUND or SHOWER (in the latter case, the coordinates are define in the showerframe, with z along the shower axis) 
- ArrayXavg...ArrayDR: define the position of theshower impact w.r.t. the array, with random variation if SampleNRepeat > 1)
- SampleNRepeat: number of repetitions with the same shower
- SampleRmin...SampleTimeSpread: parameters o do the resampling (unthinning) 
- TankProperties: name of the file describing the properies of the tank
- OutputName: prefix to build the name of output files .summary , .trace
- PrintTrace: number of slots written per trace

General description of the objects
- GroundParticle: in a shower (ground file)  or a single particle depending on the 'SampleMode'option
- ParticleInTank: one a particle is assigned to a given tank, he coordinates are expressed in the tank frame
- Station: a tank with 3 PMTs. The particles are followed until they are absorbed, decay, stop, of hit a PfMT
- PMT: it produces a low-gain and a high-gain trace with a baseline at level 50 
-------------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string.h>

using namespace std;

#include "Global.h"
#include "Shower.h"
#include "Array.h"
#include "Sample.h"
#include "Tank.h"
#include "Utilities.cc"
#include "Trigger.cc"

#define EPSILON 1e-14

#define RESCALE_WEIGHT 0
#define PRINT_WALLHIT 0
#define PRINT_NPE 1
#define PRINT_RISETIME 0

/***********************************************************************************************/
void FractionalTimes(int nslot, int *trace, double base, int nfrac, double *time)
{
  // computes the times where the integrated signal (-baseline) reaches some fractions of the total
  // (nfrac fractions between 0 and 1)
  // time[k] is the time for fraction k/nfrac 
  double sum = 0;
  for(int i=0; i<nslot; i++) sum += (trace[i]-base);
  for(int k=0; k<nfrac; k++) time[k] = 0;
  double sump = 0;
  double dsum = sum/(nfrac+1);
  for(int i=0; i<nslot; i++) {
    sump += trace[i]-base;
    for(int k=0; k<nfrac; k++) {
      double level = (k+1)*dsum;
      if(sump>level&&!time[k]) time[k] = i-(sump-level)/(trace[i]-base);
    }
  }
}
/*--------------------------------------------------------------------------------------*/ 
void FractionalTimes(int nslot, double *trace, int nfrac, double *time)
{
  // same as before, without subtracting a baseline
  double sum = 0;
  for(int i=0; i<nslot; i++) sum += trace[i];
  for(int k=0; k<nfrac; k++) time[k] = 0;
  double sump = 0;
  double dsum = sum/(nfrac+1);
  for(int i=0; i<nslot; i++) {
    sump += trace[i];
    for(int k=0; k<nfrac; k++) {
      double level = (k+1)*dsum;
      if(sump>level&&!time[k]) time[k] = i-(sump-level)/trace[i];
    }
  }
} 
/**************************************************************************************************************/
class GroundParticle
{
  // particle arriving at ground level (either from a shower or injected by hand)
public:
  int Code;
  double E,X,Y,R,Azi,CX,CY,CZ,CXY,Time,Weight;  // CX,CY,CZ are the direction cosines; CXY = sqrt(CX^2+CY^2)
  GroundParticle(int code,double e,double x,double y,double cx,double cy,double time,double weight);
  ~GroundParticle() {};
};
/*------------------------------------------*/
GroundParticle::GroundParticle(int code,double e,double x,double y,double cx,double cy,double time,double weight)
{
  Code = code; E = e; X = x; Y = y;
  double xs = X*ShowerUX+Y*ShowerUY;
  double ys = X*ShowerVX+Y*ShowerVY;
  // R,Azi are the polar coordinates of the particle 
  R = sqrt(xs*xs+ys*ys);
  Azi = atan2(ys,xs);
  CX = cx; CY = cy; CXY = sqrt(cx*cx+cy*cy); CZ = -sqrt(1-(cx*cx+cy*cy));
  Time = time; Weight = weight; 
}
/*--------------------------------------------------------------------------------------------*/
class ParticleInTank
{
public:
  int Entry,Code;
  double E,X,Y,Z,CX,CY,CZ,RelTime,Weight;
  ParticleInTank() {};
  ParticleInTank(int entry,int code,double e,double x,double y,double z,double cx,double cy,double cz,double relt,double w)
  { Entry = entry; Code = code; E = e; X = x; Y = y; Z = z; CX = cx; CY = cy; CZ = cz; RelTime = relt; Weight = w; };
  ~ParticleInTank() {};
};
/*--------------------------------------------------------------------------------------------*/
class PMT
{
public:
  double XCenter, YCenter, HCenter, Radius;
  int Npe_tot,Npe_em,Npe_mu,Npe[NTIMEMAX], LowGainTrace[NSLOTMAX], HighGainTrace[NSLOTMAX];
  int Npe_em_slots[NTIMEMAX], Npe_mu_slots[NTIMEMAX];
  int Npe_tot_FADC, Npe_em_FADC, Npe_mu_FADC;
  double LowGainBase, HighGainBase, LowGain[NSLOTMAX], HighGain[NSLOTMAX];
  double TfracLowGain[9], TfracHighGain[9], RiseTime;
  int SaturHighGain,SaturLowGain;
  PMT() {};
  PMT(double x, double y, double h, double r) { XCenter = x; YCenter = y; HCenter = h; Radius = r; };
  ~PMT() {}; 
};
/*--------------------------------------------------------------------------------------------*/
class Station
{
public:
  int Id,Trigger,Npe_tot;
  double X,Y,R,Azi,Tfront;
  vector <ParticleInTank> Part;
  ParticleInTank *Last_regenerated, *Last_created;
  double smu_rough,sem_rough;
  PMT Pmt[NPMTMAX];
  double LowGainAvg[NSLOTMAX],HighGainAvg[NSLOTMAX],TfracLowGainAvg[9],TfracHighGainAvg[9],RiseTime;
  Station(int id, double x, double y);
  ~Station() {};
  void Inject(GroundParticle gp);
  void FollowParticles();
  void FollowGamma(ParticleInTank tpart);
  void FollowElectron(ParticleInTank tpart);
  void FollowMuon(ParticleInTank tpart);
  void MuonDecay(ParticleInTank tpart, double x, double y, double z, double timr);
  void CherenkovPhotons(double x,double y,double z,double cx,double cy,double cz,double beta,double step,double weight, double timr);
  void FollowPhoton(double x,double y,double z,double cx,double cy,double cz,double t,double prb_refdif,double abs_length);
  int HitPMT(int k, double x, double y, double z, double cx, double cy, double cz, double z_top, double& length);
  void DeltaRays(double x,double y,double z,double cx,double cy,double cz,double beta,double step,double weight,double timr);
  void BuildFADC();
};
/*-------------------------------------------*/
bool sortByNpe(const Station s1, const Station s2) {
  int n1=0, n2=0;
  for(int k=0; k<3; k++) { n1 += s1.Pmt[k].Npe_tot; n2 += s2.Pmt[k].Npe_tot; }
  return n1<n2;
}
bool sortByR(Station &s1, Station &s2) {
  printf("%f %f\n",s1.R,s2.R);
  return s1.R<s2.R;
}
/*------------------------------------------*/
Station::Station(int id, double x, double y)
{
  Id = id;

  double xs,ys;
  if(!strcmp(ArrayGeom,"GROUND")) {
    double xx = x*cos(ArrayRotat)-y*sin(ArrayRotat);
    double yy = x*sin(ArrayRotat)+y*cos(ArrayRotat);
    X = xx+ArrayXCenter;
    Y = yy+ArrayYCenter;
    xs = X*ShowerUX+Y*ShowerUY;
    ys = X*ShowerVX+Y*ShowerVY;
  }
  else {
    xs = x;
    ys = y;
  }
  R = sqrt(xs*xs+ys*ys);
  Azi = atan2(ys,xs);
  Tfront = -xs*ShowerSinTh/ShowerCosTh/CSPEED;
  smu_rough = sem_rough = 0;
  Part.clear();
  for(int k=0; k<NPMT; k++) {
    Pmt[k].Npe_tot = Pmt[k].Npe_em = Pmt[k].Npe_mu = 0;
    Pmt[k].Npe_tot_FADC = Pmt[k].Npe_em_FADC = Pmt[k].Npe_mu_FADC = 0;
    for(int j=0; j<NTIME; j++) Pmt[k].Npe[j] = 0;
  }
  Last_regenerated = NULL;
  Last_created = NULL;
}

/*------------------------------------------*/
void Station::Inject(GroundParticle gp)
{
  // "single" mode -------------------------------------------------------
  if(strcmp(SampleMode,"SHOWER")) {
    ParticleInTank part;
    if(!strcmp(SampleMode,"SNGL_FIX"))
    part = ParticleInTank(0,gp.Code,gp.E,gp.X,gp.Y,TANK_HEIGHT,gp.CX,gp.CY,gp.CZ,gp.Time,1.);
    else {
      double top = TankTopArea*fabs(gp.CZ);
      double side = TankSideArea*gp.CXY;
      if(frand()<top/(top+side)) {
        double ri = TANK_RADIUS*sqrt(frand());
        double azi = 2*M_PI*frand();   
        double x = ri*cos(azi);
        double y = ri*sin(azi);
        part = ParticleInTank(1,gp.Code,gp.E,x,y,TANK_HEIGHT,gp.CX,gp.CY,gp.CZ,gp.Time,1.);
      }
      else {
        double yy = TANK_RADIUS*(2*frand()-1.);
        double xx = sqrt(TANK_RADIUS*TANK_RADIUS-yy*yy);
        double x = (-xx*gp.CX+yy*gp.CY)/gp.CXY;
        double y = (-xx*gp.CY-yy*gp.CX)/gp.CXY;
        double z = TANK_HEIGHT*frand();
        part = ParticleInTank(2,gp.Code,gp.E,x,y,z,gp.CX,gp.CY,gp.CZ,gp.Time,1.);
      }
    }
    Part.push_back(part);
    //printf("inject %d %6.3f\n",gp.Code,gp.E);
    return;
  }

  // "normal" mode (SHOWER) ----------------------------------------------
  // sampling region
  
  if(fabs(gp.R/R-1)>SampleDeltaR) return;
  double dazi = gp.Azi-Azi;
  if(dazi<-M_PI) dazi += 2*M_PI;
  if(dazi>M_PI) dazi -= 2*M_PI;
  if(fabs(dazi)>SampleDeltaAzi) return;
  double zs = gp.X*ShowerWX+gp.Y*ShowerWY;
  // relative time (delay front front plane)
  double reltime = gp.Time+zs/CSPEED;
  if(reltime<0) { fprintf(diag,"relative time %6.2f\n",reltime); return; }
    
  // resampling (unthinning) ------------------------
  double weight = gp.Weight;
  if(abs(gp.Code)<3) weight *= SampleEMWeightFactor;
  if(abs(gp.Code)==3) weight *= SampleMuWeightFactor;
#if RESCALE_WEIGHT
  weight *= pow(gp.R/SampleRthin,3);
#else
  // post-thinning
  if(gp.R<SampleRthin) {
    double val = sqr(gp.R/SampleRthin);
    if(frand()>val) return;
    weight /= val;
  }
#endif
  double sampling_area = 4*R*R*SampleDeltaR*SampleDeltaAzi/ShowerCosTh;
  // resampling weight
  double ntop_avg = weight*TankTopArea/sampling_area;
  int nentry_top = Poisson(ntop_avg);
  double nside_avg = weight*TankSideArea/sqrt(1/(sqr(gp.CX)+sqr(gp.CY))-1)/sampling_area;
  int nentry_side = Poisson(nside_avg);
  if(abs(gp.Code)<3) sem_rough += (nentry_top+nentry_side)*gp.E/.24; 
  if(abs(gp.Code)==3) smu_rough += (nentry_top+nentry_side);

  //if(nentry_top+nentry_side>1) printf("clone %d %6.0f %7.3f %d %f %f\n",nentry_top+nentry_side,R,Azi,gp.Code,gp.E,gp.Weight);
  for(int i=0; i<nentry_top+nentry_side; i++) {
    double x,y,z;
    // top entry
    if(i<nentry_top) {
      //double ri = (TANK_RADIUS-EPS4)*sqrt(frand());
      double ri = TANK_RADIUS*sqrt(frand());
      double azi = 2*M_PI*frand();   
      x = ri*cos(azi);
      y = ri*sin(azi);
      //z = TANK_HEIGHT-EPS4;
      z = TANK_HEIGHT;
    }
    // side entry
    else {
      //double yy = (TANK_RADIUS-EPS4)*(2*frand()-1.);
      double yy = TANK_RADIUS*(2*frand()-1.);
      double xx = sqrt(TANK_RADIUS*TANK_RADIUS-yy*yy);
      x = (-xx*gp.CX+yy*gp.CY)/gp.CXY;
      y = (-xx*gp.CY-yy*gp.CX)/gp.CXY;
      //z = (TANK_HEIGHT-EPS4)*frand()+EPS4;
      z = TANK_HEIGHT*frand();
    }
    // time spread for clones
    double timr = reltime;
    if(nentry_top+nentry_side>0) timr *= exp(grand()*SampleTimeSpread); 
    // delay of entry related to the position on the wall (0 for center at ground level)
    timr += (ShowerWX*x+ShowerWY*y+ShowerWZ*z)/CSPEED;
    ParticleInTank part(1+(i>=nentry_top),gp.Code,gp.E,x,y,z,gp.CX,gp.CY,gp.CZ,timr,1.);
    Part.push_back(part);
  }
  //printf("inject %d  %6.1f %6.3f %6.1f   %7.1f  %f %d  %f %d\n",gp.Code,R,Azi,reltime,gp.Weight,ntop_avg,nentry_top,nside_avg,nentry_side);  

}
/*------------------------------------------*/
void Station::FollowParticles()
{
  for(unsigned int ip=0; ip<Part.size(); ip++) {
#if PRINT_WALLHIT
    printf('Print wall hits')
    for(int k=0; k<3; k++) fprintf(output[k],"-1 -1 -1 %2d %7.4f %7.4f %7.4f\n",Part[ip].Code,Part[ip].E,Part[ip].CX,Part[ip].CY);
#endif
    //printf("%6.0f Follow %2d %8.5f  %6.3f %6.3f %6.3f    %6.3f %6.3f %6.3f  %6.0f %f\n",
    //	   R,Part[ip].Code,Part[ip].E,Part[ip].X,Part[ip].Y,Part[ip].Z,Part[ip].CX,Part[ip].CY,Part[ip].CZ,Part[ip].RelTime,Part[ip].Weight);
    Last_regenerated = &Part[ip];
    //if(fabs(sqr(Part[ip].CX)+sqr(Part[ip].CY)+sqr(Part[ip].CZ)-1)>.0001) printf("!!!! %d  %f %f %f   %f %f %f\n",Part[ip].Code,Part[ip].X,Part[ip].Y,Part[ip].Z,Part[ip].CX,Part[ip].CY,Part[ip].CZ);
    if(Part[ip].Code==1) FollowGamma(Part[ip]);
    if(abs(Part[ip].Code)==2) FollowElectron(Part[ip]);
    if(abs(Part[ip].Code)==3) FollowMuon(Part[ip]);
  }
}
/*------------------------------------------*/
void Station::FollowGamma(ParticleInTank tpart)
{
  double GAM_EMIN = 0.0004;
  double GAM_INT_LENGTH[] =
   {.06,.07,.08,.10,.12,.18,.20,.27,.30,.40,
    .50,.55,.59,.60,.59,.56,.55,.52,.51,.50,
    .48,.47,.46,.46,.46};
  int GAM_NSTEP = sizeof(GAM_INT_LENGTH)/sizeof(double);
  double GAM_PRB_PAIR[] =
   { 0., 0., 0., 0., 0., 0., 0.,.03,.08,.17,
    .33,.42,.59,.70,.79,.88,.91,.93,.95,.97,
    .99, 1., 1., 1., 1.};

  if(tpart.E<GAM_EMIN) return;
  double x = tpart.X; double y = tpart.Y; double z = tpart.Z;
  double cx = tpart.CX; double cy = tpart.CY; double cz = tpart.CZ;
  double t = tpart.RelTime; double e = tpart.E;

  //loop over photon interactions ---------------------------------------
  while(1) {
    int wall;
    double length = DistToWall(x,y,z,cx,cy,cz,0,TANK_HEIGHT,wall);
    // interaction length and pair production probability : from tabulation
    int itab = log10(e)*5+20.5;
    itab = minval(maxval(itab,0),GAM_NSTEP-1);
    double int_length = GAM_INT_LENGTH[itab];
    double prb_pair = GAM_PRB_PAIR[itab];
    double range = -log(frand())*int_length;
    if(range>=length) return;
    x += range*cx;
    y += range*cy;
    z += range*cz;
    t += range/CSPEED;

    // pair production ----------------------------------
    // final direction of electrons is the same as incident photon
    if(frand()<prb_pair) {
      double e1 = (e-2*EL_MASS)*frand();
      ParticleInTank el1(0,2,e1,x,y,z,cx,cy,cz,t,tpart.Weight);
      Last_created = &el1;
      FollowElectron(el1);
      double e2 = (e-2*EL_MASS)-e1;
      ParticleInTank el2(0,-2,e2,x,y,z,cx,cy,cz,t,tpart.Weight);
      Last_created = &el2;
      FollowElectron(el2);
      //printf("      range %6.3f    pair %6.3f %6.3f %6.3f    %6.3f %6.3f %6.3f   %8.5f %8.5f\n",
      //   range,x,y,z,cx,cy,cz,e1,e2);
      return;
    }
    // Compton scattering ------------------------------
    else {
      // Klein-Nishina energy distribution for final photon
      //     (rejection method) 
      double emin = 1/(1/e+2/EL_MASS);
      double smax = emin/e+e/emin+sqr(1+EL_MASS/e-EL_MASS/emin)-1;
      double efin,sfin;
      do {
        efin = emin+frand()*(e-emin);
	sfin = efin/e+e/efin;
      } while(frand()>sfin/smax);
      double e1 = e-efin;
      if(e1>EL_EMIN||efin>GAM_EMIN) {
        // diffusion angles of the photon
        double costhd = 1+EL_MASS*(1/e-1/efin);
        double cxs, cys, czs;
        Scatter(cx,cy,cz,costhd,cxs,cys,czs);
        // scattered electron ---------
        if(e1>EL_EMIN) {
          double p1 = sqrt(2*EL_MASS*e1+e1*e1);
	  double cx_el = (e*cx-efin*cxs)/p1;
          double cy_el = (e*cy-efin*cys)/p1;
	  double cz_el = (e*cz-efin*czs)/p1;
          //if(fabs(sqr(cx_el)+sqr(cy_el)+sqr(cz_el)-1)>.0001) printf("!!!compton %f %f %f\n",cx_el,cy_el,cz_el);
          ParticleInTank el(0,-2,e1,x,y,z,cx_el,cy_el,cz_el,t,tpart.Weight);
          Last_created = &el;
          FollowElectron(el);
        }
        // scattered photon ---------
        if(efin<GAM_EMIN) return;
	e = efin;
	cx = cxs;
	cy = cys;
	cz = czs;
        //printf("    range %6.3f  Compton %6.3f %6.3f %6.3f    %6.3f %6.3f %6.3f   %8.5f %8.5f\n",
	//      range,x,y,z,cxs,cys,czs,efin,e1);
      }
    }
  }
}
/*------------------------------------------*/
void Station::FollowElectron(ParticleInTank tpart)
{
  //printf("   follow %d %f      %d %d %9.7f\n",Part[Last_regenerated].Code,Part[Last_regenerated].E,tpart.Entry,tpart.Code,tpart.E);
  double x = tpart.X; double y = tpart.Y; double z = tpart.Z;
  double cx = tpart.CX; double cy = tpart.CY; double cz = tpart.CZ;
  double range = 0;
  double e = tpart.E;
  double timr = tpart.RelTime;

  int wall;
  //double eloss_brems = 0;
  for(int i=0; i<EL_NSTEP; i++) {
    if(e>EL_ESTEP[i]) {
      // geometrical length to next wall
      double glength = DistToWall(x,y,z,cx,cy,cz,0,TANK_HEIGHT,wall);
      // length to next energy step
      double elength = (e-EL_ESTEP[i])/EL_DEDX[i];
      int exit = 0;
      double step;
      if (glength<elength) { step = glength; exit = 1; }
      else {
        step = elength;
      }
      //printf("   step %2d %7.5f  pos %6.3f %6.3f %6.3f  %d  %8.5f\n",i,step,x,y,z,exit,e);
      double mom = sqrt(sqr(EL_MASS+e)-sqr(EL_MASS));
      double beta = mom/(EL_MASS+e);
      CherenkovPhotons(x,y,z,cx,cy,cz,beta,step,tpart.Weight,timr);
      // displacement
      x += cx*step;
      y += cy*step;
      z += cz*step;
      range += step;
      timr += step/(CSPEED*beta);

      // Multiple scattering
      MultipleScattering(cx,cy,cz,beta,mom,step);

      // Bremsstrahlung
      double prob_brems = BremsTable[i][0]*(step*100);
      double e_brems = 0;
      if(frand()<prob_brems) {
        double costh;
        e_brems = Bremsstrahlung(e,i,step,beta,costh);
        if(e_brems) {
          double cxs,cys,czs;
          Scatter(cx,cy,cz,costh,cxs,cys,czs);
          ParticleInTank ph(0,1,e_brems,x,y,z,cxs,cys,czs,timr,tpart.Weight);
          Last_created = &ph;
          FollowGamma(ph);
        }
        //eloss_brems += e_brems;
        //if(e_brems) printf("brems e %8.6f step %7.5f e_brems %8.6f %6.4f costh %7.5f\n",e,step,e_brems,e_brems/e,costh);
      }
      if(exit) break;
      //if(i<EL_NSTEP-1) e = EL_ESTEP[i];
      e = EL_ESTEP[i];
      e -= e_brems;
    }
  }
  //if(eloss_brems) printf("bremsloss e %f de %f\n",tpart.E,eloss_brems);
}
/*------------------------------------------*/
void Station::FollowMuon(ParticleInTank tpart)
{
  double x = tpart.X; double y = tpart.Y; double z = tpart.Z;
  double cx = tpart.CX; double cy = tpart.CY; double cz = tpart.CZ;
  double range = 0;
  double e_mu = tpart.E;
  double timr = tpart.RelTime;
  //printf("FollowMuon %6.3f   %6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f\n",e_mu,x,y,z,cx,cy,cz);
  int wall;
  
  int npe_before[4],npe_mu[4];
  for(int k=0; k<4; k++) npe_before[k] = Pmt[k].Npe_tot;
  
  double BETGAM0 = sqrt(sqr(1+MU_ESTEP[MU_NSTEP-1]/MU_MASS)-1);
  for(int i=0; i<MU_NSTEP; i++) {
    if(e_mu>MU_ESTEP[i]) {
      // geometrical length to next wall
      double glength = DistToWall(x,y,z,cx,cy,cz,0,TANK_HEIGHT,wall);
      // length to next energy step
      double elength = (e_mu-MU_ESTEP[i])/MU_DEDX[i];
      int exit = 0;
      double step;
      if (glength<elength) { step = glength; exit = 1; }
      else {
        step = elength;
      }
      double beta = sqrt(sqr(MU_MASS+e_mu)-sqr(MU_MASS))/(MU_MASS+e_mu);
      DeltaRays(x,y,z,cx,cy,cz,beta,step,tpart.Weight,timr);

      CherenkovPhotons(x,y,z,cx,cy,cz,beta,step,tpart.Weight,timr);

      // displacement
      x += cx*step;
      y += cy*step;
      z += cz*step;
      range += step;
      timr += step/(CSPEED*beta);
      //if(i<MU_NSTEP-1) e_mu = MU_ESTEP[i+1];
      //else e_mu = 0;
      e_mu = MU_ESTEP[i];
      //printf("e_mu %f  x,y,z %6.3f %6.3f %6.3f  %d\n",e_mu,x,y,z,exit);
      /*
      if(exit) {
        printf("npe ");
        for(int k=0; k<4; k++) printf("%5d   0 ",Pmt[k].Npe_tot-npe_before[k]);
        printf("\n"); 
        return;
      }
      */
      if(exit) return;
    }
  }
  // stopping and decay
  double betgam = sqrt(sqr(1+e_mu/MU_MASS)-1);
  double stop_length = .08*pow(betgam/BETGAM0,3.5);
  x += cx*stop_length;
  y += cy*stop_length;
  z += cz*stop_length;
  if(x*x+y*y<TANK_RADIUS*TANK_RADIUS && z>0 && z<TANK_HEIGHT) {
    for(int k=0; k<4; k++) npe_mu[k] = Pmt[k].Npe_tot-npe_before[k];
    MuonDecay(tpart,x,y,z,timr);
    //printf("decay %6.3f %6.3f %6.3f  %6.3f %6.3f  ",x,y,z,tpart.E,-tpart.CZ);
    //for(int k=0; k<4; k++) printf(" %4d %4d  ",npe_mu[k],Pmt[k].Npe_tot-npe_before[k]-npe_mu[k]);
    //printf("\n");
  }
}
/*---------------------------------------------------------------------*/
void Station::MuonDecay(ParticleInTank tpart, double x, double y, double z, double timr)
{
  //Michel distribution with rho=3/4, eta=0
  int try_rndm  = 1;
  double a,b;
  while(try_rndm) {
    a = frand();
    b = a*a*(3-2*a);
    if(b>frand()) try_rndm = 0;
  }
  double e_el = a*MU_MASS/2;
  double cz_el = 2*frand()-1;
  double phi = M_PI*(2*frand()-1);
  double cx_el = sqrt(1-cz_el*cz_el)*cos(phi);
  double cy_el = sqrt(1-cz_el*cz_el)*sin(phi);
  int ndec = 1;
  if(tpart.Weight>1) ndec = Poisson(tpart.Weight);
  for(int i=0; i<ndec; i++) {
    double tdec = timr;
    if(tpart.Code>0) {
      tdec -= MUPOS_TDECAY*log(frand());
      ParticleInTank el(0,2,e_el,x,y,z,cx_el,cy_el,cz_el,tdec,1.);
      Last_created = &el;
      FollowElectron(el);
    }
    else {
     tdec -= MUNEG_TDECAY*log(frand());
      ParticleInTank el(0,-2,e_el,x,y,z,cx_el,cy_el,cz_el,tdec,1.);
      Last_created = &el;
      FollowElectron(el);
    }
  }
}
/**************************************************************************************************/
void Station::CherenkovPhotons(double x,double y,double z,double cx,double cy,double cz,double beta,double step,double weight,double timr)
{
  double costh_cher = 1/(beta*WATER_INDEX);
  if(costh_cher>1) return;
  int n = Poisson(CHERENKOV_FACTOR*(1-sqr(costh_cher))*step*weight);
  //printf("      step %6.4f factor %6.3f %6.3f   %f %d\n",step,costh_cher,CHERENKOV_FACTOR,CHERENKOV_FACTOR*(1-sqr(costh_cher))*step,n);

  for(int i=0; i<n; i++) {
    // direction and point of emission
    double cx_ph,cy_ph,cz_ph;
    Scatter(cx,cy,cz,costh_cher,cx_ph,cy_ph,cz_ph);
    //if(fabs(cx_ph*cx_ph+cy_ph*cy_ph+cz_ph*cz_ph-1)>.00001) printf("scatter %f %f \n",cx*cx+cy*cy+cz*cz,cx_ph*cx_ph+cy_ph*cy_ph+cz_ph*cz_ph);
    double r_ph = step*frand();
    double x_ph = x+r_ph*cx;
    double y_ph = y+r_ph*cy;
    double z_ph = z+r_ph*cz;
    double t_ph = timr+r_ph/CSPEED;
    //if(abs(Part[Part.size()-1].Code)==3) printf("    cherphot %6.1f\n",t_ph);
    // wavelength and absorption properties at this wavelength
    double prb = frand();
    double prb_reflect = 1, abs_length = 1000;
    for(int j=0; j<NWAVELENGTH; j++) {
      if(prb<CHERENKOV_PROBABILITY[j]) {
	prb_reflect = LINER_REFLECTIVITY[j];
        abs_length = WATER_ABSORPTION_LENGTH[j];
        break;
      }
    }
    //printf("Follow %6.3f %6.3f %6.3f   %6.4f %6.4f %6.4f  %f %f\n",x_ph,y_ph,z_ph,cx_ph,cy_ph,cz_ph,prb_reflect,abs_length);
    FollowPhoton(x_ph,y_ph,z_ph,cx_ph,cy_ph,cz_ph,t_ph,prb_reflect,abs_length);
    //printf("    fini\n");
  }
}
/**************************************************************************************************/
void Station::FollowPhoton(double x, double y, double z, double cx, double cy, double cz, double t, double prb_reflect, double abs_length)
{
  // if LSD option; else 'upper' is always TRUE
  int upper = (z>TANK_ZSPLIT);
  double z_bot, z_top;
  if(upper) { z_bot = TANK_ZSPLIT; z_top = TANK_HEIGHT; }
  else { z_bot = 0; z_top = TANK_ZSPLIT; }
  double length;

  //if(fabs(cx*cx+cy*cy+cz*cz-1)>.00001) printf("entry\n");

  int nhitwall[4]={0,0,0,0};
  double t_init = t;
  while(nhitwall[1]+nhitwall[2]+nhitwall[3]<1000) {
    // PMT entry ?
    for(int k=0; k<NPMT; k++) {
      if((upper&&k<NPMT_UP)||(!upper&&k>=NPMT_UP)) {
        //printf("upper %d\n",upper);
        if(HitPMT(k,x,y,z,cx,cy,cz,z_top,length)) {
          // absorption before hitting the PMT ?
          if(length>-abs_length*log(frand())) return;
          t += length*WATER_INDEX/CSPEED;
          //printf("Path of the hit %f \n",length);
          //printf("Time of the hit %f \n",t);
          int it = int(t/TIMEUNIT);
          Pmt[k].Npe_tot++;
          if(abs(Last_regenerated->Code)<3) Pmt[k].Npe_em++;
          if(abs(Last_regenerated->Code)==3) Pmt[k].Npe_mu++;
          if(it>=0&&it<NTIME){
            Pmt[k].Npe[it]++;
            if(abs(Last_regenerated->Code)<3) Pmt[k].Npe_em_slots[it]++;
            if(abs(Last_regenerated->Code)==3) Pmt[k].Npe_mu_slots[it]++;
            //printf("Npe %d",Pmt[k].Npe[it]);
          }
#if PRINT_WALLHIT
          printf('Print wall hits')
          fprintf(output[k],"%6.1f %2d %2d %2d\n",t-t_init,nhitwall[1],nhitwall[2],nhitwall[3]);
#endif
          return;
        }
      }
    }
    // impact on wall
    if(t-t_init>3000) { 
      fprintf(diag,"tlimit %7.2f x,y,z %6.3f %6.3f %6.3f  cx,cy,cz %6.3f %6.3f %6.3f\n",t-t_init,x,y,z,cx,cy,cz);
      return;
    }
    int wall;
    length = DistToWall(x,y,z,cx,cy,cz,z_bot,z_top,wall);
    // protection against rounding overflow
    if(length<EPSILON||length>sqrt(4*TANK_RADIUS*TANK_RADIUS+sqr(z_top-z_bot))-EPSILON) {
      fprintf(diag,"FollowPhoton  start x,y,z %6.3f %6.3f %6.3f  cx,cy,cz %6.3f %6.3f %6.3f   wall %d  length %18.16f\n",x,y,z,cx,cy,cz,wall,length);
      return;
    }
    if(length>-abs_length*log(frand())) return;
    if(frand()>prb_reflect*REFLECTIVITY_FACTOR[wall-1]) return;

    nhitwall[wall]++;
    if(!wall) return;

    // protection against position outside the water
    length -= 5*EPSILON;
    x += length*cx; y += length*cy; z += length*cz;
    if(x*x+y*y>sqr(TANK_RADIUS) || z>z_top || z<z_bot) {
      fprintf(diag,"FollowPhoton length %6.4f arrival x,y,dr2,ang,z,cz %19.16f %19.16f %18.16f %19.16f  %19.16f %19.16f\n",
	      length,x,y,x*x+y*y-sqr(TANK_RADIUS),(x*cy-y*cx)/sqrt(x*x+y*y),z,cz);
      return;
    }
    t += length*WATER_INDEX/CSPEED;
    if(frand()<SPECULAR_FRACTION) {
      if(wall==2) {
        double fact = (x*cx+y*cy)/(x*x+y*y);
        cx -= 2*fact*x;
        cy -= 2*fact*y;
      }
      else cz *= -1;
      //if(fabs(cx*cx+cy*cy+cz*cz-1)>.00001) printf("spec %d\n",wall);
    }
    else {
      double c2bet = 2*frand()-1;
      double cbet = sqrt((1+c2bet)/2);
      double sbet = sqrt((1-c2bet)/2);
      double alph = 2*frand()*M_PI;
      double calph = cos(alph);
      double salph = sin(alph);
      if(wall==2) {
        cx = (-cbet*x+sbet*salph*y)/sqrt(x*x+y*y);
        cy = (-cbet*y-sbet*salph*x)/sqrt(x*x+y*y);
        cz = sbet*calph;
      }
      else {
        cx = sbet*calph;
        cy = sbet*salph;
        cz = cbet;
        if(wall==1) cz = -cz;
      }
      //if(fabs(cx*cx+cy*cy+cz*cz-1)>.00001) printf("diff %d\n",wall);
    }
  }
}
/**************************************************************************************************/
int Station::HitPMT(int k, double x, double y, double z, double cx, double cy, double cz, double z_top, double& length)
{
  double dx = x-X_PMT[k];
  double dy = y-Y_PMT[k];
  double dz = z-Z_PMT[k];
  double b = cx*dx+cy*dy+cz*dz;
  double c = sqr(RAD_PMT[k])-dx*dx-dy*dy-dz*dz;
  if(b>0||b*b+c<0) return 0;
  length = -b-sqrt(b*b+c);
  double zhit = z+length*cz;
  //if(zhit<1.1) printf("!!! %d  %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f -- %6.3f --  %f %f\n",k,dx,dy,dz,cx,cy,cz,cx*cx+cy*cy+cz*cz,b,c);
  return zhit<z_top;
}
/**************************************************************************************************/
#define DELTARAY_CONST .008445
void Station::DeltaRays(double x,double y,double z,double cx,double cy,double cz,double beta,double step,double weight, double timr)
{
  double ec_min = EL_EMIN;
  double ec_max = EL_MASS*beta*beta/(1-beta*beta);
  int nrays = Poisson(DELTARAY_CONST*step/beta*(1/ec_min-1/ec_max));
  double eloss = 0;
  for(int i=0; i<nrays; i++) {
    double ec = 1/(1/ec_min-frand()*(1/ec_min-1/ec_max));
    eloss += ec;
    //if(ec<ec_min||ec>ec_max) printf("!!! ec %8.5f\n",ec);
    double costhd = ec/ec_max*sqrt((sqr(ec_max)+2*EL_MASS*ec_max)/(sqr(ec)+2*EL_MASS*ec));
    double dist = step*frand();
    double xd = x+cx*dist;
    double yd = y+cy*dist;
    double zd = z+cz*dist;
    double timd = timr+dist/(beta*CSPEED);

    double cxs, cys, czs;
    Scatter(cx,cy,cz,costhd,cxs,cys,czs);
    ParticleInTank el(0,-2,ec,xd,yd,zd,cxs,cys,czs,timd,weight);
    FollowElectron(el);
  }
}
/*----------------------------------------------------------------*/
void Station::BuildFADC()
{
  for(int k=0; k<NPMT; k++) {
    PMT *pmt = &Pmt[k];
    for(int j=0; j<NSLOT; j++) {
      pmt->HighGain[j] = pmt->LowGain[j] = 0;
      int npe_slot = 0;
      int npe_slot_em = 0;
      int npe_slot_mu = 0;

      for(int it=j*TIMESLOT; it<(j+1)*TIMESLOT; it++){ 
        npe_slot += pmt->Npe[it];
        npe_slot_em += pmt->Npe_em_slots[it];
        npe_slot_mu += pmt->Npe_mu_slots[it];
      }
      if(npe_slot) {
        //printf("Time slot is %d",j);
        //printf("Number of PE %d",npe_slot);
        double sigma = sqrt(SIGMA_PER_PE*npe_slot);
        pmt->HighGain[j] = FADC_PER_PE*(npe_slot+sigma*grand());
        pmt->LowGain[j] = pmt->HighGain[j]/DYNODE_ANODE_RATIO;
      }
    float npe_tot_fadc = 0;
    for(int j=0; j<NSLOT; j++){npe_tot_fadc += pmt->HighGain[j];};
    pmt->Npe_tot_FADC = npe_tot_fadc;
    }
    // high gain traces ------------------------------------
    pmt->HighGainBase = 50.+frand();
    for(int j=0; j<NSLOT; j++) {
      int term = pmt->HighGainBase+SIGMA_BASE*grand()+pmt->HighGain[j]+.5;
      if(term<MAXFADC) pmt->HighGainTrace[j] = term;
      else {
        pmt->HighGainTrace[j] = MAXFADC-1;
        (pmt->SaturHighGain)++;
      }
    }
    FractionalTimes(NSLOT,pmt->HighGainTrace,pmt->HighGainBase,9,pmt->TfracHighGain);
    // low gain traces --------------------------------------
    pmt->LowGainBase = 50.+frand();
    for(int j=0; j<NSLOT; j++) {
      int term = pmt->LowGainBase+SIGMA_BASE*grand()+pmt->LowGain[j]+.5;
      if(term<MAXFADC) pmt->LowGainTrace[j] = term;
      else {
        pmt->LowGainTrace[j] = MAXFADC-1;
        (pmt->SaturLowGain)++;
      }
    }
    FractionalTimes(NSLOT,pmt->LowGainTrace,pmt->LowGainBase,9,pmt->TfracLowGain);
  }

  // average on the PMTs (baseline subtracted)
  for(int j=0; j<NSLOT; j++) {
    LowGainAvg[j] = HighGainAvg[j] = 0;
    for(int k=0; k<NPMT; k++) {
      LowGainAvg[j] += Pmt[k].LowGainTrace[j]-Pmt[k].LowGainBase;
      HighGainAvg[j] += Pmt[k].HighGainTrace[j]-Pmt[k].HighGainBase;
    }
    LowGainAvg[j] /= NPMT;
    HighGainAvg[j] /= NPMT;
  }
  FractionalTimes(NSLOT,LowGainAvg,9,TfracLowGainAvg);
  FractionalTimes(NSLOT,HighGainAvg,9,TfracHighGainAvg);
}
/**************************************************************************************************/
vector <GroundParticle> GrdPart;
vector <Station> ArrayStation;
void ReadRequest(char *request_file);
void ReadArray();
void ReadGroundParticles();
void ReadTankProperties(char* name);

/*********************************************************************************************/
int main(int argc, char *argv[])
{
  diag = fopen("diagnostics","w");
  ReadRequest(argv[1]);
  frand();
  // Shower characteristics and shower frame
  if(!strcmp(SampleMode,"SHOWER")) {
    fscanf(input,"%d %lf %lf %lf",&ShowerPrimary,&ShowerEnergy,&ShowerTheta,&ShowerPhi);
    printf("Shower  Primary %d  Energy %6.3f  Theta %6.2f  Phi %6.2f\n",ShowerPrimary,ShowerEnergy/1e18,ShowerTheta*180/M_PI,ShowerPhi*180/M_PI);
    ShowerCosTh = cos(ShowerTheta); ShowerSinTh = sin(ShowerTheta);
    ShowerWX = ShowerSinTh*cos(ShowerPhi); ShowerWY = ShowerSinTh*sin(ShowerPhi); ShowerWZ = ShowerCosTh;
    BuildFrame(ShowerWX,ShowerWY,ShowerWZ,ShowerUX,ShowerUY,ShowerUZ,ShowerVX,ShowerVY,ShowerVZ);
  }
  BuildTables();

  ReadGroundParticles();
  
  for(int irepeat=0; irepeat<SampleNRepeat; irepeat++) {
    // "single particle" mode -----------------------------------------------------
    printf("Sample mode is ");
    printf(SampleMode,"\n");
    if(strcmp(SampleMode,"SHOWER")) {
      for(unsigned int igp=0; igp<GrdPart.size(); igp++) {
        if(igp%1000==999) printf("%d\n",igp+1);
        //printf("\nCode of particle %f \n",GrdPart[igp].Code);
        //printf("Energy of particle %f \n",GrdPart[igp].E);
        //printf("X of particle %f \n",GrdPart[igp].X);
        Station stat(0,0,0);
        stat.Inject(GrdPart[igp]);
        stat.FollowParticles();
        stat.BuildFADC();

        if(igp > 1e4){break;};
        if(stat.Pmt[0].Npe_tot+stat.Pmt[1].Npe_tot+stat.Pmt[2].Npe_tot>30) {
          /* charge, et peak de sliding sum
          int maxsum[3] = {0,0,0}, charge[3] = {0,0,0};
          for(int itr=0; itr<NSLOT; itr++) {
            for(int k=0; k<3; k++) {
              charge[k] += (stat.Pmt[k].HighGainTrace[itr]-50);
              int sum = 0;
              for(int l=0; l<3; l++) sum += (stat.Pmt[k].HighGainTrace[itr+l]-50);
              if(sum>maxsum[k]) maxsum[k] = sum;
            }
	  }
          printf("%4d %4d %4d   %4d %4d %4d\n",charge[0],charge[1],charge[2],maxsum[0],maxsum[1],maxsum[2]);
          */
#if PRINT_NPE
          fprintf(output[3],"%d\t%2d\t%7.4f",igp,GrdPart[igp].Code,GrdPart[igp].E);
          for(int k=0; k<NPMT; k++) fprintf(output[3],"\t%4d",stat.Pmt[k].Npe_tot);
          for(int k=0; k<NPMT; k++) fprintf(output[3],"\t%4d", stat.Pmt[k].Npe_tot_FADC);
          fprintf(output[3],"\n");
#endif
          if(PrintTrace) {
            for(int itr=0; itr<PrintTrace; itr++) {
              for(int k=0; k<NPMT; k++) fprintf(trace," %4d",stat.Pmt[k].HighGainTrace[itr]);
              for(int k=0; k<NPMT; k++) fprintf(trace," %4d",stat.Pmt[k].LowGainTrace[itr]);
              fprintf(trace,"\n");
	    }
          }
        }
      }
    }

    // "shower" mode --------------------------------------------------------------
    else {
      ArrayStation.clear();
      ReadArray();
      printf("\n----------------------------\nsim %d  %f %f %f\n----------------------------\n",
       irepeat+1,ArrayXCenter,ArrayYCenter,ArrayRotat);
      // sort the stations by increasing R if "GROUND" geometry
      //if(strcmp(ArrayGeom,"SHOWER")) std::sort(ArrayStation.begin(),ArrayStation.end(),sortByR);
      //printf("sorted\n");
      for(unsigned int igp=0; igp<GrdPart.size(); igp++)
        for(unsigned int is=0; is<ArrayStation.size(); is++) ArrayStation[is].Inject(GrdPart[igp]);
      for(unsigned int is=0; is<ArrayStation.size(); is++) {
        Station *stat = &ArrayStation[is];
        stat->FollowParticles();
        int nmu = 0;
        double Eem = 0;
        for(unsigned int ip=0; ip<stat->Part.size(); ip++) {
          nmu += (abs(stat->Part[ip].Code)==3);
          Eem += (abs(stat->Part[ip].Code)<3)*stat->Part[ip].E;
        }
        printf("%4d %6.0f %6.3f %6.0f  %4lu %3d %6.3f  ",
         stat->Id,stat->R,stat->Azi,stat->Tfront,stat->Part.size(),nmu,Eem);
        for(int k=0; k<NPMT; k++) printf("%5d (%5d %5d) ",stat->Pmt[k].Npe_tot,stat->Pmt[k].Npe_em,stat->Pmt[k].Npe_mu);
        printf("  %7.2f %7.2f",stat->sem_rough,stat->smu_rough);
#if PRINT_NPE
        fprintf(output[3],"%6.1f %7.3f  ",stat->R,stat->Azi);
        for(int k=0; k<NPMT; k++) fprintf(output[3],"%5d %5d ",stat->Pmt[k].Npe_em,stat->Pmt[k].Npe_mu);
        fprintf(output[3],"\n");
#endif
	// build FADC traces -------------------------------------------------------------
        stat->BuildFADC();
        int trace3[3*NSLOT];
        double base[3], vempeak[3];
        for(int k=0; k<3; k++) {
          base[k] = stat->Pmt[k].HighGainBase;
          vempeak[k] = 50;
          for(int j=0; j<NSLOT; j++)  trace3[3*j+k] = stat->Pmt[k].HighGainTrace[j];
        }
        // apply trigger algorithms ---------------------------------------------------
        int thresh1 = IsThreshold(trace3,base,vempeak,1.7);
        int thresh2 = IsThreshold(trace3,base,vempeak,3.2);
        int totd = IsToTd(trace3,base,vempeak,.2,10,120);
        int mops = IsMoPS(trace3,4,30,120,5,1);
        stat->Trigger = (thresh1>0)+2*(thresh2>0)+4*(totd>0)+8*(mops>0); 
        printf(" trig %4d %4d %4d %4d",thresh2,thresh1,totd,mops);
        printf("\n");

        if(PrintTrace) {
          double signal = 0;
          for(int k=0; k<NPMT; k++) signal += (stat->Pmt[k].Npe_mu+stat->Pmt[k].Npe_em)/240.;
          fprintf(trace,"%6.1f %7.3f %6.2f\n",stat->R,stat->Azi,signal);
          for(int itr=0; itr<PrintTrace; itr++) {
            for(int k=0; k<NPMT; k++) fprintf(trace," %4d",stat->Pmt[k].HighGainTrace[itr]);
            for(int k=0; k<NPMT; k++) fprintf(trace," %4d",stat->Pmt[k].LowGainTrace[itr]);
            fprintf(trace,"\n");
	  }
        }

#if PRINT_RISETIME
        fprintf(output[0],"%7.3f %6.3f %7.3f %5.0f %6.3f %2d ",ShowerEnergy/1e18,ShowerTheta,ShowerPhi,stat->R,stat->Azi,stat->Trigger);
        // risetime -----------------
        for(int k=0; k<3; k++) {
          double risetime = stat->Pmt[k].TfracHighGain[4]-stat->Pmt[k].TfracHighGain[0];
          if(stat->Pmt[k].SaturHighGain)  risetime = -(stat->Pmt[k].TfracLowGain[4]-stat->Pmt[k].TfracLowGain[0]);
          if(stat->Pmt[k].SaturLowGain) risetime = 0;
          stat->Pmt[k].RiseTime = risetime;
          fprintf(output[0],"%6d %6.2f ",stat->Pmt[k].Npe_tot,risetime);
        }
        double risetime = stat->TfracHighGainAvg[4]-stat->TfracHighGainAvg[0];
        if(stat->Pmt[0].SaturHighGain+stat->Pmt[1].SaturHighGain+stat->Pmt[2].SaturHighGain)
          risetime = -(stat->TfracLowGainAvg[4]-stat->TfracLowGainAvg[0]);
        if(stat->Pmt[0].SaturLowGain+stat->Pmt[1].SaturLowGain+stat->Pmt[2].SaturLowGain) risetime = 0;
        stat->RiseTime = risetime;
        fprintf(output[0],"  %6.2f\n",risetime);
#endif
      }
      // one summary line per event --------------------------
      if(PrintSummary) {
        int ntrig = 0;
        for(unsigned int is=0; is<ArrayStation.size(); is++) ntrig += (ArrayStation[is].Trigger>0);
        fprintf(sum,"%7.3f %6.3f %7.3f %d",ShowerEnergy/1e18,ShowerTheta,ShowerPhi,ntrig);
        for(unsigned int is=0; is<ArrayStation.size(); is++) {
          Station *stat = &ArrayStation[is];
          fprintf(sum,"  %6.1f %6.3f %6d %6d %2d %6.2f",stat->R,stat->Azi,
	   stat->Pmt[0].Npe_em+stat->Pmt[1].Npe_em+stat->Pmt[2].Npe_em,stat->Pmt[0].Npe_mu+stat->Pmt[1].Npe_mu+stat->Pmt[2].Npe_mu,
           stat->Trigger,stat->RiseTime);
        }
        fprintf(sum,"\n");
      }
    }
  }
  return 1;
}
/*--------------------------------------------------------------------------------------------*/
void ReadRequest(char* request_file)
{
  char line[1000], item[100], value[100], outname[100], name[100];
  request = fopen(request_file,"r");

  while(fgets(line,1000,request)) {
    printf("%s",line);
    if(sscanf(line,"%s %s",item,value)<2 || !strcmp(item,"END")) break;
    if(!strcmp(item,"InputFile")) { input = fopen(value,"r"); continue; }
    if(!strcmp(item,"ArrayFile")) { array_file = fopen(value,"r"); continue; }
    if(!strcmp(item,"TankProperties")) { ReadTankProperties(value); continue; }
    if(!strcmp(item,"OutputName")) { sprintf(outname,"%s",value); continue; }  
    if(!strcmp(item,"PrintTrace")) {
      PrintTrace = atoi(value);
      sprintf(name,"%s.trace",outname);
      trace = fopen(name,"w");
      printf("trace file %s\n",name);
      continue;
    }
    if(!strcmp(item,"PrintSummary")) {
      PrintSummary = atoi(value);
      sprintf(name,"%s.summary",outname);
      sum = fopen(name,"w");
      continue;
    }
    if(!strcmp(item,"Noutput")) {
      Noutput = atoi(value);
      if(Noutput>10) Noutput=10;
      for(int k=0; k<Noutput; k++) {
        sprintf(name,"%s_%01.2f_%d",outname,FADC_PER_PE,k);
        output[k] = fopen(name,"w");
        printf("output file %s\n",name);
      }
      continue;
    }
    if(!strcmp(item,"SampleMode")) { sprintf(SampleMode,"%s",value); continue; } 
    if(!strcmp(item,"SampleNRepeat")) { SampleNRepeat = atoi(value); continue; } 
    if(!strcmp(item,"SampleRmin")) { SampleRmin = atof(value); continue; }
    if(!strcmp(item,"SampleRmax")) { SampleRmax = atof(value); continue; }
    if(!strcmp(item,"SampleRthin")) { SampleRthin = atof(value); continue; }
    if(!strcmp(item,"SampleDeltaR")) { SampleDeltaR = atof(value); continue; }
    if(!strcmp(item,"SampleDeltaAzi")) { SampleDeltaAzi = atof(value); continue; }
    if(!strcmp(item,"SampleEMWeightFactor")) { SampleEMWeightFactor = atof(value); continue; }
    if(!strcmp(item,"SampleMuWeightFactor")) { SampleMuWeightFactor = atof(value); continue; }
    if(!strcmp(item,"SampleTimeSpread")) { SampleTimeSpread = atof(value); continue; }
    if(!strcmp(item,"ArrayGeom")) { sprintf(ArrayGeom,"%s",value); continue; }
    if(!strcmp(item,"ArrayXavg")) { ArrayXavg = atof(value); continue; }
    if(!strcmp(item,"ArrayDX")) { ArrayDX = atof(value); continue; }
    if(!strcmp(item,"ArrayYavg")) { ArrayYavg = atof(value); continue; }
    if(!strcmp(item,"ArrayDY")) { ArrayDY = atof(value); continue; }
    if(!strcmp(item,"ArrayRavg")) { ArrayRavg = atof(value); continue; }
    if(!strcmp(item,"ArrayDR")) { ArrayDR = atof(value); continue; }
    printf("ReadRequest: unknown name %s\n",item);
  }
}
/*--------------------------------------------------------------------------------------------*/
void ReadGroundParticles()
{
  GrdPart.clear();
  int code;
  double x,y,e,cx,cy,t,w;
  int n = 0;
  while((fscanf(input,"%d %lf %lf %lf %lf %lf %lf %lf",&code,&e,&x,&y,&cx,&cy,&t,&w)==8)) {
    // conversion to Aires particle code
    if(code==3) code = -2;
    if(code==5) code = 3;
    if(code==6) code = -3;
    // keep only photons, e+-, mu+- between Rmin and Rmax or anywhere is not SHOWER sampling mode 
    if(abs(code)<4) {
      GroundParticle grdp(code,e,x,y,cx,cy,t,w);
      n++;
      if(strcmp(SampleMode,"SHOWER") || (grdp.R>SampleRmin*(1-SampleDeltaR)&&grdp.R<SampleRmax*(1+SampleDeltaR))) GrdPart.push_back(grdp);
    }
  }
  printf("%d particles - %lu within (Rmin,Rmax)\n",n,GrdPart.size());
}
/*--------------------------------------------------------------------------------------------*/
void ReadArray()
{
  if(!strcmp(ArrayGeom,"GROUND")) {
    ArrayXCenter = ArrayXavg + (2*frand()-1)*ArrayDX;
    ArrayYCenter = ArrayYavg + (2*frand()-1)*ArrayDY;
    ArrayRotat = ArrayRavg + (2*frand()-1)*ArrayDR;
  }
  int id;
  double x,y;
  rewind(array_file);
  while((fscanf(array_file,"%d %lf %lf",&id,&x,&y)==3)) {
    Station stat(id,x,y);
    if(stat.R>SampleRmin&&stat.R<SampleRmax) ArrayStation.push_back(stat);
  }
}
/*--------------------------------------------------------------------------------------------*/
void ReadTankProperties(char* name)
{
  FILE* tank;
  tank = fopen(name,"r");
  printf("Tank Properties -------------------------------\n");
  fscanf(tank,"%*s %lf",&TANK_RADIUS);
  fscanf(tank,"%*s %lf",&TANK_HEIGHT);
  fscanf(tank,"%*s %lf",&TANK_ZSPLIT);
  TankTopArea = M_PI*sqr(TANK_RADIUS);
  TankSideArea = 2*TANK_RADIUS*TANK_HEIGHT;
  printf("RADIUS %6.3f  HEIGHT %6.3f  ZSPLIT %6.3f\n",TANK_RADIUS,TANK_HEIGHT,TANK_ZSPLIT);
  fscanf(tank,"%*s %d",&NPMT);
  fscanf(tank,"%*s %d",&NPMT_UP);
  fscanf(tank,"%*s");
  for(int k=0; k<NPMT; k++) {
    fscanf(tank,"%lf %lf %lf %lf",X_PMT+k,Y_PMT+k,Z_PMT+k,RAD_PMT+k);
    printf("PMT %d   %7.4f %7.4f %7.4f %6.4f\n",k+1,X_PMT[k],Y_PMT[k],Z_PMT[k],RAD_PMT[k]);
  }
  fscanf(tank,"%*s %lf",&WATER_ABSORPTION_LENGTH_MAX);
  fscanf(tank,"%*s %lf",&ADD_ABSORPTION_LENGTH);
  fscanf(tank,"%*s %lf",REFLECTIVITY_FACTOR);
  fscanf(tank,"%*s %lf",REFLECTIVITY_FACTOR+1);
  fscanf(tank,"%*s %lf",REFLECTIVITY_FACTOR+2);
  fscanf(tank,"%*s %lf",&COLLECTION_EFFICIENCY);
  fscanf(tank,"%*s %lf",&SPECULAR_FRACTION);
  fscanf(tank,"%*s");
  for(int j=0; j<NWAVELENGTH; j++) {
    fscanf(tank,"%lf",WATER_ABSORPTION_LENGTH+j);
    WATER_ABSORPTION_LENGTH[j] = 1/(1/(WATER_ABSORPTION_LENGTH_MAX*WATER_ABSORPTION_LENGTH[j])+1/ADD_ABSORPTION_LENGTH);
  }
  fscanf(tank,"%*s");
  for(int j=0; j<NWAVELENGTH; j++) fscanf(tank,"%lf",LINER_REFLECTIVITY+j);
  fscanf(tank,"%*s");
  for(int j=0; j<NWAVELENGTH; j++) fscanf(tank,"%lf",PMT_QUANTUM_EFFICIENCY+j);

  for(int j=0; j<NWAVELENGTH; j++) {
    double term = PMT_QUANTUM_EFFICIENCY[j]/sqr(WAVELENGTH_MIN+j*DWAVELENGTH);
    if(j==0) CHERENKOV_PROBABILITY[j] = term;
    else CHERENKOV_PROBABILITY[j] = CHERENKOV_PROBABILITY[j-1]+term;
  }
  CHERENKOV_FACTOR = CHERENKOV_CONSTANT*CHERENKOV_PROBABILITY[NWAVELENGTH-1]*DWAVELENGTH*COLLECTION_EFFICIENCY;
  for(int j=0; j<NWAVELENGTH; j++) CHERENKOV_PROBABILITY[j] /= CHERENKOV_PROBABILITY[NWAVELENGTH-1];
  printf("CHERENKOV_FACTOR %6.1f\n",CHERENKOV_FACTOR);
  for(int j=0; j<NWAVELENGTH; j++) printf("%2d  %6.1f %6.3f %6.3f %6.3f\n",
     j,WATER_ABSORPTION_LENGTH[j],LINER_REFLECTIVITY[j],PMT_QUANTUM_EFFICIENCY[j],CHERENKOV_PROBABILITY[j]);

  fscanf(tank,"%*s %lf",&TIMEUNIT);
  fscanf(tank,"%*s %d",&NTIME);
  fscanf(tank,"%*s %d",&TIMESLOT);
  fscanf(tank,"%*s %d",&NSLOT);
  //printf("%d %lf %d\n",NTIME,TIMESLOT,NSLOT);

  fscanf(tank,"%*s %lf",&FADC_PER_PE);
  fscanf(tank,"%*s %lf",&SIGMA_PER_PE);
  fscanf(tank,"%*s %lf",&DYNODE_ANODE_RATIO);
  fscanf(tank,"%*s %lf",&SIGMA_BASE);
  fclose(tank);
}
