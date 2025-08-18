#define minval(a,b) a<b ? a : b
#define maxval(a,b) a>b ? a : b
#define sqr(x) ((x)*(x))
//#define frand() ((double)random()/RAND_MAX)
#define frand() ((double)((arc4random()%((unsigned)RAND_MAX+1)))/RAND_MAX)
#define grand() sqrt(-2*log(frand()))*cos(2*M_PI*frand())

#define LINERSAG  .05

#define EL_NSTEP 46
double EL_ESTEP[EL_NSTEP] = 
 {   1.,    .9,    .8,    .7,    .6,    .5,    .4,    .3,    .2,    .1,
    .08,   .06,   .05,   .04,  .035,  .030,  .025,   .02,  .015,  .010,
   .008,  .006,  .005,  .004, .0035, .0030, .0025, .0020, .0015, .0010,
  .0009, .0008, .0007, .0006,.00055, .0005,.00045,.00040,.00038,.00036,
 .00034,.00032,.00030,.00029,.00028,.00027};
double EL_DEDX[EL_NSTEP] = 
 { 2.75, 2.50, 2.25, 1.95, 1.70, 1.40, 1.15, 0.85, 0.60, 0.34,
    0.3, 0.26, 0.24, 0.23,0.225,0.220,0.215,0.210,0.205,0.205,
    0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2, 
    0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2,  0.2, 
    0.2,  0.2,  0.2,  0.2,  0.2,  0.2};

#define MU_NSTEP 26
double MU_ESTEP[MU_NSTEP] =  
 {100., 50.,  20., 10.,   5.,  2.,  1.,  .8,  .6,  .5,
    .4,  .3,  .25,  .2,  .18, .16, .14, .12, .10, .09,
   .08, .07, .065, .06, .055, .05};
double MU_DEDX[MU_NSTEP] =
  {.240,.235,.230,.225,.215,.200,.192,.188,.185,.183,
   .181,.180,.180,.182,.184,.186,.190,.196,.204,.210,
   .218,.228,.235,.242,.252,.263};

/*------------------------------------------------------------------------------------------------------------*/
int Poisson(double lambda)
{
  double f,s,t;
  int i,npoiss;
  int max=30;

  //  gaussian approximation for large value of a
  if(lambda>max/2.) return lambda+grand()*sqrt(lambda)+.5;
  else {
    f=exp(lambda)*frand();
    npoiss=0;
    s=0.;
    t=1.;
    for(i=0; i<max; i++) {
      s+=t;
      if(f<s) return i;
      t*=lambda/(i+1);
    }
    return max-1;
  }
}
/*------------------------------------------------------------------------------------------------------------*/
void BuildFrame(double wx,double wy,double wz,double& ux,double& uy,double& uz,double& vx,double& vy,double& vz)
{
  // from a unit vector w, builds in the plane perp. to w two orthogonal unit vectors u and v, with v in the horizontal plane 
  double wxy = sqrt(wx*wx+wy*wy);
  // if w is not vertical, unique solution for v
  if(wxy>0) { vx = -wy/wxy; vy = wx/wxy; }
  // if w is vertical, arbitrary choice
  else { vx = 0; vy = 1; }
  vz = 0;
  // builds u to make (u,v,w) orthonormal
  ux = vy*wz-vz*wy; uy = vz*wx-vx*wz; uz = vx*wy-vy*wx;
}
/*--------------------------------------------------------------------------------------------*/
#if CURVEDLINER
double DistToWall(double x, double y, double z, double cx, double cy, double cz, double z_bot, double z_top, int& wall)
{
  double aa = cx*cx+cy*cy;
  double bb = x*cx+y*cy;
  double cc = TANK_RADIUS*TANK_RADIUS-x*x-y*y;
  // distance to side cylinder
  double length;
  if(aa>0) length = (-bb+sqrt(bb*bb+aa*cc))/aa;
  // exact vertical direction ------------
  else {
    if(cz>0) { wall = 1; return z_top-z; }
    else { wall = 3; return z-z_bot; }
  }
  // distance to curved top
  if(cz>0) {
    double ltop, xtop, ytop, ztop, xtopp, ytopp;
    ztop = z_top;
    ltop = (ztop-z)/cz;
    xtopp = x + ltop*cx; ytopp = y + ltop*cy;
    //printf("cz_up\n");
    if(sqrt(xtopp*xtopp+ytopp*ytopp)<TANK_RADIUS) {
      //printf("curv  %7.4f %7.4f %7.4f    %7.4f %7.4f %7.4f\n",x,y,z,cx,cy,cz);
      int niter = 0;
      while(1) {
        ztop = z_top+LINERSAG*sin(M_PI*(sqrt(xtopp*xtopp+ytopp*ytopp)-.15)/(TANK_RADIUS-.15));
        ltop = (ztop-z)/cz;
        xtop = x + ltop*cx; ytop = y + ltop*cy;
        //printf("      %7.4f %7.4f %7.4f\n",xtop,ytop,ztop);
        niter++;
        if(niter>10||sqr(xtop-xtopp)+sqr(ytop-ytopp)<1e-6) break;
        xtopp = xtop; ytopp = ytop;
      }
      if(niter>10) {
        printf("curv  %7.4f %7.4f %7.4f    %7.4f %7.4f %7.4f\n",x,y,z,cx,cy,cz);
        ztop = z_top;
        ltop = (ztop-z)/cz;
        xtopp = x + ltop*cx; ytopp = y + ltop*cy;
        for(int k=0; k<niter; k++) {
          ztop = z_top+LINERSAG*sin(M_PI*(sqrt(xtopp*xtopp+ytopp*ytopp)-.15)/(TANK_RADIUS-.15));
          ltop = (ztop-z)/cz;
          xtop = x + ltop*cx; ytop = y + ltop*cy;
          //printf("      %7.4f %7.4f %7.4f  %7.4f %7.4f\n",xtopp,ytopp,ztop,xtop,ytop);
          xtopp = xtop; ytopp = ytop;
        }
      }
    }
  }
  if(length<-1e-14) fprintf(diag,"length %19.16f  x,y,z,r %7.4f %7.4f %7.4f %7.4f  cx,cy,cz %7.4f %7.4f %7.4f \n",length,x,y,z,sqrt(x*x+y*y),cx,cy,cz);
  // non-vertical direction --------------
  // downwards
  if(cz<0) {
    if(length>(z_bot-z)/cz) { wall = 3; return (z_bot-z)/cz; }
    else { wall = 2; return length; }
  }
  // upwards
  if(cz>0) {
    if(length>(z_top-z)/cz) { wall = 1; return (z_top-z)/cz; }
    else { wall = 2; return length; }
  }
  return 0;
}
#else
double DistToWall(double x, double y, double z, double cx, double cy, double cz, double z_bot, double z_top, int& wall)
{
  double aa = cx*cx+cy*cy;
  double bb = x*cx+y*cy;
  double cc = TANK_RADIUS*TANK_RADIUS-x*x-y*y;
  // distance to side cylinder
  double length;
  if(aa>0) length = (-bb+sqrt(bb*bb+aa*cc))/aa;
  // exact vertical direction ------------
  else {
    if(cz>0) { wall = 1; return z_top-z; }
    else { wall = 3; return z-z_bot; }
  }
  if(length<-1e-14) fprintf(diag,"length %19.16f  x,y,z,r %7.4f %7.4f %7.4f %7.4f  cx,cy,cz %7.4f %7.4f %7.4f \n",length,x,y,z,sqrt(x*x+y*y),cx,cy,cz);
  // non-vertical direction --------------
  // bottom exit
  if(cz<0 && length>(z_bot-z)/cz) { 
    wall = 3;
    return (z_bot-z)/cz;
  }
  // top exit
  if(cz>0 && length>(z_top-z)/cz) {
    wall = 1;
    return (z_top-z)/cz;
  }
  // side exit
  wall = 2;
  return length;
}
#endif
/*--------------------------------------------------------------------------------------------*/
void Scatter(double cx, double cy,double cz, double costh, double& cxs, double& cys, double& czs)
{
  double sinth = sqrt(1-costh*costh);
  double phi = 2*M_PI*frand();
  double ux = sinth*cos(phi);
  double uy = sinth*sin(phi);
  double uz = costh;
  double cxy = sqrt(cx*cx+cy*cy);
  if(cxy==0) {
    cxs = ux;
    cys = uy;
    czs = uz;
    if(cz<0) czs = -czs;
  }
  else {
    cxs = (cz*cx*ux-cy*uy)/cxy+cx*uz;
    cys = (cz*cy*ux+cx*uy)/cxy+cy*uz;
    czs =   -cxy*ux           +cz*uz;
  }
}
/*--------------------------------------------------------------------------------------------*/
void MultipleScattering(double& cx, double& cy,double& cz, double beta, double mom, double step)
{
  double th_ms,th1,th2,ux,uy,uz,uxy,vx,vy,vz;
  th_ms = CONST_MULTSCAT/(beta*mom)*sqrt(step/WATER_X0);
  th1 = th_ms*grand();
  th2 = th_ms*grand();
  ux = -cy; uy = cx; uz = 0; uxy = sqrt(ux*ux+uy*uy);
  if(uxy>0) { ux /= uxy; uy /= uxy; }
  else { ux = 1; uy = 0; }
  vx = cy*uz-cz*uy;
  vy = cz*ux-cx*uz;
  vz = cx*uy-cy*ux;
  cx = cx+th1*ux+th2*vx;
  cy = cy+th1*uy+th2*vy;
  cz = cz+th1*uz+th2*vz;
  double norm = sqrt(cx*cx+cy*cy+cz*cz);
  cx /= norm;
  cy /= norm;
  cz /= norm;
}
/*-------------------------------------------------------------------------------*/
double  CalcPhiUNudNu(double U, double Nu, double dNu, double A, double Z)
{
#define ALPHA (1/137.0)
#define NAVAGADRO 6.02e+23
#define R0 2.82e-13

  double  FuncConst,GammaScreen,PhiUNudNu=0.;
  double  term1,term2,f1gam,f2gam,cgam;

  FuncConst = 4*ALPHA*(NAVAGADRO/A)*(Z*Z+Z)*R0*R0;
  GammaScreen = 100*(EL_MASS/U)*(Nu/(1-Nu))*pow(Z,-0.33);
  
  if(GammaScreen<0.001) {
    term1 = (1+pow(1-Nu,2)-(2./3)*(1-Nu))*log(183*pow(Z,-(1./3)));
    term2 = (1./9)*(1-Nu);
    PhiUNudNu = FuncConst*(dNu/Nu)*(term1+term2);
  }
  if(GammaScreen>=0.001&&GammaScreen<2) {
    if(GammaScreen<0.8) {
      f1gam = -3.52*GammaScreen+20.79;
      f2gam = -2.86*GammaScreen+20.29;
    }
    else {
      f1gam = f2gam = -1.91*GammaScreen+19.40;
    }
    term1 = (1+pow(1-Nu,2))*((f1gam/4)-(1./3)*log(Z));
    term2 = -((2./3)*(1-Nu))*((f2gam/4)-(1./3)*log(Z));
    PhiUNudNu = FuncConst*(dNu/Nu)*(term1+term2);
  }
  if(GammaScreen>=2&&GammaScreen<=15) {
    cgam = 0.5*exp(-0.45*GammaScreen)+0.01;
    term1 = 1+pow(1-Nu,2)-(2./3)*(1-Nu);
    term2 = log((2*U/EL_MASS)*((1-Nu)/Nu))-(1./2)-pow(cgam,1);
    PhiUNudNu = FuncConst*(dNu/Nu)*term1*term2;
  }
  if(GammaScreen>15) {
    term1 = 1+pow(1-Nu,2)-(2./3)*(1-Nu);        
    term2 = log((2*U/EL_MASS)*((1-Nu)/Nu))-(1./2);
    PhiUNudNu = FuncConst*(dNu/Nu)*term1*term2;
  }
  if(GammaScreen<0)  PhiUNudNu = 0;  
  
  return(PhiUNudNu);
}
/*-------------------------------------------------------------------*/
#define NDIV 100
double BremsTable[EL_NSTEP][NDIV];
void BuildTables()
{
  for(int i=0; i<EL_NSTEP; i++) {
    double e = EL_ESTEP[i];
      //T=0.1*pow(10,(i+0.5)*0.2);
    double U = e + EL_MASS;
    double MaxNu = 1-(EL_MASS/U);
    double dNu = MaxNu/NDIV;
    double IntPhi = 0;
      
    for(int j=1; j<NDIV; j++) { 
      double Nu = j*dNu;
      double PhiOxygen = CalcPhiUNudNu(U,Nu,dNu,15.999,8.0);
      double PhiHydrogen = CalcPhiUNudNu(U,Nu,dNu, 1.008,1.0);
      double Phi = PhiOxygen*(16.0/18)+PhiHydrogen*(2.0/18);
      IntPhi += Phi;
      BremsTable[i][j] = IntPhi;
    }
    BremsTable[i][0] = IntPhi;      
  }
  
  for(int i=0; i<EL_NSTEP; i++) {
    for(int j=1; j<NDIV; j++) {
      BremsTable[i][j] /= BremsTable[i][0];
    }
  }
}
/*--------------------------------------------------------------------*/
#define EMINBREM .100
double Bremsstrahlung(double e,int istep,double step,double beta,double& costh)
{
  double prb = frand();
  int ind = 0;
  for(int i=1; i<NDIV; i++) {
    if(BremsTable[istep][i]>prb) {ind=i; break;}
  }
  double frac = 1.*ind/NDIV;
  double egam = e*frac;
  if(egam<EMINBREM*step) return 0;  

  // Determine approx RMS gamma emission angle
  // (ref Rossi page 53 with q factor=1)
  // In fact angles are small compared to scattering so no problem
  
  double mom2 = e*(e+2*EL_MASS);
  double sig2 = pow(.0136/beta,2) / mom2 * step/.361;
  double urandom = frand();
  while(urandom==0) urandom = frand();
  costh = cos(sqrt(-2*sig2*log(urandom)));
  return egam;
}
