#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include "VBBinaryLensingLibrary.h"
using namespace std;

///**************** constants *********************//
const int YZ=3615;
const int Num=2000;
const double MaxD=20.0;///kpc
const double RA=180.0/M_PI;
const double step=MaxD/(double)Num/1.0;///step in kpc
const int Ntem=int(8300);//11979);  
const int M=5;
const int Nog=988;
const int Nkm=1253; 
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
const double pi= M_PI;
const double Pc=3.0857*pow(10.0,16);///[m]
const double G= 6.67408*pow(10.0,-11);/// [m^3/kg.s^2]
const double velocity =299792458.0;///[m/s]
const double Msun=1.989*pow(10.0,30.0);///[kg]
const double Rsun=6.9551*pow(10.0,8);///[m]
const double Au =1.495978707*pow(10.0,11);///[m]
const double binary_fraction=double(2.0/3.0);
const double vro_sun=226.0;
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
const double R_sun=8.0;
const double Avks=double(8.20922); 
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]=  {3.1,2.5,3.1,3.1};
const double FWHM[M]={3.5*0.26, 3.5*0.26, 3.5*0.26, 3.5*0.26 ,3.5*0.26};//UBVRI_OGLE
const double AlAv[M]={1.555, 1.332, 1.009, 0.841, 0.600};///From besancon model[UBVI]
const double sigma[M]={0.022, 0.022, 0.02, 0.023 ,0.025};//MOAاستفاده از مقاله کاردلی 
const int N1=19744, N2=40000, N3=4404, N4=2597;///CMD_BESANCON_new, thinD, bulge, thickD, halo
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
struct pulsing{
   double peri;///period
   double phir, phit;
   double epsiR, epsiT; ///amplitude
   double omega; 
   double inc;  
   int l,m;
   double ampl; 
};
struct lens{
    double RE, tE, Vt, Dl, Ml;
    double t0, ksi, u0, pt1, pt2,tstar;
    double proj,dt1, dt2;
    double vl, vs, xls;
    double rhomaxl;
    int numl, struc;
};
struct source{
    int nums, struc, cl; 
    double lon, lat, type;
    double mass,Ds, logl,logg, Type;
    double Lstar[M], Tstar0, Rstar0, ro_star0;
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart,nstarti;
    double FI,TET;
    double Mab[M],Map[M], col;
    double Fluxb[M],magb[M], Fluxb2[M]; 
    double blend[M], blend2[M], nsbl[M]; 
    double ext[M];
    double Tem[Ntem], Mm[M][Ntem];
};
struct CMD{
    double Teff_d[N1],logl_d[N1],Mab_d[M][N1],mass_d[N1], type_d[N1], metal_d[N1], gra_d[N1],Rs_d[N1]; int cl_d[N1];//thindisk
    double Teff_b[N2],logl_b[N2],Mab_b[M][N2],mass_b[N2], type_b[N2], metal_b[N2], gra_b[N2],Rs_b[N2]; int cl_b[N2];// bulge
    double Teff_t[N3],logl_t[N3],Mab_t[M][N3],mass_t[N3], type_t[N3], metal_t[N3], gra_t[N3],Rs_t[N3]; int cl_t[N3];//thickdisk
    double Teff_h[N4],logl_h[N4],Mab_h[M][N4],mass_h[N4], type_h[N4], metal_h[N4], gra_h[N4],Rs_h[N4]; int cl_h[N4];// halo
};
struct extinc{
   double dis[100];///distance 
   double Extks[100];///ks-band extinction
   double Ai[M],Av,Aks,AI;
   double exks;
};
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_lens(source &s,lens & l);
void Func_source(source & s, CMD & cm, extinc & ex);
void Func_puls(lens & l, pulsing & p);
double RandN(double N, double sigma);
double randR(double down, double up);
double TETA(double y, double x);  
double Pfunc(int l, int m, double tets);///Legandre function
int Extinction(extinc & ex,source & s);
void read_cmd(CMD & cm);
void vrel(source & s,lens & l);
void Disk_model(source & s);
double Interpol(double ds, extinc & ex);
double ErrorCal(double mag, int sur);
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    time_t _timeNow;
    unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
///=========================================
int main()
{
//================================================
   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
///=================================================
   printf("START_TIME >>>>>>>> %s",ctime(&_timeNow));
   pulsing p;
   lens l;
   source s; 
   CMD cm; 
   extinc ex;
   read_cmd(cm);
   
   
   
  
   printf("START_TIME >>>>>>>> %s",ctime(&_timeNow));
   VBBinaryLensing *vbb=new VBBinaryLensing;
   vbb->Tol=1.e-2;
   vbb->LoadESPLTable("./files/ESPL.tbl");
   
   char filname1[40], filname2[40], filname3[40], filname4[40]; 
   FILE* ogle;     FILE* kmtnet; 
   FILE* param;    FILE* magm;    
   FILE* danish;   FILE* tems;    
   FILE* lumino;   FILE* sour;
   param=fopen("./files/MONT/NRadial1L/NRadial1L_parama.txt","a+");
   fclose(param); 




   double dto[Nog],  dtk[Nkm], dmo, magg; 
   ogle=fopen("./files/OGLE_delta.dat","r");
   for (int i=0; i<Nog; ++i) {fscanf(ogle,"%lf   %lf   %lf \n", &dto[i], &dmo, &magg); }
   kmtnet=fopen("./files/KMTNet_delta_2.dat","r");
   for (int i=0; i<Nkm; ++i) {fscanf(kmtnet,"%lf   %lf   %lf \n",&dtk[i], &dmo, &magg); }
   tems= fopen("./files/FstarT_Bessel_new.dat","r"); 
   for(int i=0; i<Ntem; ++i){
   fscanf(tems,"%lf   %lf   %lf   %lf   %lf   %lf\n",&s.Tem[i],&s.Mm[0][i],&s.Mm[1][i],&s.Mm[2][i],&s.Mm[3][i],&s.Mm[4][i]); }
   fclose(ogle);  fclose(kmtnet);   fclose(tems); 
   cout<<"The files were read *************************"<<endl; 
   
   
   
   double chi_r1, chi_r2, chi_b1, chi_b2, chi_np, chi_re; 
   double dchi_lens, dchi_puls1, dchi_puls2, sigAt;
   double AA, AA0, dm0, improv;
   int fla[2]; 
   int ndat[2];
   int nog, nkm, j;  
   double timo, timk, chi0, chi1, dchi, xlens, ylens, Tstar;
   double tt, AtotI, deltaA, deld, peak, magp; 
   double t, u, Rstar, Lstar[M];
   double tets, phis, tets0 , phis0, Legan, ux, uy, magni, area;
   double del[M], Av[M], A0[M], Bs[M], Amagx;
   double xs, ys, zs, xo, yo, zo, rhoe, dis, dmin;
   double dtet=0.99, dphi0=0.99, dphi;
   int NN=0;
   for(tets=dtet; tets<180.0; tets+= dtet){
   for(phis=0.0;  phis<360.0; phis+=double(dphi0/sin(tets/RA)))  NN+=1;}
   NN=int(NN/2.+500.0); 
   cout<<"NN:  "<<NN<<endl;
   double Yo[NN], Zo[NN], tet[NN], phi[NN]; 
   double dzmax, dymax, dmax, ddel, lonn; 
   int nt,j1,Nlens, flag[2], yye, Nperi, fgh, flagE; 
   double ymin, ymax, zmin, zmax, tp, sigA0, nw, ffbl;
   double aveL[M]; 
   double ave_error[3]={0.0}; 
   
   
   
   
   for(int icon=2000; icon<2200; ++icon){ 
   cout<<"***** Step:  "<<icon<<endl;
   
   s.lat=randR(-6.0,-4.0);
   lonn= randR(+3.0,+6.0);
   cout<<"latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
   if(lonn<=0.0)   s.lon=360.0+lonn;
   else            s.lon=lonn;
   s.TET=(360.0-s.lon)/RA;///radian
   s.FI=s.lat/RA;///radian
   Disk_model(s);
   flagE=Extinction(ex,s);
   cout<<"flagE:  "<<flagE<<endl;
   if(flagE==1){


   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
  
   do{
   Func_source(s, cm, ex);
   Func_lens(s,l);
   Func_puls(l,p); 
   }while(l.tE<=0.5 or l.tE>300.0 or s.magb[4]>22.0 or  s.Tstar0>17000.0 or s.Tstar0<400.0);
   cout<<"u0:      "<<l.u0<<"\t ro_star0:  "<<s.ro_star0<<endl;
   cout<<"magb[4]:  "<<s.magb[4]<<"\t Tstar0:  "<<s.Tstar0<<endl;  
   cout<<"tE: "<<l.tE<<"\t tstar:    "<<l.tstar<<"\t p.peri:  "<<p.peri<<endl;
   
   ffbl=randR(0.0,1.0);  
   if(ffbl<s.blend[4]){
   
  
    fgh=0;
   if(icon%1==0){
   fgh=1;
   sprintf(filname1,"./files/MONT/NRadial1L/%c%c%c%c%d.dat",'m','a','g','_',icon);
   magm=fopen(filname1,"w");
   sprintf(filname2,"./files/MONT/NRadial1L/%c%c%c%c%d.dat",'d','a','t','_',icon);
   danish=fopen(filname2,"w");
   sprintf(filname3,"./files/MONT/NRadial1L/%c%c%c%c%d.dat",'l','u','m','_',icon);
   lumino=fopen(filname3,"w");
   sprintf(filname4,"./files/MONT/NRadial1L/%c%c%c%c%d.dat",'s','o','u','_',icon);
   sour=fopen(filname4,"w"); }
  
////==================================================================


    l.dt1=double(5.0/60.0/24.0);//5 min as time interval between data
    Nlens= int(l.pt2-l.pt1)/l.dt1 + 2;
    
    Nperi=100; 
    l.dt2=double(p.peri/Nperi);
    if(l.dt2<l.dt1){
    l.dt2=l.dt1;
    Nperi= int(p.peri/l.dt2);}   
    double lumi[M][Nperi]; 
    double *flensp=new double[Nlens]; 
    int *clens=new int[Nlens]; 
    double **Mto0=new double*[M];
    double **Mto1=new double*[M];
    for(int i=0; i<M; ++i){
    Mto0[i] =new double[Nlens];
    Mto1[i] =new double[Nlens];}///Mto0[M][Nlens], Mto1[M][Nlens]
    for(int i=0; i<Nlens; ++i){
    flensp[i]=0;
    for(int j=0; j<M; ++j){
    Mto0[j][i]=Mto1[j][i]=0.0;}}
    cout<<"arrayes are made Nlens:   "<<Nlens<<endl;
   
   
   int Npr=0; 
   for(int nl=0; nl<Nlens; ++nl){
   clens[nl]=Npr; 
   Npr+=1;
   if(Npr==Nperi)  Npr=0; }
  
  
  
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   ave_error[0]=ave_error[1]=ave_error[2]=0.0; 
   aveL[0]=aveL[1]=aveL[2]=aveL[3]=aveL[4]=0.0;  p.ampl=0.0;
   for(int Npr=0;  Npr<Nperi;  Npr+=1){////time of pulsation   
   tp=double(Npr*l.dt2); 
   lumi[0][Npr]=lumi[1][Npr]=lumi[2][Npr]=lumi[3][Npr]=lumi[4][Npr]=0.0; 

   
   if(Npr==0){
   phis=0.0;  
   tets= pi/2.0-p.inc-dtet/RA;//the points with ys=0
   Legan= Pfunc(p.l,p.m,tets);
   Rstar=s.Rstar0*(1.0+p.epsiR*Legan*cos(p.omega*tp+ p.m*phis + p.phir));
   xs=Rstar*sin(tets)*cos(phis);
   zs=Rstar*cos(tets);
   dzmax=fabs(cos(p.inc)*zs-sin(p.inc)*xs); ///zo-0
   tets=pi/2.0-p.inc;
   dphi=double(dphi0/RA/sin(tets)); 
   phis += dphi;//radian
   Legan= Pfunc(p.l,p.m,tets);
   Rstar=s.Rstar0*(1.0+p.epsiR*Legan*cos(p.omega*tp+ p.m*phis + p.phir));
   dymax=fabs(Rstar*sin(tets)*sin(phis));///yo-0
   if(tets==0 or p.inc==pi/2.0) dmax= dzmax;
   else{
   if(dymax>dzmax)  dmax=dymax;
   else             dmax=dzmax;}
   ddel=fabs(double(dmax*0.9) );
   area= fabs(ddel*ddel);
   rhoe= sqrt(area/pi)*l.proj;
   cout<<"dzmax:  "<<dzmax<<"\t dymax:  "<<dymax<<"\t dmax: "<<dmax<<endl;}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      
      
   nt=0; ymin=0.0, ymax=0.0, zmin=0.0, zmax=0.0; 
   for(int i=0; i<NN; ++i){Yo[i]=Zo[i]=tet[i]=phi[i]=0.0;}
   for(tets=dtet; tets<180.0; tets+=dtet){///degree
   Legan= Pfunc(p.l,p.m,tets/RA);
   if(sin(tets/RA)>0.0)  dphi=double(dphi0/fabs(sin(tets/RA)));
   else{cout<<"ERROR  tets(degree):  "<<tets<<"\t sin(tets):  "<<sin(tets/RA)<<endl;  cin>>yye;}  
   for(phis=0.0;  phis<360.0; phis+=dphi){///degree 
   Rstar=fabs(s.Rstar0*(1.0+p.epsiR*Legan*cos(p.omega*tp + p.m*phis/RA + p.phir)));
   xs=Rstar*sin(tets/RA)*cos(phis/RA);
   ys=Rstar*sin(tets/RA)*sin(phis/RA);
   zs=Rstar*cos(tets/RA);
   xo=cos(p.inc)*xs+sin(p.inc)*zs;
   yo=ys;
   zo=cos(p.inc)*zs-sin(p.inc)*xs;
   if(xo>0.0 or xo==0.0){
   Yo[nt]=yo;  Zo[nt]=zo; tet[nt]=tets/RA;   phi[nt]= phis/RA; 
   if(yo>ymax) ymax=yo;  if(yo<ymin) ymin=yo; 
   if(zo>zmax) zmax=zo;  if(zo<zmin) zmin=zo; 
   nt+=1;}
   if(nt>(NN-1)){cout<<"ERROR nt: "<<nt<<"\t NN: "<<NN<<"\t tets:  "<<tets<<"\t phis:  "<<phis<<"\t inc: "<<p.inc<<endl; exit(0);}}}
   
   
 
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH    
   for(yo=double(ymin-3.0*ddel);  yo<=double(ymax+3.0*ddel);  yo += ddel){
   for(zo=double(zmin-3.0*ddel);  zo<=double(zmax+3.0*ddel);  zo += ddel){
   flag[0]=flag[1]=0; 
   dmin=1000000.0; 
   for(int i=0; i<nt; ++i){   
   dis= sqrt((yo-Yo[i])*(yo-Yo[i])+(zo-Zo[i])*(zo-Zo[i])); 
   if(dis<dmin){ dmin=dis;  tets=tet[i];  phis=phi[i];}}
   if(fabs(dmin)<=fabs(dmax)){
   flag[1]=1;
   Rstar=s.Rstar0*(1.0+p.epsiR*Pfunc(p.l,p.m,tets)*cos(p.omega*tp + p.m*phis + p.phir));
   Tstar= s.Tstar0 +   p.epsiT*Pfunc(p.l,p.m,tets)*cos(p.omega*tp + p.m*phis + p.phit);  
   j1=-1; 
   if(Tstar<s.Tem[0] or Tstar==s.Tem[0])  j1=0;   
   else if(Tstar>s.Tem[Ntem-1])         j1=int(Ntem-1); 
   else{
   for(int k=1; k<Ntem; ++k){
   if((Tstar-s.Tem[k])*(Tstar-s.Tem[k-1])<0.0 or Tstar==s.Tem[k]){ j1=k;   break;}}}
   for(int k=0; k<M; ++k){
   Lstar[k]= pow(10.0,-0.4*s.Mm[k][j1])/s.Lstar[k];}} 
  
  
   if(sqrt(yo*yo+zo*zo)<s.Rstar0 or sqrt(yo*yo+zo*zo)==s.Rstar0){
   if(sqrt(yo*yo+zo*zo)==s.Rstar0)     xo=0.0; 
   else    xo=sqrt(fabs(s.Rstar0*s.Rstar0-yo*yo-zo*zo));
   xs=cos(p.inc)*xo-sin(p.inc)*zo;
   ys=yo;
   zs=sin(p.inc)*xo+cos(p.inc)*zo; 
   tets0= TETA(sqrt(xs*xs+ys*ys),zs);
   phis0= TETA(ys,xs); 
   if(phis0==2.0*pi) phis0=0.0; 
   flag[0]=1;}
  
  
   if( (flag[0]+flag[1])>1 and (Tstar>12000.0 or Tstar<21.0 or fabs(Rstar-s.Rstar0)>s.Rstar0 or j1<0 or j1>=Ntem or
   tets>pi or tets<0.0 or phis<0.0 or phis>=(2.0*pi) or fabs(tets-tets0)>double(90.0/RA) or 
   tets0>pi or tets0<0.0 or phis0<0.0 or phis0>(2.0*pi) or Lstar[1]>10.0 or Lstar[1]<0.0)){
   cout<<"ERROR***   yo: "<<yo<<"\t zo:  "<<zo<<"\t xo:  "<<xo<<endl; 
   cout<<"Rstar:  "<<Rstar<<"\t Rstar0: "<<s.Rstar0<<endl; 
   cout<<"Tstar: "<<Tstar<<"\t Tstar0:  "<<s.Tstar0<<endl;
   cout<<"tets0:  "<<tets0*RA<<"\t phis0:  "<<phis0*RA<<endl;
   cout<<"tets:  "<<tets*RA<<"\t phis:  "<<phis*RA<<endl; cin>>yye;}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH    
   if( (flag[1]>0 or flag[0]>0) and  Npr==0){
   if(fgh==1){  
   sour=fopen(filname4,"a+");
   fprintf(sour,"%d  %d   %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf\n",
   flag[0],flag[1],yo,zo,Rstar,Tstar,tets*RA,tets0*RA,phis*RA,phis0*RA);
   fclose(sour);}}
     
   
   if(flag[1]>0){
   for(int k=0; k<M; ++k){
   lumi[k][Npr] +=Lstar[k]*area/(pi*s.Rstar0*s.Rstar0); 
   aveL[k] +=     Lstar[k]*area/(pi*s.Rstar0*s.Rstar0);}
   for(int nl=0;  nl<Nlens;  ++nl){
   //nw=double(nl*l.dt1);  
   //while(nw>p.peri or nw==p.peri){nw=double(nw-p.peri);}
   //if(nw==p.peri)  nw=0.0;  
   //if(nw>p.peri or nw<0.0 or p.peri<=0.0){cout<<"ERROR nw:  "<<nw<<"\t peri:   "<<p.peri<<endl;  int wwe; cin>>wwe;}
   //if((nw-tp)*(nw-tp-l.dt2)<0.0 or nw==tp){
   if(clens[nl]==Npr){
   flensp[nl]=1;
   t=double(l.pt1 + nl*l.dt1);
   xlens = (t-l.t0)/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi);
   ylens = (t-l.t0)/l.tE * sin(l.ksi) + l.u0 * cos(l.ksi);
   ux= xlens - yo*l.proj;
   uy= ylens - zo*l.proj;
   u=sqrt(ux*ux+ uy*uy); 
   if(u>double(10.0*rhoe))  magni=fabs( (u*u+2.0)/sqrt(u*u*(u*u+4.0)) ); 
   else                     magni=vbb->ESPLMag2(u,rhoe);  
   for(int k=0; k<M; ++k){
   Mto1[k][nl] += Lstar[k]*magni*area/(pi*s.Rstar0*s.Rstar0);
   Mto0[k][nl] += Lstar[k]*1.000*area/(pi*s.Rstar0*s.Rstar0);
   if(magni<1.0 or Lstar[k]<0.0 or Lstar[k]>10.0 or area<=0.0 or fabs(Mto1[k][nl])<fabs(Mto0[k][nl])){
   cout<<"Error magni:  "<<magni<<"\t Lstar[k]:  "<<Lstar[k]<<"\t area:  "<<area<<endl;
   cout<<"Mto1:  "<<Mto1[k][nl]<<"\t Mto0:  "<<Mto0[k][nl]<<endl;   } 
   }}
   }} }} 
   //if(Npr%10==0)   
   cout<<"tp:  "<<tp<<"\t Npr:   "<<Npr<<"\t lumi:  "<<lumi[4][Npr]<<endl;   
   for(int k=0;k<M; ++k) {
   cout<<"filter:  "<<k<<"\t lumi[k]:  "<<lumi[k][Npr]<<endl; }
   cout<<"***************************************************"<<endl;
   }///time loop
   

   p.ampl=-1000.0;  
   for(int k=0; k<M; ++k){ 
   aveL[k]= double(aveL[k]/Nperi/1.00);
   s.blend2[k]= fabs(aveL[k]/(aveL[k] + s.Fluxb2[k])); 
   if(s.blend2[k]<=0.0 or s.blend2[k]>1.0 or aveL[k]<0.0  or aveL[k]==0.0){
   cout<<"Error aveL[k]:  "<<aveL[k]<<"\t Fluxb2:  "<<s.Fluxb2[k]<<"\t blend2: "<<s.blend2[k]<<endl;
   cout<<"**********************************************************************************"<<endl;
   int uue; cin>>uue;} }
   
   
   for(int Npr=0;  Npr<Nperi;  Npr+=1){
   tp=double(Npr*l.dt2); 
   for(int k=0; k<M; ++k)        lumi[k][Npr]= -2.5*log10(lumi[k][Npr]/aveL[k]);
   if(fabs(lumi[4][Npr])>p.ampl) p.ampl=fabs(lumi[4][Npr]);
   if(fgh==1) 
   fprintf(lumino,"%.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf\n",tp,lumi[0][Npr],lumi[1][Npr],lumi[2][Npr],lumi[3][Npr],lumi[4][Npr]);}
    if(fgh==1) fclose(lumino);  
  
    
    
    
    
    
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH 
   timo=timk=0.0;
   chi_r1=0.0, chi_r2=0.0, chi_b1=0.0, chi_b2=0.0, chi_np=0.0, chi_re=0.0; 
   ndat[0]=ndat[1]=0;
   nog=int((double)rand()*(Nog-5)/(double)(RAND_MAX+1.))+1;
   nkm=int((double)rand()*(Nkm-5)/(double)(RAND_MAX+1.))+1;
   for(int nl=0; nl<Nlens; ++nl){
   if(flensp[nl]>0){
   t=double( l.pt1 + nl*l.dt1);
   xlens = (t-l.t0)/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi);
   ylens = (t-l.t0)/l.tE * sin(l.ksi) + l.u0 * cos(l.ksi);
   u=sqrt(xlens*xlens+ ylens*ylens); 
   
   if(u>double(10.0*s.ro_star0)) Amagx=fabs((u*u+2.0)/sqrt(u*u*(u*u+4.0)) ); 
   else                          Amagx=fabs(vbb->ESPLMag2(u,s.ro_star0) );
   for(int k=0;  k<M;  ++k){
   A0[k]=fabs(Amagx*s.blend2[k] + 1.0- s.blend2[k]);///magnification without pulsation
   Av[k]=fabs(Mto1[k][nl]/aveL[k]) * s.blend2[k] + 1.0- s.blend2[k];//Magnification with pulsation
   Bs[k]=fabs(Mto0[k][nl]/aveL[k]) * s.blend2[k] + 1.0- s.blend2[k];// Base_line 
   del[k]=-2.5*log10(fabs(Av[k]/A0[k]));
   if(A0[k]<1.0 or Amagx<1.0 or Av[k]<0.0 or Bs[k]<0.0 or s.blend2[k]<0.0 or s.blend2[k]>1.0 or fabs(Mto1[k][nl])<fabs(Mto0[k][nl]) ){
   cout<<"ERROR    Amagx:   "<<Amagx<<"\t  A0:  "<<A0[k]<<"\t blend2:  "<<s.blend2[k]<<endl;
   cout<<"Av:    "<<Av[k]<<"\t Bs:  "<<Bs[k]<<"\t del:   "<<del[k]<<endl; 
   cout<<"Mto1:   "<<Mto1[k][nl]<<"\t  Mto0:   "<<Mto0[k][nl]<<endl;  }}
   tp=double(nl*l.dt1)/p.peri;
//=================================================================   
   timo+= double(l.dt1); 
   timk+= double(l.dt1);
   fla[0]=fla[1]=-1;
   if(timo>dto[nog]){//OGLE data
   timo -= dto[nog];
   fla[0]=1; 
   nog=int((double)rand()*(Nog-3)/(double)(RAND_MAX+1.))+1; }
   if(timk>dtk[nkm]){//KMTNet data
   timk -= dtk[nkm];
   fla[1]=1; 
   nkm=int((double)rand()*(Nkm-3)/(double)(RAND_MAX+1.))+1; }
   magg=double(s.magb[4]-2.5*log10(Av[4]));
   for(int i=0; i<2; ++i){ 
   dmo= ErrorCal(magg, i); 
   if(fla[i]>0 and magg>12.0){  
   ndat[i]+=1; 
   tt= RandN(3.0,1.0);
   dm0= ErrorCal(s.magb[4],i); 
   sigAt= double(fabs(pow(10.0,-0.4*dmo)-1.0)*Av[4]);
   sigA0= double(fabs(pow(10.0,-0.4*dm0)-1.0)*1.000);
   deltaA= sigAt * tt;  
   AtotI= Av[4] + deltaA;
   AA0=  -2.5*log10(Av[4]/A0[4]);//Lstar/Lave 
   deld= -2.5*log10(AtotI/A0[4]);
   chi_b1 +=(deld-0.0)*(deld-0.0)/(dm0*dm0*2.0);//  (AA-1.0)*(AA-1.0)/(sigL0*sigL0); 
   chi_r1 +=(deld-AA0)*(deld-AA0)/(dm0*dm0*2.0);//(AA-AA0)*(AA-AA0)/(sigL0*sigL0); 
   if(fabs(u)>s.ro_star0){
   chi_b2 +=(deld-0.0)*(deld-0.0)/(dmo*dmo*2.0); 
   chi_r2 +=(deld-AA0)*(deld-AA0)/(dmo*dmo*2.0);}
   chi_np  +=  (AtotI-A0[4])*(AtotI-A0[4])/(sigAt*sigAt);////without any pulsation
   chi_re  +=  (deltaA*deltaA)/(sigAt*sigAt); ///real model 
   if(fgh==1){  fprintf(danish,"%.5lf  %.5lf  %.8lf %.5lf  %.8lf  %d\n",(t-l.t0), AtotI , sigAt, deld, dmo*sqrt(2.0), i); } 
   ave_error[0]+= fabs(dm0);
   ave_error[1]+= fabs(dmo);
   ave_error[2]+=1.0; 
   }}  
///===================================================================================   
       if(fgh==1){ fprintf(magm,"%.4lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.4lf  %.4lf  %.4lf %.4lf  %.5lf  %.5lf  %.5lf %.5lf %.7lf  %.7lf  %.7lf %.7lf %.4lf %.4lf %.5lf\n",   
   (t-l.t0),Amagx,Av[1],Av[2],Av[3],Av[4],A0[1],A0[2],A0[3],A0[4],Bs[1],Bs[2],Bs[3],Bs[4],
   del[1],del[2],del[3],del[4],xlens,ylens,tp); }
   if(nl%100==0){
   cout<<"================================================"<<endl;
   cout<<"\t **** counter:  "<<nl<<endl;
   cout<<"(t-t0)/tstar:  "<<(t-l.t0)/l.tstar<<endl;
   cout<<"xlens: "<<xlens<<"\t ylens: "<<ylens<<endl;
   cout<<"Amag0:  "<<Amagx<<"\t Av[1]: "<<Av[1]<<endl;
   cout<<"A[2]"<<Av[2]<<"\t A[3]"<<Av[3]<<"\r A[4]:  "<<Av[4]<<endl;
   cout<<"lumi[2]"<<A0[2]<<"\t lumi[3]"<<A0[3]<<"\t lumi[4]:  "<<A0[4]<<endl;
   cout<<"del[0]:  "<<del[0]<<"\t del[1]: "<<del[1]<<"del[2]:  "<<del[2]<<endl;
   cout<<"del[3]:  "<<del[3]<<"\t del[4]:  "<<del[4]<<endl; 
   cout<<"================================================"<<endl;}}
   else{cout<<"ERROR: flensp[nl]:  "<<flensp[nl]<<endl; }}
   if(fgh==1){fclose(magm);   fclose(danish); }
   
   dchi_lens= fabs(chi_re-chi_np);///How much is the pulsation perturbation? 
   dchi_puls1=fabs(chi_r1-chi_b1);///In residual without lensing how much is pulsation perturbation
   dchi_puls2=fabs(chi_r2-chi_b2);///In residual of lensing howmuch is teh pulsation perturnation 
   improv= fabs(dchi_puls1- dchi_puls2); 
   ave_error[0] = double(ave_error[0]/ave_error[2]);
   ave_error[1] = double(ave_error[1]/ave_error[2]);

   param=fopen("./files/MONT/NRadial1L/NRadial1L_parama.txt","a+");
   fprintf(param,"%d  %d  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.7lf   %.5lf  %.5lf  %d  %d  %d  %.5lf  %.5lf  %.5lf    %.5lf %.5lf  %.5lf  %.8lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %d  %d  %.5lf %.5lf %.5lf %.5lf %.5lf  %d  %d  %d  %.3lf  %.3lf %.3lf  %.3lf  %.2lf  %.5lf   %.8lf   %.8lf  %.6lf %.6lf\n",
   l.numl,l.struc,l.Ml,l.Dl,l.xls,l.RE/Au,l.Vt,l.vs,l.vl,l.tE,l.u0/s.ro_star0,l.ksi, //12
   l.tstar,s.nums,s.struc,s.cl,s.mass,s.Ds,s.logg,s.logl,s.Tstar0,s.Rstar0,s.ro_star0,  //23
   s.Lstar[1],s.Lstar[2],s.Lstar[3],s.Lstar[4],s.Map[4],s.magb[4], //29
   p.peri,p.phir*RA,p.phit*RA,p.epsiR,p.epsiT,p.inc*RA,p.l,p.m,  ///37
   lumi[0][0],lumi[1][0],lumi[2][0],lumi[3][0],lumi[4][0],icon, ///43
   ndat[0],ndat[1], dchi_lens, dchi_puls1, dchi_puls2,improv,s.type,p.ampl,ave_error[0],ave_error[1],s.blend2[2],s.blend2[4]);
   fclose(param); 
    cout<<"******************************************************"<<endl;
    cout<<"********* icon:  "<<icon<<"\t dchi:  "<<dchi<<"\t ampl:  "<<p.ampl<<endl; 
    cout<<"l:  "<<p.l<<"\t m:  "<<p.m<<"inclination:  "<<p.inc*RA<<endl;
    cout<<"LENS:   mass[Msun]: "<<l.Ml<<"\t lens_dis[kpc]:  "<<l.Dl<<"\t RE[AU]: "<<l.RE/Au<<endl;
    cout<<"xls: "<<l.xls<<"\t u0:  "<<l.u0<<"\t ksi:  "<<l.ksi<<endl;
    cout<<"tstar:  "<<l.tstar<<"\t pt1:  "<<l.pt1<<"\t pt2:   "<<l.pt2<<endl;
    cout<<"tE[days]: "<<l.tE<<"\t Vt:  "<<l.Vt<<"\t Vl:  "<<l.vl<<endl;
    cout<<"source:  mass:  "<<s.mass<<"Ds[Kpc]: "<<s.Ds<<endl;
    cout<<"mI:  "<<s.Map[4]<<"\t map: "<<s.magb[4]<<endl;
    cout<<"Rstar[Rsun]: "<<s.Rstar0<<"\t Tstar[K]:  "<<s.Tstar0<<endl;
    cout<<"ro_star: "<<s.ro_star0<<"\t u0: "<<l.u0/s.ro_star0<<endl;
    cout<<"logg:  "<<s.logg<<"\t logl:  "<<s.logl<<"\t Tstar0:  "<<s.Tstar0<<endl;
    cout<<"PULS: peri: "<<p.peri<<"\t phir: "<<p.phir*180.0/pi<<endl;
    cout<<"epsiR/Rstar: "<<p.epsiR<<"\t epsiT:  "<<p.epsiT<<endl;
    cout<<"aveL[1]:  "<<lumi[0][0]<<"\t aveL[2]:  "<<lumi[1][0]<<endl;
    cout<<"aveL[3]:  "<<lumi[2][0]<<"\t aveL[4]:  "<<lumi[3][0]<<endl; 
    cout<<"dchi_puls1:  "<<dchi_puls1<<"\t dchi_puls2:  "<<dchi_puls2<<"\t improve:  "<<improv<<endl;
    cout<<"***************************************************** "<<endl; 
    
    for(int i = 0;  i<M;  ++i){   delete [] Mto0[i], Mto1[i]; }
    delete [] flensp,  clens;}
    }}
    
    time(&_timeNow);
    printf("END_TIME >>>>>>>>  %s ",ctime(&_timeNow));
    delete vbb;
    return(0);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double ErrorCal(double mag, int sur){
   double error; 
   if(sur==0) { ///OGLE
   if(mag<15.0 ) error=0.003;   //8.62339038e-03 6.82867290e-03 2.27422710e-03 1.66494077e+01
   if(mag>=15.0) error=0.00862339038 + 0.00682867290*(mag-16.6494077) +  0.00227422710*(mag-16.6494077)*(mag-16.6494077);  
   if(error<0.003 or error >0.1) {cout<<"Error(OGLE) error:  "<<error<<"\t mag:  "<<mag<<endl;  }}
   if(sur==1) {//KMTNet
   if(mag<14.0)  error= 0.0038;//7.36447617e-02 -1.61925360e-02  9.52245559e-04  5.61612331e+00
   if(mag>=14.0) error= 0.0736447617 -0.0161925360*(mag-5.61612331) +  0.000952245559*(mag-5.61612331)*(mag-5.61612331); 
   if(error<0.0038 or error>0.1) {cout<<"Error(KMTNet)  error:  "<<error<<"\t mag:  "<<mag<<endl;   }}
   return(error);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_source                             ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_source(source & s, CMD & cm, extinc & ex)
{
    int num,struc,nums;
    double rho,rf,Nblend[M];
    double Akv,Avk,Alv;
    double Ds,Ai[M],Av;
    double Map[M];
    double Lstar[M];
    
    

    double maxnb=0.0;
    for(int i=0; i<M; ++i){
    s.Fluxb[i]=s.Fluxb2[i]=Nblend[i]=0.0;
    Nblend[i]=s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0);
    Nblend[i]=Nblend[i]+RandN(1.0,sqrt(Nblend[i]));
    if(Nblend[i]<0.0)  Nblend[i]=0.0; 
    s.nsbl[i]=  double(Nblend[i]);
    if(Nblend[i]<=1.0) Nblend[i]=1.0;
    if(Nblend[i]>maxnb) maxnb=Nblend[i]; }
    //cout<<"maxnb:  "<<maxnb<<endl;
    


    for(int k=1; k<=int(maxnb); ++k){
    do{
    num=int(fabs((double)rand()/(double)(RAND_MAX+1.)*Num*1.0));
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    }while(rho>s.Rostari[num] || num<5);///distance larger than 50.0 
    Ds=(double)(num*step);///in kpc
    nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums;}




    rf=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[nums];
         if (rf<= s.rho_disk[nums]) struc=0;///thin disk
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
    else if (rf<=s.Rostar0[nums]) struc=3;///halo
    if(k==1)    s.struc=struc;
    //cout<<"Ds:  "<<Ds<<"\t struc:  "<<struc<<endl;


   double Mab[M];   double Tstar, Rstar; 
    if(struc==0){///thin disk
    //if(k==1){
    //do{ num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-2.0)); 
     // }while(cm.rr_d[num]<1);}
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-2.0)); 
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    Tstar=cm.Teff_d[num];
    Rstar=  cm.Rs_d[num]; 
    if(k==1){ 
    s.type=cm.type_d[num];
    s.mass=cm.mass_d[num];
    s.logl=cm.logl_d[num];
    s.cl=    cm.cl_d[num]; 
    s.logg= cm.gra_d[num]; 
    s.Tstar0=cm.Teff_d[num]; 
    s.Rstar0 = cm.Rs_d[num];}}



    if(struc==1){///bulge
    //if(k==1){
    //do{ num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-2.0)); 
     /// }while(cm.rr_b[num]<1);}
     num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-2.0)); 
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num];}
    Tstar=cm.Teff_b[num];
    Rstar=  cm.Rs_b[num]; 
    if(k==1){ 
    s.type=cm.type_b[num];
    s.mass=cm.mass_b[num];
    s.logl=cm.logl_b[num];
    s.cl=    cm.cl_b[num];
    s.logg =cm.gra_b[num];
    s.Tstar0=cm.Teff_b[num]; 
    s.Rstar0 = cm.Rs_b[num]; }}
   


    if(struc==2){///thick disk
   // if(k==1){
    //do{ num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-2.0)); 
      //}while(cm.rr_t[num]<1);}
     num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-2.0)); 
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num]; }
    Tstar=cm.Teff_t[num];
    Rstar=  cm.Rs_t[num]; 
    if(k==1){ 
    s.type=cm.type_t[num];
    s.mass=cm.mass_t[num];
    s.logl=cm.logl_t[num];
    s.cl=    cm.cl_t[num];
    s.logg =cm.gra_t[num];
    s.Tstar0=cm.Teff_t[num];
    s.Rstar0 = cm.Rs_t[num];}}


    if(struc==3){///stellar halo
   // if(k==1){
   // do{ num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-2.0)); 
     // }while(cm.rr_h[num]<1);}
     num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-2.0)); 
    for(int i=0; i<M; ++i) {Mab[i]=cm.Mab_h[i][num];}
    Tstar=cm.Teff_h[num];
    Rstar=  cm.Rs_h[num]; 
    if(k==1){
    s.cl  =  cm.cl_h[num]; 
    s.type=cm.type_h[num];
    s.mass=cm.mass_h[num];
    s.logl=cm.logl_h[num];
    s.logg =cm.gra_h[num];
    s.Tstar0=cm.Teff_h[num]; 
    s.Rstar0 = cm.Rs_h[num];}}
    
    if(s.type>8.0 or s.type<2.0 or s.Rstar0<0.0 or s.Rstar0>1000.0 or s.mass<0.0 or s.Tstar0<400.0 or s.cl>=6 or s.logg<0.0){
    cout<<"ERROR:  type: "<<s.type<<"\t struc: "<<struc<<"\t num: "<<num<<endl; int rre;  cin>>rre;}
    
   
    ex.Aks=Interpol(Ds,ex);///extinction in Ks-band
    Av=ex.Aks*Avks;
    if(Av>20.0 or Av<0.0 or Ds>20.0  or Ds<0.0){cout<<"ERROR Ds:  "<<Ds<<" \t Av:  "<<Av<<endl; int yyw;  cin>>yyw; }
    if(Av<0.0)    Av=0.0;

    for(int i=0; i<M; ++i){    
    Ai[i]=fabs(Av*AlAv[i])+RandN(1.0,sigma[i]);
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i]+5.0*log10(Ds*100.0)+Ai[i];
    if(Nblend[i]>=k)  s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i]));
    if(k==1){s.ext[i]=Ai[i];  s.Map[i]=Map[i];    s.Mab[i]= Mab[i];}
    if(Ai[i]<0.0  or Ai[i]>100.0 or Map[i]<Mab[i] or Ds<0.0  or Ds>20.0  or Map[i]<0.0 or s.Fluxb[i]<0.0){
    cout<<"ERROR filter:  "<<i<<"\t extinction:  "<<Ai[i]<<"\t App_mag:  "<<Map[i]<<"\t Abso_mag:  "<<Mab[i]<<endl;
    int rre; cin>>rre;}}




   Lstar[0]=Lstar[1]=Lstar[2]=Lstar[3]=Lstar[4]=0.0; 
   for(int i=1; i<Ntem; ++i){
   if(k==1 and  ( ((s.Tstar0-s.Tem[i])*(s.Tstar0-s.Tem[i-1])<0.0) or s.Tstar0==s.Tem[i]) ){
   s.Lstar[0]= pow(10.0,-0.4*s.Mm[0][i]);
   s.Lstar[1]= pow(10.0,-0.4*s.Mm[1][i]);
   s.Lstar[2]= pow(10.0,-0.4*s.Mm[2][i]);
   s.Lstar[3]= pow(10.0,-0.4*s.Mm[3][i]);
   s.Lstar[4]= pow(10.0,-0.4*s.Mm[4][i]);}
   if( ((Tstar-s.Tem[i])*(Tstar-s.Tem[i-1])<0.0) or Tstar==s.Tem[i]){
   Lstar[0]= pow(10.0,-0.4*s.Mm[0][i]);
   Lstar[1]= pow(10.0,-0.4*s.Mm[1][i]);
   Lstar[2]= pow(10.0,-0.4*s.Mm[2][i]);
   Lstar[3]= pow(10.0,-0.4*s.Mm[3][i]);
   Lstar[4]= pow(10.0,-0.4*s.Mm[4][i]);}}
   
   for(int i=0; i<M; ++i){
   Lstar[i]= Lstar[i]/s.Lstar[i]*pow(Rstar/s.Rstar0,2.0);
   if(Nblend[i]>=k and k!=1)     s.Fluxb2[i]+=fabs(Lstar[i]); }
   if(s.Lstar[0]==0.0 or s.Lstar[1]==0.0 or s.Lstar[2]==0.0 or s.Lstar[3]==0.0 or s.Lstar[4]==0.0 or
   s.Tstar0<400.0 or s.Tstar0>17000.0 or Tstar<400.0  or Tstar>17000.0){  
   cout<<"ERROR Tstar:  "<<s.Tstar0<<"\t Lstar0:  "<<s.Lstar[0]<<"\t Lstar1:  "<<s.Lstar[1]<<"\t Lstar2:  "<<s.Lstar[2]<<endl;
   cout<<"Tstar:  "<<Tstar<<"\t Rstar:  "<<Rstar<<endl;
   int yye ; cin>> yye;} }///loop 


    for(int i=0; i<M; ++i){
    s.magb[i] = -2.5*log10(s.Fluxb[i]);
    s.col=s.magb[2]-s.magb[4];//V-I
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    //s.blend2[i]=double(1.0/(1.0 + s.Fluxb2[i])); 
    if((Nblend[i]==1.0 && s.blend[i]<1.0) or s.blend[i]>1.0  or s.blend[i]<0.0 ){ 
    cout<<"BIGG ERRROR nsbl: "<<s.nsbl[i]<<"\t Nlend: "<<Nblend[i]<<"\t s.blend  "<<s.blend[i]<<endl; int uue; cin>>uue;}}
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_lens                               ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_lens(source & s , lens & l)
{
    double rholens[Num],test, f;
    l.rhomaxl=0.0;
      
    for(int k=1;k<=s.nums;++k){
    rholens[k]=0.0;
    l.Dl=k*step;
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds) {cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  int yye; cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}



    do{
    l.numl = (int)((double)rand()*1000./((double)(RAND_MAX+1.)*1000.)*(s.nums-2.0)+1.0);
    test = ((double)rand()/(double)(RAND_MAX+1.)*l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;
     int ue; cin>>ue;}
    }while(test>rholens[l.numl]);



   double  randflag=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[l.numl];
       if (randflag<=s.rho_disk[l.numl]) l.struc=0;///thin disk
   else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1; // bulge structure
   else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2; //thick disk
   else if (randflag<= s.Rostar0[l.numl]) l.struc=3;//halo
   else {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}



  if(l.struc==0){///thin disk
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(4.5-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*57.0);
  if(l.Ml<=1.0) f=pow(l.Ml,-1.6);
  if(l.Ml>1.0) f=pow(l.Ml,-3.0);
  }while(test>f);}


  if(l.struc==1){///Galactic bulge
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(1.4-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*378.0)+0.3;
  f=pow(l.Ml,-2.35);
  }while(test>f);}


  if(l.struc==2){///thick disk
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(1.4-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*3.0)+0.8;
  f=pow(l.Ml,-0.5);
  }while(test>f);}


  if(l.struc==3){///stellar halo
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(0.8-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*1.0)+0.8;
  f=pow(l.Ml,-0.5);
  }while(test>f);}





    l.Dl=l.numl*step;///kpc
    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*Msun*s.Ds*Pc*1000.0)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s,l);
    l.tE=l.RE/fabs(l.Vt*1000.0*3600.0*24.0)/1.6;///in day

    
    l.proj=double(l.xls*Rsun/l.RE);
    s.ro_star0=fabs(s.Rstar0*l.proj);
    l.ksi=randR(0.0,360.0)*pi/180.0;///[radian
    l.u0= randR(0.0,0.2);///*s.ro_star0;///*s.ro_star0;
    l.tstar= fabs(l.tE*s.ro_star0);///days
    if(l.Ml<0.0  or l.Dl<0.0  or l.Dl>s.Ds or s.ro_star0<0.0  or l.numl>s.nums  or l.Dl>20.0  or l.tE<0.0 or l.Vt<0.0 or l.xls>1.0){ 
    cout<<"Ml:  "<<l.Ml<<"\t Dl:  "<<l.Dl<<"\t tE:  "<<l.tE<<"\t u0:  "<<l.u0<<endl;
    cout<<"numl:  "<<l.numl<<"\t strucl: "<<l.struc<<"\t Ds:  "<<s.Ds<<"\t l.Vt:  "<<l.Vt<<endl;
    int uue;  cin>>uue; }    
   
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_Pulsing                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_puls(lens & l , pulsing & p){

  p.peri= randR(3.5,7.0);
  p.phit= randR(0.0,359.0)*pi/180.0;
  p.phir= p.phit-pi/2.0;
  p.epsiT= randR(100.0,550.0); 
  p.epsiR= randR(0.1,0.45); 
  p.omega=2.0*pi/p.peri; 
 
   double w=randR(0.0,6.0);    
    
   if(w<1.0)      {p.l=0;  p.m=0; }
   else if(w<2.0) {p.l=1;  p.m=0; }
   else if(w<3.0) {p.l=1;  p.m=1; }
   else if(w<4.0) {p.l=2;  p.m=0; }
   else if(w<5.0) {p.l=2;  p.m=1; }
   else if(w<6.0) {p.l=2;  p.m=2; }
   else if(w<7.0) {p.l=3;  p.m=0; }
   else if(w<8.0) {p.l=3;  p.m=1; }
   else if(w<9.0) {p.l=3;  p.m=2; }
   else if(w<10.0){p.l=3;  p.m=3; }
   else if(w<11.0){p.l=4;  p.m=0; }
   else if(w<12.0){p.l=4;  p.m=1; }
   else if(w<13.0){p.l=4;  p.m=2; }
   else if(w<14.0){p.l=4;  p.m=3; }
   else if(w<15.0){p.l=4;  p.m=4; }
   else if(w<16.0){p.l=5;  p.m=0; }
   else if(w<17.0){p.l=5;  p.m=1; }
   else if(w<18.0){p.l=5;  p.m=2; }
   else if(w<19.0){p.l=5;  p.m=3; }
   else if(w<20.0){p.l=5;  p.m=4; }    
   else if(w<21.0){p.l=5;  p.m=5; }    
 
   p.inc=randR(5.0,85.5)*pi/180.0;    
   l.t0=0.0;//randR(-0.5,0.5)*p.peri; 
  
   l.pt1=-0.9*l.tE;//-30.0*l.tstar;// + l.t0;///days
   l.pt2= 0.9*l.tE;//+30.0*l.tstar;// + l.t0;///days
 
   if(p.m>p.l or p.m<0.0 or double(p.inc*RA)>89.9){
   cout<<"ERRPR l: "<<p.l<<"\t m:  "<<p.m<<"\t inc:  "<<p.inc*RA<<"\t dt:  "<<l.dt2<<endl;  int ue; cin>>ue; }
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
double TETA(double y, double x){
    double tan=-1;
    if(y==0 and x>0.0)        tan=0.0;
    else if(y==0.0 and x<0.0) tan=(pi); 
    else if(x==0.0 and y>0.0) tan=(pi/2.0);
    else if(x==0.0 and y<0.0) tan=(3.0*pi/2.0);  
    else if(x>0.0 and y>0.0)  tan=(atan(y/x));
    else if(x<0.0 and y>0.0)  tan=(pi-atan(y/fabs(x)));
    else if(x<0.0 and y<0.0)  tan=(pi+atan(fabs(y/x)));
    else if(x>0.0 and y<0.0)  tan=(2.0*pi-atan(fabs(y/x))); 
    if(tan<0.0 or tan>2.0*pi){cout<<"ERROR (TETA): y:  "<<y<<"\t x:  "<<x<<"\t tan:  "<<tan<<endl;  exit(0);} 
    return(tan);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double RandN(double N, double sigma){
    double p, fp, frand;
    do{
    p=(double)rand()/((double)(RAND_MAX)+(double)(1))*N*sigma; ///[-N:N]
    fp=exp(-p*p/(2.*sigma*sigma));
    frand= (double)rand()/((double)(RAND_MAX)+(double)(1));
    }while(frand>fp);
    double sign= (double)rand()/((double)(RAND_MAX)+(double)(1));
    if(sign<0.5)     return(p);
    else             return(-p);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double randR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double Pfunc(int l, int m, double tets){
     double func, si, co; 
     co=cos(tets);  si=sin(tets);
    if(l==0)  func=1.0/sqrt(4.0);
    else if(l==1){
    if(abs(m)==1)   func=si*sqrt(3.0/8.0); 
    if(abs(m)==0)   func=co*sqrt(3.0/4.0); }
    else if(l==2){
    if(abs(m)==2)   func=si*si*sqrt(15.0/32.0); 
    if(abs(m)==1)   func=si*co*sqrt(15.0/8.0) ;
    if(abs(m)==0)   func=(3.0*co*co-1.0)*sqrt(5.0/16.0); }
    else if(l==3){
    if(abs(m)==3)   func=si*si*si*sqrt(35.0/64.0);
    if(abs(m)==2)   func=si*si*co*sqrt(105.0/32.0); 
    if(abs(m)==1)   func=si*(5.0*co*co-1.0)*sqrt(21.0/64.0);
    if(abs(m)==0)   func=co*(5.0*co*co-3.0)*sqrt(7.0/16.0); }
    else if(l==4){
    if(abs(m)==4)   func=si*si*si*si*sqrt(35.0/2.0)*3.0/16.0;
    if(abs(m)==3)   func=si*si*si*co*sqrt(35.0)*3.0/8.0; 
    if(abs(m)==2)   func=si*si*(7.0*co*co-1.0)*sqrt(45.0/2.0)/8.0;
    if(abs(m)==1)   func=si*co*(7.0*co*co-3.0)*sqrt(45.0)/8.0; 
    if(abs(m)==0)   func=(co*co*(35.0*co*co-30.0)+3.0)*sqrt(9.0)/16.0;}
    else if(l==5){
    if(abs(m)==5)   func=si*si*si*si*si*sqrt(77.0)*3.0/32.0;
    if(abs(m)==4)   func=si*si*si*si*co*sqrt(385.0/2.0)*3.0/16.0; 
    if(abs(m)==3)   func=si*si*si*(9.0*co*co-1.0)*sqrt(385.0)/32.0;
    if(abs(m)==2)   func=si*si*co*(3.0*co*co-1.0)*sqrt(1155.0/2.0)/8.0; 
    if(abs(m)==1)   func=si*(co*co*(21.0*co*co-14.0)+1.0)*sqrt(165.0/2.0)/16.0; 
    if(abs(m)==0)   func=co*(63.0*co*co*co*co-70.0*co*co+15.0)*sqrt(11.0)/16.0;}
    ///https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#%7F'%22%60UNIQ--postMath-0000000F-QINU%60%22'%7F_=_5[1]
    return(func/sqrt(pi));
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
if (l.Dl==0.0) l.Dl=0.00034735;

 double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*l.Dl*cos(s.TET)*cos(s.FI));
 double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*s.Ds*cos(s.TET)*cos(s.FI));
 if(Rlc==0.0) Rlc=0.00034346123;
 if(Rsc==0.0) Rsc=0.0004762654134;  
 ///Source and Lens velocity components in Galactocentric cylindrical coordinates
 double SVT, SVR, SVZ, LVT, LVR, LVZ,SVt,LVt;
 ///Source and Lens velocity components in heliocenteric galactic coordinates
 double SVb, SVl, LVb, LVl;
 double fv, testfv;
 double VSunl,VSunt,VSunb,vls_b,vls_l;
 double betal,betas,deltal,deltas,deltao;



 double NN=3.0;
 double VSunR =-10.3;
 double VSunT =vro_sun*(1.00762+0.00712)+6.3;
 double VSunZ = 5.9;
 double sigma_R_Disk=43.0, sigma_T_Disk=27.8, sigma_Z_Disk=17.5;
 double sigma_R_TDisk= 67.0, sigma_T_TDisk= 51.0, sigma_Z_TDisk= 42.0;
 double sigma_R_halo= 131.0, sigma_T_halo= 106.0, sigma_Z_halo= 85.0;
 double sigma_R_Bulge = 113.0,sigma_T_Bulge = 115.0,sigma_Z_Bulge = 100.0;
 double Rho[8]={00.0}; double maxr=0.0;
 for(int i=0; i<8; ++i){ 
 Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}




  double test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])     {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0;}
else if(test<=(Rho[0]+Rho[1])) {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]))  {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]))  
                           {sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))  
                           {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5]))
                           {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4;}
else if(test<=maxr)        {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5;}
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}  


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
/// Generate Source velocity components in Glactocenteric cylindrical coordinates(x',y')
    if(s.struc==0){///Galactic disk
    do{
    SVR = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_Disk;
    fv = exp(-1./2.*SVR*SVR/(sigma_R_Disk*sigma_R_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.); }while(testfv>fv);
    do{
    SVT = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_Disk;
    fv = exp(-1./2.*SVT*SVT/(sigma_T_Disk*sigma_T_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    do{
    SVZ = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_Disk;
    fv = exp(-1./2.*SVZ*SVZ/(sigma_Z_Disk*sigma_Z_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    SVT =SVT +vro_sun*(1.00762 * pow(Rsc/R_sun,0.0394) + 0.00712);}

    else if(s.struc==1){///Galactic bulge
    do{
    SVZ= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_Bulge;
    fv = exp(-1./2.*SVZ*SVZ/(sigma_Z_Bulge*sigma_Z_Bulge));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
    do{
    SVR= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_Bulge;
    fv = exp(-1./2.*SVR*SVR/(sigma_R_Bulge*sigma_R_Bulge));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
    do{
    SVT= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_Bulge;
    fv = exp(-1./2.*SVT*SVT/(sigma_T_Bulge*sigma_T_Bulge));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);}

    else if(s.struc==2) {///thick disk
    do{
    SVR = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_TDisk;
    fv = exp(-1./2.*SVR*SVR/(sigma_R_TDisk*sigma_R_TDisk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    do{
    SVT = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_TDisk;
    fv = exp(-1./2.*SVT*SVT/(sigma_T_TDisk*sigma_T_TDisk));
    testfv = (double)rand()/(double)(RAND_MAX+1.); }while(testfv>fv);
    do{
    SVZ = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_TDisk;
    fv = exp(-1./2.*SVZ*SVZ/(sigma_Z_TDisk*sigma_Z_TDisk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    SVT =SVT+ vro_sun *(1.00762*pow(Rsc/R_sun,0.0394) + 0.00712); }

    else if(s.struc==3){///stellar halo
    do{
    SVZ= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_halo;
    fv = exp(-1./2.*SVZ*SVZ/(sigma_Z_halo*sigma_Z_halo));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
    do{
    SVR= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_halo;
    fv = exp(-1./2.*SVR*SVR/(sigma_R_halo*sigma_R_halo));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
    do{
    SVT= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_halo;
    fv = exp(-1./2.*SVT*SVT/(sigma_T_halo*sigma_T_halo));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);}

    l.vs=sqrt(SVR*SVR+SVT*SVT+SVZ*SVZ);


///======================================================================================
/// Generate Lens velocity components in Glactocenteric cylindrical coordinates(x',y')
if(l.struc==0){///Galactic disk
    do{
    LVR = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_Disk;
    fv = exp(-1./2.*LVR*LVR/(sigma_R_Disk*sigma_R_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    do{
    LVT = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_Disk;
    fv = exp(-1./2.*LVT*LVT/(sigma_T_Disk*sigma_T_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    do{
    LVZ =(-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_Disk;
    fv = exp(-1./2.*LVZ*LVZ/(sigma_Z_Disk*sigma_Z_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    LVT =LVT+ vro_sun *(1.00762 * pow(Rlc/R_sun,0.0394) + 0.00712);}

   else if(l.struc==1){///Galactic bulge
   do{
   LVZ= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_Bulge;
   fv = exp(-1./2.*LVZ*LVZ/(sigma_Z_Bulge*sigma_Z_Bulge));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
   do{
   LVR= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_Bulge;
   fv = exp(-1./2.*LVR*LVR/(sigma_R_Bulge*sigma_R_Bulge));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
   do{
   LVT= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_Bulge;
   fv = exp(-1./2.*LVT*LVT/(sigma_T_Bulge*sigma_T_Bulge));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);}

   else if(l.struc==2){///thick disk
   do{
   LVR = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_TDisk;
   fv = exp(-1./2.*LVR*LVR/(sigma_R_TDisk*sigma_R_TDisk));
   testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
   do{
   LVT = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_TDisk;
   fv = exp(-1./2.*LVT*LVT/(sigma_T_TDisk*sigma_T_TDisk));
   testfv = (double)rand()/(double)(RAND_MAX+1.); }while(testfv>fv);
   do{
   LVZ = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_TDisk;
   fv = exp(-1./2.*LVZ*LVZ/(sigma_Z_TDisk*sigma_Z_TDisk));
   testfv = (double)rand()/(double)(RAND_MAX+1.); }while(testfv>fv);
   LVT =LVT+ vro_sun *(1.00762*pow(Rlc/R_sun,0.0394) + 0.00712); }

   else if(l.struc==3){///stellar halo
   do{
   LVZ= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_halo;
   fv = exp(-1./2.*LVZ*LVZ/(sigma_Z_halo*sigma_Z_halo));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
   do{
   LVR= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_halo;
   fv = exp(-1./2.*LVR*LVR/(sigma_R_halo*sigma_R_halo));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
   do{
   LVT= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_halo;
   fv = exp(-1./2.*LVT*LVT/(sigma_T_halo*sigma_T_halo));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);}
   l.vl=sqrt(LVT*LVT+LVZ*LVZ+LVR*LVR);
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

   if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc-1.0)<0.01) betal=pi/2.0; 
   else if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc+1.0)<0.01) betal=-pi/2.0; 
   else  betal=asin(l.Dl*cos(s.FI)*sin(s.TET)/Rlc);///lens[-pi/2,pi/2]
   if(fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc-1.0)<0.01) betas=pi/2.0; 
   else if(fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc+1.0)<0.01) betas=-pi/2.0; 
   else  betas=asin(s.Ds*cos(s.FI)*sin(s.TET)/Rsc);///lens[-pi/2,pi/2]
    if(fabs(l.Dl*cos(s.FI)*sin(s.TET))>Rlc || fabs(s.Ds*cos(s.FI)*sin(s.TET))>Rsc || Rlc==0.0 || Rsc==0.0){
    cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
    cout<<"FI: "<<s.FI<<"\t TET: "<<s.TET<<endl; 
    cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<endl;
    cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(s.TET)/Rlc<<"\t sin(s): "<<s.Ds*cos(s.FI)*sin(s.TET)/Rsc<<endl;
     //int ew; cin>>ew;
      }
      
       //betao=0.0; ///observer
    if(fabs(l.Dl*cos(s.FI))>sqrt(pow(R_sun,2.0)+pow(l.Dl*cos(s.FI)*sin(s.TET),2.0)) ) betal= pi-betal;
    if(fabs(s.Ds*cos(s.FI))>sqrt(pow(R_sun,2.0)+pow(s.Ds*cos(s.FI)*sin(s.TET),2.0)) ) betas= pi-betas;



    if(fabs((Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI))-1.0)<0.01)   deltal=0.0; 
    else if (fabs((Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI))+1.0)<0.01) deltal=pi; 
    else    deltal=acos((Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI)));
    if(fabs((Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI))-1.0)<0.01)   deltas=0.0; 
    else if (fabs((Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI))+1.0)<0.01) deltas=pi; 
    else    deltas=acos((Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI)));
   if(fabs((Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI)))>1.002 || fabs((Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI)))>1.002  || l.Dl==0.0 || s.Ds==0.0 || fabs(s.FI)==pi/2.0){
    cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
    cout<<"betal: "<<betal<<"\t betas: "<<betas<<endl; 
    cout<<"FI: "<<s.FI<<"\t TET: "<<s.TET<<endl; 
    cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<endl;
    cout<<"cos(dl): "<<(Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI))<<"\t cos(ds): "<<(Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI))<<endl;
    // int ew; cin>>ew; 
     }


    deltao=pi/2.0;
    SVl=-SVR*    sin(deltas)+ SVT* cos(deltas);
    LVl=-LVR*    sin(deltal)+ LVT* cos(deltal);
    VSunl=-VSunR*sin(deltao)+VSunT*cos(deltao);

    SVt=  1.0*SVR*cos(deltas)+  SVT*sin(deltas);
    LVt=  1.0*LVR*cos(deltal)+  LVT*sin(deltal);
    VSunt=1.0*VSunR*cos(deltao)+VSunT*sin(deltao);

    SVb=-sin(s.FI)*(SVt) + cos(s.FI)* (SVZ);
    LVb=-sin(s.FI)*(LVt) + cos(s.FI)* (LVZ);
    VSunb=-sin(s.FI)*(VSunt)+cos(s.FI)*(VSunZ);

    vls_l= LVl-l.xls*SVl -(1.0-l.xls)*VSunl;
    vls_b= LVb-l.xls*SVb -(1.0-l.xls)*VSunb;
    l.Vt=sqrt(fabs(vls_l*vls_l+ vls_b*vls_b));

if (l.Vt<0.0 || l.Vt>1.0e6 ){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;   int yee; cin>>yee;}
//cout<<"Vt: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;
}
///==================================================================
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nn=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   double fd=1.0; ///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;///No limitation 
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;
/*
I assume that the up-limit of the mass is indicated by the simulation. Because it depends on the structre, .... all details can be found in mass_averaged.cpp  code. */

   /*  
    char filename[40];
    FILE *fill;
    sprintf(filename,"./files/density/%c%.2lf%c%.2lf.dat",'D',s.lat,'_',s.lon);
    fill=fopen(filename,"w");
    if(!fill){cout<<"cannot open file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl;  exit(0);}
   */



for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;

   x=i*step;
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = R_sun-x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///M_sun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nn)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*exp(nn)*exp(-fabs(zb)/0.8)/(1.0+0.5*nn);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(0.5/R_sun,-2.44);
   else            s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(rdi/R_sun,-2.44);///M_sun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///M_sun/pc^3
///=================================================
///در اینجا اینکه تعداد ستاره ه
///ا به این بستگی دارد که ما  نمودار قدر رنگ مطلق ستاره ها را چگونه درست کرده باشیم.اگر هیچ گونه محدودیتی برای
///درست کردن آن  در نظر نگرفته ایم،  پس تعداد کل ستاره ها را نظر  میگیریم.
/// ولی بهتر است که ما رابطه بین قدر مطلق و جرم را تعیین کنیم. در این صورت می توانیم  ورودی قدر رنگ
///وارد شده به کد را خودمان محدود به ستاره های روشن کنیم تا سرعت اجرای برنامه بالارود.
///averaged mass are the same as the previous work!!! because we did not change the besancon model


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]);///[M_sun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[M_sun/deg^2]
s.Nstari[i]=binary_fraction*(s.rho_disk[i]*fd/0.403445+s.rho_ThD[i]*fh/0.4542+s.rho_halo[i]*fh/0.4542+s.rho_bulge[i]*fb/0.308571);////[Nt/pc^3] 

s.Nstari[i]= s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]

s.Nstart  +=s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
//fprintf(fill,"%e   %e   %e   %e   %e  %e   %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]);
   }
 // fclose(fill);
}
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex){
  double F=-1.0;
  if(ds<ex.dis[0])        F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{ 
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] && ds<ex.dis[i+1]){
  F = ex.Extks[i]+(double)(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0||F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; int uut; cin>>uut;}
  return(F);
}
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     if(s.lon==0.0 or Lon<0.0 or Lon<0.1  or s.lon<0.1) Lon=360.00;

     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     if(Lon>360.000 || Lon<0.25 || fabs(Lat)>10.0 || (Lon>100 && Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}


     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");
     cout<<"Lat:  "<<Lat<<"\t Lon:    "<<Lon<<endl;

     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<64881; ++i) {
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 && fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     //cout<<"distance: "<<ex.dis[i]<<"\t Extks: "<<ex.Extks[i]<<endl;
     if(ex.dis[i]<0.2  || ex.dis[i]>50.0 || ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; } }     
     fclose(fpd);}
     cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<flag<<endl;
     return(flag);
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm)
{
    int nrr[4]={0}; 
    int yye, uui,h, k1, k2, g; double metal,age, mk;   
    char filename[40];
    FILE *fp2;

    int number[70]={0};   int count[70]={0};   double Metal[70]={0.0}; 
    FILE *meta; 
    meta=fopen("./files/CMD_BESANCON/metal.txt","r"); 
    for(int i=0; i<70; ++i){
    fscanf(meta,"%lf %d  %d\n",&Metal[i],&count[i],&number[i]);    
    if((Metal[i]<Metal[i-1] and i>0) or float(Metal[i])<-0.001 or number[i]==0 or count[i]>YZ or 
       (abs(count[i-1]+number[i-1]-count[i])>2 and i>0)){
    cout<<"ERROR Metal[i]: "<<Metal[i]<<"\t count[i]: "<<count[i]<<"\t number[i]: "<<number[i]<<"\t i: "<<i<<endl; cin>>uui;} }
    fclose(meta);
    double Age[YZ]={0.0}; double B1[YZ]={0.0};  double M1[YZ]={0.0};   double mm[YZ]={0.0}; 
    FILE *ji; 
    ji=fopen("./files/CMD_BESANCON/RVI.txt", "r"); 
    for(int i=0; i<YZ; ++i){
    fscanf(ji,"%lf   %lf   %lf  %lf\n",&Age[i],&mm[i],&B1[i],&M1[i]); 
    if(Age[i]<0.0 or mm[i]<0.0 or fabs(B1[i])>1.7 or M1[i]<0.5 or Age[i]>18.0){   
    cout<<"ERROR Age(JI): "<<Age[i]<<"\t metal: "<<mm[i]<<"\t B[i]"<<B1[i]<<"\t M[i]: "<<M1[i]<<"\t i: "<<i<<endl;
    cin>>uui;}}
    fclose(ji); 

////=================================== THIN DISK ==============================

    int j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c.dat",'C','M','D','T','i');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTi.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.Teff_d[j],&age,&cm.logl_d[j],&cm.gra_d[j],&cm.metal_d[j],&cm.Rs_d[j],
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.Mab_d[4][j],&mk,&cm.cl_d[j],&cm.type_d[j]);
///*******************************************************
    metal= cm.metal_d[j];
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0])         h=0; 
    else if(metal>Metal[69] or  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm[k1]!=mm[k2-1] or number[h]==0 or fabs(metal-Metal[h])>0.05){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  ||  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_d[3][j]= double(B1[g]+M1[g]*(cm.Mab_d[2][j]+cm.Mab_d[4][j])*0.5); ///R-band versus (V+I)
    if(fabs(cm.Mab_d[4][j]-cm.Mab_d[3][j])>2.5 or fabs(age-Age[g])>3.0){
    cout<<"ERROR: Mab_d(z-band): "<<cm.Mab_d[4][j]<<"\t Mab[3]:  "<<cm.Mab_d[3][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
///*******************************************************   
    if(cm.mass_d[j]<0.0 or cm.mass_d[j]==0.0 or cm.gra_d[j]>6.0 or cm.Teff_d[j]<0.0 or cm.metal_d[j]>0.12 or age>10.0 or age<0.0 or cm.cl_d[j]>5 or cm.type_d[j]>=9.0 or cm.type_d[j]<2.0 or (cm.cl_d[j]==5 and int(cm.type_d[j])>7) or (cm.cl_d[j]==6 and int(cm.type_d[j])!=9) or 
    (cm.cl_d[j]<5 and int(cm.type_d[j])==9)){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; 
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<endl;
    cin>>yye;}
    j++;  } fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;




////=================================== Galactic Bulge ===========================================================
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c.dat",'C','M','D','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDb.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.Teff_b[j],&age,&cm.logl_b[j],&cm.gra_b[j],&cm.metal_b[j],&cm.Rs_b[j],
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.Mab_b[4][j],&mk,&cm.cl_b[j],&cm.type_b[j]);
///*******************************************************
    metal= cm.metal_b[j];
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0])         h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm[k1]!=mm[k2-1] or number[h]==0 or fabs(metal-Metal[h])>0.4){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  ||  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_b[3][j]= double(B1[g]+M1[g]*(cm.Mab_b[2][j]+cm.Mab_b[4][j])*0.5); ///R-band versus (V+I)
    if(fabs(cm.Mab_b[4][j]-cm.Mab_b[3][j])>2.5 or fabs(age-Age[g])>3.0){
    cout<<"ERROR:  Mab_b(z-band): "<<cm.Mab_b[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
///*******************************************************   
    if(cm.mass_b[j]<0.0 or cm.mass_b[j]==0.0 or cm.Teff_b[j]<0.0 or age>10.0 or cm.metal_b[j]>0.2 or cm.cl_b[j]>5 or cm.type_b[j]>8.0 or (cm.cl_b[j]==5 and int(cm.type_b[j])>7) or (cm.cl_b[j]==6 and int(cm.type_b[j])!=9) or (cm.cl_b[j]<5 and int(cm.type_b[j])==9)){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<endl; cin>>yye;}
    //if(cm.cl_b[j]>2  and  cm.cl_b[j]<5){cm.rr_b[j]=1;  nrr[1]+=1;}
    //else cm.rr_b[j]=0;  
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;










////==================================== THICK DISK ==========================================================
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c.dat",'C','M','D','T','k');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTk.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.Teff_t[j],&age,&cm.logl_t[j],&cm.gra_t[j],&cm.metal_t[j],&cm.Rs_t[j],
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.Mab_t[4][j],&mk,&cm.cl_t[j],&cm.type_t[j]);
///*******************************************************
    metal= cm.metal_t[j];
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0])         h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm[k1]!=mm[k2-1] or number[h]==0 or fabs(metal-Metal[h])>0.4){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  ||  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_t[3][j]= double(B1[g]+M1[g]*(cm.Mab_t[2][j]+cm.Mab_t[4][j])*0.5); ///R-band versus (V+I)
    if(fabs(cm.Mab_t[4][j]-cm.Mab_t[3][j])>2.5 or fabs(age-Age[g])>3.0){
    cout<<"ERROR: Mab_t(z-band): "<<cm.Mab_t[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
///*******************************************************   
    if(cm.mass_t[j]<0.0||  cm.mass_t[j]==0.0 or cm.Teff_t[j]<0.0 or cm.metal_t[j]>0.06 || cm.cl_t[j]>5 || cm.type_t[j]>8.0 or 
    (cm.cl_t[j]==5 and float(cm.type_t[j])>8.0) or (cm.cl_t[j]==6 and int(cm.type_t[j])!=9) or (cm.cl_t[j]<5 and int(cm.type_t[j])==9)){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"mass: "<<cm.mass_t[j]<<"\t TefF:  "<<cm.Teff_t[j]<<"\t metal:  "<<cm.metal_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;






////=================================== STELLAR HALO ============================================================ 
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c.dat",'C','M','D','h');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDh.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&age,&cm.logl_h[j],&cm.gra_h[j],&cm.metal_h[j],&cm.Rs_h[j],
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.Mab_h[4][j],&mk,&cm.cl_h[j],&cm.type_h[j]);
///*******************************************************
    metal= cm.metal_h[j];
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0])         h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm[k1]!=mm[k2-1] or number[h]==0 or fabs(metal-Metal[h])>0.4){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  ||  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_h[3][j]= double(B1[g]+M1[g]*(cm.Mab_h[2][j]+cm.Mab_h[4][j])*0.5); ///R-band versus (V+I)
    if(fabs(cm.Mab_h[4][j]-cm.Mab_h[3][j])>2.5 or fabs(age-Age[g])>3.0 ){
    cout<<"ERROR: Mab_h(z-band): "<<cm.Mab_h[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
///*******************************************************       
    if(cm.mass_h[j]<0.0 || cm.mass_h[j]==0.0 || age<13.0 or age>15.0 or cm.cl_h[j]<0  or cm.cl_h[j]>5  or  cm.Teff_h[j]<0.0 ||
    cm.metal_h[j]>0.04 || cm.cl_h[j]>8 || cm.type_h[j]>9 or (cm.cl_h[j]==5 and int(cm.type_h[j])>7) or (cm.cl_h[j]==6 and int(cm.type_h[j])!=9) or   (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
    cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
    cout<<"nrr[0]:  "<<nrr[0]<<"\t nrr[1]:   "<<nrr[1]<<"\t nrr[2]:   "<<nrr[2]<<"\t nrr[3]:    "<<nrr[3]<<endl;
    cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

