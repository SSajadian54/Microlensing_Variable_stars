#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include "VBBinaryLensingLibrary.h"
using namespace std;

///**************** constants *********************//


const int Num=2000;
const double MaxD=20.0;///kpc
const double RA=180.0/M_PI;
const double step=MaxD/(double)Num/1.0;///step in kpc
const int Ntem=int(11979);  
const int Nperi=100; 

const int M=4;
const int nf=int(5);//number of filter 
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
const double hplanck= 6.62607004*pow(10.0,-34.0);//## m2 kg / s
const double kbol=1.380649*pow(10.0,-23); 
const double BB= double(8.0*pi*hplanck*velocity*velocity*pow(10.0,30.0));
const double AA= double(hplanck*velocity*1000000.0/kbol); 

const double wave[nf]={0.35,0.5,0.65,0.8,1.0};

///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
const double R_sun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]=  {3.1,2.5,3.1,3.1};
const int N1=27004,N2=36558, N3=2612, N4=492;///CMD_BESANCON, thinD, bulge, thickD, halo
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
struct pulsing{
   double peri;///period
   double phir, phit;
   double epsiR, epsiT; ///amplitude
   double omega; 
   double inc;  
   int l,m;
};
struct lens{
    double RE, tE, Vt, Dl, Ml;
    double t0, ksi, u0, pt1, pt2,tstar;
    double proj,dt;
    double vl, vs, xls;
    double rhomaxl;
    int numl, struc;
    double dis, q;

};
struct source{
    double lon, lat; 
    double mass,Ds; 
    double Lstar[nf], Tstar0, Rstar0, ro_star0;
    double logl,logt,logg;
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart,nstarti;
    double FI, TET;
    int nums, struc, cl; 
    double mI,map;  
};
struct CMD{
    double logt_d[N1],logl_d[N1],Mab_d[M][N1],mass_d[N1], type_d[N1], gra_d[N1]; int cl_d[N1];  ///thin disk
    double logt_b[N2],logl_b[N2],Mab_b[M][N2],mass_b[N2], type_b[N2], gra_b[N2]; int cl_b[N2];  /// bulge
    double logt_t[N3],logl_t[N3],Mab_t[M][N3],mass_t[N3], type_t[N3], gra_t[N3]; int cl_t[N3];  ///thick disk
    double logt_h[N4],logl_h[N4],Mab_h[M][N4],mass_h[N4], type_h[N4], gra_h[N4]; int cl_h[N4];  /// halo
};
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_lens(source &s,lens & l);
void Func_source(source & s, CMD & cm);
void Func_puls(lens & l, pulsing & p);
double randN(double N, double sigma);
double randR(double down, double up);
double TETA(double y, double x);  
double Pfunc(int l, int m, double tets);///Legandre function
void read_cmd(CMD & cm);
void vrel(source & s,lens & l);
void Disk_model(source & s);
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

   pulsing p;
   lens l;
   source s; 
   CMD cm; 
   read_cmd(cm);
   
   
  
   printf("START_TIME >>>>>>>> %s",ctime(&_timeNow));
   VBBinaryLensing *vbb=new VBBinaryLensing;
   vbb->Tol=1.e-4;
   vbb->LoadESPLTable("./files/ESPL.tbl");
   char filnam1[40]; char filnam2[40]; 
   FILE* param;   FILE* magm;    FILE* sour;     
   FILE * tems;
   tems= fopen("./files/FstarT_Bessel.dat","r"); 
   double Tem[Ntem], M[nf][Ntem];
   for(int i=0; i<Ntem; ++i){
   fscanf(tems,"%lf   %lf   %lf   %lf   %lf   %lf\n",&Tem[i],&M[0][i],&M[1][i],&M[2][i],&M[3][i],&M[4][i]);}
   fclose(tems); 
   s.lat=-1.5;
   double lonn=0.5;
   cout<<"latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
   if(lonn<=0.0)   s.lon=360.0+lonn;
   else            s.lon=lonn;
   s.TET=(360.0-s.lon)/RA;///radian
   s.FI=s.lat/RA;///radian
   Disk_model(s);
   
   
   
   double t, xlens, ylens, tets, phis, tets0 , phis0, Legan, ux, uy, u,magni,Tstar, area;
   double delA, del[nf], A[nf], B[nf], A0, A1, Lstar[nf];
   double Rstar, xs, ys, zs, xo, yo, zo, rhoe, dis, dmin;
   double dtet=0.55, dphi0=0.55, dphi;
   int NN=0;
   for(tets=dtet; tets<180.0; tets+= dtet){
   for(phis=0.0;  phis<360.0; phis+=double(dphi0/sin(tets/RA)))  NN+=1;}
   NN=int(NN/2.+500.0); 
   cout<<"NN:  "<<NN<<endl;
   double Yo[NN], Zo[NN], tet[NN], phi[NN]; 
   double dzmax, dymax, dmax, ddel; 
   int nt,j1,Nlens, flag[2],yye; 
   double ymin, ymax, zmin, zmax, tp, ratio_a;
   double aveL[nf], aveLR, aveR0, Area0, conl,L[nf]; 
  
 
 
   param=fopen("./files/param_all_b1.txt","a+");
   fclose(param); 
   for(int icon=1; icon<100; ++icon){    
   
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   


   do{
   Func_source(s,cm);
   Func_lens(s,l);
   Func_puls(l,p);
   Nlens= int(l.pt2-l.pt1)/l.dt + 2; 
   }while(l.tE<=1.0 or l.tE>300.0 or Nlens<80 or Nlens>5000 or l.dt<=0.0 or s.map>21.0 or s.map<12.5); 
   cout<<"Nlens:  "<<Nlens<<"\t tE: "<<l.tE<<"\t tstar:  "<<l.tstar<<"\t p.peri:  "<<p.peri<<endl;
   
   
  
   
   int count[Nlens]={0}; int j=0; 
   for(int i=0; i<Nlens; ++i){
   if(j>=Nperi)  j=j-Nperi;
   count[i]=j;   j+=1;}
   //double M0[2][Nlens]={0.0}; 
   double Msh[2][Nlens]={0.0};  double Mto[2][nf][Nlens]={0.0};
   for(int i=0; i<Nlens; ++i){
   //M0[0][i]=   M0[1][i]=0.0; 
   Msh[0][i]= Msh[1][i]=0.0; 
   for(int j=0; j<nf; ++j){
   Mto[0][j][i]=Mto[1][j][i]=0.0;}}
   conl=0.0; 
   for(int i=0;  i<Ntem;  ++i){
   if(( (s.Tstar0-Tem[i])*(s.Tstar0-Tem[i-1])<0.0  and  i>0) or s.Tstar0==Tem[i]){
   for(int j=0; j<nf; ++j)  s.Lstar[j]= pow(10.0,-0.4*M[j][i]);
   break; }}
   if(s.Lstar[0]==0.0  or s.Lstar[1]==0.0  or s.Lstar[2]==0.0 or s.Lstar[3]==0.0 or s.Lstar[4]==0.0 or s.Tstar0<21.0  or s.Tstar0>12000){  
   cout<<"ERROR Tstar: "<<s.Tstar0<<"\t Lstar0: "<<s.Lstar[0]<<"\t Lstar1:  "<<s.Lstar[1]<<"\t Lstar2:  "<<s.Lstar[2]<<endl;
   cin>>yye; }
   conl= double(s.Lstar[1]+ s.Lstar[2] + s.Lstar[3] + s.Lstar[4])/4.0 ; 
   for(int k=0; k<nf; ++k) L[k]=double(s.Lstar[k]*Area0/conl); 
   
   
   
   sprintf(filnam1,"./files/Bins/%c%c%c%c%d.dat",'M','a','g','_',icon);
   magm=fopen(filnam1,"w");
   sprintf(filnam2,"./files/Bins/%c%c%c%c%d.dat",'S','o','u','_',icon);
   sour=fopen(filnam2,"w");
   fclose(sour);
   cout<<"Files were made"<<endl;
  
  
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   aveL[0]=aveL[1]=aveL[2]=aveL[3]=aveL[4]=aveLR=aveR0=0.0;  
   Area0=pi*s.Rstar0*s.Rstar0;
   for(int Npr=0;  Npr<Nperi;  Npr +=1){////time of pulsation   
   tp=double(Npr*l.dt); 
   
   
   if(Npr==0){
   phis=0.0;  
   tets=pi/2.0-p.inc-dtet/RA;//the points with ys=0
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
   ddel=fabs(double(dmax*0.6) );
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
   if(Tstar<Tem[0] or Tstar==Tem[0])  j1=0;   
   else if(Tstar>Tem[Ntem-1])         j1=int(Ntem-1); 
   else{
   for(int k=1; k<Ntem; ++k){
   if( (Tstar-Tem[k])*(Tstar-Tem[k-1])<0.0 or Tstar==Tem[k]){ j1=k;   break;}}}
   for(int k=0; k<nf; ++k){
   Lstar[k]= pow(10.0,-0.4*M[k][j1]);}} 
  
  
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
   tets0>pi or tets0<0.0 or phis0<0.0 or phis0>(2.0*pi) )){
   cout<<"ERROR***   yo: "<<yo<<"\t zo:  "<<zo<<"\t xo:  "<<xo<<endl; 
   cout<<"Rstar:  "<<Rstar<<"\t Rstar0: "<<s.Rstar0<<endl; 
   cout<<"Tstar: "<<Tstar<<"\t Tstar0:  "<<s.Tstar0<<endl;
   cout<<"tets0:  "<<tets0*RA<<"\t phis0:  "<<phis0*RA<<endl;
   cout<<"tets:  "<<tets*RA<<"\t phis:  "<<phis*RA<<endl; cin>>yye;}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH    
   if( flag[1]>0 or flag[0]>0 ){  
   if( Npr==0 ){ 
   sour=fopen(filnam2,"a+");
   fprintf(sour,"%d  %d   %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf\n",
   flag[0],flag[1],yo,zo,Rstar,Tstar,tets*RA,tets0*RA,phis*RA,phis0*RA);
   fclose(sour);}


  
  
  if(flag[1]>0){
  aveLR += area;
  for(int k=0; k<nf; ++k){
  aveL[k] += double(Lstar[k]/conl)*area;} }
  if(flag[0]>0)  aveR0 += area; 
  
   
   for(int nl=0;  nl<Nlens;  ++nl){
   if(int(count[nl]) == Npr){
   t=double(l.pt1 + nl*l.dt);
   xlens = (t-l.t0)/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi);
   ylens = (t-l.t0)/l.tE * sin(l.ksi) + l.u0 * cos(l.ksi);
   ux=xlens - yo*l.proj;
   uy=ylens - zo*l.proj;
   u=sqrt(ux*ux+ uy*uy); 
   //magni=vbb->ESPLMag2(u,rhoe);
     
   magni=vbb->BinaryMag2(l.dis,l.q, ux, uy, rhoe );
   
   //if(flag[0]>0){
   //M0[1][nl]+= area*magni; 
   //M0[0][nl]+= area*1.000;}
   if(flag[1]>0){
   Msh[1][nl]+= area*magni; 
   Msh[0][nl]+= area*1.000;
   for(int k=0; k<nf; ++k){
   Mto[1][k][nl] += double(Lstar[k]/conl)*magni*area;
   Mto[0][k][nl] += double(Lstar[k]/conl)*1.000*area;}}}}} 
   }} 
   cout<<"tp:  "<<tp<<"\t Npr:   "<<Npr<<"\t Nperiod:  "<<Nperi<<endl;
   }///time loop
    
   for(int k=0; k<nf; ++k){ aveL[k]= double(aveL[k]/Nperi/1.00);}
   aveLR=double(aveLR/Nperi/1.0);
   aveR0=double(aveR0/Nperi/1.0);
    
    
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH    
   for(int nl=0; nl<Nlens; ++nl){
   t=double( l.pt1 + nl*l.dt);
   xlens = (t-l.t0)/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi);
   ylens = (t-l.t0)/l.tE * sin(l.ksi) + l.u0 * cos(l.ksi);
   
   u=sqrt(xlens*xlens+ ylens*ylens); 
//   A0=vbb->ESPLMag2(u,s.ro_star0);
   
   A0=vbb->BinaryMag2(l.dis,l.q, xlens, ylens, s.ro_star0 );
   
   
   //A0=double( M0[1][nl])/aveR0;//magnification factor for a non-palsating star
   A1=double(Msh[1][nl])/aveLR;//magnification factor for a palsating star with change in R
   delA=double(A1-A0)/A0; 
   for(int k=0;  k<nf;  ++k){
   A[k]=double(Mto[1][k][nl])/aveL[k] ;////double(M0[0][nl]);
   B[k]=double(Mto[0][k][nl])/aveL[k];///  /double(M0[0][nl]);///total luminosity 
   del[k]=(A[k]-A0)/A0;}
   ratio_a=double(Msh[0][nl]-aveLR)*10.0/double(aveLR);
   tp=double(count[nl]*l.dt)/p.peri;
   fprintf(magm,"%.4lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.4lf  %.4lf  %.4lf %.4lf  %.7lf  %.7lf  %.7lf %.7lf  %.7lf %.4lf %.4lf %.7lf %.5lf\n",(t-l.t0)/l.tstar,A0,A1,A[1],A[2],A[3],A[4],B[1],B[2],B[3],B[4],del[1],del[2],del[3],del[4],delA,xlens,ylens,ratio_a,tp);//20
   if(nl%10==0){
   cout<<"================================================"<<endl;
   cout<<"\t **** counter:  "<<nl<<endl;
   cout<<"(t-t0)/tstar:  "<<(t-l.t0)/l.tstar<<endl;
   cout<<"xlens: "<<xlens<<"\t ylens: "<<ylens<<endl;
   cout<<"A0:  "<<A0<<"\t A1:  "<<A1<<"\t A[1]: "<<A[1]<<endl;
   cout<<"A[2]"<<A[2]<<"\t A[3]"<<A[3]<<"\r A[4]:  "<<A[4]<<endl;
   cout<<"lumi[2]"<<B[2]<<"\t lumi[3]"<<B[3]<<"\t lumi[4]:  "<<B[4]<<endl;
   cout<<"delA:  "<<delA<<"\t del[1]: "<<del[1]<<"del[2]:  "<<del[2]<<endl;
   cout<<"del[3]:  "<<del[3]<<"\t del[4]:  "<<del[4]<<endl; 
   cout<<"ratio_area[%]:  "<<ratio_a<<endl;
   cout<<"================================================"<<endl;} }



   param=fopen("./files/param_all_b1.txt","a+");   
   fprintf(param,"%d  %d  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.7lf   %.5lf  %.5lf  %d  %d  %d  %.5lf  %.5lf  %.5lf    %.5lf %.5lf  %.5lf  %.8lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %d  %d  %.5lf %.5lf %.5lf %.5lf %.5lf  %d  %.5lf  %.5lf \n",
   l.numl,l.struc,l.Ml,l.Dl,l.xls,l.RE/Au,l.Vt,l.vs,l.vl,l.tE,l.u0/s.ro_star0,l.ksi,l.tstar,   
   s.nums,s.struc,s.cl,s.mass,s.Ds,s.logg,s.logl,s.Tstar0,s.Rstar0,s.ro_star0,s.Lstar[1],s.Lstar[2],s.Lstar[3],s.Lstar[4],s.mI,s.map, 
   p.peri,p.phir*RA,p.phit*RA,p.epsiR,p.epsiT,p.inc*RA,p.l,p.m,aveLR,aveL[1],aveL[2],aveL[3],aveL[4],icon, l.dis, l.q);//43
   fclose(param); 
   
    cout<<"******************************************************"<<endl;
    cout<<"********* icon:  "<<icon<<endl; 
    cout<<"l:  "<<p.l<<"\t m:  "<<p.m<<"inclination:  "<<p.inc*RA<<endl;
    cout<<"LENS:   mass[Msun]: "<<l.Ml<<"\t lens_dis[kpc]:  "<<l.Dl<<"\t RE[AU]: "<<l.RE/Au<<endl;
    cout<<"xls: "<<l.xls<<"\t u0:  "<<l.u0<<"\t ksi:  "<<l.ksi<<endl;
    cout<<"tstar:  "<<l.tstar<<"\t pt1:  "<<l.pt1<<"\t pt2:   "<<l.pt2<<endl;
    cout<<"tE[days]: "<<l.tE<<"\t Vt:  "<<l.Vt<<"\t Vl:  "<<l.vl<<endl;
    cout<<"source:  mass:  "<<s.mass<<"Ds[Kpc]: "<<s.Ds<<endl;
    cout<<"mI:  "<<s.mI<<"\t map: "<<s.map<<endl;
    cout<<"Rstar[Rsun]: "<<s.Rstar0<<"\t Tstar[K]:  "<<s.Tstar0<<endl;
    cout<<"ro_star: "<<s.ro_star0<<"\t u0: "<<l.u0/s.ro_star0<<endl;
    cout<<"logg:  "<<s.logg<<"\t logl:  "<<s.logl<<"\t logt:  "<<s.logt<<endl;
    cout<<"PULS: peri: "<<p.peri<<"\t phir: "<<p.phir*180.0/pi<<endl;
    cout<<"epsiR/Rstar: "<<p.epsiR<<"\t epsiT:  "<<p.epsiT<<endl;
    cout<<">>>>>:  Area0:  "<<Area0<<"\t <area0>_P:   "<<aveR0<<"\t <area1>_P: "<<aveLR<<endl;
    cout<<">>aveL[1]:  "<<aveL[1]<<"\t aveL[2]:  "<<aveL[2]<<"\t aveL[3]:  "<<aveL[3]<<"\t aveL[4]:  "<<aveL[4]<<endl; 
    cout<<">>L(R0,T0)[1]:  "<<L[1]<<"\tL(R0,T0)[2]:  "<<L[2]<<"\tL(R0,T0)[3]:  "<<L[3]<<"\tL(R0,T0)[4]:  "<<L[4]<<endl;   
    cout<<"***************************************************** "<<endl;     
    fclose(magm);}
    
   // fclose(param);  
    time(&_timeNow);
    printf("END_TIME >>>>>>>>  %s ",ctime(&_timeNow));
    delete vbb;
    return(0);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_source                             ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_source(source & s, CMD & cm)
{
    double rho, rf, Av;   
    int num, yye; 
    do{
    s.nums=int(randR(1.0,Num-1));
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    }while(rho>s.Rostari[s.nums] || s.nums<5);///distance larger than 50.0
    s.Ds=(double)(s.nums*step);
    if(s.Ds>MaxD or s.Ds<0.1){cout<<"ERROR (1): Ds: "<<s.Ds<<"\t MaxD: "<<MaxD<<"\t step: "<<step<<"\t num: "<<s.nums<<endl; cin>>yye;}


    
    rf=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[s.nums];
         if (rf<= s.rho_disk[s.nums])                    s.struc=0;///thin disk
    else if (rf<=(s.rho_disk[s.nums]+s.rho_bulge[s.nums])) s.struc=1;/// bulge structure
    else if (rf<=(s.rho_disk[s.nums]+s.rho_bulge[s.nums]+s.rho_ThD[s.nums])) s.struc=2;///thick disk
    else if (rf<=s.Rostar0[s.nums]) s.struc=3;///halo




    if(s.struc==0){
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-1.0));
    s.logt=cm.logt_d[num]; 
    s.logl=cm.logl_d[num];
    s.cl=    cm.cl_d[num];
    s.logg=cm.gra_d[num];
    s.mass=cm.mass_d[num];
    s.mI=cm.Mab_d[3][num]; }


    if(s.struc==1){///bulge
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-1.0));
    s.logt=cm.logt_b[num];
    s.logl=cm.logl_b[num];
    s.cl=    cm.cl_b[num];
    s.logg=cm.gra_b[num]; 
    s.mass=cm.mass_b[num];
    s.mI=cm.Mab_b[3][num]; }



    if(s.struc==2){///thick disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-1.0));
    s.logt=cm.logt_t[num];
    s.logl=cm.logl_t[num];
    s.cl=    cm.cl_t[num];
    s.logg= cm.gra_t[num]; 
    s.mass=cm.mass_t[num];
    s.mI=cm.Mab_t[3][num]; }

    if(s.struc==3){///stellar halo
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-1.0));
    s.logt=cm.logt_h[num];
    s.logl=cm.logl_h[num];
    s.cl   = cm.cl_h[num];
    s.logg=cm.gra_h[num];  
    s.mass=cm.mass_h[num];
    s.mI=cm.Mab_h[3][num]; }
    s.Tstar0= pow(10.0,s.logt);
    s.Rstar0= sqrt(s.mass*27400.0/pow(10.0,s.logg));     
    Av= 4.381e-06*pow(s.Ds,6.0)-0.0002941*pow(s.Ds,5.0)+0.007594*pow(s.Ds,4.0)-0.09202*pow(s.Ds,3.0)+0.4775*pow(s.Ds,2.0)-0.2454*s.Ds+0.6476;
    s.map= s.mI + 5.0*log10(s.Ds*100.0) + fabs(Av*0.6); 
    
    if(Av>6.0 or Av<0.0 or s.Ds>20.0  or s.Ds<0.0 or s.mass<0.0 or s.cl>9 or s.Tstar0<0.0 or s.nums>Num or s.Rstar0<0.0 or s.Rstar0>1000.0){
    cout<<"Ds;  "<<s.Ds<<"\t nums:  "<<s.nums<<endl; 
    cout<<"struc:  "<<s.struc<<endl;
    cout<<"ERROR Ds:  "<<s.Ds<<" \t Av:  "<<Av<<endl; 
    cout<<"Source_mass:  "<<s.mass<<"\t logt:  "<<s.logt<<"\t logl:  "<<s.logl<<endl;
    cout<<"cl:  "<<s.cl<<"\t logg:  "<<s.logg<<"\t Tstar0:  "<<s.Tstar0<<"\t Rstar0:  "<<s.Rstar0<<endl;
    int yyw;  cin>>yyw;}  
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


l.dis=randR(0.6,1.6); 
l.q=randR(0.1,1.0); 


    l.Dl=l.numl*step;///kpc
    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*(1.0+l.q)*Msun*s.Ds*Pc*1000.0)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s,l);
    l.tE=l.RE/fabs(l.Vt*1000.0*3600.0*24.0);///in day
 //l.tE=24.0;    
    
    l.proj=double(l.xls*Rsun/l.RE);
    s.ro_star0=fabs(s.Rstar0*l.proj);
    l.ksi=randR(0.0,30.0)*pi/180.0;///[radian
    l.u0=randR(0.0,0.5);//*s.ro_star0;///*s.ro_star0;
    l.tstar= fabs(l.tE*s.ro_star0);///days
    if(l.Ml<0.0  or l.Dl<0.0  or l.Dl>s.Ds or s.ro_star0<0.0  or l.numl>s.nums  or l.Dl>20.0  or l.tE<0.0 or l.Vt<0.0 or l.xls>1.0){ 
    cout<<"Ml:  "<<l.Ml<<"\t Dl:  "<<l.Dl<<"\t tE:  "<<l.tE<<"\t u0:  "<<l.u0<<endl;
    cout<<"numl:  "<<l.numl<<"\t strucl: "<<l.struc<<"\t Ds:  "<<s.Ds<<"\t l.Vt:  "<<l.Vt<<endl;
    int uue;  cin>>uue; }    
    
    cout<<"dis:  "<<l.dis<<"\t q:  "<<l.q<<endl;
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_Pulsing                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_puls(lens & l , pulsing & p){
  p.peri=1.8;/// randR(0.5,7.0);
  p.phit=pi/2.0;///randR(0.0,359.0)*pi/180.0;
  p.phir=0.0;///p.phit-pi/2.0;
  p.epsiT=350.0;//randR(250.0,600.0); 
  p.epsiR=0.25; //randR(0.1,0.4); //
  p.omega=2.0*pi/p.peri; 
  //p.l= int(randR(0.0,5.5)); 
  //p.m= int(randR(0.0,p.l+0.01));
 
   double w=randR(0.0,21.0);

w=4.5;    
     
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
 
 
 
 
  p.inc=randR(1.0,89.9)*pi/180.0;    
  l.t0=0.0;//randR(-0.5,0.5)*p.peri; 
  
 // l.pt1=-40.0*l.tstar;// + l.t0;///days
 // l.pt2=+40.0*l.tstar;// + l.t0;///days
  
  l.pt1=-0.9*l.tE;// l.tstar;///days
  l.pt2=+0.9*l.tE;//l.tstar;///days
  
  l.dt=double(p.peri/Nperi);    
  if(p.m>p.l or p.m<0.0 or double(p.inc*RA)>89.9  or l.dt==0.0){
  cout<<"ERRPR l: "<<p.l<<"\t m:  "<<p.m<<"\t inc:  "<<p.inc*RA<<"\t dt:  "<<l.dt<<endl;  int ue; cin>>ue; }
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
double randN(double N, double sigma){
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
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm){
////mass log10(T) log10(age) log10(L) log10(g) metal U G R I Z Bj Vj Rj Ij CL TYP
    int yye; int age;  
    double Mu,Mg,Mr,Mi,Mz,metal; 
    char filename[40];
    FILE *fp2;
///ابتدا سعی کردیم به کمک روابط درون مقاله Fernie 1983 اینها را به فیلتر های کرون کویزن تبدیل کنیم و لی بعد ما متوجه شدیم که اینها همان فیلتر های استاندارد ///هستند.  

////=================================== THIN DISK ==============================

    int j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c%c.dat",'C','M','D','T','i','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTi.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %d %lf\n",
    &cm.mass_d[j],&cm.logt_d[j],&age,&cm.logl_d[j],&cm.gra_d[j],&metal,&Mu,&Mg,&Mr,&Mi,&Mz,
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.Mab_d[3][j],&cm.cl_d[j],&cm.type_d[j]);
   if(cm.mass_d[j]<0.0 or cm.mass_d[j]==0.0 or cm.logt_d[j]<0.0 or Mr>20.0 or metal>0.12 or age>10 or cm.cl_d[j]>7 or cm.type_d[j]>9 or 
    cm.type_d[j]<2.0 or (cm.cl_d[j]==5 and int(cm.type_d[j])>7) or (cm.cl_d[j]==6 and int(cm.type_d[j])!=9) or 
    (cm.cl_d[j]<5 and int(cm.type_d[j])==9)){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; 
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<endl;
    cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;


////=================================== BULGE ================================== 
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c%d.dat",'C','M','D','B','b',2);
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDB.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %d %lf \n",
    &cm.mass_b[j],&cm.logt_b[j],&age,&cm.logl_b[j],&cm.gra_b[j],&metal,&Mu,&Mg,&Mr,&Mi,&Mz,
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.Mab_b[3][j],&cm.cl_b[j],&cm.type_b[j]);
    if(cm.mass_b[j]<0.0|| cm.mass_b[j]==0.0 ||cm.logt_b[j]<0.0||Mr>18.0||age>10 || metal>0.9||cm.cl_b[j]>7||cm.type_b[j]>9 or 
    (cm.cl_b[j]==5 and int(cm.type_b[j])>7) or (cm.cl_b[j]==6 and int(cm.type_b[j])!=9) or (cm.cl_b[j]<5 and int(cm.type_b[j])==9)){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<endl;
    cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;



////=================================== THICK DISK =============================
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c%c.dat",'C','M','D','T','k','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTk.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %d %lf \n",
    &cm.mass_t[j],&cm.logt_t[j],&age,&cm.logl_t[j],&cm.gra_t[j],&metal,&Mu,&Mg,&Mr,&Mi,&Mz,
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.Mab_t[3][j],&cm.cl_t[j],&cm.type_t[j]);
    if(cm.mass_t[j]<0.0||  cm.mass_t[j]==0.0 or cm.logt_t[j]<0.0 or Mr>20.0|| metal>0.025||cm.cl_t[j]>7|| cm.type_t[j]>9 or 
    (cm.cl_t[j]==5 and int(cm.type_t[j])>7) or (cm.cl_t[j]==6 and int(cm.type_t[j])!=9) or (cm.cl_t[j]<5 and int(cm.type_t[j])==9)){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;



////=================================== STELLAR HALO =========================== 
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON/%c%c%c%c%c.dat",'C','M','D','H','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDH.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf %d  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.logt_h[j],&age,&cm.logl_h[j],&cm.gra_h[j],&metal,&Mu,&Mg,&Mr,&Mi,&Mz,
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.Mab_h[3][j],&cm.cl_h[j],&cm.type_h[j]);
    if(cm.mass_h[j]<0.0 || cm.mass_h[j]==0.0 || age<0 or cm.cl_h[j]<0  or cm.cl_h[j]>7  or  cm.logt_h[j]<0.0 || Mr>20.0 ||
    metal>0.01 || cm.cl_h[j]>7|| cm.type_h[j]>9 or (cm.cl_h[j]==5 and int(cm.type_h[j])>7) or (cm.cl_h[j]==6 and int(cm.type_h[j])!=9) or
    (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
    cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
    cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
   
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
