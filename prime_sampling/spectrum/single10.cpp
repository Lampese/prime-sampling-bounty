#include <algorithm>
#include <cmath>
#include <vector>
#include <utility>
struct Model{double w,phi,w1,p1,p2,divs;std::vector<std::pair<int,double>> terms;};
extern "C" void choose(double a,double b,double v,double s,double k,double span,double pm,double p,double md,double current,int* id){
 static const auto rows=[](){
  std::vector<Model> v;
  const int primes[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
  double w=1,phi=1;
  for(int d=1;d<=15;d++){
   w*=primes[d-1];phi*=primes[d-1]-1;
   double w1=1,p1=1;std::vector<std::pair<int,double>> t{{1,1}};
   for(int i=0;i<=std::min(d,5);i++){
    v.push_back({w,phi,w1,p1,phi/p1,std::pow(2.,d-i),t});
    if(i==std::min(d,5))break;
    const int z=primes[i];std::vector<std::pair<int,double>> n;
    for(auto [x,k]:t){if(z!=2)n.push_back({x,k*(z-2)/(z-1)});n.push_back({x*z,k/(z-1)});}t=n;w1*=z;p1*=z-1;
   }
  }
  return v;
 }();
 const double sq=std::sqrt(p);
 static const double ratios[]={1./128,1./32,1./8,1./3,.7,1.5,3.,6.,12.,24.};
 double best[2]={current,current};id[0]=id[1]=-1;
 for(int i=0;i<rows.size();++i){
  const auto&r=rows[i];if(r.w>=a)continue;
  double d0=std::ceil(span/r.w)+5;
  for(int j=0;j<10;j++){
   double m=std::max(1.,std::floor(d0*ratios[j])),h=d0+2*m,sz=r.p2*m;
   if(sz<1.e6||2*r.w*h>=pm)continue;
   double J=r.p1*h,K=r.w1*h;
   double near=0;
   for(auto [d,z]:r.terms){
    double bestd=v+1;
    for(int aux=1;aux<=2;aux++){
     double den=std::floor(std::sqrt(v*pm*d/(2*b*K*aux)));if(den<1)continue;
     double U=b*(std::ceil(p/den)+std::ceil((den+1)*aux*K/d));
     double f=std::max(1.,std::floor(std::log(U)/std::log(a)));
     double rt=std::floor(std::sqrt(2*v*b*K/(pm*d*aux)))+1;
     bestd=std::min(bestd,(1+1/h)*(f-1)+f*(4./3+2/h)*rt);
    }
    near+=z*bestd;
   }
   double rr=(2+1/(m*m))/(3*sz),co=(1+md*near/s)/(J*s*s),front=r.phi*h/k;
   double l1=std::min(1./std::sqrt(r.p2),r.divs*(1+std::sqrt(2*r.w/(r.w1*m)))/r.p2),c=16.*sq*std::pow(p/m,2)*rr*l1*l1;
   for(int q=0;q<2;q++){
    double var=c+(q+1)*std::pow(p*rr,2)+2*p*rr*(1+2*rr);
    double value=std::sqrt(front*(1/s+std::sqrt(rr+std::sqrt(co*var))));
    for(int r:{3,4,5}){
     double mult=(r==3?15.:r==4?105.:945.)*std::pow(q?2.:1.,r);
     double main=q?(r==3?15.:r==4?105.:945.):(r==3?6.:r==4?24.:120.);
     double gen=(std::pow(4.,r)-2)*p*sq*rr*std::pow(p*l1*l1/(m*m),r-1)+main*p*p*std::pow(rr,r)+(mult-main)*p*std::pow(rr,r-1);
     double nv=std::sqrt(front*(1/s+std::pow(co*gen,1./(2*r))));
     if(nv<value)value=nv;
    }
    if(value<best[q]){best[q]=value;id[q]=10*i+j;}
   }
  }
 }
}
