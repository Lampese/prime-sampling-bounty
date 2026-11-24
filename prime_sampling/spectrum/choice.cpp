#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
struct Wheel {int w,phi; std::vector<std::pair<int,double>> terms;};
static Wheel make_wheel(int w,int phi){
    Wheel out{w,phi,{{1,1.0}}};
    for(int q:{2,3,5,7,11})if(w%q==0){
        std::vector<std::pair<int,double>> next;
        for(auto [d,k]:out.terms){
            if(q!=2)next.push_back({d,k*(q-2)/(q-1)});
            next.push_back({d*q,k/(q-1)});
        }
        out.terms=next;
    }
    return out;
}
extern "C" void choose(double a,double b,double V,double s,double k,double span,
                       double pmin,double p,double md,double best,int* ids){
    static const std::vector<Wheel> wheels={
        make_wheel(1,1),make_wheel(2,1),make_wheel(3,2),make_wheel(6,2),
        make_wheel(15,8),make_wheel(30,8),make_wheel(210,48),make_wheel(2310,480)};
    static const std::array<double,8> sizes={1./64,1./16,1./4,1./2,1.,2.,4.,8.};
    double order[2][2]={{best,best},{best,best}};
    for(int j=0;j<4;++j)ids[j]=-1;
    const double L=std::log(a), sq=std::sqrt(p);
    for(int wi=0;wi<(int)wheels.size();++wi){
        const auto& wh=wheels[wi]; int w=wh.w, phi=wh.phi;
        double D=std::ceil(span/w)+3;
        for(int mode=0;mode<2;++mode){
            if(!mode && phi==1)continue;
            const auto& supports=mode?wheels[0].terms:wh.terms;
            int factor=mode?1:w;
            for(int l=0;l<8;++l){
                double M=std::max(1.,std::floor(D*sizes[l])),H=D+M, mass=M*(mode?phi:1.);
                if(mass<1e6||2*w*H>=pmin)continue;
                double J=H*(mode?1.:phi),K=factor*H;
                double den=std::max(1.,std::floor(std::sqrt(V*pmin/(2*b*K))));
                double size=b*(std::ceil(p/den)+den*K);
                double f=std::max(1.,std::floor(std::log(size)/L));
                double near=f*(1+std::ceil(V/den)+std::ceil(2*b*den*K/pmin));
                if(V<=b){
                    double weighted=0;
                    for(auto [d,coef]:supports){
                        double best_d=V+1;
                        for(int aux=1;aux<=2;++aux){
                            double d0=std::floor(std::sqrt(V*pmin*d/(2*b*K*aux)));
                            if(d0<1)continue;
                            double Z=b*(std::ceil(p/d0)+std::ceil((d0+1)*aux*K/d));
                            double freq=std::max(1.,std::floor(std::log(Z)/L));
                            double t=std::floor(std::sqrt(2*V*b*K/(pmin*d*aux)))+1;
                            best_d=std::min(best_d,(1+1/H)*(freq-1)+freq*(4./3+2./H)*t);
                        }
                        weighted+=coef*best_d;
                    }
                    near=std::min(near,weighted);
                }
                double rr=1/mass,co=(1+md*near/s)/(J*s*s),front=phi*H/k;
                double base=std::min(3.0001*p,p*(2*sq+1)*rr+4*sq+5+6*rr*std::pow(2*sq+6,2));
                for(int j=0;j<=1;++j){
                    double var=base+(j+1)*std::pow(p*rr,2)+2*p*rr*(1+2*rr);
                    double val=std::sqrt(front*(1/s+std::sqrt(rr+std::sqrt(co*var))));
                    int id=(wi*2+mode)*8+l;
                    if(val<order[j][0]){order[j][1]=order[j][0];ids[j*2+1]=ids[j*2];order[j][0]=val;ids[j*2]=id;}
                    else if(val<order[j][1]){order[j][1]=val;ids[j*2+1]=id;}
                }
            }
        }
    }
}
