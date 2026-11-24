#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>

int main() {
    int levels, primes, queries;
    int64_t unit;
    if (!(std::cin>>levels>>primes>>queries>>unit)) return 1;
    std::vector<int64_t> penalty(levels),upper(levels),lower(levels);
    for(int j=0;j<levels;++j){
        std::cin>>penalty[j];
        upper[j]=unit+penalty[j]; lower[j]=unit-penalty[j];
    }
    std::vector<int> query_prime(queries), query_col(queries);
    std::vector<std::vector<int>> at(primes+1);
    std::vector<int64_t> dl(queries),du(queries);
    for (int k=0;k<queries;++k) {
        std::cin >> query_prime[k] >> query_col[k];
        if(query_prime[k]<=0 || query_prime[k]>primes || query_col[k]<0 || query_col[k]>=levels) return 3;
        at[query_prime[k]].push_back(k);
    }
    for(int i=0; i<primes;++i){
        int p,step; std::cin>>p>>step;
        for(int j=levels-1;j>=0;--j){
            if(j<step){
                lower[j]=std::max<int64_t>(0,lower[j]-penalty[j]-(unit+p-1)/p);
            }else{
                const auto u=upper[j-step],l=lower[j-step];
                upper[j]-=l/p;
                lower[j]=std::max<int64_t>(0,lower[j]-(u+p-1)/p);
            }
        }
        for(int k : at[i+1]) { dl[k]=lower[query_col[k]]; du[k]=upper[query_col[k]]; }
    }
    for(int k=0;k<queries;++k) std::cout<<query_prime[k]<<" "<<query_col[k]<<" "<<dl[k]<<" "<<du[k]<<"\n";
}
