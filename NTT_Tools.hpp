#include <algorithm>

#ifndef NTT_MOD
#define NTT_MOD 998244353
#endif

#ifndef NTT_MAXN
#define NTT_MAXN 800005
#endif

// define start
#define NTT_MUL
#define NTT_INV
#define NTT_SQRT
#define NTT_LN
#define NTT_EXP
// define end

namespace NTT_Tools{
    using namespace std;

    const int MOD=NTT_MOD,INV2=(MOD+1)/2;
    const int MAXN=NTT_MAXN;
    int c[MAXN];
    int ksm(int p,int q){
        int r=1;
        while(q){
            if(q&1) r=1ll*r*p%MOD;
            p=1ll*p*p%MOD,q>>=1;
        }
        return r;
    }
    void NTT(int *a,int N,bool invFlag){
        int n=__lg(N);
        c[0]=0;
        for(int i=1;i<N;i++){
            int lbit=i&-i;
            c[i]=c[i^lbit]|(1<<(n-1-__lg(lbit)));
            if(i<c[i]) swap(a[i],a[c[i]]);
        }
        for(int k=1;k<N;k<<=1){
            int w=ksm(3,(MOD-1)/(k<<1));
            if(invFlag) w=ksm(w,MOD-2);
            for(int i=0;i<N;i+=k<<1){
                int I=1;
                for(int j=0;j<k;j++){
                    a[i|j|k]=(a[i|j]-1ll*I*a[i|j|k]%MOD+MOD)%MOD;
                    a[i|j]=(2ll*a[i|j]-a[i|j|k]+MOD)%MOD;
                    if(invFlag) a[i|j]=1ll*a[i|j]*INV2%MOD,a[i|j|k]=1ll*a[i|j|k]*INV2%MOD;
                    I=1ll*I*w%MOD;
                }
            }
        }
    }
    
    #ifdef NTT_MUL
    int multmp[MAXN],multmp2[MAXN];
    void mul(int *a,int *b,int *res,int N1,int N2,int NN){
        int N=1;
        while(N<N1+N2-1) N<<=1;
        for(int i=0;i<N;i++) multmp[i]=a[i],multmp2[i]=b[i];
        for(int i=N1;i<N;i++) multmp[i]=0;
        for(int i=N2;i<N;i++) multmp2[i]=0;
        NTT(multmp,N,false),NTT(multmp2,N,false);
        for(int i=0;i<N;i++) res[i]=1ll*multmp[i]*multmp2[i]%MOD;
        NTT(res,N,true);
        for(int i=N;i<NN;i++) res[i]=0;
        for(int i=NN;i<N;i++) res[i]=0;
    }

    #ifdef NTT_INV
    int invtmp[MAXN],invtmp2[MAXN];
    void getInv(int *a,int *res,int N){
        int N1=1,NN1;
        res[0]=ksm(a[0],MOD-2);
        for(int i=1;i<N;i++) res[i]=0;
        while(N1<N){
            NN1=N1<<1;
            mul(res,res,invtmp,N1,N1,NN1);
            mul(a,invtmp,invtmp2,NN1,NN1,NN1);
            for(int i=0;i<NN1;i++) res[i]=(2ll*res[i]-invtmp2[i]+MOD)%MOD;
            N1=NN1;
        }
    }

    #ifdef NTT_SQRT
    int sqrttmp[MAXN],sqrttmp2[MAXN];
    void getSqrt(int *a,int *res,int N){
        int N1=1,NN1;
        res[0]=1;
        while(N1<N){
            NN1=N1<<1;
            mul(res,res,sqrttmp,N1,N1,NN1);
            for(int i=0;i<NN1;i++) sqrttmp[i]=(sqrttmp[i]+a[i])%MOD,res[i]=2*res[i]%MOD;
            getInv(res,sqrttmp2,NN1);
            mul(sqrttmp,sqrttmp2,res,NN1,NN1,NN1);
            N1=NN1;
        }
    }
    #endif

    #ifdef NTT_LN
    bool invInited=false;
    int inv[MAXN],lntmp[MAXN],lntmp2[MAXN];
    void initInv(){
        inv[0]=inv[1]=1;
        for(int i=2;i<MAXN;i++) inv[i]=1ll*inv[MOD%i]*(MOD-MOD/i)%MOD;
        invInited=true;
    }
    void getLn(int *a,int *res,int N){
        getInv(a,lntmp,N);
        for(int i=0;i<N-1;i++) lntmp2[i]=1ll*(i+1)*a[i+1]%MOD;
        mul(lntmp,lntmp2,res,N,N,N);
        inv[1]=1;
        if(!invInited) initInv();
        for(int i=N-1;i>=1;i--) res[i]=1ll*res[i-1]*inv[i]%MOD;
        res[0]=0;
    }

    #ifdef NTT_EXP
    int exptmp[MAXN],exptmp2[MAXN];
    void getExp(int *a,int *res,int N){
        res[0]=1;
        int N1=1,NN1;
        while(N1<2*N){
            NN1=N1<<1;
            getLn(res,exptmp,N1);
            for(int i=0;i<NN1;i++) exptmp[i]=(a[i]-exptmp[i]+MOD)%MOD,exptmp2[i]=res[i];
            exptmp[0]=(exptmp[0]+1)%MOD;
            mul(exptmp,exptmp2,res,NN1,N1,NN1);
            N1=NN1;
        }
    }
    #endif
    #endif
    #endif
    #endif
}