#include<NTL/ZZ.h>
#include <NTL/vector.h>
#include <hash_map>
NTL_CLIENT;
using namespace std;
struct sp_f{
    ZZ *q;
    ZZ *e;
};


ZZ gxx(ZZ g,ZZ x,ZZ p){
    int l;
	ZZ z;
    int str[50];
    if (x==0) 
		return to_ZZ(1);
    for(int i=0;i<50;i++) 
		str[i]=2;
    int i=0;
    while (x!=0){
        str[i]=x%2;
        x=x/2;
        i++;
        l=i;
    }
    z=g;
    for(i=l-2;i>=0;i--){
        z=(z*z)%p;
         if (str[i]==1) 
			z=(z*g)%p;
    }
	z=z%p;
    return z;
}

ZZ e_gcd(ZZ a,ZZ b,ZZ &x,ZZ &y)
{
    if (b == 0)
    {
        x = 1;
        y = 0;
        return a;
    }
    ZZ ans = e_gcd(b, a%b, x, y);
    ZZ temp = x;
    x = y;
	y = temp - (a/b)*y;
    return ans;
}

ZZ cal(ZZ a,ZZ m)
{
    ZZ x, y, ans;
    ZZ gcd = e_gcd(a, m, x, y);
	if (gcd != 1) 
		return to_ZZ("-1");
	
    ans = x % m;
    if (ans<=0)
		ans += m;
    return ans;
}

ZZ exinv(ZZ a,ZZ m)
{
    ZZ s;
    s = cal(a, m);
    if (s==-1) 
		return to_ZZ("0");
    else return s;
}
bool is_p(ZZ n){
	if(n==to_ZZ("5853296629055571413"))
		return true;
	
	else{
		if(n==to_ZZ("185378401392658780760176993633305263162028554127893025955620165648027002231849"))
			return false;
		for(ZZ x=to_ZZ(2);x<SqrRoot(n)+1;x++)
			if(n%x==0)
				return false;
	}
	return true;
}

sp_f p_f(ZZ n){ //素数分解
	struct sp_f ans;
	const ZZ s2=n;
	ZZ i,s=to_ZZ(0);
	int len=10000;
	ans.q=new ZZ[len];
	ans.e=new ZZ[len];
	ZZ s1=n;
	i=1;

	for(int k=0;k<len;k++){
		ans.e[k]=1;
	}
	for(int k=0;k<len;k++){
		ans.q[k]=1;
	}
	int j=0;
	while(1){
		if((s1/2)%2==0){
			s1=s1/2;
			ans.e[j]=ans.e[j]+1;
			ans.q[j]=2;
		}
		else{
			if(n%2==0){
				ans.q[j]=2;
				j++;
				s1=s1/2;
			}
			break;
		}
	}
	while(!is_p(s1)){
		while(1){
			if(ans.q[j]==1 ||ans.q[j]==n){
				ans.q[to_int(j)]=GCD(gxx(to_ZZ("2"),i,s1)-1,s1);
			}
			else
				break;
			i++;
			}
		while((s1/ans.q[j])%ans.q[j]==0){
				ans.e[j]=ans.e[j]+1;
				s1=s1/ans.q[j];	
			}
	/*	cout<<"第"<<j<<"个位的指数为 "<<ans.e[j]<<endl;
		cout<<ans.q[j]<<" ";
		cout<<ans.e[j]<<endl;*/
		s1=s1/ans.q[j];
		j++;
	}
	ans.q[j]=s1;
	return ans;
}
/*sp_f p_f(ZZ s0){
	struct sp_f a;
	const ZZ s2=s0;
	ZZ i,s=to_ZZ(0);
	int j=0;
	int len=10000;
	a.q=new ZZ[len];
	a.e=new ZZ[len];
	for(int k=0;k<len;k++){
		a.e[k]=1;
	}
    s=SqrRoot(s0)+1;
	ZZ s1=s0;
    for (i=2;i<s;i++) 
		if (s0%i==0){
			while((s1/i)%i==0){
				a.e[j]=a.e[j]+1;
				s1=s1/i;
			}
			s1=s0;
			a.q[j]=i;
			s0=s0/(gxx(i,a.e[j],s2));
			j++;
		}
		return a;
}*/
ZZ order(ZZ g,ZZ p){
	ZZ *q=p_f(p-1).q;
	ZZ *e=p_f(p-1).e;
	int i=0;
	ZZ ans;
	while(q[i]!=0){
		ZZ N=p-1;
		/*cout<<q[i]<<endl;
		cout<<e[i]<<endl<<endl;*/
		ans=gxx(g,N/q[i],p);
		if(ans==1){
			return N/q[i];
		}
		i++;
	}
	return p-1;
	
}
ZZ BSGS(ZZ g,ZZ h,ZZ p){
		ZZ N,n,x,gi;
		N=order(g,p);
		n= SqrRoot(N)+1;
		cout<<n<<" "<<to_long(n+1)<<endl;
		int length=2147483647;
		int i=0;
		int j=0;
		cout<<length<<endl;
		ZZ *a=new ZZ[length];
		for(i=0;i<n+1;i++) {
			a[i]=gxx(g,conv<ZZ>(i),p);
		}
		ZZ *b=new ZZ[to_int(n+1)];
		gi=exinv(a[to_int(n)],p);
		for(i=0;i<=n;i++){
			b[i]=(gxx(gi,conv<ZZ>(i),p)*h)%p;
		}
		for(i=0;i<=n;i++)
			for(j=0;j<n+1;j++)
				if (a[i]==b[j]){
					return to_ZZ(i)+to_ZZ(j)*n;
					break;
				}
}

ZZ Chinese_Remainder(ZZ r[],ZZ prime[],ZZ len)
{
    ZZ i;
	ZZ m=to_ZZ(1);
	ZZ sum=to_ZZ(0);
	ZZ d,x,y;
	ZZ n=to_ZZ(1);
    //计算所以除数的积n，也是所以除数的最小公倍数
    for(i=0;i<len;i++)
        n=n*prime[to_int(i)];
    //计算符合所以条件的数
    for(i=0;i<len;i++)
    {
        m=n/prime[to_int(i)];//计算除去本身的所有除数的积m
        d=e_gcd(prime[to_int(i)],m,x,y);//计算w[i]*x+m*y=gcd(w[i],m)的一个解y
        //累加整数解y的同并不断对n取余，其利用公式:(a+b)%c=(a%c+b%c)%c
        sum=(sum+y*m*r[to_int(i)])%n;
    }
    return (n+sum%n)%n;//满足所以方程的最小解
}
ZZ ppf(ZZ g,ZZ h,ZZ p){
	int i=0;
	ZZ *e=p_f(p-1).e;
	ZZ *prime=p_f(p-1).q;
	int len=10000;
	ZZ *qei=new ZZ[len];
	ZZ *x=new ZZ[len];
	while(prime[i]!=1){
		qei[i]=gxx(prime[i],e[i],p);
		cout<<"qe"<<i<<"="<<qei[i]<<endl;
		x[i]=0;
		i++;
	}
	for(;i<len;i++){
		qei[i]=0;
		x[i]=0;
	}
	i=0;
	while(qei[i]!=0){
		x[i]=BSGS(gxx(g,(p-1)/qei[i],p),gxx(h,(p-1)/qei[i],p),p); //计算g^((p-1)x/qiei)=h^((p-1)/qiei) solve for x
		cout<<gxx(g,(p-1)/qei[i],p)<<" "<<gxx(h,(p-1)/qei[i],p)<<endl;
		cout<<"x["<<i<<"]="<<x[i]<<endl<<endl;
		i++;
	}
	
	return Chinese_Remainder(x,qei,to_ZZ(i));
}



int main(){
//	cout<<BSGS(to_ZZ("4"),to_ZZ("11844727"),to_ZZ("41022299"))<<endl<<BSGS(to_ZZ("1558416"),to_ZZ("32438768"),to_ZZ("41022299"))<<endl;
//	cout<<ppf(to_ZZ("2"),to_ZZ("3415549615389079368"),to_ZZ("11706593258111142827"))<<endl;
	ZZ g,gi,h,p;
	sp_f ans=p_f(to_ZZ("37075680278531756152035398726661052632405a7108255786051911240331296054004463698"));
	for(int i=0;ans.q[i]!=1;i++){
		cout<<ans.q[i]<<" "<<ans.e[i]<<endl;
	};
	//cout<<endl<<gxx(to_ZZ("1054592380"),to_ZZ("1021763679"),to_ZZ("1889570071"))<<endl;
	//cout<<endl<<gxx(to_ZZ("1477875652"),to_ZZ("1519424709"),to_ZZ("1889570071"))<<endl;
	cout<<BSGS(to_ZZ("2"),to_ZZ("3415549615389079368"),to_ZZ("11706593258111142827"))<<endl;
//	cout<<gxx(qq,to_ZZ(1),to_ZZ("61063"))<<endl;
	cout<<gxx(to_ZZ("2"),to_ZZ("120"),to_ZZ("1739"))<<endl;
	cout<<GCD(gxx(to_ZZ("2"),to_ZZ("12"),to_ZZ("1739"))-1,to_ZZ("1739"));
	cin>>g>>h>>p;
	cout<<endl<<gxx(to_ZZ("1244183534"),to_ZZ("1671690479"),to_ZZ("1889570071"))<<endl;
	BSGS(g,h,p);
	
	
    system("pause");
    return 0;
}