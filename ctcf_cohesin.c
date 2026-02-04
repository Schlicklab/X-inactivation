/***************************************************************

METROPOLIS Monte Carlo Simulation for self-avoiding bead-spring
polymer in 3D

***************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<mpi.h>

#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))
#define Sqr(xx)     ((xx) * (xx))
#define Cube(xx)    ((xx) * (xx) * (xx))
#define AllocMem(a, n, t)  a = (t *) malloc ((n) * sizeof (t))

//Random number generator
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}

double Energy(double xo,double yo,double zo, int i);
void BindCohesin(int id);
void BuildNbr();
void UnbindCohesin(int id);
void SlideRightCohesin(int id);
void SlideLeftCohesin(int id);

/********************************************************************
VARIABLES:
N                       Number of beads
m                       mass
k                       spring constant
x,y,z                   Position of particle
xn,yn,zn                new position of particle
xo,yo,zo                old position of particle
dr                      random displacement factor
E,Eo,En                 current, old, new potential energy
ks                      spring constant
r0                      equilibrium bond length
rg                      radius of gyration

********************************************************************/

double *x,*y,*z;
int *ctcf;
int *fwd,*bwd;
int nfwd,nbwd;
int **cohesin;
int **nbr;
int *cnbr;
int *s;
double *pf,*pb;
double sigma,eps,ecut,rc,rc2,rca2,Rs;
double epsa,ecuta;
double epsb,ecutb;
double epsab,ecutab;
int N,ncoh,nsite,nctcf;
int maxnbr;
double J;
double ks,kc,kappa;
double r0;
//temp
double en_lj,en_sp;
int world_rank=0;								//rank of the process in MPI_COMM_WORLD
int world_size=0;								//size of MPI_COMM_WORLD
long idum;


int main()
{
MPI_Init(NULL,NULL);
char file[64];
MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
MPI_Comm_size(MPI_COMM_WORLD,&world_size);

double m=1.0;
int i,j,k,id,count=0;
int n,nstep,naccept=0;
int nf;
int sum;
int nn,t;
double xr,yr,zr,rr;
double xn,yn,zn;
double xo,yo,zo;
double *x0,*y0,*z0;
double cx,cy,cz;
double dr,d,dd;
double E,Eo,En;
double r_slide,r_bind;
double p_slide,p_bind,p_unbind;
double dE,dE_b;
double T;
double beta;
double eta,mu;
double r2i,r6i;
double rg;
double rnum;
double s0,A0;
double maxr;
double scale;
int nc0,ns0;
int last,diff;
FILE *ip,*ip1,*ip2,*op;

N=500;
maxnbr=500;
nstep=11000000;
nfwd=7;
nbwd=17;


nc0=5;
eps = 0.3;
epsab = 0.2;
srand(time(NULL));
epsa = 0.5;
epsb = 0.5;
mu=10.0;
r_bind=0.000005;
r_slide=0.0003;
sigma = 1.0;
dr = sigma/2.0;
//rc = pow(2.0,1.0/6.0)*sigma;
rc = 2.5*sigma;
rc2 = Sqr(rc);
Rs = 3.0;
r0 = sigma;
ks = 100.0;
T = 1.0;
beta  = 1.0/T;
//srand(12345);
idum=-955431-world_rank-rand()%100;
r2i=Sqr(sigma)/rc2;
r6i=r2i*r2i*r2i;
ecut=4.0*eps*r6i*(r6i-1.0);
ecuta=4.0*epsa*r6i*(r6i-1.0);
ecutab=4.0*epsab*r6i*(r6i-1.0);
ecutb=4.0*epsb*r6i*(r6i-1.0);

//Memory allocation
AllocMem(x,N,double);
AllocMem(y,N,double);
AllocMem(z,N,double);
AllocMem(x0,N,double);
AllocMem(y0,N,double);
AllocMem(z0,N,double);
AllocMem(ctcf,N,int);
AllocMem(s,N,int);
AllocMem(fwd,N,int);
AllocMem(bwd,N,int);
AllocMem(pf,N,double);
AllocMem(pb,N,double);
AllocMem(cnbr,N,int);
AllocMem(cohesin,N,int*);
  for(i=0;i<N;i++)
    AllocMem(cohesin[i],3,int);
AllocMem(nbr,N,int*);
  for(i=0;i<N;i++)
    AllocMem(nbr[i],maxnbr,int);

//Generate initial random configuration 
x[0]=0.0;
y[0]=0.0;
z[0]=0.0;
s[0]=0;
    for(i=1;i<N;i++)
    {
    s[i]=0;
    d=0.0;
        while(d<sigma)
        {
        xr=1.9*ran2(&idum)-1.0;
        yr=1.9*ran2(&idum)-1.0;
        zr=1.9*ran2(&idum)-1.0;
        dd=sqrt(Sqr(xr)+Sqr(yr)+Sqr(zr));
        xr=xr/dd;
        yr=yr/dd;
        zr=zr/dd;
        x[i]=x[i-1]+xr*r0;
        y[i]=y[i-1]+yr*r0;
        z[i]=z[i-1]+zr*r0;
            for(j=0;j<i;j++)
            {
                if(j==0)
                {
                xr=x[i]-x[j];
                yr=y[i]-y[j];
                zr=z[i]-z[j];
                d=Sqr(xr)+Sqr(yr)+Sqr(zr);
                }
                else
                {
                xr=x[i]-x[j];
                yr=y[i]-y[j];
                zr=z[i]-z[j];
                dd=Sqr(xr)+Sqr(yr)+Sqr(zr);
                }
                if(dd<d)d=dd;
            }
        }
    }

    for(i=0;i<N;i++)
    {
        x0[i]=x[i];
        y0[i]=y[i];
        z0[i]=z[i];
    }
    for(i=0;i<N;i++)
    {
    ctcf[i]=0;
    fwd[i]=0;
    bwd[i]=0;
    pf[i]=0.0;
    pb[i]=0.0;
    cohesin[i][0]=0;
    cohesin[i][0]=0;
    cohesin[i][0]=0;
    }
    //ctcf value 0  => ctcf absent
    //ctcf value 1  => ctcf present
    //ctcf value 2  => cohesin present
    //fwd value 1   => forward CTCF
    //bwd value 1   => backward CTCF
    
    nctcf=0;
    // Read the positions of CTCF with forward orientation
    ip=fopen("fwd.txt","r");
    //fscanf(ip,"%d\n",&nfwd);
    for(i=0;i<nfwd;i++)
    {
    double buf;
        fscanf(ip,"%d %lf\n",&k,&buf);
        if(ran2(&idum)<buf)
        {   
            k+=rand()%3-1;
            pf[k]=1.0;
            fwd[k]=1;
            ctcf[k]=1;
            nctcf++;
        }
    }
    fclose(ip);
    // Read the positions of CTCF with backward orientation
    ip=fopen("bwd.txt","r");
    //fscanf(ip,"%d\n",&nbwd);
    for(i=0;i<nbwd;i++)
    {
    double buf;
        fscanf(ip,"%d %lf\n",&k,&buf);
        if(ran2(&idum)<buf)
        {
            k+=rand()%3-1;
            pb[k]=1.0;
            bwd[k]=1;
            ctcf[k]=1;
            nctcf++;
        }
    }
    fclose(ip);


    // Read the positions of CTCF with forward orientation
    double p27,p36;
    scale = 1.0;
    ip1=fopen("./h3k27me3.txt","r");
    ip2=fopen("./h3k36me3.txt","r");
    for(i=0;i<N;i++)
    {
        fscanf(ip1,"%*d %lf\n",&p27);
        fscanf(ip2,"%*d %lf\n",&p36);
        rnum = ran2(&idum);
        p27 = p27 * scale;
        p36 = p36 * scale;
        if(rnum<p36)
            s[i]=2;     //K36 mark
        else if(rnum<p27+p36)
            s[i]=1;     //K27 mark
        if(s[i]>0)
        {
        count++;
        }
    }
    fclose(ip1);
    fclose(ip2);

    ncoh=0;
    ns0=N-2*nctcf-3*nc0;
    A0=(double)(2*nc0-ns0);
    
BuildNbr();

last=0;
diff=0;
    //Monte Carlo steps
    for(n=0;n<nstep;n++)
    {
        //Dump the positions
        if(n%5000==0  && n>=5000000)
        {
        sprintf(file,"r%d.xyz",world_rank);
        op=fopen(file,"a");
        if(op==NULL)printf("ERROR: File %s not found\n",file);
        fprintf(op,"%d\n%d\n",N,n);
            for(i=0;i<N;i++)
            {
            fprintf(op,"%d %lf %lf %lf\n",ctcf[i],x[i],y[i],z[i]);
            }
        fclose(op);
        }

        for(i=0;i<N;i++)
        {
        maxr=0.0;
        //Move 1: Change position of the bead
        id=(int)(ran2(&idum)*N)%N;
        xo=x[id];
        yo=y[id];
        zo=z[id];
        Eo=Energy(xo,yo,zo,id);
        xr=(2.0*ran2(&idum)-1.0);
        yr=(2.0*ran2(&idum)-1.0);
        zr=(2.0*ran2(&idum)-1.0);
        xn=x[id]+xr*dr;
        yn=y[id]+yr*dr;
        zn=z[id]+zr*dr;
        En=Energy(xn,yn,zn,id);
        dE=En-Eo;
        dE_b=dE*beta;
            if(dE<=0.0)
            {
            E=E+dE;
            x[id]=xn;
            y[id]=yn;
            z[id]=zn;
            // Compute maximum displacement
            rr=Sqr(xn-x0[id])+Sqr(yn-y0[id])+Sqr(zn-z0[id]);
            rr=sqrt(rr);
                if(maxr<rr)
                maxr=rr;
            }
            else if(ran2(&idum)<exp(-dE_b))
            {
            E=E+dE;
            x[id]=xn;
            y[id]=yn;
            z[id]=zn;
            // Compute maximum displacement
            rr=Sqr(xn-x0[id])+Sqr(yn-y0[id])+Sqr(zn-z0[id]);
            rr=sqrt(rr);
                if(maxr<rr)
                maxr=rr;
            }

        // Update neighbor list if maxr>Rs/2
            if(maxr>(Rs/2.0))
            {
            BuildNbr();
                for(k=0;k<N;k++)
                {
                    x0[k]=x[k];
                    y0[k]=y[k];
                    z0[k]=z[k];
                }
            }
        //average number of available sites
        nsite=N-2*nctcf-3*ncoh;
        eta=(double)ns0/(double)nc0-1.0;
        s0=(double)(2*ncoh-nsite)-A0;
            if(eta+exp(mu*s0)==0.0)printf("ERROR: denominator is zero in p_bind");
        p_bind = r_bind/(eta+exp(mu*s0)) ; //binding
        //Move 2: Bind new cohesin
        rnum=ran2(&idum);
            if(rnum<p_bind)		// Bind new cohesin
            {
            id=rand()%(N-1);
                if(ctcf[id]==0 && ctcf[id+1]==0)
                {
                BindCohesin(id);
                }
            }
        //compute number of available sites
        nsite=N-2*nctcf-3*ncoh;
        s0=(double)(2*ncoh-nsite)-A0;
            if(1.0+exp(-mu*s0)==0.0)printf("ERROR: denominator is zero in p_unbind");
        p_unbind = r_bind/(1.0+exp(-mu*s0)) ; //dissociation
        //Move 3: Dissociate bound cohesin
        rnum=ran2(&idum);
            if(rnum<p_unbind)
            {
            id=rand()%(N-1);
                if(ctcf[id]==2)
                {
                    UnbindCohesin(id);
                }
            }
        }

	    //Move 4: slide the cohesin
	    for(i=0;i<ncoh;i++)
	    {
            //slide left head of cohesin
	    	if(ncoh>0)
		    {
            	id=rand()%ncoh;
            	rnum=ran2(&idum);
	            if(rnum<r_slide)
	            {
            	    SlideLeftCohesin(id);
            	}
		    }
            //slide right head of cohesin
	        if(ncoh>0)
	        {
                id=rand()%ncoh;
                rnum=ran2(&idum);
	            if(rnum<r_slide)
	            {
            	    SlideRightCohesin(id);
            	}
	        }
        }
        if(n%50000==0 && world_rank==0)printf("%d\n",n);
    }
//Free memory
free(x); 
free(y); 
free(z); 
free(x0); 
free(y0); 
free(z0); 
free(ctcf); 
free(cnbr); 
    for(i=0;i<N;i++)
        free(cohesin[i]);
    for(i=0;i<N;i++)
        free(nbr[i]);
free(cohesin);
free(nbr);

MPI_Finalize();
return 0;
}

double Energy(double xo,double yo,double zo,int i)
{
double en=0.0;
double xr,yr,zr,r2,r2i,r6i;
double tx1,tx2;
double ty1,ty2;
double tz1,tz2;
double dot;
double ep,ec;
int j,k;
int left,right;
en_sp=0.0;
en_lj=0.0;

    for(k=0;k<cnbr[i];k++)
    {
        j=nbr[i][k];
        if(j!=i)
        {
        xr=xo-x[j];
        yr=yo-y[j];
        zr=zo-z[j];
        r2=Sqr(xr)+Sqr(yr)+Sqr(zr);
            if(r2<rc2 && abs(i-j)>1)
            {
                r2i=Sqr(sigma)/r2;
                r6i=r2i*r2i*r2i;
                //LJ potential
                if(s[i]+s[j]==3)
                {
                    ep=epsab;
                    ec=ecutab;
                    en=en+4.0*ep*r6i*(r6i-1.0)-ec;
                }
                else if(s[i]==1 && s[j]==1)
                {
                    ep=epsa;
                    ec=ecuta;
                    en=en+4.0*ep*r6i*(r6i-1.0)-ec;
                }
                else if(s[i]==2 && s[j]==2)
                {
                    ep=epsb;
                    ec=ecutb;
                    en=en+4.0*ep*r6i*(r6i-1.0)-ec;
                }
                else
                {
                    ep=eps;
                    ec=ecut;
                    en=en+4.0*ep*r6i*(r6i-1.0)-ec;
                }
            }
            //Spring potential
            if(j==i-1)
            {
            en=en+0.5*ks*Sqr((sqrt(r2)-r0));
            }
            if(j==i+1)
            {
            en=en+0.5*ks*Sqr((sqrt(r2)-r0));
            }
        }
    }
    //Spring between cohesin heads
    for(j=0;j<ncoh;j++)
    {
    left=cohesin[j][1];
    right=cohesin[j][2];
        if(i==left && left>=0 && right >=0 && left<N && right<N)
	    {
	    xr=xo-x[right];
	    yr=yo-y[right];
	    zr=zo-z[right];
	    r2=Sqr(xr)+Sqr(yr)+Sqr(zr);
        en_sp=en;
	    en=en+0.5*ks*Sqr((sqrt(r2)-r0));
        en_sp=en-en_sp;
            if(en_sp>100000)printf("Cohesin spring\n");
	    }
        if(i==right && left>=0 && right >=0 && left<N && right<N)
	    {
	    xr=xo-x[left];
	    yr=yo-y[left];
	    zr=zo-z[left];
	    r2=Sqr(xr)+Sqr(yr)+Sqr(zr);
        en_sp=en;
	    en=en+0.5*ks*Sqr((sqrt(r2)-r0));
        en_sp=en-en_sp;
            if(en_sp>100000)printf("Cohesin spring\n");
	    }
    }
return en;
}


void BindCohesin(int id)
{
ctcf[id]=2;
ctcf[id+1]=2;
cohesin[ncoh][0]=1;
cohesin[ncoh][1]=id;
cohesin[ncoh][2]=id+1;
ncoh++;
}
void UnbindCohesin(int id)
{
int j;
int left,right;
    for(j=0;j<ncoh;j++)
    {
        //If id is head of LEF
        if(cohesin[j][1]==id || cohesin[j][2]==id)
        {
	    left=cohesin[j][1];
	    right=cohesin[j][2];
	    ctcf[left]=0;
	    ctcf[right]=0;
        if(fwd[left]==1 || bwd[left]==1)
            ctcf[left]=1;
        if(fwd[right]==1 || bwd[right]==1)
            ctcf[right]=1;
	    cohesin[j][1]=cohesin[ncoh-1][1];
	    cohesin[j][2]=cohesin[ncoh-1][2];
	    cohesin[ncoh-1][0]=0;
	    ncoh--;
        }
    }
}
void SlideRightCohesin(int id)
{
int left,right;
right=cohesin[id][2];
    //dissociate if reaches end of the polymer
    if(right==N-1)
    {
        UnbindCohesin(right);
    }
    else if(ctcf[right+1]!=2)
        if(bwd[right+1]==0 || ran2(&idum)>pb[right+1])
        {
        ctcf[right]=0;
            if(fwd[right]==1 || bwd[right]==1)
                ctcf[right]=1;
        right++;
        ctcf[right]=2;
        cohesin[id][2]=right;
        }
}
void SlideLeftCohesin(int id)
{
int left,right;
left=cohesin[id][1];
    //dissociate if reaches end of the polymer
    if(left==0)
    {
        UnbindCohesin(left);
    }
    else if(ctcf[left-1]!=2)
        if(fwd[left-1]==0 || ran2(&idum)>pf[left-1])
        {
        ctcf[left]=0;
            if(fwd[left]==1 || bwd[left]==1)
                ctcf[left]=1;
        left--;
        ctcf[left]=2;
        cohesin[id][1]=left;
        }
}


void BuildNbr()
{
int i,j,k;
double Rc,Rc2;
double xr,yr,zr,rr;
Rc=rc+Rs;
Rc2=Sqr(Rc);
    for(i=0;i<N;i++)
        cnbr[i]=0;
    for(i=0;i<N-1;i++)
    {
        for(j=i+1;j<N;j++)
        {
            xr=x[i]-x[j];
            yr=y[i]-y[j];
            zr=z[i]-z[j];
            rr=Sqr(xr)+Sqr(yr)+Sqr(zr);
            if(rr<Rc2)
            {
                if(cnbr[i]>maxnbr || cnbr[j]>maxnbr)
                {
                    printf("ERROR: number of neighbors %d, %d greater than %d\n",cnbr[i],cnbr[j],maxnbr);
                    exit(EXIT_FAILURE);
                }
                nbr[i][cnbr[i]]=j;
                nbr[j][cnbr[j]]=i;
                cnbr[i]++;
                cnbr[j]++;
            }
        }
    }
}


