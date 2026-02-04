#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))	//Nearest integer

int LinesIn();

char file[100];
int main()
{
int i,j,k;
int ii;
int n;
int N;
N=500;
int s[N];
int t,nn;
int n_ens[2000],ncount;
int count;
int tmax,tmin;
int nlines,nframes;
int n_sys;
int nupdate;
int flag;
int sum;
int bin;
int s0;
int ss,ee;
double pp;
double x[N],y[N],z[N];
double xl[3000],yl[3000],zl[3000];
double xr,yr,zr,rr,r2;
double cx,cy,cz;
double rg[2000];
double rg2[2000];
double sigma,rc,rc2;
int nuc_pos[N];
int ll_pos[3000];
int num_cores, num_dna[N],cum_dna[N];
FILE *ip,*op;
sigma=1.0;
bin=200.0;
rc=1.5*sigma;
rc2=rc*rc;
n_sys=300;
nframes=1000;
double mlink=6.0;
double mcore=206.0;

    ip=fopen("../../src/gene1_pos.txt","r");
    fscanf(ip,"%d %d",&ss,&ee);
    fclose(ip);
    op=fopen("gene1_dist.txt","w");
    fclose(op);
  for(i=0;i<2000;i++) 
    {
    rg[i]=0.0;
    rg2[i]=0.0;
    n_ens[i]=0;
    }
  for(ii=0;ii<n_sys;ii++)
  {
    sprintf(file,"run%d/dim.in",ii+1);
    FILE *dimFile = fopen(file, "r");
    if(dimFile==NULL)continue;
    nuc_pos[0]=103425000+147/2;
    int nlink=0;
    fscanf(dimFile, "#cores\n%d\n#LH\n0 6 22\n#DNA\n", &num_cores);
    N=num_cores;
    for (int i = 0; i < num_cores; i++) {
        fscanf(dimFile, "%d\n", &num_dna[i]);
        if(i!=0)
        {
            nuc_pos[i]=nuc_pos[i-1]+147+9*(num_dna[i]+1);
        }
        //nlink+=num_dna[i];
        for(int k=0;k<num_dna[i];k++)
        {
            if(k==0)
            {
            ll_pos[nlink]=nuc_pos[i]+147/2+9;
            nlink++;
            }
            else
            {
            ll_pos[nlink]=nuc_pos[i]+9;
            nlink++;
            }
        }
    }
    fclose(dimFile);


    sprintf(file,"xyz/traj%d.dat",ii+1);
    nlines=LinesIn();
    nframes=nlines/(N*82+4*nlink);
    printf("%s frames=%d Cores=%d Linkers=%d \n",file,nframes,N,nlink);
    if(nframes<100)continue;
    //nframes=100;
    ip=fopen(file,"r");
    ncount=0;
        if(ip==NULL)
        {
            printf("ERROR: file %s not found\n",file);
            continue;
        }
    ncount=0;
    count=0;
    while(ncount!=nframes) 
    {
        int start,end,size;
        int lstart,lend; //linker start and end
        start=0;
        end=0;
    ncount++;
        for(i=0;i<N;i++)
        {
            if(nuc_pos[i]>=ss && start==0) start=i;
            if(nuc_pos[i]>=ss && nuc_pos[i]<=ee) end=i;
        fscanf(ip,"%lf %lf %lf\n",&x[i],&y[i],&z[i]); //core positions
            for(j=0;j<3;j++)
            fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore rotation
        }
        for(j=0;j<(nlink);j++)
        {
            fscanf(ip,"%lf %lf %lf\n",&xl[j],&yl[j],&zl[j]); // linker positions 
            if(ll_pos[j]>=ss && lstart==0) lstart=j;
            if(ll_pos[j]>=ss && ll_pos[j]<=ee) lend=j;
            for(k=0;k<3;k++)
            fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore linker ratation
        }
        for(j=0;j<(N*78);j++)
            fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore histone tails and LH
    rr=0.0;
    cx=0.0;
    cy=0.0;
    cz=0.0;
    double mass=0.0;
        //printf("ncount %d\n",ncount);
        for(i=start;i<=end;i++)
        {
        cx=cx+x[i]*mcore;
        cy=cy+y[i]*mcore;
        cz=cz+z[i]*mcore;
        mass+=mcore;
        }
        for(i=lstart;i<=lend;i++)
        {
        cx=cx+xl[i]*mlink;
        cy=cy+yl[i]*mlink;
        cz=cz+zl[i]*mlink;
        mass+=mlink;
        }
    cx=cx/mass;
    cy=cy/mass;
    cz=cz/mass;
        for(i=start;i<=end;i++)
        {
	    xr=x[i]-cx;
	    yr=y[i]-cy;
	    zr=z[i]-cz;
	    rr=rr+mcore*(xr*xr+yr*yr+zr*zr);
        }
        for(i=lstart;i<=lend;i++)
        {
	    xr=xl[i]-cx;
	    yr=yl[i]-cy;
	    zr=zl[i]-cz;
	    rr=rr+mlink*(xr*xr+yr*yr+zr*zr);
        }
    rr=rr/(double)mass;
    rg2[count]+=rr;
    rr=sqrt(rr);
    rg[count]+=rr;
    n_ens[count]++;
    count++;

    op=fopen("gene1_dist.txt","a");
        if(op==NULL)printf("file not found\n");
    fprintf(op,"35 %lf\n",rr);
    fclose(op);
    }
  fclose(ip);
  }

return 0;
}


int LinesIn()
{
char c;
int count=0;
FILE *fp;
    // Check if file exists
    fp=fopen(file,"r");
    if (fp == NULL)
    {
        printf("Could not open file\n");
        return 0;
    }
    for (c = getc(fp); c != EOF; c = getc(fp))
    if (c == '\n') // Increment count if this character is newline
    count = count + 1;
    fclose(fp);
return count;
}

