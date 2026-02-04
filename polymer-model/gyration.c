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
int t,nn,clust_id[N];
int clust_size[N];
int size_dist[N];
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
double pp;
double x[N],y[N],z[N];
double xr,yr,zr,rr,r2;
double cx,cy,cz;
double rg[2000];
double rg2[2000];
double sigma,rc,rc2;
FILE *ip,*op;
sigma=1.0;
rc=1.5*sigma;
rc2=rc*rc;
n_sys=120;
nframes=1000;
  for(i=0;i<2000;i++) 
    {
    rg[i]=0.0;
    rg2[i]=0.0;
    n_ens[i]=0;
    }
  for(ii=0;ii<n_sys;ii++)
  {
    sprintf(file,"r%d.xyz",ii);
    nlines=LinesIn();
    nframes=nlines/(N+2);
        //Do not read file if frames less than 1000
        if(nframes<1000)
        {
            printf("file %s nlines %d %d nframes\n",file,nlines,nframes);
            continue;
        }
        nframes=1000;
    ip=fopen(file,"r");
        if(ip==NULL)continue;
    ncount=0;
    count=0;
    while(ncount<(N+2)*nframes)
    {
    fscanf(ip,"%d\n%d\n",&nn,&t);
    ncount=ncount+N+2;
        for(i=0;i<N;i++)
        {
        fscanf(ip,"%d %lf %lf %lf\n",&s[i],&x[i],&y[i],&z[i]);
        clust_id[i]=i;
        }
    rr=0.0;
    cx=0.0;
    cy=0.0;
    cz=0.0;
        //printf("ncount %d\n",ncount);
        for(i=0;i<N;i++)
        {
        cx=cx+x[i];
        cy=cy+y[i];
        cz=cz+z[i];
        }
    cx=cx/(double)N;
    cy=cy/(double)N;
    cz=cz/(double)N;
        for(i=0;i<N;i++)
        {
	    xr=x[i]-cx;
	    yr=y[i]-cy;
	    zr=z[i]-cz;
	    rr=rr+xr*xr+yr*yr+zr*zr;
        }
    rr=rr/(double)N;
    rg2[count]+=rr;
    rr=sqrt(rr);
    rg[count]+=rr;
    n_ens[count]++;
    count++;

    sprintf(file,"rg_dist.txt");
    op=fopen(file,"a");
    if(ncount>(N+2)*500)
    fprintf(op,"35 %lf\n",rr);
    fclose(op);
    }
  fclose(ip);
  }

    sprintf(file,"rg_time.txt");
    op=fopen(file,"w");
    for(i=0;i<nframes;i++)
    {
    rg[i]=rg[i]/(double)n_ens[i];
    rg2[i]=rg2[i]/(double)n_ens[i];
    fprintf(op,"%d %lf\n",i*50000,rg[i]);
    }
    fclose(op);
  rr=0.0;
  r2=0.0;
  count=0;
  for(i=100;i<nframes;i++)
  {
  count++;
  rr+=rg[i];
  r2+=rg2[i];
  }
  rr=rr/(double)(count);
  r2=r2/(double)(count);
double se;
//se = sqrt(r2-rr*rr);
se = sqrt(r2-rr*rr)/sqrt(count);
  op=fopen("rg_avg.txt","w");
  fprintf(op,"%lf %lf\n",rr,se);
  fclose(op);
return 0;
}


int LinesIn()
{
char c;
int count;
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

