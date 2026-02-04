#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))	//Nearest integer
#define Sqr(x)     ((x) * (x))


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
int count,ncount;
int s_count[N];
int n_sys;
int n_ens;
int nlines;
int n_flip,n_bind,n_slide,ncoh;
int J,mu;
int nupdate;
int nframes;
double P[N][N];
double Pe[N][N];
double bin;
double sum;
double pp;
double s_dist[N];
double s_P[N];
double x[N],y[N],z[N];
double xr,yr,zr,rr;
double dist[N][N];
double cx,cy,cz,rg[5000];
double sigma,rc,rc2;
double buf;
double dt=0.005*1000.0;
FILE *ip,*op;
sigma=1.0;
rc=1.2*sigma;
rc2=rc*rc;
n_sys=120;  //Number of trajectories
bin=0.1;
nframes=1000;
n_ens=0;
count=0;


  for(i=0;i<5000;i++) 
  rg[i]=0;


  for(i=0;i<N;i++) 
  {
    for(j=0;j<N;j++)
    {
    P[i][j]=0.0;
    dist[i][j]=0;
    s_count[i]=0;
    s_P[i]=0.0;
    s_dist[i]=0.0;
    }
  }
    
  for(ii=0;ii<n_sys;ii++)
  {
    sprintf(file,"r%d.xyz",ii);
    nlines=0;
    nlines=LinesIn();
    nframes=nlines/(N+2);
    if(nframes<500)
    {
        printf("file %d frames %d\n",ii+1,nframes);
        continue;
    }
    ip=fopen(file,"r");
    ncount=0;
    if(ip==NULL)
    {
        printf("file %s not found\n",file);
        continue;
    }
    while(ncount!=(N+2)*nframes) //Total 500 snapshots of the system in 1 trajectory
    {
    ncount=ncount+N+2;
    fscanf(ip,"%d\n",&nn);
    fscanf(ip,"%d\n",&t);
        for(i=0;i<N;i++)
        {
        fscanf(ip,"%d %lf %lf %lf\n",&s[i],&x[i],&y[i],&z[i]);
        clust_id[i]=i;
        }
        if(ncount>(N+2)*0 ) //For averaging system equilibrates after 200 snapshots
        {
            for(i=0;i<N-1;i++)
            {
    	        for(j=i+1;j<N;j++)
	            {
	            xr=x[i]-x[j];
	            yr=y[i]-y[j];
	            zr=z[i]-z[j];
	            rr=xr*xr+yr*yr+zr*zr;
                rr=sqrt(rr);
                    if(i==0 && j==1)
                    n_ens++;
                    //contact probability matrix
                    if(rr<rc)
                    {
                    P[i][j]+=1.0;
                    P[j][i]+=1.0;
                    }
                    //if(rr<0.6)printf("frame%d i %d j %d rr %lf\n",ncount/(N+2),i,j,rr);
                // Distance matrix
                dist[i][j]=dist[i][j]+rr;
                dist[j][i]=dist[i][j];
                }
            }
        j=0;
        }
    }
  fclose(ip);
  }

  sprintf(file,"pc.txt");
  op=fopen(file,"w");
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
        P[i][j]=P[i][j]/(double)n_ens;
        fprintf(op,"%d %d %lf\n",i,j,P[i][j]);
    }
  }
  fclose(op);

/*
  sprintf(file,"pd_nconf.dat");
  op=fopen(file,"w");
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
        fprintf(op,"%d %d %lf\n",i,j,dist[i][j]/(double)n_ens);
  }
  fclose(op);

  for(i=0;i<N-1;i++)
  {
    for(j=i+1;j<N;j++)
    {
        s_P[j-i]+=P[i][j];
        s_dist[j-i]+=dist[i][j]/(double)n_ens;
        s_count[j-i]++;
    }
  }
  sprintf(file,"pc_S.dat");
  op=fopen(file,"w");
  for(i=0;i<N;i++)
    fprintf(op,"%d %lf\n",i,s_P[i]/(double)s_count[i]);
  fclose(op);


  sprintf(file,"dist_S.dat");
  op=fopen(file,"w");
  for(i=0;i<N;i++)
    fprintf(op,"%d %lf\n",i,s_dist[i]/(double)s_count[i]);
  fclose(op);

*/

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



