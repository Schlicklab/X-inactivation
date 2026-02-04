#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#//include <gsl/gsl_fit.h>

#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))	//Nearest integer
#define Sqr(x)     ((x) * (x))

//void computeDiffusionCoefficient(int N, int num_samples, double time_step);
int LinesIn(FILE** fp)
{
char c;
int count;
    // Check if file exists
    if (*fp == NULL)
    {
        printf("Could not open file\n");
        return 0;
    }
    for (c = getc(*fp); c != EOF; c = getc(*fp))
    if (c == '\n') // Increment count if this character is newline
    count = count + 1;
return count;
}

int world_rank=0;								//rank of the process in MPI_COMM_WORLD
int world_size=0;								//size of MPI_COMM_WORLD
int main()
{
MPI_Init(NULL,NULL);
MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
MPI_Comm_size(MPI_COMM_WORLD,&world_size);
int i,j,k;
int ii;
int n;
int N;
N=500;
int s[N];
int t,nn;
int count,ncount;
int nlines;
int n_sys;
int n_ens[4000];
int nframes,maxframe;
int sk;
double old;
double bin;
double sum;
double pp;
//double x0[N],y0[N],z0[N];
double x[4000][2],y[4000][2],z[4000][2];
int dist[1000];
int gdist[1000];
int cdist;
double xr,yr,zr,rr;
double xi,yi,zi;
double xj,yj,zj;
double sigma,rc,rc2;
int start, end;
FILE *ip,*op;
char file[100];

sigma=1.0;
rc=1.5*sigma;
rc2=rc*rc;
n_sys=120;  //Number of trajectories
//Xist Ftx
//start=189;
//end=260;
//Xist Jpx
//start=189;
//end=195;
//Tsix linx
start=100;
end=170;
sk=end-start;
bin=0.2;

count=0;
cdist=0;

nframes=1000; 
maxframe=1;
  //for(ii=0;ii<n_sys;ii++)
  ii=world_rank;
  {
    sprintf(file,"r%d.xyz",ii);
    ip=fopen(file,"r");
    if(ip==NULL)printf("ERROR: File %s not found\n",file);
    nlines=LinesIn(&ip);
    nframes=nlines/(N+2);
        if(nframes<1000)printf("ERROR: File %s incomplete\n",file);
        else nframes=1000;
    fclose(ip);
    ip=fopen(file,"r");
    ncount=0;
    while(ncount!=nframes) //Total number of snapshots in 1 trajectory
    {
    fscanf(ip,"%d\n%d\n",&nn,&t);
        for(i=0;i<N;i++)
        {
            if(i==start)
            fscanf(ip,"%d %lf %lf %lf\n",&s[i],&x[ncount][0],&y[ncount][0],&z[ncount][0]);
            else if(i==end)
            fscanf(ip,"%d %lf %lf %lf\n",&s[i],&x[ncount][1],&y[ncount][1],&z[ncount][1]);
            else
            fscanf(ip,"%d %*lf %*lf %*lf\n",&s[i]);
        }
    ncount++;
    }
    fclose(ip);
    for(int ti=0;ti<nframes-1;ti+=1)
    {
        {
            i=0;
            {
                j=1;
                {
                xi=x[ti][j]-x[ti][i];
                yi=y[ti][j]-y[ti][i];
                zi=z[ti][j]-z[ti][i];

                rr=Sqr(xi)+Sqr(yi)+Sqr(zi);
                rr=sqrt(rr);
                dist[(int)(rr/bin)]++;
                cdist++;

                }
            }
        }
    }
    if(maxframe<nframes)maxframe=nframes;
  }


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&dist, &gdist, 1000, MPI_INT, MPI_SUM, 0,MPI_COMM_WORLD);
    if(ii==0)
    {
    sprintf(file,"dist_tsix_linx.txt");
    op=fopen(file,"w");
        for(i=0;i<1000;i++)
        {
            double tcount;
            tcount = (double)(world_size*cdist);
            fprintf(op,"%lf %lf\n",bin*(double)i,(double)gdist[i]/(tcount*bin));
        }
    fclose(op);
    }

MPI_Barrier(MPI_COMM_WORLD);
MPI_Finalize();
return 0;
}

