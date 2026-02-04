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
int count,ncount;
int s_count[N];
int n_sys;
int n_ens;
int nlines;
int J,mu;
int nupdate;
int nframes;
int ndist[450][450];
double P[N][N];
double dist[450][450];
double Pe[450][450];
double bin;
double sum;
double pp;
double s_P[N];
double x[N],y[N],z[N];
double xr,yr,zr,rr;
double sigma,rc,rc2;
double buf;
double dt=0.005*1000.0;
int nuc_pos[N];
FILE *ip,*op;
rc=25.0;
rc2=rc*rc;
bin=200.0;
n_sys=1;
nframes=200;
//printf("nlines %d\n",nlines);
n_ens=0;
count=0;
int num_cores, num_dna[N],cum_dna[N];


    for(i=0;i<N;i++) 
    {
        for(j=0;j<N;j++)
        {
        P[i][j]=0.0;
        s_count[i]=0;
        s_P[i]=0.0;
        }
    }
    for(i=0;i<450;i++) 
        for(j=0;j<450;j++)
        {
        Pe[i][j]=0.0;
        dist[i][j]=0.0;
        ndist[i][j]=0;
        }
    for(ii=0;ii<n_sys;ii++)
    {
        sprintf(file,"run%d/dim.in",ii+1);
        FILE *dimFile = fopen(file, "r");
        if(dimFile==NULL)continue;
        nuc_pos[0]=0;
        int nlink=0;
        fscanf(dimFile, "#cores\n%d\n#LH\n0 6 22\n#DNA\n", &num_cores);
        N=num_cores;
        for (int i = 0; i < num_cores; i++) {
            fscanf(dimFile, "%d\n", &num_dna[i]);
            if(i==0)
                cum_dna[i]=num_dna[i]+1;
            else 
                cum_dna[i]=num_dna[i]+cum_dna[i-1]+1;
            nlink+=num_dna[i];
            if(i<N-1)
            {
                xr=147.0*((double)i+1.5) + 9.0*((double)cum_dna[i]);
                nuc_pos[i+1]= (int)(xr/bin);
                if(nuc_pos[i+1]>=450)printf("Error: nucleosome position %d larger than region size\n",nuc_pos[i+1]);
            }
        }
        fclose(dimFile);

        
        //sprintf(file,"traj.dat");
        sprintf(file,"xyz/traj%d.dat",ii+1);
        nlines=LinesIn();
        nframes=nlines/(N*82+4*nlink);
        printf("%s %d\n",file,nframes);
        if(nframes<100)continue;
        //nframes=100;
        ip=fopen(file,"r");
        ncount=0;
        if(ip==NULL)
        {
            printf("ERROR: file %s not found\n",file);
            continue;
        }
        while(ncount!=nframes) 
        {
        ncount++;
            for(i=0;i<N;i++)
            {
            fscanf(ip,"%lf %lf %lf\n",&x[i],&y[i],&z[i]);
                for(j=0;j<3;j++)
                fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore rotation
            }
            for(j=0;j<(nlink*4+N*78);j++)
                fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore linker beads, histone tails and LH
            if(ncount>0) //For averaging system equilibrates after 100 snapshots
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
                        Pe[nuc_pos[i]][nuc_pos[j]]+=1.0;
                        Pe[nuc_pos[j]][nuc_pos[i]]+=1.0;
                        }
                        dist[nuc_pos[i]][nuc_pos[j]]+=rr;
                        dist[nuc_pos[j]][nuc_pos[i]]+=rr;
                        ndist[nuc_pos[i]][nuc_pos[j]]++;
                        ndist[nuc_pos[j]][nuc_pos[i]]++;
                    }
                }
            j=0;
            }
        }
    fclose(ip);
    }

  sprintf(file,"pc.txt");
  op=fopen(file,"w");
  for(i=0;i<450;i++)
  {
    for(j=0;j<450;j++)
    {
        //P[i][j]=P[i][j]/(double)n_ens;
        //fprintf(op,"%d %d %lf\n",i,j,P[i][j]);
        Pe[i][j]=Pe[i][j]/(double)n_ens;
        fprintf(op,"%d %d %lf\n",i,j,Pe[i][j]);
    }
  }
  fclose(op);

  sprintf(file,"pd.txt");
  op=fopen(file,"w");
  for(i=0;i<450;i++)
  {
    for(j=0;j<450;j++)
    {
        dist[i][j]=dist[i][j]/(double)ndist[i][j];
        fprintf(op,"%d %d %lf\n",i,j,dist[i][j]);
    }
  }
  fclose(op);



  for(i=0;i<450-1;i++)
  {
    for(j=i+1;j<450;j++)
    {
        s_P[j-i]+=Pe[i][j];
        s_count[j-i]++;
    }
  }
  sprintf(file,"pc_S.txt");
  op=fopen(file,"w");
  for(i=0;i<450;i++)
    fprintf(op,"%d %lf\n",i,s_P[i]/(double)s_count[i]);
  fclose(op);

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



