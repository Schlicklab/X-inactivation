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
double pp;
double x[N],y[N],z[N];
double xr,yr,zr,rr,r2;
double cx,cy,cz;
double rg[2000];
double rg2[2000];
double sigma,rc,rc2;
int nuc_pos[N];
int num_cores, num_dna[N],cum_dna[N];
FILE *ip,*op;
sigma=1.0;
bin=200.0;
rc=1.5*sigma;
rc2=rc*rc;
n_sys=300;
nframes=1000;
  for(i=0;i<2000;i++) 
    {
    rg[i]=0.0;
    rg2[i]=0.0;
    n_ens[i]=0;
    }

    sprintf(file,"run%d/rg.txt",ii+1);
    op=fopen(file,"w");
    fclose(op);
    op=fopen("rg_dist.txt","w");
    fclose(op);
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
            //if(nuc_pos[i+1]>=450)printf("Error: nucleosome position %d larger than region size\n",nuc_pos[i+1]);
        }
    }
    fclose(dimFile);


    sprintf(file,"xyz/traj%d.dat",ii+1);
    nlines=LinesIn();
    nframes=nlines/(N*82+4*nlink);
    printf("%s frames=%d Cores=%d Linkers=%d \n",file,nframes,N,nlink);
    if(nframes<50)continue;
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
    ncount++;
        for(i=0;i<N;i++)
        {
        fscanf(ip,"%lf %lf %lf\n",&x[i],&y[i],&z[i]);
            for(j=0;j<3;j++)
            fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore rotation
        }
        for(j=0;j<(nlink*4+N*78);j++)
            fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore linker beads, histone tails and LH
    rr=0.0;
    cx=0.0;
    cy=0.0;
    cz=0.0;
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

    sprintf(file,"run%d/rg.txt",ii+1);
    op=fopen(file,"a");
        if(op==NULL)printf("file not found\n");
    fprintf(op,"%d %lf\n",ncount*100000,rr);
    fclose(op);
    op=fopen("rg_dist.txt","a");
        if(op==NULL)printf("file not found\n");
        //if(ncount>0)
    fprintf(op,"35 %lf\n",rr);
    fclose(op);
    }
  fclose(ip);
  }

    sprintf(file,"rg_time.txt");
    op=fopen(file,"w");
        if(op==NULL)printf("file not found\n");
    for(i=0;i<nframes;i++)
    {
    rg[i]=rg[i]/(double)n_ens[i];
    rg2[i]=rg2[i]/(double)n_ens[i];
    fprintf(op,"%d %lf\n",i*100000,rg[i]);
    }
    fclose(op);
  rr=0.0;
  r2=0.0;
  count=0;
  for(i=0;i<nframes;i++)
  {
  count++;
  rr+=rg[i];
  r2+=rg2[i];
  }
  rr=rr/(double)(count);
  r2=r2/(double)(count);
double se;
se = sqrt(r2-rr*rr)/sqrt(count);
  op=fopen("rg_avg.txt","w");
        if(op==NULL)printf("file not found\n");
  fprintf(op,"%lf %lf\n",rr,se);
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

