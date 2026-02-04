#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))    //Nearest integer
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
    int n_ens=0,ncount;
    int count,nconf=0;
    int tmax,tmin;
    int nlines,nframes;
    int n_sys;
    int flag;
    int sum;
    int bin;
    int s0=103425000;
    double pp;
    int size_lim=3;
    double x[N],y[N],z[N];
    int P[N];
    int Pc[N];
    double xl[3000],yl[3000],zl[3000];
    double xr,yr,zr,rr,r2;
    double cx1,cy1,cz1;
    double cx2,cy2,cz2;
    double rg[2000];
    double rg2[2000];
    double sigma,rc,rc2;
    int nuc_pos[N];
    int ll_pos[3000];
    int num_cores, num_dna[N],cum_dna[N];
    int ss1,ee1;
    int ss2,ee2;
    int ss3,ee3;
    int nupdate;
    FILE *ip,*op;
    
    // Arrays to count overlap for majority rule
    int overlap1[N];
    int overlap2[N];
    int overlap3[N];

    bin=4;
    sigma=1.0;
    rc=20.0;
    rc2=rc*rc;
    n_sys=300;
    nframes=1000;
    int gene1_dist[N];
    int gene2_dist[N];
    int gene3_dist[N];
    int type[N];

    for(i=0;i<N;i++) 
    {
        size_dist[i]=0;
        gene1_dist[i]=0;
        gene2_dist[i]=0;
        gene3_dist[i]=0;
        type[N]=0;
        P[i]=0;
        Pc[i]=0;
    }

    // Read gene locus positions
    ip=fopen("../../src/lxite_pos.txt","r");
    fscanf(ip,"%d %d",&ss1,&ee1);
    fclose(ip);
    ip=fopen("../../src/lxist_pos.txt","r");
    fscanf(ip,"%d %d",&ss2,&ee2);
    fclose(ip);
    ip=fopen("../../src/ljpx_pos.txt","r");
    fscanf(ip,"%d %d",&ss3,&ee3);
    fclose(ip);

    for(i=0;i<2000;i++) 
    {
        rg[i]=0.0;
        rg2[i]=0.0;
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
        
        nframes=100;
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
                clust_id[i]=i;
                fscanf(ip,"%lf %lf %lf\n",&x[i],&y[i],&z[i]); //core positions
                for(j=0;j<3;j++)
                    fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore rotation
            }
            for(j=0;j<(nlink);j++)
            {
                fscanf(ip,"%lf %lf %lf\n",&xl[j],&yl[j],&zl[j]); // linker positions 
                for(k=0;k<3;k++)
                    fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore linker rotation
            }
            for(j=0;j<(N*78);j++)
                fscanf(ip,"%lf %lf %lf\n",&xr,&yr,&zr); //Ignore histone tails and LH

            if(ncount>0 && ncount<100)
            {
                for(i=0;i<N;i++)
                {
                    clust_size[i]=0;
                }
                n_ens++;
                nupdate=1; //non-zero update
                
                // --- CLUSTERING LOGIC ---
                while(nupdate!=0)
                {
                    nupdate=0;
                    for(i=0;i<N-1;i++)
                    {
                        flag=0;
                        for(j=i+1;j<N;j++)
                        {
                            xr=x[i]-x[j];
                            yr=y[i]-y[j];
                            zr=z[i]-z[j];
                            rr=xr*xr+yr*yr+zr*zr;
                            if(rr<rc2)
                            {
                                flag=1;
                                if(clust_id[i]<clust_id[j])
                                {
                                    clust_id[j]=clust_id[i];
                                    nupdate++;
                                }
                                else if(clust_id[i]>clust_id[j])
                                {
                                    clust_id[i]=clust_id[j];
                                    nupdate++;
                                }
                            }
                        }
                    }
                }
                // Calculate Cluster Sizes
                for(i=0;i<N;i++)
                {
                    clust_size[clust_id[i]]++;
                }
                // General Distribution
                for(i=0;i<N;i++)
                {
                    if(clust_size[i]>2)
                        size_dist[(clust_size[i]-1)/bin]++;
                    if(clust_size[i]>=20)Pc[(nuc_pos[i]-s0)/200]++;
                }

                // --- NEW GENE CLASSIFICATION LOGIC (Majority Rule) ---
                
                // 1. Reset overlap counters
                for(i=0; i<N; i++) {
                    overlap1[i]=0;
                    overlap2[i]=0;
                    overlap3[i]=0;
                }

                // 2. Count how many nucleosomes of each cluster are inside gene boundaries
                for(i=0; i<N; i++) {
                    int cid = clust_id[i];
                    if(nuc_pos[i] >= ss1 && nuc_pos[i] <= ee1) overlap1[cid]++;
                    if(nuc_pos[i] >= ss2 && nuc_pos[i] <= ee2) overlap2[cid]++;
                    if(nuc_pos[i] >= ss3 && nuc_pos[i] <= ee3) overlap3[cid]++;
                }

                // 3. Apply Threshold and Update Distributions
                // We iterate i from 0 to N. If clust_size[i] > 0, 'i' is a valid cluster ID.
                for(i=0; i<N; i++) 
                {
                    if(clust_size[i] > 2) // Maintain the size limit check from your original code
                    {
                        double threshold = (double)clust_size[i] * 0.5;
                        int bin_idx = (clust_size[i]-1)/bin;

                        // Check Gene 1
                        if((double)overlap1[i] > threshold) {
                            gene1_dist[bin_idx]++;
                        }
                        // Check Gene 2
                        if((double)overlap2[i] > threshold) {
                            gene2_dist[bin_idx]++;
                        }
                        // Check Gene 3
                        if((double)overlap3[i] > threshold) {
                            gene3_dist[bin_idx]++;
                        }
                    }
                }

                // --- Cluster Probability along chain (Existing Logic) ---
                for(i=0;i<N-1;i++)
                {
                    if(clust_id[i]!=clust_id[i+1])
                    {
                        int xb;
                        xb=(nuc_pos[i]+nuc_pos[i+1])/2;
                        xb=(xb-s0)/200;
                        P[xb]++;
                    }
                }
            }
        }
        fclose(ip);
    }
    
    // Output writing
    sprintf(file,"clust_pb.txt");
    op=fopen(file,"w");
    for(i=0;i<450;i++)
    {
        fprintf(op,"%d %lf\n",i,(double)P[i]/(double)n_ens);
    }
    fclose(op);
    
    sprintf(file,"pclust.txt");
    op=fopen(file,"w");
    for(i=0;i<450;i++)
    {
        fprintf(op,"%d %lf\n",i,(double)Pc[i]/(double)n_ens);
    }
    fclose(op);
    
    sum=0;
    sprintf(file,"size_clust.txt");
    op=fopen(file,"w");
    for(i=0;i<N/bin;i++)
    {
        fprintf(op,"%d %lf\n",i,(double)size_dist[i]/(double)(bin*n_ens));
        sum=0;
    }
    fclose(op);

    sprintf(file,"lxite_clust.txt");
    op=fopen(file,"w");
    for(i=0;i<N/bin;i++)
    {
        fprintf(op,"%d %lf\n",i,(double)gene1_dist[i]/(double)(bin*n_ens));
    }
    fclose(op);
    
    sprintf(file,"lxist_clust.txt");
    op=fopen(file,"w");
    for(i=0;i<N/bin;i++)
    {
        fprintf(op,"%d %lf\n",i,(double)gene2_dist[i]/(double)(bin*n_ens));
    }
    fclose(op);
    
    sprintf(file,"ljpx_clust.txt");
    op=fopen(file,"w");
    for(i=0;i<N/bin;i++)
    {
        fprintf(op,"%d %lf\n",i,(double)gene3_dist[i]/(double)(bin*n_ens));
    }
    fclose(op);
    
    return 0;
}

int LinesIn()
{
    char c;
    int count=0;
    FILE *fp;
    fp=fopen(file,"r");
    if (fp == NULL)
    {
        printf("Could not open file\n");
        return 0;
    }
    for (c = getc(fp); c != EOF; c = getc(fp))
        if (c == '\n') 
            count = count + 1;
    fclose(fp);
    return count;
}
