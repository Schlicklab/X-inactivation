#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CORES 500
#define MAX_BEADS 300
#define LINE_LEN 256

int LinesIn();
// Structure to store a 3D vector
typedef struct {
    double x, y, z;
} Vector3D;

// Function to transform core beads based on orientation
Vector3D transform_bead(Vector3D core_pos, Vector3D a, Vector3D b, Vector3D c, Vector3D rel_bead) {
    Vector3D result;
    result.x = core_pos.x + a.x * rel_bead.x + b.x * rel_bead.y + c.x * rel_bead.z;
    result.y = core_pos.y + a.y * rel_bead.x + b.y * rel_bead.y + c.y * rel_bead.z;
    result.z = core_pos.z + a.z * rel_bead.x + b.z * rel_bead.y + c.z * rel_bead.z;
    return result;
}

char file[100];
int main() {
    sprintf(file,"fort.10");
    FILE *fortFile = fopen(file, "r");
    FILE *dimFile = fopen("dim.in", "r");
    FILE *trajFile = fopen("1.dat", "w");
    int nlines, nframes = 350;
	int freq = 50;
    double cx,cy,cz;
    
    if (!dimFile  || !trajFile || !fortFile) {
        printf("Error opening input/output files.\n");
        return 1;
    }
    
    int num_cores, num_dna[MAX_CORES];
    int MAX_DNA=0;
    fscanf(dimFile, "#cores\n%d\n#LH\n0 6 22\n#DNA\n", &num_cores);
    for (int i = 0; i < num_cores; i++) {
        fscanf(dimFile, "%d\n", &num_dna[i]);
        MAX_DNA+=num_dna[i];
    }
    fclose(dimFile);
    
    nlines=LinesIn();
    nframes=nlines/(num_cores*82+4*MAX_DNA);
    
    Vector3D core_positions[MAX_CORES], orientations[MAX_CORES][3];
    Vector3D dna_positions[MAX_DNA], dna_orientations[MAX_DNA][3];
    Vector3D tail_positions[MAX_CORES*50], lh_positions[MAX_CORES*28];
    
    for(int n = 0; n < 1; n++)
    {
        printf("frame %d\n",n);
        cx=0.0;
        cy=0.0;
        cz=0.0;
        for (int i = 0; i < num_cores; i++) {
            fscanf(fortFile, "%lf %lf %lf\n", &core_positions[i].x, &core_positions[i].y, &core_positions[i].z);
            for (int j = 0; j < 3; j++) {
                fscanf(fortFile, "%lf %lf %lf\n", &orientations[i][j].x, &orientations[i][j].y, &orientations[i][j].z);
            }
        }
        int atom_id = 0;
        for (int i = 0; i < num_cores; i++) {
            for (int k = 0; k < num_dna[i]; k++) {
                fscanf(fortFile, "%lf %lf %lf\n", &dna_positions[atom_id].x, &dna_positions[atom_id].y, &dna_positions[atom_id].z);
                fprintf(trajFile,"%e %e %e\n", dna_positions[atom_id].x, dna_positions[atom_id].y, dna_positions[atom_id].z);
                for (int j = 0; j < 3; j++) {
                    fscanf(fortFile, "%lf %lf %lf\n", &dna_orientations[atom_id][j].x, &dna_orientations[atom_id][j].y, &dna_orientations[atom_id][j].z);
                }
                atom_id++;
            }
        }
        for (int i = 0; i < num_cores; i++) {
            fprintf(trajFile, "%e %e %e\n", core_positions[i].x, core_positions[i].y, core_positions[i].z);
            for (int j = 0; j < 3; j++) {
                fprintf(trajFile, "%e %e %e\n", orientations[i][j].x, orientations[i][j].y, orientations[i][j].z);
            }
        }
        
        for (int i = 0; i < num_cores*50; i++) {
            fscanf(fortFile, "%lf %lf %lf\n", &tail_positions[i].x, &tail_positions[i].y, &tail_positions[i].z);
        }
        for (int i = 0; i < num_cores*28; i++) {
            fscanf(fortFile, "%lf %lf %lf\n", &lh_positions[i].x, &lh_positions[i].y, &lh_positions[i].z);
        }
    }
    fclose(fortFile);
    fclose(trajFile);
     
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
