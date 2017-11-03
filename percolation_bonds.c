#include<stdio.h>
#include<stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>

int LATTICE_SIZE = 10; //default
int **LATTICE;
int largestCluster = 0;
int highestColumn = 0;
int highestRow = 0;
int row_percolates = 0;
int column_percolates = 0;

void printLattice(int** lattice){
    if(LATTICE_SIZE > 50) return;
    printf("   ");
    for(int i = 0; i < LATTICE_SIZE; i++){
        printf("%2i ", i);
    }
    printf("\n");
    for(int r = 0; r < LATTICE_SIZE; r++){
        printf("%2i ", r);
        for(int c = 0; c < LATTICE_SIZE; c++){
            switch(lattice[c][r]){
                case 1: printf("\x1B[0m"" █ "); break; //occupied
                case 2: printf("\x1B[36m"" █ ""\x1B[0m"); break;   // checked
                case 3: printf("\x1B[31m"" █ ""\x1B[0m");  break; //start
                default: printf("\x1B[0m"" ☐ ");
            }
        }
        printf("\n");
    }
}

// fill up randomly allocated lattice sites

void createSiteLattice(double p_seed){

    for(int c = 0; c < LATTICE_SIZE; c++){

        for(int r = 0; r < LATTICE_SIZE; r++){
            double random = ((double) rand())/((double) RAND_MAX);
            if(random <= p_seed){
                LATTICE[c][r] = 1; // occupied
            }
            else{
                LATTICE[c][r] = 0; // not occupied
            }
        }
    }
}

int checkSiteRowPerc(int r, int c, int **lattice_check){
    int clusterSize = 0;
    if(lattice_check[c][r] == 0 ) return clusterSize; // empty site
    if(lattice_check[c][r] == 3){
        //printf("Checked 3: [%i][%i]\n", r, c);
        clusterSize++;
    }
    if(lattice_check[c][r] == 1){ // unchecked, mark as checked
        clusterSize++;
        lattice_check[c][r] = 2; // checked
        //printf("Checked: [%i][%i]]\n", r, c);
    }
    if(lattice_check[c][r] == 2 || lattice_check[c][r] == 3) // checked & start point, check outer sites
    {
        int northRow    = (r-1);
        int southRow    = (r+1);
        int eastColumn  = (c+1)%(LATTICE_SIZE);
        int westColumn  = (c-1+LATTICE_SIZE)%(LATTICE_SIZE);
        
        if(southRow < LATTICE_SIZE && lattice_check[c][southRow] == 1){
            if(southRow > highestRow) highestRow = southRow;
            if(southRow == LATTICE_SIZE - 1) row_percolates = 1;
            clusterSize += checkSiteRowPerc(southRow, c, lattice_check);
        }
        if(northRow > 0 && lattice_check[c][northRow] == 1){
            clusterSize += checkSiteRowPerc(northRow, c, lattice_check);
        }
        if(lattice_check[eastColumn][r] == 1){
            clusterSize += checkSiteRowPerc(r, eastColumn, lattice_check);
        }
        if(lattice_check[westColumn][r] == 1){
            clusterSize += checkSiteRowPerc(r, westColumn, lattice_check);
        }
    }
    return clusterSize;
}

int checkSiteColPerc(int r, int c, int **lattice_check){
    int clusterSize = 0;
    if(lattice_check[c][r] == 0 ) return clusterSize; // empty site
    if(lattice_check[c][r] == 3){
        //printf("Checked 3: [%i][%i]\n", r, c);
        clusterSize++;
    }
    if(lattice_check[c][r] == 1){ // unchecked, mark as checked
        clusterSize++;
        lattice_check[c][r] = 2; // checked
        //printf("Checked: [%i][%i]]\n", r, c);
    }
    if(lattice_check[c][r] == 2 || lattice_check[c][r] == 3) // checked & start point, check outer sites
    {
        int southRow    = (r-1+LATTICE_SIZE)%(LATTICE_SIZE);
        int northRow    = (r+1)%(LATTICE_SIZE);
        int eastColumn  = (c+1);
        int westColumn  = (c-1);
        
        if(lattice_check[c][southRow] == 1){
            clusterSize += checkSiteColPerc(southRow, c, lattice_check);
        }
        if(lattice_check[c][northRow] == 1){
            clusterSize += checkSiteColPerc(northRow, c, lattice_check);
        }
        if(eastColumn < LATTICE_SIZE && lattice_check[eastColumn][r] == 1){
            if(eastColumn > highestColumn) highestColumn = eastColumn;
            if(eastColumn == (LATTICE_SIZE - 1)) column_percolates = 1;
            clusterSize += checkSiteColPerc(r, eastColumn, lattice_check);
        }
        if(westColumn >= 0 && lattice_check[westColumn][r] == 1){
            clusterSize += checkSiteColPerc(r, westColumn, lattice_check);
        }
    }
    return clusterSize;
}

void clearLattice(int** dest){
    for(int r = 0; r < LATTICE_SIZE; r++){
        for(int c = 0; c < LATTICE_SIZE; c++){
            dest[c][r] = LATTICE[c][r];
        }
    }
}

void checkSiteLattice(){
    int **lattice_check;
    int **lattice_final;
    
    lattice_check = (int**) malloc(LATTICE_SIZE * sizeof(int*));
    lattice_final = (int**) malloc(LATTICE_SIZE * sizeof(int*));
    for (int i = 0; i < LATTICE_SIZE; i++) {
        lattice_check[i] = (int*) malloc(LATTICE_SIZE * sizeof(int));
        lattice_final[i] = (int*) malloc(LATTICE_SIZE * sizeof(int));
    }
    
    //copy lattice into lattice_copy
    clearLattice(lattice_check);
    printLattice(lattice_check);
    //check lattice for percolation starting from each site in top row
    largestCluster = checkSiteRowPerc(0, 0, lattice_check);
    lattice_check[0][0] = 3;
    for(int r = 0; r < LATTICE_SIZE; r++){
        for(int c = 0; c < LATTICE_SIZE; c++){
            lattice_final[c][r] = lattice_check[c][r];
        }
    }
    clearLattice(lattice_check);
    
    for(int i = 1; i < LATTICE_SIZE; i++){
        clearLattice(lattice_check);
        if(lattice_check[i][0] == 1) lattice_check[i][0] = 3;
        int clusterSize = checkSiteRowPerc(0, i, lattice_check);
        if(clusterSize > largestCluster){
            for(int r = 0; r < LATTICE_SIZE; r++){
                for(int c = 0; c < LATTICE_SIZE; c++){
                    lattice_final[c][r] = lattice_check[c][r];
                }
            }
            largestCluster = clusterSize;
        }
    }
    for(int i = 0; i < LATTICE_SIZE; i++){
        clearLattice(lattice_check);
        if(lattice_check[0][i] == 1) lattice_check[0][i] = 3;
        int clusterSize = checkSiteColPerc(i, 0, lattice_check);
        if(clusterSize > largestCluster){
            //printLattice(lattice_check);
            for(int r = 0; r < LATTICE_SIZE; r++){
                for(int c = 0; c < LATTICE_SIZE; c++){
                    lattice_final[c][r] = lattice_check[c][r];
                }
            }
            largestCluster = clusterSize;
        }
    }
    printLattice(lattice_final);
    free(lattice_check);
    free(lattice_final);
}

void getBondLattice(double p_seed){
    
    // create bond lattice, horizontal and vertical
    int bond_lattice_size = LATTICE_SIZE;
    int hoz_bond_lattice[bond_lattice_size][bond_lattice_size - 1];
    int ver_bond_lattice[bond_lattice_size - 1][bond_lattice_size];
    
    // allocate occupied and vacant bonds
    for(int c = 0; c < bond_lattice_size; c++){
        for(int r = 0; r < bond_lattice_size - 1; r++){
            // randomly allocate horizontal bonds
            double random1 = ((double) rand())/((double) RAND_MAX);
            if(random1 <= p_seed){
                hoz_bond_lattice[c][r] = 1; // occupied
            }
            else{
                hoz_bond_lattice[c][r] = 0; // not occupied
            }
            
            // randomly allocate vertical bonds
            double random2 = ((double) rand())/((double) RAND_MAX);
            if(random2 <= p_seed){
                ver_bond_lattice[r][c] = 1; // occupied
            }
            else{
                ver_bond_lattice[r][c] = 0; // not occupied
            }
        }
    }
    
    // fill lattice with empty sites
    for(int c = 0; c < LATTICE_SIZE; c++) {
        for (int r = 0; r < LATTICE_SIZE; r++) {
            LATTICE[c][r] = 0; // not occupied
        }
    }
    // process horizontal bonds into site lattice
    for(int c = 0; c < bond_lattice_size - 1; c++){
        for(int r = 0; r < bond_lattice_size ; r++){
            if(hoz_bond_lattice[r][c] == 1){
                LATTICE[r][c] = 1;
                LATTICE[r][c+1] = 1;
            }
        }
    }
    // process vertical bonds into site lattice
    for(int c = 0; c < bond_lattice_size; c++){
        for(int r = 0; r < bond_lattice_size - 1; r++){
            if(hoz_bond_lattice[r][c] == 1){
                LATTICE[r][c] = 1;
                LATTICE[r+1][c] = 1;
            }
        }
    }
    // check new site lattice
}

int main(int argc, char* argv[]){

    double p_seed;
    
    // get start time
    struct timeval start, end;
    gettimeofday(&start, NULL);
    
    if(argc >= 5){
        
        LATTICE_SIZE = atoi(argv[1]);
        
        if(LATTICE_SIZE < 1){
            printf("LATTICE SIZE MUST BE GREATER THAN 0\n");
        }
        p_seed = atof(argv[2]);
        int percolation_type = atoi(argv[4]);
        int is_site_perc;
        if(!strcmp(argv[3], "b") || !strcmp(argv[3], "B")){
            is_site_perc = 0;
        }
        else if(!strcmp(argv[3], "s") || !strcmp(argv[3], "S")) is_site_perc = 1;
        else {
            fprintf(stderr, "Input either 'S' or 'B' for site or bond percolation\n");
            return 0;
        }
        
        //print simulation constants
        printf("Lattice Size: %i\nProbability: %f\n", LATTICE_SIZE, p_seed);
        printf("Type of Lattice: ");
        
        if(is_site_perc == 1) printf("Site\n");
        else printf("Bond\n");
        
        if(percolation_type == 0) printf("Percolation type: Row\n");
        else if(percolation_type == 1) printf("Percolation type: Column\n");
        else if(percolation_type == 2) printf("Percolation type: Row & Column\n");
        
        //dynamically allocate lattice size
        LATTICE = (int **) malloc(LATTICE_SIZE * sizeof(int*));
        for(int i = 0; i < LATTICE_SIZE; i++){
            LATTICE[i] = (int *) malloc(LATTICE_SIZE * sizeof(int));
        }
        
        // site lattice
        if(is_site_perc) {
            createSiteLattice(p_seed);
            // check created lattice for Site Percolation
            checkSiteLattice();
        }
        else {
            getBondLattice(p_seed);
            checkSiteLattice();
        }
        
        //get end time and calculate
        gettimeofday(&end, NULL);
        double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
        printf("Time =%f\n", delta);
        free(LATTICE);
        
        switch(percolation_type){
            case 0: if(row_percolates) printf("Row Percolation: true\n");
            else { printf("Row Percolation: false\n");}
                break;
            case 1:if(column_percolates) printf("Column Percolation: true\n");
            else {printf("Column Percolation: false\n");}
                break;
            case 2: if(column_percolates && row_percolates) printf("Row & Column Percolation: true\n");
            else {
                printf("Row & Column Percolation: false\n");
                if(column_percolates) printf(" Column Percolation: true\n");
                else printf(" Column Percolation: false\n");
                if(row_percolates) printf(" Row Percolation: true\n");
                else printf(" Row Percolation: false\n");
            }
                break;
            default: fprintf(stderr, "Percolation Type Invalid\n");
        }
        printf("largestCluster: %i\n", largestCluster);
    }
    else {
        fprintf(stderr, "USAGE: [grid size] [p of seed] [s (site) or b (bond)] [type of percolation]\n Type of percolation: \n 0 - cluster must span all rows \n 1 - cluster must span all columns \n 2 - cluster must span all rows & columns\n");
        return 0;
    }
    return 0;
}
