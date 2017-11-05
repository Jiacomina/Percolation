#include <mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>

int LATTICE_SIZE = 10; //default
char **SITE_LATTICE;
char **BOND_LATTICE;
int thelargestCluster = 0;
int highestColumn = 0;
int highestRow = 0;
int row_percolates = 0;
int column_percolates = 0;
int num_threads = 1; 
int pid;
int numProcess;

void printLattice(char** lattice){
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
                case '1': printf("\x1B[0m"" █ "); break; //occupied
                case '2': printf("\x1B[36m"" █ ""\x1B[0m"); break;   // checked
                case '3': printf("\x1B[31m"" █ ""\x1B[0m");  break; //start
                default: printf("\x1B[0m"" ☐ ");
            }
        }
        printf("\n");
    }
}
// fill up randomly allocated lattice sites
void createSiteLattice(float p_seed){
    #pragma omp parallel for collapse(2)
        for(int c = 0; c < LATTICE_SIZE; c++){
            
            for(int r = 0; r < LATTICE_SIZE; r++){
                float random = ((float) rand())/((float) RAND_MAX);
                if(random <= p_seed){
                    SITE_LATTICE[c][r] = '1'; // occupied
                }
                else{
                    SITE_LATTICE[c][r] = '0'; // not occupied
                }
            }
        }
}

int checkSiteRowPerc(int r, int c, char **lattice_check){
    int clusterSize = 0;
    if(lattice_check[c][r] == '0' ) return clusterSize; // empty site
    if(lattice_check[c][r] == '3'){
        clusterSize++;
    }
    if(lattice_check[c][r] == '1'){ // unchecked, mark as checked
        clusterSize++;
        lattice_check[c][r] = '2'; // checked
    }
    if(lattice_check[c][r] == '2' || lattice_check[c][r] == '3') // checked & start point, check outer sites
    {
        int northRow    = (r-1);
        int southRow    = (r+1);
        int eastColumn  = (c+1)%(LATTICE_SIZE);
        int westColumn  = (c-1+LATTICE_SIZE)%(LATTICE_SIZE);
        
        if(southRow < LATTICE_SIZE && lattice_check[c][southRow] == '1'){
            if(southRow > highestRow) highestRow = southRow;
            if(southRow == LATTICE_SIZE - 1) row_percolates = 1;
            clusterSize += checkSiteRowPerc(southRow, c, lattice_check);
        }
        if(northRow > 0 && lattice_check[c][northRow] == '1'){
            clusterSize += checkSiteRowPerc(northRow, c, lattice_check);
        }
        if(lattice_check[eastColumn][r] == '1'){
            clusterSize += checkSiteRowPerc(r, eastColumn, lattice_check);
        }
        if(lattice_check[westColumn][r] == '1'){
            clusterSize += checkSiteRowPerc(r, westColumn, lattice_check);
        }
    }
    return clusterSize;
}

int checkSiteColPerc(int r, int c, char **lattice_check){
    int clusterSize = 0;
    if(lattice_check[c][r] == '0' ) return clusterSize; // empty site
    if(lattice_check[c][r] == '3'){
        clusterSize++;
    }
    if(lattice_check[c][r] == '2'){
        return clusterSize;
    }
    if(lattice_check[c][r] == '1'){ // unchecked, mark as checked
        clusterSize++;
        lattice_check[c][r] = '2'; // checked
    }
    if(lattice_check[c][r] == '2' || lattice_check[c][r] == '3') // checked & start point, check outer sites
    {
        int southRow    = (r-1+LATTICE_SIZE)%(LATTICE_SIZE);
        int northRow    = (r+1)%(LATTICE_SIZE);
        int eastColumn  = (c+1);
        int westColumn  = (c-1);
        
        if(lattice_check[c][southRow] == '1'){
            clusterSize += checkSiteColPerc(southRow, c, lattice_check);
        }
        if(lattice_check[c][northRow] == '1'){
            clusterSize += checkSiteColPerc(northRow, c, lattice_check);
        }
        if(eastColumn < LATTICE_SIZE && lattice_check[eastColumn][r] == '1'){
            if(eastColumn > highestColumn) highestColumn = eastColumn;
            if(eastColumn == (LATTICE_SIZE - 1)) column_percolates = 1;
            clusterSize += checkSiteColPerc(r, eastColumn, lattice_check);
        }
        if(westColumn >= 0 && lattice_check[westColumn][r] == '1'){
            clusterSize += checkSiteColPerc(r, westColumn, lattice_check);
        }
    }
    return clusterSize;
}

void clearLattice(char** dest){
    #pragma omp parallel
    {
        #pragma omp for
        for(int r = 0; r < LATTICE_SIZE; r++){
            for(int c = 0; c < LATTICE_SIZE; c++){
                dest[c][r] = SITE_LATTICE[c][r];
            }
        }
    }
}

void checkSiteLattice(){
    char **lattice_check;
    int largestCluster = 0;
    lattice_check = (char**) malloc(LATTICE_SIZE * sizeof(char*));
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < LATTICE_SIZE; i++) {
            lattice_check[i] = (char*) malloc(LATTICE_SIZE * sizeof(char));
        }
    }
    //copy lattice into lattice_copy
    clearLattice(lattice_check);
    printLattice(lattice_check);
    //check lattice for percolation starting from each site in top row
    largestCluster = checkSiteRowPerc(0, 0, lattice_check);
    lattice_check[0][0] = 3;
    
    clearLattice(lattice_check);

    for(int i = 1; i < LATTICE_SIZE; i++){
        if((i + numProcess) % numProcess == pid){
            clearLattice(lattice_check);
            if(lattice_check[i][0] == 1) lattice_check[i][0] = 3;
            int clusterSize = checkSiteRowPerc(0, i, lattice_check);
            if(clusterSize > largestCluster){
                largestCluster = clusterSize;
            }
        }
        
    }

    for(int i = 0; i < LATTICE_SIZE; i++){
        if((i + numProcess) % numProcess == pid){
            clearLattice(lattice_check);
            if(lattice_check[0][i] == 1) lattice_check[0][i] = 3;
            int clusterSize = checkSiteColPerc(i, 0, lattice_check);
            if(clusterSize > largestCluster){                
                largestCluster = clusterSize;
            }
        }
    }
    
    free(lattice_check);

    printf("%i/%i: Largest Lattice %i\n", pid, numProcess, largestCluster);

    MPI_Reduce(&largestCluster, &thelargestCluster, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
}

void getBondLattice(){
    BOND_LATTICE = (char**) malloc(LATTICE_SIZE * sizeof(char*));
    for (int i = 0; i < LATTICE_SIZE; i++) {
        BOND_LATTICE[i] = (char*) malloc(LATTICE_SIZE * sizeof(char));
    }
}
int checkBond(int r, int c, float p_seed){
    int clusterSize = 0;
    short int next_row = (r+1) % LATTICE_SIZE;
    
    // check bottom bond, goes to next row
    if(((float) rand())/((float) RAND_MAX) <= p_seed){
        //BOND_LATTICE[r][c].bottom = '1'; // bond exists, check connecting lattice
        if(BOND_LATTICE[next_row][c] != '1'){
            if(!row_percolates && next_row == LATTICE_SIZE - 1) row_percolates = 1;
            clusterSize++;
            clusterSize += checkBond(next_row, c, p_seed);
        }
    }
    else{
    }
    
    // check right bond, goes to next column
    if(((float) rand())/((float) RAND_MAX) <= p_seed){
        short int next_col = (c+1) % LATTICE_SIZE;
        if(BOND_LATTICE[r][next_col] != '1'){
            if(!column_percolates && next_col == LATTICE_SIZE - 1) column_percolates = 1;
            clusterSize++;
            clusterSize += checkBond(r, next_col, p_seed);
        }
    }
    else{
    }
    BOND_LATTICE[r][c] = '1';
    return clusterSize;
}
void checkBondLattice(float p_seed){
    int largestCluster;
    // divide lattice into parts
    // int num_seg = numProcess;
    // int sub_lattice_size = LATTICE_SIZE / num_seg;
    
    // for(int i = 0; i < num_seg; i++){
    //     int start = i*sub_lattice_size;
    //     int end = (i+1)*sub_lattice_size;
    //     for(int c = start; c < end; c++)
    //     {
    //         int clusterSize = 0;
    //         clusterSize = checkBond(0, c, p_seed);
    //         if(clusterSize > largestCluster) largestCluster = clusterSize;
    //     }
    // }

    for(int i = 1; i < LATTICE_SIZE; i++){
        if((i + numProcess) % numProcess == pid){
            clusterSize = checkBond(0, i, p_seed);
            if(clusterSize > largestCluster){
                largestCluster = clusterSize;
            }
        }
    }

    for(int i = 1; i < LATTICE_SIZE; i++){
        if((i + numProcess) % numProcess == pid){
            clusterSize = checkBond(i, 0, p_seed);
            if(clusterSize > largestCluster){
                largestCluster = clusterSize;
            }
        }
    }

    MPI_Reduce(&largestCluster, &thelargestCluster, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    // // starting from each row of leftmost column (excluding top left site)
    // for(int i = 0; i < num_seg; i++){
    //     int start = i*sub_lattice_size;
    //     int end = (i+1)*sub_lattice_size;
    //     for(int r = start; r < end; r++)
    //     {
    //         int clusterSize = 0;
    //         clusterSize = checkBond(r, 0, p_seed);
    //         if(clusterSize > largestCluster[threadNum]) largestCluster = clusterSize;
    //     }
    // }

}

int main(int argc, char* argv[]){
    
    // OPEN MPI
    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &pid); // reports number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &numProcess); // reports the rank, a number between 0 and size-1 identifying the calling 

    float p_seed;
    srand(10);
    // get start time
    struct timeval start, end;
    gettimeofday(&start, NULL);
    int percolation_type = atoi(argv[4]);
    if(argc >= 5){

        LATTICE_SIZE = atoi(argv[1]);
            
        if(LATTICE_SIZE < 1){
            if(pid == 0) printf("LATTICE SIZE MUST BE GREATER THAN 0\n");
        }
        p_seed = atof(argv[2]);
        int is_site_perc;
        if(!strcmp(argv[3], "b") || !strcmp(argv[3], "B")){
            is_site_perc = 0;
        }
        else if(!strcmp(argv[3], "s") || !strcmp(argv[3], "S")){
            is_site_perc = 1;
        }
        else {
            if(pid == 0)
                fprintf(stderr, "Input either 'S' or 'B' for site or bond percolation\n");
            return 0;
        }

        if(pid == 0){

            //print simulation constants
            printf("Lattice Size: %i\nProbability: %f\n", LATTICE_SIZE, p_seed);
            printf("Type of Lattice: ");
            
            if(is_site_perc == 1) printf("Site\n");
            else printf("Bond\n");
            
            if(percolation_type == 0) printf("Percolation type: Row\n");
            else if(percolation_type == 1) printf("Percolation type: Column\n");
            else if(percolation_type == 2) printf("Percolation type: Row & Column\n");

        }

        num_threads = atoi(argv[5]);
        omp_set_num_threads(num_threads); 
        
        MPI_Barrier(MPI_COMM_WORLD);
        // site lattice
        if(is_site_perc){
            //dynamically allocate lattice size
            SITE_LATTICE = (char **) malloc(LATTICE_SIZE * sizeof(char*));
            for(int i = 0; i < LATTICE_SIZE; i++){
                SITE_LATTICE[i] = (char *) malloc(LATTICE_SIZE * sizeof(char));
            }
            createSiteLattice(p_seed);
            // check created lattice for Site Percolation
            checkSiteLattice();
            free(SITE_LATTICE);
        }
        else {
            getBondLattice();
            //checkBondLattice(p_seed);
            free(BOND_LATTICE);
        }
        
        if(pid == 0){
            //get end time and calculate
            gettimeofday(&end, NULL);
            double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
            printf("Time =%f\n", delta);
            
            switch(percolation_type){
                case 0: if(highestRow) printf("Row Percolation: true\n");
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
                MPI_Barrier(MPI_COMM_WORLD);
                printf("Largest Cluster: %i\n", thelargestCluster);
        }
    }
    else {
        if(pid == 0){
            fprintf(stderr, "USAGE: [grid size] [p of seed] [s (site) or b (bond)] [type of percolation]\n Type of percolation: \n 0 - cluster must span all rows \n 1 - cluster must span all columns \n 2 - cluster must span all rows & columns\n");
        }
        return 0;
    }
    MPI_Finalize();
    return 0;
}
