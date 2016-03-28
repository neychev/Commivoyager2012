#include <stdio.h>
#include <malloc.h>
#include "head.h"

int main(int argc, char** argv){
//int main(){
    int i;
    if(argc < 3) {printf("wrong input of arguments!\n format: input_file output_file\n"); return 2;}
    freopen(argv[2], "w", stdout);
    int n = get_length(argv[1]);
    //int n = get_length("input.txt");
    int **M = (int**)malloc(n*sizeof(int*));
    int **D = (int**)malloc(n*sizeof(int*));
    int **save = (int**)malloc(n*sizeof(int*));
    for(i = 0; i < n; i++){
        M[i] = (int*)malloc(n*sizeof(int));
        D[i] = (int*)malloc(n*sizeof(int));
        save[i] = (int*)malloc(n*sizeof(int));
    }
    int* R = (int*)malloc(n*sizeof(int));
    int* S = (int*)malloc(n*sizeof(int));
    
    input_(argv[1], save, n);
    //input_("input.txt", save, n);
    int al = check_data(save, n);
    //output_(1,save,n);
    if(al == 1) {printf("wrong data input\n"); free_all(save, M, D, R, S, n); return 2;}
    int H[1];
    //get_data_from_within(save, "output.txt");
    int cur_pos[1];
    int record[1];
    record[0] = get_first_record(save, M, D, R, S, n, cur_pos, H);  
    if(record[0] == -1) {
        return 1;
    }
    backtracking(save, M, D, R, S, n, record, H);
    //printf("Got here!\n");
    printf("%d\n", H[0] + V[record[0]].L);
    make_cycle(record[0]);
    free_all(save, M, D, R, S, n);
    return 0;
}

//Где-то не увеличивает длину пути или делает это неверно :(
//Замечено только при использовании backtracking'а

//в тесте с единицами - не работает вывод 0_о