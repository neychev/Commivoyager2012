#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
int v_quant = 0;
int* fork = (int*)malloc(sizeof(int));
int f_q = 0;
typedef struct edge{
        int s;
        int t;
    } edge;

edge *eq_w = (edge*)malloc(sizeof(edge));

typedef struct vertex{
    edge* w;//ребра в пути
    edge* b;//плохие ребра
    int l;//кол-во ребер в пути
    int q;//кол-во плохих ребер
    int L;//нижняя грань длины пути
    char been;
} vertex;
vertex* V = (vertex*)malloc(sizeof(vertex));


int get_length(char* filename){
    freopen(filename, "r", stdin);
    int n = 0;
    unsigned int i;
    char s[4096];
    gets(s);
    while(strlen(s) == 0)
        gets(s);
    int len = strlen(s);
    for(i = 0; i < strlen(s);i++)
        if (((s[i] == ' ')||(s[i] == '\t'))&&((s[i+1] != ' ')&&(s[i+1]!= '\t')&&(s[i+1]!=0))) 
            n++;
    if (s[0] != ' ') n++;
    return n;
}

void output_(int a, int **M, int n){
    int i,j;
    printf("after %d\n", a);
    for(j = 0; j < n; j++){
        for(i = 0; i < n; i++)
            printf("%d\t", M[j][i]);
    printf("\n");
    }
}
void input_(char *filename, int **M, int n){
    freopen(filename, "r", stdin);
    int i, k, j;
    char l[15];
    for (k = 0; k < n; k++)
        for(i = 0; i < n; i++){
            scanf("%s", l);
            M[k][i] = atoi(l);
            if (M[k][i] == 0){
                j = 0;
                while(j < 15)
                    if (l[j] == '-') {M[k][i] = -1;break;}
                    else j++;
            }
        }    
}


int min_(int *mas, int l){
    int i = 0;
    int flag = 0;
    while(mas[i] == -1)
        i++;
    int m;
    if (i < l) m = mas[i];
    else return -1;
    for(; i < l; i++)
        if(mas[i] != -1){
            flag = 1;
            if(mas[i] < m) m = mas[i];
            else;
        }
        if (flag == 0) m = -1;     
   return m;
}    

int check_data(int **M, int n){
    int i;
    for(i = 0; i < n; i++)
        if(M[i][i] != -1) return 1;
    return 0;
}

int priv_s(int *S, int **M, int n){
    int k = 0;
    int i,j;
    for (i = 0; i<n;i++)
        S[i] = 0;
    for (i = 0; i<n;i++){
        k = min_(M[i], n);
        if (k == -1) {S[i] = 0; continue;}
        S[i] = k;
        for (j = 0; j < n; j++)
                if(M[i][j] != -1)M[i][j] -= k;
    }
    return 0;
}

int priv_r(int *R, int **M, int n){
    int* tM = (int*)malloc(n*sizeof(int));
    int k = 0;
    int i,j;
    for(i = 0; i < n; i++)
        R[i] = 0;
    for(i = 0; i < n; i++){
        for(j = 0; j<n; j++)
            tM[j] = M[j][i];
        k = min_(tM,n);
        if (k == -1) {R[i] = 0; continue;}
        R[i] = k;
        for(j = 0; j<n; j++)
            if(M[j][i] != -1)M[j][i] -= k;
    }
    free(tM);
    return 0;
}

int inc_H(int *R, int *S, int n){
    int i;
    int tmp = 0;
    for (i = 0; i < n; i++){
        tmp += R[i];
        tmp += S[i];
    }
    return tmp;
}

int matrix_priv(int **M, int *R, int *S, int n){
    int tmp;
    priv_s(S,M,n);
    priv_r(R,M,n);
    tmp = inc_H(R,S,n);
    return tmp;
}

void matrix_red(int **M, int str, int row, int n){
    int i,j;
    for (i = 0; i < n; i++)
        M[i][row] = -1;
    for (j = 0; j < n; j++)
        M[str][j] = -1;
}


void matrix_copy(int **M, int **D, int n){
    int i,j;
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            D[i][j] = M[i][j];
}

int best_match(int **D, int *R, int *S, int n, int **M){
     int k = 1;
     int i,j, tmp;
     int L = -1;
     for (i = 0; i<n; i++){
         S[i] = 0;
         R[i] = 0;
     }
     for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if (D[i][j] == 0){
                D[i][j] = -1;
                priv_s(S, D, n);
                priv_r(R, D, n);
                tmp = inc_H(R,S,n);

                if (tmp > L){
                    k = 1;
                    eq_w = (edge*)realloc(eq_w, k*sizeof(edge));
                    
                    L = tmp;
                    eq_w[0].s = i;
                    eq_w[0].t = j;

                }
                else if(tmp == L){
                    k++;
                    eq_w = (edge*)realloc(eq_w, k*sizeof(edge));
                    eq_w[k-1].s = i;
                    eq_w[k-1].t = j;//набираем массив равноценных путей :)
                }
            }
            matrix_copy(M, D, n);
        }
     }
     return L;
}

int rec_search(edge *ways, int l, int tmpb, int tmpd, int *ret, int n){
    int i;
    int begin = tmpb;
    int end = tmpd;
    int flag = 0;
    //printf("1 %d %d\n", begin, end);
    for (i = 0; i < l; i++)
        if (ways[i].t == begin){
            begin = ways[i].s;
            if (begin == end) {printf("FATAL ERROR - INVALID CYCLE\n"); exit(EXIT_FAILURE);}
            flag = 1;
        }
    for (i = 0; i < l; i++)
        if (ways[i].s == end){
            end = ways[i].t;
            if (begin == end) {printf("FATAL ERROR - INVALID CYCLE\n"); exit(EXIT_FAILURE);}
            flag = 1;
        }
    
    //printf("%d %d\n", begin, end);
    //if (begin == end) {printf("FATAL ERROR - INVALID CYCLE\n"); exit(EXIT_FAILURE);}
    ret[0] = begin;
    ret[1] = end;
    if(flag == 1) rec_search(ways, l, ret[0], ret[1], ret, n);//Если в графе есть цикл - зависнет нафиг
    return 0;
}

int if_we_go_0(int **D, int* S, int *R, int n, edge a){
    D[a.t][a.s] = -1;
    matrix_red(D, a.s, a.t, n);
    return matrix_priv(D, R, S, n);
}

int if_we_go(int **D, int *S, int *R, int n, edge *ways, int l, int *bver, edge a, int n_l_s){
    int tmpret[2];
    if (n_l_s == 1) {
        rec_search(ways, l, a.s, a.t, tmpret, n);
        D[tmpret[1]][tmpret[0]] = -1;
    }
    matrix_red(D, a.s, a.t, n);
    bver[0] = tmpret[1];
    bver[1] = tmpret[0];
    return matrix_priv(D, R, S, n);
}

void copy_from_parent(int cur_pos, int p_pos){
    int i;
    int l = V[p_pos].l;
    int q = V[p_pos].q;
    for(i = 0; i < l; i++)
        V[cur_pos].w[i] = V[p_pos].w[i];
    for(i = 0; i < q; i++)
        V[cur_pos].b[i] = V[p_pos].b[i];
    V[cur_pos].L = V[p_pos].L;
    V[cur_pos].been = 0;
}

void block_way(int **M, int pos){
    int a,b;
    a = V[pos].b[V[pos].q - 1].s;
    b = V[pos].b[V[pos].q - 1].t;
    M[a][b] = -1;
}


int init_root(int **M, int **D, int *S, int *R, int n, int *cur_pos){
    int i;
    int H = 0;
    matrix_copy(M,D,n);
    int tmpH = best_match(D, S, R, n, M);
    if(tmpH != -1) H+=tmpH;
    else return -1;
    v_quant += 2;
    V = (vertex*)realloc(V, v_quant*sizeof(vertex));
    
        V[0].b = (edge*)malloc(sizeof(edge));
        V[0].w = (edge*)malloc(sizeof(edge));
        V[0].b[0] = eq_w[0];
        V[0].L = H;
        V[0].l = 0;
        V[0].q = 1;
        V[1].b = (edge*)malloc(sizeof(edge));
        V[1].w = (edge*)malloc(sizeof(edge));
        V[1].w[0] = eq_w[0];
        V[1].b[0].s = eq_w[0].t;
        V[1].b[0].t = eq_w[0].s;
        V[1].l = 1;
        V[1].q = 1;
        matrix_copy(M,D,n);
        int tmpL = if_we_go_0(D, S, R, n, eq_w[0]);
        if (tmpL != -1)
            V[1].L = tmpL;
        else return -1;
    for(i = 0; i < v_quant; i++)
            V[i].been = 0;
    int min_0 = V[0].L;
    int pos = 0;
    if (V[1].L < min_0){
            min_0 = V[1].L;
            pos = 1;
    }
    cur_pos[0] = pos;
    if(V[pos].L == V[pos+1].L){
            //printf("forking!\n");
            f_q++;
            fork = (int*)realloc(fork, f_q*sizeof(int));
            fork[f_q -1] = pos + 1;
        }   
    if (V[pos].l > 0)
        matrix_red(M, V[pos].w[V[pos].l-1].s, V[pos].w[V[pos].l-1].t, n);
    block_way(M, pos);
    int checker = matrix_priv(M, R, S, n);
    int quq = V[pos].L;
    if (checker != quq) printf("\terror type 2\n");
    V[pos].been = 1;
    return 0;
}

void get_to_pos(int **save, int **M, int *R, int *S, int n, int pos, int *cur_pos){
    int i;
    matrix_copy(save, M, n);
    for(i = 0; i < V[pos].q; i++)
        M[V[pos].b[i].s][V[pos].b[i].t] = -1;
    for( i = 0; i < V[pos].l; i++)
        matrix_red(M, V[pos].w[i].s, V[pos].w[i].t, n);
    matrix_priv(M, R, S, n);
    cur_pos[0] = pos;
}



int rec_work(int **M, int **D, int *S, int *R, int n, int *cur_p){
    int cur_pos = cur_p[0];
    matrix_copy(M,D,n);
    int H = 0;
    int n_l_s = 0;//not last step flag
    int tmpH = best_match(D, S, R, n, M);
    if(tmpH != -1) H+=tmpH;
    else return -1;
    V = (vertex*)realloc(V, (v_quant+2)*sizeof(vertex));
    if(V[cur_pos].l < (n - 2)) n_l_s = 1;
        V[v_quant].b = (edge*)malloc((V[cur_pos].q + 1)*sizeof(edge));//не идем по пути 
        V[v_quant].w = (edge*)malloc((V[cur_pos].l)*sizeof(edge));
        copy_from_parent(v_quant, cur_pos);
        V[v_quant].L += H;
        V[v_quant].b[V[cur_pos].q] = eq_w[0];
        V[v_quant].l = V[cur_pos].l;
        V[v_quant].q = V[cur_pos].q + 1;
        V[v_quant+1].b = (edge*)malloc((V[cur_pos].q + n_l_s)*sizeof(edge));//или идем
        V[v_quant+1].w = (edge*)malloc((V[cur_pos].l + 1)*sizeof(edge));
        copy_from_parent(v_quant+1, cur_pos);
        matrix_copy(M,D,n);
        int bver[2];
        int tmpL = if_we_go(D, S, R, n, V[cur_pos].w, V[cur_pos].l, bver, eq_w[0], n_l_s);
        if (tmpL != -1)
            V[v_quant+1].L += tmpL;
        else return -1;
        V[v_quant+1].w[V[cur_pos].l] = eq_w[0];
        V[v_quant+1].l = V[cur_pos].l + 1;
        V[v_quant+1].q = V[cur_pos].q + n_l_s;
        if(n_l_s == 1){
            V[v_quant+1].b[V[cur_pos].q].s = bver[0];
            V[v_quant+1].b[V[cur_pos].q].t = bver[1];
        }
    int pos = v_quant;
    if(V[cur_pos].l >= (n - 2)){
        pos = v_quant + 1;
    }
    else{
        int min_0 = V[v_quant].L;
        pos = v_quant;
        if (V[v_quant+1].L < min_0){
            min_0 = V[v_quant+1].L;
            pos = v_quant + 1;
        }
        if((V[pos].l == V[cur_pos].l)&&(V[pos].L == V[pos+1].L)){
            //printf("\t\t\tforking!\n");
            f_q++;
            fork = (int*)realloc(fork, f_q*sizeof(int));
            fork[f_q -1] = pos + 1;
        }   
    }
    if (V[pos].l > V[cur_pos].l)
        matrix_red(M, V[pos].w[V[pos].l-1].s, V[pos].w[V[pos].l-1].t, n);
    block_way(M, pos);
    cur_p[0] = pos;
    int checker = matrix_priv(M, R, S, n);
    int quq = V[pos].L - V[cur_pos].L;
    if (checker != quq) printf("\t error type 2\n");
    V[pos].been = 1;
    v_quant += 2;
    return 0;
}

int action(int **save, int **M, int **D, int *R, int *S, int n, int *cur_pos, int frH, int mode){
    int al = 0;
    while(1){
        al = rec_work(M, D, R, S, n, cur_pos);//Ищем лучшее ребро - идем или не идем
        if (mode == 1)
            if (V[cur_pos[0]].L >= frH)
                return -1;
        if (V[cur_pos[0]].l == n) break;
        if(al == -1){
            if(f_q > 0){
                get_to_pos(save, M, R, S, n, fork[f_q - 1], cur_pos);
                f_q --;
                //fork = (int*)realloc(fork, f_q*sizeof(int));
                al = 0;
            }
            else return -1;
        }
    }
    return 0;
}

void make_cycle(int pos){
    int i, j;
    int c = V[pos].l;
    edge *path = (edge*)malloc(c*sizeof(edge));
    edge temp;
    for(i = 0; i < c; i++)
        path[i] = V[pos].w[i];
    for(i = 0; i < c - 1; i++){
        if (path[i+1].s != path[i].t)
            for(j = i+2; j < c; j++){
                if (path[j].s == path[i].t){
                    temp = path[j];
                    path[j] = path[i + 1];
                    path[i + 1] = temp;
                }
            }
    }
    for (i = 0; i < c ; i++){
        printf("%d->", path[i].s+1);
        if(i == c-1) printf("%d\n", path[i].t+1);
    }
}

void free_all(int **save, int **M, int **D, int *R, int *S, int n){
    int i;
    for(i = 0; i < n; i++){
        free(save[i]);
        free(M[i]);
        free(D[i]);
    }
    free(save);
    free(M);
    free(D);
    free(S);
    free(R);
}

int get_first_record(int **save, int **M, int **D, int *R, int *S, int n, int *cur_pos, int *H){
    int al;
    int h = 0;
    al = matrix_priv(save, R, S, n);      //Привели матрицу - сумма констант приведения идет в H
    if(al!=-1)
        h+=al;
    else {
        printf("no cycle\n");
        free_all(save, M, D, R, S, n);
        exit(1);
    }
    matrix_copy(save, M, n);
    if (n == 2) {
        if((M[0][1]!=-1)&&(M[1][0]!=-1)){
            printf("%d\n1->2->1\n", h);
            free_all(save, M, D, R, S, n);
            exit(0);
        }
        else{
            printf("no cycle\n");
            free_all(save, M, D, R, S, n);
            exit(1);
        }
    }

    //printf("h = %d\n", h);
    //output_(0, M, n);
    //output_(1, M, n);
    al = init_root(M, D, S, R, n, cur_pos);// Инициализация корней дерева - порождаем две вершины, одна значит идем, другая - не идем
    //output_(2, M, n);
    if (al == -1) {printf("no cycle\n"); free_all(save, M, D, R, S, n); return -1;}
    //output_(0, M, n);
    //printf("position %d\n", cur_pos[0]);
    al = action(save, M, D, S, R, n, cur_pos, 0, 0);// Основное действо - из текущей вершины ищем рекорд
    if (al == -1) {printf("no cycle\n"); free_all(save, M, D, R, S, n); return -1;}
    int record = cur_pos[0];
    /*for(i = 0; i < v_quant; i++)
        if(V[i].been == 0)printf("not been - %d\n", i);*/
    //printf("first record H = %d, place %d, incL = %d\n", h + V[record].L, record, V[record].L);
    H[0] = h;
    return record;
}
    
int backtracking(int **save, int **M, int **D, int *R, int *S, int n, int *rec, int *H){
    int i;
    int pos = v_quant;
    int record = rec[0];
    int al;
    int cur_pos[1];
    int prev_quant = v_quant;
    while(1){
        pos--;
        if (pos < 0)
            break;
        prev_quant = v_quant;
        if((V[pos].L < V[record].L)&&(V[pos].been == 0)){
            V[pos].been = 1;
            get_to_pos(save, M, R, S, n, pos, cur_pos);
            al = action(save, M, D, S, R, n, cur_pos, V[record].L, 1);
            if (al == 0){
                //printf("%d - new record!\n cur_pos %d, v_quant %d, incL = %d\n", V[cur_pos[0]].L + H[0], cur_pos[0], v_quant, V[cur_pos[0]].L);
                record = cur_pos[0];
            }
            if (prev_quant != v_quant) pos = v_quant;
        }
    }
    for(i = 0; i < v_quant; i++)
        if((V[i].L < V[record].L)&&(V[i].been == 0))printf("\t\t\terror type 3%d\n", i);
    rec[0] = record;
    return 0;
}

void check_sum(int **save, int record, int *H){
    int i;
    int Y = 0;
    for(i = 0; i < V[record].l; i++)
        Y += save[V[record].w[i].s][V[record].w[i].t];
    
    printf("%d %d\n", H + V[record].L, Y + H[0]);
}

    
int get_data_from_within(int **M, char *path){
    freopen(path, "r", stdin);
    int n;
    int i;
    scanf("%d", &n);
    char l[16];
    int* tmp =(int*)malloc(2*n*sizeof(int));
    edge *way = (edge*)malloc(n * sizeof(edge));
    for(i = 0; i <= n; i++){
        scanf("%s", l);
        tmp[i] = atoi(l);
    }
    for(i = 0; i < n; i++){
        way[i].s = tmp[i]-1;
        way[i].t = tmp[i+1]-1;
    }
    int SUM = 0;
    for(i = 0; i < n; i++){
      //  printf("%d %d %d\n", way[i].s, way[i].t, M[way[i].s][way[i].t]);
        SUM += M[way[i].s][way[i].t];
    }
    
    printf("sum = %d\n", SUM);
    return n;
}
        

        

        