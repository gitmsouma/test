#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// プロトタイプ宣言は略

void multiAx(double *y, double **a, double *x, int m, int n)
{
    // y = Ax の計算 A in R^{m×n}, x in R^{n}, y in R^{m}
    int i, j;
    for (i = 0; i < m; i++)
    {
        y[i] = 0.0;
        for (j = 0; j < n; j++)
        {
            y[i] = y[i] + a[i][j] * x[j];
        }
    }
}

void multiAB(double **c, double **a, double **b, int m, int l, int n)
{
    // C = AB の計算 A \in R^{m×l}, B \in R^{l×n}, C \in R^{m×n}
    int i, j, k;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            c[i][j] = 0.0;
            for (k = 0; k < l; k++)
            {
                c[i][j] = c[i][j] + a[i][k] * b[k][j];
            }
        }
    }
}

double **matrix(int m, int n)
{ // m X n 行列の行列確保
    int i;
    double **x;
    x = malloc(m * sizeof(double *));
    if (x == NULL)
    {
        printf("allocation error 1 in matrix() \n");
        exit(-1);
    }
    for (i = 0; i < m; i++)
    {
        x[i] = malloc(n * sizeof(double));
        if (x[i] == NULL)
        {
            printf("allocation error 2 in matrix() \n");
            exit(-1);
        }
    }
    return (x);
}

double *vector(int m)
{ // m成分のベクトル
    double *x;
    x = malloc(m * sizeof(double));
    if (x == NULL)
    {
        printf("allocation error in vector()\n");
        exit(-1);
    }
    return (x);
}

void freemat(double **x, int m)
{ //行列の解放
    int i;
    for (i = 0; i < m; i++)
    {
        free(x[i]);
    }
    free(x);
}

double GS_SDL(double **a, double *b, double *x, int size)
{
    double *buf;
    int i, j, k;
    double error, z;

    buf = vector(size);

    for (i = 0; i < 1000; i++)
    {
        printf("ite=%d ", i);

        for (j = 0; j < size; j++)
        {
            z = 0;
            for (k = 0; k < size; k++)
            {
                if (k != j)
                {
                    z = z + a[j][k] * x[k]; // z = Aのj行目（A[j][j]以外）× x
                }
            }
            x[j] = (b[j] - z) / a[j][j];
        }

        // error チェック error = \| A x -b \|^2
        multiAx(buf, a, x, size, size); // buf = A * x

        for (error = 0.0, j = 0; j < size; j++)
        {
            error = error + (buf[j] - b[j]) * (buf[j] - b[j]);
        }

        if (error < 1.0E-6)
        {
            free(buf);
            printf(" v=(");
            for (j = 0; j < size; j++)
            {
                printf("%.6f ", x[j]);
            }
            printf(") error= %.6f \n", error);
            break;
        }
        else
        {
            printf(" v=(");
            for (j = 0; j < size; j++)
            {
                printf("%.6f ", x[j]);
            }
            printf(") error= %.6f \n", error);
        }
    }
    return error;
}

void printWb(double *v, int dim, double error)
{ // w, theta, error の表示
    int i;

    printf("w=("); // vの1～dim番目の要素(=w)を表示
    for (i = 0; i < dim; i++)
    {
        printf("%.6lf ", v[i]);
    }
    printf(") ");
    printf("theta=%.6lf error=%.6lf\n", v[dim], error); // (dim+1)番目の要素(=theta)，誤差表示
}

int main(int argc, char *argv[])
{
    double **a, **x, **tx, *b, *v, *y, error;
    int i, k, num, dim; // num データ数 , dim 次元
    FILE *fp;

    // コマンドライン引数として得たデータファイルを開く
    if (NULL == (fp = fopen("data01.txt", "r")))
    {
        exit(-1);
    }

    fscanf(fp, "%d %d", &dim, &num); // 次元，データ数の読み込み
    printf("%d %d\n", dim, num);

    x = matrix(num, dim + 1);     // x の領域確保：num × dim+1
    tx = matrix(dim + 1, num);    // 行列 X の転置行列 X^t の領域確保 dim+1 × num
    a = matrix(dim + 1, dim + 1); // X^t X を格納する 2 次元配列 a の領域確保 dim+1 × dim+1
    y = vector(num);              // y を格納する配列 y の動的確保 num × 1
    b = vector(dim + 1);          // X^t y を格納する配列 b の動的確保 dim+1 × 1
    v = vector(dim + 1);          // 最小2乗法で求める重みw 閾値θ からなるベクトル v の領域確保 dim+1 × 1

    for (k = 0; k < num; k++)
    {
        // X, X^t へのデータ入力
        for (i = 0; i < dim; i++)
        {
            fscanf(fp, "%lf", &x[k][i]); // X のデータ読み込み
            /* 確認debug用` */ printf("%f ", x[k][i]);
            tx[i][k] = x[k][i]; // X から X^t へのデータ入力
        }
        x[k][dim] = 1;            // 行列 x の (k,dim+1) 要素は常に 1 を代入
        tx[dim][k] = 1;           // 対応する tx の (dim+1,k) 要素に 1 を代入
        fscanf(fp, "%lf", &y[k]); // y のデータ入力
        /* 確認debug用 */ printf("%f \n", y[k]);
    }
    multiAB(a, tx, x, dim + 1, num, dim + 1); // A = X^t X の計算
    /* 確認debug用 */ for (i = 0; i < dim + 1; i++)
    {
        for (k = 0; k < num; k++)
        {
            printf("%f ", a[i][k]);
        }
        printf("\n");
    }
    /* 確認debug用 */ printf("\n");
    multiAx(b, tx, y, dim + 1, num); // b = X^t y の計算
    /* 確認debug用 */ for (i = 0; i < dim + 1; i++)
    {
        printf("%f ", b[i]);
    }
    printf("\n");

    srand((unsigned int)time(NULL));
    for (i = 0; i < dim + 1; i++)
    {
        v[i] = GS_SDL(a, b, *x, dim + 1);
    }

    /* 確認debug用 */ printf("\n");
    error = GS_SDL(a, b, *x, dim + 1);
    printWb(v, dim, error); // (w, theta), error を表示

    freemat(x, num);
    freemat(tx, dim + 1);
    freemat(a, dim + 1);
    free(y);
    free(b);
    free(v);

    return 0;
}