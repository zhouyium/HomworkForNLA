/*
*功能: 列选主元消元法 
*@Language: C++
*@File Name: Gaussian_wp.cpp
*编码：GBK
*/
#include <bits/stdc++.h>

using namespace std;

//#define __OUTPUT

const int MAXN = 1e2+2;
double matrix[MAXN][MAXN];// augmented matrix

int row;//the row of augmented matrix
int col;//the column of augmented matrix

//test output
void printM() {
    for (int i=1; i<=row; i++) {
        for (int j=1; j<=col; j++) {
            printf("%12.4lf", matrix[i][j]);
        }
        printf("\n");
    }
}

//select pivoting
void SelectColE() {
    double temp;
    for (int i=1; i<row; i++) {
        int r = i;
        for (int j=i+1; j<=col; j++) {
            if (fabs(matrix[j][i]) > fabs(matrix[r][i])) {
                r = j;
            }
        }
        if (r != i) {
            for (int j=i; j<=col; j++) {
                swap(matrix[i][j], matrix[r][j]);//exchange
            }
        }
        for (int j=i+1; j<=row; j++) {//Gaussian elimination
            temp=matrix[j][i]/matrix[i][i];
            for (int k=i; k<=col; k++) {
                matrix[j][k] -= matrix[i][k]*temp;
            }
        }

#if defined(__OUTPUT)
        printf("%d Gaussian elimination：\n", i);
        printM();
#endif
    }
}

//Gaussian elimination with column pivoting
void Gauss() {
    SelectColE();
    printf("Upper triangular matrix：\n");
    printM();
    for (int i=row; i>=1; i--) {
        for (int j=i+1; j<=row; j++) {
            matrix[i][col] -= matrix[i][j] * matrix[j][col];
        }
        matrix[i][col] /= matrix[i][i];
    }
}

int main() {
    //freopen("gaussian.out", "w", stdout);

#if 1
    //定义作业的行列
    row=84;
    col=85;

    //第1行
    matrix[1][1] = 6;
    matrix[1][2] = 1;
    matrix[1][85] = 7;
    //第84行
    matrix[84][83] = 8;
    matrix[84][84] = 6;
    matrix[84][85] = 14;
    int t=1;
    for (int i=2; i<84; i++) {
        matrix[i][t+0] = 8;
        matrix[i][t+1] = 6;
        matrix[i][t+2] = 1;
        matrix[i][85]  = 15;
        t++;
    }
#else
    cin>>row;
    col=row+1;
    for (int i=1; i<=row; i++) {
        for (int j=1; j<=col; j++) {
            scanf("%lf", &matrix[i][j]);
        }
    }
#endif

    Gauss();
    
    printf("The answer is:\n");
    for (int i=1; i<=row; i++) {
        printf("X%d = %14.10lf\n", i, matrix[i][col]);
    }
    system("pause");

    return 0;
}