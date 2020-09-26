/*
*功能: 列选主元消元法 
*@Language: C++
*@File Name: Gaussian_wp.cpp
*编码：GBK
*/
#include <bits/stdc++.h>

using namespace std;

#define __OUTPUT

const int MAXN = 6;//1e2+2;//增广矩阵最大值
double matrix[MAXN][MAXN];//增广矩阵，注意最后一列答案

int row;//增广矩阵的实际行
int col;//增广矩阵的实际列

//输出矩阵
void printM() {
    for (int i=1; i<=row; i++) {
        for (int j=1; j<=col; j++) {
            printf("%12.4lf", matrix[i][j]);
        }
        printf("\n");
    }
}

//选择列主元并进行消元
void SelectColE() {
    double temp; //用于记录消元时的因数
    for (int i=1; i<row; i++) {
        int r = i;
        for (int j=i+1; j<=col; j++) {
            if (fabs(matrix[j][i]) > fabs(matrix[r][i])) {
                r = j;
            }
        }
        if (r != i) {
            for (int j=i; j<=col; j++) {
                swap(matrix[i][j], matrix[r][j]);//与最大主元所在行交换
            }
        }
        for (int j=i+1; j<=row; j++) {//消元
            temp=matrix[j][i]/matrix[i][i];
            for (int k=i; k<=col; k++) {
                matrix[j][k] -= matrix[i][k]*temp;
            }
        }

#if defined(__OUTPUT)
        printf("第%d列消元后：\n", i);
        printM();
#endif
    }
}

//高斯消元法(列选主元)
void Gauss() {
    SelectColE();//列选主元并消元成上三角
    printf("上三角的结果：\n");
    printM();
    for (int i=row; i>=1; i--) {//回代求解
        for (int j=i+1; j<=row; j++) {
            matrix[i][col] -= matrix[i][j] * matrix[j][col];
        }
        matrix[i][col] /= matrix[i][i];
    }
}

//测试函数
int main() {
    //freopen("gaussian.out", "w", stdout);

#if 0
    row=3;
    col=4;

    //测试用
    matrix[1][1] = 1;
    matrix[1][2] = 3;
    matrix[1][3] = 4;
    matrix[1][4] = 5;
    matrix[2][1] = 1;
    matrix[2][2] = 4;
    matrix[2][3] = 7;
    matrix[2][4] = 3;
    matrix[3][1] = 9;
    matrix[3][2] = 3;
    matrix[3][3] = 2;
    matrix[3][4] = 2;
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
    for (int i=1; i<=row; i++) {
        printf("X%d = %10.2f\n", i, matrix[i][col]);
    }

    return 0;
}