/*
*����: ��ѡ��Ԫ��Ԫ�� 
*@Language: C++
*@File Name: Gaussian_wp.cpp
*���룺GBK
*/
#include <bits/stdc++.h>

using namespace std;

#define __OUTPUT

const int MAXN = 6;//1e2+2;//����������ֵ
double matrix[MAXN][MAXN];//�������ע�����һ�д�

int row;//��������ʵ����
int col;//��������ʵ����

//�������
void printM() {
    for (int i=1; i<=row; i++) {
        for (int j=1; j<=col; j++) {
            printf("%12.4lf", matrix[i][j]);
        }
        printf("\n");
    }
}

//ѡ������Ԫ��������Ԫ
void SelectColE() {
    double temp; //���ڼ�¼��Ԫʱ������
    for (int i=1; i<row; i++) {
        int r = i;
        for (int j=i+1; j<=col; j++) {
            if (fabs(matrix[j][i]) > fabs(matrix[r][i])) {
                r = j;
            }
        }
        if (r != i) {
            for (int j=i; j<=col; j++) {
                swap(matrix[i][j], matrix[r][j]);//�������Ԫ�����н���
            }
        }
        for (int j=i+1; j<=row; j++) {//��Ԫ
            temp=matrix[j][i]/matrix[i][i];
            for (int k=i; k<=col; k++) {
                matrix[j][k] -= matrix[i][k]*temp;
            }
        }

#if defined(__OUTPUT)
        printf("��%d����Ԫ��\n", i);
        printM();
#endif
    }
}

//��˹��Ԫ��(��ѡ��Ԫ)
void Gauss() {
    SelectColE();//��ѡ��Ԫ����Ԫ��������
    printf("�����ǵĽ����\n");
    printM();
    for (int i=row; i>=1; i--) {//�ش����
        for (int j=i+1; j<=row; j++) {
            matrix[i][col] -= matrix[i][j] * matrix[j][col];
        }
        matrix[i][col] /= matrix[i][i];
    }
}

//���Ժ���
int main() {
    //freopen("gaussian.out", "w", stdout);

#if 0
    row=3;
    col=4;

    //������
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