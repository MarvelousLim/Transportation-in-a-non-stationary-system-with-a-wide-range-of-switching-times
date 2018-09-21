#include <math.h>
#include <malloc.h>
#include "mt64.h"
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <limits.h>
#include <time.h>
#include <stddef.h>

#define L 2000 //Lattice size

// ��������� ������������� - ����� � ����������� � ��� �����, ������ � �������
struct triangle{
    int flag;// were we here?

    int SV;  // state for vert
    int SH;  //       for hor

};

//���������� 1 � ���-��� p � 0 � ���-��� 1 - p
int genSZ(double p) 
{
    double r;

    r = genrand64_real3();

    if (r < p) return 1;
    else return 0;

}


struct triangle Lattice[L][L];

int corr_size;
int corr_array[L / 2];


// ����������� ����� � ������
int DFS(int X, int Y){

    if (Lattice[X][Y].flag == 1)
        return 0;
	
	Lattice[X][Y].flag = 1; // �� ��� ����
	corr_size++;


    if ((X < L - 2) && (Lattice[X][Y].SH == 1)) //����� �� �� ������ ������?
        DFS(X + 1, Y);

    if ((Y < L - 2) && (Lattice[X][Y].SV == 1)) // �����?
        DFS(X, Y + 1);

    if ((X > 0) && (Lattice[X - 1][Y].SH == 1)) // �����?
        DFS(X - 1, Y);

    if ((Y > 0) && (Lattice[X][Y - 1].SV == 1)) //����?
        DFS(X, Y - 1);

    return 0;

}


//��������� ��������� ��� ���. ������� � ���������� �� �� ������ ������ 
//��� �� ��������� ����� �����, ��� �������� - �� ����, �� ���� ��������

int corr_function(){

    int i, center = L / 2;

    DFS(center, center); //����� ���� ����������� �������, ���� ��� ������� ����� ����� ���� 1

    for(i = 0; i < L / 2; i++){
        if(Lattice[center + i][center].flag == 1)   
            corr_array[i]++;
    }

    return 0;

}


int main()
{
    FILE * corr_output1, * corr_output2;

	//������������� ������� ���������� ��������
    srand((unsigned int) time (NULL)/2); 


	//������������� merssen twister 64bit 
	unsigned long long init[4] = { rand(), rand(), rand(), rand() }, length = 4;
	init_by_array64(init, length);	

    int u, t0, i, j, k, t, Rmax, d, x, y, z, rSquare, center = L / 2, Mean, b, realization_number, Z;
    double p;
	char s1[40], s2[40];

	p = 0.490;
	realization_number = 10000;

	sprintf(s1, "corr_size_L%dp%lf.txt", L, p); //����� ��� n_s
	sprintf(s2, "corr_lenght_L%dp%lf.txt", L, p); //����� ��� ������� \xi
	 
    corr_output1 = fopen(s1, "a"); 
	corr_output2 = fopen(s2, "a");


    for(i = 1; i < L / 2; i++)
        corr_array[i] = 0;
	
    for(b = 0; b < realization_number; b++){

        printf("%d\n", b);

        //initialize arrays
		for (i = 0; i < L; i++) {
			for (j = 0; j < L; j++) {
				Lattice[i][j].flag = 0;

				Lattice[i][j].SV = genSZ(p);
				Lattice[i][j].SH = genSZ(p);
			}
		}
			
		corr_size = 0;

		corr_function();

		fprintf(corr_output1, "%d ", corr_size);
    }
	

    for(i = 0; i < L / 2; i++)
        fprintf(corr_output2, "%d %d ", i, corr_array[i]);

	fclose(corr_output1);
	fclose(corr_output2);

    printf("\7 \7 \7 \7 \7 \7 \7 \7 \7 \7");

    return 0;
}




















