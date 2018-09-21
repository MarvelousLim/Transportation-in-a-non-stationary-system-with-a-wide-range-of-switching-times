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
#define PN 4000 //Particle number
#define OT 3000000 //Observation time

// структура треугольничек - точка и примыкающие к ней связи, правая и верхняя
struct triangle{
    double CV;  // characteristic time for vert
    double CH;  //                     for hor
	double AV;  // actual time for vert
	double AH;  //             for hor
    int SV;  // state for vert
    int SH;  //       for hor
};

struct coord{
    int X;
    int Y;
};

//генерирует число, распределенное как t_0 * e^{ur} при r равномерно на [0,1]
double gendouble(int t0, int u) 
{
    double r;
    r = genrand64_real3();
	return ((double)t0 * exp((double) u * r));
} 

//генерирует число, распределенное как ~ e^{-t\tau}
double genlog(double tau)
{
    double r;
    r = genrand64_real3();
    return (double) (- tau * log(r));

}

//генерирует 1 с вер-тью p и 0 с вер-тью 1 - p
int genSZ(double p) 
{
    double r;

    r = genrand64_real3();

    if (r < p) return 1;
    else return 0;

}


struct triangle Lattice[L][L];
struct coord Particles[PN];


int main()
{
    FILE * histogramm_output;

	//инициализация родного генератора временем
    srand((unsigned int) time (NULL)/2); 


	//инициализация merssen twister 64bit 
	unsigned long long init[4] = { rand(), rand(), rand(), rand() }, length = 4;
	init_by_array64(init, length);	

    int u, t0, i, j, k, t, Rmax, d, x, y, z, rSquare, center = L / 2, Mean, b, realization_number, Z, d1, d2;
    double p;
	char s[40];

	t0 = 100;
	u = 120;

    p = 0.500;
    realization_number = 1;

    for(b = 0; b < realization_number; b++){

        t = 0;
        Rmax = 0;

        printf("%d\n", b);

        //initialize arrays
        for(i = 0; i < L; i++)
            for(j = 0; j < L; j++){
                Lattice[i][j].CV = gendouble(t0, u);
                Lattice[i][j].CH = gendouble(t0, u);
                Lattice[i][j].AV = genlog(Lattice[i][j].CV);
                Lattice[i][j].AH = genlog(Lattice[i][j].CH);
                Lattice[i][j].SV = genSZ(p);
                Lattice[i][j].SH = genSZ(p);
            }

        //place particles
        for(k = 0; k < PN; k++){
            Particles[k].Y = center;
            Particles[k].X = center;
        }


        while((Rmax < center * center) && (t < OT)){
            t++;
            //Mean = 0;

            for(k = 0; k < PN; k++){
                //обрабатываем все частички
                x = Particles[k].X;
                y = Particles[k].Y;

				//кусок кода, выбирающий направление прыжка
				d1 = Lattice[x][y].SH + Lattice[x - 1][y].SH + Lattice[x][y].SV + Lattice[x][y - 1].SV;
				if (d1 != 0) {
					d2 = (int)floor(d1 * genrand64_real3());

					d2 -= Lattice[x][y].SH;
					if (d2 < 0) Particles[k].X++;
					else {
						d2 -= Lattice[x - 1][y].SH;
						if (d2 < 0) Particles[k].X--;
						else {
							d2 -= Lattice[x][y].SV;
							if (d2 < 0) Particles[k].Y++;
							else Particles[k].Y--;
						}
					}
				}

				//проверяем чтобы частицы не вылетели за решетку
                if (((Particles[k].X - center) * (Particles[k].X - center) + (Particles[k].Y - center) * (Particles[k].Y - center)) > Rmax * Rmax)
					Rmax = (Particles[k].X - center) * (Particles[k].X - center) + (Particles[k].Y - center) * (Particles[k].Y - center);

                //делаем нужные вычисления
                //Mean += (Particles[k].X - center) * (Particles[k].X - center) + (Particles[k].Y - center) * (Particles[k].Y - center);

            }

			//прогресс
			if (100 * t % OT == 0)
				printf("%d%%\n", (int)(100 * t / OT));

            for (i = 0; i < L; i++){
                for(j = 0; j < L; j++){
                    Lattice[i][j].AV -= 1;
                    Lattice[i][j].AH -= 1;

                    if (Lattice[i][j].AV < 1){
                        Lattice[i][j].SV = 1 - Lattice[i][j].SV;
                        Lattice[i][j].AV = genlog(Lattice[i][j].CV);
                    }

                    if (Lattice[i][j].AH < 1){
                       Lattice[i][j].SH = 1 - Lattice[i][j].SH;
                       Lattice[i][j].AH = genlog(Lattice[i][j].CH);
                    }
                }
            }

			//главный вывод
			
			if (t == 100) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}
			
			if (t == 300) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}

			if (t == 1000) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}

			if (t == 3000) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}

			if (t == 10000) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}
			
			if (t == 30000) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}

			if (t == 100000) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}

			if (t == 300000) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}

			if (t == 1000000) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}


			if (t == 3000000) {

				sprintf(s, "P_L%dOT%d_u%dt0%d.txt", L, t, u, t0);
				histogramm_output = fopen(s, "a");
				for (i = 0; i < PN; i++)
					fprintf(histogramm_output, "%lf ", sqrt((double)((Particles[i].X - center) * (Particles[i].X - center) + (Particles[i].Y - center) * (Particles[i].Y - center))));
				fclose(histogramm_output);

			}
			

        }


    }

    printf("\7");

    return 0;
}










