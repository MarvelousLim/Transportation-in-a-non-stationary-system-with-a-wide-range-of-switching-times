#include <math.h>
#include "mt64.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <stddef.h>

#define L 1500 //Lattice size
#define PN 1000 //Particle number
#define OT 1000000 //Observation time

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
	return (double)t0 * exp((double)u * r);
} 

//генерирует число, распределенное как ~ e^{-t\tau}
double genlog(double tau)
{
    double r;
    r = genrand64_real3();
    return (double) (-tau * log(r));
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
    FILE * D_output;

	//инициализация родного генератора временем
    srand((unsigned int) time (NULL)/2); 

	//инициализация merssen twister 64bit 
	unsigned long long init[4] = { rand(), rand(), rand(), rand() }, length = 4;
	init_by_array64(init, length);	

    int u, t0, i, j, k, t, Rmax, d1, d2, x, y, center = L / 2, Mean, b, realization_number;
    double p;
	char s[40];

	t0 = 100;
	u = 30;

	sprintf(s, "MaxPower_L%dOT%d_u%dt%d.txt", L, OT, u, t0);
    D_output = fopen(s, "w");

    p = 0.500;
    realization_number = 1;

    for(b = 0; b < realization_number; b++){

        t = 0;
        Rmax = 0;

        printf("%d\n", b);

        //initialize arrays
		for (i = 0; i < L; i++) {
			for (j = 0; j < L; j++) {
				Lattice[i][j].CV = gendouble(t0, u);
				Lattice[i][j].CH = gendouble(t0, u);
				Lattice[i][j].AV = genlog(Lattice[i][j].CV);
				Lattice[i][j].AH = genlog(Lattice[i][j].CH);
				Lattice[i][j].SV = genSZ(p);
				Lattice[i][j].SH = genSZ(p);
			}
		}

        //place particles
        for(k = 0; k < PN; k++){
            Particles[k].Y = center;
            Particles[k].X = center;
        }


        while((Rmax < center) && (t < OT)){

            for(k = 0; k < PN; k++){
                //обрабатываем все частички
                x = Particles[k].X;
                y = Particles[k].Y;

				//кусок кода, выбирающий направление прыжка
				d1 = Lattice[x][y].SH + Lattice[x - 1][y].SH + Lattice[x][y].SV + Lattice[x][y - 1].SV;
				if (d1 != 0) {
					d2 = (int) floor(d1 * genrand64_real3());

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
                if (abs(Particles[k].X - center) >= Rmax)
                    Rmax = abs(Particles[k].X - center);
                if (abs(Particles[k].Y - center) >= Rmax)
                    Rmax = abs(Particles[k].Y - center);
            }

			for (i = 0; i < L; i++) {
				for (j = 0; j < L; j++) {
					Lattice[i][j].AV -= 1;
					Lattice[i][j].AH -= 1;

					if (Lattice[i][j].AV < 1) {
						Lattice[i][j].SV = 1 - Lattice[i][j].SV;
						Lattice[i][j].AV = genlog(Lattice[i][j].CV);
					}

					if (Lattice[i][j].AH < 1) {
						Lattice[i][j].SH = 1 - Lattice[i][j].SH;
						Lattice[i][j].AH = genlog(Lattice[i][j].CH);
					}
				}
			}

			//главный вывод
			if (t % 1000 == 0) {
				Mean = 0;

				//делаем нужные вычисления
				for (k = 0; k < PN; k++)
					Mean += (Particles[k].X - center) * (Particles[k].X - center) + (Particles[k].Y - center) * (Particles[k].Y - center);


				fprintf(D_output, "%d %lf ", t, ((double)Mean) / PN);
				//for (k = 0; k < PN; k++)
				//	fprintf(D_output, "%d %d %d ", t, Particles[k].X - center, Particles[k].Y - center);


				printf("%d%%\n", (100 * t) / OT);

			}


			t++;
        }

    
    }

    //printf("\7");

    return 0;
}