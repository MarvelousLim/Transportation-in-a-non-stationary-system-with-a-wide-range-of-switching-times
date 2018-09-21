#include <math.h>
#include <malloc.h>
#include "mt64.h"
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <limits.h>
#include <time.h>
#include <stddef.h>

#define L 500 //Lattice size
#define PN 1 //Particle number
#define OT 1000000 //Observation time
#define TN 1000000//Transition Number for bonds

// структура треугольничек - точка и примыкающие к ней связи, правая и верхняя
struct triangle{

    double CV;  // characteristic time for vert
    double CH;  //                     for hor

	double AV;  // actual time for vert
	double AH;  //             for hor

    int SV;  // state for vert
    int SH;  //       for hor

	int SVI;  // state for vert //or transition number
	int SHI;  //       for hor

};

struct coord{
    int X;
    int Y;
};

//генерирует число, распределенное как t_0 * e^{u r} при r равномерно на [0,1] с обрезкой на больших числах
double genCdouble(int t0, int u) 
{
    double r;
    r = genrand64_real3();

	return ((double)t0 * exp((double)u * r));
} 

//генерирует число, распределенное как ~ e^{-t\tau} с обрезкой на больших числах
double genlog(double tau)
{
    double r;
    r = genrand64_real3();

	return  (- tau * log(r));
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
    FILE * bond_freq_output;

	//инициализация родного генератора временем
    srand((unsigned int) time (NULL)/2); 


	//инициализация merssen twister 64bit 
	unsigned long long init[4] = { rand(), rand(), rand(), rand() }, length = 4;
	init_by_array64(init, length);	

	int u, t0, i, j, k, t, Rmax, d, x, y, z, rSquare, center = L / 2, Mean, b, realization_number, Z, d1, d2;
    double p;
	char s[50];
	
	p = 0.500;
	t0 = 1;
	u = 60;
	realization_number = 1;

	sprintf(s, "bond_freq_TNtransition%d_L%dT%d_u%dt%d.txt", TN, L, OT, u, t0);
	bond_freq_output = fopen(s, "w");

    for(b = 0; b < realization_number; b++){

        t = 0;
        Rmax = 0;

        printf("%d\n", b);

        //initialize arrays
        for(i = 0; i < L; i++)
            for(j = 0; j < L; j++){

                Lattice[i][j].CV = genCdouble(t0, u);
                Lattice[i][j].CH = genCdouble(t0, u);

                Lattice[i][j].AV = genlog(Lattice[i][j].CV);
                Lattice[i][j].AH = genlog(Lattice[i][j].CH);

                Lattice[i][j].SV = genSZ(p);
                Lattice[i][j].SH = genSZ(p);

				Lattice[i][j].SVI = 0;
				Lattice[i][j].SHI = 0;

            }
		
        //place particles
        for(k = 0; k < PN; k++){
            Particles[k].Y = center;
            Particles[k].X = center;
        }


        while((Rmax < (center - 1)) && (t < OT)){

            for(k = 0; k < PN; k++){
                //обрабатываем все частички

				//делаем нужные вычисления
				if (t % 1000 == 0)
					fprintf(bond_freq_output, "%d  ", t, x, y);

                d = 0;
                x = Particles[k].X;
                y = Particles[k].Y;

				//кусок кода, выбирающий направление прыжка
				d1 = Lattice[x][y].SH + Lattice[x - 1][y].SH + Lattice[x][y].SV + Lattice[x][y - 1].SV;
				if (d1 != 0) {
					d2 = (int)floor(d1 * genrand64_real3());

					d2 -= Lattice[x][y].SH;
					if (d2 < 0) { 
						Particles[k].X++;
						if (Lattice[x][y].SHI < TN) {
							Lattice[x][y].SHI++;
							//fprintf(bond_freq_output, "%lf ", Lattice[x][y].CH);
						}
					}
					else {
						d2 -= Lattice[x - 1][y].SH;
						if (d2 < 0) {
							Particles[k].X--;
							if (Lattice[x - 1][y].SHI < TN) {
								Lattice[x - 1][y].SHI++;
								//fprintf(bond_freq_output, "%lf ", Lattice[x - 1][y].CH);
							}
					}
						else {
							d2 -= Lattice[x][y].SV;
							if (d2 < 0) {
								Particles[k].Y++;
								if (Lattice[x][y].SVI < TN) {
									Lattice[x][y].SVI++;
									//fprintf(bond_freq_output, "%lf ", Lattice[x][y].CV);
								}
							}
							else { 
								Particles[k].Y--; 
								if (Lattice[x][y - 1].SVI < TN) {
									Lattice[x][y - 1].SVI++;
									//fprintf(bond_freq_output, "%lf ", Lattice[x][y - 1].CV);
								}
							}
						}
					}
				}

				//проверяем чтобы частицы не вылетели за решетку
                if (abs(Particles[k].X - center) > Rmax)
                    Rmax = abs(Particles[k].X - center);
                if (abs(Particles[k].Y - center) > Rmax)
                    Rmax = abs(Particles[k].Y - center);

            }


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

			//прогресс
            if (t % (OT / 20) == 0)
                printf("%d%%\n", (int) (100 * t / OT));

			t++;
			//Mean = 0;

        }
		
    }

	fclose(bond_freq_output);

    return 0;

}










