
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <omp.h>

int start = GetTickCount();

int n_cars = 50, j;
FILE * track[100];
typedef struct via {
	double dist;
	double raio;
	double incl;
} Via;
double dx[400];
double pont[400];


Via ant[200];
Via deps[200];

double foract(double Fact, double pos, int j) {
	double Fd, Fu;
	double a;
	double b;
	double z;
	double Fdrag, Fcar, Fdes, fe;
	double deltax;
	double mu = 0.084667;
	double sigma = 0.052839;
	if (fabs(pos) < 150e-3)
		deltax = 0.00005;
	else
		deltax = 0.005;

	z = fabs(fabs(pos) - 76.2e-3);
	a = 1.5857e16* pow(z, 9) - 1.1107e16 * pow(z, 8) + 3.2446e15*pow(z, 7) - 5.099e14*pow(z, 6);
	b = 4.6255e13*pow(z, 5) - 2.3864e12*pow(z, 4) + 6.1436e10 * pow(z, 3) - 4.1609e8*pow(z, 2) + 1.2011e6*z - 401.59;
	Fu = a + b;
	Fd = -4.0875e8*pow(z / 2, 3) + 5.0802e7*pow(z / 2, 2) + 2.6861e5*z / 2 - 560.57;
	if (fabs(pos) <= 76.2e-3)
		Fdrag = 0;
	else {
		//if (pos > 0 && dvel >= 0) {
		//	Fdrag = (Fact + Fu) / 2;
		//}
		//else if (pos > 0 && dvel < 0) {
		//	Fdrag = (Fact + Fd) / 2;
		//}
		//else if (pos < 0 && dvel <= 0) {
		//	Fdrag = -(-Fact + Fu) / 2;
		//}
		//else if (pos < 0 && dvel >0) {
		//	Fdrag = -(-Fact + Fd) / 2;
		//}
		if (pos > 0) {
			Fcar = (Fact + Fu) / 2;
			Fdes = (Fact + Fd) / 2;
		}
		else if (pos < 0) {
			Fcar = (Fact - Fd) / 2;
			Fdes = (Fact - Fu) / 2;
		}
		fe = 1 / (1 + exp(-(pos - dx[j]) / deltax));
		Fdrag = Fcar*fe + (1 - fe)*Fdes;
	}
	//printf("%f\n", Fu);
	return Fdrag;
}

Via interpol(double x, int k) {
	Via w;
	char buff[200];
	int t = 0;
	//printf("%f %f\n", z.dist, y.dist);
	while (t != 1) {
		if (x == ant[k].dist) {
			t = 1;
			w = ant[k];
		}
		else if (x == deps[k].dist) {
			t = 1;
			//fclose(track);
			w = deps[k];
		}
		else if (x > deps[k].dist) {
			ant[k].dist = deps[k].dist;
			ant[k].raio = deps[k].raio;
			ant[k].incl = deps[k].incl;
			fscanf(track[k], "%s", buff);
			deps[k].dist = atof(buff);
			fscanf(track[k], "%s", buff);
			deps[k].raio = atof(buff);
			fscanf(track[k], "%s", buff);
			deps[k].incl = atof(buff);
			//printf("%f %f\n", z.dist, y.dist);

		}
		else if (x > ant[k].dist && x < deps[k].dist) {
			w.dist = x;
			w.raio = deps[k].raio - (deps[k].raio - ant[k].raio)*(deps[k].dist - x) / (deps[k].dist - ant[k].dist);
			w.incl = deps[k].incl - (deps[k].incl - ant[k].incl)*(deps[k].dist - x) / (deps[k].dist - ant[k].dist);
			//fclose(track);
			t = 1;
		}
		else if (x < 0) {
			w.dist = x;
			w.raio = 0;
			w.incl = 0;
		}
	}
	//printf("%f %f\n", ant[k].dist, deps[k].dist);
	return w;
}
double fortrat(int p, double v) {
	double Ft;
	if (p == 1) {
		if (v <= 7.17 / 3.6)
			Ft = 5.7623e3*9.81;
		else
			Ft = -0.0026129*pow(v, 7) + 0.24504*pow(v, 6) - 9.532*pow(v, 5) + 199.17*pow(v, 4) - 2420.2*pow(v, 3) + 17301 * pow(v, 2) - 70041 * v + 1.4299e5;
	}
	else if (p == 2) {
		if (v <= 5.83 / 3.6)
			Ft = 1.624e4*9.81;
		else
			Ft = -0.0081568*pow(v, 7) + 0.75099*pow(v, 6) - 28.57 * pow(v, 5) + 580.98*pow(v, 4) - 6828.1*pow(v, 3) + 46866 * pow(v, 2) - 1.8083e5*v + 3.5164e5;
		//printf("%f\n", Ft);
	}
	else if (p == 3) {
		if (v <= 7.4 / 3.6)
			Ft = 2.7869e4*9.81;
		else
			Ft = -0.0087022*pow(v, 7) + 0.82398*pow(v, 6) - 32.458*pow(v, 5) + 689.87*pow(v, 4) - 8588.9 * pow(v, 3) + 63714 * pow(v, 2) - 2.7433e5 * v + 6.2916e5;
	}
	else if (p == 4) {
		if (v <= 8.2 / 3.6)
			Ft = 3.7076e4*9.81;
		else
			Ft = -0.0096343*pow(v, 7) + 0.92262 * pow(v, 6) - 36.851*pow(v, 5) + 796.97*pow(v, 4) - 10147 * pow(v, 3) + 77544 * pow(v, 2) - 3.4816e5 * v + 8.5269e5;
	}
	else if (p == 5) {
		if (v <= 10.36 / 3.6)
			Ft = 4.2791e4*9.81;
		else
			Ft = -0.0075549 * pow(v, 7) + 0.74542 * pow(v, 6) - 30.87*pow(v, 5) + 697.98*pow(v, 4) - 9393.6 * pow(v, 3) + 77046 * pow(v, 2) - 3.7944e5 * v + 1.0542e6;
	}
	else if (p == 6) {
		if (v <= 11.89 / 3.6)
			Ft = 4.8795e4*9.81;
		else
			Ft = -0.0061512 * pow(v, 7) + 0.62111 * pow(v, 6) - 26.445*pow(v, 5) + 618.43*pow(v, 4) - 8677.8 * pow(v, 3) + 75044 * pow(v, 2) - 3.964e5 * v + 1.2175e6;
	}
	else if (p == 7) {
		if (v <= 13.7119 / 3.6)
			Ft = 5.3118e4*9.81;
		else
			Ft = -0.0043548* pow(v, 7) + 0.45424 * pow(v, 6) - 20.105*pow(v, 5) + 492.73*pow(v, 4) - 7320.5 * pow(v, 3) + 67916 * pow(v, 2) - 3.9148e5 * v + 1.3421e6;
	}
	else if (p == 8) {
		if (v <= 14.5 / 3.6)
			Ft = 5.9e4*9.81;
		else
			Ft = -0.0033802 * pow(v, 7) + 0.36038 * pow(v, 6) - 16.375*pow(v, 5) + 414.35*pow(v, 4) - 6403.7 * pow(v, 3) + 62434 * pow(v, 2) - 3.8367e5 * v + 1.4319e6;
	}
	else
		Ft = 0;

	return Ft;
}

double forres(double mass, double vel, double xdis, int b) {
	double A1, A2, B, C, Fresis;
	int Na, signal;
	Via local;
	if (vel >= 0)
		signal = 1;
	else
		signal = -1;

	vel = fabs(vel);
	if (mass == 162e3) {
		if (vel <= 0.1) {
			A1 = 5;
			A2 = 425;
			B = 0;
			C = 0;
			Na = 8;
		}
		else {
			A1 = 0.1;
			A2 = 8.5;
			B = 0.00938;
			C = 0.0046;
			Na = 8;
		}
	}
	else {
		if (vel <= 0.1) {
			A1 = 15;
			A2 = 425;
			B = 0;
			C = 0;
			Na = 4;
		}
		else {
			A1 = 0.3 / 2;
			A2 = 8.5 / 2;
			B = 0.01398 / 2;
			C = 0.00065 / 2;
			Na = 4;
		}
	}
	if (xdis < 0) {
		Fresis = signal*mass*9.81*(A1*A2*Na / mass + B*vel + C*pow(vel, 2));
	}
	else {
		local = interpol(xdis, b);
		Fresis = signal*(mass*9.81*(A1*A2*Na / mass + B*vel + C*pow(vel, 2)) + mass*9.81*sin(local.incl) + mass*6.116*local.raio);
	}
	return Fresis;
}

int main()
{
	FILE * pFile_Pos;
	FILE * pFile_Vel;
	FILE * pFile_Acel;
	FILE * pFile_Fdg;
	FILE * pFile_Ft;
	FILE * Via;

	Via = fopen("Via.dat", "r");
	pFile_Pos = fopen("Pos.dat", "w");
	pFile_Vel = fopen("Vel.dat", "w");
	pFile_Acel = fopen("Acel.dat", "w");
	pFile_Fdg = fopen("Fdg.dat", "w");
	pFile_Ft = fopen("Ft.dat", "w");
	double x[400], x1[400], x2[400];
	int notch;
	double t = 0;
	int i, t_final;
	double Ftr, Fdg[400], Fr[400], m[400];
	double dt, contador;
	int write_count = 0;
	dt = 0.0001;
	t_final = 100;
#pragma omp parallel for
	for (i = 0; i < n_cars; i++) {
		track[i] = fopen("Via.txt", "r");
		char buff[200];
		fscanf(track[i], "%s", buff);
		ant[i].dist = atof(buff);
		fscanf(track[i], "%s", buff);
		ant[i].raio = atof(buff);
		fscanf(track[i], "%s", buff);
		ant[i].incl = atof(buff); fscanf(track[i], "%s", buff);
		deps[i].dist = atof(buff);
		fscanf(track[i], "%s", buff);
		deps[i].raio = atof(buff);
		fscanf(track[i], "%s", buff);
		deps[i].incl = atof(buff);
		//printf("%f %f", ant[i].dist, deps[i].dist);
	}
	contador = 0.1 / dt;
	//entrar com vetor posição
	//definir velocidades iniciais
	//definir valores das forças

	for (i = 0; i < n_cars; i++) {
		if (i == 0)
			m[i] = 162e3;
		else
			m[i] = 100e3;

		x[i] = -i * 9.825;
		x1[i] = 0;
		Fdg[i] = 0;
		if (i != n_cars - 1)
			dx[i] = x[i] - x[i + 1];
	}
	while (t < t_final) {
		if (t < 200)
			notch = 3;
		else if (t < 400)
			notch = 4;
		else if (t < 600)
			notch = 5;
		Ftr = fortrat(notch, x1[0]);
		//printf("%f\n", Ftr);
#pragma omp parallel for
		for (i = 0; i < n_cars; i++) {
			Fr[i] = forres(m[i], x1[i], x[i], i);
			if (i == n_cars - 1) {
				Fdg[i] = 0;
			}
			else {
				Fdg[i] = foract(Fdg[i], (fabs(x[i + 1] - x[i])) - 9.825, i);
				//printf("%f %f %f\n", Fdg[i], x1[i]-x1[i+1], (x[i + 1] - x[i]) - 9.825);
			}
		}
#pragma omp parallel for shared (Fr[400])
		for (i = 0; i < n_cars; i++) {
			//printf("%f\n", Fdg[i]);
			if (i == 0)
				x2[i] = (Ftr - Fr[i] - Fdg[i]) / m[i];
			else if (i == n_cars - 1) {
				//if (x1[i] == 0 && Fr[i] > Fdg[i])
				//	x2[i] = 0;
				//else
				x2[i] = (-Fr[i] + Fdg[i - 1]) / m[i];
			}
			else {
				//if (x1[i] == 0 && (Fr[i] > Fdg[i - 1] - Fdg[i]))
				//	x2[i] = 0;
				//else
				x2[i] = (-Fr[i] + Fdg[i - 1] - Fdg[i]) / m[i];
			}
			//printf("%f\n", Fr[i]);
			dx[i] = (fabs(x[i] - x[i + 1]) - 9.825);

			x1[i] = x1[i] + dt*x2[i];
			x[i] = x[i] + dt*x1[i];
		}
		if (write_count == 0 || write_count == contador) {
			write_count = 0;
			int n_int;
			fprintf(pFile_Pos, "%f ", t);
			for (n_int = 0; n_int < n_cars; n_int++) {
				fprintf(pFile_Pos, "%f ", x[n_int]);
			}
			fprintf(pFile_Pos, "\n");

			fprintf(pFile_Vel, "%f ", t);
			for (n_int = 0; n_int < n_cars; n_int++) {
				fprintf(pFile_Vel, "%f ", x1[n_int]);
			}
			fprintf(pFile_Vel, "\n");

			fprintf(pFile_Acel, "%f ", t);
			for (n_int = 0; n_int < n_cars; n_int++) {
				fprintf(pFile_Acel, "%f ", Fr[n_int]);
			}
			fprintf(pFile_Acel, "\n");

			fprintf(pFile_Fdg, "%f ", t);
			for (n_int = 0; n_int < n_cars; n_int++) {
				fprintf(pFile_Fdg, "%f ", Fdg[n_int]);
			}
			fprintf(pFile_Fdg, "\n");

			fprintf(pFile_Ft, "%f ", t);
			for (n_int = 0; n_int < n_cars; n_int++) {
				fprintf(pFile_Ft, "%f ", Ftr);
			}
			fprintf(pFile_Ft, "\n");
		}
		write_count++;
		t = t + dt;
	}
	int End = GetTickCount() - start;
	printf("%d\n", End);
	getchar();

	return 0;
}

