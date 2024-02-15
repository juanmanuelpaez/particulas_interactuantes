#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <array>
#include <vector>
#include <algorithm>
#include <cassert>
#include <time.h>

using namespace std;

//Sistema de N átomos que interactúan mediante un potencial de Lennard-Jones
//dentro de una caja en dos dimensiones de tamaño L × L

//Para el nro de particulas
const int Nc=30;
//numero de partículas
const int N=Nc*Nc;

//paso temporal
const double h=0.005;

//Longitud de la caja
const double L=sqrt(N/0.3);

struct Ranq1 { //Recommended generator for everyday use. The period is ~1.8E19, so it should not be used by an application that makes more than 1E12 calls.
	unsigned long long int v;
	Ranq1(unsigned long long int j) : v(4101842887655102017LL) {
		v ^= j;
		v = int64();
	}
	inline unsigned long long int int64() {//Return 64-bit random integer
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}
	inline double doub() { return 5.42101086242752217E-20 * int64(); } //Return random double-precision floating value in the range 0. to 1
	inline unsigned int int32() { return (unsigned int)int64(); }  //Return 32-bit random integer
};

struct particula
{
    double x, y, vx, vy; //posición y velocidad en los dos ejes
    particula()
    {
        x=0;
        y=0;
        vx=0;
        vy=0;
    }
    particula(double x_, double y_, double v_x, double v_y)
    {
        x=x_;
        y=y_;
        vx=v_x;
        vy=v_y;
    }
    double fuerza_x(particula p)
    {
        //distancia relativa al cuadrado
        double r2=(x-p.x)*(x-p.x)+(y-p.y)*(y-p.y);
        double a=1/r2;
        //fuerza en x
        double fx=4*(12*a*a*a*a*a*a*a*(x-p.x)-6*a*a*a*a*(x-p.x));
        return fx;
    }
    double fuerza_y(particula p)
    {
        //distancia relativa al cuadrado
        double r2=(x-p.x)*(x-p.x)+(y-p.y)*(y-p.y);
        double a=1/r2;
        //fuerza en y
        double fy=4*(12*a*a*a*a*a*a*a*(y-p.y)-6*a*a*a*a*(y-p.y));
        return fy;
    }
};

//función que hace un paso de tiempo, actualizando posiciones y velocidades de cada partícula
void Verlet(particula gas[], double aceleracion[N][4])
{
    //actualizo las posiciones
    for(int i=0; i<N; i++)
    {
        gas[i].x+=h*gas[i].vx+0.5*h*h*aceleracion[i][0];
        gas[i].y+=h*gas[i].vy+0.5*h*h*aceleracion[i][1];
        //Verifico que las partículas no se salgan de la caja rígida
        if(gas[i].x<0)
        {
            gas[i].x*=(-1);
            gas[i].vx*=(-1);
        }
        if(gas[i].x>L)
        {
            gas[i].x-=(2*(gas[i].x-L));
            gas[i].vx*=(-1);
        }
        if(gas[i].y<0)
        {
            gas[i].y*=(-1);
            gas[i].vy*=(-1);
        }
        if(gas[i].y>L)
        {
            gas[i].y-=(2*(gas[i].y-L));
            gas[i].vy*=(-1);
        }
    }
    //calculo las nuevas aceleraciones
    for(int i=0; i<N; i++)
    {
        double ax=0, ay=0;
        for(int j=0; j<N; j++)
        {
            if(j==i)
                continue;
            ax+=gas[i].fuerza_x(gas[j]);
            ay+=gas[i].fuerza_y(gas[j]);
        }
        aceleracion[i][2]=ax;
        aceleracion[i][3]=ay;
    }
    //calculo las nuevas velocidades
    for(int i=0; i<N; i++)
    {
        gas[i].vx+=0.5*h*(aceleracion[i][0]+aceleracion[i][2]);
        gas[i].vy+=0.5*h*(aceleracion[i][1]+aceleracion[i][3]);
        //pongo las aceleraciones nuevas en el lugar de las viejas
        aceleracion[i][0]=aceleracion[i][2];
        aceleracion[i][1]=aceleracion[i][3];
    }
}

//Energía potencial total
double Potencial(particula gas[])
{
    double E_pot=0;
    for(int i=0; i<N; i++)
    {
        for(int j=i+1; j<N; j++)
        {
            //distancia al cuadrado
            double r=(gas[j].x-gas[i].x)*(gas[j].x-gas[i].x)+(gas[j].y-gas[i].y)*(gas[j].y-gas[i].y);
            double a=1/r;
            E_pot+=(4*(a*a*a*a*a*a-a*a*a));
        }
    }
    return E_pot;
}

//Energía cinética total
double Cinetica(particula gas[])
{
    double E_cin=0;
    for(int i=0; i<N; i++)
        E_cin+=(0.5*((gas[i].vx*gas[i].vx)+(gas[i].vy*gas[i].vy)));
    return E_cin;
}

//número de velocidades que quiero guardar
const int nv=5;
//función para guardar las velocidades
void Guarda_V(particula gas[], double velocidad_x[nv][N+1], double velocidad_y[nv][N+1], int nro , int tiempo)
{
    velocidad_x[nro][0]=tiempo;
    velocidad_y[nro][0]=tiempo;
    for(int i=0; i<N; i++)
    {
        velocidad_x[nro][i+1]=gas[i].vx;
        velocidad_y[nro][i+1]=gas[i].vy;
    }
    return;
}

int main()
{
    //Tiempo final
    int T=2000;
    //vector de partículas
    particula gas[N];
    //inicializo las partículas
    Ranq1 myran(time(0));
    double a=L/(Nc+1);
    for(int i=1; i<=Nc; i++)
    {
        for(int j=1; j<=Nc; j++)
        {
            double v_x=1.1*2*((myran.int32()%2)-0.5); //v en x: +1.1 � -1.1
            gas[(Nc*(i-1))+(j-1)]=particula(i*a,j*a,v_x,0);
        }
    }
    //vector de aceleraciones: por cada partícula a_x y a_y a tiempo actual y al tiempo siguiente
    double aceleracion[N][4];
    //calculo las aceleraciones iniciales
    for(int i=0; i<N; i++)
    {
        double a_x=0, a_y=0;
        for(int j=0; j<N; j++)
        {
            if(j==i)
                continue;
            a_x+=gas[i].fuerza_x(gas[j]);
            a_y+=gas[i].fuerza_y(gas[j]);
        }
        aceleracion[i][0]=a_x;
        aceleracion[i][1]=a_y;
    }
    double E_cinetica[T+1];
    double E_potencial[T+1];
    double velocidad_x[nv][N+1];
    double velocidad_y[nv][N+1];
    int contador=0; //para contar cuántas velocidades guardé
    //guardo los valores iniciales
    E_cinetica[0]=Cinetica(gas);
    E_potencial[0]=Potencial(gas);
    Guarda_V(gas,velocidad_x,velocidad_y,contador,0);
    contador++;
    //transcurro el tiempo
    for(int t=1; t<=T; t++)
    {
        if(!(t%50))
            cout << "t = " << t << endl;
        Verlet(gas,aceleracion);
        E_cinetica[t]=Cinetica(gas);
        E_potencial[t]=Potencial(gas);
        if(!(t%(T/(nv-1)))) //guardo algunas distribuciones de velocidad
        {
            Guarda_V(gas,velocidad_x,velocidad_y,contador,t);
            contador++;
        }
    }
    //guardo
    ofstream datos_energia;
    datos_energia.open("11_datos_E_t_2000.dat");
    for(int i=0; i<T+1; i++)
        datos_energia << i << '\t' << E_cinetica[i] << '\t' << E_potencial[i] << '\n';
    datos_energia.close();
    ofstream datos_v;
    datos_v.open("11_datos_velocidad_t_2000.dat");
    for(int j=0; j<N+1; j++)
    {
        for(int i=0; i<nv; i++)
            datos_v << velocidad_x[i][j] << '\t' << velocidad_y[i][j] << '\t';
        datos_v << '\n';
    }
    datos_v.close();
    ofstream mapa;
    mapa.open("11_mapa_final_t_2000.dat");
    for(int i=0; i<N; i++)
        mapa << gas[i].x << '\t' << gas[i].y << '\n';
    mapa.close();
    return 0;
}
