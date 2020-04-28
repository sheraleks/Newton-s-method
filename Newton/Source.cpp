#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <conio.h>
using namespace std;
const double eps = 1e-7;
static int f_calc;

struct point
{
	double x, y;
	point operator+ (const point& X)
	{
		point res;
		res.x = this->x + X.x;
		res.y = this->y + X.y;
		return res;
	}
	point operator- (const point& X)
	{
		point res;
		res.x = this->x - X.x;
		res.y = this->y - X.y;
		return res;
	}
	point operator/ (int x)
	{
		point res;
		res.x = this->x / x;
		res.y = this->y / x;
		return res;
	}
	point operator* (int x)
	{
		point res;
		res.x = this->x * x;
		res.y = this->y * x;
		return res;
	}
	point operator-()
	{
		this->x = -this->x;
		this->y = -this->y;
		return *this;
	}
};

double Function(point x)
{
	f_calc++;
	return (100 * (x.y - x.x) *(x.y - x.x) + (1 - x.x) *(1 - x.x)); 
	//return 100 * (x.y - x.x * x.x) *(x.y - x.x*x.x) + (1 - x.x) *(1 - x.x); 
}

double* Gradient(point x) // градиент 
{
	double* res = new double[2];


	res[0] = (-2.0 * 100 * (x.y - x.x) - 2 * (1 - x.x)); // по х 
	res[1] = 100 * 2 * (x.y - x.x); // по у 

	//Розенброк
	//res[0] = 2 * (200 * x.x*x.x*x.x - 200 * x.x*x.y + x.x - 1); 
	//res[1] = 200 * (x.y - x.x*x.x); 

	return res;
}

double** Gesse(point x) // считаем вторые производные 
{
	double** res = new double* [2];
	res[0] = new double[2];
	res[1] = new double[2];
	double t;


	res[0][0] = 202; // xx 
	res[0][1] = -200; // xy 
	res[1][0] = -200; // yx 
	res[1][1] = 200; // yy 

	//res[0][0] = -400 * (x.y - x.x*x.x) + 800 * x.x*x.x + 2; 
	//res[0][1] = -400 * x.x; 
	//res[1][0] = -400 * x.x; 
	//res[1][1] = 200; 

	return res;
}

double** matrix_H(point x)
{
	double** res = Gesse(x);
	double opred = 0.0;
	double t;
	t = res[0][0];
	res[0][0] = res[1][1];
	res[1][1] = t;
	t = res[0][1];
	res[0][1] = -res[1][0];
	res[1][0] = -t;
	opred += res[0][0] * res[1][1] - res[1][0] * res[0][1];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			res[i][j] /= opred;

	return res;
}

double** multipl_grad(double** H, point x)
{
	double** res = new double* [2];
	res[0] = new double;
	res[1] = new double;

	double* grad = Gradient(x);
	res[0][0] = H[0][0] * grad[0] + H[0][1] * grad[1];
	res[1][0] = H[1][0] * grad[0] + H[1][1] * grad[1];
	return res;
}

double lymbda(double** H, point x)
{
	double h = 1;
	double delta = 0.01;
	point res1, res2, res3, m;
	res1.x = x.x - h * H[0][0];
	res1.y = x.y - h * H[1][0];

	res2.x = x.x - (delta + h) * H[0][0];
	res2.y = x.y - (h + delta) * H[1][0];

	res3.x = x.x - (h - delta) * H[0][0];
	res3.y = x.y - (h - delta) * H[1][0];

	if (Function(res1) < Function(res2) && Function(res3) > Function(res1))
		return h;
	else
	{
		if (Function(res2) > Function(res3))
		{
			delta = -delta;
			m = res1;
		}
		else
		{
			m = res1;
		}

	}
	h += delta;
	res1.x = x.x - h * H[0][0];
	res1.y = x.y - h * H[1][0];
	while (Function(res1) < Function(m))
	{

		h += delta;
		m = res1;
		res1.x = x.x - h * H[0][0];
		res1.y = x.y - h * H[1][0];
	}
	return h;
}

void multipl_H(double alpha, double** H)
{
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 1; j++)
			H[i][j] *= alpha;
}

point NewPoint(point x, int iter, double& h)
{
	point res = x;
	double** H = matrix_H(x);
	double** grad_f;
	grad_f = multipl_grad(H, x);
	h = lymbda(grad_f, x);
	multipl_H(h, grad_f);
	res.x = x.x - grad_f[0][0];
	res.y = x.y - grad_f[1][0];
	return res;
}

double Norm(point X)//норма
{
	double res = 0.0;
	double* g = Gradient(X);
	res += sqrt(g[0] * g[0] + g[1] * g[1]);
	return res;
}

point Newton(point x0, double h, double eps)
{
	point x = x0;
	int iter = 0;
	do
	{
		point mem = x;
		iter++;
		x = NewPoint(x, iter, h);

		point s;
		s.x = x.x - mem.x;
		s.y = x.y - mem.y;
		double deltaF;
		deltaF = fabs(Function(x) - Function(mem));

		double* grad = Gradient(x);
		double** H = Gesse(x);

		printf("(%f,%f) %f (%f,%f) %f %f %f %f (%f,%f) \n",
			x.x, x.y, Function(x), s.x, s.y, h, fabs(s.x), fabs(s.y), deltaF, grad[0], grad[1]);
		printf("(%f,%f,%f,%f)\n", -H[0][0], -H[0][1], -H[1][0], -H[1][1]);
		h = 0.5;
	} while (fabs(Norm(x)) > eps);
	printf("\niter = %d\n", iter);
	printf("count = %d\n", f_calc);
	return x;
}

void main()
{
	double h = 1;
	point X;
	point p0;
	p0.x = 2;
	p0.y = 2;
	X = Newton(p0, h, eps);

	printf("X is: %.12lf, %.12lf\n", X.x, X.y);
	printf("func(x) = %.12lf\n", Function(X));

	_getch();
}
