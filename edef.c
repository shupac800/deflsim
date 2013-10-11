#include <stdio.h>
#include <math.h>

#define PI	3.14159265359
#define MU	1.2566370614e-6 * 4
#define VACC	19500
#define CRT_LEN	0.27
#define ME	9.10938215e-31
#define QE	1.60217646e-19

typedef struct
{
	double x,y,z;
} Vector;

typedef struct
{
	Vector origin, slope;
} Wire;

Vector Vector_Create(double x, double y, double z);

// globals
static int loops = 40;
static Wire w[40];  // this syntax requires a constant in []
static double Minor_radius = 0.06 / 2.0;
static double Major_radius = 0.127 / 2.0;
static double Yoke_length = 0.07;
static double curr = 4.0;
static double K = MU / (4.0 * PI);

Vector Vector_Create(double x, double y, double z) {
	Vector vc;

	vc.x = x;
	vc.y = y;
	vc.z = z;

	return vc;
}

Vector Vector_Mult(Vector v1, double num) {
	Vector vm;

	vm.x = v1.x * num;
	vm.y = v1.y * num;
	vm.z = v1.z * num;

	return vm;
}

Vector Cross_Product(Vector v1, Vector v2) {
	Vector cp;

	cp.x = v1.y * v2.z - v1.z * v2.y;
	cp.y = v1.z * v2.x - v1.x * v2.z;
	cp.z = v1.x * v2.y - v1.y * v2.x;

	return cp;
}

Vector Vector_Add(Vector v1, Vector v2) {
	Vector va;

	va.x = v1.x + v2.x;
	va.y = v1.y + v2.y;
	va.z = v1.z + v2.z;

	return va;
}

Vector Vector_Sub(Vector v1, Vector v2) {
	Vector va;

	va.x = v1.x - v2.x;
	va.y = v1.y - v2.y;
	va.z = v1.z - v2.z;

	return va;
}

double Magnitude(Vector v) {
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

void generate_wire_formulas() {
	int i;
	// for some reason, defining wire angles as INT's gives wrong answer
	double first_wire_angle = 25.0, last_wire_angle = 80.0;
	double angle_inc = (last_wire_angle - first_wire_angle) / (loops-1);
	double theta;
	Vector origin, endpoint, slope;

	for(i=0; i<loops; i++) {
		theta = (first_wire_angle + i * angle_inc) * PI/180;
		origin.x = Minor_radius * cos(theta);
		origin.y = Minor_radius * sin(theta);
		origin.z = 0.0;

		endpoint.x = Major_radius * cos(theta);
		endpoint.y = Major_radius * sin(theta);
		endpoint.z = Yoke_length;

		slope = Vector_Mult(Vector_Sub(endpoint, origin),1.0/Yoke_length);

		w[i].origin = origin;
		w[i].slope = slope;
	}
}

Vector Calculate_B(Vector p) {
	int i, j, segments = 1000;
	double z = 0, dz = Yoke_length / segments;
	Vector lt, rt, dlt, dbt, lb, rb, dlb, dbb, b;
	b = Vector_Create(0,0,0);

	for(j=0; j < segments; j++) {
		for(i=0; i < loops; i++) {
			lt = Vector_Add(w[i].origin,Vector_Mult(w[i].slope,z));
			rt = Vector_Sub(p,lt);
			dlt = Vector_Mult(w[i].slope,dz);
			dbt = Vector_Mult(
				Cross_Product(dlt,rt),
				curr * K / pow(Magnitude(rt),3)
			);
	
			lb = lt;
			lb.y = -lb.y;
			rb = Vector_Sub(p,lb);
			dlb = dlt;
			dlb.y = -dlb.y;

			dbb = Vector_Mult(
				Cross_Product(dlb,rb),
				-curr * K / pow(Magnitude(rb),3)
			);
	
			// left-right symmetry will double x-component
			// left-right symmetry will cancel y-component
			// z-components will cancel because of top-bottom symmetry
			b.x += 2 * (dbt.x + dbb.x);
		}
		z += dz;
	}
	return b;
}

Vector Calculate_F(Vector b, Vector v) {
	Vector f = Cross_Product(Vector_Mult(v,-QE),b);

	// add contribution of electric field
	f.z += VACC * QE / CRT_LEN;

	return f;
}

Vector Calculate_A(Vector f, Vector v) {
	// calculate relativistic mass of electron
	double em = ME / sqrt(1 - (pow(Magnitude(v),2) / 9e16));

	Vector a = Vector_Mult(f,1.0/em);

	return a;
}

Vector Calculate_dv(Vector a, double dt) {
	return Vector_Mult(a,dt);
}

Vector Calculate_dp(Vector v, double dt) {
	return Vector_Mult(v,dt);
}
		
int main(int argc, char *argv[]) {
	generate_wire_formulas();

	Vector b, f, a, v, p;
	double dt = 5e-11;
	int points_calculated = 0;

	v = Vector_Create(0,0, sqrt(2 * 1000 * 1.60e-19 / ME));
	p = Vector_Create(0,0,-0.01);

	while (p.z <= 0.26) {
		b = Calculate_B(p);
		f = Calculate_F(b,v);
		a = Calculate_A(f,v);
		v = Vector_Add(v,Calculate_dv(a,dt));
		p = Vector_Add(p,Calculate_dp(v,dt));
		points_calculated++;
	}

	printf("%d points calculated\n",points_calculated);
	printf("Final v: %15.15f = %2.0f%% c\n",v.z,100*v.z/3e8);
	printf("%15.15f inches deflection at %f A\n",p.y*100/2.54,curr);

	return 0;
}
