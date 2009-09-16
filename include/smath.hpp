inline double sine(double& x)
{	
    const double B = 4/M_PI;
    const double C = -4/(M_PI*M_PI);

    double y = B * x + C * x * abs(x);

    //  const float Q = 0.775;
    const double P = 0.225;

    y = P * (y * abs(y) - y) + y;   // Q * y + P * y * abs(y)
	return y;
}