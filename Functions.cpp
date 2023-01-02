#include "Header.h";

//function for reading files into MatrixXd
void read(const string& f, MatrixXd& x)
{
	int row;
	int col;
	vector <vector<double>> data;
	col = 0;
	row = 0;
	ifstream in(f);
	if (in.fail())
	{
		cout << "couldn't open the file" << endl;
		system("pause");
		exit(1);
	}

	string line;
	if (!getline(in, line))
	{
		cout << "Couldn't read the first line" << endl;
		system("pause");
		exit(1);
	}
	istringstream stream(line);
	stream.str(line);
	// creates a vector of rows based on the number of rows in the file
	vector <double> rows;
	double value;
	while (stream >> value)
		rows.push_back(value);
	int cols = rows.size();
	data.push_back(rows);
	// if the number of columns isn't the same, the code breaks
	while (!in.eof())
	{
		if (!getline(in, line) && !in.eof())
		{
			cout << "couldn't read the file" << endl;
			system("pause");
			exit(1);
		}

		stream.clear(); //clears the stream for the next row
		stream.str(line);

		rows.clear();
		for (int k = 0; k < cols; ++k)
		{
			if (stream >> value)
			{
				rows.push_back(value);
			}
		}
		data.push_back(rows); // adds the line to the matrix
	}
	in.close();
	row = data.size();
	col = cols;
	x.resize(row, col);
	for (int j = 0; j < row; j++)
	{
		for (int k = 0; k < col; k++)
		{
			x(j, k) = (data[j])[k];// takes the data matrix and puts its data into the Eigen matrix
		}
	}
	return;
}
//function for writing matrices into text files
void write(const string & f, MatrixXd & x)
{
	ofstream out(f, ios::out);
	if (out.fail())
	{
		cout << "can't open the file" << endl;
		system("pause");
		exit(1);
	}
	for (int i = 0; i < x.rows(); i++) // prints the matrix line by line
	{
		for (int j = 0; j < x.cols(); j++)
		{
			out << fixed;
			out << setprecision(6);
			out << x(i, j);
			if (j != x.cols() - 1)
			{

				out << " ";
			}

		}
		// The file adds a new line unless it's at the last line
		if (i != (x.rows() - 1))
		{
			out << endl;
		}
	}
	out.close();
	return;
}

//function to create icf files
void writeicf(const string& f, MatrixXd& x)
{
	ofstream out(f, ios::out);
	if (out.fail())
	{
		cout << "can't open the file" << endl;
		system("pause");
		exit(1);
	}
	out << fixed;
	for (int a = 0; a < x.rows(); a++)
	{
		if (x(a, 3) < 100)
		{
			out << setprecision(0) << x(a, 0) << ".tif" << "\t" << "MD" << setprecision(0) << x(a, 3) << "\t" << setprecision(6) << x(a,1) << "\t" << x(a,2) << "\t" << 1.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 1.0;
		}
		else
		{
			out << setprecision(0) << x(a, 0) << ".tif" << "\t" << setprecision(0) << x(a, 3) << "\t" << setprecision(6) << x(a, 1) << "\t" << x(a, 2) << "\t" <<  1.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 1.0;
		}
		if (a != (x.rows() - 1))
		{
			out << endl;
		}
	}
	out.close();
	return;
}

//function to create gcp file
void writegcp(const string& f, MatrixXd& x2)
{
	ofstream out(f, ios::out);
	if (out.fail())
	{
		cout << "can't open the file" << endl;
		system("pause");
		exit(1);
	}
	out << fixed;
	for (int a = 0; a < x2.rows(); a++)
	{
		
			
				if (x2(a, 3) < 100)
				{
					out << "MD" << setprecision(0) << x2(a, 3) << "\t" << setprecision(6) << x2(a, 0) << "\t" << x2(a, 1) << "\t"  << x2(a,2) << "\t" << 10000 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 10000 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 <<"\t" << 10000;
				}
				else if (x2(a,3) > 800)
				{
					out << setprecision(0) << x2(a, 3) << "\t" << setprecision(6) << x2(a, 0) << "\t" << x2(a, 1) << "\t" << x2(a,2) << "\t" << 10 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 10 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 10;
				}
				else
				{
					out << setprecision(0) << x2(a, 3) << "\t" << setprecision(6) << x2(a, 0) << "\t" << x2(a, 1) << "\t" << x2(a,2) << "\t" << 1 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 1 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 1;
				}
				if (a != (x2.rows() - 1))
				{
					out << endl;
				}
			}
	out.close();
}

// function to write orientation file
void writeori(const string& f, MatrixXd& eop)
{
	ofstream out(f, ios::out);
	if (out.fail())
	{
		cout << "can't open the file" << endl;
		system("pause");
		exit(1);
	}
	out << fixed;
	for (int i = 0; i < eop.rows(); i++)
	{
		out << setprecision(0) << eop(i, 0) << ".tif" << "\t" << "CAMERA" << "\t" << setprecision(3) << 0.005 << setprecision(1) << "\t" << 0.0 << "\t" << 0 << "\t" << 0 << endl;
		out << 0 << endl;
		out << 0 << endl;
		out << 1 << endl;
		out << 1 << "\t" << setprecision(1) << i+1 << "\t" <<setprecision(6) << eop(i, 1) << "\t" << eop(i, 2) << "\t" << eop(i, 3) << "\t" << eop(i, 4) << "\t" << eop(i, 5) << "\t" << eop(i, 6);

		if (i != (eop.rows() - 1))
		{
			out << endl;
		}
	}
	return;
}

//function to write a .cam file
void writecam(const string& f)
{
	ofstream out(f, ios::out);
	if (out.fail())
	{
		cout << "can't open the file" << endl;
		system("pause");
		exit(1);
	}
	out << fixed;
	out << "! Camera ID Type \t xp \t yp \t c" << endl;
	out << "CAMERA \t FRAME \t" << -0.069 << "\t" << -0.104 << "\t" << 55.046 << endl;
	out << "! Variance-covariance matrix of xp, yp, c" << endl;
	out << "1.0e-8 \t 0.0 \t 0.0" << endl << "0.0 \t 1.0e-8 \t 0.0" << endl << "0.0 \t 0.0 \t 1.0e-8" << endl;
	out << "! Number of fiducials (if more than 0, include xy image coord. of the fid. marks)" << endl;
	out << "0" << endl << "! User defined Ro" << endl << "1" << endl << "! Distortion Model used: 1 = UofC, 2 = SMAC" << endl << "1" << endl;
	out << "! Number of distortion parameters and array elements" << endl << "6" << endl << "0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0" << endl;
	out << "! Variance-covariance matrix of distortion parameters" << endl;
	out << "1.0e-8 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0" << endl;
	out << "0.0 \t 1.0e-8 \t 0.0 \t 0.0 \t 0.0 \t 0.0" << endl;
	out << "0.0 \t 0.0 \t 1.0e-8 \t 0.0 \t 0.0 \t 0.0" << endl;
	out << "0.0 \t 0.t \t 0.0 \t 1.0e-8 \t 0.0 \t 0.0" << endl;
	out << "0.0 \t 0.0 \t 0.0 \t 0.0 \t 1.0e-8 \t 0.0" << endl;
	out << "0.0 \t 0.t \t 0.0 \t 0.0 \t 0.0 \t 1.0e-8" << endl;
	out << "! GPS offset: dx, dy, dz and variance-covariance matrix" << endl;
	out << "0.0 \t 0.0 \t 0.0" << endl;
	out << "1.0e-8 \t 0.0 \t 0.0" << endl;
	out << "0.0 \t 1.0e-8 \t 0.0" << endl;
	out << "0.0 \t 0.0 \t 1.0e-8";
	return;
}

//interpolation for radial distortion correction
MatrixXd interpolation(MatrixXd& distortion, MatrixXd& gcp, MatrixXd& x, double& xp, double& yp)
{
	// initialize matrices for radial distance, delta r, delta x and delta y
	MatrixXd radd(gcp.rows(), 1);
	MatrixXd correction(gcp.rows(), 1);
	MatrixXd dx(gcp.rows(), 1), dy(gcp.rows(), 1);
	for (int i = 0; i < gcp.rows(); i++)
	{
		// get radial distance for each point
		radd(i, 0) = sqrt(pow((gcp(i, 1) - xp), 2) + pow((gcp(i, 2) - yp), 2));
		for (int j = 0; j < distortion.rows(); j++)
		{
			if (distortion(j, 0) == floor(radd(i, 0)))
			{
				// find the interpolated value of radial distortion for all points
				correction(i, 0) = (radd(i, 0) - distortion(j, 0)) * (distortion(j + 1, 1) - distortion(j, 1)) / (distortion(j + 1, 0) - distortion(j, 0)) + distortion(j, 1);
				
			}
		}
		// get correction to x and y using the delta r
		dx(i, 0) = correction(i, 0) * (gcp(i, 1) - xp) / radd(i, 0);
		dy(i, 0) = correction(i, 0) * (gcp(i, 2) - yp) / radd(i, 0);
		// assign the values into matrix x for use in the similarity transformation
		x(i, 0) = gcp(i, 0);
		x(i, 1) = gcp(i, 1) - dx(i, 0);
		x(i, 2) = gcp(i, 2) - dy(i, 0);
		x(i, 3) = gcp(i, 3);
	}
	
	return x;
}

MatrixXd Amat(MatrixXd& gcp)
{
	MatrixXd A(gcp.rows() * 2, 4);
	for (int i = 0; i < gcp.rows(); i++)
	{
		A(2 * i+1, 0) = gcp(i, 1);
		A(2 * i+1, 1) = -gcp(i, 0);
		A(2 * i+1, 2) = 0;
		A(2 * i+1, 3) = 1;
		A(2 * i, 0) = gcp(i, 0);
		A(2 * i, 1) = gcp(i, 1);
		A(2 * i, 2) = 1;
		A(2 * i, 3) = 0;
	}
	return A;
}

MatrixXd lmat(MatrixXd& meas)
{
	MatrixXd l(meas.rows() * 2, 1);
	for (int i = 0; i < meas.rows(); i++)
	{
		l(2 * i + 1, 0) = meas(i, 1);
		l(2 * i, 0) = meas(i, 0);
	}
	return l;
}

double z(MatrixXd& meas)
{
	double z = 0;
	for (int i = 0; i < meas.rows(); i++)
	{
		z = z + meas(i, 2);
		
	}
	z = z / meas.rows();
	return z;
}