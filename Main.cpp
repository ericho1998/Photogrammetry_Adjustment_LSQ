#include "header.h"
#include <set>
#include <functional>
#include <algorithm>
int main()
{
	MatrixXd dist, ctrl, x, meas; // where x is the unknowns
	double xp, yp;
	xp = -0.069;
	yp = -0.104;
	double flyh = 1520;
	read("all_p.txt", ctrl); // file containing Ground Control Point coordinates in mm
	read("distortion.txt", dist); // file containing radial distortion corrections in micrometers
	read("obs.txt", meas);
	for (int i = 0; i < dist.rows(); i++)
	{
		dist(i, 1) = dist(i, 1) / 1000;
	}
	x.resize(ctrl.rows(), ctrl.cols());
	x = interpolation(dist, ctrl, x, xp, yp);
	vector<double> ID(x.rows(), 1);
	for (int i = 0; i < ID.size(); ++i)
	{
		ID.at(i) = x(i, 0);
	}
	vector<double>::iterator y;
	y = unique(ID.begin(), ID.end());
	ID.resize(distance(ID.begin(), y));
	MatrixXd obs(1,3), a, l,m(1,4), eop(ID.size(),7), param(ID.size(),4);
	
	for (int k = 0; k < ID.size(); k++)
	{
		double z0 = 0;
		for (int j = 0; j < x.rows(); j++)
		{
			if (x(j, 0) == ID.at(k) && x(j, 3) < 100)
			{
				//assign control points to a temporary matrix for inputs to the A matrix calculation
				obs(obs.rows() - 1, 0) = x(j, 1);
				obs(obs.rows() - 1, 1) = x(j, 2);
				obs(obs.rows() - 1, 2) = x(j, 3);
				obs.conservativeResize(obs.rows() + 1, obs.cols());
			}	
			
			
		}
		obs.conservativeResize(obs.rows() - 1, obs.cols());
		
		for (int j = 0; j < obs.rows(); j++)
		{
			for (int l = 0; l < meas.rows(); l++)
			{
				if (meas(l, 0) == obs(j, 2))
				{
					m(m.rows() - 1, 0) = meas(l, 1);
					m(m.rows() - 1, 1) = meas(l, 2);
					m(m.rows() - 1, 2) = meas(l, 3);
					m(m.rows() - 1, 3) = meas(l, 0);
					m.conservativeResize(m.rows() + 1, m.cols());
					if (m.rows() > obs.rows())
					{
						break;
					}

				}

			}
		}
		m.conservativeResize(m.rows() - 1, m.cols());
		
		a = Amat(obs);
		l = lmat(m);
		
		z0 = z(m);
		MatrixXd xhat = (a.transpose() * a);
		xhat = xhat.inverse() * a.transpose() * l;
		// get EOPs (External Orientation Parameters) for the image
		eop(k, 0) = ID.at(k);
		eop(k, 1) = 0;
		eop(k, 2) = 0;
		eop(k, 3) = atan2(xhat(1, 0), xhat(0, 0)) * 180 / acos(-1);
		/*if (eop(k, 3) < 0)
		{
			eop(k, 3) = eop(k, 3) + 360;
		}
		else if (eop(k, 3) > 360)
		{
			eop(k, 3) = eop(k, 3) - 360;
		}*/
		eop(k, 4) = xhat(2, 0);
		eop(k, 5) = xhat(3, 0);
		eop(k, 6) = z0 + flyh;
		// put parameters into a single matrix
		param(k, 0) = xhat(0, 0);
		param(k, 1) = xhat(1, 0);
		param(k, 2) = xhat(2, 0);
		param(k, 3) = xhat(3, 0);
		for (int b = 0; b < x.rows(); b++)
		obs.resize(1, 3);
		m.resize(1, 4);

	}
	// perform similarity transform to all points
	MatrixXd x2(x.rows(), x.cols()+1);
	x2.block(0,0,x.rows(), x.cols())= x;
	
	for (int i = 0; i < ID.size(); i++)
	{
		for (int j = 0; j < x.rows(); j++)
		{
			if (ID.at(i) == x(j, 0))
			{
				x2(j, 1) = param(i, 0) * x(j, 1) + param(i, 1) * x(j, 2) + param(i, 2);
				x2(j, 2) = -param(i, 1) * x(j, 1) + param(i, 0) * x(j, 2) + param(i, 3);
				x2(j, 4) = eop(i, 6) - flyh;
			}
		}
	}
	
	// sort the points and average them for the gcp file (condense into a function)
	vector<double> tie(x2.rows(),1);
	for (int i = 0; i < x2.rows(); i++)
	{
		tie.at(i) = x2(i, 3);
	}
	set<double> duplicates;
	unsigned size = tie.size();
	for (unsigned int i = 0; i < size; i++)
	{
		duplicates.insert(tie[i]);
	}
	tie.assign(duplicates.begin(), duplicates.end());

	MatrixXd x3(tie.size(),4);
	
	for (int i = 0; i < tie.size(); i++)
	{
		double x_avg = 0, y_avg = 0, z_avg = 0, c = 0, n = 0;
		for (int j = 0; j < x2.rows(); j++)
		{
			if (tie.at(i) == x2(j, 3) && x2(j,3) > 100)
			{
				x_avg = x_avg + x2(j, 1);
				y_avg = y_avg + x2(j, 2);
				z_avg = z_avg + x2(j, 4);
				n = x2(j, 3);
				c = c + 1;
				continue;
			}
			else if (tie.at(i) == x2(j,3) && x2(j,3) < 100)
			{
				for (int k = 0; k < meas.rows(); k++)
				{
					if (meas(k, 0) == x2(j, 3))
					{
						x_avg = meas(k, 1);
						y_avg = meas(k, 2);
						z_avg = meas(k, 3);
						n = x2(j, 3);
						continue;
					}
				}
			}
		}	
		if (tie.at(i) > 100)
		{
			x_avg = x_avg / c;
			y_avg = y_avg / c;
			z_avg = z_avg / c;
		}
		x3(i, 0) = x_avg;
		x3(i, 1) = y_avg;
		x3(i, 2) = z_avg;
		x3(i, 3) = n;
	}
	
	//write the files for MSAT
	writeicf("Lab1.icf", x);
	writegcp("Lab1.gcp", x3);
	writeori("Lab1.ori", eop);
	writecam("Lab1.cam");
	cout << "No compilation errors here! #bless" << endl;
	system("pause");
	return 0;
}