#include "lrm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

int main(int argc, char** argv) {

	ifstream infile("salary_data.csv");
	vector<vector<string>> rows;

    string line;
    while (getline(infile, line)) {
        vector<string> row;

        string field;
        stringstream ss(line);

        while (getline(ss, field, ',')) {
            row.push_back(field);
        }

        rows.push_back(row);
    }
    vector<double> x;
    vector<double> y;
    for (int i = 1; i < rows.size(); i++) {
        x.push_back(stod(rows[i][0]));
        y.push_back(stod(rows[i][1]));
    }

    darray xD(x);
    darray yD(y);

    vector<double> results = slrm(xD, yD);
    cout << "\nBeta Value: " << results[0] << '\n';
    cout << "Intercept: " << results[1] << '\n';

	return 0;
}