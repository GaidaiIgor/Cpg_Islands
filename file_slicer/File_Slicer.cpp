#include <string>
#include <iostream>
#include <fstream>

using namespace std;


int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cout << "Not enough parameters" << endl;
		return 1;
	}

	ifstream file_for_slicing(argv[1]);
	if (!file_for_slicing.good())
	{
		cout << "Bad sliced file" << endl;
		return 2;
	}

	string start(argv[2]);
	string end(argv[3]);
	string line;
	
	bool out = false;
	while (file_for_slicing)
	{
		getline(file_for_slicing, line);

		if (line == end)
		{
			break;
		}

		if (out)
		{
			cout << line << endl;
		}

		if (line == start)
		{
			out = true;
		}	
	}

	return 0;
}

