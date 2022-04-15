#include <iostream>
#include <vector>

#include "Array.h"
#include "ArrayUtils.h"

int main(int argc, char *argv[]){

	std::vector<int> size;
	size.push_back(1);
	DArray m(size);	

	print(m);

	return 0;
}
