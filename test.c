#include <iostream>
int x = 20;

namespace outer {
	int x = 10;
	namespace inner {
		int z = x;
	}
}

int main(){
std::cout << outer::inner::z;
getchar();
return 0;
}
