#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <time.h>


int main(int argc, char** args)
{
	char line [1000];
	srand(time(NULL));
	float subs;

	long int cnt=0;
	while (fgets(line, 1000, stdin) != NULL)
	{
		fprintf(stdout, ">seq%li\n", cnt++);
		fprintf(stdout, "%s", line);
	}
return 0;
}
