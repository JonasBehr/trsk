#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <time.h>


int main(int argc, char** args)
{
	char line [1000];
	srand(time(NULL));
	float subs;
	if (argc>1)
		subs = atof(args[1]);
	else
	{
		fprintf(stderr, "expected subsample factor as arg1\n");
		exit(1);
	}

	while (fgets(line, 1000, stdin) != NULL)
	{
		float r = 1.0*rand()/RAND_MAX;

		if (r<subs)
			fprintf(stdout, "%s", line);
	}
return 0;
}
