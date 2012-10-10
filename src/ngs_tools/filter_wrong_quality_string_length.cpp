#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/mman.h>
#include <time.h>
#include <assert.h>
#include <signal.h>
#define MAXLINE 1000


int main(int argc, char* argv[])
{
	FILE* pFile = stdin;

	char* line = new char[MAXLINE];
	char* seq = new char[MAXLINE];
	char* qual = new char[MAXLINE];
	char* ret = new char[MAXLINE];

	int i;

	while (!feof(pFile))
	{
		i++;
		ret = fgets(line, MAXLINE, pFile);
	
		if (ret==NULL)
		{
			fprintf(stderr, "number of lines read: %d\n", i);
			break;	
		}

		int num_read = sscanf(line, "%*s%*d%*s%*d%*d%*s%*s%*d%*d%s%s", seq, qual);

		if (num_read!=2)
		{
			fprintf(stderr, "error: could not parse line: \n%s\n", line);
			continue;
		}

		if (strlen(seq)==strlen(qual) && num_read==2)
		{
			printf("%s", line);
		}
		else
		{
			//printf("filter_out %s", line);
			fprintf(stderr, "num_read: %d\n", num_read);
			fprintf(stderr, " seq: %s len: %d\nqual: %s len: %d\n", seq, strlen(seq), qual, strlen(qual));
			fprintf(stderr, "%s\n", line);
		}
	}
}
