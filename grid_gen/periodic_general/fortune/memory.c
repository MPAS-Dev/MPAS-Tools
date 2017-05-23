
/*** MEMORY.C ***/

#include <stdio.h>
#include <stdlib.h>  /* malloc(), exit() */

#include "vdefs.h"

extern int sqrt_nsites, siteidx ;
char** memory_map;
int nallocs = 0;

void
freeinit(Freelist * fl, int size)
{
	fl->head = (Freenode *)NULL ;
	fl->nodesize = size ;
}

char *
getfree(Freelist * fl)
{
	int i ;
	Freenode * t ;
	if (fl->head == (Freenode *)NULL)
		{
		t =  (Freenode *) myalloc(sqrt_nsites * fl->nodesize) ;
		for(i = 0 ; i < sqrt_nsites ; i++)
			{
			makefree((Freenode *)((char *)t+i*fl->nodesize), fl) ;
			}
		}
	t = fl->head ;
	fl->head = (fl->head)->nextfree ;
	return ((char *)t) ;
}

void
makefree(Freenode * curr, Freelist * fl)
{
	curr->nextfree = fl->head ;
	fl->head = curr ;
}

int total_alloc;

char *
myalloc(unsigned n)
{
	char * t ;
	if ((t=(char*)malloc(n)) == (char *) 0)
		{
		fprintf(stderr,"Insufficient memory processing site %d (%d bytes in use)\n",
		siteidx, total_alloc) ;
		exit(0) ;
		}
	total_alloc += n ;

	if (nallocs % 20000 == 0)
	{
		if (nallocs == 0)
			memory_map = (char **)malloc((nallocs+20000)*sizeof(char*));
		else
			memory_map = (char **)realloc(memory_map,(nallocs+20000)*sizeof(char*));
	}
	memory_map[nallocs++] = t;
	return (t) ;
}

void free_all(void)
{
	int i;

	for (i=0; i<nallocs; i++)
	{
		if (memory_map[i] != (char*)0)
		{
			free(memory_map[i]);
			memory_map[i] = (char*)0;
		}
	}

	free(memory_map);
	nallocs = 0;
}
