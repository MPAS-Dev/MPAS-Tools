/*** MAIN.C ***/

#include <stdio.h>
#include <stdlib.h>  /* realloc(), qsort() */
#include "vdefs.h"
using namespace std;

Site * readone(void), * nextone(void) ;
void readsites(PointSet * p) ;

int sorted, triangulate, plot, debug, nsites, siteidx ;
double xmin, xmax, ymin, ymax ;
Site * sites ;
Freelist sfl ;

void
voronoi_main(PointSet * p)
{
	int c ;
	Site *(*next)() ;

	sorted = plot = debug = 0 ;
	triangulate = 1 ;

	freeinit(&sfl, sizeof(Site)) ;
	readsites(p) ;
	next = nextone ;
	siteidx = 0 ;
	geominit() ;
	voronoi(p, next) ;
	free_all();
}

/*** sort sites on y, then x, coord ***/

int
scomp(const void * vs1, const void * vs2)
{
	VPoint * s1 = (VPoint *)vs1 ;
	VPoint * s2 = (VPoint *)vs2 ;

	if (s1->y < s2->y)
		{
		return (-1) ;
		}
	if (s1->y > s2->y)
		{
		return (1) ;
		}
	if (s1->x < s2->x)
		{
		return (-1) ;
		}
	if (s1->x > s2->x)
		{
		return (1) ;
		}
	return (0) ;
}

/*** return a single in-storage site ***/

Site *
nextone(void)
{
	Site * s ;

	if (siteidx < nsites)
		{
		s = &sites[siteidx++];
		return (s) ;
		}
	else
		{
		return ((Site *)NULL) ;
		}
}

/*** read all sites, sort, and compute xmin, xmax, ymin, ymax ***/

void
readsites(PointSet * p)
{
	int i ;
	int j ;

	nsites = 0 ;
	sites = (Site *) myalloc(20000 * sizeof(Site));
	for(j=0; j<p->nPoints; j++) {
	   sites[nsites].coord.x = p->points[j]->getX();
	   sites[nsites].coord.y = p->points[j]->getY();
	   sites[nsites].sitenbr = p->points[j]->getNum() ;
	   sites[nsites++].refcnt = 0 ;
	   if (nsites % 20000 == 0) {
		  sites = (Site *)realloc(sites,(nsites+20000)*sizeof(Site));
	   }
	}

	qsort((void *)sites, nsites, sizeof(Site), scomp) ;
	xmin = sites[0].coord.x ;
	xmax = sites[0].coord.x ;
	for (i = 1 ; i < nsites ; ++i)
		{
		if(sites[i].coord.x < xmin)
			{
			xmin = sites[i].coord.x ;
			}
		if (sites[i].coord.x > xmax)
			{
			xmax = sites[i].coord.x ;
			}
		}
	ymin = sites[0].coord.y ;
	ymax = sites[nsites-1].coord.y ;
}

/*** read one site ***/

Site *
readone(void)
{
	Site * s ;

	s = (Site *)getfree(&sfl) ;
	s->refcnt = 0 ;
	s->sitenbr = siteidx++ ;
	if (scanf("%lf %lf", &(s->coord.x), &(s->coord.y)) == EOF)
		{
		return ((Site *)NULL ) ;
		}
	return (s) ;
}

