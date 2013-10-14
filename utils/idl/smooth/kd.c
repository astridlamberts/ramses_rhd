#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include "kd.h"
#include "tipsydefs.h"


#define MAX_ROOT_ITTR	32


void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;

	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}


int kdInit(KD *pkd,int nBucket)
{
	KD kd;

	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	kd->p = NULL;
	kd->kdNodes = NULL;
	*pkd = kd;
	return(1);
	}


void kdascii2binary(FILE *infile,FILE *outfile)
{
    int ndim ;
    int nbodies ;
    int ngas ;
    int ndark ;
    int nstar ;
    int count ;
    struct gas_particle *gp, *lastgp;
    struct dark_particle *dp, *lastdp ;
    struct star_particle *sp, *lastsp ;

    forever {
	count=fscanf(infile, "%d%*[, \t\n]%d%*[, \t\n]%d"
		    ,&header.nbodies, &header.nsph, &header.nstar) ;
	if ( (count == EOF) || (count==0) ){
	    break ;
	}
	fscanf(infile,"%d",&header.ndim) ;
	fscanf(infile,"%lf",&header.time) ;
	ndim=header.ndim;
	nbodies=header.nbodies;
	ngas=header.nsph;
	nstar = header.nstar ;
	ndark = header.ndark = nbodies - nstar - ngas ;
	if(gas_particles != NULL) free(gas_particles);
	if(ngas != 0) {
	    gas_particles =
		(struct gas_particle *) malloc(ngas*sizeof(*gas_particles));
	    if(gas_particles == NULL) {
		printf("<sorry, no memory for gas particles, %s>\n") ;
		return ;
	    }
	}
	else
	  gas_particles = NULL;
	if(dark_particles != NULL) free(dark_particles);
	if(ndark != 0) {
	    dark_particles =
		(struct dark_particle *) malloc(ndark*sizeof(*dark_particles));
	    if(dark_particles == NULL) {
		printf("<sorry, no memory for dark particles, %s>\n") ;
		return ;
	    }
	}
	else
	  dark_particles = NULL;
	if(star_particles != NULL) free(star_particles);
	if(nstar != 0) {
	    star_particles =
		 (struct star_particle *)malloc(nstar*sizeof(*star_particles));
	    if(star_particles == NULL) {
		printf("<sorry, no memory for star particles, %s>\n") ;
		return ;
	    }
	}
	else
	  star_particles = NULL;

	lastgp = gas_particles + ngas ;
	lastdp = dark_particles + ndark ;
	lastsp = star_particles + nstar ;

	for(gp=gas_particles; gp < lastgp ; gp++){
	    fscanf(infile,"%f%*[, \t\n]",&gp->mass);
	}
	for(dp=dark_particles; dp < lastdp;  dp++){
	    fscanf(infile,"%f%*[, \t\n]",&dp->mass);
	}
	for(sp=star_particles; sp < lastsp; sp++){
	    fscanf(infile,"%f%*[, \t\n]",&sp->mass);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(infile,"%f%*[, \t\n]",&gp->pos[0]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fscanf(infile,"%f%*[, \t\n]",&dp->pos[0]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(infile,"%f%*[, \t\n]",&sp->pos[0]);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(infile,"%f%*[, \t\n]",&gp->pos[1]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fscanf(infile,"%f%*[, \t\n]",&dp->pos[1]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(infile,"%f%*[, \t\n]",&sp->pos[1]);
	}
	if (ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++) {
		fscanf(infile,"%f%*[, \t\n]",&gp->pos[2]);
	    }
	    for(dp=dark_particles; dp < lastdp ; dp++) {
		fscanf(infile,"%f%*[, \t\n]",&dp->pos[2]);
	    }
	    for(sp=star_particles; sp < lastsp ; sp++) {
		fscanf(infile,"%f%*[, \t\n]",&sp->pos[2]);
	    }
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(infile,"%f%*[, \t\n]",&gp->vel[0]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
		fscanf(infile,"%f%*[, \t\n]",&dp->vel[0]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(infile,"%f%*[, \t\n]",&sp->vel[0]);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(infile,"%f%*[, \t\n]",&gp->vel[1]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fscanf(infile,"%f%*[, \t\n]",&dp->vel[1]);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(infile,"%f%*[, \t\n]",&sp->vel[1]);
	}
	if (ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++) {
		fscanf(infile,"%f%*[, \t\n]",&gp->vel[2]);
	    }
	    for(dp=dark_particles; dp < lastdp ; dp++) {
		fscanf(infile,"%f%*[, \t\n]",&dp->vel[2]);
	    }
	    for(sp=star_particles; sp < lastsp ; sp++) {
		fscanf(infile,"%f%*[, \t\n]",&sp->vel[2]);
	    }
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    fscanf(infile,"%f%*[, \t\n]",&dp->eps);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(infile,"%f%*[, \t\n]",&sp->eps);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(infile,"%f%*[, \t\n]",&gp->rho);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(infile,"%f%*[, \t\n]",&gp->temp);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(infile,"%f%*[, \t\n]",&gp->hsmooth);
	}
	for(gp=gas_particles; gp < lastgp ; gp++) {
	    fscanf(infile,"%f%*[, \t\n]",&gp->metals);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(infile,"%f%*[, \t\n]",&sp->metals);
	}
	for(sp=star_particles; sp < lastsp ; sp++) {
	    fscanf(infile,"%f%*[, \t\n]",&sp->tform);
	}
	for(gp=gas_particles; gp < lastgp ; gp++){
	    fscanf(infile,"%f%*[, \t\n]",&gp->phi);
	}
	for(dp=dark_particles; dp < lastdp;  dp++){
	    fscanf(infile,"%f%*[, \t\n]",&dp->phi);
	}
	for(sp=star_particles; sp < lastsp; sp++){
	    fscanf(infile,"%f%*[, \t\n]",&sp->phi);
	}
	fwrite((char *)&header,sizeof(header),1,outfile) ;
	fwrite((char *)gas_particles,sizeof(struct gas_particle),ngas,outfile) ;
	fwrite((char *)dark_particles,sizeof(struct dark_particle),
	       ndark,outfile) ;
	fwrite((char *)star_particles,sizeof(struct star_particle),
	   nstar,outfile) ;
	printf("read time %lf\n",header.time) ;
    }
}

int kdReadTipsy(KD kd,FILE *fp,int bDark,int bGas,int bStar)
{
	int i,j,nCnt;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;

	fread(&h,sizeof(struct dump),1,fp);
	kd->nParticles = h.nbodies;
	kd->nDark = h.ndark;
	kd->nGas = h.nsph;
	kd->nStar = h.nstar;
	kd->fTime = h.time;
	kd->nActive = 0;
	if (bDark) kd->nActive += kd->nDark;
	if (bGas) kd->nActive += kd->nGas;
	if (bStar) kd->nActive += kd->nStar;
	kd->bDark = bDark;
	kd->bGas = bGas;
	kd->bStar = bStar;
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);
	/*
	 ** Read Stuff!
	 */
	nCnt = 0;
	for (i=0;i<h.nsph;++i) {
		fread(&gp,sizeof(struct gas_particle),1,fp);
		if (bGas) {
			kd->p[nCnt].fMass = gp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = gp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = gp.vel[j];
			++nCnt;
			}
		}
	for (i=0;i<h.ndark;++i) {
		fread(&dp,sizeof(struct dark_particle),1,fp);
		if (bDark) {
			kd->p[nCnt].fMass = dp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = dp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = dp.vel[j];
			++nCnt;
			}
		}
	for (i=0;i<h.nstar;++i) {
		fread(&sp,sizeof(struct star_particle),1,fp);
		if (bStar) {
			kd->p[nCnt].fMass = sp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = sp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = sp.vel[j];
			++nCnt;
			}
		}
	return(kd->nParticles);
	}


void kdInMark(KD kd,char *pszFile)
{
	FILE *fp;
	char ach[80];
	int i,iCnt,iDum;

	fp = fopen(pszFile,"r");
	if (!fp) {
		fprintf(stderr,"Could not open mark array, %s\n",pszFile);
		exit(1);
		}
	fgets(ach,80,fp);	/* ignore the array header! */
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) fscanf(fp,"%d",&kd->p[iCnt++].iMark);
		else fscanf(fp,"%d",&iDum);
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) fscanf(fp,"%d",&kd->p[iCnt++].iMark);
		else fscanf(fp,"%d",&iDum);
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) fscanf(fp,"%d",&kd->p[iCnt++].iMark);
		else fscanf(fp,"%d",&iDum);
		}
	fclose(fp);
	}


void kdSelect(KD kd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	double v;
	int i,j;

	p = kd->p;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void kdCombine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	}


void kdUpPass(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		kdUpPass(kd,l);
		kdUpPass(kd,u);
		kdCombine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->p[u].r[j];
			c[iCell].bnd.fMax[j] = kd->p[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->p[pj].r[j];
				if (kd->p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->p[pj].r[j];
				}
			}
		}
	}


void kdBuildTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	n = kd->nActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes != NULL) free(kd->kdNodes);
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->p[0].r[j];
		bnd.fMax[j] = kd->p[0].r[j];
		}
	for (i=1;i<kd->nActive;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->p[i].r[j]) 
				bnd.fMin[j] = kd->p[i].r[j];
			else if (bnd.fMax[j] < kd->p[i].r[j])
				bnd.fMax[j] = kd->p[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelect(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->p[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m+1;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	kdUpPass(kd,ROOT);
	}


int cmpParticles(const void *v1,const void *v2)
{
	PARTICLE *p1=(PARTICLE *)v1,*p2=(PARTICLE *)v2;
	
	return(p1->iOrder - p2->iOrder);
	}


void kdOrder(KD kd)
{
	qsort(kd->p,kd->nActive,sizeof(PARTICLE),cmpParticles);
	}


void kdFinish(KD kd)
{
	free(kd->p);
	free(kd->kdNodes);
	free(kd);
	}

