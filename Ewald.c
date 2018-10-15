
/*
 * zone = 0 if atom is in cluster.
 * zone = 1 if atom is in a non-moveable unit cell.
 * zone = 2 if atom is in a moveable unit cell.
 * zone = 3 if atom is in an initially.
 * empty unit cell.
 *
 */

/* General Ewald Potential program */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef double doublereal;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int integer;

#define IARITHIF(x)  (x > 0 ? 1 : x < 0 ? -1 : 0)

struct t_fgercm {
        integer ierr, ierct;
}       fgercm;

#define MAXPT 120000
#define MAXUC 10000
#define MAXMAD 10000
#define MAXTYPE 2100
#define MAX2N 20
#define MXPAR 42000 /* 35000 */
#define MXEQU 2100  /* 1100 */

FILE *ptsfilenw,*ptsfile,*ptsfilefro,*outfile,*ucellinfile,*rmspotfile,*seedfile,*pts,*defectfiler;
FILE *clustfile,*potjamfil,*ucelloutfile,*ucellstdfile,*listfile,*pubfile;

static double avepotuc2, rmspotuc2;
static double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z;
static double b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z;
static double uc1[MAXUC], uc2[MAXUC], uc3[MAXUC];
static double quc[MAXUC], ucx[MAXUC], ucy[MAXUC], ucz[MAXUC];
static integer itemp, nsites;
static double uctmp1, uctmp2, uctmp3,radmom;
static double energy, energy1;
static double xmad[MAXMAD], ymad[MAXMAD], zmad[MAXMAD], qmad[MAXMAD];
static double dspot[MAXMAD], dspotS[MAXMAD], dspotFA[MAXMAD];
static double rmsewalda[MAXMAD];
static double dmin,scalemad,r2, ran;
static double rmspot,avepot;
static integer ipt, iptmax, ifill, imad, imadmax;
static double eshellm[20], eshellr[20], ewaldmad, evconv;
static integer ir1, ir2, ir3, irmax, ishell;
static double xptr, yptr, zptr, qptr,rmspot2;
static integer in1, in2,maxpts;
static integer nrandpts,nrandpts_ck,nomad, nave;
static integer it1, it2, it3, itr1, itr2, itr3;
static integer in, inmax,itr11, itr12, itr13;
static integer it, itmax, itype[MAXUC], ntype[MAXTYPE], nuctype[MAXTYPE];
static char typetab[MAXTYPE][6], typuc[MAXUC][6];
static integer input_flag, ifatom,ifmatch;
static double qsum, d,temp, temp1, temp2, temp3;
static double xmin, ymin, zmin, xmax, ymax, zmax;
static integer m1, m2, m3, mmax;
static double ux, uy, uz, uq, rx, ry, rz, dot, x;
static double alpha, volume, vx, vy, vz, pi, tpi;
static double ewaldr[MAXMAD], ewaldm[MAXMAD];
static double ewaldc[MAXMAD], ewald[MAXMAD];
static integer iseed1, iseed2;
double xpt[MAXPT], ypt[MAXPT], zpt[MAXPT], qpt[MAXPT], ept[MAXPT];
double tept[MAXPT], sept[MAXPT], bept[MAXPT], delept[MAXPT];
integer zone[MAXPT], tzone[MAXPT], bzone[MAXPT], typept[MAXPT];
static double dx, dy, dz, qq, d2;
static double dxx, dyy, dz,dzz;
static double xsgrid[50], ysgrid[50], zsgrid[50];
static integer ipti[1001], ixyzs, ixs, iys, izs, nseg;
static double xsmin, xsmax, ysmin, ysmax, zsmin, zsmax, dsx, dsy, dsz;
static double xval, yval, zval;
static double xlmin, xlmax, ylmin, ylmax, zlmin, zlmax;

static integer nr1,nr2,nr3,nr4,nr5,mxzon2,zone2list[MAXPT],z2list[MAXPT];
integer nr_equ,nr_par,nr_e_eq,i;
static double qlist[100];
integer M,LWORK=10,j,jj,jp,RANK;
double *A,*B,*WORK,xx[10][10],dq,q0,ddq,*S,RCOND,*JPVT;
integer NRHS,LDA,LDB,INFO=0,N,*IPIV;
static double qtmp[11000],dspottmp[MAXMAD],dipx,dipy,dipz,qmin,qmax;
static integer varycharges,ignore,vs,nrclusteratoms,DEBUG,nrclusteratoms_rmax_rnd;
double itera,itermax,itera_rms,itera_dip,pot_clust_small;
static double dipr,diprtmp,dipx_org,dipxtmp,dipy_org,dipytmp,dipz_org,dipztmp,dipfact2,rmsfact2;
double *xxg,*xgstep,*dg,q,r;
double xmadr[MAXMAD], ymadr[MAXMAD], zmadr[MAXMAD], qmadr[MAXMAD];
double xval_tmp[11000],yval_tmp[11000],zval_tmp[11000],y,z,wtmadr[MAXMAD];
double dspotr[MAXMAD],rrms,rrms2,atom_r_rnd,rr;
char sss[80],labatm[MAXMAD][80];
integer imadmaxr,rnomad,nrchk,listan[MAXMAD],listan2[MAXMAD],n0,n1,n2,n3,n4,n5;
FILE *rmspotfiler;

static time_t now;
struct tm *date, *date1;
static double pot,xyz[MAXMAD][4];
static integer typxyz[MAXMAD], typmad[MAXMAD], typmada[MAXMAD];
static integer nn,n,nnn,cnr1,cnr2,nrformulaunitsincell,do_lapack5;
static double antala[MAXMAD],ewalda[MAXMAD];
static double contin,qmada[MAXMAD];
static double dspotSa[MAXMAD], dspotFAa[MAXMAD], dspota[MAXMAD];
static double rmsmadSa[MAXMAD], rmsmad1a[MAXMAD], rmsmad2a[MAXMAD];
static double rmscS, rmsc1, rmsc2, avecS, avec1, avec2;
static double rmsS, rms1, rms2, aveS, ave1, ave2;
static double aveSa[MAXMAD], ave1a[MAXMAD], ave2a[MAXMAD];
char fileroot[60], filermspot[60], fileucell[60], fileout[60], filepts[60],fileptsfro[60],fileptsnw[60],
     fileatoms[60], filexch[60], resname[60], fileclust[60], filelist[60], strdate1[80];

void randome (integer *, integer *, double *);
double erfc(double);


int main(){

        DEBUG = 0;  /*
                     * 0 -
                     * 1 - debug output
                     * 2 - debug output including parameter minimisation routin, 1 - lite debug output
                     */
        pi = 2.*asin(1.);
        tpi = 2.*pi;
        evconv = 14.3998;

        /* zone = 0 if atom is in cluster; zone = 1 if atom is in a non-moveable unit cell
           zone = 2 if atom is in a moveable unit cell; zone = 3 if atom is in an initially
           empty unit cell */


        printf("Enter root filename\n ");
        printf("     (press return to read default 'Epotfit.listx' file): ");
        scanf("%s",fileroot);
        /*printf("Enter number of random points among central 50+ atoms:");
           scanf("%ld", &nrandpts);gets(sss);*/
        nrandpts=0;
        printf("Enter number of random points among central 50+ atoms for check:");
        scanf("%ld", &nrandpts_ck); gets(sss);
        /*
           printf("Enter 1 to minimize rms using '4-equations' method. Best if starting with dip=0.Conserves dip.\n");
           printf("      2 to minimize rms using 'ladder' method. Alternates minimization of rms and dip.\n");
           printf("      3 first minimize dipole moment and then method 1.\n");
           scanf("%ld", &varycharges);gets(sss);
           printf(" varycharges = %ld\n", varycharges);
         */
        varycharges=3;
        /*
           printf("Enter max nr random atttempts for(1-3) above: ");
           scanf("%le", &itermax);gets(sss);
           printf("itermax = %e\n",itermax);
         */
        itermax=3.0e4;

        printf("Enter nr clusteratoms: ");
        scanf("%ld", &nrclusteratoms); gets(sss);
        printf("nrclusteratoms = %ld\n",nrclusteratoms);
        printf("Enter nr clusteratoms-rnd-chk: ");
        scanf("%ld", &nrclusteratoms_rmax_rnd); gets(sss);
        printf("nrclusteratoms_rmax_rnd = %ld\n",nrclusteratoms_rmax_rnd);
        /*
           printf("Enter max dist. from rmspot atoms for random check point selection :");
           scanf("%lf", &atom_r_rnd);gets(sss);
         */
        atom_r_rnd=2.5;

        strncpy(fileout,fileroot,50);
        strncat(fileout,".out",20);
        outfile=fopen(fileout, "w");
        if (outfile == NULL) {
                fprintf(outfile, " unable to open output file %s- exit\n", fileout);
                exit(1);
        }
        strncpy(fileucell,fileroot,50);
        strncat(fileucell,".uc",20);
        ucellinfile=fopen(fileucell, "r");
        if (ucellinfile == NULL) {
                fprintf(outfile, " unable to open unit cell input file %s- exit\n",fileucell);
                exit(1);
        }

        fprintf(outfile, "output file opened\n");
        fflush(outfile);

        time(&now);
        date1=localtime(&now);
        strncpy(strdate1, asctime(date1), 60);
        fprintf(outfile, "\n%s  %s",fileroot, asctime(date1));

        strncpy(filermspot,fileroot,50);
        strncat(filermspot,".qc",20);

        strncpy(filepts,fileroot,50);
        strncat(filepts,".pts-jag",20);
        strncpy(fileptsfro,fileroot,50);
        strncat(fileptsfro,".pts-fro",20);
        strncpy(fileptsnw,fileroot,50);
        strncat(fileptsnw,".pts-nw",20);
        ptsfilenw=fopen(fileptsnw,"w");
        ptsfile=fopen(filepts, "w");
        ptsfilefro=fopen(fileptsfro,"w");
        if (ptsfile == NULL) {
                fprintf(outfile, " unable to open output point charge file %s- exit\n", filepts);
                exit(1);
        }

        maxpts = MAXPT;

        seedfile = fopen("seedfile","r");
        fscanf(seedfile, "%ld %ld", &iseed1, &iseed2);
        fclose(seedfile);

        xxg    = (double*) calloc(100, sizeof(double));
        xgstep = (double*) calloc(100, sizeof(double));
        dg     = (double*) calloc(100, sizeof(double));

        time(&now);
        date=localtime(&now);

        fprintf(outfile, " atom cluster for rmspot:\n");
        xmin=1.E20;
        ymin=1.E20;
        zmin=1.E20;
        xmax=-1.E20;
        ymax=-1.E20;
        zmax=-1.E20;
        for (imad=0; imad<imadmax; imad++) {
                if (xmad[imad] < xmin) xmin = xmad[imad];
                if (ymad[imad] < ymin) ymin = ymad[imad];
                if (zmad[imad] < zmin) zmin = zmad[imad];
                if (xmad[imad] > xmax) xmax = xmad[imad];
                if (ymad[imad] > ymax) ymax = ymad[imad];
                if (zmad[imad] > zmax) zmax = zmad[imad];
        }
        if (DEBUG > 0) {
                for (imad=0; imad<imadmax; imad++) {
                        fprintf(outfile, " imad %ld x %.5f y %.5f z %.5f q %.4f\n",
                                imad, xmad[imad], ymad[imad], zmad[imad], qmad[imad]);
                }
        } /* end if DEBUG > 0 */

        fscanf(ucellinfile, " %lf %lf %lf %ld", &a1x, &a1y, &a1z, &itr1);
        fscanf(ucellinfile, " %lf %lf %lf %ld", &a2x, &a2y, &a2z, &itr2);
        fscanf(ucellinfile, " %lf %lf %lf %ld", &a3x, &a3y, &a3z, &itr3);
        fflush(ucellinfile);

        /* read unit cell coordinates in principal crystal axis units */
        for (in=0; in < MAXUC-2; in++) {
                input_flag=fscanf(ucellinfile," %lf %lf %lf %lf %s",
                                  &uc1[in], &uc2[in], &uc3[in], &quc[in], typuc[in]);
                if (input_flag < 4) break;
        }
        inmax = in;

        itr11 = itr1;
        itr12 = itr2;
        itr13 = itr3;

        for (it=0; it<MAXTYPE; it++) {
                nuctype[it]=0;
                for (nnn=0; nnn<5; nnn++) typetab[it][nnn]='\0';
        }
        for (nnn=0; nnn<5; nnn++) typetab[0][nnn]=typuc[0][nnn];
        itype[0]=0;
        nuctype[0]=1;
        itmax=1;
        for (in=1; in<inmax; in++) {
                for (it=0; it<itmax; it++) {
                        ifmatch=1;
                        for (nnn=0; nnn<5; nnn++) {
                                if (typetab[it][nnn] != typuc[in][nnn]) {
                                        ifmatch = 0;
                                        break;
                                }
                                if ( (typetab[it][nnn] == '\0') && (typuc[in][nnn]=='\0') ) break;
                        } /* end nnn loop over string characters*/
                        if (ifmatch == 1) break;
                } /* end loop over previously found types */
                if (ifmatch == 0) {
                        for (nnn=0; nnn<5; nnn++) typetab[itmax][nnn]=typuc[in][nnn];
                        nuctype[itmax]=1;
                        itype[in]=itmax;
                        itmax++;
                }
                if (ifmatch ==1)
                        nuctype[it] += 1;
                itype[in] = it;
        } /* end in loop over unit cell atoms */
        for (nnn=0; nnn<4; nnn++) typetab[itmax][nnn]='?';
        typetab[itmax][5]='\0';

        for(it=0; it<itmax; it++) fprintf(outfile, "Unit cell has %ld atoms of type %s\n",
                                          nuctype[it], typetab[it]);

        fprintf(outfile, " Translation vectors:\n");
        fprintf(outfile, " a1x %.5f a1y %.5f a1z %.5f\n", a1x, a1y, a1z);
        fprintf(outfile, " a2x %.5f a2y %.5f a2z %.5f\n", a2x, a2y, a2z);
        fprintf(outfile, " a3x %.5f a3y %.5f a3z %.5f\n", a3x, a3y, a3z);

        fprintf(outfile, " Unit cell coordinates in crystal principal directions and charges:\n");
        for (in =0; in < inmax; in++)
                fprintf(outfile, " in %ld uc1 %f uc2 %f uc3 %f q %.4f %s\n",in, uc1[in], uc2[in], uc3[in], quc[in], typuc[in]);
        fflush(outfile);
        fprintf(outfile, "\n Unit cell x, y, z coordinates and charges:\n");

        dx = 0.;  dy = 0.;  dz = 0.; qq =0.;
        dxx =0.;  dyy =0.;  dzz =0.;
        radmom = 0.;
        for (in = 0; in<inmax; in++) {
                ucx[in] = uc1[in] * a1x + uc2[in] * a2x +uc3[in] * a3x;
                ucy[in] =uc1[in] * a1y + uc2[in] * a2y + uc3[in] * a3y;
                ucz[in] = uc1[in] * a1z + uc2[in] * a2z + uc3[in] * a3z;
                dx += ucx[in] * quc[in];
                dxx += ucx[in] * ucx[in] * quc[in];
                dy += ucy[in] * quc[in];
                dyy += ucy[in] * ucy[in] * quc[in];
                dz += ucz[in] * quc[in];
                dzz += ucz[in] * ucz[in] * quc[in];
                qq += quc[in];
                fprintf(outfile, " in %ld ucx %.5f ucy %.5f ucz %.5f q %.4f %s\n",
                        in, ucx[in], ucy[in], ucz[in], quc[in], typuc[in]);
        }
        d2 = dx*dx + dy*dy + dz*dz;
        radmom = sqrt(dxx*dxx + dyy*dyy + dzz*dzz);
        fflush(outfile);
        energy1 =0.;
        for (in1 =1; in1<inmax; in1++) {
                for (in2=0; in2<in1; in2++) {
                        r2 = (ucx[in1]-ucx[in2]) * (ucx[in1]-ucx[in2]) +
                             (ucy[in1]-ucy[in2]) * (ucy[in1]-ucy[in2]) +
                             (ucz[in1]-ucz[in2]) * (ucz[in1]-ucz[in2]);
                        if (r2 > 0.001)
                                energy1 += quc[in1]*quc[in2]/sqrt(r2);
                }
        }
        energy1 *= evconv;


        fprintf(outfile, " initial unit cell moments:\n");
        fprintf(outfile, "   dx %.8f  dy %.8f  dz %.8f  d2 %.8f\n",dx,dy,dz,d2);
        fprintf(outfile, "   dxx %.8f  dyy %.8f  dzz %.8f  radmom %.8f\n",dxx,dyy,dzz,radmom);
        fprintf(outfile, "  unit cell: energy %f\n",energy1);
        fprintf(outfile, " initial unit cell charge %f\n", qq);
        fflush(outfile);
        vx = a2y*a3z - a2z*a3y;
        vy = a2z*a3x - a2x*a3z;
        vz = a2x*a3y - a2y*a3y;
        volume = a1x*vx + a1y*vy + a1z*vz;
        fprintf(outfile, "volume of unit cell = %f �^3\n", volume);

        /* compute reciprical lattice vectors */
        b1x = a2y*a3z - a2z*a3y;
        b1y = a2z*a3x - a2x*a3z;
        b1z = a2x*a3y - a2y*a3x;
        temp = a1x*b1x + a1y*b1y + a1z*b1z;
        b1x /= temp;
        b1y /= temp;
        b1z /= temp;

        b2x = a1y*a3z - a1z*a3y;
        b2y = a1z*a3x - a1x*a3z;
        b2z = a1x*a3y - a1y*a3x;
        temp = a2x*b2x + a2y*b2y + a2z*b2z;
        b2x /= temp;
        b2y /= temp;
        b2z /= temp;

        b3x = a1y*a2z - a1z*a2y;
        b3y = a1z*a2x - a1x*a2z;
        b3z = a1x*a2y - a1y*a2x;
        temp = a3x*b3x + a3y*b3y + a3z*b3z;
        b3x /= temp;
        b3y /= temp;
        b3z /= temp;

        fprintf(outfile, "Reciprical lattice vectors\n");
        fprintf(outfile, " b1x %.5f b1y %.5f b1z %.5f\n", b1x, b1y, b1z);
        fprintf(outfile, " b2x %.5f b2y %.5f b2z %.5f\n", b2x, b2y, b2z);
        fprintf(outfile, " b3x %.5f b3y %.5f b3z %.5f\n", b3x, b3y, b3z);

        /* compute minimum inter atomic distance in unit cell*/
        dmin = 1000000.;
        for (in1=0; in1<inmax; in1++) {
                for (in2=in1+1; in2<inmax; in2++) {
                        d = (ucx[in1]-ucx[in2]) *  (ucx[in1]-ucx[in2]) +
                            (ucy[in1]-ucy[in2]) *  (ucy[in1]-ucy[in2]) +
                            (ucz[in1]-ucz[in2]) *  (ucz[in1]-ucz[in2]);
                        d = sqrt(d);
                        if (d < dmin) dmin = d;
                } /* end in2 loop */
        } /* end in1 loop */

        for (in1=0; in1<inmax; in1++) {
                for (in2=in1+1; in2<inmax; in2++) {
                        d = (ucx[in1]-ucx[in2]) *  (ucx[in1]-ucx[in2]) +
                            (ucy[in1]-ucy[in2]) *  (ucy[in1]-ucy[in2]) +
                            (ucz[in1]-ucz[in2]) *  (ucz[in1]-ucz[in2]);
                        d = sqrt(d);
                        if (d < 1.05*dmin) {
                                fprintf(outfile, " in1 %ld %s in2 %ld %s distance %f\n",
                                        in1, typuc[in1], in2,typuc[in2],  d);
                        }
                } /* end in2 loop */
        } /* end in1 loop */
        fflush(outfile);
        nsites = 2 * itr11 * 2 * itr12 * 2 * itr13 * inmax;
        fprintf(outfile, " compute xyz coordinates of lattice of %ld sites\n",nsites);

        /* find minimum and maximum x,y,z values of entire lattice */
        xlmin = 1.E20;
        xlmax = -1.E20;
        ylmin = 1.E20;
        ylmax = -1.E20;
        zlmin = 1.E20;
        zlmax = -1.E20;
        for (it1 = -itr11; it1 < itr11; it1++) { /* loop over first translation vector */
                for (it2 =-itr12; it2 < itr12; it2++) { /* loop over second translation vector */
                        for (it3 = -itr13; it3 < itr13; it3++) { /* loop over third translation vector */
                                for (in=0; in<inmax; in++) { /* loop over atoms in unit cell */
                                        uctmp1 = uc1[in] + it1;
                                        uctmp2 = uc2[in] + it2;
                                        uctmp3 = uc3[in] + it3;
                                        xval = uctmp1 * a1x + uctmp2 * a2x + uctmp3 * a3x;
                                        yval = uctmp1 * a1y + uctmp2 * a2y + uctmp3 * a3y;
                                        zval = uctmp1 * a1z + uctmp2 * a2z + uctmp3 * a3z;
                                        if (xval < xlmin) xlmin = xval;
                                        if (xval > xlmax) xlmax = xval;
                                        if (yval < ylmin) ylmin = yval;
                                        if (yval > ylmax) ylmax = yval;
                                        if (zval < zlmin) zlmin = zval;
                                        if (zval > zlmax) zlmax = zval;
                                } /* end in loop */
                        } /* end it3 loop */
                } /* end it2 loop */
        } /* end it1 loop */

        nseg = 8;
        dsx = 1.01*(xlmax-xlmin)/nseg;
        dsy = 1.01*(ylmax-ylmin)/nseg;
        dsz = 1.01*(zlmax-zlmin)/nseg;

        if(DEBUG > 0) {
                fprintf(outfile, "\n xlmin  %f   xlmax %f     dsx %f\n", xlmin, xlmax, dsx);
                fprintf(outfile, " ylmin  %f   ylmax %f    dsy %f\n", ylmin, ylmax, dsy);
                fprintf(outfile, " zlmin  %f   zlmax %f    dsz %f\n", zlmin, zlmax, dsz);
        }
        ipt=0;
        ipti[0]=0;
        for (ixs =0; ixs<=nseg; ixs++)
                xsgrid[ixs] = 1.01*xlmin + dsx*ixs;
        for (iys =0; iys<=nseg; iys++)
                ysgrid[iys] = 1.01*ylmin + dsy*iys;
        for (izs =0; izs<=nseg; izs++)
                zsgrid[izs] = 1.01*zlmin + dsz*izs;

        /* loop over x,y,z volume segments */
        for (ixs =0; ixs<nseg; ixs++) {
                xsmin = xsgrid[ixs];
                xsmax = xsgrid[ixs+1];
                for (iys =0; iys<nseg; iys++) {
                        ysmin = ysgrid[iys];
                        ysmax = ysgrid[iys+1];
                        for (izs =0; izs<nseg; izs++) {
                                zsmin = zsgrid[izs];
                                zsmax =zsgrid[izs+1];
                                ixyzs = ixs*nseg*nseg + iys*nseg +izs;
                                if (ixyzs > 1000) {
                                        fprintf(outfile, "error -EXIT  ixys %ld\n",ixyzs);
                                        exit(1);
                                }
                                /* loop over all unit cells- determine which atoms are in the
                                   volume segment being considered */
                                for (it1 = -itr11; it1 < itr11; it1++) { /* loop over first translation vector */
                                        for (it2 =-itr12; it2 < itr12; it2++) { /* loop over second translation vector */
                                                for (it3 = -itr13; it3 < itr13; it3++) { /* loop over third translation vector */
                                                        for (in=0; in<inmax; in++) { /* loop over atoms in unit cell */
                                                                uctmp1 = uc1[in] + it1;
                                                                uctmp2 = uc2[in] + it2;
                                                                uctmp3 = uc3[in] + it3;
                                                                xval = uctmp1 * a1x + uctmp2 * a2x + uctmp3 * a3x;
                                                                yval = uctmp1 * a1y + uctmp2 * a2y + uctmp3 * a3y;
                                                                zval = uctmp1 * a1z + uctmp2 * a2z + uctmp3 * a3z;
                                                                if ( (xval >= xsmin) && (xval < xsmax) &&
                                                                     (yval >= ysmin) && (yval < ysmax) &&
                                                                     (zval >= zsmin) && (zval < zsmax) ) { /*this atom is in volume segment ixs, iys, izs */
                                                                        xpt[ipt] = xval;
                                                                        ypt[ipt] = yval;
                                                                        zpt[ipt] = zval;
                                                                        qpt[ipt] = quc[in];
                                                                        typept[ipt] = itype[in];
                                                                        zone[ipt] =2; /*next to last layer, initially occupied sites */
                                                                        if ( (it1 > -itr11+1) && (it1 < itr11-2) &&
                                                                             (it2 > -itr12+1) && (it2 < itr12-2) && (it3 > -itr13+1) && (it3 < itr13-2) )
                                                                                zone[ipt] = 1;  /* interior, perminantly occupied sites */
                                                                        if ( (it1 == -itr11) || (it1 == itr11-1) ||
                                                                             (it2 == -itr12) || (it2 == itr12-1) || (it3 == -itr13) || (it3 == itr13-1) )
                                                                                zone[ipt] = 2;  /* org 3... outermost layer, initially vacant sites */
                                                                        if (DEBUG > 0) {
                                                                                fprintf(outfile, "ipt %ld x %.5f y %.5f z %.5f q %.4f zone %ld\n",
                                                                                        ipt, xpt[ipt], ypt[ipt], zpt[ipt], qpt[ipt], zone[ipt]);
                                                                        }
                                                                        ipt++;
                                                                        if (ipt >= maxpts) {
                                                                                fprintf(outfile, " error- ipt %ld    maxpts %ld - exit\n",
                                                                                        ipt, maxpts);
                                                                                exit(1);
                                                                        }
                                                                } /*end if x,y,z in volume segment test */
                                                        } /* end in loop */
                                                } /* end it3 loop */
                                        } /* end it2 loop */
                                } /* end it1 loop */
                                ixyzs = ixs*nseg*nseg + iys*nseg +izs;
                                ipti[ixyzs+1] = ipt;
                                if (DEBUG > 0) {
                                        fprintf(outfile, "ixs %ld  iys %ld  izs %ld  ixyzs %ld  ipti %ld\n",
                                                ixs, iys, izs, ixyzs, ipti[ixyzs]);
                                        fprintf(outfile,
                                                "    xsmin %f  xsmax %f  ysmin %f  ysmax %f  zsmin %f  zsmax %f\n",
                                                xsmin, xsmax, ysmin, ysmax, zsmin, zsmax);
                                } /* end case DEBUG > 0 */
                        } /* end izs loop */
                } /* end iys loop */
        } /* end ixs loop */

        iptmax = ipt;
        fflush(outfile);

        qsum=0.;
        for (ipt=0; ipt<iptmax; ipt++) {
                if (zone[ipt] < 3)
                        qsum += qpt[ipt];
        }

        if (fabs(qsum)>10e-7) {
                fprintf (outfile, "FATAL error- LATTICE CHARGE %f\n", qsum);
                fflush(outfile);
                exit(1);
        }

        /* Fixa .rmspot har...*/

        /* If nrclusteratoms_rmax_rnd != 0 generate a spherical cluster with at least nrclusteratoms inside. */
        /* Otherwise use .rmspot file as given. */
        /* Note that this is only important for the check since the fitting is done on */
        /* the nrclusteratom atoms. */

        rmspotfile=fopen("tmpfile", "w");
        for(pot=1.0; pot<=100.0; pot+=0.01) {                            /* Make rmspot file.*/
                nn=0;
                for (ipt=0; ipt< iptmax; ipt++) {
                        r2 = xpt[ipt]*xpt[ipt] + ypt[ipt]*ypt[ipt] + zpt[ipt]*zpt[ipt];
                        if (r2 > 0.000001) r2 = sqrt(r2);
                        if(r2<pot) nn++;
                }
                if(nn>=nrclusteratoms) break;
        }

        nn=0;
        for (ipt=0; ipt< iptmax; ipt++) {
                r2 = xpt[ipt]*xpt[ipt] + ypt[ipt]*ypt[ipt] + zpt[ipt]*zpt[ipt];
                if (r2 > 0.000001) r2 = sqrt(r2);
                if (r2 < pot) fprintf(rmspotfile,"%f\t%f\t%f\t%f\n",xpt[ipt],ypt[ipt],zpt[ipt],qpt[ipt]);
        }
        fclose(rmspotfile);
        rmspotfile=fopen("tmpfile", "r");                        /* .rmspot ok */

        /*
           Make  xmadr[] ... for the random chk...
           Only if no quantum qluster is given. (.qc file)
         */
        if(nrclusteratoms_rmax_rnd != 0) {                        /* mad for chk points...*/
                for(pot=1.0; pot<=100.0; pot+=0.01) {
                        nn=0;
                        for (ipt=0; ipt< iptmax; ipt++) {
                                r2 = xpt[ipt]*xpt[ipt] + ypt[ipt]*ypt[ipt] + zpt[ipt]*zpt[ipt];
                                if (r2 > 0.000001) r2 = sqrt(r2);
                                if(r2<pot) nn++;
                        }
                        if(nn<=nrclusteratoms_rmax_rnd) {
                                pot_clust_small=pot;
                                imadmaxr=nn;
                                if(DEBUG >0) fprintf(outfile,"-->imadmaxr=%ld   %f  %ld\n",
                                                     imadmaxr,pot_clust_small,nrclusteratoms_rmax_rnd); fflush(outfile);
                        }
                        if(nn>=nrclusteratoms) break;
                }

                fprintf(outfile, "\n %ld atoms within %f (%f)� of 0, 0, 0\n", imadmaxr, pot,pot_clust_small);
                nn=0;
                for (ipt=0; ipt< iptmax; ipt++) {
                        r2 = xpt[ipt]*xpt[ipt] + ypt[ipt]*ypt[ipt] + zpt[ipt]*zpt[ipt];
                        if (r2 > 0.000001) r2 = sqrt(r2);

                        if (r2 < pot_clust_small) {
                                xmadr[nn]=xpt[ipt];
                                ymadr[nn]=ypt[ipt];
                                zmadr[nn]=zpt[ipt];
                                qmadr[nn]=qpt[ipt];

                                nn++;
                                /*      imadmaxr=nn;*/
                        }
                }
        }
        else{ /* Read .qc file and put in xmadr[] ...*/
                rmspotfiler=fopen(filermspot, "r");

                for (imad=0; imad<MAXMAD; imad++) {
                        input_flag = fscanf(rmspotfiler,"%s %lf %lf %lf %lf",labatm[imad],
                                            &xmadr[imad], &ymadr[imad], &zmadr[imad], &qmadr[imad]);
                        if (input_flag < 4) break;
                        fprintf(outfile,"rmspot:  %s%ld %f %f %f\n",
                                labatm[imad],imad+1,xmadr[imad], ymadr[imad], zmadr[imad]);

                }
                imadmaxr=imad;
                fprintf(outfile,"Efter inlasning fran .rmspot imadmaxr=%ld",imadmaxr); fflush(outfile);
                fclose(rmspotfiler);
        }
        for (imad=0; imad<MAXMAD; imad++) {
                input_flag = fscanf(rmspotfile,"%lf %lf %lf %lf",
                                    &xmad[imad], &ymad[imad], &zmad[imad], &qmad[imad]);

                if (input_flag < 4) break;
        } /* end imad loop */
        imadmax = imad;

        if(imadmax > MAXMAD) fprintf(outfile,"imadmax > MAXMAD  (%ld > %d\n",imadmax,MAXMAD);
        for(pot=1.0; pot<=25.0; pot+=0.01) {
                nn=0;
                for (ipt=0; ipt< iptmax; ipt++) {
                        r2 = xpt[ipt]*xpt[ipt] + ypt[ipt]*ypt[ipt] + zpt[ipt]*zpt[ipt];
                        if (r2 > 0.000001) r2 = sqrt(r2);
                        if (r2 < pot) nn++;
                }
                if(nn<=nrclusteratoms_rmax_rnd) pot_clust_small=pot;
        }

        pot_clust_small = 0.0;
        /*
           satt istallet lika med max av cluxt xyz.
         */

        for(nn=0; nn<imadmaxr; nn++) {
                x=xmadr[nn];
                y=ymadr[nn];
                z=zmadr[nn];
                if(sqrt(x*x+y*y+z*z)>pot_clust_small) pot_clust_small=sqrt(x*x+y*y+z*z);
        }
        pot_clust_small += 1.0;

        /*
         *  Calculate the Madelung distance.
         */
        dmin = 1.0e10;
        for (in1=0; in1<inmax; in1++) {
                for (in2=in1+1; in2<inmax; in2++) {
                        d = (ucx[in1]-ucx[in2]) *  (ucx[in1]-ucx[in2]) +
                            (ucy[in1]-ucy[in2]) *  (ucy[in1]-ucy[in2]) +
                            (ucz[in1]-ucz[in2]) *  (ucz[in1]-ucz[in2]);
                        d = sqrt(d);
                        if (d < dmin) { dmin = d; cnr1=in1; cnr2=in2; }
                }
        }
        scalemad=dmin;
        fflush(outfile);
        /*
         *  End Madelung scalemad calc.
         */

        fprintf(outfile, " Madelung scale distance %f\n", scalemad);
        fflush(outfile);


        fprintf(outfile,"\n\nGenerate random points among central  imadmaxr=%ld+ mad atoms  pot_clust_small=%f\n",
                imadmaxr,pot_clust_small);

        /*
         * Check points with charge = 0.0 !!!
         */

        for (n=0; n<nrandpts_ck; n++) {
                r2 = 1.0e10;
                while (r2 > pot_clust_small*pot_clust_small) {
                        randome(&iseed1, &iseed2, &ran);
                        xval = 2.* (ran-0.5) *pot;
                        randome(&iseed1, &iseed2, &ran);
                        yval = 2.* (ran-0.5) *pot;
                        randome(&iseed1, &iseed2, &ran);
                        zval = 2.* (ran-0.5) *pot;
                        r2 = xval*xval + yval*yval +zval*zval;

                        r=100.0;
                        for(nn=0; nn<imadmaxr; nn++) {
                                x=xval-xmadr[nn];
                                y=yval-ymadr[nn];
                                z=zval-zmadr[nn];
                                if(sqrt(x*x+y*y+z*z)<r) r=sqrt(x*x+y*y+z*z);
                        }
                        if(r<0.1 || r>atom_r_rnd) r2 = 1.0e10;

                } /* end while loop */
                xval_tmp[n]=xval;
                yval_tmp[n]=yval;
                zval_tmp[n]=zval;
                xmad[imadmax] = xval;
                ymad[imadmax] = yval;
                zmad[imadmax] = zval;
                qmad[imadmax] = 0.0;

                imadmax++;
        } /* end n loop */

        fclose(rmspotfile);
        fprintf(outfile," compute xyz coordinates of central %ld+ atoms\n",nrclusteratoms);

        for(pot=1.0; pot<=100.0; pot+=0.01) {
                nn=0;
                for (ipt=0; ipt< iptmax; ipt++) {
                        r2 = xpt[ipt]*xpt[ipt] + ypt[ipt]*ypt[ipt] + zpt[ipt]*zpt[ipt];
                        if (r2 > 0.000001) r2 = sqrt(r2);
                        if(r2<pot) nn++;
                } /* end ipt loop */
                if(nn>=nrclusteratoms && pot>=pot_clust_small) break;
        } /* end pot loop */

        fprintf(outfile, "\n %ld atoms within %f � of 0, 0, 0\n", nn, pot);
        fprintf(outfile,"Note: 2 criteria must be met 1) #zone2 atoms and 2) radius must be > then cluster\n");
        if(nn<=imadmaxr) {
                fprintf(outfile,"nr constraint atoms defining sphere containing zones 1 and 2\n");
                fprintf(outfile,"(%d) must > nr quantum atoms (%d)\n",nn,imadmaxr);
                exit(0);
        }
        nn=0;
        for (ipt=0; ipt< iptmax; ipt++) {
                r2 = xpt[ipt]*xpt[ipt] + ypt[ipt]*ypt[ipt] + zpt[ipt]*zpt[ipt];
                if (r2 > 0.000001) r2 = sqrt(r2);
                if (r2 < pot) {
                        itemp= typept[ipt];
                        if(DEBUG > 0) fprintf(outfile, " ipt %ld x %.5f y %.5f z %.5f q %.4f r %f %s\n",
                                              ipt, xpt[ipt], ypt[ipt], zpt[ipt], qpt[ipt], r2, typetab[itemp]);
                }

                if(r2< pot) {
                        xyz[nn][0]=xpt[ipt];
                        xyz[nn][1]=ypt[ipt];
                        xyz[nn][2]=zpt[ipt];
                        xyz[nn][3]=qpt[ipt];
                        typxyz[nn]=typept[ipt];
                        nn++;
                }
        } /* end ipt loop */

        strncpy(fileclust,fileroot,50);
        strncat(fileclust,".clust",20);
        clustfile = fopen(fileclust,"w");

        for(it=0; it<itmax; it++) {
                ntype[it]=0;
                for (n=0; n<nn; n++) {
                        if (typxyz[n] == it)
                                ntype[it] += 1;
                } /* end n loop over central 50+ atoms */
        } /* end it loop over atom types */

        xmin=1.E20;
        ymin=1.E20;
        zmin=1.E20;
        xmax=-1.E20;
        ymax=-1.E20;
        zmax=-1.E20;

        for(n=0; n<nn; n++) {
                itemp = typxyz[n];
                fprintf(clustfile, "%7.5f  %7.5f  %7.5f  %5.1f\n",xyz[n][0],
                        xyz[n][1],xyz[n][2], xyz[n][3]);
                if (xyz[n][0] < xmin) xmin = xyz[n][0];
                if (xyz[n][1] < ymin) ymin = xyz[n][1];
                if (xyz[n][2] < zmin) zmin = xyz[n][2];
                if (xyz[n][0] > xmax) xmax = xyz[n][0];
                if (xyz[n][1] > ymax) ymax = xyz[n][1];
                if (xyz[n][2] > zmax) zmax = xyz[n][2];
        } /* end n loop over central cluster */


        /*
         * generate random points among central 50+ atoms and write to the
         * end of the .clust file for use as a .rmspot file  for subsequent runs
         */

        for (n=0; n<nrandpts; n++) {
                r2 = 2. * pot*pot;
                while (r2 > pot*pot) {
                        randome(&iseed1, &iseed2, &ran);
                        xval = 2.* (ran-0.5) *pot;
                        randome(&iseed1, &iseed2, &ran);
                        yval = 2.* (ran-0.5) *pot;
                        randome(&iseed1, &iseed2, &ran);
                        zval = 2.* (ran-0.5) *pot;
                        r2 = xval*xval + yval*yval +zval*zval;
                } /* end while loop */
                fprintf(clustfile, "%7.5f  %7.5f  %7.5f    0.0   0.000000\n",xval, yval, zval);
        } /* end n loop */
        fclose(clustfile);

        fprintf(outfile, " central cluster xmin %f  xmax %f\n", xmin, xmax);
        fprintf(outfile, " central cluster ymin %f  ymax %f\n", ymin, ymax);
        fprintf(outfile, " central cluster zmin %f  zmax %f\n", zmin, zmax);
        fflush(outfile);
        for (imad=0; imad<imadmax; imad++) {
                for (ipt=0; ipt<iptmax; ipt++) { /* find atom type for this rmspot atom*/
                        ifmatch =1;
                        if (fabs(xmad[imad]-xpt[ipt]) > 0.001) ifmatch=0;
                        if (fabs(ymad[imad]-ypt[ipt]) > 0.001) ifmatch=0;
                        if (fabs(zmad[imad]-zpt[ipt]) > 0.001) ifmatch=0;
                        if (ifmatch == 0) continue;
                        typmad[imad] = typept[ipt];
                        ifmatch=1;
                        break;
                }
                if (ifmatch ==0) typmad[imad] = itmax;
        }

        /*
         *  Compute Ewald sums.
         */

        fprintf(outfile," Compute Ewald sums\n");

        alpha = 0.2;                                   /* Ewald convergence parameter!!!*/
        fprintf(outfile, " alpha %f  pi %f\n", alpha, pi);
        irmax =8; /* 8 */                                  /* Max shells in real space sum  */
        for (imad=0; imad<imadmax; imad++) { /* loop to find Ewald sums for each atom of interest */
                ewaldr[imad] =0.; /* real space term */
                ewaldm[imad] =0.; /* reciprical lattice term */
                for (ishell = 0; ishell <10; ishell++) {
                        eshellr[ishell] = 0.;
                        eshellm[ishell] = 0.;
                }
                for (ir1 = -irmax; ir1 <= irmax; ir1++) { /* loop over first translation vector */
                        for (ir2 = -irmax; ir2 <= irmax; ir2++) { /* loop over second translation vector */
                                for (ir3 = -irmax; ir3 <= irmax; ir3++) { /* loop over third translation vector */
                                        for (in=0; in<inmax; in++) { /* loop over atoms in unit cell */
                                                uctmp1 = uc1[in] + ir1;
                                                uctmp2 = uc2[in] + ir2;
                                                uctmp3 = uc3[in] + ir3;
                                                xptr = uctmp1 * a1x + uctmp2 * a2x + uctmp3 * a3x;
                                                yptr = uctmp1 * a1y + uctmp2 * a2y + uctmp3 * a3y;
                                                zptr = uctmp1 * a1z + uctmp2 * a2z + uctmp3 * a3z;
                                                qptr = quc[in];
                                                r2 = (xmad[imad]-xptr) * (xmad[imad]-xptr) +
                                                     (ymad[imad]-yptr) * (ymad[imad]-yptr) +
                                                     (zmad[imad]-zptr) * (zmad[imad]-zptr);
                                                if (r2 > 0.000001) {
                                                        r2 = sqrt(r2);
                                                        temp = evconv * qptr * erfc(alpha*r2)/r2;
                                                        ewaldr[imad] += temp;
                                                        ishell = abs(ir1);
                                                        if (abs(ir2) > ishell) ishell = abs(ir2);
                                                        if (abs(ir3) > ishell) ishell = abs(ir3);
                                                        eshellr[ishell] += temp;
                                                }
                                        } /* end in loop */
                                } /* end it3 loop */
                        } /* end it2 loop */
                } /* end it1 loop */

                for (in=0; in<inmax; in++) { /* loop over atoms in unit cell */
                        ux = uc1[in]*a1x +uc2[in]*a2x + uc3[in]*a3x;
                        uy = uc1[in]*a1y +uc2[in]*a2y + uc3[in]*a3y;
                        uz = uc1[in]*a1z +uc2[in]*a2z + uc3[in]*a3z;
                        uq = quc[in];
                        /* loop over three reciprical lattice vectors */
                        mmax = 5;                    /* Max shells in reciprocal space  */
                        for (m1=-mmax; m1<=mmax; m1++) {
                                for (m2=-mmax; m2<=mmax; m2++) {
                                        for (m3=-mmax; m3<=mmax; m3++) {
                                                if ( (m1==0) && (m2==0) && (m3 == 0)) continue;
                                                rx = m1*b1x + m2*b2x + m3*b3x;
                                                ry = m1*b1y + m2*b2y + m3*b3y;
                                                rz = m1*b1z + m2*b2z + m3*b3z;
                                                r2 = rx*rx + ry*ry + rz*rz;
                                                dot = rx*(xmad[imad]-ux) + ry*(ymad[imad]-uy) +
                                                      rz*(zmad[imad]-uz);
                                                temp1 = evconv * uq/pi/volume;
                                                temp2 = exp(-pi*pi*r2/alpha/alpha);
                                                temp3 = cos(tpi*dot)/r2;
                                                temp =  temp1 * temp2 * temp3;
                                                ewaldm[imad] +=  temp;
                                                ishell = abs(m1);
                                                if (abs(m2) > ishell) ishell = abs(m2);
                                                if (abs(m3) > ishell) ishell = abs(m3);
                                                eshellm[ishell] +=  temp;
                                        } /* end loop over m3 */
                                } /* end loop over m2 */
                        } /* end loop over m1 */
                } /* end loop over unit cell atoms to compute Um*/
                ewaldc[imad] = -evconv * 2. * alpha*qmad[imad]/sqrt(pi);
                ewald[imad] = ewaldr[imad] + ewaldm[imad] + ewaldc[imad];
                itemp = typmad[imad];
                fprintf(outfile, "\n imad %ld xmad %f ymad %f zmad %f  %s\n",
                        imad, xmad[imad], ymad[imad], zmad[imad], typetab[itemp]);
                if(DEBUG > 0)
                        for (ishell =0; ishell<10; ishell+=1) {
                                fprintf(outfile, "ishell %ld eshellr %e eshellm %e\n",
                                        ishell, eshellr[ishell], eshellm[ishell]);
                        }
                itemp=typmad[imad];
                if (qmad[imad] != 0.0) {
                        fprintf(outfile, " Ewald pot. (V): r %.8f m %.8f c %.8f tot %.8f\n",
                                ewaldr[imad], ewaldm[imad], ewaldc[imad], ewald[imad]);
                }
                else{
                        fprintf(outfile, "  ewald pot (V): r %.8f m %.8f c %.8f tot %.8f\n",
                                ewaldr[imad], ewaldm[imad], ewaldc[imad], ewald[imad]);

                }
        } /* end imad loop */

        ewaldmad = 0.;


        /*
         *  Moments...
         */

        ifill =0.;
        dx=0.;   dy=0.;   dz=0.;
        dxx=0.;  dyy=0.;  dzz=0.;
        xmin = 1.E20;
        ymin = 1.E20;
        zmin = 1.E20;
        xmax = -1.E20;
        ymax = -1.E20;
        zmax = -1.E20;
        qsum =0.;
        for (ipt=0; ipt<iptmax; ipt++) {
                if (zone[ipt] < 3) {
                        ifill++;
                        dx   += xpt[ipt]*qpt[ipt];
                        dxx  += xpt[ipt]*xpt[ipt]*qpt[ipt];
                        dy   += ypt[ipt]*qpt[ipt];
                        dyy  += ypt[ipt]*ypt[ipt]*qpt[ipt];
                        dz   += zpt[ipt]*qpt[ipt];
                        dzz  += zpt[ipt]*zpt[ipt]*qpt[ipt];
                        qsum += qpt[ipt];
                        if (xpt[ipt] < xmin) xmin = xpt[ipt];
                        if (ypt[ipt] < ymin) ymin = ypt[ipt];
                        if (zpt[ipt] < zmin) zmin = zpt[ipt];
                        if (xpt[ipt] > xmax) xmax = xpt[ipt];
                        if (ypt[ipt] > ymax) ymax = ypt[ipt];
                        if (zpt[ipt] > zmax) zmax = zpt[ipt];
                }
        } /* end ipt loop */
        radmom = sqrt(dxx*dxx + dyy*dyy + dzz*dzz);
        fprintf(outfile,"dipole: dx % .8f dy %.8f dz %.8f qsum %f\n",
                dx, dy, dz, qsum);
        fprintf(outfile, "    dxx %.8f  dyy %.8f  dzz %.8f   radmad %f\n",
                dxx, dyy, dzz, radmom);
        if (qsum != 0.0) fprintf(outfile, " error- qsum %f\n", qsum);

        if(DEBUG > 0) {
                fprintf(outfile, "lattice xmin %f xmax %.5f\n",xmin, xmax);
                fprintf(outfile, "lattice ymin %f ymax %.5f\n",ymin, ymax);
                fprintf(outfile, "lattice zmin %f zmax %.5f\n",zmin, zmax);
        }

        /*
         * compute initial cluster potentials
         */

        fprintf(outfile, "\ncluster atoms and initial potentials:\n");


        for (imad = 0; imad < imadmax; imad++) { /* loop over all cluster atoms */
                dspot[imad] = 0.;
                itemp=typmad[imad];
                for (ipt =0; ipt<iptmax; ipt++) { /* loop over all occupied sites */
                        if (zone[ipt] >= 3) continue;
                        r2 = (xpt[ipt]-xmad[imad]) * (xpt[ipt]-xmad[imad]) +
                             (ypt[ipt]-ymad[imad]) * (ypt[ipt]-ymad[imad]) +
                             (zpt[ipt]-zmad[imad]) * (zpt[ipt]-zmad[imad]);
                        if (r2 < 0.000001) zone[ipt] = 0;
                        if (r2 > 0.000001) dspot[imad] += evconv*qpt[ipt]/sqrt(r2);
                }
                if (qmad[imad] != 0.0) {
                        fprintf(outfile, "imad %ld x %.5f y %.5f z %.5f q %.4f %s\n",
                                imad, xmad[imad], ymad[imad], zmad[imad], qmad[imad],  typetab[itemp]);
                        fprintf(outfile, "     dspot %f eV  ewald %f eV  � %f V\n",
                                dspot[imad], ewald[imad], dspot[imad]-ewald[imad]);
                }
                else{
                        fprintf(outfile, "imad %ld x %.5f y %.5f z %.5f q %.4f %s\n",
                                imad, xmad[imad], ymad[imad], zmad[imad], qmad[imad],typetab[itemp]);
                        fprintf(outfile, "     dspot %f V  ewald %f V  � %f V\n",
                                dspot[imad], ewald[imad], dspot[imad]-ewald[imad]);
                }

                dspotFA[imad]=dspot[imad];
        } /* end imad loop over cluster atoms */

        /*
         * Moments...
         */

        dx = 0.;  dy = 0.;  dz = 0.; qq =0.;
        dxx = 0.;  dyy = 0.;  dzz = 0.;
        radmom = 0.;
        for (in1 = 0; in1<inmax; in1++) {
                dx += ucx[in1] * quc[in1];
                dxx += ucx[in1] * ucx[in1] * quc[in1];
                dy += ucy[in1] * quc[in1];
                dyy += ucy[in1] * ucy[in1] * quc[in1];
                dz += ucz[in1] * quc[in1];
                dzz += ucz[in1] * ucz[in1] * quc[in1];
                qq += quc[in1];
        }
        d2 = dx*dx + dy*dy + dz*dz;
        radmom = sqrt(dxx*dxx + dyy*dyy + dzz*dzz);
        fflush(outfile);


        rmspot = 0.;
        avepot =0.;
        nomad =0;
        for (imad=0; imad<imadmax; imad++) {
                if (qmad[imad] == 0.0) continue;
                nomad++;
                avepot += ewald[imad] - dspot[imad];
                rmspot += (ewald[imad] - dspot[imad] ) * (ewald[imad] - dspot[imad] );
        }
        rmspot = sqrt (rmspot/nomad);
        avepot = avepot/nomad;
        fprintf(outfile, "\nResults before solving the system of linear eq.\n");
        avepotuc2 = avepot;
        rmspotuc2 = rmspot;
        fprintf(outfile, " (Ewald - direct sum) potential average %f rms %f\n",
                avepot, rmspot);


        rmspot2=rmspot;

        for (in1 =0; in1<inmax; in1++) {
                ucx[in1] = uc1[in1] * a1x + uc2[in1] * a2x + uc3[in1] * a3x;
                ucy[in1] = uc1[in1] * a1y + uc2[in1] * a2y + uc3[in1] * a3y;
                ucz[in1] = uc1[in1] * a1z + uc2[in1] * a2z + uc3[in1] * a3z;
        }

        dx = 0.;  dy = 0.;  dz = 0.; qq =0.;
        dxx = 0.;  dyy = 0.;  dzz = 0.;
        radmom = 0.;
        for (in1 = 0; in1<inmax; in1++) {
                dx  += ucx[in1] * quc[in1];
                dxx +=ucx[in1] * ucx[in1] * quc[in1];
                dy  += ucy[in1] * quc[in1];
                dyy +=ucy[in1] * ucy[in1] * quc[in1];
                dz  += ucz[in1] * quc[in1];
                dzz +=ucz[in1] * ucz[in1] * quc[in1];
                qq  += quc[in1];
        }
        d2 = dx*dx + dy*dy + dz*dz;
        radmom = sqrt(dxx*dxx + dyy*dyy + dzz*dzz);
        fprintf(outfile, " dipole dx %f  dy %f  dz %f  d2 %f  radmom %f\n",dx, dy, dz, d2, radmom);
        fflush(outfile);

        energy =0.;
        for (in1 =1; in1<inmax; in1++) {
                for (in2=0; in2<in1; in2++) {
                        r2 = (ucx[in1]-ucx[in2]) * (ucx[in1]-ucx[in2]) +
                             (ucy[in1]-ucy[in2]) * (ucy[in1]-ucy[in2]) +
                             (ucz[in1]-ucz[in2]) * (ucz[in1]-ucz[in2]);
                        if (r2 > 0.001) energy += quc[in1]*quc[in2]/sqrt(r2);
                        else{
                                fprintf(outfile, "Detect zero r2 in final energy \n");
                                fprintf(outfile, "  r2 %f   in1 %ld  ucx %f  uxy   %f  ucz  %f\n",
                                        r2, in1, ucx[in1],  ucy[in1], ucz[in1]);
                                fprintf(outfile, "  r2 %f   in %ld  ucx %f  uxy   %f  ucz  %f\n",
                                        r2, in, ucx[in],  ucy[in], ucz[in]);
                        }
                }
        }
        energy *= evconv;

        fprintf(outfile, " unit cell: energy %f\n",energy);

        /*
         *
         *  Printouts for Jaguar and NWChem.
         *
         */

        q=0.0;
        for (ipt=0; ipt<iptmax; ipt++) {q += qpt[ipt]; }
        fprintf(ptsfilenw,"Total charge of this file (madelung atoms + point charges) is %5.4e\n",q);
        fprintf(ptsfilenw,"The first %ld atoms in this file are the central madelung atoms\n",imadmax-1);
        fprintf(ptsfilenw,"The madelung atoms are sorted with increasing r (distance to origin)\n");
        fprintf(ptsfilenw,"Replace the appropriate madelung bq's with quantum atoms (the central cluster)\n\n");

        nrchk=0;
        /*MIGUEL WAS HERE*/
        fprintf(outfile,"START SAMPLE SPHERICAL POINTS\n");
        for (imad=0; imad<imadmax; imad++) {
                int ityptmp1 = typmad[imad];
                fprintf(outfile,"%s  %f  %f  %f  \n", typetab[ityptmp1], xmad[imad], ymad[imad], zmad[imad]);
        }
        fprintf(outfile,"END SAMPLE SPHERICAL POINTS\n");

        /*MIGUEL WAS HERE*/

        for (imad=0; imad< imadmax; imad++) {
                if(fabs(qmad[imad])!= 0.0) {
                        r=sqrt(xmad[imad]*xmad[imad]+ymad[imad]*ymad[imad]+zmad[imad]*zmad[imad]);
                        if(nrclusteratoms_rmax_rnd!=0) {
                                fprintf(ptsfile, "  %s%ld\t%25.20f \t%25.20f \t%25.20f \n",
                                        typetab[typmad[imad]],imad+1,xmad[imad],  ymad[imad], zmad[imad]);
                                fprintf(ptsfilenw, "  %s%ld\t%25.20f\t%25.20f\t%25.20f\n",typetab[typmad[imad]],
                                        imad,xmad[imad],  ymad[imad], zmad[imad]);

                        }
                        if(nrclusteratoms_rmax_rnd==0) {
                                nn=0;
                                for(n=0; n<imadmaxr; n++) { /*imadmaxr ???*/
                                        if((fabs(xmad[imad]-xmadr[n])<0.01)&&(fabs(ymad[imad]-ymadr[n])<0.01)&&
                                           (fabs(zmad[imad]-zmadr[n])<0.01)&&(fabs(qmad[imad]-qmadr[n])<0.01)) {nn=1;
                                                                                                                /*MIGUEL WAS HERE*/
                                                                                                                int ityptmp = typmad[imad];
                                                                                                                fprintf(outfile,"%s  %f  %f  %f  \n", typetab[ityptmp], xmad[imad], ymad[imad], zmad[imad]);
                                                                                                                /*MIGUEL WAS HERE*/
                                                                                                                break; }
                                }

                                if(nn==1 && imad<imadmax) {
                                        listan[n] = imad; /* n om xmadr[]... is used*/
                                        listan2[n] = n;
                                        nrchk+=1;
                                }
                        }
                }
        }

        if(nrclusteratoms_rmax_rnd==0) {
                for(nn=0; nn<nrchk; nn++) {
                        n=listan[nn];
                        fprintf(ptsfile, "  %s%ld\t%25.20f\t%25.20f\t%25.20f \n",
                                labatm[listan2[nn]],listan2[nn]+1,xmad[n],  ymad[n], zmad[n]);
                        fprintf(ptsfilenw, "  %s\t%25.20f\t%25.20f\t%25.20f\n",labatm[listan2[nn]],
                                xmad[n],  ymad[n], zmad[n]);
                }
        }

        if(nrclusteratoms_rmax_rnd==0) {
                if(nrchk!=imadmaxr) {
                        fprintf(outfile,"\nAll rmspot are not on grid points nrchk=%ld imadmaxr=%ld)!!\n",nrchk,imadmaxr);
                        fprintf(outfile,"A probable cause is too few constraint points or .uc and .qc are not compatible.\n\n\n");
                        exit(0);
                }
        }
        if(nrclusteratoms_rmax_rnd==0)
                for (imad=0; imad< imadmax; imad++) {
                        if(fabs(qmad[imad])>0.01) {
                                nn=0;
                                for(n=0; n<imadmaxr; n++) {
                                        if((fabs(xmad[imad]-xmadr[n])<0.01)&&(fabs(ymad[imad]-ymadr[n])<0.01)&&
                                           (fabs(zmad[imad]-zmadr[n])<0.01)&&(fabs(qmad[imad]-qmadr[n])<0.01)) nn=1;
                                }
                                if(nn!=1) {
                                        fprintf(ptsfile, "  %25.20f \t%25.20f \t%25.20f \t%25.20f\n",
                                                qmad[imad],xmad[imad],  ymad[imad], zmad[imad]);
                                        fprintf(ptsfilenw, "  bq%ld %25.20f %25.20f %25.20f charge %25.20f\n",
                                                imad+1,xmad[imad],ymad[imad],zmad[imad],qmad[imad]);
                                        fprintf(ptsfilefro,"%25.20f %25.20f %25.20f %25.20f\n",
                                                xmad[imad],ymad[imad],zmad[imad],qmad[imad]);
                                }
                        }
                }


        /*
         *  Fix the dip.mom.
         */

        if((varycharges==1)||(varycharges==2)||(varycharges==3)) {
                rmspot = 0.;
                avepot = 0.;
                nomad  = 0;


                for (imad = 0; imad < imadmax; imad++) { /* loop over all cluster atoms to compute direct sum potentials*/
                        dspot[imad] = 0.;
                        dipx_org = 0.0;
                        dipy_org = 0.0;
                        dipz_org = 0.0;
                        for (ipt =0; ipt<iptmax; ipt++) {
                                if (zone[ipt] >= 3) continue;

                                dipx_org += xpt[ipt]*qpt[ipt];
                                dipy_org += ypt[ipt]*qpt[ipt];
                                dipz_org += zpt[ipt]*qpt[ipt];

                                r2 = (xpt[ipt]-xmad[imad]) * (xpt[ipt]-xmad[imad]) +
                                     (ypt[ipt]-ymad[imad]) * (ypt[ipt]-ymad[imad]) +
                                     (zpt[ipt]-zmad[imad]) * (zpt[ipt]-zmad[imad]);
                                if (r2 > 0.000001) dspot[imad] += evconv*qpt[ipt]/sqrt(r2);
                        }
                } /* end imad loop */
                dipr=sqrt(dipx_org*dipx_org+dipy_org*dipy_org+dipz_org*dipz_org);

                for (imad = 0; imad < imadmax; imad++) { /* loop over all cluster atoms to compute direct sum potentials*/
                        if (qmad[imad] != 0.0) {


                                nomad++;
                                avepot += ewald[imad] - dspot[imad];
                                rmspot += (ewald[imad] - dspot[imad] ) * (ewald[imad] - dspot[imad] );
                        }
                } /* end loop over cluster atoms */
                rmspot2 = sqrt (rmspot/nomad);
                avepot  = avepot/nomad;


                rmspot = rmspot2;
                xxg[4] = rmspot;


                mxzon2 = 0;
                for (ipt =0; ipt<iptmax; ipt++) if (zone[ipt] >0) { /*==2*/
                                zone2list[mxzon2]=ipt;
                                mxzon2++;
                        }

                qmax      = -1.0e2;
                qmin      = 1.0e2;
                dipfact2  = 1.0;
                rmsfact2  = 1.0;
                ddq       = 0.05; /* 0.01 orig*/

                itera     = 0.0;
                itera_rms = 0.0;
                itera_dip = 0.0;
                /* The number that is important is the number of rms trials. Dip trials goes very fast. */

                while((itera_rms<=itermax)&&(itera_dip<=2.0e7)) {
                        ignore = 0;
                        itera+=1.0;
                        if(rmsfact2==1.0) itera_rms += 1.0;
                        if(dipfact2==1.0) itera_dip += 1.0;
                        nr1=1; nr2=1; /* Maste vara lika fran borjan */
                        while((nr1==nr2)||(nr1==nr3)||(nr1==nr4)||(nr1==nr5)||
                              (nr2==nr3)||(nr2==nr4)||(nr2==nr5)||
                              (nr3==nr4)||(nr3==nr5)||(nr4==nr5)) {
                                randome(&iseed1, &iseed2, &ran);
                                nr1 = (int)(ran*mxzon2);
                                randome(&iseed1, &iseed2, &ran);
                                nr2 = (int)(ran*mxzon2);
                                randome(&iseed1, &iseed2, &ran);
                                nr3 = (int)(ran*mxzon2);
                                randome(&iseed1, &iseed2, &ran);
                                nr4 = (int)(ran*mxzon2);
                                randome(&iseed1, &iseed2, &ran);
                                nr5 = (int)(ran*mxzon2);
                        }
                        z2list[0]=zone2list[nr1];
                        z2list[1]=zone2list[nr2];
                        z2list[2]=zone2list[nr3];
                        z2list[3]=zone2list[nr4];
                        z2list[4]=zone2list[nr5];

                        /*
                         *  Lapack stuff if varycharges==1.
                         */

                        if(varycharges==1) {

                                fprintf(outfile,"\n\nlapack...\n\n"); fflush(outfile);
                                fprintf(outfile,"Results after solving the system of linear eq.\n"); fflush(outfile);

                                itera_rms=itermax+1;

                                nr_e_eq = 4;
                                imadmax += nr_e_eq;
                                nr_equ = imadmax-nrandpts_ck;
                                nr_par = (int)(mxzon2);

                                M=(integer)(imadmax-nrandpts_ck);
                                N=(integer)(nr_par);
                                NRHS=(integer)(1);
                                LDA=(integer)(imadmax-nrandpts_ck);
                                LDB=(integer)(nr_par);
                                LWORK=(integer)(N+N+N+N);
                                A = malloc(LDA*LDB*sizeof(double));
                                B = malloc(LDB*sizeof(double));
                                WORK = malloc(LWORK*sizeof(double));

                                S = malloc(M*sizeof(double));
                                RCOND = (double)(-1.0);

                                JPVT = malloc(LDB*sizeof(double));

                                if(DEBUG > 0) {
                                        fprintf(outfile,"M=%ld\n",M);
                                        fprintf(outfile,"N=%ld\n",N);
                                        fprintf(outfile,"NRHS=%ld\n",NRHS);
                                        fprintf(outfile,"LDA=%ld\n",LDA);
                                        fprintf(outfile,"LDB=%ld\n",LDB);
                                        fprintf(outfile,"LWORK=%ld\n",LWORK);
                                        fprintf(outfile,"\n");
                                        fprintf(outfile,"imadmax=%ld\n",imadmax);
                                        fprintf(outfile,"nr_equ=%ld\n",nr_equ);
                                        fprintf(outfile,"nr_e_eq=%ld\n",nr_e_eq);
                                        fprintf(outfile,"nr_par=%ld\n",nr_par);
                                        fprintf(outfile,"\n");
                                        fprintf(outfile,"qsum=%f\n",qsum);
                                        fprintf(outfile,"dipx_org=%f\n",dipx_org);
                                        fprintf(outfile,"dipy_org=%f\n",dipy_org);
                                        fprintf(outfile,"dipz_org=%f\n",dipz_org);
                                        for(i=0; i<nr_equ-nr_e_eq; i++) {
                                                if(DEBUG==1) fprintf(outfile,"ewald[i]=%f\tdspot[i]=%f\tewald[i]-dspot[i]=%f\n",ewald[i],dspot[i],ewald[i]-dspot[i]);
                                        }
                                }
                                fflush(outfile);

                                for(n=0; n<nr_par; n++) {z2list[n]=zone2list[n]; JPVT[n]=0.0; }

                                jj=0;
                                for(j=0; j<nr_par; j++) {

                                        n=nr_equ-nr_e_eq;
                                        for(i=0; i<n; i++) {

                                                x = xmad[i] - xpt[z2list[j]];
                                                y = ymad[i] - ypt[z2list[j]];
                                                z = zmad[i] - zpt[z2list[j]];
                                                r=sqrt(x*x + y*y + z*z);

                                                A[jj] = (1.0/r)*evconv;
                                                jj++;
                                                B[i] = ewald[i]-dspot[i]; /* I solve for dq, not q!!!*/

                                                if(j==0 && DEBUG >0) fprintf(outfile,"%ld\tx=%f,%f,%f\tB[%ld]=%f\t%ld=%ld-%ld\n",
                                                                             i,x,y,z,i,B[i],n,nr_equ,nr_e_eq);
                                        }
                                        fflush(outfile);
                                        A[jj] = 1.0;
                                        jj++;
                                        B[i+0] = -qsum; /* q if only parameters should be 0 */

                                        x = xpt[z2list[j]];
                                        y = ypt[z2list[j]];
                                        z = zpt[z2list[j]];
                                        r2=x*x + y*y + z*z;

                                        A[jj] = x;
                                        jj++;
                                        B[i+1] = -dipx_org;
                                        A[jj] = y;
                                        jj++;
                                        B[i+2] = -dipy_org;
                                        A[jj] = z;
                                        jj++;
                                        B[i+3] = -dipz_org;
                                }

                                if(DEBUG > 0) fprintf(outfile,"\njapp imadmax-nrandpts_ck=%ld\n\n",imadmax-nrandpts_ck);

                                dgelsy_(&M,&N,&NRHS,A,&LDA,B,&LDB,JPVT,&RCOND,&RANK,WORK,&LWORK,&INFO);
                                /* printf("RCOND=%lf\n",RCOND);
                                 */
                                if(INFO!=0) {fprintf(outfile,"\nInversion problem in lapack rountine  %ld\n\n\n",INFO); ignore=1; }
                                fflush(outfile);
                                /*-- hmmmmmmmm kolla upp mer.
                                   for(j=0;j<nr_par;j++){
                                   qtmp[j]=qpt[z2list[j]];
                                   }
                                 */

                                for(j=0; j<nr_par; j++) {
                                        qpt[z2list[j]] += B[j];
                                        if(DEBUG > 0) fprintf(outfile,"j=%ld\tqpt[z2list[j]]=%f\n",j,qpt[z2list[j]]);
                                }

                                imadmax -= nr_e_eq;

                                if(DEBUG > 0)
                                        for(j=0; j<nr_par; j++) {
                                                fprintf(outfile,"--- %5.3f + %5.3f = %5.3f\n",B[j],qpt[z2list[j]]-B[j],qpt[z2list[j]]);
                                        }

                                free(A);
                                free(B);
                                free(WORK);
                                free(S);
                                free(JPVT);

                                rmspot = 0.;
                                avepot = 0.;
                                nomad  = 0;


                                for (imad = 0; imad < imadmax-nrandpts_ck; imad++) {
                                        if (qmad[imad] != 0.0) {

                                                nomad++;
                                                avepot += ewald[imad] - dspottmp[imad];
                                                rmspot += (ewald[imad] - dspottmp[imad] ) * (ewald[imad] - dspottmp[imad] );
                                        }
                                }
                                rmspot = sqrt (rmspot/nomad);
                                avepot = avepot/nomad;



                                /* printf("%lf    ",rmspot);*/

                                for (imad = 0; imad < imadmax-nrandpts_ck; imad++) dspot[imad]=dspottmp[imad];
                                for(n=0; n<nr_par; n++) {
                                        if(qpt[z2list[n]]>qmax) qmax=qpt[z2list[n]];
                                        if(qpt[z2list[n]]<qmin) qmin=qpt[z2list[n]];
                                }

                                qmin=1.0e10;
                                qmax=-1.0e10;
                                for (ipt =0; ipt<iptmax; ipt++) {
                                        if(qpt[ipt]<qmin) qmin=qpt[ipt];
                                        if(qpt[ipt]>qmax) qmax=qpt[ipt];
                                }

                                fprintf(outfile,"q ranges:   %f  to  %f\n",qmin,qmax);

                                rmspot = 0.;
                                avepot = 0.;
                                nomad  = 0;


                                for (imad = 0; imad < imadmax-nrandpts_ck; imad++) { /* loop over all cluster atoms to compute direct sum potentials*/
                                        dspot[imad] = 0.;
                                        for (ipt =0; ipt<iptmax; ipt++) { /* loop over all occupied sites */
                                                if (zone[ipt] >= 3) continue;
                                                r2 = (xpt[ipt]-xmad[imad]) * (xpt[ipt]-xmad[imad]) +
                                                     (ypt[ipt]-ymad[imad]) * (ypt[ipt]-ymad[imad]) +
                                                     (zpt[ipt]-zmad[imad]) * (zpt[ipt]-zmad[imad]);
                                                if (r2 > 0.000001) dspot[imad] += evconv*qpt[ipt]/sqrt(r2);
                                        }
                                } /* end imad loop */

                                for (imad = 0; imad < imadmax-nrandpts_ck; imad++) { /* loop over all cluster atoms to compute direct sum potentials*/
                                        if (qmad[imad] != 0.0) {

                                                nomad++;
                                                avepot += ewald[imad] - dspot[imad];
                                                rmspot += (ewald[imad] - dspot[imad] ) * (ewald[imad] - dspot[imad] );
                                        }
                                } /* end loop over cluster atoms */

                                rmspot = sqrt (rmspot/nomad);
                                avepot = avepot/nomad;

                                fprintf(outfile,"rmspot = %f\n",rmspot);
                        } /* varycharges==1 */


                        /*
                         *  The random approach. Default.
                         */

                        if((varycharges==2)||(varycharges==3)) {

                                ignore=0;

                                qlist[0]=qpt[z2list[0]];
                                qlist[1]=qpt[z2list[1]];
                                q0= qlist[0];
                                randome(&iseed1, &iseed2, &ran);

                                qlist[0]+=2.0*(ran-0.5)*ddq;
                                qlist[1]-=2.0*(ran-0.5)*ddq;
                                dq=-q0+qlist[0];

                                for(j=0; j<2; j++) {
                                        qtmp[j]=qpt[z2list[j]];
                                        qpt[z2list[j]]=qlist[j];
                                }

                                for (imad = 0; imad < imadmax-nrandpts_ck; imad++) {
                                        dspottmp[imad] = dspot[imad];
                                        dipxtmp=dipx_org;
                                        dipytmp=dipy_org;
                                        dipztmp=dipz_org;
                                        for (j =0; j<2; j++) {
                                                ipt=z2list[j];
                                                if (zone[ipt] >= 3) continue;
                                                r2 = (xpt[ipt]-xmad[imad]) * (xpt[ipt]-xmad[imad]) +
                                                     (ypt[ipt]-ymad[imad]) * (ypt[ipt]-ymad[imad]) +
                                                     (zpt[ipt]-zmad[imad]) * (zpt[ipt]-zmad[imad]);
                                                if (r2 > 0.000001) {
                                                        dspottmp[imad] -= evconv*qtmp[j]/sqrt(r2);
                                                        dspottmp[imad] += evconv*qpt[ipt]/sqrt(r2);
                                                }
                                                dipxtmp += xpt[ipt]*(qpt[ipt]-qtmp[j]);
                                                dipytmp += ypt[ipt]*(qpt[ipt]-qtmp[j]);
                                                dipztmp += zpt[ipt]*(qpt[ipt]-qtmp[j]);
                                        }
                                }
                                if(varycharges==2) diprtmp=sqrt(dipxtmp*dipxtmp+dipytmp*dipytmp+dipztmp*dipztmp)/10.0;
                                if(varycharges==3) diprtmp=sqrt(dipxtmp*dipxtmp+dipytmp*dipytmp+dipztmp*dipztmp);
                        } /* varycharges==2 || 3 */

                        if((fmod(itera,1000)==0.0)||(itera==itermax-1)) {
                                rrms = 0.;
                                rnomad  = 0;

                                rmspot=0.0;
                                nomad =0;

                                for (imad = 0; imad < imadmax; imad++) {
                                        dspotr[imad] = 0.;
                                        for (ipt =0; ipt<iptmax; ipt++) {
                                                if (zone[ipt] >= 3) continue;
                                                r2 = (xpt[ipt]-xmad[imad]) * (xpt[ipt]-xmad[imad]) +
                                                     (ypt[ipt]-ymad[imad]) * (ypt[ipt]-ymad[imad]) +
                                                     (zpt[ipt]-zmad[imad]) * (zpt[ipt]-zmad[imad]);
                                                if (r2 > 0.000001) dspotr[imad] += evconv*qpt[ipt]/sqrt(r2);
                                        }
                                } /* end imad loop */

                                for (imad = 0; imad < imadmax; imad++) {
                                        if (qmad[imad] <= 1.0e-16) {
                                                rnomad++;
                                                rrms += (ewald[imad] - dspotr[imad] ) * (ewald[imad] - dspotr[imad] );
                                        }
                                        else{
                                                nomad++;
                                                rmspot += (ewald[imad] - dspotr[imad] ) * (ewald[imad] - dspotr[imad] );
                                        }
                                } /* end loop over cluster atoms */
                                rrms2 = sqrt (rrms/rnomad);
                                avepot = avepot/rnomad;


                                rmspot = sqrt (rmspot/rnomad);

                                if(DEBUG >0) fprintf(outfile,"%10.0f (d %8.0f, r %8.0f)\trmspot_old   new :   %f\t%f  (%f)",
                                                     itera,itera_dip,itera_rms,rmspot2,rmspot,rrms2);
                                fflush(outfile);
                        } /* fmod(,1000) */

                        if(varycharges==2) {
                                if(fabs(rmspot/diprtmp)<=0.5) {rmsfact2=0.0; dipfact2=1.0; }
                                else if(fabs(rmspot/diprtmp)>2.0) {rmsfact2=1.0; dipfact2=0.0; }
                                else{rmsfact2=1.0; dipfact2=1.0; }
                        }

                        if(varycharges==3) {rmsfact2=0.0; dipfact2=1.0; }

                        if(((rmspot*rmsfact2+diprtmp*dipfact2)<=(rmspot2*rmsfact2+dipr*dipfact2))&&(ignore==0)) {
                                if((fmod(itera,1000)==0.0)||(itera==itermax-1)) {if(DEBUG >0) fprintf(outfile," * ");
                                                                                 fflush(outfile); }
                                rmspot2=rmspot;
                                dipx_org=dipxtmp;
                                dipy_org=dipytmp;
                                dipz_org=dipztmp;
                                dipr=diprtmp;
                                /*   if(dipr<1.0e-3)dipfact2=0.0;
                                     else dipfact2=1.0;*/
                                for (imad = 0; imad < imadmax; imad++) dspot[imad]=dspottmp[imad];
                                if(varycharges==1) for(n=0; n<5; n++) {
                                                if(qpt[z2list[n]]>qmax) qmax=qpt[z2list[n]];
                                                if(qpt[z2list[n]]<qmin) qmin=qpt[z2list[n]];
                                        }
                                if(varycharges==2) for(n=0; n<2; n++) {
                                                if(qpt[z2list[n]]>qmax) qmax=qpt[z2list[n]];
                                                if(qpt[z2list[n]]<qmin) qmin=qpt[z2list[n]];
                                        }
                        }
                        else{ /* No improvement or ignore*/
                                if((fmod(itera,1000)==0.0)||(itera==itermax-1)) {
                                        if(DEBUG >0) fprintf(outfile,"   "); fflush(outfile);
                                }
                                /*
                                   if(varycharges==1)for(n=0;n<5;n++)qpt[z2list[n]]=qtmp[n];
                                 */
                                if(varycharges==2) for(n=0; n<2; n++) qpt[z2list[n]]=qtmp[n];
                                if(varycharges==3) for(n=0; n<2; n++) qpt[z2list[n]]=qtmp[n];
                                /*
                                   if(varycharges==1){
                                   xpt[z2list[0]]=xx0[0];
                                   ypt[z2list[0]]=xx0[1];
                                   zpt[z2list[0]]=xx0[2];
                                   }
                                 */
                                /* om dr maste fixas sa gor det har */

                        }

                        if(varycharges==2) vs=1;
                        if(varycharges==3) {
                                rmsfact2=0.0;
                                dipfact2=1.0;
                                vs=1;
                                if(diprtmp<1000000.0001 ) { /* Mattias har stalls dipgransen in!!!!*/
                                        varycharges=1; rmsfact2=1.0; dipfact2=0.0; ddq=0.001; vs=4;

                                }

                        }
                        if(do_lapack5 !=1) {
                                if(((vs==4)&&(diprtmp>=1000000.0001))) { /* Mattias har stalls dipgransen in!!!!*/
                                        varycharges=3; rmsfact2=0.0; dipfact2=1.0; ddq=0.01; vs=1;

                                }
                        }
                        if((fmod(itera,1000)==0.0)||(itera==itermax-1)) {
                                dq=0.0;
                                dipx=0.0; dipy=0.0; dipz=0.0;
                                for (ipt =0; ipt<iptmax; ipt++) {
                                        if (zone[ipt] >= 3) continue;
                                        dq+=qpt[ipt];
                                        dipx+=qpt[ipt]*xpt[ipt];
                                        dipy+=qpt[ipt]*ypt[ipt];
                                        dipz+=qpt[ipt]*zpt[ipt];
                                }
                                if(DEBUG >0)
                                        fprintf(outfile,"\tq=%4.2e (%4.3f,%4.3f)\trq=(%4.2e,%4.2e,%4.2e) (%f) (%1.0f %1.0f)\n",
                                                dq,qmin,qmax,dipx,dipy,dipz,sqrt(dipx*dipx+dipy*dipy+dipz*dipz),rmsfact2,dipfact2);
                                fflush(outfile);
                        }
                } /* iter */

                qmin=1.0e10;
                qmax=-1.0e10;
                for (ipt =0; ipt<iptmax; ipt++) {
                        if(qpt[ipt]<qmin) qmin=qpt[ipt];
                        if(qpt[ipt]>qmax) qmax=qpt[ipt];
                }

                if(DEBUG >0) fprintf(outfile,"q ranges:   %f  to  %f\n",qmin,qmax);

                rmspot = 0.;
                avepot = 0.;
                nomad  = 0;


                for (imad = 0; imad < imadmax; imad++) { /* loop over all cluster atoms to compute direct sum potentials*/
                        dspot[imad] = 0.;
                        for (ipt =0; ipt<iptmax; ipt++) { /* loop over all occupied sites */
                                if (zone[ipt] >= 3) continue;
                                r2 = (xpt[ipt]-xmad[imad]) * (xpt[ipt]-xmad[imad]) +
                                     (ypt[ipt]-ymad[imad]) * (ypt[ipt]-ymad[imad]) +
                                     (zpt[ipt]-zmad[imad]) * (zpt[ipt]-zmad[imad]);
                                if (r2 > 0.000001) dspot[imad] += evconv*qpt[ipt]/sqrt(r2);
                        }
                } /* end imad loop */

                for (imad = 0; imad < imadmax; imad++) { /* loop over all cluster atoms to compute direct sum potentials*/
                        if (qmad[imad] != 0.0) {

                                nomad++;
                                avepot += ewald[imad] - dspot[imad];
                                rmspot += (ewald[imad] - dspot[imad] ) * (ewald[imad] - dspot[imad] );
                        }
                } /* end loop over cluster atoms */
                rmspot = sqrt (rmspot/nomad);
                avepot = avepot/nomad;

        } /* varycharges == 1 || varycharges == 2 */

        /*
           pts=fopen("q.pts", "w");
           if(DEBUG > 0)
           for (ipt =0; ipt<iptmax; ipt++){
           fprintf(pts,"%25.20f  %25.20f  %25.20f  %25.20f\n",
           xpt[ipt],ypt[ipt],zpt[ipt],qpt[ipt]);
           }
           fclose(pts);
         */
        /*
         *
         *  End fix dip.mom..
         *
         */

        seedfile = fopen("seedfile","w");
        fprintf(seedfile, "%ld %ld", iseed1, iseed2);
        fclose(seedfile);
        if(DEBUG > 0)
                fprintf(outfile, "\nfinal results- cluster atoms:\n");

        rmspot = 0.;
        avepot =0.;
        nomad =0;

        if(DEBUG > 0) {
                fprintf(outfile, " Unit cell optimization results \n");
                for (imad=0; imad<imadmax; imad++) {
                        fprintf(outfile, "***     dspot %f V  ewald %f V  � %f (%e)V\n",
                                dspot[imad], ewald[imad], dspot[imad]-ewald[imad],dspot[imad]-ewald[imad]);
                }
        }

        /*
         *  Print outs.
         */

        for (imad = 0; imad < imadmax; imad++)
        { /* loop over all cluster atoms to compute direct sum potentials*/
                if (qmad[imad] != 0.0) {

                        nomad++;
                        avepot += ewald[imad] - dspot[imad];
                        rmspot += (ewald[imad] - dspot[imad] ) * (ewald[imad] - dspot[imad] );
                        if(DEBUG > 0) {
                                fprintf(outfile, " imad %ld x %.5f y %.5f z %.5f q %.4f\n",
                                        imad, xmad[imad], ymad[imad], zmad[imad], qmad[imad]);
                                fprintf(outfile, "     dspot %f eV  ewald %f eV  � %f eV\n",
                                        dspot[imad], ewald[imad], dspot[imad]-ewald[imad]);
                        }
                }
                else{
                        if(DEBUG > 0) {
                                fprintf(outfile,
                                        " imad %ld x %.5f y %.5f z %.5f q %.4f dspot %f V\n",
                                        imad, xmad[imad], ymad[imad], zmad[imad],
                                        qmad[imad], dspot[imad]);
                                fprintf(outfile, "     dspot %f V  ewald %f V  � %f V\n",
                                        dspot[imad], ewald[imad], dspot[imad]-ewald[imad]);
                        }

                        nomad++;
                        avepot += ewald[imad] - dspot[imad];
                        rmspot += (ewald[imad] - dspot[imad] ) * (ewald[imad] - dspot[imad] );
                }
        } /* end loop over cluster atoms */
        rmspot = sqrt (rmspot/nomad);
        avepot = avepot/nomad;



        fprintf(outfile,"\nNot only the random chk points:  rmspot=%e\n",rmspot);
        time(&now);
        date=localtime(&now);

        dx = 0.;  dy = 0.;  dz = 0.; qq =0.;
        dxx =0.;  dyy=0.;  dzz=0.;  radmom = 0.;
        for (in = 0; in<inmax; in++) {
                ucx[in] = uc1[in] * a1x + uc2[in] * a2x +uc3[in] * a3x;
                ucy[in] = uc1[in] * a1y + uc2[in] * a2y + uc3[in] * a3y;
                ucz[in] = uc1[in] * a1z + uc2[in] * a2z + uc3[in] * a3z;
                dx += ucx[in] * quc[in];
                dy += ucy[in] * quc[in];
                dz += ucz[in] * quc[in];
                dxx += ucx[in] * ucx[in] * quc[in];
                dyy += ucy[in] * ucy[in] * quc[in];
                dzz += ucz[in] *ucz[in] * quc[in];
                qq += quc[in];
        } /* end in loop */
        d2 = dx*dx + dy*dy + dz*dz;
        radmom = sqrt(dxx*dxx + dyy*dyy + dzz*dzz);

        for(in= 0; in<imadmax; in++) {
                antala[in] = 0.0;
                ewalda[in] = 0.0;
                dspotFAa[in] = 0.0;
                dspota[in] = 0.0;
                dspotSa[in] = 0.0;
                qmada[in] = 0.0;
                rmsmad1a[in] = 0.0;
                rmsmad2a[in] = 0.0;
                rmsmadSa[in] = 0.0;
                rmsewalda[in] = 0.0;
                ave1a[in] = 0.0;
                ave2a[in] = 0.0;
                aveSa[in] = 0.0;
        }

        avec1=0.0;
        avec2=0.0;
        avecS=0.0;
        rmsc1=0.0;
        rmsc2=0.0;
        rmscS=0.0;
        nave = 0;
        for (in = 0; in<imadmax; in++) {
                if(qmad[in] == 0.0) continue;
                nave++;
                itemp=typmad[in];
                avec1 += dspotFA[in]-ewald[in];
                avec2 += dspot[in]-ewald[in];
                avecS += dspotS[in]-ewald[in];
        } /* end in loop */
        if (nave > 0) {
                avec1 = avec1/nave;
                avec2 = avec2/nave;
                avecS = avecS/nave;
        } /* end nave > 0 case */

        if (nave > 1) {
                for (in = 0; in<imadmax; in++) {
                        if(qmad[in] == 0.0) continue;
                        temp = avec1 - (dspotFA[in]-ewald[in]);
                        rmsc1 += temp * temp;
                        temp = avec2 - (dspot[in]-ewald[in]);
                        rmsc2 += temp * temp;
                        temp = avecS - (dspotS[in]-ewald[in]);
                        rmscS += temp * temp;
                } /* end in loop */
                rmsc1 = sqrt(rmsc1/(nave-1.0));
                rmsc2 = sqrt(rmsc2/(nave-1.0));
                rmscS = sqrt(rmscS/(nave-1.0));
        } /* end nave > 1 case */

        ave1=0.0;
        ave2=0.0;
        aveS=0.0;
        rms1=0.0;
        rms2=0.0;
        rmsS=0.0;
        nave =0;
        for (in = 0; in<imadmax; in++) {
                if(qmad[in] != 0.0) continue;
                nave++;
                itemp=typmad[in];
                ave1 += dspotFA[in]-ewald[in];
                ave2 += dspot[in]-ewald[in];
                aveS += dspotS[in]-ewald[in];
        } /* end in loop */
        if (nave > 0) {
                ave1 = ave1/nave;
                ave2 = ave2/nave;
                aveS = aveS/nave;
        }

        if (nave > 1) {
                for (in = 0; in<imadmax; in++) {
                        if(qmad[in] != 0.0) continue;
                        temp = ave1 - (dspotFA[in]-ewald[in]);
                        rms1 += temp * temp;
                        temp = ave2 - (dspot[in]-ewald[in]);
                        rms2 += temp * temp;
                        temp = aveS - (dspotS[in]-ewald[in]);
                        rmsS += temp * temp;
                } /* end in loop */
                rms1 = sqrt(rms1/(nave-1.0));
                rms2 = sqrt(rms2/(nave-1.0));
                rmsS = sqrt(rmsS/(nave-1.0));
        } /* end if nave >1 case */

        typmada[0] = typmad[0];
        antala[0] = 1.0;
        ewalda[0] = ewald[0];
        dspotFAa[0] = dspotFA[0];
        dspotSa[0] = dspotS[0];
        dspota[0] = dspot[0];
        aveSa[0] = dspotS[0];
        ave1a[0] = dspotFA[0];
        ave2a[0] = dspot[0];

        nn=0;
        for(in=1; in<imadmax; in++) {
                contin = 0.0;
                for(n = 0; n<=nn; n++) {
                        if(typmad[in] == typmada[n]) { /*loop over atom types already seen */
                                antala[n] += 1.0;
                                ewalda[n] += ewald[in];
                                dspotSa[n] += dspotS[in];
                                dspotFAa[n] += dspotFA[in];
                                dspota[n] += dspot[in];
                                qmada[n] = qmad[in];
                                aveSa[n] += dspotS[in];
                                ave1a[n] += dspotFA[in];
                                ave2a[n] += dspot[in];
                                contin = 1.0;
                        }
                } /* end n loop */
                if(contin == 0.0) { /* new atom type encountered */
                        nn++;
                        typmada[nn] = typmad[in];
                        antala[nn] += 1.0;
                        ewalda[nn] += ewald[in];
                        dspotSa[nn] += dspotS[in];
                        dspotFAa[nn] += dspotFA[in];
                        dspota[nn] += dspot[in];
                        qmada[nn] = qmad[in];
                        aveSa[nn] += dspotS[in];
                        ave1a[nn] += dspotFA[in];
                        ave2a[nn] += dspot[in];
                } /* end case of new atom type */
        } /* end in loop */

        for(n=0; n<=nn; n++) { /* loop over atom types */
                aveSa[n] = aveSa[n]/antala[n];
                ave1a[n] = ave1a[n]/antala[n];
                ave2a[n] = ave2a[n]/antala[n];
                ewalda[n] = ewalda[n]/antala[n];
                dspotFAa[n] = dspotFAa[n]/antala[n];
                dspotSa[n] = dspotSa[n]/antala[n];
                dspota[n] = dspota[n]/antala[n];
        }

        for(in=0; in<imadmax; in++) {
                for(n = 0; n<=nn; n++) {
                        if(typmad[in]==typmada[n]) {
                                rmsmadSa[n] += (dspotS[in]-ave1a[n])*(dspotS[in]-ave1a[n]);
                                rmsmad1a[n] += (dspotFA[in]-ave1a[n])*(dspotFA[in]-ave1a[n]);
                                rmsmad2a[n] += (dspot[in]-ave2a[n])*(dspot[in]-ave2a[n]);
                                rmsewalda[n] += (ewald[in]-ewalda[n])*(ewald[in]-ewalda[n]);
                        }
                } /* end n loop */
        } /* end in loop */

        for(n=0; n<=nn; n++) {
                if(antala[n]>1.0) {
                        rmsmadSa[n] = sqrt(rmsmadSa[n]/(antala[n]-1.0));
                        rmsmad1a[n] = sqrt(rmsmad1a[n]/(antala[n]-1.0));
                        rmsmad2a[n] = sqrt(rmsmad2a[n]/(antala[n]-1.0));
                        rmsewalda[n] = sqrt(rmsewalda[n]/(antala[n]-1.0));
                }
                else{
                        rmsmadSa[n] = 0.0;
                        rmsmad1a[n] = 0.0;
                        rmsmad2a[n] = 0.0;
                        rmsewalda[n] = 0.0;
                }
        } /* end n loop */

        /* output section for atom type and Ewald */

        for (n = 0; n<nn; n++) {
                if(rmsewalda[n] > 0.0001) {
                        itemp=typmada[n];
                        fprintf(outfile, "WARNING rms Ewald for %s = %f V\n",
                                typetab[itemp], rmsewalda[n]);
                }
        }

        ifill =0.;
        dx=0.;
        dy=0.;
        dz=0.;
        qsum =0.;
        xmin = 10000000000.;
        ymin = 10000000000.;
        zmin = 10000000000.;
        xmax = -10000000000.;
        ymax = -10000000000.;
        zmax = -10000000000.;
        for (ipt=0; ipt<iptmax; ipt++)
                if (zone[ipt] < 3) { /* loop over all occupied sites */
                        ifill++;
                        dx += xpt[ipt]*qpt[ipt];
                        dy += ypt[ipt]*qpt[ipt];
                        dz += zpt[ipt]*qpt[ipt];
                        qsum += qpt[ipt];
                        if (xpt[ipt] < xmin) xmin = xpt[ipt];
                        if (ypt[ipt] < ymin) ymin = ypt[ipt];
                        if (zpt[ipt] < zmin) zmin = zpt[ipt];
                        if (xpt[ipt] > xmax) xmax = xpt[ipt];
                        if (ypt[ipt] > ymax) ymax = ypt[ipt];
                        if (zpt[ipt] > zmax) zmax = zpt[ipt];
                }
        if(DEBUG > 0) {
                fprintf(outfile,"dipole: x %.8f y %.8f z %.8f qsum %f\n",dx, dy, dz, qsum);
                fprintf(outfile, "lattice xmin %f xmax %f\n",xmin, xmax);
                fprintf(outfile, "lattice ymin %f ymax %f\n",ymin, ymax);
                fprintf(outfile, "lattice zmin %f zmax %f\n",zmin, zmax);
                fprintf(outfile," %ld sites %ld atoms\n",iptmax, ifill);
        }

        for (ipt=0; ipt<iptmax; ipt++) { /* loop over all occupited sites to match with madelung rms atoms */
                if (zone[ipt] > 2) continue;
                ifatom = 0;
                for (imad=0; imad< imadmax; imad++) {
                        if ( (fabs(xpt[ipt]-xmad[imad]) < 0.001) &&
                             (fabs(ypt[ipt]-ymad[imad]) < 0.001) &&
                             (fabs(zpt[ipt]-zmad[imad]) < 0.001) ) {
                                ifatom = 1;
                                break;
                        }
                } /* end imad loop */

                if (ifatom == 0) {
                        n++;
                        fprintf(ptsfilenw, "  bq%ld %25.20f %25.20f %25.20f charge %25.20f\n",
                                n+1+imadmax,xpt[ipt],ypt[ipt],zpt[ipt],qpt[ipt]);
                }
                if (ifatom != 1) fprintf(ptsfilefro, "%25.20f %25.20f %25.20f %25.20f\n",
                                         xpt[ipt],  ypt[ipt], zpt[ipt], qpt[ipt]);
                if (ifatom != 1) fprintf(ptsfile, "%25.20f\t%25.20f\t%25.20f\t%25.20f\n", qpt[ipt],
                                         xpt[ipt],  ypt[ipt], zpt[ipt]);
        } /* end ipt loop */

        fclose(ptsfile);
        fclose(ptsfilenw);
        fclose(ptsfilefro);

        rmspot = 0.;
        nomad =0;
        avepot =0.;
        for (imad=0; imad<imadmax; imad++) {
                if (qmad[imad] != 0.0) continue;
                nomad++;
                avepot += ewald[imad] - dspot[imad];
                rmspot += (ewald[imad] - dspot[imad] ) * (ewald[imad] - dspot[imad] );
        }
        rmspot = sqrt (rmspot/nomad);
        avepot = avepot/nomad;

        fprintf(outfile,"Only the random chk points:      rmspot=%e\n\n",rmspot);

        seedfile = fopen("seedfile","w");
        fprintf(seedfile, "%ld %ld", iseed1, iseed2);
        fclose(seedfile);

        if(DEBUG > 0) fprintf(outfile, "\n\nSVAR:  %f\t%f\t%f\t%f (%fmV)\n\n\n",
                              avepotuc2,rmspotuc2,avepot,rmspot,rmspot*1000.0);
        /*
           strncpy(resname,fileroot,50);
           strncat(resname,".pub2",20);
           pubfile = fopen(resname,"w");
           fprintf(pubfile,"\nTable B.X.1 Non-equivalent atoms\n");
           fprintf(pubfile,"Atom\tq\tu1\tu2\tu3\tVe(V)\n");
         */
        for(it=0; it<itmax; it++) {
                for(n=0; n<=nn; n++) {
                        if(qmada[n]!=0.0) {
                                itemp=typmada[n];
                                /*
                                   if(strncmp(typetab[itemp],typetab[it],strlen(typetab[it]))==0.0 &&
                                   strlen(typetab[it]) == strlen(typetab[itemp]))
                                   fprintf(pubfile,"%s\t%2.0f\t\t\t\t%8.4f\n",
                                    typetab[itemp],qmada[n], ewalda[n]);*/
                        }
                }
        }

        fprintf(outfile,"\n");
        for(it=0; it<itmax; it++) {
                fprintf(outfile, "Unit cell has %ld atoms of type %s\n",
                        nuctype[it], typetab[it]);
        }

        for(r=10.0; r>0.0; r-=1.0) {
                rr=0.0;
                for(it=0; it<itmax; it++) rr+= fabs(fmod(nuctype[it],r));
                if(rr==0.0) {
                        nrformulaunitsincell=(integer)(r);
                        r=0;
                }
        }

        fprintf(outfile, "This gives Z=%ld\n",nrformulaunitsincell);

        ewaldmad=0.0;
        for(n=0; n<=nn; n++) {
                if(qmada[n]!=0.0) {
                        itemp=typmada[n];
                        ewaldmad += nuctype[itemp]*qmada[n]*ewalda[n]/nrformulaunitsincell;

                }
        }
        ewaldmad = -ewaldmad*scalemad/0.529177/2.0/27.2116;

        fprintf(outfile,"\nMadelung constant %.6f for ",ewaldmad);
        fprintf(outfile,"R(%s-%s) %.6f �\n", typuc[cnr1],typuc[cnr2],scalemad);


        fprintf(outfile, "\n\nseed1 = %ld, seed2 = %ld\n", iseed1, iseed2);
        fprintf(outfile, "\nStarted:  %s  %s",fileroot, strdate1);
        time(&now);
        date=localtime(&now);
        fprintf(outfile, "Finished: %s  %s",fileroot, asctime(date));
        fprintf(outfile, "END\n");

        fclose(outfile);

        return n;
} /* end of main */
