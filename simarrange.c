/*Copyright (c) 2013 Kliment Yanev

This file is part of simarrange.

Simarrange is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
#include <admesh/stl.h>
#include <stdio.h>
#include <sys/stat.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "utlist.h"
#include <dirent.h>
#include <argtable2.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <omp.h>


//#define PARALLEL
#if defined PARALLEL
#include <omp.h>
#endif
#define FILENAME_LEN 350

inline int max ( int a, int b ) { return a > b ? a : b; }
inline int min ( int a, int b ) { return a < b ? a : b; }

typedef struct img_list{
    IplImage *image;
    long area;
    int count;
    int done; // counts the number of copies already plated
    int countinplate; // counts the number of copies in the latest plate
    int *x;
    int *y;
    int *rotangle;
    int *plate;
    char filename[FILENAME_LEN];
    stl_file *stl;
    struct img_list *prev,*next;
} img_list;


void sqspiral(int n, int *i, int *j)
{
  double x = sqrt(n+1);
  int xi = floor(x);
  int xisq = xi * xi;
  if (xi % 2 == 0)  // even
  { 
    int xi2 = xi / 2;
    if (x == xi)
    {
      *i = -xi2 + 1;
      *j = xi2;
    }
    // top left number is xisq - xi - 1, bottom left n - xisq - 1
    else if (n <= xisq + xi + 1)
    {
      *i = -xi2;
      *j = xi2 - (n - xisq - 1); 
    }
    else 
    {
      *i = -xi2 + (n - xisq - xi - 1);
      *j = -xi2;
    }
  }
  else  // odd 
  {
    if (x == xi)
    {
      int xi2 = (xi-1)/2;
      *i = xi2;
      *j = -xi2;
    } 
    else 
    {
      int xi2 = (xi + 1) / 2;
      if (n <= xisq + xi + 1)
      {
        *i = xi2;
        *j = -xi2 + (n - xisq);
      } 
      else
      {
        *i = xi2 - (n - xisq - xi - 1);
        *j = xi2;
      }
    }
  }
}

static void
stl_put_little_int(FILE *fp, int value_in)
{
  int new_value;
  union 
    {
      int  int_value;
      char char_value[4];
    } value;
  
  value.int_value = value_in;
  
  new_value  = value.char_value[0] & 0xFF;
  new_value |= (value.char_value[1] & 0xFF) << 0x08;
  new_value |= (value.char_value[2] & 0xFF) << 0x10;
  new_value |= (value.char_value[3] & 0xFF) << 0x18;
  fwrite(&new_value, sizeof(int), 1, fp);
}

static void
stl_put_little_float(FILE *fp, float value_in)
{
  int new_value;
  union 
    {
      float float_value;
      char  char_value[4];
    } value;
  
  value.float_value = value_in;
  
  new_value  = value.char_value[0] & 0xFF;
  new_value |= (value.char_value[1] & 0xFF) << 0x08;
  new_value |= (value.char_value[2] & 0xFF) << 0x10;
  new_value |= (value.char_value[3] & 0xFF) << 0x18;
  fwrite(&new_value, sizeof(int), 1, fp);
}


void
stl_write_binary_block(stl_file *stl, FILE *fp)
{
  int i;
  for(i = 0; i < stl->stats.number_of_facets; i++)
    {
      stl_put_little_float(fp, stl->facet_start[i].normal.x);
      stl_put_little_float(fp, stl->facet_start[i].normal.y);
      stl_put_little_float(fp, stl->facet_start[i].normal.z);
      stl_put_little_float(fp, stl->facet_start[i].vertex[0].x);
      stl_put_little_float(fp, stl->facet_start[i].vertex[0].y);
      stl_put_little_float(fp, stl->facet_start[i].vertex[0].z);
      stl_put_little_float(fp, stl->facet_start[i].vertex[1].x);
      stl_put_little_float(fp, stl->facet_start[i].vertex[1].y);
      stl_put_little_float(fp, stl->facet_start[i].vertex[1].z);
      stl_put_little_float(fp, stl->facet_start[i].vertex[2].x);
      stl_put_little_float(fp, stl->facet_start[i].vertex[2].y);
      stl_put_little_float(fp, stl->facet_start[i].vertex[2].z);
      fputc(stl->facet_start[i].extra[0], fp);
      fputc(stl->facet_start[i].extra[1], fp);
    }
}

void
stl_translate_rel(stl_file *stl, float x, float y, float z)
{
  int i;
  int j;
  
  for(i = 0; i < stl->stats.number_of_facets; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  stl->facet_start[i].vertex[j].x += x;
	  stl->facet_start[i].vertex[j].y += y;
	  stl->facet_start[i].vertex[j].z += z;
	}
    }
  stl->stats.max.x += x;
  stl->stats.max.y += y;
  stl->stats.max.z += z;
  stl->stats.min.x += x;
  stl->stats.min.y += y;
  stl->stats.min.z += z;
}



void add_stl(char *filename, int count, int width, int height, img_list **shapes){
    img_list *elt;
    DL_FOREACH(*shapes,elt) {
        if (strncmp(filename, elt->filename, FILENAME_LEN) == 0) {
            elt->count += count;
            elt->x = (int*) realloc(elt->x, sizeof(int)*elt->count);
            elt->y = (int*) realloc(elt->y, sizeof(int)*elt->count);
            elt->rotangle = (int*) realloc(elt->rotangle, sizeof(int)*elt->count);
            elt->plate = (int*) realloc(elt->plate, sizeof(int)*elt->count);
            return;
        }
    }
    stl_file *s=(stl_file *)malloc(sizeof(stl_file));
    memset(s,0,sizeof(stl_file));
    stl_open(s,filename);
    stl_check_facets_exact(s);
    s->stats.facets_w_1_bad_edge = (s->stats.connected_facets_2_edge - s->stats.connected_facets_3_edge);
    s->stats.facets_w_2_bad_edge = (s->stats.connected_facets_1_edge - s->stats.connected_facets_2_edge);
    s->stats.facets_w_3_bad_edge = (s->stats.number_of_facets - s->stats.connected_facets_1_edge);
    height=(int)(2*height);
    width=(int)(2*width);
    stl_translate(s, ((float)width/2)-(s->stats.max.x-s->stats.min.x)/2.0, ((float)height/2)-(s->stats.max.y-s->stats.min.y)/2.0, 0.0);
    
    IplImage* img = NULL;
    img=cvCreateImage(cvSize(width,height), IPL_DEPTH_8U, 1);
    cvZero(img);
    unsigned int i;
    for (i = 0; i < s->stats.number_of_facets; i++) {
        CvPoint points[3]={cvPoint((int)(s->facet_start[i].vertex[0].x),(int)(s->facet_start[i].vertex[0].y)),
                           cvPoint((int)(s->facet_start[i].vertex[1].x),(int)(s->facet_start[i].vertex[1].y)),
                           cvPoint((int)(s->facet_start[i].vertex[2].x),(int)(s->facet_start[i].vertex[2].y))
                       };
        
        cvFillConvexPoly(img,points,3,cvScalarAll(127),8,0);
    }

    img_list *e=(img_list*)malloc(sizeof(img_list));
    e->image=img;
    e->x=(int*)malloc(sizeof(int)*count);
    e->y=(int*)malloc(sizeof(int)*count);
    e->rotangle=(int*)malloc(sizeof(int)*count);
    e->plate=(int*)malloc(sizeof(int)*count);
    for (i = 0; i < count; ++i)
      e->plate[i] = -1;
    e->stl=s;
    e->count=count;
    e->done = 0;
    e->area=cvCountNonZero(e->image);
    strcpy(e->filename,filename);
    DL_APPEND(*shapes,e);
    
    
}

int areacmp(img_list *a, img_list *b) {
    return (b->area)-(a->area);
}

int dl_count(img_list *l){
    int c=0;
    img_list *elt;
    DL_FOREACH(l,elt){
        c += elt->count - elt->done;
    }
    return c;
}

int add_files(struct arg_file *arg, int w, int h, img_list **shapes, int withrepeat) {
    int ifile;
    char *indir=NULL;
    struct stat filestat;
    for(ifile = 0; ifile < arg->count; ifile++){
        char filename[FILENAME_LEN];
        strncpy(filename, arg->filename[ifile], FILENAME_LEN);
        int copies = 1;
        if (withrepeat) {
            int filename_len = strlen(filename);
            int k;
            for (k = filename_len - 1; k >= 0; --k) {
                if (filename[k] == '+') {
                    copies = strtol(filename + k + 1, NULL, 10);
                    filename[k] = '\0';
                    break;
                }
            }
        }
        if( stat(filename, &filestat) == 0 ){
            if( filestat.st_mode & S_IFREG )
                add_stl(filename, copies, w, h, shapes);
            if( filestat.st_mode & S_IFDIR ){
                indir= filename;
                DIR *dir;
                struct dirent *ent;
                if ((dir = opendir (indir)) != NULL) {
                    while ((ent = readdir (dir)) != NULL) {
                        int i;
                        char d[FILENAME_LEN];
                        strcpy(d,ent->d_name);
                        for(i = 0; d[i]; i++)
                            d[i] = tolower(d[i]);
                        if(strstr(d,".stl")!=0){
                            char f[FILENAME_LEN];
                            f[0]=0;
                            strcat(f,indir);
                            strcat(f+strlen(f),"/");
                            strcat(f+strlen(f),ent->d_name);
                            
                            add_stl(f, copies, w, h, shapes);
                        }
                            //printf ("%s\n", ent->d_name);
                    }
                    closedir (dir);
                } else {
                    fprintf(stderr, "Input directory not found %s\n",indir);
                    return EXIT_FAILURE;
                }        
            }
        }else{
            fprintf(stderr, "Could not access %s\n",filename);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}



int search(int rotangle, int posstep, int w, int h, int firstpassed, int middle, int *minxpos, int *minypos, int *minrotangle, int *mincentricords,CvPoint2D32f center, IplImage *itmp, IplImage *rpatch, IplImage *img,IplImage *testfit,CvMat *rot){
    int xpos,ypos,placed=0;
    cvWarpAffine(itmp,rpatch,cv2DRotationMatrix(center, rotangle, 1.0, rot),CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS, cvScalarAll(0) );
    if(firstpassed && middle){
        int centricords=0;
        for(centricords=0;centricords<*mincentricords;centricords++){
            xpos=ypos=0;
            sqspiral(centricords, &xpos, &ypos);
            xpos*=posstep;
            ypos*=posstep;
            
            xpos=max(0,min(xpos+w/2,w-1));
            ypos=max(0,min(ypos+h/2,h-1));
            cvSetImageROI(rpatch, cvRect(w-xpos,h-ypos,w,h));
            cvAnd(rpatch, img, testfit, NULL);
            
            if(!cvCountNonZero(testfit)){
                int prec=cvCountNonZero(img);
                cvAdd(rpatch, img, testfit, NULL);
                if(prec!=cvCountNonZero(testfit)){
                    #if defined PARALLEL
                    #pragma omp critical
                    #endif
                    {
                        if(centricords<*mincentricords){
                            *minxpos=xpos;
                            *minypos=ypos;
                            *minrotangle=rotangle;
                            *mincentricords=centricords;
                            placed=1;
                        }
                    }

                    
                }
            }
            cvResetImageROI(rpatch);
        }
    }else{
        for(ypos=1;ypos<*minypos;ypos+=posstep){
            for(xpos=1;xpos<*minxpos;xpos+=posstep){
                cvSetImageROI(rpatch, cvRect(w-xpos,h-ypos,w,h));
                cvAnd(rpatch, img, testfit, NULL);
                if(!cvCountNonZero(testfit)){
                    int prec=cvCountNonZero(img);
                    cvAdd(rpatch, img, testfit, NULL);
                    if(prec!=cvCountNonZero(testfit)){
                        #if defined PARALLEL
                        #pragma omp critical
                        #endif
                        {
                            if (xpos+ypos < (*minxpos+*minypos))
                            {
                                *minxpos=xpos;
                                *minypos=ypos;
                                *minrotangle=rotangle;
                                placed=1;
                            }
                        }
                    }
                }
                cvResetImageROI(rpatch);
                
            }
        }
    }
    return placed;
}

int main(int argc, char** argv){
    int w=200,h=200;
    int spacing=1;
    int rotstep=10;
    int posstep=5;
    int quiet=0;
    int c;
    struct arg_int  *aw  = arg_int0("x","width",NULL,              "plate width in mm (default is 200)");
    struct arg_int  *ah  = arg_int0("y","height",NULL,              "plate height in mm (default is 200)");
    struct arg_int  *as  = arg_int0("s","spacing",NULL,              "spacing between parts in mm (default is 1)");
    struct arg_int  *ar  = arg_int0("r","rotstep",NULL,              "rotation step when searching (default 10 degrees)");
    struct arg_int  *ap  = arg_int0("p","posstep",NULL,              "positional step when searching (default 5mm)");
    struct arg_lit  *ac  = arg_lit0("c","circle",              "circular print area with diameter given by -x");
    struct arg_lit  *acorigin  = arg_lit0("m","middle",              "place objects from middle of build area out");
    #ifdef PARALLEL
    struct arg_int  *athreads  = arg_int0("j","threads",NULL,         "number of threads (default is to use as much as possible, set to 1 to disable multithreading)");
    #endif
    struct arg_lit  *adryrun  = arg_lit0("d","dryrun",              "only do a dry run, computing placement but not producing any output file");
    struct arg_lit  *ahelp  = arg_lit0("h","help",              "display this help message");
    struct arg_lit  *aquiet  = arg_lit0("q","quiet",              "silence information messages");
    struct arg_str  *aodir = arg_str0("o","outputdir",NULL,  "output directory (default .)");
    struct arg_file  *arepeat = arg_filen(NULL,"repeat",NULL,0,argc+2,  "add a given number of copies of the input file or dir by specifying filepath+count");
    struct arg_file  *ainfile = arg_filen(NULL,NULL,NULL,0,argc+2,  "input file or dir (any number allowed)");
    struct arg_end  *end      = arg_end(20);
    void* argtable[] = {aw,ah,as,ar,ap,ac,acorigin,aodir,ainfile,arepeat,
                        #ifdef PARALLEL
                        athreads,
                        #endif
                        adryrun,aquiet,ahelp,end};
    
    int nerrors;
    nerrors = arg_parse(argc,argv,argtable);
    if(ahelp->count){
        printf("Usage: %s", argv[0]);
        arg_print_syntax(stdout, argtable, "\n");
        printf("Options:\n");
        arg_print_glossary(stdout, argtable, " %-25s %s\n");
        return EXIT_SUCCESS;
    }
    if (nerrors > 0)
        {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stderr,end,argv[0]);
        return EXIT_FAILURE;
        }
   if (arepeat->count == 0 && ainfile->count == 0) {
        fprintf(stderr, "%s: please specify one or more input file or directory\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    if(aw->count){
        w=aw->ival[0];
    }
    if(ah->count){
        h=ah->ival[0];
    }
    if(ac->count){
        h=w;
    }
    if(as->count){
        spacing=as->ival[0];
    }
    if(ar->count){
        rotstep=ar->ival[0];
    }
    if(ap->count){
        posstep=ap->ival[0];
    }
    if(aquiet->count){
        quiet = 1;
    }
    if(adryrun->count && !quiet){
        printf("Running in dry run mode (no output file will be produced)\n");
    }
    #ifdef PARALLEL
    if(athreads->count) {
        if (athreads->ival[0] > 0){
            int numthreads = athreads->ival[0];
            omp_set_dynamic(0);
            omp_set_num_threads(numthreads);
        } else {
            fprintf(stderr, "%s: number of threads cannot be lower than 0\n", argv[0]);
            return EXIT_FAILURE;
        }
    }
    #endif
    
    
    
    img_list *shapes=NULL;
    img_list *curplate=NULL;
    int platecount=0;
    int ret;
    if(ainfile->count){
        ret = add_files(ainfile, w, h, &shapes, 0);
        if (ret != EXIT_SUCCESS)
            return ret;
    }
    if(arepeat->count){
        ret = add_files(arepeat, w, h, &shapes, 1);
        if (ret != EXIT_SUCCESS)
            return ret;
    }
    
    DL_SORT(shapes, areacmp);
    img_list *elt,*tmp;
    
    int plate=0;
    
    IplImage* img = NULL;
    img=cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 1);
    CvMat* rot = NULL;
    rot=cvCreateMat(2,3, CV_32FC1);
    IplImage* rpatch = NULL;
    rpatch=cvCreateImage(cvSize((w*2),(h*2)), IPL_DEPTH_8U, 1);
    IplImage* testfit = NULL;
    testfit=cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 1);
    IplImage* itmp = NULL;
    itmp=cvCreateImage(cvSize(2*w,2*h), IPL_DEPTH_8U, 1);
    CvPoint2D32f center=cvPoint2D32f((w),(h));
    
    unsigned copy;    
    while(dl_count(shapes)){
        if (!quiet) printf("Generating plate %d\n",plate);
        cvZero(img);
        cvLine(img, cvPoint(0,0), cvPoint(w-1,0), cvScalarAll(127), 1, 8, 0);
        cvLine(img, cvPoint(w-1,0), cvPoint(w-1,h-1), cvScalarAll(127), 1, 8, 0);
        cvLine(img, cvPoint(w-1,h-1), cvPoint(0,h-1), cvScalarAll(127), 1, 8, 0);
        cvLine(img, cvPoint(0,h-1), cvPoint(0,0), cvScalarAll(127), 1, 8, 0);
        if(ac->count){
            cvRectangle(img, cvPoint(0,0), cvPoint(w-1,h-1), cvScalarAll(127), CV_FILLED, 8, 0);
            cvCircle(img, cvPoint(w/2,h/2), (w/2)-1, cvScalarAll(0), CV_FILLED, 8, 0);
        }
        cvZero(rpatch);
        cvZero(testfit);
        cvZero(itmp);
        int firstpassed=0, placed=0;
        DL_FOREACH(shapes,elt) {
            for (copy = elt->done; copy < elt->count; copy++) {
                placed=0;
                //printf("File: %s\n",elt->filename);
                cvCopy(img, testfit, NULL);
                int xpos=1, ypos=1, rotangle=0, minxpos=w-1, minypos=h-1, minrotangle=0, mincentricords=max(w-1,h-1)*max(w-1,h-1)/(posstep*posstep);
                cvDilate(elt->image,itmp,NULL,spacing);
                #ifdef PARALLEL
                #pragma omp parallel for
                for(rotangle=0;rotangle<360;rotangle+=rotstep){
                    //rpatch,rot,testfit
                    CvMat* rot = NULL;
                    rot=cvCreateMat(2,3, CV_32FC1);
                    IplImage* rpatch = NULL;
                    rpatch=cvCreateImage(cvSize((w*2),(h*2)), IPL_DEPTH_8U, 1);
                    IplImage* testfit = NULL;
                    testfit=cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 1);
                    placed|=search(rotangle, posstep, w, h, firstpassed, acorigin->count, &minxpos, &minypos, &minrotangle, &mincentricords,center, itmp, rpatch, img,testfit,rot);
                    cvReleaseImage(&testfit);
                    cvReleaseImage(&rpatch);
                    cvReleaseMat(&rot);
                }
                #else
                for(rotangle=0;rotangle<360;rotangle+=rotstep){
                    placed|=search(rotangle, posstep, w, h, firstpassed, acorigin->count, &minxpos, &minypos, &minrotangle, &mincentricords,center, itmp, rpatch, img,testfit,rot);
                }
                #endif
                if(!firstpassed && acorigin->count){
                    int centricords=0;
                    for(centricords=0;centricords<max(w-1,h-1)*max(w-1,h-1);centricords++){
                        xpos=ypos=0;
                        sqspiral(centricords, &xpos, &ypos);
                        xpos=max(0,min(xpos+w/2,w-1));
                        ypos=max(0,min(ypos+h/2,h-1));
                        cvSetImageROI(rpatch, cvRect(w-xpos,h-ypos,w,h));
                        cvAnd(rpatch, img, testfit, NULL);
                        if(!cvCountNonZero(testfit)){
                            int prec=cvCountNonZero(img);
                            cvAdd(rpatch, img, testfit, NULL);
                            if(prec!=cvCountNonZero(testfit)){
                                minxpos=xpos;
                                minypos=ypos;
                                placed=1;
                                cvResetImageROI(rpatch);
                                break;
                            }
                        }
                        
                    }
                }
                if(placed){
                    if (!quiet) printf("File: %s minx: %d, miny: %d, minrot: %d\n",elt->filename, minxpos, minypos, minrotangle);
                    cvWarpAffine(elt->image,rpatch,cv2DRotationMatrix(center, minrotangle, 1.0, rot),CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS, cvScalarAll(0) );
                    cvSetImageROI(rpatch, cvRect(w-minxpos,h-minypos,w,h));
                    cvAdd(rpatch,img,testfit,NULL);
                    cvCopy(testfit,img,NULL);
                    cvResetImageROI(rpatch);
                    elt->x[elt->done]=minxpos;
                    elt->y[elt->done]=minypos;
                    elt->rotangle[elt->done]=minrotangle;
                    elt->plate[elt->done]=plate;
                    elt->done++;
                    platecount++;
                }else{
                    if (!quiet) printf("SKIP: %s skipped for this plate\n",elt->filename);
                    if (!firstpassed) {
                        fprintf(stderr, "Could not fit this file as the first item of the plate in any tested orientation! It might be too large for the print area.\n");
                        return EXIT_FAILURE;
                    }
                    break;
                }
                firstpassed=1;
            }
        } // end of the DL_FOREACH(shapes,elt) { loop

        if(!platecount){
            fprintf(stderr, "The files skipped in the last stage do not fit on plate in any tested orientation! They might be too large for the print area.\n");
            return EXIT_FAILURE;
        }

        int maxextraplates = INT_MAX;
        DL_FOREACH(shapes,elt){
            elt->countinplate = 0;
            for (copy = 0; copy < elt->done; ++copy)
            {
                if (elt->plate[copy] == plate) {
                    elt->countinplate++;
                }
            }
            if (elt->countinplate > 0) {
                int extraplates = (elt->count - elt->done) / elt->countinplate;
                maxextraplates = min(maxextraplates, extraplates);
            }
        }
        if (maxextraplates > 0) {
            if (!quiet) printf("Making %d duplicates of plate %d\n", maxextraplates, plate);
            DL_FOREACH(shapes,elt){
                elt->done += elt->countinplate * maxextraplates;
            }
        }
        
        if (!adryrun->count)
        {
            char tmpfn[512], outdir[512];
            char imagefn[1024], stlfn[1024];
            outdir[0]=0;
            if(aodir->count){
                mkdir(aodir->sval[0],0777);
                strcpy(imagefn,aodir->sval[0]);
                strcat(imagefn,"/");
            }
 
            sprintf(tmpfn,"plate%02d.png",plate);
            imagefn[0]=0;
            strcat(imagefn,outdir);
            strcat(imagefn,tmpfn);
            cvFlip(img,NULL,0);
            cvSaveImage(imagefn,img,0);
            sprintf(tmpfn,"plate%02d.stl",plate);
            stlfn[0]=0;
            strcat(stlfn,outdir);
            strcat(stlfn,tmpfn);
            FILE *fp;
            int i;
            fp = fopen(stlfn, "w");
            if(fp == NULL){
              fprintf(stderr, "Could not open output file %s\n",stlfn);
              return EXIT_FAILURE;
            }
            fprintf(fp, "%s", stlfn);
            for(i = strlen(stlfn); i < LABEL_SIZE; i++) putc(0, fp);
            fseek(fp, LABEL_SIZE, SEEK_SET);
            
            int totalfacets=0;
            DL_FOREACH(shapes,elt){
                for (copy = 0; copy < elt->done; ++copy)
                {
                    if(elt->plate[copy]==plate){
                        totalfacets+=elt->stl->stats.number_of_facets;
                    }
                }
            }
            stl_put_little_int(fp, totalfacets);
            DL_FOREACH(shapes,elt){
                if (elt->countinplate == 0)
                    continue;
                for (copy = 0; copy < elt->done; ++copy)
                {
                    if(elt->plate[copy]==plate){
                        stl_file *s=elt->stl;
                        stl_translate(s,0-(s->stats.max.x-s->stats.min.x)/2.0,0-(s->stats.max.y-s->stats.min.y)/2.0,0);
                        stl_rotate_z(s,-elt->rotangle[copy]);
                        stl_translate_rel(s,elt->x[copy],elt->y[copy], 0);
                        stl_write_binary_block(s,fp);
                        if (elt->count > 1) { // Reset object state
                            stl_translate_rel(s,-elt->x[copy],-elt->y[copy], 0);
                            stl_rotate_z(s,elt->rotangle[copy]);
                        }
                    }
                }
            }
            
            fclose(fp);
            stl_file stl_in;
            int last_edges_fixed = 0;
            int iterations = 2;
            stl_open(&stl_in, stlfn);
            stl_check_facets_exact(&stl_in);
            stl_in.stats.facets_w_1_bad_edge = 
                (stl_in.stats.connected_facets_2_edge -
                stl_in.stats.connected_facets_3_edge);
            stl_in.stats.facets_w_2_bad_edge = 
                (stl_in.stats.connected_facets_1_edge -
                stl_in.stats.connected_facets_2_edge);
            stl_in.stats.facets_w_3_bad_edge = 
                (stl_in.stats.number_of_facets -
                stl_in.stats.connected_facets_1_edge);
            float tolerance = stl_in.stats.shortest_edge;
            float increment = stl_in.stats.bounding_diameter / 10000.0;
            if(stl_in.stats.connected_facets_3_edge < stl_in.stats.number_of_facets)
            {
              for(i = 0; i < iterations; i++)
                {
                  if(stl_in.stats.connected_facets_3_edge < 
                 stl_in.stats.number_of_facets)
                {
                  stl_check_facets_nearby(&stl_in, tolerance);
                  last_edges_fixed = stl_in.stats.edges_fixed;
                  tolerance += increment;
                }
              }
              stl_remove_unconnected_facets(&stl_in);
              stl_fill_holes(&stl_in);
            }
            stl_fix_normal_directions(&stl_in);
            stl_fix_normal_values(&stl_in);
            stl_verify_neighbors(&stl_in);
            stl_generate_shared_vertices(&stl_in);
            stl_write_binary(&stl_in, stlfn, stlfn);

            if (maxextraplates > 0)
            {
                int k;
                for (k = 0; k < maxextraplates; ++k) {
                    int platek = plate + k + 1;
                    sprintf(tmpfn,"plate%02d.png", platek);
                    imagefn[0]=0;
                    strcat(imagefn,outdir);
                    strcat(imagefn,tmpfn);
                    cvSaveImage(imagefn,img,0);
                    sprintf(tmpfn,"plate%02d.stl", platek);
                    stlfn[0]=0;
                    strcat(stlfn,outdir);
                    strcat(stlfn,tmpfn);
                    stl_write_binary(&stl_in, stlfn, stlfn);
                }
            }
            stl_close(&stl_in);
        }
            
        plate += 1 + maxextraplates;
    
        platecount=0;
    }    
    
    
    DL_FOREACH_SAFE(shapes,elt,tmp) {
      cvReleaseImage(&(elt->image)); 
      free(elt->stl);
      free(elt->x);
      free(elt->y);
      free(elt->rotangle);
      free(elt->plate);
      DL_DELETE(shapes,elt);
      free(elt);
    }
    cvReleaseImage(&itmp);
    cvReleaseImage(&img);
    cvReleaseImage(&rpatch);
    cvReleaseImage(&testfit);

    return EXIT_SUCCESS;
}

