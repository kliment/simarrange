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
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "utlist.h"
#include <dirent.h>
#include <argtable2.h>

typedef struct img_list{
    IplImage *image;
    long area;
    int count;
    int x;
    int y;
    int rotangle;
    int plate;
    char filename[350];
    stl_file *stl;
    struct img_list *prev,*next;
} img_list;


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
    stl_file *s=(stl_file *)malloc(sizeof(stl_file));
    memset(s,0,sizeof(stl_file));
    stl_open(s,filename);
    stl_check_facets_exact(s);
    s->stats.facets_w_1_bad_edge = (s->stats.connected_facets_2_edge - s->stats.connected_facets_3_edge);
    s->stats.facets_w_2_bad_edge = (s->stats.connected_facets_1_edge - s->stats.connected_facets_2_edge);
    s->stats.facets_w_3_bad_edge = (s->stats.number_of_facets - s->stats.connected_facets_1_edge);
    stl_fix_normal_directions(s);
    stl_fix_normal_values(s);
    stl_fill_holes(s);
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
    e->plate=-1;
    e->stl=s;
    e->count=count;
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
        if(elt->plate==-1)
            ++c;
    }
    return c;
}

int main(int argc, char** argv){
    
    int w=200,h=200;
    int spacing=1;
    int rotstep=10;
    int posstep=5;
    int c;
    char defdir[]="test";
    char *indir=NULL;
    struct arg_int  *aw  = arg_int0("x","width",NULL,              "plate width in mm (default is 200)");
    struct arg_int  *ah  = arg_int0("y","height",NULL,              "plate height in mm (default is 200)");
    struct arg_int  *as  = arg_int0("s","spacing",NULL,              "spacing between parts in mm (default is 1)");
    struct arg_int  *ar  = arg_int0("r","rotstep",NULL,              "rotation step when searching (default 10 degrees)");
    struct arg_int  *ap  = arg_int0("p","posstep",NULL,              "positional step when searching (default 5mm)");
    struct arg_str  *aindir = arg_str0("i","inputdir",NULL,  "input directory (default \"test\")");
    struct arg_str  *ainfile = arg_strn("f","inputfile",NULL,0,argc+2,  "input file (any number allowed)");
    struct arg_end  *end      = arg_end(20);
    void* argtable[] = {aw,ah,as,ar,ap,aindir,ainfile,end};
    
    int nerrors;
    nerrors = arg_parse(argc,argv,argtable);
    if (nerrors > 0)
        {
        /* Display the error details contained in the arg_end struct.*/
        arg_print_errors(stdout,end,"");
        return EXIT_FAILURE;
        }
    
    if(aindir->count){
        indir=(char *)(aindir->sval[0]);
    }
    if(aw->count){
        w=aw->ival[0];
    }
    if(ah->count){
        h=ah->ival[0];
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
    
    if(!indir){
        indir=defdir;
    }
    
    
    
    img_list *shapes=NULL;
    img_list *curplate=NULL;
    int platecount=0;
    if(ainfile->count){
        int i;
        for(i=0;i<ainfile->count;i++){
            add_stl((char *)(ainfile->sval[i]),1, w, h, &shapes);
        }
    }
    if (aindir->count || !ainfile->count){
        DIR *dir;
        struct dirent *ent;
        if ((dir = opendir (indir)) != NULL) {
            while ((ent = readdir (dir)) != NULL) {
                int i;
                char d[350];
                strcpy(d,ent->d_name);
                for(i = 0; d[i]; i++)
                    d[i] = tolower(d[i]);
                if(strstr(d,".stl")!=0){
                    char f[350];
                    f[0]=0;
                    strcat(f,indir);
                    strcat(f+strlen(f),"/");
                    strcat(f+strlen(f),ent->d_name);
                    
                    add_stl(f, 1, w, h, &shapes);
                }
                    //printf ("%s\n", ent->d_name);
            }
            closedir (dir);
        } else {
            printf("Input directory not found\n");
            return EXIT_FAILURE;
        }
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
    
            
    while(dl_count(shapes)){
        printf("Generating plate %d\n",plate);
        cvZero(img);
        cvLine(img, cvPoint(0,0), cvPoint(w-1,0), cvScalarAll(127), 1, 8, 0);
        cvLine(img, cvPoint(w-1,0), cvPoint(w-1,h-1), cvScalarAll(127), 1, 8, 0);
        cvLine(img, cvPoint(w-1,h-1), cvPoint(0,h-1), cvScalarAll(127), 1, 8, 0);
        cvLine(img, cvPoint(0,h-1), cvPoint(0,0), cvScalarAll(127), 1, 8, 0);
        cvZero(rpatch);
        cvZero(testfit);
        cvZero(itmp);
        
        DL_FOREACH(shapes,elt) {
            if(elt->plate>-1)
                continue;
            //printf("File: %s\n",elt->filename);
            cvCopy(img, testfit, NULL);
            int xpos=1, ypos=1, rotangle=0, minxpos=w-1, minypos=h-1, minrotangle=0;
            cvDilate(elt->image,itmp,NULL,spacing);
            for(rotangle=0;rotangle<360;rotangle+=rotstep){
                cvWarpAffine(itmp,rpatch,cv2DRotationMatrix(center, rotangle, 1.0, rot),CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS, cvScalarAll(0) );
                for(ypos=5;ypos<minypos;ypos+=posstep){
                    for(xpos=5;xpos<minxpos;xpos+=posstep){
                        cvSetImageROI(rpatch, cvRect(w-xpos,h-ypos,w,h));
                        cvAnd(rpatch, img, testfit, NULL);
                        if(!cvCountNonZero(testfit)){
                            int prec=cvCountNonZero(img);
                            cvAdd(rpatch, img, testfit, NULL);
                            if(prec!=cvCountNonZero(testfit)){
                                minxpos=xpos;
                                minypos=ypos;
                                minrotangle=rotangle;
                            }
                        }
                        cvResetImageROI(rpatch);
                        
                    }
                }
            }
            
            if(minxpos<w-1 || minypos<h-1){
                printf("File: %s minx: %d, miny: %d, minrot: %d\n",elt->filename, minxpos, minypos, minrotangle);
                cvWarpAffine(elt->image,rpatch,cv2DRotationMatrix(center, minrotangle, 1.0, rot),CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS, cvScalarAll(0) );
                cvSetImageROI(rpatch, cvRect(w-minxpos,h-minypos,w,h));
                cvAdd(rpatch,img,testfit,NULL);
                cvCopy(testfit,img,NULL);
                cvResetImageROI(rpatch);
                elt->x=minxpos;
                elt->y=minypos;
                elt->rotangle=minrotangle;
                elt->plate=plate;
                platecount++;
            }else{
                printf("SKIP: %s skipped for this plate\n",elt->filename);
            }
            
        
            
        }
        
        if(!platecount){
            printf("The files skipped in the last stage do not fit on plate in any tested orientation! They might be too large for the print area.\n");
            return EXIT_FAILURE;
        }
        char fn[20];
        sprintf(fn,"plate%02d.png",plate);
        cvFlip(img,NULL,0);
        cvSaveImage(fn,img,0);
        sprintf(fn,"plate%02d.stl",plate);
        FILE *fp;
        int i;
        fp = fopen(fn, "w");
        if(fp == NULL){
          printf("Could not open output file %s\n",fn);
          return EXIT_FAILURE;
        }
        fprintf(fp, "%s", fn);
        for(i = strlen(fn); i < LABEL_SIZE; i++) putc(0, fp);
        fseek(fp, LABEL_SIZE, SEEK_SET);
        
        int totalfacets=0;
        DL_FOREACH(shapes,elt){
            if(elt->plate==plate){
                totalfacets+=elt->stl->stats.number_of_facets;
            }
        }
        stl_put_little_int(fp, totalfacets);
        DL_FOREACH(shapes,elt){
            if(elt->plate==plate){
                stl_file *s=elt->stl;
                stl_translate(s,0-(s->stats.max.x-s->stats.min.x)/2.0,0-(s->stats.max.y-s->stats.min.y)/2.0,0);
                stl_rotate_z(s,-elt->rotangle);
                stl_translate_rel(s,elt->x,elt->y,0 );
                stl_write_binary_block(s,fp);
            }
        }
        
        fclose(fp);
        plate+=1;
    
        platecount=0;
    }    
    
    
    DL_FOREACH_SAFE(shapes,elt,tmp) {
      cvReleaseImage(&(elt->image)); 
      free(elt->stl);
      DL_DELETE(shapes,elt);
      free(elt);
    }
    cvReleaseImage(&itmp);
    cvReleaseImage(&img);
    cvReleaseImage(&rpatch);
    cvReleaseImage(&testfit);
    
}

