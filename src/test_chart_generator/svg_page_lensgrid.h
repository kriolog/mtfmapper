/*
Copyright 2011 Frans van den Bergh. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Frans van den Bergh ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Frans van den Bergh OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of the Council for Scientific and Industrial Research (CSIR).
*/
#ifndef SVG_PAGE_LENSGRID_H
#define SVG_PAGE_LENSGRID_H

#include "svg_page.h"
#include "include/mtf50_edge_quality_rating.h"
#include <vector>
#include <cmath>
using std::vector;
using std::pair;
using std::make_pair;

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

class Svg_page_lensgrid : public Svg_page {
  public:
    Svg_page_lensgrid(const string& page_spec, const string& fname) 
    : Svg_page(page_spec, fname), centre(0.5, 0.5*sqrt(2.0)),
      off(0), scale(10.0/100) {
    
    }
    
    void render(void) {
        grid(0.04, 17, 11); 
    }
    
  protected:
    virtual iPoint project(double x, double y) {
        
        x = floor(x*width);
        y = floor(y*width);
        
        return iPoint(int(x), int(y));
    }

    double quant_angle(double x) {
        double quad1 = fabs(fmod(x, M_PI/2.0));
        if (quad1 > M_PI/4.0) {
            quad1 = M_PI/2.0 - quad1;
        }
        quad1 = quad1 / M_PI * 180;
        return quad1;
    }
    
    void triangle(dPoint a, dPoint b, dPoint c) {
       
        iPoint p = project(a.x, a.y);
        fprintf(fout, "  <polygon points=\"%d,%d ", p.x, p.y);
        p = project(b.x, b.y);
        fprintf(fout, "%d,%d ", p.x, p.y);
        p = project(c.x, c.y);
        fprintf(fout, "%d,%d\" style=\"%s\"/>\n", p.x, p.y, style.c_str());
    }
    
    void place_trapezoid(double xpos, double ypos, double width, double angle_thresh=1.5) {
        vector<iPoint> coords(4);
        
        static double direction = 1;
        
        //direction *= -1;
        
        double ewidth = 0;
        
        bool done = false;
        int tries = 0;
        int last_change = 0;
        do {
            dPoint delta = dPoint(xpos, ypos) - centre;
            
            //#define SQUARES
            
            double norm = sqrt(delta.x*delta.x + delta.y*delta.y);
            delta.x /= norm;
            delta.y /= norm;
            
            dPoint omx = dPoint(xpos + delta.x*width*sqrt(0.5), ypos + delta.y*width*sqrt(0.5));
            dPoint tang = dPoint(-delta.y, delta.x);
            
            dPoint tl(omx.x - 0.5*width*tang.x, omx.y - 0.5*width*tang.y);
            iPoint p = project(tl.x, tl.y);
            coords[0] = p;
            
            dPoint rad = tl - centre;
            norm = sqrt(rad.x*rad.x + rad.y*rad.y);
            rad.x /= norm;
            rad.y /= norm;
            double radang1 = atan2(fabs(rad.y), fabs(rad.x)) * 180.0/M_PI;
            if (radang1 > 45) radang1 = 90 - radang1;
            #ifdef SQUARES
            rad = delta;
            #endif
            
            dPoint bl(tl.x - rad.x*width + width*ewidth*tang.x, tl.y - rad.y*width + width*ewidth*tang.y);
            p = project(bl.x, bl.y);
            coords[1] = p;
            
            dPoint tr(omx.x + 0.5*width*tang.x, omx.y + 0.5*width*tang.y);
            
            rad = tr - centre;
            norm = sqrt(rad.x*rad.x + rad.y*rad.y);
            rad.x /= norm;
            rad.y /= norm;
            double radang2 = atan2(fabs(rad.y), fabs(rad.x)) * 180.0/M_PI;
            if (radang2 > 45) radang2 = 90 - radang2;
            #ifdef SQUARES
            rad = delta;
            #endif
            
            dPoint br(tr.x - rad.x*width - width*ewidth*tang.x, tr.y - rad.y*width - width*ewidth*tang.y);
            p = project(br.x, br.y);
            coords[2] = p;
            
            p = project(tr.x, tr.y);
            coords[3] = p;
            
            double angles[4];
            for (int i=0; i < 4; i++) {
                int next = (i + 1) % 4;
                double angle = atan2(fabs(double(coords[i].y - coords[next].y)), fabs(double(coords[i].x - coords[next].x)));
                angle *= 180.0/M_PI;
                if (angle > 45) {
                    angle = 90 - angle;
                }
                // identify which edge this is
                int j = -1;
                double l2 = sqrt(SQR(coords[i].x - coords[next].x) + SQR(coords[i].y - coords[next].y));
                double dot = (coords[i].x - coords[next].x)*delta.x + (coords[i].y - coords[next].y)*delta.y;
                dot /= l2;
                if (fabs(dot) < 0.707) {
                    j = (dot < 0) ? 0 : 1;
                } else {
                    j = (dot < 0) ? 2 : 3;
                }
                if (j < 0) {
                    printf("error. not a valid orientation\n");
                    exit(1);
                }
                angles[j] = angle;
            }
            
            if (tries % 2 == 0) {
                if (angles[0] < angle_thresh || angles[0] > 43.5 || fabs(angles[0] - 26.565) < 2 ||
                    angles[1] < angle_thresh || angles[1] > 43.5 || fabs(angles[1] - 26.565) < 2 ) { 
                    xpos -= 0.001*delta.y*direction;
                    ypos += 0.001*delta.x*direction;
                    last_change = tries;
                }
                
            } else {
                if (angles[2] < 2 || angles[2] > 43 || angles[3] < 2 || angles[3] > 43 ||
                    fabs(angles[2] - 26.565) < 2 || fabs(angles[3] - 26.565) < 2) {
                    ewidth += 0.01;
                    last_change = tries;
                }
            }
            tries++;
            if (tries > last_change+3) done = true;
            
            // estimate radial error
            printf("%lf %lf\n", 
                std::min(fabs(radang1-angles[2]), fabs(radang2-angles[2])),
                std::min(fabs(radang1-angles[3]), fabs(radang2-angles[3]))
            );
            
        } while (!done);
        
        
        fprintf(fout, "  <polygon points=\"%d,%d ", coords[0].x, coords[0].y);
        fprintf(stderr, "4\n%lf %lf\n", coords[0].y*scale + off, coords[0].x*scale + off);
        fprintf(fout, "%d,%d ", coords[1].x, coords[1].y);
        fprintf(stderr, "%lf %lf\n", coords[1].y*scale + off, coords[1].x*scale + off);
        fprintf(fout, "%d,%d ", coords[2].x, coords[2].y);
        fprintf(stderr, "%lf %lf\n", coords[2].y*scale + off, coords[2].x*scale + off);
        fprintf(fout, "%d,%d\" style=\"%s\"/>\n", coords[3].x, coords[3].y, style.c_str());
        fprintf(stderr, "%lf %lf\n", coords[3].y*scale + off, coords[3].x*scale + off);
    }

    void grid(double swidth, size_t nrows, size_t ncols) {
        
        const double phi = 1.9/180.0*M_PI;
        for (int row=0; row < (int)nrows; row++) {
            for (int col=0; col < (int)ncols; col++) {
                if (row == (int)nrows/2 && col == (int)ncols/2) continue;
                double sx = (col-int(ncols)/2)/double(ncols)*0.93;
                double sy = sqrt(2.0)*(row-int(nrows)/2)/double(nrows)*0.93;
                double rx = sx*cos(phi) - sy*sin(phi);
                double ry = sx*sin(phi) + sy*cos(phi);
                place_trapezoid(
                    0.5 + rx,
                    sqrt(0.5) + ry,
                    swidth,
                    (col == (int)ncols/2) || (row == (int)nrows/2) ? 3 : 1.5
                );
            }
        }
        
        
        
        iPoint p = project(0.985, sqrt(2)*0.985);
        fprintf(stderr, "3\n%lf %lf\n", p.y*scale + off, p.x*scale + off);
        p = project(0.987, sqrt(2)*0.985);
        fprintf(stderr, "%lf %lf\n", p.y*scale + off, p.x*scale + off);
        p = project(0.987, sqrt(2)*0.9865);
        fprintf(stderr, "%lf %lf\n", p.y*scale + off, p.x*scale + off);
        
        
        triangle(
          dPoint(centre.x, centre.y),
          dPoint(centre.x - 0.04, centre.y + 0.04),
          dPoint(centre.x - 0.04, centre.y - 0.04)
        );
        
        
        triangle(
          dPoint(centre.x, centre.y),
          dPoint(centre.x + 0.04, centre.y + 0.04),
          dPoint(centre.x + 0.04, centre.y - 0.04)
        );
        
        triangle(
          dPoint(1, 0),
          dPoint(1, 0.05),
          dPoint(0.5,  0)
        );
        
        triangle(
          dPoint(1-0.04, 0.05 - 0.04*0.1),
          dPoint(1, 0.05),
          dPoint(1,  0.1)
        );
        
        triangle(
          dPoint(0, sqrt(2.0)),
          dPoint(0, sqrt(2.0)-0.05),
          dPoint(0.5,  sqrt(2.0))
        );
        
        triangle(
          dPoint(0.04, sqrt(2.0) - 0.05 + 0.04*0.1),
          dPoint(0,  sqrt(2.0) - 0.05),
          dPoint(0,  sqrt(2.0) - 0.1)
        );
        
        triangle(
          dPoint(0.05, 0),
          dPoint(0, 0.5),
          dPoint(0, 0.0)
        );
        
        triangle(
          dPoint(0.05, 0),
          dPoint(0.1, 0),
          dPoint(0.05 - 0.04*0.1, 0.04)
        );
        
        triangle(
          dPoint(1 - 0.05, sqrt(2.0)),
          dPoint(1, sqrt(2.0) - 0.5),
          dPoint(1, sqrt(2.0))
        );
        
        triangle(
          dPoint(1 - 0.05, sqrt(2.0)),
          dPoint(1 - 0.1, sqrt(2.0)),
          dPoint(1 - 0.05 + 0.04*0.1, sqrt(2.0) - 0.04)
        );
    }
    
    dPoint centre;
    double off;
    double scale;
};

#endif
