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
#ifndef MTF_CORE_H
#define MTF_CORE_H

#include "include/component_labelling.h"
#include "include/gradient.h"
#include "include/block.h"
#include "include/rectangle.h"
#include "include/ellipse.h"
#include "include/edge_record.h"
#include "include/loess_fit.h"
#include "include/afft.h"
#include "include/mtf_profile_sample.h"
#include "include/bayer.h"

#include <map>
using std::map;

typedef vector<Block> block_vector;

// global constants for ESF-fourier MTF method
// TODO: these can be dynamic parameters, with some effort
const double max_dot = 28;
const int SAMPLES_PER_PIXEL = 32;
const size_t FFT_SIZE = int(16)*SAMPLES_PER_PIXEL;
const int NYQUIST_FREQ = FFT_SIZE/16;
const double max_edge_length = 200;

class Mtf_core {
  public:
    

    Mtf_core(const Component_labeller& in_cl, const Gradient& in_g, 
             const cv::Mat& in_img, const cv::Mat& in_bayer_img, std::string bayer_subset)
      : cl(in_cl), g(in_g), img(in_img), bayer_img(in_bayer_img), absolute_sfr(false),
        snap_to(false), snap_to_angle(0), sfr_smoothing(true),
        sliding(false), samples_per_edge(0), border_width(0) {

        bayer = Bayer::from_string(bayer_subset);
        printf("bayer subset is %d\n", bayer);
      
        for (Boundarylist::const_iterator it=cl.get_boundaries().begin(); it != cl.get_boundaries().end(); ++it) {
            valid_obj.push_back(it->first);
        }
        
        cv::Mat temp;
        in_img.convertTo(temp, CV_8U, 256.0/16384.0);
        cv::cvtColor(temp, od_img, CV_GRAY2RGB);
    }
    
    ~Mtf_core(void) {
        //cv::imwrite(string("detections.png"), od_img);
    }
    
    size_t num_objects(void) {
        return valid_obj.size();
    }
    
    void search_borders(const Point2d& cent, int label);
    bool extract_rectangle(const Point2d& cent, int label, Mrectangle& rect);
    double compute_mtf(const Point2d& in_cent, const map<int, scanline>& scanset, 
                       Edge_record& er,
                       double& poor, 
                       vector<double>& sfr, vector<double>& esf);
    
    vector<Block>& get_blocks(void) {
        // make a copy into an STL container if necessary
        if (detected_blocks.size() == 0) {
            for (map<int,Block>::const_iterator it=shared_blocks_map.begin();
                 it != shared_blocks_map.end(); ++it) {
                
                bool allzero = true;
                for (int k=0; k < 4 && allzero; k++) {
                    if (fabs(it->second.get_mtf50_value(k)) > 1e-6) {
                        allzero = false;
                    }
                }
                
                if (it->second.valid && !allzero) { 
                    detected_blocks.push_back(it->second);
                }
            }
        }
        return detected_blocks;
    }
    
    vector<Mtf_profile_sample>& get_samples(void) {
        return samples;
    }
    
    void set_absolute_sfr(bool val) {
        absolute_sfr = val;
    }
    
    void set_sfr_smoothing(bool val) {
        sfr_smoothing = val;
    }
    
    void set_sliding(bool val) {
        sliding = val;
    }
    
    void set_samples_per_edge(int s) {
        samples_per_edge = s;
    }
    
    void set_snap_angle(double angle) {
        snap_to = true;
        snap_to_angle = angle;
    }
    
    void set_border(int in_border_width) {
        border_width = in_border_width;
    }
    
    const Component_labeller& cl;
    const Gradient&           g;
    const cv::Mat&            img;
    const cv::Mat&            bayer_img;
    Bayer::bayer_t bayer;
    
    AFFT<512> afft; // FFT_SIZE = 512 ??
    vector<int> valid_obj;
    
    vector<Block> detected_blocks;  
    map<int, Block> shared_blocks_map;
    vector<Point2d> solid_ellipses;
    vector<Ellipse_detector> ellipses;
    
    vector<Mtf_profile_sample> samples;
    
    cv::Mat od_img;
  private:
    bool absolute_sfr;
    bool snap_to;
    double snap_to_angle;
    bool sfr_smoothing;
    bool sliding;
    int samples_per_edge;
    int border_width;
    
    void process_with_sliding_window(Mrectangle& rrect);
  
    void sample_at_angle(double ea, vector<Ordered_point>& local_ordered, 
        const map<int, scanline>& scanset, const Point2d& cent,
        double& edge_length) {

        double max_along_edge = -1e50;
        double min_along_edge = 1e50;
        
        if (snap_to) {
            
            double max_dot_angle = snap_to_angle;
            double max_dot = 0;
            
            double angles[4] = {snap_to_angle, -snap_to_angle, M_PI/2 - snap_to_angle, snap_to_angle - M_PI/2};
            
            for (int k=0; k < 4; k++) {
            
                double sa = angles[k];
                
                double dot = cos(ea)*cos(sa) + sin(ea)*sin(sa);
                if (dot > max_dot) {
                    max_dot = dot;
                    max_dot_angle = sa;
                }
            }
            
            ea = max_dot_angle;
        }

        Point2d mean_grad(cos(ea), sin(ea));
        Point2d edge_direction(-sin(ea), cos(ea));

        if (bayer == Bayer::NONE) {
            for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
                int y = it->first;
                if (y < border_width || y > img.rows-1-border_width) continue;
                
                for (int x=it->second.start; x <= it->second.end; ++x) {
                    
                    if (x < border_width || x > img.cols-1-border_width) continue;
                    
                    Point2d d((x) - cent.x, (y) - cent.y);
                    double dot = d.ddot(mean_grad); 
                    double dist_along_edge = d.ddot(edge_direction);
                    if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                        local_ordered.push_back(Ordered_point(dot, img.at<uint16_t>(y,x) ));
                        max_along_edge = max(max_along_edge, dist_along_edge);
                        min_along_edge = min(min_along_edge, dist_along_edge);
                    }
                }
            }
        } else {
            for (map<int, scanline>::const_iterator it=scanset.begin(); it != scanset.end(); ++it) {
                int y = it->first;
                int rowcode = (y & 1) << 1;
                if (y < border_width || y > img.rows-1-border_width) continue;
                
                for (int x=it->second.start; x <= it->second.end; ++x) {
                    int code = rowcode | (x & 1);
                    
                    if (x < border_width || x > img.cols-1-border_width) continue;

                    // skip the appropriate sites if we are operating only on a subset
                    if (bayer == Bayer::RED && code != 0) {
                        continue;
                    } 
                    if (bayer == Bayer::BLUE && code != 3) {
                        continue;
                    } 
                    if (bayer == Bayer::GREEN && (code == 0 || code == 3)) {
                        continue;
                    }

                    Point2d d((x) - cent.x, (y) - cent.y);
                    double dot = d.ddot(mean_grad); 
                    double dist_along_edge = d.ddot(edge_direction);
                    if (fabs(dot) < max_dot && fabs(dist_along_edge) < max_edge_length) {
                        local_ordered.push_back(Ordered_point(dot, bayer_img.at<uint16_t>(y,x) ));
                        max_along_edge = max(max_along_edge, dist_along_edge);
                        min_along_edge = min(min_along_edge, dist_along_edge);
                    }
                }
            }
        }
    
        
        edge_length = max_along_edge - min_along_edge;
        //printf("edge length = %lf\n", edge_length);
    }

};

#endif
