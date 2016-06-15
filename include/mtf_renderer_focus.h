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
#ifndef MTF_RENDERER_FOCUS_H
#define MTF_RENDERER_FOCUS_H

#include "mtf_renderer.h"
#include "common_types.h"
#include "loess_fit.h"
#include "mtf_profile_sample.h"

#include <stdlib.h>
#include <vector>
using std::vector;
using std::pair;
using std::make_pair;

#include "focus_surface.h"
#include "distance_scale.h"
#include "mtf50_edge_quality_rating.h"
#include "bayer.h"
#include "ellipse.h"

class Mtf_renderer_focus : public Mtf_renderer {
  public:
    Mtf_renderer_focus(Distance_scale& distance_scale,
        const std::string& wdir, const std::string& prof_fname, 
        const cv::Mat& img, 
        bool lpmm_mode=false, double pixel_size=1.0) 
      :  zero(distance_scale.zero), 
         transverse(distance_scale.transverse), 
         longitudinal(distance_scale.longitudinal),
         wdir(wdir), prname(prof_fname),
         img(img), 
         lpmm_mode(lpmm_mode), pixel_size(pixel_size),
         distance_scale(distance_scale) {
         
    }
    
    void render(const vector<Block>&) {
        printf("Fatal error. This function should not be used. Aborting\n");
        exit(1);
        return;
    }
    
    void render(const vector<Mtf_profile_sample>& samples, Bayer::bayer_t bayer = Bayer::NONE, vector<Ellipse_detector>* ellipses = NULL) {
        Point2d centroid(0,0);
        
        if (samples.size() < 10) { // probably not a valid image for profiles
            return;
        }
        
        if (!distance_scale.fiducials_found) {
            printf("Error: No valid fiducials found, refusing to generate Focus chart output.\nAre you using the right chart type?\n");
            return;
        }
        
        double mean_y = 0;
        vector<Sample> data;
        double min_long = 1e50;
        double max_long = -1e50;
        double min_pix_long = 1e50;
        double max_pix_long = -1e50;
        for (size_t i=0; i < samples.size(); i++) {
            if (samples[i].quality > very_poor_quality) {
                Point2d d = samples[i].p - zero;
                Point2d coord(
                    d.x*longitudinal.x + d.y*longitudinal.y,
                    d.x*transverse.x + d.y*transverse.y
                );
                
                Point2d wc = distance_scale.estimate_world_coords(samples[i].p.x, samples[i].p.y);
                
                if (fabs(wc.y) < 20 && fabs(wc.x) < 180) { 
                
                    mean_y += coord.y;
                    
                    double depth = 0;
                    distance_scale.estimate_depth_img_coords(samples[i].p.x, samples[i].p.y, depth);
                    coord.x = depth;
                    
                    // coord is now projected into the surface
                    data.push_back(Sample(coord.x, samples[i].mtf, 1.0, 1.0));
                    min_long = min(coord.x, min_long);
                    max_long = max(coord.x, max_long);
                    min_pix_long = min(min_pix_long, coord.x);
                    max_pix_long = max(max_pix_long, coord.x);
                }
            }
        }
        mean_y /= double(data.size());
        
        const int sh = 4;
        const double sgw[] = {-21/231.0, 14/231.0, 39/231.0, 54/231.0, 59/231.0, 54/231.0, 39/231.0, 14/231.0, -21/231.0};
        sort(data.begin(), data.end());
        
        // just pretend our samples are equally spaced
        vector<Sample> ndata;
        for (size_t i=sh; i < data.size() - sh; i++) {
            double val = 0;
            for (int w=-sh; w <= sh; w++) {
                val += sgw[w+sh] * data[i+w].y;    
            }
            
            
            if (fabs(val - data[i].y)/val < 0.04) {
                ndata.push_back(data[i]);
            }
            
        }
        printf("dropped %lu samples, %lu remain\n", data.size() - ndata.size(), ndata.size());
        data = ndata;
        
        Ratpoly_fit cf(data, 4, 2);
        VectorXd sol = rpfit(cf);
        
        // perform a few steps of IRLS
        double prev_err = 1e50;
        int dccount;
        for (int iter=0; iter < 50; iter++) {
            double errsum = 0;
            double wsum = 0;
            dccount = 0;
            for (size_t k=0; k < data.size(); k++) {
                double y = cf.rpeval(sol, data[k].x*cf.xsf)/cf.ysf;
                double e = fabs(y - data[k].y);
                errsum += e * data[k].yweight;
                wsum += data[k].yweight;
                data[k].yweight = 1.0; 
                if (e/y > 0.05) { // kill really poor outliers
                    data[k].yweight = 0;
                    dccount++;
                }
                
            }
            errsum /= wsum;
            printf("iter %d err: %lg dc=%lg\n", iter, errsum, dccount/double(data.size()));
            sol = rpfit(cf);
            if (iter > 0 && (prev_err - errsum)/prev_err < 0.0001) {
                printf("bailing out at iter %d\n", iter);
                break;
            }
            prev_err = errsum;
        }
        double errsum = 0;
        double wsum = 0;
        for (size_t k=0; k < data.size(); k++) {
            double y = cf.rpeval(sol, data[k].x*cf.xsf)/cf.ysf;
            double e = fabs(y - data[k].y);
            errsum += e * data[k].yweight;
            wsum += data[k].yweight;
        }
        errsum /= wsum;
        printf("final model fit error (weighted): %lg\n", errsum);
        
        // do a quick bootstrap to estimate some bounds on the focus peak
        vector<double> mc_peaks;
        const int n_mc = 30;
        for (int i=0; i < n_mc; i++) {
            vector<Sample> mc_data;
            for (size_t j=i; j < data.size(); j += 10) {
                mc_data.push_back(data[j]);
            }
            Ratpoly_fit mc_cf(mc_data, 4, 2);
            VectorXd mc_sol = rpfit(mc_cf);
            
            double mc_peak = mc_cf.peak(mc_sol);
            mc_peaks.push_back(mc_peak);
        }
        sort(mc_peaks.begin(), mc_peaks.end());
        double mc_p5 = mc_peaks[lrint(n_mc*0.05)];
        double mc_p95 = mc_peaks[lrint(n_mc*0.95)];
        
        string rp_name = wdir + ((bayer != Bayer::NONE) ? Bayer::to_string(bayer) + "_" : "") + "profile_points.txt";
        FILE* fraw = fopen(rp_name.c_str(), "wt");
        for (size_t i=0; i < data.size(); i++) {
            fprintf(fraw, "%lf %lf\n", data[i].x, data[i].y);
        }
        fclose(fraw);
        
        // now we can plot the reconstruction?
        double peak_y = 0;
        string cp_name = wdir + ((bayer != Bayer::NONE) ? Bayer::to_string(bayer) + "_" : "") + "profile_curve.txt";
        FILE* fout = fopen(cp_name.c_str(), "wt");
        for (double x=min_long; x < max_long; x += 1) {
            double y = cf.rpeval(sol, x*cf.xsf)/cf.ysf;
            fprintf(fout, "%lf %lf\n", x, y);
            
            if (y > peak_y) {
                peak_y = y;
            }
        }
        fclose(fout);
        
        double focus_peak = cf.peak(sol); // works on pre-stretched depth
        printf("focus_plane %lg\n", focus_peak);
        
        cv::Mat channel(img.rows, img.cols, CV_8UC1);
        double imin;
        double imax;
        cv::minMaxLoc(img, &imin, &imax);
        img.convertTo(channel, CV_8U, 255.0/(imax - imin), 0);
        
        vector<cv::Mat> channels;
        channels.push_back(channel);
        channels.push_back(channel);
        channels.push_back(channel);
        cv::Mat merged;
        merge(channels, merged);
        int initial_rows = merged.rows;
        merged.resize(merged.rows + 100);
        
        if (ellipses) {
            for (auto e: *ellipses) {
                if (!e.valid) continue;
                for (double theta=0; theta < 2*M_PI; theta += M_PI/720.0) {
                    double synth_x = e.major_axis * cos(theta);
                    double synth_y = e.minor_axis * sin(theta);
                    double rot_x = cos(e.angle)*synth_x - sin(e.angle)*synth_y + e.centroid_x;
                    double rot_y = sin(e.angle)*synth_x + cos(e.angle)*synth_y + e.centroid_y;

                    // clip to image size, just in case
                    rot_x = max(rot_x, 0.0);
                    rot_x = min(rot_x, (double)(merged.cols-1));
                    rot_y = max(rot_y, 0.0);
                    rot_y = min(rot_y, (double)(merged.rows-1));

                    cv::Vec3b& color = merged.at<cv::Vec3b>(lrint(rot_y), lrint(rot_x));
                    color[0] = 255;
                    color[1] = 255;
                    color[2] = 0;
                }
            }
        }
        
        vector<Point2d> curve;
        min_pix_long = -merged.cols/2;
        max_pix_long = merged.cols/2;
        double mtf_peak_value = 0;
        double peak_wx = 0;
        for (double x=min_pix_long; x < max_pix_long; x += 1) {
            double px = longitudinal.x * x  + longitudinal.y * mean_y + zero.x;
            double py = transverse.x * x + transverse.y *  mean_y + zero.y;
            
            double depth = 0;
            distance_scale.estimate_depth_img_coords(px, py, depth);
            
            if (depth > min_long && depth < max_long) {
                double mtf = cf.rpeval(sol, depth*cf.xsf)/cf.ysf;
                
                double world_y = (mtf / peak_y) * 130 - 130;
                    
                Point2d wc = distance_scale.estimate_world_coords(px, py);
                
                if (mtf > mtf_peak_value) {
                    mtf_peak_value = mtf;
                    peak_wx = wc.x;
                }
                    
                Point2d proj_p = distance_scale.world_to_image(wc.x, world_y);
                curve.push_back(proj_p);
            }
        }
        draw_curve(merged, curve, cv::Scalar(128, 128, 128), 6);
        draw_curve(merged, curve, cv::Scalar(40, 90, 40), 3, cv::Scalar(40, 255, 40));
        
        curve.clear();
        for (double wy=-130; wy < 130; wy += 1) {
            Point2d proj_p = distance_scale.world_to_image(peak_wx, wy);
            curve.push_back(proj_p);
        }
        draw_curve(merged, curve, cv::Scalar(255, 30, 30), 3);
        
        
        cv::Scalar axisdark = cv::Scalar(30, 140, 190);
        cv::Scalar axislight = cv::Scalar(40, 187, 255);
        
        curve.clear();
        for (double ystep=-135; ystep <= 135; ystep += 2) {
            curve.push_back(distance_scale.world_to_image(-180, ystep));
        }
        draw_curve(merged, curve, axisdark, 2, axisdark);
        curve.clear();
        for (double xstep=-180; xstep <= 180; xstep += 2) {
            curve.push_back(distance_scale.world_to_image(xstep, -135));
        }
        draw_curve(merged, curve, axisdark, 2, axislight);
        curve.clear();
        curve.push_back(distance_scale.world_to_image(-180, -135));
        curve.push_back(distance_scale.world_to_image(-180, -135, -100));
        draw_curve(merged, curve, axisdark, 2, axisdark);
        
        int font = cv::FONT_HERSHEY_DUPLEX; 
        char tbuffer[1024];
        
        sprintf(tbuffer, "Y");
        Point2d textpos = distance_scale.world_to_image(-178, 135);
        cv::putText(merged, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2.5, CV_AA);
        cv::putText(merged, tbuffer, textpos, font, 1, axisdark, 1, CV_AA);
        
        sprintf(tbuffer, "X");
        textpos = distance_scale.world_to_image(180, -133);
        cv::putText(merged, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2.5, CV_AA);
        cv::putText(merged, tbuffer, textpos, font, 1, axislight, 1, CV_AA);
        
        sprintf(tbuffer, "Z");
        textpos = distance_scale.world_to_image(-180, -132, -100);
        cv::putText(merged, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2.5, CV_AA);
        cv::putText(merged, tbuffer, textpos, font, 1, axisdark, 1, CV_AA);
        
        cv::Scalar resultcolour(255, 255, 30);
        sprintf(tbuffer, "%.3lf (c/p)", mtf_peak_value);
        int tbaseline=0;
        cv::Size tsize = cv::getTextSize(tbuffer, font, 1.0, 2.5, &tbaseline);
        textpos = distance_scale.world_to_image(peak_wx, -20);
        textpos -= Point2d(tsize.width/2.0, 0);
        cv::putText(merged, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2.5, CV_AA);
        cv::putText(merged, tbuffer, textpos, font, 1, resultcolour, 1, CV_AA);
        double prev_line_height = tsize.height;
        
        sprintf(tbuffer, "%.1lf mm", focus_peak);
        tsize = cv::getTextSize(tbuffer, font, 1.0, 2.5, &tbaseline);
        textpos = distance_scale.world_to_image(peak_wx, -20);
        textpos -= Point2d(tsize.width/2.0, -prev_line_height*1.5);
        cv::putText(merged, tbuffer, textpos, font, 1, CV_RGB(20, 20, 20), 2.5, CV_AA);
        cv::putText(merged, tbuffer, textpos, font, 1, resultcolour, 1, CV_AA);
        
        
        sprintf(tbuffer, "Focus peak at depth %.1lf mm [%.1lf,%.1lf] relative to chart origin."
            " Estimated chart distance=%.2lf mm. Bundle adjustment RMSE=%.3lf pixels", 
            focus_peak, mc_p5, mc_p95, distance_scale.centre_depth, distance_scale.bundle_rmse
        );
        
        rectangle(merged, Point2d(0, initial_rows), Point2d(merged.cols, merged.rows), cv::Scalar::all(255), CV_FILLED);
        cv::putText(merged, tbuffer, Point2d(50, initial_rows + (merged.rows-initial_rows)/2), font, 1, cv::Scalar::all(0), 1, CV_AA);
        
        imwrite(wdir + prname, merged);
        
    }

  private:
  
    void draw_curve(cv::Mat& image, const vector<Point2d>& data, cv::Scalar col, double width, cv::Scalar col2=cv::Scalar(0, 0, 0)) {
        int prevx = 0;
        int prevy = 0;
        
        double total_l = 0;
        for (size_t i=1; i < data.size(); i++) {
            double dx = data[i].x - data[i-1].x;
            double dy = data[i].y - data[i-1].y;
            total_l += sqrt(dx*dx + dy*dy);
        }
        
        bool shade = col2[0] != 0 || col2[1] != 0 || col2[1] != 0;
        cv::Scalar blended_col = col;
        
        double running_l = 0;
        for (size_t i=0; i < data.size(); i++) {
            if (i > 1) {
                double dx = data[i].x - data[i-1].x;
                double dy = data[i].y - data[i-1].y;
                running_l += sqrt(dx*dx + dy*dy);
            }
            
            double progress = running_l / total_l;
            
            if (shade) {
                for (int k=0; k < 3; k++) {
                    blended_col[k] = col[k] + (col2[k] - col[k])*progress;
                }
            }
            
            int ix = lrint(data[i].x);
            int iy = lrint(data[i].y);
            if (ix >= 0 && ix < image.cols &&
                iy >= 0 && iy < image.rows && i > 0) {
                
                cv::line(image, Point2d(prevx, prevy), Point2d(ix, iy), blended_col, width, CV_AA);
            }
            
            prevx = ix;
            prevy = iy;
        }
    }
    
    VectorXd rpfit(Ratpoly_fit& cf, bool scale=true, bool refine=true) {
        const vector<Sample>& pts_row = cf.get_data();
        
        if (scale) {
            double xsf=0;
            double ysf=0;
            for (size_t i=0; i < pts_row.size(); i++) {
                xsf = max(xsf, fabs(pts_row[i].x));
                ysf = max(ysf, fabs(pts_row[i].y));
            }
            cf.xsf = xsf = 1.0/xsf;
            cf.ysf = ysf = 1.0/ysf;
        }
        
        int tdim = cf.dimension();
        MatrixXd cov = MatrixXd::Zero(tdim, tdim);
        VectorXd b = VectorXd::Zero(tdim);
        VectorXd a = VectorXd::Zero(tdim);
        
        VectorXd sol;
        
        for (int iter=0; iter < 1; iter++) {
            cov.setZero();
            b.setZero();
            a.setZero();
            
            for (size_t i=0; i < pts_row.size(); i++) {
                const Sample& sp = pts_row[i];
                double w = sp.weight * sp.yweight;
                a[0] = 1*w;
                double prod = sp.x * cf.xsf; // top poly
                for (int j=1; j <= cf.order_n; j++) {
                    a[j] = prod*w;
                    prod *= sp.x * cf.xsf;
                }
                prod = sp.x*cf.xsf; // bottom poly
                for (int j=1; j <= cf.order_m; j++) {
                    a[j+cf.order_n] = prod*w*sp.y*cf.ysf;
                    prod *= sp.x*cf.xsf;
                }
                
                for (int col=0; col < tdim; col++) { 
                    for (int icol=0; icol < tdim; icol++) { // build covariance of design matrix : A'*A
                        cov(col, icol) += a[col]*a[icol];
                    }
                    b[col] += cf.base_value*a[col]*sp.y*cf.ysf*w; // build rhs of system : A'*b
                }
            }
            
            for (int col=cf.order_n+1; col < cov.cols(); col++) {
                cov.col(col) = -cov.col(col);
            }
            
            sol = cov.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
        }
        
        // now perform non-linear optimization
        
        if (refine) {
            sol = cf.gauss_newton_armijo(sol);
        }
        
        return sol;
    }
    
    double softclamp(double x, double lower, double upper, double p=0.98) {
        double s = (x - lower) / (upper - lower);
        if (s > p) {
            return 1.0/(1.0 + exp(-3.89182*s));
        }
        return s < 0 ? 0 : s;
    }
  
    Point2d& zero;
    Point2d& transverse;
    Point2d& longitudinal;
    std::string wdir;
    std::string prname;
    const cv::Mat& img;
    bool    lpmm_mode;
    double  pixel_size;
    Distance_scale& distance_scale;
};

#endif
