#include "include/ellipse.h"
#include "include/component_labelling.h"
#include "include/point_helpers.h"

#include <Eigen/Cholesky>
#include <Eigen/Dense>
using namespace Eigen;

#include <limits>

typedef Matrix<double, 5, 5> Matrix5d;
typedef Matrix<double, 5, 1> Vector5d;


int Ellipse_detector::fit(const Component_labeller& cl, const Gradient& gradient,
    const Pointlist& raw_points, int tl_x, int tl_y, int dilate) {
    
    
    int width  = gradient.width();
    int height = gradient.height();
    
    const cv::Mat& grad_x = gradient.grad_x();
    const cv::Mat& grad_y = gradient.grad_y();
    const cv::Mat& grad_m = gradient.grad_magnitude();

    set<iPoint> boundary;

    const int border = 1;
    bool edge_touched = false;
    double mx = 0;
    double my = 0;
    for (size_t i=0; i < raw_points.size(); i++) {
        mx += raw_points[i].x;
        my += raw_points[i].y;
        boundary.insert(iPoint(lrint(raw_points[i].x), lrint(raw_points[i].y)) );
        
        if (raw_points[i].x <= border || raw_points[i].x >= width - 1 - border ||
            raw_points[i].y <= border || raw_points[i].y >= height - 1 - border) {
            
            edge_touched = true;
        }
    }
    
    mx /= raw_points.size();
    my /= raw_points.size();
    
    if (edge_touched) return 0; // skip objects that touch the edge of the image (ellipses not really allowed there)

    if (isnan(mx) || isnan(my)) {
        return 0; // 0 -> not a circle
    }

    _dilate(boundary, width, height, dilate);
    _dilate_outer_only(boundary, width, height);

    Matrix<double, 5, 5> wK;
    Matrix<double, 5, 1> K; 
    K.setZero();
    Vector3d l;
    l.setZero();
    Matrix<double, 5, 1> rhs;
    rhs.setZero();

    double sum_c4 = 0;

    wK.setZero();
    int count = (int)boundary.size();
    Matrix<double, Eigen::Dynamic, 3> L(count, 3);
    L.setZero();

    double mean_dist = 0;
    int counter = 0;
    for (set<iPoint>::const_iterator it=boundary.begin(); it != boundary.end(); it++) {
        const int& x_pos = it->first;
        const int& y_pos = it->second;

        if (grad_m.at<float>(y_pos, x_pos) > 1e-7) {
            Vector3d v;
            v << x_pos,y_pos,1;
            double dist = sqrt( SQR(x_pos-mx) + SQR(y_pos-my) );
            mean_dist += dist;
            counter++;
        }
    }
    mean_dist /= counter;
    double iso_scale = sqrt(2.0) / mean_dist;

    size_t idx = 0;
    for (set<iPoint>::const_iterator it=boundary.begin(); it != boundary.end(); it++) {
        const int& x_pos = it->first;
        const int& y_pos = it->second;
        
        if (grad_m.at<float>(y_pos, x_pos) > 1e-7) {


            l[0] = grad_x.at<float>(y_pos, x_pos);
            l[1] = grad_y.at<float>(y_pos, x_pos);
            l[2] = -(grad_x.at<float>(y_pos, x_pos) * (x_pos-mx) * iso_scale + grad_y.at<float>(y_pos, x_pos) * (y_pos-my) * iso_scale);

            L.row(idx++) = l;

        }
    }

    for (size_t r=0; r < (size_t)L.rows(); r++) {

        Vector3d l = L.row(r);

        sum_c4 += SQR(SQR(l[2]));

        K[0] = l[0]*l[0];
        K[1] = l[0]*l[1];
        K[2] = l[1]*l[1];
        K[3] = l[0]*l[2];
        K[4] = l[1]*l[2];

        double weight = (l[0]*l[0] + l[1]*l[1]);

        wK += weight * K * K.transpose();

        rhs -= weight * K*(l[2]*l[2]);
    }
    
    if (rhs.rows() != 5 || rhs.cols() != 1) {
        printf("rhs undefined\n");
        exit(1);
    }
    
    bool has_nans = false;
    for (size_t ri=0; ri < rhs.rows(); ri++) {
        if (isnan(rhs(ri,0))) {
            has_nans = true;   
            return 0;
        }
    }

    Vector5d sol;
    sol = wK.fullPivHouseholderQr().solve(rhs);

    Matrix3d Cstar;
    Cstar(0, 0) = sol[0];
    Cstar(0, 1) = sol[1]*0.5;
    Cstar(0, 2) = sol[3]*0.5;
    Cstar(1, 0) = sol[1]*0.5;
    Cstar(1, 1) = sol[2];
    Cstar(1, 2) = sol[4]*0.5;
    Cstar(2, 0) = sol[3]*0.5;
    Cstar(2, 1) = sol[4]*0.5;
    Cstar(2, 2) = 1;      // F* == 1


    Matrix3d C = Cstar.inverse();
    C *= 1.0/C(2,2);
    _C = C;

    Matrix<double, 1, 5> s;
    s.row(0) = sol;
    double sAAs = (s * wK * (s.transpose()))(0,0);
    double R = (sAAs - 2*(sol.dot(rhs)) + sum_c4) / double(count - 5);
    Matrix<double, 5, 5> cov = wK.inverse();
    cov = cov * R;
    Matrix2d cov2;
    cov2(0, 0) = cov(3, 3);
    cov2(0, 1) = cov(3, 4);
    cov2(1, 0) = cov(4, 3);
    cov2(1, 1) = cov(4, 4);
    JacobiSVD<Matrix2d> svd(cov2, ComputeFullU | ComputeFullV);
    Matrix2d Vs = svd.matrixV();
    Vector2d ws = svd.singularValues();

    Matrix2d V;
    V.row(0) = Vs.row(1);
    V.row(1) = Vs.row(0);

    Vector2d w;
    w[0] = ws[1];
    w[1] = ws[0];

    Vector2d centre_uncertainty;
    centre_uncertainty[0] = sqrt(w[0]) * 0.25;
    centre_uncertainty[1] = sqrt(w[1]) * 0.25;

    // shift the ellipse back to the original pixel coordinate system
    Matrix3d S;
    S.setIdentity();
    S(0, 0) = iso_scale;
    S(1, 1) = iso_scale;
    C = _C = S.transpose()*C*S;
    S.setIdentity();
    S(0, 2) = -(tl_x + mx);
    S(1, 2) = -(tl_y + my);
    C = _C = S.transpose()*C*S;

    int result = _matrix_to_ellipse(C);

    if (minor_axis > major_axis) {
        _C = -_C;
        _matrix_to_ellipse(C);
    }

    bool gradient_ok = gradient_check(cl, gradient, raw_points);

    int is_circle = 0;

    if (result == 0 &&
        major_axis >= min_major_axis &&
        minor_axis >= min_minor_axis &&
        gradient_ok) {

        is_circle = 1;
    }

    printf("centre (%.2lf, %.2lf), major = %lf, minor = %lf, angle = %lf, is_circle = %d\n",
        centroid_x, centroid_y, major_axis, minor_axis, angle/M_PI*180.0, is_circle);

    if (isnan(centroid_x) || isnan(centroid_y) || isnan(major_axis) || isnan(minor_axis)) {
        is_circle = 0;
    }
    
    if (is_circle) {
        scanset.clear();
        for (size_t i=0; i < raw_points.size(); i++) {
            int ix = lrint(raw_points[i].x);
            int iy = lrint(raw_points[i].y);
            
            map<int, scanline>::iterator it = scanset.find(iy);
            if (it == scanset.end()) {
                scanline sl(ix,ix);
                scanset.insert(make_pair(iy, sl));
            }
            if (ix < scanset[iy].start) {
                scanset[iy].start = ix;
            }
            if (ix > scanset[iy].end) {
                scanset[iy].end = ix;
            }
        }
        int clabel = cl(lrint(raw_points[0].x), lrint(raw_points[0].y));
        int total = 0;
        int foreground = 0;
        for (map<int, scanline>::iterator it=scanset.begin(); it != scanset.end(); it++) {
            int y=it->first;
            for (int x=it->second.start; x <= it->second.end; x++) {
                if (cl(x,y) == clabel) {
                    foreground++;
                }
                total++;
            }
        }
        if (foreground >= total - 1) {
            solid = true;
        }
    }
    

    return is_circle;
}


void Ellipse_detector::_dilate(set<iPoint>& s, int width, int height, int iters) {
    const int border = 1;

    if (iters > 0) {

        set<iPoint> gen_s;

        for (int k=0; k < iters; k++) {
            for (set<iPoint>::const_iterator it=s.begin(); it != s.end(); it++) {

                gen_s.insert(*it);

                int left = std::max(border, it->first-1);
                int right = std::min(width-1-border, it->first+1);
                int top = std::max(border, it->second-1);
                int bottom = std::min(height-1-border, it->second+1);

                gen_s.insert(iPoint(right, it->second));
                gen_s.insert(iPoint(right, bottom));
                gen_s.insert(iPoint(it->first, bottom));
                gen_s.insert(iPoint(left, bottom));
                gen_s.insert(iPoint(left, it->second));
                gen_s.insert(iPoint(left, top));
                gen_s.insert(iPoint(it->first, top));
                gen_s.insert(iPoint(right, top));

            }
            s = gen_s; // copy it back
        }
    }
}

void Ellipse_detector::_dilate_outer_only(set<iPoint>& s, int width, int height) {
    const int border = 1;
    
    double cx = 0;
    double cy = 0;
    for (set<iPoint>::const_iterator it=s.begin(); it != s.end(); it++) {
        cx += it->first;
        cy += it->second;
    }
    cx /= s.size();
    cy /= s.size();

    set<iPoint> gen_s;

    for (set<iPoint>::const_iterator it=s.begin(); it != s.end(); it++) {

        gen_s.insert(*it);
        
        Point2d dir(it->first - cx, it->second - cy); // current radial direction

        int left = std::max(border, it->first-1);
        int right = std::min(width-1-border, it->first+1);
        int top = std::max(border, it->second-1);
        int bottom = std::min(height-1-border, it->second+1);

        if ((right - it->first)*dir.x >= 0) gen_s.insert(iPoint(right, it->second));
        if ((right - it->first)*dir.x >= 0 && (bottom - it->second)*dir.y >= 0) gen_s.insert(iPoint(right, bottom));
        if ((bottom - it->second)*dir.y >= 0) gen_s.insert(iPoint(it->first, bottom));
        if ((left - it->first)*dir.x >= 0 && (bottom - it->second)*dir.y >= 0) gen_s.insert(iPoint(left, bottom));
        if ((left - it->first)*dir.x >= 0) gen_s.insert(iPoint(left, it->second));
        if ((left - it->first)*dir.x >= 0 && (top - it->second)*dir.y >= 0) gen_s.insert(iPoint(left, top));
        if ((top - it->second)*dir.y >= 0) gen_s.insert(iPoint(it->first, top));
        if ((right - it->first)*dir.x >= 0 && (top - it->second)*dir.y >= 0) gen_s.insert(iPoint(right, top));

    }
    s = gen_s; // copy it back
}

int Ellipse_detector::_matrix_to_ellipse(Matrix3d& C) {

    double a = C(0,0);
    double b = C(0,1)*2;
    double c = C(1,1);  
    double d = C(0,2)*2;
    double e = C(1,2)*2;
    double f = C(2,2);  


    double thetarad = 0.5*atan2(b, a - c);
    double cost = cos(thetarad);
    double sint = sin(thetarad);
    double sin_squared = sint*sint;
    double cos_squared = cost*cost;
    double cos_sin = sint*cost;

    double Ao = f;
    double Au =   d * cost + e * sint;
    double Av = - d * sint + e * cost;
    double Auu = a * cos_squared + c * sin_squared + b * cos_sin;
    double Avv = a * sin_squared + c * cos_squared - b * cos_sin;

    if(Auu==0 || Avv==0) {
        // invalid ellipse
        return -1;
    }

    double tuCentre = - Au/(2*Auu);
    double tvCentre = - Av/(2*Avv);
    double wCentre = Ao - Auu*tuCentre*tuCentre - Avv*tvCentre*tvCentre;

    double uCentre = tuCentre * cost - tvCentre * sint;
    double vCentre = tuCentre * sint + tvCentre * cost;

    double Ru = -wCentre/Auu;
    double Rv = -wCentre/Avv;

    Ru = sqrt(fabs(Ru))*(Ru < 0 ? -1 : 1);
    Rv = sqrt(fabs(Rv))*(Rv < 0 ? -1 : 1);

    centroid_x = uCentre;
    centroid_y = vCentre;
    major_axis = Ru;
    minor_axis = Rv;
    angle = thetarad;

    return 0;
}

double Ellipse_detector::calculate_curve_length(const Pointlist& points) {

    int n = points.size();

    // calculate the curve length of the  boundary
    double curve_len = 0;
    double prev_x = points[0].x - centroid_x;
    double prev_y = points[0].y - centroid_y;
    for (int i=1; i < n; i++) {
        double x = points[i].x - centroid_x;
        double y = points[i].y - centroid_y;
        curve_len += sqrt(SQR(x - prev_x) + SQR(y - prev_y));
        prev_x = x;
        prev_y = y;
    }
    curve_len += sqrt(SQR(points[0].x - centroid_x - prev_x) +
        SQR(points[0].y - centroid_y - prev_y));

    return curve_len;
}

bool Ellipse_detector::gradient_check(const Component_labeller& cl, const Gradient& gradient, const Pointlist& raw_points) {
    
    // gradient just outside ellipse must be perpendicular to ellipse tangent
    double cosa = cos(-angle);
    double sina = sin(-angle);
    int neighbours[8][2] = {
        {-1, -1}, {0, -1}, {1, -1}, 
        {-1, 0}, {1, 0},
        {-1,  1}, {0,  1}, {1, 1}
    };
    
    // TODO: we need improved handling for objects falling on the edge of the scene?
    
    int not_fg_count = 0;
    vector<double> phi_diff;
    for (size_t i=0; i < raw_points.size(); i++) {
        // now generate a point just outside the ellipse ....
        int ox = lrint(raw_points[i].x);
        int oy = lrint(raw_points[i].y);
        int px = ox;
        int py = oy;
        double maxdist = sqrt((px - centroid_x)*(px - centroid_x) + (py - centroid_y)*(py - centroid_y));
        for (int n=0; n < 8; n++) {
            int lx = px + neighbours[n][0];
            int ly = py + neighbours[n][1];
            
            if (lx >= 5 && lx < gradient.width() - 5 &&
                ly >= 5 && ly < gradient.height() - 5) {
            
                double dist = sqrt((lx - centroid_x)*(lx - centroid_x) + (ly - centroid_y)*(ly - centroid_y));
                if (dist > maxdist) {
                    px = lx;
                    py = ly;
                    maxdist = dist;
                }
            }
        }
        
        // points just outside ellipse should have labels of 0 or -1 (not foreground)
        not_fg_count += cl(px, py) <= 0 ? 1 : 0;
        
        Point2d d(raw_points[i].x - centroid_x, raw_points[i].y - centroid_y);
        double rx = cosa*d.x - sina*d.y;
        double ry = sina*d.x + cosa*d.y;
        double theta = atan2(ry, rx); // ellipse curve parameter theta
        Point2d tangent = normalize(Point2d(-major_axis*sin(theta), minor_axis*cos(theta))); 
        // rotate the tangent vector back to image coordinates
        rx = cosa*tangent.x + sina*tangent.y;
        ry = (-sina)*tangent.x + cosa*tangent.y;
        tangent.x = rx; 
        tangent.y = ry; 
        
        Point2d grad = normalize(Point2d(gradient.grad_x().at<float>(py, px), gradient.grad_y().at<float>(py, px)));
        
        double dot = tangent.x*grad.x + tangent.y*grad.y;
        double phi = acos(dot);
        
        phi_diff.push_back(phi/M_PI*180 - 90);
    }
    if ((raw_points.size() - not_fg_count) > 1) {
        return false; // we can leave early here ...
    }
    
    sort(phi_diff.begin(), phi_diff.end());
    const double phi_percentile = 0.9;
    
    double phi_delta = phi_diff[phi_percentile*phi_diff.size()];
    
    return (phi_delta < max_ellipse_gradient_error) && phi_delta >= 0;
}

Vector3d Ellipse_detector::pose(int imgrows, int imgcols, double circle_diameter, double f, Eigen::Vector3d ext_norm) {
    Matrix3d A(_C);
     
    double fc_x = f;
    double fc_y = f;
    double cc_x = imgcols/2;
    double cc_y = imgrows/2;
    
                     
    Matrix3d K;
    K.setZero();
    K(0, 0) = fc_x;
    K(1, 1) = fc_y;
    K(2, 2) = 1;   
    K(0, 2) = cc_x;
    K(1, 2) = cc_y;

    Matrix3d Cn = (K.transpose()*A*K);

    Cn = -Cn; // ?

    SelfAdjointEigenSolver<Matrix3d> Eeig(Cn);
    Matrix3d Ae = Eeig.eigenvectors();
    Vector3d Aeigs = Eeig.eigenvalues();

    //printf("eigenvector matrix\n");
    //cout << Ae << endl;
    //printf("\neigvals: %lg %lg %lg\n", Aeigs[0], Aeigs[1], Aeigs[2]);

    double l1 = Aeigs[0];
    double l2 = Aeigs[1];
    double l3 = Aeigs[2];

    Matrix3d eigenvectors;
    Vector3d z_axis;
    z_axis << 0.0, 0.0, 1.0;

    int signs[3] = {
        Aeigs[0] < 0 ? 1 : 0,
        Aeigs[1] < 0 ? 1 : 0,
        Aeigs[2] < 0 ? 1 : 0,
    };
    

    if (signs[0] != signs[1]) {
        // last two eigenvalues have same sign
        l3 = Aeigs[0];
        if (z_axis.dot(Ae.col(0)) < 0) {
            eigenvectors.col(2) = -Ae.col(0);
        } else {
            eigenvectors.col(2) = Ae.col(0);
        }
        if (fabs(Aeigs[1]) > fabs(Aeigs[2])) {
            l1 = Aeigs[1];
            l2 = Aeigs[2];
            eigenvectors.col(1) = Ae.col(2);
            eigenvectors.col(0) = eigenvectors.col(1).cross(eigenvectors.col(2));
        } else {
            l1 = Aeigs[2];
            l2 = Aeigs[1];
            eigenvectors.col(1) = Ae.col(1);
            eigenvectors.col(0) = eigenvectors.col(1).cross(eigenvectors.col(2));
        }
    } else {
        // first two eigenvalues have same sign
        l3 = Aeigs[2];
        if (z_axis.dot(Ae.col(2)) < 0) {
            eigenvectors.col(2) = -Ae.col(2);
        } else {
            eigenvectors.col(2) = Ae.col(2);
        }
        if (fabs(Aeigs[0]) > fabs(Aeigs[1])) {
            l1 = Aeigs[0];
            l2 = Aeigs[1];
            eigenvectors.col(1) = Ae.col(1);
            eigenvectors.col(0) = eigenvectors.col(1).cross(eigenvectors.col(2));
        } else {
            l1 = Aeigs[1];
            l2 = Aeigs[0];
            eigenvectors.col(1) = Ae.col(0);
            eigenvectors.col(0) = eigenvectors.col(1).cross(eigenvectors.col(2));
        }
    }
    // l3 contains smallest eigenvalue, eigenvectors.col(2) contains corresponding z-axis
    // l2 contains the middle eigenvalue, plus ev(1) as its eigenvector
    // l1 contains the largest eigenvalue, plut ev(0) as the cross of ev(1) and ev(2)

    Matrix3d R1 = eigenvectors;

    double al1 = fabs(l1);
    double al2 = fabs(l2);
    double al3 = fabs(l3);

    Vector3d s1, v1;
    s1 << sqrt( (al3*(al1 - al2)) / (al1*(al1 + al3)) ),
          0,
          sqrt( (al1*(al2 + al3)) / (al3*(al1 + al3)) );
    v1 << sqrt( (al1 - al2) / (al1 + al3) ),
          0,
          -sqrt( (al2 + al3) / (al1 + al3) );

    Vector3d s2, v2;
    s2 << -sqrt( (al3*(al1 - al2)) / (al1*(al1 + al3)) ),
          0,
          sqrt( (al1*(al2 + al3)) / (al3*(al1 + al3)) );
    v2 << -sqrt( (al1 - al2) / (al1 + al3) ),
          0,
          -sqrt( (al2 + al3) / (al1 + al3) );


    s1 = (0.5 * circle_diameter * s1.transpose()) * R1.transpose();
    v1 = v1.transpose() * R1.transpose();
    s2 = (0.5 * circle_diameter * s2.transpose()) * R1.transpose();
    v2 = v2.transpose() * R1.transpose();
    
    if (ext_norm.norm() > 0.1) { // valid external norm was provided
        if (fabs(v2.dot(ext_norm)) > fabs(v1.dot(ext_norm))) {
            std::swap(s1, s2);
            std::swap(v1, v2);
        }
    } else {
        // use a fixed heuristic instead
        if (v1[0] > 0 && v2[0] < 0) {
            std::swap(s1, s2);
            std::swap(v1, v2);
        }
    }

    //cout << "pos 1 : " << s1.transpose() << " norm 1 : " << v1.transpose() << endl;
    //cout << "pos 2 : " << s2.transpose() << " norm 2 : " << v2.transpose() << endl;

    pos1 = s1;
    pos2 = s2;
    
    n1 = v1;
    n2 = v2;

    return s1;
}
