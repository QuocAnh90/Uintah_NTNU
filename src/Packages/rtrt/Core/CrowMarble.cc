#include <Packages/rtrt/Core/CrowMarble.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Packages/rtrt/Core/Ray.h>
#include <Packages/rtrt/Core/HitInfo.h>
#include <Packages/rtrt/Core/Scene.h>
#include <Packages/rtrt/Core/Stats.h>
#include <Packages/rtrt/Core/Object.h>
#include <Packages/rtrt/Core/Worker.h>
#include <Packages/rtrt/Core/Context.h>
#include <Packages/rtrt/Core/Light.h>

using namespace rtrt;
using namespace SCIRun;

CrowMarble::CrowMarble(double scale,
               const Vector& direction,
               const Color&  c1,
               const Color&  c2,
               const Color&  c3,
               double R0, 
               double phong_exponent)
    : direction(direction), scale(scale), c1(c1), c2(c2), c3(c3),
     spline(c1,c1,c1,c1,c2,c2,c2,c3,c3,c3), R0(R0), phong_exponent(phong_exponent)
{
}

CrowMarble::~CrowMarble()
{
}

void CrowMarble::shade(Color& result, const Ray& ray,
		  const HitInfo& hit, int depth, 
		  double atten, const Color& accumcolor,
		  Context* cx)
{


    double nearest=hit.min_t;
    Object* obj=hit.hit_obj;
    Point hitpos(ray.origin()+ray.direction()*nearest);
    Vector normal(obj->normal(hitpos, hit));
    double cos_prime=-Dot(normal, ray.direction());
    if(cos_prime<0){
	cos_prime=-cos_prime;
	normal=-normal;
    }

    //double pixelsize=0.001;
    Point p(scale*hitpos.x(), scale*hitpos.y(), scale*hitpos.z());
    //double noise=turbulence(p, pixelsize);
    double noise=turbulence(p);
    double w=sin(direction.x()*p.x()+direction.y()*p.y()+
                 direction.z()*p.z()+8*noise);
    //double csp=.45*w+.55;
    double csp = 0.5*(1+w);
    csp*=csp;
    Color Rd=spline(csp);
    Color difflight(0,0,0);
    Color speclight(0,0,0);
    int nlights=cx->scene->nlights();
    cx->stats->ds[depth].nshadow+=nlights;
    double k1 = (1-cos_prime);
    k1 *= k1*k1*k1*k1;
    double ray_objnormal_dot(Dot(ray.direction(),normal));
    for(int i=0;i<nlights;i++){
	Light* light=cx->scene->light(i);
	Vector light_dir=light->get_pos()-hitpos;
	if (ray_objnormal_dot*Dot(normal,light_dir)>0) continue;
//	if (Dot(ray.direction(),light_dir)>0) continue;
	double dist=light_dir.normalize();
	Color shadowfactor(1,1,1);
	if(cx->scene->lit(hitpos, light, light_dir, dist, shadowfactor, depth, cx) ){
            double cos_theta=Dot(light_dir, normal);
	    if(cos_theta < 0){
		cos_theta=-cos_theta;
		light_dir=-light_dir;
	    }
            double k2 = (1-cos_theta);
            k2 *= k2*k2*k2*k2;
	    difflight+=light->get_color()*((1-k1)*(1-k2)*shadowfactor);
            speclight+=light->get_color() * shadowfactor * phong_term( ray.direction(), light_dir,
                                 normal, phong_exponent);
	} else {
	    cx->stats->ds[depth].inshadow++;
	}
    }
    
    Color surfcolor=Rd * (difflight + ambient_hack(cx->scene, normal)*(1-k1))
	+speclight;

    double spec_refl = (R0 + (1-R0)*k1);
    atten *= spec_refl;
    if (depth < cx->scene->maxdepth && atten > 0.02){
            Vector refl_dir = reflection( ray.direction(), normal );
            Ray rray(hitpos, refl_dir);
            Color rcolor;
            cx->worker->traceRay(rcolor, rray, depth+1,  atten,
                                 accumcolor+difflight+speclight, cx);
            surfcolor+= rcolor * spec_refl;
            cx->stats->ds[depth].nrefl++;
    }
    

    result=surfcolor;
}
