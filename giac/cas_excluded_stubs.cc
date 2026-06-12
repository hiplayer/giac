// Stubs for giac modules excluded from the CAS build (plot3d.cc, cocoa.cc).

#include "giacPCH.h"
#include "cocoa.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif

  // --- cocoa.cc (Groebner / CoCoA integration) ---

  bool f5(vectpoly &, const gen &) { return false; }

  bool cocoa_gbasis(vectpoly &, const gen &) { return false; }

  vecteur cocoa_in_ideal(const vectpoly & r, const vectpoly &, const gen &) {
    return vecteur(r.size(), -1);
  }

  bool cocoa_greduce(const vectpoly &, const vectpoly &, const gen &, vectpoly &) {
    return false;
  }

  bool gbasis8(const vectpoly &, order_t &, vectpoly &, environment *, bool, bool,
               int &, GIAC_CONTEXT, gbasis_param_t) {
    return false;
  }

  bool greduce8(const vectpoly &, const vectpoly &, order_t &, vectpoly &,
                environment *, GIAC_CONTEXT) {
    return false;
  }

  longlong memory_usage() { return 0; }

  // --- plot3d.cc (3D geometry / plotting) ---

  static gen plot3d_unavailable(GIAC_CONTEXT) {
    return gensizeerr(gettext("3D plotting/geometry not available in CAS build"));
  }

  gen _plan(const gen & args, GIAC_CONTEXT) {
    if (args.type == _STRNG && args.subtype == -1) return args;
    return plot3d_unavailable(contextptr);
  }

  gen _sphere(const gen & args, GIAC_CONTEXT) {
    if (args.type == _STRNG && args.subtype == -1) return args;
    return plot3d_unavailable(contextptr);
  }

  gen _cylindre(const gen & args, GIAC_CONTEXT) {
    if (args.type == _STRNG && args.subtype == -1) return args;
    return plot3d_unavailable(contextptr);
  }

  static const char _plan_s[] = "plan";
  static define_unary_function_eval(__plan, &_plan, _plan_s);
  define_unary_function_ptr5(at_plan, alias_at_plan, &__plan, 0, true);

  static const char _sphere_s[] = "sphere";
  static define_unary_function_eval(__sphere, &_sphere, _sphere_s);
  define_unary_function_ptr5(at_sphere, alias_at_sphere, &__sphere, 0, true);

  static const char _cylindre_s[] = "cylindre";
  static define_unary_function_eval(__cylindre, &_cylindre, _cylindre_s);
  define_unary_function_ptr5(at_cylindre, alias_at_cylindre, &__cylindre, 0, true);

  bool is3d(const gen &) { return false; }

  gen do_point3d(const gen & g) { return g; }

  vecteur hyperplan_normal(const gen &) { return vecteur(0); }

  bool hyperplan_normal_point(const gen &, vecteur &, vecteur &) { return false; }

  bool normal3d(const gen &, vecteur &, vecteur &) { return false; }

  bool perpendiculaire_commune(const gen &, const gen &, gen &, gen &, vecteur &,
                               GIAC_CONTEXT) {
    return false;
  }

  gen similitude3d(const vecteur &, const gen &, const gen &, const gen &, int,
                   GIAC_CONTEXT) {
    return undef;
  }

  gen hyperplan2hypersurface(const gen & g) { return g; }

  gen hypersphere2hypersurface(const gen & g) { return g; }

  gen hypersphere_equation(const gen & g, const vecteur &) { return g; }

  gen hypersurface_equation(const gen & g, const vecteur &, GIAC_CONTEXT) {
    return g;
  }

  vecteur interpolyedre(const vecteur &, const gen &, GIAC_CONTEXT) {
    return vecteur(0);
  }

  vecteur interdroitehyperplan(const gen &, const gen &, GIAC_CONTEXT) {
    return vecteur(0);
  }

  vecteur interhyperplan(const gen &, const gen &, GIAC_CONTEXT) {
    return vecteur(0);
  }

  vecteur interhypersurfacecurve(const gen &, const gen &, GIAC_CONTEXT) {
    return vecteur(0);
  }

  vecteur inter2hypersurface(const gen &, const gen &, GIAC_CONTEXT) {
    return vecteur(0);
  }

  vecteur interplansphere(const gen &, const gen &, GIAC_CONTEXT) {
    return vecteur(0);
  }

  vecteur rand_3d() { return vecteur(0); }

  gen plotparam3d(const gen &, const vecteur &, double, double, double, double,
                  double, double, double, double, double, double, bool, bool,
                  const vecteur &, double, double, const gen &, const vecteur &,
                  GIAC_CONTEXT) {
    return undef;
  }

  gen plotimplicit(const gen &, const gen &, const gen &, const gen &, double,
                   double, double, double, double, double, int, int, int, double,
                   const vecteur &, bool, const context *) {
    return undef;
  }

#ifndef NO_NAMESPACE_GIAC
}
#endif
