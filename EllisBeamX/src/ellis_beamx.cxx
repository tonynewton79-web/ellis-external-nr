#include <loop.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

// Exact homogeneous Bianchi-I two-stream matter reduction:
// comoving dust plus two oppositely directed null dust beams (x direction).
namespace EllisBeamX {
using namespace Loop;
extern "C" void EllisBeamX_InitialGeometry(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_EllisBeamX_InitialGeometry;
  DECLARE_CCTK_PARAMETERS;
  GF3D<CCTK_REAL,0,0,0> gx(cctkGH,gxx),xy(cctkGH,gxy),xz(cctkGH,gxz),gy(cctkGH,gyy),yz(cctkGH,gyz),gz(cctkGH,gzz);
  GF3D<CCTK_REAL,0,0,0> kx(cctkGH,kxx),kxy_(cctkGH,kxy),kxz_(cctkGH,kxz),ky(cctkGH,kyy),kyz_(cctkGH,kyz),kz(cctkGH,kzz);
  GF3D<CCTK_REAL,0,0,0> rd(cctkGH,rho8pi),ep(cctkGH,beamplus8pi),em(cctkGH,beamminus8pi),R(cctkGH,R_Ellis);
  loop_all<0,0,0>(cctkGH,[&](const PointDesc&p){
    gx(p.I)=1;xy(p.I)=0;xz(p.I)=0;gy(p.I)=1;yz(p.I)=0;gz(p.I)=1;
    kx(p.I)=-H0;kxy_(p.I)=0;kxz_(p.I)=0;ky(p.I)=-H0;kyz_(p.I)=0;kz(p.I)=-H0;
    // H constraint: 8 pi (rho_d + E_+ + E_-) = 3 H0^2 for flat metric.
    rd(p.I)=3*H0*H0-rad8pi;
    ep(p.I)=rad8pi/2;em(p.I)=rad8pi/2;R(p.I)=1;
  });
}
extern "C" void EllisBeamX_InitialGauge(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_EllisBeamX_InitialGauge;
  GF3D<CCTK_REAL,0,0,0> a(cctkGH,alp),da(cctkGH,dtalp),bx(cctkGH,betax),by(cctkGH,betay),bz(cctkGH,betaz),dbx(cctkGH,dtbetax),dby(cctkGH,dtbetay),dbz(cctkGH,dtbetaz);
  loop_all<0,0,0>(cctkGH,[&](const PointDesc&p){
    a(p.I)=1;da(p.I)=0;bx(p.I)=0;by(p.I)=0;bz(p.I)=0;
    dbx(p.I)=0;dby(p.I)=0;dbz(p.I)=0;
  });
}
extern "C" void EllisBeamX_AddTmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_EllisBeamX_AddTmunu;
  GF3D<const CCTK_REAL,0,0,0> rd(cctkGH,rho8pi),ep(cctkGH,beamplus8pi),em(cctkGH,beamminus8pi),gx(cctkGH,gxx);
  GF3D<CCTK_REAL,0,0,0> tt(cctkGH,eTtt),tx(cctkGH,eTtx),ty(cctkGH,eTty),tz(cctkGH,eTtz),xx(cctkGH,eTxx),xy(cctkGH,eTxy),xz(cctkGH,eTxz),yy(cctkGH,eTyy),yz(cctkGH,eTyz),zz(cctkGH,eTzz);
  const CCTK_REAL inv8pi=1/(8*std::acos(CCTK_REAL(-1)));
  loop_int<0,0,0>(cctkGH,[&](const PointDesc&p){
    const CCTK_REAL E=ep(p.I)+em(p.I);
    tt(p.I)+=(rd(p.I)+E)*inv8pi;
    tx(p.I)-=(ep(p.I)-em(p.I))*std::sqrt(gx(p.I))*inv8pi;
    ty(p.I)+=0;tz(p.I)+=0;
    xx(p.I)+=E*gx(p.I)*inv8pi;
    xy(p.I)+=0;xz(p.I)+=0;yy(p.I)+=0;yz(p.I)+=0;zz(p.I)+=0;
  });
}
extern "C" void EllisBeamX_SetTmunuRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_EllisBeamX_SetTmunuRHS;
  GF3D<const CCTK_REAL,0,0,0> rd(cctkGH,rho8pi),ep(cctkGH,beamplus8pi),em(cctkGH,beamminus8pi),gx(cctkGH,gxx);
  GF3D<CCTK_REAL,0,0,0> tt(cctkGH,eTtt),tx(cctkGH,eTtx),ty(cctkGH,eTty),tz(cctkGH,eTtz),xx(cctkGH,eTxx),xy(cctkGH,eTxy),xz(cctkGH,eTxz),yy(cctkGH,eTyy),yz(cctkGH,eTyz),zz(cctkGH,eTzz);
  const CCTK_REAL inv8pi=1/(8*std::acos(CCTK_REAL(-1)));
  loop_int<0,0,0>(cctkGH,[&](const PointDesc&p){
    const CCTK_REAL E=ep(p.I)+em(p.I);
    tt(p.I)=(rd(p.I)+E)*inv8pi;
    tx(p.I)=-(ep(p.I)-em(p.I))*std::sqrt(gx(p.I))*inv8pi;
    ty(p.I)=0;tz(p.I)=0;
    xx(p.I)=E*gx(p.I)*inv8pi;
    xy(p.I)=0;xz(p.I)=0;yy(p.I)=0;yz(p.I)=0;zz(p.I)=0;
  });
}
extern "C" void EllisBeamX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_EllisBeamX_RHS;
  GF3D<const CCTK_REAL,0,0,0> rd(cctkGH,rho8pi),ep(cctkGH,beamplus8pi),em(cctkGH,beamminus8pi),R(cctkGH,R_Ellis);
  GF3D<CCTK_REAL,0,0,0> drd(cctkGH,rho8pi_rhs),dep(cctkGH,beamplus8pi_rhs),dem(cctkGH,beamminus8pi_rhs),dR(cctkGH,R_Ellis_rhs);
  GF3D<const CCTK_REAL,0,0,0> gx(cctkGH,gxx),xy(cctkGH,gxy),xz(cctkGH,gxz),gy(cctkGH,gyy),yz(cctkGH,gyz),gz(cctkGH,gzz);
  GF3D<const CCTK_REAL,0,0,0> kx(cctkGH,kxx),kxy_(cctkGH,kxy),kxz_(cctkGH,kxz),ky(cctkGH,kyy),kyz_(cctkGH,kyz),kz(cctkGH,kzz),a(cctkGH,alp);
  loop_int<0,0,0>(cctkGH,[&](const PointDesc&p){
    const CCTK_REAL x=gx(p.I),b=xy(p.I),c=xz(p.I),d=gy(p.I),e=yz(p.I),f=gz(p.I);
    const CCTK_REAL det=x*(d*f-e*e)-b*(b*f-c*e)+c*(b*e-c*d);
    const CCTK_REAL ixx=(d*f-e*e)/det,iyy=(x*f-c*c)/det,izz=(x*d-b*b)/det,
      ixy=(c*e-b*f)/det,ixz=(b*e-c*d)/det,iyz=(b*c-x*e)/det;
    const CCTK_REAL tr=ixx*kx(p.I)+iyy*ky(p.I)+izz*kz(p.I)
      +2*(ixy*kxy_(p.I)+ixz*kxz_(p.I)+iyz*kyz_(p.I));
    const CCTK_REAL kxxMixed=ixx*kx(p.I)+ixy*kxy_(p.I)+ixz*kxz_(p.I);
    // 8pi rho_d R^3 and 8pi E_+/- R^3 sqrt(g_xx) are first integrals.
    drd(p.I)=a(p.I)*tr*rd(p.I);
    dep(p.I)=a(p.I)*(tr+kxxMixed)*ep(p.I);
    dem(p.I)=a(p.I)*(tr+kxxMixed)*em(p.I);
    dR(p.I)=-a(p.I)*tr*R(p.I)/3;
  });
}
extern "C" void EllisBeamX_Boundaries(CCTK_ARGUMENTS) {}
} // namespace EllisBeamX
