#!/usr/bin/env python3
"""Fail-closed Phase C Bianchi-I dust + 2 null-stream external Einstein gate.

This test does not include inhomogeneous photon streaming, nonzero beam
momentum, scattering, or full Einstein-Boltzmann transport.
"""
import json
import math
import sys
from pathlib import Path

NS=(16,24,32)
END=.18
TOL=2e-8

def load(path):
    header=None
    rows=[]
    for line in path.read_text(errors="replace").splitlines():
        s=line.strip()
        if s.startswith("#") and "1:iteration" in s:
            header=[c.strip().split(":",1)[1] for c in s.lstrip("# ").split("\t")]
        elif s and not s.startswith("#") and header:
            ss=s.split()
            if len(ss)!=len(header): continue
            try: record=dict(zip(header,map(float,ss)))
            except ValueError: continue
            if not all(math.isfinite(v) for v in record.values()):
                raise ValueError("Nonfinite CarpetX diagnostics in "+str(path))
            rows.append(record)
    if not rows: raise ValueError("No CarpetX norm rows: "+str(path))
    return rows

def only(base,name):
    matches=[p for p in base.rglob("*.tsv") if p.name.lower()==name.lower()]
    if len(matches)!=1:
        raise ValueError("Expected exactly one "+name+" in "+str(base)+", got "+str(len(matches)))
    return load(matches[0])

def entry(row,name):
    keys=[k for k in row if k.lower()==name.lower()]
    if len(keys)!=1: raise ValueError("Missing CarpetX column "+name)
    return row[keys[0]]

def at(rows,t):
    item=min(rows,key=lambda r:abs(r["time"]-t))
    if abs(item["time"]-t)>TOL: raise ValueError("Time mismatch at "+str(t))
    return item

def collect(root):
    data={}
    paths={
        "matter":"ellisbeamx-state.tsv",
        "metric":"admbasex-metric.tsv",
        "hamiltonian":"z4c-hc.tsv",
        "momentum":"z4c-mtc.tsv",
        "theta":"z4c-theta.tsv",
        "ztc":"z4c-ztc.tsv",
        "combined":"z4c-allc.tsv",
    }
    for n in NS:
        base=root/("N"+str(n))
        groups={k:only(base,f) for k,f in paths.items()}
        times=[r["time"] for r in groups["matter"]]
        if len(times)<10 or abs(times[0])>TOL or abs(times[-1]-END)>TOL:
            raise ValueError("Incomplete history N="+str(n))
        if any(a>=b for a,b in zip(times,times[1:])):
            raise ValueError("Nonmonotone times")
        hist=[]
        for t in times:
            r={k:at(g,t) for k,g in groups.items()}
            m=r["matter"];g=r["metric"]
            get=lambda v,part:entry(m,"EllisBeamX::"+v+"."+part)
            gv=lambda v,part:entry(g,"ADMBaseX::"+v+"."+part)
            radp=get("beamplus8pi","min")
            radm=get("beamminus8pi","min")
            h=entry(r["hamiltonian"],"Z4c::HC.maxabs")
            z=max(entry(r["ztc"],"Z4c::ZtC"+c+".maxabs") for c in "xyz")
            mo=max(entry(r["momentum"],"Z4c::MtC"+c+".maxabs") for c in "xyz")
            hist.append(dict(time=t, dust=get("rho8pi","min"),
                 radplus=radp,radminus=radm,
                 R=get("R_Ellis","min"),
                 gxx=gv("gxx","min"),gyy=gv("gyy","min"),gzz=gv("gzz","min"),
                 gxxmax=gv("gxx","max"),gyymax=gv("gyy","max"),
                 H=h,M=mo,Theta=entry(r["theta"],"Z4c::Theta.maxabs"),
                 ZtC=z,allC=entry(r["combined"],"Z4c::allC.maxabs")))
        data[str(n)]=hist
    ref=[r["time"] for r in data["16"]]
    for n in NS[1:]:
        arr=[r["time"] for r in data[str(n)]]
        if len(ref)!=len(arr) or any(abs(a-b)>TOL for a,b in zip(ref,arr)):
            raise ValueError("Cross-resolution time grid mismatch")
    return data

def check(histories):
    outcome={}
    for n in NS:
        h=histories[str(n)]
        first=h[0]
        m0=first["dust"]*first["R"]**3
        ep0=first["radplus"]*first["R"]**3*math.sqrt(first["gxx"])
        em0=first["radminus"]*first["R"]**3*math.sqrt(first["gxx"])
        if not (m0>0 and ep0>0 and em0>0): raise ValueError("Radiation/dust absent")
        if abs(ep0-em0)>1e-10: raise ValueError("Initial counterstreams imbalanced")
        if abs(first["gxx"]-1)>1e-7 or abs(first["gyy"]-1)>1e-7:
            raise ValueError("Unexpected Cauchy geometry")
        if abs(m0+ep0+em0-3)>1e-5: raise ValueError("Bad initial Hamiltonian source")
        errors={"dust":0.0,"radplus":0.0,"radminus":0.0}
        residuals={"H":0.0,"M":0.0,"Theta":0.0,"ZtC":0.0,"allC":0.0}
        for row in h:
            d,eplus,emin,R,gx,gy,gz=(row[k] for k in
                ("dust","radplus","radminus","R","gxx","gyy","gzz"))
            if not all(math.isfinite(v) and v>0 for v in (d,eplus,emin,R,gx,gy,gz)):
                raise ValueError("Nonpositive or nonfinite evolved matter/metric N="+str(n))
            if abs(row["gxxmax"]-gx)>1e-6 or abs(row["gyymax"]-gy)>1e-6:
                raise ValueError("Lost homogeneous Bianchi-I control")
            if abs(eplus-emin)>1e-8:
                raise ValueError("Counterpropagating stream imbalance")
            quantities={
                "dust":d*R**3,
                "radplus":eplus*R**3*math.sqrt(gx),
                "radminus":emin*R**3*math.sqrt(gx),
            }
            for key,expected in (("dust",m0),("radplus",ep0),("radminus",em0)):
                errors[key]=max(errors[key],abs(quantities[key]/expected-1))
            for k in residuals:
                value=row[k]
                if not math.isfinite(value) or value<0:
                    raise ValueError("Invalid Z4c constraint norm")
                residuals[k]=max(residuals[k],value)
        if max(errors.values())>=5e-4:
            raise ValueError("Covariant conservation invariant drift N="+str(n)+" "+str(errors))
        if max(residuals.values())>=1e-2:
            raise ValueError("Large Einstein constraint norm N="+str(n))
        last=h[-1]
        if last["radplus"] >= first["radplus"]:
            raise ValueError("Missing radiation redshift")
        if abs(last["gxx"]-last["gyy"])<1e-4:
            raise ValueError("No demonstrated anisotropic gravitational backreaction")
        if abs(last["R"]**3-math.sqrt(last["gxx"]*last["gyy"]*last["gzz"]))>5e-4:
            raise ValueError("Volume scale inconsistent with physical spatial metric")
        outcome[str(n)]={"max_conservation_relative_error":errors,
            "maximum_Z4c_norms":residuals,"gxx_end":last["gxx"],
            "gyy_end":last["gyy"],"anisotropy_end":last["gxx"]-last["gyy"],
            "rho_end":last["dust"],"radplus_end":last["radplus"],
            "radminus_end":last["radminus"]}
    return outcome

def selftest():
    demo={}
    for n in NS:
        rows=[]
        for i in range(19):
            t=.01*i
            R=1+0.7*t
            gx=(1+t)**2
            gy=gz=R**3/math.sqrt(gx)
            rows.append({"time":t,"dust":2.88/R**3,
                "radplus":.06/(R**3*math.sqrt(gx)),
                "radminus":.06/(R**3*math.sqrt(gx)),
                "R":R,"gxx":gx,"gyy":gy,"gzz":gz,
                "gxxmax":gx,"gyymax":gy,
                **{key:1e-6 for key in ("H","M","Theta","ZtC","allC")}})
        demo[str(n)]=rows
    assert len(check(demo))==3
    d=json.loads(json.dumps(demo));d["32"][-1]["radplus"]*=1.03
    try: check(d)
    except ValueError: pass
    else: raise AssertionError("Accepted violated photon conservation")
    d=json.loads(json.dumps(demo))
    for row in d["24"]: row["radplus"]=row["radminus"]=0
    try: check(d)
    except ValueError: pass
    else: raise AssertionError("Accepted absent photons")
    print("PASS_PHASEC_ADVERSARIAL_SELFTEST")

def main():
    if len(sys.argv)>1 and sys.argv[1]=="--selftest":return selftest()
    root=Path(sys.argv[1] if len(sys.argv)>1 else "workspace/Cactus/external_results")
    data=collect(root)
    result={"phase":"C_homogeneous_BianchiI_two_null_streams",
        "resolutions":list(NS),"history":data}
    try:
        result["certificate"]=check(data)
        result["status"]="PASS"
    except ValueError as e:
        result["status"]="REFUSE"
        result["refusal_reason"]=str(e)
    Path("external_phaseC_metrics.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],
        "certificate":result.get("certificate"),
        "refusal_reason":result.get("refusal_reason")},indent=2))
    if result["status"]!="PASS":
        raise SystemExit("REFUSE_PHASE_C: "+result["refusal_reason"])
    print("PASS_PHASE_C_EXTERNAL_TWO_NULL_STREAMS")
    print("Only homogeneous two-stream reduction; not full Einstein-Boltzmann.")

if __name__=="__main__":main()
