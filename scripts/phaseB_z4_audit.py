#!/usr/bin/env python3
"""Fail-closed full-time constraint audit for free-Theta SpacetimeX Z4c."""
import json
import math
import sys
from pathlib import Path

NS=(16,24,32)
FIELDS=("hamiltonian","momentum","theta","ztc","allc")
TOL=2e-8
FLOOR=1e-13

def norm(path):
    rows=[]; header=None
    for line in path.read_text(errors="replace").splitlines():
        line=line.strip()
        if line.startswith("#") and "1:iteration" in line:
            header=[v.strip().split(":",1)[1] for v in line.lstrip("# ").split("\t")]
        elif line and not line.startswith("#") and header:
            words=line.split()
            if len(words)==len(header):
                try: row=dict(zip(header,map(float,words)))
                except ValueError: continue
                if not all(map(math.isfinite,row.values())):
                    raise ValueError(f"Non-finite data: {path}")
                rows.append(row)
    if not rows:
        raise ValueError(f"Missing or unparseable norm file {path}")
    return rows

def only(directory,filename):
    matches=[p for p in directory.rglob("*.tsv") if p.name.lower()==filename.lower()]
    if len(matches)!=1:
        raise ValueError(f"Expected exactly one {filename} under {directory}, got {len(matches)}")
    return norm(matches[0])

def at(rows,t):
    row=min(rows,key=lambda r:abs(r["time"]-t))
    if abs(row["time"]-t)>TOL:
        raise ValueError(f"Missing synchronized time={t}; nearest={row['time']}")
    return row

def col(row,key):
    keys=[k for k in row if k.lower()==key.lower()]
    if len(keys)!=1: raise ValueError(f"Expected one norm column {key}, found {list(row)}")
    return row[keys[0]]

def collect(root):
    results={}
    names={"state":"ellisdustx-state.tsv","hc":"z4c-hc.tsv",
           "mtc":"z4c-mtc.tsv","theta":"z4c-theta.tsv",
           "ztc":"z4c-ztc.tsv","allc":"z4c-allc.tsv"}
    for n in NS:
        files={k:only(root/f"N{n}",name) for k,name in names.items()}
        ts=[r["time"] for r in files["theta"]]
        if len(ts)<10 or any(a>=b for a,b in zip(ts,ts[1:])) or abs(ts[-1]-.18)>TOL:
            raise ValueError(f"Incomplete or unordered free-Z4c history at N={n}")
        history=[]
        for t in ts:
            rows={k:at(v,t) for k,v in files.items()}
            s=rows["state"]
            vals={
                "time":t,
                "rho_min":col(s,"EllisDustX::rho8pi.min"),
                "R_min":col(s,"EllisDustX::R_Ellis.min"),
                "hamiltonian":col(rows["hc"],"Z4c::HC.maxabs"),
                "momentum":max(col(rows["mtc"],f"Z4c::MtC{c}.maxabs") for c in "xyz"),
                "theta":col(rows["theta"],"Z4c::Theta.maxabs"),
                "ztc":max(col(rows["ztc"],f"Z4c::ZtC{c}.maxabs") for c in "xyz"),
                "allc":col(rows["allc"],"Z4c::allC.maxabs")}
            if vals["rho_min"]<=0 or vals["R_min"]<=0:
                raise ValueError(f"Nonpositive state N={n} t={t}")
            if any(vals[k]<0 or not math.isfinite(vals[k]) for k in FIELDS):
                raise ValueError(f"Invalid norm N={n} t={t}")
            history.append(vals)
        results[str(n)]=history
    reference=[r["time"] for r in results["16"]]
    for n in (24,32):
        other=[r["time"] for r in results[str(n)]]
        if len(other)!=len(reference) or any(abs(a-b)>TOL for a,b in zip(reference,other)):
            raise ValueError(f"Nonmatching physical time samples N={n}")
    return results

def judge(data):
    samples=data["16"]
    if max(r["theta"] for r in samples)<1e-16:
        raise ValueError("Theta identically zero: no free-Z4c evidence")
    times=[r["time"] for r in samples]
    summary={}
    for field in FIELDS:
        entries=[]
        for i,t in enumerate(times):
            if t>=.02-TOL and data["16"][i][field]>FLOOR:
                entries.append((t,data["16"][i][field],data["32"][i][field]))
        if field in ("theta","ztc","allc") and not entries:
            raise ValueError(f"{field}: no above-floor data; cannot certify")
        rate=sum(f<c for _,c,f in entries)/len(entries) if entries else None
        if rate is not None and rate < .90:
            raise ValueError(f"{field}: nonconvergent full-time samples {rate:.1%}")
        checkpoints={}
        for t in (.09,.18):
            ii=min(range(len(times)),key=lambda i:abs(times[i]-t))
            if abs(times[ii]-t)>TOL: raise ValueError(f"Missing checkpoint t={t}")
            a,b,c=[data[str(n)][ii][field] for n in NS]
            if field in ("theta","ztc","allc") and a<=FLOOR:
                raise ValueError(f"{field} unresolved at checkpoint t={t}")
            if a>FLOOR and not c<a:
                raise ValueError(f"{field} fails refinement at t={t}")
            checkpoints[f"{t:.2f}"]={"N16":a,"N24":b,"N32":c,
                "order16to32":math.log2(a/c) if c>0 and a>0 else None}
        summary[field]={"comparable_samples":len(entries),
            "improvement_fraction":rate,"checkpoints":checkpoints}
    return summary

def selftest():
    ts=[i*.01 for i in range(19)]
    d={str(n):[{"time":t,**{k:(16/n)**4*(1+t) for k in FIELDS},
        "rho_min":1,"R_min":1} for t in ts] for n in NS}
    assert judge(d)["theta"]["improvement_fraction"]==1
    d["32"]=[dict(r,theta=5.0) for r in d["32"]]
    try: judge(d)
    except ValueError: pass
    else: raise AssertionError("Nonconvergent injected Theta accepted")
    d["16"]=[dict(r,theta=0.0) for r in d["16"]]
    try: judge(d)
    except ValueError: pass
    else: raise AssertionError("Zero injected Theta accepted")
    print("PASS_PHASEB_ADVERSARIAL_SELFTEST")

def main():
    if len(sys.argv)>1 and sys.argv[1]=="--selftest": return selftest()
    root=Path(sys.argv[1] if len(sys.argv)>1 else "workspace/Cactus/external_results")
    data=collect(root)
    report={"phase":"B_dust_free_Z4c_fulltime","resolutions":list(NS),
        "samples_per_resolution":{str(n):len(data[str(n)]) for n in NS},
        "history":data}
    try:
        report["certificate"]=judge(data);report["status"]="PASS"
    except ValueError as exc:
        report["status"]="REFUSE";report["refusal_reason"]=str(exc)
    Path("external_phaseB_metrics.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],
        "samples":report["samples_per_resolution"],
        "certificate":report.get("certificate"),
        "refusal_reason":report.get("refusal_reason")},indent=2))
    if report["status"]!="PASS":
        raise SystemExit("REFUSE_PHASE_B: "+report["refusal_reason"])
    print("PASS_PHASE_B_EXTERNAL_FREE_Z4C_FULLTIME")
    print("Dust + Ellis scale only; no Boltzmann transport or memory theorem.")

if __name__=="__main__":
    main()
