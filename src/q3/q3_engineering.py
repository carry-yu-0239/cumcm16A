#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CUMCM 2016 A题 Q3 工程计算原型；不写论文。

对外符号按 04_symbols.tex：theta_i, theta_d, varphi_a, delta, r。
手稿张力方向角内部改用水平/竖直张力分量，避免与 theta_i 冲突。

工程闭合（尚未写入 model card）：
- 钢管/钢桶姿态由手稿力矩方程求，不使用与其冲突的孤立 tan(phi_i) 式；
- F_cb 严格按手稿 374*delta*v_c^2；
- F_cs 用钢球球体投影面积；
- 锚链受流特征直径由单位长度钢质量换算等面积圆杆直径。
"""
from dataclasses import dataclass, replace
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "q3"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
LOGS = OUT / "logs"
PROCESSED = ROOT / "data" / "processed"
for d in (TABLES, FIGURES, LOGS, PROCESSED): d.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.family"] = ["Noto Serif CJK JP", "DejaVu Serif", "serif"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.dpi"] = 220

CHAIN = {
    "I": (0.078, 3.2), "II": (0.105, 7.0), "III": (0.120, 12.5),
    "IV": (0.150, 19.5), "V": (0.180, 28.12),
}

@dataclass(frozen=True)
class P:
    g: float = 9.8; rho: float = 1025.0; rho_steel: float = 7850.0
    H: float = 18.0; vw: float = 36.0; vc: float = 0.0
    ms: float = 2250.0; chain_type: str = "II"; n: int = 210
    buoy_r: float = 1.0; buoy_h: float = 2.0; buoy_m: float = 1000.0
    pipe_L: float = 1.0; pipe_D: float = 0.05; pipe_m: float = 10.0
    barrel_L: float = 1.0; barrel_D: float = 0.30; barrel_m: float = 100.0

    @property
    def ell(self): return CHAIN[self.chain_type][0]
    @property
    def mu(self): return CHAIN[self.chain_type][1]
    @property
    def L(self): return self.n*self.ell
    @property
    def w_chain(self): return self.mu*self.g*(1-self.rho/self.rho_steel)
    @property
    def d_chain(self): return math.sqrt(4*(self.mu/self.rho_steel)/math.pi)
    @property
    def W_pipe(self):
        V = math.pi*(self.pipe_D/2)**2*self.pipe_L
        return self.pipe_m*self.g-self.rho*self.g*V
    @property
    def W_barrel(self):
        V = math.pi*(self.barrel_D/2)**2*self.barrel_L
        return self.barrel_m*self.g-self.rho*self.g*V
    @property
    def W_ball(self): return self.ms*self.g*(1-self.rho/self.rho_steel)
    @property
    def r_ball(self): return (3*(self.ms/self.rho_steel)/(4*math.pi))**(1/3)

class BadState(RuntimeError): pass

def member_angle(Th, Tv, W, k):
    """(Tv-W/2)sin(beta)=(Th+Fc/2)cos(beta), Fc=k cos(beta)."""
    vm = Tv-W/2
    if vm <= 0: raise BadState("构件中点竖向张力非正")
    if abs(k) < 1e-15: return math.atan2(Th, vm), 0.0
    beta = math.atan2(Th+k/2, vm)
    for _ in range(30):
        b2 = math.atan2(Th+k*math.cos(beta)/2, vm)
        if abs(b2-beta) < 1e-13: beta=b2; break
        beta=b2
    else:
        f=lambda b: vm*math.sin(b)-(Th+k*math.cos(b)/2)*math.cos(b)
        beta=brentq(f,0,math.pi/2-1e-10)
    return beta, k*math.cos(beta)

def solve_chain(Th, Tv, p):
    if Tv <= 0: raise BadState("锚链顶部竖向张力非正")
    x=[0.0]; z=[0.0]; seabed=0.0; suspended=0.0; touchdown=False
    for _ in range(p.n):
        if touchdown or Tv <= 1e-12:
            dx,dz=p.ell,0.0; seabed+=p.ell; touchdown=True
        else:
            W=p.w_chain*p.ell
            if Tv-W >= 0:
                beta,Fc=member_angle(Th,Tv,W,374*p.d_chain*p.vc**2*p.ell)
                dx=p.ell*math.sin(beta); dz=p.ell*math.cos(beta)
                Th+=Fc; Tv-=W; suspended+=p.ell
            else:
                ls=max(0.0,Tv/p.w_chain); lb=p.ell-ls
                if ls>1e-12:
                    beta,Fc=member_angle(Th,Tv,p.w_chain*ls,374*p.d_chain*p.vc**2*ls)
                    dx=ls*math.sin(beta)+lb; dz=ls*math.cos(beta); Th+=Fc
                    suspended+=ls; seabed+=lb
                else: dx,dz=p.ell,0.0; seabed+=p.ell
                Tv=0.0; touchdown=True
        x.append(x[-1]+dx); z.append(z[-1]+dz)
    phi=0.0 if seabed>1e-9 else math.atan2(max(Tv,0),Th)
    return dict(x=np.array(x),z=np.array(z),span=x[-1],drop=z[-1],
                suspended=suspended,seabed=seabed,phi=phi,
                regime="存在卧底段" if seabed>1e-9 else "全部悬空")

def eval_state(delta,p):
    if not 0<delta<p.buoy_h: raise BadState("delta越界")
    Fb=p.rho*p.g*math.pi*p.buoy_r**2*delta; Gb=p.buoy_m*p.g
    Fw=0.625*((p.buoy_h-delta)*2*p.buoy_r)*p.vw**2
    Fcb=374*delta*p.vc**2
    Th=Fw+Fcb; Tv=Fb-Gb
    if Th<=0 or Tv<=0: raise BadState("顶部张力非正")
    angles=[]; rv=0.0; rh=0.0
    for _ in range(4):
        b,Fc=member_angle(Th,Tv,p.W_pipe,374*p.pipe_D*p.vc**2*p.pipe_L)
        angles.append(b); rv+=p.pipe_L*math.cos(b); rh+=p.pipe_L*math.sin(b)
        Th+=Fc; Tv-=p.W_pipe
    bd,Fcd=member_angle(Th,Tv,p.W_barrel,374*p.barrel_D*p.vc**2*p.barrel_L)
    rv+=p.barrel_L*math.cos(bd); rh+=p.barrel_L*math.sin(bd)
    Fcs=374*(math.pi*p.r_ball**2)*p.vc**2
    chain=solve_chain(Th+Fcd+Fcs, Tv-p.W_barrel-p.W_ball, p)
    return dict(delta=delta,H=delta+rv+chain["drop"],r=rh+chain["span"],
                theta=np.degrees(angles),theta_d=math.degrees(bd),
                phi_a=math.degrees(chain["phi"]),chain=chain)

def solve(p):
    grid=np.linspace(0.30,1.95,120); prev=None; br=None
    for d in grid:
        try: f=eval_state(float(d),p)["H"]-p.H
        except BadState: continue
        if prev and prev[1]*f<=0: br=(prev[0],d); break
        prev=(d,f)
    if br is None: raise BadState("未找到H(delta)=目标水深的根")
    root=brentq(lambda d: eval_state(d,p)["H"]-p.H,*br,xtol=1e-12)
    return eval_state(root,p)

def row(p,s):
    return {
        "锚链型号":p.chain_type,"链环长度(mm)":p.ell*1000,"单位长度质量(kg/m)":p.mu,
        "锚链节数 n":p.n,"锚链总长度 L(m)":p.L,"海面风速 v_w(m/s)":p.vw,
        "海水流速 v_c(m/s)":p.vc,"海水深度 H(m)":p.H,"重物球质量 m_s(kg)":p.ms,
        "钢桶倾角 θ_d(°)":s["theta_d"],"锚端夹角 ϕ_a(°)":s["phi_a"],
        "浮标吃水深度 δ(m)":s["delta"],"浮标游动半径 r(m)":s["r"],
        "锚链状态":s["chain"]["regime"],"θ_d≤5°":"是" if s["theta_d"]<=5 else "否",
        "ϕ_a≤16°":"是" if s["phi_a"]<=16 else "否",
    }

def save_table(rows,name,keep):
    df=pd.DataFrame(rows); df[keep].to_csv(TABLES/name,index=False,encoding="utf-8-sig")

def chain_coords(s):
    c=s["chain"]; X=c["span"]; Z=c["drop"]
    return pd.DataFrame({"距锚点水平距离 x(m)":(X-c["x"])[::-1],
                         "距海床竖直高度 z(m)":(Z-c["z"])[::-1]})

def plot_required(base):
    rr=[]
    for d in np.linspace(0.35,1.70,280):
        try: rr.append((d,eval_state(float(d),base)["H"]))
        except BadState: pass
    hd=pd.DataFrame(rr,columns=["浮标吃水深度 δ(m)","模型海水深度 H(m)"])
    hd.to_csv(PROCESSED/"q3_H_delta_curve.csv",index=False,encoding="utf-8-sig")
    fig,ax=plt.subplots(figsize=(7.4,5.0)); ax.plot(hd.iloc[:,0],hd.iloc[:,1],lw=1.8)
    ax.set_xlabel(r"浮标吃水深度 $\delta$（m）"); ax.set_ylabel(r"模型海水深度 $H(\delta)$（m）")
    ax.set_title(r"吃水深度 $\delta$ 与海水深度 $H$ 的关系"+"\n"+
                 f"{base.chain_type}型，{base.n}节，"+rf"$m_s$={base.ms:.0f} kg，$v_w$={base.vw:.0f} m/s，$v_c$={base.vc:.1f} m/s")
    ax.grid(True,alpha=.25); fig.tight_layout()
    for ext in ("png","svg"): fig.savefig(FIGURES/f"q3_H_delta_curve.{ext}",bbox_inches="tight")
    plt.close(fig)
    states={}
    for H in (16.0,20.0):
        p=replace(base,H=H); s=solve(p); states[H]=(p,s)
        xy=chain_coords(s); xy.to_csv(PROCESSED/f"q3_chain_shape_H{H:.0f}.csv",index=False,encoding="utf-8-sig")
        fig,ax=plt.subplots(figsize=(8,5.2)); ax.plot(xy.iloc[:,0],xy.iloc[:,1],lw=2,label="锚链")
        ax.axhline(0,ls="--",lw=1,label="海床"); ax.axhline(H,ls=":",lw=1,label="海面")
        ax.set_xlabel("距锚点水平距离（m）"); ax.set_ylabel("距海床竖直高度（m）")
        ax.set_title(rf"海水深度 $H={H:.0f}$ m 时的锚链形状"+"\n"+
                     f"{p.chain_type}型，{p.n}节，"+rf"$m_s$={p.ms:.0f} kg，$v_w$={p.vw:.0f} m/s，$v_c$={p.vc:.1f} m/s")
        ax.grid(True,alpha=.25); ax.legend(); ax.set_aspect("equal",adjustable="datalim"); fig.tight_layout()
        for ext in ("png","svg"): fig.savefig(FIGURES/f"q3_chain_shape_H{H:.0f}.{ext}",bbox_inches="tight")
        plt.close(fig)
    return states

def main():
    base=P(H=18,vw=36,vc=1.5,ms=2250,chain_type="II",n=210)
    q2=solve(replace(base,vc=0.0))
    ref=(4.42514012168273,15.8259582197017,0.992846090286962,18.531479680426)
    err=max(abs(a-b) for a,b in zip((q2["theta_d"],q2["phi_a"],q2["delta"],q2["r"]),ref))
    if err>1e-8: raise RuntimeError(f"Q2回归失败:{err}")

    common=["钢桶倾角 θ_d(°)","锚端夹角 ϕ_a(°)","浮标吃水深度 δ(m)","浮标游动半径 r(m)","锚链状态","θ_d≤5°","ϕ_a≤16°"]
    rows=[]
    for t,(ell,_) in CHAIN.items():
        p=replace(base,chain_type=t,n=round(22.05/ell),vc=0.0); rows.append(row(p,solve(p)))
    save_table(rows,"q3_海水静止_锚链型号与参数.csv",["锚链型号","链环长度(mm)","单位长度质量(kg/m)"]+common)

    rows=[]
    for vc in (0,.3,.6,.9,1.2,1.5):
        for t,(ell,_) in CHAIN.items():
            p=replace(base,chain_type=t,n=round(22.05/ell),vc=vc); rows.append(row(p,solve(p)))
    save_table(rows,"q3_水流速度与锚链型号.csv",["海水流速 v_c(m/s)","锚链型号","链环长度(mm)","单位长度质量(kg/m)"]+common)

    rows=[row(replace(base,n=n),solve(replace(base,n=n))) for n in range(160,261,10)]
    save_table(rows,"q3_链条节数与参数.csv",["锚链节数 n","锚链总长度 L(m)"]+common)
    rows=[row(replace(base,vw=v),solve(replace(base,vw=v))) for v in (12.,24.,36.)]
    save_table(rows,"q3_风速与参数.csv",["海面风速 v_w(m/s)"]+common)
    rows=[row(replace(base,vc=float(v)),solve(replace(base,vc=float(v)))) for v in np.linspace(0,1.5,7)]
    save_table(rows,"q3_水流速度与参数.csv",["海水流速 v_c(m/s)"]+common)
    rows=[row(replace(base,ms=float(m)),solve(replace(base,ms=float(m)))) for m in range(1200,3601,200)]
    save_table(rows,"q3_重物球质量与参数.csv",["重物球质量 m_s(kg)"]+common)
    rows=[row(replace(base,H=float(H)),solve(replace(base,H=float(H)))) for H in np.linspace(16,20,9)]
    save_table(rows,"q3_海水深度与参数.csv",["海水深度 H(m)"]+common)

    states=plot_required(base)
    with (LOGS/"q3_engineering_checks.txt").open("w",encoding="utf-8") as f:
        f.write(f"Q2静水回归最大绝对误差={err:.3e}\n")
        f.write("工程基准：II型210节，m_s=2250kg，v_w=36m/s，v_c=1.5m/s。\n")
        for H,(p,s) in states.items():
            f.write(f"H={H:.0f}m: theta_d={s['theta_d']:.6f}deg, phi_a={s['phi_a']:.6f}deg, delta={s['delta']:.6f}m, r={s['r']:.6f}m\n")
        f.write("注意：链条等面积受流直径、F_cs球体投影面积属于工程闭合，尚未写入model card/论文。\n")
    print("Q3工程完成；Q2回归误差",err)

if __name__=="__main__": main()
