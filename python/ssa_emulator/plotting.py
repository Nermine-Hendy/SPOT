from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def overlay(ax,E,c,angular=False):
    for j,T in enumerate(E.Targets):
        h=T.Samples;ok=h.Visible&(h.RadiometricRate>0);col=c.Plot.Colors[j%len(c.Plot.Colors)]
        x=h.Xi_deg if angular else h.x;y=h.Eta_deg if angular else h.y
        ax.plot(np.where(ok,x,np.nan),np.where(ok,y,np.nan),color=col,lw=1.5)
        if ok.any():
            i=np.flatnonzero(ok)[-1]
            if T.Appearance=='object':ax.plot(x[i],y[i],'o',color=col,ms=3)
            ax.text(x[i],y[i],f'T{T.ID}',color=col,clip_on=True)

def plot_product(out,kind='detector',exposure=0):
    E=out.Exposures[exposure];c=out.Configuration
    if kind=='dataset':
        fig,axes=plt.subplots(2,3,figsize=(13,8));kinds=['detector','binary','instance','semantic','boxes','metadata']
        for ax,name in zip(axes.flat,kinds):draw(ax,out,name,exposure)
    elif kind=='stages':
        t=E.Truth;o=E.OpticalStages;q=E.DetectorStages
        blocked=t.Active&(t.Elevation_deg<c.Visibility.MinElevation_deg) if t.GeometryRegime=='ground_sky' else t.EarthMask
        sky=E.SkyRate/np.maximum(t.SolidAngle_sr,np.finfo(float).tiny);earth=E.EarthRate/np.maximum(t.SolidAngle_sr,np.finfo(float).tiny)
        bias=c.Detector.Bias_DN if c.Detector.BiasEnabled else 0
        titles=['In-FoV mask','Earth / horizon blocked','Visible sky mask','Earth-in-FoV mask','Stars flux only','Sky background radiance','Earth radiance only','Targets flux only','Point-source flux input','Diffuse radiance input','Radiometric input to optics','Vignetting map','After vignetting','After optics / PSF','Electron mean image','Dark current mean','Read noise realization','Pre-ADC signal','Bias/dark-subtracted image','Final detector image','After shot noise','After dark current','Saturation mask','Detector truth map']
        data=[t.Active,blocked,t.SkyMask,t.EarthMask,E.StarRate,sky,earth,E.TargetRate,E.StarRate+E.TargetRate,sky+earth,E.SceneRate,o.Vignetting,o.Vignetted,o.PSF,q.ExpectedSignal_e,q.ExpectedDark_e,q.ReadNoise_e,q.Analog_DN,E.Image_DN-bias-q.ExpectedDark_e/c.Detector.ConversionGain_e_per_DN,E.Image_DN,q.SignalElectrons,q.BeforeFullWell_e,q.FullWellSaturated|q.ADCSaturated,np.zeros_like(E.Image_DN)]
        for i in list(range(4,11))+[12,13]:
            a=data[i]
            if np.any(a>0):data[i]=np.log10(1+np.maximum(a,0)/max(a.max()*1e-5,np.finfo(float).tiny))
        products=list(zip(titles,data))
        fig,axes=plt.subplots(4,6,figsize=(19,12))
        for i,(ax,(name,a)) in enumerate(zip(axes.flat,products)):
            limits=dict(vmin=0,vmax=1) if i in (0,1,2,3,22) else {}
            ax.imshow(a,cmap='gray',origin='upper',extent=[.5,c.Camera.Nx+.5,c.Camera.Ny+.5,.5],**limits);ax.set_title(name,fontsize=9);ax.set_xticks([]);ax.set_yticks([])
        overlay(axes.flat[-1],E,c)
    elif kind=='photometry':
        fig,axes=plt.subplots(2,1,figsize=(11,7))
        for j,T in enumerate(E.Targets):
            h=T.Samples;t=h.Time_s-h.Time_s[0];col=c.Plot.Colors[j%len(c.Plot.Colors)]
            if 'Reflection' in h:
                axes[0].plot(t,h.Reflection.DiffuseIrradiance_W_m2,color=col,label=T.Name+' diffuse');axes[0].plot(t,h.Reflection.SpecularIrradiance_W_m2,'--',color=col,label=T.Name+' specular')
            axes[1].plot(t,h.ReceivedIrradiance_W_m2*h.Visible,color=col,label=T.Name)
        for ax in axes:ax.legend(fontsize=8);ax.set_xlabel('Time within exposure (s)');ax.set_ylabel('Irradiance (W m$^{-2}$)');ax.grid(alpha=.3)
    else:
        fig,ax=plt.subplots(figsize=(8,7));draw(ax,out,kind,exposure)
        if kind=='detector' and c.Plot.get('DetectorColorbar',False):
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            cax=make_axes_locatable(ax).append_axes('right',size='2.5%',pad=.1)
            fig.colorbar(ax.images[0],cax=cax,label='Digital number (DN)')
            ax.set_title(f'Detector image plane | {E.Truth.GeometryRegime} | exposure {E.Index}')
    fig.tight_layout();return fig

def draw(ax,out,kind,k=0):
    from matplotlib.patches import Rectangle
    E=out.Exposures[k];c=out.Configuration;extent=[.5,c.Camera.Nx+.5,c.Camera.Ny+.5,.5]
    if kind=='metadata':
        lines=[f'Exposure {E.Index}; duration {c.Time.Exposure_s:g} s',f'Scene: {E.Truth.GeometryRegime}',f'Image: {c.Camera.Nx} × {c.Camera.Ny}',f'Annotated targets: {sum(t.HasAnnotation for t in E.Targets)} / {len(E.Targets)}','', 'ID   appearance   smear / threshold (px)']
        lines += [f'{t.ID:<4} {t.Appearance:<11} {t.Morphology.EquivalentSmear_pix:.3g} / {t.Morphology.Threshold_pix:.3g}' for t in E.Targets]
        lines += ['', 'One exposure-integrated renderer for all targets.']
        ax.axis('off');ax.text(0,1,'\n'.join(lines),transform=ax.transAxes,va='top',fontsize=9,family='monospace');ax.set_title('Label and metadata summary');return
    if kind=='camera':
        f=c.Camera.FOV_deg;extent=[-f[0]/2,f[0]/2,-f[1]/2,f[1]/2]
        ax.imshow((.07*E.Truth.Active+.13*E.Truth.EarthMask)[::-1],cmap='gray',origin='lower',extent=extent,vmin=0,vmax=1)
        s=E.StarTruth;ax.scatter(s.Xi_deg,s.Eta_deg,s=2,c='.7');overlay(ax,E,c,True);ax.set_xlabel('Camera angle ξ (deg)');ax.set_ylabel('Camera angle η (deg)')
        if c.Camera.Shape=='circular':
            th=np.linspace(0,2*np.pi,361);q=np.deg2rad(f[0]/2);ax.plot(np.rad2deg(np.arctan(np.tan(q)*np.cos(th))),np.rad2deg(np.arctan(np.tan(q)*np.sin(th))),'g-',lw=2)
        else:ax.plot(np.array([-1,1,1,-1,-1])*f[0]/2,np.array([-1,-1,1,1,-1])*f[1]/2,'g-',lw=2)
        ax.plot(0,0,'wx',ms=10,mew=2)
    elif kind=='radec':
        middle=np.argmin(abs(E.TargetTime_s-(E.TargetTime_s[0]+c.Time.Exposure_s/2)))
        b=E.CameraToGCRF[middle,:,2];centre=np.degrees(np.arctan2(b[1],b[0]))%360;wrap=lambda r:(r-centre+180)%360-180
        ra=wrap(E.Truth.RA_deg);dec=E.Truth.DEC_deg
        ax.set_facecolor('black');ax.scatter(ra[::8,::8],dec[::8,::8],c=.08+.13*E.Truth.EarthMask[::8,::8],s=2,cmap='gray',vmin=0,vmax=1)
        s=E.StarTruth;ax.scatter(wrap(s.RA_deg),s.DEC_deg,s=8,c='.8',alpha=.35)
        for j,t in enumerate(E.Targets):
            h=t.Samples;ok=h.Visible&(h.RadiometricRate>0);x=wrap(h.RA_deg).copy();y=h.DEC_deg.copy();seam=np.r_[False,abs(np.diff(x))>180];x[~ok|seam]=np.nan;y[~ok|seam]=np.nan;col=c.Plot.Colors[j%len(c.Plot.Colors)]
            ax.plot(x,y,color=col,lw=1.5)
            if ok.any():
                last=np.flatnonzero(ok)[-1];ax.text(wrap(h.RA_deg[last]),h.DEC_deg[last],f'T{t.ID}',color=col)
        ax.plot(0,np.rad2deg(np.arcsin(b[2])),'wx',ms=10,mew=2);ax.set_xlim(-180,180);ax.set_ylim(-90,90)
        ax.set_xlabel(f'RA offset from {centre:.2f}°');ax.set_ylabel('DEC (deg)')
    else:
        data={'binary':E.Labels.Binary,'instance':E.Labels.Instance,'semantic':E.Labels.Semantic,'target':E.TargetRate}.get(kind,E.Image_DN);kwargs={}
        if kind in ('detector','boxes'):lo,hi=np.percentile(data,c.Plot.DisplayPercentiles);kwargs=dict(vmin=lo,vmax=max(hi,lo+1e-9))
        from matplotlib.colors import ListedColormap,BoundaryNorm
        cmap='gray'
        if kind=='instance':
            im=np.zeros_like(data);colors=['black']
            for j,T in enumerate(E.Targets):im[data==T.ID]=j+1;colors.append(c.Plot.Colors[j%len(c.Plot.Colors)])
            data=im;cmap=ListedColormap(colors if len(colors)>1 else ['black','white']);kwargs=dict(vmin=0,vmax=max(1,len(E.Targets)))
        elif kind=='semantic':cmap=ListedColormap(['black','#8c8c8c','white']);kwargs=dict(vmin=0,vmax=2)
        ax.imshow(data,cmap=cmap,extent=extent,**kwargs)
        if kind in ('detector','boxes'):overlay(ax,E,c)
        if kind=='boxes':
            for T in E.Targets:
                if T.HasAnnotation:x,y,w,h=T.BoundingBox_xywh;ax.add_patch(Rectangle((x-.5,y-.5),w,h,fill=False,edgecolor='yellow',lw=1))
        ax.set_xlabel('x pixel (one-based)');ax.set_ylabel('y pixel (one-based)')
    ax.set_title(f'{kind.replace("_"," ").title()} | {E.Truth.GeometryRegime}')
