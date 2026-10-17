"""Single Tk GUI: shared SPOT 7 engine; worker never accesses Tk objects."""
import json,queue,threading,traceback
from copy import deepcopy
from pathlib import Path
import tkinter as tk
from tkinter import ttk,filedialog,messagebox
from ssa_emulator.config import default_config,plain,struct,load_config,save_config,NAMES
from ssa_emulator.pipeline import run
from ssa_emulator.examples import example_config,EXAMPLES

class SPOTApp:
    def __init__(self,root):
        self.root=root;root.title('SPOT 7 — optical SDA emulator');root.geometry('1280x850');self.events=queue.Queue();self.result=None;self.busy=False;self.canvas=None;self.figure=None
        top=ttk.Frame(root,padding=8);top.pack(fill='x');self.preset=tk.StringVar(value='ground');self.quality=tk.StringVar(value='preview')
        ttk.Combobox(top,textvariable=self.preset,values=list(EXAMPLES)+[n for n in NAMES if n not in EXAMPLES],state='readonly',width=25).pack(side='left');ttk.Combobox(top,textvariable=self.quality,values=['preview','full'],state='readonly',width=10).pack(side='left')
        ttk.Button(top,text='Load example',command=self.load_example).pack(side='left');ttk.Button(top,text='Open configuration',command=self.open_config).pack(side='left');ttk.Button(top,text='Save configuration',command=self.save_config).pack(side='left');self.runbutton=ttk.Button(top,text='Run simulation',command=self.start);self.runbutton.pack(side='left');ttk.Button(top,text='Export results',command=self.export).pack(side='left')
        pane=ttk.Panedwindow(root,orient='horizontal');pane.pack(fill='both',expand=True);left=ttk.Frame(pane,width=440);right=ttk.Frame(pane);pane.add(left,weight=1);pane.add(right,weight=2)
        tabs=ttk.Notebook(left);tabs.pack(fill='both',expand=True);self.form=ttk.Frame(tabs,padding=10);advanced=ttk.Frame(tabs);tabs.add(self.form,text='Common controls');tabs.add(advanced,text='Full configuration')
        self.variables={};self.form_values={};fields=[('Exposure (s)',('Time','Exposure_s')),('Cadence (s)',('Time','Cadence_s')),('Exposures',('Time','NumExposures')),('Image columns',('Camera','Nx')),('Image rows',('Camera','Ny')),('Horizontal FoV (deg)',('Camera','FOV_deg',0)),('Vertical FoV (deg)',('Camera','FOV_deg',1)),('Sensor CSV enabled',('Input','UseSensorCSV')),('Target CSV enabled',('Input','UseTargetCSV')),('Targets from pool',('Input','NumTargets')),('Selection: ordered/random/ids',('Input','TargetSelection','Mode')),('First row (one-based)',('Input','TargetSelection','StartRow')),('Selection seed',('Input','TargetSelection','Seed'))]
        for i,(label,path) in enumerate(fields):
            ttk.Label(self.form,text=label).grid(row=i,column=0,sticky='w',pady=3);v=tk.StringVar();ttk.Entry(self.form,textvariable=v,width=22).grid(row=i,column=1);self.variables[path]=v
        ttk.Button(self.form,text='Apply common controls to configuration',command=self.apply_controls).grid(row=len(fields),column=0,columnspan=2,pady=10)
        ttk.Label(self.form,text='Target materials and attitude: edit the selected\ntarget below. All state/pointing/instrument options\nare available in Full configuration.').grid(row=len(fields)+1,column=0,columnspan=2,sticky='w')
        self.target_index=tk.IntVar(value=1);ttk.Spinbox(self.form,from_=1,to=1000,textvariable=self.target_index,width=10).grid(row=len(fields)+2,column=0);ttk.Button(self.form,text='Edit target',command=self.edit_target).grid(row=len(fields)+2,column=1)
        instrument=ttk.Frame(tabs,padding=10);tabs.add(instrument,text='Instrument')
        optical_fields=[('Transmission',('Optics','Transmission')),('Vignetting enabled',('Optics','VignettingEnabled')),('Edge transmission',('Optics','EdgeTransmission')),('PSF centre sigma (px)',('Optics','SigmaCenter_pix')),('PSF edge addition (px)',('Optics','SigmaEdgeAdd_pix')),('Aperture diameter (m)',('Radiometry','ApertureDiameter_m')),('Quantum efficiency',('Radiometry','QuantumEfficiency')),('Read noise (electrons)',('Detector','ReadNoise_e')),('Gain (electrons/DN)',('Detector','ConversionGain_e_per_DN')),('Bias (DN)',('Detector','Bias_DN')),('Full well (electrons)',('Detector','FullWell_e')),('ADC enabled',('Detector','ADCEnabled')),('ADC bits',('Detector','ADCBits'))]
        for i,(label,path) in enumerate(optical_fields):
            ttk.Label(instrument,text=label).grid(row=i,column=0,sticky='w',pady=4);v=tk.StringVar();ttk.Entry(instrument,textvariable=v,width=20).grid(row=i,column=1);self.variables[path]=v
        motion=ttk.Frame(tabs,padding=10);tabs.add(motion,text='Motion / pointing')
        motion_fields=[('Angle to LOS (deg)',('ExampleTarget','AngleToLOS_deg')),('Relative speed (m/s)',('ExampleTarget','RelativeSpeed_mps')),('Motion direction (deg)',('ExampleTarget','Direction_deg')),('Initial angle pairs (deg)',('ExampleTarget','InitialAngles_deg')),('Range (m)',('ExampleTarget','Range_m')),('Azimuth / elevation (deg)',('ExamplePointing','AzEl_deg')),('Off nadir (deg)',('ExamplePointing','OffNadir_deg')),('Orbital azimuth (deg)',('ExamplePointing','SpaceAzimuth_deg')),('Roll (deg)',('ExamplePointing','Roll_deg'))]
        for i,(label,path) in enumerate(motion_fields):
            ttk.Label(motion,text=label).grid(row=i,column=0,sticky='w',pady=4);v=tk.StringVar();ttk.Entry(motion,textvariable=v,width=25).grid(row=i,column=1);self.variables[path]=v
        ttk.Label(motion,text='[] preserves the configured / selected CSV value.\nUse a number or JSON list. Target overrides apply\nafter CSV selection and require seeded targets.').grid(row=len(motion_fields),column=0,columnspan=2,pady=12,sticky='w')
        self.editor=tk.Text(advanced,wrap='none',undo=True);self.editor.pack(fill='both',expand=True);ttk.Button(advanced,text='Refresh common controls from JSON',command=self.refresh).pack(fill='x')
        plotbar=ttk.Frame(right);plotbar.pack(fill='x');self.product=tk.StringVar(value='detector');self.exp=tk.IntVar(value=1)
        ttk.Combobox(plotbar,textvariable=self.product,values=['detector','camera','radec','dataset','stages','photometry','binary','instance','semantic'],state='readonly').pack(side='left');ttk.Spinbox(plotbar,from_=1,to=1000,textvariable=self.exp,width=6).pack(side='left');ttk.Button(plotbar,text='Show plot',command=self.show).pack(side='left');ttk.Button(plotbar,text='Save plot',command=self.save_plot).pack(side='left');self.plotframe=ttk.Frame(right);self.plotframe.pack(fill='both',expand=True)
        self.status=tk.StringVar(value='Ready');ttk.Label(root,textvariable=self.status,padding=8).pack(fill='x');self.set_config(example_config('ground','preview'));root.after(100,self.poll)
    def set_config(self,c):
        c=deepcopy(c)
        c.setdefault('ExampleTarget',struct({k:[] for k in ['AngleToLOS_deg','RelativeSpeed_mps','Direction_deg','InitialAngles_deg','Range_m']}))
        c.setdefault('ExamplePointing',struct({k:[] for k in ['AzEl_deg','OffNadir_deg','SpaceAzimuth_deg','Roll_deg']}))
        self.editor.delete('1.0','end');self.editor.insert('1.0',json.dumps(plain(c),indent=2));self.refresh()
    def get_config(self):return struct(json.loads(self.editor.get('1.0','end')))
    def refresh(self):
        try:
            c=self.get_config()
            for path,var in self.variables.items():
                value=c
                for key in path:value=value[key]
                var.set(json.dumps(value) if not isinstance(value,str) else value)
                self.form_values[path]=var.get()
        except Exception as e:messagebox.showerror('Configuration',str(e))
    def apply_controls(self):
        c=self.get_config()
        for path,var in self.variables.items():
            if var.get()==self.form_values.get(path):continue
            target=c
            for key in path[:-1]:target=target[key]
            old=target[path[-1]];text=var.get()
            if isinstance(old,str):value=text
            else:value=json.loads(text)
            target[path[-1]]=value
        self.set_config(c);return c
    def load_example(self):self.set_config(example_config(self.preset.get(),self.quality.get()))
    def open_config(self):
        p=filedialog.askopenfilename(filetypes=[('JSON','*.json')])
        if p:
            try:self.set_config(load_config(p))
            except Exception as e:messagebox.showerror('Open configuration',str(e))
    def save_config(self):
        p=filedialog.asksaveasfilename(defaultextension='.json')
        if p:
            try:save_config(self.apply_controls(),p)
            except Exception as e:messagebox.showerror('Save configuration',str(e))
    def edit_target(self):
        try:c=self.get_config();index=self.target_index.get()-1;target=c.Targets[index]
        except Exception as e:messagebox.showerror('Target',str(e));return
        win=tk.Toplevel(self.root);win.title(f'Target {index+1}: geometry, surface and attitude');win.geometry('750x650');text=tk.Text(win,undo=True);text.pack(fill='both',expand=True);text.insert('1.0',json.dumps(plain(target),indent=2))
        def apply():
            try:
                fresh=self.get_config();fresh.Targets[index]=struct(json.loads(text.get('1.0','end')));self.set_config(fresh);win.destroy()
            except Exception as e:messagebox.showerror('Target',str(e))
        ttk.Button(win,text='Apply target',command=apply).pack()
    def start(self):
        if self.busy:return
        try:c=self.apply_controls() # Merge changed common controls into the current JSON.
        except Exception as e:messagebox.showerror('Configuration',str(e));return
        self.busy=True;self.runbutton.configure(state='disabled');self.status.set('Running; configuration snapshot captured')
        def worker():
            try:self.events.put(('result',run(c,progress=lambda s:self.events.put(('status',s)))))
            except Exception:self.events.put(('error',traceback.format_exc()))
        threading.Thread(target=worker,daemon=True).start()
    def poll(self):
        try:
            while True:
                kind,value=self.events.get_nowait()
                if kind=='status':self.status.set(value)
                else:
                    self.busy=False;self.runbutton.configure(state='normal')
                    if kind=='error':self.status.set('Run failed');messagebox.showerror('SPOT',value)
                    else:self.result=value;self.status.set('Completed — plots and export use this same acquisition');self.show()
        except queue.Empty:pass
        self.root.after(100,self.poll)
    def show(self):
        if self.result is None:return
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
        import matplotlib.pyplot as plt
        from ssa_emulator.plotting import plot_product
        try:f=plot_product(self.result,self.product.get(),self.exp.get()-1)
        except Exception as e:messagebox.showerror('Plot',str(e));return
        if self.figure:plt.close(self.figure)
        for child in self.plotframe.winfo_children():child.destroy()
        self.figure=f;self.canvas=FigureCanvasTkAgg(f,master=self.plotframe);self.canvas.draw();self.canvas.get_tk_widget().pack(fill='both',expand=True);NavigationToolbar2Tk(self.canvas,self.plotframe)
    def save_plot(self):
        if self.figure:
            p=filedialog.asksaveasfilename(defaultextension='.png',filetypes=[('PNG','*.png'),('PDF','*.pdf')])
            if p:self.figure.savefig(p,dpi=300)
    def export(self):
        if self.result is None:return
        from ssa_emulator.export import export_run
        base=filedialog.askdirectory()
        if base:
            try:
                clone=struct(dict(self.result));clone.Configuration=struct(plain(self.result.Configuration));clone.Configuration.Output.Directory=base;clone.Configuration.Output.SaveResults=True;p=export_run(clone);self.status.set('Exported: '+str(p))
            except Exception as e:messagebox.showerror('Export',str(e))

def main():
    root=tk.Tk();SPOTApp(root);root.mainloop()
if __name__=='__main__':main()
