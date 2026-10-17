"""Controller/worker tests do not require a native window server."""
import queue
import time
from gui import SPOTApp
from ssa_emulator.examples import example_config


class Widget:
    def __init__(self):self.value=None
    def configure(self,**kwargs):self.value=kwargs
    def set(self,value):self.value=value


def controller(c):
    app=object.__new__(SPOTApp)
    app.busy=False;app.result=None;app.events=queue.Queue()
    app.runbutton=Widget();app.status=Widget();app.apply_controls=lambda:c
    return app


def test_gui_worker_returns_real_acquisition():
    c=example_config('ground','preview');c.Camera.Nx=c.Camera.Ny=16
    c.Time.Exposure_s=.02;c.Stars.Enabled=False
    app=controller(c);app.start();assert app.busy
    deadline=time.monotonic()+30
    while time.monotonic()<deadline:
        kind,value=app.events.get(timeout=30)
        if kind=='status':continue
        assert kind=='result',value
        assert value.Exposures[0].Image_DN.shape==(16,16)
        assert value.Configuration.Time.Exposure_s==.02
        break
    else:raise AssertionError('GUI worker timed out')


def test_gui_worker_reports_model_errors():
    c=example_config('ground','preview');c.Time.Exposure_s=-1
    app=controller(c);app.start()
    kind,value=app.events.get(timeout=30)
    assert kind=='error' and 'positive finite' in value
