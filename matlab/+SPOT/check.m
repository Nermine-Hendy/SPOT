function report=check
%CHECK Focused MATLAB regression checks for preset defaults and configuration controls.
% Uses small explicit test fixtures for rendering; examples stay full quality.
% Temporary files/figures are removed after checking. No benchmark is implied.
SPOT_Setup;
names={'original_defaults','catalogue_reader','empty_overrides_preserve_states', ...
    'csv_override_precedence','save_switches_and_multi_exposure','original_plot_layouts'};
tests={@defaults,@catalogue,@unchanged,@pools,@saving,@layouts};
status=strings(numel(tests),1); detail=status;
for k=1:numel(tests)
    try
        tests{k}(); status(k)="PASS";
    catch err
        status(k)="FAIL"; detail(k)=string(err.message);
    end
    fprintf('%s  %s  %s\n',status(k),names{k},detail(k));
end
report=table(string(names(:)),status,detail,'VariableNames',{'Test','Status','Detail'});
assert(~any(status=="FAIL"),'SPOT:Check','One or more checks failed; see printed details.');

    function defaults()
        presets={'ground','space_sky','space_limb','space_earth'};
        for name=presets
            C=SPOT_Config(name{1},'full');
            assert(C.Camera.Nx==1024 && C.Camera.Ny==1024);
            assert(strcmp(C.Stars.Source,'csv') && isfile(C.Stars.File));
            if strcmp(name{1},'space_limb')
                assert(C.Stars.MaxStars==200000 && C.Time.Exposure_s==1 && C.Time.Cadence_s==.8);
                assert(C.Time.StartUTC==datetime(2026,5,1,12,0,0,'TimeZone','UTC'));
                assert(isequal(C.Observer.Trajectory.Elements,[6878137 .001 97.5 40 0 0]));
                assert(strcmp(C.Radiometry.DetectorCalibration,'legacy_gain'));
                assert(strcmp(C.Seeding.Propagation,'two_body') && strcmp(C.Integration.AttitudeMode,'dynamic'));
                assert(C.Plot.DetectorColorbar);
            else
                assert(C.Stars.MaxStars==648000 && C.Time.Exposure_s==1.5 && C.Time.Cadence_s==1.5);
                assert(strcmp(C.Radiometry.DetectorCalibration,'band_photon'));
            end
            assert(strcmp(C.Frames.Model,'iau2006') && C.Detector.PoissonEnabled);
            assert(all(arrayfun(@(t) strcmp(t.Photometry.Model,'surface_brdf'),C.Targets)));
        end
        C=SPOT_Config('multi_exposure','full');
        assert(C.Time.NumExposures==3 && C.Time.Exposure_s==1.5 && C.Time.Cadence_s==1);
        C=SPOT_Config('current_ground','full');
        assert(C.Time.Exposure_s==1 && C.Time.Cadence_s==.8 && C.Stars.MaxStars==200000);
        assert(strcmp(C.Seeding.Propagation,'linear') && strcmp(C.Radiometry.DetectorCalibration,'legacy_gain'));
    end
    function catalogue()
        C=SPOT_Config('ground','full'); C.Stars.MaxStars=16;
        S=SPOT.loadStars(C);
        assert(numel(S.ID)==16 && strcmp(S.Metadata.Source,'csv'));
        assert(all(abs(vecnorm(S.Direction_GCRF)-1)<1e-12));
        assert(~any(startsWith(S.ID,'synthetic_')));
    end
    function unchanged()
        for name={'ground','space_sky','space_limb','space_earth','current_ground'}
            C=SPOT_Config(name{1},'full'); D=C;
            D.ExampleTarget=emptyTarget(); D.ExamplePointing=struct('Roll_deg',[]);
            D=SPOT.applyExampleOverrides(D);
            assert(isequaln(C.Targets,D.Targets) && isequaln(C.Pointing,D.Pointing));
        end
    end
    function pools()
        C=SPOT_Config('current_ground','full');
        C.Input.UseSensorCSV=true; C.Input.UseTargetCSV=true;
        C.Input.NumTargets=5; C.Input.TargetSelection.StartRow=10;
        C.ExampleTarget=emptyTarget(); C.ExampleTarget.AngleToLOS_deg=[90 .2 90 .2 90];
        C.ExamplePointing=struct('Roll_deg',17);
        A=SPOT.resolvePools(C); B=SPOT.applyExampleOverrides(A);
        assert(isequal(B.Input.SelectedTargetRows,(10:14)') && B.Pointing.Roll_deg==17);
        for j=1:5
            before=A.Targets(j).Seed; after=B.Targets(j).Seed;
            speed=hypot(before.TransverseSpeed_mps,before.RadialSpeed_mps);
            angle=C.ExampleTarget.AngleToLOS_deg(j);
            assert(abs(after.TransverseSpeed_mps-speed*sind(angle))<1e-8);
            assert(abs(after.RadialSpeed_mps-speed*cosd(angle))<1e-8);
        end
        C.ExampleTarget.AngleToLOS_deg=[0 90];
        rejected=false;
        try, SPOT.applyExampleOverrides(SPOT.resolvePools(C));
        catch err, rejected=strcmp(err.identifier,'SPOT:Example'); end
        assert(rejected);
    end
    function saving()
        base=tempname; mkdir(base); cleanup=onCleanup(@() removeTestFolder(base));
        for saveResults=[false true]
            for savePlots=[false true]
                C=small(); C.Output.Directory=fullfile(base,sprintf('r%d_p%d',saveResults,savePlots));
                C.Output.SaveResults=saveResults; C.Output.Plot=false; C.Plot.Save=savePlots;
                C.Time.NumExposures=2; C.Plot.Exposures='all';
                O=RUN_SPOT(C); folder=fullfile(O.Configuration.Output.Directory,C.Dataset.SampleID);
                assert(isfile(fullfile(folder,'metadata.json'))==saveResults);
                assert(isfile(fullfile(folder,'acquisition.mat'))==saveResults);
                assert(isfolder(fullfile(folder,'Figures'))==savePlots);
                assert(isempty(O.FigureHandles));
                if savePlots
                    images=dir(fullfile(folder,'Figures','*.png'));
                    assert(numel(images)==2 && numel(O.PlotFiles)==2);
                end
                if ~saveResults && ~savePlots, assert(~isfolder(C.Output.Directory)); end
            end
        end
        % Display is independent of saving. Close only this call's handles.
        C=small(); C.Output.Directory=fullfile(base,'display_only');
        C.Output.SaveResults=false; C.Output.Plot=true; C.Plot.Save=false;
        O=RUN_SPOT(C); figures=onCleanup(@() close(O.FigureHandles));
        assert(~isempty(O.FigureHandles) && ~isfolder(C.Output.Directory));
        clear figures cleanup
    end
    function layouts()
        C=small(); C.Plot.Save=false; C.Plot.Visible='off';
        O=SPOT.run(C,false);
        [f,files]=SPOT_Plots(O,1,{'camera','detector','radec','stages','dataset'});
        cleanup=onCleanup(@() close(f));
        assert(numel(f)==5 && isempty(files));
        assert(numel(findall(f(4),'Type','axes'))==24);
        assert(numel(findall(f(5),'Type','axes'))==6);
        clear cleanup
    end
end
function T=emptyTarget()
T=struct('AngleToLOS_deg',[],'RelativeSpeed_mps',[],'Direction_deg',[], ...
    'InitialAngles_deg',[],'Range_m',[]);
end
function C=small()
C=SPOT_Config('ground','full');
% Explicit changes apply only to this test fixture.
C.Frames.Model='gmst_precession'; C.Camera.Nx=16; C.Camera.Ny=12;
C.Stars.Enabled=false; C.Targets=C.Targets([]);
C.Background.Enabled=false; C.Earth.Enabled=false; C.Atmosphere.Enabled=false;
C.Detector.PoissonEnabled=false; C.Detector.ReadEnabled=false;
C.Time.Exposure_s=.02; C.Time.Cadence_s=.02; C.Time.NumExposures=1;
C.Integration.Step_s=.01; C.Output.Verbose=false; C.Output.Plot=false;
C.Plot.Camera=false; C.Plot.Detector=true; C.Plot.RADEC=false;
C.Plot.Stages=false; C.Plot.Dataset=false; C.Plot.Photometry=false;
C.Plot.ExportRawDetector=false; C.Plot.Formats={'png'}; C.Plot.DPI=72;
C.Plot.Exposures=1; C.Plot.Visible='off'; C.Dataset.SampleID='check_sample';
end
function removeTestFolder(folder)
assert(startsWith(string(folder),string(tempdir),'IgnoreCase',ispc) && ...
    ~strcmp(folder,tempdir),'SPOT:Check','Refusing cleanup outside the temporary test directory.');
if isfolder(folder), rmdir(folder,'s'); end
end
