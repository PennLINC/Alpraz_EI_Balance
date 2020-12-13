#!/usr/bin/env python
# coding: utf-8

import flywheel
import datetime
import pytz

# Get project and sessions
fw = flywheel.Client()
proj = fw.lookup('bbl/PNC_CS_810336')
sessions = proj.sessions()
sessions = [fw.get(x.id) for x in sessions]

# Get analysis gear (fmriprep)
fmriprep = fw.lookup('gears/fmriprep-hpc')

# Analysis name
now = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M")
analysis_label = 'IDEMO_fMRIPrep_{}_{}_{}'.format(now, fmriprep.gear.name, fmriprep.gear.version)

analysis_ids = []
fails = []
#Launch gear on identified sessions
for i, x in enumerate(sessions):
    ses = fw.get(x.id)
    
    analysis_input = {'freesurfer_license':  proj.files[5]}
    analysis_config = {
        'force_syn': False,
        'bold2t1w_dof': 9,
        'intermediate_folders': '',
        'sge-cpu': '8',
        'skull_strip_fixed_seed': False,
        'use_syn_sdc': False,
        'save_intermediate_work': False,
        'singularity-writable': False,
        'no_submm_recon': False,
        'sge-short': False,
        'longitudinal': False,
        'fs_no_reconall': False,
        'template': 'MNI152NLin2009cAsym',
        'anat_only': False,
        'force_bbr': False,
        'output_space': 'template fsaverage5',
        'template_resampling_grid': 'native',
        'intermediate_files': '',
        't2s_coreg': False, 
        'medial_surface_nan': False,
        'save_outputs': False,
        'sloppy_mode': False,
        'fmap_no_demean': False,
        'low_mem': False,
        'singularity-debug': False,
        'ignore': '',
        'skull_strip_template': 'OASIS30ANTs',
        'cifti_output': '91k',
        'use_aroma': False,
        'use_all_sessions': False,
        'fmap_bspline': False,
        'force_no_bbr': False,
        'sge-ram': '64G'
    }
    try:
        jobid = fmriprep.run(analysis_label=analysis_label, destination=ses, inputs=analysis_input, config=analysis_config)
        analysis_ids.append(ses.label)
        print(ses.label)
    except Exception as e:
        print(e)
        fails.append(ses.label)


#Write analysis IDs and failed sessions to files
with open('{}_{}_{}_analysisIDS.txt'.format(fmriprep.gear.name,fmriprep.gear.version,now), 'w') as f: 
    for id in analysis_ids:
        f.write("%s\n" % id)

with open('{}_{}_{}_failSES.txt'.format(fmriprep.gear.name,fmriprep.gear.version,now), 'w') as a:
    for ses in fails:
        a.write("%s\n" % ses) 
