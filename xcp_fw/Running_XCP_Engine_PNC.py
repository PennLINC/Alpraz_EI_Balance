#!/usr/bin/env python
# coding: utf-8


import flywheel
import datetime
import pytz

fw = flywheel.Client()
proj = fw.lookup('bbl/PNC_CS_810336')
sessions = proj.sessions()
# ids_to_get = ['8469']
# sessions = [fw.get(x.id) for x in sessions if x.label in ids_to_get]
sessions = [fw.get(x.id) for x in sessions]


def get_latest_fmriprep(session):
    
    if session.analyses:
        
        timezone = pytz.timezone("UTC")
        init_date = datetime.datetime(2018, 1, 1)
        latest_date = timezone.localize(init_date)

        latest_run = None

        for i in session.analyses:
            gear_name = i.gear_info['name']
            state = i.job.state
            date = i.created
            if 'fmriprep' in gear_name and date > latest_date and state =='complete':
                latest_date = date
                latest_run = i
        
        if latest_run is not None:
            fmriprep_out = [x for x in latest_run.files if 'fmriprep' in x.name and 'html' not in x.name and x.name.endswith(".zip")][0]
            fmriprep_out
            return(fmriprep_out)

        else:
            return None

    else:
        return None


xcp = fw.lookup('gears/xcpengine-fw')

# Get task file and design file
designFile = [x for x in proj.files if "task_pnc_idemo_acompcor_GSR.dsn" in x.name][0]
designFile = proj.get_file(designFile.name)

taskFile = [x for x in proj.files if "task.zip" in x.name][0]


xcp_runs = {}

for i, x in enumerate(sessions):
    ses = fw.get(x.id)

    #struc = get_latest_struct(ses)
    fmriprep = get_latest_fmriprep(ses)
    
    if  fmriprep:
        myconfig = {
            'analysis_type': 'acompcor_gsr',
            'task_name': 'idemo',
            'session': 'PNC1'
        }

        myinput = {
            'taskfile': taskFile,
            'fmriprepdir': fmriprep,
            'designfile': designFile
        }
        jobid = xcp.run(analysis_label="IDEMO_ACOMPCOR_GSR{}".format(datetime.datetime.now()), destination=ses, inputs=myinput, config=myconfig)
        xcp_runs[ses.label] = jobid
    else:
        xcp_runs[ses.label] = None


len(xcp_runs)
