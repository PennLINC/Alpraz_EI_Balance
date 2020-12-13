#!/usr/bin/env python
# coding: utf-8


import flywheel
import datetime
import pytz

fw = flywheel.Client()
proj = fw.lookup('bbl/ALPRAZ_805556')
sessions = proj.sessions()
ids_to_get = ['001571']
sessions = [fw.get(x.id) for x in sessions if x.label in ids_to_get]


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
            fmriprep_out = [x for x in latest_run.files if 'fmriprep' in x.name][0]
            fmriprep_out
            return(fmriprep_out)

        else:
            return None

    else:
        return None



#get task file 
def get_taskfile(session):
    for i in range(len(session.files)):
        if session.files[i].name.endswith('task.zip'):
            taskfile=session.files[i]
        else:
            taskfile=None
    return taskfile

xcp = fw.lookup('gears/xcpengine-fw')


designFile = [x for x in proj.files if "task_acompcor_GSR_alff_reho_fcon" in x.name][0]
designFile = proj.get_file(designFile.name)


xcp_runs = {}

for i, x in enumerate(sessions):
    ses = fw.get(x.id)

    #struc = get_latest_struct(ses)
    fmriprep = get_latest_fmriprep(ses)
    taskfile1 = get_taskfile(ses)
    
    if  fmriprep and taskfile1 :
        myconfig = {
            'analysis_type': 'task_acompcor_GSR_alff_reho_fcon',
            'task_name' :'emotionid'
        }

        myinput = {
            'fmriprepdir': fmriprep,
            'designfile': designFile,
            'taskfile': taskfile1,
        }
        jobid = xcp.run(analysis_label="XCP_SDK_TASK_GSR{}".format(datetime.datetime.now()), destination=ses, inputs=myinput, config=myconfig)
        xcp_runs[ses.label] = jobid
    else:
        xcp_runs[ses.label] = None


len(xcp_runs)
