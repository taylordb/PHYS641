import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
plt.ion()

pathname = 'LOSC_Event_tutorial/'

import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
plt.ion()

def read_template(filename):
    dataFile=h5py.File(filename,'r')
    template=dataFile['template']
    tp=template[0]
    tx=template[1]
    return tp,tx
def read_file(filename):
    dataFile=h5py.File(filename,'r')
    dqInfo = dataFile['quality']['simple']
    qmask=dqInfo['DQmask'][...]

    meta=dataFile['meta']
    #gpsStart=meta['GPSstart'].value
    gpsStart=meta['GPSstart'][()]
    #print meta.keys()
    #utc=meta['UTCstart'].value
    utc=meta['UTCstart'][()]
    #duration=meta['Duration'].value
    duration=meta['Duration'][()]
    #strain=dataFile['strain']['Strain'].value
    strain=dataFile['strain']['Strain'][()]
    dt=(1.0*duration)/len(strain)

    dataFile.close()
    return strain,dt,utc



#fnames=glob.glob("[HL]-*.hdf5")
#fname=fnames[0]

# first event

# Hanford detector:
fname='H-H1_LOSC_4_V2-1126259446-32.hdf5'
print('reading file ',pathname+fname)
strain,dt,utc=read_file(pathname+fname)
np.savez('H-H1_LOSC_4_V2-1126259446-32.npz', strain=strain, dt=dt,utc=utc, allow_pickle=True)
# template:
template_name='GW150914_4_template.hdf5'
tp,tx=read_template(pathname+template_name)
np.savez('GW150914_4_template.npz', tp=tp, tx=tx, allow_pickle=True)
# Livingston detector:
fname='L-L1_LOSC_4_V2-1126259446-32.hdf5'
print('reading file ',pathname+fname)
strain,dt,utc=read_file(pathname+fname)
np.savez('L-L1_LOSC_4_V2-1126259446-32.npz', strain=strain, dt=dt,utc=utc, allow_pickle=True)

# second event

# Hanford detector:
fname='H-H1_LOSC_4_V2-1128678884-32.hdf5'
print('reading file ',pathname+fname)
strain,dt,utc=read_file(pathname+fname)
np.savez('H-H1_LOSC_4_V2-1128678884-32.npz', strain=strain, dt=dt,utc=utc, allow_pickle=True)
# template: 
template_name='LVT151012_4_template.hdf5'
tp,tx=read_template(pathname+template_name)
np.savez('LVT151012_4_template.npz', tp=tp, tx=tx, allow_pickle=True)
# Livingston detector:
fname='L-L1_LOSC_4_V2-1128678884-32.hdf5'
print('reading file ',pathname+fname)
strain,dt,utc=read_file(pathname+fname)
np.savez('L-L1_LOSC_4_V2-1128678884-32.npz', strain=strain, dt=dt,utc=utc, allow_pickle=True)

# third event

# Hanford detector:
fname='H-H1_LOSC_4_V2-1135136334-32.hdf5'
print('reading file ',pathname+fname)
strain,dt,utc=read_file(pathname+fname)
np.savez('H-H1_LOSC_4_V2-1135136334-32.npz', strain=strain, dt=dt,utc=utc, allow_pickle=True)
# template: 
template_name='GW151226_4_template.hdf5'
tp,tx=read_template(pathname+template_name)
np.savez('GW151226_4_template.npz', tp=tp, tx=tx, allow_pickle=True)
# Livingston detector:
fname='L-L1_LOSC_4_V2-1135136334-32.hdf5'
print('reading file ',pathname+fname)
strain,dt,utc=read_file(pathname+fname)
np.savez('L-L1_LOSC_4_V2-1135136334-32.npz', strain=strain, dt=dt,utc=utc, allow_pickle=True)

# fourth event

# Hanford detector:
fname='H-H1_LOSC_4_V1-1167559920-32.hdf5'
print('reading file ',pathname+fname)
strain,dt,utc=read_file(pathname+fname)
np.savez('H-H1_LOSC_4_V1-1167559920-32.npz', strain=strain, dt=dt,utc=utc, allow_pickle=True)
# template: 
template_name='GW170104_4_template.hdf5'
tp,tx=read_template(pathname+template_name)
np.savez('GW170104_4_template.npz', tp=tp, tx=tx, allow_pickle=True)
# Livingston detector:
fname='L-L1_LOSC_4_V1-1167559920-32.hdf5'
print('reading file ',pathname+fname)
strain,dt,utc=read_file(pathname+fname)
np.savez('L-L1_LOSC_4_V1-1167559920-32.npz', strain=strain, dt=dt,utc=utc, allow_pickle=True)


#np.savez('GW150914_4_template.npz', tp=tp, tx=tx, allow_pickle=True)