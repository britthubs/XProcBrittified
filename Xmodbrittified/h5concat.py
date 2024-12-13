import h5py
import numpy as np

def h5_combine(files, newfile):
    ax = 0
    for index, file in enumerate(files):
        f = h5py.File(file,'r', locking=True)
        if index == 0:
            cmd = ''
            mot1 = []
            mot2 = []
            i0 = []
            i1 = []
            tm = []
            icr0 = []
            icr1 = []
            ocr0 = []
            ocr1 = []
            spectra0 = []
            spectra1 = []
        try:
            cmd += f['cmd'][()].decode('utf8')
        except AttributeError:
            cmd += f['cmd'][()]

        mot1.append(np.array(f['mot1']))
        mot1_name = str(f['mot1'].attrs["Name"])
        mot2.append(np.array(f['mot2']))
        mot2_name = str(f['mot2'].attrs["Name"])
        i0.append(np.array(f['raw/I0']))
        i1.append(np.array(f['raw/I1']))
        tm.append(np.array(f['raw/acquisition_time']))
        icr0.append(np.array(f['raw/channel00/icr']))
        icr1.append(np.array(f['raw/channel01/icr']))
        ocr0.append(np.array(f['raw/channel00/ocr']))
        ocr1.append(np.array(f['raw/channel01/ocr']))
        spectra0.append(np.array(f['raw/channel00/spectra']))
        spectra1.append(np.array(f['raw/channel01/spectra']))
        f.close()

    # add in one array
    spectra0 = np.concatenate(spectra0, axis=ax)
    icr0 = np.concatenate(icr0, axis=ax)
    ocr0 = np.concatenate(ocr0, axis=ax)
    spectra1 = np.concatenate(spectra1, axis=ax)
    icr1 = np.concatenate(icr1, axis=ax)
    ocr1 = np.concatenate(ocr1, axis=ax)
    i0 = np.concatenate(i0, axis=ax)
    i1 = np.concatenate(i1, axis=ax)
    mot1 = np.concatenate(mot1, axis=ax)
    mot2 = np.concatenate(mot2, axis=ax)
    tm = np.concatenate(tm, axis=ax)

    if len(spectra0.shape) == 3:
        sumspec0 = np.sum(spectra0, axis=(0,1))
        maxspec0 = np.zeros(sumspec0.shape[0])
        for i in range(sumspec0.shape[0]):
            maxspec0[i] = spectra0[:,:,i].max()
    else:
        sumspec0 = np.sum(spectra0, axis=0)
        maxspec0 = np.zeros(sumspec0.shape[0])
        for i in range(sumspec0.shape[0]):
            maxspec0[i] = spectra0[:,i].max()
            
    if len(spectra1.shape) == 3:
        sumspec1 = np.sum(spectra1, axis=(0,1))
        maxspec1 = np.zeros(sumspec1.shape[0])
        for i in range(sumspec1.shape[0]):
                maxspec1[i] = spectra1[:,:,i].max()
    else:
        sumspec1 = np.sum(spectra1, axis=0)
        maxspec1 = np.zeros(sumspec1.shape[0])
        for i in range(sumspec1.shape[0]):
            maxspec1[i] = spectra1[:,i].max()
        
    # write the new file
    f = h5py.File(newfile, 'w', locking=True)
    f.create_dataset('cmd', data=cmd)
    dset = f.create_dataset('mot1', data=mot1, compression='gzip', compression_opts=4)
    dset.attrs['Name'] = mot1_name
    dset = f.create_dataset('mot2', data=mot2, compression='gzip', compression_opts=4)
    dset.attrs['Name'] = mot2_name
    f.create_dataset('raw/acquisition_time', data=tm, compression='gzip', compression_opts=4)
    f.create_dataset('raw/channel00/icr', data=icr0, compression='gzip', compression_opts=4)
    f.create_dataset('raw/channel00/maxspec', data=maxspec0, compression='gzip', compression_opts=4)
    f.create_dataset('raw/channel00/ocr', data=ocr0, compression='gzip', compression_opts=4)
    f.create_dataset('raw/channel00/spectra', data=spectra0, compression='gzip', compression_opts=4)
    f.create_dataset('raw/channel00/sumspec', data=sumspec0, compression='gzip', compression_opts=4)
    f.create_dataset('raw/channel01/icr', data=icr1, compression='gzip', compression_opts=4)
    f.create_dataset('raw/channel01/maxspec', data=maxspec1, compression='gzip', compression_opts=4)
    f.create_dataset('raw/channel01/ocr', data=ocr1, compression='gzip', compression_opts=4)
    f.create_dataset('raw/channel01/spectra', data=spectra1, compression='gzip', compression_opts=4)    
    f.create_dataset('raw/channel01/sumspec', data=sumspec1, compression='gzip', compression_opts=4)
    f.create_dataset('raw/I0', data=i0, compression='gzip', compression_opts=4)
    f.create_dataset('raw/I1', data=i1, compression='gzip', compression_opts=4)
    f.close()                   
    
    print("FINISHED")
    
    
path1 = '/Users/burrito/Library/CloudStorage/OneDrive-UGent/Thesis/scripts/or1_90_0005_scan1.h5'
path2 = '/Users/burrito/Library/CloudStorage/OneDrive-UGent/Thesis/scripts/or1_90_0006_scan1.h5'
newpath = '/Users/burrito/Library/CloudStorage/OneDrive-UGent/Thesis/scripts/or1_90_56comb_scan1.h5'
h5_combine([path1,path2], newpath)