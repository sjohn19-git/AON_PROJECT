import os
import glob
import shutil

filed=[]
files=[]
for ele in glob.glob("//scratch/sjohn/spectrogram/*/*.npy"):
    filed.append(os.path.join("/home/sjohn/Data",ele.split("/")[-2]))
    files.append(ele)
for ele in glob.glob("//scratch/sjohn/spectrogram/*/*.db"):
    filed.append(os.path.join("/home/sjohn/Data",ele.split("/")[-2]))
    files.append(ele)




for i in range(len(files)):
    isExist = os.path.exists(filed[i])
    if not isExist:
        os.makedirs(filed[i])
        
    shutil.copy(files[i],filed[i])


    
isExist
