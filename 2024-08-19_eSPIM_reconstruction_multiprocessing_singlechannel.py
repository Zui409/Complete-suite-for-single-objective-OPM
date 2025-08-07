# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:19:54 2024

@author: zz409
"""
from skimage import io
import numpy as np
from scipy import ndimage
from scipy.ndimage import gaussian_filter1d
import tifffile
import multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
from functools import partial
import os 
import re
import imagej

# Main directory of the data where your positions in (eg. folder where the Pos0, Pos1 belong)
path = r'D:\Example\cell_image'

# Where your fiji is. For cropping & projection 
path_fiji = r'D:\Example\Fiji\Fiji.app'              

# Have you done these processing steps?
renamed = False                    # Rename and reorganise the file
cropped = False                    # Fiji cropping + subtract background + reverse
processed = False                  # Deskew
resliced = False                   # Reslice
projectioned = False               # Fiji 3D projection

# crop image
roi = "134,558,672,94"   

# Rolling ball background subtraction
sub_background = "50"

# Parameter for deskewing
angle = 35
spacing_factor = 2.4           # no binning 6.5

# Parameter for 3D projection:
direction = "Y"       # X or Y. Must be CAPITAL letter
angle_increment = 10

# Normally no need to change
# For Multi/Parallel Processing
param = [angle, spacing_factor]
num_cpu = multiprocessing.cpu_count()-2          
num_workers = num_cpu

#%% Rename

# Find all directories of the positions in the data_path folder in a list
def findPos(path):
    pos_path_list = [f.path for f in os.scandir(path) if f.is_dir()]
    return pos_path_list

# Find all directories of the slices in a list
def generate_slice(pos_path):
    cycle_path_list = [f.path for f in os.scandir(pos_path) if f.is_dir()]
    
    return cycle_path_list

# Find all directories of the TIFF files in a list
def getListFiles(path):
    filelist = [] 
    for root, dirs, files in os.walk(path):  
        for filespath in files: 
            if filespath[-4:] == '.tif':
                filelist.append(os.path.join(root,filespath)) 
    return filelist

# Rename and save images
def concat_img(filelist,con_name):
    num_imgs = len(filelist)
    for i in range(0, num_imgs):
        im = io.imread(filelist[i])
        tifffile.imsave(con_name, im)

def rename_img(main_path, img_name, i):
    im_newname = main_path + r'\\' + str(i) + '.tif'
    if not os.path.exists(im_newname):
        os.rename(img_name[0], im_newname)
    return im_newname
        
def sort_order(a_list):
    index_list = []
    for element in a_list:
        index = re.findall(r'\d+', element)
        index_list.append(int(index[-1]))
    sorted_list = [x for _,x in sorted(zip(index_list, a_list))] 
    return sorted_list

def rename(path):
    pos_path_list = findPos(path)
    for main_path in pos_path_list:
        slice_list = generate_slice(main_path)
        for i in range(0,len(slice_list)):
            #con_name = main_path + '/' + str(i) +'.tif'
            img_name = getListFiles(slice_list[i])
    return 0

def get_img_list(path):
    pos_path_list = findPos(path)
    pos_name_list = []
    for main_path in pos_path_list:
        img_name = getListFiles(main_path)
        sorted_name_list = sort_order(img_name)
        pos_name_list.append(sorted_name_list)
    return pos_name_list

    
#%% Create analysis folder

# Create sub folder(s) with a given name
def create_sub_folder(path, new_folder_name):
    if isinstance(path, list):
        new_dir_list = []
        for element in path:
            new_dir_name = element + '/' + new_folder_name
            new_dir_list.append(new_dir_name)
            if not os.path.exists(new_dir_name):
                os.makedirs(new_dir_name)
    elif isinstance(path, str):
        new_dir_name = path + '\\' + new_folder_name
        if not os.path.exists(new_dir_name):
            os.makedirs(new_dir_name)
        new_dir_list = [new_dir_name]
    else:
        print("Error: Invalid input path.")
    return new_dir_list

def get_last_name(path):
    if isinstance(path, list):
        name_list = []
        for path_k in path:
            name_k = path_k.split("\\")[-1]
            name_list.append(name_k)
    elif isinstance(path, str):
        name_list = [path.split("\\")[-1]]
    return name_list

# Get a list of direct subfolder of the path folder
def get_sub_folder(path):
    sub_folder_list = [f.path for f in os.scandir(path) if f.is_dir()]
    return sub_folder_list

def create_analysis_folder(path):
    pos_list = get_sub_folder(path)
    pos_name = get_last_name(pos_list)
    ana_pos_list = []
    for i in range(len(pos_list)):
        new_analysis_path = path + "_Analysis" + "/" + pos_name[i]
        ana_pos_list.append(new_analysis_path)
        if not os.path.exists(new_analysis_path):
            os.makedirs(new_analysis_path)
    folder_sb = create_sub_folder(ana_pos_list, 'Sub_bg')
    folder_process = create_sub_folder(ana_pos_list, 'Processed')
    folder_reslice = create_sub_folder(ana_pos_list, 'Resliced')
    folder_projection = create_sub_folder(ana_pos_list, 'Projection')
    
    return ana_pos_list, folder_sb, folder_process, folder_reslice, folder_projection

#%% Crop images in Fiji

def crop_images(path_fiji, pos_name_list, roi, sub_background, folder_sb):
    IJ = imagej.init(path_fiji, mode = 'interactive')
    IJ.ui().showUI()
    for i in range(len(pos_name_list)):
        for j in range(len(pos_name_list[i])):
            image = pos_name_list[i][j]
            image = image.replace("\\", "/")
            #imp =  IJ.io().open(image)     
            macro_new1 = '''
            open("''' +image + '''");
            
            
            '''
            IJ.py.run_macro(macro_new1)
            process_img = folder_sb[i] + "/" + str(j) + ".tif"
            process_img = process_img.replace("\\", "/")
            macro_new = '''
            selectWindow("''' +image.split('/')[-1] +'''");
            makeRectangle(''' + roi + ''');
            run("Duplicate...", "duplicate");      
            run("Subtract Background...", "rolling='''+ sub_background + ''' stack");
            run("Reverse");
            saveAs("Tiff", "'''+ process_img + '''");
            close();
            '''
            IJ.py.run_macro(macro_new)

            macro_new2 = '''
            close();
            '''
            IJ.py.run_macro(macro_new2)
    #del IJ
    return IJ

def fake_ij(cropped):
    if cropped == True:
        IJ = "No"
    return IJ
            
#%% Read cropped images

def read_cropped_images(ana_pos_list):    
    SB_img_list = [] 
    for main_path in ana_pos_list:
        folder_img = main_path + r'/Sub_bg'
        img_name = getListFiles(folder_img)
        sorted_name_list = sort_order(img_name)
        SB_img_list.append(sorted_name_list)
    return SB_img_list

def sort_order_third(a_list):
    index_list = []
    for element in a_list:
        index = re.findall(r'\d+', element)
        index_list.append(int(index[-4]))
    sorted_list = [x for _,x in sorted(zip(index_list, a_list))] 
    return sorted_list

def read_processed_images(ana_pos_list):    
    process_img_list = [] 
    process_img_last_list = []
    for main_path in ana_pos_list:
        folder_img = main_path + r'/Processed'
        img_name = getListFiles(folder_img)
        sorted_name_list = sort_order_third(img_name)
        last_img = sorted_name_list[0]
        process_img_last_list.append(last_img)
        process_img_list.append(sorted_name_list)
    return process_img_list, process_img_last_list

def read_projection_images(ana_pos_list):    
    projection_img_list = [] 
    for main_path in ana_pos_list:
        folder_img = main_path + r'/Projection'
        img_name = getListFiles(folder_img)
        sorted_name_list = sort_order(img_name)
        projection_img_list.append(sorted_name_list)
    return projection_img_list


def get_intensity_range(process_img_list, process_img_last_list):
    im_range = []
    for imname in process_img_last_list:
        img = io.imread(imname)
        im_max = int(np.max(img))
        im_min = np.min(img)
        range_pair = "(" + str(im_min) +", "+ str(im_max) + ")"
        im_range.append(range_pair)
    return im_range
            
            
#%%
def get_all_sb_image_list(SB_img_list):
    total_img_list = []
    for i in range(len(SB_img_list)):
        total_img_list += SB_img_list[i]
    return total_img_list

def get_save_path(img_list):
    name_list = get_last_name(img_list)
    path_name_list = []
    new_save_path_list = []
    for path_k in img_list:
        name_k = path_k.split("\\")[-1]
        length =len(name_k) +1
        name_folder = path_k[:-length]
        new_path_name = name_folder + '_Processed'
        path_name_list.append(new_path_name)
    for i in range(len(name_list)):
        save_path = path_name_list[i] + "\\" +  name_list[i]
        new_save_path_list.append(save_path)
    return new_save_path_list

#%% Multiprocessing Deskew

def process_img(imname, param):
    angle = param[0]
    spacing_factor = param[1]
    
    #imname = path+'\\'+str(img_index)+".tif"
    name_k = imname.split("\\")[-1]
    length =len(name_k) +1
    name_folder = imname[:-length]
    new_path_name = name_folder[:-7]
    save_path = new_path_name + "\\Processed\\" + name_k
    
    if os.path.exists(imname):
        img0 = io.imread(imname)
        shear_factor = 1/np.tan(np.radians(angle))
        img = ndimage.zoom(img0, (spacing_factor, np.sin(np.radians(angle)), 1))
        
        (frame, x, y) = img.shape
        transform_matrix = [[1, shear_factor, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
        sheared_img = ndimage.affine_transform(img, transform_matrix, offset = (-x*shear_factor, 0, 0), output_shape = (int(frame+x*shear_factor), x, y))

        sheared_imname = save_path[:-4]+'_angle_'+str(angle)+'_sf_'+str(spacing_factor)+'.tif'
        tifffile.imsave(sheared_imname, sheared_img)

#%% Reslice
def reslice_img(imname, imname_reslice):
    im = io.imread(imname)
    im_reslice = np.swapaxes(im, 0, 1)  # swaps axis 0 and 2
    tifffile.imsave(imname_reslice, im_reslice)

def save_resliceimg(process_img_list, folder_reslice):
    for i in range(len(process_img_list)):
        for j in range(len(process_img_list[i])):
            imname = process_img_list[i][j]
            imname_reslice = folder_reslice[i] + '\\'+imname.split('\\')[-1]
            reslice_img(imname, imname_reslice)
   
#%%  3D Projection
def projection(cropped, IJ, path_fiji, process_img_list, folder_projection, direction, angle_increment, im_range):
    if cropped == True:
        IJ = imagej.init(path_fiji, mode = 'interactive')
        IJ.ui().showUI()

    for i in range(len(process_img_list)):
        for j in range(len(process_img_list[i])):
            image = process_img_list[i][j]
            image = image.replace("\\", "/")
            #imp =  IJ.io().open(image)         
            process_img = folder_projection[i] + "/" + direction + "-" + str(j) + ".tif"
            process_img = process_img.replace("\\", "/")
            macro_new2 = '''
            open("''' +image + '''");
            selectWindow("''' +image.split('/')[-1] +'''");
            setMinAndMax'''+im_range[i]+''';
            run("3D Project...", "projection=[Brightest Point] axis=''' + direction + '''-Axis slice=1 initial=0 total=360 rotation=''' + str(angle_increment) +''' lower=1 upper=255 opacity=0 surface=100 interior=50 interpolate");
            saveAs("Tiff", "'''+ process_img + '''");
            close();
            close();
            '''
            IJ.py.run_macro(macro_new2)

    return 0

#%%


if __name__ == '__main__':
    print("Rename all files..")
    if renamed == False:
        rename(path)
        pos_name_list = get_img_list(path)
    else:
        pos_name_list = get_img_list(path)
    print("Cropping in Fiji..")
    ana_pos_list, folder_sb, folder_process, folder_reslice, folder_projection = create_analysis_folder(path)
    if cropped == False:
        IJ = crop_images(path_fiji, pos_name_list, roi, sub_background, folder_sb)
    else:
        IJ = fake_ij(cropped)
    print("Deskewing..")
    SB_img_list = read_cropped_images(ana_pos_list)
    img_list = get_all_sb_image_list(SB_img_list)
    
    if processed == False:
        partial_func = partial(process_img, param = param)
        if __name__ == '__main__':
            pool = Pool(num_workers)
            pool.map(partial_func, img_list)
            pool.close()
            pool.join()
    else:
        pass
    
    print("3D projection")   
    
    process_img_list, process_img_last_list = read_processed_images(ana_pos_list)
    im_range = get_intensity_range(process_img_list, process_img_last_list)
    
    if resliced == False:
        save_resliceimg(process_img_list, folder_reslice)
    
    if projectioned == False:
        projection(cropped, IJ, path_fiji, process_img_list, folder_projection, direction, angle_increment, im_range)

    print("Combining frame")
    projection_img_list = read_projection_images(ana_pos_list)

    print("Completed")
    
    
        
        
    
    
