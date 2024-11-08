# import os  
# import rasterio  
# import xarray as xr  
# import numpy as np  
# from datetime import datetime  
  
# # 文件夹路径，包含所有的TIF文件  
# folder_path = r'I:\data\ERA5_Daily\temperature'  
  
# # 读取所有TIF文件  
# tif_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.tif')]  
# tif_files.sort()  # 确保文件按日期排序  
  
# # 读取第一个文件以获取维度和坐标信息  
# with rasterio.open(tif_files[0]) as src: 
    # print(tif_files[0])
    # transform = src.transform  
    # print(transform)
    # crs = src.crs  
    # height, width = src.shape  
    # x_coords = np.linspace(transform.c, transform.c + transform.a * width, width)  
    # y_coords = np.linspace(transform.f + transform.e * height, transform.f, height)  
    # print(x_coords)
    # print(y_coords)
      
# # # 初始化xarray DataArray列表  
# # data_arrays = []  
  
# # # 遍历所有TIF文件  
# # for tif_file in tif_files:  
    # # with rasterio.open(tif_file) as src:  
        # # # 提取日期信息（假设文件名包含日期）  
        # # date_str = os.path.splitext(os.path.basename(tif_file))[0].split('_')[-1] 
        # # print(date_str)
        # # date = datetime.strptime(date_str, '%Y%m%d')  
          
        # # # 读取数据  
        # # data = src.read(1)  # 假设我们只对第一个波段感兴趣  
          
        # # # 创建DataArray  
        # # da = xr.DataArray(data,  
                          # # dims=('y', 'x', 'time'),  
                          # # coords={'y': (('y',), y_coords),  
                                  # # 'x': (('x',), x_coords),  
                                  # # 'time': (('time',), [date])})  
          
        # # # 添加到列表中  
        # # data_arrays.append(da.expand_dims('time'))  
  
# # # 合并DataArray  
# # combined_da = xr.concat(data_arrays, dim='time')  
  
# # # 保存为NetCDF  
# # output_nc_path = 'I:\data\ERA5_Daily\temperature\output.nc'  \

# # 读取栅格数据并合并
# stack = []
# for tif in tif_files:
    # with rasterio.open(tif) as src:
        # stack.append(src.read())
 
# # 创建xarray DataArray
# da = xr.DataArray(
    # data=np.dstack(stack),
    # dims=['y', 'x', 'time'],  # 确保时间维度在最外层
    # coords={'time': [f'{i+1:03d}' for i in range(len(tif_files))]})

# combined_da.to_netcdf(output_nc_path)  
  
# # 验证  
# ds = xr.open_dataset(output_nc_path)  
# print(ds)


import os
import xarray as xr
import numpy as np
from netCDF4 import Dataset
from datetime import datetime 

# 设置要合并的文件名列表，例如： ['tiff_data_2019_01.tif', 'tiff_data_2019_02.tif', 'tiff_data_2019_03.tif']
# 读取所有TIF文件  
folder_path = r'J:\data\ERA5_Daily\india\winy'  
file_list = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.tif')]  
file_list.sort()  # 确保文件按日期排序  
print(file_list)
# file_list = []
# 读取第一个tiff文件来获取地理信息
ds = xr.open_rasterio(file_list[0])
lon = np.linspace(ds.coords['x'][0], ds.coords['x'][-1], len(ds.coords['x']))
lat = np.linspace(ds.coords['y'][0], ds.coords['y'][-1], len(ds.coords['y']))
# 创建nc文件
nc_file = Dataset(r'J:\data\ERA5_Daily\india\winy\winy2020.nc', 'w', format='NETCDF4')
# 创建nc文件的维度信息
nc_file.createDimension('time', len(file_list))
nc_file.createDimension('lat', len(lat))
nc_file.createDimension('lon', len(lon))
# 创建nc文件的变量信息
time_var = nc_file.createVariable('time', 'str', ('time',))
lat_var = nc_file.createVariable('lat', 'f8', ('lat',))
lon_var = nc_file.createVariable('lon', 'f8', ('lon',))
data_var = nc_file.createVariable('wind', 'f4', ('time', 'lat', 'lon'))
# 将地理信息写入nc文件
lat_var[:] = lat
lon_var[:] = lon
# 遍历tiff文件列表，将数据逐步合并到nc文件中
for i in range(len(file_list)):
    ds = xr.open_rasterio(file_list[i])
    data_var[i,:,:] = ds.values.squeeze()

    filename = os.path.basename(file_list[i]) 
    filename, extension = os.path.splitext(filename)
    # time = filename.split('_')[2]
    time = filename[-8:] 
    print(time)
    # time = time + "0101"
    # 转换为datetime对象（假设是YYYYMMDD格式）  
    original_date = datetime.strptime(time, "%Y%m%d")   
    new_datetime = original_date.replace(hour=0, minute=0, second=0)    
    formatted_datetime_str = new_datetime.strftime("%Y-%m-%d %H:%M:%S") 
    print(formatted_datetime_str)
    time_var[i] = formatted_datetime_str  # 使用数字表示时间戳也可以.2002/1/1 10:30:00  2004-12-31 10:30:00
nc_file.close()
