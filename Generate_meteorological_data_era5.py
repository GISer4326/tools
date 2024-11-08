import pandas as pd
import xarray as xr
import numpy as np
import tkinter as tk
from tkinter import filedialog
from collections import namedtuple
import glob
import csv
from datetime import datetime  
import sys
import math

# 函数提取月最高温度以及最低温度
def add_monthly_extreme_temperatures(data_array):
    # 将数组转换为Pandas DataFrame
    df = pd.DataFrame(data_array, columns=['StationID', 'DATETIME', 'TMEAN'])

    # 将DATETIME列转换为日期时间格式
    df['DATETIME'] = pd.to_datetime(df['DATETIME'])

    # 提取年份和月份作为新列
    df['Year'] = df['DATETIME'].dt.year
    df['Month'] = df['DATETIME'].dt.month

    # 计算每个站点每个月的最大和最小平均温度
    monthly_extremes = df.groupby(['StationID', 'Year', 'Month'])['TMEAN'].agg(['max', 'min']).reset_index()

    # 将最大和最小平均温度合并到原始数组
    merged_data = pd.merge(df, monthly_extremes, on=['StationID', 'Year', 'Month'], how='left').drop(
        columns=['Year', 'Month'])

    # 将最大和最小平均温度添加到原始数组的每个记录
    data_array_with_extremes = merged_data.values

    return data_array_with_extremes
# 函数用于获取站点信息----------------------------------------------------------------------------
def find_valid_coordinates(csv_file_path):
    valid_coordinates = []  # 初始化一个空列表，用于存储有效的经纬度坐标
    with open(csv_file_path, mode='r') as file:
        reader = csv.DictReader(file, delimiter=',')  # 指定CSV文件的分隔符为制表符
        for row in reader:
            try:
                latitude = float(row['Lat'])  # 从CSV文件中读取纬度信息
                longitude = float(row['Lon'])  # 从CSV文件中读取经度信息
                valid_coordinates.append((latitude, longitude))  # 将有效的经纬度坐标添加到列表中
            except ValueError:
                pass  # 处理可能出现的值错误
    return valid_coordinates

# 函数用于提取每天每个站点的气象数据-------------------------------------------------------------------------
def ExtractvaribleBycoordinates(ds,Sites_file,new_array,Variablekey):
    # 提取流域范围内的格网点经纬度
    latitudes = ds['lat'][:]
    longitudes = ds['lon'][:]
    variable_values = ds[Variablekey].values
    time = ds['time'].values
    coordinates = find_valid_coordinates(Sites_file)
    grid_points = []
    for coordinate in coordinates:
        latitude, longitude = coordinate
        lat_index = np.where(latitudes == latitude)[0][0]
        lon_index = np.where(longitudes == longitude)[0][0]
        for t, current_time in enumerate(time):
            current_value = variable_values[t, lat_index, lon_index]
            if not np.isnan(current_value):
                grid_point = (current_time, latitude, longitude, current_value)
                grid_points.append(grid_point)
    # 定义命名元组
    Point = namedtuple('Station', ['latitude', 'longitude'])
    # 创建一个空字典
    points_dict = {}
    # 使用循环将经纬度信息映射为字典的键
    for index, coordinate in enumerate(coordinates, start=0):
        # 创建命名元组对象
        point = Point(coordinate[0], coordinate[1])  # 直接使用元组的值
        # 使用键名0,1,2...等来存储点的信息
        key = index
        points_dict[key] = point
    for point in grid_points:
        # 从点的经纬度信息中创建命名元组
        point_tuple = Point(point[1], point[2])
        # 遍历字典，查找匹配的ID名
        for key, value in points_dict.items():
            if value == point_tuple:
                # 如果点的经纬度信息与字典中的某个值匹配，则将该ID名、current_time和current_value添加到新数组中
                new_array.append((key, point[0], point[3]))
                break
    return new_array
    
def ExtractvaribleBycoordinates2(ds,Sites_file,new_array,Variablekey):
    # 提取流域范围内的格网点经纬度
    latitudes = ds['lat'][:]
    longitudes = ds['lon'][:]
    variable_values = ds[Variablekey].values
    time = ds['time'].values
    coordinates = find_valid_coordinates(Sites_file)
    grid_points = []
    for coordinate in coordinates:
        latitude, longitude = coordinate
        lat_index = np.where(latitudes == latitude)[0][0]
        lon_index = np.where(longitudes == longitude)[0][0]
        for t, current_time in enumerate(time):
            current_value = variable_values[t, lat_index, lon_index]
            if not np.isnan(current_value):
                grid_point = (current_time, latitude, longitude, current_value)
                grid_points.append(grid_point)
    # 定义命名元组
    Point = namedtuple('Station', ['latitude', 'longitude'])
    # 创建一个空字典
    points_dict = {}
    # 使用循环将经纬度信息映射为字典的键
    for index, coordinate in enumerate(coordinates, start=0):
        # 创建命名元组对象
        point = Point(coordinate[0], coordinate[1])  # 直接使用元组的值
        # 使用键名0,1,2...等来存储点的信息
        key = index
        points_dict[key] = point
    for point in grid_points:
        # 从点的经纬度信息中创建命名元组
        point_tuple = Point(point[1], point[2])
        # 遍历字典，查找匹配的ID名
        for key, value in points_dict.items():
            if value == point_tuple:
                # 如果点的经纬度信息与字典中的某个值匹配，则将该ID名、current_time和current_value添加到新数组中
                # new_array[key] = [point[0], point[3]]
                new_array.setdefault(key,new_array.setdefault(key,[]).append([point[0], point[3]]))
                break
    return new_array    

# 合并同一类型的nc数据根据时间序列--------------------------------------------------------------------------
def MergeNcData(CMFDFiles):
    all_files = glob.glob(CMFDFiles + '\*.nc')
    all_files.sort()
    print(all_files)
    data = None
    for i in range(len(all_files)):
        file = xr.open_dataset(all_files[i])
        if data is None:
            data = file
        else:
            data = xr.concat([data, file], dim='time')
    return data

def MergeNcData99_08(CMFDFiles):
    all_files = glob.glob(CMFDFiles + '\*.nc')
    all_files.sort()
    print(all_files)
    print(len(all_files))
    data = None
    for i in range(int(len(all_files)/2)):
        print(i)
        file = xr.open_dataset(all_files[i])
        if data is None:
            data = file
        else:
            data = xr.concat([data, file], dim='time')
    return data

def MergeNcData08_18(CMFDFiles):
    all_files = glob.glob(CMFDFiles + '\*.nc')
    all_files.sort()
    print(all_files)
    
    data = None
    
    for i in range(int(len(all_files)/2)):
        file = xr.open_dataset(all_files[i + int(len(all_files)/2)])
        if data is None:
            data = file
        else:
            data = xr.concat([data, file], dim='time')
    return data


if __name__ == "__main__":

    project_path = sys.argv[1]
    climate_path = project_path + r"\data_prepare\climate"
    # 选择文件输入输出路径
    Sites_file = climate_path + r"\Sites_M.csv" # 请替换为您的格网点CSV文件路径

    # 1.生成平均温度以及当月最高温度时间站点信息
    # 生成一个空数组用来存放所有平均温度气象数据
    # TMEAN_array_all = []
    TMEAN_array = []
    # TMEAN_CMFDFiles = r"H:\huanjing\CMFD\temp"
    TMEAN_CMFDFiles = r"C:\Data\CMFD\1\temp"
    
    TMEAN_data = MergeNcData(TMEAN_CMFDFiles)
    # print(TMEAN_data)
    # 生成每日平均温度
    TMEAN_array =ExtractvaribleBycoordinates(TMEAN_data,Sites_file,TMEAN_array,'temp')
    TMEAN_data.close()
 
    
    # 处理：提取出每月最高最低温度
    TMEAN_array = add_monthly_extreme_temperatures(TMEAN_array)
    print('1.平均温度气象数据读取完成！')

    # 2.生成最高温度气象数据
    # 生成一个空数组用来存放所有最高气温气象数据
    # TMAX_array_all = []
    TMAX_array = []
    # TMAX_CMFDFiles = r"H:\huanjing\CMFD\tmax"
    TMAX_CMFDFiles = r"C:\Data\CMFD\1\tmax"
    
    TMAX_data = MergeNcData(TMAX_CMFDFiles)
    # 生成每日最高温度
    TMAX_array =ExtractvaribleBycoordinates(TMAX_data,Sites_file,TMAX_array,'temp')
    TMAX_data.close()

    print('2.最高温度气象数据读取完成！')

    # 3.生成最低温度气象数据
    # 生成一个空数组用来存放所有最低气温气象数据
    # TMIN_array_all = []
    TMIN_array = []
    # TMIN_CMFDFiles = r"H:\huanjing\CMFD\tmin"
    TMIN_CMFDFiles = r"C:\Data\CMFD\1\tmin"
    
    TMIN_data = MergeNcData(TMIN_CMFDFiles)
    # 生成每日最低温度
    TMIN_array =ExtractvaribleBycoordinates(TMIN_data,Sites_file,TMIN_array,'temp')
    TMIN_data.close()

    print('3.最低温度气象数据读取完成！')

    # # 4.生成相对湿度气象数据
    # # 生成一个空数组用来存放所有湿度气象数据
    # # RM_array_all = []
    # RM_array = []
    # # RM_CMFDFiles = r"H:\huanjing\CMFD\shum"
    # RM_CMFDFiles = r"C:\Data\CMFD\1\shum"
    
    # RM_data = MergeNcData(RM_CMFDFiles)
    # # 生成湿度数据
    # RM_array = ExtractvaribleBycoordinates(RM_data,Sites_file,RM_array,'shum')
    # RM_data.close()

    # print('4.湿度气象数据读取完成！')
    
    # 4.生成露点温度气象数据
    # RM_array_all = []
    TDAW_array = []
    # RM_CMFDFiles = r"H:\huanjing\CMFD\shum"
    TDAW_CMFDFiles = r"C:\Data\CMFD\1\tdaw"
    
    TDAW_data = MergeNcData(TDAW_CMFDFiles)
    # 生成湿度数据
    TDAW_array = ExtractvaribleBycoordinates(TDAW_data,Sites_file,TDAW_array,'temp')
    TDAW_data.close()

    print('4.露点温度气象数据读取完成！')

    # 5.生成风速气象数据
    # 生成一个空数组用来存放所有湿度气象数据
    # WS_array_all = []
    WSx_array = []
    # WS_CMFDFiles = r"H:\huanjing\CMFD\wind"
    WSx_CMFDFiles = r"C:\Data\CMFD\1\winx"
    
    WSx_data = MergeNcData(WSx_CMFDFiles)
    # 生成风速气象数据
    WSx_array = ExtractvaribleBycoordinates(WSx_data,Sites_file,WSx_array,'wind')
    WSx_data.close()

    print('5.1风速气象数据读取完成！')
    
    # 5.生成风速气象数据
    # 生成一个空数组用来存放所有湿度气象数据
    # WS_array_all = []
    WSy_array = []
    # WS_CMFDFiles = r"H:\huanjing\CMFD\wind"
    WSy_CMFDFiles = r"C:\Data\CMFD\1\winy"
    
    WSy_data = MergeNcData(WSy_CMFDFiles)
    # 生成风速气象数据
    WSy_array = ExtractvaribleBycoordinates(WSy_data,Sites_file,WSy_array,'wind')
    WSy_data.close()

    print('5.2风速气象数据读取完成！')

    # 6.生成太阳辐射气象数据
    # 生成一个空数组用来存放所有太阳辐射气象数据
    # SR_array_all = []
    SR_array = []
    # SR_CMFDFiles = r"H:\huanjing\CMFD\srad"
    SR_CMFDFiles = r"C:\Data\CMFD\1\srad"
    
    SR_data = MergeNcData(SR_CMFDFiles)
    # 生成太阳辐射气象数据
    SR_array = ExtractvaribleBycoordinates(SR_data,Sites_file,SR_array,'srad')
    SR_data.close()

    print('6.太阳辐射气象数据读取完成！')
    
    # 7.生成气压气象数据
    # 生成一个空数组用来存放所有气压气象数据
    PRES_array = []
    # PRES_array_all = []
    # PRES_CMFDFiles = r"H:\huanjing\CMFD\pres"
    PRES_CMFDFiles = r"C:\Data\CMFD\1\pres"
    
    PRES_data = MergeNcData(PRES_CMFDFiles)
    # 生成气压气象数据
    PRES_array = ExtractvaribleBycoordinates(PRES_data,Sites_file,PRES_array,'pres')
    PRES_data.close()

    print('7.气压气象数据读取完成！')
    
    # 8.生成降雨气象数据
    # 生成一个空数组用来存放所有降雨气象数据
    # PREC_array_all = []
    PREC_array = {}
    # PREC_CMFDFiles = r"H:\huanjing\CMFD\prec"
    PREC_CMFDFiles = r"C:\Data\CMFD\1\prec"
    
    PREC_data = MergeNcData(PREC_CMFDFiles)
    # 生成降雨气象数据
    PREC_array = ExtractvaribleBycoordinates2(PREC_data,Sites_file,PREC_array,'prec')
    PREC_data.close()

    print('8.降雨气象数据读取完成！')
    
    # 9.定义输出文件CSV文件路径
    M_Outputfile_path = climate_path + r"\meteo_daily_CMFD.csv"
    P_Outputfile_path = climate_path + r"\pcp_daily_CMFD.csv"
    # #写入CSV文件
    # with open(P_Outputfile_path, mode='w', newline='') as file:
        # writer = csv.writer(file)
        # # 写入表头
        # #UTCTIME
        # writer.writerow(["#UTCTIME"])
        # writer.writerow(['DATETIME', 0])
        
        # # 写入数据
        # for TMEAN, PREC in zip(TMEAN_array , PREC_array):
            # # print(TMEAN[1])
            # # 将字符串转换为datetime对象  
            # date_obj = datetime.strptime(str(TMEAN[1]), "%Y-%m-%d %H:%M:%S")  
            # # 将时间部分设置为00:00:00  
            # new_date_obj = date_obj.replace(hour=0, minute=0, second=0, microsecond=0)
            # writer.writerow([str(new_date_obj), PREC[2] * 24])

    # print("降雨数据文件已生成:", P_Outputfile_path)
    
    #写入CSV文件
    with open(P_Outputfile_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        # 写入表头
        #UTCTIME
        writer.writerow(["#UTCTIME"])
        
        num = len(PREC_array.keys())
        headarray = []
        headarray.append('DATETIME')
        for i in range(num):
            headarray.append(str(i))
        
        writer.writerow(headarray)
        
        # # 写入数据
        # for TMEAN, PREC in zip(TMEAN_array , PREC_array):
            # # print(TMEAN[1])
            # # 将字符串转换为datetime对象  
            # date_obj = datetime.strptime(str(TMEAN[1]), "%Y-%m-%d %H:%M:%S")  
            # # 将时间部分设置为00:00:00  
            # new_date_obj = date_obj.replace(hour=0, minute=0, second=0, microsecond=0)
            # writer.writerow([str(new_date_obj), PREC[2] * 24])
        
        timenum = len(PREC_array[0])
        
        for index in range(timenum):
            row = []
            for id in range(num):
                if id == 0:
                    print("time", str(PREC_array[id][index][0]))
                    time = str(PREC_array[id][index][0]).replace("T"," ").split('.')[0]
                    print("time", time)
                    # 将字符串转换为datetime对象  
                    date_obj = datetime.strptime(time, "%Y-%m-%d %H:%M:%S")  
                    # 将时间部分设置为00:00:00  
                    new_date_obj = date_obj.replace(hour=0, minute=0, second=0, microsecond=0)
                    row.append(str(new_date_obj))
                row.append(PREC_array[id][index][1] * 1000)   
            writer.writerow(row)

    print("降雨数据文件已生成:", P_Outputfile_path)
    
    
    #写入CSV文件
    with open(M_Outputfile_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        # 写入表头
        writer.writerow(["#UTCTIME"])
        writer.writerow(['StationID', 'DATETIME', 'TMEAN','TMAX','TMIN','RM','WS','SR','MAXMONT','MINMONT'])

        # 写入数据
        for TMEAN , TMAX , TMIN , TDAW , Wx , Wy , SR , PRES in zip(TMEAN_array, TMAX_array, TMIN_array, TDAW_array, WSx_array, WSy_array, SR_array, PRES_array):
            # print(TMEAN[1])
            # 将字符串转换为datetime对象  
            date_obj = datetime.strptime(str(TMEAN[1]), "%Y-%m-%d %H:%M:%S")  
            # 将时间部分设置为00:00:00  
            new_date_obj = date_obj.replace(hour=0, minute=0, second=0, microsecond=0)
            
            # # CMFD
            # temp = round(TMEAN[2]-273.15,2)
            # sr = SR[2] * 24 * 3600 / 1000000
            # shum = RM[2]
            # pres = PRES[2]
            # rm = 0.236 * pres * shum * np.exp((17.67 * (TMEAN[2] - 273.16))/(TMEAN[2] - 29.65)) ** (-1)
            # writer.writerow([TMEAN[0],str(new_date_obj), temp, round(TMAX[2]-273.15,2),round(TMIN[2]-273.15,2), rm ,WS[2], sr ,round(TMEAN[3]-273.15,2),round(TMEAN[4]-273.15,2)])
            
            # ERA5
            # 平均气温 摄氏度
            temp = round(TMEAN[2], 2)
            # 太阳辐射  J/m^2/天 -- MJ/m^2/天
            sr = SR[2] / 1000000
            # 相对湿度  %
            e = 6.11 * ( 10 ** ( ( 7.45 * round(TDAW[2], 2) ) / ( 235 + round(TDAW[2], 2) ) ) )
            E = 6.11 * ( 10 ** ( ( 7.45 * temp ) / ( 235 + temp ) ) )
            rm = ( e / E ) * 100
            # 气压 暂时没用
            pres = PRES[2]
            # 风速 m/s
            wind = math.sqrt(round(Wx[2], 2) ** 2 + round(Wy[2], 2) ** 2)
            
            #                  id           date         平均气温   最高气温        最低气温    相对湿度  风速  太阳辐射  月最高气温     月最低气温
            writer.writerow([TMEAN[0], str(new_date_obj), temp, round(TMAX[2],2), round(TMIN[2],2) , rm , wind , sr , round(TMEAN[3],2), round(TMEAN[4],2)])

        # # 写入数据
        # for TMEAN in zip(TMEAN_array):
            # # print(TMEAN[1])
            # # 将字符串转换为datetime对象  
            # date_obj = datetime.strptime(str(TMEAN[1]), "%Y-%m-%d %H:%M:%S")  
            # # 将时间部分设置为00:00:00  
            # new_date_obj = date_obj.replace(hour=0, minute=0, second=0, microsecond=0)
            # temp = round(TMEAN[2]-273.15,2)
            
            # rm = 0.236 * pres * shum * np.exp((17.67 * (TMEAN[2] - 273.16))/(TMEAN[2] - 29.65)) ** (-1)
            # writer.writerow([TMEAN[0],str(new_date_obj), temp])

    print("气象数据文件已生成:", M_Outputfile_path)








