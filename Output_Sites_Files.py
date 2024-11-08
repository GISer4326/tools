import xarray as xr
import numpy as np
from osgeo import osr, ogr
import pandas as pd
from osgeo import gdal
import geopandas as gpd
import sys
import os
import math

os.environ['PROJ_LIB'] = r'D:\Anaconda3\envs\python36_2\Lib\site-packages\osgeo\data\proj'

# 提取高程函数-------------------------------------------------------------------------------------
def get_elevation_from_tif(longitude, latitude, dem_file):
    dataset = gdal.Open(dem_file)
    geo_transform = dataset.GetGeoTransform()
    x = int((longitude - geo_transform[0]) / geo_transform[1])
    y = int((latitude - geo_transform[3]) / geo_transform[5])
    band = dataset.GetRasterBand(1)
    
    if band.ReadAsArray(x, y, 1, 1):
        elevation = band.ReadAsArray(x, y, 1, 1)[0][0]
    else:
        elevation = -9999
    return elevation

# 提取坐标函数-------------------------------------------------------------------------------------
def find_valid_coordinates_with_values(latitudes, longitudes, variable_values):
    print(len(variable_values))
    valid_coordinates_with_values = []  # 初始化一个空列表，用于存储有效的经纬度坐标及对应的变量值
    for lat_index, latitude in enumerate(latitudes):
        print("lat_index",lat_index)
        for lon_index, longitude in enumerate(longitudes):
            try:
                # value = variable_values[lat_index, lon_index]  # 获取当前经纬度坐标对应的变量值
                value = variable_values[lon_index, lat_index]  # 获取当前经纬度坐标对应的变量值
                print("value",value)
                if not np.isnan(value).any():  # 检查变量值是否为NaN
                    valid_coordinates_with_values.append(((latitude.values.item(), longitude.values.item()), value))
                    # 将有效的经纬度坐标及对应的变量值添加到列表中
            except IndexError:
                pass  # 处理可能出现的索引错误
    return valid_coordinates_with_values

def create_rectangle_shp(min_lon, min_lat, max_lon, max_lat, output_shp):  
    # 创建一个新的Shapefile  
    driver = ogr.GetDriverByName('ESRI Shapefile')  
    if os.path.exists(output_shp):  
        driver.DeleteDataSource(output_shp)  
    data_source = driver.CreateDataSource(output_shp)  
      

    # 定义空间参考系统为WGS84  
    srs = osr.SpatialReference()  
    srs.ImportFromEPSG(4326)  
    # 创建一个新的图层  
    out_layer = data_source.CreateLayer('polygon_layer', geom_type=ogr.wkbPolygon, srs=srs)    
      
    # 添加字段  
    field_name = 'id'  
    field_defn = ogr.FieldDefn(field_name, ogr.OFTString)  
    field_defn.SetWidth(255)  
    if out_layer.CreateField(field_defn) != 0:  
        print(f"Failed to create field {field_name}")    
          
    # 创建一个矩形面  
    ring = ogr.Geometry(ogr.wkbLinearRing)  
    ring.AddPoint(min_lon, min_lat)  
    ring.AddPoint(max_lon, min_lat)  
    ring.AddPoint(max_lon, max_lat)  
    ring.AddPoint(min_lon, max_lat)  
    ring.AddPoint(min_lon, min_lat)  # 闭合多边形  
    polygon = ogr.Geometry(ogr.wkbPolygon)  
    polygon.AddGeometry(ring)  
      
    # 创建一个新的要素并设置其几何形状  
    feature = ogr.Feature(out_layer.GetLayerDefn())  
    feature.SetGeometry(polygon)  
    
    # 设置字段值  
    feature.SetField(field_name, '0')  
      
    # 将要素添加到图层  
    out_layer.CreateFeature(feature)  
      
    # 清理并关闭数据源  
    feature = None  
    data_source = None  

if __name__ == "__main__":

    project_path = sys.argv[1]
    # 设置输入流域范围shp文件路径
    shp_file_path = project_path + r"/workspace/spatial_shp/basin.shp"
    # 设置输出Sites_M文件路径
    sites_file_path1 = project_path + r"/data_prepare/climate/Sites_M.csv"
    sites_file_path2 = project_path + r"/data_prepare/climate/Sites_P.csv"
    # 设置高程文件路径-------------------------------------------------------------------------------------------------------------------------------
    dem_file = project_path + r"/data_prepare/spatial/basin_elv.tif"
    # 设置nc数据文件路径
    # nc_file_path = r"H:\huanjing\CMFD\tmax\tmax_2016.nc"
    nc_file_path = r"C:\Data\CMFD\1\tmax\tmax_2002.nc"
    # nc_file_path = r"C:\Data\CMFD\1\tmax\tmax2019.nc"

    # 程序开始运行-----------------------------------------------------------------------------------------------------------------------------
    # 读取shp文件
    gdf = gpd.read_file(shp_file_path)
    # 获取裁剪区域的边界框
    bbox = gdf.bounds.iloc[0]
    # 分别获取边界框的最小经度、最小纬度、最大经度、最大纬度
    min_lon, min_lat, max_lon, max_lat = bbox['minx'], bbox['miny'], bbox['maxx'], bbox['maxy']
    
    # # 生成外界矩形作为泰森多边形.在java中生成了。
    # TS_m_path = project_path + r"\data_prepare\spatial\TS_M_WGS84.shp"
    # TS_p_path = project_path + r"\data_prepare\spatial\TS_P_WGS84.shp"
    # create_rectangle_shp(min_lon, min_lat, max_lon, max_lat,TS_m_path)
    # create_rectangle_shp(min_lon, min_lat, max_lon, max_lat,TS_p_path)

    # # 读取NetCDF数据/打开NetCDF文件
    ds = xr.open_dataset(nc_file_path)
    # 根据裁剪区域的边界选择NetCDF数据的子集
    subset_ds = ds.sel(lon=slice(min_lon, max_lon), lat=slice(min_lat, max_lat))  # 选择NetCDF数据集的子集，根据裁剪区域的边界
    # 打印当前nc数据的信息
    # print(subset_ds.keys())

    # 提取流域范围内的格网点经纬度
    latitudes = subset_ds['lat'][:]  # 获取纬度数据
    longitudes = subset_ds['lon'][:]  # 获取经度数据

    if len(latitudes) == 0:
        min_lat = min_lat + 0.44915764
        max_lat = max_lat - 0.44915764
        subset_ds = ds.sel(lon=slice(min_lon, max_lon), lat=slice(min_lat, max_lat))
        # 提取流域范围内的格网点经纬度
        latitudes = subset_ds['lat'][:]  # 获取纬度数据
        longitudes = subset_ds['lon'][:]  # 获取经度数据
    
    if len(longitudes) == 0:
        print("longitudes")
        min_lon = min_lon - 0.44915764
        max_lon = max_lon + 0.44915764
        subset_ds = ds.sel(lon=slice(min_lon, max_lon), lat=slice(min_lat, max_lat)) 
        # 提取流域范围内的格网点经纬度
        latitudes = subset_ds['lat'][:]  # 获取纬度数据
        longitudes = subset_ds['lon'][:]  # 获取经度数据

    variable_values = subset_ds['temp'].values  # 提取温度变量数据

    # 设置间隔范围阈值设置1间隔经度为0.1度（纬度为0.001），2为0.2度以此类推---------------------------------------------------------------
    # 方法一：目前只取一个站点：位于流域中间
    # lat_Difference = max_lat - min_lat
    # lon_Difference = max_lon - min_lon 
    # nlat = math.floor(lat_Difference / 0.01)
    # nlon = math.floor(lon_Difference / 0.1)
    ilat = math.floor(len(latitudes) / 2)
    ilon = math.floor(len(longitudes) / 2)

    # latitudes = latitudes[ilat:ilat + 1]
    # longitudes =longitudes[ilon:ilon + 1]
    # print(latitudes)
    # print(longitudes)
    
    # 方法二：设置间隔范围阈值取多个站点：
    lat_interval = 1
    lon_interval = 1
    latitudes = latitudes[::lat_interval]
    longitudes =longitudes[::lon_interval]
    print(latitudes)
    print(longitudes)

    # 获取每个坐标索引到的变量值
    coordinates_with_values = find_valid_coordinates_with_values(latitudes, longitudes, variable_values)

    # print(len(coordinates_with_values))

    # 创建地理坐标系和投影坐标系的转换对象
    # 开始投影坐标
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)  # 设置WGS 84坐标系
    beijing_gauss_kruger = osr.SpatialReference()
    beijing_gauss_kruger.ImportFromEPSG(2412)  # 设置北京 1954 年的 3 度 Gauss-Kruger 坐标系
    transform = osr.CoordinateTransformation(wgs84, beijing_gauss_kruger)
    # 遍历经纬度坐标并转换为投影坐标
    projected_coordinates = []
    for coord in coordinates_with_values:
        latitude, longitude = coord[0]  # 获取经纬度坐标
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(latitude,longitude)
        point.Transform(transform)  # 转换为投影坐标
        projected_coord = (point.GetX(), point.GetY())
        projected_coordinates.append(((latitude, longitude), projected_coord))

    # 写入数据
    data = {
        'StationID': [i for i in range(len(projected_coordinates))],
        'Name': [f'Station {i+1}' for i in range(len(projected_coordinates))],
        'LocalX': [coord[1][0] for coord in projected_coordinates],
        'LocalY': [coord[1][1] for coord in projected_coordinates],
        'Lon': [coord[0][1] for coord in projected_coordinates],
        'Lat': [coord[0][0] for coord in projected_coordinates],
        'Elevation': [get_elevation_from_tif(coord[0][1], coord[0][0], dem_file) for coord in projected_coordinates]
        }
    df = pd.DataFrame(data)

    df.to_csv(sites_file_path1, index=False)
    df.to_csv(sites_file_path2, index=False)
    print(f"CSV文件已保存至: {sites_file_path1}")
    print(f"CSV文件已保存至: {sites_file_path2}")
