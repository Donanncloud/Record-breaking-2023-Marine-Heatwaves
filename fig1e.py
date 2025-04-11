# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 15:33:05 2024

@author: dty
"""

import netCDF4 as nc
import xarray as xr
import numpy as np
import os 
import pandas as pd
import matplotlib.pyplot as plt 

from scipy.stats import genextreme as gev

#%%
file_path = 'C:/code-py/globalmhw/gev_v2/bar/'
file_name1 = os.path.join(file_path,'cum-ann-mhw-time-90_northatlantic.nc')
file_name2 = os.path.join(file_path,'cum-ann-mhw-time-90_tropicalpacific.nc')
file_name3 = os.path.join(file_path,'cum-ann-mhw-time-90_northpacific.nc')
file_name4 = os.path.join(file_path,'cum-ann-mhw-time-90_southwestpacific.nc')


year = np.arange(1982,2024,1)

ds = xr.open_dataset(file_name1)
cum_1 = np.array(ds.cum[:])
cum_1 = cum_1-cum_1.mean()
ds = xr.open_dataset(file_name2)
cum_2 = np.array(ds.cum[:])

ds = xr.open_dataset(file_name3)
cum_3 = np.array(ds.cum[:])

ds = xr.open_dataset(file_name4)
cum_4 = np.array(ds.cum[:])



#%% 
#1 area=================================================   
data = cum_1[0:42] 
shape, loc, scale = gev.fit(data)   
# 假设 cum_1[42] 是您感兴趣的观测值
observed_value = data[41]
# 使用 gev.ppf 计算对应的累积概率
probability = gev.cdf(observed_value, shape, loc, scale)
# 计算与该累积概率对应的回归周期
return_period_1 = 1 / (1-probability)

#2 area=================================================   
data = cum_2[0:42] 
shape, loc, scale = gev.fit(data)   
# 假设 cum_1[42] 是您感兴趣的观测值
observed_value = data[41]
# 使用 gev.ppf 计算对应的累积概率
probability = gev.cdf(observed_value, shape, loc, scale)
# 计算与该累积概率对应的回归周期
return_period_2 = 1 / (1-probability)
#3 area=================================================   
data = cum_3[0:42] 
shape, loc, scale = gev.fit(data)   
# 假设 cum_1[42] 是您感兴趣的观测值
observed_value = data[41]
# 使用 gev.ppf 计算对应的累积概率
probability = gev.cdf(observed_value, shape, loc, scale)
# 计算与该累积概率对应的回归周期
return_period_3 = 1 / (1-probability)
#4 area=================================================   
data = cum_4[0:42] 
shape, loc, scale = gev.fit(data)   
# 假设 cum_1[42] 是您感兴趣的观测值
observed_value = data[41]
# 使用 gev.ppf 计算对应的累积概率
probability = gev.cdf(observed_value, shape, loc, scale)
# 计算与该累积概率对应的回归周期
return_period_4 = 1 / (1-probability)


# Create a data frame
df = pd.DataFrame ({
        'Group':  ['North Atlantic', 'Tropical Eastern Pacific', 'North Pacific', \
                   'Southwest Pacific'],
        'Value': [return_period_1,return_period_2,return_period_3,return_period_4]
})

df.head()
df = df.sort_values(by=['Value'])


#%%



max_value = df['Value'].max()

new_lefts = [max_value - value for value in df['Value']]

# 创建figure和axes对象
fig, ax = plt.subplots(figsize=(6,8))

# 绘制水平条形图，从y轴右侧开始
colors = [ 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey'] 
ax.barh(y=df['Group'], width=df['Value'],height=0.5, left=new_lefts, color=colors)

# 设置标题和轴标签
# ax.set_title('Return Levels')
# ax.set_xlabel('Values')  # 这里实际上是x轴的标签，对应的是数值
# ax.set_ylabel('Categories')

# 将y轴刻度和标签移动到右侧
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

# 设置x轴范围
ax.set_xlim([0, 276])
# 假设 x 轴刻度间隔
xticks =[  1,  76, 126, 176, 226, 276]

# 自定义的 x 轴刻度标签
aa = [275,  200, 150, 100, 50, 0]
 # 确保 xticks 和 aa 长度一致，若不一致则需要调整其中任一数组的长度或步长
assert len(xticks) == len(aa), "Length of xticks and labels should match."

# 将自定义标签赋给 xlabels
xlabels = aa
# 设置 x 轴刻度及对应的标签
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels,size=20,fontname="Arial")
ax.set_yticklabels(df['Group'],size=20,fontname="Arial")
ax.axvline(176, 0.0, 0.95,color='black')

ax.set_xlabel('Return Levels (year)',size=20,fontname="Arial")

# 关闭y轴左侧轴线和刻度
ax.spines['left'].set_visible(False)  # 隐藏左侧边框
ax.yaxis.tick_left=False  # 取消左侧刻度线

# 关闭x轴上侧轴线和刻度
ax.spines['top'].set_visible(False)  # 隐藏上侧边框
ax.xaxis.tick_top=False  # 取消上侧刻度线
# 显示图形
plt.savefig( 'Figure1e.eps', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.show()

